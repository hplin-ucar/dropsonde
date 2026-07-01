# Synthetic test of differ.py's portable MAM constituent realignment.
# Builds fake manifests for a routine called from two different constituent
# index conventions (CAM: gas_pcnst vmr, loffset>0; SIMA: packed array,
# loffset=0) and checks that differ.py:
#   - realigns the interstitial constituent axis per species via each model's
#     own loffset and finds an injected per-species difference with the right
#     cam/sima index labels
#   - leaves matching species (and out-of-range species like Q) alone
#   - treats index-space/convention args (loffset, num_q, lmassptr_amode) as
#     notes, not divergences
#   - compares mode/species constant tables over valid slots only, ignoring
#     unused padding that differs by construction
#   - defers cloud-borne (qqcw) args with a note
#
# Run: python3 tests/portable_differ_test.py

import io
import json
import os
import shutil
import struct
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import differ  # noqa: E402


def f8(*v):
    return struct.pack("<{}d".format(len(v)), *v)


def i4(*v):
    return struct.pack("<{}i".format(len(v)), *v)


class Builder(object):
    def __init__(self, outdir, role):
        self.dir = os.path.join(outdir, role)
        os.makedirs(self.dir)
        self.man = {"role": role, "breakpoints": {}, "constituents": None,
                    "hits": [], "notes": []}
        self.idx = 0

    def blobfile(self, data):
        self.idx += 1
        fn = "{:05d}.bin".format(self.idx)
        with open(os.path.join(self.dir, fn), "wb") as f:
            f.write(data)
        return fn

    def arr(self, extents, data, dtype="f8"):
        return {"kind": "array", "dtype": dtype, "los": [1] * len(extents),
                "extents": extents, "strides": [0] * len(extents),
                "entry_file": self.blobfile(data),
                "exit_file": self.blobfile(data)}

    def scal(self, val):
        return {"kind": "scalar", "dtype": "i4",
                "entry_value": val, "exit_value": val}

    def hit(self, scheme, hit, step, caller, args):
        self.man["hits"].append({"scheme": scheme, "hit": hit, "step": step,
                                 "caller": caller, "args": args,
                                 "complete": True})

    def write(self):
        with open(os.path.join(self.dir, "manifest.json"), "w") as f:
            json.dump(self.man, f)


SCHEME = "modal_aero_rename_ccpp"

# One aerosol mode with one mass species (so4) plus a number, each with an
# interstitial (*_a*) and a cloud-borne (*_c*) form. Pointer arrays hold
# 1-based pcnst-space indices; the routine reads array(:,:, ptr - loffset).
#
#   CAM  : vmr/vmrcw are gas_pcnst=3 wide, loffset=1. Cloud-borne lives in a
#          separate vmrcw at the SAME index as its interstitial partner and is
#          not a registered constituent.  full 1-based: Q,so4_a1,num_a1,H2SO4
#   SIMA : one packed array, loffset=0, all species (incl. *_c*) registered.
#          full 1-based: H2SO4,so4_a1,num_a1,Q,so4_c1,num_c1
CAM_CNST = ["Q", "so4_a1", "num_a1", "H2SO4"]
SIMA_CNST = ["H2SO4", "so4_a1", "num_a1", "Q", "so4_c1", "num_c1"]


def args_for(b, is_cam):
    if is_cam:
        # q local (full-1): 0=so4_a1, 1=num_a1, 2=H2SO4  (Q at full 1 excluded)
        q = f8(1.0, 1.0) + f8(2.0, 2.0) + f8(3.0, 3.0)
        # qqcw local: 0=so4_c1, 1=num_c1, 2=unused
        qqcw = f8(10.0, 10.0) + f8(20.0, 20.0) + f8(0.0, 0.0)
        specmw = f8(115.0, 999.0)                 # slot 2 = unused padding
        lmassptr = i4(2, -999999)                 # so4_a1 @ full 2; slot 2 pad
        lmassptrcw = i4(2, 0)                      # cloud uses interstitial idx
        numptr = i4(3)                             # num_a1 @ full 3
        numptrcw = i4(3)                           # cloud uses interstitial idx
        loff, nq = 1, 3
    else:
        # q local: 0=H2SO4,1=so4_a1,2=num_a1,3=Q,4=so4_c1,5=num_c1
        q = (f8(3.0, 3.0) + f8(1.0, 99.0) + f8(2.0, 2.0) +
             f8(7.0, 7.0) + f8(5.0, 5.0) + f8(6.0, 6.0))     # so4_a1 differs
        qqcw = (f8(0.0, 0.0) + f8(0.0, 0.0) + f8(0.0, 0.0) +
                f8(0.0, 0.0) + f8(10.0, 10.0) + f8(20.0, 77.0))  # num_c1 differs
        specmw = f8(115.0, 888.0)                 # slot 2 = unused padding
        lmassptr = i4(2, -1)
        lmassptrcw = i4(5, -1)                     # so4_c1 @ full 5
        numptr = i4(3)
        numptrcw = i4(6)                           # num_c1 @ full 6
        loff, nq = 0, 6
    nspec = nq
    return {
        "loffset": b.scal(loff),
        "num_q": b.scal(nq),
        "q": b.arr([2, 1, nspec], q),
        "qqcw": b.arr([2, 1, nspec], qqcw),
        "specmw_amode": b.arr([2, 1], specmw),
        "lmassptr_amode": b.arr([2, 1], lmassptr, "i4"),
        "lmassptrcw_amode": b.arr([2, 1], lmassptrcw, "i4"),
        "numptr_amode": b.arr([1], numptr, "i4"),
        "numptrcw_amode": b.arr([1], numptrcw, "i4"),
    }


def main():
    outdir = tempfile.mkdtemp(prefix="dropsonde_port_")
    cam = Builder(outdir, "cam")
    sima = Builder(outdir, "sima")
    for b in (cam, sima):
        b.man["breakpoints"] = {SCHEME: SCHEME + "_run"}
    cam.man["constituents"] = CAM_CNST
    sima.man["constituents"] = SIMA_CNST

    cam.hit(SCHEME, 0, 1, "tphysac", args_for(cam, True))    # step1 = garbage
    cam.hit(SCHEME, 0, 2, "tphysac", args_for(cam, True))    # step2 = real
    sima.hit(SCHEME, 0, 1, "cap", args_for(sima, False))
    cam.write()
    sima.write()

    buf = io.StringIO()
    old = sys.stdout
    sys.stdout = buf
    try:
        differ.report(outdir, [SCHEME], 1)
    finally:
        sys.stdout = old
    out = buf.getvalue()
    shutil.rmtree(outdir)

    def need(s):
        assert s in out, "expected in report:\n  {!r}\n---\n{}".format(s, out)

    def absent(s):
        assert s not in out, "unexpected in report:\n  {!r}\n---\n{}".format(
            s, out)

    # interstitial realignment: only so4_a1 differs, labelled with both indices
    need("(interstitial, per-model local index), 1/3 mapped species differ: "
         "so4_a1")
    need("so4_a1 (cam idx 1, sima idx 2) <-> so4_a1:")
    # cloud-borne realignment via pointer arrays: only num_c1 differs
    need("(cloud-borne, per-model local index), 1/2 mapped species differ: "
         "num_c1")
    need("num_c1 (cam idx 2, sima idx 6) <-> num_c1:")
    # matching / out-of-range species must not show up as differences
    absent("H2SO4 (cam")
    absent(" Q (cam")
    absent("so4_c1 (cam")     # cloud-borne match: no per-species diff line
    # convention args -> a note, never a divergence
    need("index-space/convention args not compared")
    for a in ("lmassptr_amode", "lmassptrcw_amode", "numptr_amode",
              "numptrcw_amode", "loffset", "num_q"):
        need(a)
        absent("arg {}:".format(a))
    # constant table compared over valid slots only (padding slot 2 ignored)
    absent("specmw_amode")

    print("portable test passed: interstitial + cloud-borne per-species "
          "realignment, convention notes, valid-slot metadata")


if __name__ == "__main__":
    main()
