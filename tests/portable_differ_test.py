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


def args_for(b, is_cam):
    # CAM: 3-wide vmr (gas_pcnst), loffset=1, full indices 2,3,4 -> local 1,2,3
    #      = bc_a1, so4_a1, H2SO4 (Q at full idx 1 is out of the vmr range).
    # SIMA: 4-wide packed array, loffset=0, permuted order.
    if is_cam:
        q = f8(1.0, 1.0) + f8(2.0, 2.0) + f8(3.0, 3.0)        # bc, so4, H2SO4
        specmw = f8(10.0, 20.0, 30.0, 111.0)                  # slot 4 = padding
        lmassptr = i4(1, 3, 2, -999999)                       # slot 4 invalid
        qqcw = f8(0.0, 0.0) + f8(0.0, 0.0) + f8(0.0, 0.0)
        loff, nq = 1, 3
        nspec = 3
    else:
        # sima order: H2SO4(0), so4_a1(1), bc_a1(2), Q(3)
        q = f8(3.0, 3.0) + f8(2.0, 99.0) + f8(1.0, 1.0) + f8(7.0, 7.0)
        specmw = f8(10.0, 20.0, 30.0, 222.0)                  # slot 4 = padding
        lmassptr = i4(1, 3, 2, -1)                            # slot 4 invalid
        qqcw = f8(0.0, 0.0) + f8(0.0, 0.0) + f8(0.0, 0.0) + f8(0.0, 0.0)
        loff, nq = 0, 4
        nspec = 4
    return {
        "loffset": b.scal(loff),
        "num_q": b.scal(nq),
        "q": b.arr([2, 1, nspec], q),
        "qqcw": b.arr([2, 1, nspec], qqcw),
        "specmw_amode": b.arr([2, 2], specmw),
        "lmassptr_amode": b.arr([2, 2], lmassptr, "i4"),
    }


def main():
    outdir = tempfile.mkdtemp(prefix="dropsonde_port_")
    cam = Builder(outdir, "cam")
    sima = Builder(outdir, "sima")
    for b in (cam, sima):
        b.man["breakpoints"] = {SCHEME: SCHEME + "_run"}
    cam.man["constituents"] = ["Q", "bc_a1", "so4_a1", "H2SO4"]
    sima.man["constituents"] = ["H2SO4", "so4_a1", "bc_a1", "Q"]

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
    need("1/3 mapped species differ: so4_a1")
    need("so4_a1 (cam idx 2, sima idx 2) <-> so4_a1:")
    # matching / out-of-range species must not show up as differences
    absent("differ: bc_a1")
    absent("H2SO4 (cam")
    absent(" Q (cam")
    # convention args -> a note, never a divergence
    need("index-space/convention args not compared")
    for a in ("lmassptr_amode", "loffset", "num_q"):
        need(a)
        absent("arg {}:".format(a))
    # constant table compared over valid slots only (padding slot 4 ignored)
    absent("specmw_amode")
    # cloud-borne deferred
    need("cloud-borne args not yet compared")
    need("qqcw")

    print("portable test passed: interstitial per-species realignment, "
          "convention notes, valid-slot metadata, cloud-borne deferral")


if __name__ == "__main__":
    main()
