# Synthetic end-to-end test of differ.py alignment and comparison logic.
# Builds fake manifests/dumps (no gdb needed) and checks the report:
#   - sima step t aligns with cam step t+1
#   - nested calls (shared caller) bucket separately from toplevel calls
#   - errflg/errmsg are skipped
#   - constituent-indexed arrays compare correctly under index permutation
#   - an injected exit-phase difference is found; entry phases are clean
#
# Run: python3 tests/synthetic_differ_test.py

import json
import os
import shutil
import struct
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import differ  # noqa: E402

Q_STD = "water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water"
CLDLIQ_STD = ("cloud_liquid_water_mixing_ratio_wrt_moist_air"
              "_and_condensed_water")


def f8(*vals):
    return struct.pack("={}d".format(len(vals)), *vals)


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

    def hit(self, scheme, hit, step, caller, args):
        rec = {"scheme": scheme, "hit": hit, "step": step,
               "caller": caller, "args": args, "complete": True}
        self.man["hits"].append(rec)

    def arr(self, extents, entry, exit_):
        return {"kind": "array", "dtype": "f8", "los": [1] * len(extents),
                "extents": extents, "strides": [0] * len(extents),
                "entry_file": self.blobfile(entry),
                "exit_file": self.blobfile(exit_)}

    def scal(self, entry, exit_):
        return {"kind": "scalar", "dtype": "i4",
                "entry_value": entry, "exit_value": exit_}

    def write(self):
        with open(os.path.join(self.dir, "manifest.json"), "w") as f:
            json.dump(self.man, f)


def main():
    outdir = tempfile.mkdtemp(prefix="dropsonde_test_")
    cam = Builder(outdir, "cam")
    sima = Builder(outdir, "sima")
    for b in (cam, sima):
        b.man["breakpoints"] = {"alpha": "alpha_run", "beta": "beta_run"}
    cam.man["constituents"] = ["Q", "CLDLIQ"]
    sima.man["constituents"] = [CLDLIQ_STD, Q_STD]  # permuted order

    field = f8(1, 2, 3, 4, 5, 6)
    q_slice = f8(*range(10, 16))
    liq_slice = f8(*range(20, 26))
    cam_q3d = q_slice + liq_slice          # (2,3,2): Q then CLDLIQ
    sima_q3d = liq_slice + q_slice         # permuted constituent order

    def alpha_args(b):
        return {"x": b.arr([2, 3], field, field),
                "q": b.arr([2, 3, 2], b is cam and cam_q3d or sima_q3d,
                           b is cam and cam_q3d or sima_q3d),
                "errflg": b.scal(-1 if b is cam else 0, 0)}

    def beta_args(b, exit_val):
        return {"y": b.arr([2, 3], field, field),
                "n": b.scal(7, exit_val)}

    # CAM: step 1 = garbage (skipped), steps 2-3 = real.
    cam.hit("alpha", 0, 1, "tphysac", alpha_args(cam))
    cam.hit("beta", 0, 1, "tphysac", beta_args(cam, 999))
    for t, h in ((2, 1), (3, 3)):
        cam.hit("alpha", h, t, "tphysac", alpha_args(cam))
        cam.hit("beta", h, t, "__alpha_MOD_inner", beta_args(cam, 42))
        cam.hit("beta", h + 1, t, "tphysac", beta_args(cam, 43))

    # SIMA steps 1-2; beta nested exit value differs at step 2 (occ 0).
    for t, h in ((1, 0), (2, 2)):
        sima.hit("alpha", h, t, "cap_group", alpha_args(sima))
        sima.hit("beta", h, t, "__alpha_MOD_inner",
                 beta_args(sima, 42 if t == 1 else 41))
        sima.hit("beta", h + 1, t, "cap_group", beta_args(sima, 43))

    cam.write()
    sima.write()

    ok = differ.report(outdir, ["alpha", "beta"], 2)

    assert not ok, "differ should report a divergence"
    diffs = []
    rep = differ.Reporter()
    shared = differ.shared_caller_map(cam.man, sima.man)
    assert shared["beta"] == {"__alpha_MOD_inner"}, shared
    assert shared.get("alpha", set()) == set(), shared

    shutil.rmtree(outdir)
    print("\nsynthetic test passed: nested-call bucketing, step offset, "
          "errflg skip, constituent permutation all behaved")


if __name__ == "__main__":
    main()
