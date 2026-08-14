# Synthetic end-to-end test of same-model (base/test) mode.
# Builds fake manifests/dumps (no gdb needed) plus a suite.json carrying
# roles/base-test, step_offset 0, and a non-default sentinel, and checks:
#   - report labels are base/test everywhere (no cam/sima leakage)
#   - step-offset 0 pairs step t with step t, including in the alignment
#     verdict (which hardcoded t+1 before same-model mode)
#   - identical constituent registries identity-map
#   - an injected exit-phase difference is found
# Also unit-tests the driver's dropsonde-targets-v1 loader.
#
# Run: python3 tests/same_model_differ_test.py

import importlib.machinery
import json
import os
import re
import shutil
import struct
import sys
import tempfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import differ  # noqa: E402


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
        self.man["hits"].append({"scheme": scheme, "hit": hit, "step": step,
                                 "caller": caller, "args": args,
                                 "complete": True})

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


def targets_loader_test():
    drv = importlib.machinery.SourceFileLoader(
        "dropsonde_drv", os.path.join(HERE, "dropsonde")).load_module()
    d = tempfile.mkdtemp(prefix="dropsonde_targets_")

    def spec_file(obj):
        p = os.path.join(d, "t.json")
        with open(p, "w") as f:
            json.dump(obj, f)
        return p

    targets, sentinel, capture = drv.load_targets(spec_file(
        {"spec": "dropsonde-targets-v1", "step_sentinel": "clm_drv",
         "capture": {"element_t": ["state%v", "State%T "]},
         "targets": ["foo", {"sub": "bar", "test": "bar_new"}]}))
    assert targets == [{"sub": "foo", "base": "foo", "test": "foo"},
                       {"sub": "bar", "base": "bar", "test": "bar_new"}], \
        targets
    assert sentinel == "clm_drv", sentinel
    assert capture == {"element_t": ["state%v", "state%t"]}, capture

    # defaults: sentinel cam_run1, no capture
    _t, sentinel, capture = drv.load_targets(spec_file(
        {"spec": "dropsonde-targets-v1", "targets": [{"sub": "foo"}]}))
    assert sentinel == "cam_run1" and capture == {}

    for bad in (
            {"targets": ["foo"]},                              # no spec tag
            {"spec": "dropsonde-targets-v1", "targets": []},   # empty
            {"spec": "dropsonde-targets-v1",
             "targets": ["foo", "foo"]},                       # duplicate
            {"spec": "dropsonde-targets-v1",
             "targets": [{"sub": "foo", "sym": "x"}]},         # unknown key
    ):
        try:
            drv.load_targets(spec_file(bad))
        except SystemExit:
            pass
        else:
            raise AssertionError("accepted bad spec: {}".format(bad))
    shutil.rmtree(d)


def main():
    targets_loader_test()

    outdir = tempfile.mkdtemp(prefix="dropsonde_same_")
    base = Builder(outdir, "base")
    test = Builder(outdir, "test")
    for b in (base, test):
        b.man["breakpoints"] = {"alpha": "alpha", "beta": "beta"}
        b.man["constituents"] = ["Q", "CLDLIQ"]   # identical registries

    field = f8(1, 2, 3, 4, 5, 6)
    field_perturbed = f8(9, 2, 3, 4, 5, 6)   # 1/6 elements: field-specific

    def alpha_args(b, t):
        # x: entry diff on the test side (both steps); q: identity-mapped,
        # equal across sides but step-dependent (so step 2 is not bitwise
        # identical to step 1, which would trip the repeated-snapshot note)
        x_entry = field if b is base else field_perturbed
        q3d = f8(*range(10 + t, 22 + t))
        return {"x": b.arr([2, 3], x_entry, x_entry),
                "q": b.arr([2, 3, 2], q3d, q3d),
                "errflg": b.scal(0, 0)}

    def beta_args(b, exit_val):
        return {"y": b.arr([2, 3], field, field),
                "n": b.scal(7, exit_val)}

    # both sides run steps 1-2 (offset 0); beta exit diverges at step 2
    for t in (1, 2):
        base.hit("alpha", t - 1, t, "drv", alpha_args(base, t))
        base.hit("beta", t - 1, t, "drv", beta_args(base, 42))
        test.hit("alpha", t - 1, t, "drv", alpha_args(test, t))
        test.hit("beta", t - 1, t, "drv",
                 beta_args(test, 42 if t == 1 else 41))
    base.write()
    test.write()

    with open(os.path.join(outdir, "suite.json"), "w") as f:
        json.dump({"schemes": ["alpha", "beta"], "steps": 2,
                   "step_offset": 0,
                   "roles": {"c": "base", "s": "test"},
                   "step_sentinel": "phys_step",
                   "intents": {}, "portable": {},
                   "version": "testver",
                   "base_case": "/scratch/one", "test_case": "/scratch/two"},
                  f)

    ok = differ.report(outdir, ["alpha", "beta"], 2)
    assert not ok, "differ should report a divergence"

    with open(os.path.join(outdir, "report.txt")) as f:
        text = f.read()

    # role labels everywhere, no historical-label leakage
    assert re.search(r"\bsima\b", text, re.I) is None, "sima label leaked"
    assert re.search(r"\bcam\b", text, re.I) is None, "cam label leaked"
    assert "test=9 base=1" in text, "entry diff not labeled base/test"
    assert "test=41 base=42" in text, "exit diff not labeled base/test"
    assert "[OUTPUTS DIFFER]" in text and "[INPUTS DIFFER]" in text
    assert "constituent mapping (BASE 2 species, TEST 2):" in text
    assert "base q(:,:, 1) Q" in text, "identity constituent map missing"
    assert "test step t <-> base step t+0" in text

    # offset-0 alignment verdict must score the SAME step, not t+1
    assert "the paired base step 1" in text, \
        "alignment verdict ignored step_offset"
    assert "timestep alignment is\ncorrect" in text
    assert "*** WARNING" not in text and "*** note" not in text, \
        "offset-0 pairing raised a spurious alignment warning"

    shutil.rmtree(outdir)
    print("\nsame-model test passed: base/test labels, offset-0 pairing + "
          "alignment verdict, identity constituent map, targets-v1 loader")


if __name__ == "__main__":
    main()
