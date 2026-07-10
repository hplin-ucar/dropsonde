# Synthetic test of differ.py for pseudo-SDF (SE dycore) captures:
#   - step_offset 0 (two free-running models): sima step t pairs with
#     cam step t, driven by suite.json's step_offset
#   - derived-type component pseudo-args (elem%state%t) inherit the intent
#     of their root dummy: intent(in) roots not compared at exit,
#     intent(out) roots not compared at entry
#   - portable=self (the pseudo-SDF convention) displays as a plain name
#
# Run: python3 tests/se_differ_test.py

import json
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import differ  # noqa: E402
from synthetic_differ_test import Builder, f8  # noqa: E402


def main():
    outdir = tempfile.mkdtemp(prefix="dropsonde_se_test_")
    with open(os.path.join(outdir, "suite.json"), "w") as f:
        json.dump({"schemes": ["prim_step"], "steps": 1,
                   "step_offset": 0,
                   "portable": {"prim_step": "prim_step"}}, f)

    cam = Builder(outdir, "cam")
    sima = Builder(outdir, "sima")
    for b in (cam, sima):
        b.man["breakpoints"] = {"prim_step": "prim_step"}

    field = f8(1, 2, 3, 4)

    def args(b, t_exit, tl_exit, junk_entry):
        return {"elem": {"kind": "skipped",
                         "why": "derived-type array element_t: expanded "
                                "into 1 component pseudo-args"},
                "elem%state%t": b.arr([2, 2], field, t_exit),
                "tl%n0": b.scal(1, tl_exit),
                "junk%x": b.scal(junk_entry, 9),
                "nets": b.scal(1, 1)}

    # step 1 in BOTH models (offset 0). Injected differences:
    #   elem%state%t exit  -> must be reported (root elem is inout)
    #   tl%n0 exit         -> must NOT be (root tl is intent(in))
    #   junk%x entry       -> must NOT be (root junk is intent(out))
    cam.hit("prim_step", 0, 1, "dyn_run", args(cam, f8(9, 9, 9, 9), 1, 5))
    sima.hit("prim_step", 0, 1, "dyn_run", args(sima, f8(9, 9, 9, 8), 2, 7))
    cam.write()
    sima.write()

    intents = {"prim_step": {"elem": "inout", "tl": "in", "junk": "out",
                             "nets": "in"}}
    ok = differ.report(outdir, ["prim_step"], 1, intents)
    assert not ok, "differ should report the elem%state%t divergence"

    with open(os.path.join(outdir, "report.txt")) as f:
        rep = f.read()
    assert "cam step t+0" in rep, "step_offset 0 not applied:\n" + rep
    # dycore callers are the same symbol in both models, so hits bucket
    # under the shared caller and the label carries "via <caller>"
    assert "[OUTPUTS DIFFER] prim_step via dyn_run [step 1]" in rep and \
        "elem%state%t" in rep, "component diff not reported:\n" + rep
    assert "arg tl%n0" not in rep, \
        "intent(in) root leaked into exit comparison:\n" + rep
    assert "arg junk%x" not in rep, \
        "intent(out) root leaked into entry comparison:\n" + rep
    assert "prim_step -> prim_step" not in rep, \
        "portable=self should display plain:\n" + rep
    assert "no CAM hits" not in rep, "offset-0 pairing failed:\n" + rep

    shutil.rmtree(outdir)
    print("\nse differ test passed: offset-0 pairing, pseudo-arg intent "
          "fallback, portable=self label all behaved")


if __name__ == "__main__":
    main()
