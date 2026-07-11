# Synthetic test of differ.py for pseudo-SDF (SE dycore) captures:
#   - step_offset 0 (two free-running models): sima step t pairs with
#     cam step t, driven by suite.json's step_offset
#   - derived-type component pseudo-args (elem%state%t) inherit the intent
#     of their root dummy: intent(in) roots not compared at exit,
#     intent(out) roots not compared at entry
#   - portable=self (the pseudo-SDF convention) displays as a plain name
#   - constituent-count mismatch (CAM 4 species vs SIMA 3): qdp-style args
#     compare the mapped species along the (interior) constituent axis
#     instead of reporting a shape mismatch
#   - uninitialized intent(inout) scratch (omega_cn-style): entry garbage
#     downgraded to a note, exit compared over written elements only
#
# Run: python3 tests/se_differ_test.py

import json
import os
import shutil
import struct
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import differ  # noqa: E402
from synthetic_differ_test import Builder, f8, Q_STD, CLDLIQ_STD  # noqa: E402

RAIN_STD = "rain_mixing_ratio_wrt_moist_air_and_condensed_water"


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


def qdp_bytes(species_vals, ni=2, nj=2):
    """Fortran-order bytes of a (ni, nspecies, nj) f8 array whose species
    slab s is the constant species_vals[s] (constituent axis INTERIOR, like
    qdp's)."""
    ns = len(species_vals)
    out = []
    for j in range(nj):
        for s in range(ns):
            out.extend([species_vals[s]] * ni)
    return struct.pack("={}d".format(ni * ns * nj), *out)


def mapped_and_uninit():
    """Constituent-count-mismatch mapping + uninitialized-entry handling."""
    outdir = tempfile.mkdtemp(prefix="dropsonde_se_test2_")
    with open(os.path.join(outdir, "suite.json"), "w") as f:
        json.dump({"schemes": ["prim_step"], "steps": 1,
                   "step_offset": 0,
                   "portable": {"prim_step": "prim_step"}}, f)

    cam = Builder(outdir, "cam")
    sima = Builder(outdir, "sima")
    for b in (cam, sima):
        b.man["breakpoints"] = {"prim_step": "prim_step"}
    # CAM registers one extra (test-tracer) constituent; orders differ
    cam.man["constituents"] = ["Q", "CLDLIQ", "RAINQM", "CL"]
    sima.man["constituents"] = [CLDLIQ_STD, RAIN_STD, Q_STD]

    # mapped entry values agree; at exit RAINQM differs; CL (unmatched)
    # is wildly different throughout and must never be compared
    cam_qdp_in = qdp_bytes([1.0, 2.0, 3.0, 99.0])
    sima_qdp_in = qdp_bytes([2.0, 3.0, 1.0])
    cam_qdp_out = qdp_bytes([1.0, 2.0, 3.5, 88.0])
    sima_qdp_out = qdp_bytes([2.0, 3.0, 1.0])
    # fq: CLDLIQ already differs at entry -> INPUTS DIFFER + exit suppressed
    cam_fq_in = qdp_bytes([4.0, 5.0, 6.0, 77.0])
    sima_fq_in = qdp_bytes([5.5, 6.0, 4.0])

    # omega_cn: entry garbage on both sides; exits agree on every written
    # element; the third element is never written (entry == exit bitwise in
    # each model) and still differs across models
    om_in_s = f8(1.1285e277, 3.9e-315, 6.66e-313)
    om_in_c = f8(7e-320, 2.2e280, 1.2e290)
    om_out_s = f8(0.1, 0.2, 6.66e-313)
    om_out_c = f8(0.1, 0.2, 1.2e290)
    # omega_bad: same entry garbage, but a WRITTEN element differs at exit
    ob_out_s = f8(0.1, 0.2, 0.3)
    ob_out_c = f8(0.1, 0.25, 0.3)

    def args(b, qdp_in, qdp_out, fq_in, om_in, om_out, ob_out):
        na = 4 if b.man["role"] == "cam" else 3
        return {"elem": {"kind": "skipped",
                         "why": "derived-type array element_t: expanded "
                                "into 2 component pseudo-args"},
                "elem%state%qdp": b.arr([2, na, 2], qdp_in, qdp_out),
                "elem%derived%fq": b.arr([2, na, 2], fq_in, fq_in),
                "omega_cn": b.arr([3], om_in, om_out),
                "omega_bad": b.arr([3], om_in, ob_out),
                "nets": b.scal(1, 1)}

    cam.hit("prim_step", 0, 1, "dyn_run",
            args(cam, cam_qdp_in, cam_qdp_out, cam_fq_in, om_in_c,
                 om_out_c, ob_out_c))
    sima.hit("prim_step", 0, 1, "dyn_run",
             args(sima, sima_qdp_in, sima_qdp_out, sima_fq_in, om_in_s,
                  om_out_s, ob_out_s))
    cam.write()
    sima.write()

    intents = {"prim_step": {"elem": "inout", "nets": "in"}}
    ok = differ.report(outdir, ["prim_step"], 1, intents)
    assert not ok, "differ should report divergences"

    with open(os.path.join(outdir, "report.txt")) as f:
        rep = f.read()
    # qdp: no shape-mismatch bail; exit diff isolated to the mapped RAINQM
    assert "shape mismatch" not in rep, \
        "count-mismatched qdp fell back to shape mismatch:\n" + rep
    assert ("[OUTPUTS DIFFER]" in rep and "elem%state%qdp" in rep and
            "1/3 mapped species differ: RAINQM" in rep), \
        "mapped qdp exit diff not isolated to RAINQM:\n" + rep
    assert "CL (cam idx" not in rep, \
        "unmatched CL species leaked into a comparison:\n" + rep
    # fq: entry diff on CLDLIQ, exit suppressed
    assert ("[INPUTS DIFFER]" in rep and "elem%derived%fq" in rep and
            "CLDLIQ" in rep), "mapped fq entry diff missing:\n" + rep
    assert "suppressed" in rep and "elem%derived%fq" in rep, \
        "fq exit not suppressed after entry diff:\n" + rep
    # omega_cn: entry garbage -> note, not INPUTS DIFFER; exit clean
    assert "look uninitialized" in rep, \
        "uninit entry note missing:\n" + rep
    assert "arg omega_cn: shape" not in rep and \
        "[INPUTS DIFFER] prim_step via dyn_run [step 1] (hit 0) " \
        "arg omega_cn" not in rep, \
        "uninit omega_cn counted as an input diff:\n" + rep
    assert "[OUTPUTS DIFFER] prim_step via dyn_run [step 1] (hit 0) " \
        "arg omega_cn" not in rep, \
        "never-written omega_cn element compared at exit:\n" + rep
    # omega_bad: exit must still be compared (not suppressed) and differ
    assert "[OUTPUTS DIFFER] prim_step via dyn_run [step 1] (hit 0) " \
        "arg omega_bad" in rep, \
        "uninit-entry arg's real exit diff was lost:\n" + rep

    shutil.rmtree(outdir)
    print("\nse mapped/uninit test passed: count-mismatched constituent "
          "axis mapped per species, uninit entry downgraded to note, "
          "written-only exit comparison behaved")


if __name__ == "__main__":
    main()
    mapped_and_uninit()
