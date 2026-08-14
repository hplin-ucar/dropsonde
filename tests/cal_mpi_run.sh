#!/bin/bash
# End-to-end calibration of same-model MPMD one-rank mode, driven through
# the real ./dropsonde driver (unlike cal_run.sh/cal_se_run.sh, which drive
# dropsonde_gdb.py directly): targets spec, --source-root intents, --launch
# MPMD template, fake base/test case dirs, differ report.
#
# Scenario A (money test): an output-side divergence planted inside
#   flux_calc on rank 0 at step 2 must be localized as [OUTPUTS DIFFER]
#   flux_calc [step 2] arg b at the planted index, with base/test labels
#   and inputs bit-for-bit. The program runs 3 steps but --steps 2 kills
#   both jobs at step 3's sentinel; since the ranks synchronize in
#   MPI_Allreduce every step, mpiexec exiting at all (not the `timeout`
#   backstop) IS the kill-propagation check. The driver's "gdb exited with
#   rc" warning is expected -- mpiexec reports the killed job as abnormal.
#
# Scenario B (documented blind spot): the same divergence planted on rank
#   2 -- which MPMD one-rank mode does NOT instrument -- must produce a
#   CLEAN report. That is the accepted limitation of instrumenting rank 0
#   only; this scenario asserts it so it stays documented. Pick the
#   instrumented rank to own the columns cprnc flagged, or fall back to
#   NTASKS=1, when the divergence may be spatially confined.
#
# Requirements: mpif90/mpiexec (override via $MPIFC/$MPIEXEC), gdb.
# On Izumi: module load compiler/intel/20.0.1, then bash cal_mpi_run.sh.
set -e
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MPIFC=${MPIFC:-mpif90}
MPIEXEC=${MPIEXEC:-mpiexec}
FCFLAGS=${FCFLAGS:--g -O0}
NRANKS=${NRANKS:-4}
echo "=== compiling cal_mpi.f90 with $MPIFC $FCFLAGS ==="
$MPIFC $FCFLAGS -o cal_mpi "$SCRIPT_DIR/cal_mpi.f90"

cat > cal_mpi_targets.json <<EOF
{"spec": "dropsonde-targets-v1",
 "step_sentinel": "step_begin",
 "targets": ["flux_calc"]}
EOF

setup_cases() {  # $1 = perturb.txt line for the test case ('' = none)
  rm -rf case_base case_test
  mkdir -p case_base/bld case_base/run case_test/bld case_test/run
  cp cal_mpi case_base/bld/cesm.exe
  cp cal_mpi case_test/bld/cesm.exe
  if [ -n "$1" ]; then echo "$1" > case_test/run/perturb.txt; fi
}

LAUNCH="$MPIEXEC -n 1 {gdb} : -n $((NRANKS-1)) {exe}"

run_driver() {  # $1 = out dir
  set +e
  timeout 300 "$SCRIPT_DIR/../dropsonde" --base case_base --test case_test \
    --targets cal_mpi_targets.json --source-root "$SCRIPT_DIR" \
    --steps 2 --launch "$LAUNCH" --out "$1" > "$1.log" 2>&1
  rc=$?
  set -e
}

echo "=== scenario A: divergence on the instrumented rank (must be found) ==="
setup_cases "out 0 2 5"
rm -rf out_mpi_a out_mpi_b
run_driver out_mpi_a
rc_a=$rc
echo "(driver rc=$rc_a)"

echo "=== scenario B: divergence on a non-instrumented rank (blind spot) ==="
setup_cases "out 2 2 5"
run_driver out_mpi_b
rc_b=$rc
echo "(driver rc=$rc_b)"

RC_A=$rc_a RC_B=$rc_b python3 - <<'EOF'
import json
import os

rc_a = int(os.environ["RC_A"])
rc_b = int(os.environ["RC_B"])
assert rc_a not in (124, 137), \
    "scenario A hit the timeout backstop: kill did not propagate " \
    "through mpiexec"
assert rc_b not in (124, 137), \
    "scenario B hit the timeout backstop: kill did not propagate"
assert rc_a == 1, "scenario A: driver rc {} (want 1 = divergence)".format(rc_a)
assert rc_b == 0, "scenario B: driver rc {} (want 0 = clean)".format(rc_b)

a = open("out_mpi_a/report.txt").read()
assert "[OUTPUTS DIFFER]" in a and "flux_calc" in a, a[-500:]
assert "[step 2]" in a and "arg b" in a, "wrong step/arg localized"
assert "at (5)" in a, "wrong element localized (perturb planted at b(5))"
assert "test=" in a and "base=" in a, "labels wrong"
# finding lines start "  [TAG] ..."; the reading guide's mention of the
# tag starts "  3. [INPUTS DIFFER]" and must not trip this
assert "\n  [INPUTS DIFFER]" not in a, "inputs must be bit-for-bit in A"

for role in ("base", "test"):
    m = json.load(open("out_mpi_a/{}/manifest.json".format(role)))
    print("{}: mpi_rank={} hits={}".format(
        role, m["mpi_rank"],
        [(r["scheme"], r["step"]) for r in m["hits"]]))
    assert m["mpi_rank"] == 0, \
        "instrumented rank not recorded as 0 (launcher env not seen?)"
    steps = set(r["step"] for r in m["hits"])
    assert steps == {1, 2}, "hits not confined to compared steps: " + \
        str(steps)
    log = open("out_mpi_a/{}/gdb.log".format(role)).read()
    assert "collected 2 timesteps; terminating run" in log, \
        "step-sentinel kill did not fire for " + role

b = open("out_mpi_b/report.txt").read()
assert "No differences found" in b, \
    "blind-spot scenario unexpectedly found a divergence"

print("cal_mpi: ALL CHECKS PASSED")
EOF
