#!/bin/bash
# Run the dropsonde_gdb.py micro-calibration on a machine with gdb + a
# Fortran compiler. Usage: bash cal_run.sh   (run from anywhere; cal.f90 is
# located next to this script and dropsonde_gdb.py in the repo root above it).
# Compiler selectable via $FC (default gfortran); both gfortran and ifort
# accept -g -O0. Used to calibrate per machine (Derecho gfortran, Izumi
# Intel, ...).
#
# Flags selectable via $FCFLAGS (default "-g -O0"). To calibrate the
# optimized-binary capture mode (production builds), use the production
# flags plus -g and -fno-inline, e.g.:
#   FCFLAGS="-O2 -ffp-contract=off -g -fno-inline" bash cal_run.sh
# (-fno-inline stands in for the cross-TU calls of the real model, keeping
# the nested inner hit out-of-line.)
#
# Expected when everything works:
#   cam role  (kill_after_steps=1): outer+inner step-1 hits only, then
#             "collected 1 timesteps; terminating run"; 5 constituent names
#   sima role (kill_after_steps=0): 3 hits; nested inner has caller
#             phys::outer_run and step 1; strided inner hit has step 2,
#             ext=[3,3] with a non-contiguous stride, and
#             CHECK a = (1,2,3,5,6,7,9,10,11); all CHECK b == 2a: True;
#             3 constituent names std_name_1..3
# Optimized (-O2) runs differ by design: breakpoints show as *0x<addr>, and
# b (explicit-shape, runtime bounds) is skipped, so CHECK b lines vanish.
set -e
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
FC=${FC:-gfortran}
FCFLAGS=${FCFLAGS:--g -O0}
echo "=== compiling cal.f90 with $FC $FCFLAGS ==="
if echo "$FCFLAGS" | grep -qE -- "-O[123sz]"; then
  # Optimized calibration: split the program unit into its own translation
  # unit. In the real model every SDF scheme call is cross-TU (no LTO), so
  # callers populate full array descriptors; single-TU -O2 lets gcc IPA
  # dead-store-eliminate descriptor fields the callee provably never reads,
  # which is unrepresentative and breaks the raw-descriptor decode.
  awk "/^program cal/{exit} {print}" "$SCRIPT_DIR/cal.f90" > cal_mods.f90
  awk "/^program cal/{f=1} f{print}" "$SCRIPT_DIR/cal.f90" > cal_main.f90
  $FC $FCFLAGS -c cal_mods.f90
  $FC $FCFLAGS cal_main.f90 cal_mods.o -o cal
else
  $FC $FCFLAGS -o cal "$SCRIPT_DIR/cal.f90"
fi
KILL_cam=1
KILL_sima=0
for role in cam sima; do
  eval kill_steps=\$KILL_$role
  rm -rf out_$role && mkdir -p out_$role
  cat > out_$role/config.json <<EOF
{"role": "$role", "out_dir": "$PWD/out_$role",
 "schemes": ["inner", "outer"], "kill_after_steps": $kill_steps}
EOF
  # -u PYTHONHOME/PYTHONPATH: some sites (e.g. Izumi) point these at an
  # anaconda whose stdlib version mismatches gdb's embedded Python, which
  # then aborts with "No module named 'encodings'". dropsonde's gdb side is
  # stdlib-only, so stripping them is always safe.
  DROPSONDE_CONFIG=$PWD/out_$role/config.json \
    env -u PYTHONHOME -u PYTHONPATH \
    gdb --batch -x "$SCRIPT_DIR/../dropsonde_gdb.py" ./cal > out_$role/gdb.log 2>&1 || true
done
echo "=== cam gdb.log ==="; cat out_cam/gdb.log
echo "=== sima gdb.log ==="; cat out_sima/gdb.log
python3 - <<'EOF'
import json, os, struct

def vals(d, fn):
    raw = open(os.path.join(d, fn), 'rb').read()
    return struct.unpack('={}d'.format(len(raw) // 8), raw)

for role in ('cam', 'sima'):
    d = 'out_' + role
    m = json.load(open(os.path.join(d, 'manifest.json')))
    print('=== {} manifest ==='.format(role))
    print('constituents:', m['constituents'])
    print('breakpoints:', m['breakpoints'])
    for r in m['hits']:
        print('{} hit {} step {} caller {} complete={}'.format(
            r['scheme'], r['hit'], r.get('step'), r.get('caller'),
            r.get('complete', False)))
        for a, i in sorted(r['args'].items()):
            print('   {:8s} {:7s} {:4s} ext={} str={} in={} out={} '
                  'files={}/{} {}'.format(
                      a, str(i.get('kind')), str(i.get('dtype')),
                      i.get('extents'), i.get('strides'),
                      i.get('entry_value'), i.get('exit_value'),
                      'entry_file' in i, 'exit_file' in i,
                      i.get('why', '')))
        ai = r['args'].get('a', {})
        bi = r['args'].get('b', {})
        if r['scheme'] == 'inner' and 'entry_file' in ai:
            av = vals(d, ai['entry_file'])
            print('   CHECK a:', tuple(int(x) for x in av))
            if 'exit_file' in bi:
                bv = vals(d, bi['exit_file'])
                ok = (len(av) == len(bv) and
                      all(abs(y - 2 * x) < 1e-12 for x, y in zip(av, bv)))
                print('   CHECK b == 2a:', ok)
EOF
