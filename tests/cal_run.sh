#!/bin/bash
# Run the dropsonde_gdb.py micro-calibration on a machine with gdb+gfortran.
# Usage: bash cal_run.sh   (expects cal.f90 and dropsonde_gdb.py in cwd)
#
# Expected when everything works:
#   cam role  (kill_after_steps=1): outer+inner step-1 hits only, then
#             "collected 1 timesteps; terminating run"; 5 constituent names
#   sima role (kill_after_steps=0): 3 hits; nested inner has caller
#             phys::outer_run and step 1; strided inner hit has step 2,
#             ext=[3,3] with a non-contiguous stride, and
#             CHECK a = (1,2,3,5,6,7,9,10,11); all CHECK b == 2a: True;
#             3 constituent names std_name_1..3
set -e
gfortran -g -O0 -o cal cal.f90
KILL_cam=1
KILL_sima=0
for role in cam sima; do
  eval kill_steps=\$KILL_$role
  rm -rf out_$role && mkdir -p out_$role
  cat > out_$role/config.json <<EOF
{"role": "$role", "out_dir": "$PWD/out_$role",
 "schemes": ["inner", "outer"], "kill_after_steps": $kill_steps}
EOF
  DROPSONDE_CONFIG=$PWD/out_$role/config.json \
    gdb --batch -x dropsonde_gdb.py ./cal > out_$role/gdb.log 2>&1 || true
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
        if (r['scheme'] == 'inner' and 'entry_file' in ai
                and 'exit_file' in bi):
            av = vals(d, ai['entry_file'])
            bv = vals(d, bi['exit_file'])
            ok = (len(av) == len(bv) and
                  all(abs(y - 2 * x) < 1e-12 for x, y in zip(av, bv)))
            print('   CHECK a:', tuple(int(x) for x in av))
            print('   CHECK b == 2a:', ok)
EOF
