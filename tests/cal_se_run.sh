#!/bin/bash
# Run the derived-type (--capture) micro-calibration on a machine with gdb +
# a Fortran compiler. Usage: bash cal_se_run.sh   (run from anywhere;
# cal_se.f90 sits next to this script, dropsonde_gdb.py in the repo root
# above it). Compiler selectable via $FC (default gfortran), flags via
# $FCFLAGS (default "-g -O0"; derived-type capture is a debug-build path,
# so there is no optimized variant of this calibration).
#
# Expected when everything works (details asserted by the python block):
#   cam role  (kill_after_steps=1): the 4 step-1 hits only, then
#             "collected 1 timesteps; terminating run"; 3 constituent names
#   sima role (kill_after_steps=0): 6 hits over 2 steps; element_t dummies
#             expand into elem%state%v [2,2,2,3,2,4] etc.; the allocatable
#             qdp goes through the per-element address path; tl%n0/np1
#             rotate between steps; fvm (null pointer) skips gracefully;
#             entry capture lands on the first executable statement
#             ("no errmsg idiom" note in gdb.log)
set -e
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
FC=${FC:-gfortran}
FCFLAGS=${FCFLAGS:--g -O0}
echo "=== compiling cal_se.f90 with $FC $FCFLAGS ==="
$FC $FCFLAGS -o cal_se "$SCRIPT_DIR/cal_se.f90"
KILL_cam=1
KILL_sima=0
for role in cam sima; do
  eval kill_steps=\$KILL_$role
  rm -rf out_se_$role && mkdir -p out_se_$role
  cat > out_se_$role/config.json <<EOF
{"role": "$role", "out_dir": "$PWD/out_se_$role",
 "schemes": ["prim_step_like", "compute_like"],
 "portable": {"prim_step_like": "prim_step_like",
              "compute_like": "compute_like"},
 "capture": {"element_t": ["state%v", "state%t", "state%qdp",
                           "derived%ft", "spheremp", "localid", "tag"],
             "timelevel_t": ["n0", "np1"],
             "hvcoord_t": ["hyai", "ps0"]},
 "kill_after_steps": $kill_steps}
EOF
  # -u PYTHONHOME/PYTHONPATH: some sites (e.g. Izumi) point these at an
  # anaconda whose stdlib version mismatches gdb's embedded Python.
  DROPSONDE_CONFIG=$PWD/out_se_$role/config.json \
    env -u PYTHONHOME -u PYTHONPATH \
    gdb --batch -x "$SCRIPT_DIR/../dropsonde_gdb.py" ./cal_se \
    > out_se_$role/gdb.log 2>&1 || true
done
echo "=== cam gdb.log ==="; cat out_se_cam/gdb.log
echo "=== sima gdb.log ==="; cat out_se_sima/gdb.log
python3 - <<'EOF'
import json, os, struct

NP, NLEV, QSIZE, TLV, NELEM = 2, 3, 2, 2, 4


def vals(d, fn):
    raw = open(os.path.join(d, fn), 'rb').read()
    return struct.unpack('={}d'.format(len(raw) // 8), raw)


def ivals(d, fn, width):
    raw = open(os.path.join(d, fn), 'rb').read()
    code = {4: 'i', 8: 'q'}[width]
    return struct.unpack('={}{}'.format(len(raw) // width, code), raw)


def li(subs, extents):
    """0-based linear index of 1-based subscripts, first dim fastest."""
    idx, mult = 0, 1
    for s, e in zip(subs, extents):
        idx += (s - 1) * mult
        mult *= e
    return idx


def check(cond, what):
    if not cond:
        raise SystemExit('FAIL: ' + what)
    print('   ok: ' + what)


for role in ('cam', 'sima'):
    d = 'out_se_' + role
    m = json.load(open(os.path.join(d, 'manifest.json')))
    print('=== {} manifest ==='.format(role))
    print('constituents:', m['constituents'])
    print('breakpoints:', m['breakpoints'])
    for r in m['hits']:
        print('{} hit {} step {} caller {} complete={}'.format(
            r['scheme'], r['hit'], r.get('step'), r.get('caller'),
            r.get('complete', False)))
        for a, i in sorted(r['args'].items()):
            print('   {:16s} {:7s} {:4s} ext={} in={} {}'.format(
                a, str(i.get('kind')), str(i.get('dtype')),
                i.get('extents'), i.get('entry_value'),
                i.get('why', '')))

    check(all(v != 'missing' for v in m['breakpoints'].values()),
          role + ': all breakpoints resolved')
    psl = [r for r in m['hits'] if r['scheme'] == 'prim_step_like']
    cl = [r for r in m['hits'] if r['scheme'] == 'compute_like']
    if role == 'cam':
        check(len(psl) == 2 and len(cl) == 2,
              'cam: 4 step-1 hits (kill after 1 step)')
        check(all(r.get('step') == 1 for r in m['hits']),
              'cam: all hits tagged step 1')
        check(len(m['constituents']) == 3, 'cam: 3 constituent names')
        continue

    check(len(psl) == 3 and len(cl) == 3, 'sima: 3+3 hits')
    check([r.get('step') for r in psl] == [1, 1, 2],
          'sima: prim_step_like steps 1,1,2')
    check(all(r.get('caller', '').endswith('prim_step_like')
              for r in cl), 'sima: compute_like caller bucketed')
    check(m['constituents'] == ['std_name_1', 'std_name_2'],
          'sima: constituent names')

    h0 = psl[0]['args']
    # --- expansion shapes -------------------------------------------------
    check('expanded into 7 component' in h0['elem'].get('why', ''),
          'elem parent marked expanded')
    exp_ext = {'elem%state%v': [NP, NP, 2, NLEV, TLV, NELEM],
               'elem%state%t': [NP, NP, NLEV, TLV, NELEM],
               'elem%state%qdp': [NP, NP, NLEV, QSIZE, 2, NELEM],
               'elem%derived%ft': [NP, NP, NLEV, NELEM],
               'elem%spheremp': [NP, NP, NELEM],
               'elem%localid': [NELEM]}
    for k, e in sorted(exp_ext.items()):
        check(h0[k].get('kind') == 'array' and h0[k].get('extents') == e,
              '{} extents {}'.format(k, e))
    check(h0['elem%state%qdp'].get('plan', {}).get('addrs') is not None and
          len(h0['elem%state%qdp']['plan']['addrs']) == NELEM,
          'qdp (allocatable) took the per-element address path')
    check(h0['elem%state%v'].get('plan', {}).get('addrs') is None,
          'v (fixed shape) took the strided single-plan path')
    check(h0['elem%tag'].get('kind') == 'skipped',
          'character component skipped')
    check(h0['fvm'].get('kind') in ('skipped', 'error'),
          'null fvm pointer array handled gracefully')
    # --- scalar-struct dummies ---------------------------------------------
    check(h0['tl%n0'].get('entry_value') == 1 and
          h0['tl%np1'].get('entry_value') == 2, 'tl indices at step 1')
    h2 = psl[2]['args']
    check(h2['tl%n0'].get('entry_value') == 2 and
          h2['tl%np1'].get('entry_value') == 1, 'tl rotation at step 2')
    hy = vals(d, h0['hvcoord%hyai']['entry_file'])
    check(all(abs(x - 0.1 * (i + 1)) < 1e-12 for i, x in enumerate(hy)),
          'hvcoord%hyai values')
    check(h0['hvcoord%ps0'].get('entry_value') == 1000.0, 'hvcoord%ps0')
    # --- values: dim order, per-element strides, entry/exit relations ------
    tve = vals(d, h0['elem%state%t']['entry_file'])
    text = exp_ext['elem%state%t']
    for ie in range(1, NELEM + 1):
        check(abs(tve[li((1, 1, 1, 1, ie), text)] - (100 * ie + 20)) < 1e-12,
              't entry base value elem {}'.format(ie))
        check(abs(tve[li((1, 2, 3, 1, ie), text)] -
                  (100 * ie + 20.5)) < 1e-12,
              't interior bump at (1,2,3,1) elem {}'.format(ie))
    vv = vals(d, h0['elem%state%v']['entry_file'])
    vext = exp_ext['elem%state%v']
    check(abs(vv[li((2, 1, 2, 3, 1, 1), vext)] - 110.25) < 1e-12,
          'v interior bump at (2,1,2,3,1) elem 1')
    qv = vals(d, h0['elem%state%qdp']['entry_file'])
    qext = exp_ext['elem%state%qdp']
    check(abs(qv[li((2, 1, 3, 2, 1, 1), qext)] - 130.125) < 1e-12,
          'qdp interior bump at (2,1,3,2,1) elem 1')
    tvx = vals(d, h0['elem%state%t']['exit_file'])
    ok = all(abs(tvx[li((i, j, k, 2, ie), text)] -
                 (tve[li((i, j, k, 1, ie), text)] + 0.5)) < 1e-12
             for ie in range(1, NELEM + 1) for k in range(1, NLEV + 1)
             for j in range(1, NP + 1) for i in range(1, NP + 1))
    check(ok, 't exit(np1) == entry(n0) + dt everywhere')
    qx = vals(d, h0['elem%state%qdp']['exit_file'])
    ok = all(abs(qx[li((i, j, k, q, 2, ie), qext)] -
                 2 * qv[li((i, j, k, q, 1, ie), qext)]) < 1e-12
             for ie in range(1, NELEM + 1) for q in range(1, QSIZE + 1)
             for k in range(1, NLEV + 1) for j in range(1, NP + 1)
             for i in range(1, NP + 1))
    check(ok, 'qdp exit(:,2) == 2 * entry(:,1) everywhere')
    lid_w = int(h0['elem%localid']['dtype'][1:])
    lid_in = ivals(d, h0['elem%localid']['entry_file'], lid_w)
    lid_out = ivals(d, h0['elem%localid']['exit_file'], lid_w)
    check(lid_in == (1, 2, 3, 4) and lid_out == (2, 3, 4, 5),
          'localid scalar component per element, exit == entry + 1')
    vx = vals(d, h0['elem%state%v']['exit_file'])
    check(abs(vx[li((1, 1, 1, 1, 2, 1), vext)] -
              (vv[li((1, 1, 1, 1, 2, 1), vext)] + 250.0)) < 1e-12,
          'v exit reflects the nested compute_like update (+dt2*ps0)')
    log = open(os.path.join(d, 'gdb.log')).read()
    check('no errmsg idiom; using first executable statement' in log,
          'entry advance used the first-executable-statement fallback')

print('cal_se: ALL CHECKS PASSED')
EOF
