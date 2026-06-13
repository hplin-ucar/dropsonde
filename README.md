# dropsonde

Find the first CCPP scheme where CAM and CAM-SIMA physics answers diverge,
without modifying or recompiling either model. Runs both executables under
gdb with breakpoints at every scheme `_run` entry/exit listed in the SDF,
dumps all scheme arguments to disk, and diffs them byte-for-byte.

## Usage

```
./dropsonde --cam  /glade/derecho/scratch/.../cam_case  \
            --sima /glade/derecho/scratch/.../sima_case \
            --sdf  suite_park_macrop.xml \
            [--meta-root /path/to/atmospheric_physics] [--steps 1]
```

`--steps` defaults to 1: in FPHYStest snapshot runs every timestep is
independent (the model state is re-read from `ncdata` each step), so
comparing CAM timestep 2 against SIMA timestep 1 is decisive — if one
timestep matches, the rest will too (barring `is_first_timestep` logic,
for which `--steps 2` compares an extra pair).

`--meta-root` (default: three levels up from the SDF) locates the schemes'
CCPP `.meta` files; their intents let the differ skip `intent(out)` args
at entry (caller-side garbage, e.g. 67 of bretherton_park_diff's 113
args) and `intent(in)` args at exit.

Run it where the models can run (i.e., a compute node). Both gdb sessions
launch in parallel; the report prints to stdout when they finish.
Re-run just the comparison with `python3 differ.py <out_dir>`.

WARNING: any direct `cesm.exe` execution (dropsonde, manual gdb) reopens
and TRUNCATES the component logs named in `nuopc.runconfig` — i.e. it
clobbers `atm.log.<last-jobid>.*` etc. from the most recent batch run.
Copy aside any logs you want to keep before running.

## Prerequisites

- Case dirs contain `bld/cesm.exe` and a populated `run/` directory
  (`preview_namelists` already run).
- CAM: GNU, `DEBUG=TRUE`, `NTASKS=1`, `NTHRDS=1`, dechunked
  (`pcols` = total columns). Must call the same CCPPized
  `atmospheric_physics` routines as the SIMA build.
- CAM-SIMA: GNU, `DEBUG=TRUE`, single rank.
- `STOP_N` only needs to be *long enough* (CAM's 6-step coupling minimum
  is fine): each gdb session counts timesteps via a sentinel breakpoint on
  `cam_run1` and kills its run once it has what it needs (CAM:
  `--steps + 1` timesteps, SIMA: `--steps`).
- `gdb` on PATH. No numpy needed anywhere — gdb side and differ are
  stdlib-only, Python >= 3.6 (works with gdb's embedded interpreter and
  `/usr/bin/python3` on Derecho).

## How it works

1. `dropsonde` parses the SDF for the ordered scheme list and launches
   `gdb --batch -x dropsonde_gdb.py` on each executable (config passed via
   `$DROPSONDE_CONFIG`).
2. `dropsonde_gdb.py` breaks on every `<scheme>_run`. At entry it first
   steps over the declaration line (gfortran materializes explicit-shape
   dummy bounds there, past gdb's post-prologue stop), then walks the
   frame's dummy arguments, probes each array's base address and
   per-subscript byte strides empirically (element-address arithmetic, so
   strided actual args like `state%q(:ncol,:,m)` are handled), dumps raw
   bytes, and plants a `FinishBreakpoint` that re-reads the same addresses
   at scheme exit. Every hit is tagged with the current timestep (sentinel
   breakpoint on `cam_run1`) and its caller, and the run is killed after
   enough timesteps.
3. At the first hit it also captures the ordered constituent name lists:
   CAM `constituents::cnst_name`, SIMA
   `cam_constituents::const_props(i)%prop%var_std_name`.
4. `differ.py` aligns offline: SIMA timestep t pairs with CAM timestep
   t+1 (CAM's first dumped step is skipped). Within a step, hits are
   bucketed by caller: a scheme called from inside another scheme (e.g.
   the diffusion solver inside `bretherton_park`'s iteration, or
   `compute_cloud_fraction` inside RK) has the *same shared-code caller
   symbol* in both binaries and is compared against its nested
   counterpart; model-specific call sites (CCPP cap vs CAM driver) bucket
   together as toplevel. `errmsg`/`errflg`/`iulog` are never compared.
   Intent filtering (from `.meta`) drops entry-garbage and
   input-repeated-at-exit noise; exit reports are suppressed for args
   whose entry already differed, and `intent(out)` exit elements that
   BOTH models left untouched (entry == exit bitwise) are ignored —
   partially-defined outputs (e.g. `dtk`'s surface row) otherwise echo
   each caller's unrelated buffer contents. Constituent-indexed arrays —
   rank-3 `q(ncol,pver,pcnst)`, rank-2 fluxes `cflx(ncol,pcnst)`, and
   rank-1 per-species config like `qmincg(pcnst)` — are compared per
   mapped species via a short-name ↔ standard-name table
   (`SPECIAL_CONST_MAP`, then exact-name, then `cnst_<short name>` —
   constituents auto-registered from a snapshot file keep their netcdf
   variable name as their standard name). When hit
   counts differ within a bucket (CAM hosts call `geopotential_temp`
   from ~20 sites per step), each SIMA hit is paired with the CAM hit
   whose comparable inputs ALL match bitwise (constituent arrays
   per-species), or reported unpaired with the closest candidate. The
   report, in execution order:
   - inputs differ → wiring problem upstream of the scheme;
   - inputs match, outputs differ → problem inside the scheme;
   - shape/dtype/kind mismatches, hit-count mismatches, SIMA-only schemes;
   - the constituent index mapping (useful for manual gdb sessions);
   - if the very first compared scheme already has input differences, an
     alignment warning plus an *offset scan*: the first hit's entry args
     are checked bitwise against every dumped CAM step (wrong offset →
     some other step matches; no step matches → the two runs are
     different trajectories, e.g. the snapshot driving CAM-SIMA came
     from a different CAM build) and against the other SIMA steps
     (identical → the snapshot time record is being re-read);
   - bit-for-bit: `No differences found in any subroutines!`

Entry-time addresses of every argument are kept in
`<out>/<role>/manifest.json` — handy for setting watchpoints manually
afterwards.

## Calibration status

Micro-calibration (`tests/cal.f90` + `tests/cal_run.sh`, run on a machine
with gdb/gfortran — Derecho login node has gdb 16.2 + gfortran 12.5 in the
default `bash -l` environment) verified on 2026-06-12:

- plain `<scheme>_run` breakpoints resolve; caller capture works
  (`phys::outer_run` form, identical across binaries for shared code)
- array extents/strides pair correctly (`ext=[4,3] str=[8,32]`,
  3-D `[8,32,96]`) with `DIM_ORDER = "reversed"` + per-array sizeof check
- CAM `cnst_name` capture works

Full-model run (UW PBL suite, 2026-06-12) verified end-to-end mechanics:
step tagging + auto-kill ("collected 3 timesteps; terminating run"),
explicit-shape capture, exits, caller bucketing (nested
`bretherton_park_diff_run` calls of the diffusion routines pair
correctly), CAM constituent names, pending-breakpoint detection
(13 CAM-missing schemes correctly classified). The report correctly
flagged alignment as suspect; manual offset scanning showed the
snapshot driving the SIMA case came from a different CAM trajectory
(no CAM step matched bitwise) — that scan is now automated.

A second full run verified intent filtering end-to-end and exposed that
deferred-length strings need the hidden-length path: gdb prints SIMA's
`const_props(i)%prop%var_std_name` as `(PTR TO -> character*0) 0x...`
because DWARF carries no dynamic length on the component. The capture
now reads the hidden sibling member `_var_std_name_length` (gfortran
appends it after the type's own visible components; verified against
gfortran DWARF output) and then reads exactly that many bytes at the
data pointer.

The second run also validated the differ on real data: with the
constituent permutation recovered bitwise from the q dumps themselves,
the decisive timestep pair (CAM nstep 1 vs SIMA nstep 0) was bit-for-bit
across the whole UW PBL suite except for genuine, separately-explained
findings (SIMA reading NUMLIQ/NUMICE as zero from the snapshot; the
known `do_diffusion_const` NUMLIQ flag mismatch; SIMA passing `qmin`
where CAM passes `qmincg`; `update_dry_static_energy` invoked from a
different host phase — `dp_coupling` in CAM).

Still to verify on the next run:

1. SIMA constituent names via the hidden `_<name>_length` member path.
2. A true BFB pairing (with the SIMA-side input findings fixed) must
   report zero diffs.

## Known limitations

- Derived-type and pointer dummy arguments are skipped (recorded in the
  manifest with reason). CCPP-compliant schemes should not have them,
  except `ccpp_constituent_prop_ptr_t` arrays, which are skipped too.
- NAG/Izumi unsupported for now: the f2c-style lowering mangles names and
  loses descriptors. The stride-probe design degrades to needing explicit
  shapes (e.g. from `.meta` files) — future work.
- Negative-stride (reversed) array slices are skipped with a note.
- One subcycle scheme called with *different arguments* per sub-step
  aligns fine (hits are sequence-matched), but the report labels hits by
  per-scheme index, not subcycle position.
