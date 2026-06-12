# dropsonde

Find the first CCPP scheme where CAM and CAM-SIMA physics answers diverge,
without modifying or recompiling either model. Runs both executables under
gdb with breakpoints at every scheme `_run` entry/exit listed in the SDF,
dumps all scheme arguments to disk, and diffs them byte-for-byte.

## Usage

```
./dropsonde --cam  /glade/derecho/scratch/.../cam_case  \
            --sima /glade/derecho/scratch/.../sima_case \
            --sdf  suite_park_macrop.xml
```

Run it where the models can run (i.e., a compute node). Both gdb sessions
launch in parallel; the report prints to stdout when they finish.
Re-run just the comparison with `python3 differ.py <out_dir>`.

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
   together as toplevel. `errmsg`/`errflg` are not compared (uninitialized
   at entry on the CAM side). Constituent indices are matched via a
   short-name ↔ standard-name table (`SPECIAL_CONST_MAP` — extend it as
   ports grow). The report, in execution order:
   - inputs differ → wiring problem upstream of the scheme;
   - inputs match, outputs differ → problem inside the scheme;
   - shape/dtype/kind mismatches, hit-count mismatches, SIMA-only schemes;
   - the constituent index mapping (useful for manual gdb sessions);
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

Still to verify on re-run (fixes applied since the first attempt):

1. Explicit-shape dummies (`b(ncol,pver)`) capture after the
   declaration-line `next` (was: `no such vector element`).
2. SIMA constituent chase with `cam_constituents::` spelling (linker-name
   spelling hit "unknown type" on gdb 16.2).
3. FinishBreakpoint exits: hit records have `"complete": true` and
   `exit_file` entries; `CHECK b == 2a: True` in cal_run.sh.
4. Pending-breakpoint detection (the Python API ignores
   `set breakpoint pending off`): missing symbols must report `missing`,
   not resolve silently.
5. Step tagging + auto-kill: cam role stops after 1 "timestep" in
   cal_run.sh; hits carry correct `step`.
6. Then the real known-BFB case must report zero diffs end-to-end; if the
   first scheme's inputs differ the report flags alignment as suspect.

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
