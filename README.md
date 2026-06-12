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
  (`pcols` = total columns), `STOP_N` covering at least `--steps + 1`
  physics timesteps. Must call the same CCPPized `atmospheric_physics`
  routines as the SIMA build.
- CAM-SIMA: GNU, `DEBUG=TRUE`, single rank, `STOP_N` covering exactly
  `--steps` timesteps (the differ infers calls-per-timestep from the SIMA
  hit count).
- `gdb` on PATH. No numpy needed anywhere — gdb side and differ are
  stdlib-only, Python >= 3.6 (works with gdb's embedded interpreter and
  `/usr/bin/python3` on Derecho).

## How it works

1. `dropsonde` parses the SDF for the ordered scheme list and launches
   `gdb --batch -x dropsonde_gdb.py` on each executable (config passed via
   `$DROPSONDE_CONFIG`).
2. `dropsonde_gdb.py` breaks on every `<scheme>_run`. At entry it walks the
   frame's dummy arguments, probes each array's base address and
   per-subscript byte strides empirically (element-address arithmetic, so
   strided actual args like `state%q(:ncol,:,m)` are handled), dumps raw
   bytes, and plants a `FinishBreakpoint` that re-reads the same addresses
   at scheme exit. Every hit is recorded — no timestep bookkeeping in gdb.
3. At the first hit it also captures the ordered constituent name lists:
   CAM `constituents::cnst_name`, SIMA
   `cam_constituents::const_props(i)%prop%var_std_name`.
4. `differ.py` aligns hit streams offline (CAM's first timestep is skipped;
   calls-per-step inferred per scheme), matches constituent indices via a
   short-name ↔ standard-name table (`SPECIAL_CONST_MAP` — extend it as
   ports grow), and reports in execution order:
   - inputs differ → wiring problem upstream of the scheme;
   - inputs match, outputs differ → problem inside the scheme;
   - shape/dtype/kind mismatches, hit-count mismatches, SIMA-only schemes;
   - the constituent index mapping (useful for manual gdb sessions);
   - bit-for-bit: `No differences found in any subroutines!`

Entry-time addresses of every argument are kept in
`<out>/<role>/manifest.json` — handy for setting watchpoints manually
afterwards.

## Calibration checklist (first run against a known-BFB case)

Not yet validated on a real case. Verify, in roughly this order:

1. Breakpoints resolve: `gdb.log` shows `N/N scheme entry points resolved`
   (spellings tried: `<scheme>_run`, `__<scheme>_MOD_<scheme>_run`).
2. Arg walking: each hit record in `manifest.json` lists the expected
   dummy names; arrays have sane `extents`/`strides` (contiguous arrays
   should show strides `[8, 8*ncol, ...]` for f8).
3. Dimension order: if report indices print as `(lev, col)` instead of
   `(col, lev)`, flip `DIM_ORDER` in `dropsonde_gdb.py` (cosmetic only —
   both runs use the same convention, so diffs are valid regardless).
4. Constituent capture: both manifests have a `constituents` list; the
   SIMA chase through `const_props(i)%prop%var_std_name` and
   `num_constituents` is the most fragile part (deferred-length allocatable
   chars). Failure degrades gracefully: q-arrays compare element-wise with
   a warning.
5. FinishBreakpoint exits: hit records have `"complete": true` and
   `exit_file` entries.
6. Alignment: BFB case must report zero diffs; if the first scheme's
   inputs differ the report flags alignment as suspect.

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
