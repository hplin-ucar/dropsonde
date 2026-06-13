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

Everything lands in `--out` (default `./dropsonde_out`):

```
report.txt    the printed report, archived verbatim
suite.json    scheme order, --steps, parsed .meta intents
cam/, sima/   manifest.json + gdb.log + one .bin dump per argument
              per phase (entry/exit)
```

Re-run just the comparison — to re-read the report, regenerate it after
editing `differ.py`, or compare with different code without re-running
the models — with:

```
python3 differ.py dropsonde_out
```

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
   the diffusion solver inside `bretherton_park`'s iteration) has the
   *same shared-code caller symbol* in both binaries and is compared
   against its nested counterpart; model-specific call sites (CCPP cap vs
   CAM driver) bucket together as toplevel. When hit counts still differ
   within a bucket (CAM hosts call `geopotential_temp` from ~20 sites per
   step), each SIMA hit is paired with the CAM hit whose comparable
   inputs ALL match bitwise, or reported unpaired with the closest
   candidate. Constituent-indexed arrays — rank-3 `q(ncol,pver,pcnst)`,
   rank-2 fluxes `cflx(ncol,pcnst)`, rank-1 per-species config like
   `qmincg(pcnst)` — are compared per mapped species via a short-name ↔
   standard-name table (`SPECIAL_CONST_MAP`, then exact name, then
   `cnst_<short name>`: constituents auto-registered from a snapshot file
   keep their netcdf variable name as their standard name). Intent
   filtering (from `.meta`) drops entry-garbage and
   input-repeated-at-exit noise; `errmsg`/`errflg`/`iulog` are never
   compared.

gdb/gfortran facts the capture relies on (gdb 16.2, gfortran 12, Derecho):

- Array dims appear REVERSED through the gdb Python API; a per-array
  sizeof cross-check guards the convention.
- The Python API ignores `set breakpoint pending off`; unresolved
  breakpoints must be detected via `bp.pending`.
- gdb's post-prologue stop is before explicit-shape dummy bounds exist;
  one `next` materializes them.
- Deferred-length character components print as
  `(PTR TO -> character*0)` (DWARF carries no dynamic length); the
  length lives in a hidden member `_<name>_length` that gfortran appends
  after each type extension level's visible components, readable like
  any other component.
- Prefer `module::var` spelling for globals.

## Reading the report

In order of appearance:

- **schemes compared** — SDF schemes whose `_run` symbol resolved in both
  binaries. SIMA-only schemes (no CAM symbol — usually CCPP-only glue
  like `*_stub`, `*_default`) are listed and skipped; a CAM-only entry
  would mean the SIMA build is missing a scheme.
- **constituent mapping** — the recovered `q(:,:,i) <-> q(:,:,j)` index
  permutation with both names. Unmatched species are listed and excluded
  from per-species comparison (extend `SPECIAL_CONST_MAP` in `differ.py`
  if a name pair should match but doesn't). Also handy as a lookup table
  during manual gdb sessions.
- **call alignment** — one row per call site with per-side hit counts
  (`sima x4 cam x4`); indented rows are calls made from inside the
  parent scheme (e.g. the diffusion routines inside `bretherton_park`'s
  iteration). `<--` markers flag count mismatches and say how they were
  handled.
- **comparison (execution order)** — the heart of the report. The first
  divergence is boxed. Each line is one argument of one hit pair:
  - `[INPUTS DIFFER]`: the two models handed the scheme different data —
    the problem is *upstream* (wiring, initialization, an earlier
    scheme), not in the scheme itself.
  - `[OUTPUTS DIFFER]`: identical inputs, different results — the
    problem is *inside* the scheme (or a build difference). Exit
    comparison is suppressed for args whose entry already differed.
  - Constituent-indexed args report per species:
    `NUMLIQ (cam idx 4, sima idx 1) <-> mass_number_... : ...` with
    element stats, or just `sima=… cam=…` for rank-1 (one value per
    species).
  - `[note]` lines are non-verdicts: hits paired by input match, args
    absent on one side, and `intent(out)` args where every element the
    scheme wrote matches but some elements were never written in either
    model (exit == entry bitwise) — those still hold each caller's
    pre-call memory and are not compared. A scheme that fills its
    outputs only partially (e.g. `dtk` rows above the surface layer)
    shows up this way; persistent notes of this kind can also point at
    an oversized actual argument on one side.
- **per-scheme summary** (on failure) — diff counts per scheme, schemes
  that were bit-for-bit.
- **alignment warning + offset scan** — printed when the very *first*
  compared scheme already has input diffs. The first hit's entry args
  are then matched bitwise against every dumped CAM step (some other
  step matches → step offset is wrong) and the other SIMA steps
  (identical → the snapshot time record is being re-read). Distrust
  everything below the first scheme until this is resolved.

## Digging deeper after a run

The comparison is rerunnable and the raw data self-describing; nothing
about a finding requires re-running the models.

- `<out>/<role>/manifest.json` — per hit: `scheme`, `step`, `caller`,
  `hit`, and per argument: `kind`, `dtype` (`f8`, `i4`, ...), `extents`
  and `los` (Fortran order, first subscript fastest; `los` = lower
  bounds), `addr` (entry-time base address), `entry_file`/`exit_file`
  (dump names), or `entry_value`/`exit_value` for scalars. Skipped args
  carry a `why`.
- Dumps are named `NNNNN_<scheme>_h<hit>_<arg>_<in|out>.bin`: raw
  native-endian bytes in Fortran memory order, no header. Bitwise
  questions need no unpacking: `cmp a.bin b.bin`. To look at values:

  ```python
  from array import array
  a = array("d")            # "d" for dtype f8, "i" for i4, ...
  a.frombytes(open("00123_..._q1_out.bin", "rb").read())
  # element (i,j,k), extents [ni,nj,nk], lower bounds los:
  # a[(i-los[0]) + ni*((j-los[1]) + nj*(k-los[2]))]
  ```

- Argument storage is caller-owned, so the manifest `addr` of an
  argument is valid for the whole scheme call: break at the scheme in a
  manual gdb session and `watch *(double *)0x<addr>` to catch the exact
  statement that writes a differing element.
- `<out>/<role>/gdb.log` — full gdb transcript (breakpoint resolution,
  capture notes, the kill).

## Validation

Micro-calibration: `tests/cal.f90` + `tests/cal_run.sh` (any machine
with gdb + gfortran; Derecho login nodes work) checks breakpoint
resolution, caller capture, extent/stride probing against known arrays,
and constituent-name capture.

Full-model calibration (2026-06, CAM5 UW PBL suite
`suite_vdiff_bretherton_park`, ne3 = 486 columns, CAM FHIST_C5 vs
CAM-SIMA FPHYStest): end-to-end mechanics verified — step tagging and
auto-kill, explicit-shape capture, exits, nested-call bucketing,
constituent capture on both sides, pending-breakpoint classification,
best-match hit pairing, and the offset scan (which correctly identified
a snapshot from a mismatched CAM trajectory before clean snapshots were
regenerated). With clean snapshots the decisive timestep pair was
bit-for-bit across the whole suite except real, separately-confirmed
SIMA-side issues the tool surfaced (constituents silently reading as
zero when `ic_file_input_names` is missing; a `diffuse_tracers_init`
constituent flag mismatch; `qmin` passed where CAM passes `qmincg`) —
divergences it reports are real wiring differences, not tool noise.

Not yet exercised: a fully-fixed model pair reporting
`No differences found in any subroutines!`.

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
