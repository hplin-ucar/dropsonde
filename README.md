# dropsonde

![a lovely squid representing the SIMA project flying with a dropsonde](dropsonde_logo.png)

Find the first CCPP scheme where CAM and CAM-SIMA physics answers diverge,
without modifying or recompiling either model. Runs both executables under
gdb with breakpoints at every scheme `_run` entry/exit listed in the SDF,
dumps all scheme arguments to disk, and diffs them byte-for-byte.

## Usage

1. Build the models: `DEBUG=TRUE,NTASKS=1,NTHRDS=1`, chunks disabled for CAM.
   **Important:** CAM and CAM-SIMA must be running the same underlying atmos_phys code.
2. Run `./preview_namelists`.
3. Run `dropsonde` where it is checked out: point it to the case directories on scratch with the SDF file:

```
./dropsonde --cam  /glade/derecho/scratch/.../cam_case_dir  \
            --sima /glade/derecho/scratch/.../sima_case_dir \
            --sdf  suite_park_macrop.xml \
            [--meta-root /path/to/atmospheric_physics] [--steps 1]
```

`--steps` defaults to 1: in FPHYStest snapshot runs each timestep is
independent, so a single matching pair is decisive. Use `--steps 2` to
also cover `is_first_timestep` logic.

`--meta-root` (default: three levels up from the SDF) locates the schemes'
CCPP `.meta` files; their intents let the differ skip `intent(out)` args
at entry (caller-side garbage) and `intent(in)` args at exit.

Run it where the models can run (i.e., a compute node). Both gdb sessions
launch in parallel; the report prints to stdout when they finish.

Everything is written out in `--out` (default `./dropsonde_out`):

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
2. `dropsonde_gdb.py` breaks on every `<scheme>_run`. gdb's prologue skip
   stops inside gfortran's dummy-argument descriptor setup, where the
   dummies aren't readable yet and single-stepping doesn't reliably escape.
   So it instead `advance`s to the first body statement (`errmsg=''`/
   `errflg=0`, located by scanning the source from the subroutine's
   definition line), where every dummy is live and no input has been
   touched. It then walks the frame's dummy arguments, probes each array's
   base address and per-subscript byte strides empirically (element-address
   arithmetic, so strided actual args like `state%q(:ncol,:,m)` are
   handled), dumps raw bytes, and plants a `FinishBreakpoint` that re-reads
   the same addresses at scheme exit. Dummies that gfortran's `-O0` debug
   info fails to locate at all — the stack-passed arguments of very large
   subroutines — are recovered directly from the System V AMD64 argument
   slots instead. Every hit is tagged with the current timestep (sentinel
   breakpoint on `cam_run1`) and its caller, and the run is killed after
   enough timesteps.
3. At the first hit it also captures the ordered constituent name lists:
   CAM `constituents::cnst_name`, SIMA
   `cam_constituents::const_props(i)%prop%var_std_name`.
4. `differ.py` aligns offline: SIMA timestep t pairs with CAM timestep
   t+1 (CAM's first dumped step is skipped). Hits are bucketed by
   caller — nested calls (e.g. diffusion inside `bretherton_park`) share
   the same caller symbol in both binaries and are compared against their
   counterparts; toplevel calls (CCPP cap vs CAM driver) bucket together.
   When hit counts differ within a bucket, SIMA hits are paired with the
   CAM hit whose comparable inputs all match bitwise, or reported unpaired
   with the closest candidate. Constituent-indexed arrays are compared per
   mapped species via `SPECIAL_CONST_MAP`, then exact-name match, then
   `cnst_<short name>`. Intent filtering (from `.meta`) drops
   entry-garbage and input-repeated-at-exit noise;
   `errmsg`/`errflg`/`iulog` are never compared.

gdb/gfortran assumptions (gdb 16.2, gfortran 12, Derecho):

- Array dims are REVERSED in the gdb Python API.
- `set breakpoint pending off` is ignored by the Python API; detect via
  `bp.pending`.
- gdb's prologue skip lands the scheme breakpoint inside the argument-
  descriptor setup, where dummies raise "Location address is not set" and
  single-stepping doesn't reliably escape. `advance` to the body's first
  statement instead (found by scanning the source for `errmsg=''`/`errflg=0`
  from the subroutine's definition line); dummies are live there, with
  inputs untouched.
- Very large subroutines (`park_macrophysics_run`: 66 dummies, ~1500 lines,
  a ~30 KB automatic-array frame) defeat gfortran `-O0` debug info outright:
  it emits an *empty* `DW_AT_location` for every stack-passed dummy (the
  first 6, register-passed, still resolve), so gdb can't read them by name
  anywhere in the routine — not a stopping-point problem. The capture falls
  back to the System V AMD64 ABI: dummy `i` (declaration order, `i>6`) is the
  pointer at `[rbp + 16 + 8*(i-7)]` — a gfortran array descriptor for
  assumed-shape args, a value pointer for scalars. The descriptor (`base_addr`
  `+0`, `elem_len` `+16`, `rank` `+28`, `dim[]` `+40`) yields the same base/
  extents/strides the normal path probes. Assumes all-by-reference dummies;
  by-value float args (which shift the integer-register accounting) aren't
  handled.
- Deferred-length character components: length is in a hidden
  `_<name>_length` member, not the DWARF type.
- Use `module::var` spelling for globals.

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
    absent on one side, and `intent(out)` elements never written by
    either model (exit == entry bitwise, still holding each caller's
    pre-call memory). Partially-filled outputs or oversized actual
    arguments show up this way.
- **per-scheme summary** (on failure) — diff counts per scheme, schemes
  that were bit-for-bit.
- **alignment warning + offset scan** — printed when the first compared
  scheme already has input diffs. Scans all dumped steps on both sides
  to detect a wrong step offset or a re-read snapshot record. Distrust
  everything below the first scheme until resolved.

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
resolution, caller capture, extent/stride probing, and constituent-name
capture.

Full-model calibration (2026-06, CAM5 UW PBL suite
`suite_vdiff_bretherton_park`, ne3, FHIST_C5 vs FPHYStest): all
end-to-end mechanics verified. With clean snapshots the suite was
bit-for-bit except for real, separately-confirmed SIMA-side wiring
issues the tool surfaced — divergences it reports are real, not tool
noise. A fully-clean model pair has not yet been tested.

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

## Credits

Dropsonde was created in a session with Claude Mythos/Fable 5 in a discussion
about programmatically operating gdb to automate the discovery of answer
differences between CAM and CAM-SIMA.

Later improvements were added by Claude Opus 4.6.

The supervising human was Haipeng Lin <hplin@ucar.edu>, NSF NCAR/CGD/AMP.
