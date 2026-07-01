# dropsonde

![a lovely squid representing the SIMA project flying with a dropsonde](dropsonde_logo.png)

Find the first CCPP scheme where CAM and CAM-SIMA physics answers diverge,
without needing to modify either model to add instrumentation or manually
operating the debugger.

`dropsonde` runs both models under `gdb` with breakpoints at every scheme's
`run` phase (as specified in the SDF) and dumps all arguments entering/exiting
these subroutines to disk, and diffs them byte-for-byte.

## When to use `dropsonde`?

`dropsonde` can be used during the final stretch of converting CAM physics to CCPP,
where CAM already runs the CCPPized subroutines, but despite the same code,
CAM-SIMA driven with the CAM snapshots have answer differences compared to CAM.

## Usage

1. Build the models: `DEBUG=TRUE,NTASKS=1,NTHRDS=1`, with chunks disabled for CAM.
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
suite.json    scheme order, --steps, parsed intents, portable overrides
realign.json  copy of the --realign spec, if one was given
cam/, sima/   manifest.json + gdb.log + one .bin dump per argument
              per phase (entry/exit)
```

Re-run just the comparison: to re-read the report, regenerate it after
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
- CAM: GNU or Intel, `DEBUG=TRUE`, `NTASKS=1`, `NTHRDS=1`, dechunked
  (`pcols` = total columns). Must call the same CCPPized
  `atmospheric_physics` routines as the SIMA build.
- CAM-SIMA: GNU or Intel, `DEBUG=TRUE`, single rank. Both runs must use
  the same compiler (caller symbols and constituent capture are matched
  across the pair). The compiler is auto-detected from the DWARF producer
  string; no flag needed. NAG is not yet supported (see Known limitations).
- `STOP_N` only needs to be *long enough* (CAM's 6-step coupling minimum
  is fine): each gdb session counts timesteps via a sentinel breakpoint on
  `cam_run1` and kills its run once it has what it needs (CAM:
  `--steps + 1` timesteps, SIMA: `--steps`).
- `gdb` on PATH. No numpy needed anywhere — gdb side and differ are
  stdlib-only, Python >= 3.6 (works with gdb's embedded interpreter and
  `/usr/bin/python3` on Derecho). Calibrated on gdb 16.2 (Derecho) and
  gdb 8.2 (Izumi). dropsonde strips `PYTHONHOME`/`PYTHONPATH` when it
  launches gdb, so a site anaconda on those vars (e.g. Izumi) does not
  crash gdb's embedded interpreter with `No module named 'encodings'`.

## How it works

1. `dropsonde` parses the SDF for the ordered scheme list and launches
   `gdb --batch -x dropsonde_gdb.py` on each executable (config passed via
   `$DROPSONDE_CONFIG`).
2. `dropsonde_gdb.py` breaks on every `<scheme>_run` (or, for schemes with a
   `dropsonde:portable=` SDF annotation, on the named portable subroutine
   instead — see [Comparing shared portable subroutines](#comparing-shared-portable-subroutines)). gdb's prologue skip
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

The compiler is auto-detected at the first scheme hit from the DWARF
producer string (`symtab.producer`: "GNU Fortran ..." / "Intel(R)
Fortran ...") and selects the memory-layout-specific paths below; the
by-name element-address probing is identical for both.

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

Intel assumptions (gdb 8.2, ifort 19.1, Izumi):

- Symbols are mangled `module_mp_<proc>_`; gdb demangles them to
  `module::<proc>`, which the bare name does not resolve to. Breakpoints
  fall back to the qualified spelling found by listing the symbol table
  (`info functions \b<proc>\b` → `module::<proc>`) — module-agnostic, so a
  scheme whose `_run` lives in a differently named module still resolves.
- Dummies are live right at the entry breakpoint (no descriptor-setup
  prologue to skip); the body `advance` is harmless and still used.
- gdb 8.2 misreads Intel's array descriptor for the module pointer-array
  `cam_constituents::const_props` (bogus bounds), so the SIMA constituent
  names are decoded from raw memory: descriptor `base_addr` `+0`,
  `elem_len` `+8`; each element's `%prop` points to a properties type whose
  `var_std_name` (`character(len=:),allocatable`) is laid out inline as
  `{char *data; size_t len}` (data at the component offset, length `+8`).
- The System V AMD64 ABI fallback is gfortran-descriptor-specific and so
  gfortran-only; under Intel it bails with a note rather than misreading a
  descriptor (Intel has not been observed to drop dummy locations).

## Comparing shared portable subroutines

Some CCPP schemes wrap *portable* science code that CAM and CAM-SIMA share
verbatim but reach through different wrappers — e.g. `modal_aero_calcsize_ccpp`
(SIMA) and `modal_aero_calcsize_sub` (CAM) both call the portable
`modal_aero_calcsize_run`. CAM never calls the `<scheme>_run` CCPP wrapper, so
the wrapper symbol is absent from the CAM binary and the scheme would be
reported "SIMA-only (no CAM symbol)". The real comparison point is the shared
portable subroutine, which has an identical signature in both binaries.

Retarget such a scheme by annotating its `<scheme>` line in the SDF with a
**same-line XML comment** naming the portable subroutine:

```xml
<scheme>modal_aero_calcsize_ccpp</scheme>  <!-- dropsonde:portable=modal_aero_calcsize_run -->
```

The comment is invisible to the CCPP framework (no schema change). dropsonde
then plants the breakpoint on the portable subroutine in *both* binaries
(internally still keyed by the SDF scheme name, so pairing and suite ordering
are unchanged) and sources arg intents from the portable `.F90` declarations
(portable subroutines have no `.meta`). The report shows such schemes as
`scheme -> portable_sub` (e.g. `modal_aero_calcsize_ccpp -> modal_aero_calcsize_run`)
so it is clear what is actually being compared. Non-annotated schemes are
unaffected and keep using `<scheme>_run` and `.meta` intents.

Notes and limits: this fits *driver-level* portable subroutines called once
per (dechunked) timestep, not per-column kernels called in a loop. Derived-type
dummy args (e.g. `class(aerosol_state)`) are recorded but not captured (gdb
can't read them portably). Wrapper-only post-processing (e.g. mapping
tendencies into the constituent array) is *not* compared at the portable
boundary — only the shared science inputs/outputs are.

### Realigning constituent index conventions (`--realign`)

Some portable subroutines are additionally called from two different
*constituent index conventions*. The MAM aerosol routines
(`modal_aero_gasaerexch`/`_rename`) are the archetype: CAM passes the mozart
`vmr`/`vmrcw` arrays (`gas_pcnst` wide, `loffset = imozart-1`, mode/species
pointer arrays indexed in full `pcnst` space) while CAM-SIMA passes one packed
constituent array (`loffset = 0`). The routine reaches every constituent as
`array(:,:, ptr - loffset)`, so the constituent axis of these args cannot be
compared byte-for-byte and would otherwise drown the report in false
shape/value diffs.

`--realign spec.json` teaches the differ to realign such args per species.
The spec (see `docs/realign_mam.json`, which covers MAM — a CARMA spec would
follow the same shape) names:

- `offset_arg` — the per-model index-offset dummy (`loffset`);
- `spaces` — groups of constituent-indexed args and how species pair across
  the models: `constituent-map` (shift the name-matched constituent map by
  each model's offset) or `pointer-arrays` (pair via captured mode/species
  pointer arrays such as `lmassptrcw_amode`/`numptrcw_amode`, for species
  with no registered constituent name — e.g. MAM cloud-borne, which CAM
  keeps in a separate unregistered `vmrcw` array);
- `convention_args` — index-space/sentinel args expected to differ between
  the models (reported as a note, never counted as divergences);
- `metadata_args` + `valid_slot_mask` — physical-constant tables compared
  over real (species,mode) slots only, ignoring per-model padding fill.

The spec is validated before the models launch and copied to
`<out>/realign.json`, where the differ reads it — edit that copy and re-run
`python3 differ.py <out>` to iterate without re-running the models (drop a
spec there by hand to reprocess an existing dump). Realignment applies only
to hits of portable-annotated schemes that carry the offset arg on both
sides; everything else compares exactly as before.

## Reading the report

In order of appearance:

- **schemes compared** — SDF schemes whose `_run` symbol resolved in both
  binaries. SIMA-only schemes (no CAM symbol — usually CCPP-only glue
  like `*_stub`, `*_default`) are listed and skipped; a CAM-only entry
  would mean the SIMA build is missing a scheme. A scheme is also SIMA-only
  when CAM and CAM-SIMA share the *science* but call it through different
  wrappers — CAM never calls the `<scheme>_run` CCPP wrapper, so it has no
  matching symbol. See [Comparing shared portable subroutines](#comparing-shared-portable-subroutines)
  to retarget such a scheme to the portable subroutine both models call.
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
with gdb + a Fortran compiler; Derecho login nodes work) checks breakpoint
resolution, caller capture, extent/stride probing, and constituent-name
capture. The compiler is selectable via `$FC` (default `gfortran`):
`FC=ifort bash cal_run.sh`. Passes clean under gfortran (Derecho, gdb
16.2) and ifort 19.1 (Izumi, gdb 8.2). Under gfortran on gdb 8.2 the SIMA
constituent names read as `<unreadable>` (that gdb's gfortran deferred-
length-character support is too old); the gdb 16.2 gfortran path is
unaffected.

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
- NAG unsupported for now: it lowers Fortran to C (DWARF producer is
  `GNU C...`), mangles names `module_MP_<proc>` with C-mangled,
  pointer-typed dummies (`a_Dummy`, `ncol_`), and the Fortran name `a` is
  gone, so by-name element addressing can't be used. Assumed-shape dummies
  are still recoverable — gdb reads NAG's `__NAGf90_Dope2` descriptor as a
  struct (`addr`, `offset`, `dim[].{extent,mult,lower}`) — but explicit-
  shape dummies (`b(ncol,pver)`) lower to a bare data pointer with no
  shape, which would have to come from the `.meta` files. A dedicated NAG
  backend (arg demangling + dope-vector decode + `.meta`-sourced shapes)
  is future work. (Intel on Izumi, by contrast, is supported.)
- Negative-stride (reversed) array slices are skipped with a note.
- One subcycle scheme called with *different arguments* per sub-step
  aligns fine (hits are sequence-matched), but the report labels hits by
  per-scheme index, not subcycle position.

## Credits

Dropsonde was created in a session with Claude Mythos/Fable 5 in a discussion
about programmatically operating gdb to automate the discovery of answer
differences between CAM and CAM-SIMA.

Later improvements were added by Claude Opus 4.6 and Claude Opus 4.8.

The supervising human was Haipeng Lin <hplin@ucar.edu>, NSF NCAR/CGD/AMP.

(c) 2026 University Corporation for Atmospheric Research.
