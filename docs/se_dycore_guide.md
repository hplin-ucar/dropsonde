# Comparing the SE dycore between CAM and CAM-SIMA

The SE dycore is not CCPP: its state travels in derived types
(`elem(:)%state%v` ...) through subroutines that exist under the same
names, modules, and signatures in both models. dropsonde compares it with
three pieces on top of the normal workflow:

1. a **pseudo-SDF** (`docs/suite_se_dycore.xml`): one `<scheme>` per dycore
   subroutine, each annotated `dropsonde:portable=<itself>` so the
   breakpoint lands on the bare symbol and intents come from the Fortran
   source (`--meta-root` pointing at the dycore);
2. a **capture spec** (`--capture docs/capture_se_dycore.json`): maps
   derived-type names to the component paths worth comparing; matching
   dummies expand into pseudo-args (`elem%state%v`) that the differ treats
   like any array. Components at a fixed offset in the element are read as
   one strided plan; allocatable components (all of CAM-SIMA's
   `elem_state_t`, CAM's `Qdp`/`FQ`) are read per element;
3. `--step-offset 0`: both models run freely from identical (analytic)
   initial conditions, so SIMA step t pairs with CAM step t (the default 1
   is for FPHYStest snapshot runs).

```
./dropsonde --cam  /scratch/.../cam_fkessler   \
            --sima /scratch/.../sima_fkessler  \
            --sdf  docs/suite_se_dycore.xml    \
            --capture docs/capture_se_dycore.json \
            --step-offset 0 --steps 1          \
            --meta-root $CAM/src/dynamics/se
```

Both cases: DEBUG (`-O0 -g`), GNU, single rank (the capture spec path is
debug-only; the optimized raw-entry mode skips derived types).

## Before blaming the code

The dycore's control flow comes from namelist knobs held in module state
(`control_mod`), which dropsonde does not capture. Confirm parity of
`se_nsplit`, `se_rsplit`, `se_qsplit`, `se_tstep_type`, `se_ftype`,
`nu*`, `se_hypervis_subcycle*`, and the vertical grid between the two
cases first — a mismatch changes call counts and branches, and shows up in
the report's alignment table as differing hit counts rather than as a
localized answer difference.

## Reading the report

- Pseudo-args are named `dummy%path` (`elem%state%t`); the element index
  is the LAST (slowest) dimension, so a worst-diff subscript
  `(i,j,k,tl,ie)` reads GLL point (i,j), level k, timelevel tl, element ie.
- Timelevel indices (`tl%n0`, `np1`, `nm1`, and the plain integer
  `np1/n0/qn0` args) are compared as scalars; the state arrays carry every
  timelevel, so first check WHICH timelevel slice differs against the
  captured indices.
- Repeated calls (subcycles, RK stages) pair by occurrence within
  (timestep, caller); nested rows show as `via <caller>` exactly like
  nested physics calls.
- `elem%state%qdp` / `elem%derived%fq`: when the two models register
  DIFFERENT constituent counts (e.g. CAM carries test tracers SIMA does
  not), the constituent axis is compared through the report's
  constituent-mapping table, one mapped species at a time; unmatched
  species are skipped, so CAM-only tracers no longer hide the water
  species behind a shape-mismatch line. With EQUAL counts the axis is
  still compared raw (index m vs m) -- check the mapping table for
  reordering artifacts in that case.
- Entry "diffs" whose values are NaN/denormal/huge on both sides
  (omega_cn-style intent(inout) diagnostics the dycore only writes) are
  downgraded to a [note] and do NOT suppress the exit comparison; the
  exit is then compared over the elements the scheme actually wrote.
- Stack-passed args of big-frame subroutines (compute_and_apply_rhs's
  qwater/inv_cp_full/kappa/qidx, deriv%dvv) lose their gfortran -O0 DWARF
  locations and are recovered from the ABI argument slots plus the source
  declaration's bounds; a remaining "capture errored" [note] there means
  that recovery failed too -- check gdb.log.

## Blind spots

- Module state carried BETWEEN calls (`qmin`/`qmax` limiter bounds, edge
  buffers, `ur_weights`) is not captured: a divergence entering there
  surfaces one routine later than its true origin.
- Physics<->dynamics coupling: in CAM-SIMA `d_p_coupling` runs in
  `cam_timestep_init`, BEFORE the `cam_run1` sentinel (in CAM it runs
  inside it), so coupling routines would be tagged with different steps.
  Keep pseudo-SDF entries inside `dyn_run` (both models run it in
  `cam_run3`); `prim_run_subcycle` entry state already reflects the
  coupled-in forcing.
- `dyn_run` itself takes only `dyn_state` and reads `TimeLevel`/`hvcoord`
  from module scope — start the SDF at `prim_run_subcycle`, where
  everything is an argument.

## Cost

At ne16np4/L30 with the shipped capture spec, one `elem` expansion is
~100 MB; a full step (~40 hits, entry+exit) is roughly 8 GB per model.
Start with `--steps 1` and the coarse suite, then zoom: generate a finer
pseudo-SDF for the implicated region (externally — a script or an LLM can
emit `<scheme>` entries with `portable=<itself>` for that region's call
tree) and re-run. `--skip-cam`/`--skip-sima` reuse captures while
iterating on the spec.

## SE-CSLAM (physgrid, e.g. ne16pg3)

On `pg` grids CSLAM advects the tracers; use
`docs/suite_se_dycore_cslam.xml` + `docs/capture_se_dycore_cslam.json`
(GLL dynamics chain unchanged; `prim_advec_tracers_remap`/`euler_step`/
`advance_hypervis_scalar` replaced by `prim_advec_tracers_fvm` →
`run_consistent_se_cslam`; `fvm_struct` state captured). pg-specific
gotchas:

- **Hit counts**: CSLAM is supercycled — tracer advection runs every
  `fvm_supercycling` (and `fvm_supercycling_jet`) rsteps, so add those to
  the namelist-parity checklist; a mismatch shows as differing
  `prim_advec_tracers_fvm` hit counts.
- **Tracer split**: GLL `elem%state%qdp` carries only the thermodynamic
  active species (`qsize = thermodynamic_active_species_num`, CAM: from
  air_composition; CAM-SIMA: its port) while `fvm%c` carries all advected
  tracers (`ntrac`). A count mismatch between the models is compared per
  mapped species (see above) when the counts equal each model's
  registered-constituent count; a `qdp` axis that matches NEITHER model's
  registry (wrong `qsize`) still reports as a shape mismatch — that is
  itself the finding, not tool noise.
- **Halos**: `fvm%c`/`dp_fvm`/`se_flux` include halo rings (lbounds
  `1-nhc` etc., visible in the reported subscripts). Halo cells hold
  exchange leftovers; treat a worst-diff at a subscript outside `1..nc`
  with suspicion and look for the first *interior* diff.
- The fvm dummies are fixed-shape in CAM but allocatable in CAM-SIMA —
  handled transparently (same as `elem_state_t`), shapes still compare.
- Zooming further: the stages inside `run_consistent_se_cslam`
  (reconstruction, swept-area remap in `fvm_consistent_se_cslam.F90`)
  can be added as pseudo-SDF entries the same way if the divergence lands
  there.

## Old-gdb sites

gdb <= 8.x (e.g. Izumi's 8.2) cannot read the DWARF5 that gcc >= 11
emits by default: dummy expansion then probes garbage addresses. Build the
model (and calibrations) with `-gdwarf-4` there. Derecho's gdb 16 needs
nothing. Calibrate a machine with `tests/cal_se_run.sh` (gfortran needs
`FCFLAGS="-g -gdwarf-4 -O0"` on such sites; ifort's default DWARF is fine).
