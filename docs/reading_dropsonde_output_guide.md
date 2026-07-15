# Reading dropsonde output: a field guide

To read the `dropsonde` output, we must understand what `dropsonde` does.
In short:
* `dropsonde` runs CAM and CAM-SIMA copies that point to the same underlying
  `atmospheric_physics` code (i.e., CAM must already be calling the same
  CCPPized schemes under a CAM interface) under the `gdb` debugger;
* It looks at every subroutine in the suite definition file (SDF) and captures
  then matches their input and output arguments in CAM and CAM-SIMA.

This saves hours of debugging time bisecting differences, because the output
will tell you exactly when the two runs have started diverging from each other.

The primary output of `dropsonde` is the report. It is meant to be read from
top to bottom. It contains five sections:
1. the schemes being compared: because the SIMA SDF may contain SIMA-specific
   interstitials, this section will tell you which ones are being compared
   and which ones were not found in CAM.
2. constituent mapping, because SIMA has no guaranteed ordering of constituents.
3. call alignment: how dropsonde analyzed the CAM and SIMA runtimes and aligned
   them with each other.
4. **comparison**: the heart of the report, which tells you which arguments
   do not match.
5. per-scheme summary: which schemes differed and which were bit-for-bit?

The three `example_*.log` files in this folder are the running examples. Open
one alongside this guide.

Note: you want to look at the **first** (non-roundoff) divergence, because
everything below it often is cascading differences downstream. The workflow with
`dropsonde` is intended to be:
* run `dropsonde`
* find the first divergence, triage it, and fix it in code
* rinse and repeat

## Going through one report, one section at a time

**1. schemes compared.** Which `_run` symbols resolved in both binaries.
It is expected that SIMA-only schemes will not match (of course), and schemes
that do not have a `_run` phase will also not match.

```
schemes: 21 compared
  SIMA-only (no CAM symbol, not compared): initialize_constituents, ...
```

**2. constituent mapping.** The recovered `q(:,:,i) <-> q(:,:,j)` permutation
between the two models' tracer arrays.

**Don't skip this easily!** If you are missing constituents on either side, it
may mean a standard name mismatch (see example below).

**3. call alignment.** One row per call site, with hit counts per side.
Indented rows are calls made from *inside* a parent scheme. `<--` markers flag
anything unusual and tell you how it was handled.

This should roughly resemble your suite definition file's order, and it tells
you exactly how many times SIMA and CAM called each run phase and how `dropsonde`
resolved the alignment of these calls. If the counts look off, there may be a loop
issue, or CAM is running with chunks (don't! make sure chunking is off by setting pcols
in `CAM_CONFIG_OPTS`).

```
call alignment (hits per compared step; sima step t pairs with cam step t+1;
indented rows are calls made from inside the parent scheme):
  beljaars_zero_stub                                     sima x1   cam x0    <-- no CAM hits: not compared
  vertical_diffusion_prepare_inputs                      sima x1   cam x1
  hb_diff_prepare_vertical_diffusion_inputs              sima x1   cam x1
  vertical_diffusion_set_temperature_at_toa_default      sima x1   cam x0    <-- no CAM hits: not compared
  vertical_diffusion_interpolate_to_interfaces           sima x1   cam x1
  vertical_diffusion_set_total_surface_stress            sima x1   cam x0    <-- no CAM hits: not compared
  bretherton_park_diff                                   sima x1   cam x1
      implicit_surface_stress_add_drag_coefficient       sima x1   cam x1
      turbulent_mountain_stress_add_drag_coefficient     sima x1   cam x1
      vertical_diffusion_wind_damping_rate               sima x1   cam x1
      vertical_diffusion_diffuse_horizontal_momentum     sima x4   cam x4
      vertical_diffusion_diffuse_dry_static_energy       sima x4   cam x4
      vertical_diffusion_diffuse_tracers                 sima x4   cam x4
  compute_kinematic_fluxes_and_obklen                    sima x1   cam x1
  eddy_diffusivity_adjustment_above_pbl                  sima x1   cam x1
  vertical_diffusion_sponge_layer                        sima x1   cam x1
  implicit_surface_stress_add_drag_coefficient           sima x1   cam x1
  turbulent_mountain_stress_add_drag_coefficient         sima x1   cam x1
  vertical_diffusion_wind_damping_rate                   sima x1   cam x1
  vertical_diffusion_diffuse_horizontal_momentum         sima x1   cam x1
  vertical_diffusion_set_dry_static_energy_at_toa_zero   sima x1   cam x0    <-- no CAM hits: not compared
  vertical_diffusion_diffuse_dry_static_energy           sima x1   cam x1
  vertical_diffusion_diffuse_tracers                     sima x1   cam x1
  dropmixnuc_apply_surface_fluxes                        sima x1   cam x0    <-- no CAM hits: not compared
  vertical_diffusion_tendencies                          sima x1   cam x1
  geopotential_temp                                      sima x1   cam x20   <-- hit counts differ: paired by bitwise input match (see notes)
  update_dry_static_energy                               sima x1   cam x1
```

Note: `geopotential_temp` having much more counts in CAM is expected, because
CAM is a full run (i.e., includes all other physics), but SIMA is pure test SDF,
so it is totally expected that CAM will call it many more times. `dropsonde` will
match the one that naturally follows the parameterization being debugged.

**4. comparison (execution order).** The heart of the report. One line per
argument per hit pair. The first divergence is boxed:

```
FIRST DIVERGENCE -----------------------------------------------
  [INPUTS DIFFER] vertical_diffusion_diffuse_tracers [step 1] (hit 0) arg do_diffusion_const: constituent-indexed array, 1/25 mapped species differ: NUMLIQ
    NUMLIQ (cam idx 4, sima idx 12) <-> mass_number_concentration_of_cloud_liquid...: sima=1 cam=0
----------------------------------------------------------------
```

The example above is extremely high quality signal, as it tells you that:
* the input argument `do_diffusion_const` to `vertical_diffusion_diffuse_tracers` diverged;
* that argument is an array with a constituent dimension,
* the divergent entry corresponds to the index for `NUMLIQ`
  (yes, `dropsonde` will map CAM to SIMA! it compared SIMA index 12 to CAM index 4 to find the difference)
* the divergent value is that SIMA has `1` and CAM has `0`.

But it may be possible that the first difference is roundoff (`1e-14`). In that case, you have to look further.

**5. per-scheme summary.** On failure, a tally of which schemes diffed and which
were bit-for-bit. Good for a quick "how bad is it."

```
per-scheme: 18/21 compared schemes bit-for-bit; differing:
  vertical_diffusion_diffuse_tracers: 2 input / 1 output diffs
```

**Now you know how to read `dropsonde` outputs!**

Below is some additional guidance:

## Generally, the first divergence is the lead

`[INPUTS DIFFER]` vs `[OUTPUTS DIFFER]` is the single most important
distinction in the whole report:

- **`[INPUTS DIFFER]`** — the two models handed the scheme *different data*. The
  bug is **upstream**: wiring, initialization, or an earlier scheme. Not here.
- **`[OUTPUTS DIFFER]`** — identical inputs, different results. The bug is
  **inside this scheme** (or a build difference).

Example 1 is a clean illustration. The first divergence is an input flag:

```
[INPUTS DIFFER] ... arg do_diffusion_const: 1/25 mapped species differ: NUMLIQ
    NUMLIQ ...: sima=1 cam=0
```

SIMA is applying diffusion to NUMLIQ where CAM isn't. The effects of this
cascade downstream:

```
[OUTPUTS DIFFER] ... arg q1: 1/25 mapped species differ: NUMLIQ        <- the flag changed the result
[INPUTS DIFFER]  vertical_diffusion_tendencies ... arg q1: ... NUMLIQ   <- next scheme inherits the bad q1
[OUTPUTS DIFFER] vertical_diffusion_tendencies ... arg tend_q: ... NUMLIQ
```

The same root cause controlled by this per-constituent flag shows up in four places.
Once the root cause is fixed, all four disappear.

## Signal vs noise: a catalog

Not every line is a finding that will help you with your answer differences,
but it may uncover a different bug.

**`[note]` lines are just that, notes**. For example:
```
[note] ... arg dtk: every element the scheme wrote matches; 486 elements were never written in either model ... so they are not compared
```
this uncovers that the `dtk` array was sized one vertical layer too big
(486 is the number of horizontal columns, so it works out to exactly one layer),
but as the note says, it "(was) never written in either model". But it is a real bug to fix.

**Near-roundoff values are real but rarely the culprit.** Example 1's `qmincg`
differs for 18 species, but their magnitudes are tiny:

```
[INPUTS DIFFER] ... arg qmincg: 18/25 mapped species differ: Q, H2O2, ...
    Q ...: sima=9.9999999999999998e-13 cam=0
    H2O2 ...: sima=0 cam=9.9999999999999994e-37
```

These are constituent minimum thresholds: `1e-37` vs `0`, `1e-13` vs `0`. A
genuine config difference, worth a note, but it isn't what's moving the answer
(none of those species' *outputs* diverge).

**Per-species counts tell you the blast radius.** `1/25 mapped species differ`
is a precise, localized problem (one tracer). `18/25` or `28188/28188 elements
differ` is broad — and broad-everything, especially starting at the very first
scheme, points at wiring rather than physics (next section).

## Another example: input differences

Looking at example 2 (PUMAS):

```
FIRST DIVERGENCE -----------------------------------------------
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_numliq: 28188/28188 elements differ ...: sima=0 cam=64978883.8...
----------------------------------------------------------------
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_numice: 28188/28188 elements differ, max |diff| 5.955e+05 at (201,28): sima=0 cam=595530.6096265011
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_numrain: 28188/28188 elements differ, max |diff| 1.115e+06 at (70,54): sima=0 cam=1114690.4985518879
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_numsnow: 28188/28188 elements differ, max |diff| 3.389e+03 at (415,55): sima=0 cam=3388.657623429518
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_numgraup: 28188/28188 elements differ, max |diff| 4.427e+02 at (139,28): sima=0 cam=442.67857812209826
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_strat_liq_cldfrc: 13464/28188 elements differ, max |diff| 9.990e-01 at (199,12): sima=0 cam=0.999
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_strat_ice_cldfrc: 1768/28188 elements differ, max |diff| 1.000e+00 at (70,57): sima=0 cam=1
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_qsatfac: 28188/28188 elements differ, max |diff| 1.000e+00 at (1,1): sima=0 cam=1
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_naai: 812/28188 elements differ, max |diff| 1.031e+02 at (289,29): sima=0 cam=103.07385026541536
  [INPUTS DIFFER] micro_pumas_ccpp [step 1] (hit 0) arg pumas_npccn: 28188/28188 elements differ, max |diff| 3.610e+04 at (208,36): sima=0 cam=-36099.379900846732
...
```

In this example, the run phase to `micro_pumas_ccpp` has diverged on its inputs: all of the input fields
`pumas_numliq`, `pumas_numice`, ... are *completely* different between CAM and SIMA. And note
that CAM has real values (e.g., `sima=0 cam=64978883.8`) but SIMA has zero. This could either mean
* SIMA is not passing the values into the scheme correctly
* SIMA is passing correctly, but not read from snapshot: check `atm.log` at `debug_output = 2`
* SIMA is reading the snapshot, but something is overwriting it (check interstitials)

When the *first* compared scheme already has input diffs like this, the report
ends with an `alignment check:` verdict. `dropsonde` bitwise-matches that
scheme's entry args against every dumped CAM step and states the conclusion:
if most args match the paired step, alignment is correct and the differing
inputs are field-specific upstream issues; it only warns when another step
matches better (wrong step offset) or no step matches well (wrong snapshot).
(The example logs in this folder predate this and show the older
`*** WARNING` + raw offset scan output.)

**`nan` on one side** is uninitialized / not-provided-in-snapshot memory:

```
[INPUTS DIFFER] ... arg pumas_frzimm: ... sima=0 cam=nan
[INPUTS DIFFER] ... arg pumas_effi_external: ... sima=25 cam=nan
```

This means that CAM might have uninitialized pbuf fields when it ran.
This is entirely possible, e.g., QPC (aquaplanet) compsets lack aerosols, generally.
The report prints a `[hint]` after such args; the decisive check is whether the
scheme's *outputs* still match — if they do, the nan inputs were never read and
the whole block is a red herring. (If a run crashed inside the scheme, a
`[note]` will say `no exit capture ... outputs NOT compared` instead — then you
cannot draw that conclusion.)

## Don't skip the constituent map

The mapping table catches problems before you ever reach the diffs. Compare the
two PUMAS runs.

Example 2 is clean — 11 species each, every one matched:

```
constituent mapping (CAM 11 species, SIMA 11):
  cam q(:,:, 4) NUMLIQ   <-> sima q(:,:, 6) mass_number_concentration_of_cloud_liquid_wrt_moist_air_and_condensed_water
```

Example 3 is the *same model pair* but the count is off — 11 vs **13** — with
two unmatched leftovers:

```
constituent mapping (CAM 11 species, SIMA 13):
  cam q(:,:, 4) NUMLIQ   <-> sima q(:,:, 7) mass_number_concentration_of_cloud_liquid_wrt_moist_air_and_condensed_water
  ...
  sima q(:,:, 5) mass_number_concentration_of_cloud_liquid_water_wrt_moist_air_and_condensed_water (unmatched, not compared)
  sima q(:,:,12) mass_number_concentration_of_cloud_ice_wrt_moist_air_and_condensed_water (unmatched, not compared)
```

Look closely at the names: SIMA registered a `cloud_liquid_water` *and* a
`cloud_liquid`, an `ice` *and* a `cloud_ice` — near-duplicate standard names for
the same physical quantity. That's a registry bug (two schemes each registering
their own number concentration under slightly different names). The matched
indices even shifted because of it (NUMLIQ went from idx 6 to idx 7). Catch this
in the map and you've found a real problem in ten seconds, before reading a
single diff.

## TL;DR checklist

1. **Check the constituent map** — right species count, no unmatched
   duplicates?
2. **Is the first compared scheme already diffing?** Read the `alignment
   check` verdict near the bottom: it says whether the step pairing/snapshot
   is at fault (nearly every arg differs) or the differing inputs are
   field-specific upstream issues (most args bitwise match the paired step).
3. **Go to the FIRST DIVERGENCE box.** `[INPUTS DIFFER]` = upstream;
   `[OUTPUTS DIFFER]` = in-scheme.
4. **Tag every diff by field and trace it to its first appearance.** That first
   appearance is the fix; the rest is usually propagation.
5. **Filter the noise:** `[note]` lines, tiny floor values (`1e-37` vs `0`),
   utility-routine hit-count mismatches. Note them, don't chase them.
6. **`nan` or `sima=0`-everywhere = setup/wiring**, not a physics answer.

Everything is rerunnable and the raw dumps are self-describing — once you've
spotted the lead, the README's "Digging deeper after a run" section shows how to
go from a line in the report to the exact statement that wrote the bad value.
