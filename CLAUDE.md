# helios_postprocess -- Project Guide

## Repository
- **GitHub**: `tamehlhorn/helios_postprocess`
- **Package**: `helios_postprocess/`
- **Version**: 3.0.0 (March 2026)
- **Dev machine**: MacBook (`~/Codes/helios_postprocess/`) -- editing, pushing; use `python` not `python3`
- **Run machine**: Mac Studio (`tommehlhorn`, `~/helios_postprocess/`) -- use `python3`
- **NOTE (June 2026)**: Both the GitHub repo and both filesystem checkouts use
  the short name `helios_postprocess` (no `-or` suffix). The doc previously
  said `helios_postprocessor` everywhere — that was wrong on every machine.
  Only the legacy docx filenames `helios_postprocessor_guide.docx` and
  `helios_postprocessor_report.docx` still carry the old name (referenced
  below in "Diagnostic Scripts"); leave those as-is unless renamed on disk.
- **MacBook python**: Anaconda (`python` command); Mac Studio uses `python3`
- **Package install**: `pip install -e . --user` (MacBook); `pip install -e .` (Mac Studio)

## Architecture

### Class-Based Pipeline (primary workflow)

```
HeliosRun(exo_path)  ->  build_run_data(run)  ->  ICFRunData (dataclass, ~68 attributes)
                                                    |
                                             ICFAnalyzer(data)
                                               .analyze_drive_phase()
                                               .analyze_stagnation_phase()   <- also calls analyze_implosion_phase()
                                               .analyze_burn_phase()
                                               .compute_performance_metrics()
                                                    |
                                             ICFPlotter(data, config)
                                               .create_full_report(pdf_path)
                                                    |
                                             ICFOutputGenerator(data)
                                               .write_all(base_path)   -> _summary.txt + _history.csv
```

### Runner Script (recommended entry point)

```bash
# Mac Studio
python3 ~/helios_postprocess/examples/run_analysis.py <base_path>

# MacBook
python ~/Codes/helios_postprocess/examples/run_analysis.py <base_path>
```

Takes a path WITHOUT extension and derives all filenames:
- `<base>.exo` -- EXODUS input (required)
- `<base>.rhw` -- RHW input (optional, auto-detected)
- `<base>_published.json` -- published data for comparison (optional)
- `<base>_report.pdf` -- output: diagnostic plots
- `<base>_comparison.pdf` -- output: comparison table + burn history plots
- `<base>_summary.txt` -- output: text summary (includes burn-avg + comparison)
- `<base>_history.csv` -- output: time histories

### Burn-Averaged Metrics (published-data comparison)

```python
from helios_postprocess.burn_averaged_metrics import (
    extract_histories_from_run_data,
    calculate_burn_averaged_metrics,
    compare_with_published,
)

histories = extract_histories_from_run_data(data)
metrics = calculate_burn_averaged_metrics(histories)
print(compare_with_published(metrics, published_data, laser_energy_MJ=4.0))
```

### Additional Physics Modules (functional style)

| Module | Purpose |
|--------|---------|
| `energetics` | Kinetic energy, hydro efficiency, PdV work |
| `neutron_downscatter` | Down-scatter ratio (DSR) diagnostics |
| `pressure_gradients` | Shock identification, RT instability assessment |

These modules are standalone and not yet wired into ICFAnalyzer -- see Open Items.

## Active Source Files

| File | Lines | Role |
|------|-------|------|
| `core.py` | ~965 | `HeliosRun` -- EXODUS/netCDF4 reader |
| `data_builder.py` | ~555 | Bridge: `HeliosRun` -> `ICFRunData` dataclass |
| `icf_analysis.py` | ~1340 | `ICFAnalyzer` -- drive, stagnation, burn, implosion, IFAR, mass fractions, burn propagation, convergence ratios |
| `icf_plotting.py` | ~1580 | `ICFPlotter` -- full PDF report |
| `icf_output.py` | ~435 | `ICFOutputGenerator` -- summary text + CSV time histories |
| `burn_averaged_metrics.py` | ~560 | Temporal burn-averaging, implosion metrics, published-data comparison |
| `energetics.py` | -- | Kinetic energy, hydro efficiency, PdV work |
| `neutron_downscatter.py` | -- | Neutron down-scatter ratio diagnostics |
| `pressure_gradients.py` | -- | Pressure gradient analysis, shock ID, RT assessment |
| `rhw_parser.py` | ~240 | `RHWParser` -- reads `.rhw` input files |
| `__init__.py` | ~90 | Package exports |
| `examples/run_analysis.py` | ~350 | CLI runner with auto file derivation and comparison PDF |

## Published Data Comparison (JSON format)

Place a `<name>_published.json` file next to the `.exo` file. Format:

```json
{
    "laser_energy_MJ": 2.15,
    "T_hs":    [22.5, 2.0],
    "P_hs":    [193, 20],
    "rhoR_cf": [0.52, 0.05],
    "CR_max":  [29.0, 3.0],
    "yield":   [20.6, 1.0],
    "gain":    [9.6, 0.5],
    "peak_velocity_kms":    [410, 0.0],
    "adiabat":              [0.0, 0.0],
    "ifar":                 [0.0, 0.0],
    "hydro_efficiency_pct": [0.0, 0.0],
    "imploded_DT_mass_mg":  [0.0, 0.0],
    "inflight_KE_kJ":       [0.0, 0.0],
    "fraction_absorbed_pct":[97.0, 0.0],
    "P_hs_ignition_Gbar":   [75.0, 0.0],
    "hs_radius_ignition_um":[120.0, 0.0],
    "T_ion_onaxis_ignition_keV":[15.0, 3.0]
}
```

Each entry is `[value, uncertainty]`. Entries with `[0.0, 0.0]` are skipped.

## Unit Conventions

| Quantity | Internal Unit | Conversion from Helios |
|----------|--------------|----------------------|
| Temperature | eV (arrays), keV (scalars/reporting) | Helios stores eV; /1000 for keV |
| Pressure | J/cm3 (arrays), Gbar (reporting) | x1e-8 from J/cm3 |
| Velocity | cm/s (arrays), km/s (reporting) | x1e-5 from cm/s |
| Areal density | g/cm2 | rhoR = integral(rho dr) |
| Laser energy | MJ | x1e-6 from J |
| Drive temperature | eV | Stays in eV (radiation drive) |
| Time | ns | Helios may store in seconds -> x1e9 |

## Physics Conventions

1. **Pressure aliasing**: `rad_pressure = elec_pressure` (or `elec + actual_rad` when 3-component).
   Total pressure = `ion_pressure + rad_pressure` everywhere. Documented in data_builder.py.

2. **Stagnation** = minimum hot-spot outer radius (NOT peak density).
   For igniting capsules, peak rho occurs after alpha heating.

3. **Ablation front** = steepest negative drho/dr outside the hot spot boundary.
   Smoothed with iterative algorithm anchored at stagnation.
   Stored as zone indices (`ablation_front_indices`) -- use these, NOT the smoothed radius,
   for mass fraction and imploded mass calculations.

4. **Region boundaries**: `region_interfaces_indices` are NODE indices (not zone indices).
   - Hot-spot boundary: `ri[:, 0]`
   - Cold fuel / ablated fuel boundary: `ri[:, 1]`
   - Fuel/ablator boundary: `ri[:, -2]`

5. **Fusion yield**: Preferred method uses `TimeIntFusionProd_n_1406 x 17.6 MeV x 1.602e-19 MJ/MeV`.
   Fallback: time-integrate `FusionRate_DT_nHe4` (less precise).

5b. **FusionRate_DT_nHe4 units (CORRECTED May 27 2026)**: Helios emits this
   variable as **reactions/s per gram of plasma** in each zone (a
   mass-specific intensive rate), NOT reactions/s/zone (volume-integrated
   extensive rate). Earlier comments in this file and in
   `icf_analysis.py` describing it as "reactions/s per zone,
   volume-integrated" are wrong and need updating.

   Confirmed by zero-D static verification (May 27 2026,
   `examples/verify_zero_d.py`): uniform sphere at n_DT=1e22 cm⁻³,
   T=2.982 keV, 100 zones, ρ=0.0418 g/cc. Each zone reports
   FusionRate_DT_nHe4 = 1.090e+27, which equals the Bosch-Hale per-gram
   analytic to 0.04%. Summing across zones and dividing by total mass
   gave 4.5×10⁶× the correct value — a definitive disproof of the
   per-zone-per-second interpretation.

   **Correct usage patterns:**
   - **Per-gram bulk rate (uniform conditions)**: any zone value gives it
     directly; for non-uniform: `np.average(fusion_power[t], weights=zone_mass[t])`
   - **Total reactions/s in the run**: `np.sum(fusion_power[t] * zone_mass[t])`
   - **Cumulative reactions to time t**: time-integrate the total
     reactions/s above (use `TimeIntFusionProd_n_1406` preferentially)
   - **Bang time** (`np.argmax(np.sum(fusion_power, axis=1))`): works
     either way since only the time of peak matters; the absolute value
     is summed-mass-specific and not physically meaningful by itself

   `icf_analysis.py` line ~2145 uses `np.sum(fusion_rate[t]) * dt` for
   the neutron-yield fallback path — this is incorrect under the
   corrected interpretation and should be `np.sum(fusion_rate[t] * zone_mass[t]) * dt`.
   The preferred yield path (TimeIntFusionProd_n_1406) is unaffected.
   Production runs validated against Olson_PDD_9 (yield 20.6 MJ) use
   the preferred path, so this fallback-path bug is academic until a
   run lacks the cumulative-product field.

6. **Mass fractions**: Use zone-index method exclusively (no temperature mask).
   - `unablated_fuel_mass` = fuel zones (0..fuel_bnd) inside ablation front / initial_fuel_mass
   - `unablated_ablatar_mass` = ablator zones (fuel_bnd..) inside ablation front / initial_ablator_mass
   - `stagnated_fuel_mass` = fuel zones between hs_bnd and ablation front / initial_fuel_mass
   - Typo `unablated_ablatar_mass` (ablatar not ablator) is intentional -- do not fix without
     coordinated find-replace across all files.

7. **Burn propagation** (Olson et al. convention): Tracks hot-spot rhoR vs total rhoR over time.
   Ignition identified when hot-spot rhoR (absolute) >= 0.3 g/cm2. The T_ion > 4.5 keV mask
   is used only for the burn propagation plot. Scalars ignition_time, ignition_hs_pressure,
   ignition_hs_radius all use the 0.3 g/cm2 threshold.

8. **EXODUS sampling principle**: EXODUS files contain only a fraction of the actual simulation
   timesteps. Helios's own time-integrated quantities are authoritative. The burn rate from
   EXODUS is used only as a weighting function for burn-averaging, never for absolute yield.

9. **Adiabat**: alpha = P / P_Fermi where P_Fermi = 2.17 (rho/rho_0)^(5/3) Mbar,
   rho_0 = 0.205 g/cc (equimolar DT ice, Lindl convention). Evaluated at peak
   implosion velocity (pre-stagnation) in the cold unablated fuel only
   (ri[t, 0] to ri[t, 1]), excluding the ablated fuel region.

10. **Hydrodynamic efficiency**: eta_hydro = max(KE_inward) / E_absorbed.
    KE_inward = sum of 0.5 m v^2 for zones with v < 0 (imploding), maximized
    over all timesteps. E_absorbed = max(laser_energy_deposited) in Joules.

11. **IFAR (In-Flight Aspect Ratio)**: IFAR = R_shell / Delta_R at peak
    implosion velocity. Shell boundaries determined from density profile using
    rho > rho_peak / e threshold (NOT Lagrangian region interfaces).

12. **Convergence ratios** -- two distinct quantities, both computed and reported:
    - `comp_ratio` = CR_stag = R0 / R_hs_at_stagnation
      R0 = initial inner shell radius = zbnd[0, ri[0,0]]
      R_hs = hot-spot boundary radius at stagnation = zbnd[stag_idx, ri[stag_idx,0]]
      Validated: Olson_PDD_9 gives 29.6 vs published 29.0 (2.1%)
    - `cr_inflight` = R0 / R_ablfront_at_peak_velocity. Uses ablation_front_radius[pv_idx].
      Validated ~4.5-4.6 for TM/26-series (physically correct). FIXED April 2026.
    - DO NOT use density ratio for CR -- that was the old incorrect implementation,
      now fixed as of April 2026.

13. **Hot-spot radius**: Use region interface `ri[stag_idx, 0]` (node index) to get
    hot-spot boundary radius from `zone_boundaries[stag_idx, hs_node]`.
    DO NOT use temperature mask -- alpha heating warms dense shell above 1 keV.
    NOTE: `_compute_hot_spot_properties()` currently returns inconsistent values
    (0.1786 cm vs correct 0.0068 cm for PDD_9) -- known bug, see Open Items.

14. **Imploded DT mass**: Use zone-index method in burn_averaged_metrics.py --
    sum `zone_mass[stag_idx, :min(abl_idx+1, fuel_bnd)]`.

15. **Peak velocity index**: Stored as `data.peak_velocity_index` (integer timestep index).
    Search restricted to shell zones only (outside ri[:, 0]) and strictly pre-bang time
    (t < bang_time, strict less-than). This prevents the stagnation/ignition velocity
    spike from being incorrectly identified as the peak implosion velocity.
    Fixed April 2026.

16. **Flux-limiter convention + regime (May 2026, updated May 29 2026)**:

    **CORRECTED GLOBAL VALUE (May 29 2026): FL_prism = 0.012** (f_standard
    = 0.048) is the appropriate global flux-limiter value for the present
    version of Helios, replacing the per-region 0.06 DT / 0.007 foam-CH
    split that was used in the fab007 production calibration. This
    correction came from a separate Prism-side calibration update and is
    treated here as the standing default for new design-study work.

    Note: applying FL_prism=0.012 to all 4 regions changes the production
    fab007 picture meaningfully — DT regions go 0.06 → 0.012 (5× less
    conduction allowed), foam/CH regions go 0.007 → 0.012 (~1.7× more).
    See `notes/foam_vs_ice_investigation.md` §3 dated 2026-05-29 for the
    full design study (FL=0.012 + geometric defocus reduction recovers
    72% of the Helios-vs-LILAC yield gap with 31.5% foam yield share at
    `wf_fl012_c20_burn`).

    **KNOWN ISSUE (May 29 2026): FL_prism = 0.012 + burn-on configurations
    hit recurring SIGSEGV/malloc-error crashes during Helios shutdown.**
    Two of three burn-on FL=0.012 runs crashed (baseline_burn on first
    attempt, c20_burn on re-extraction attempt). The .exo files are
    substantially complete at the time of crash and remain recoverable for
    postprocessing — the metrics extract correctly. Worth flagging for
    larger scans at this FL: budget time for at least one retry per run.

    a) **Per-region read**: Helios reads the thermal-conductivity flux
       limiter `f` from the .rhw file **per region** (not a single global
       default). The values are parsed into
       `data.flux_limiter_per_region` and reported uniformly or per-
       region in the summary `Flux limiter (f)` block.

    b) **Convention factor**: Prism's Helios reports `f` as a value that
       is **4× smaller than the ICF-community standard convention**
       (Spitzer–Harm cap). To convert:

            f_standard  =  4 × f_prism

       So `f_prism = 0.015` is equivalent to `f_standard = 0.060` (the
       simulation default for ICF). A run reported as `fab02` /
       `f = 0.020` is actually at `f_standard = 0.080`. **This is the
       Prism team's working hypothesis as of May 2026 — pending formal
       confirmation. NOT directly verified from EXODUS diagnostics
       because the FL doesn't engage in the calibration-range runs
       (see (c) below).**

    c) **Regime sensitivity (institutional knowledge, May 2026)**: the
       flux-limiter knob has **two distinct regimes** in Helios behavior:

       - **Small-f regime (FL-dominated)**: heat conduction strongly
         capped, T_corona spikes, ablation pressure and implosion
         velocity drop sharply with `f`. The simulation is sensitive
         to `f` here.
       - **Large-f regime (geometry-dominated)**: the cap is loose
         enough that q_actual rarely hits f·q_FS in the conduction
         zone; the simulation behavior is set by beam geometry
         (cone, spot, focus), absorption physics, and adiabat — not
         by `f`. Changing `f` between, say, 0.015 and 0.06 (Prism)
         produces minimal change in V, T_crit, or HS ρR.

       The current Olson_PDD_20 calibration runs span `f_prism ∈
       [0.015, 0.06]` which is **all in the large-f regime** — the FL
       knob does not engage there. `examples/check_flux_limiter.py`
       confirmed this empirically: T_crit ≈ 2.4–2.8 keV regardless of
       `f_file`, with `q_SH/q_FS` saturation-indicator implying the
       limiter *should* engage but the observed temperature ratio
       doesn't show the expected `f^(-2/3)` scaling.

       To verify the convention factor and locate the FL-dominated
       regime, run a slab/sphere problem at `f_prism = 0.003`
       (`f_standard = 0.012` under ×4) — if FL is real and the ×4
       hypothesis holds, this should produce a sharp T_corona rise
       (~25 keV) and a substantially slower implosion. If it doesn't,
       either (i) the FL value isn't being plumbed end-to-end through
       Helios's electron conductivity, or (ii) the convention factor
       is even larger than ×4. Either result is calibration-relevant.

    `examples/scan_summary.py` exposes both columns (`flux_limiter_prism`,
    `flux_limiter_standard`) so cross-run comparison against published
    `f` values is unambiguous; `check_flux_limiter.py` plots
    `q_SH/q_FS` saturation indicator and prints `T_crit` for inferred
    regime classification.

17. **Alpha-deposition mode (May 2026)**: `data.alpha_deposition_local` and
    `data.alpha_deposition_nonlocal` flags from the RHW parser determine the
    alpha-transport model. **Local deposition** (alpha energy deposited in
    the zone of birth instantaneously) over-estimates burn / hot-spot
    pressure / yield by typical factors of 2-5× relative to non-local
    transport (which respects mean free path).

    For LILAC / xRAGE / HYDRA comparisons, use **non-local transport only**
    (`alpha_deposition_nonlocal = True`, `alpha_deposition_local = False`).
    Local-deposition runs are diagnostic curiosities and should be flagged
    by their `alpha_deposition_mode = 'local'` column in `scan_summary.csv`
    before being included in any calibration comparison.

    Example diagnosis: `PDD_TM_4nb` (local α) gives gain = 92, yield = 131 MJ
    even though its cold-implosion diagnostics (coupling 68.5%, HS ρR 1.01,
    adiabat 1.98) match LILAC well — the burn numbers are inflated by the
    local-α model, not by physics.

## Laser Deposition Model (Helios vs reference codes)

Helios uses a 1D spherical ray-trace with refraction (Snell's law, geometrical optics).
Parameters: focus position d (cm from plasma origin), half-cone angle, spot size s,
number of rays = number of zones (default). Rays originate from the focal plane,
treated as parallel at that plane.

**Key calibration issue (April 2026):** Helios 1D model cannot capture the 3D geometric
effect of NIF PDD beams increasingly missing the imploding capsule as it shrinks. In 3D
codes (HYDRA/LILAC/Xrage), fixed-pointing NIF beams illuminate a smaller fraction of
the shrinking capsule at late times -- this is captured automatically by their 3D ray trace.
Helios's 1D model maintains spherical symmetry throughout, so all cone rays continue to
interact with the target regardless of radius.

**Practical consequence:** Helios over-drives the capsule at late times relative to the
reference codes, producing peak implosion velocities ~40% above the LILAC reference (587
vs 410 km/s for PDD_22). The hydro efficiency is locked at ~10.5% regardless of laser
geometry parameters, indicating the excess drive is in the ablation physics rather than
purely geometric coupling.

**Calibration approach:** The correct approach requires empirical comparison of cumulative
absorbed energy histories between Helios and the reference codes (HYDRA/LILAC/Xrage).
A time-dependent power multiplier applied to the RHW power table can compensate, but must
be anchored to reference absorbed-energy data rather than a geometric model.
The focus position parameter d offers natural geometric defocusing: when d ≈ R_initial,
coupling is maximised early and decreases as R(t) < d. This is the most physically
motivated knob available in Helios without modifying the power table.

**CBET:** Cross-Beam Energy Transfer is absent from Helios. In NIF PDD geometry, CBET
reduces drive primarily at the equatorial region as rays cross during implosion. The
Olson 2021 HYDRA simulations excluded CBET explicitly because the quarter-critical radius
does not shrink below the original capsule radius during their pulse. So CBET is not the
primary calibration target for this design.

## Helios Laser Source Parameters (RHW format)

Key fields in `[Laser Source Data]` block (beam 1):
- `Focus position` (cm) -- distance from plasma origin to focal plane; set to ~R_initial
  (~0.23 cm for Olson PDD target) to get natural geometric defocusing
- `Half cone angle` (deg) -- NIF PDD rings at 23.5°, 30°, 44.5°, 50°; use 35° for average
- `Spot size` (cm) -- 1/e radius for Gaussian profile; keep ≥ 0.02 to avoid singularity
- `Laser spatial profile model` -- 0=Uniform, 1=Gaussian
- `Number of points at focus` -- must be > 1 for non-zero spot size
- `Beam energy normalization` -- MUST be 1.0 for single active beam; 0.5 caused 2× energy
  deficit in PDD_27 (critical bug found April 2026)
- `Number of laser beams` -- if set to 3 but only beam 1 is active, beams 2/3 must have
  `Laser power model is on = 0`

## PDD Design-Study Anchors: c33_burn + c25_burn (May 30 2026)

The full FL=0.012 burn-on cone-angle sweep revealed a **two-plateau staircase**, not a
smooth ramp. Two recommended anchors for two different objectives:

**`wf_fl012_c33_burn`** — closest-to-LILAC anchor (minimum geometric intervention):
- cone 33°, spot 0.16 cm, focus 0.20 cm, FL_prism = 0.012 (global)
- Yield: **55 MJ** with **25.3% foam yield share** (matches fab02's 26.5% baseline)
- V_peak 447 km/s (+9% vs LILAC's 410), coupling 79.6%, adiabat 1.54, HS ρR 0.53

**`wf_fl012_c25_burn`** — max-yield anchor (engineering ceiling):
- cone 25°, spot 0.14 cm, focus 0.16 cm, FL_prism = 0.012 (global)
- Yield: **69 MJ** with **31.1% foam yield share** (exceeds fab02's 26.5%)
- V_peak 470 km/s (+15% vs LILAC), coupling 85.3%, adiabat 1.01, HS ρR 0.64

**For HDD calibration transfer: use c33-class as the working geometric anchor.** Same
strong foam burn as c25/c20, less drive overshoot, lower kinematic deviation from
LILAC. Use c25-class only when maximum yield is the objective.

**Staircase plateaus (all at FL=0.012, burn on):**

| Cone | Yield (MJ) | Foam share | Note |
|---:|---:|---:|---|
| 37° baseline | 33 | 14% | marginal (FL fix alone) |
| 33° | 55 | 25% | **Plateau 1** — bootstrap baseline |
| 28° | 55 | 26% | Plateau 1 (foam over-compresses to 62 g/cc but burns no more) |
| 25° | 69 | 31% | **Plateau 2** — unlocks at (spot 0.16→0.14, focus 0.18→0.16) |
| 20° | 69 | 32% | Plateau 2 — saturated |

**c20 is superseded by c25**: same yield, same foam share, marginally less aggressive
geometry — use c25 as the cleaner choice when max yield is needed.

**Sub-finding:** c28's foam over-compression (61.6 g/cc, 2× ice's 30 g/cc) at unchanged
foam burn productivity confirms **foam burn is bootstrap-temperature-limited at this
drive, not compression-limited.** Consistent with May 24 fab007 diagnostic.

**KNOWN ISSUE:** FL=0.012 + burn-on hits recurring SIGSEGV/malloc-error crashes during
Helios shutdown (~4 of 6 runs in this scan). The .exo files are substantially complete
and metrics extract correctly; budget one retry per run for larger scans.

See `notes/foam_vs_ice_investigation.md` §3 dated 2026-05-29-30 for the design study
narrative; `docs/Xcimer_foam_burn_deficit_report.md` §5 for the external write-up;
`comparisons/pdd_design_comparison.png` for the 10-row decision-matrix figure.

## PDD Calibration Scan (April 2026)

Goal: reproduce LILAC peak implosion velocity ~410 km/s for Olson PDD target at 1.4×
drive multiplier. Reference pulse: foot ~23-25 TW (0-5 ns), ramp 5-9 ns, peak ~329 TW
(9-12.7 ns). LILAC reference energy: 2.150 MJ.

| Run | Spot (cm) | Cone (°) | Profile | d (cm) | Norm | v (km/s) | Abs (%) | Yield (MJ) | Adiabat | Notes |
|-----|-----------|----------|---------|--------|------|----------|---------|------------|---------|-------|
| PDD_9  | 0.00 | 1  | — | 0.00 | — | 504 | 100 | 20.6 | 3.03 | Reference (wrong pulse shape) |
| PDD_22 | 0.12 | 35 | Gaussian | 0.20 | — | 587 | 87.5 | 75.6 | 1.22 | Ignites |
| PDD_24 | 0.12 | 35 | Gaussian | 0.20 | — | 587 | 87.5 | 75.7 | 1.22 | Same as PDD_22 |
| PDD_25 | 0.16 | 35 | Gaussian | 0.20 | — | 563 | 82.4 | 58.4 | 1.05 | Larger spot |
| PDD_26 | 0.16 | 15 | Gaussian | 0.20 | — | 556 | 87.2 | 59.5 | 1.09 | Narrow cone |
| PDD_26a| 0.16 | 15 | Gaussian | 0.20 | — | 1226 | 87.2 | 130 | 102 | Local burn -- discard |
| PDD_27 | 0.16 | 11 | Uniform | 5.00 | 0.5 | 73 | 10 | 0 | — | Focus too far + norm=0.5 bug |
| PDD_28 | 0.02 | 35 | Uniform | 0.23 | 1.0 | ~1226 | 87.2 | 127 | 111 | Over-driven (peak v bug was active) |
| DM_01b | 0.00 | 3.5 | Uniform | -2.50 | 1.0 | 405 | 90.9 | 1.58 | 1.56 | Inverted focus + narrow cone; matches LILAC v within 10%; 1.365 MJ drive |
| fab015_foot25_s018_c37 | 0.18 | 37 | Uniform | 0.22 | 1.0 | 429 | 74.4 | 29.3 | 1.75 | First geometry+FL closure attempt (May 2026); superseded by fab007 |
| **fab007_foot25_s018_c37** | **0.18** | **37** | **Uniform** | **0.22** | **1.0** | **421** | **73.1** | **26.0** | **1.95** | **★ PDD production calibration (May 2026 v2), see closure note below** |

**Key findings:**
- Hydro efficiency locked at ~10.5% regardless of geometry -- excess drive is in ablation physics
- Absorption insensitive to cone angle/spot (87±3%) -- refraction controls deposition, not geometry
- Narrower cone → lower adiabat (deeper deposition) but minimal velocity effect
- PDD_27 failure: d=5.0 cm (rays miss entirely) + norm=0.5 (half power) -- two compounding bugs
- PDD_28 peak velocity was artifactually high due to velocity bug (now fixed); rerun needed

### PDD calibration closure (May 2026 v2)

Production calibration point: **`Olson_PDD_20_fab007_foot25_s018_c37_burn`**

Supersedes the earlier `fab015_foot25_s018_c37` closure (May 21). The
FL Prism convention test (fab003) demonstrated the FL-dominated
regime exists and that f_prism between 0.003 and 0.015 brackets the
ignition cliff. fab007 sits in the saturation knee, giving the best
LILAC thermodynamic-state match while retaining ignition.

| Parameter | Value | Note |
|-----------|-------|------|
| cone half-angle | 37° | wider than baseline 35° for geometric defocus |
| spot radius | 0.18 cm | larger than baseline 0.16 cm for late-time defocus |
| focus position d | 0.22 cm | ≈ R_initial |
| FL (Prism) | **0.007** | f_standard = 0.028, at the FL saturation knee |
| foot power | 25 TW | unchanged from baseline pulse table |
| α deposition | non-local | clean (no local-α yield inflation) |
| burn | ON | |
| **eff_avg_coupling** | **73.1%** | from `data.eff_avg_coupling_pct` |

**LILAC kinematics + thermodynamic state matched** (within 7% on every non-mass-deficit metric):

| Metric | This run (fab007) | LILAC | Δ |
|--------|------------------:|------:|---|
| V_peak (km/s) | 421 | 410 | **+2.7%** ✓ |
| CR_max | 30.7 | 29.0 | +5.9% |
| Peak total ρR (g/cm²) | 1.01 | 1.05 | −3.6% ✓ |
| **⟨T_hs⟩ (keV)** | **23.1** | **22.5** | **+2.7%** ✓ |
| **⟨P_hs⟩ (Gbar)** | **206** | **193** | **+6.6%** ✓ |
| Adiabat | 1.95 | 3.0 | −35% |
| Foot shock (ns) | 6.05 | 7.5 | −19% |
| Ramp shock (ns) | 8.15 | 10 | −19% |
| Peak shock (ns) | (lost in consolidator) | 13 | — |
| Ignition | YES | YES | ✓ |
| HS pressure at ignition (Gbar) | 237 | 90 | +163% (over) |
| HS ρR T>4.5 (g/cm²) | 0.35 | 0.85 | −58% (mass deficit) |
| Yield (MJ) | 26 | 87 | −70% (mass deficit) |

**Corrected diagnosis (May 24 2026 final, supersedes May 23 burn-front-arrest framing): the remaining HS ρR / yield gap is anemic alpha bootstrap at the alpha-ignition cliff, limiting bootstrap penetration into the foam VOLUME. Not interface physics, not EOS, not transport at the Z̄ jump.**

Four earlier framings of the residual were wrong, each retired in turn:

1. **Target-mass deficit (April 2026, retired May 22):** the 0.60 mg "imploded DT" tag was misread. The radial-at-ignition CSV showed Helios's hot spot (r_hs = 148 µm) extends to within 6 µm of the ice/foam boundary (r₂ = 153.8 µm) — nearly the entire DT ice layer is *inside* the hot region. The "cold fuel" tag was capturing DT-CH foam, not missing DT mass. The DT ice is conserved.

2. **Early-shock heating (May 22, retired May 23 first revision):** the ~29 keV on-axis T_ion at ignition was *not* coming from shock convergence heating of the central gas. It's alpha self-heating. Confirmed by running fab007 with burn flag OFF (same geometry, same pulse, no alpha deposition) — the hydrodynamic state is essentially identical to fab007 burn-on (V_peak 421, peak ρ 165 g/cc, peak total ρR 1.25, foot shock 6.05 ns, ramp shock 8.15 ns) but ⟨T_hs⟩ falls to 5 keV. Energy-ledger comparison (`energy_balance_diagnostic.py`) confirms the pre-alpha hydrodynamic state is identical between burn-ON and burn-OFF runs to <0.1% on every tracked channel at peak velocity.

3. **Generic "alpha-amplified hot-spot expansion" (May 23 first revision, retired May 23 final):** too vague.

4. **Burn-front arrest at the DT-ice / DT-CH-foam interface (May 23 final, retired May 24):** the interface ITSELF is benign in fab007 — 3T quantities and α-deposition rate are essentially continuous across zone 190 (T_e ratio 1.017, T_rad 1.000, α-dep ratio 1.06). The Z̄ jump is only 8% (0.999 → 1.084), because the "DT-CH foam" is in fact 0.222 g/cc DT in a 0.020 g/cc CH structural binder — ~92% DT by mass, not a high-Z mixed material. There's no significant Z̄ discontinuity to drive a brems sink or alpha-MFP cliff. The phenomenon is real (foam contributes only 10% of fab007 yield, vs 26.5% in fab02) but is BODY-level (anemic bootstrap penetration through foam VOLUME), not interface-level. The ice→foam swap experiment that motivated this framing still stands but is consistent with the body-level picture: replacing the ice with foam removes the burn SUBSTRATE that sustains a hot enough alpha source to fire any foam mass at all — not the same as "burn arrests at the interface".

**The correct picture: alpha-bootstrap strength sets foam-burn penetration depth.**

Per-zone diagnostics (`examples/dump_per_zone_burn_share.py`, `examples/dump_burn_rate_timing.py`) on fab007 (yield-limited) vs fab02 (high-yield over-drive) give the conclusive side-by-side:

|  | fab02 (59 MJ) | fab007 (26 MJ) |
|---|---|---|
| Ice CR_max | 1390 | 1080 |
| **Foam CR_max** | **1458** | **979** |
| Foam CR / Ice CR | **1.05** | **0.91** |
| Ice ⟨T_ion⟩ at bang | 34 keV | 17 keV |
| Ice P_avg at bang | 408 Gbar | 116 Gbar |
| **Foam ⟨T_ion⟩ at bang** | **4.4 keV** | **1.85 keV** |
| Foam zones crossing 4.5 keV | 68% | 21% |
| % yield from foam | **26.5%** | **10.1%** |
| Foam yield / DT-mass (kJ/mg) | 3,811 | 642 |
| Ice yield / DT-mass (kJ/mg) | 70,064 | 36,786 |

What the data refutes:
- **PROPACEOS foam EOS is NOT broken.** Foam reaches peak CR equal to or higher than ice in both runs at compression times within 30-40 ps of the ice peak. Foam compresses fine.
- **The interface is NOT a brems/transport barrier.** Continuity across zone 190 in fab007 at peak burn (`dump_burn_propagation_profile.py`).

What the data supports:
- The foam-burn efficiency is set by the **alpha-bootstrap STRENGTH at bang time** — how hot/dense the ice hot spot becomes, which sets how much alpha-deposited power reaches the foam VOLUME during the ~0.3 ns burn FWHM. fab02's ice hot spot is 3.5× hotter and 3.5× higher-pressure than fab007's; fab02's bootstrap fires 68% of foam zones, fab007's fires 21%.
- fab007's foam is actually DENSER at bang time (12.9 g/cc avg) than fab02's (7.4 g/cc) — fab02's stronger bootstrap has already started expanding the foam by bang. The foam compresses fine in both cases; only its heating rate differs.
- The fab007 residual is **bootstrap reach** (compatible with what was previously labeled "mechanism 2 — burn duration too short"). The foam is well-compressed but cold (1.85 keV avg) because the ice bootstrap is anemic.
- Even fab02's foam burns ~18× less per unit DT mass than its ice (3811 vs 70064 kJ/mg). LILAC's 87 MJ implies near-parity. Residual likely (i) Helios 1D non-local α-transport details at Z̄≈3, (ii) ~0.3 ns burn FWHM (would need longer compression hold), (iii) fab02's over-drive produces a hot ice hot spot but not LILAC-like compression history.

**Restated calibration residual:** Helios and LILAC compress the target similarly (V_peak, CR within ~5%, foam total compression matches ice). Both can ignite the DT ice layer. The yield gap (26 vs 87 MJ for fab007) is set by **how strongly the ice hot spot ignites** (sub-cliff in fab007, over-driven in fab02, calibrated drive in LILAC) and **how much alpha power that bootstrap delivers to the foam volume** during burn FWHM.

**Implications:**
- fab007 is the best Helios produces *at LILAC's calibrated thermodynamic state*. It sits on the alpha-ignition cliff and that limits foam-burn penetration. Its drive-phase calibration remains correct; the foam-burn deficit is downstream.
- `Olson_PDD_20_fab02_foot25_s016_burn` (coupling 84%, V=463, α=1.05) is the bootstrap-strength reference. It over-drives but produces robust foam-penetrating burn (26.5% of yield from foam).
- **HDD calibration risk: PARTIALLY RELAXED.** Helios CAN fire foam DT mass given a strong alpha bootstrap (proven by fab02). The risk is conditional on landing HDD in a robust-bootstrap regime (ice-equivalent hot spot ≥30 keV / ≥400 Gbar at bang) rather than a marginal-ignition regime like fab007. See HDD transfer plan below.

**Carry both runs as HDD reference baselines:**
- `Olson_PDD_20_fab02_foot25_s016_burn` — bootstrap-strength baseline. cone=20°, spot=0.16 cm, FL_DT=0.06, FL_foam+CH=0.02. Comparison anchor for what Helios CAN achieve when the alpha bootstrap fires cleanly.
- `Olson_PDD_20_fab007_foot25_s018_c37_burn` — thermodynamic-state baseline. Sensitivity-bracketing point at the alpha-ignition cliff with calibrated thermo state. Not the primary HDD transfer baseline.

**Diagnostic scripts (added May 24 2026, examples/):**
- `dump_burn_propagation_profile.py` — radial T_e / T_rad / T_ion / ρ / α-dep across the ice/foam interface at peak burn. Refuted the Z̄-jump brems-sink hypothesis.
- `dump_per_zone_burn_share.py` — cumulative DT neutron count per zone at end-of-run, region-aggregated. Showed fab007 fires 10% of yield from foam vs fab02's 26.5%.
- `dump_burn_rate_timing.py` — per-zone burn-rate timing + compression-state by region. Refuted the PROPACEOS foam EOS hypothesis (foam CR ≥ ice CR in both runs) and isolated the discriminator (foam ⟨T_ion⟩ at bang).

### May 22-23 2026 perturbation study (retired)

A series of perturbations to fab007 was run to test the (then-current) early-shock-heating hypothesis: foot-power reduction (foot=15/18/20 TW, single-beam, burn ON) and foot-geometry defocus (cone 39/42/45° + spot 0.20/0.22/0.25 cm via `apply_foot_beam_split.py`, 2-beam recycling a placeholder beam, burn OFF). Results in three lines:

- **Foot-power scan (burn ON):** foot shock arrival slipped monotonically toward LILAC's 7.5 ns target (6.05 → 6.75 → 7.30 → 7.35 ns at foot=25 → 20 → 18 → 15). foot=18 nailed the timing target. **But all three perturbed runs failed to ignite** — yields collapsed to 0.25–0.30 MJ. Below the alpha-ignition cliff.

- **Foot-geometry defocus (burn OFF):** foot shock slipped 6.05 → 6.30 → 6.50 → 6.70 ns at cone=37 → 39 → 42 → 45°. Ramp shock moved the *wrong* way (8.15 → 7.95 → 7.80 → 7.60 ns) — weaker foot → less ice pre-compression → faster ramp shock. Foot-ramp gap shrank from 2.10 to 0.90 ns vs LILAC's 2.50 ns.

- **fab007 burn-OFF reference:** identical hydro state to burn-ON fab007 (V_peak 421, peak ρ 166, peak ρR 1.25, foot 6.05 ns, ramp 8.15 ns) but ⟨T_hs⟩ falls from 23 to 5 keV. **This is what closed the early-shock hypothesis** — the 23 keV in burn-on fab007 is alpha amplification of a 5 keV hydrodynamic floor, not shock convergence heating.

**Lessons:**
1. Foot intensity (geometric or power) IS a real lever for foot-shock arrival, ~+0.7 ns slip per 25-percent foot-intensity reduction. Useful telemetry, not a calibration knob.
2. Ramp arrival is anti-correlated with foot intensity (weaker foot → less ice pre-compression → faster ramp). Foot defocus alone *can't* land both shocks at LILAC values; would need an independent ramp-defocus knob, which we didn't pursue.
3. fab007 sits exactly on the alpha-ignition cliff. **Every** perturbation we tried failed to ignite. Don't perturb the production calibration speculatively — pre-test with a no-burn equivalent first.

These perturbation runs (`foot25_s018_c37_nb`, `foot{20,18,15}_s018_c37_burn`, `footgeo_c{39s20,42s22,45s25}_nb`) are retained for diagnostic value but are not candidate calibration points. **fab007 (`Olson_PDD_20_fab007_foot25_s018_c37_burn`) is the production PDD calibration.**

**Tuning sensitivities observed during closure:**
- `cone+spot` is the dominant lever — coupling 84% → 74% via cone 35°→37° + spot 0.16→0.18
- `foot22` (vs `foot25`) loses ignition at this coupling. Foot drive sets a narrow ignition band.
- `s020/c38` over-defocuses → coupling 71.6%, ignition fails. `s018/c37` is at the ignition cliff.
- FL changes in `f_prism ∈ [0.015, 0.060]` produce minimal change — geometry-dominated regime.
- **FL saturation knee at `f_prism ≈ 0.005–0.010` (f_standard ≈ 0.02–0.04)** — established by
  bracketing tests (fab003 fails ignition, fab015 over-compresses, fab007 sits at the knee
  with best ⟨T_hs⟩/⟨P_hs⟩ match). The full FL saturation curve maps as:

  | f_prism | f_standard | coupling | V_peak | ⟨T_hs⟩ | Ignition |
  |--------:|-----------:|---------:|-------:|-------:|---------:|
  | 0.060   | 0.240      | ~85%     | 463    | 38     | yes (over-driven) |
  | 0.015   | 0.060      | 74.4%    | 429    | 25.0   | yes |
  | **0.007** | **0.028**  | **73.1%** | **421** | **23.1** | **yes (production)** |
  | 0.003   | 0.012      | 70.0%    | 412    | 16.3   | NO (under-driven) |

**Next runs in the cycle (May 22):**
- `Olson_PDD_20_fab007_foot25_s018_c37_burn` IS the calibration sweet spot
  (this is the production point above).
- Open question: the central T_ion plateau extending to ~25 µm at ignition
  (vs LILAC's smooth profile) -- could be (a) Lagrangian first-zone
  artifact, (b) AV multiplier effect, (c) actual physics. Test via:
  - Pull latest `examples/run_analysis.py` (commit 3429903+) to generate
    `<base>_radial_at_ignition.csv`. Zone-count of T_ion > 25 keV in
    that file resolves (a) directly.
  - AV scan at fab007 if .rhw exposes the coefficient.
- HDD calibration transfer (see "HDD calibration transfer plan" section below).

### HDD calibration transfer plan (May 2026, revised May 23)

The real production target is HDD (Vulcan / VI_6-class), not PDD. The PDD calibration above
establishes the *method* -- now port it to HDD.

**⚠️ CALIBRATION FLAG (updated May 24 2026 — supersedes May 23 "fundamentally limited"
framing): the PDD closeout established that fab007's foam-yield deficit is from anemic
alpha bootstrap at the ignition cliff, NOT a Helios-fundamental limit on burning foam
mass. fab02 (over-driven) achieves 26.5% of yield from foam zones, proving Helios CAN
propagate burn into wetted foam given a strong-enough alpha source. The PROPACEOS foam
EOS, the Z̄-jump physics, and the ice/foam interface are not the binding limits. See
the "Corrected diagnosis (May 24 2026 final)" PDD section above for the fab02-vs-fab007
table and refutation evidence.

For HDD calibration: aim for a robust-bootstrap regime (ice-equivalent hot spot
≥30 keV / ≥400 Gbar at bang time) rather than the marginal-ignition regime fab007 sits
in. HDD's wetted-foam cold fuel will burn proportionally to the alpha-source strength
delivered to it. If HDD ends up at fab007-like marginal ignition the foam yield will
under-perform; if it lands at fab02-like over-drive the foam will fire. The geometric
defocus knob from PDD applies, but the target coupling % may need to be HIGHER than
the PDD calibration value (74%) to land in the bootstrap-strong regime. Use fab02
settings (cone=20°, spot=0.16 cm, FL_DT=0.06, FL_foam+CH=0.02) as the primary
HDD-transfer geometric baseline, fab007 as a sensitivity-bracketing reference only.**

**What carries over directly:**

1. **Geometric defocus is the primary calibration knob.** PDD closed at cone 35°→37°,
   spot 0.16→0.18 cm with coupling dropping 84% → 74%. Expect HDD to need a similar
   relative defocus, anchored on its own LILAC-equivalent target velocity.
2. **`eff_avg_coupling_pct` is the headline diagnostic.** Pipeline reports it directly;
   `compare_absorbed_energy.py` and `scan_summary.csv` make cross-run comparison trivial.
3. **3-shock-train detection** confirms structural similarity to the reference code.
4. **Foot drive is a secondary knob with a narrow ignition band.** Don't change it before
   geometry is set.
5. **FL is in the geometry-dominated regime** for current f_prism values. Don't bother
   tuning it unless geometry-only closure plateaus.

**What's DIFFERENT for HDD (don't blindly copy PDD values):**

1. **Effective coupling target.** PDD aims at ~74% to match LILAC's ~68%. HDD reference
   (e.g., Vulcan publications, or whichever code Xcimer compares against) likely has a
   different absorbed-energy fraction -- ascertain it first.
2. **Geometry baseline.** VI_6 starts at cone=20° (per CLAUDE.md table). Going wider
   to defocus may require more aggressive cone changes (e.g. 25–30°) than PDD's 2° step.
3. **Hot-spot mass / cold-shell mass.** VI_6 reported imploded DT mass is currently
   the order ~1.6 mg per CLAUDE.md notes (verify against the digitized HDD reference
   when added to `Aux_Standalone_files/`). HDD likely does NOT have the target-mass
   deficit that limited PDD's HS ρR (which we now know was ~37% of the
   Olson-2021-Fig-7-integrated 1.5–1.6 mg reference). If VI_6's Helios imploded mass
   matches its reference within 20%, full HS ρR closure should be possible. If
   not, the same structural-deficit ceiling applies.
4. **Over-drive magnitude.** Per CLAUDE.md PDD scan, VI_6 was at V=763 km/s vs reference
   410 km/s -- 86% over, vs PDD_22's 43% over. HDD over-drive is roughly 2× worse than
   PDD's was. Plan for a larger coupling reduction (target ~50–60%?) to land on reference V.
5. **HDD pulse shape** (610 TW peak vs PDD's 329 TW) and **target geometry** (~0.5 cm R
   vs 0.23 cm) put it in a different ablation regime. Critical-surface intensity is
   ~2.8e15 W/cm² vs PDD's 1.5e15 -- 2× higher. Resonant / SBS / parametric thresholds
   are different.

**Suggested HDD calibration sweep (drop-in template):**

```bash
# After identifying VI_6 baseline geometry from the .rhw file, sweep:
#
#   {cone in [20, 25, 30, 35]} × {spot in [base, base+0.05, base+0.10] cm}
#
# Target: eff_avg_coupling lands in the same fractional position relative
# to HDD-reference-coupling as PDD's 74%/68% ratio implied for PDD.
#
# Run scan_summary after each batch:
python3 ~/helios_postprocess/examples/scan_summary.py \
    ~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/ \
    --out ~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/HDD_calibration_scan.csv \
    --lilac-velocity 410   # adjust to HDD reference velocity
```

**Diagnostic priorities for HDD (same toolset, retargeted):**
- `examples/compare_absorbed_energy.py` -- compare candidate HDD geometry to VI_6 baseline.
- `examples/plot_shock_trajectories.py` -- expect 3-shock train if the HDD pulse design is
  multi-shock (most are); verify it survives at the calibrated geometry.
- `examples/check_flux_limiter.py` -- if HDD's `r_crit` sits in a different conduction
  regime than PDD's, the FL knob may engage differently. Worth one diagnostic run.
- `examples/scan_summary.py` -- once 4+ HDD runs exist, the same composite-distance
  ranking will surface the production calibration point.

**Pre-existing HDD scan in CLAUDE.md** (the "HDD Calibration Setup" / VI_6 sections
elsewhere in this doc) covers flux-limiter sensitivity in slab geometry. The PDD-closure
work suggests revisiting whether the flux-limiter scan was in the FL-dominated or
geometry-dominated regime, and whether HDD's calibration mismatch is also primarily
geometric.

## Test Data

| Case | Regions | Zones | Key Features |
|------|---------|-------|-------------|
| **Olson_PDD_9** | 4 | 351 | 3-component pressure, PDD, igniting, primary validation target |
| **VI_6** | 4 | 350 | Vulcan HDD target, over-ablation issue, work in progress |

### Olson_PDD_9 Region Structure
- Region 1: DT Vapor (zones 0-150) -- hot spot
- Region 2: DT Solid (zones 151-190) -- cryo DT ice
- Region 3: DT-CH foam (zones 191-320) -- wetted foam ablator
- Region 4: CH Skin (zones 321-350) -- outer ablator

### Validated Reference Values (Olson_PDD_9)
- Stagnation time: 12.599 ns (min HS radius = 0.0068 cm)
- Peak density: 181.17 g/cc (at t = 12.670 ns)
- Bang time: 12.681 ns
- Hot spot pressure: 106.97 Gbar
- Hot spot internal energy: 597.76 kJ
- Core radius: 0.0068 cm
- Stagnation CR: 29.6 (R0=0.2008 cm, R_hs=0.0068 cm) -- validated April 2026
- In-flight CR: reported as 6.3 (known bug -- see Open Items)
- Peak implosion velocity: 504 km/s (shell zones, pre-bang) -- validated April 2026
  NOTE: LILAC reference is 410 km/s; Helios over-drives by ~23% even for PDD_9
  (different pulse shape than calibration runs -- PDD_9 uses lower-energy pulse)
- Adiabat: 3.03 (cold fuel at peak v_imp)
- Fusion yield: 20.594 MJ, Target gain: 9.574
- <Ti>_n: 22.46 keV, <P>_n: 193.40 Gbar
- <rhoR>_n (fuel): 0.5189 g/cm2
- Burn-averaged: T=23.75 keV, P=208.4 Gbar, rhoR=0.5515 g/cm2
- Initial DT mass: 4.718 mg, Initial ablator mass: 0.349 mg
- Unablated fuel: 0.1262 (12.6%), Unablated ablator: 0.000
- Stagnated fuel: 0.1219, Imploded DT: ~0.60 mg

### VI_6 Reference Values (work in progress)
- Laser energy: 3.383 MJ (absorbed); published design is 4.0 MJ
- Stagnation time: 14.630 ns
- Peak implosion velocity: 763.1 km/s (published 410 km/s -- over-ablation)
- In-flight KE: 311.6 kJ, Hydro efficiency: 9.2%
- IFAR: 18.1 (density-based)
- Fusion yield: 28.811 MJ, Target gain: 8.516

#### Pulse shape (measured from EXODUS, April 2026)
- Peak power: 610 TW at t = 10.40 ns
- Laser off (P < 1% of peak): t ~= 14.30 ns -- 0.33 ns before stagnation
- Simulation time range: 0.00 - 16.00 ns over 296 EXODUS timesteps
- Corona R_outer at peak power: 0.91 cm (5x initial target radius of 0.19 cm)
- Incident intensity at peak power: I_outer = 58.7 TW/cm^2 at R_outer = 0.91 cm
- Peak intensity at critical surface: ~6.3 x 10^14 W/cm^2 at R ~= 0.22 cm
  (amplified by (R_out/r_crit)^2 ~= 17 from geometric convergence)
- Optical depth to peak-I location: tau ~= 0.44 at peak power
- **NOTE**: EXODUS timesteps are clustered around stagnation/burn. Temporal midpoint
  of EXODUS output (idx = nt/2) is ~0.1 ns before stagnation, when laser is already
  off. Always use `np.argmax(P_on_target)` to pick peak-drive timestep for
  laser-coupling diagnostics.

#### Critical-surface trajectory (measured via plot_laser_intensity.py, April 2026)
Using `elec_density` crossing (primary method, confirmed April 2026):

| t (ns) | P (TW) | r_crit (um) | R_outer (cm) | R_out/r_crit |
|---|---|---|---|---|
| 4.50 | 65 | 214 | 0.35 | 16 |
| 6.40 | 65 | 209 | 0.45 | 22 |
| 8.40 | 295 | 198 | 0.60 | 30 |
| 10.30 | 610 | 176 | 0.80 | 45 |
| 12.30 | 610 | 129 | 1.15 | 89 |
| 14.30 | 560 | 64 | 1.50 | 234 |
| 14.63 (stag) | 0 | ~46 | -- | -- |

Critical surface shrinks 214 -> 64 um during drive (3.3x), then converges another
1.4x to stagnation radius. Peak coronal intensity occurs at **laser turn-off**, not
at peak power, because (R_out/r_crit)^2 grows faster than P_delivered falls.

#### Confirmed intensity metrics (plot_laser_intensity.py, April 2026)
- Peak I_outer                = 6.63e13 W/cm^2
- Peak I (Method 2, max over t,r) = 2.83e15 W/cm^2 (at late-time r_crit)
- Peak I at r_crit (time history peak) = 2.01e15 W/cm^2
- These match the 0.22 cm snapshot analysis (6.3e14 at peak power) after
  accounting for the late-time geometric convergence boost.

## EXODUS Variable Reference

Real variable names as they appear in Helios EXODUS output, with array shapes and
conventions. Internal postprocessor names (used after `data_builder.py`) differ --
this section documents the *raw EXODUS* names needed when writing new diagnostic
scripts that go straight to the .exo file.

### Geometry
| EXODUS name | Shape | Units | Notes |
|---|---|---|---|
| `time_whole` | (nt,) | s | Multiply by 1e9 for ns |
| `zone_boundaries` | (nt, nzone+1) | cm | Zone boundary radii. Internal alias: `zbnd` |
| `radius_squared` | (nt, nzone+1) | cm^2 | Redundant with `zone_boundaries` |
| `zone_mass` | (nt, nzone) | g | Per-zone mass (Lagrangian, time-constant) |

### Plasma state
| EXODUS name | Shape | Units | Notes |
|---|---|---|---|
| `elec_density` | (nt, nzone) | # per cm^3 | Electron number density. Use for n_crit crossing. |
| `ion_density` | (nt, nzone) | # per cm^3 | Ion number density |
| `mean_charge` | (nt, nzone) | dimensionless | Zbar; self-consistency: ne ~= Zbar * ni |
| `mass_density` | (nt, nzone) | g per cm^3 | Internal alias: `rho` |
| `elec_temperature` | (nt, nzone) | eV | Electron temperature |
| `elec_pressure` | (nt, nzone) | J per cm^3 | Electron pressure (x 1e-8 for Gbar) |

NOTE: the variable name is lowercase `elec_density`, NOT `NumElecDensity`.
An older CLAUDE.md note claimed `NumElecDensity` was absent from Xcimer EXODUS
files; that name was never correct. `elec_density` is present in all Helios
outputs verified to date (PDD_9, VI_6, PDD_26b).

### Laser
| EXODUS name | Shape | Units | Notes |
|---|---|---|---|
| `LaserPwrOnTargetForBeam` | (nt, nbeam) | W | Power reaching outer target surface |
| `LaserPwrDeliveredForBeam` | (nt, nbeam) | W | Power delivered by laser system |
| `LaserPwrDensForBeamAtZoneBd` | (nt, nbeam, nzone+2) | W/cm^3 | Deposition at zone boundaries |
| `LaserPwrDensTotForBeamAtZoneBd` | (nt, nbeam, nzone+2) | W/cm^3 | All-beam sum |
| `LaserPwrSrc` | (nt, nzone) | W/cm^3 | Volumetric deposition at zone centers |
| `LaserEnTimeIntg` | (nt, nzone) | J | Cumulative deposited energy per zone |
| `LaserEnDeliveredTimeInt` | (nt,) | J | Cumulative delivered |
| `EnLaserDepositedTimeIntg` | (nt,) | J | Cumulative absorbed |
| `laserAttinuationCoeff` | (nt, nbeam, nzone+2) | 1/cm | IB absorption coefficient; see conventions below |

### laserAttinuationCoeff conventions (CRITICAL -- get these wrong and diagnostics silently fail)

1. **Layout**: shape is `(nt, nbeam, nzone+2)` = 1 leading ghost at r=0 + (nzone+1) physical
   boundaries. To match `zone_boundaries` shape `(nt, nzone+1)`, drop the LEADING point:
   ```python
   if att.shape[-1] == zone_boundaries.shape[-1] + 1:
       att = att[:, :, 1:]          # NOT att[:, :, :-1]
   ```
   Verified by inspecting `att[idx, 0, 0]` which is 0.0 (unphysical ghost at r=0),
   with trailing boundary values being real physical attenuation in the corona.

2. **Opaque-zone sentinel**: Helios uses `alpha = 1e30` to flag zones the laser
   cannot propagate into (past the critical surface, or outside the first beam's
   reach working inward). These are NOT numerical junk -- they're deliberate flags.
   Clip before use in `exp(-tau)` or any sum:
   ```python
   att = np.where(att > 1e20, 1.0e6, att)   # 1e6 cm^-1 saturates tau cleanly
   ```
   Threshold 1e20 is well above any physical alpha (typical range 1e-7 to ~1e3 cm^-1)
   and well below the 1e30 flag.

3. **Center-averaging**: EXODUS stores alpha at zone *boundaries*, but `LaserPwrSrc`
   is at zone *centers*. For Method 1 (`I = P_src / alpha`) use zone-centered alpha:
   ```python
   alpha_zone = 0.5 * (att[t, 0, :-1] + att[t, 0, 1:])
   ```

### Derived quantities
| Derived quantity | Formula |
|---|---|
| Zone centers | `0.5 * (zbnd[:, :-1] + zbnd[:, 1:])` |
| Zone widths | `np.diff(zbnd, axis=1)` |
| Critical density | `n_crit = 1.115e21 / lambda_um^2`  [cm^-3] |
| Incident intensity | `I_outer = LaserPwrOnTargetForBeam / (4*pi*R_outer^2)` |
| Critical radius | Outermost r where `elec_density >= n_crit` (primary); OR first zone outside `att >= 1e20` opaque sentinel (fallback) |

### Timestep selection for laser diagnostics
EXODUS output is **not uniformly spaced in time**. Helios concentrates timesteps
around stagnation and burn (highest-gradient physics). For a typical implosion
simulation:
- First ~1/3 of EXODUS timesteps: drive phase (laser on)
- Last ~2/3: stagnation, burn, disassembly (laser OFF)

**Rule for laser-coupling diagnostics**: pick `idx = int(np.argmax(P_on_target))`.
NEVER use `idx = nt // 2` -- the temporal midpoint is almost always past laser-off.

## Simulation Paths

### Mac Studio
- Olson PDD: `~/Sims/Xcimer/Olson_PDD/Olson_PDD_<N>/Olson_PDD_<N>`
- VI_6: `~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6/VI_6`

### MacBook
- Olson PDD_9: `~/Sims/Helios_Sims/Xcimer_Sims/Olson_PDD/Olson_PDD_9/Olson_PDD_9`

## Workflow

1. Edit on MacBook (`python`) -> push to GitHub
2. Pull on Mac Studio (`python3`) -> run against simulation data
3. Terminal commands run **one line at a time** (avoid multi-line paste errors)
4. Validate against known reference values before moving to next feature
5. Output goes to `<sim>_summary.txt` -- use `grep` on the file, not stdout
6. Mac Studio edits: use Python file manipulation, not `sed` (macOS sed newline issues)

## Quick Tests

```bash
# Mac Studio -- Olson (primary validation target)
python3 ~/helios_postprocess/examples/run_analysis.py \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9

# Validation checks
grep "Stagnation CR\|In-flight CR\|Peak implosion" \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9_summary.txt
grep -A 6 "MASS FRACTIONS" \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9_summary.txt
grep -A 8 "LASER CONFIGURATION" \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9_summary.txt
grep -A 25 "COMPARISON WITH PUBLISHED" \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9_summary.txt
```

## Laser Intensity Diagnostics

Reconstructing I(r, t) from EXODUS via two independent methods. Cross-checking them
is the main way we diagnose whether Helios's absorption behaves as expected.

### Method 1 (exact, local)
```
I(r) = LaserPwrSrc(r) / alpha_zone(r)
```
Derived from the definition P_src = alpha * I. Undefined where alpha <= 0 or where
alpha is at the opaque-sentinel ceiling. Only meaningful in zones where absorption
is physically active -- filter on `alpha_zone * dr > 1e-4` (zones contributing at
least ~0.01% to total tau). Without this filter, Method 1 blows up in the tenuous
outer corona where alpha is ~1e-7 cm^-1 and the ratio becomes numerical noise.

### Method 2 (Beer-Lambert, spherical)
```
I(r, t) = I_outer(t) * (R_out(t) / r)^2 * exp(-tau(r, t))
tau(r, t) = integral from r to R_out of alpha(r') dr'
I_outer(t) = LaserPwrOnTargetForBeam(t) / (4 * pi * R_out(t)^2)
```
Defined everywhere. The `(R_out/r)^2` term is geometric convergence (critical for
spherical geometry -- omitting it silently undercounts peak intensity by ~20x in
HDD-scale coronas). Integrate alpha*dr inward from the outer boundary using
`np.cumsum(...[::-1])[::-1]`.

### Expected agreement
In zones where alpha > 0 and IB is the only absorption mechanism, Method 1 and
Method 2 should agree to within zone-discretization accuracy (~10-20%). Systematic
disagreement is diagnostic:
- **M1 > M2**: `LaserPwrSrc` contains deposition beyond what IB alpha captures
  (e.g. resonant absorption at n_crit, or a dump of unabsorbed power at critical
  surface that Helios attributes to a specific zone).
- **M1 < M2**: typically a center-vs-boundary averaging inconsistency; usually fixable.

### VI_6 at peak power (verified April 2026, t = 10.4 ns, P = 610 TW)
- R_outer = 0.91 cm, I_outer = 58.7 TW/cm^2
- Critical surface at r ~= 0.22 cm (first transparent zone working inward)
- Peak Method-2 intensity = 633 TW/cm^2 = 6.3e14 W/cm^2 at r ~= 0.22 cm
- tau at peak-I location ~= 0.44 -- ~35% absorbed before reaching critical surface
- Geometric amplification (R_out/r_crit)^2 ~= 17

### Peak coronal intensity across targets (ne-based r_crit, April 2026)

First three targets run with validated `elec_density` r_crit (not sentinel fallback):

| Target | P_peak (TW) | I_outer (W/cm^2) | I_M2 peak (W/cm^2) | I at r_crit (W/cm^2) | I/TW at r_crit |
|---|---|---|---|---|---|
| Olson_PDD_9 | 315 | 5.07e13 | 2.05e15 | 2.05e15 | 6.50e12 |
| Olson_PDD_26b | 329 | 4.03e13 | 1.55e15 | 1.46e15 | 4.45e12 |
| VI_6 (HDD) | 610 | 6.63e13 | 2.83e15 | 2.65e15 | 4.35e12 |

Two effects visible:

**Beam-geometry effect (PDD_9 vs PDD_26b):** Both are Olson PDD targets at
essentially the same peak power, but PDD_9 (pencil beam: cone=1 deg, spot=0)
delivers 1.46x higher I/TW at r_crit than PDD_26b (realistic NIF PDD: cone=20 deg,
spot=0.16 cm Gaussian). Normal-incidence pencil-beam rays deposit all their power
concentrated at the strict critical surface; oblique cone rays spread deposition
and turn at ne = n_crit * cos^2(theta). The diagnostic is sensitive to a real
geometric difference between the two laser-source specs, validating the
ne-based method.

**Turning-point offset (M2-peak vs at-r_crit):** Gap is 0% for PDD_9 (pencil
beam, M2 peak coincides with strict r_crit), 6% for PDD_26b (cone=20 deg),
7% for VI_6 (realistic cone). Consistent with oblique-incidence rays turning
slightly outside strict critical density; M2-peak sits just outside r_crit
where (R_out/r)^2 geometric amplification hasn't yet been killed by exp(-tau)
at the turning region.

**HDD coronal deposition null result:** PDD_26b and VI_6 have statistically
identical I/TW at r_crit (4.45e12 vs 4.35e12 W/cm^2/TW, 2% agreement).
Per-incident-watt coronal coupling is the same for PDD and HDD geometries at
lambda=0.351 um. HDD over-drive is therefore NOT in the coronal deposition
channel -- absorption per delivered watt is consistent between drive types.
The over-drive has to be downstream, in ablation-pressure conversion or
thermal conduction structure. Consistent with the existing note that hydro
efficiency is locked near 10% regardless of geometry.

NOTE (refined April 2026): the 4.45e12 vs 4.35e12 "cluster" is NOT a
universal Helios property -- both PDD_26b and VI_6 use forward-focused
realistic-cone geometry. DM_01b (cone=3.5 deg, d=-2.5 cm, inverted focus)
gives I/TW at r_crit = 2.77e12 W/cm^2/TW, 38% below the cluster. Per-watt
coronal coupling depends on beam geometry, not just drive type; the null
result applies between targets with similar effective incidence angle.

### Two paths to reference velocity (April 2026)

DM_01b (inverted focus, d=-2.5 cm, cone=3.5 deg, 230 TW, 1.365 MJ) matches
LILAC reference velocity (405 vs 370 km/s, +9.4%) via a different mechanism
than PDD_26b (forward focus, d=0.20 cm, cone=20 deg, 329 TW, 2.15 MJ, matches
470 km/s within 2%). DM_01b has LOWER coronal I/TW but HIGHER fractional
absorption (91% vs 87%) -- diverging-beam geometry spreads absorption over
a larger coronal volume, with longer ray paths at lower local intensity.

Both runs hit reference velocity; both MISS reference ablation-physics
metrics by the same rough amount (adiabat ~1.5 vs 7.4; IFAR ~3 vs 18;
imploded DT ~0.60 vs 1.68 mg). The shared residual across two different
geometric paths to velocity match is the strongest evidence to date that
the ablation-physics deficit is not a laser-coupling issue. Flux limiter
and EOS remain the primary calibration levers per existing HDD Calibration
Setup notes.

DM_01b also has a much larger turning-point offset (M2 peak / at-r_crit gap
= 29%) than the forward-focused targets (0-7%). Effective local incidence
angle -- not nominal cone -- drives the offset, consistent with rays
fanning out over 2.5 cm before hitting plasma.

Peak-coronal-intensity table (updated with DM_01b):

| Target | Cone | Focus d | P_peak (TW) | I at r_crit (W/cm^2) | I/TW at r_crit | M2/r_crit gap |
|---|---|---|---|---|---|---|
| PDD_9 | 1 deg | 0 | 315 | 2.05e15 | 6.50e12 | 0% |
| PDD_26b | 20 deg | 0.20 | 329 | 1.46e15 | 4.45e12 | 6% |
| VI_6 | realistic | - | 610 | 2.65e15 | 4.35e12 | 7% |
| DM_01b | 3.5 deg | -2.50 | 230 | 6.36e14 | 2.77e12 | 29% |

### Diagnostic Scripts (repo root)
| Script | Purpose |
|---|---|
| `plot_adiabat_shock.py` | Adiabat history + first shock plots (any target with DT ice layer) |
| `plot_laser_deposition.py` | Laser energy deposition spatial profile at multiple timesteps |
| `helios_postprocessor_guide.docx` | User guide for collaborators |

*(`plot_laser_intensity.py` retired April 2026 -- fully integrated into pipeline; see `helios_postprocess/laser_intensity.py` and the `analyze_laser_intensity` phase in `ICFAnalyzer`.)*

Standard usage (all scripts accept .exo path as first positional arg):
```bash
# Standalone diagnostic scripts (still supported):
python3 ~/helios_postprocess/plot_adiabat_shock.py \
  ~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6/VI_6.exo
python3 ~/helios_postprocess/plot_laser_deposition.py \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_2021_01a/Olson_PDD_2021_01a.exo

# Laser intensity now part of the main pipeline -- no separate script needed:
python3 ~/helios_postprocess/examples/run_analysis.py \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_26b_burn/Olson_PDD_26b_burn
# produces <base>_report.pdf with 3 intensity pages,
#          <base>_summary.txt with LASER INTENSITY section,
#          <base>_comparison.pdf with intensity rows (if _published.json
#            includes "I_at_crit_peak_Wcm2")
```

### Critical-surface location as cross-check
Two independent ways to locate the critical surface:
1. Electron density (primary): outermost r where `elec_density >= n_crit`
2. Attenuation sentinel (fallback): first zone outside opaque-sentinel region (`att >= 1e20`)

Sentinel fallback validated against `elec_density` method on Olson_PDD_9:
agreement within 1 zone during laser-on (median offset -0.50 zones,
max 2.31 zones at t=2.2 ns during early corona formation before a clean
coronal shelf is established). Safe to use as a fallback when
`elec_density` is not available, though no Helios output verified to date
lacks it. Post-laser-off (t > ~13 ns for PDD_9), r_crit from EITHER
method is NOT physically meaningful -- sentinel field is stale and ne
tracks disassembling stagnation region. Mask laser-on in any
r_crit-trajectory plot.

### Script implementation lessons (April 2026)
Hard-won during debugging the April 2026 VI_6 run:

1. **r_crit primary method is `elec_density`**; sentinel fallback (`att < 1e20`)
   is a safe substitute. In the fallback, use `opaque.max() + 1` (outermost
   opaque zone, step one outward), NOT `transparent.min()`. The 1e30 sentinel
   tags only the turning-point region, not every inner zone, so the
   "innermost transparent zone" is some arbitrary central zone well below
   r_crit. VI_6 observed failure (before `elec_density` was loaded correctly):
   `transparent.min()` returned 0.026 cm vs. actual r_crit of 0.2 cm (8x off).

2. **Method 1 (P_src/alpha) filter**: alpha > 1% of per-timestep max alpha.
   Without this filter, M1 returns 10^18+ W/cm^2 in the tenuous outer corona
   where both alpha (~1e-7 cm^-1) and P_src are noise-dominated. The 1% cutoff
   keeps the absorbing-layer signal (where IB alpha is significant) and masks
   the numerical-noise region.

3. **Y-axis auto-scaling**: matplotlib's semilogy will span 300+ decades when
   exp(-tau) underflows inside the opaque region (giving I ~= 10^-300 values).
   Always clip with `ax.set_ylim(y_hi/1e8, y_hi)` using finite data only.

4. **I-at-r_crit lookup**: use direct zone-index lookup, NOT np.interp across
   the exp(-tau) cliff. np.interp produces garbage when interpolating a curve
   that drops 280 decades in one zone.

5. **Peak coronal intensity grows through the pulse**, not at peak power. The
   (R_out/r_crit)^2 geometric factor accelerates as r_crit shrinks, outpacing
   the drop in delivered power after peak. Relevant for parametric-instability
   threshold analysis.

### Pipeline integration (April 2026 -- Task 2)

All logic above lives in the automated pipeline. The standalone
`plot_laser_intensity.py` has been retired; reconstruction runs as part
of `run_analysis.py`.

| Layer | Contribution |
|---|---|
| `data_builder.py` | Loads `laser_attenuation_coeff` + `laser_power_on_target` via `_VARIABLE_MAP`; squeezes beam axis with labeled log line |
| `helios_postprocess/laser_intensity.py` | Module with `clean_attenuation`, `compute_method1`, `compute_method2`, `find_critical_radius_*`, `analyze_laser_intensity` entry point |
| `ICFAnalyzer.analyze_laser_intensity()` | Called after `analyze_drive_phase`; populates scalars and histories on `data`; caches 2D arrays in `_laser_intensity_arrays` for the plotter |
| `ICFPlotter._plot_laser_intensity()` | 3 PDF pages (log10 I(r,t) heatmap with r_crit overlay; P_laser + I histories; Method 1 vs 2 cross-check at peak power) |
| `ICFOutputGenerator` | `LASER INTENSITY` section in `<name>_summary.txt`, between `LASER CONFIGURATION` and `EOS MODELS` |
| `burn_averaged_metrics.py` | `I_at_crit_peak_Wcm2` and `I_grid_outer_peak_Wcm2` keys plumbed through `histories` -> `sim_metrics` -> `compare_with_published` |

**Renamed attribute:** `I_outer_*` -> `I_grid_outer_*` throughout the
pipeline code. The simulation grid extends well into vacuum (e.g., 0.8
cm while the capsule outer radius is 0.2 cm), so `I = P / (4*pi*R_grid^2)`
is much lower than intensity at the capsule / critical surface.
`I at critical surface` remains the primary physical quantity;
grid-outer is reported for audit purposes only. Inside
`laser_intensity.py` the local variable is still named `I_outer` since
it is the incident ray at the grid outer boundary (Method 2
Beer-Lambert's starting point); all attributes, log messages, and
summary output use `I_grid_outer`.

**M1 filter change:** the old standalone script's threshold
(`alpha > 0.01 * alpha_max_t`) cut at ~100 cm^-1 at the critical
surface, masking the entire absorbing layer. Replaced with an absolute
floor `ALPHA_MIN_M1 = 1e-2 cm^-1` in
`helios_postprocess/laser_intensity.py`. Coronal-noise exclusion is
preserved without losing the absorbing layer itself.

**Published-JSON keys:**

    "I_at_crit_peak_Wcm2":     [value, unc]   # primary comparison target
    "I_grid_outer_peak_Wcm2":  [0.0, 0.0]     # leave at 0; geometry-dependent,
                                              # not cross-code comparable

**Graceful degradation:** if a simulation lacks `LaserPwrOnTargetForBeam`
or `laserAttinuationCoeff` (pre-upgrade `.exo`), `analyze_laser_intensity`
logs a warning and `None`-populates attributes. Plotter emits a
placeholder page; summary block is omitted; compare rows are filtered by
the zero-skip rule. No crashes.

**Cross-check:** on `Olson_PDD_26b_burn`, the pipeline's peak-power
r_crit (0.099 cm / 992 um) sits inside the drive-phase formula-method
1-sigma band (0.079 +/- 0.039 cm), consistent with r_crit expanding
outward at peak drive relative to the all-timesteps mean.

### Shock train pipeline (May 2026 — Task 3 Stage 3 multi-shock)

Foot/ramp/peak gas/ice breakouts extracted automatically and propagated
through summary, comparison, and report PDF. Used to diagnose Helios
over-drive vs LILAC's published 7.5 / 10 / 13 ns arrival pattern, and
to inform foot-pulse calibration. **Shock timing drives the HS ρR
closure work; treat the SHOCK TRAIN block as a first-class diagnostic.**

| Layer | Contribution |
|---|---|
| `helios_postprocess/pressure_gradients.py` | `track_shock_trajectories()` (greedy nearest-radius linker with inward-only motion, velocity bound, max_gap_steps, terminate-at-breakout) and `consolidate_breakouts()` (merges raw breakouts within `min_separation_ns` into foot/ramp/peak events). |
| `ICFAnalyzer._compute_shock_train()` | Runs inside `analyze_implosion_phase`. Populates `data.shock_trajectories`, `data.shock_coalescence_events`, `data.shock_breakouts`, `data.shock_events` (consolidated) and class-keyed scalars `t_foot_shock_ns`, `t_ramp_shock_ns`, `t_peak_shock_ns` (NaN when not detected). Prints `[shock_train]` diagnostic table sampled uniformly in time. |
| `ICFPlotter._plot_shock_train()` | 1-page R-T overlay (top, 2/3) + inward-velocity panel (bottom, 1/3). Trajectories colored by class (`foot=blue, ramp=orange, peak=red`); pile-up gray. Headline summary box in bottom-left. |
| `ICFOutputGenerator` | New `SHOCK TRAIN` block between `LASER INTENSITY` and `EOS MODELS` — trajectory count, consolidated event count, per-event table (class, t [ns], r [µm], P_post [Mbar], P_ratio, merged-raw count). |
| `burn_averaged_metrics.py` | `t_foot_shock_breakout_ns`, `t_ramp_shock_breakout_ns`, `t_peak_shock_breakout_ns` plumbed through `histories → sim_metrics → compare_with_published`. NaN sim values mapped to `-1.0` so the existing `<= 0 == missing` skip logic in the comparison loop kicks in. |
| `examples/plot_shock_trajectories.py` | Standalone reproduction (R-T + velocity + CSV emission + `[shock_summary]` grep tag) for ad-hoc runs or scans not going through `run_analysis.py`. |

**Published-JSON keys** (place next to `<base>.exo` as `<base>_published.json`):

    "t_foot_shock_breakout_ns": [7.5,  0.0],   # LILAC reference for foot shock
    "t_ramp_shock_breakout_ns": [10.0, 0.0],
    "t_peak_shock_breakout_ns": [13.0, 0.0]    # leave as [0.0, 0.0] to skip

The comparison row prints "—" instead of a number when the sim didn't
detect that shock (NaN → −1.0 internal sentinel).

**Tuning knobs (via `ICFAnalyzer` config dict):**

    shock_train_dP_dr_threshold            5e7      # J/cm⁴; raise to filter weak shocks
    shock_train_min_P_ratio                1.5      # compression-ratio gate per detection
    shock_train_group_separation           1e-2     # cm; identify_shocks grouping width
    shock_train_max_gap_steps              20       # trajectory persistence across detection dropouts
    shock_train_min_traj_len               5        # post-filter (bypassed for breakout-ended)
    shock_train_min_traj_span              5e-3     # cm; drop stationary "trajectories"
    shock_train_min_breakout_separation_ns 0.3      # consolidate raw breakouts within this window

**Headline measurements (PDD_20_fab02_foot25_s016_burn, α=1.05):**

    foot:  Helios 5.65 ns  vs LILAC 7.5 ns   (−1.85 ns over-drive)
    ramp:  Helios 8.25 ns  vs LILAC 10  ns   (−1.75 ns over-drive)
    peak:  Helios —        vs LILAC 13  ns   (absorbed into convergence)

Consistent ~1.8 ns over-drive across both detected shocks → single
ablation-physics mechanism, not a per-shock issue. Same lever that
fixes shock timing should also reduce HS pressure (270 → 90 Gbar).

**Graceful degradation:** if the search window is degenerate (gas
cavity too thin or n_capsule_regions < 2), `_compute_shock_train`
logs a skip and leaves all four `data.shock_*` containers at their
empty defaults; the SHOCK TRAIN summary block is omitted and the
plotter emits no page (rather than crashing on empty arrays).

### Priority 1 -- cr_inflight (FIXED April 2026)
- Now uses ablation_front_radius[pv_idx]. Value ~4.5-4.6 is physically correct.
- Old reference of 2.16 was from wrong (pre-alpha-onset-fix) timestep.

- Currently uses HS boundary radius at peak velocity: R0/R_hs = 6.3 for PDD_9
- CLAUDE.md convention 12 specifies R0 / R_ablfront_at_peak_velocity = 2.16
- The ablation front radius at peak v (~0.0930 cm) is the correct denominator
- Fix: change Rf line in compute_performance_metrics() to use
  `ablation_front_radius[pv_idx]` not `zone_boundaries[pv_idx, ri[pv_idx,0]]`

### Priority 2 -- Hot-spot properties (FIXED April 2026)
- _compute_hot_spot_properties() now uses ri[stag_idx,0] instead of T mask.
- Validated: PDD_9 radius=0.0068 cm, pressure=208 Gbar, IE=40.75 kJ.

### Priority 3 -- PDD calibration (active, cone=20 deg geometry)
- 26b: CR=29.6, imploded DT=0.60 mg, v=478 km/s -- geometry matched
- Adiabat 1.43 vs 3.0 is sole remaining deficit
- 26c running: cone=20, foot=45-50 TW -- mapping foot/adiabat lever
- If foot lever insufficient: time-dependent peak power reduction (~0.75x after 9 ns)
- Longer term: obtain LILAC absorbed-energy history to anchor empirical correction

### Priority 4 -- Physics Module Integration
- `energetics`, `neutron_downscatter`, `pressure_gradients` work standalone
  but are not yet wired into `ICFAnalyzer` or `ICFPlotter`.
- Add `ICFAnalyzer.integrate_physics_modules()`, plotter pages, output sections.

### Priority 5 -- data_builder.py cleanup
- Duplicate laser wiring block exists (~lines 286-305).
  Both blocks assign same values -- harmless but should be cleaned up.

### Priority 6 -- Laser diagnostic workstream (ACTIVE April 2026)
Ongoing analysis using `plot_laser_intensity.py` and EXODUS direct access.

**COMPLETED April 2026**:
- VI_6 full pipeline: r_crit trajectory extracted, peak intensities confirmed
  (I_outer 6.6e13, peak coronal 2.8e15, at r_crit 2.0e15 W/cm^2).
- Script gotchas identified and fixed: r_crit fallback logic, M1 filter,
  y-axis clipping. See Laser Intensity Diagnostics > Script lessons.

**NEXT**:
1. **Olson PDD comparison suite**: run plot_laser_intensity.py on three
   targets to contrast with VI_6:
   - Olson_PDD_9 (4-region igniting baseline, LILAC comparison target).
     First check whether NumElecDensity is present -- if so, cross-check
     the `att < 1e20` fallback against the direct ne-based r_crit.
   - Olson_PDD_2021_01a (230 TW peak, 1.44 MJ, Olson 2021 paper baseline).
     Lower intensity; direct comparison of deposition profile shape.
   - Olson_PDD_26b (BEST MATCH PDD calibration, cone=20 deg).
     Key comparison for HDD over-drive diagnosis: is peak I at r_crit
     similar to VI_6 (=> flux-limiter signature) or different
     (=> other physics in play)?

2. **M1 vs M2 comparison across targets**: on VI_6 pages 1/3, Method 1 is
   largely masked by the 1% filter leaving only a narrow band at r_crit.
   Check whether M1/M2 ~= 1 holds for Olson targets too, and whether the
   outer-corona M1 explosion observed in early VI_6 diagnostics (10^19
   W/cm^2) is VI_6-specific or universal to Helios.

3. **Critical-surface cross-check (if ne available)**: compare
   `att < 1e20` r_crit vs `ne >= n_crit` r_crit on any target that has
   NumElecDensity in EXODUS. Agreement to ~1 zone validates the fallback.

4. **Flux limiter scan (DT ice slab test problem)**: set up planar-equivalent
   DT slab (R=0.5 cm, 100 um DT ice, 25 TW flat foot, burn off) at f = 0.06,
   0.08, 0.10, 0.12. Decision rule: if d(adiabat)/d(f) < 0.2 over f=[0.06, 0.10],
   flux limiter alone cannot close HDD adiabat gap and we pivot to EOS/empirical.

5. **VI_6 at f=0.08**: HDD calibration run. Expect velocity DROP (not rise) and
   modest adiabat rise if flux limiter is a useful lever for HDD geometry.

Longer-term:
- Extend plot_laser_intensity.py with time-dependent view: I at r_crit(t), not
  just spatial snapshots. (DONE April 2026 -- added peak coronal I and I-at-r_crit
  to Page 2 of the PDF.)
- Integrate I-profile output into ICFAnalyzer as an optional page in the main PDF. (DONE April 2026 -- Task 2 Stage C.2, see `_plot_laser_intensity` in icf_plotting.py.)

## Implosion & Ablation Diagnostics (Task 3.A, April 2026)

Timing milestones and ablation-physics scalars. This request made by David Montgomery, an Xcimer collaborator
in this priority order:

1. First shock breakout from fuel to gas
2. Second shock breakout (same interface)
3. Subsequent shock breakouts (3+ shock designs)
4. Shock flash at r=0
5. Time of peak velocity
6. Ablation pressure in ablator (CH / foam phase)
7. Ablation pressure in fuel (when ablator fully consumed)

### Stage 1 (shipped April 2026)

Three additions wired through `ICFAnalyzer` and summary:

| Attribute on ICFRunData | Source | Exposure |
|---|---|---|
| `t_peak_velocity_ns` | `time[peak_velocity_index]` during `analyze_implosion_phase` | Augmented "Peak implosion velocity" log line; new "Peak velocity time" row in TIMING summary |
| `ablation_pressure_Gbar` | `(ion + rad pressure)[t, ablation_front_indices[t]] * 1e-8` for each timestep | History on `data`, available to plotters |
| `P_abl_peak_Gbar` | `max(ablation_pressure_Gbar)` over valid timesteps | New "Peak ablation pressure" log line at end of `_track_ablation_front` |

**PDD_26b_burn reference values (Stage 1 validation):**

    t_peak_velocity_ns          12.000 ns
    P_abl_peak_Gbar             176.69 Gbar

### Important semantic note on `P_abl_peak_Gbar`

Stage 1 measures **total pressure at the ablation-front zone**, where
the ablation front is defined as the steepest negative density gradient
outside the hot spot. This tracks the outer boundary of the DENSE shell,
which sits inside the compressed region and is closer to shell pressure
than to the classical Lindl momentum-balance ablation pressure.

Lindl's scaling `P_abl = 57 (I/1e14)^(2/3) (lambda/um)^(-2/3)` Mbar
at PDD_26b's peak critical-surface intensity (I_crit = 1.46e15 W/cm^2,
lambda = 0.35 um) predicts 6.77 Gbar -- about 26x lower than the
176.69 Gbar reported here. The difference is physical, not a bug:
shell pressure > ablation drive pressure by the time the shell is
converging. Stage 4 (material-split) will clarify the interpretation
by separating ablator-phase and fuel-phase averages.

### Stage 2-10 queued

- **Stage 2:** Debug `_compute_shock_breakout` -- currently returns
  zeros for PDD_26b due to `self.data.time * 1e9` double-multiplication
  (`data.time` is already in ns per `data_builder.py:344`).
- **Stage 3:** Extend shock breakout to N-shock detection, producing
  `shock_breakout_times_ns` list + per-shock scalars.
  **(SHIPPED May 2026 — see "Shock train pipeline" below.)**
- **Stage 4:** Material-split `P_abl` via
  `material_index[ablation_front_indices[t]]`; adds
  `P_abl_ablator_peak_Gbar`, `P_abl_fuel_peak_Gbar`,
  `t_fuel_ablation_start_ns`, `fuel_mass_ablated_mg`.
- **Stage 5:** Shock flash at r=0 -- peak pressure at innermost zone
  before stagnation -> `t_shock_flash_ns`.
- **Stage 6:** TIMING MILESTONES and ABLATION sections in summary text;
  comparison JSON keys for the new scalars.
- **Stage 7:** New ABLATION plotter page with material-phase shading.
- **Stage 8:** Wire `adiabat_history.py` module into `ICFAnalyzer`.
- **Stage 9:** Integrate `plot_adiabat_shock.py` into pipeline
  (overlaps with Stages 2-5).
- **Stage 10:** Integrate `plot_laser_deposition.py` into pipeline.

### Note on shock breakout bug

The `analyze_first_shock` call in `_compute_shock_breakout` converts
`self.data.time * 1e9`, which would be correct if `data.time` were in
seconds. But `data_builder.py` line 344 already multiplies by 1e9
during load, so `data.time` is in nanoseconds. The double-multiplication
makes every timestep ~1e10, which fails the foot-pulse mask
`(t_ns >= 4.0) & (t_ns <= 6.0)` trivially. Result: NaN breakout values,
which the error handler converts to 0.0. Stage 2 will fix this with a
one-line change.

## Dependencies

**Required**: numpy, scipy, matplotlib, netCDF4
**Optional**: scikit-learn (RANSAC shock fitting in icf_plotting.py -- guarded with `_HAS_SKLEARN`)

## Laser Configuration Parsing

`RHWParser._parse_laser_geometry()` extracts beam-1 ray-trace parameters and pulse shape.

### Fields on RHWConfiguration and ICFRunData:
| Field | Unit | Description |
|-------|------|-------------|
| `laser_wavelength_um` | um | Laser wavelength |
| `laser_spot_size_cm` | cm | Focal spot radius |
| `laser_half_cone_angle_deg` | deg | Half cone angle |
| `laser_focus_position_cm` | cm | Focus position from plasma origin |
| `laser_power_multiplier` | — | Power table multiplier |
| `laser_spatial_profile` | str | "Gaussian" or "Uniform" |
| `laser_foot_power_TW` | TW | Mean foot power (< 50% of peak) |
| `laser_peak_power_TW` | TW | Peak power |
| `laser_foot_start_ns` | ns | Foot start time |
| `laser_foot_end_ns` | ns | Foot end time |
| `laser_peak_start_ns` | ns | Peak start time |
| `laser_peak_end_ns` | ns | Peak end time |
| `laser_pulse_duration_ns` | ns | Total pulse duration |

## Comparison Framework

### Keys supported in `_published.json`:
All standard keys plus:
- `P_hs_ignition_Gbar` -- hot-spot pressure at ignition (rhoR_hs = 0.3 g/cm²)
- `hs_radius_ignition_um` -- hot-spot radius at ignition (micrometers)
- `peak_velocity_kms` -- peak implosion velocity
- `fraction_absorbed_pct` -- fraction of delivered laser energy absorbed

### Ignition criterion:
Hot-spot rhoR (absolute) >= 0.3 g/cm² (Olson et al. convention).
T_ion > 4.5 keV mask used only for burn propagation plot, not for scalar thresholds.

### PDF report improvements (April 12, 2026)
- path.simplify + threshold=0.5 added to icf_plotting.py (reduces vector complexity)
- DPI reduced 150->100 for all savefig calls
- Contour plots now off by default (include_contours=False in config)
- --contours CLI flag added to run_analysis.py via argparse (default off)
- MacBook must use "python" not "python3" (Anaconda vs system Python)

### run_analysis.py usage (updated April 2026)
    # Mac Studio
    python3 ~/helios_postprocess/examples/run_analysis.py <base_path>
    python3 ~/helios_postprocess/examples/run_analysis.py <base_path> --contours

    # MacBook
    python ~/Codes/helios_postprocess/examples/run_analysis.py <base_path>
    python ~/Codes/helios_postprocess/examples/run_analysis.py <base_path> --contours
## Session Update — April 26, 2026

This appendix captures changes from the April 2026 calibration session. Older
sections of this document may contain stale information that this section
supersedes; an integration pass to fold these notes into the main body is a
future TODO.

### New conventions (apply to all future analysis)

**DT ice adiabat is evaluated over Region 2 only.** For multi-region capsule
targets, both `_compute_adiabat` (peak-velocity) and `_compute_adiabat_at_breakout`
(base) take their zone selection as `ri[:, 0]:ri[:, 1]` — gas/fuel interface
through ice/foam interface — with NO density mask. Mass-weighted average over
the full ice layer regardless of partial-ablation state. This matches Olson 2021's
adiabat convention and is now apples-to-apples with published references. Earlier
foam-inclusive adiabat numbers (e.g., the WfCDT_01b 0.21 base adiabat) are
superseded.

**Plasma pressure = ion + electron pressure throughout.** All standard ICF
diagnostics — shock breakout, ablation pressure, adiabat (peak-vel and base),
hot-spot pressure (stagnation and ignition), neutron-averaged pressure,
burn-averaged pressure — use `data.plasma_pressure = ion_pressure + elec_pressure`
computed once in `data_builder.py`. Radiation pressure is excluded. The legacy
`data.rad_pressure` field still contains `elec_pressure + true_radiation_pressure`
for backward compatibility but should be considered deprecated; new code uses
`data.plasma_pressure` and `data.rad_pressure_true` instead.

**Ablation pressure = spatial peak at breakout instant** (not peak-over-time).
`shock_foot_pressure_Gbar` reports `np.max(plasma_pressure[i_breakout, :])`.
This matches the typical ICF reporting convention and reproduces Tom's published
Helios reference plot (107 Mbar at 0.83 ns for the CH sphere baseline).

**Shock breakout detection** is now target-geometry-agnostic. For capsule
targets (multi-region), a 5-zone gas-cavity probe inside the gas/fuel interface
detects the moment the shock crosses into the gas. For solid-shell targets
(single region, e.g. CH sphere/slab), a 5-zone rear-face probe at zones 0..4
detects shock arrival at the inner face. Threshold: `max(100 × P_initial, 1e-3 Gbar)`
i.e. floor of 1 Mbar to filter preheat. Both target types report
`shock_breakout_time_ns`, `shock_breakout_P_gas_Gbar` (rear face for solid shells),
`shock_breakout_P_ice_Gbar` (drive side for solid shells), `shock_foot_pressure_Gbar`
(spatial peak at breakout = ablation pressure).

**Pressure unit display.** Pre-stagnation pressures (shock breakout, ablation,
foot) are reported in **Mbar**. Stagnation-era pressures (hot-spot, burn-averaged)
remain in **Gbar**. Conversion is display-only; stored attributes keep their
`*_Gbar` suffix throughout.

### Single-region target support

The pipeline now handles CH-only slab and sphere targets that were previously
unsupported. Guards added in `analyze_stagnation_phase`, `_compute_adiabat`,
`_compute_adiabat_at_breakout`, `_compute_hs_radius_vs_time`,
`_compute_fuel_rhoR_vs_time`, peak-velocity search (in implosion analysis),
`extract_histories_from_run_data` (returns NaN-stub dict with the keys that
`calculate_burn_averaged_metrics` expects), and `compute_performance_metrics`
(CR section). Skip messages are logged when relevant analysis cannot apply.
Non-applicable summary-output sections (mass fractions, hot-spot block, ignition
block) are guarded by attribute presence so the summary is partial but clean.

This is structurally messy. A future TODO is the **target_class attribute
refactor** — adding `ICFRunData.target_class ∈ {"capsule", "slab", "sphere_test"}`
and routing the pipeline through dedicated analysis paths.

### Calibration baseline: PDD_26b is the best Helios match

`Olson_PDD_26b_burn` is the validated best-match calibration to the published
Olson 2021 LILAC reference at 2.15 MJ. Located at
`~/Sims/Xcimer/Olson_PDD/Olson_PDD_26b_burn/` on Mac Studio.

**Geometry:** half-cone 20°, spot radius 0.16 cm (Gaussian), focus d = +0.20 cm,
wavelength 0.351 μm, power multiplier 1.0.

**Drive:** foot 30 TW (0.10–5.00 ns), ramp 5.00–9.00 ns, peak 329 TW
(9.00–12.70 ns), total duration 12.70 ns, 2.15 MJ delivered.

**Other settings:** flux limiter f = 0.06; non-local alpha transport;
3T burn enabled.

**Target structure (4 regions, 350 zones total):**
- DT Vapor: zones 0–150 (ρ₀ = 6e-4 g/cc)
- DT Solid (ice): zones **151–190** (ρ₀ = 0.222 g/cc) — only 40 zones
- DT-CH foam: zones 191–320 (ρ₀ = 0.242 g/cc)
- CH Skin: zones 321–350 (ρ₀ = 1.049 g/cc)

`region_interfaces_indices` shape is `(n_times, 4)` with boundaries at 151, 191,
321, 351. Adiabat uses `ri[:, 0]:ri[:, 1]` = zones 151–190.

**Headline results vs Olson 2021 LILAC reference:**

| Metric | PDD_26b | Reference | Δ |
|---|---|---|---|
| Peak implosion velocity | 478.4 km/s | 470 km/s | +1.8% ✓ |
| Stagnation time | 13.25 ns | ~13.47 ns | −1.6% ✓ |
| Imploded DT mass | 0.60 mg | ~0.60 mg | ≈ 0 ✓ |
| Stagnation CR | 29.6 | 29.0 ± 3 | +2.1% ✓ |
| Fraction absorbed | 85.6% | 65.1 ± 9.3% (3D) | +31.5% |
| **DT ice adiabat (peak v)** | **1.43** | **3.0 ± 0.5** | **−52%** |
| DT ice base adiabat (at breakout) | 0.63 | — | — |
| Shock breakout time | 9.00 ns | — | — |
| Ablation pressure at breakout | 67.7 Mbar | — | — |
| ρR (cold fuel) | 0.61 g/cm² | 1.10 ± 0.05 g/cm² | −44% |
| Hot-spot pressure (stagnation) | 161.7 Gbar | — | — |
| ⟨T_hs⟩ (burn-averaged) | 19.4 keV | 22.5 ± 2.0 keV | −13.6% |
| ⟨P_hs⟩ (burn-averaged) | 165 Gbar | 193 ± 20 Gbar | −14.5% |
| Yield | 13.8 MJ | 87.4 MJ | −84% |
| Gain | 6.4 | 40.6 | −84% |

Geometric and timing channels match within a few percent. Burn-averaged temperature
and pressure are within ~14% of LILAC. Persistent residuals are in adiabat (−52%),
ρR_cf (−44%), absorbed fraction (+32%), and yield (−84%). These all link to the
1D-vs-3D limitation: Helios maintains full spherical coupling throughout the
implosion while the fixed-pointing NIF PDD beams illuminate a decreasing solid
angle as the capsule shrinks.

A handoff document for cross-code comparison with MULTI-IFE has been generated
at `~/Sims/Xcimer/Olson_PDD/MULTI comparison/PDD_26b_MULTI-IFE_Comparison.docx`.

### Key calibration finding: ablation residual is not laser-coupling

Two geometrically distinct beam configurations both reach the LILAC reference
peak velocity:
- **PDD_26b**: forward focus, cone 20°, d = +0.20 cm
- **DM_01b**: inverted focus, cone 3.5°, d = −2.5 cm, 230 TW, 1.365 MJ

Both leave the same residuals on adiabat, IFAR, and imploded-fuel-mass channels.
This rules out laser-coupling geometry as the cause of the residual. The
remaining suspects are **flux limiter implementation** and **EOS choices** for
the cold fuel and ablator.

### CH sphere flux-limiter scan

Validated test case (CH-only 100-μm shell at R=2 mm, 502.7 TW, 2.5 ns square
pulse, equivalent to 1e15 W/cm² slab intensity, half-cone 1° focus d=−2.0 cm).
Reference values: shock breakout 0.83 ns, ablation pressure ~107 Mbar, 100%
absorbed.

**Validation result:** Postprocessor reproduces 0.84 ns breakout, 106.6 Mbar
ablation pressure (essentially exact match) at f = 0.06.

**Flux-limiter scan (f = 0.06, 0.08, 0.10, 0.12):**

| f | Breakout (ns) | Ablation (Mbar) | P_rear (Mbar) |
|---|---|---|---|
| 0.06 | 0.840 | 106.56 | 10.03 |
| 0.08 | 0.840 | 106.75 | 13.41 |
| 0.10 | 0.840 | 106.53 | 15.58 |
| 0.12 | 0.830 | 108.12 | 1.10 |

**Conclusion (provisional):** flux limiter is a weak lever for ablation pressure
across the standard range. dα/df ≈ 0.23 per unit f (marginally above the
CLAUDE.md 0.2 threshold). MULTI-IFE feedback subsequently indicated Helios's
flux-limiter implementation may not be running properly — motivating a much
tighter f = 0.005 test on WfCDT_01b as the highest-priority next-session item.

### Background TODOs (when convenient)

- `target_class` attribute refactor (replaces ad-hoc single-region guards)
- `_track_ablation_front` for single-region targets (currently fails silently
  with "all zero" history; not blocking since `shock_foot_pressure_Gbar` gives
  the headline number)
- `data_builder.build_run_data()` `time_unit` default — should auto-detect
  rather than passthrough
- Hot-spot pressure jumped from 103 → 162 Gbar after plasma-pressure +
  stagnation-time changes; worth understanding before next big calibration
  round, not a blocker
- PDD_26b Stagnation CR rerun under current conventions (the doc footnotes
  this; geometric definition shouldn't move it but verification is cheap)
- VI_6 (Vulcan HDD) restart — set aside during PDD work, several open items
  on over-ablation physics and CR convention mismatch

### Retired items (no longer pursued)

- Geometry sweep with cone = 35°, d = 0.22–0.23 cm (laser coupling not the
  residual)
- Foot power / adiabat counter-intuitive trend (likely artifact of old
  foam-inclusive adiabat convention)

### Patches landed this session

In chronological order, all on origin/main as of session close:

- Stage 2d shock breakout (gas-cavity / rear-face pressure monitor with
  geometry-agnostic probe selection; density-based shell adiabat;
  base-adiabat-at-breakout method)
- Single-region target guards across pipeline + summary output
- Mbar/Gbar display convention for pre-stag/stag pressures
- Ablation-pressure = spatial peak at breakout instant
- Plasma-pressure convention end-to-end (data_builder + 8 call sites)
- DT ice-only adiabat zone selection (both `_compute_adiabat` and
  `_compute_adiabat_at_breakout`)
- Single-region stub-dict fix for `extract_histories_from_run_data`
- CR computation guard against single-region targets

## Session Update — May 3, 2026

Supersedes earlier next-session priorities, the EOS-variation suggestion,
and all WfCDT items.

### WfCDT retired

WfCDT was a David Montgomery test problem with no published benchmark to
calibrate against — a dead end for code validation. All WfCDT next-session
items (f=0.005 rerun, ice-only baseline, full FL scan) are dropped.

### Calibration strategy: PDD → HDD

The active path is now (1) close the PDD calibration against the Olson
2021 LILAC / xRAGE / HYDRA 1.4× cluster on the PDD_26 series, then
(2) transfer the calibrated knobs to the Vulcan HDD target (VI_6). PDD
is the testbench because the published reference is the strongest
cross-code anchor we have; HDD is the deliverable per Christopherson
11/3/2025 (v=410 km/s, α=6, KE=300 kJ, stag DT=3 mg, P_n-avg=210 Gbar,
ρR_n-avg=2.1 g/cm², T_n-avg=4.8 keV, yield=0.6 MJ).

### EOS is not a knob

PROPACEOS and SESAME are the only EOS options available, and both are
already in use (SESAME 5271 for DT regions, PROPACEOS for DT-CH foam
and CH skin). EOS variation is removed from the priority list. Earlier
notes treating EOS as a calibration lever (e.g. lines ~590, ~1063,
priority 4 of the April-26 appendix) are superseded.

### Calibration knobs going forward

1. **Per-material flux limiter.** Each RHW region has its own FL block;
   the laser-coupling-relevant ones are CH Skin (outermost) and DT-CH
   foam (inner ablator). Tightening FL on those two while leaving DT
   regions at f=0.06 is the standard move.
2. **Time-dependent power scaling on the RHW pulse.** Mimics the 3D
   ray-trace miss as the capsule shrinks below the NIF PDD beam-pointing
   solid angle. Empirical target is a time-resolved absorbed-energy
   curve from Olson 2021 (HYDRA preferred); first-cut piecewise
   multiplier on `LaserPwrDeliveredForBeam` past ~9 ns is acceptable
   until the time-resolved curve is extracted.
3. **Helios source fixes (escalation only).** If FL + power scaling
   can't close the residual on PDD, identify and propose fixes to
   Helios's coupling/transport. Not pursued speculatively.

### PDD_26 flux-limiter scan results

| f      | Stag (ns) | v_imp (km/s) | α_ice | Abs % | ρ_max (g/cc) | P_hs (Gbar) | Yield (MJ) | Imploded DT (mg) |
|--------|-----------|--------------|-------|-------|--------------|-------------|------------|------------------|
| 0.06   | 13.43     | 479          | 1.09  | 87.2  | 359          | 390         | 59.5       | 1.65             |
| 0.005  | 13.83     | 438          | 1.37  | 78.1  | 300          | 214         | 49.3       | 0.60             |
| 0.001  | 14.25     | 417          | 1.58  | 70.4  | 197          | 154         |  9.8       | 0.60             |

Olson 1.4× cluster reference: v=470, α=3.0, abs=65±9%, P_hs=193±20 Gbar,
yield=87 MJ.

**Findings:**
- Helios FL implementation IS responsive on the PDD capsule, contra
  the CH-sphere flat result (which was the degenerate case — solid
  shell + square pulse has no shock-train or compressional adiabat
  dynamics). The MULTI-IFE-flagged "FL not running" worry is retired
  for this target class.
- All channels move monotonically and in the textbook direction with
  decreasing f: lower absorption, lower drive, higher α, lower yield.
- Adiabat slope is small: dα/d(log₁₀ f) ≈ 0.4. Closing the 1.09→3.0
  gap by FL alone would require f ~ 1e-7, well past extinction.
- f = 0.001 extinguishes ignition (HS ρR never reaches 0.3 g/cm²).
  f = 0.005 sits at the edge — ignition fires (13.87 ns) but no
  complete-propagation timestamp.
- Imploded DT mass floors at 0.60 mg between f=0.005 and f=0.001 →
  ablation-channel residual is NOT in the laser-coupling channel.
  Confirms the "two paths to velocity" April finding.
- Conclusion: FL is a real lever but insufficient alone. Combined
  FL + late-time laser-power scaling is the natural next try.

### Per-material FL parser (FIXED May 3, 2026)

`RHWParser._parse_flux_limiter` previously kept only the *first*
region's value (DT Vapor, always 0.06). When foam/CH-skin FL was
changed, summary output stayed at 0.06 — making it invisible that a
knob was being moved. The warning that fired was easy to miss.

Fix: parser now returns a per-region list,
`[{'region', 'enabled', 'value'}, ...]`, exposed on
`RHWConfiguration.flux_limiters` and `ICFRunData.flux_limiters`. The
scalar `flux_limiter` field now holds the *outermost* region's value
(CH skin / ablator — laser-coupling-relevant), not the first.
`icf_output.py` prints all regions in the LASER CONFIGURATION block
when per-region data is available.

### Published-JSON pointer fix

`Olson_PDD_26af001_burn_published.json` and
`Olson_PDD_26af005_burn_published.json` currently reference the 1.0×
drive Olson values (v=370, α=7.4, yield=62, abs=99%) — both runs are
1.4× drive variants and should compare against the cluster reference
(v=470, α=3.0, yield=87, abs=65±9%). Action: drop in
`Olson_PDD_26c_burn_published.json` (created May 3, 2026 with per-code
LILAC/xRAGE/HYDRA values) under both filenames.

### Next-session priorities (ordered)

1. **Combined FL + late-time power scaling on PDD_26.** Foam + CH-skin
   FL = 0.005, plus piecewise multiplier on the RHW power table
   (initial guess: ~0.75× past ~9 ns; calibrate against absorbed-energy
   target ≈ 65%). Target metrics: absorbed → 65%, α → 3, yield → 30+ MJ.
2. **Backfill 26af001 / 26af005 published JSONs** with the 26c reference.
3. **Restore notebook RANSAC shock-tracking algorithm.** Iterative
   `sklearn.linear_model.RANSACRegressor` + `LinearRegression`, fits
   straight lines to shock-front points in (t, r) space, removes
   inliers, computes pairwise intersections for coalescence points.
   Critical for the Christopherson "simultaneous shock breakout"
   tuning target on Vulcan HDD.
4. **Extend comparison framework with per-code Δ columns** (LILAC,
   xRAGE, HYDRA), per Tom 5/3: "different algorithms and databases in
   the three codes." Use the underscore-prefixed per-code keys in the
   new 26c JSON as the data source.
5. **HDD transfer.** Once PDD calibration is closed, port the
   FL + power-scaling settings to VI_6 and validate against
   Christopherson 11/3/2025 deliverables.

### Background TODOs (unchanged from April 26 appendix)

- `target_class` attribute refactor
- `_track_ablation_front` for single-region targets
- `data_builder.build_run_data()` `time_unit` auto-detect
- Hot-spot pressure 103→162 Gbar root-cause after plasma-pressure +
  stagnation-time changes
- PDD_26b stagnation CR rerun under current conventions

### Retired (no longer pursued)

- WfCDT_01b f=0.005 rerun
- Full WfCDT FL scan
- WfCDT clean ice-only baseline
- EOS variation as a calibration lever (no remaining options)
- Geometry sweep cone=35°, d=0.22-0.23 cm (April 2026; laser-coupling
  not the residual)

### Session Update — May 4, 2026

Continues PDD calibration work from May 3. Per-region FL parser fix
landed and validated; FL/geometry/foot calibration scan executed; best
1D calibration point identified and burn run submitted.

### Per-region FL parser landed (committed May 3 evening)

`RHWParser._parse_flux_limiter` rewritten to return per-region list
`[{'region', 'enabled', 'value'}, ...]`. Surfaced on
`RHWConfiguration.flux_limiters` and mirrored on
`ICFRunData.flux_limiters`. Scalar `flux_limiter` field now reports
the *outermost* region (CH skin / ablator — laser-coupling-relevant).
`icf_output.py` prints all regions in the LASER CONFIGURATION block
when per-region data is present. Validated end-to-end on Studio
against PDD_26a/26af005/26af001: all three RHW files parse correctly,
showing CH Skin and DT-CH foam at the run-specific f, DT Solid and DT
Vapor at canonical 0.06.

### Published JSON expanded with per-code reference values

`Olson_PDD_26c_burn_published.json` extended from 38 → 79 lines.
Per-code blocks (LILAC, xRAGE, HYDRA) extracted from Olson 2021
Figs 6, 7, 8 covering stagnation time, ignition time, bang time,
peak total ρR, peak HS ρR, peak compressed-shell density at ignition,
T_ion HS at ignition, HS radius at ignition. Per-code keys carry
leading underscore so the comparison framework currently ignores them
(activate by stripping underscore once Tom decides to add per-code Δ
columns; per his May 3 decision option 1 — show three Δ columns,
one per reference code).

Three codes are *close* in shape but differ meaningfully:
- ~1 ns spread in stagnation time (HYDRA earliest at 12.5 ns,
  LILAC at 13.5 ns, xRAGE at 13.7 ns)
- ~30% spread in peak compressed-shell density at ignition
  (xRAGE 170, LILAC 220, HYDRA 250 g/cm³)
- HS radius at ignition pinned at 120 μm in all three (annotated)

Per-code differences should inform "Helios within the cluster spread"
framing rather than "Helios matches one number."

### Published JSON pointer fix on 26af001 / 26af005

`Olson_PDD_26af001_burn_published.json` and `_26af005_published.json`
were pointing at the 1.0× drive Olson reference (v=370, α=7.4, abs=99,
yield=62). Both runs are 1.4× drive variants and should compare against
the cluster reference (v=470, α=3.0, abs=65±9.3, yield=87.4). Fixed
May 3 by copying the corrected 26c JSON over both. Comparison columns
now read sensibly: f=0.005 → v=438 vs 470 (−6.7%) instead of the
misleading +18.5% against the wrong reference.

### PDD calibration scan: FL × geometry × foot at the Olson pulse

Six no-burn runs at FL=0.02 on foam+CH skin (DT regions held at 0.06):

| Run                         | Cone  | Spot  | Foot | v_imp | α    | Abs % | Bang  |
|-----------------------------|-------|-------|------|-------|------|-------|-------|
| PDD_28_fab005 (FL=0.005)    | 20°   | 0.16  | 30   | 375   | 1.87 | 62.3  | 15.00 |
| PDD_35_fab005 (FL=0.005)    | 35°   | 0.25  | 25   | 375   | 1.86 | 62.3  | 15.00 |
| PDD_35_fab02_nb             | 35°   | 0.25  | 25   | 397   | 2.02 | 67.7  | 14.50 |
| PDD_35_fab02_foot40_nb      | 35°   | 0.25  | 40   | 397   | 2.40 | 66.2  | 14.13 |
| PDD_20_fab02_foot40_nb      | 20°   | 0.20  | 40   | 409   | 3.25 | 75.5  | 13.47 |
| **PDD_20_fab02_foot25_nb**  | **20°** | **0.20** | **25** | **437** | **1.78** | **77.6** | **13.80** |

Cluster envelope reference: v=470, α=3.0±0.3, abs=65±9.3,
bang=13.47 ns.

### Four physics findings on the PDD target

**1. Per-region FL parsing is real and required.** Yesterday's
inferred conclusions about "FL response" silently treated the foam+CH
skin pair as the lever, but parsed FL still showed 0.06 due to the
parser bug. With the fix, FL response on coupling materials only is
now an explicit, traceable knob.

**2. FL=0.005 saturates the heat-flow channel; geometry becomes
irrelevant.** PDD_28_fab005 (cone=20°, spot=0.16) and PDD_35_fab005
(cone=35°, spot=0.25) gave *identical* downstream metrics — same v,
same α, same abs, same bang time, same yield, same imploded mass —
despite very different ray-trace deposition profiles (peak coronal
intensity, r_crit, intensity at r_crit all differed). The flux limiter
was the rate-limiting bottleneck; geometry within the explored range
couldn't propagate to downstream observables.

**3. FL=0.02 restores geometry sensitivity.** Same comparison at
FL=0.02 (cone=20° vs cone=35°, foot=40) gave dramatically different
results: v 397 → 409 (+3%), α 2.40 → 3.25 (+35%), abs 66.2 → 75.5
(+14%), bang 14.13 → 13.47 (essentially exact match to cluster), peak
coronal I 1.32e15 → 1.77e15 W/cm², r_crit 1417 → 841 μm, I at r_crit
4.43e14 → 1.51e15 W/cm² (3.4×). At looser FL the system has heat-flow
headroom; geometry-driven differences in deposition concentration now
propagate to the implosion.

**4. Foot/peak orthogonality is geometry-dependent, not universal.**
At cone=35°/spot=0.25, foot 25 → 40 TW raised α by 0.38 with v_imp
unchanged (foot is the α knob, peak is the v knob, independent). At
cone=20°/spot=0.20, foot 25 → 40 TW raised α by 1.47 (much steeper)
*and* dropped v_imp by 28 km/s (foot couples to v). At the tighter
geometry the foot's deeper deposition (smaller r_crit, higher I at
critical) drives velocity-relevant ablation, breaking the
orthogonality.

### Calibration ceiling: α ≈ 2 at the Olson pulse, FL=0.02

At foot=25 (Olson reference), FL=0.02, no geometry combination
explored produced α ≥ 2.7 (lower edge of cluster envelope 3.0±0.3).
Confirms the April 2026 finding ("adiabat cannot be raised from 2.0
to 3.0 by FL alone without crossing into catastrophic over-drive").
The α=3.25 point at PDD_20_foot40 required foot=40 TW, 60% above
Olson's nominal foot — off-spec from the reference pulse and not a
defensible calibration.

### Selected calibration point: PDD_20_fab02_foot25_burn

Final geometry and FL settings for PDD calibration:

| Parameter      | Value                            |
|----------------|----------------------------------|
| Wavelength     | 0.351 µm                         |
| Cone (γ)       | 20°                              |
| Spot radius    | 0.20 cm (Gaussian)               |
| Focus position | 0.22 cm                          |
| Foot power     | 25 TW (0.10 – 5.00 ns)           |
| Peak power     | 329 TW (9.00 – 12.70 ns)         |
| Pulse duration | 12.70 ns (Olson 2021 1.4× drive) |
| FL CH Skin     | 0.020                            |
| FL DT-CH foam  | 0.020                            |
| FL DT Solid    | 0.060                            |
| FL DT Vapor    | 0.060                            |
| Alpha transp.  | Non-local                        |

No-burn results: v=437 (−7% vs 470), α=1.78 (−41% vs 3.0±0.3),
abs=77.6% (just outside cluster upper edge 74.4%), bang=13.80 ns
(+2.4% vs 13.47), CR=37 (no-burn — will drop with α-pressure),
imploded DT 0.61 mg.

Burn run submitted; preliminary observation ⟨T_hs⟩ reaching 20+ keV
during burn, consistent with cluster ⟨T_hs⟩=22.5 keV. Final results
pending.

### Persistent residuals (1D Helios limits at this calibration)

- **v_imp short by 7%** (437 vs 470). Smallest velocity gap achieved
  at any FL in any geometry combination explored today, and at the
  Olson reference pulse. Attributable to 1D ablation-channel
  efficiency vs 3D codes; consistent with April "two paths to
  velocity" finding (Helios reaches reference velocity only at
  geometries that miss ablation physics).
- **α short by 41%** (1.78 vs 3.0±0.3). The α=2 ceiling at the Olson
  pulse stands. Closing this gap would require off-spec foot power
  (foot=40+) or non-Olson pulse modifications.
- **abs slightly out of envelope** (77.6 vs 74.4 upper edge, +1.1pp).
  Absorbed energy 1.61 MJ vs cluster 1.40 MJ (+15%). Real but minor.

### Naming convention adopted for the PDD calibration scan
Olson_PDD_<base><variant><run-type>
base       = integer or integer+letter, geometry baseline
(e.g. 28 = first schematic-geometry test,
35 = full schematic geometry,
20 = cone=20°/spot=0.20 hybrid,
26b = April 2026 best-velocity match)
variant    = compact descriptor of knob settings beyond base
fab005, fab02 = foam+CH skin FL = 0.005, 0.020
foot40        = foot pulse 40 TW
p110, p115    = power multiplier 1.10, 1.15
(omit if vanilla — uniform f=0.06, foot=25,
power×1.0)
run-type   = burn = full simulation with fusion enabled
nb   = no-burn / hydro-only

Example: `Olson_PDD_20_fab02_foot25_burn` = cone=20° geometry,
foam+skin FL=0.02, foot 25 TW, burn enabled.

Base number changes only when geometry changes. Knob variations live
in variant suffix, not as new base numbers. Avoids the 26a/26af005
ambiguity where parallel knob settings looked like sequential
geometry. Variant suffix sorts naturally in `ls`.

### Output schedule for new no-burn runs (PDD)

454-timestep staged schedule, EXODUS time-point list in RHW:

| Phase                  | Range          | Δt       |
|------------------------|----------------|----------|
| Pre-drive              | 0 – 4 ns       | 0.25 ns  |
| Foot + 1st shock       | 4 – 10 ns      | 0.05 ns  |
| Main drive ramp        | 10 – 13 ns     | 0.05 ns  |
| Peak velocity / late   | 13 – 13.4 ns   | 0.01 ns  |
| Stagnation + burn      | 13.4 – 14.4 ns | 0.005 ns |
| Post-burn coast        | 14.4 – 16 ns   | 0.1 ns   |

Resolves bang time to 0.005 ns, peak velocity to 0.01 ns, and shock
breakouts to 0.05 ns. Used for all PDD_28/35/20 runs today; results
reproducible.

### Next-session priorities (ordered)

1. **Verify burn run lands cleanly** at PDD_20_fab02_foot25_burn.
   Decision criteria: ⟨T_hs⟩ in 20–25 keV range, yield > 5 MJ,
   bang time within ±2% of no-burn, CR drops to 28–32 envelope under
   α-pressure. If all four → calibration is closed. If yield far
   below cluster despite T_hs match → spot-size scan at cone=20°
   becomes priority (PDD_20 with s=0.16, 0.18, 0.22 to bracket
   between PDD_26b s=0.16 and today's s=0.20).
2. **HDD transfer**, contingent on (1). Apply Christopherson
   11/3/2025 Euler scaling to take PDD_20 calibration → VI_6 HDD
   geometry. Targets: v=410 km/s, α=6, KE=300 kJ, stag DT=3 mg,
   P_n-avg=210 Gbar, ρR_n-avg=2.1 g/cm², T_n-avg=4.8 keV,
   yield=0.6 MJ. Note: HDD α=6 is achievable in 1D since the higher
   foot pulse for HDD targets directly drives the higher adiabat —
   the α=2 ceiling in PDD calibration was Olson-pulse-specific.
3. **Per-code Δ columns** in the comparison framework (Tom's option 1
   from May 3). Activate underscore-prefixed per-code keys in the
   26c/26af001/26af005 JSONs; extend `compare_with_published()` to
   render three additional Δ columns vs LILAC, xRAGE, HYDRA. Lower
   priority than (1) and (2) but useful diagnostic for the
   "Helios within cluster spread" framing.
4. **Restore notebook RANSAC shock-tracker** (carried from May 3).
   `sklearn.linear_model.RANSACRegressor` + `LinearRegression`,
   iterative line-fitting with inlier removal, pairwise intersections
   = coalescence points. Critical for Christopherson's "simultaneous
   shock breakout" tuning target on Vulcan HDD.
5. **Background TODOs** from April 26 / May 3 appendices unchanged:
   `target_class` attribute refactor, `_track_ablation_front` for
   single-region targets, `data_builder.build_run_data()` `time_unit`
   auto-detect.

### Retired (no longer pursued)

Carried forward from May 3 unchanged: WfCDT_01b f=0.005 rerun, full
WfCDT FL scan, WfCDT clean ice-only baseline, EOS variation as
calibration lever. Plus retired May 4:

- **Time-dependent power scaling** as a calibration knob. The April
  motivation (mock 3D beam-miss as R(t) shrinks) was sound, but
  FL=0.02 alone brings absorbed energy close enough to the cluster
  envelope (1.61 vs 1.40 MJ; abs 77.6 vs 74.4 upper edge) that the
  added complexity of a time-dependent power table is not justified
  for the current calibration. Keep as a backup if HDD transfer
  reveals a gap that FL alone can't close.
- **Geometry exploration past cone=20° / spot=0.20** at the Olson
  pulse with FL=0.02. Today's scan covered the relevant range; the
  v/α/abs/bang trades are well-characterized and adding intermediate
  cone angles would refine gradients without changing conclusions.

  ### Spot-size endpoint at s=0.16 (May 5)

Burn run at PDD_20 calibration with spot=0.16 (the April PDD_26b
spot value) instead of 0.20: `Olson_PDD_20_fab02_foot25_s016_burn`.
v=462.8 km/s (−1.5% vs 470), bang=13.645 ns (DEAD on LILAC's 13.5),
HS radius at ignition 125 μm (vs cluster 120, +4%), yield 59.1 MJ
(−32% vs cluster 87.4 — best yield match achieved), imploded DT
1.77 mg (3× the s=0.20 value of 0.60 mg). However: ⟨T_hs⟩ 38.1 keV
(+69%), ⟨P_hs⟩ 429 Gbar (+122%), abs 83.9% (+29%, well out of
envelope), α=1.05 (−65%, very low).

The s=0.16 vs s=0.20 comparison brackets cluster physics:
- s=0.20 undercouples compression: T_hs/P_hs/CR/density in envelope,
  but yield 70% short and ρR_cf 39% short.
- s=0.16 overcouples burn: v/bang/HS radius at ignition match cluster,
  yield only 32% short, but T_hs/P_hs overshoot and abs out of
  envelope.
- Cluster physics sits between these two regimes. Helios cannot
  produce a single point that lands all eight metrics simultaneously
  at the Olson reference pulse with FL=0.02 — confirming that the
  residual is a 1D vs 3D physics difference (specifically: 1D
  Helios's ablation channel cannot reproduce the cluster's combined
  efficient mass-coupling AND moderate absorption).

Both runs valid as calibration anchors. For HDD transfer, both
reference points carried forward; spot-size choice in HDD will be
informed by the HDD coupling intent rather than direct PDD transfer.

# Session updates — 2026-05-08 (HDD halfraum design + postprocessor halfraum support)

This block can be appended to the existing CLAUDE.md (or merged into the
relevant existing sections by hand). Items are grouped by topic so they
slot into the document structure already in use.

---

## HDD halfraum target architecture (1D Helios approximation)

5-region target representing the Vulcan HDD halfraum:

| Region | Material | Initial r (cm) | Thickness | Initial ρ | Role |
|---|---|---|---|---|---|
| 1 | Gas fill (DT vapor) | 0 – 0.1875 | — | 3×10⁻⁴ g/cc | Hot spot |
| 2 | Wetted foam (DT-CH) | 0.1875 – 0.2145 | 270 µm | 0.30 g/cc | Cold fuel |
| 3 | CD shell | 0.2145 – 0.2185 | 40 µm | 1.09 g/cc | Ablator |
| 4 | Pseudo void (He) | 0.2185 – 0.9185 | 7 mm | 3×10⁻⁴ g/cc | Hohlraum interior |
| 5 | Cu shell | 0.9185 – 0.919 | 5 µm | 8.93 g/cc | Radiation converter |

3-beam laser splits the published P1–P7 pulse table:

| Beam | Cone | Spot | Target | Pulse segments |
|---|---|---|---|---|
| 1 | 1.0° | 0.22 cm | Cu shell (r=0.92 cm) | P1 (220 TW, 0–1 ns), P2 (65 TW, 4.5–7.75 ns), P3 (295 TW, 7.75–10.3 ns), P4 (614 TW, 10.31–11.7 ns) |
| 2 | 0.7° | 0.16 cm | CD shell (r=0.22 cm) | P5 (604 TW, 11.71–13.0 ns), P6 (564 TW, 13.0–14.41 ns) |
| 3 | 0.5° | 0.11 cm | CD shell (r=0.22 cm) | P7 (424 TW, 14.41–15.6 ns) |

Total delivered energy 4.0 MJ. Prescribed `[Rad Source Data]` is OFF in
this configuration; rad drive is generated self-consistently from the
beam-1/Cu interaction rather than being prescribed at Rmax. Beams 2 and
3 must penetrate the expanded Cu corona to reach the CD shell — that
geometric coupling is the open question.

Reference deck: `~/Sims/Xcimer/HDD_26/HDD26_Cu_He_FL02_Zoom_1_nb/`

---

## Vulcan published pulse (P1–P7), read from Thomas figure

| Segment | Time (ns) | Power (TW) | Energy (kJ) |
|---|---|---|---|
| P1 | 0 – 2 | 220 | 440 |
| (off) | 2 – 4 | 0 | 0 |
| P2 | 4 – 8 | 70 | 280 |
| P3 | 8 – 10 | 290 | 580 |
| P4 | 10 – 11 | 610 | 610 |
| P5 | 11 – 12 | 600 | 600 |
| P6 | 12 – 13 | 560 | 560 |
| P7 | 13 – 15 | 420 | 840 |

Total ≈ 3.9 MJ, matches published Thomas burn-off delivered energy
within reading error. Values uncertain to ±10% on each segment until
cross-checked against the underlying paper data.

---

## First halfraum run result (HDD26_Cu_He_FL02_Zoom_1_nb)

Sim ran 20 ns; capsule barely imploded:
- CR = 1.19 (effectively no compression)
- Peak velocity = 370 km/s at t = 19.1 ns (capsule still accelerating at end of sim)
- No bang time within 20 ns window
- ⟨T_hs⟩ = 0.03 keV, yield = 3.5×10⁶ reactions
- Frac absorbed = 99% — but to Cu corona, not capsule

Three competing culprits requiring diagnosis:
1. **5 µm Cu too thin** for sufficient x-ray emission at the 100 eV target T_rad
2. **7 mm He gap absorbing more than expected** (low-density He should be
   transparent to soft x-rays; not transparent if rad emission is too soft)
3. **Sim too short** relative to drive geometry (capsule may bang at 25–30 ns)

Resolution path: look at T_rad(t) at zone 425 (Cu inner surface) vs zone
376 (CD outer surface). Strong T_rad at Cu inner but weak at CD outer →
He gap eats radiation. Weak T_rad at Cu inner → Cu thickness issue.
Both look healthy → just need longer sim.

Sensitivity sweep suggestions: Cu thickness 5→10→15 µm; He gap 7→3→1 mm;
max sim time 20→30 ns. Run one variation per culprit to triage.

---

## Persistent calibration findings (HDD26)

**Adiabat is unmoved by radiation-drive details.** Three-run sweep of
HDD26 with progressively more rad-drive energy:

| Run | Rad config | Adiabat | Yield (MJ) | Imploded DT (mg) |
|---|---|---|---|---|
| CA1 | cutoff at 2.17 ns | 2.12 | 0.037 | 2.55 |
| CA1x | extended to 4.5 ns, mult=0.8 | 1.99 | 0.018 | 3.13 |
| CA1x1 | extended to 4.5 ns, mult=1.0 | 1.93 | 0.021 | 3.18 |

Adiabat is flat (2.12 → 1.93). Closing the rad-laser gap improves shell
preconditioning — imploded DT moves from 2.55 mg (15% under published)
to 3.13 mg (4% over published), unablated fuel rises from 0.61 to 0.75 —
but *degrades* laser coupling because the corona is more diffuse. Net
absorbed energy drops 2.22 → 1.76 MJ (-21%). Yield drops because
T_hs⁴ scaling dominates over the modest ρR_cf gain.

Burn-on version of CA1x: α onset at 14.3 ns (1 ns before bang), yield
0.018 → 0.023 MJ (+28%), ⟨T_hs⟩ 3.07 → 3.35 keV. No propagating burn
(ρR_hs = 0.059 g/cm² vs threshold 0.3). The HDD26 calibration gap is
"several factors short on everything" — not "near-miss on ignition."

**Confirmed null absorption levers in 1D Helios:**
- Reflection at r_crit (rp0/rp1/rp2/rp3 sweep on Olson_PDD_20: all four
  runs identical, frac absorbed stays 83.9% across 0–30% reflection)
- Flux limiter f=0.02 vs f=0.06 (PDD)
- Flux multiplier 0.8 vs 1.0 (radiation drive)
- Rad-drive timing extension

In 1D Helios, reflected rays go back through the same corona and are
re-absorbed there — energy cannot escape the system. The only knob that
actually reduces absorbed energy is **Power multiplier** (which just
delivers less laser).

---

## Postprocessor halfraum support (May 2026)

`ICFRunData` now auto-detects target geometry from `region_names`:
- `target_class = "capsule"` (default — current behaviour)
- `target_class = "halfraum_capsule"` (when external structure detected)

Detection keywords (case-insensitive substring, scanned from outside inward):
`pseudo void`, `void`, `cu shell`, `pb shell`, `au shell`, `u shell`,
`ta shell`, `hohlraum`, `halfraum`, `he fill`, `helium`, `hohlraum gas`.

Two new helper properties on `ICFRunData`:
- `capsule_outer_idx` — column in `region_interfaces_indices` for the
  capsule outer surface. For halfraum_capsule, points to ablator outer
  (= inner edge of first external region) rather than grid edge.
- `fuel_ablator_idx` — column for the fuel/ablator interface, one column
  inside `capsule_outer_idx`.
- Plus a method `capsule_outer_node(t)` returning the node index at
  timestep `t` for capsule-bounded zone-index searches.

For 3- and 4-region standard capsule targets these reduce to the
historical `ri[:, -1]` and `ri[:, -2]` — zero regression on existing
runs. For 5-region halfraum targets the analyzer indexes into the
capsule outer surface and fuel/ablator interface within the capsule
rather than the grid edge.

Analyzer call sites updated to use the new properties:
- `_compute_adiabat`, `_compute_adiabat_at_breakout` (fuel-velocity range,
  no-ablator branches)
- `_compute_ifar` (fuel-velocity range, density-shell search)
- `_track_ablation_front` (search bounded at capsule outer node — keeps
  the He/Cu density step from registering as ablation front)
- `_compute_areal_densities` (Total ρR restricted to capsule outer for
  halfraum)
- Neutron-averaged fuel ρR
- Initial fuel/ablator mass split

Two halfraum-specific corrections:
- **Peak-density search bounded to capsule** in `_find_stagnation_time`
  — without this, the initial Cu shell compression at simulation start
  registers as the imploded-fuel peak (Cu at 8.93 g/cc compressing
  transiently to 23 g/cc at t=0.1 ns).
- **`I at r_crit at peak laser power` underflow guard** in
  `analyze_laser_intensity` — for halfraum, peak total power and r_crit
  can refer to different absorbing regions; if the lookup yields
  <1 MW/cm² the diagnostic is set to NaN and reported as `n/a (halfraum:
  peak total power not at capsule r_crit)`.

Also for the first-shock analyzer (`_analyze_first_shock`), the
`ablator_outer_zone` and `fuel_inner_zone` variables now use
`getattr(self.data, 'capsule_outer_idx', -1)` and `getattr(self.data,
'fuel_ablator_idx', -2)` respectively, so the first-shock tracker
bounds at the capsule outer rather than the grid edge.

---

## Multi-beam pulse visibility (May 2026)

`data_builder._collapse_beam_axis` sums `laser_power_delivered` and
`laser_power_on_target` across beams while preserving per-beam tables on
`*_per_beam` attributes. The earlier silent `np.squeeze` bug (no-op on
multi-beam shape, hidden bug for any HDD-style focal-zoom deck) is
fixed. Canary log line: `summed across N beams → (n_t,), per-beam
preserved → (n_t, N)`.

`icf_output` adds a `PULSE SHAPE (per beam)` block for multi-beam runs
showing per-beam on/off time, peak power, and integrated energy. Total
energy across beams sanity-checks against `laser_energy_delivered`. The
existing `(beam 1)` hardcoded section is now correctly labeled as
beam-1-specific only when multiple beams are present. Single-beam runs
unchanged.

`icf_plotting._plot_laser_power` overlays per-beam contributions on the
laser-power PDF page (tab10 colors, summed total in thick black). Makes
the focal-zoom hand-off between beams visible at a glance — this would
have caught the "is the post-peak missing?" confusion immediately.

---

## Convention update — region indexing

For halfraum_capsule (5-region), the convention table is:

| ri column | Boundary | Role |
|---|---|---|
| `ri[:, 0]` | Gas/foam | Hot-spot boundary |
| `ri[:, 1]` | Foam/CD | Fuel/ablator interface (= `fuel_ablator_idx`) |
| `ri[:, 2]` | CD/He | Capsule outer surface (= `capsule_outer_idx`) |
| `ri[:, 3]` | He/Cu | Hohlraum wall inner |
| `ri[:, 4]` | Cu outer | Grid edge |

Capsule analyses use columns 0–2 only. External structure (columns 3, 4)
is excluded from peak velocity, IFAR, adiabat, ablation front,
areal density, mass fractions, and peak density searches.

Old convention for standard 3- and 4-region capsules unchanged.

---

## Open priority list for next session

1. **Halfraum sensitivity sweep** — Cu thickness 5→10→15 µm; He gap
   7→3→1 mm; sim duration 20→30 ns. T_rad(t) at Cu inner (zone 425)
   and CD outer (zone 376) is the diagnostic that triages between the
   three culprits.
2. **Once halfraum produces a viable T_rad(t)**, drop the trace into
   the HDD26 production deck as the new self-consistent rad source
   (replaces extended-rad-drive placeholder).
3. **Power_mult sweep** on best HDD config (1.0 / 1.5 / 2.0) — tests
   drive-limited vs partition-limited hypotheses. If yield scales
   linearly with absorbed energy, drive is the bottleneck. If yield
   saturates with α flat, partition is the bottleneck.
4. **Apply the three pending patches** (icf_analysis_halfraum,
   icf_output_per_beam, icf_plotting_per_beam) and verify regression
   on the standard 3-region (HDD26_CA1x) and 4-region (PDD_20) cases.

---

## Files added or modified this session

| File | Change |
|---|---|
| `helios_postprocess/data_builder.py` | Multi-beam summing + per-beam preservation; `target_class` auto-detection; `capsule_outer_idx` / `fuel_ablator_idx` / `capsule_outer_node()` helpers |
| `helios_postprocess/rhw_parser.py` | Radiation drive temperature parser; `drive_location`, `drive_flux_multiplier` fields |
| `helios_postprocess/energetics.py` | `compute_radiation_drive_energy()` |
| `helios_postprocess/icf_analysis.py` | First-shock analyzer integration; halfraum-aware boundary indexing; bounded peak-density search; LPI underflow guard |
| `helios_postprocess/icf_output.py` | Radiation drive + total drive energy lines; first-shock summary block; per-beam pulse summary block |
| `helios_postprocess/icf_plotting.py` | First-shock 3-panel PDF page; per-beam laser power overlay |
| Project knowledge (this document) | Halfraum target architecture, P1–P7 pulse table, halfraum coupling diagnosis path, postprocessor halfraum convention table |

## Boundary-tally and global cumulative quantities

These are scalar-per-timestep `(n_times,)` quantities that measure energy
**leaving the grid** or **summed over the whole grid**. They're the
right-hand-side terms in the global energy-balance closure equation
and are preferred over inferred fractions when computing residuals.

### Currently loaded (wired into `data_builder`)

| Pipeline attribute | EXODUS variable | Shape | Units | Description |
|---|---|---|---|---|
| `radiation_energy_at_boundary_cum` | `TimeIntRadiationLossAtBds` | (n_t,) | J cumulative | Total radiation energy that has crossed the grid boundaries (Rmin + Rmax combined) up to time t. Used directly in the energy-ledger closure as the rad-escape channel. |
| `particle_energy_escaped_cum` | `particle_time_int_energy_escaped` | (n_t,) | J cumulative | Total particle (predominantly DT neutrons) kinetic energy that has escaped the grid up to time t. Used as a direct measurement of neutron-loss energy, replacing the nominal 14.1/17.6 fraction of fusion yield. |

### Available but not yet wired — wishlist for further accounting

| EXODUS variable | Shape | Likely units | Use case |
|---|---|---|---|
| `FreqIntgRadEnLossRmax` | (n_t,) | J/s? | Frequency-integrated radiation power loss at outer boundary only. Lets us split the combined `TimeIntRadiationLossAtBds` into outer-vs-inner contributions. |
| `FreqIntgRadEnLossRmin` | (n_t,) | J/s? | Same, inner boundary. For most direct-drive runs this should be ~0 (no inner sink); useful as a sanity check. |
| `EnTotRadiation` | (n_t,) | J | Total radiation energy resident in the grid. Cross-check on our `U_rad = 3·P_rad·V` ideal-blackbody calculation. |
| `EnRadSinkTimeIntg` | (n_t,) | J cumulative | Cumulative radiation sink energy in zones. Useful for tracking absorbed-rad / re-emitted-rad balance. |
| `EnExchEleToRadTimeIntg` | (n_t,) | J cumulative | Cumulative electron-to-radiation energy exchange. Indicates how much plasma thermal energy converts to in-grid radiation (vs. radiation that just escapes). |
| `EnJouleHeatingTimeIntg` | (n_t,) | J cumulative | Cumulative Joule heating. Source term if magnetic effects are non-negligible (usually small in ICF). |
| `EnMagTimeIntg` | (n_t,) | J cumulative | Cumulative magnetic-field energy. Same — usually negligible but should close if non-zero. |
| `HeatFluxAtRegionBdEle` | (n_t, 4) | J/(s·cm²)? | Electron heat flux at the inter-region boundaries. Lets the energy ledger be done **per region** rather than globally — would isolate whether energy is stuck in corona, ablator, or fuel. |
| `HeatFluxAtRegionBdIon` | (n_t, 4) | J/(s·cm²)? | Ion heat flux at inter-region boundaries (same use case). |
| `RadNetCoolingRateRegion` | (n_t, 4) | J/(s·cm²) | Net radiation cooling rate per region. Useful for understanding where the trapped coronal thermal energy is converted to radiation. |
| `TimeIntFusionProd_He4_0352` | (n_t,) | reactions | Cumulative DT alpha count (3.52 MeV alphas). Direct cross-check on our DT-neutron-derived alpha-energy estimate. |
| `TimeIntFusionProd_He4_0367` | (n_t,) | reactions | Cumulative DD alpha count (3.67 MeV branch). |
| `TimeIntFusionProd_*_zone` | (n_t, n_z) | reactions | Per-zone fusion-product counts — enable spatial fusion-burn diagnostics. |

### Closure equation with current wiring

After patches 1–3:
## Session 2026-05-12 — Energy-balance diagnostic + lrm4 HDD milestone

### New tooling

- **examples/energy_balance_diagnostic.py** (standalone, 387 lines, uses HeliosRun
  + build_run_data, NOT ICFRunBuilder). Side-by-side ledger for two runs.
  Produces console snapshot tables at peak-v / stagnation / end-of-run plus
  3-page PDF. Tracks: KE inward, KE outward, plasma thermal (ideal gas
  (3/2)(P_i+P_e)V), in-grid radiation (3 P_rad_true V), absorbed laser,
  fusion released (split into alpha-deposited and neutron-escaped), direct
  rad escape from EXODUS boundary tally. Closure: absorbed + alpha_deposited
  = Σ in-plasma channels + rad_boundary + residual_gap. Residual <1% for
  all tested runs.

- **data_builder** now loads boundary tallies:
  - `radiation_energy_at_boundary_cum` ← `TimeIntRadiationLossAtBds` (J)
  - `particle_energy_escaped_cum` ← `particle_time_int_energy_escaped` (J)
    Note: particle escape variable does NOT include neutrons; for DT
    neutron escape we use the 14.1/17.6 fractional split. Logic prefers
    direct tally only when it's >30% of fusion (sanity check).

- **rhw_parser** cleanups:
  - Legacy `_parse_flux_limiter` retired (now thin wrapper over
    per-region parser; no more "FL varies" UserWarning).
  - `_parse_drive_temp_table` trims sentinel times (Helios pads with
    t≈1e18 s meaning "extend last value forever"). Anything ≥ 1 ms dropped.

- **helios_exodus_variable_reference.md** gained a "Boundary-tally and
  global cumulative quantities" section documenting both currently-wired
  variables and an expansion wishlist (FreqIntgRadEnLossRmax/Rmin per-boundary
  split, EnTotRadiation cross-check on in-grid U_rad, EnExchEleToRadTimeIntg
  e→rad coupling, HeatFluxAtRegionBdEle/Ion for per-region accounting,
  TimeIntFusionProd_* alpha cross-check).

### Key scientific findings (energy-ledger derived)

After all patches the ledger closes to <1% residual for VI_6 and all HDD26
variants. Direct measurements (not inferred):

| Channel (peak v, % of absorbed) | VI_6 | lrm1 | lrm2 | lrm4 |
|---|---|---|---|---|
| Absorbed (kJ)                   | 3075 | 1537 | 1571 | 1580 |
| KE inward                       | 10.1% | 14.8% | 14.9% | 15.2% |
| KE outward (blowoff)            | 51.3% | 0.0% | 0.0% | 0.0% |
| Plasma thermal                  | 32.4% | 57.5% | 56.1% | 55.6% |
| Rad escape (boundary tally)     | 5.5% | 26.7% | 28.0% | 28.2% |
| Residual closure                | 0.6% | 0.9% | 0.9% | 0.9% |

- HDD radiates 5× harder per unit absorbed than VI_6. The 27–28% rad
  escape is structurally floored for current target architecture +
  laser geometry; FL changes within 0.02–0.04 don't move it.
- HDD's hydro efficiency per absorbed (15%) is actually BETTER than
  VI_6's (10%). The HDD residual against published is in the absorbed-
  energy budget (50% vs 97%), not in coupling-efficiency.
- Energy ledger does NOT see mass-partitioning changes — same shell KE
  and v_peak with very different cold-fuel composition shows up as
  identical ledger but very different implosion physics.
- Beam 3 in HDD26 (fires 14.5–15.6 ns) contributes nothing to the
  implosion in 1D — fires after peak velocity. Geometry-fixing it
  (focus -10.22 → -0.22) changed absorbed energy +5% with zero change
  to any implosion metric.
- Indirect-drive prepulse is NOT the gap: HDD has hotter prepulse
  (130 eV / 2 ns) than VI_6 (104 eV / 2 ns). Hypothesis retired.

### Milestone: HDD26 lrm4 calibration state

lrm4 = `HDD26_DTI40_1ns130_FL04_lrm4_nb`. First HDD config to match
Imploded DT and ⟨T_hs⟩ to published within error bars.

| Metric | lrm4 | Published | Status |
|---|---|---|---|
| Imploded DT (mg) | 2.91 | 3.00 ± 0.15 | ✓ |
| ⟨T_hs⟩ (keV) | 4.62 | 4.80 ± 0.2 | ✓ |
| In-flight KE (kJ) | 239 | 300 ± 15 | -20% |
| v_peak (km/s) | 503 | 410 ± 20 | +23% |
| Adiabat | 2.12 | 6.0 ± 0.5 | -65% |
| ρR_cf (g/cm²) | 0.62 | 2.06 ± 0.10 | -70% |
| ⟨P_hs⟩ (Gbar) | 42 | 212 ± 11 | -80% |
| Frac absorbed (%) | 50.1 | 97 ± 5 | -48% |
| Yield (MJ) | 0.057 | 0.6 | -90% |

Configuration:
- FL: gas=0.06, DT ice/foam/CD=0.04
- Beam 1 spot=0.25 cm, beam 2 spot=0.22, beam 3 spot=0.18, all focus +0.22, cone 1.0°
- Pulse: foot 180 TW (4.5–10.3 ns), peak 614 TW (10.31–11.7 ns)
- Prepulse: 130 eV peak, 0–2 ns
- DT ice 40 µm

lrm progression (held FL=0.04 for lrm2/lrm4, varied spot):
- lrm1 (FL=0.02, spot 0.30): over-thermal-trap, cold fuel preserved but poor compression
- lrm2 (FL=0.04, spot 0.22): over-ablation, hot spot good but cold fuel destroyed (0.77 mg)
- lrm4 (FL=0.04, spot 0.25): SWEET SPOT — cold fuel preserved AND good compression

### Remaining residuals decompose to two independent levers

1. **Absorbed-fraction deficit (50% vs 97%)** — geometric beam-target overlap.
   Candidates: spot=0.27/0.28 (further into intermediate), cone=0.5°
   (more parallel beam to handle long focal column).
2. **Adiabat too low (2.1 vs 6.0)** — independent of beam geometry; needs
   prepulse + foot work. Candidates: Tr 130 → 150 eV, prepulse duration
   2.0 → 2.5 ns, foot power 180 → 220 TW.

Cascading residuals (P_hs, ρR_cf, yield) will follow from these two.
P_hs scales like ρR_cf², so closing α and ρR_cf is the leverage.

### Limitations of the energy-ledger diagnostic

- Sees energy distribution, NOT mass partitioning. lrm1→lrm2 and
  lrm2→lrm4 had ~zero shift in ledger but very different cold-fuel
  composition.
- Stagnation snapshot falls back to t_end when data.stag_time isn't
  set on the loaded ICFRunData (we don't run icf_analysis). For lrm2
  this samples at 16.0 ns instead of actual 14.8 ns — visible in
  output but doesn't break peak-v interpretation.
- Ideal-gas plasma thermal slightly overstates in dense ignited
  states (VI_6 stagnation residual = -0.8%, consistent with missed
  EOS-internal energy). Tolerable; could refine with EOS-derived U.
- particle_time_int_energy_escaped tally semantically ambiguous;
  not used for neutron escape (fractional 14.1/17.6 split used instead).

### Open items / minor bugs

- lrm2 fusion-yield discrepancy: EXODUS dt_neutron_count gives 48 kJ
  but summary reports 81 kJ. May be stale EXODUS file or different
  yield-computation path in icf_analysis. Worth checking which is
  authoritative.
- Filename convention drift: lrm2_nb and lrm4_nb have Burn ON in rhw.
  Either rename files or document the convention.

  **PDD anchor (validated): PDD_20_s016**
Geometry: cone 20°, spot 0.16 cm, f=0.06
Per May 2026 briefing — wins on every priority-1/2 metric:
  Yield:        59.1 MJ        (cluster 87.4, −32%)
  Neutron:      2.10×10¹⁹      (cluster 3.10, −32%)
  Peak v:       463 km/s       (cluster 470, −1.5%)
  Bang:         13.65 ns       (LILAC 13.5, +1%)
  HS ρR:        ~0.11          (1D Helios universal residual)
  ⟨ρR_cf⟩:      0.67 g/cm²     (cluster 1.10, geometry-independent)
  CR_max:       31.7           (cluster 29±3, ✓)
  Adiabat:      1.05           (cluster 3.0±0.3, low — 1D residual)
  Frac abs:     83.9%          (cluster 65±9.3%, slightly over)
  T_hs/P_hs:    overshoot (de-prioritized — Helios burn physics may differ
                from LILAC/xRAGE/HYDRA)

Note: prior CLAUDE.md / memory references to "PDD_26b" as the calibrated
state refer to an April 2026 vintage that was superseded by the s016
calibration in May 2026. PDD_20_s016 is the current authoritative anchor.

## Session 2026-05-12 — Energy-balance diagnostic, HDD lrm4 milestone, PDD anchor correction

### PDD anchor (validated): PDD_20_s016 — supersedes PDD_26b
Geometry: cone 20°, spot 0.16 cm, f=0.06 (May 2026 briefing)
Wins on every priority-1/2 metric vs LILAC/xRAGE/HYDRA cluster:
  Yield 59.1 MJ (cluster 87.4, -32%) | Neutron 2.10e19 (-32%)
  Peak v 463 km/s (-1.5%) | Bang 13.65 ns (+1%)
  HS ρR ~0.11 | ρR_cf 0.67 g/cm² (cluster 1.10) | CR_max 31.7 (cluster 29±3) ✓
  Adiabat 1.05 (cluster 3.0±0.3, low — 1D residual)
  Frac abs 83.9% (cluster 65±9.3%, slightly over)
  T_hs/P_hs overshoot — deprioritized; Helios burn physics may differ from LILAC.
NOTE: HS ρR ~0.11 and ρR_cf 0.67 stated as geometry-independent 1D Helios
residuals — implies a structural 1D floor below cluster values. Worth
cross-checking against HDD residuals.

### HDD26 lrm4 milestone
File: HDD26_DTI40_1ns130_FL04_lrm4_nb
Path: ~/Sims/Xcimer/HDD_26/HDD26_DTI40_1ns130_FL04_lrm4_nb/
Config: FL gas=0.06, ice/foam/CD=0.04. Beam 1 spot=0.25, beam 2 spot=0.22,
beam 3 spot=0.18, all focus +0.22, cone 1.0°. Pulse foot 180 TW (4.5-10.3),
peak 614 TW (10.31-11.7). Prepulse 130 eV / 0-2 ns. DT ice 40 µm.

FIRST HDD config to match published within error bars:
  Imploded DT: 2.91 mg vs 3.00 ± 0.15 ✓
  ⟨T_hs⟩: 4.62 keV vs 4.8 ± 0.2 ✓

Remaining residuals:
  v_peak 503 vs 410 (+23%)
  Adiabat 2.12 vs 6.0 (-65%) — independent of beam geometry
  ρR_cf 0.62 vs 2.06 (-70%)
  ⟨P_hs⟩ 42 vs 212 Gbar (-80%)
  Frac absorbed 50% vs 97% (-48%)
  Yield 0.057 vs 0.6 MJ (-90%)

Two independent levers for remaining work:
  1. Absorbed-fraction (geometry — spot/cone/focus)
  2. Adiabat (prepulse Tr, prepulse duration, foot power)

### Energy-balance diagnostic (examples/energy_balance_diagnostic.py)
Standalone side-by-side ledger using build_run_data (not ICFRunBuilder).
Tracked channels: KE inward, KE outward, plasma thermal (ideal-gas
(3/2)(P_i+P_e)V), in-grid radiation (3·P_rad_true·V — uses
rad_pressure_true not 3-component rad_pressure), absorbed laser
(cumulative), fusion released (split into alpha-deposited + neutron-escaped
via 14.1/17.6 fractional split since EXODUS particle_time_int_energy_escaped
does NOT include neutrons), direct rad escape from boundary tally.
Closure equation: E_absorbed + E_alpha = Σ_in-plasma + E_rad_boundary + residual.
Closes to <1% residual across VI_6 and all HDD26 variants. Runs as:
  python3 examples/energy_balance_diagnostic.py <run1_base> <run2_base>
Output: console snapshot tables at peak-v/stag/end + 3-page PDF.

KNOWN LIMITATIONS:
- Sees energy distribution, NOT mass partitioning. lrm1→lrm2→lrm4 had
  ~zero shift in peak-v ledger but very different cold-fuel composition.
- Stagnation snapshot falls back to t_end when data.stag_time isn't on
  loaded ICFRunData. For runs with stag < t_end (e.g. lrm2 stag=14.8,
  t_end=16.0), snapshot is mislabeled. Peak-v snapshot is correct.
- Plasma thermal uses ideal gas — slightly overstates in dense ignited
  states (VI_6 stag residual -0.8%, consistent with EOS-internal miss).

### Key HDD vs VI_6 ledger findings (energy distribution at peak v, % of absorbed)
                                VI_6   lrm1   lrm2   lrm4
Absorbed (kJ)                   3075   1537   1571   1580
KE inward                       10.1%  14.8%  14.9%  15.2%
KE outward (blowoff)            51.3%   0.0%   0.0%   0.0%
Plasma thermal                  32.4%  57.5%  56.1%  55.6%
Rad escape (boundary, measured)  5.5%  26.7%  28.0%  28.2%
Residual closure                 0.6%   0.9%   0.9%   0.9%

- HDD hydro efficiency per absorbed (15%) is BETTER than VI_6 (10%).
- HDD's residual against published is in ABSORBED ENERGY (50% vs 97%),
  not coupling-efficiency.
- HDD corona radiates 5× harder per absorbed than VI_6's. Structural
  floor for current target+geometry+FL class — FL changes within
  0.02-0.04 don't move it.
- Beam 3 contributes nothing in 1D (fires 14.5-15.6 ns, after peak-v).
- Indirect drive prepulse NOT the lever — HDD has 130 eV vs VI_6's 104 eV.

### Code changes this session (all committed/pushed)
- data_builder loads radiation_energy_at_boundary_cum and
  particle_energy_escaped_cum from EXODUS.
- rhw_parser: legacy _parse_flux_limiter retired (thin wrapper);
  drive-temperature sentinel filter (trims t ≥ 1 ms).
- helios_exodus_variable_reference.md: new "Boundary-tally" section
  with wishlist of future accounting variables.
- examples/energy_balance_diagnostic.py: new standalone tool.

## May 2026 Appendix — 1D Geometric Coupling Correction (CLOSED)

### Motivation

Helios is a 1D spherical hydrodynamics code with a laser ray-trace that
models the beam as a single cone of half-angle θ from a focal point at
distance d.  By construction the cone always intercepts the spherical
target — Helios cannot reproduce the geometric coupling reduction that
occurs in 3D PDD/HDD configurations, where fixed-pointing beams
increasingly miss the imploding capsule as it shrinks.  This appendix
documents the empirical investigation of whether that 1D limitation is
the binding constraint blocking the cluster match.

### Tools added

- `examples/pdd_geometric_coupling.py` — pre-processor that computes
  `f_geom(t) = eta(t) / eta(anchor)` from a baseline no-burn EXODUS run.
  Uses view-factor formula `(1 − √(1 − (R/d)²)) / (1 − cos θ)` in the
  partial-coupling regime and f = 1 when R/d ≥ sin θ (sphere fills cone)
  or R/d ≥ 1 (focal point inside corona).  Three R(t) source options:
  `r_crit`, `r_outer`, `r_ablation`.  Use `--r-source r_ablation` for
  PDD calibration — it tracks the cold-shell trajectory rather than the
  laser-coupling surface.  Writes `<base>_geomcorr_power.csv` and a
  diagnostic PDF.

- `examples/apply_geomcorr_to_rhw.py` — rhw patcher that injects the
  corrected pulse table back into the rhw.  Locates `[Laser Source
  Data]:` block, walks `Parameters for beam:  = N` headers, finds each
  `[table format=2]:` block, **resamples the pulse table from its
  original 6 control points to 200** (configurable via `--n-points`) and
  rewrites three rhw lines per beam: `# table rows = N`, times row, and
  powers row.  Resampling is critical because the original 6-point
  table uses piecewise-linear interpolation between control points; a
  single-point per-control modification creates artificial linear ramps
  across the peak plateau that don't reflect the actual f_geom shape.
  Resampling correctly captures the f_geom variation under Helios's
  linear-interp model.  Round-trip energy matches the analytical
  `∫(P_orig × f_geom) dt` to within 0.03%.

### rhw file format (mapped this session)

The Olson_PDD rhw layout is (line numbers approximate for an
unmodified PDD_35 rhw):

```
[Header Data]:
[Geometry Data]:
[Spatial Grid Data]:        <-- many [table format=2] sub-blocks
[Hydro Data]:               <-- "Allow free expansion at Rmax" flag
[Rad Trans Data]:
[Rad Source Data]:          <-- prepulse Tr table for X-ray drive
[Laser Source Data]:        <-- starts ~line 808 for PDD_35
   Number of laser beams = 3
   Parameters for beam:          = 1     <-- note colon after "beam"
     Wavelength = ...
     Spot size = ...
     Half cone angle = ...
     Focus position = ...
     [table format=2]:    Time-dependent laser powers:
       table is 3D  = 0
       # table rows = K      <-- number of (time, power) pairs
       # table cols = 2
       <K times in seconds, space-separated>
       <K powers in TW, space-separated>
   Parameters for beam:          = 2
   Parameters for beam:          = 3
[End Laser Source Data]
[Ion Beam]:
[MHD Data]:
[Time Control Data]:
[Output Control Data]:
[Atomic Processes Parameters]:
[Dialog Settings Data]:
[End of Workspace File]
```

The beam-header colon syntax (`Parameters for beam:  = N`, with colon
*after* "beam") is non-obvious — original assumption was no colon, which
caused the patcher to silently locate zero pulse tables on its first
run.  Both rhw_parser and apply_geomcorr_to_rhw now accept either form.

### Helios CLI on macOS

Helios is installed at `~/Codes/Prism/Helios_11.0.0/`.  On macOS the
`.app` is a bundle — the actual binary lives inside.  Invocation:

```bash
~/Codes/Prism/Helios_11.0.0/Helios.app/Contents/MacOS/Helios -b \
    -i /path/to/run.rhw \
    -d /path/to/sims_root \
    -o run_name \
    -x
```

`-b` runs headless, `-d` sets the results-directory root, `-o` sets
the run folder name, `-x` overwrites if it already exists.

### PDD_20_s016 (narrow cone) — pre-processor only

Geometry: θ = 20°, d = 0.22 cm, spot = 0.16 cm, λ = 0.350 µm.
Miss threshold d·sin(θ) = 753 µm.

| Metric                          | Value   |
|---------------------------------|---------|
| f_geom at peak power (t = 12 ns)| 1.000   |
| Mean f_geom over laser-on       | 0.971   |
| Min f_geom over laser-on        | 0.384   |
| Integrated ΔE_delivered         | −3.7 %  |

The ablation front sits well above 753 µm during peak power — the
correction kicks in only on the late tail of the main pulse.  Pre-
processor only; no Helios round-trip needed (too small to matter).

### PDD_35_s016 (wide cone) — full Helios round-trip

Geometry: θ = 35°, d = 0.22 cm, spot = 0.16 cm, λ = 0.350 µm.
Miss threshold d·sin(θ) = 1262 µm.

Pre-processor: f_geom at peak power = 1.000 (t = 10.7 ns), mean f_geom
on-pulse = 0.897, integrated ΔE_delivered = −14.9 %.

Helios round-trip (corrected rhw, no-burn, Helios 11.0.0):

| Metric                       | Baseline   | Corrected  | Δ       | Cluster   |
|------------------------------|-----------:|-----------:|--------:|----------:|
| Laser energy on target (MJ)  | 1.656      | 1.452      | −12.3 % | 1.40      |
| Spot–sphere coupling (%)     | 79.6       | 82.0       | +2.4 pp | 65 ± 9    |
| Peak velocity (km/s)         | 450.6      | 443.3      | −1.6 %  | 470       |
| Mass-avg adiabat at peak v   | 1.79       | **1.07**   | **−40%**| 3.0 ± 0.3 |
| Base adiabat at breakout     | 0.39       | 0.39       | 0 %     | —         |
| ⟨ρR_cf⟩ (g/cm²)              | 0.96       | 0.76       | −21 %   | 1.10      |
| ⟨P_hs⟩ (Gbar)                | 90.6       | 56.3       | −38 %   | 193 ± 20  |
| CR_max                       | 38.3       | 33.4       | −13 %   | 29 ± 3    |
| Imploded DT mass (mg)        | **0.60**   | **0.60**   | **0 %** | ~3.0      |
| Hydro efficiency (%)         | 10.8       | 11.4       | +6 %    | —         |
| Yield no-burn (MJ)           | 0.316      | 0.170      | −46 %   | —         |

### Key findings

1. **Pre-processor prediction validated.**  Predicted −14.9 %, Helios
   delivered −12.3 %.  The 2.6 pp slippage is a second-order shift in
   Helios's spot–sphere coupling factor (79.6 → 82.0 %): the slightly
   less violent implosion keeps the target larger longer and captures
   more of the Gaussian spot.  Real second-order coupling effect, not
   a numerical artifact.

2. **Foot-pulse-set metrics are structurally outside the correction's
   reach.**  Base adiabat at breakout (0.39 → 0.39 exactly), imploded
   DT mass (0.60 → 0.60 mg exactly), and shock breakout time (9.55 ns
   both) are all set during the foot pulse where f_geom = 1, so they
   are unaffected by the correction by construction.

3. **Peak velocity barely moves (−1.6 %).**  Confirms the structural
   argument: f_geom = 1.000 at peak laser power for both anchors, so
   the rocket-equation work integral on the imploding shell is
   essentially unchanged.  Residual comes from reduced late-tail
   coasting drive.

4. **Compression channel absorbs the full reduction.**  ⟨P_hs⟩ −38 %,
   ⟨ρR_cf⟩ −21 %, yield −46 %, peak density −25 %.  **Adiabat falls
   from 1.79 to 1.07 — moving the wrong way relative to cluster 3.0
   — because the correction reduces late-pulse re-heating of the
   assembling cold fuel, leaving it softer than baseline.**

5. **The correction lands delivered energy on the cluster mark.**
   1.452 MJ vs cluster 1.40 MJ.  Despite this, the binding calibration
   residuals (adiabat shortfall, ρR_cf shortfall, imploded-mass
   shortfall, yield gap) survive intact or worsen.

### Disposition

The 1D ray-trace is not the limiting physics for matching the 3D-code
cluster.  The compression-deficit chain — shell-fuel mass distribution
at stagnation, EOS, drive-symmetry effects not captured in 1D —
remains the binding constraint.

**Options 2 and 3 (in-Helios source-code modifications) are not
pursued.**  The pre-processor + rhw patcher pipeline is retained as
an end-to-end diagnostic for future targets where the geometric miss
could be larger (wider cones, shorter focal distances, more aggressive
implosions that compress the ablation surface below the miss threshold
earlier in the pulse).  The Option 1 test can be re-run in minutes on
any new candidate geometry.

**Anchor-selection circular logic resolved.**  The natural concern
("the anchor was selected without the correction in the loop, biasing
selection toward geometries where the correction is small") was
addressed by running the correction on a wide-cone anchor where the
geometric miss is 6× larger than at the narrow-cone PDD_20 anchor.
The corrected wide-cone result does not unlock a better cluster
match — it just shifts which subset of metrics agrees.  Re-running
the calibration loop with the correction applied would shift the
chosen anchor but would not close the compression residuals.

### Reference document

Full memo with derivation, four-options analysis, validation
criteria, and disposition: `docs/Helios_1D_Coupling_Correction_Memo.docx`.

## Open items / Next priorities (post-May 23 2026 PDD closeout)

In order of expected information yield:

0. **DONE (May 23): Energy-ledger comparison fab007 burn-ON vs fab007 burn-OFF.**
   Confirmed pre-alpha hydro state is identical at peak velocity to <0.1%. Alpha
   bootstrap when it catches deposits 5.2 MJ (3.4× the absorbed laser); when it
   doesn't, deposits 50 kJ (3% of absorbed). 86% of alpha-deposited energy goes
   to outward blowoff KE, only 5% sustained as plasma thermal at stagnation.
   Closure 0.4% at peak v, 0.5% at stagnation. Validates the diagnostic
   toolset for HDD work. PDF and console table at:
   `~/helios_postprocess/Olson_PDD_20_fab007_foot25_s018_c37_burn_vs_Olson_PDD_20_fab007_foot25_s018_c37_nb_energy_balance.pdf`

0a. **DONE (May 24 2026): foam burn-propagation diagnosis.** Three-part
    diagnostic series resolved priority 0a: the May 23 "burn-front arrest
    at the ice/foam interface" framing is wrong. Findings (full table and
    disposition in the "Corrected diagnosis (May 24 2026 final)" PDD
    section above):

    - The ice/foam interface is essentially continuous at peak burn in
      fab007 (T_e ratio 1.017, T_rad 1.000, α-dep ratio 1.06).
      `examples/dump_burn_propagation_profile.py`.
    - The Z̄ jump is only 8% (foam is 92% DT by mass — 0.222 g/cc DT in a
      0.020 g/cc CH structural binder, not a high-Z mixed material).
    - PROPACEOS foam EOS is NOT broken — foam peak CR equals or exceeds
      ice peak CR in both fab02 (1.05×) and fab007 (0.91×) at peak
      compression times within 30-40 ps. `dump_burn_rate_timing.py`.
    - The residual IS alpha-bootstrap STRENGTH at bang time. fab02
      ignites the ice strongly (Ti 34 keV / P 408 Gbar at bang),
      bootstrap fires 68% of foam zones (4.4 keV avg), foam contributes
      26.5% of yield. fab007 marginally ignites (Ti 17 keV / P 116 Gbar),
      bootstrap fires 21% of foam (1.85 keV avg), foam contributes 10%.
      `dump_per_zone_burn_share.py`.
    - HDD calibration risk flag updated from "fundamentally limited" to
      "conditional on landing in a robust-bootstrap regime." Use fab02
      settings (cone=20°, spot=0.16, FL_DT=0.06, FL_foam+CH=0.02) as the
      primary HDD-transfer geometric baseline.

    Residual open question (deferred, NOT a calibration-blocker): even
    fab02's foam burns ~18× less per unit DT mass than its ice (3811 vs
    70064 kJ/mg), vs LILAC's implied near-parity. Likely (i) Helios 1D
    non-local α-transport details at Z̄≈3, (ii) ~0.3 ns burn FWHM too
    short for full foam burnup, (iii) over-drive gives hot ice rather
    than LILAC-like compression history. Not worth pursuing inside
    helios_postprocess until HDD calibration converges.

1. **HDD lrm4 with 150 eV prepulse** (Tr 130 → 150 eV at lrm4
   geometry, no-burn).  The cleanest isolated adiabat-lever test we
   can run on HDD.  The lrm4 free-Rmax baseline already matches
   imploded DT mass and ⟨T_hs⟩; the residual is α = 2.12 vs target 6.
   If stiffer prepulse lifts α toward the target without breaking the
   other matched metrics, the prepulse strength is a real knob.  If
   not, the adiabat shortfall is structural to Helios's foot-shock
   physics.  Path: `~/Sims/Xcimer/HDD_26/HDD26_DTI40_1ns150_FL06_lrm4_nb/`
   (run name pattern).  Modify the prepulse temperature in the rhw's
   `[Rad Source Data]:` block — specifically the time-dependent drive
   temperature table at Rmax.

2. **HDD lrm4-at-FL-0.06 inverse-FL test.**  Revert the FL=0.04 on
   ice/foam/CD shell back to 0.06 (the PDD baseline value) via
   `examples/edit_rhw_flux_limiter.py`.  Prediction: yield and energy
   ledger essentially unchanged, confirming FL is a non-lever in HDD
   as it was in PDD.  Quick test (1 hr of compute, 5 min of analysis).

3. **PDD compression-channel investigation.**  EOS sensitivity (SESAME
   vs PROPACEOS vs alternative tables), foot-pulse adiabat physics,
   shock-timing sensitivity to the foot/ramp transition.  The
   compression deficit is the binding residual across all calibration
   anchors — both PDD_20_s016 (geomcorr-confirmed-clean drive
   delivery, still α = 1.05) and PDD_35_s016 corrected (delivery on
   cluster mark, still α = 1.07).  The lever is somewhere in the
   foot/peak EOS-shock chain that's identical at both anchors.

4. **HDD lrm5 spot=0.27 at FL=0.04** — absorbed-fraction monotonicity
   test.  Deprioritized pending (1) and (2); will pick up later if
   the prepulse lever does not close the adiabat gap.

## Session Update — June 2 2026 (WT_cthomas HDD calibration breakthrough + retractions)

### Production HDD calibration: `WT_cthomas_baseline_picket_012`

`WT_cthomas_baseline_picket_012` (FL=0.012, EOS-fixed CD ablator, picket
modification) matches Thomas Vulcan HDD reference on the four headline
structural metrics:

| Metric | picket_012 | Thomas | Δ |
|---|---:|---:|---:|
| **Yield (MJ)** | **256.4** | **256** | **+0.1%** ✓ |
| **Hydro efficiency (%)** | **8.3** | **8.0** | **+4%** ✓ |
| In-flight KE (kJ) | 258 | 300 ± 15 | -14% ✓ |
| ⟨ρR_cf⟩ (g/cm²) | 1.17 | 1.60 ± 0.08 | -27% close |
| RHINO V_impl (km/s) | 389 | 410 | -5% ✓ |
| Foot shock (ns) | 7.88 | — | 2-shock train emerges |
| Ramp shock (ns) | 11.45 | — | — |
| Base adiabat at breakout | 1.08 | ~6 | tripled vs picket-OFF (0.35) |
| CR_max | 53 | ~20 | persistent 1D residual (still over) |
| Imploded DT mass (mg) | 0.05 | 3.0 | persistent 1D residual (-98%) |
| ⟨T_hs⟩ (keV) | 120 | 47 | persistent 1D over-shoot (+157%) |

This is the HDD analog of `Olson_PDD_20_fab007` for PDD — structural
metrics within 15% of reference, persistent 1D-vs-3D residuals in the
hot-spot fine structure that aren't fixable inside Helios. Production
calibration anchor going forward.

Run path (Studio): `~/Sims/Xcimer/HDD_26/HDD_scan/WT_cthomas_baseline_picket_012/`

Three settings define this anchor:
1. **CD ablator EOS path resolved to Studio-local PROPACEOS table** (NOT
   Will's hard-coded Windows OneDrive path — see CD EOS bug below)
2. **Flux limiter = 0.012 in all 5 regions**
3. **Picket modification on the Trad table at Rmax** (peak still 135 eV;
   the picket SHAPE / TIMING was modified to triple the base adiabat at
   breakout — exact .rhw diff to be documented after one more validation
   pass)

### CD EOS Windows-path fallback: ±50% yield bug

The original `WT_cthomas_baseline.rhw` from Will Trickey contains a
hard-coded Windows OneDrive path for the CD ablator EOS:

```
EOS filepath = C:/Users/WilliamTrickey/OneDrive - Xcimer Energy Corporation/Documents/tables/TMelhorn_PROPACEOS/CD.prp
```

On Mac Studio (and likely on any non-Windows host), Helios prints:

```
The EOS file C:/Users/WilliamTrickey/...CD.prp does not exist.
Check spatial regions for valid EOS filepath.
```

…then **silently falls back to an ideal-gas-like default**. The
simulation completes without error. The downstream impact is structural:

| | Will baseline (fallback) | baseline_tm (proper EOS) | Δ |
|---|---:|---:|---:|
| RHINO V_impl (km/s) | 396.7 | 395.3 | -0.4% (identical) |
| Shock breakout, base adiabat, ablation P | — | — | all identical |
| Peak density (g/cc) | 102,971 | 31,062 | **-70%** |
| Burn FWHM (ns) | 0.065 | 0.082 | +27% |
| **Yield (MJ)** | **196** | **295** | **+50%** |

In-flight phase is identical; the entire effect lives at stagnation,
where the CD ablator EOS controls the bounce. Soft (ideal-gas-fallback)
ablator → easy compression → high peak ρ but short burn → low integrated
yield. Stiff (tabulated) ablator → moderate peak ρ but sustained burn →
high integrated yield.

**Implications:**

- Will's published `WT_cthomas` yield of 196 MJ likely used the fallback
  EOS on his own machine too — his RHINO postprocess of his own .exo
  also reports 196 MJ, matching what fallback-Helios produces.
- Before quoting any yield number for a WT-derived run, grep the .rhw
  for `.prp` paths and confirm Helios resolves all of them. The script:

  ```bash
  grep -nE "EOS filepath|Opacity filepath" <run>.rhw | grep -v '^\s*#'
  ```

  Any path starting with `C:/` will fall back silently on macOS/Linux.
  Remap to a Studio-local table.

### RHINO convention validation status

Implosion velocity convention (W. Trickey) is **validated and production-ready**:

- Per-timestep shell = zones with `rho > rho_peak/e`
- v_shell = sqrt(2 KE_shell / m_shell)
- Implosion velocity = first significant local maximum of v_shell pre-stagnation

Helios reproduces Will's reported V_impl to within 3–5% across every
WT_cthomas configuration tested (baseline, baseline_tm, baseline_012,
picket_012). Available as `data.implosion_velocity_rhino_kms` and in the
comparison table as `Implosion velocity RHINO (km/s)`.

Min shell adiabat convention is **VALIDATED (June 2 2026, commit 5704bca)**:

Clean-room re-implemented from RHINO source (wtrickey27/RHINO,
private; TM has read access). The previous implementation was off by
~27× from Will's reported 4.13 because of three independent bugs:

1. **Fermi pressure formula.** Old code used Lindl convention
   `P_F = 2.17 × (ρ/0.205)^(5/3) Mbar` — that's a *normalized*
   adiabat (α=1 at cryo DT solid density), not the actual degenerate
   electron-gas Fermi pressure. At ICF densities it's ~14× larger
   than the proper formula. New code uses
   `P_F = (3π²)^(2/3)/5 × ℏ²/m_e × n_e^(5/3)` with `n_e` from
   `data.electron_density` per zone (RHINO's `partially_ionized`
   mode — correct for multi-material shells).

2. **Shell inner threshold.** Old code used 1/e always. New code
   uses RHINO's time-dependent default — 1% of peak density
   *pre-breakout*, 1/e *post-breakout*. Critical for the t=0
   shell anchor: at t=0 in WT_cthomas, only the dense CH ablator
   passes 1/e of peak; with the 1% threshold, the foam + ablator
   together comprise the shell, giving R_inner_0 at the gas/foam
   interface (0.1875 cm) as Will intends.

3. **CR=1.5 trigger.** Old code used "innermost zone of
   rho>peak/e mask" radius. New code uses Will's
   `cr_inner = shell_inner(t=0) / shell_inner(t)` reaching 1.5, with
   the time-dependent shell_inner above. Breakout time itself is
   computed per Will's definition: the time when peak-density
   position first crosses shell_inner_0 derived from a constant 1%
   threshold.

Validation on `WT_cthomas_baseline_tm`: **5.64 vs Thomas 6.00 (-6%)**.

Will's reported 4.13 was almost certainly computed on his original
fallback-EOS run (the WT_cthomas_baseline.rhw with the Windows path
that doesn't resolve on Mac). Pending validation that our same
algorithm on the fallback-EOS .exo gives ~4.13 too -- expected
behaviour since the fallback EOS gives -50% on yield and similar
shifts on other compression-channel metrics.

**Other adiabat methods in this file** (`_compute_adiabat`,
`_compute_adiabat_at_breakout`, `_compute_cr15_metrics`) retain
the Lindl convention they always used. Changing them would shift
the documented PDD/HDD calibration history by ~14×. The RHINO
method is the only cross-tool-comparable one; the Lindl-convention
methods stay for historical continuity. Future analyses should
prefer the RHINO method (`data.adiabat_min_rhino`) for any
comparison to RHINO-postprocessed reference simulations.

Audit attributes on data: `t_breakout_rhino_ns`,
`r_inner_initial_cm`, `r_inner_at_cr15_cm` — appear in the
`Adiabat min (RHINO convention, ...)` log line for sanity-checking.

### Retired diagnoses (incorrect — do not pursue)

1. **"FL=0.012 + zoning mismatch couples to crash Helios hydro"** —
   WRONG. Heap-corruption / malloc abort errors on `WT_cthomas_fl012.rhw`
   and `WT_cthomas_picket_v1.rhw` were caused by the CD EOS path issue,
   not by the FL value or the zoning. FL=0.012 runs cleanly when CD EOS
   is properly resolved. User pointed out the FL-causes-malloc story made
   no physical sense; they were right.

2. **"Foam_2 / CD_2 bisections introduce dMass discontinuities"** —
   WRONG. Computed dMass(i)/dMass(i-1) ratios at every region boundary:

   | Boundary | Type | Ratio |
   |---|---|---:|
   | 0 → 1 (gas → gas) | spherical-geometry first zone | **7.0×** (the MAX Helios complained about) |
   | 58 → 59 (foam_1 → foam_2 bisection) | clean | 0.92× |
   | 125 → 126 (foam_2 → CD ablator) | material interface | **0.577×** (the MIN) |
   | 155 → 156 (CD_1 → CD_2 bisection) | clean | 0.885× |

   The 7× max is intrinsic to spherical geometry at r=0 (full sphere of
   width Δr vs first shell of similar Δr; ratio = 7 regardless of FL or
   bisection). The bisections introduced by Will for diagnostic purposes
   are essentially perfect (~0.9× ratios). The "Poor mass matching" warning
   was a red herring; it's present in baseline too and doesn't correlate
   with the crash.

### Validation run on tap: `WT_cthomas_baseline_picket_015`

Same picket as picket_012, FL=0.015 instead of 0.012. Tests whether the
picket alone (without tighter FL) is sufficient for the Thomas yield
match, or whether the FL=0.012 + picket combination is what's needed.

Expected results:
- If yield ~256 MJ: picket is the primary lever, FL secondary
- If yield ~295 MJ (like baseline_tm): FL is essential to reach Thomas

Either result is informative. User to run.
