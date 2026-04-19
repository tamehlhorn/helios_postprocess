# helios_postprocessor -- Project Guide

## Repository
- **GitHub**: `tamehlhorn/helios_postprocessor`
- **Package**: `helios_postprocess/`
- **Version**: 3.0.0 (March 2026)
- **Dev machine**: MacBook (`~/Codes/helios_postprocessor`) -- editing, pushing; use `python` not `python3`
- **Run machine**: Mac Studio (`tommehlhorn`, `~/helios_postprocessor`) -- use `python3`
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
python3 ~/helios_postprocessor/examples/run_analysis.py <base_path>

# MacBook
python ~/Codes/helios_postprocessor/examples/run_analysis.py <base_path>
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
    "hs_radius_ignition_um":[120.0, 0.0]
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

**Key findings:**
- Hydro efficiency locked at ~10.5% regardless of geometry -- excess drive is in ablation physics
- Absorption insensitive to cone angle/spot (87±3%) -- refraction controls deposition, not geometry
- Narrower cone → lower adiabat (deeper deposition) but minimal velocity effect
- PDD_27 failure: d=5.0 cm (rays miss entirely) + norm=0.5 (half power) -- two compounding bugs
- PDD_28 peak velocity was artifactually high due to velocity bug (now fixed); rerun needed

**Next runs recommended:**
- PDD_29: spot=0.02, cone=35°, uniform, d=0.23 cm, norm=1.0 -- rerun PDD_28 with fixed postprocessor
- PDD_30: spot=0.02, cone=35°, uniform, d=0.20 cm, norm=1.0 -- compare focus position effect
- If geometric approach insufficient: reduce peak power toward 200 TW, or increase flux limiter
  above f=0.06 to reduce ablation efficiency

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
Using the `att < 9e5` fallback (NumElecDensity not in VI_6.exo):

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
| Critical radius | Innermost r where `NumElecDensity >= n_crit`, OR innermost r where `att < 1e20` |

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
python3 ~/helios_postprocessor/examples/run_analysis.py \
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

### Diagnostic Scripts (repo root)
| Script | Purpose |
|---|---|
| `plot_adiabat_shock.py` | Adiabat history + first shock plots (any target with DT ice layer) |
| `plot_laser_deposition.py` | Laser energy deposition spatial profile at multiple timesteps |
| `plot_laser_intensity.py` | I(r,t) reconstruction using Methods 1 & 2; 3-page PDF with cross-check |
| `helios_postprocessor_guide.docx` | User guide for collaborators |

Standard usage (all scripts accept .exo path as first positional arg):
```bash
python3 ~/helios_postprocessor/plot_laser_intensity.py \
  ~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6/VI_6.exo
python3 ~/helios_postprocessor/plot_laser_intensity.py \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_2021_01a/Olson_PDD_2021_01a.exo \
  --wavelength_um 0.351 --ntimes 6
```

### Critical-surface location as cross-check
Two independent ways to locate the critical surface:
1. Electron density: innermost r where `NumElecDensity >= n_crit`
2. Attenuation: innermost r where `att < 1e20` (i.e. not an opaque-sentinel zone)

These should agree to within one zone. If they don't, that points at either (a) a
subtle issue in how Helios handles the critical-surface turning point, or (b) a
non-LTE electron population where the simple plasma-frequency cutoff doesn't apply.

### Script implementation lessons (April 2026)
Hard-won during debugging the April 2026 VI_6 run:

1. **r_crit fallback when NumElecDensity is absent**: use `opaque.max() + 1`
   (outermost opaque zone, step one outward), NOT `transparent.min()`. The 1e30
   sentinel tags only the turning-point region, not every inner zone, so the
   "innermost transparent zone" is some arbitrary central zone well below r_crit.
   VI_6 observed failure: `transparent.min()` returned 0.026 cm vs. actual
   r_crit of 0.2 cm (8x off).

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
- Integrate I-profile output into ICFAnalyzer as an optional page in the main PDF.

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
    python3 ~/helios_postprocessor/examples/run_analysis.py <base_path>
    python3 ~/helios_postprocessor/examples/run_analysis.py <base_path> --contours

    # MacBook
    python ~/Codes/helios_postprocessor/examples/run_analysis.py <base_path>
    python ~/Codes/helios_postprocessor/examples/run_analysis.py <base_path> --contours
