# helios_postprocessor -- Project Guide

## Repository
- **GitHub**: `tamehlhorn/helios_postprocessor`
- **Package**: `helios_postprocess/`
- **Version**: 3.0.0 (March 2026)
- **Dev machine**: MacBook (`~/Codes/helios_postprocessor`) -- editing, pushing; use `python` not `python3`
- **Run machine**: Mac Studio (`tommehlhorn`, `~/helios_postprocessor`) -- use `python3`
- **MacBook python**: `/Users/mehlhorn/anaconda3/bin/python` (NOT `/usr/local/bin/python3`)
- **Package install on MacBook**: `pip install -e . --user` (requires `--user` due to Anaconda permissions)

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

After the pipeline runs, burn-averaged metrics can be computed for comparison
with published ICF target designs:

```python
from helios_postprocess.burn_averaged_metrics import (
    extract_histories_from_run_data,
    calculate_burn_averaged_metrics,
    compare_with_published,
)

histories = extract_histories_from_run_data(data)   # computes from 2D arrays
metrics = calculate_burn_averaged_metrics(histories)
print(compare_with_published(metrics, published_data, laser_energy_MJ=4.0))
```

`extract_histories_from_run_data()` computes per-timestep hot-spot averages
from the 2D arrays on ICFRunData using `region_interfaces_indices`:
- Mass-averaged ion temperature (keV) over hot-spot zones
- Mass-averaged total pressure (Gbar) over hot-spot zones
- Mass-averaged density (g/cm3) over hot-spot zones
- Hot-spot outer radius from `zone_boundaries[:, ri[:, 0]]`
- Cold-fuel rhoR from `areal_density_vs_time` (requires `analyze_burn_phase()`)

Also computes implosion metrics:
- In-flight KE (inward-moving shell only, max over time)
- Hydrodynamic efficiency (max KE_inward / E_absorbed)
- Fraction absorbed (laser_energy_deposited / integrated laser_power_delivered)
- Imploded DT mass at stagnation (zone-index method, consistent with mass fractions)
- IFAR (density-based shell boundaries)
- CR_max (stagnation convergence ratio, passed through from icf_analysis.py)

Yield, laser energy, and gain come from Helios's own time-integrated
quantities -- NOT re-integrated from sampled EXODUS data. See Physics
Convention #8 below.

### Additional Physics Modules (functional style)

| Module | Purpose |
|--------|---------|
| `energetics` | Kinetic energy, hydro efficiency, PdV work |
| `neutron_downscatter` | Down-scatter ratio (DSR) diagnostics |
| `pressure_gradients` | Shock identification, RT instability assessment |

These modules are standalone and not yet wired into ICFAnalyzer -- see Open Items.

### Optional

`RHWParser` reads `.rhw` input files for drive configuration (direct/indirect,
burn on/off, drive temperature profile).

## Active Source Files

| File | Lines | Role |
|------|-------|------|
| **Pipeline** | | |
| `core.py` | ~965 | `HeliosRun` -- EXODUS/netCDF4 reader |
| `data_builder.py` | ~555 | Bridge: `HeliosRun` -> `ICFRunData` dataclass |
| `icf_analysis.py` | ~1340 | `ICFAnalyzer` -- drive, stagnation, burn, implosion, IFAR, mass fractions, burn propagation, convergence ratios |
| `icf_plotting.py` | ~1580 | `ICFPlotter` -- full PDF report |
| `icf_output.py` | ~435 | `ICFOutputGenerator` -- summary text + CSV time histories |
| **Physics modules** | | |
| `burn_averaged_metrics.py` | ~560 | Temporal burn-averaging, implosion metrics, published-data comparison |
| `energetics.py` | -- | Kinetic energy, hydro efficiency, PdV work |
| `neutron_downscatter.py` | -- | Neutron down-scatter ratio diagnostics |
| `pressure_gradients.py` | -- | Pressure gradient analysis, shock ID, RT assessment |
| **Support** | | |
| `rhw_parser.py` | ~240 | `RHWParser` -- reads `.rhw` input files |
| `__init__.py` | ~90 | Package exports |
| **Runner** | | |
| `examples/run_analysis.py` | ~350 | CLI runner with auto file derivation and comparison PDF |

## Archived Files

`archive/legacy_modules/` contains 3 standalone modules from v2.0 whose physics
is fully incorporated into `icf_analysis.py`:

- `areal_density.py` -- superseded by `ICFAnalyzer._compute_areal_densities()`
- `burn.py` -- superseded by `ICFAnalyzer.analyze_burn_phase()`
- `hot_spot.py` -- superseded by `ICFAnalyzer._compute_hot_spot_properties()`

## Published Data Comparison (JSON format)

Place a `<name>_published.json` file next to the `.exo` file. Format:

```json
{
    "laser_energy_MJ": 4.0,
    "T_hs":    [46.7, 4.8],
    "P_hs":    [2720, 212],
    "rhoR_cf": [1.60, 0.46],
    "CR_max":  [20.1, 0.0],
    "yield":   [256, 0.6],
    "gain":    [65, 0],
    "peak_velocity_kms":    [410, 0.0],
    "adiabat":              [6, 0.0],
    "ifar":                 [20, 0.0],
    "hydro_efficiency_pct": [8, 0.0],
    "imploded_DT_mass_mg":  [3, 0.0],
    "inflight_KE_kJ":      [300, 0.0],
    "fraction_absorbed_pct":[97, 0.0]
}
```

Each entry is `[value, uncertainty]`. Entries with `[0.0, 0.0]` are skipped.
Keys starting with `_` are treated as comments.

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
   Temperature mask fails for igniting capsules where alpha heating warms dense shell above 1 keV.
   - `unablated_fuel_mass` = fuel zones (0..fuel_bnd) inside ablation front / initial_fuel_mass
   - `unablated_ablatar_mass` = ablator zones (fuel_bnd..) inside ablation front / initial_ablator_mass
   - `stagnated_fuel_mass` = fuel zones between hs_bnd and ablation front / initial_fuel_mass
   - `initial_fuel_mass_mg`, `initial_ablator_mass_mg` stored on data and written to summary
   - Typo `unablated_ablatar_mass` (ablatar not ablator) is intentional -- do not fix without
     coordinated find-replace across all files.

7. **Burn propagation** (Olson et al. convention): Tracks hot-spot rhoR vs total rhoR over time.
   Ignition identified when hot-spot rhoR fraction exceeds 50%.

8. **EXODUS sampling principle**: EXODUS files contain only a fraction of the
   actual simulation timesteps. Therefore Helios's own time-integrated quantities
   (neutron count, deposited energy, etc.) are far more accurate than anything
   re-integrated from EXODUS data. The burn rate from EXODUS is used only as a
   *weighting function* for burn-averaging, never for absolute yield calculation.

9. **Adiabat**: alpha = P / P_Fermi where P_Fermi = 2.17 (rho/rho_0)^(5/3) Mbar,
   rho_0 = 0.205 g/cc (equimolar DT ice, Lindl convention). Evaluated at peak
   implosion velocity (pre-stagnation) in the cold unablated fuel only
   (ri[t, 0] to ri[t, 1]), excluding the ablated fuel region.

10. **Hydrodynamic efficiency**: eta_hydro = max(KE_inward) / E_absorbed.
    KE_inward = sum of 0.5 m v^2 for zones with v < 0 (imploding), maximized
    over all timesteps. E_absorbed = max(laser_energy_deposited) in Joules.

11. **IFAR (In-Flight Aspect Ratio)**: IFAR = R_shell / Delta_R at peak
    implosion velocity. Shell boundaries determined from density profile using
    rho > rho_peak / e threshold (NOT Lagrangian region interfaces, which
    include uncompressed vapor). Validated: VI_6 gives IFAR=18.1 vs published 20.

12. **Convergence ratios** -- two distinct quantities, both computed and reported:
    - `comp_ratio` = CR_stag = R0 / R_hs_at_stagnation (matches published data convention)
      R0 = initial inner shell radius = zbnd[0, ri[0,0]]
      R_hs = hot-spot boundary radius at stagnation = zbnd[stag_idx, ri[stag_idx,0]]
      Validated: Olson_PDD_9 gives 29.6 vs published 29.0 (2.1%)
    - `cr_inflight` = R0 / R_ablfront_at_peak_velocity (in-flight shell diagnostic)
      Reported in implosion section of summary; NOT used in published-data comparison.
    - DO NOT use density ratio for CR -- that was the old incorrect implementation.

13. **Hot-spot radius**: Use region interface `ri[stag_idx, 0]` (node index) to get
    hot-spot boundary radius from `zone_boundaries[stag_idx, hs_node]`.
    DO NOT use temperature mask (`T > threshold`) -- alpha heating in igniting capsules
    warms the dense shell above 1 keV, causing the mask to extend far outside the true
    hot spot, giving unphysically large radii (e.g. 0.68 cm > R0 for VI_6).

14. **Imploded DT mass**: Use zone-index method in burn_averaged_metrics.py --
    sum `zone_mass[stag_idx, :min(abl_idx+1, fuel_bnd)]` where `abl_idx` comes from
    `ablation_front_indices[stag_idx]`. DO NOT use radius-based method -- the smoothed
    ablation front radius is unreliable at stagnation and gave 1.55 mg vs correct 2.14 mg for VI_6.

15. **Peak velocity index**: Stored as `data.peak_velocity_index` (integer timestep index).
    Used for adiabat, IFAR, in-flight CR, and ablation front at peak velocity.

## Test Data

| Case | Regions | Zones | Key Features |
|------|---------|-------|-------------|
| **Olson_PDD_9** | 4 | 351 | 3-component pressure (ion + elec + rad), PDD, igniting |
| **VI_6** | 4 | 350 | Vulcan HDD target, 3-component pressure, work in progress |

### Olson_PDD_9 Region Structure
- Region 1: DT Vapor (zones 0-150) -- hot spot
- Region 2: DT Solid (zones 151-190) -- cryo DT ice
- Region 3: DT-CH foam (zones 191-320) -- wetted foam ablator
- Region 4: CH Skin (zones 321-350) -- outer ablator

### VI_6 Region Structure (Vulcan HDD)
- Region 1: DT gas (zones 0-50) -- hot spot
- Region 2: DT fuel (zones 51-100) -- cold DT fuel
- Region 3: DT ablated (zones 101-150) -- ablated DT
- Region 4: CD ablator (zones 151-350) -- outer ablator

### Validated Reference Values (Olson_PDD_9)
- Stagnation time: 12.599 ns (min HS radius = 0.0068 cm)
- Peak density: 181.17 g/cc (at t = 12.670 ns)
- Bang time: 12.681 ns
- Hot spot pressure: 106.97 Gbar
- Hot spot internal energy: 597.76 kJ
- Core radius: 0.0068 cm
- Stagnation CR: 29.6 (R0=0.2008 cm, R_hs=0.0068 cm)
- In-flight CR: 2.16 (R0=0.2008 cm, Rf=0.0930 cm at peak v)
- Fusion yield: 20.594 MJ, Target gain: 9.574
- <Ti>_n: 22.46 keV, <P>_n: 193.40 Gbar
- <rhoR>_n (fuel): 0.5189 g/cm2
- Burn-averaged: T=23.75 keV, P=208.4 Gbar, rhoR=0.5515 g/cm2
- Initial DT mass: 4.718 mg, Initial ablator mass: 0.349 mg
- Unablated fuel: 0.1262 (12.6%), Unablated ablator: 0.000 (all CH ablated)
- Stagnated fuel: 0.1219, Imploded DT: ~0.60 mg

### VI_6 Reference Values (work in progress -- laser deposition over-ablating target)
- Laser energy: 3.383 MJ (absorbed); published design is 4.0 MJ
- Stagnation time: 14.630 ns (min HS radius = 0.0046 cm)
- Peak density: 221.84 g/cc (at t = 14.700 ns)
- Bang time: 14.670 ns
- Hot spot pressure: 88.20 Gbar
- Stagnation CR: 41.1 (R0=0.1890 cm, R_hs=0.0046 cm)
- In-flight CR: 2.49 (R0=0.1890 cm, Rf=0.0759 cm at peak v)
- Peak implosion velocity: 763.1 km/s (published 410 km/s -- over-ablation)
- In-flight KE: 311.6 kJ, Hydro efficiency: 9.2%
- IFAR: 18.1 (density-based, rho > rho_peak/e)
- Adiabat: 1.13 (cold fuel at peak v_imp); published 6.0 -- different target variant
- Fusion yield: 28.811 MJ, Target gain: 8.516
- Burn-averaged: T=29.9 keV, P=340 Gbar, rhoR=0.88 g/cm2
- Initial DT mass: 3.709 mg, Initial ablator mass: 2.548 mg
- Unablated fuel: 0.5771 (57.7%), Unablated ablator: 0.000 (all CD ablated)
- Imploded DT: 2.14 mg (published 3.00 mg -- energy/ablation difference)

## Simulation Paths

### Mac Studio
- Olson: `~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9`
- VI_6:  `~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6/VI_6`

### MacBook
- Olson: `~/Sims/Helios_Sims/Xcimer_Sims/Olson_PDD/Olson_PDD_9/Olson_PDD_9`

## Workflow

1. Edit on MacBook (`python`) -> push to GitHub
2. Pull on Mac Studio (`python3`) -> run against simulation data
3. Terminal commands run one line at a time (avoid multi-line paste errors)
4. Validate against known reference values before moving to next feature
5. Output goes to `<sim>_summary.txt` -- use `grep` on the file, not stdout

## Quick Tests

```bash
# Mac Studio -- Olson (primary validation target)
python3 ~/helios_postprocessor/examples/run_analysis.py \
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9
grep -A 6 "MASS FRACTIONS" ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9_summary.txt

# Mac Studio -- VI_6
python3 ~/helios_postprocessor/examples/run_analysis.py \
  ~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6/VI_6
```

## Open Items

### Priority 1 -- Physics Module Integration (next session)
- `energetics`, `neutron_downscatter`, `pressure_gradients` work standalone
  but are not yet wired into `ICFAnalyzer` or `ICFPlotter`.
- Start by inspecting `def` signatures in each module, then wire into
  `ICFAnalyzer.integrate_physics_modules()` (new method), then add plotter pages
  and output sections.
- Cross-check: `energetics.py` KE_inward should match `burn_averaged_metrics.py`
  hydro efficiency (both give 9.2% for VI_6).

### Priority 2 -- VI_6 Laser Deposition Fix
- Over-ablation issue: simulated velocity 763 vs published 410 km/s.
- Low adiabat (1.13 vs published 6.0) and low remaining mass (2.14 vs 3.00 mg)
  are consistent with over-ablation, not just energy scaling.
- Published values (peak_velocity=410, adiabat=6.0) may be from a different
  target variant -- verify which published table VI_6_published.json references.

### Priority 3 -- CR_max Definition for VI_6
- Stagnation CR = 41.1 vs published 20.1 for VI_6.
- Published value likely uses shell mid-radius at stagnation, not hot-spot boundary.
- Olson matches well (29.6 vs 29.0) because hot-spot and shell nearly coincide
  at stagnation for high-convergence igniting targets.
- May need a separate `cr_shell_stag` metric using the ablation front radius
  at stagnation rather than the hot-spot boundary.

### Priority 4 -- `region_interfaces_indices` Robustness
- Currently relies on EXODUS region data; could add fallback based on
  density/composition gradients for targets without explicit region data.

## Dependencies

**Required**: numpy, scipy, matplotlib, netCDF4
**Optional**: scikit-learn (RANSAC shock fitting in icf_plotting.py -- guarded with `_HAS_SKLEARN`)
