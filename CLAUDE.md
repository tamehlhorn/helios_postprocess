# helios_postprocessor -- Project Guide

## Repository
- **GitHub**: `tamehlhorn/helios_postprocessor`
- **Package**: `helios_postprocess/`
- **Version**: 3.0.0 (March 2026)
- **Dev machine**: MacBook (`~/Codes/helios_postprocessor`) -- editing, pushing
- **Run machine**: Mac Studio (`tommehlhorn`, `~/helios_postprocessor`) -- running against simulation data

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
python3 ~/helios_postprocessor/examples/run_analysis.py <base_path>
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
- Imploded DT mass at stagnation
- IFAR (density-based shell boundaries)

Yield, laser energy, and gain come from Helios's own time-integrated
quantities -- NOT re-integrated from sampled EXODUS data. See Physics
Convention #8 below.

### Additional Physics Modules (functional style)

| Module | Purpose |
|--------|---------|
| `energetics` | Kinetic energy, hydro efficiency, PdV work |
| `neutron_downscatter` | Down-scatter ratio (DSR) diagnostics |
| `pressure_gradients` | Shock identification, RT instability assessment |

### Optional

`RHWParser` reads `.rhw` input files for drive configuration (direct/indirect,
burn on/off, drive temperature profile).

## Active Source Files

| File | Lines | Role |
|------|-------|------|
| **Pipeline** | | |
| `core.py` | ~965 | `HeliosRun` -- EXODUS/netCDF4 reader |
| `data_builder.py` | ~555 | Bridge: `HeliosRun` -> `ICFRunData` dataclass |
| `icf_analysis.py` | ~1300 | `ICFAnalyzer` -- drive, stagnation, burn, implosion, IFAR, mass fractions, burn propagation |
| `icf_plotting.py` | ~1580 | `ICFPlotter` -- full PDF report |
| `icf_output.py` | ~430 | `ICFOutputGenerator` -- summary text + CSV time histories |
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

4. **Region boundaries**: `region_interfaces_indices` are NODE indices (not zone indices).
   - Hot-spot boundary: `ri[:, 0]`
   - Cold fuel / ablated fuel boundary: `ri[:, 1]`
   - Fuel/ablator boundary: `ri[:, -2]`

5. **Fusion yield**: Preferred method uses `TimeIntFusionProd_n_1406 x 17.6 MeV x 1.602e-19 MJ/MeV`.
   Fallback: time-integrate `FusionRate_DT_nHe4` (less precise).

6. **Mass fractions**: Stagnated fuel uses zone-boundary definition --
   fuel zones between hot-spot boundary and ablation front at stagnation.
   Temperature-based cold mask was too aggressive for igniting targets.

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

## Test Data

| Case | Regions | Zones | Key Features |
|------|---------|-------|-------------|
| **Olson_PDD_9** | 4 | 351 | 3-component pressure (ion + elec + rad), PDD |
| **VI_6** | 4 | 350 | Vulcan HDD target, 3-component pressure |

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
- Core radius: 0.1853 cm
- Fusion yield: 20.594 MJ, Target gain: 9.574
- <Ti>_n: 22.46 keV, <P>_n: 193.40 Gbar
- <rhoR>_n (fuel): 0.5189 g/cm2
- Burn-averaged: T=23.75 keV, P=208.4 Gbar, rhoR=0.5515 g/cm2

### Validated Reference Values (VI_6)
- Laser energy: 3.383 MJ (absorbed)
- Stagnation time: 14.630 ns
- Peak density: 221.84 g/cc (at t = 14.700 ns)
- Bang time: 14.670 ns
- Peak implosion velocity: 763.1 km/s
- In-flight KE: 311.6 kJ, Hydro efficiency: 9.2%
- IFAR: 18.1 (density-based, rho > rho_peak/e)
- Adiabat: 1.13 (cold fuel at peak v_imp)
- Fusion yield: 28.811 MJ, Target gain: 8.516
- Burn-averaged: T=29.9 keV, P=340 Gbar, rhoR=0.88 g/cm2

## Workflow

1. Edit on MacBook -> push to GitHub
2. Pull on Mac Studio -> run against simulation data at `~/Sims/Xcimer/`
3. Terminal commands run one line at a time (avoid multi-line paste errors)
4. Validate against known reference values before moving to next feature

## Quick Test

```bash
python3 ~/helios_postprocessor/examples/run_analysis.py ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9
```

## Open Items

### Priority 1 -- Mass Fractions Refinement
- Unablated fuel, ablator, and stagnated fuel fractions compute but may need
  further validation against published values.
- Current approach: zone-boundary definition (no temperature mask).

### Priority 2 -- Physics Module Integration
- `energetics`, `neutron_downscatter`, `pressure_gradients` work standalone
  with `HeliosRun` data but are not yet wired into `ICFAnalyzer` or `ICFPlotter`.

### Priority 3 -- `region_interfaces_indices` Robustness
- Material boundary identification between hot spot, fuel, and ablator.
- Currently relies on EXODUS region data; could add fallback based on
  density/composition gradients.
- Note: for adiabat, cold fuel = ri[t,0] to ri[t,1] only.
  For IFAR, density-based boundaries (not Lagrangian) are required.

## Dependencies

**Required**: numpy, scipy, matplotlib, netCDF4
**Optional**: scikit-learn (RANSAC shock fitting in icf_plotting.py -- guarded with `_HAS_SKLEARN`)
