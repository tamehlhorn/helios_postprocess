# helios_postprocessor — Project Guide

## Repository
- **GitHub**: `tamehlhorn/helios_postprocessor`
- **Package**: `helios_postprocess/`
- **Dev machine**: MacBook (`~/Codes/helios_postprocessor`) — editing, pushing
- **Run machine**: Mac Studio (`tommehlhorn`, `~/helios_postprocessor`) — running against simulation data

## Architecture

### Class-Based Pipeline (primary workflow)

```
HeliosRun(exo_path)  →  build_run_data(run)  →  ICFRunData (dataclass, ~68 attributes)
                                                    ↓
                                             ICFAnalyzer(data)
                                               .analyze_drive_phase()
                                               .analyze_stagnation_phase()   ← also calls analyze_implosion_phase()
                                               .analyze_burn_phase()
                                               .compute_performance_metrics()
                                                    ↓
                                             ICFPlotter(data, config)
                                               .create_full_report(pdf_path)
                                                    ↓
                                             ICFOutputGenerator(data)
                                               .write_all(base_path)   → _summary.txt + _history.csv
```

### Burn-Averaged Metrics (published-data comparison)

After the pipeline runs, burn-averaged metrics can be computed for comparison
with published ICF target designs:

```python
from helios_postprocess.burn_averaged_metrics import (
    extract_histories_from_run_data,
    calculate_burn_averaged_metrics,
    compare_with_published,
)

histories = extract_histories_from_run_data(data)   # pulls from ICFRunData
metrics = calculate_burn_averaged_metrics(histories)
print(compare_with_published(metrics, published_data, laser_energy_MJ=4.0))
```

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
| `core.py` | ~965 | `HeliosRun` — EXODUS/netCDF4 reader |
| `data_builder.py` | ~555 | Bridge: `HeliosRun` → `ICFRunData` dataclass |
| `icf_analysis.py` | ~1100 | `ICFAnalyzer` — drive, stagnation, burn, implosion, mass fractions, burn propagation |
| `icf_plotting.py` | ~1580 | `ICFPlotter` — full PDF report (contours, histories, trajectories, burn propagation, radial lineouts) |
| `icf_output.py` | ~430 | `ICFOutputGenerator` — summary text + CSV time histories |
| **Physics modules** | | |
| `burn_averaged_metrics.py` | ~310 | Temporal burn-averaging, published-data comparison (refactored for ICFRunData) |
| `energetics.py` | — | Kinetic energy, hydro efficiency, PdV work |
| `neutron_downscatter.py` | — | Neutron down-scatter ratio diagnostics |
| `pressure_gradients.py` | — | Pressure gradient analysis, shock ID, RT assessment |
| **Support** | | |
| `rhw_parser.py` | ~240 | `RHWParser` — reads `.rhw` input files |
| `__init__.py` | ~90 | Package exports |

## Archived Files

`archive/legacy_modules/` contains 3 standalone modules from v2.0 whose physics
is fully incorporated into `icf_analysis.py`:

- `areal_density.py` — superseded by `ICFAnalyzer._compute_areal_densities()`
- `burn.py` — superseded by `ICFAnalyzer.analyze_burn_phase()`
- `hot_spot.py` — superseded by `ICFAnalyzer._compute_hot_spot_properties()`

Note: `burn_averaged_metrics.py` formerly imported these three modules. It was
refactored in v3.0 to read from `ICFRunData` instead, eliminating those dependencies.

## Unit Conventions

| Quantity | Internal Unit | Conversion from Helios |
|----------|--------------|----------------------|
| Temperature | eV (arrays), keV (scalars/reporting) | Helios stores eV; ÷1000 for keV |
| Pressure | J/cm³ (arrays), Gbar (reporting) | ×1e-8 from J/cm³ |
| Velocity | cm/s (arrays), km/s (reporting) | ×1e-5 from cm/s |
| Areal density | g/cm² | ρR = ∫ρ dr |
| Laser energy | MJ | ×1e-6 from J |
| Drive temperature | eV | Stays in eV (radiation drive) |
| Time | ns | Helios may store in seconds → ×1e9 |

## Physics Conventions

1. **Pressure aliasing**: `rad_pressure = elec_pressure` (or `elec + actual_rad` when 3-component).
   Total pressure = `ion_pressure + rad_pressure` everywhere. Documented in data_builder.py.

2. **Stagnation** = minimum hot-spot outer radius (NOT peak density).
   For igniting capsules, peak ρ occurs after alpha heating.

3. **Ablation front** = steepest negative dρ/dr outside the hot spot boundary.
   Smoothed with iterative algorithm anchored at stagnation.

4. **Region boundaries**: `region_interfaces_indices` are NODE indices (not zone indices).
   - Hot-spot boundary: `ri[:, 0]`
   - Fuel/ablator boundary: `ri[:, -2]`

5. **Fusion yield**: Preferred method uses `TimeIntFusionProd_n_1406 × 17.6 MeV × 1.602e-19 MJ/MeV`.
   Fallback: time-integrate `FusionRate_DT_nHe4` (less precise).

6. **Mass fractions**: Stagnated fuel uses zone-boundary definition —
   fuel zones between hot-spot boundary and ablation front at stagnation.
   Temperature-based cold mask was too aggressive for igniting targets.

7. **Burn propagation** (Olson et al. convention): Tracks hot-spot ρR vs total ρR over time.
   Ignition identified when hot-spot ρR fraction exceeds 50%.

## Test Data

| Case | Regions | Zones | Key Features |
|------|---------|-------|-------------|
| **Olson_PDD_8** | 2 (DT gas + DT ice) | 350 | 2-component pressure (ion + elec) |
| **Olson_PDD_9** | 4 | 351 | 3-component pressure (ion + elec + rad) |

### Olson_PDD_9 Region Structure
- Region 1: DT Vapor (zones 0–150) — hot spot
- Region 2: DT Solid (zones 151–190) — cryo DT ice
- Region 3: DT-CH foam (zones 191–320) — wetted foam ablator
- Region 4: CH Skin (zones 321–350) — outer ablator

### Validated Reference Values (Olson_PDD_8)
- 1095 timesteps × 350 zones
- Peak density: ~112 g/cc
- Bang time: ~12.68 ns
- Laser energy: ~2.15 MJ
- Fusion yield: ~20.6 MJ
- Target gain: ~9.57

### Validated Reference Values (Olson_PDD_9)
- Stagnation time: 12.599 ns (min HS radius = 0.0068 cm)
- Peak density: 181.17 g/cc (at t = 12.670 ns)
- Bang time: 12.681 ns
- Hot spot pressure: 106.97 Gbar
- Hot spot internal energy: 597.76 kJ
- Core radius: 0.1853 cm
- Fusion yield: 20.594 MJ, Target gain: 9.574
- ⟨Ti⟩_n: 22.46 keV, ⟨P⟩_n: 193.40 Gbar
- ⟨ρR⟩_n (fuel): 0.5189 g/cm²

## Workflow

1. Edit on MacBook → push to GitHub
2. Pull on Mac Studio → run against simulation data at `~/Sims/Xcimer/`
3. Terminal commands run one line at a time (avoid multi-line paste errors)
4. Validate against known reference values before moving to next feature

## Test Script

```python
from helios_postprocess import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.icf_analysis import ICFAnalyzer
from helios_postprocess.icf_plotting import ICFPlotter
from helios_postprocess.icf_output import ICFOutputGenerator
import logging
logging.basicConfig(level=logging.INFO)

run = HeliosRun('/Users/tommehlhorn/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9.exo', verbose=True)
data = build_run_data(run, time_unit='s')
run.close()

analyzer = ICFAnalyzer(data)
analyzer.analyze_drive_phase()
analyzer.analyze_stagnation_phase()
analyzer.analyze_burn_phase()
analyzer.compute_performance_metrics()

plotter = ICFPlotter(data, {})
plotter.create_full_report('/Users/tommehlhorn/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9_report.pdf')

output = ICFOutputGenerator(data)
output.write_all('/Users/tommehlhorn/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9')
print('Done — report + summary + history written')
```

## Open Items

### Priority 1 — Mass Fractions Refinement
- Unablated fuel, ablator, and stagnated fuel fractions compute but may need
  further validation against published values.
- Current approach: zone-boundary definition (no temperature mask).

### Priority 2 — Burn-Averaged Histories Validation
- `extract_histories_from_run_data()` uses a flexible attribute lookup to
  pull hot-spot time histories from ICFRunData. The attribute names it looks
  for (e.g. `hot_spot_ion_temperature`, `fuel_areal_density`) need to be
  verified against the actual ICFRunData attributes populated by ICFAnalyzer.
- May need to add per-timestep hot-spot history arrays to ICFAnalyzer if
  they don't exist yet (currently many quantities are stored as scalars at
  stagnation/bang time rather than full time histories).

### Priority 3 — Physics Module Integration
- `energetics`, `neutron_downscatter`, `pressure_gradients` work standalone
  with `HeliosRun` data but are not yet wired into `ICFAnalyzer` or `ICFPlotter`.

### Priority 4 — `region_interfaces_indices` Robustness
- Material boundary identification between hot spot, fuel, and ablator.
- Currently relies on EXODUS region data; could add fallback based on
  density/composition gradients.

## Dependencies

**Required**: numpy, scipy, matplotlib, netCDF4
**Optional**: scikit-learn (RANSAC shock fitting in icf_plotting.py — guarded with `_HAS_SKLEARN`)
