# CLAUDE.md - Helios Postprocessor Project Context

## Project Overview
Production-ready Python package for analyzing Helios ICF (Inertial Confinement Fusion) rad-hydro simulation outputs. Transforms raw EXODUS/netCDF simulation files into comprehensive performance analyses including PDF visualizations, JSON metrics, text summaries, and CSV time series data.

## Repository
- **GitHub:** https://github.com/tamehlhorn/helios_postprocess
- **Package name:** `helios_postprocess`
- **Repo directory:** `helios_postprocessor`
- **Version:** 2.0.0
- **Author:** Prof T (Xcimer, ICF research)

## Development Environment
- **MacBook Pro (M3):** Development/editing machine, repo at `~/Codes/helios_postprocessor/`
- **Mac Studio:** Data processing machine, repo at `~/helios_postprocessor/`, data at `~/Sims/Xcimer/`
- **Workflow:** Edit on MacBook → git push → git pull on Mac Studio → run against data
- **Editor:** VS Code (local for editing, Remote-SSH for Mac Studio execution)
- **Mac Studio Python:** `/usr/bin/python3` (3.9.6, system Python, no conda)
- **Package installed:** `python3 -m pip install -e .` (editable mode on Mac Studio)

## Package Structure
```
helios_postprocessor/           ← repo root
    setup.py
    README.md
    SETUP_GUIDE.md
    CLAUDE.md                   ← this file
    examples/
        compare_with_published_integrated.py
    helios_postprocess/         ← Python package
        __init__.py
        core.py                 ← HeliosRun class (main interface, ~1100 lines)
        burn.py                 ← Burn diagnostics
        hot_spot.py             ← Hot spot identification
        areal_density.py        ← Areal density calculations
        neutron_downscatter.py  ← DSR diagnostics
        pressure_gradients.py   ← Shock/RT analysis
        energetics.py           ← Kinetic energy/efficiency
        burn_averaged_metrics.py ← Burn-averaged quantities
```

## Key Class: HeliosRun
```python
from helios_postprocess import HeliosRun
run = HeliosRun('simulation.exo')

# Core attributes
run.times          # 1D array of simulation times (seconds)
run.n_times        # number of time steps
run.coords         # dict (currently empty - use zone_boundaries instead)

# Core methods
run.list_variables()                        # list all available variable names
run.get_variable('mass_density', time_idx=0) # get zone data at time index
run.get_burn_diagnostics()                  # bang time, burn width, yield
run.calculate_areal_density(time_idx=-1)    # returns dict with 'rhoR'
run.check_ignition_criteria()               # Lindl criteria check
run.calculate_dsr(time_idx=-1)              # downscatter ratio
run.close()
```

## Helios EXODUS Variable Names
The EXODUS files use Helios-specific variable names, NOT generic ones. Common mappings:
- **Density:** `mass_density` (NOT `dens` or `density`)
- **Ion temperature:** `ion_temperature` (NOT `temp` or `temperature`)
- **Electron temperature:** `elec_temperature`
- **Radiation temperature:** `radiation_temperature`
- **Ion pressure:** `ion_pressure` (units: J/cm³; divide by 1e8 for Gbar)
- **Electron pressure:** `elec_pressure`
- **Velocity:** `fluid_velocity`
- **Zone boundaries:** `zone_boundaries` (N+1 nodes for N zones)
- **Zone mass:** `zone_mass`
- **DT fusion rate:** `FusionRate_DT_nHe4` (zone-level, sum spatially for global rate)
- **DD fusion rate:** `FusionRate_DD_nHe3`, `FusionRate_DD_pT`
- **DT neutron count:** `TimeIntFusionProd_n_1406` (14.06 MeV, time-integrated)
- **DD neutron count:** `TimeIntFusProd_n_0245` (2.45 MeV, time-integrated)
- **Laser energy deposited:** `EnLaserDepositedTimeIntg`
- **Laser power:** `LaserPwrSrc`, `LaserEnDelivered`

## Physics Conventions
- **Pressure:** Report in Gbar. 1 Gbar = 1e8 J/cm³
- **Total pressure:** ion_pressure + elec_pressure (NOT ion + rad_pressure for reporting)
- **Fusion yield:** Calculate from `TimeIntFusionProd_n_1406` × 17.6 MeV × 1.602e-19 MJ/MeV
  - Do NOT integrate fusion rate in time (imprecise due to time grid)
- **Laser energy:** Helios stores in Joules, convert to MJ (× 1e-6)
- **Target gain:** yield_MJ / laser_energy_MJ
- **DT fusion energy:** 17.6 MeV per reaction
- **Variable centering:** Temperature and density are zone-centered; velocity is node-centered (N+1 points)
- **Indexing:** Helios uses Fortran 1-indexed conventions; Python arrays are 0-indexed

## Lindl Ignition Criteria
- Temperature: > 5 keV
- Pressure: > 100 Gbar
- Areal density (ρR): > 0.3 g/cm²
- Convergence ratio: > 15

## Known Issues / Recent Fixes
1. **FIXED:** `get_burn_diagnostics()` now searches for `FusionRate_DT_nHe4` first, builds global rate vs time by summing over zones. Falls back to T² proxy only if no fusion rate variable exists.
2. **FIXED:** `calculate_areal_density()` uses `mass_density` and `zone_boundaries` directly. Removed broken import of nonexistent `.physics` subpackage.
3. **WARNING:** `get_burn_diagnostics_corrected()` standalone function (line ~190) - receives neutron_rate as parameter, does NOT use `self`. The class method at line ~995 finds the right data and passes it in.
4. **ISSUE:** `run.coords` returns empty dict `{}`. Use `run.get_variable('zone_boundaries', time_idx=i)` instead.
5. **ISSUE:** Some methods in `core.py` still reference `'dens'` or `'temp'` instead of `'mass_density'` / `'ion_temperature'`. These need auditing.
6. **ISSUE:** FWHM calculation sometimes fails, falls back to RMS width. Not critical but should be investigated.

## What's Missing / Next Steps
### Priority 1: Visualization Module
- No plotting module in current repo
- Best candidate: `~/Sims/Helios_Sims/helios_processor/icf_plotting.py` (1295 lines, most complete)
- Also available: `~/Sims/Helios_Sims/Xcimer_Sims/Helios postprocess 3/icf_plotting.py` (1023 lines)
- These depend on `ICFRunData` dataclass (old interface), need adaptation to `HeliosRun`
- **Recommended approach:** Create `data_builder.py` that loads all time steps from HeliosRun into 2D arrays and scalar metrics, then feed to existing plotter

### Priority 2: Output Generation Module
- `icf_output.py` exists in old versions at:
  - `~/Sims/Helios_Sims/helios_processor/icf_output.py`
  - `~/Sims/Helios_Sims/Xcimer_Sims/Helios postprocess 3/icf_output.py`
- Generates JSON, CSV, TXT, NPZ outputs

### Priority 3: Analysis Module
- `icf_analysis.py` exists in old versions
- Computes all derived metrics (30+): timing, compression, hot spot, energy, etc.

### Priority 4: RHW Config Parser
- Parser exists: `~/Sims/Helios_Sims/helios_processor/parsers/rhw_parser.py`
- Parses Helios .rhw input files for drive type, burn status, drive temperature profile

### Priority 5: Consolidation
- Multiple old copies exist on MacBook (see below), should be archived after integration
- Variable name audit throughout all modules

## Old Code Locations (MacBook)
These contain code to mine for the integration:
```
~/Sims/Helios_Sims/helios_processor/          ← icf_* modules + RHW parser
~/Sims/Helios_Sims/Xcimer_Sims/Helios postprocess 3/  ← icf_* modules
~/Sims/Helios_Sims/Postprocess archives/helios_postprocess/     ← archive copy
~/Sims/Helios_Sims/Postprocess archives/helios_postprocess_dir/ ← older with physics/ subdir
```

## Validated Test Case
**Olson PDD 8/9** (`~/Sims/Xcimer/Olson_PDD/Olson_PDD_8/Olson_PDD_8.exo` on Mac Studio):
- 1095 time steps, 0 to 16.018 ns
- 350 zones, 351 boundaries
- Peak density: 112 g/cc at stagnation (12.5 ns)
- Bang time: 12.68 ns
- Peak Ti: 39 keV, Peak Te: 25 keV
- Total pressure: 266 Gbar
- ρR: 0.64 g/cm² at bang time
- DT yield: 20.6 MJ from 2.15 MJ laser → gain 9.57
- 248 nm KrF laser, direct drive

## For Starting a New Chat
Upload these files:
1. **CLAUDE.md** (this file) - complete project context
2. **helios_postprocess/core.py** - current main module
3. **icf_plotting.py** (from `~/Sims/Helios_Sims/helios_processor/`) - plotting to integrate
4. **icf_analysis.py** (from same directory) - analysis to integrate
5. **icf_output.py** (from same directory) - output generation to integrate
