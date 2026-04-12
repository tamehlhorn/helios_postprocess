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
HeliosRun(exo_path)  ->  build_run_data(run)  ->  ICFRunData (dataclass, ~80 attributes)
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
| `data_builder.py` | ~580 | Bridge: `HeliosRun` -> `ICFRunData` dataclass |
| `icf_analysis.py` | ~1340 | `ICFAnalyzer` -- drive, stagnation, burn, implosion, IFAR, mass fractions, burn propagation, convergence ratios |
| `icf_plotting.py` | ~1580 | `ICFPlotter` -- full PDF report |
| `icf_output.py` | ~470 | `ICFOutputGenerator` -- summary text + CSV time histories |
| `burn_averaged_metrics.py` | ~560 | Temporal burn-averaging, implosion metrics, published-data comparison |
| `energetics.py` | -- | Kinetic energy, hydro efficiency, PdV work |
| `neutron_downscatter.py` | -- | Neutron down-scatter ratio diagnostics |
| `pressure_gradients.py` | -- | Pressure gradient analysis, shock ID, RT assessment |
| `rhw_parser.py` | ~300 | `RHWParser` -- reads `.rhw` input files including EOS and burn model |
| `__init__.py` | ~90 | Package exports |
| `examples/run_analysis.py` | ~350 | CLI runner with auto file derivation and comparison PDF |

## Published Data Comparison (JSON format)

Place a `<name>_published.json` file next to the `.exo` file. Format:

```json
{
    "_note": "Olson LILAC/Xrage/Hydra reference at 1.4x drive multiplier, ~2.15 MJ",
    "laser_energy_MJ": 2.15,
    "T_hs":    [22.5, 2.0],
    "P_hs":    [193, 20],
    "rhoR_cf": [1.10, 0.05],
    "CR_max":  [29.0, 3.0],
    "yield":   [87.4, 0.0],
    "gain":    [40.6, 0.0],
    "peak_velocity_kms":    [470, 0.0],
    "adiabat":              [3.0, 0.5],
    "ifar":                 [0.0, 0.0],
    "hydro_efficiency_pct": [0.0, 0.0],
    "imploded_DT_mass_mg":  [0.0, 0.0],
    "inflight_KE_kJ":       [0.0, 0.0],
    "fraction_absorbed_pct":[65.1, 9.3],
    "P_hs_ignition_Gbar":   [90.0, 0.0],
    "hs_radius_ignition_um":[120.0, 0.0],
    "_stagnation_time_ns":      [13.47, 0.0],
    "_hs_internal_energy_kJ":   [70.0, 0.0],
    "_ignition_time_ns":        [13.3, 0.0],
    "_peak_rhoR_cf_time_ns":    [13.36, 0.0],
    "_hs_rhoR_at_13p24ns":      [0.942, 0.0]
}
```

Keys starting with `_` are comments (skipped). Each entry is `[value, uncertainty]`.
`[0.0, 0.0]` entries are skipped in comparison.

**CRITICAL**: peak_velocity reference is **470 km/s** (not 410 km/s -- that was an error).
**CRITICAL**: fraction_absorbed is **65.1 ± 9.3%** (not 97% -- that was symmetric DD).

## Unit Conventions

| Quantity | Internal Unit | Conversion from Helios |
|----------|--------------|----------------------|
| Temperature | eV (arrays), keV (scalars/reporting) | Helios stores eV; /1000 for keV |
| Pressure | J/cm3 (arrays), Gbar (reporting) | x1e-8 from J/cm3 |
| Velocity | cm/s (arrays), km/s (reporting) | x1e-5 from cm/s |
| Areal density | g/cm2 | rhoR = integral(rho dr) |
| Laser energy | MJ | x1e-6 from J |
| Time | ns | Helios may store in seconds -> x1e9 |

## Physics Conventions

1. **Pressure aliasing**: Total pressure = `ion_pressure + rad_pressure` everywhere.
2. **Stagnation** = minimum hot-spot outer radius (NOT peak density).
3. **Ablation front** = steepest negative drho/dr outside hot spot. Use zone indices.
4. **Region boundaries**: `region_interfaces_indices` are NODE indices.
   - Hot-spot boundary: `ri[:, 0]`; Cold fuel boundary: `ri[:, 1]`; Fuel/ablator: `ri[:, -2]`
5. **Fusion yield**: `TimeIntFusionProd_n_1406 x 17.6 MeV x 1.602e-19 MJ/MeV`.
6. **Mass fractions**: Zone-index method. Typo `unablated_ablatar_mass` is intentional.
7. **Burn propagation**: Ignition when hot-spot rhoR >= 0.3 g/cm2 (Olson et al.).
8. **EXODUS sampling**: Helios time-integrated quantities authoritative. EXODUS burn rate
   used only as weighting function.
9. **Adiabat**: P/P_Fermi, P_Fermi=2.17(rho/0.205)^(5/3) Mbar. At peak velocity,
   cold unablated fuel only (ri[t,0] to ri[t,1]).
10. **Hydrodynamic efficiency**: eta = max(KE_inward) / E_absorbed.
11. **IFAR**: rho > rho_peak/e threshold (NOT Lagrangian region interfaces).
12. **Convergence ratios**:
    - `comp_ratio` = CR_stag = R0/R_hs_at_stagnation
    - `cr_inflight` = R0/R_ablfront_at_peak_velocity (known bug: currently uses HS boundary)
13. **Hot-spot radius**: Use region interface ri[stag_idx,0]. NOT temperature mask.
14. **Peak velocity index**: `data.peak_velocity_index`. Shell zones only, strictly pre-bang.

## Burn Model Classification (RHW flags)

| Fusion transport on | Use alpha deposition | Label |
|--------------------|---------------------|-------|
| 0 | any | Disabled (no alpha heating) |
| 1 | 1 | Local only (instantaneous — use with caution) |
| 1 | 0 | Non-local transport (time/space dependent) |

**WARNING**: Local alpha model is non-physical. Causes adiabat >70, velocity spikes.
NEVER use local alpha for calibration runs. Use non-local or disabled.

## EOS Model Reporting

Per-region EOS parsed by `_parse_eos_models()`. Shown in summary EOS MODELS section.
- `EOS data type = 1` → SESAME; `EOS data type = 0` → PROPACEOS
- EOS choice has NO effect on adiabat -- confirmed by controlled test April 2026.
- Current TM runs: SESAME DT Vapor+Solid, PROPACEOS DT-CH foam, SESAME or PROPACEOS CH Skin.

## PDD Calibration Status (April 2026)

### Reference values (Olson et al. Phys. Plasmas 28, 122704 (2021))
- Peak velocity: **470 km/s**
- Absorbed energy: **1.2-1.6 MJ** = 65.1 ± 9.3% of 2.15 MJ delivered
- Adiabat: **~3.0 ± 0.5**
- CR: **29.0 ± 3.0**
- Stagnation time: **~13.47 ns**
- ρR_cf: **1.10 ± 0.05 g/cm²**
- Yield: **87.4 MJ**, Gain: **40.6**
- Ignition time: ~13.3 ns; HS pressure at ignition: 90 Gbar; HS radius: 120 μm
- HS rhoR at 13.24 ns: 0.942 g/cm²; HS internal energy: ~70 kJ

### Best geometry: spot=0.25 cm, cone=40°, d=0.21 cm, Gaussian, f=0.06, peak=329 TW

### Calibration scan results

| Run | Spot | d | Foot | f | Burn | v (km/s) | Abs% | Adiabat | CR | ρR_cf | Notes |
|-----|------|---|------|---|------|----------|------|---------|-----|-------|-------|
| TM_1nb | 0.30 | 0.20 | 25 | 0.06 | None | 379 | 61.5 | 2.01 | 30.8 | 0.46 | |
| TM_2nb | 0.25 | 0.23 | 25 | 0.06 | None | 402 | 68.0 | 2.04 | 46.1 | 0.70 | |
| TM_3nb | 0.27 | 0.21 | 25 | 0.06 | None | 394 | 65.8 | 1.94 | 38.7 | 0.58 | |
| TM_4nb | 0.25 | 0.21 | 25 | 0.06 | None | 403 | 68.6 | 1.98 | 45.1 | 0.75 | Best baseline |
| TM_4f8nb | 0.25 | 0.21 | 25 | 0.08 | None | 404 | 68.8 | 1.95 | 33.8 | 0.80 | |
| TM_4f8alpha | 0.25 | 0.21 | 25 | 0.08 | Non-local | 404 | 68.8 | 1.96 | **29.1** | 0.72 | CR=29.1 ✓ |
| TM_4f8local | 0.25 | 0.21 | 25 | 0.08 | Local | 644 | 68.8 | 75 | 50.7 | 0.71 | DISCARD |
| TM_5f6nb | 0.25 | 0.21 | 26 | 0.06 | None | ~408 | ~68 | ~1.85 | ~34 | ~0.82 | foot +1 |
| TM_6f6nb | 0.25 | 0.21 | 27 | 0.06 | None | ~406 | ~68 | ~1.79 | ~34 | ~0.82 | foot +2 |
| TM_7f6nb | 0.25 | 0.21 | 28 | 0.06 | None | ~410 | ~68 | ~1.69 | ~34 | ~0.83 | foot +3 |
| TM_8f6nb | 0.25 | 0.21 | 26.5 | 0.06 | None | ~408 | ~68 | ~1.85 | ~34 | ~0.82 | foot 26.5 |
| TM_9f6nb | 0.25 | 0.21 | 23 | 0.06 | None | ~405 | ~69 | ~2.20 | ~33 | ~0.77 | foot -2 |
| **Target** | — | — | — | — | — | **470** | **65±9** | **3.0** | **29** | **1.10** | |

### Key calibration findings

1. **Absorption solved**: 65-69% matches published 65.1 ± 9.3% ✓
2. **CR matches with non-local alpha**: TM_4f8alpha CR=29.1 vs published 29.0 ✓
3. **Imploded DT matches**: 0.60 mg consistent across valid runs ✓
4. **Adiabat stubbornly ~1.7-2.2** regardless of geometry, EOS, or flux limiter
5. **Counter-intuitive foot trend**: higher foot power → LOWER adiabat at f=0.06
   (23 TW → adiabat 2.20; 28 TW → adiabat 1.69) -- under investigation
6. **Velocity gap**: ~404-410 vs 470 km/s (-13%) -- primary remaining target
7. **Stagnation ~1.1 ns late**: 14.5-14.7 ns vs 13.47 ns
8. **Local alpha model catastrophic**: NEVER use for calibration
9. **EOS has no effect on adiabat**: confirmed by controlled test
10. **Flux limiter cliff**: f=0.06 safe, f=0.08 borderline (non-local alpha ok,
    local alpha catastrophic)

### Active calibration direction

After understanding foot/adiabat relationship, planned geometry sweep:
- cone=35°, d=0.22 cm (TM_10f6nb) -- reduce late-time geometric miss
- cone=35°, d=0.23 cm (TM_11f6nb)
Then switch on non-local alpha transport for ignition check.

## Simulation Paths

### Mac Studio
- Olson PDD: `~/Sims/Xcimer/Olson_PDD/PDD_<name>/PDD_<name>`
- VI_6: `~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6/VI_6`

### MacBook
- Olson PDD_9: `~/Sims/Helios_Sims/Xcimer_Sims/Olson_PDD/Olson_PDD_9/Olson_PDD_9`

## Workflow

1. Edit on MacBook (`python`) -> push to GitHub
2. Pull on Mac Studio (`python3`) -> run against simulation data
3. Terminal commands run **one line at a time** (no multi-line paste)
4. Output goes to `<sim>_summary.txt` -- use `grep`, not stdout
5. Mac Studio direct edits: push from Studio, pull on MacBook to resync
6. `git config --global pull.rebase false` on MacBook (prevents divergence errors)
7. When MacBook has uncommitted changes blocking pull: `git stash; git pull; git stash drop`

## Quick Validation

```bash
# Run and check key metrics
python3 ~/helios_postprocessor/examples/run_analysis.py \
  ~/Sims/Xcimer/Olson_PDD/PDD_TM_4f8nb/PDD_TM_4f8nb
grep "Alpha burn model\|Adiabat\|Peak implosion\|Fraction absorbed\|CR_max" \
  ~/Sims/Xcimer/Olson_PDD/PDD_TM_4f8nb/PDD_TM_4f8nb_summary.txt
grep -A 8 "EOS MODELS" ~/Sims/Xcimer/Olson_PDD/PDD_TM_4f8nb/PDD_TM_4f8nb_summary.txt
grep -A 25 "COMPARISON WITH PUBLISHED" \
  ~/Sims/Xcimer/Olson_PDD/PDD_TM_4f8nb/PDD_TM_4f8nb_summary.txt
```

## Open Items

### Priority 1 -- Foot power / adiabat mystery (active investigation)
- Increasing foot 23→28 TW decreases adiabat 2.20→1.69 -- counter-intuitive
- Hypothesis: foot timing shifts which timestep is identified as peak velocity,
  changing the zone where adiabat is evaluated
- Check: plot adiabat vs time for TM_9f6nb (foot=23) vs TM_7f6nb (foot=28)
- Also check: is the mass-averaged adiabat region consistent across runs?

### Priority 2 -- Velocity gap (-13%)
- 404-410 vs 470 km/s target
- Next lever: geometry sweep, cone=35°, d=0.22-0.23 cm
- Do NOT increase peak power -- reference codes use same 329 TW peak

### Priority 3 -- cr_inflight bug
- Currently uses HS boundary radius at peak velocity (~6.6-6.8)
- Should use ablation front radius (~2.16 for PDD_9)
- Fix: `Rf = data.ablation_front_radius[pv_idx]` in compute_performance_metrics()

### Priority 4 -- Hot-spot radius inconsistency in _compute_hot_spot_properties()
- Returns ~0.08-0.22 cm at stagnation for TM runs (should be ~0.01 cm)
- Hot-spot pressure and internal energy at stagnation may be wrong
- Stagnation finder correctly reports 0.004-0.007 cm -- inconsistency needs fixing

### Priority 5 -- stagnation_time_ns as active comparison metric
- Currently `_stagnation_time_ns` in JSON (commented out, not active)
- Add as comparison key reading `data.stag_time`

### Priority 6 -- Physics Module Integration
- `energetics`, `neutron_downscatter`, `pressure_gradients` not yet wired into ICFAnalyzer

## New ICFRunData Fields (added April 2026)

```python
# Laser (all from RHW, default to 0 if RHW missing)
self.laser_wavelength_um: float = 0.0
self.laser_spot_size_cm: float = 0.0
self.laser_half_cone_angle_deg: float = 0.0
self.laser_focus_position_cm: float = 0.0
self.laser_power_multiplier: float = 1.0
self.laser_spatial_profile: str = "Uniform"
self.laser_foot_power_TW: float = 0.0
self.laser_peak_power_TW: float = 0.0
self.laser_foot_start_ns, laser_foot_end_ns: float = 0.0
self.laser_peak_start_ns, laser_peak_end_ns: float = 0.0
self.laser_pulse_duration_ns: float = 0.0
# EOS and burn model
self.eos_models: list = None          # [{'region', 'type', 'file'}]
self.alpha_deposition_local: bool = False
self.alpha_deposition_nonlocal: bool = False
```

## Summary Output Sections (added April 2026)

- **Alpha burn model** line in SIMULATION OVERVIEW
- **EOS MODELS** section after LASER CONFIGURATION
- Both degrade gracefully when RHW is absent

## Dependencies

**Required**: numpy, scipy, matplotlib, netCDF4
**Optional**: scikit-learn (RANSAC shock fitting -- guarded with `_HAS_SKLEARN`)
