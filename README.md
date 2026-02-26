# Helios Postprocess - Integrated Package

Complete Python package for analyzing Helios ICF simulations with burn-averaged metrics for comparison with published target designs.

## 📦 Package Structure

```
helios_postprocess/
├── __init__.py                         # Main package interface
├── core.py                            # HeliosRun class & ExodusII reading
├── hot_spot.py                        # Hot spot identification & metrics
├── burn.py                            # Burn diagnostics (bang time, yield)
├── areal_density.py                   # Areal density (ρR) calculations
├── neutron_downscatter.py             # Neutron DSR diagnostics
├── pressure_gradients.py              # Shock & RT instability analysis
├── energetics.py                      # Kinetic energy & hydro efficiency
├── burn_averaged_metrics.py           # ⭐ NEW: Burn-averaged comparison
└── examples/
    └── compare_with_published_integrated.py  # Complete workflow example
```

## 🎯 What's New: Burn-Averaged Metrics

**burn_averaged_metrics.py** adds temporal burn-weighted averaging to complement your existing spatial neutron-averaging:

- **Existing modules**: Neutron-average T, P, ρ at each time (spatial weighting)
- **New module**: Burn-average over time (temporal weighting)

This enables direct comparison with published ICF performance metrics like your target table!

## 🚀 Quick Start

### Installation

```bash
# Option 1: Install as editable package (recommended for development)
cd /path/to/helios_postprocess
pip install -e .

# Option 2: Add to Python path
export PYTHONPATH=/path/to/helios_postprocess:$PYTHONPATH
```

### Basic Usage

```python
from helios_postprocess import HeliosRun

# Load simulation
run = HeliosRun('Vulcan_HDD_NB_111125_AI_2.exo')

# Get data at specific time
temperature = run.get_variable('temp', time_idx=-1)
density = run.get_variable('dens', time_idx=-1)

# Check what's available
print(run.list_variables())
```

### Burn-Averaged Comparison Workflow

```python
from helios_postprocess import (
    HeliosRun,
    extract_hot_spot_histories,
    calculate_burn_averaged_metrics,
    compare_with_published
)

# 1. Load simulation
run = HeliosRun('sim.exo')

# 2. Extract time histories
histories = extract_hot_spot_histories(run, T_threshold=1000.0)

# 3. Calculate burn-averaged metrics
metrics = calculate_burn_averaged_metrics(histories)

# 4. Compare with published data
published = {
    'T_hs': (46.7, 4.8),      # keV
    'P_hs': (2720, 212),      # Gbar
    'rhoR_cf': (1.60, 0.46),  # g/cm²
    'CR_max': (20.1, 9.6),
    'yield': (256, 0.6),      # MJ
    'gain': (65, 0)
}

comparison = compare_with_published(metrics, published, laser_energy_MJ=4.0)
print(comparison)
```

## 📊 Complete Example

Run the integrated comparison script:

```bash
cd helios_postprocess/examples
python compare_with_published_integrated.py path/to/your/file.exo
```

This will:
1. Load your Helios simulation
2. Extract hot spot time histories
3. Calculate burn-averaged metrics
4. Generate comparison table with published data
5. Check ignition criteria
6. Create diagnostic plots

## 🔬 Physics Modules Overview

### Spatial Analysis (at each time)

**hot_spot.py**
```python
from helios_postprocess import hot_spot

# Identify hot spot region
hot_mask, r_hs = hot_spot.identify_hot_spot(temp, zone_boundaries, T_threshold=1000)

# Calculate hot spot mass
m_hs = hot_spot.hot_spot_mass(density, zone_mass, temp, T_threshold=1000)

# Get mass-averaged temperature
T_avg = hot_spot.mass_averaged_temperature(temp, zone_mass, hot_mask)
```

**areal_density.py**
```python
from helios_postprocess import areal_density

# Calculate shell areal density
rhoR_shell = areal_density.calculate_shell_rhoR(density, zone_boundaries, T_threshold=1000)

# Line-of-sight ρR
rhoR_los = areal_density.calculate_line_rhoR(density, radii)
```

**neutron_downscatter.py**
```python
# Calculate neutron downscatter ratio
dsr = run.calculate_dsr(time_idx=-1)

# Convert to areal density
rhoR = run.dsr_to_areal_density(dsr['DSR'], calibration='NIF')
```

### Temporal Analysis (over time)

**burn.py**
```python
from helios_postprocess import burn as burn_module

# Calculate total fusion rate
rate = burn_module.calculate_total_fusion_rate(fusion_rate_density, zone_mass)

# Calculate neutron yield
yield_n = burn_module.calculate_neutron_yield(rate, time)

# Find bang time
bang_time = burn_module.find_bang_time(rate, time)
```

**burn_averaged_metrics.py** ⭐
```python
from helios_postprocess.burn_averaged_metrics import (
    extract_hot_spot_histories,
    calculate_burn_averaged_metrics
)

# Extract time histories
histories = extract_hot_spot_histories(run)

# Calculate burn-averaged quantities
metrics = calculate_burn_averaged_metrics(histories)

print(f"⟨T_hs⟩ = {metrics['T_burn_avg']:.1f} keV")
print(f"⟨P_hs⟩ = {metrics['P_burn_avg']:.0f} Gbar")
print(f"Yield = {metrics['yield_MJ']:.1f} MJ")
```

**energetics.py**
```python
from helios_postprocess import energetics

# Calculate kinetic energy
KE = energetics.calculate_kinetic_energy(mass, velocity)

# Hydrodynamic efficiency
eta = energetics.calculate_hydro_efficiency(KE, absorbed_energy)
```

## 📐 Key Metrics Explained

### Burn-Averaged Quantities

Burn-averaged values weight each time by fusion reaction rate:

```
<Q> = ∫ Q(t) · Ṙ(t) dt / ∫ Ṙ(t) dt
```

where Ṙ(t) ∝ ρ² · <σv>(T) is the fusion rate.

### Your Published Target Table

| Metric | Symbol | Published Value | Units |
|--------|--------|----------------|-------|
| Temperature | ⟨T_hs⟩ | 46.7 ± 4.8 | keV |
| Pressure | ⟨P_hs⟩ | 2720 ± 212 | Gbar |
| Areal Density | ⟨ρR_cf⟩ | 1.60 ± 0.46 | g/cm² |
| Convergence | CR_max | 20.1 ± 9.6 | - |
| Yield | Y | 256 ± 0.6 | MJ |
| Gain | G | 65 | - |

### Lindl Ignition Criteria

```python
from helios_postprocess import print_ignition_criteria
print_ignition_criteria()
```

- Temperature: > 5 keV
- Pressure: > 100 Gbar
- Areal Density: > 0.3 g/cm²
- Convergence Ratio: > 15

## ⚙️ Integration Notes

### How It Works Together

1. **HeliosRun** (core.py) reads ExodusII files
2. **hot_spot.py** identifies hot spot at each time → neutron-averages T, P, ρ
3. **burn_averaged_metrics.py** extracts time histories → burn-averages over time
4. **compare_with_published()** generates comparison table

### Data Flow

```
ExodusII file
    ↓
HeliosRun.get_variable()
    ↓
hot_spot.identify_hot_spot() at each time
    ↓
extract_hot_spot_histories() → time arrays
    ↓
calculate_burn_averaged_metrics() → <T>, <P>, <ρR>
    ↓
compare_with_published() → comparison table
```

## 🔧 Customization

### Adjust Hot Spot Threshold

```python
# Use 2 keV threshold instead of 1 keV
histories = extract_hot_spot_histories(run, T_threshold=2000.0)  # eV
```

### Select Specific Time Range

```python
# Analyze only times near peak burn
bang_idx = np.argmax(temperature_history)
time_window = slice(bang_idx-20, bang_idx+20)
histories = extract_hot_spot_histories(run, time_indices=time_window)
```

### Modify Comparison Metrics

```python
# Add your own metrics
published = {
    'T_hs': (46.7, 4.8),
    'P_hs': (2720, 212),
    'rhoR_cf': (1.60, 0.46),
    # ... add more as needed
}
```

## ⚠️ Important Notes

### Units
- **Temperature**: ExodusII stores in eV, convert to keV for comparison
- **Pressure**: Check if stored as J/cm³ or dyne/cm² in your files
- **Time**: ExodusII in seconds, convert to ns for display
- **Radius**: Convert cm → μm for display

### Hot Spot vs Cold Fuel
- **Hot spot**: High T, low ρ (where fusion occurs)
- **Cold fuel**: Low T, high ρ (surrounding shell)
- Use **cold fuel ρR** for comparison, not hot spot ρR!

### Coordinate Systems
The `extract_hot_spot_histories()` function assumes:
- Spherical symmetry for volume calculations
- Radial coordinates available as 'coordr' or can be calculated
- Adjust if using Cartesian coordinates

## 🐛 Troubleshooting

### "ModuleNotFoundError: No module named 'helios_postprocess'"
- Install package: `pip install -e /path/to/helios_postprocess`
- Or add to path: `export PYTHONPATH=/path/to/helios_postprocess:$PYTHONPATH`

### "Variable 'temp' not found"
- Check available variables: `run.list_variables()`
- Your file might use different names (e.g., 'tev', 'tion')

### Very different from published values
- Verify hot spot identification (T_threshold)
- Check that you're using cold fuel ρR
- Confirm time range includes full burn
- Validate unit conversions

### NaN or Inf in results
- Check for zero/negative temperatures
- Verify time array is monotonic
- Ensure densities are positive

## 📚 References

1. **Lindl et al.** (2004) - Ignition criteria  
   *Phys. Plasmas* 11, 339

2. **Bosch & Hale** (1992) - DT reactivity  
   *Nucl. Fusion* 32, 611

3. **Hurricane et al.** (2019) - Burning plasmas  
   *Phys. Plasmas* 26, 052704

## 👤 Author

Prof T  
Version: 2.0 (November 2025)

## 📝 License

Part of the helios_postprocess research package
