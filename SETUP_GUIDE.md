# Setup Guide for Integrated Helios Postprocess Package

## What You Have

A complete, integrated Python package that combines:
- Your existing 7 physics modules (core, hot_spot, burn, areal_density, neutron_downscatter, pressure_gradients, energetics)
- New burn_averaged_metrics module for comparing with published data
- Example scripts and comprehensive documentation

## Quick Setup (3 Steps)

### Step 1: Download and Extract

Download the complete `helios_postprocess` folder to your Mac.

Your directory structure should look like:
```
helios_postprocess/
├── __init__.py
├── setup.py
├── README.md
├── core.py
├── hot_spot.py
├── burn.py
├── areal_density.py
├── neutron_downscatter.py
├── pressure_gradients.py
├── energetics.py
├── burn_averaged_metrics.py
└── examples/
    └── compare_with_published_integrated.py
```

### Step 2: Install the Package

Open Terminal on your Mac and navigate to the package directory:

```bash
cd /path/to/helios_postprocess
```

Install in editable mode (recommended for development):

```bash
pip install -e .
```

This will:
- Install required dependencies (numpy, scipy, matplotlib, netCDF4)
- Make `helios_postprocess` available system-wide
- Allow you to edit files and see changes immediately

### Step 3: Test the Installation

```bash
python -c "from helios_postprocess import HeliosRun; print('✓ Installation successful!')"
```

## Usage

### Basic Test

```python
from helios_postprocess import HeliosRun

# Load your simulation
run = HeliosRun('path/to/your/Vulcan_file.exo')

# Check what's available
print(run.list_variables())
print(f"Time steps: {run.n_times}")
```

### Run Complete Comparison

```bash
cd helios_postprocess/examples
python compare_with_published_integrated.py /path/to/your/Vulcan_HDD_NB_111125_AI_2.exo
```

This will:
1. Extract hot spot time histories
2. Calculate burn-averaged metrics
3. Compare with published table
4. Generate diagnostic plots

## Customization Needed

The `extract_hot_spot_histories()` function in `burn_averaged_metrics.py` makes some assumptions about your ExodusII file structure. You may need to adjust:

### Variable Names

Check what variables your ExodusII file uses:
```python
run = HeliosRun('your_file.exo')
print(run.list_variables())
```

Common variations:
- Temperature: 'temp', 'tev', 'tion', 'tele'
- Density: 'dens', 'density', 'rho'
- Pressure: 'pres', 'pressure', 'p'

Update in `burn_averaged_metrics.py` line ~165:
```python
temp = run.get_variable('temp', time_idx=t_idx)  # Change 'temp' if needed
dens = run.get_variable('dens', time_idx=t_idx)  # Change 'dens' if needed
pres = run.get_variable('pres', time_idx=t_idx)  # Change 'pres' if needed
```

### Coordinate System

The code assumes spherical coordinates. If you have Cartesian:

Update in `burn_averaged_metrics.py` around line ~170:
```python
if 'coordr' in run.coords:
    zone_boundaries = run.coords['coordr']
elif 'coordx' in run.coords:
    # Calculate radii from Cartesian coordinates
    x = run.coords['coordx']
    y = run.coords['coordy'] if 'coordy' in run.coords else np.zeros_like(x)
    z = run.coords['coordz'] if 'coordz' in run.coords else np.zeros_like(x)
    zone_boundaries = np.sqrt(x**2 + y**2 + z**2)
```

### Areal Density Function

The code tries to call `areal_density.calculate_shell_rhoR()`. Check your actual function signature in `areal_density.py` and update the call around line ~210 in `burn_averaged_metrics.py`.

## Integration with Your Workflow

### Option 1: Use HeliosRun directly

```python
from helios_postprocess import (
    HeliosRun,
    extract_hot_spot_histories,
    calculate_burn_averaged_metrics,
    compare_with_published
)

# Your workflow
run = HeliosRun('sim.exo')
histories = extract_hot_spot_histories(run)
metrics = calculate_burn_averaged_metrics(histories)
# ... analyze
```

### Option 2: Import individual modules

```python
from helios_postprocess import core, hot_spot, burn, areal_density
from helios_postprocess.burn_averaged_metrics import calculate_burn_averaged_metrics

# Use modules separately
# ...
```

### Option 3: Jupyter Notebook

```python
# In your Jupyter notebook
%load_ext autoreload
%autoreload 2

from helios_postprocess import *

# Interactive analysis
run = HeliosRun('sim.exo')
# ...
```

## Testing Checklist

- [ ] Package installed successfully
- [ ] Can import HeliosRun
- [ ] Can load your ExodusII file
- [ ] Can list variables in file
- [ ] Can extract temperature/density at one time
- [ ] Hot spot identification works
- [ ] Time history extraction completes
- [ ] Burn-averaged metrics calculated
- [ ] Comparison table generated

## Common Issues

### Import Error
```
ModuleNotFoundError: No module named 'helios_postprocess'
```
**Solution**: Run `pip install -e .` from the package directory

### NetCDF4 not found
```
ModuleNotFoundError: No module named 'netCDF4'
```
**Solution**: Install netCDF4: `pip install netCDF4`

### Variable not found
```
KeyError: "Variable 'temp' not found"
```
**Solution**: Check `run.list_variables()` and update variable names

### Extraction fails
```
Error extracting data at time index X
```
**Solution**: 
1. Check coordinate system assumptions
2. Verify hot spot threshold is appropriate
3. Add debug prints to see what's failing

## Next Steps

1. **Test on one time step first**: Modify example to only process one time to debug faster
2. **Validate hot spot identification**: Plot temperature profiles to ensure hot spot is correctly identified
3. **Check units**: Verify temperature in eV, pressure in J/cm³ or Gbar
4. **Compare metrics**: Run full comparison and check if values are physically reasonable

## Getting Help

If you encounter issues:

1. Check what's in your ExodusII file:
   ```python
   run = HeliosRun('file.exo', verbose=True)
   print(run.list_variables())
   print(run.coords.keys())
   ```

2. Test hot spot ID at one time:
   ```python
   from helios_postprocess import hot_spot
   temp = run.get_variable('temp', time_idx=-1)
   # Check what temperature values you have
   print(f"T range: {temp.min():.1f} - {temp.max():.1f} eV")
   ```

3. Share error messages with full traceback

## Success Criteria

You'll know it's working when you see output like:

```
================================================================================
COMPARISON WITH PUBLISHED DATA
================================================================================
Metric                         Simulation           Published            Δ (%)
--------------------------------------------------------------------------------
⟨T_hs⟩ (keV)                         48.2        46.7±4.8           +3.2
⟨P_hs⟩ (Gbar)                        2650       2720±212            -2.6
⟨ρR_cf⟩ (g/cm²)                      1.55        1.60±0.46          -3.1
CR_max                               19.8        20.1±9.6           -1.5
Yield (MJ)                          245.0       256±0.6            -4.3
================================================================================
```

Good luck! 🚀
