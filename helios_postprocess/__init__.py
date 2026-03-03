"""
Helios Postprocessing Package

A comprehensive Python package for analyzing Helios ICF simulation outputs
with integrated diagnostics including neutron downscatter ratio (DSR),
areal density, hot spot metrics, burn diagnostics, and burn-averaged
performance metrics for comparison with published data.

Main Classes
-----------
HeliosRun : Primary interface for simulation analysis
    Provides integrated access to all physics modules

Core Modules
------------
- core : ExodusII reading and HeliosRun class
- hot_spot : Hot spot identification and metrics
- burn : Burn diagnostics (bang time, yield, burn width)
- areal_density : Areal density (ρR) calculations
- neutron_downscatter : Neutron DSR diagnostics
- pressure_gradients : Shock and RT instability analysis
- energetics : Kinetic energy and hydro efficiency
- burn_averaged_metrics : Burn-averaged quantities for published comparisons

Quick Start - Basic Analysis
-----------------------------
>>> from helios_postprocess import HeliosRun
>>> 
>>> # Load simulation
>>> run = HeliosRun('Vulcan_HDD_NB_111125_AI_2.exo')
>>> 
>>> # Calculate neutron downscatter ratio
>>> dsr = run.calculate_dsr(time_idx=-1)
>>> print(f"DSR = {dsr['DSR']:.3f}")
>>> 
>>> # Convert to areal density
>>> rho_R = run.dsr_to_areal_density(dsr['DSR'], calibration='NIF')
>>> print(f"ρR = {rho_R['rho_R']:.2f} ± {rho_R['uncertainty']:.2f} g/cm²")
>>> 
>>> # Check ignition
>>> ignition = run.check_ignition_criteria()
>>> print(f"Ignition: {ignition['ignition_achieved']}")

Quick Start - Burn-Averaged Comparison
---------------------------------------
>>> from helios_postprocess import HeliosRun
>>> from helios_postprocess.burn_averaged_metrics import (
...     extract_hot_spot_histories,
...     calculate_burn_averaged_metrics,
...     compare_with_published
... )
>>> 
>>> # Load simulation
>>> run = HeliosRun('sim.exo')
>>> 
>>> # Extract time histories
>>> histories = extract_hot_spot_histories(run, T_threshold=1000.0)
>>> 
>>> # Calculate burn-averaged metrics
>>> metrics = calculate_burn_averaged_metrics(histories)
>>> 
>>> # Define published values for comparison
>>> published = {
...     'T_hs': (46.7, 4.8),      # ⟨T_hs⟩ in keV
...     'P_hs': (2720, 212),      # ⟨P_hs⟩ in Gbar
...     'rhoR_cf': (1.60, 0.46),  # ⟨ρR_cf⟩ in g/cm²
...     'CR_max': (20.1, 9.6),    # Maximum convergence ratio
...     'yield': (256, 0.6),      # Yield in MJ
...     'gain': (65, 0)           # Fusion gain
... }
>>> 
>>> # Generate comparison table
>>> comparison = compare_with_published(metrics, published, laser_energy_MJ=4.0)
>>> print(comparison)

Module Organization
-------------------

CORE FUNCTIONALITY:
    core.py
        - HeliosRun : Main analysis class
        - ExodusII file reading
        - Variable extraction
        - Coordinate handling

SPATIAL ANALYSIS (at each time):
    hot_spot.py
        - identify_hot_spot() : Find hot spot region
        - hot_spot_mass() : Calculate hot spot mass
        - mass_averaged_temperature() : Mass-averaged T
        - neutron_averaged_temperature() : Neutron-averaged T
        
    areal_density.py
        - calculate_shell_rhoR() : Shell areal density
        - calculate_line_rhoR() : Line-of-sight ρR
        - calculate_hot_spot_rhoR() : Hot spot ρR
        
    pressure_gradients.py
        - identify_shocks() : Shock identification
        - calculate_RT_growth() : Rayleigh-Taylor instability
        
    neutron_downscatter.py
        - calculate_dsr() : Downscatter ratio
        - dsr_to_areal_density() : DSR → ρR conversion

TEMPORAL ANALYSIS (over time):
    burn.py
        - calculate_neutron_yield() : Total yield
        - find_bang_time() : Peak burn time
        - calculate_burn_width() : Burn pulse width
        
    burn_averaged_metrics.py  ← NEW!
        - extract_hot_spot_histories() : Get time histories
        - calculate_burn_averaged_metrics() : ⟨T⟩, ⟨P⟩, ⟨ρR⟩
        - compare_with_published() : Comparison tables
        
    energetics.py
        - calculate_kinetic_energy() : Implosion KE
        - calculate_hydro_efficiency() : Drive efficiency

Physics Constants
-----------------
- DT fusion energy: 17.6 MeV per reaction
- Ignition criteria (Lindl):
    * Temperature: >5 keV
    * Pressure: >100 Gbar
    * Areal density: >0.3 g/cm²
    * Convergence ratio: >15

Unit Conventions
----------------
- Energy: J or MJ (use *1e6 for erg conversion)
- Temperature: eV or keV
- Pressure: J/cm³ (use /1e6 for Gbar)
- Density: g/cm³
- Areal density: g/cm²
- Time: s or ns
- Length: cm or μm
- Mass: g or mg

Author: Prof T
Version: 2.0 (November 2025) - With Burn-Averaged Metrics
"""

# Core functionality
from .core import HeliosRun

# Make key functions easily accessible
from .burn_averaged_metrics import (
    extract_hot_spot_histories,
    calculate_burn_averaged_metrics,
    compare_with_published
)

__version__ = '2.0.0'
__author__ = 'Prof T'

__all__ = [
    'HeliosRun',
    'extract_hot_spot_histories',
    'calculate_burn_averaged_metrics',
    'compare_with_published'
]

# Module information for programmatic access
MODULES = {
    'core': 'ExodusII reading and HeliosRun class',
    'hot_spot': 'Hot spot identification and metrics',
    'burn': 'Burn diagnostics (bang time, yield, burn width)',
    'areal_density': 'Areal density (ρR) calculations',
    'neutron_downscatter': 'Neutron DSR diagnostics',
    'pressure_gradients': 'Shock and RT instability analysis',
    'energetics': 'Kinetic energy and hydro efficiency',
    'burn_averaged_metrics': 'Burn-averaged quantities (NEW!)'
}

# Ignition criteria (Lindl)
LINDL_CRITERIA = {
    'temperature_keV': 5.0,
    'pressure_Gbar': 100.0,
    'areal_density_gcm2': 0.3,
    'convergence_ratio': 15.0
}

def print_module_info():
    """Print information about available modules."""
    print("=" * 70)
    print("HELIOS POSTPROCESS - Available Modules")
    print("=" * 70)
    for module, description in MODULES.items():
        print(f"  {module:<25} : {description}")
    print("=" * 70)
    print(f"\nVersion: {__version__}")
    print(f"Author: {__author__}")

def print_ignition_criteria():
    """Print Lindl ignition criteria."""
    print("=" * 70)
    print("LINDL IGNITION CRITERIA")
    print("=" * 70)
    print(f"  Temperature:      > {LINDL_CRITERIA['temperature_keV']:.1f} keV")
    print(f"  Pressure:         > {LINDL_CRITERIA['pressure_Gbar']:.0f} Gbar")
    print(f"  Areal Density:    > {LINDL_CRITERIA['areal_density_gcm2']:.1f} g/cm²")
    print(f"  Convergence Ratio: > {LINDL_CRITERIA['convergence_ratio']:.0f}")
    print("=" * 70)
from .data_builder import build_run_data, ICFRunData
