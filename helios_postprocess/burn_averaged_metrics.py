"""
Burn-Averaged Metrics Module (INTEGRATED VERSION)

Calculates burn-averaged (temporal) quantities for comparison with published
ICF performance metrics. Integrates with existing hot_spot, burn, and 
areal_density modules.

This module provides TIME-HISTORY burn-averaging, complementing the existing
SPATIAL neutron-averaging in the core module.

Key Functions
-------------
extract_hot_spot_histories : Extract time histories from HeliosRun
calculate_burn_averaged_metrics : Calculate all burn-averaged quantities
compare_with_published : Generate comparison tables

Author: Prof T
Date: November 2025
"""

import numpy as np
from scipy.integrate import simpson
from typing import Dict, Optional, Tuple
import warnings

# Import existing modules
try:
    from . import hot_spot
    from . import burn as burn_module
    from . import areal_density
    from .core import sigma_v_DT
except ImportError:
    # For standalone use
    import hot_spot
    import burn as burn_module
    import areal_density
    from core import sigma_v_DT


def calculate_burn_rate_from_sim(temperature_keV: np.ndarray,
                                  density_gcc: np.ndarray,
                                  ion_fraction: float = 0.5) -> np.ndarray:
    """
    Calculate DT fusion reaction rate from temperature and density.
    
    Uses existing sigma_v_DT function from core module.
    
    Parameters
    ----------
    temperature_keV : np.ndarray
        Ion temperature in keV at each time
    density_gcc : np.ndarray
        Mass density in g/cm³ at each time
    ion_fraction : float, optional
        Fraction of each ion species (default 0.5 for 50-50 DT)
        
    Returns
    -------
    np.ndarray
        Fusion reaction rate in reactions/(cm³·s)
    """
    # Constants
    amu_to_g = 1.66054e-24  # g
    m_DT = 2.5 * amu_to_g    # Average mass for 50-50 DT
    
    # Number densities (ions/cm³)
    n_total = density_gcc / m_DT
    n_D = n_total * ion_fraction
    n_T = n_total * ion_fraction
    
    # Reactivity using existing function
    sigma_v = sigma_v_DT(temperature_keV)
    
    # Reaction rate
    R_dot = n_D * n_T * sigma_v
    
    return R_dot


def extract_hot_spot_histories(run,
                               T_threshold: float = 1000.0,
                               time_indices: Optional[np.ndarray] = None) -> Dict:
    """
    Extract hot spot time histories from HeliosRun object.
    
    This function extracts all quantities needed for burn-averaged calculations
    using your existing hot_spot and areal_density modules.
    
    Parameters
    ----------
    run : HeliosRun
        HeliosRun object with loaded ExodusII data
    T_threshold : float, optional
        Hot spot temperature threshold in eV (default: 1000 eV = 1 keV)
    time_indices : np.ndarray, optional
        Specific time indices to extract. If None, uses all times.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'time_ns' : Time array in nanoseconds
        - 'temperature_keV' : Hot spot temperature (neutron-averaged)
        - 'pressure_Gbar' : Hot spot pressure (neutron-averaged)
        - 'density_gcc' : Hot spot density (neutron-averaged)
        - 'areal_density_gcm2' : Cold fuel areal density
        - 'radius_um' : Hot spot radius
        - 'volume_cm3' : Hot spot volume
        - 'mass_mg' : Hot spot mass
        - 'initial_radius_um' : Initial hot spot radius
        
    Examples
    --------
    >>> from helios_postprocess import HeliosRun
    >>> run = HeliosRun('sim.exo')
    >>> histories = extract_hot_spot_histories(run, T_threshold=1000.0)
    >>> print(f"Extracted {len(histories['time_ns'])} time steps")
    """
    # Get time array
    times_s = run.times  # Should be in seconds from ExodusII
    
    if time_indices is None:
        time_indices = np.arange(len(times_s))
    
    times_ns = times_s[time_indices] * 1e9  # Convert to ns
    n_times = len(time_indices)
    
    # Initialize arrays
    temperatures = np.zeros(n_times)
    pressures = np.zeros(n_times)
    densities = np.zeros(n_times)
    radii = np.zeros(n_times)
    volumes = np.zeros(n_times)
    masses = np.zeros(n_times)
    rhoR_cf = np.zeros(n_times)
    
    print(f"Extracting hot spot histories ({n_times} time steps)...")
    
    for i, t_idx in enumerate(time_indices):
        try:
            # Get data at this time
            temp = run.get_variable('temp', time_idx=t_idx)  # eV
            dens = run.get_variable('dens', time_idx=t_idx)  # g/cm³
            pres = run.get_variable('pres', time_idx=t_idx)  # J/cm³ (check units!)
            
            # Get zone boundaries for hot spot ID
            # This depends on your coordinate system - adjust as needed
            if 'coordr' in run.coords:
                zone_boundaries = run.coords['coordr']
            elif 'coordx' in run.coords:
                # For Cartesian, would need to calculate radii
                # Placeholder - adjust based on your geometry
                zone_boundaries = np.linspace(0, 1, len(temp)+1)
            else:
                # Fallback
                zone_boundaries = np.arange(len(temp)+1) * 1e-4  # cm
            
            # Identify hot spot using existing function
            hot_mask, r_hs = hot_spot.identify_hot_spot(
                temp, zone_boundaries, T_threshold=T_threshold
            )
            
            radii[i] = r_hs * 1e4  # Convert cm to μm
            
            # Calculate hot spot volume (assuming spherical)
            volumes[i] = (4.0/3.0) * np.pi * r_hs**3  # cm³
            
            # Get hot spot properties using neutron-averaging
            # Need to extract hot spot regions
            temp_hs = temp[hot_mask]
            dens_hs = dens[hot_mask]
            pres_hs = pres[hot_mask]
            
            if len(temp_hs) > 0:
                # Get zone masses (may need to calculate from density and volume)
                # Placeholder - you may need to extract this differently
                zone_mass = dens[hot_mask]  # Simplified
                
                # For DT, assume 50-50 mix
                m_DT = 2.5 * 1.66054e-24  # g
                n_total = dens_hs / m_DT
                n_D = n_total * 0.5
                n_T = n_total * 0.5
                
                # Import neutron averaging function from core
                from .core import get_neutron_averaged_conditions
                
                # Calculate neutron-averaged conditions in hot spot
                avg_conditions = get_neutron_averaged_conditions(
                    temp_hs, dens_hs, pres_hs / 1e6,  # Convert J/cm³ to Gbar
                    n_D, n_T, zone_mass, fuel_type='DT'
                )
                
                temperatures[i] = avg_conditions['T_n_avg'] / 1000.0  # eV to keV
                densities[i] = avg_conditions['rho_n_avg']  # g/cm³
                pressures[i] = avg_conditions['P_n_avg']  # Gbar
                
                # Hot spot mass
                masses[i] = np.sum(dens_hs * zone_mass) * 1000  # Convert g to mg
            else:
                # No hot spot identified
                temperatures[i] = 0.0
                densities[i] = 0.0
                pressures[i] = 0.0
                masses[i] = 0.0
            
            # Calculate cold fuel areal density
            # Use your existing areal_density module
            # This is a simplified version - adjust based on your actual function signatures
            try:
                # Assuming you have a function to calculate ρR of cold fuel
                # Adjust this call based on your actual areal_density API
                rhoR_cf[i] = areal_density.calculate_shell_rhoR(
                    dens, zone_boundaries, T_threshold=T_threshold
                )
            except Exception as e:
                warnings.warn(f"Could not calculate cold fuel ρR at t_idx={t_idx}: {e}")
                rhoR_cf[i] = 0.0
            
            if (i + 1) % 10 == 0:
                print(f"  Progress: {i+1}/{n_times}")
                
        except Exception as e:
            warnings.warn(f"Error extracting data at time index {t_idx}: {e}")
            continue
    
    print("✓ Extraction complete")
    
    return {
        'time_ns': times_ns,
        'temperature_keV': temperatures,
        'pressure_Gbar': pressures,
        'density_gcc': densities,
        'areal_density_gcm2': rhoR_cf,
        'radius_um': radii,
        'volume_cm3': volumes,
        'mass_mg': masses,
        'initial_radius_um': radii[0]
    }


def calculate_burn_averaged_metrics(histories: Dict,
                                    ion_fraction: float = 0.5) -> Dict:
    """
    Calculate burn-averaged metrics from time histories.
    
    Parameters
    ----------
    histories : dict
        Dictionary from extract_hot_spot_histories()
    ion_fraction : float, optional
        Ion fraction for burn rate calculation (default 0.5)
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'T_burn_avg' : Burn-averaged temperature (keV)
        - 'P_burn_avg' : Burn-averaged pressure (Gbar)
        - 'rho_burn_avg' : Burn-averaged density (g/cm³)
        - 'rhoR_burn_avg' : Burn-averaged areal density (g/cm²)
        - 'CR_max' : Maximum convergence ratio
        - 'yield_MJ' : Total fusion yield (MJ)
        - 'T_peak', 'P_peak', etc. : Peak values
        
    Examples
    --------
    >>> metrics = calculate_burn_averaged_metrics(histories)
    >>> print(f"⟨T_hs⟩ = {metrics['T_burn_avg']:.1f} keV")
    """
    # Extract time histories
    time_ns = histories['time_ns']
    temp_keV = histories['temperature_keV']
    pres_Gbar = histories['pressure_Gbar']
    dens_gcc = histories['density_gcc']
    rhoR_gcm2 = histories['areal_density_gcm2']
    radius_um = histories['radius_um']
    volume_cm3 = histories['volume_cm3']
    
    # Convert time to seconds for integration
    time_s = time_ns * 1e-9
    
    # Calculate burn rate at each time
    burn_rate = calculate_burn_rate_from_sim(temp_keV, dens_gcc, ion_fraction)
    
    # Calculate burn-averaged temperature
    T_num = simpson(temp_keV * burn_rate, x=time_s)
    T_den = simpson(burn_rate, x=time_s)
    T_burn_avg = T_num / T_den if T_den > 0 else 0.0
    
    # Calculate burn-averaged pressure
    P_num = simpson(pres_Gbar * burn_rate, x=time_s)
    P_den = simpson(burn_rate, x=time_s)
    P_burn_avg = P_num / P_den if P_den > 0 else 0.0
    
    # Calculate burn-averaged density
    rho_num = simpson(dens_gcc * burn_rate, x=time_s)
    rho_den = simpson(burn_rate, x=time_s)
    rho_burn_avg = rho_num / rho_den if rho_den > 0 else 0.0
    
    # Calculate burn-averaged areal density
    rhoR_num = simpson(rhoR_gcm2 * burn_rate, x=time_s)
    rhoR_den = simpson(burn_rate, x=time_s)
    rhoR_burn_avg = rhoR_num / rhoR_den if rhoR_den > 0 else 0.0
    
    # Calculate convergence ratio
    if histories['initial_radius_um'] > 0:
        CR_history = histories['initial_radius_um'] / radius_um
        CR_max = np.max(CR_history[radius_um > 0])
    else:
        CR_max = 0.0
    
    # Calculate total yield
    # Each DT fusion releases 17.6 MeV = 2.82e-12 J
    E_fusion = 17.6 * 1.60218e-13  # MeV to J
    
    # Total reactions = ∫ Ṙ(t) · V(t) dt
    reaction_rate_total = burn_rate * volume_cm3
    total_reactions = simpson(reaction_rate_total, x=time_s)
    yield_J = total_reactions * E_fusion
    yield_MJ = yield_J * 1e-6
    
    # Burn fraction (normalized)
    burn_fraction = burn_rate / np.sum(burn_rate) if np.sum(burn_rate) > 0 else burn_rate
    
    return {
        # Burn-averaged values
        'T_burn_avg': T_burn_avg,
        'P_burn_avg': P_burn_avg,
        'rho_burn_avg': rho_burn_avg,
        'rhoR_burn_avg': rhoR_burn_avg,
        
        # Peak values
        'T_peak': np.max(temp_keV),
        'P_peak': np.max(pres_Gbar),
        'rho_peak': np.max(dens_gcc),
        'rhoR_peak': np.max(rhoR_gcm2),
        
        # Convergence and yield
        'CR_max': CR_max,
        'min_radius_um': np.min(radius_um[radius_um > 0]),
        'yield_MJ': yield_MJ,
        'yield_J': yield_J,
        'total_reactions': total_reactions,
        
        # Diagnostic info
        'burn_fraction': burn_fraction,
        'burn_rate': burn_rate
    }


def compare_with_published(sim_metrics: Dict,
                          published_metrics: Dict,
                          laser_energy_MJ: float = 4.0) -> str:
    """
    Generate comparison table between simulation and published results.
    
    Parameters
    ----------
    sim_metrics : dict
        Calculated metrics from simulation
    published_metrics : dict
        Published values with uncertainties
        Format: {metric_name: (value, uncertainty)}
    laser_energy_MJ : float
        Laser energy in MJ
        
    Returns
    -------
    str
        Formatted comparison table
        
    Examples
    --------
    >>> published = {
    ...     'T_hs': (46.7, 4.8),
    ...     'P_hs': (2720, 212),
    ...     'rhoR_cf': (1.60, 0.46),
    ...     'CR_max': (20.1, 9.6),
    ...     'yield': (256, 0.6),
    ...     'gain': (65, 0)
    ... }
    >>> table = compare_with_published(metrics, published, laser_energy_MJ=4.0)
    >>> print(table)
    """
    # Calculate fusion gain
    sim_gain = sim_metrics['yield_MJ'] / laser_energy_MJ if laser_energy_MJ > 0 else 0.0
    pub_gain, pub_gain_unc = published_metrics.get('gain', (0.0, 0.0))
    
    # Build table
    lines = []
    lines.append("=" * 80)
    lines.append("COMPARISON WITH PUBLISHED DATA")
    lines.append("=" * 80)
    lines.append(f"{'Metric':<30} {'Simulation':<20} {'Published':<20} {'Δ (%)':<10}")
    lines.append("-" * 80)
    
    # Temperature
    sim_T = sim_metrics['T_burn_avg']
    pub_T, pub_T_unc = published_metrics.get('T_hs', (0.0, 0.0))
    if pub_T > 0:
        delta_T = 100 * (sim_T - pub_T) / pub_T
        lines.append(f"{'⟨T_hs⟩ (keV)':<30} {sim_T:>19.1f} {pub_T:>12.1f}±{pub_T_unc:<6.1f} {delta_T:>9.1f}")
    
    # Pressure
    sim_P = sim_metrics['P_burn_avg']
    pub_P, pub_P_unc = published_metrics.get('P_hs', (0.0, 0.0))
    if pub_P > 0:
        delta_P = 100 * (sim_P - pub_P) / pub_P
        lines.append(f"{'⟨P_hs⟩ (Gbar)':<30} {sim_P:>19.0f} {pub_P:>12.0f}±{pub_P_unc:<6.0f} {delta_P:>9.1f}")
    
    # Areal density
    sim_rhoR = sim_metrics['rhoR_burn_avg']
    pub_rhoR, pub_rhoR_unc = published_metrics.get('rhoR_cf', (0.0, 0.0))
    if pub_rhoR > 0:
        delta_rhoR = 100 * (sim_rhoR - pub_rhoR) / pub_rhoR
        lines.append(f"{'⟨ρR_cf⟩ (g/cm²)':<30} {sim_rhoR:>19.2f} {pub_rhoR:>12.2f}±{pub_rhoR_unc:<6.2f} {delta_rhoR:>9.1f}")
    
    # Convergence ratio
    sim_CR = sim_metrics['CR_max']
    pub_CR, pub_CR_unc = published_metrics.get('CR_max', (0.0, 0.0))
    if pub_CR > 0:
        delta_CR = 100 * (sim_CR - pub_CR) / pub_CR
        lines.append(f"{'CR_max':<30} {sim_CR:>19.1f} {pub_CR:>12.1f}±{pub_CR_unc:<6.1f} {delta_CR:>9.1f}")
    
    # Yield
    sim_Y = sim_metrics['yield_MJ']
    pub_Y, pub_Y_unc = published_metrics.get('yield', (0.0, 0.0))
    if pub_Y > 0:
        delta_Y = 100 * (sim_Y - pub_Y) / pub_Y
        lines.append(f"{'Yield (MJ)':<30} {sim_Y:>19.1f} {pub_Y:>12.1f}±{pub_Y_unc:<6.1f} {delta_Y:>9.1f}")
    
    # Gain
    if pub_gain > 0:
        delta_gain = 100 * (sim_gain - pub_gain) / pub_gain
        lines.append(f"{'Fusion Gain':<30} {sim_gain:>19.1f} {pub_gain:>12.1f}±{pub_gain_unc:<6.1f} {delta_gain:>9.1f}")
    
    lines.append("=" * 80)
    
    return "\n".join(lines)
