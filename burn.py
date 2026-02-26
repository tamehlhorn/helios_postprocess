"""
Burn diagnostics for ICF simulations

Functions for analyzing fusion burn, neutron yield, bang time, burn width,
and stagnation timing.
"""

import numpy as np
from typing import Tuple, Optional, Dict


def calculate_total_fusion_rate(fusion_rate_density: np.ndarray,
                                 zone_mass: np.ndarray) -> np.ndarray:
    """
    Calculate total fusion rate by integrating over all zones
    
    Converts fusion rate density [neutrons/(s·g)] to total rate [neutrons/s]
    
    Parameters
    ----------
    fusion_rate_density : np.ndarray
        Fusion rate per unit mass [neutrons/(s·g)], shape (n_times, n_zones)
    zone_mass : np.ndarray
        Zone masses [g], shape (n_times, n_zones)
        
    Returns
    -------
    np.ndarray
        Total fusion rate [neutrons/s], shape (n_times,)
        
    Examples
    --------
    >>> fusion_rate = calculate_total_fusion_rate(
    ...     run.reader.get_variable('fusion_rate'),
    ...     run.zone_mass
    ... )
    >>> print(f"Peak fusion rate: {fusion_rate.max():.2e} neutrons/s")
    
    Notes
    -----
    This is the instantaneous neutron production rate at each time.
    Integrate over time to get total neutron yield.
    """
    # Vectorized: sum over zones for each time
    total_rate = np.sum(fusion_rate_density * zone_mass, axis=1)
    return total_rate


def calculate_neutron_yield(fusion_rate: np.ndarray,
                            time: np.ndarray) -> float:
    """
    Calculate total neutron yield by integrating fusion rate over time
    
    Yield = ∫ fusion_rate(t) dt
    
    Parameters
    ----------
    fusion_rate : np.ndarray
        Total fusion rate [neutrons/s], shape (n_times,)
    time : np.ndarray
        Time [s], shape (n_times,)
        
    Returns
    -------
    float
        Total neutron yield [neutrons]
        
    Examples
    --------
    >>> yield_n = calculate_neutron_yield(fusion_rate, time_s)
    >>> print(f"Neutron yield: {yield_n:.2e} neutrons")
    
    Notes
    -----
    Uses trapezoidal integration for accuracy.
    """
    # Trapezoidal rule: more accurate than simple sum
    yield_total = np.trapz(fusion_rate, time)
    return yield_total


def calculate_fusion_energy_output(neutron_yield: float,
                                   energy_per_neutron: float = 17.6) -> float:
    """
    Calculate total fusion energy output from neutron yield
    
    For D-T fusion: Each reaction releases 17.6 MeV
    (14.1 MeV neutron + 3.5 MeV alpha)
    
    Parameters
    ----------
    neutron_yield : float
        Total number of fusion reactions [neutrons]
    energy_per_neutron : float, optional
        Energy per fusion reaction [MeV], default 17.6 for D-T
        
    Returns
    -------
    float
        Total fusion energy [J]
        
    Examples
    --------
    >>> energy_J = calculate_fusion_energy_output(1e15)
    >>> print(f"Fusion energy: {energy_J * 1e-6:.2f} MJ")
    
    Notes
    -----
    Default assumes D-T fusion (17.6 MeV per reaction).
    For D-D: use energy_per_neutron = 3.27 MeV
    For D-He3: use energy_per_neutron = 18.3 MeV
    """
    MeV_to_J = 1.602e-13  # Conversion factor
    energy_J = neutron_yield * energy_per_neutron * MeV_to_J
    return energy_J


def find_bang_time(fusion_rate: np.ndarray,
                  time: np.ndarray) -> Tuple[float, int, float]:
    """
    Find bang time (time of peak fusion rate)
    
    Parameters
    ----------
    fusion_rate : np.ndarray
        Total fusion rate [neutrons/s], shape (n_times,)
    time : np.ndarray
        Time [ns], shape (n_times,)
        
    Returns
    -------
    bang_time : float
        Time of peak fusion rate [ns]
    bang_idx : int
        Index of bang time
    peak_rate : float
        Peak fusion rate [neutrons/s]
        
    Examples
    --------
    >>> t_bang, idx, rate = find_bang_time(fusion_rate, time)
    >>> print(f"Bang time: {t_bang:.2f} ns")
    >>> print(f"Peak rate: {rate:.2e} neutrons/s")
    
    Notes
    -----
    Bang time is the moment of peak fusion power, typically near or
    slightly after stagnation in ICF. If bang time is before stagnation,
    it may indicate missing physics (e.g., alpha transport).
    """
    bang_idx = np.argmax(fusion_rate)
    bang_time = time[bang_idx]
    peak_rate = fusion_rate[bang_idx]
    
    return bang_time, bang_idx, peak_rate


def calculate_burn_width(fusion_rate: np.ndarray,
                         time: np.ndarray,
                         fraction: float = 0.5) -> float:
    """
    Calculate burn width (FWHM of fusion rate curve)
    
    Full Width at Half Maximum (FWHM) of the fusion rate time history.
    
    Parameters
    ----------
    fusion_rate : np.ndarray
        Total fusion rate [neutrons/s], shape (n_times,)
    time : np.ndarray
        Time [ns], shape (n_times,)
    fraction : float, optional
        Fraction of peak for width calculation, default 0.5 (FWHM)
        Use 0.1 for Full Width at Tenth Maximum (FWTM)
        
    Returns
    -------
    float
        Burn width [ns]
        
    Examples
    --------
    >>> burn_width = calculate_burn_width(fusion_rate, time)
    >>> print(f"Burn width (FWHM): {burn_width:.2f} ns")
    >>> 
    >>> # Or use tenth-maximum for broader pulse
    >>> fwtm = calculate_burn_width(fusion_rate, time, fraction=0.1)
    
    Notes
    -----
    Burn width indicates the duration of significant fusion activity.
    Shorter burn width generally indicates better compression and timing.
    Typical ICF: 50-200 ps for direct drive, 100-300 ps for indirect drive.
    """
    max_rate = np.max(fusion_rate)
    threshold = fraction * max_rate
    
    # Find indices where rate exceeds threshold
    above_threshold = fusion_rate >= threshold
    
    if not np.any(above_threshold):
        return 0.0  # No burn
    
    # Find first crossing (rising edge)
    rising_idx = np.where(above_threshold)[0][0]
    
    # Find last crossing (falling edge)
    falling_idx = np.where(above_threshold)[0][-1]
    
    # If we need to interpolate for more accuracy
    if rising_idx > 0 and falling_idx < len(time) - 1:
        # Linear interpolation at threshold crossing
        # Rising edge
        t0, t1 = time[rising_idx-1], time[rising_idx]
        r0, r1 = fusion_rate[rising_idx-1], fusion_rate[rising_idx]
        t_rise = t0 + (t1 - t0) * (threshold - r0) / (r1 - r0)
        
        # Falling edge
        t0, t1 = time[falling_idx], time[falling_idx+1]
        r0, r1 = fusion_rate[falling_idx], fusion_rate[falling_idx+1]
        t_fall = t0 + (t1 - t0) * (threshold - r0) / (r1 - r0)
        
        burn_width = t_fall - t_rise
    else:
        # Simple calculation without interpolation
        burn_width = time[falling_idx] - time[rising_idx]
    
    return burn_width


def calculate_target_gain(fusion_energy: float,
                         laser_energy: float) -> float:
    """
    Calculate target gain (fusion energy out / laser energy in)
    
    Parameters
    ----------
    fusion_energy : float
        Total fusion energy output [J]
    laser_energy : float
        Total laser energy delivered [J]
        
    Returns
    -------
    float
        Target gain (dimensionless)
        
    Examples
    --------
    >>> gain = calculate_target_gain(energy_out, laser_in)
    >>> print(f"Target gain: {gain:.2f}")
    
    Notes
    -----
    Target gain = E_fusion / E_laser
    
    Gain regimes:
    - Gain < 1: No net energy gain
    - Gain = 1: Breakeven
    - Gain > 1: Net energy gain
    - Gain > 10: Scientific breakeven (compensates for driver efficiency)
    - Gain > 100: Potential for energy production
    
    NIF achieved gain > 1 in December 2022.
    """
    if laser_energy <= 0:
        return 0.0
    
    gain = fusion_energy / laser_energy
    return gain


def find_stagnation_time(density: np.ndarray,
                         time: np.ndarray,
                         method: str = 'max_density') -> Tuple[float, int]:
    """
    Find stagnation time (moment of peak compression)
    
    Stagnation is when the implosion reaches maximum compression,
    typically identified by peak density.
    
    Parameters
    ----------
    density : np.ndarray
        Density array [g/cm³], shape (n_times, n_zones)
    time : np.ndarray
        Time [ns], shape (n_times,)
    method : str, optional
        Method for finding stagnation:
        - 'max_density': Time of maximum density anywhere (default)
        - 'max_avg_density': Time of maximum mass-averaged density
        - 'min_radius': Time of minimum inner radius (requires radius data)
        
    Returns
    -------
    stag_time : float
        Stagnation time [ns]
    stag_idx : int
        Index of stagnation time
        
    Examples
    --------
    >>> stag_time, idx = find_stagnation_time(run.density, run.time)
    >>> print(f"Stagnation at {stag_time:.2f} ns")
    >>> print(f"Peak density: {run.density[idx].max():.2f} g/cm³")
    
    Notes
    -----
    Stagnation time is critical for ICF analysis. The fuel should ideally
    bang (peak fusion rate) at or slightly after stagnation for maximum
    compression and heating.
    """
    if method == 'max_density':
        # Find time of maximum density anywhere in the simulation
        max_densities = np.max(density, axis=1)
        stag_idx = np.argmax(max_densities)
        
    elif method == 'max_avg_density':
        # Find time of maximum spatially-averaged density
        avg_densities = np.mean(density, axis=1)
        stag_idx = np.argmax(avg_densities)
        
    else:
        raise ValueError(f"Unknown method: {method}")
    
    stag_time = time[stag_idx]
    return stag_time, stag_idx


def find_stagnation_time_mass_averaged(density: np.ndarray,
                                       zone_mass: np.ndarray,
                                       time: np.ndarray,
                                       region_mask: Optional[np.ndarray] = None) -> Tuple[float, int, float]:
    """
    Find stagnation time using mass-averaged density
    
    More sophisticated than simple maximum - uses mass weighting.
    Can optionally focus on specific region (e.g., fuel only).
    
    Parameters
    ----------
    density : np.ndarray
        Density [g/cm³], shape (n_times, n_zones)
    zone_mass : np.ndarray
        Zone masses [g], shape (n_times, n_zones)
    time : np.ndarray
        Time [ns], shape (n_times,)
    region_mask : np.ndarray, optional
        Boolean mask for region of interest, shape (n_zones,)
        If provided, only considers these zones
        
    Returns
    -------
    stag_time : float
        Stagnation time [ns]
    stag_idx : int
        Index of stagnation
    max_density : float
        Peak mass-averaged density [g/cm³]
        
    Examples
    --------
    >>> # Stagnation of full implosion
    >>> t_stag, idx, rho_max = find_stagnation_time_mass_averaged(
    ...     run.density,
    ...     run.zone_mass,
    ...     run.time
    ... )
    >>> 
    >>> # Stagnation of fuel region only
    >>> fuel_mask = (run.zone_centers[0] < 0.05)  # r < 0.05 cm
    >>> t_stag, idx, rho_max = find_stagnation_time_mass_averaged(
    ...     run.density,
    ...     run.zone_mass,
    ...     run.time,
    ...     region_mask=fuel_mask
    ... )
    
    Notes
    -----
    Mass-averaged density = Σ(ρ_i * m_i) / Σ(m_i)
    
    This is more physically meaningful than simple maximum density,
    especially when tracking compression of specific regions.
    """
    n_times = density.shape[0]
    avg_density = np.zeros(n_times)
    
    for t in range(n_times):
        rho = density[t]
        mass = zone_mass[t]
        
        if region_mask is not None:
            rho = rho[region_mask]
            mass = mass[region_mask]
        
        # Mass-weighted average
        avg_density[t] = np.sum(rho * mass) / np.sum(mass)
    
    stag_idx = np.argmax(avg_density)
    stag_time = time[stag_idx]
    max_density = avg_density[stag_idx]
    
    return stag_time, stag_idx, max_density


def find_stagnation_by_radius(zone_boundaries: np.ndarray,
                              time: np.ndarray,
                              density_threshold: float = 1.0) -> Tuple[float, int, float]:
    """
    Find stagnation time by minimum radius of dense region
    
    Tracks the innermost boundary of a region exceeding density threshold.
    Stagnation is when this radius is minimum (maximum compression).
    
    Parameters
    ----------
    zone_boundaries : np.ndarray
        Zone boundaries [cm], shape (n_times, n_zones+1)
    time : np.ndarray
        Time [ns], shape (n_times,)
    density_threshold : float, optional
        Density threshold [g/cm³] defining "dense region", default 1.0
        
    Returns
    -------
    stag_time : float
        Stagnation time [ns]
    stag_idx : int
        Index of stagnation
    min_radius : float
        Minimum inner radius [cm]
        
    Examples
    --------
    >>> t_stag, idx, r_min = find_stagnation_by_radius(
    ...     run.zone_boundaries,
    ...     run.time,
    ...     density_threshold=10.0
    ... )
    >>> print(f"Minimum radius: {r_min:.4f} cm at {t_stag:.2f} ns")
    
    Notes
    -----
    This method is useful for tracking the actual geometric compression
    of the fuel region, rather than just density magnitudes.
    """
    # This is a placeholder - full implementation would need density data too
    # For now, just track minimum inner radius
    inner_radius = zone_boundaries[:, 0]
    stag_idx = np.argmin(inner_radius)
    stag_time = time[stag_idx]
    min_radius = inner_radius[stag_idx]
    
    return stag_time, stag_idx, min_radius


def calculate_convergence_ratio(zone_boundaries: np.ndarray,
                                stagnation_idx: int,
                                initial_idx: int = 0) -> float:
    """
    Calculate convergence ratio (initial radius / stagnation radius)
    
    Key metric for ICF performance. Higher convergence typically means
    better compression.
    
    Parameters
    ----------
    zone_boundaries : np.ndarray
        Zone boundaries [cm], shape (n_times, n_zones+1)
    stagnation_idx : int
        Time index of stagnation
    initial_idx : int, optional
        Initial time index, default 0
        
    Returns
    -------
    float
        Convergence ratio (dimensionless)
        
    Examples
    --------
    >>> stag_time, stag_idx = find_stagnation_time(run.density, run.time)
    >>> CR = calculate_convergence_ratio(run.zone_boundaries, stag_idx)
    >>> print(f"Convergence ratio: {CR:.1f}")
    
    Notes
    -----
    Typical ICF convergence ratios:
    - Direct drive: 25-40
    - Indirect drive: 30-45
    - Higher CR generally better, but stability issues increase
    """
    r_initial = zone_boundaries[initial_idx, 0]  # Inner radius at start
    r_stagnation = zone_boundaries[stagnation_idx, 0]  # At stagnation
    
    if r_stagnation > 0:
        return r_initial / r_stagnation
    else:
        return np.inf


def comprehensive_stagnation_analysis(density: np.ndarray,
                                     zone_mass: np.ndarray,
                                     zone_boundaries: np.ndarray,
                                     time: np.ndarray) -> Dict:
    """
    Perform complete stagnation analysis
    
    Parameters
    ----------
    density : np.ndarray
        Density [g/cm³], shape (n_times, n_zones)
    zone_mass : np.ndarray
        Zone masses [g], shape (n_times, n_zones)
    zone_boundaries : np.ndarray
        Zone boundaries [cm], shape (n_times, n_zones+1)
    time : np.ndarray
        Time [ns], shape (n_times,)
        
    Returns
    -------
    dict
        Dictionary with:
        - 'stagnation_time': Time of stagnation [ns]
        - 'stagnation_idx': Index of stagnation
        - 'max_density': Peak density at stagnation [g/cm³]
        - 'max_avg_density': Peak mass-averaged density [g/cm³]
        - 'min_radius': Minimum inner radius [cm]
        - 'convergence_ratio': Initial/stagnation radius ratio
        
    Examples
    --------
    >>> results = comprehensive_stagnation_analysis(
    ...     run.density,
    ...     run.zone_mass,
    ...     run.zone_boundaries,
    ...     run.time
    ... )
    >>> print(f"Stagnation: {results['stagnation_time']:.2f} ns")
    >>> print(f"Convergence ratio: {results['convergence_ratio']:.1f}")
    """
    results = {}
    
    # Find stagnation by max density
    stag_time, stag_idx = find_stagnation_time(density, time, method='max_density')
    results['stagnation_time'] = stag_time
    results['stagnation_idx'] = stag_idx
    
    # Peak density at stagnation
    results['max_density'] = density[stag_idx].max()
    
    # Mass-averaged density at stagnation
    _, _, max_avg_dens = find_stagnation_time_mass_averaged(
        density, zone_mass, time
    )
    results['max_avg_density'] = max_avg_dens
    
    # Minimum radius
    results['min_radius'] = zone_boundaries[stag_idx, 0]
    
    # Convergence ratio
    results['convergence_ratio'] = calculate_convergence_ratio(
        zone_boundaries, stag_idx
    )
    
    return results


def calculate_target_gain(fusion_energy: float,
                         laser_energy: float) -> float:
    """
    Calculate target gain (fusion energy out / laser energy in)
    
    Parameters
    ----------
    fusion_energy : float
        Total fusion energy output [J]
    laser_energy : float
        Total laser energy delivered [J]
        
    Returns
    -------
    float
        Target gain (dimensionless)
        
    Examples
    --------
    >>> gain = calculate_target_gain(energy_out, laser_in)
    >>> print(f"Target gain: {gain:.2f}")
    
    Notes
    -----
    Target gain = E_fusion / E_laser
    
    Gain regimes:
    - Gain < 1: No net energy gain
    - Gain = 1: Breakeven
    - Gain > 1: Net energy gain
    - Gain > 10: Scientific breakeven (compensates for driver efficiency)
    - Gain > 100: Potential for energy production
    
    NIF achieved gain > 1 in December 2022.
    """
    if laser_energy <= 0:
        return 0.0
    
    gain = fusion_energy / laser_energy
    return gain


def comprehensive_burn_analysis(fusion_rate_density: np.ndarray,
                                zone_mass: np.ndarray,
                                time_s: np.ndarray,
                                time_ns: np.ndarray,
                                laser_energy: float = None) -> dict:
    """
    Perform complete burn analysis and return all metrics
    
    Convenience function that calculates all burn diagnostics at once.
    
    Parameters
    ----------
    fusion_rate_density : np.ndarray
        Fusion rate per mass [neutrons/(s·g)], shape (n_times, n_zones)
    zone_mass : np.ndarray
        Zone masses [g], shape (n_times, n_zones)
    time_s : np.ndarray
        Time in seconds for integration, shape (n_times,)
    time_ns : np.ndarray
        Time in nanoseconds for reporting, shape (n_times,)
    laser_energy : float, optional
        Total laser energy [J] for gain calculation
        
    Returns
    -------
    dict
        Dictionary containing all burn metrics:
        - 'fusion_rate': Total fusion rate vs time [neutrons/s]
        - 'neutron_yield': Total neutron yield [neutrons]
        - 'fusion_energy': Fusion energy output [J]
        - 'bang_time': Bang time [ns]
        - 'bang_idx': Index of bang time
        - 'peak_rate': Peak fusion rate [neutrons/s]
        - 'burn_width': Burn width FWHM [ns]
        - 'target_gain': Fusion/laser energy ratio (if laser_energy provided)
        
    Examples
    --------
    >>> results = comprehensive_burn_analysis(
    ...     fusion_rate_dens,
    ...     zone_mass,
    ...     time_s,
    ...     time_ns,
    ...     laser_energy=1e6
    ... )
    >>> print(f"Bang time: {results['bang_time']:.2f} ns")
    >>> print(f"Neutron yield: {results['neutron_yield']:.2e}")
    >>> print(f"Target gain: {results['target_gain']:.3f}")
    """
    results = {}
    
    # Total fusion rate
    results['fusion_rate'] = calculate_total_fusion_rate(
        fusion_rate_density, zone_mass
    )
    
    # Neutron yield
    results['neutron_yield'] = calculate_neutron_yield(
        results['fusion_rate'], time_s
    )
    
    # Fusion energy
    results['fusion_energy'] = calculate_fusion_energy_output(
        results['neutron_yield']
    )
    
    # Bang time
    bang_time, bang_idx, peak_rate = find_bang_time(
        results['fusion_rate'], time_ns
    )
    results['bang_time'] = bang_time
    results['bang_idx'] = bang_idx
    results['peak_rate'] = peak_rate
    
    # Burn width
    results['burn_width'] = calculate_burn_width(
        results['fusion_rate'], time_ns
    )
    
    # Target gain (if laser energy provided)
    if laser_energy is not None:
        results['target_gain'] = calculate_target_gain(
            results['fusion_energy'], laser_energy
        )
    
    return results


def print_burn_summary(results: dict, stagnation_time: float = None):
    """
    Print formatted summary of burn analysis results
    
    Parameters
    ----------
    results : dict
        Results from comprehensive_burn_analysis()
    stagnation_time : float, optional
        Stagnation time [ns] for comparison with bang time
    """
    print("\n" + "="*60)
    print("BURN ANALYSIS SUMMARY")
    print("="*60)
    
    print(f"\nNeutron Yield:")
    print(f"  Total yield: {results['neutron_yield']:.2e} neutrons")
    
    print(f"\nFusion Energy:")
    print(f"  Total output: {results['fusion_energy']*1e-6:.3f} MJ")
    
    print(f"\nBang Time:")
    print(f"  Occurs at: {results['bang_time']:.2f} ns")
    print(f"  Peak rate: {results['peak_rate']:.2e} neutrons/s")
    
    if stagnation_time is not None:
        dt = results['bang_time'] - stagnation_time
        if dt > 0:
            print(f"  {dt:.2f} ns AFTER stagnation")
        else:
            print(f"  {-dt:.2f} ns BEFORE stagnation")
            print("  ⚠️  Bang before stagnation may indicate missing alpha transport")
    
    print(f"\nBurn Width:")
    print(f"  FWHM: {results['burn_width']:.2f} ns")
    print(f"       ({results['burn_width']*1000:.0f} ps)")
    
    if 'target_gain' in results:
        print(f"\nTarget Gain:")
        print(f"  G = {results['target_gain']:.3f}")
        if results['target_gain'] >= 1:
            print(f"  ✓ Net energy gain achieved!")
        else:
            print(f"  Below breakeven")
    
    print("="*60 + "\n")
