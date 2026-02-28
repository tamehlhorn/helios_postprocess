"""
Energetics: Energy calculations for ICF analysis

Functions for computing kinetic energy, hydro efficiency, PdV work, etc.
"""

import numpy as np
from typing import Optional, Dict, Tuple


def calculate_kinetic_energy(mass: np.ndarray, 
                             velocity: np.ndarray,
                             by_region: bool = False) -> np.ndarray:
    """
    Calculate total kinetic energy vs time
    
    KE(t) = Σ (1/2) * m_i * v_i²
    
    Parameters
    ----------
    mass : np.ndarray
        Zone masses [g], shape (n_times, n_zones)
    velocity : np.ndarray
        Zone velocities [cm/s], shape (n_times, n_zones)
    by_region : bool, optional
        If True, return dict with regional breakdown (not yet implemented)
        
    Returns
    -------
    np.ndarray
        Total kinetic energy [J], shape (n_times,)
        
    Examples
    --------
    >>> ke = calculate_kinetic_energy(mass, velocity)
    >>> peak_ke = ke.max()
    >>> t_peak = np.argmax(ke)
    
    Notes
    -----
    Conversion: 1 erg = 1e-7 J
    Typical ICF implosion: Peak KE ~ 1-10 kJ for direct drive
    """
    # Calculate KE per zone [erg]
    ke_zones = 0.5 * mass * velocity**2
    
    # Sum over all zones at each time [erg]
    ke_total = np.sum(ke_zones, axis=1)
    
    # Convert to Joules
    return ke_total * 1e-7


def calculate_hydro_efficiency(kinetic_energy: np.ndarray,
                               absorbed_laser_energy: np.ndarray) -> float:
    """
    Calculate hydrodynamic coupling efficiency
    
    η_hydro = max(KE) / E_absorbed_total
    
    Key metric for drive efficiency in ICF. Measures what fraction
    of absorbed laser energy is converted to implosion kinetic energy.
    
    Parameters
    ----------
    kinetic_energy : np.ndarray
        Kinetic energy vs time [J]
    absorbed_laser_energy : np.ndarray
        Absorbed laser energy vs time [J]
        
    Returns
    -------
    float
        Hydrodynamic efficiency (dimensionless, 0-1)
        
    Examples
    --------
    >>> efficiency = calculate_hydro_efficiency(ke, laser_energy)
    >>> print(f"Hydro efficiency: {efficiency:.1%}")
    
    Notes
    -----
    Typical values:
    - Direct drive: 5-7%
    - Indirect drive: 1-2%
    - Shock ignition: 3-5%
    
    Higher efficiency generally means better implosion performance,
    though other factors (symmetry, timing) are equally critical.
    """
    max_ke = np.max(kinetic_energy)
    total_absorbed = absorbed_laser_energy[-1]
    
    if total_absorbed > 0:
        return max_ke / total_absorbed
    else:
        return 0.0


def calculate_pdv_work(pressure: np.ndarray,
                       volume: np.ndarray, 
                       time: np.ndarray) -> np.ndarray:
    """
    Compute PdV work done by/on each zone
    
    W = ∫ P dV ≈ Σ P_i * ΔV_i
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure [dyne/cm²], shape (n_times, n_zones)
    volume : np.ndarray
        Zone volumes [cm³], shape (n_times, n_zones)
    time : np.ndarray
        Time [s], shape (n_times,)
        
    Returns
    -------
    np.ndarray
        PdV work [erg], shape (n_times, n_zones)
        
    Notes
    -----
    Positive work: Zone doing work on surroundings (expansion)
    Negative work: Work done on zone (compression)
    
    This will be fully implemented in Phase 3
    """
    # Placeholder - to be fully implemented
    dV = np.diff(volume, axis=0)
    P_avg = 0.5 * (pressure[:-1] + pressure[1:])
    
    work = P_avg * dV
    
    # Pad to match time dimension
    work_full = np.vstack([np.zeros_like(work[0]), work])
    
    return work_full


def find_peak_kinetic_energy_time(kinetic_energy: np.ndarray,
                                   time: np.ndarray) -> Tuple[float, int, float]:
    """
    Find time of peak kinetic energy
    
    Parameters
    ----------
    kinetic_energy : np.ndarray
        KE vs time [J]
    time : np.ndarray
        Time [ns]
        
    Returns
    -------
    peak_time : float
        Time of peak KE [ns]
    peak_idx : int
        Index of peak
    peak_ke : float
        Peak KE value [J]
        
    Examples
    --------
    >>> t_peak, idx, ke_peak = find_peak_kinetic_energy_time(ke, time)
    >>> print(f"Peak KE of {ke_peak*1e-6:.2f} MJ at {t_peak:.2f} ns")
    """
    peak_idx = np.argmax(kinetic_energy)
    peak_time = time[peak_idx]
    peak_ke = kinetic_energy[peak_idx]
    
    return peak_time, peak_idx, peak_ke


def calculate_implosion_velocity(zone_boundaries: np.ndarray,
                                 time: np.ndarray,
                                 surface_type: str = 'shell_inner') -> np.ndarray:
    """
    Calculate implosion velocity of a surface
    
    v = dr/dt
    
    Parameters
    ----------
    zone_boundaries : np.ndarray
        Radial positions [cm], shape (n_times, n_zones+1)
    time : np.ndarray
        Time [ns]
    surface_type : str
        'shell_inner', 'shell_outer', or 'interface'
        
    Returns
    -------
    np.ndarray
        Velocity [cm/s], shape (n_times-1,)
        
    Notes
    -----
    To be fully implemented with surface tracking in Phase 3
    """
    # Placeholder
    # For now, track innermost boundary (approximation for shell inner surface)
    radius = zone_boundaries[:, 0]
    dt = np.diff(time) * 1e-9  # Convert ns to s
    dr = np.diff(radius)
    
    velocity = dr / dt  # cm/s
    
    return velocity
