"""
Energetics: Energy calculations for ICF analysis

Functions for computing kinetic energy, hydro efficiency, PdV work, etc.
"""

import numpy as np
from typing import Optional, Dict, Tuple
import logging

logger = logging.getLogger(__name__)


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


# ============================================================================
# Radiation drive energetics
# ============================================================================

# Stefan-Boltzmann constant in W cm^-2 K^-4
_SIGMA_SB_W_CM2_K4 = 5.6704e-12
# eV -> K conversion (1 eV = 11604.518 K, exact within machine precision)
_EV_TO_K = 11604.518


def compute_radiation_drive_energy(data) -> Optional[Dict]:
    """
    Integrate the radiation drive energy delivered to the simulation boundary.

    Helios applies a time-dependent blackbody drive at Rmin or Rmax based on
    the [Rad Source Data] block of the .rhw file (parsed by RHWParser). The
    incident flux is

        F(t) = flux_multiplier × σ_SB × T_rad(t)^4         [W/cm^2]

    and the total energy delivered to the boundary is

        E_rad = ∫ F(t) × 4π R_drive(t)^2 dt                [J]

    where R_drive(t) is the inner (Rmin) or outer (Rmax) zone-boundary
    radius at the time of the drive, interpolated onto the drive-time grid
    from the simulation's zone_boundaries history.

    Parameters
    ----------
    data : ICFRunData
        Must have populated `drive_time`, `drive_temperature`, `time`,
        `zone_boundaries`, and `rhw_config` (with `drive_location` and
        `drive_flux_multiplier` set by the RHW parser).

    Returns
    -------
    dict or None
        None if no radiation drive is present. Otherwise, a dict with:
            E_rad_J              : total radiation energy at boundary [J]
            peak_flux_Wcm2       : peak σT⁴ × flux_mult              [W/cm^2]
            peak_T_eV            : max T_rad                          [eV]
            peak_T_time_s        : time of peak T_rad                 [s]
            fluence_Jcm2         : time-integrated flux               [J/cm^2]
            R_outer_at_peak_cm   : R_drive at t_peak_T                [cm]
            location             : "Rmin" or "Rmax"
            flux_multiplier      : from RHW
            n_points             : drive-table row count
            drive_start_s        : first drive_time entry             [s]
            drive_end_s          : last drive_time entry              [s]
    """
    if data.drive_temperature is None or data.drive_time is None:
        return None
    if len(data.drive_time) == 0:
        return None

    rhw = getattr(data, "rhw_config", None)
    flux_mult = float(getattr(rhw, "drive_flux_multiplier", 1.0)) if rhw else 1.0
    location = str(getattr(rhw, "drive_location", "")) if rhw else ""

    if location not in ("Rmin", "Rmax"):
        # Unknown boundary — can't pick the right zone-boundary column
        logger.warning("Radiation drive present but location is not "
                       "'Rmin' or 'Rmax'; cannot compute energy")
        return None

    if data.zone_boundaries is None or data.time is None:
        logger.warning("Cannot compute radiation drive energy without "
                       "zone_boundaries and time arrays")
        return None

    T_eV = np.asarray(data.drive_temperature, dtype=float)
    t_drive = np.asarray(data.drive_time, dtype=float)

    # Flux
    T_K = T_eV * _EV_TO_K
    flux = flux_mult * _SIGMA_SB_W_CM2_K4 * T_K**4    # W/cm^2

    # R_drive at each drive sample, interpolated from sim zone boundaries.
    # zone_boundaries shape: (n_t, n_zones+1); column [-1] = outer, [0] = inner.
    R_t = (data.zone_boundaries[:, -1]
           if location == "Rmax"
           else data.zone_boundaries[:, 0])    # cm
    R_drive = np.interp(t_drive, data.time, R_t)

    # Fluence and total energy
    fluence_Jcm2 = float(np.trapezoid(flux, t_drive))
    dE_dt_W = flux * 4.0 * np.pi * R_drive**2          # W
    E_rad_J = float(np.trapezoid(dE_dt_W, t_drive))

    peak_idx = int(np.argmax(T_eV))
    peak_T_eV = float(T_eV[peak_idx])
    peak_T_time_s = float(t_drive[peak_idx])
    R_at_peak_cm = float(np.interp(peak_T_time_s, data.time, R_t))

    return {
        "E_rad_J":            E_rad_J,
        "peak_flux_Wcm2":     float(flux.max()),
        "peak_T_eV":          peak_T_eV,
        "peak_T_time_s":      peak_T_time_s,
        "fluence_Jcm2":       fluence_Jcm2,
        "R_outer_at_peak_cm": R_at_peak_cm,
        "location":           location,
        "flux_multiplier":    flux_mult,
        "n_points":           int(len(t_drive)),
        "drive_start_s":      float(t_drive[0]),
        "drive_end_s":        float(t_drive[-1]),
    }
