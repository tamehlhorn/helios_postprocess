"""
Pressure Gradient Analysis for ICF Simulations

Pressure gradients are critical for:
- Shock identification and characterization
- Rayleigh-Taylor instability assessment
- Isobaric core identification
- Drive uniformity analysis

Key metrics:
- dP/dr: Radial pressure gradient [J/cm⁴]
- L_P = P / |dP/dr|: Pressure scale length [cm]
- (1/P) dP/dr: Normalized gradient [1/cm]
- Isobaric core: Region where |dP/dr| ≈ 0

Physical significance:
- Large dP/dr → Shock location
- Small L_P → Steep gradient (unstable)
- Large L_P → Gentle gradient (stable)
- Isobaric core → Good hot spot formation

Units:
- Pressure: J/cm³
- Radius: cm
- dP/dr: J/cm⁴
- L_P: cm

Author: Prof T
"""

import numpy as np
from typing import Optional, Tuple, Dict, List
from scipy.ndimage import gaussian_filter1d


def calculate_pressure_gradient(
    pressure: np.ndarray,
    zone_boundaries: np.ndarray,
    smoothing_sigma: float = 0.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate radial pressure gradient dP/dr
    
    Uses centered finite differences for interior points,
    forward/backward differences at boundaries.
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    smoothing_sigma : float, default=0.0
        Gaussian smoothing width (zones). 0 = no smoothing.
        Recommended: 1-3 for noisy data
        
    Returns
    -------
    dP_dr : np.ndarray
        Pressure gradient [J/cm⁴], shape (n_zones,)
    zone_centers : np.ndarray
        Zone center radii [cm], shape (n_zones,)
        
    Examples
    --------
    >>> dP_dr, r = calculate_pressure_gradient(
    ...     run.pressure[idx],
    ...     run.zone_boundaries[idx],
    ...     smoothing_sigma=2.0
    ... )
    >>> 
    >>> # Find shocks (large dP/dr)
    >>> shock_mask = np.abs(dP_dr) > 1e10  # Threshold depends on target
    >>> print(f"Shock at r = {r[shock_mask]} cm")
    
    Notes
    -----
    - Positive dP/dr: Pressure increases outward (compression wave)
    - Negative dP/dr: Pressure decreases outward (rarefaction)
    - Large |dP/dr|: Shock location
    """
    # Zone centers
    zone_centers = 0.5 * (zone_boundaries[:-1] + zone_boundaries[1:])
    
    # Optional smoothing
    if smoothing_sigma > 0:
        pressure_smooth = gaussian_filter1d(pressure, smoothing_sigma)
    else:
        pressure_smooth = pressure
    
    # Calculate gradient using centered differences
    dP_dr = np.zeros_like(pressure)
    
    # Interior points: centered difference
    for i in range(1, len(pressure) - 1):
        dr = zone_centers[i+1] - zone_centers[i-1]
        dP = pressure_smooth[i+1] - pressure_smooth[i-1]
        dP_dr[i] = dP / dr
    
    # Boundary points: forward/backward difference
    # First point: forward difference
    dr = zone_centers[1] - zone_centers[0]
    dP = pressure_smooth[1] - pressure_smooth[0]
    dP_dr[0] = dP / dr
    
    # Last point: backward difference
    dr = zone_centers[-1] - zone_centers[-2]
    dP = pressure_smooth[-1] - pressure_smooth[-2]
    dP_dr[-1] = dP / dr
    
    return dP_dr, zone_centers


def pressure_scale_length(
    pressure: np.ndarray,
    zone_boundaries: np.ndarray,
    smoothing_sigma: float = 0.0,
    floor: float = 1e-10
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate pressure scale length L_P = P / |dP/dr|
    
    The scale length represents the distance over which pressure
    changes by a factor of e. Large L_P indicates gentle gradients.
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    smoothing_sigma : float, default=0.0
        Gaussian smoothing width (zones)
    floor : float, default=1e-10
        Minimum |dP/dr| to avoid division by zero [J/cm⁴]
        
    Returns
    -------
    L_P : np.ndarray
        Pressure scale length [cm], shape (n_zones,)
    zone_centers : np.ndarray
        Zone center radii [cm], shape (n_zones,)
        
    Examples
    --------
    >>> L_P, r = pressure_scale_length(
    ...     run.pressure[idx],
    ...     run.zone_boundaries[idx],
    ...     smoothing_sigma=2.0
    ... )
    >>> 
    >>> # Identify isobaric core (large L_P)
    >>> isobaric_mask = L_P > 0.001  # 10 μm
    >>> print(f"Isobaric core radius: {r[isobaric_mask].max()*1e4:.1f} μm")
    
    Notes
    -----
    - Large L_P (> 100 μm): Gentle gradient, stable
    - Small L_P (< 10 μm): Steep gradient, unstable
    - Very large L_P: Nearly isobaric (flat pressure)
    """
    dP_dr, zone_centers = calculate_pressure_gradient(
        pressure, zone_boundaries, smoothing_sigma
    )
    
    # Avoid division by zero
    dP_dr_safe = np.where(np.abs(dP_dr) < floor, floor, dP_dr)
    
    # L_P = P / |dP/dr|
    L_P = np.abs(pressure / dP_dr_safe)
    
    return L_P, zone_centers


def normalized_pressure_gradient(
    pressure: np.ndarray,
    zone_boundaries: np.ndarray,
    smoothing_sigma: float = 0.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate normalized pressure gradient (1/P) dP/dr
    
    Dimensionless gradient useful for comparing different simulations.
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    smoothing_sigma : float, default=0.0
        Gaussian smoothing width (zones)
        
    Returns
    -------
    normalized_gradient : np.ndarray
        (1/P) dP/dr [1/cm], shape (n_zones,)
    zone_centers : np.ndarray
        Zone center radii [cm], shape (n_zones,)
        
    Examples
    --------
    >>> grad_norm, r = normalized_pressure_gradient(
    ...     run.pressure[idx],
    ...     run.zone_boundaries[idx]
    ... )
    >>> 
    >>> # Plot normalized gradient
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(r*1e4, grad_norm)
    >>> plt.xlabel('Radius [μm]')
    >>> plt.ylabel('(1/P) dP/dr [1/cm]')
    >>> plt.show()
    """
    dP_dr, zone_centers = calculate_pressure_gradient(
        pressure, zone_boundaries, smoothing_sigma
    )
    
    # Avoid division by zero
    pressure_safe = np.where(pressure < 1e-10, 1e-10, pressure)
    
    normalized_gradient = dP_dr / pressure_safe
    
    return normalized_gradient, zone_centers


def identify_isobaric_core(
    pressure: np.ndarray,
    zone_boundaries: np.ndarray,
    L_P_threshold: float = 1e-3,
    smoothing_sigma: float = 2.0
) -> Tuple[np.ndarray, float]:
    """
    Identify isobaric core region (nearly flat pressure)
    
    The isobaric core is the central region where pressure is
    approximately constant. Important indicator of hot spot quality.
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    L_P_threshold : float, default=1e-3
        Minimum scale length for isobaric region [cm]
        Default: 10 μm
    smoothing_sigma : float, default=2.0
        Gaussian smoothing width (zones)
        
    Returns
    -------
    isobaric_mask : np.ndarray
        Boolean mask for isobaric zones
    r_isobaric : float
        Outer radius of isobaric core [cm]
        
    Examples
    --------
    >>> isobaric_mask, r_iso = identify_isobaric_core(
    ...     run.pressure[idx],
    ...     run.zone_boundaries[idx],
    ...     L_P_threshold=1e-3  # 10 μm
    ... )
    >>> 
    >>> print(f"Isobaric core radius: {r_iso*1e4:.1f} μm")
    >>> print(f"Isobaric core has {np.sum(isobaric_mask)} zones")
    
    Notes
    -----
    - Large isobaric core: Good hot spot formation
    - Small/no isobaric core: Poor hot spot quality
    - Typical size: 10-50 μm for ICF
    """
    # Calculate pressure scale length
    L_P, zone_centers = pressure_scale_length(
        pressure, zone_boundaries, smoothing_sigma
    )
    
    # Identify isobaric region
    isobaric_mask = L_P > L_P_threshold
    
    # Find outer radius of isobaric core
    if np.any(isobaric_mask):
        isobaric_indices = np.where(isobaric_mask)[0]
        r_isobaric = zone_centers[isobaric_indices[-1]]
    else:
        r_isobaric = 0.0
    
    return isobaric_mask, r_isobaric


def identify_shocks(
    pressure: np.ndarray,
    zone_boundaries: np.ndarray,
    dP_dr_threshold: float = 1e10,
    smoothing_sigma: float = 1.0,
    min_separation: float = 5e-4
) -> List[Dict[str, float]]:
    """
    Identify shock locations from pressure gradient
    
    Shocks appear as large pressure gradients (jumps).
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    dP_dr_threshold : float, default=1e10
        Minimum |dP/dr| to identify shock [J/cm⁴]
    smoothing_sigma : float, default=1.0
        Gaussian smoothing width (zones)
    min_separation : float, default=5e-4
        Minimum separation between shocks [cm] = 5 μm
        
    Returns
    -------
    list of dict
        Each dict contains:
        - 'radius': Shock radius [cm]
        - 'dP_dr': Pressure gradient at shock [J/cm⁴]
        - 'P_jump': Pressure jump across shock [J/cm³]
        - 'P_ratio': Post-shock / pre-shock pressure ratio
        
    Examples
    --------
    >>> shocks = identify_shocks(
    ...     run.pressure[idx],
    ...     run.zone_boundaries[idx],
    ...     dP_dr_threshold=1e10
    ... )
    >>> 
    >>> for i, shock in enumerate(shocks):
    ...     print(f"Shock {i+1}:")
    ...     print(f"  Radius: {shock['radius']*1e4:.1f} μm")
    ...     print(f"  P jump: {shock['P_jump']/1e6:.1f} Gbar")
    ...     print(f"  P ratio: {shock['P_ratio']:.2f}")
    
    Notes
    -----
    - Primary shock: First (innermost) shock
    - Secondary shock: Reshock from shell
    - Tertiary shock: Additional shocks from drive
    """
    dP_dr, zone_centers = calculate_pressure_gradient(
        pressure, zone_boundaries, smoothing_sigma
    )
    
    # Find zones with large gradient
    shock_candidates = np.abs(dP_dr) > dP_dr_threshold
    
    if not np.any(shock_candidates):
        return []
    
    # Get indices of shock candidates
    shock_indices = np.where(shock_candidates)[0]
    
    # Group nearby shocks
    shocks = []
    current_group = [shock_indices[0]]
    
    for idx in shock_indices[1:]:
        if zone_centers[idx] - zone_centers[current_group[-1]] < min_separation:
            current_group.append(idx)
        else:
            # Process current group
            shocks.append(_process_shock_group(
                current_group, pressure, dP_dr, zone_centers
            ))
            current_group = [idx]
    
    # Process last group
    if current_group:
        shocks.append(_process_shock_group(
            current_group, pressure, dP_dr, zone_centers
        ))
    
    return shocks


def _process_shock_group(
    indices: List[int],
    pressure: np.ndarray,
    dP_dr: np.ndarray,
    zone_centers: np.ndarray
) -> Dict[str, float]:
    """
    Process a group of shock candidate zones
    
    Returns shock properties at peak gradient location
    """
    # Find peak gradient in group
    peak_idx = indices[np.argmax(np.abs(dP_dr[indices]))]
    
    # Estimate pressure jump
    # Look at pressure before and after shock
    if peak_idx > 0 and peak_idx < len(pressure) - 1:
        P_pre = pressure[peak_idx - 1]
        P_post = pressure[peak_idx + 1]
        P_jump = P_post - P_pre
        P_ratio = P_post / P_pre if P_pre > 1e-10 else 0.0
    else:
        P_jump = 0.0
        P_ratio = 1.0
    
    return {
        'radius': zone_centers[peak_idx],
        'dP_dr': dP_dr[peak_idx],
        'P_jump': P_jump,
        'P_ratio': P_ratio
    }


def rayleigh_taylor_criterion(
    pressure: np.ndarray,
    density: np.ndarray,
    zone_boundaries: np.ndarray,
    smoothing_sigma: float = 2.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate Rayleigh-Taylor instability criterion
    
    RT instability occurs when: (∇P · ∇ρ) < 0
    I.e., pressure gradient opposes density gradient
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    smoothing_sigma : float, default=2.0
        Gaussian smoothing width (zones)
        
    Returns
    -------
    RT_unstable : np.ndarray
        Boolean mask: True where RT unstable
    zone_centers : np.ndarray
        Zone center radii [cm], shape (n_zones,)
        
    Examples
    --------
    >>> RT_mask, r = rayleigh_taylor_criterion(
    ...     run.pressure[idx],
    ...     run.density[idx],
    ...     run.zone_boundaries[idx]
    ... )
    >>> 
    >>> if np.any(RT_mask):
    ...     print("RT unstable regions:")
    ...     print(f"  Radii: {r[RT_mask]*1e4} μm")
    ... else:
    ...     print("No RT instability")
    
    Notes
    -----
    - RT instability grows fastest at ablation front
    - Can seed mix, reducing hot spot performance
    - Mitigation: smoother drive, better symmetry
    """
    # Calculate gradients
    dP_dr, zone_centers = calculate_pressure_gradient(
        pressure, zone_boundaries, smoothing_sigma
    )
    
    drho_dr, _ = calculate_pressure_gradient(
        density, zone_boundaries, smoothing_sigma
    )
    
    # RT unstable when gradients oppose
    # (dP/dr > 0 and drho/dr < 0) OR (dP/dr < 0 and drho/dr > 0)
    RT_unstable = (dP_dr * drho_dr) < 0
    
    return RT_unstable, zone_centers


def pressure_gradient_evolution(
    pressure: np.ndarray,
    zone_boundaries: np.ndarray,
    time: np.ndarray,
    smoothing_sigma: float = 2.0
) -> Dict[str, np.ndarray]:
    """
    Track pressure gradient metrics vs time
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_times, n_zones)
    zone_boundaries : np.ndarray
        Zone boundaries [cm], shape (n_times, n_zones+1)
    time : np.ndarray
        Time array [ns], shape (n_times,)
    smoothing_sigma : float, default=2.0
        Gaussian smoothing width
        
    Returns
    -------
    dict
        Dictionary with time series:
        - 'time': Time array [ns]
        - 'max_dP_dr': Maximum |dP/dr| [J/cm⁴]
        - 'shock_radius': Primary shock radius [cm]
        - 'isobaric_radius': Isobaric core radius [cm]
        
    Examples
    --------
    >>> evolution = pressure_gradient_evolution(
    ...     run.pressure,
    ...     run.zone_boundaries,
    ...     run.time,
    ...     smoothing_sigma=2.0
    ... )
    >>> 
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(evolution['time'], evolution['shock_radius']*1e4)
    >>> plt.xlabel('Time [ns]')
    >>> plt.ylabel('Shock Radius [μm]')
    >>> plt.show()
    """
    n_times = len(time)
    
    results = {
        'time': time,
        'max_dP_dr': np.zeros(n_times),
        'shock_radius': np.zeros(n_times),
        'isobaric_radius': np.zeros(n_times)
    }
    
    for i in range(n_times):
        # Maximum gradient
        dP_dr, zone_centers = calculate_pressure_gradient(
            pressure[i], zone_boundaries[i], smoothing_sigma
        )
        results['max_dP_dr'][i] = np.abs(dP_dr).max()
        
        # Shock radius (location of max gradient)
        shock_idx = np.argmax(np.abs(dP_dr))
        results['shock_radius'][i] = zone_centers[shock_idx]
        
        # Isobaric core radius
        _, r_iso = identify_isobaric_core(
            pressure[i], zone_boundaries[i],
            L_P_threshold=1e-3, smoothing_sigma=smoothing_sigma
        )
        results['isobaric_radius'][i] = r_iso
    
    return results


def comprehensive_pressure_gradient_analysis(
    pressure: np.ndarray,
    density: np.ndarray,
    zone_boundaries: np.ndarray,
    smoothing_sigma: float = 2.0
) -> Dict:
    """
    Comprehensive pressure gradient analysis at single time
    
    Calculates all gradient metrics for complete assessment.
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundaries [cm], shape (n_zones+1,)
    smoothing_sigma : float, default=2.0
        Gaussian smoothing width
        
    Returns
    -------
    dict
        Dictionary with:
        - 'dP_dr': Pressure gradient array [J/cm⁴]
        - 'L_P': Pressure scale length array [cm]
        - 'zone_centers': Radii [cm]
        - 'max_dP_dr': Maximum |dP/dr| [J/cm⁴]
        - 'isobaric_radius': Isobaric core radius [cm]
        - 'shocks': List of shock properties
        - 'RT_unstable_fraction': Fraction of mass RT unstable
        
    Examples
    --------
    >>> results = comprehensive_pressure_gradient_analysis(
    ...     run.pressure[idx],
    ...     run.density[idx],
    ...     run.zone_boundaries[idx],
    ...     smoothing_sigma=2.0
    ... )
    >>> 
    >>> print(f"Isobaric core: {results['isobaric_radius']*1e4:.1f} μm")
    >>> print(f"Number of shocks: {len(results['shocks'])}")
    >>> for i, shock in enumerate(results['shocks']):
    ...     print(f"  Shock {i+1} at {shock['radius']*1e4:.1f} μm")
    """
    results = {}
    
    # Basic gradients
    dP_dr, zone_centers = calculate_pressure_gradient(
        pressure, zone_boundaries, smoothing_sigma
    )
    results['dP_dr'] = dP_dr
    results['zone_centers'] = zone_centers
    results['max_dP_dr'] = np.abs(dP_dr).max()
    
    # Pressure scale length
    L_P, _ = pressure_scale_length(
        pressure, zone_boundaries, smoothing_sigma
    )
    results['L_P'] = L_P
    
    # Isobaric core
    isobaric_mask, r_iso = identify_isobaric_core(
        pressure, zone_boundaries,
        L_P_threshold=1e-3, smoothing_sigma=smoothing_sigma
    )
    results['isobaric_radius'] = r_iso
    
    # Shocks
    shocks = identify_shocks(
        pressure, zone_boundaries,
        dP_dr_threshold=1e10, smoothing_sigma=smoothing_sigma
    )
    results['shocks'] = shocks
    
    # Rayleigh-Taylor stability
    RT_unstable, _ = rayleigh_taylor_criterion(
        pressure, density, zone_boundaries, smoothing_sigma
    )
    results['RT_unstable_fraction'] = np.sum(RT_unstable) / len(RT_unstable)
    
    return results


def analyze_first_shock(
    pressure: np.ndarray,
    zone_boundaries: np.ndarray,
    density: np.ndarray,
    time: np.ndarray,
    ablator_outer_zone: int,
    fuel_inner_zone: int,
    search_inner_zone: int = 0,
    gamma: float = 5.0 / 3.0,
    dP_dr_threshold: float = 1e10,
    smoothing_sigma: float = 1.5,
    min_time_ns: Optional[float] = None,
    min_P_ratio: float = 1.5,
    lookback_steps: int = 5,
) -> Dict[str, np.ndarray]:
    """
    Track first (outermost) shock through the shell and compute strength metrics.

    The first shock is the outermost pressure discontinuity in the ablator at
    each timestep. It is identified as the shock with the largest radius among
    all shocks found by identify_shocks(). Shock velocity is computed by
    differentiating the shock radius history. Rankine-Hugoniot relations for a
    gamma-law gas give Mach number and density jump from the pressure ratio.

    Parameters
    ----------
    pressure : np.ndarray
        Total pressure [J/cm³], shape (n_times, n_zones).
        Use ion_pressure + rad_pressure from ICFRunData.
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_times, n_zones+1)
    density : np.ndarray
        Mass density [g/cm³], shape (n_times, n_zones)
    time : np.ndarray
        Time array [ns], shape (n_times,)
    ablator_outer_zone : int
        Zone index of outer ablator boundary (outermost zone to search for shock).
        Typically ri[0, -2] - 1 (fuel/ablator interface node index minus 1).
    fuel_inner_zone : int
        Zone index of fuel/ablator interface (inner ablator boundary).
        Used to detect shock breakout into fuel. Typically ri[0, -2].
    gamma : float, default=5/3
        Adiabatic index for Rankine-Hugoniot relations. Use 5/3 for DT.
    dP_dr_threshold : float, default=1e10
        Minimum |dP/dr| to identify a shock [J/cm⁴].
    smoothing_sigma : float, default=1.5
        Gaussian smoothing width [zones] for pressure gradient.

    Returns
    -------
    dict with keys:
        'time' : np.ndarray [ns]
            Time array (same as input).
        'shock_radius' : np.ndarray [cm]
            Outermost shock radius at each timestep. NaN if no shock found.
        'shock_velocity' : np.ndarray [cm/ns]
            Shock propagation velocity (d r_shock / dt). NaN where undefined.
        'P_ratio' : np.ndarray
            Post-shock / pre-shock pressure ratio (compression ratio proxy).
        'P_jump_Gbar' : np.ndarray
            Pressure jump across shock [Gbar].
        'mach_number' : np.ndarray
            Shock Mach number from Rankine-Hugoniot (gamma-law).
        'density_jump' : np.ndarray
            Post/pre density jump ratio from Rankine-Hugoniot.
        'shock_pressure_Gbar' : np.ndarray
            Absolute post-shock pressure [Gbar].
        'breakout_time_ns' : float
            Time when shock first crosses fuel_inner_zone [ns].
            NaN if breakout not observed in simulation window.
        'breakout_pressure_Gbar' : float
            Post-shock pressure at breakout time [Gbar].
        'breakout_mach' : float
            Mach number at breakout time.

    Notes
    -----
    Rankine-Hugoniot for gamma-law gas:
        M^2 = [(gamma+1)*P_ratio + (gamma-1)] / (2*gamma)
        rho2/rho1 = (gamma+1)*M^2 / [(gamma-1)*M^2 + 2]

    For strong shocks (M >> 1):
        P_ratio -> (2*gamma*M^2) / (gamma+1)
        rho2/rho1 -> (gamma+1) / (gamma-1) = 4 for gamma=5/3

    Shock velocity is computed using central differences on the shock radius
    history, with forward/backward differences at the endpoints. Time is in ns
    so velocity is in cm/ns = 10^7 cm/s = 100 km/s.

    Examples
    --------
    >>> from helios_postprocess.pressure_gradients import analyze_first_shock
    >>> result = analyze_first_shock(
    ...     pressure=data.ion_pressure + data.rad_pressure,
    ...     zone_boundaries=data.zone_boundaries,
    ...     density=data.density,
    ...     time=data.time,
    ...     ablator_outer_zone=data.region_interfaces_indices[0, -1] - 1,
    ...     fuel_inner_zone=data.region_interfaces_indices[0, -2],
    ... )
    >>> print(f"Breakout time:     {result['breakout_time_ns']:.3f} ns")
    >>> print(f"Breakout pressure: {result['breakout_pressure_Gbar']:.1f} Gbar")
    >>> print(f"Breakout Mach:     {result['breakout_mach']:.1f}")
    """
    n_times = len(time)

    shock_radius   = np.full(n_times, np.nan)
    P_ratio        = np.full(n_times, np.nan)
    P_jump_Gbar    = np.full(n_times, np.nan)
    mach_number    = np.full(n_times, np.nan)
    density_jump   = np.full(n_times, np.nan)
    shock_pressure_Gbar = np.full(n_times, np.nan)

    # Zone center radii at each timestep
    zone_centers_all = 0.5 * (zone_boundaries[:, :-1] + zone_boundaries[:, 1:])

    for i in range(n_times):
        # Search only within the ablator (zones up to ablator_outer_zone)
        p_slice  = pressure[i, search_inner_zone:ablator_outer_zone + 1]
        zb_slice = zone_boundaries[i, search_inner_zone:ablator_outer_zone + 2]

        shocks = identify_shocks(
            p_slice, zb_slice,
            dP_dr_threshold=dP_dr_threshold,
            smoothing_sigma=smoothing_sigma,
        )

        if not shocks:
            continue

        # First shock = outermost (largest radius)
        first = max(shocks, key=lambda s: s['radius'])

        shock_radius[i] = first['radius']
        pr = first['P_ratio']
        P_ratio[i] = pr
        P_jump_Gbar[i] = first['P_jump'] * 1e-8  # J/cm³ -> Gbar

        # Post-shock absolute pressure [Gbar]
        # find zone nearest shock radius
        zc = zone_centers_all[i, :ablator_outer_zone + 1]
        nearest = int(np.argmin(np.abs(zc - first['radius'])))
        shock_pressure_Gbar[i] = pressure[i, nearest] * 1e-8

        # Rankine-Hugoniot
        pr = max(pr, 1.0 / pr) if pr > 0 else pr  # inward shock: compressed side at lower radius
        P_ratio[i] = pr
        if pr > 1.0:
            m2 = ((gamma + 1.0) * pr + (gamma - 1.0)) / (2.0 * gamma)
            if m2 > 0:
                mach_number[i] = np.sqrt(m2)
                density_jump[i] = ((gamma + 1.0) * m2) / ((gamma - 1.0) * m2 + 2.0)

    # Shock velocity: d(r_shock)/dt using np.gradient (central differences)
    # Only where shock_radius is valid
    valid = ~np.isnan(shock_radius)
    shock_velocity = np.full(n_times, np.nan)
    if np.sum(valid) >= 2:
        t_valid = time[valid]
        r_valid = shock_radius[valid]
        v_valid = np.gradient(r_valid, t_valid)  # cm/ns
        shock_velocity[valid] = v_valid

    # Breakout detection: first shock crosses the gas-cavity interface (node fuel_inner_zone).
    # ICF definition: the foot-launched shock traversing ablator + foams + ice/wetted-foam
    # breaks out into the central gas cavity, beginning cavity compression and heating.
    # This breakout time sets the base adiabat for the cold shell.
    #
    # Two gates protect against spurious detections at t=0:
    #   (1) min_time_ns - skip initial-condition discontinuities; defaults to time[0]
    #       + 1e-3 ns so the very first step is always excluded.
    #   (2) min_P_ratio - a real foot shock has P_post/P_pre >> 1 (typically 5-50);
    #       static IC steps have P_ratio ~ 1.
    breakout_time_ns       = np.nan
    breakout_pressure_Gbar = np.nan
    breakout_mach          = np.nan

    t_floor = min_time_ns if min_time_ns is not None else (time[0] + 1e-3)

    # Breakout criterion: first timestep where the tracked shock has crossed the
    # gas-cavity interface (r_sh <= r_interface), AND a real shock was present
    # in the shell within the preceding `lookback_steps` timesteps.
    #
    # Look-back on P_ratio is required because at the moment of crossing the
    # transmitted wave into the low-density gas cavity is weak (P_ratio ~ 1)
    # due to impedance mismatch — the *incoming* shock strength (just before
    # crossing) is the relevant diagnostic. This matches the ICF convention
    # and handles weak-first-shock designs (e.g., CH foam ablators) where the
    # peak P_ratio in the shell is modest (~2-4), not the 5-50 of HDC-style
    # strong-first-shock designs.
    for i in range(n_times):
        if time[i] < t_floor:
            continue
        if np.isnan(shock_radius[i]):
            continue
        interface_radius = zone_boundaries[i, fuel_inner_zone]
        if shock_radius[i] > interface_radius:
            continue
        # Crossing candidate: check that a real shock was present just before.
        lo = max(0, i - lookback_steps)
        pr_window = P_ratio[lo:i + 1]
        pr_valid = pr_window[~np.isnan(pr_window)]
        if pr_valid.size == 0 or np.max(pr_valid) < min_P_ratio:
            continue
        breakout_time_ns       = time[i]
        breakout_pressure_Gbar = shock_pressure_Gbar[i]
        breakout_mach          = mach_number[i] if not np.isnan(mach_number[i]) else np.nan
        break

    return {
        'time':                  time,
        'shock_radius':          shock_radius,
        'shock_velocity':        shock_velocity,
        'P_ratio':               P_ratio,
        'P_jump_Gbar':           P_jump_Gbar,
        'mach_number':           mach_number,
        'density_jump':          density_jump,
        'shock_pressure_Gbar':   shock_pressure_Gbar,
        'breakout_time_ns':      breakout_time_ns,
        'breakout_pressure_Gbar': breakout_pressure_Gbar,
        'breakout_mach':         breakout_mach,
    }


def track_shock_trajectories(
    pressure: np.ndarray,
    zone_boundaries: np.ndarray,
    time: np.ndarray,
    interface_radius: np.ndarray,
    search_inner_zone: int = 0,
    search_outer_zone: Optional[int] = None,
    dP_dr_threshold: float = 1e9,
    smoothing_sigma: float = 1.5,
    min_P_ratio: float = 1.5,
    max_shock_velocity: float = 0.05,
    min_separation: float = 5e-4,
    group_separation: Optional[float] = None,
    breakout_tolerance: float = 0.0,
    max_gap_steps: int = 0,
    min_time_ns: Optional[float] = None,
    lookback_steps: int = 5,
) -> Dict:
    """
    Multi-shock tracker: link per-snapshot identify_shocks() detections into
    trajectories via greedy nearest-radius matching with inward-only motion
    and a velocity bound.

    Parameters
    ----------
    pressure : np.ndarray
        Total pressure [J/cm³], shape (n_times, n_zones).
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_times, n_zones+1).
    time : np.ndarray
        Time array [ns], shape (n_times,).
    interface_radius : np.ndarray
        Gas/ice interface radius [cm], shape (n_times,). Used for breakout
        detection (e.g. zone_boundaries[t, ri[t, 0]] for capsules).
    search_inner_zone, search_outer_zone : int
        Inclusive zone-index window in which identify_shocks is run.
        Default outer zone is n_zones - 1.
    dP_dr_threshold, smoothing_sigma, min_separation : floats
        Passed through to identify_shocks. min_separation is also the
        coalescence-detection threshold.
    min_P_ratio : float, default 1.5
        Minimum compression ratio (max(pr, 1/pr)) for a detection to seed or
        extend a trajectory.
    max_shock_velocity : float, default 0.05
        Magnitude bound on inward shock motion [cm/ns] (~500 km/s).
    min_time_ns : float or None
        Earliest time to start tracking. Default = time[0] + 1e-3 (skips IC).
    lookback_steps : int, default 5
        Window over which a trajectory must have shown P_ratio >= min_P_ratio
        for a radius crossing to count as breakout.

    Returns
    -------
    dict with keys
        'trajectories'        : list of per-trajectory dicts with arrays
                                (indices, radius, time_ns, P_ratio,
                                P_post_Gbar) plus scalars (started, ended,
                                reason_ended).
        'coalescence_events'  : list of {t_idx, time_ns, radius,
                                trajectory_pair} where two active
                                trajectories converged within
                                min_separation.
        'breakouts'           : list of per-trajectory breakout records
                                (trajectory_id, t_idx, time_ns, radius,
                                P_post_Gbar, P_ratio), sorted by time.
        'time', 'interface_radius' : inputs echoed for plotting convenience.

    Notes
    -----
    identify_shocks reports P_ratio = P[peak+1]/P[peak-1], which is < 1 for
    inward shocks (compressed side is inner). The tracker flips this to a
    compression ratio via max(pr, 1/pr) before applying min_P_ratio.

    Greedy matching: each active trajectory (iterated in deterministic order)
    claims the closest unmatched candidate at smaller radius whose implied
    velocity is within the bound. Unmatched trajectories deactivate;
    unmatched candidates seed new trajectories.

    Coalescence is post-hoc -- a pair of trajectories alive at the same
    timestep whose radii are within min_separation triggers one event;
    the younger trajectory is tagged 'coalesced_into:<id>'.
    """
    n_times = len(time)
    n_zones = pressure.shape[1]
    if search_outer_zone is None:
        search_outer_zone = n_zones - 1
    if group_separation is None:
        group_separation = min_separation

    t_floor = min_time_ns if min_time_ns is not None else (time[0] + 1e-3)

    # ---- 1. Per-timestep shock detection (above min_P_ratio compression) ----
    detections: List[List[Dict]] = [[] for _ in range(n_times)]
    for i in range(n_times):
        if time[i] < t_floor:
            continue
        p_slice  = pressure[i, search_inner_zone:search_outer_zone + 1]
        zb_slice = zone_boundaries[i, search_inner_zone:search_outer_zone + 2]
        if p_slice.size < 5:
            continue
        shocks = identify_shocks(
            p_slice, zb_slice,
            dP_dr_threshold=dP_dr_threshold,
            smoothing_sigma=smoothing_sigma,
            min_separation=group_separation,
        )
        if not shocks:
            continue
        zc = 0.5 * (zb_slice[:-1] + zb_slice[1:])
        for sh in shocks:
            pr_raw = sh['P_ratio']
            if pr_raw <= 0 or not np.isfinite(pr_raw):
                continue
            strength = pr_raw if pr_raw >= 1.0 else 1.0 / pr_raw
            if strength < min_P_ratio:
                continue
            nearest = int(np.argmin(np.abs(zc - sh['radius'])))
            # P_post = local peak pressure in a 5-zone window centred on the
            # shock zone. Using the nearest-zone value alone picks up the
            # low-density gas zone when the shock has just crossed an
            # interface, masking the actual post-shock pressure of the
            # transmitted/reflected wave on the compressed side.
            abs_zone = search_inner_zone + nearest
            lo = max(0, abs_zone - 2)
            hi = min(pressure.shape[1], abs_zone + 3)
            P_post_Gbar = float(np.max(pressure[i, lo:hi])) * 1e-8
            detections[i].append({
                'radius':      float(sh['radius']),
                'P_ratio':     float(strength),
                'P_post_Gbar': float(P_post_Gbar),
            })

    # ---- 2. Greedy linking into trajectories ----
    trajectories: List[Dict] = []

    def _seed(i: int, det: Dict) -> Dict:
        return {
            'indices':     [i],
            'radius':      [det['radius']],
            'time_ns':     [time[i]],
            'P_ratio':     [det['P_ratio']],
            'P_post_Gbar': [det['P_post_Gbar']],
            'started':     i,
            'ended':       i,
            'reason_ended': 'active',
            '_active':     True,
            '_misses':     0,
        }

    for i in range(n_times):
        dets = detections[i]
        # Iterate active trajectories in deterministic order (innermost first)
        active = sorted(
            [tr for tr in trajectories if tr['_active']],
            key=lambda tr: tr['radius'][-1],
        )
        unmatched = list(range(len(dets)))
        for tr in active:
            r_last = tr['radius'][-1]
            i_last = tr['indices'][-1]
            dt = time[i] - time[i_last]
            if dt <= 0.0:
                continue
            best_j = -1
            best_dr = np.inf
            for j in unmatched:
                r_j = dets[j]['radius']
                if r_j > r_last:
                    continue
                v = (r_last - r_j) / dt
                if v > max_shock_velocity:
                    continue
                dr = r_last - r_j
                if dr < best_dr:
                    best_dr = dr
                    best_j = j
            if best_j >= 0:
                d = dets[best_j]
                tr['indices'].append(i)
                tr['radius'].append(d['radius'])
                tr['time_ns'].append(time[i])
                tr['P_ratio'].append(d['P_ratio'])
                tr['P_post_Gbar'].append(d['P_post_Gbar'])
                tr['ended'] = i
                tr['_misses'] = 0
                unmatched.remove(best_j)
                # Terminate trajectory at the first gas/ice breakout so
                # subsequent shocks have to seed new trajectories rather
                # than chaining onto the post-breakout Lagrangian
                # interface continuation. Matches LILAC R-T convention.
                if d['radius'] <= interface_radius[i] + breakout_tolerance:
                    recent_pr = tr['P_ratio'][-(lookback_steps + 1):]
                    if recent_pr and max(recent_pr) >= min_P_ratio:
                        tr['_breakout_record'] = {
                            't_idx':       i,
                            'time_ns':     float(time[i]),
                            'radius':      float(d['radius']),
                            'P_post_Gbar': float(d['P_post_Gbar']),
                            'P_ratio':     float(d['P_ratio']),
                        }
                        tr['_active'] = False
                        tr['reason_ended'] = 'breakout'
            else:
                tr['_misses'] += 1
                if tr['_misses'] > max_gap_steps:
                    tr['_active'] = False
                    tr['reason_ended'] = 'no_match'

        for j in unmatched:
            trajectories.append(_seed(i, dets[j]))

    for tr in trajectories:
        if tr['_active']:
            tr['reason_ended'] = 'end_of_sim'
        del tr['_active']
        del tr['_misses']
        tr['indices']     = np.asarray(tr['indices'], dtype=int)
        tr['radius']      = np.asarray(tr['radius'])
        tr['time_ns']     = np.asarray(tr['time_ns'])
        tr['P_ratio']     = np.asarray(tr['P_ratio'])
        tr['P_post_Gbar'] = np.asarray(tr['P_post_Gbar'])

    # ---- 3. Coalescence: post-hoc pairwise convergence scan ----
    coalescence_events: List[Dict] = []
    n_tr = len(trajectories)
    for a in range(n_tr):
        ta = trajectories[a]
        ia_map = {int(idx): k for k, idx in enumerate(ta['indices'])}
        for b in range(a + 1, n_tr):
            tb = trajectories[b]
            common = sorted(set(ia_map.keys()) & {int(idx) for idx in tb['indices']})
            for i_t in common:
                ka = ia_map[i_t]
                kb = int(np.where(tb['indices'] == i_t)[0][0])
                ra = float(ta['radius'][ka])
                rb = float(tb['radius'][kb])
                if abs(ra - rb) < min_separation:
                    coalescence_events.append({
                        't_idx':            int(i_t),
                        'time_ns':          float(time[i_t]),
                        'radius':           0.5 * (ra + rb),
                        'trajectory_pair':  (a, b),
                    })
                    older   = a if ta['started'] <= tb['started'] else b
                    younger = b if older == a else a
                    if not trajectories[younger]['reason_ended'].startswith('coalesced_into'):
                        trajectories[younger]['reason_ended'] = f'coalesced_into:{older}'
                    break  # one event per pair (earliest convergence)

    # ---- 4. Breakout records: collected during linking ----
    # Trajectories self-record their breakout at the first interface
    # crossing in step 2 (then terminate). Harvest those records here.
    breakouts: List[Dict] = []
    for tr_id, tr in enumerate(trajectories):
        rec = tr.pop('_breakout_record', None)
        if rec is not None:
            rec = dict(rec)
            rec['trajectory_id'] = tr_id
            breakouts.append(rec)
    breakouts.sort(key=lambda b: b['time_ns'])

    return {
        'trajectories':       trajectories,
        'coalescence_events': coalescence_events,
        'breakouts':          breakouts,
        'time':               time,
        'interface_radius':   interface_radius,
    }
