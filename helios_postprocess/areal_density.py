"""
Areal Density (ρR) Calculations for ICF Simulations

Areal density is the integral of density along a radial path: ρR = ∫ρ dr

This is one of the most critical ICF metrics because it determines:
- Neutron confinement time (τ ∝ ρR)
- Alpha particle energy deposition
- Hot spot tamping
- Burn efficiency

Units:
- Density: g/cm³
- Radius: cm
- Areal density: g/cm²

Functions:
- calculate_areal_density: Basic ∫ρ dr from center outward
- calculate_areal_density_region: ρR in specific radial region
- calculate_neutron_averaged_areal_density: Weighted by DT reaction rate
- find_peak_areal_density: Find max ρR and when it occurs
- areal_density_evolution: Track ρR(t) at specific radius or region
- shell_areal_density: ρR in shell/pusher region
- hot_spot_areal_density: ρR in hot spot core
"""

import numpy as np
from typing import Optional, Tuple, Dict
from scipy.integrate import cumulative_trapezoid


def calculate_areal_density(
    density: np.ndarray,
    zone_boundaries: np.ndarray,
    max_radius: Optional[float] = None
) -> np.ndarray:
    """
    Calculate cumulative areal density from center outward: ρR = ∫ρ dr
    
    Integrates from r=0 to each zone boundary.
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
        Assumes boundaries[0] is inner edge, boundaries[-1] is outer edge
    max_radius : float, optional
        Maximum radius to integrate to. If None, uses full domain.
        
    Returns
    -------
    np.ndarray
        Cumulative areal density at each zone boundary [g/cm²]
        Shape (n_zones+1,)
        
    Examples
    --------
    >>> # Get areal density profile
    >>> rho_r = calculate_areal_density(run.density[idx], run.zone_boundaries[idx])
    >>> # Total areal density to outer radius
    >>> rho_r_total = rho_r[-1]
    >>> print(f"Total ρR: {rho_r_total:.3f} g/cm²")
    
    Notes
    -----
    - Returns cumulative integral, so rho_r[i] = ∫₀^(r[i]) ρ dr
    - For areal density in a shell, use rho_r[outer] - rho_r[inner]
    - High ρR (>0.3 g/cm²) indicates good confinement for ICF
    """
    n_zones = len(density)
    
    if max_radius is not None:
        # Find zones within max_radius
        mask = zone_boundaries[:-1] <= max_radius
        density = density[mask]
        zone_boundaries = zone_boundaries[:len(density)+1]
    
    # Zone centers (arithmetic mean of boundaries)
    zone_centers = 0.5 * (zone_boundaries[:-1] + zone_boundaries[1:])
    
    # Cumulative integral using trapezoidal rule
    # We integrate ρ(r) dr from center outward
    rho_r_zones = cumulative_trapezoid(
        density,
        zone_centers,
        initial=0.0
    )
    
    # Interpolate to zone boundaries
    # rho_r_zones is at zone centers, we want at boundaries
    rho_r_boundaries = np.zeros(len(zone_boundaries))
    rho_r_boundaries[0] = 0.0  # Zero at center
    
    # Linear interpolation to boundaries
    for i in range(1, len(zone_boundaries)-1):
        # Boundary i is between zones i-1 and i
        r_boundary = zone_boundaries[i]
        r_left = zone_centers[i-1]
        r_right = zone_centers[i]
        
        # Linear interpolation weight
        w = (r_boundary - r_left) / (r_right - r_left)
        rho_r_boundaries[i] = (1-w) * rho_r_zones[i-1] + w * rho_r_zones[i]
    
    # Last boundary gets last zone value
    rho_r_boundaries[-1] = rho_r_zones[-1] + \
        0.5 * density[-1] * (zone_boundaries[-1] - zone_centers[-1])
    
    return rho_r_boundaries


def calculate_areal_density_region(
    density: np.ndarray,
    zone_boundaries: np.ndarray,
    r_inner: float,
    r_outer: float
) -> float:
    """
    Calculate areal density in a specific radial region
    
    Computes ρR = ∫(r_inner to r_outer) ρ dr
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    r_inner : float
        Inner radius of region [cm]
    r_outer : float
        Outer radius of region [cm]
        
    Returns
    -------
    float
        Areal density in region [g/cm²]
        
    Examples
    --------
    >>> # Shell areal density
    >>> rho_r_shell = calculate_areal_density_region(
    ...     run.density[idx], 
    ...     run.zone_boundaries[idx],
    ...     r_inner=0.005,  # 50 μm
    ...     r_outer=0.010   # 100 μm
    ... )
    """
    # Get cumulative areal density
    rho_r = calculate_areal_density(density, zone_boundaries)
    
    # Interpolate to desired radii
    rho_r_inner = np.interp(r_inner, zone_boundaries, rho_r)
    rho_r_outer = np.interp(r_outer, zone_boundaries, rho_r)
    
    return rho_r_outer - rho_r_inner


def calculate_neutron_averaged_areal_density(
    density: np.ndarray,
    temperature: np.ndarray,
    zone_boundaries: np.ndarray,
    fuel_fraction: Optional[np.ndarray] = None
) -> float:
    """
    Calculate neutron-averaged areal density (ρR)ₙ
    
    Weighted by DT reaction rate: (ρR)ₙ = ∫ρ²<σv> dr / ∫ρ<σv> dr
    
    This gives the effective areal density seen by neutrons born in the fuel.
    More relevant for understanding neutron confinement than simple ρR.
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    temperature : np.ndarray
        Ion temperature [eV], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    fuel_fraction : np.ndarray, optional
        Fraction of DT fuel in each zone (0-1). If None, assumes pure fuel.
        
    Returns
    -------
    float
        Neutron-averaged areal density [g/cm²]
        
    Examples
    --------
    >>> rho_r_neutron = calculate_neutron_averaged_areal_density(
    ...     run.density[idx],
    ...     run.temperature[idx],
    ...     run.zone_boundaries[idx]
    ... )
    >>> print(f"Neutron-averaged ρR: {rho_r_neutron:.3f} g/cm²")
    
    Notes
    -----
    Uses Bosch-Hale DT reactivity fit (1992)
    Valid for 0.2 keV < T < 100 keV
    """
    # Zone centers
    zone_centers = 0.5 * (zone_boundaries[:-1] + zone_boundaries[1:])
    dr = zone_boundaries[1:] - zone_boundaries[:-1]
    
    # DT reactivity <σv> as function of temperature
    # Bosch-Hale 1992 fit for DT
    reactivity = _dt_reactivity(temperature)
    
    if fuel_fraction is not None:
        # Only count fuel regions
        reactivity *= fuel_fraction
    
    # Numerator: ∫ρ²<σv> dr (neutron-weighted areal density)
    numerator = np.sum(density**2 * reactivity * dr)
    
    # Denominator: ∫ρ<σv> dr (total reaction rate)
    denominator = np.sum(density * reactivity * dr)
    
    if denominator < 1e-30:
        return 0.0
    
    return numerator / denominator


def _dt_reactivity(T_eV: np.ndarray) -> np.ndarray:
    """
    DT fusion reactivity <σv> using Bosch-Hale 1992 fit
    
    Parameters
    ----------
    T_eV : np.ndarray
        Ion temperature [eV]
        
    Returns
    -------
    np.ndarray
        Reactivity [cm³/s]
        
    Notes
    -----
    Valid for 0.2 keV < T < 100 keV
    For T < 0.2 keV, returns very small value
    For T > 100 keV, extrapolates (less accurate)
    """
    # Convert to keV
    T_keV = T_eV / 1000.0
    
    # Bosch-Hale parameters for DT
    C1 = 1.17302e-9
    C2 = 1.51361e-2
    C3 = 7.51886e-2
    C4 = 4.60643e-3
    C5 = 1.35000e-2
    C6 = -1.06750e-4
    C7 = 1.36600e-5
    
    # Protect against very low temperatures
    T_keV = np.maximum(T_keV, 0.1)
    
    theta = T_keV / (1 - (C2*T_keV + C4*T_keV**2 + C6*T_keV**3) / 
                          (1 + C3*T_keV + C5*T_keV**2 + C7*T_keV**3))
    
    xi = (C1 / T_keV**(1/3))**2
    
    reactivity = C1 * theta * np.sqrt(xi / (np.exp(xi) * T_keV**3))
    
    return reactivity


def find_peak_areal_density(
    density: np.ndarray,
    zone_boundaries: np.ndarray,
    time: np.ndarray,
    region: str = 'total'
) -> Tuple[float, float, int]:
    """
    Find peak areal density and when it occurs
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_times, n_zones)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_times, n_zones+1)
    time : np.ndarray
        Time array [ns], shape (n_times,)
    region : str, default='total'
        Region to analyze:
        - 'total': Entire domain
        - 'core': Inner 50% by mass
        - 'shell': Outer 50% by mass
        
    Returns
    -------
    peak_rho_r : float
        Peak areal density [g/cm²]
    peak_time : float
        Time of peak [ns]
    peak_idx : int
        Time index of peak
        
    Examples
    --------
    >>> rho_r_peak, t_peak, idx = find_peak_areal_density(
    ...     run.density,
    ...     run.zone_boundaries,
    ...     run.time,
    ...     region='total'
    ... )
    >>> print(f"Peak ρR: {rho_r_peak:.3f} g/cm² at t={t_peak:.2f} ns")
    """
    n_times = len(time)
    rho_r_max = np.zeros(n_times)
    
    for i in range(n_times):
        rho_r = calculate_areal_density(density[i], zone_boundaries[i])
        
        if region == 'total':
            rho_r_max[i] = rho_r[-1]
        elif region == 'core':
            # Inner 50% by radius
            r_mid = zone_boundaries[i, -1] * 0.5
            rho_r_max[i] = np.interp(r_mid, zone_boundaries[i], rho_r)
        elif region == 'shell':
            # Outer 50% by radius
            r_mid = zone_boundaries[i, -1] * 0.5
            rho_r_inner = np.interp(r_mid, zone_boundaries[i], rho_r)
            rho_r_max[i] = rho_r[-1] - rho_r_inner
        else:
            raise ValueError(f"Unknown region: {region}")
    
    peak_idx = np.argmax(rho_r_max)
    peak_rho_r = rho_r_max[peak_idx]
    peak_time = time[peak_idx]
    
    return peak_rho_r, peak_time, peak_idx


def areal_density_evolution(
    density: np.ndarray,
    zone_boundaries: np.ndarray,
    time: np.ndarray,
    radius: Optional[float] = None
) -> np.ndarray:
    """
    Track areal density evolution ρR(t)
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_times, n_zones)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_times, n_zones+1)
    time : np.ndarray
        Time array [ns], shape (n_times,)
    radius : float, optional
        Radius to evaluate ρR at [cm]. If None, uses outer boundary.
        
    Returns
    -------
    np.ndarray
        Areal density vs time [g/cm²], shape (n_times,)
        
    Examples
    --------
    >>> # Track total areal density
    >>> rho_r_t = areal_density_evolution(
    ...     run.density, run.zone_boundaries, run.time
    ... )
    >>> 
    >>> # Track at specific radius (e.g., initial outer radius)
    >>> r_initial = run.zone_boundaries[0, -1]
    >>> rho_r_shell_t = areal_density_evolution(
    ...     run.density, run.zone_boundaries, run.time,
    ...     radius=r_initial
    ... )
    """
    n_times = len(time)
    rho_r_t = np.zeros(n_times)
    
    for i in range(n_times):
        rho_r = calculate_areal_density(density[i], zone_boundaries[i])
        
        if radius is None:
            # Use outer boundary
            rho_r_t[i] = rho_r[-1]
        else:
            # Interpolate to desired radius
            rho_r_t[i] = np.interp(radius, zone_boundaries[i], rho_r)
    
    return rho_r_t


def hot_spot_areal_density(
    density: np.ndarray,
    temperature: np.ndarray,
    zone_boundaries: np.ndarray,
    T_threshold: float = 1000.0
) -> float:
    """
    Calculate areal density in hot spot (T > threshold)
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    T_threshold : float, default=1000.0
        Temperature threshold for hot spot [eV]
        
    Returns
    -------
    float
        Hot spot areal density [g/cm²]
        
    Examples
    --------
    >>> rho_r_hs = hot_spot_areal_density(
    ...     run.density[idx],
    ...     run.temperature[idx],
    ...     run.zone_boundaries[idx],
    ...     T_threshold=1000.0  # 1 keV
    ... )
    """
    # Find hot spot region
    hot_mask = temperature > T_threshold
    
    if not np.any(hot_mask):
        return 0.0
    
    # Find outermost hot zone
    hot_indices = np.where(hot_mask)[0]
    r_outer = zone_boundaries[hot_indices[-1] + 1]
    
    # Calculate areal density from center to hot spot edge
    rho_r = calculate_areal_density(density, zone_boundaries, max_radius=r_outer)
    
    return rho_r[-1]


def shell_areal_density(
    density: np.ndarray,
    zone_boundaries: np.ndarray,
    zone_mass: np.ndarray,
    mass_fraction: float = 0.5
) -> float:
    """
    Calculate areal density in shell (outer fraction by mass)
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_zones,)
    mass_fraction : float, default=0.5
        Fraction of mass to include in shell (0-1)
        0.5 = outer half by mass
        
    Returns
    -------
    float
        Shell areal density [g/cm²]
        
    Examples
    --------
    >>> # Outer half by mass
    >>> rho_r_shell = shell_areal_density(
    ...     run.density[idx],
    ...     run.zone_boundaries[idx],
    ...     run.zone_mass[idx],
    ...     mass_fraction=0.5
    ... )
    """
    # Find boundary between inner and outer regions
    total_mass = np.sum(zone_mass)
    cumulative_mass = np.cumsum(zone_mass)
    
    # Find zone where cumulative mass = (1 - mass_fraction) * total_mass
    target_mass = (1 - mass_fraction) * total_mass
    idx_boundary = np.searchsorted(cumulative_mass, target_mass)
    
    if idx_boundary >= len(zone_mass):
        idx_boundary = len(zone_mass) - 1
    
    r_inner = zone_boundaries[idx_boundary]
    r_outer = zone_boundaries[-1]
    
    return calculate_areal_density_region(
        density, zone_boundaries, r_inner, r_outer
    )


def comprehensive_areal_density_analysis(
    density: np.ndarray,
    temperature: np.ndarray,
    zone_boundaries: np.ndarray,
    zone_mass: np.ndarray,
    time_idx: Optional[int] = None
) -> Dict[str, float]:
    """
    Comprehensive areal density analysis at a single time
    
    Calculates multiple ρR metrics:
    - Total areal density
    - Hot spot areal density  
    - Shell areal density
    - Core areal density
    - Neutron-averaged areal density
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_zones,)
    time_idx : int, optional
        Time index for labeling (for printing)
        
    Returns
    -------
    dict
        Dictionary with keys:
        - 'rho_r_total': Total areal density [g/cm²]
        - 'rho_r_hot_spot': Hot spot areal density [g/cm²]
        - 'rho_r_shell': Shell areal density [g/cm²]
        - 'rho_r_core': Core areal density [g/cm²]
        - 'rho_r_neutron': Neutron-averaged areal density [g/cm²]
        - 'r_outer': Outer radius [cm]
        
    Examples
    --------
    >>> results = comprehensive_areal_density_analysis(
    ...     run.density[idx],
    ...     run.temperature[idx],
    ...     run.zone_boundaries[idx],
    ...     run.zone_mass[idx]
    ... )
    >>> print(f"Total ρR: {results['rho_r_total']:.3f} g/cm²")
    >>> print(f"Hot spot ρR: {results['rho_r_hot_spot']:.3f} g/cm²")
    """
    results = {}
    
    # Total areal density
    rho_r = calculate_areal_density(density, zone_boundaries)
    results['rho_r_total'] = rho_r[-1]
    results['r_outer'] = zone_boundaries[-1]
    
    # Hot spot areal density (T > 1 keV)
    results['rho_r_hot_spot'] = hot_spot_areal_density(
        density, temperature, zone_boundaries, T_threshold=1000.0
    )
    
    # Shell areal density (outer 50% by mass)
    results['rho_r_shell'] = shell_areal_density(
        density, zone_boundaries, zone_mass, mass_fraction=0.5
    )
    
    # Core areal density (inner 50% by mass)
    results['rho_r_core'] = shell_areal_density(
        density, zone_boundaries, zone_mass, mass_fraction=0.5
    )
    
    # Neutron-averaged areal density
    results['rho_r_neutron'] = calculate_neutron_averaged_areal_density(
        density, temperature, zone_boundaries
    )
    
    return results
