"""
Hot Spot Metrics for ICF Simulations

The hot spot is the central high-temperature region where most fusion occurs.
Understanding hot spot properties is critical for predicting ignition.

Key metrics:
- Mass-averaged temperature: <T> = Σ(m_i * T_i) / Σm_i
- Mass-averaged pressure: <P> = Σ(m_i * P_i) / Σm_i
- Burn-averaged temperature: Weighted by reaction rate
- Hot spot mass, volume, radius
- Hot spot ρR (see areal_density module)

Ignition criteria (Lindl):
- <T> > 5 keV (mass-averaged)
- <P> > 100 Gbar
- Hot spot mass ~ 0.1-1 mg
- ρR > 0.3 g/cm²

Units:
- Temperature: eV or keV
- Pressure: J/cm³ (divide by 1e5 for Mbar, 1e6 for Gbar)
- Mass: g
- Volume: cm³
- Radius: cm

Author: Prof T
"""

import numpy as np
from typing import Optional, Tuple, Dict
from scipy.integrate import trapezoid


def identify_hot_spot(
    temperature: np.ndarray,
    zone_boundaries: np.ndarray,
    T_threshold: float = 1000.0
) -> Tuple[np.ndarray, float]:
    """
    Identify hot spot region (T > threshold)
    
    Parameters
    ----------
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    T_threshold : float, default=1000.0
        Temperature threshold [eV]
        
    Returns
    -------
    hot_mask : np.ndarray
        Boolean mask for hot zones
    r_hot_spot : float
        Hot spot outer radius [cm]
        
    Examples
    --------
    >>> hot_mask, r_hs = identify_hot_spot(
    ...     run.temperature[idx],
    ...     run.zone_boundaries[idx],
    ...     T_threshold=1000.0  # 1 keV
    ... )
    >>> print(f"Hot spot radius: {r_hs*1e4:.1f} μm")
    """
    hot_mask = temperature > T_threshold
    
    if not np.any(hot_mask):
        return hot_mask, 0.0
    
    # Find outermost hot zone
    hot_indices = np.where(hot_mask)[0]
    r_hot_spot = zone_boundaries[hot_indices[-1] + 1]
    
    return hot_mask, r_hot_spot


def hot_spot_mass(
    density: np.ndarray,
    zone_mass: np.ndarray,
    temperature: np.ndarray,
    T_threshold: float = 1000.0
) -> float:
    """
    Calculate total mass in hot spot
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_zones,)
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    T_threshold : float, default=1000.0
        Temperature threshold [eV]
        
    Returns
    -------
    float
        Hot spot mass [g]
        
    Examples
    --------
    >>> m_hs = hot_spot_mass(
    ...     run.density[idx],
    ...     run.zone_mass[idx],
    ...     run.temperature[idx],
    ...     T_threshold=1000.0
    ... )
    >>> print(f"Hot spot mass: {m_hs*1e3:.2f} mg")
    
    Notes
    -----
    Typical hot spot masses for ignition: 0.1-1 mg
    """
    hot_mask = temperature > T_threshold
    
    if not np.any(hot_mask):
        return 0.0
    
    return np.sum(zone_mass[hot_mask])


def hot_spot_volume(
    zone_boundaries: np.ndarray,
    temperature: np.ndarray,
    T_threshold: float = 1000.0
) -> float:
    """
    Calculate hot spot volume
    
    Assumes spherical geometry: V = (4/3)πr³
    
    Parameters
    ----------
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    T_threshold : float, default=1000.0
        Temperature threshold [eV]
        
    Returns
    -------
    float
        Hot spot volume [cm³]
        
    Examples
    --------
    >>> V_hs = hot_spot_volume(
    ...     run.zone_boundaries[idx],
    ...     run.temperature[idx]
    ... )
    >>> print(f"Hot spot volume: {V_hs:.2e} cm³")
    """
    hot_mask = temperature > T_threshold
    
    if not np.any(hot_mask):
        return 0.0
    
    # Find outer radius
    hot_indices = np.where(hot_mask)[0]
    r_hot_spot = zone_boundaries[hot_indices[-1] + 1]
    
    # Spherical volume
    volume = (4.0/3.0) * np.pi * r_hot_spot**3
    
    return volume


def mass_averaged_temperature(
    temperature: np.ndarray,
    zone_mass: np.ndarray,
    hot_spot_only: bool = False,
    T_threshold: float = 1000.0
) -> float:
    """
    Calculate mass-averaged temperature
    
    <T> = Σ(m_i * T_i) / Σm_i
    
    This is THE key temperature metric for ignition assessment.
    
    Parameters
    ----------
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_zones,)
    hot_spot_only : bool, default=False
        If True, only average over hot spot (T > threshold)
    T_threshold : float, default=1000.0
        Temperature threshold for hot spot [eV]
        
    Returns
    -------
    float
        Mass-averaged temperature [eV]
        
    Examples
    --------
    >>> # Whole domain
    >>> T_avg = mass_averaged_temperature(
    ...     run.temperature[idx],
    ...     run.zone_mass[idx]
    ... )
    >>> print(f"<T>: {T_avg/1000:.2f} keV")
    >>> 
    >>> # Hot spot only
    >>> T_hs_avg = mass_averaged_temperature(
    ...     run.temperature[idx],
    ...     run.zone_mass[idx],
    ...     hot_spot_only=True,
    ...     T_threshold=1000.0
    ... )
    >>> print(f"<T>_hs: {T_hs_avg/1000:.2f} keV")
    
    Notes
    -----
    Ignition criterion: <T> > 5 keV in hot spot
    """
    if hot_spot_only:
        mask = temperature > T_threshold
        if not np.any(mask):
            return 0.0
        temperature = temperature[mask]
        zone_mass = zone_mass[mask]
    
    total_mass = np.sum(zone_mass)
    
    if total_mass < 1e-30:
        return 0.0
    
    return np.sum(temperature * zone_mass) / total_mass


def mass_averaged_pressure(
    pressure: np.ndarray,
    zone_mass: np.ndarray,
    hot_spot_only: bool = False,
    T_threshold: float = 1000.0,
    temperature: Optional[np.ndarray] = None
) -> float:
    """
    Calculate mass-averaged pressure
    
    <P> = Σ(m_i * P_i) / Σm_i
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_zones,)
    hot_spot_only : bool, default=False
        If True, only average over hot spot
    T_threshold : float, default=1000.0
        Temperature threshold for hot spot [eV]
    temperature : np.ndarray, optional
        Temperature field [eV], required if hot_spot_only=True
        
    Returns
    -------
    float
        Mass-averaged pressure [J/cm³]
        Divide by 1e5 for Mbar, 1e6 for Gbar
        
    Examples
    --------
    >>> P_avg = mass_averaged_pressure(
    ...     run.pressure[idx],
    ...     run.zone_mass[idx],
    ...     hot_spot_only=True,
    ...     T_threshold=1000.0,
    ...     temperature=run.temperature[idx]
    ... )
    >>> print(f"<P>_hs: {P_avg/1e6:.1f} Gbar")
    
    Notes
    -----
    Ignition criterion: <P> > 100 Gbar in hot spot
    """
    if hot_spot_only:
        if temperature is None:
            raise ValueError("temperature required when hot_spot_only=True")
        mask = temperature > T_threshold
        if not np.any(mask):
            return 0.0
        pressure = pressure[mask]
        zone_mass = zone_mass[mask]
    
    total_mass = np.sum(zone_mass)
    
    if total_mass < 1e-30:
        return 0.0
    
    return np.sum(pressure * zone_mass) / total_mass


def mass_averaged_density(
    density: np.ndarray,
    zone_mass: np.ndarray,
    hot_spot_only: bool = False,
    T_threshold: float = 1000.0,
    temperature: Optional[np.ndarray] = None
) -> float:
    """
    Calculate mass-averaged density
    
    <ρ> = Σ(m_i * ρ_i) / Σm_i
    
    Parameters
    ----------
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_zones,)
    hot_spot_only : bool, default=False
        If True, only average over hot spot
    T_threshold : float, default=1000.0
        Temperature threshold for hot spot [eV]
    temperature : np.ndarray, optional
        Temperature field [eV], required if hot_spot_only=True
        
    Returns
    -------
    float
        Mass-averaged density [g/cm³]
    """
    if hot_spot_only:
        if temperature is None:
            raise ValueError("temperature required when hot_spot_only=True")
        mask = temperature > T_threshold
        if not np.any(mask):
            return 0.0
        density = density[mask]
        zone_mass = zone_mass[mask]
    
    total_mass = np.sum(zone_mass)
    
    if total_mass < 1e-30:
        return 0.0
    
    return np.sum(density * zone_mass) / total_mass


def burn_averaged_temperature(
    temperature: np.ndarray,
    density: np.ndarray,
    zone_mass: np.ndarray,
    zone_boundaries: np.ndarray
) -> float:
    """
    Calculate burn-averaged (reactivity-weighted) temperature
    
    Weighted by DT reaction rate: <T>_burn = Σ(ρ²<σv> V_i T_i) / Σ(ρ²<σv> V_i)
    
    More representative than mass-averaged for understanding burn conditions.
    
    Parameters
    ----------
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
        
    Returns
    -------
    float
        Burn-averaged temperature [eV]
        
    Examples
    --------
    >>> T_burn = burn_averaged_temperature(
    ...     run.temperature[idx],
    ...     run.density[idx],
    ...     run.zone_mass[idx],
    ...     run.zone_boundaries[idx]
    ... )
    >>> print(f"<T>_burn: {T_burn/1000:.2f} keV")
    
    Notes
    -----
    Uses Bosch-Hale DT reactivity fit
    """
    # Get zone volumes
    r_inner = zone_boundaries[:-1]
    r_outer = zone_boundaries[1:]
    volumes = (4.0/3.0) * np.pi * (r_outer**3 - r_inner**3)
    
    # DT reactivity
    reactivity = _dt_reactivity(temperature)
    
    # Reaction rate weight: ρ² <σv> V
    weights = density**2 * reactivity * volumes
    
    total_weight = np.sum(weights)
    
    if total_weight < 1e-30:
        return 0.0
    
    return np.sum(weights * temperature) / total_weight


def burn_averaged_pressure(
    pressure: np.ndarray,
    temperature: np.ndarray,
    density: np.ndarray,
    zone_boundaries: np.ndarray
) -> float:
    """
    Calculate burn-averaged pressure
    
    Weighted by DT reaction rate
    
    Parameters
    ----------
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundary radii [cm], shape (n_zones+1,)
        
    Returns
    -------
    float
        Burn-averaged pressure [J/cm³]
    """
    # Zone volumes
    r_inner = zone_boundaries[:-1]
    r_outer = zone_boundaries[1:]
    volumes = (4.0/3.0) * np.pi * (r_outer**3 - r_inner**3)
    
    # DT reactivity
    reactivity = _dt_reactivity(temperature)
    
    # Weights
    weights = density**2 * reactivity * volumes
    
    total_weight = np.sum(weights)
    
    if total_weight < 1e-30:
        return 0.0
    
    return np.sum(weights * pressure) / total_weight


def _dt_reactivity(T_eV: np.ndarray) -> np.ndarray:
    """
    DT fusion reactivity using Bosch-Hale 1992 fit
    
    Parameters
    ----------
    T_eV : np.ndarray
        Temperature [eV]
        
    Returns
    -------
    np.ndarray
        Reactivity [cm³/s]
    """
    # Convert to keV
    T_keV = T_eV / 1000.0
    
    # Bosch-Hale parameters
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


def hot_spot_evolution(
    temperature: np.ndarray,
    density: np.ndarray,
    pressure: np.ndarray,
    zone_mass: np.ndarray,
    zone_boundaries: np.ndarray,
    time: np.ndarray,
    T_threshold: float = 1000.0
) -> Dict[str, np.ndarray]:
    """
    Track hot spot metrics vs time
    
    Parameters
    ----------
    temperature : np.ndarray
        Temperature field [eV], shape (n_times, n_zones)
    density : np.ndarray
        Density field [g/cm³], shape (n_times, n_zones)
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_times, n_zones)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_times, n_zones)
    zone_boundaries : np.ndarray
        Zone boundaries [cm], shape (n_times, n_zones+1)
    time : np.ndarray
        Time array [ns], shape (n_times,)
    T_threshold : float, default=1000.0
        Temperature threshold [eV]
        
    Returns
    -------
    dict
        Dictionary with time series:
        - 'time': Time array [ns]
        - 'T_avg': Mass-averaged T [eV]
        - 'P_avg': Mass-averaged P [J/cm³]
        - 'rho_avg': Mass-averaged ρ [g/cm³]
        - 'mass': Hot spot mass [g]
        - 'radius': Hot spot radius [cm]
        - 'volume': Hot spot volume [cm³]
        
    Examples
    --------
    >>> evolution = hot_spot_evolution(
    ...     run.temperature, run.density, run.pressure,
    ...     run.zone_mass, run.zone_boundaries, run.time,
    ...     T_threshold=1000.0
    ... )
    >>> 
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(evolution['time'], evolution['T_avg']/1000)
    >>> plt.xlabel('Time [ns]')
    >>> plt.ylabel('<T> [keV]')
    >>> plt.show()
    """
    n_times = len(time)
    
    results = {
        'time': time,
        'T_avg': np.zeros(n_times),
        'P_avg': np.zeros(n_times),
        'rho_avg': np.zeros(n_times),
        'mass': np.zeros(n_times),
        'radius': np.zeros(n_times),
        'volume': np.zeros(n_times)
    }
    
    for i in range(n_times):
        # Temperature
        results['T_avg'][i] = mass_averaged_temperature(
            temperature[i], zone_mass[i], 
            hot_spot_only=True, T_threshold=T_threshold
        )
        
        # Pressure
        results['P_avg'][i] = mass_averaged_pressure(
            pressure[i], zone_mass[i],
            hot_spot_only=True, T_threshold=T_threshold,
            temperature=temperature[i]
        )
        
        # Density
        results['rho_avg'][i] = mass_averaged_density(
            density[i], zone_mass[i],
            hot_spot_only=True, T_threshold=T_threshold,
            temperature=temperature[i]
        )
        
        # Mass
        results['mass'][i] = hot_spot_mass(
            density[i], zone_mass[i], temperature[i], T_threshold
        )
        
        # Volume and radius
        results['volume'][i] = hot_spot_volume(
            zone_boundaries[i], temperature[i], T_threshold
        )
        
        if results['volume'][i] > 0:
            results['radius'][i] = (3*results['volume'][i]/(4*np.pi))**(1/3)
        else:
            results['radius'][i] = 0.0
    
    return results


def comprehensive_hot_spot_analysis(
    temperature: np.ndarray,
    density: np.ndarray,
    pressure: np.ndarray,
    zone_mass: np.ndarray,
    zone_boundaries: np.ndarray,
    T_threshold: float = 1000.0
) -> Dict[str, float]:
    """
    Comprehensive hot spot analysis at single time
    
    Calculates all key hot spot metrics for ignition assessment.
    
    Parameters
    ----------
    temperature : np.ndarray
        Temperature field [eV], shape (n_zones,)
    density : np.ndarray
        Density field [g/cm³], shape (n_zones,)
    pressure : np.ndarray
        Pressure field [J/cm³], shape (n_zones,)
    zone_mass : np.ndarray
        Mass in each zone [g], shape (n_zones,)
    zone_boundaries : np.ndarray
        Zone boundaries [cm], shape (n_zones+1,)
    T_threshold : float, default=1000.0
        Temperature threshold [eV]
        
    Returns
    -------
    dict
        Dictionary with:
        - 'T_mass_avg': Mass-averaged T [eV]
        - 'T_burn_avg': Burn-averaged T [eV]
        - 'T_peak': Peak T [eV]
        - 'P_mass_avg': Mass-averaged P [J/cm³]
        - 'P_burn_avg': Burn-averaged P [J/cm³]
        - 'P_peak': Peak P [J/cm³]
        - 'rho_avg': Mass-averaged ρ [g/cm³]
        - 'mass': Hot spot mass [g]
        - 'radius': Hot spot radius [cm]
        - 'volume': Hot spot volume [cm³]
        
    Examples
    --------
    >>> idx_stag = np.argmax(np.max(run.density, axis=1))
    >>> results = comprehensive_hot_spot_analysis(
    ...     run.temperature[idx_stag],
    ...     run.density[idx_stag],
    ...     run.pressure[idx_stag],
    ...     run.zone_mass[idx_stag],
    ...     run.zone_boundaries[idx_stag]
    ... )
    >>> 
    >>> print(f"<T>: {results['T_mass_avg']/1000:.2f} keV")
    >>> print(f"<P>: {results['P_mass_avg']/1e6:.1f} Gbar")
    >>> print(f"Mass: {results['mass']*1e3:.2f} mg")
    >>> 
    >>> # Check ignition criteria
    >>> if results['T_mass_avg'] > 5000 and results['P_mass_avg'] > 1e8:
    ...     print("✓ Ignition conditions met!")
    """
    results = {}
    
    # Temperature metrics
    results['T_mass_avg'] = mass_averaged_temperature(
        temperature, zone_mass, hot_spot_only=True, T_threshold=T_threshold
    )
    results['T_burn_avg'] = burn_averaged_temperature(
        temperature, density, zone_mass, zone_boundaries
    )
    
    hot_mask = temperature > T_threshold
    if np.any(hot_mask):
        results['T_peak'] = temperature[hot_mask].max()
    else:
        results['T_peak'] = 0.0
    
    # Pressure metrics
    results['P_mass_avg'] = mass_averaged_pressure(
        pressure, zone_mass, hot_spot_only=True, 
        T_threshold=T_threshold, temperature=temperature
    )
    results['P_burn_avg'] = burn_averaged_pressure(
        pressure, temperature, density, zone_boundaries
    )
    
    if np.any(hot_mask):
        results['P_peak'] = pressure[hot_mask].max()
    else:
        results['P_peak'] = 0.0
    
    # Density
    results['rho_avg'] = mass_averaged_density(
        density, zone_mass, hot_spot_only=True,
        T_threshold=T_threshold, temperature=temperature
    )
    
    # Geometric properties
    results['mass'] = hot_spot_mass(
        density, zone_mass, temperature, T_threshold
    )
    results['volume'] = hot_spot_volume(
        zone_boundaries, temperature, T_threshold
    )
    
    if results['volume'] > 0:
        results['radius'] = (3*results['volume']/(4*np.pi))**(1/3)
    else:
        results['radius'] = 0.0
    
    return results
