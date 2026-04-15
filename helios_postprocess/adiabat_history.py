"""
Adiabat History for ICF Implosion Simulations

Computes the mass-averaged fuel adiabat alpha(t) = P / P_Fermi at every
timestep in the cold DT fuel region, using the same physics as
ICFAnalyzer._compute_adiabat() but returning the full time history rather
than a single scalar at peak velocity.

Useful for:
- Diagnosing when and how entropy is deposited in the fuel
- Correlating first-shock strength with final adiabat
- Comparing adiabat evolution across flux limiter / EOS scan runs

Author: Prof T
"""
import numpy as np
from typing import Optional, Dict


def compute_adiabat_history(
    ion_pressure: np.ndarray,
    rad_pressure: np.ndarray,
    mass_density: np.ndarray,
    zone_mass: np.ndarray,
    region_interfaces_indices: np.ndarray,
    time: np.ndarray,
    rho0: float = 0.205,
    fuel_region_col: int = 0,
    fuel_outer_col: int = 1,
) -> Dict[str, np.ndarray]:
    """
    Compute mass-averaged adiabat in cold DT fuel at every timestep.

    Uses the Lindl convention:
        alpha = P / P_Fermi
        P_Fermi = 2.17 * (rho / rho0)^(5/3)  [Mbar]
    for equimolar DT ice (rho0 = 0.205 g/cc).

    The cold fuel zone range is taken from region_interfaces_indices:
        z_start = ri[t, fuel_region_col]   (hot-spot outer boundary)
        z_end   = ri[t, fuel_outer_col]    (cold fuel outer boundary)

    Parameters
    ----------
    ion_pressure : np.ndarray
        Ion pressure [J/cm³], shape (n_times, n_zones)
    rad_pressure : np.ndarray
        Radiation+electron pressure [J/cm³], shape (n_times, n_zones).
        Added to ion_pressure for total pressure.
    mass_density : np.ndarray
        Mass density [g/cm³], shape (n_times, n_zones)
    zone_mass : np.ndarray
        Zone mass [g], shape (n_times, n_zones)
    region_interfaces_indices : np.ndarray
        Region interface NODE indices, shape (n_times, n_regions+1).
        ri[:, fuel_region_col] = hot-spot boundary node (inner cold fuel edge)
        ri[:, fuel_outer_col]  = cold fuel outer boundary node
    time : np.ndarray
        Time array [ns], shape (n_times,)
    rho0 : float, default=0.205
        Reference DT ice density [g/cm³] (Lindl convention).
    fuel_region_col : int, default=0
        Column index in ri for inner cold fuel boundary (hot-spot edge).
        PDD_9: 0 (node 151); VI_6: 0 (node 51)
    fuel_outer_col : int, default=1
        Column index in ri for outer cold fuel boundary.
        PDD_9: 1 (node 191); VI_6: 1 (node 101)

    Returns
    -------
    dict with keys:
        'time' : np.ndarray [ns]
            Time array (same as input).
        'adiabat' : np.ndarray
            Mass-averaged adiabat alpha(t) at each timestep.
            NaN where cold fuel region is empty or data invalid.
        'adiabat_min' : np.ndarray
            Minimum zone adiabat in cold fuel at each timestep.
        'adiabat_max' : np.ndarray
            Maximum zone adiabat in cold fuel at each timestep.
        'P_fuel_mean_Gbar' : np.ndarray
            Mass-averaged total pressure in cold fuel [Gbar].
        'rho_fuel_mean' : np.ndarray
            Mass-averaged density in cold fuel [g/cm³].
        'z_start' : np.ndarray (int)
            Inner cold fuel zone index at each timestep.
        'z_end' : np.ndarray (int)
            Outer cold fuel zone index at each timestep.

    Notes
    -----
    The adiabat is meaningful in the cold fuel before stagnation. After
    stagnation, alpha rises sharply as the fuel is reheated by alpha
    particles and compression -- this is expected and physical.

    For calibration comparisons across flux limiter runs, evaluate alpha
    at the time of peak implosion velocity (data.peak_velocity_index) or
    at a fixed reference time during the foot pulse to assess first-shock
    entropy deposition.

    Examples
    --------
    >>> from helios_postprocess.adiabat_history import compute_adiabat_history
    >>> result = compute_adiabat_history(
    ...     ion_pressure=data.ion_pressure,
    ...     rad_pressure=data.rad_pressure,
    ...     mass_density=data.mass_density,
    ...     zone_mass=data.zone_mass,
    ...     region_interfaces_indices=data.region_interfaces_indices,
    ...     time=data.time * 1e9,
    ... )
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(result['time'], result['adiabat'])
    >>> plt.xlabel('Time [ns]')
    >>> plt.ylabel('Mass-averaged adiabat')
    >>> plt.title('Cold fuel adiabat history')
    >>> plt.show()
    """
    n_times = len(time)
    ri = region_interfaces_indices

    adiabat         = np.full(n_times, np.nan)
    adiabat_min     = np.full(n_times, np.nan)
    adiabat_max     = np.full(n_times, np.nan)
    P_mean_Gbar     = np.full(n_times, np.nan)
    rho_mean        = np.full(n_times, np.nan)
    z_start_arr     = np.zeros(n_times, dtype=int)
    z_end_arr       = np.zeros(n_times, dtype=int)

    for i in range(n_times):
        try:
            z_start = int(ri[i, fuel_region_col])
            z_end   = int(ri[i, fuel_outer_col])
            z_start_arr[i] = z_start
            z_end_arr[i]   = z_end

            if z_end <= z_start:
                continue

            # Total pressure in Mbar
            p_tot = ion_pressure[i, z_start:z_end].copy()
            if rad_pressure is not None:
                p_tot += rad_pressure[i, z_start:z_end]
            p_Mbar = p_tot * 1e-5  # J/cm³ -> Mbar

            rho  = mass_density[i, z_start:z_end]
            mass = zone_mass[i, z_start:z_end]

            if np.all(rho <= 0) or np.all(mass <= 0):
                continue

            # Fermi pressure
            with np.errstate(divide='ignore', invalid='ignore'):
                p_fermi = 2.17 * (rho / rho0) ** (5.0 / 3.0)  # Mbar
                alpha_zones = p_Mbar / p_fermi
                alpha_zones[~np.isfinite(alpha_zones)] = np.nan

            valid = np.isfinite(alpha_zones) & (mass > 0)
            if not np.any(valid):
                continue

            adiabat[i]     = np.average(alpha_zones[valid], weights=mass[valid])
            adiabat_min[i] = np.nanmin(alpha_zones)
            adiabat_max[i] = np.nanmax(alpha_zones)
            P_mean_Gbar[i] = np.average(p_tot[valid] * 1e-8, weights=mass[valid])
            rho_mean[i]    = np.average(rho[valid], weights=mass[valid])

        except Exception:
            continue

    return {
        'time':           time,
        'adiabat':        adiabat,
        'adiabat_min':    adiabat_min,
        'adiabat_max':    adiabat_max,
        'P_fuel_mean_Gbar': P_mean_Gbar,
        'rho_fuel_mean':  rho_mean,
        'z_start':        z_start_arr,
        'z_end':          z_end_arr,
    }


def print_adiabat_history(result: Dict[str, np.ndarray],
                          t_start: float = 0.0,
                          t_end: float = 999.0,
                          stride: int = 1) -> None:
    """
    Print adiabat history table to stdout.

    Parameters
    ----------
    result : dict
        Output of compute_adiabat_history().
    t_start, t_end : float
        Time window [ns] to print.
    stride : int
        Print every Nth timestep.
    """
    time = result['time']
    print(f"{'Time(ns)':>10} {'Adiabat':>9} {'Alpha_min':>10} {'Alpha_max':>10} "
          f"{'P_mean(Gbar)':>13} {'rho_mean':>10}")
    print("-" * 67)
    for i in range(0, len(time), stride):
        t = time[i]
        if t < t_start or t > t_end:
            continue
        a  = result['adiabat'][i]
        mn = result['adiabat_min'][i]
        mx = result['adiabat_max'][i]
        p  = result['P_fuel_mean_Gbar'][i]
        rh = result['rho_fuel_mean'][i]
        if np.isnan(a):
            continue
        print(f"{t:10.3f} {a:9.3f} {mn:10.3f} {mx:10.3f} {p:13.4f} {rh:10.4f}")
