"""
Burn-Averaged Metrics Module
=============================

Calculates burn-averaged (temporal) quantities for comparison with published
ICF performance metrics.

This module provides TIME-HISTORY burn-averaging using the 2D arrays stored
on ICFRunData (after build_run_data + ICFAnalyzer have run), complementing
the SPATIAL neutron-averaging in icf_analysis.py.

Important
---------
EXODUS files sample only a fraction of the actual simulation timesteps.
Therefore Helios's own time-integrated quantities (neutron count, deposited
energy, etc.) are far more accurate than anything re-integrated from EXODUS
data.  This module uses the Helios-computed fusion yield and laser energy
for absolute values, and only uses the EXODUS-sampled burn rate profile as
a *weighting function* for burn-averaging.

Workflow
-------
    >>> from helios_postprocess import HeliosRun
    >>> from helios_postprocess.data_builder import build_run_data
    >>> from helios_postprocess.icf_analysis import ICFAnalyzer
    >>> from helios_postprocess.burn_averaged_metrics import (
    ...     extract_histories_from_run_data,
    ...     calculate_burn_averaged_metrics,
    ...     compare_with_published,
    ... )
    >>>
    >>> run = HeliosRun('sim.exo', verbose=True)
    >>> data = build_run_data(run, time_unit='s')
    >>> run.close()
    >>>
    >>> analyzer = ICFAnalyzer(data)
    >>> analyzer.analyze_drive_phase()
    >>> analyzer.analyze_stagnation_phase()
    >>> analyzer.analyze_burn_phase()
    >>> analyzer.compute_performance_metrics()
    >>>
    >>> histories = extract_histories_from_run_data(data)
    >>> metrics = calculate_burn_averaged_metrics(histories)
    >>> print(compare_with_published(metrics, published_metrics, laser_energy_MJ=4.0))

Author: Prof T
Date: November 2025 (original), March 2026 (refactored for ICFRunData pipeline)
"""

import numpy as np
from scipy.integrate import simpson
from typing import Dict, Optional
import warnings
import logging

logger = logging.getLogger(__name__)

try:
    from .core import sigma_v_DT
except ImportError:
    from core import sigma_v_DT


# ── Burn rate calculation ────────────────────────────────────────────────────

def calculate_burn_rate_from_sim(temperature_keV: np.ndarray,
                                density_gcc: np.ndarray,
                                ion_fraction: float = 0.5) -> np.ndarray:
    """
    Calculate DT fusion reaction rate from temperature and density.

    This is used as a *weighting function* for burn-averaging, NOT for
    computing absolute yield (which comes from Helios's own integrals).

    Parameters
    ----------
    temperature_keV : np.ndarray
        Ion temperature in keV at each time step, shape (n_times,)
    density_gcc : np.ndarray
        Mass density in g/cm³ at each time step, shape (n_times,)
    ion_fraction : float, optional
        Fraction of each ion species (default 0.5 for 50-50 DT)

    Returns
    -------
    np.ndarray
        Fusion reaction rate in reactions/(cm³·s), shape (n_times,)
    """
    amu_to_g = 1.66054e-24
    m_DT = 2.5 * amu_to_g   # Average DT ion mass

    n_total = density_gcc / m_DT
    n_D = n_total * ion_fraction
    n_T = n_total * ion_fraction

    sv = sigma_v_DT(temperature_keV)

    return n_D * n_T * sv


# ── History extraction from ICFRunData ───────────────────────────────────────

def extract_histories_from_run_data(data) -> Dict:
    """
    Compute per-timestep hot-spot averaged quantities from ICFRunData 2D arrays.

    Uses region_interfaces_indices to identify the hot-spot boundary at each
    timestep, then computes mass-averaged temperature, pressure, and density
    within the hot spot.  Also extracts cold-fuel areal density from the
    cumulative areal_density_vs_time array.

    Also passes through Helios-computed integral quantities and scalar
    implosion metrics from ICFAnalyzer.

    Parameters
    ----------
    data : ICFRunData
        Populated dataclass — must have had build_run_data() and
        ICFAnalyzer (at least through compute_performance_metrics()) run on it.

    Returns
    -------
    dict
        Keys for calculate_burn_averaged_metrics().
    """
    # ── Validate required data ──
    _check_required(data, [
        'time', 'ion_temperature', 'mass_density',
        'ion_pressure', 'rad_pressure',
        'zone_boundaries', 'zone_mass',
        'region_interfaces_indices',
    ])

    time_ns = np.asarray(data.time)
    n_times = len(time_ns)

    ri = data.region_interfaces_indices          # (n_times, n_interfaces)

    # Hot-spot averaging requires a defined gas/fuel interface (>=2 regions).
    # Single-region targets (e.g. CH-only slab/sphere flux-limiter tests) have
    # no hot spot to extract histories for -- return empty result.
    if ri is None or ri.shape[1] < 2:
        # Single-region target (e.g. CH-only slab/sphere flux-limiter test).
        # Return a dict with the keys calculate_burn_averaged_metrics() expects,
        # filled with NaN/zero so its Simpson-weighted integrals evaluate to
        # 0 / NaN cleanly without crashing.
        logger.info("Hot-spot history extraction: single-region target -- returning NaN stub")
        nan_arr  = np.full(n_times, np.nan)
        zero_arr = np.zeros(n_times)
        return {
            'time_ns':            np.asarray(data.time),
            'temperature_keV':    zero_arr,    # zero -> burn_rate = 0 -> clean zero integral
            'pressure_Gbar':      nan_arr,
            'density_gcc':        zero_arr,
            'areal_density_gcm2': nan_arr,
            'radius_um':          nan_arr,
            'CR_max':             0.0,
            'energy_output_MJ':   0.0,
            'laser_energy_MJ':    getattr(data, 'laser_energy_delivered_MJ', 0.0),
            'target_gain':        0.0,
            'stagnation_time_ns': 0.0,
            'single_region':      True,
        }

    T_ion = data.ion_temperature                 # (n_times, n_zones) eV
    rho   = data.mass_density                    # (n_times, n_zones) g/cm³
    P_ion = data.ion_pressure                    # (n_times, n_zones) J/cm³
    P_rad = data.rad_pressure                    # (n_times, n_zones) J/cm³
    P_tot = P_ion + P_rad                        # total pressure
    zmass = data.zone_mass                       # (n_times, n_zones) g
    zbnd  = data.zone_boundaries                 # (n_times, n_nodes) cm
    vel   = data.velocity                        # (n_times, n_nodes) cm/s

    # Cumulative areal density (optional — for cold fuel ρR)
    cumul_rhoR = data.areal_density_vs_time      # (n_times, n_zones+1) or None

    # Output arrays
    temperatures = np.zeros(n_times)
    pressures    = np.zeros(n_times)
    densities    = np.zeros(n_times)
    radii_um     = np.zeros(n_times)
    volumes      = np.zeros(n_times)
    masses_mg    = np.zeros(n_times)
    rhoR_cf      = np.zeros(n_times)

    logger.info(f"Extracting hot-spot histories ({n_times} timesteps)...")

    for t in range(n_times):
        # Hot-spot outer boundary = first region interface (node index)
        hs_bnd = int(ri[t, 0])

        if hs_bnd <= 0:
            continue   # no hot spot at this timestep

        # Hot-spot zones: 0 to hs_bnd-1
        hs = slice(0, hs_bnd)

        m = zmass[t, hs]                # zone masses in hot spot
        m_total = np.sum(m)

        if m_total <= 0:
            continue

        # Mass-averaged ion temperature → keV
        temperatures[t] = np.sum(T_ion[t, hs] * m) / m_total / 1000.0

        # Mass-averaged total pressure → Gbar
        pressures[t] = np.sum(P_tot[t, hs] * m) / m_total * 1e-8

        # Mass-averaged density
        densities[t] = np.sum(rho[t, hs] * m) / m_total

        # Hot-spot outer radius (node position at hs_bnd)
        r_cm = zbnd[t, hs_bnd]
        radii_um[t] = r_cm * 1e4            # cm → μm

        # Spherical volume
        volumes[t] = (4.0 / 3.0) * np.pi * r_cm**3

        # Hot-spot mass
        masses_mg[t] = m_total * 1e3         # g → mg

        # Cold-fuel areal density
        if cumul_rhoR is not None:
            # fuel/ablator boundary is the second-to-last interface
            fuel_bnd = int(ri[t, -2]) if ri.shape[1] >= 2 else hs_bnd
            rhoR_cf[t] = cumul_rhoR[t, fuel_bnd] - cumul_rhoR[t, hs_bnd]

    logger.info("Hot-spot history extraction complete")

    # ── Helios-computed integral quantities (authoritative) ──
    energy_output_MJ = getattr(data, 'energy_output', 0.0)
    laser_energy_MJ  = getattr(data, 'laser_energy', 0.0)
    target_gain      = getattr(data, 'target_gain', 0.0)

    # ── Scalar implosion metrics from ICFAnalyzer ──
    peak_vimp_kms       = abs(getattr(data, 'peak_implosion_velocity', 0.0))
    peak_vimp_kms_cr15  = abs(getattr(data, 'peak_implosion_velocity_at_cr15', 0.0))
    impl_vel_rhino_kms  = abs(getattr(data, 'implosion_velocity_rhino_kms', 0.0))
    adiabat             = getattr(data, 'adiabat_mass_averaged_ice', 0.0)
    adiabat_cr15        = getattr(data, 'adiabat_mass_averaged_ice_cr15', 0.0)
    adiabat_min_rhino   = getattr(data, 'adiabat_min_rhino', 0.0)
    # RHINO-formula parallel adiabats (proper degenerate electron gas
    # Fermi pressure with actual n_e per zone). Same zone-selection and
    # time-selection as the Lindl-convention adiabats above; just a
    # different Fermi-pressure denominator. Typically ~14x larger at
    # ICF densities.
    adiabat_rhino_formula        = getattr(data, 'adiabat_mass_averaged_ice_rhino_formula', 0.0)
    adiabat_cr15_rhino_formula   = getattr(data, 'adiabat_mass_averaged_ice_cr15_rhino_formula', 0.0)
    adiabat_breakout_rhino_formula = getattr(data, 'adiabat_at_breakout_rhino_formula', 0.0)
    adiabat_at_breakout            = getattr(data, 'adiabat_at_breakout', 0.0)
    # RHINO diagnostics suite (June 2026)
    t_max_shell_velocity_rhino_ns   = getattr(data, 't_max_shell_velocity_rhino_ns', 0.0)
    stag_time_rhino_ns              = getattr(data, 'stag_time_rhino_ns', 0.0)
    assembled_mass_rhino_mg         = getattr(data, 'assembled_mass_rhino_mg', 0.0)
    burn_fraction_rhino             = getattr(data, 'burn_fraction_rhino', 0.0)
    ablation_pressure_at_cr_3p5_Mbar = getattr(data, 'ablation_pressure_at_cr_3p5_Mbar', 0.0)
    t_at_cr_3p5_ns                  = getattr(data, 't_at_cr_3p5_ns', 0.0)
# IFAR variants (May 2026 — multi-variant patch)
    # - data.ifar          : compound shell at peak v_imp (legacy default)
    # - data.ifar_ice      : DT ice only at peak v_imp
    # - data.ifar_cr15     : compound shell at CR=1.5  (Thomas/Vulcan convention)
    # - data.ifar_ice_cr15 : DT ice only at CR=1.5     (artifact — see _compute_ifar)
    ifar_compound_pv = getattr(data, 'ifar', 0.0)
    ifar_ice_pv      = getattr(data, 'ifar_ice', None)
    ifar_cr15        = getattr(data, 'ifar_cr15', None)
    ifar_ice_cr15    = getattr(data, 'ifar_ice_cr15', None)

    # For published comparison: prefer compound at CR=1.5 (Thomas convention).
    # Fall back to legacy peak-v compound for older runs without the variants.
    ifar = ifar_cr15 if ifar_cr15 is not None else ifar_compound_pv

    # ── Derived implosion metrics ──

    # Absorbed energy from laser_energy_deposited (cumulative, in J)
    # This is the total energy coupled into the target — the denominator
    # for hydrodynamic efficiency.
    led = getattr(data, 'laser_energy_deposited', None)
    E_absorbed_J = 0.0
    if led is not None:
        E_absorbed_J = np.max(np.sum(led, axis=-1)) if led.ndim > 1 else np.max(led)

    # Fraction absorbed: deposited / delivered
    # laser_power_delivered integrated gives total incident energy;
    # laser_energy_deposited gives absorbed energy.
    fraction_absorbed = 0.0
    lpd = getattr(data, 'laser_power_delivered', None)
    if lpd is not None and E_absorbed_J > 0:
        # Integrate delivered power over time
        E_delivered_J = np.trapezoid(lpd, x=time_ns * 1e-9)
        if E_delivered_J > 0:
            fraction_absorbed = E_absorbed_J / E_delivered_J * 100.0

    # In-flight kinetic energy: INWARD-MOVING shell only, max over time.
    # η_hydro = max(KE_inward) / E_absorbed
    # KE computed at every timestep; only zones with v < 0 (imploding) count.
    inflight_KE_kJ = 0.0
    if vel is not None and zmass is not None:
        n_zones = zmass.shape[1]
        if vel.shape[1] > n_zones:
            # Node-centered velocity → zone-center average
            v_zone = 0.5 * (vel[:, :n_zones] + vel[:, 1:n_zones+1])
        else:
            v_zone = vel[:, :n_zones]

        # At each timestep, sum KE of inward-moving zones only (v < 0)
        KE_inward = np.zeros(n_times)
        for t in range(n_times):
            inward = v_zone[t, :] < 0           # imploding zones
            if np.any(inward):
                KE_inward[t] = 0.5 * np.sum(
                    zmass[t, inward] * v_zone[t, inward]**2
                )  # g·(cm/s)² = erg

        # Convert max KE from erg to kJ: 1 erg = 1e-7 J, 1 kJ = 1e3 J
        inflight_KE_kJ = np.max(KE_inward) * 1e-10  # erg → kJ

    # Hydrodynamic efficiency = max(KE_inward) / E_absorbed
    hydro_efficiency_pct = 0.0
    if inflight_KE_kJ > 0 and E_absorbed_J > 0:
        E_absorbed_kJ = E_absorbed_J * 1e-3  # J → kJ
        hydro_efficiency_pct = inflight_KE_kJ / E_absorbed_kJ * 100.0

    # Imploded DT mass at stagnation
    imploded_DT_mass_mg = 0.0
    stag_idx = getattr(data, '_stag_index', None)
    if stag_idx is None:
        # Find stagnation index from stag_time
        stag_time = getattr(data, 'stag_time', 0.0)
        if stag_time > 0:
            stag_idx = np.argmin(np.abs(time_ns - stag_time))
    if stag_idx is not None and stag_idx > 0:
        # Imploded mass = all DT zones inside the ablation front at stagnation.
        # Prefer zone-index method (consistent with icf_analysis.py mass fractions)
        # over radius method, which is sensitive to smoothing artefacts.
        abl_indices = getattr(data, 'ablation_front_indices', None)
        ri = getattr(data, 'region_interfaces_indices', None)
        fuel_bnd = int(ri[0, -2]) if ri is not None and ri.shape[1] >= 2 else zmass.shape[1]
        if abl_indices is not None and abl_indices[stag_idx] > 0:
            abl_idx = int(abl_indices[stag_idx])
            imploded_DT_mass_mg = float(np.sum(zmass[stag_idx, :min(abl_idx+1, fuel_bnd)])) * 1e3
        else:
            # Fallback: radius-based (less reliable at stagnation)
            abl_front = getattr(data, 'ablation_front_radius', None)
            if abl_front is not None:
                r_abl = abl_front[stag_idx]
                for z in range(zmass.shape[1]):
                    if zbnd[stag_idx, z + 1] <= r_abl * 1.01:
                        imploded_DT_mass_mg += zmass[stag_idx, z]
                imploded_DT_mass_mg *= 1e3  # g → mg

    return {
        # Time histories
        'time_ns':            time_ns,
        'temperature_keV':    temperatures,
        'pressure_Gbar':      pressures,
        'density_gcc':        densities,
        'areal_density_gcm2': rhoR_cf,
        'radius_um':          radii_um,
        'volume_cm3':         volumes,
        'mass_mg':            masses_mg,
        'initial_radius_um':  radii_um[0],
        # Helios-computed integrals (authoritative)
        'energy_output_MJ':   energy_output_MJ,
        'stagnation_time_ns':  getattr(data, 'stag_time', 0.0),
        'laser_energy_MJ':    laser_energy_MJ,
        'target_gain':        target_gain,
        # Implosion metrics
        'peak_velocity_kms':     peak_vimp_kms,
        'peak_velocity_kms_cr15': peak_vimp_kms_cr15,
        'implosion_velocity_rhino_kms': impl_vel_rhino_kms,
        'adiabat':               adiabat,
        'adiabat_cr15':          adiabat_cr15,
        'adiabat_min_rhino':     adiabat_min_rhino,
        'adiabat_rhino_formula': adiabat_rhino_formula,
        'adiabat_cr15_rhino_formula': adiabat_cr15_rhino_formula,
        'adiabat_breakout_rhino_formula': adiabat_breakout_rhino_formula,
        'adiabat_at_breakout':           adiabat_at_breakout,
        # RHINO diagnostics suite
        't_max_shell_velocity_rhino_ns': t_max_shell_velocity_rhino_ns,
        'stag_time_rhino_ns':            stag_time_rhino_ns,
        'assembled_mass_rhino_mg':       assembled_mass_rhino_mg,
        'burn_fraction_rhino':           burn_fraction_rhino,
        'ablation_pressure_at_cr_3p5_Mbar': ablation_pressure_at_cr_3p5_Mbar,
        't_at_cr_3p5_ns':                t_at_cr_3p5_ns,
        'ifar':                  ifar,                # comparison-ready (CR=1.5 if available)
        'ifar_compound_peak_v':  ifar_compound_pv,    # legacy
        'ifar_ice_peak_v':       ifar_ice_pv,         # diagnostic
        'ifar_compound_cr15':    ifar_cr15,           # Thomas convention
        'ifar_ice_cr15':         ifar_ice_cr15,       # artifact (see code comments)
        'fraction_absorbed_pct': fraction_absorbed,
        'P_hs_ignition_Gbar':    getattr(data, 'ignition_hs_pressure', 0.0),
        'hs_radius_ignition_um': getattr(data, 'ignition_hs_radius',   0.0) * 1e4,
        'T_ion_onaxis_ignition_keV': getattr(data, 'ignition_T_ion_onaxis_keV', 0.0),
        'T_ion_hs_avg_ignition_keV': getattr(data, 'ignition_T_ion_hs_avg_keV', 0.0),
        'peak_total_rhoR':       getattr(data, 'peak_total_rhoR',       0.0),
        'peak_hs_rhoR_T_mask':   getattr(data, 'peak_hs_rhoR_T_mask',   0.0),
        'peak_density_at_ignition': getattr(data, 'peak_density_at_ignition', 0.0),
        'peak_density_at_ignition_is_stagnation':
            getattr(data, 'peak_density_at_ignition_is_stagnation', False),
        'inflight_KE_kJ':       inflight_KE_kJ,
        'hydro_efficiency_pct':  hydro_efficiency_pct,
        'imploded_DT_mass_mg':   imploded_DT_mass_mg,
        'CR_max':                getattr(data, 'comp_ratio', 0.0),
        # Laser intensity metrics (from analyze_laser_intensity); 0.0 if skipped
        'I_at_crit_peak_Wcm2':    getattr(data, 'I_at_crit_peak', None) or 0.0,
        'I_grid_outer_peak_Wcm2': getattr(data, 'I_grid_outer_peak', None) or 0.0,
        # Shock-train breakouts at gas/ice interface (from _compute_shock_train).
        # NaN sentinel means "not detected"; comparison code treats it as a
        # missing value rather than a real 0.0 ns measurement.
        't_foot_shock_breakout_ns': float(getattr(data, 't_foot_shock_ns', float('nan'))),
        't_ramp_shock_breakout_ns': float(getattr(data, 't_ramp_shock_ns', float('nan'))),
        't_peak_shock_breakout_ns': float(getattr(data, 't_peak_shock_ns', float('nan'))),
        # ── Will list extensions (June 4 2026) ───────────────────
        'cr_outer':                       getattr(data, 'cr_outer', 0.0),
        'laser_overlapped_intensity_Wcm2': getattr(data, 'laser_overlapped_intensity_Wcm2', 0.0),
        'shell_mass_at_stagnation_mg':    getattr(data, 'shell_mass_at_stagnation_mg', 0.0),
        'hot_spot_mass_at_stagnation_mg': getattr(data, 'hot_spot_mass_at_stagnation_mg', 0.0),
        'stag_time_areal_density':        getattr(data, 'stag_time_areal_density', 0.0),
        'stag_time_fuel_areal_density':   getattr(data, 'stag_time_fuel_areal_density', 0.0),
        'neutron_ave_electron_temperature': getattr(data, 'neutron_ave_electron_temperature', 0.0),
        # ── Will shell-convention scalars (June 2026, Phys. Conv. §18b) ──
        'shell_mass_will_at_stagnation_mg': getattr(data, 'shell_mass_will_at_stagnation_mg', 0.0),
        'adiabat_mass_avg_will_cr15':       getattr(data, 'adiabat_mass_avg_will_cr15', 0.0),
        'sound_speed_shell_will_cr15_kms':  getattr(data, 'sound_speed_shell_will_cr15_kms', 0.0),
        't_will_shell_cr15_ns':             getattr(data, 't_will_shell_cr15_ns', 0.0),
        # ── Will ignition-product timing (June 2026, Phys. Conv. §18d/e) ──
        't_peak_rhoR_Ti_ns':                       getattr(data, 't_peak_rhoR_Ti_ns', 0.0),
        'peak_rhoR_Ti_gcm2_keV':                   getattr(data, 'peak_rhoR_Ti_gcm2_keV', 0.0),
        'rhoR_hs_at_peak_rhoR_Ti_gcm2':            getattr(data, 'rhoR_hs_at_peak_rhoR_Ti_gcm2', 0.0),
        'T_i_volume_avg_at_peak_rhoR_Ti_keV':      getattr(data, 'T_i_volume_avg_at_peak_rhoR_Ti_keV', 0.0),
        'T_i_mass_avg_at_peak_rhoR_Ti_keV':        getattr(data, 'T_i_mass_avg_at_peak_rhoR_Ti_keV', 0.0),
        't_peak_Phs_Rhs_ns':                       getattr(data, 't_peak_Phs_Rhs_ns', 0.0),
        'peak_Phs_Rhs_Gbar_um':                    getattr(data, 'peak_Phs_Rhs_Gbar_um', 0.0),
        'R_hs_at_peak_Phs_Rhs_um':                 getattr(data, 'R_hs_at_peak_Phs_Rhs_um', 0.0),
        'P_hs_volume_avg_at_peak_Phs_Rhs_Gbar':    getattr(data, 'P_hs_volume_avg_at_peak_Phs_Rhs_Gbar', 0.0),
        'P_hs_mass_avg_at_peak_Phs_Rhs_Gbar':      getattr(data, 'P_hs_mass_avg_at_peak_Phs_Rhs_Gbar', 0.0),
    }


def _check_required(data, attr_names):
    """Raise ValueError if any required attribute is None."""
    missing = [a for a in attr_names if getattr(data, a, None) is None]
    if missing:
        raise ValueError(
            f"ICFRunData is missing required attributes: {missing}\n"
            f"Make sure build_run_data() and ICFAnalyzer have been run."
        )


# ── Legacy compatibility wrapper ─────────────────────────────────────────────

def extract_hot_spot_histories(run_or_data, T_threshold=1000.0,
                               time_indices=None) -> Dict:
    """
    Compatibility wrapper — routes to extract_histories_from_run_data().
    """
    if hasattr(run_or_data, 'ion_temperature') and hasattr(run_or_data, 'zone_mass'):
        return extract_histories_from_run_data(run_or_data)

    raise TypeError(
        "extract_hot_spot_histories() no longer accepts HeliosRun directly.\n"
        "Use the pipeline instead:\n"
        "  data = build_run_data(run, time_unit='s')\n"
        "  analyzer = ICFAnalyzer(data)\n"
        "  analyzer.analyze_stagnation_phase()\n"
        "  analyzer.analyze_burn_phase()\n"
        "  histories = extract_histories_from_run_data(data)\n"
    )


# ── Burn-averaged metrics ────────────────────────────────────────────────────

def calculate_burn_averaged_metrics(histories: Dict,
                                    ion_fraction: float = 0.5) -> Dict:
    """
    Calculate burn-averaged metrics from time histories.

    Burn-averaging: <Q> = ∫ Q(t) · w(t) dt / ∫ w(t) dt
    where w(t) is a burn-rate weighting function ∝ ρ² · <σv>(T).

    The weighting function is computed from the EXODUS-sampled data and is
    used only for averaging.  Absolute yield and gain come from Helios's own
    time-integrated quantities, which use every simulation timestep and are
    therefore far more accurate than re-integrating from the sampled EXODUS
    output.

    Parameters
    ----------
    histories : dict
        Dictionary from extract_histories_from_run_data().
    ion_fraction : float, optional
        Ion fraction for burn rate weighting (default 0.5 for 50-50 DT).

    Returns
    -------
    dict
        Burn-averaged quantities, peak values, Helios integrals,
        implosion metrics, and burn-rate weighting profile.
    """
    time_ns    = histories['time_ns']
    temp_keV   = histories['temperature_keV']
    pres_Gbar  = histories['pressure_Gbar']
    dens_gcc   = histories['density_gcc']
    rhoR_gcm2  = histories['areal_density_gcm2']
    radius_um  = histories['radius_um']

    time_s = time_ns * 1e-9

    # ── Burn-rate weighting function ──
    burn_rate = calculate_burn_rate_from_sim(temp_keV, dens_gcc, ion_fraction)

    # ── Burn-averaged quantities ──
    def _burn_avg(quantity):
        num = simpson(quantity * burn_rate, x=time_s)
        den = simpson(burn_rate, x=time_s)
        return num / den if den > 0 else 0.0

    T_burn_avg    = _burn_avg(temp_keV)
    P_burn_avg    = _burn_avg(pres_Gbar)
    rho_burn_avg  = _burn_avg(dens_gcc)
    rhoR_burn_avg = _burn_avg(rhoR_gcm2)

    # ── Convergence ratio ──
    # Use CR computed in icf_analysis.py: R0 / Rf where R0 = initial inner
    # shell radius and Rf = ablation front radius at peak implosion velocity.
    CR_max = histories.get('CR_max', 0.0)

    # ── Helios-computed integrals ──
    yield_MJ        = histories.get('energy_output_MJ', 0.0)
    laser_energy_MJ = histories.get('laser_energy_MJ', 0.0)
    target_gain     = histories.get('target_gain', 0.0)

    # ── Normalized burn fraction ──
    br_sum = np.sum(burn_rate)
    burn_fraction = burn_rate / br_sum if br_sum > 0 else burn_rate

    valid_r = radius_um[radius_um > 0]
    min_radius = np.min(valid_r) if len(valid_r) > 0 else 0.0

    return {
        # Burn-averaged values
        'stagnation_time_ns': histories.get('stagnation_time_ns', 0.0),
        'T_burn_avg':       T_burn_avg,
        'P_burn_avg':       P_burn_avg,
        'rho_burn_avg':     rho_burn_avg,
        'rhoR_burn_avg':    rhoR_burn_avg,
        # Peak values
        'T_peak':           np.max(temp_keV),
        'P_peak':           np.max(pres_Gbar),
        'rho_peak':         np.max(dens_gcc),
        'rhoR_peak':        np.max(rhoR_gcm2),
        # Convergence
        'CR_max':           CR_max,
        'min_radius_um':    min_radius,
        # Helios-computed integrals (authoritative)
        'yield_MJ':         yield_MJ,
        'laser_energy_MJ':  laser_energy_MJ,
        'target_gain':      target_gain,
        # Implosion metrics (pass through from histories)
        'peak_velocity_kms':     histories.get('peak_velocity_kms', 0.0),
        'peak_velocity_kms_cr15': histories.get('peak_velocity_kms_cr15', 0.0),
        'implosion_velocity_rhino_kms': histories.get('implosion_velocity_rhino_kms', 0.0),
        'adiabat':               histories.get('adiabat', 0.0),
        'adiabat_cr15':          histories.get('adiabat_cr15', 0.0),
        'adiabat_min_rhino':     histories.get('adiabat_min_rhino', 0.0),
        'adiabat_rhino_formula': histories.get('adiabat_rhino_formula', 0.0),
        'adiabat_cr15_rhino_formula': histories.get('adiabat_cr15_rhino_formula', 0.0),
        'adiabat_breakout_rhino_formula': histories.get('adiabat_breakout_rhino_formula', 0.0),
        'adiabat_at_breakout':           histories.get('adiabat_at_breakout', 0.0),
        # RHINO diagnostics suite
        't_max_shell_velocity_rhino_ns': histories.get('t_max_shell_velocity_rhino_ns', 0.0),
        'stag_time_rhino_ns':            histories.get('stag_time_rhino_ns', 0.0),
        'assembled_mass_rhino_mg':       histories.get('assembled_mass_rhino_mg', 0.0),
        'burn_fraction_rhino':           histories.get('burn_fraction_rhino', 0.0),
        'ablation_pressure_at_cr_3p5_Mbar': histories.get('ablation_pressure_at_cr_3p5_Mbar', 0.0),
        't_at_cr_3p5_ns':                histories.get('t_at_cr_3p5_ns', 0.0),
        'ifar':                  histories.get('ifar', 0.0),
        'fraction_absorbed_pct': histories.get('fraction_absorbed_pct', 0.0),
        'P_hs_ignition_Gbar':    histories.get('P_hs_ignition_Gbar',    0.0),
        'hs_radius_ignition_um': histories.get('hs_radius_ignition_um', 0.0),
        'T_ion_onaxis_ignition_keV': histories.get('T_ion_onaxis_ignition_keV', 0.0),
        'T_ion_hs_avg_ignition_keV': histories.get('T_ion_hs_avg_ignition_keV', 0.0),
        'peak_total_rhoR':       histories.get('peak_total_rhoR',       0.0),
        'peak_hs_rhoR_T_mask':   histories.get('peak_hs_rhoR_T_mask',   0.0),
        'peak_density_at_ignition': histories.get('peak_density_at_ignition', 0.0),
        'peak_density_at_ignition_is_stagnation':
            histories.get('peak_density_at_ignition_is_stagnation', False),
        'inflight_KE_kJ':       histories.get('inflight_KE_kJ', 0.0),
        'hydro_efficiency_pct':  histories.get('hydro_efficiency_pct', 0.0),
        'imploded_DT_mass_mg':   histories.get('imploded_DT_mass_mg', 0.0),
        # Laser intensity metrics (pass through from histories)
        'I_at_crit_peak_Wcm2':    histories.get('I_at_crit_peak_Wcm2',    0.0),
        'I_grid_outer_peak_Wcm2': histories.get('I_grid_outer_peak_Wcm2', 0.0),
        # Shock-train breakouts (pass through; NaN if not detected)
        't_foot_shock_breakout_ns': histories.get('t_foot_shock_breakout_ns', float('nan')),
        't_ramp_shock_breakout_ns': histories.get('t_ramp_shock_breakout_ns', float('nan')),
        't_peak_shock_breakout_ns': histories.get('t_peak_shock_breakout_ns', float('nan')),
        # Will list extensions (pass through)
        'cr_outer':                          histories.get('cr_outer', 0.0),
        'laser_overlapped_intensity_Wcm2':   histories.get('laser_overlapped_intensity_Wcm2', 0.0),
        'shell_mass_at_stagnation_mg':       histories.get('shell_mass_at_stagnation_mg', 0.0),
        'hot_spot_mass_at_stagnation_mg':    histories.get('hot_spot_mass_at_stagnation_mg', 0.0),
        'stag_time_areal_density':           histories.get('stag_time_areal_density', 0.0),
        'stag_time_fuel_areal_density':      histories.get('stag_time_fuel_areal_density', 0.0),
        'neutron_ave_electron_temperature':  histories.get('neutron_ave_electron_temperature', 0.0),
        # Will shell-convention (June 2026)
        'shell_mass_will_at_stagnation_mg':  histories.get('shell_mass_will_at_stagnation_mg', 0.0),
        'adiabat_mass_avg_will_cr15':        histories.get('adiabat_mass_avg_will_cr15', 0.0),
        'sound_speed_shell_will_cr15_kms':   histories.get('sound_speed_shell_will_cr15_kms', 0.0),
        # Will ignition-product timing (pass through)
        't_peak_rhoR_Ti_ns':                    histories.get('t_peak_rhoR_Ti_ns', 0.0),
        'peak_rhoR_Ti_gcm2_keV':                histories.get('peak_rhoR_Ti_gcm2_keV', 0.0),
        'rhoR_hs_at_peak_rhoR_Ti_gcm2':         histories.get('rhoR_hs_at_peak_rhoR_Ti_gcm2', 0.0),
        'T_i_volume_avg_at_peak_rhoR_Ti_keV':   histories.get('T_i_volume_avg_at_peak_rhoR_Ti_keV', 0.0),
        'T_i_mass_avg_at_peak_rhoR_Ti_keV':     histories.get('T_i_mass_avg_at_peak_rhoR_Ti_keV', 0.0),
        't_peak_Phs_Rhs_ns':                    histories.get('t_peak_Phs_Rhs_ns', 0.0),
        'peak_Phs_Rhs_Gbar_um':                 histories.get('peak_Phs_Rhs_Gbar_um', 0.0),
        'R_hs_at_peak_Phs_Rhs_um':              histories.get('R_hs_at_peak_Phs_Rhs_um', 0.0),
        'P_hs_volume_avg_at_peak_Phs_Rhs_Gbar': histories.get('P_hs_volume_avg_at_peak_Phs_Rhs_Gbar', 0.0),
        'P_hs_mass_avg_at_peak_Phs_Rhs_Gbar':   histories.get('P_hs_mass_avg_at_peak_Phs_Rhs_Gbar', 0.0),
        # Burn-rate weighting profile
        'burn_fraction':    burn_fraction,
        'burn_rate':        burn_rate,
    }


# ── Comparison with published data ───────────────────────────────────────────

def compare_with_published(sim_metrics: Dict,
                           published_metrics: Dict,
                           laser_energy_MJ: Optional[float] = None) -> str:
    """
    Generate comparison table between simulation and published results.

    Renders six columns:
        Metric | Sim | Cluster ± unc | Δ_cluster | HYDRA ± unc | Δ_HYDRA

    The HYDRA column comes from per-code keys in the published JSON
    (e.g. 'stagnation_time_HYDRA_ns', 'peak_rhoR_total_HYDRA_gcm2').
    Rows whose per-code key is absent render '—' in the HYDRA cell.
    This is backwards-compatible: JSONs with no HYDRA keys render
    the cluster columns exactly as before.

    Parameters
    ----------
    sim_metrics : dict
        From calculate_burn_averaged_metrics().
    published_metrics : dict
        Published values: {key: (value, uncertainty)}.
        Supported cluster-level keys:
          Burn-averaged: 'T_hs', 'P_hs', 'rhoR_cf', 'CR_max', 'yield', 'gain'
          Implosion: 'peak_velocity_kms', 'adiabat', 'ifar',
                     'hydro_efficiency_pct', 'imploded_DT_mass_mg',
                     'inflight_KE_kJ', 'fraction_absorbed_pct'
        Per-code HYDRA keys (optional): '{stagnation,ignition,bang}_time_HYDRA_ns',
          'peak_rhoR_total_HYDRA_gcm2', 'peak_rhoR_hs_HYDRA_gcm2',
          'peak_density_at_ignition_HYDRA_gcm3', 'T_ion_hs_at_ignition_HYDRA_keV',
          'hs_radius_ignition_HYDRA_um'.
    laser_energy_MJ : float, optional
        Override laser energy for gain calculation.

    Returns
    -------
    str
        Formatted comparison table.
    """
    if laser_energy_MJ is None:
        laser_energy_MJ = sim_metrics.get('laser_energy_MJ', 0.0)

    sim_gain = sim_metrics['yield_MJ'] / laser_energy_MJ if laser_energy_MJ > 0 else 0.0

    # Column widths -- keep cluster column at original width for readability
    HDR_W   = 30
    SIM_W   = 15
    PUB_W   = 15
    DELTA_W = 10

    lines = []
    lines.append("=" * (HDR_W + SIM_W + 2 * (PUB_W + DELTA_W) + 5))
    lines.append("COMPARISON WITH PUBLISHED DATA")
    lines.append("=" * (HDR_W + SIM_W + 2 * (PUB_W + DELTA_W) + 5))
    lines.append(f"  Laser energy:  Sim = {sim_metrics.get('laser_energy_MJ', 0.0):.3f} MJ"
                 f"   Published = {laser_energy_MJ:.3f} MJ")
    lines.append("")
    hdr = (f"{'Metric':<{HDR_W}} {'Simulation':>{SIM_W}} "
           f"{'Cluster':>{PUB_W}} {'Δ_cluster':>{DELTA_W}} "
           f"{'HYDRA':>{PUB_W}} {'Δ_HYDRA':>{DELTA_W}}")
    lines.append(hdr)
    lines.append("-" * len(hdr))

    def _render_cell(pub_val, pub_unc, fmt):
        """Format a 'Published ± unc' cell."""
        if pub_val is None or pub_val <= 0:
            return "—"
        out = format(pub_val, fmt)
        if pub_unc and pub_unc > 0:
            out += f"±{format(pub_unc, fmt)}"
        return out

    def _render_delta(sim_val, pub_val):
        """Format a 'Δ (%)' cell as a signed percent of pub_val."""
        if (sim_val is None or sim_val <= 0
                or pub_val is None or pub_val <= 0):
            return "—"
        d = 100.0 * (sim_val - pub_val) / pub_val
        if not np.isfinite(d):
            return "—"
        return f"{d:>+.1f}"

    def _emit_row(label, sim_val, cluster_entry, hydra_entry, fmt, label_suffix=""):
        """Emit one fully-formatted row with all 5 numeric cells."""
        cluster_val, cluster_unc = (_to_tuple(cluster_entry)
                                     if cluster_entry is not None else (None, None))
        hydra_val,   hydra_unc   = (_to_tuple(hydra_entry)
                                     if hydra_entry is not None else (None, None))
        # Skip the whole row only if sim_val and BOTH refs are missing
        if (sim_val <= 0
                and (cluster_val is None or cluster_val <= 0)
                and (hydra_val is None or hydra_val <= 0)):
            return False
        sv = format(sim_val, fmt) if sim_val > 0 else "—"
        cv = _render_cell(cluster_val, cluster_unc, fmt)
        dc = _render_delta(sim_val, cluster_val)
        hv = _render_cell(hydra_val,   hydra_unc,   fmt)
        dh = _render_delta(sim_val, hydra_val)
        lbl = f"{label}{label_suffix}"
        lines.append(
            f"{lbl:<{HDR_W}} {sv:>{SIM_W}} "
            f"{cv:>{PUB_W}} {dc:>{DELTA_W}} "
            f"{hv:>{PUB_W}} {dh:>{DELTA_W}}"
        )
        return True

    # ─── Adiabat convention-aware reference-key selection ───
    # Each adiabat row in the implosion_rows table compares against a
    # published-JSON key that matches its OWN convention. Conventions
    # in the published-data JSON (only fill the relevant ones per
    # publication):
    #   'adiabat'                          — Lindl peak v mass-avg
    #                                        (Olson 2021 / classical
    #                                        ICF convention, default).
    #   'adiabat_lindl_cr15'               — Lindl mass-avg at CR=1.5.
    #   'adiabat_rhino_min_cr15'           — RHINO min shell adiabat at
    #                                        CR=1.5 (Thomas Vulcan HDD;
    #                                        Will Trickey postprocessor).
    #   'adiabat_rhino_formula_peak_v'     — proper-Fermi mass-avg at peak v.
    #   'adiabat_rhino_formula_cr15'       — proper-Fermi mass-avg at CR=1.5.
    #   'adiabat_rhino_formula_breakout'   — proper-Fermi mass-avg at breakout.
    # When a convention-specific key is absent in the published data,
    # the corresponding row's cluster_key is set to a sentinel name not
    # present in published_metrics so the comparison helper renders '—'
    # for the published column instead of showing a misleading
    # +860%-type mismatch against the wrong-convention reference.
    def _adi_ref(key):
        """Return key if it's in published_metrics, else a sentinel that
        won't be found (which suppresses the published-side comparison)."""
        if published_metrics and key in published_metrics:
            return key
        return '__adi_ref_missing__'
    # Lindl peak v (standard ICF / Olson convention): try 'adiabat'.
    adi_lindl_pv_key   = _adi_ref('adiabat')
    # Lindl CR=1.5: prefer dedicated key, no fallback to 'adiabat'
    # (those two are different physical quantities).
    adi_lindl_cr15_key = _adi_ref('adiabat_lindl_cr15')
    # RHINO min CR=1.5 (Thomas Vulcan HDD): prefer the convention-
    # specific key; fall back to 'adiabat' (compares to Olson Lindl,
    # which is a different convention -- value shown but Δ should be
    # interpreted carefully).
    adi_rhino_min_key  = ('adiabat_rhino_min_cr15'
                          if (published_metrics and 'adiabat_rhino_min_cr15' in published_metrics)
                          else _adi_ref('adiabat'))
    # RHINO-formula adiabats (proper Fermi): no fallback. Rare in
    # publications; suppress Δ when absent.
    adi_rhino_pv_key   = _adi_ref('adiabat_rhino_formula_peak_v')
    adi_rhino_cr15_key = _adi_ref('adiabat_rhino_formula_cr15')
    adi_rhino_bo_key   = _adi_ref('adiabat_rhino_formula_breakout')

    # ── Implosion metrics section ──
    # Tuple: (label, sim_val, cluster_key, hydra_key, fmt)
    # hydra_key=None means no per-code HYDRA value for this metric (cluster only).
    implosion_rows = [
        ('Stagnation time (ns)',
         sim_metrics.get('stagnation_time_ns', 0.0),
         'stagnation_time_ns', 'stagnation_time_HYDRA_ns', '.3f'),
        ('Peak velocity (km/s)',
         sim_metrics.get('peak_velocity_kms', 0.0),
         'peak_velocity_kms', None, '.1f'),
        ('Peak velocity at CR=1.5 (km/s)',
         sim_metrics.get('peak_velocity_kms_cr15', 0.0),
         'peak_velocity_kms', None, '.1f'),
        ('Implosion velocity RHINO (km/s)',
         sim_metrics.get('implosion_velocity_rhino_kms', 0.0),
         'peak_velocity_kms', None, '.1f'),
        ('Min shell adiabat RHINO',
         sim_metrics.get('adiabat_min_rhino', 0.0),
         adi_rhino_min_key, None, '.2f'),
        ('Adiabat (Lindl peak v)',
         sim_metrics.get('adiabat', 0.0),
         adi_lindl_pv_key, None, '.2f'),
        ('Adiabat at CR=1.5 (Lindl)',
         sim_metrics.get('adiabat_cr15', 0.0),
         adi_lindl_cr15_key, None, '.2f'),
        ('Adiabat at peak v (RHINO formula)',
         sim_metrics.get('adiabat_rhino_formula', 0.0),
         adi_rhino_pv_key, None, '.2f'),
        ('Adiabat at CR=1.5 (RHINO formula)',
         sim_metrics.get('adiabat_cr15_rhino_formula', 0.0),
         adi_rhino_cr15_key, None, '.2f'),
        ('Adiabat at breakout (Lindl)',
         sim_metrics.get('adiabat_at_breakout', 0.0),
         _adi_ref('adiabat_lindl_breakout'), None, '.2f'),
        # --- RHINO diagnostics suite (June 2026) ---
        ('Peak shell velocity time RHINO (ns)',
         sim_metrics.get('t_max_shell_velocity_rhino_ns', 0.0),
         _adi_ref('t_max_shell_velocity_rhino_ns'), None, '.3f'),
        ('Stagnation time RHINO (ns)',
         sim_metrics.get('stag_time_rhino_ns', 0.0),
         _adi_ref('stag_time_rhino_ns'), None, '.3f'),
        ('Assembled mass RHINO (mg)',
         sim_metrics.get('assembled_mass_rhino_mg', 0.0),
         _adi_ref('assembled_mass_rhino_mg'), None, '.3f'),
        ('Burn fraction RHINO',
         sim_metrics.get('burn_fraction_rhino', 0.0),
         _adi_ref('burn_fraction_rhino'), None, '.3f'),
        ('Ablation P at CR=3.5 (Mbar)',
         sim_metrics.get('ablation_pressure_at_cr_3p5_Mbar', 0.0),
         _adi_ref('ablation_pressure_at_cr_3p5_Mbar'), None, '.1f'),
        ('Adiabat at breakout (RHINO formula)',
         sim_metrics.get('adiabat_breakout_rhino_formula', 0.0),
         adi_rhino_bo_key, None, '.2f'),
        ('Fraction absorbed (%)',
         sim_metrics.get('fraction_absorbed_pct', 0.0),
         'fraction_absorbed_pct', None, '.1f'),
        ('In-flight KE (kJ)',
         sim_metrics.get('inflight_KE_kJ', 0.0),
         'inflight_KE_kJ', None, '.1f'),
        ('IFAR (CR=1.5)',
         sim_metrics.get('ifar', 0.0),
         'ifar', None, '.1f'),
        ('Hydro efficiency (%)',
         sim_metrics.get('hydro_efficiency_pct', 0.0),
         'hydro_efficiency_pct', None, '.1f'),
        ('Imploded DT mass (mg)',
         sim_metrics.get('imploded_DT_mass_mg', 0.0),
         'imploded_DT_mass_mg', None, '.2f'),
        ('HS pressure at ignition (Gbar)',
         sim_metrics.get('P_hs_ignition_Gbar', 0.0),
         'P_hs_ignition_Gbar', None, '.1f'),
        ('HS radius at ignition (μm)',
         sim_metrics.get('hs_radius_ignition_um', 0.0),
         'hs_radius_ignition_um', 'hs_radius_ignition_HYDRA_um', '.1f'),
        ('On-axis T_ion at ignition (keV)',
         sim_metrics.get('T_ion_onaxis_ignition_keV', 0.0),
         'T_ion_onaxis_ignition_keV', None, '.2f'),
        ('⟨T_ion⟩_HS at ignition (keV)',
         sim_metrics.get('T_ion_hs_avg_ignition_keV', 0.0),
         None, 'T_ion_hs_at_ignition_HYDRA_keV', '.2f'),
        # --- Laser intensity diagnostics (Task 2 Stage C) ---
        ('I at r_crit peak (W/cm²)',
         sim_metrics.get('I_at_crit_peak_Wcm2', 0.0),
         'I_at_crit_peak_Wcm2', None, '.2e'),
        ('I grid outer peak (W/cm²)',
         sim_metrics.get('I_grid_outer_peak_Wcm2', 0.0),
         'I_grid_outer_peak_Wcm2', None, '.2e'),
        # --- Shock train breakouts (Task 3 Stage 3 multi-shock tracker) ---
        # NaN sentinel from extract_histories -> mapped to -1.0 so the
        # "<= 0 == missing" skip logic in _emit_row kicks in.
        ('Foot shock breakout (ns)',
         (sim_metrics.get('t_foot_shock_breakout_ns', float('nan'))
          if np.isfinite(sim_metrics.get('t_foot_shock_breakout_ns', float('nan')))
          else -1.0),
         't_foot_shock_breakout_ns', None, '.2f'),
        ('Ramp shock breakout (ns)',
         (sim_metrics.get('t_ramp_shock_breakout_ns', float('nan'))
          if np.isfinite(sim_metrics.get('t_ramp_shock_breakout_ns', float('nan')))
          else -1.0),
         't_ramp_shock_breakout_ns', None, '.2f'),
        ('Peak shock breakout (ns)',
         (sim_metrics.get('t_peak_shock_breakout_ns', float('nan'))
          if np.isfinite(sim_metrics.get('t_peak_shock_breakout_ns', float('nan')))
          else -1.0),
         't_peak_shock_breakout_ns', None, '.2f'),
        # --- Will list extensions (June 4 2026) ---
        ('Outer convergence ratio',
         sim_metrics.get('cr_outer', 0.0),
         'cr_outer', None, '.2f'),
        ('Laser overlapped intensity (W/cm²)',
         sim_metrics.get('laser_overlapped_intensity_Wcm2', 0.0),
         'laser_overlapped_intensity_Wcm2', None, '.2e'),
        ('Shell mass at stagnation (mg)',
         sim_metrics.get('shell_mass_at_stagnation_mg', 0.0),
         'shell_mass_at_stagnation_mg', None, '.3f'),
        ('Hot-spot mass at stagnation (mg)',
         sim_metrics.get('hot_spot_mass_at_stagnation_mg', 0.0),
         'hot_spot_mass_at_stagnation_mg', None, '.4f'),
        ('Total ρR at stagnation (g/cm²)',
         sim_metrics.get('stag_time_areal_density', 0.0),
         'stag_time_areal_density', None, '.3f'),
        ('Cold-fuel ρR at stagnation (g/cm²)',
         sim_metrics.get('stag_time_fuel_areal_density', 0.0),
         'stag_time_fuel_areal_density', None, '.3f'),
        ('Neutron-averaged T_e (keV)',
         sim_metrics.get('neutron_ave_electron_temperature', 0.0),
         'neutron_ave_electron_temperature', None, '.2f'),
        # --- Will shell convention (June 2026, Phys. Conv. §18b) ---
        ('Shell mass at stag (Will, mg)',
         sim_metrics.get('shell_mass_will_at_stagnation_mg', 0.0),
         'shell_mass_will_at_stagnation_mg', None, '.3f'),
        ('Mass-avg shell adiabat at CR=1.5 (Will)',
         sim_metrics.get('adiabat_mass_avg_will_cr15', 0.0),
         'adiabat_mass_avg_will_cr15', None, '.2f'),
        ('Sound speed in shell at CR=1.5 (km/s, Will)',
         sim_metrics.get('sound_speed_shell_will_cr15_kms', 0.0),
         'sound_speed_shell_will_cr15_kms', None, '.1f'),
        # --- Will ignition-product timing (June 2026, Phys. Conv. §18d/e) ---
        ('t_peak(rhoR x Ti) (ns)',
         sim_metrics.get('t_peak_rhoR_Ti_ns', 0.0),
         't_peak_rhoR_Ti_ns', None, '.3f'),
        ('rhoR_hs at peak(rhoR x Ti) (g/cm²)',
         sim_metrics.get('rhoR_hs_at_peak_rhoR_Ti_gcm2', 0.0),
         'rhoR_hs_at_peak_rhoR_Ti_gcm2', None, '.3f'),
        ('T_i HS vol-avg at peak(rhoR x Ti) (keV)',
         sim_metrics.get('T_i_volume_avg_at_peak_rhoR_Ti_keV', 0.0),
         'T_i_volume_avg_at_peak_rhoR_Ti_keV', None, '.2f'),
        ('T_i HS mass-avg at peak(rhoR x Ti) (keV)',
         sim_metrics.get('T_i_mass_avg_at_peak_rhoR_Ti_keV', 0.0),
         'T_i_mass_avg_at_peak_rhoR_Ti_keV', None, '.2f'),
        ('t_peak(P x R) (ns)',
         sim_metrics.get('t_peak_Phs_Rhs_ns', 0.0),
         't_peak_Phs_Rhs_ns', None, '.3f'),
        ('R_hs at peak(P x R) (μm)',
         sim_metrics.get('R_hs_at_peak_Phs_Rhs_um', 0.0),
         'R_hs_at_peak_Phs_Rhs_um', None, '.1f'),
        ('P_hs vol-avg at peak(P x R) (Gbar)',
         sim_metrics.get('P_hs_volume_avg_at_peak_Phs_Rhs_Gbar', 0.0),
         'P_hs_volume_avg_at_peak_Phs_Rhs_Gbar', None, '.1f'),
        ('P_hs mass-avg at peak(P x R) (Gbar)',
         sim_metrics.get('P_hs_mass_avg_at_peak_Phs_Rhs_Gbar', 0.0),
         'P_hs_mass_avg_at_peak_Phs_Rhs_Gbar', None, '.1f'),
    ]

    has_implosion = False
    for label, sim_val, cluster_key, hydra_key, fmt in implosion_rows:
        cluster_entry = published_metrics.get(cluster_key) if cluster_key else None
        hydra_entry   = published_metrics.get(hydra_key) if hydra_key else None
        emitted = _emit_row(label, sim_val, cluster_entry, hydra_entry, fmt)
        if emitted:
            has_implosion = True

    if has_implosion:
        lines.append("-" * len(hdr))

    # ── Burn-averaged metrics section ──
    # Cluster-only rows: T_hs, P_hs, rhoR_cf, CR_max, yield, gain (no per-code
    # HYDRA values for these in Olson 2021).
    burn_rows = [
        ('⟨T_hs⟩ (keV)',
         sim_metrics['T_burn_avg'],     'T_hs',    None, '.1f'),
        ('⟨P_hs⟩ (Gbar)',
         sim_metrics['P_burn_avg'],     'P_hs',    None, '.0f'),
        ('⟨ρR_cf⟩ (g/cm²)',
         sim_metrics['rhoR_burn_avg'],  'rhoR_cf', None, '.2f'),
        ('CR_max',
         sim_metrics['CR_max'],         'CR_max',  None, '.1f'),
        ('Peak total ρR (g/cm²)',
         sim_metrics.get('peak_total_rhoR', 0.0),
         'peak_total_rhoR_gcm2',  'peak_rhoR_total_HYDRA_gcm2', '.2f'),
        ('Peak HS ρR T>4.5 (g/cm²)',
         sim_metrics.get('peak_hs_rhoR_T_mask', 0.0),
         'peak_hs_rhoR_T_mask_gcm2', 'peak_rhoR_hs_HYDRA_gcm2', '.2f'),
        ('Yield (MJ)',
         sim_metrics['yield_MJ'],       'yield',   None, '.1f'),
        ('Fusion Gain',
         sim_gain,                      'gain',    None, '.1f'),
    ]

    for label, sim_val, cluster_key, hydra_key, fmt in burn_rows:
        cluster_entry = published_metrics.get(cluster_key) if cluster_key else None
        hydra_entry   = published_metrics.get(hydra_key) if hydra_key else None
        _emit_row(label, sim_val, cluster_entry, hydra_entry, fmt)

    # ── Peak ρ at ignition (HYDRA-specific, no cluster value) ──
    # Headline no-burn hydro metric. For no-burn runs the value is at
    # stagnation (no ignition crossing); label the row accordingly.
    _rho_ign_val = sim_metrics.get('peak_density_at_ignition', 0.0)
    _is_stag_fb  = bool(sim_metrics.get('peak_density_at_ignition_is_stagnation', False))
    _rho_label   = 'Peak ρ at ignition (g/cm³)'
    if _is_stag_fb:
        _rho_label = 'Peak ρ at stagnation (g/cm³)'   # no-burn fallback
    _emit_row(_rho_label, _rho_ign_val,
              None,
              published_metrics.get('peak_density_at_ignition_HYDRA_gcm3'),
              '.0f')

    # ── Over-amplified / oversized hot-spot diagnostic (May 23 2026) ──
    # Two prior framings of this residual were wrong: the "target-mass
    # mismatch" of May 2026 (Helios's cold-fuel tag was misread; the DT
    # ice is actually inside the hot spot), and the "early-shock heating"
    # of May 22 2026 (the high on-axis T at ignition is alpha self-
    # heating, not shock convergence — confirmed by running fab007 with
    # burn OFF: identical hydro, ⟨T_hs⟩ drops from 23 to 5 keV).
    #
    # The corrected picture: when the alpha bootstrap catches, it
    # amplifies T_hs from the ~5 keV hydrodynamic floor up to ignition-
    # range values AND spreads the T>4.5 keV mask outward into less-
    # dense ice. The result is an oversized, over-heated hot spot at
    # ignition with a confinement deficit (Helios HS ρR T>4.5 = 0.35 vs
    # LILAC's 0.85). This is downstream of the calibrated drive — not a
    # drive-phase knob. The warning surfaces the symptom; the cause is
    # alpha-amplification-driven hot-spot bloat at the ignition cliff.
    _T_axis     = float(sim_metrics.get('T_ion_onaxis_ignition_keV', 0.0))
    _R_hs_ign   = float(sim_metrics.get('hs_radius_ignition_um', 0.0))
    _pub_T_e    = published_metrics.get('T_ion_onaxis_ignition_keV', None)
    _pub_R_e    = published_metrics.get('hs_radius_ignition_um', None)
    _pub_foot_e = published_metrics.get('t_foot_shock_breakout_ns', None)
    _pub_ramp_e = published_metrics.get('t_ramp_shock_breakout_ns', None)
    _sim_foot   = float(sim_metrics.get('t_foot_shock_breakout_ns', float('nan')))
    _sim_ramp   = float(sim_metrics.get('t_ramp_shock_breakout_ns', float('nan')))

    _T_hot = False
    _R_big = False
    if _pub_T_e is not None and _T_axis > 0:
        _pub_T, _ = _to_tuple(_pub_T_e)
        if _pub_T > 0 and _T_axis > 1.5 * _pub_T:
            _T_hot = True
    if _pub_R_e is not None and _R_hs_ign > 0:
        _pub_R, _ = _to_tuple(_pub_R_e)
        if _pub_R > 0 and _R_hs_ign > 1.2 * _pub_R:
            _R_big = True

    if _T_hot and _R_big:
        lines.append("")
        lines.append("-" * 80)
        lines.append(
            "WARNING: HOT SPOT OVER-HEATED AND OVERSIZED AT IGNITION.")
        if _pub_T_e is not None:
            _pub_T, _ = _to_tuple(_pub_T_e)
            lines.append(
                f"  On-axis T_ion = {_T_axis:.1f} keV vs ref {_pub_T:.1f} keV "
                f"({_T_axis/_pub_T:.1f}x).")
        if _pub_R_e is not None:
            _pub_R, _ = _to_tuple(_pub_R_e)
            lines.append(
                f"  Hot-spot radius = {_R_hs_ign:.1f} um vs ref {_pub_R:.1f} um "
                f"({_R_hs_ign/_pub_R:.2f}x).")
        # Shock-timing context (if available)
        if _pub_foot_e is not None and np.isfinite(_sim_foot):
            _pub_foot, _ = _to_tuple(_pub_foot_e)
            if _pub_foot > 0:
                lines.append(
                    f"  Foot shock arrival: sim {_sim_foot:.2f} ns vs ref "
                    f"{_pub_foot:.2f} ns ({_sim_foot - _pub_foot:+.2f} ns).")
        if _pub_ramp_e is not None and np.isfinite(_sim_ramp):
            _pub_ramp, _ = _to_tuple(_pub_ramp_e)
            if _pub_ramp > 0:
                lines.append(
                    f"  Ramp shock arrival: sim {_sim_ramp:.2f} ns vs ref "
                    f"{_pub_ramp:.2f} ns ({_sim_ramp - _pub_ramp:+.2f} ns).")
        lines.append(
            "  Likely cause: alpha self-heating amplifying the central hot")
        lines.append(
            "  region beyond LILAC's confined state. In 1D Helios the alpha-")
        lines.append(
            "  amplified T>4.5 keV mask spreads outward into less-dense ice,")
        lines.append(
            "  producing an oversized & under-confined hot spot relative to")
        lines.append(
            "  the reference codes. This is downstream of the drive — not a")
        lines.append(
            "  drive-phase calibration knob. Run a burn-OFF equivalent to see")
        lines.append(
            "  the hydrodynamic floor (~5 keV for fab007-class PDD configs).")

    lines.append("=" * 80)
    return "\n".join(lines)


def _to_tuple(val):
    """Convert [value, unc] list or (value, unc) tuple to (float, float)."""
    if isinstance(val, (list, tuple)) and len(val) >= 2:
        return (float(val[0]), float(val[1]))
    elif isinstance(val, (int, float)):
        return (float(val), 0.0)
    return (0.0, 0.0)
