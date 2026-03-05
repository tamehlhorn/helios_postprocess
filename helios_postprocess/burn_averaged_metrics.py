"""
Burn-Averaged Metrics Module
=============================

Calculates burn-averaged (temporal) quantities for comparison with published
ICF performance metrics.

This module provides TIME-HISTORY burn-averaging using the 2D arrays stored
on ICFRunData (after build_run_data + ICFAnalyzer have run), complementing
the SPATIAL neutron-averaging in icf_analysis.py.

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
    >>> print(compare_with_published(metrics, published_data, laser_energy_MJ=4.0))

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

    This replaces the old extract_hot_spot_histories() which looped over
    timesteps calling the legacy hot_spot/areal_density modules.

    Parameters
    ----------
    data : ICFRunData
        Populated dataclass — must have had build_run_data() and
        ICFAnalyzer (at least analyze_stagnation_phase) run on it.
        Required attributes:
            time, ion_temperature, mass_density, ion_pressure, rad_pressure,
            zone_boundaries, zone_mass, region_interfaces_indices
        Optional (for areal density):
            areal_density_vs_time

    Returns
    -------
    dict
        Dictionary with keys expected by calculate_burn_averaged_metrics():
        - 'time_ns'            : (n_times,) time in nanoseconds
        - 'temperature_keV'    : (n_times,) mass-avg hot-spot ion T in keV
        - 'pressure_Gbar'      : (n_times,) mass-avg hot-spot total P in Gbar
        - 'density_gcc'        : (n_times,) mass-avg hot-spot density in g/cm³
        - 'areal_density_gcm2' : (n_times,) cold-fuel ρR in g/cm²
        - 'radius_um'          : (n_times,) hot-spot outer radius in μm
        - 'volume_cm3'         : (n_times,) hot-spot volume in cm³ (sphere)
        - 'mass_mg'            : (n_times,) hot-spot mass in mg
        - 'initial_radius_um'  : scalar, radius at first timestep

    Raises
    ------
    ValueError
        If required arrays are missing from ICFRunData.
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
    T_ion = data.ion_temperature                 # (n_times, n_zones) eV
    rho   = data.mass_density                    # (n_times, n_zones) g/cm³
    P_ion = data.ion_pressure                    # (n_times, n_zones) J/cm³
    P_rad = data.rad_pressure                    # (n_times, n_zones) J/cm³
    P_tot = P_ion + P_rad                        # total pressure
    zmass = data.zone_mass                       # (n_times, n_zones) g
    zbnd  = data.zone_boundaries                 # (n_times, n_nodes) cm

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

    return {
        'time_ns':            time_ns,
        'temperature_keV':    temperatures,
        'pressure_Gbar':      pressures,
        'density_gcc':        densities,
        'areal_density_gcm2': rhoR_cf,
        'radius_um':          radii_um,
        'volume_cm3':         volumes,
        'mass_mg':            masses_mg,
        'initial_radius_um':  radii_um[0],
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

    If called with an ICFRunData object (has 'ion_temperature' attribute),
    extracts histories directly.  If called with a HeliosRun object, raises
    an informative error directing the caller to use the pipeline.

    Parameters
    ----------
    run_or_data : ICFRunData or HeliosRun
        Populated ICFRunData (preferred) or HeliosRun (raises error).
    T_threshold : float, optional
        Ignored (kept for API compatibility).
    time_indices : array-like, optional
        Ignored (kept for API compatibility).

    Returns
    -------
    dict
        Histories dictionary for calculate_burn_averaged_metrics().
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

    Burn-averaging: <Q> = ∫ Q(t) · Ṙ(t) dt / ∫ Ṙ(t) dt
    where Ṙ is the volumetric fusion reaction rate ∝ ρ² · <σv>(T).

    Parameters
    ----------
    histories : dict
        Dictionary from extract_histories_from_run_data().
    ion_fraction : float, optional
        Ion fraction for burn rate calculation (default 0.5 for 50-50 DT).

    Returns
    -------
    dict
        - 'T_burn_avg'    : Burn-averaged temperature (keV)
        - 'P_burn_avg'    : Burn-averaged pressure (Gbar)
        - 'rho_burn_avg'  : Burn-averaged density (g/cm³)
        - 'rhoR_burn_avg' : Burn-averaged cold-fuel areal density (g/cm²)
        - 'CR_max'        : Maximum convergence ratio
        - 'yield_MJ'      : Total fusion yield (MJ)
        - 'T_peak', 'P_peak', 'rho_peak', 'rhoR_peak' : Peak values
        - 'burn_fraction'  : Normalized burn rate profile
        - 'burn_rate'      : Raw burn rate array
    """
    time_ns    = histories['time_ns']
    temp_keV   = histories['temperature_keV']
    pres_Gbar  = histories['pressure_Gbar']
    dens_gcc   = histories['density_gcc']
    rhoR_gcm2  = histories['areal_density_gcm2']
    radius_um  = histories['radius_um']
    volume_cm3 = histories['volume_cm3']

    time_s = time_ns * 1e-9

    # Burn rate at each timestep
    burn_rate = calculate_burn_rate_from_sim(temp_keV, dens_gcc, ion_fraction)

    # Burn-averaged quantities
    def _burn_avg(quantity):
        num = simpson(quantity * burn_rate, x=time_s)
        den = simpson(burn_rate, x=time_s)
        return num / den if den > 0 else 0.0

    T_burn_avg    = _burn_avg(temp_keV)
    P_burn_avg    = _burn_avg(pres_Gbar)
    rho_burn_avg  = _burn_avg(dens_gcc)
    rhoR_burn_avg = _burn_avg(rhoR_gcm2)

    # Convergence ratio
    r0 = histories['initial_radius_um']
    valid = radius_um > 0
    if r0 > 0 and np.any(valid):
        CR_max = np.max(r0 / radius_um[valid])
    else:
        CR_max = 0.0

    # Total yield from volumetric burn rate
    E_fusion = 17.6 * 1.60218e-13   # MeV → J
    reaction_rate_total = burn_rate * volume_cm3
    total_reactions = simpson(reaction_rate_total, x=time_s)
    yield_J  = total_reactions * E_fusion
    yield_MJ = yield_J * 1e-6

    # Normalized burn fraction (for plotting)
    br_sum = np.sum(burn_rate)
    burn_fraction = burn_rate / br_sum if br_sum > 0 else burn_rate

    # Min radius (skip zeros)
    valid_r = radius_um[radius_um > 0]
    min_radius = np.min(valid_r) if len(valid_r) > 0 else 0.0

    return {
        'T_burn_avg':       T_burn_avg,
        'P_burn_avg':       P_burn_avg,
        'rho_burn_avg':     rho_burn_avg,
        'rhoR_burn_avg':    rhoR_burn_avg,
        'T_peak':           np.max(temp_keV),
        'P_peak':           np.max(pres_Gbar),
        'rho_peak':         np.max(dens_gcc),
        'rhoR_peak':        np.max(rhoR_gcm2),
        'CR_max':           CR_max,
        'min_radius_um':    min_radius,
        'yield_MJ':         yield_MJ,
        'yield_J':          yield_J,
        'total_reactions':  total_reactions,
        'burn_fraction':    burn_fraction,
        'burn_rate':        burn_rate,
    }


# ── Comparison with published data ───────────────────────────────────────────

def compare_with_published(sim_metrics: Dict,
                           published_metrics: Dict,
                           laser_energy_MJ: float = 4.0) -> str:
    """
    Generate comparison table between simulation and published results.

    Parameters
    ----------
    sim_metrics : dict
        From calculate_burn_averaged_metrics().
    published_metrics : dict
        Published values: {key: (value, uncertainty)}.
        Supported keys: 'T_hs', 'P_hs', 'rhoR_cf', 'CR_max', 'yield', 'gain'
    laser_energy_MJ : float
        Laser energy in MJ for gain calculation.

    Returns
    -------
    str
        Formatted comparison table.
    """
    sim_gain = sim_metrics['yield_MJ'] / laser_energy_MJ if laser_energy_MJ > 0 else 0.0

    lines = []
    lines.append("=" * 80)
    lines.append("COMPARISON WITH PUBLISHED DATA")
    lines.append("=" * 80)
    lines.append(f"{'Metric':<30} {'Simulation':<20} {'Published':<20} {'Δ (%)':<10}")
    lines.append("-" * 80)

    rows = [
        ('⟨T_hs⟩ (keV)',    sim_metrics['T_burn_avg'],    'T_hs',    '.1f'),
        ('⟨P_hs⟩ (Gbar)',   sim_metrics['P_burn_avg'],    'P_hs',    '.0f'),
        ('⟨ρR_cf⟩ (g/cm²)', sim_metrics['rhoR_burn_avg'], 'rhoR_cf', '.2f'),
        ('CR_max',           sim_metrics['CR_max'],         'CR_max',  '.1f'),
        ('Yield (MJ)',       sim_metrics['yield_MJ'],       'yield',   '.1f'),
        ('Fusion Gain',      sim_gain,                      'gain',    '.1f'),
    ]

    for label, sim_val, pub_key, fmt in rows:
        pub_val, pub_unc = published_metrics.get(pub_key, (0.0, 0.0))
        if pub_val > 0:
            delta = 100 * (sim_val - pub_val) / pub_val
            sv = format(sim_val, fmt)
            pv = f"{format(pub_val, fmt)}±{format(pub_unc, fmt)}"
            lines.append(f"{label:<30} {sv:>19} {pv:>19} {delta:>9.1f}")

    lines.append("=" * 80)
    return "\n".join(lines)
