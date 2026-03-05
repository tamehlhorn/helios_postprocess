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
    peak_vimp_kms = abs(getattr(data, 'peak_implosion_velocity', 0.0))
    adiabat       = getattr(data, 'adiabat_mass_averaged_ice', 0.0)

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
        E_delivered_J = np.trapz(lpd, x=time_ns * 1e-9)
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
        # Imploded mass = all mass inside ablation front or fuel/ablator boundary
        abl_front = getattr(data, 'ablation_front_radius', None)
        if abl_front is not None:
            r_abl = abl_front[stag_idx]
            # Sum mass of zones with outer boundary <= ablation front
            for z in range(zmass.shape[1]):
                if zbnd[stag_idx, z + 1] <= r_abl * 1.01:  # small tolerance
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
        'laser_energy_MJ':    laser_energy_MJ,
        'target_gain':        target_gain,
        # Implosion metrics
        'peak_velocity_kms':     peak_vimp_kms,
        'adiabat':               adiabat,
        'fraction_absorbed_pct': fraction_absorbed,
        'inflight_KE_kJ':       inflight_KE_kJ,
        'hydro_efficiency_pct':  hydro_efficiency_pct,
        'imploded_DT_mass_mg':   imploded_DT_mass_mg,
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
    r0 = histories['initial_radius_um']
    valid = radius_um > 0
    if r0 > 0 and np.any(valid):
        CR_max = np.max(r0 / radius_um[valid])
    else:
        CR_max = 0.0

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
        'adiabat':               histories.get('adiabat', 0.0),
        'fraction_absorbed_pct': histories.get('fraction_absorbed_pct', 0.0),
        'inflight_KE_kJ':       histories.get('inflight_KE_kJ', 0.0),
        'hydro_efficiency_pct':  histories.get('hydro_efficiency_pct', 0.0),
        'imploded_DT_mass_mg':   histories.get('imploded_DT_mass_mg', 0.0),
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

    Parameters
    ----------
    sim_metrics : dict
        From calculate_burn_averaged_metrics().
    published_metrics : dict
        Published values: {key: (value, uncertainty)}.
        Supported keys:
          Burn-averaged: 'T_hs', 'P_hs', 'rhoR_cf', 'CR_max', 'yield', 'gain'
          Implosion: 'peak_velocity_kms', 'adiabat', 'ifar',
                     'hydro_efficiency_pct', 'imploded_DT_mass_mg',
                     'inflight_KE_kJ', 'fraction_absorbed_pct'
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

    lines = []
    lines.append("=" * 80)
    lines.append("COMPARISON WITH PUBLISHED DATA")
    lines.append("=" * 80)
    lines.append(f"  Laser energy:  Sim = {sim_metrics.get('laser_energy_MJ', 0.0):.3f} MJ"
                 f"   Published = {laser_energy_MJ:.3f} MJ")
    lines.append("")
    lines.append(f"{'Metric':<30} {'Simulation':>15} {'Published':>15} {'Δ (%)':>10}")
    lines.append("-" * 80)

    # ── Implosion metrics section ──
    implosion_rows = [
        ('Peak velocity (km/s)',      sim_metrics.get('peak_velocity_kms', 0.0),
         'peak_velocity_kms',         '.1f'),
        ('Adiabat',                   sim_metrics.get('adiabat', 0.0),
         'adiabat',                   '.2f'),
        ('Fraction absorbed (%)',     sim_metrics.get('fraction_absorbed_pct', 0.0),
         'fraction_absorbed_pct',     '.1f'),
        ('In-flight KE (kJ)',        sim_metrics.get('inflight_KE_kJ', 0.0),
         'inflight_KE_kJ',           '.1f'),
        ('Hydro efficiency (%)',     sim_metrics.get('hydro_efficiency_pct', 0.0),
         'hydro_efficiency_pct',     '.1f'),
        ('Imploded DT mass (mg)',    sim_metrics.get('imploded_DT_mass_mg', 0.0),
         'imploded_DT_mass_mg',      '.2f'),
    ]

    has_implosion = False
    for label, sim_val, pub_key, fmt in implosion_rows:
        pub_entry = published_metrics.get(pub_key, None)
        if pub_entry is None:
            continue
        pub_val, pub_unc = _to_tuple(pub_entry)
        if pub_val <= 0 and sim_val <= 0:
            continue
        has_implosion = True
        sv = format(sim_val, fmt) if sim_val > 0 else "—"
        if pub_val > 0:
            delta = 100 * (sim_val - pub_val) / pub_val if sim_val > 0 else float('nan')
            pv = format(pub_val, fmt)
            if pub_unc > 0:
                pv += f"±{format(pub_unc, fmt)}"
            delta_str = f"{delta:>9.1f}" if np.isfinite(delta) else "      —"
        else:
            pv = "—"
            delta_str = "      —"
        lines.append(f"{label:<30} {sv:>15} {pv:>15} {delta_str}")

    if has_implosion:
        lines.append("-" * 80)

    # ── Burn-averaged metrics section ──
    burn_rows = [
        ('⟨T_hs⟩ (keV)',    sim_metrics['T_burn_avg'],    'T_hs',    '.1f'),
        ('⟨P_hs⟩ (Gbar)',   sim_metrics['P_burn_avg'],    'P_hs',    '.0f'),
        ('⟨ρR_cf⟩ (g/cm²)', sim_metrics['rhoR_burn_avg'], 'rhoR_cf', '.2f'),
        ('CR_max',           sim_metrics['CR_max'],         'CR_max',  '.1f'),
        ('Yield (MJ)',       sim_metrics['yield_MJ'],       'yield',   '.1f'),
        ('Fusion Gain',      sim_gain,                      'gain',    '.1f'),
    ]

    for label, sim_val, pub_key, fmt in burn_rows:
        pub_entry = published_metrics.get(pub_key, None)
        if pub_entry is None:
            continue
        pub_val, pub_unc = _to_tuple(pub_entry)
        if pub_val <= 0:
            continue
        delta = 100 * (sim_val - pub_val) / pub_val
        sv = format(sim_val, fmt)
        pv = format(pub_val, fmt)
        if pub_unc > 0:
            pv += f"±{format(pub_unc, fmt)}"
        lines.append(f"{label:<30} {sv:>15} {pv:>15} {delta:>9.1f}")

    lines.append("=" * 80)
    return "\n".join(lines)


def _to_tuple(val):
    """Convert [value, unc] list or (value, unc) tuple to (float, float)."""
    if isinstance(val, (list, tuple)) and len(val) >= 2:
        return (float(val[0]), float(val[1]))
    elif isinstance(val, (int, float)):
        return (float(val), 0.0)
    return (0.0, 0.0)
