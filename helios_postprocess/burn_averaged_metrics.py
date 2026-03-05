"""
Burn-Averaged Metrics Module
=============================

Calculates burn-averaged (temporal) quantities for comparison with published
ICF performance metrics.

This module provides TIME-HISTORY burn-averaging using data already computed
by the ICFAnalyzer pipeline, complementing the SPATIAL neutron-averaging
in icf_analysis.py.

Workflow
-------
    >>> from helios_postprocess import HeliosRun
    >>> from helios_postprocess.data_builder import build_run_data
    >>> from helios_postprocess.icf_analysis import ICFAnalyzer
    >>> from helios_postprocess.burn_averaged_metrics import (
    ...     extract_histories_from_run_data,
    ...     calculate_burn_averaged_metrics,
    ...     compare_with_published
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
    >>> print(compare_with_published(metrics, published_data))

Key Functions
-------------
extract_histories_from_run_data : Build histories dict from ICFRunData
calculate_burn_averaged_metrics : Calculate all burn-averaged quantities
compare_with_published : Generate comparison tables

Author: Prof T
Date: November 2025 (original), March 2026 (refactored for ICFRunData pipeline)
"""

import numpy as np
from scipy.integrate import simpson
from typing import Dict, Optional
import warnings

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

    Uses the sigma_v_DT reactivity function from core module.

    Parameters
    ----------
    temperature_keV : np.ndarray
        Ion temperature in keV at each time step
    density_gcc : np.ndarray
        Mass density in g/cm³ at each time step
    ion_fraction : float, optional
        Fraction of each ion species (default 0.5 for 50-50 DT)

    Returns
    -------
    np.ndarray
        Fusion reaction rate in reactions/(cm³·s)
    """
    amu_to_g = 1.66054e-24  # g
    m_DT = 2.5 * amu_to_g   # Average mass for 50-50 DT

    n_total = density_gcc / m_DT
    n_D = n_total * ion_fraction
    n_T = n_total * ion_fraction

    sv = sigma_v_DT(temperature_keV)

    return n_D * n_T * sv


# ── History extraction from ICFRunData ───────────────────────────────────────

def extract_histories_from_run_data(data) -> Dict:
    """
    Build histories dict from an ICFRunData object (after ICFAnalyzer has run).

    This replaces the old extract_hot_spot_histories() which looped over
    timesteps and called hot_spot/areal_density modules directly.  All the
    needed arrays are already present on the ICFRunData dataclass after
    ICFAnalyzer.analyze_stagnation_phase() and analyze_burn_phase().

    Parameters
    ----------
    data : ICFRunData
        Populated dataclass — must have had ICFAnalyzer run on it so that
        attributes like hot_spot_temperature, hot_spot_pressure, etc. exist.

    Returns
    -------
    dict
        Dictionary with keys expected by calculate_burn_averaged_metrics():
        - 'time_ns'
        - 'temperature_keV'
        - 'pressure_Gbar'
        - 'density_gcc'
        - 'areal_density_gcm2'
        - 'radius_um'
        - 'volume_cm3'
        - 'mass_mg'
        - 'initial_radius_um'
    """
    # Time — ICFRunData stores time in ns
    time_ns = np.asarray(data.time)

    # ── Hot-spot temperature (keV) ──
    # ICFAnalyzer stores neutron-averaged T in keV on data as a scalar;
    # for a time history we need the per-timestep ion temperature in the
    # hot spot.  Several possible attribute names depending on what
    # ICFAnalyzer has computed.
    temp_keV = _get_history(data, [
        'hot_spot_ion_temperature',      # per-timestep array (keV)
        'neutron_avg_temperature',       # per-timestep array (keV)
    ], fallback_scalar='neutron_avg_Ti_keV', time_ns=time_ns)

    # ── Hot-spot pressure (Gbar) ──
    pres_Gbar = _get_history(data, [
        'hot_spot_pressure',             # per-timestep array (Gbar)
        'neutron_avg_pressure',
    ], fallback_scalar='neutron_avg_P_Gbar', time_ns=time_ns)

    # ── Hot-spot density (g/cm³) ──
    dens_gcc = _get_history(data, [
        'hot_spot_density',
        'neutron_avg_density',
    ], fallback_scalar='neutron_avg_rho', time_ns=time_ns)

    # ── Cold fuel areal density (g/cm²) ──
    rhoR = _get_history(data, [
        'fuel_areal_density',
        'cold_fuel_rhoR',
        'areal_density',
    ], fallback_scalar='neutron_avg_rhoR_fuel', time_ns=time_ns)

    # ── Hot-spot radius (μm) ──
    radius_um = _get_history(data, [
        'hot_spot_radius',               # per-timestep (already μm or cm?)
    ], fallback_scalar=None, time_ns=time_ns)
    # If it looks like cm (values < 1), convert to μm
    if radius_um is not None and np.nanmax(radius_um) < 1.0:
        radius_um = radius_um * 1e4

    # ── Hot-spot volume (cm³) ──
    volume_cm3 = _get_history(data, [
        'hot_spot_volume',
    ], fallback_scalar=None, time_ns=time_ns)
    # If we have radius but no volume, compute assuming sphere
    if (volume_cm3 is None or np.all(volume_cm3 == 0)) and radius_um is not None:
        r_cm = radius_um * 1e-4
        volume_cm3 = (4.0 / 3.0) * np.pi * r_cm ** 3

    # ── Hot-spot mass (mg) ──
    mass_mg = _get_history(data, [
        'hot_spot_mass',
    ], fallback_scalar=None, time_ns=time_ns)

    # ── Initial radius ──
    initial_radius_um = radius_um[0] if radius_um is not None else 0.0

    # Package into dict
    n = len(time_ns)
    return {
        'time_ns':            time_ns,
        'temperature_keV':    _ensure_array(temp_keV, n),
        'pressure_Gbar':      _ensure_array(pres_Gbar, n),
        'density_gcc':        _ensure_array(dens_gcc, n),
        'areal_density_gcm2': _ensure_array(rhoR, n),
        'radius_um':          _ensure_array(radius_um, n),
        'volume_cm3':         _ensure_array(volume_cm3, n),
        'mass_mg':            _ensure_array(mass_mg, n),
        'initial_radius_um':  initial_radius_um,
    }


def _get_history(data, attr_names, fallback_scalar=None, time_ns=None):
    """Try multiple attribute names; return the first that exists and is an array."""
    for name in attr_names:
        val = getattr(data, name, None)
        if val is not None:
            arr = np.asarray(val)
            if arr.ndim >= 1 and len(arr) > 1:
                return arr
    # Try scalar fallback — expand to constant array
    if fallback_scalar is not None and time_ns is not None:
        val = getattr(data, fallback_scalar, None)
        if val is not None and np.isfinite(val):
            warnings.warn(f"Using scalar {fallback_scalar}={val:.4g} as constant history")
            return np.full(len(time_ns), float(val))
    return None


def _ensure_array(arr, n):
    """Return arr if valid, else zeros."""
    if arr is not None and len(arr) == n:
        return arr
    return np.zeros(n)


# ── Legacy compatibility wrapper ─────────────────────────────────────────────

def extract_hot_spot_histories(run_or_data, T_threshold=1000.0,
                               time_indices=None) -> Dict:
    """
    Compatibility wrapper — routes to extract_histories_from_run_data().

    If called with an ICFRunData object (has attribute 'time'), extracts
    histories directly.  If called with a HeliosRun object, raises an
    informative error directing the caller to use the pipeline.

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
    # Check if this is an ICFRunData (has 'time' and 'density' arrays)
    if hasattr(run_or_data, 'time') and hasattr(run_or_data, 'density'):
        return extract_histories_from_run_data(run_or_data)

    # Otherwise assume it's a HeliosRun — guide the user to the new workflow
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

    Parameters
    ----------
    histories : dict
        Dictionary from extract_histories_from_run_data() or
        extract_hot_spot_histories().
    ion_fraction : float, optional
        Ion fraction for burn rate calculation (default 0.5)

    Returns
    -------
    dict
        Dictionary containing:
        - 'T_burn_avg' : Burn-averaged temperature (keV)
        - 'P_burn_avg' : Burn-averaged pressure (Gbar)
        - 'rho_burn_avg' : Burn-averaged density (g/cm³)
        - 'rhoR_burn_avg' : Burn-averaged areal density (g/cm²)
        - 'CR_max' : Maximum convergence ratio
        - 'yield_MJ' : Total fusion yield (MJ)
        - 'T_peak', 'P_peak', etc. : Peak values
        - 'burn_fraction' : Normalized burn rate profile
        - 'burn_rate' : Raw burn rate array
    """
    time_ns   = histories['time_ns']
    temp_keV  = histories['temperature_keV']
    pres_Gbar = histories['pressure_Gbar']
    dens_gcc  = histories['density_gcc']
    rhoR_gcm2 = histories['areal_density_gcm2']
    radius_um = histories['radius_um']
    volume_cm3 = histories['volume_cm3']

    time_s = time_ns * 1e-9

    # Burn rate at each time step
    burn_rate = calculate_burn_rate_from_sim(temp_keV, dens_gcc, ion_fraction)

    # Burn-averaged quantities: <Q> = ∫ Q·R dt / ∫ R dt
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
    if r0 > 0:
        valid = radius_um > 0
        CR_max = np.max(r0 / radius_um[valid]) if np.any(valid) else 0.0
    else:
        CR_max = 0.0

    # Total yield from volumetric burn rate
    E_fusion = 17.6 * 1.60218e-13   # MeV → J
    reaction_rate_total = burn_rate * volume_cm3
    total_reactions = simpson(reaction_rate_total, x=time_s)
    yield_J = total_reactions * E_fusion
    yield_MJ = yield_J * 1e-6

    # Normalized burn fraction (for plotting)
    br_sum = np.sum(burn_rate)
    burn_fraction = burn_rate / br_sum if br_sum > 0 else burn_rate

    # Min radius (skip zeros)
    valid_r = radius_um[radius_um > 0]
    min_radius = np.min(valid_r) if len(valid_r) > 0 else 0.0

    return {
        # Burn-averaged values
        'T_burn_avg':    T_burn_avg,
        'P_burn_avg':    P_burn_avg,
        'rho_burn_avg':  rho_burn_avg,
        'rhoR_burn_avg': rhoR_burn_avg,
        # Peak values
        'T_peak':    np.max(temp_keV),
        'P_peak':    np.max(pres_Gbar),
        'rho_peak':  np.max(dens_gcc),
        'rhoR_peak': np.max(rhoR_gcm2),
        # Convergence and yield
        'CR_max':          CR_max,
        'min_radius_um':   min_radius,
        'yield_MJ':        yield_MJ,
        'yield_J':         yield_J,
        'total_reactions':  total_reactions,
        # Diagnostic
        'burn_fraction': burn_fraction,
        'burn_rate':     burn_rate,
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
        Calculated metrics from calculate_burn_averaged_metrics().
    published_metrics : dict
        Published values with uncertainties.
        Format: {metric_name: (value, uncertainty)}
        Supported keys: 'T_hs', 'P_hs', 'rhoR_cf', 'CR_max', 'yield', 'gain'
    laser_energy_MJ : float
        Laser energy in MJ for gain calculation.

    Returns
    -------
    str
        Formatted comparison table.
    """
    sim_gain = sim_metrics['yield_MJ'] / laser_energy_MJ if laser_energy_MJ > 0 else 0.0
    pub_gain, pub_gain_unc = published_metrics.get('gain', (0.0, 0.0))

    lines = []
    lines.append("=" * 80)
    lines.append("COMPARISON WITH PUBLISHED DATA")
    lines.append("=" * 80)
    lines.append(f"{'Metric':<30} {'Simulation':<20} {'Published':<20} {'Δ (%)':<10}")
    lines.append("-" * 80)

    rows = [
        ('⟨T_hs⟩ (keV)',    sim_metrics['T_burn_avg'],    'T_hs',   '.1f', '.1f'),
        ('⟨P_hs⟩ (Gbar)',   sim_metrics['P_burn_avg'],    'P_hs',   '.0f', '.0f'),
        ('⟨ρR_cf⟩ (g/cm²)', sim_metrics['rhoR_burn_avg'], 'rhoR_cf','.2f', '.2f'),
        ('CR_max',           sim_metrics['CR_max'],         'CR_max', '.1f', '.1f'),
        ('Yield (MJ)',       sim_metrics['yield_MJ'],       'yield',  '.1f', '.1f'),
        ('Fusion Gain',      sim_gain,                      'gain',   '.1f', '.1f'),
    ]

    for label, sim_val, pub_key, sfmt, pfmt in rows:
        pub_val, pub_unc = published_metrics.get(pub_key, (0.0, 0.0))
        if pub_val > 0:
            delta = 100 * (sim_val - pub_val) / pub_val
            sv = format(sim_val, sfmt)
            pv = f"{format(pub_val, pfmt)}±{format(pub_unc, pfmt)}"
            lines.append(f"{label:<30} {sv:>19} {pv:>19} {delta:>9.1f}")

    lines.append("=" * 80)
    return "\n".join(lines)
