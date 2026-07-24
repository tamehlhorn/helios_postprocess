"""
Synthetic neutron time-of-flight (nTOF) from a Helios birth spectrum.
====================================================================

Thin, analytic forward model that turns the synthetic *birth* spectrum produced
by :mod:`helios_postprocess.neutron_spectrum` into the signal an nTOF detector at
distance ``d`` would record. This is the display/observable layer: the primary
DT peak gives the yield (area) and the ion temperature (width); the
down-scattered neutrons — once added by the scatter model — land on the same
time axis as the late-time tail.

Physics (relativistic; matches Gopalaswamy et al. 2025, Eqs. 2-4)
----------------------------------------------------------------
- neutron speed   v(E) = c·β,  β = sqrt(1 − 1/γ²),  γ = 1 + E/(m_n c²)
- time of flight  t(E) = t0 + d / v(E)
- signal          dN/dt = dN/dE · |dE/dt|      (change of variables; conserves counts)

A 14.06 MeV neutron moves at 5.13 cm/ns → 58.5 ns at 3 m. The primary-peak
*time* width is the nTOF ion temperature: for a narrow peak Δt/t ≈ ½ ΔE/E, and
ΔE_FWHM[keV] = 177·√(Ti[keV]) for DT (Brysk).

Author: Prof T (helios_postprocess)
"""

from __future__ import annotations

import logging
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.integrate import trapezoid

logger = logging.getLogger(__name__)

#: neutron rest energy (MeV) and speed of light (cm/ns).
M_N_C2_MEV = 939.56542
C_CM_NS = 29.9792458

# Brysk coefficients reused from neutron_spectrum (kept local to avoid a hard
# import cycle); FWHM[keV] = coeff · sqrt(Ti[keV]).
_BRYSK = {"DT": 177.0, "DD": 82.5}
_E0 = {"DT": 14.06, "DD": 2.45}
_FWHM_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


# ---------------------------------------------------------------------------
# Kinematics
# ---------------------------------------------------------------------------

def neutron_velocity(E_MeV: np.ndarray) -> np.ndarray:
    """Relativistic neutron speed (cm/ns) for kinetic energy ``E_MeV``."""
    E = np.asarray(E_MeV, dtype=float)
    gamma = 1.0 + E / M_N_C2_MEV
    beta = np.sqrt(np.clip(1.0 - 1.0 / gamma ** 2, 0.0, None))
    return beta * C_CM_NS


def energy_to_tof(E_MeV: np.ndarray, distance_m: float, t0_ns: float = 0.0) -> np.ndarray:
    """Arrival time (ns) of energy ``E_MeV`` neutrons at ``distance_m`` metres."""
    d_cm = distance_m * 100.0
    return t0_ns + d_cm / neutron_velocity(E_MeV)


def tof_to_energy(t_ns: np.ndarray, distance_m: float, t0_ns: float = 0.0) -> np.ndarray:
    """Kinetic energy (MeV) implied by arrival time ``t_ns`` (inverse of
    :func:`energy_to_tof`). NaN for ``t <= t0``."""
    t = np.asarray(t_ns, dtype=float)
    d_cm = distance_m * 100.0
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(t > t0_ns, d_cm / (t - t0_ns), np.nan)      # cm/ns
        beta = v / C_CM_NS
        beta = np.where(beta < 1.0, beta, np.nan)
        gamma = 1.0 / np.sqrt(1.0 - beta ** 2)
    return (gamma - 1.0) * M_N_C2_MEV


def tion_from_tof_width(t_peak_ns: float, fwhm_t_ns: float,
                        E0_MeV: float = 14.06, reaction: str = "DT") -> float:
    """Ion temperature (keV) from the primary-peak time width.

    Δt/t ≈ ½ ΔE/E  ⇒  ΔE_FWHM = 2·E·(Δt/t);  then invert Brysk
    Ti = (ΔE_FWHM[keV] / coeff)². Valid for a narrow, isolated primary peak.
    """
    if t_peak_ns <= 0:
        return float("nan")
    dE_keV = 2.0 * E0_MeV * (fwhm_t_ns / t_peak_ns) * 1000.0
    return (dE_keV / _BRYSK[reaction]) ** 2


# ---------------------------------------------------------------------------
# Spectrum -> TOF signal
# ---------------------------------------------------------------------------

def spectrum_to_tof(energy_MeV: np.ndarray, dNdE: np.ndarray,
                    distance_m: float, t0_ns: float = 0.0
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """Transform an energy spectrum ``dN/dE`` into the TOF signal ``dN/dt``
    (counts/ns) at ``distance_m``. Returns ``(t_ns, dNdt)`` sorted ascending in
    time. The change of variables conserves counts:
    ``∫ dN/dt dt == ∫ dN/dE dE``.
    """
    E = np.asarray(energy_MeV, dtype=float)
    dNdE = np.asarray(dNdE, dtype=float)
    t = energy_to_tof(E, distance_m, t0_ns)
    dtdE = np.gradient(t, E)                                     # < 0 (t falls as E rises)
    with np.errstate(divide="ignore", invalid="ignore"):
        dNdt = np.where(np.abs(dtdE) > 0, dNdE / np.abs(dtdE), 0.0)
    order = np.argsort(t)
    return t[order], dNdt[order]


def resample_uniform(t_ns: np.ndarray, y: np.ndarray, n_t: int = 4000
                     ) -> Tuple[np.ndarray, np.ndarray]:
    """Interpolate ``(t, y)`` onto a uniform time grid of ``n_t`` points
    (needed for IRF convolution / plotting)."""
    t = np.asarray(t_ns, dtype=float)
    grid = np.linspace(t.min(), t.max(), n_t)
    return grid, np.interp(grid, t, y, left=0.0, right=0.0)


def apply_irf(t_uniform_ns: np.ndarray, signal: np.ndarray,
              irf_fwhm_ns: float) -> np.ndarray:
    """Convolve a signal on a *uniform* time grid with a Gaussian instrument
    response of FWHM ``irf_fwhm_ns`` (area-preserving)."""
    t = np.asarray(t_uniform_ns, dtype=float)
    if irf_fwhm_ns <= 0 or t.size < 2:
        return np.asarray(signal, dtype=float)
    dt = t[1] - t[0]
    sigma = irf_fwhm_ns / _FWHM_SIGMA
    half = int(np.ceil(4 * sigma / dt))
    k = np.arange(-half, half + 1) * dt
    kern = np.exp(-0.5 * (k / sigma) ** 2)
    kern /= kern.sum()
    return np.convolve(np.asarray(signal, dtype=float), kern, mode="same")


def primary_peak_metrics(t_ns: np.ndarray, dNdt: np.ndarray, distance_m: float,
                         reaction: str = "DT") -> Dict:
    """Arrival time, time-FWHM, area, and nTOF ion temperature of the primary
    peak. The search window is centred on the expected primary arrival."""
    t = np.asarray(t_ns, dtype=float)
    y = np.asarray(dNdt, dtype=float)
    t_expect = float(energy_to_tof(_E0[reaction], distance_m))
    # window: +/-40% around the expected primary arrival
    win = (t > 0.6 * t_expect) & (t < 1.4 * t_expect)
    if not win.any() or y[win].max() <= 0:
        return {"t_peak_ns": float("nan"), "fwhm_ns": float("nan"),
                "area": float("nan"), "Ti_ntof_keV": float("nan")}
    tw, yw = t[win], y[win]
    ymax = yw.max()
    # centroid (robust peak position)
    t_peak = float(trapezoid(tw * yw, tw) / trapezoid(yw, tw))
    fwhm = _fwhm(tw, yw)
    area = float(trapezoid(yw, tw))
    # Ti_tof_width is the linear-approx (Δt/t ≈ ½ΔE/E) quick-look; it under-reads
    # for a non-Gaussian composite peak. The authoritative nTOF Ti comes from the
    # spectral Doppler width (see synthetic_ntof), not this.
    Ti_approx = tion_from_tof_width(t_peak, fwhm, _E0[reaction], reaction)
    return {"t_peak_ns": t_peak, "fwhm_ns": fwhm, "area": area,
            "Ti_tof_width_keV": Ti_approx,
            "t_expected_ns": t_expect, "peak_height": float(ymax)}


def _fwhm(x: np.ndarray, y: np.ndarray) -> float:
    ymax = y.max()
    if ymax <= 0:
        return float("nan")
    half = 0.5 * ymax
    idx = np.where(y >= half)[0]
    if idx.size == 0:
        return float("nan")
    i0, i1 = idx[0], idx[-1]
    xl = np.interp(half, [y[i0 - 1], y[i0]], [x[i0 - 1], x[i0]]) if i0 > 0 else x[i0]
    xr = np.interp(half, [y[i1 + 1], y[i1]], [x[i1 + 1], x[i1]]) if i1 < len(x) - 1 else x[i1]
    return float(abs(xr - xl))


# ---------------------------------------------------------------------------
# High-level entry point
# ---------------------------------------------------------------------------

def synthetic_ntof(data=None, energy_MeV: Optional[np.ndarray] = None,
                   spectrum: Optional[np.ndarray] = None,
                   distance_m: float = 3.0, t0_ns: float = 0.0,
                   irf_fwhm_ns: float = 0.0, n_t: int = 4000,
                   reaction: str = "DT") -> Optional[Dict]:
    """Synthetic nTOF trace from either a birth spectrum (``energy_MeV`` +
    ``spectrum``) or an ``ICFRunData`` (``data`` — the birth spectrum is built
    via :mod:`neutron_spectrum`).

    Returns a dict with the ideal and (optionally) IRF-convolved signals on a
    uniform time grid, plus primary-peak metrics. ``None`` on a no-burn run.
    """
    from . import neutron_spectrum as ns
    if energy_MeV is None or spectrum is None:
        if data is None:
            raise ValueError("Provide either (energy_MeV, spectrum) or data.")
        nd = ns.extract_neutronics(data=data, use_rhino=False)
        if nd is None:
            logger.warning("synthetic_ntof: no birth spectrum (no-burn run).")
            return None
        energy_MeV = nd.quicklook["energy_MeV"]
        spectrum = nd.quicklook["birth_spectrum"]

    energy_MeV = np.asarray(energy_MeV, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    t_raw, dNdt_raw = spectrum_to_tof(energy_MeV, spectrum, distance_m, t0_ns)
    t, sig_ideal = resample_uniform(t_raw, dNdt_raw, n_t)
    sig_det = apply_irf(t, sig_ideal, irf_fwhm_ns) if irf_fwhm_ns > 0 else sig_ideal
    metrics = primary_peak_metrics(t_raw, dNdt_raw, distance_m, reaction)
    # Authoritative nTOF ion temperature = spectral Doppler width (second moment),
    # the quantity a proper peak fit recovers; the TOF-width value is only a
    # linear-approx quick-look.
    metrics["Ti_ntof_keV"] = ns.infer_ntof_tion(energy_MeV, spectrum, reaction)["Ti_ntof_keV"]

    logger.info("nTOF @ %.1f m: primary peak %.2f ns (expected %.2f), FWHM %.3f ns, "
                "nTOF Ti %.2f keV (TOF-width approx %.2f keV)", distance_m,
                metrics["t_peak_ns"], metrics.get("t_expected_ns", float("nan")),
                metrics["fwhm_ns"], metrics["Ti_ntof_keV"], metrics["Ti_tof_width_keV"])
    return {"time_ns": t, "signal_ideal": sig_ideal, "signal_detector": sig_det,
            "distance_m": distance_m, "irf_fwhm_ns": irf_fwhm_ns,
            "primary": metrics}
