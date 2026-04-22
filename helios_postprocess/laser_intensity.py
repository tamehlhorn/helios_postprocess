"""
laser_intensity.py
==================
Reconstruct laser intensity I(r, t) [W/cm^2] from quantities on ICFRunData.

Two independent methods, both computed for cross-check:

  Method 1 (exact, local):
      I(r) = LaserPwrSrc(r) / laserAttinuationCoeff(r)
  Derived from P_src = alpha * I by definition. Undefined where alpha is
  below the coronal-noise floor (ALPHA_MIN_M1) or above the opaque cap.

  Method 2 (Beer-Lambert, 1D spherical):
      I(r) = I_outer * (R_out / r)^2 * exp(-tau(r))
  with
      I_outer = LaserPwrOnTargetForBeam / (4 pi R_out^2)
      tau(r)  = integral_r^{R_out} alpha(r') dr' at zone centers

Extracted from the standalone plot_laser_intensity.py (Apr 2026).

Key change from the standalone:
  - M1 filter: the per-timestep relative threshold `alpha > 0.01 * alpha_max(t)`
    was over-aggressive. alpha_max at the critical surface reaches 10^4-10^5
    cm^-1, so 1% of that is 100+ cm^-1 -- far above the physical coronal
    absorption range (10^-2 to 10^1 cm^-1). Replaced with an absolute physical
    floor ALPHA_MIN_M1 = 1e-2 cm^-1 that excludes only tenuous-vacuum noise
    while keeping the full absorbing corona.

Author: Prof T / Xcimer ICF Analysis
Date:   2026
"""
from __future__ import annotations
import numpy as np
import warnings
from typing import Optional

# 3w Nd:glass baseline; ncr = N_CR_COEFF / lambda_um^2
LAMBDA_UM_DEFAULT = 0.351
N_CR_COEFF = 1.115e21   # cm^-3 * um^2

# Opaque-sentinel handling
SENTINEL_THRESHOLD = 1.0e20
ALPHA_MAX = 1.0e6       # cap that saturates tau on any physical path length
OPAQUE_FLOOR = 0.9 * ALPHA_MAX

# Coronal-noise floor for Method 1 (see module docstring)
ALPHA_MIN_M1 = 1.0e-2   # cm^-1


# ---------------------------------------------------------------------------
# Attenuation cleanup
# ---------------------------------------------------------------------------

def clean_attenuation(atten_raw: np.ndarray, nzone: int) -> np.ndarray:
    """
    Convert raw `laserAttinuationCoeff` to zone-centered, cleaned alpha.

    Handles:
      - Beam axis (takes beam 0 if shape is (nt, nbeam, nbnd))
      - Leading ghost at r=0 if nbnd == nzone+2
      - Opaque sentinels (~1e30) capped at ALPHA_MAX
      - Boundary -> zone-center averaging

    Returns
    -------
    alpha_zone : (nt, nzone) float, [1/cm]
    """
    atten_bnd = atten_raw[:, 0, :] if atten_raw.ndim == 3 else np.asarray(atten_raw)

    nbnd = atten_bnd.shape[1]
    if nbnd == nzone + 2:
        atten_bnd = atten_bnd[:, 1:]
    elif nbnd == nzone + 1:
        pass
    else:
        warnings.warn(f'laserAttinuationCoeff has {nbnd} points along last axis; '
                      f'expected {nzone+1}. Truncating.')
        atten_bnd = atten_bnd[:, :nzone + 1]

    opaque_mask = atten_bnd > SENTINEL_THRESHOLD
    atten_bnd = np.where(opaque_mask, ALPHA_MAX, atten_bnd)
    return 0.5 * (atten_bnd[:, :-1] + atten_bnd[:, 1:])


# ---------------------------------------------------------------------------
# Intensity reconstruction
# ---------------------------------------------------------------------------

def compute_method1(pwr_src: np.ndarray, alpha_zone: np.ndarray) -> np.ndarray:
    """
    I = P_src / alpha, NaN where not physically meaningful.

    Invalid zones: alpha <= ALPHA_MIN_M1 (coronal noise) or alpha >= OPAQUE_FLOOR
    (laser cannot propagate). See module docstring for threshold rationale.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        valid = (alpha_zone > ALPHA_MIN_M1) & (alpha_zone < OPAQUE_FLOOR)
        I = np.where(valid, pwr_src / alpha_zone, np.nan)
    return I


def compute_method2(zbnd: np.ndarray,
                    alpha_zone: np.ndarray,
                    I_outer: np.ndarray) -> np.ndarray:
    """
    I(r, t) = I_outer(t) * (R_out(t) / r)^2 * exp(-tau(r, t))
    tau evaluated at zone centers via inward cumulative sum of alpha*dr.
    """
    dr = np.diff(zbnd, axis=1)
    zcen = 0.5 * (zbnd[:, :-1] + zbnd[:, 1:])
    R_out = zbnd[:, -1]

    alpha_dr = alpha_zone * dr
    # Cumulative optical depth from outer boundary inward, evaluated at zone
    # outer edges, then shifted + corrected to give tau at zone centers.
    cum_out = np.cumsum(alpha_dr[:, ::-1], axis=1)[:, ::-1]
    tau_outer_edge = np.zeros_like(cum_out)
    tau_outer_edge[:, :-1] = cum_out[:, 1:]
    tau_center = tau_outer_edge + 0.5 * alpha_dr

    r_safe = np.where(zcen > 1e-12, zcen, 1e-12)
    geom = (R_out[:, None] / r_safe) ** 2
    return I_outer[:, None] * geom * np.exp(-tau_center)


# ---------------------------------------------------------------------------
# Critical surface
# ---------------------------------------------------------------------------

def find_critical_radius_from_ne(ne: np.ndarray,
                                 zcen: np.ndarray,
                                 lam_um: float) -> tuple:
    """r_crit = outermost zone where ne >= ncr. Returns (r_crit[t] cm, ncr cm^-3)."""
    ncr = N_CR_COEFF / (lam_um ** 2)
    nt = ne.shape[0]
    r_crit = np.full(nt, np.nan)
    for t in range(nt):
        idx = np.where(ne[t] >= ncr)[0]
        if idx.size > 0:
            r_crit[t] = zcen[t, int(idx.max())]
    return r_crit, ncr


def find_critical_radius_fallback(alpha_zone: np.ndarray,
                                  pwr_src: np.ndarray,
                                  zcen: np.ndarray,
                                  lam_um: float) -> tuple:
    """
    Fallback when elec_density is absent. Strategies:
      (a) outermost opaque-sentinel zone + 1
      (b) innermost zone with significant P_src deposition
    """
    nt, nz = alpha_zone.shape
    r_crit = np.full(nt, np.nan)
    for t in range(nt):
        opaque = np.where(alpha_zone[t] > OPAQUE_FLOOR)[0]
        if opaque.size > 0:
            idx = int(opaque.max()) + 1
            if idx < nz:
                r_crit[t] = zcen[t, idx]
                continue
        p = pwr_src[t]
        if p.max() > 0:
            active = np.where(p > 0.01 * p.max())[0]
            if active.size > 0:
                r_crit[t] = zcen[t, int(active.min())]
    return r_crit, N_CR_COEFF / (lam_um ** 2)


def _I_at_first_transparent_zone(I2: np.ndarray,
                                 zcen: np.ndarray,
                                 r_crit: np.ndarray) -> np.ndarray:
    """
    I2 evaluated at first transparent zone outside r_crit[t].
    Avoids np.interp across the exp(-tau) cliff.
    """
    nt, nz = I2.shape
    out = np.full(nt, np.nan)
    for t in range(nt):
        if not np.isfinite(r_crit[t]):
            continue
        idx = int(np.argmin(np.abs(zcen[t] - r_crit[t])))
        while idx < nz - 1 and not (np.isfinite(I2[t, idx]) and I2[t, idx] > 0):
            idx += 1
        if np.isfinite(I2[t, idx]) and I2[t, idx] > 0:
            out[t] = I2[t, idx]
    return out


# ---------------------------------------------------------------------------
# Top-level driver
# ---------------------------------------------------------------------------

def analyze_laser_intensity(data, wavelength_um: float = LAMBDA_UM_DEFAULT
                            ) -> Optional[dict]:
    """
    Compute laser intensity diagnostics from an ICFRunData instance.

    Required fields on `data`:
      - laser_power_source           (nt, nz)                 [W/cm^3]
      - laser_attenuation_coeff      raw (nt, nbeam, nbnd)    [1/cm]
      - laser_power_on_target        (nt,)                    [W]
      - zone_boundaries              (nt, nz+1)               [cm]
    Optional:
      - electron_density             (nt, nz)                 [1/cm^3]
        (fallback used if absent)

    Returns None if required inputs are missing.

    Returns dict with:
      I1, I2                     (nt, nz)    Method-1 and Method-2 intensities
      I_grid_outer               (nt,)       Incident intensity at grid outer boundary
      alpha_zone                 (nt, nz)    Cleaned, zone-centered attenuation
      r_crit                     (nt,)       Critical-surface radius [cm]
      ncr                        float       [1/cm^3]
      I_at_crit_vs_t             (nt,)       I at first transparent zone outside r_crit
      I_peak_coronal_vs_t        (nt,)       max I2(r) over zones
      peak_I_grid_outer          float       max over t (at grid outer boundary)
      peak_I_at_crit             float       max over t
      peak_I_coronal             float       max over (t, r)
      I_at_crit_at_peak_power    float       I_at_crit at t = argmax(P_on_target)
      t_peak_power_ns            float       time of peak laser power
      wavelength_um              float
    """
    pwr_src = getattr(data, 'laser_power_source', None)
    atten_raw = getattr(data, 'laser_attenuation_coeff', None)
    P_on_target = getattr(data, 'laser_power_on_target', None)
    zbnd = getattr(data, 'zone_boundaries', None)
    if pwr_src is None or atten_raw is None or P_on_target is None or zbnd is None:
        return None

    nt, nz = pwr_src.shape
    zcen = 0.5 * (zbnd[:, :-1] + zbnd[:, 1:])
    R_out = zbnd[:, -1]

    alpha_zone = clean_attenuation(atten_raw, nz)

    with np.errstate(divide='ignore', invalid='ignore'):
        I_outer = np.where(R_out > 0,
                           P_on_target / (4.0 * np.pi * R_out ** 2),
                           0.0)

    I1 = compute_method1(pwr_src, alpha_zone)
    I2 = compute_method2(zbnd, alpha_zone, I_outer)

    ne = getattr(data, 'electron_density', None)
    if ne is not None:
        r_crit, ncr = find_critical_radius_from_ne(ne, zcen, wavelength_um)
    else:
        r_crit, ncr = find_critical_radius_fallback(
            alpha_zone, pwr_src, zcen, wavelength_um)

    I_at_crit_vs_t = _I_at_first_transparent_zone(I2, zcen, r_crit)

    I_peak_coronal_vs_t = np.full(nt, np.nan)
    with np.errstate(invalid='ignore'):
        for t in range(nt):
            row = I2[t]
            mask = np.isfinite(row) & (row > 0)
            if mask.any():
                I_peak_coronal_vs_t[t] = np.nanmax(row[mask])

    peak_I_outer = float(np.nanmax(I_outer)) if I_outer.size else np.nan
    peak_I_at_crit = (float(np.nanmax(I_at_crit_vs_t))
                      if np.any(np.isfinite(I_at_crit_vs_t)) else np.nan)
    peak_I_coronal = (float(np.nanmax(I_peak_coronal_vs_t))
                      if np.any(np.isfinite(I_peak_coronal_vs_t)) else np.nan)

    t_peak = int(np.argmax(P_on_target)) if P_on_target.size else 0
    I_at_crit_at_peak_power = (float(I_at_crit_vs_t[t_peak])
                               if np.isfinite(I_at_crit_vs_t[t_peak]) else np.nan)
    t_peak_ns = (float(data.time[t_peak])
                 if hasattr(data, 'time') and data.time is not None else np.nan)

    return dict(
        I1=I1, I2=I2, I_grid_outer=I_outer, alpha_zone=alpha_zone,
        r_crit=r_crit, ncr=ncr,
        I_at_crit_vs_t=I_at_crit_vs_t,
        I_peak_coronal_vs_t=I_peak_coronal_vs_t,
        peak_I_grid_outer=peak_I_outer,
        peak_I_at_crit=peak_I_at_crit,
        peak_I_coronal=peak_I_coronal,
        I_at_crit_at_peak_power=I_at_crit_at_peak_power,
        t_peak_power_ns=t_peak_ns,
        wavelength_um=wavelength_um,
    )
