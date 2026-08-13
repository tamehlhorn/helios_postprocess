"""
Spherical chord integration for 1D Helios runs.

A 1D spherically-symmetric snapshot is imaged by integrating along parallel
chords of impact parameter p.  Because the geometry is analytic, no ray
tracer is needed: the path length of a chord through the shell bounded by
r_k and r_{k+1} is

    ds_k = sqrt(r_{k+1}^2 - p^2) - sqrt(max(r_k^2 - p^2, 0))

appearing twice (far side and near side of the tangent point).

Two solvers:

  * ``integrate_thin``  -- pure optically-thin sum, I = 2 * sum(ds_k * j_k).
    Fast, exactly what "Level 0" means, and the reference the absorption
    correction is measured against.

  * ``integrate_formal`` -- ordered formal solution of the transfer equation
    with LTE continuum absorption.  Returns the emergent intensity AND the
    chord optical depth, so we can map *where in (p, E, t)* the optically-thin
    assumption fails.  That map is the quantitative hand-off argument to
    SPECT3D.

Author: Prof T / helios_postprocess
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass
class ChordResult:
    """
    Emergent intensity along a set of impact parameters.

    Attributes
    ----------
    p_cm : (n_p,)
        Impact parameters, cm.
    E_eV : (n_E,)
        Photon energy grid, eV.
    I : (n_p, n_E)
        Emergent specific intensity, erg s^-1 cm^-2 sr^-1 eV^-1.
    tau : (n_p, n_E)
        Total chord optical depth (zero when run optically thin).
    I_thin : (n_p, n_E) or None
        Optically-thin intensity, retained for the absorption diagnostic.
    """
    p_cm: np.ndarray
    E_eV: np.ndarray
    I: np.ndarray
    tau: np.ndarray
    I_thin: Optional[np.ndarray] = None


def make_impact_grid(r_max: float, n_p: int = 128,
                     oversample_limb: bool = True) -> np.ndarray:
    """
    Impact-parameter grid from 0 to ``r_max``.

    The emission profile of an imploding capsule is limb-brightened and has a
    sharp edge at the ablation front, so a uniform grid wastes samples in the
    interior and under-resolves the very feature the streak camera measures.
    ``oversample_limb`` bunches points toward large p with a sqrt stretch.
    """
    u = np.linspace(0.0, 1.0, int(n_p))
    if oversample_limb:
        u = 1.0 - np.sqrt(1.0 - u ** 2)      # dense near u -> 1
        u = u / u[-1]
    return u * float(r_max)


def _segment_lengths(r_bnd: np.ndarray, p: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Half-chord path lengths through each zone for impact parameter p.

    Returns
    -------
    ds : (n_zone,)
        Path length through zone k on ONE side of the tangent point (cm).
        Zero for zones entirely inside p.
    active : (n_zone,) bool
        Zones the chord actually crosses.
    """
    s = np.sqrt(np.maximum(r_bnd ** 2 - p ** 2, 0.0))
    ds = s[1:] - s[:-1]
    active = r_bnd[1:] > p
    ds = np.where(active, ds, 0.0)
    return ds, active


def integrate_thin(cube, p_cm: np.ndarray) -> ChordResult:
    """
    Optically-thin chord integration.

    Parameters
    ----------
    cube : EmissivityCube
        From ``helios_postprocess.xray.emissivity.build_emissivity``.
    p_cm : (n_p,) array
        Impact parameters, cm.
    """
    p_cm = np.atleast_1d(np.asarray(p_cm, float))
    r_bnd = cube.r_bnd
    n_p = p_cm.size

    ds_mat = np.empty((n_p, cube.j_E.shape[0]), dtype=float)
    for i, p in enumerate(p_cm):
        ds_mat[i], _ = _segment_lengths(r_bnd, p)

    I = 2.0 * (ds_mat @ cube.j_E)                # (n_p, n_E)
    return ChordResult(p_cm=p_cm, E_eV=cube.E_eV, I=I,
                       tau=np.zeros_like(I), I_thin=I.copy())


def integrate_formal(cube, p_cm: np.ndarray) -> ChordResult:
    """
    Formal solution with LTE continuum absorption.

    Traversal order is far side inward to the tangent zone, then back outward
    to the observer, so the exponential attenuation weights near-side emission
    correctly.  Within each segment the source function is taken constant
    (short-characteristic, first order):

        I <- I * exp(-dtau) + S * (1 - exp(-dtau))

    Also accumulates the total chord optical depth.
    """
    p_cm = np.atleast_1d(np.asarray(p_cm, float))
    r_bnd = cube.r_bnd
    n_E = cube.E_eV.size

    I_out = np.zeros((p_cm.size, n_E), dtype=float)
    tau_out = np.zeros((p_cm.size, n_E), dtype=float)
    I_thin = np.zeros((p_cm.size, n_E), dtype=float)

    for i, p in enumerate(p_cm):
        ds, active = _segment_lengths(r_bnd, p)
        idx = np.flatnonzero(active)
        if idx.size == 0:
            continue

        ds_a = ds[idx]                                   # (m,)
        kap = cube.kappa_E[idx]                          # (m, n_E)
        src = cube.source_E[idx]                         # (m, n_E)
        jem = cube.j_E[idx]                              # (m, n_E)

        I_thin[i] = 2.0 * (ds_a @ jem)

        # Ordered segments: outer -> inner (far side), then inner -> outer.
        order = np.concatenate([idx.size - 1 - np.arange(idx.size),
                                np.arange(idx.size)])
        dtau = (kap * ds_a[:, None])[order]              # (2m, n_E)
        S = src[order]

        dtau = np.clip(dtau, 0.0, 7.0e2)
        tau_tot = dtau.sum(axis=0)
        tau_out[i] = tau_tot

        # Closed form of the sequential march, vectorized over segments:
        #   I = sum_k S_k (1 - e^{-dtau_k}) e^{-tau_ahead_k}
        # where tau_ahead_k is the optical depth between segment k and the
        # observer.  Identical to marching segment by segment, but without
        # the Python loop over ~2*n_zone segments -- the difference is a
        # factor of several hundred in wall time on a 200-zone run, which
        # matters once the sweep window covers a few hundred time steps.
        # -expm1(-dt) avoids the cancellation that costs ~3 digits at
        # dt < 1e-12.
        tau_ahead = tau_tot[None, :] - np.cumsum(dtau, axis=0)
        I_out[i] = np.sum(S * (-np.expm1(-dtau))
                          * np.exp(-np.clip(tau_ahead, 0.0, 7.0e2)), axis=0)

    return ChordResult(p_cm=p_cm, E_eV=cube.E_eV, I=I_out,
                       tau=tau_out, I_thin=I_thin)


def integrate(cube, p_cm: np.ndarray, optically_thin: bool = True) -> ChordResult:
    """Dispatch to the thin or formal solver."""
    if optically_thin:
        return integrate_thin(cube, p_cm)
    return integrate_formal(cube, p_cm)


def spatially_integrate(result: ChordResult) -> np.ndarray:
    """
    Convert an intensity profile I(p, E) into the total emitted spectral power
    per unit photon energy over 4*pi:

        P_E = 4 * pi * integral( I(p) * 2*pi*p dp )

    Returns (n_E,) in erg s^-1 eV^-1.  This is the quantity a streaked
    spectrometer with no spatial resolution records (up to instrument
    solid angle and response).
    """
    p = result.p_cm
    weight = 2.0 * np.pi * p
    return 4.0 * np.pi * np.trapezoid(result.I * weight[:, None], p, axis=0)
