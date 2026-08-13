"""
Level 0 x-ray emissivity from Helios hydro state.

Physics content (deliberately minimal and fully inspectable):

  * thermal free-free (bremsstrahlung) continuum emission
  * thermally-averaged non-relativistic Born/Elwert Gaunt factor
  * LTE continuum absorption via Kirchhoff's law (kappa = j / B_E)

Deliberately ABSENT (these are exactly what the SPECT3D rungs of the ladder
supply, and what the Level 0 <-> Level 1 comparison is meant to quantify):

  * free-bound (radiative recombination) continuum and its edges
  * bound-bound line emission and absorption
  * non-LTE / collisional-radiative ionization balance
  * opacity-table (rather than pure-continuum) absorption

Conventions
-----------
Photon energy grid ``E`` is in eV.
Spectral emissivity ``j_E`` is per unit photon energy per steradian:
    [erg s^-1 cm^-3 eV^-1 sr^-1]
Absorption coefficient ``kappa_E`` is [cm^-1].
Specific intensity downstream is [erg s^-1 cm^-2 sr^-1 eV^-1].

Author: Prof T / helios_postprocess
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.special import kv

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants (CGS unless noted)
# ---------------------------------------------------------------------------
H_PLANCK_EVS = 4.135667696e-15      # eV s
H_PLANCK_ERGS = 6.62607015e-27      # erg s
C_LIGHT = 2.99792458e10             # cm/s
EV_PER_KELVIN = 8.617333262e-5      # eV/K
K_TO_EV = EV_PER_KELVIN
EV_TO_K = 1.0 / EV_PER_KELVIN       # 11604.5 K per eV
ERG_PER_EV = 1.602176634e-12

# Rybicki & Lightman (5.14b): 4pi-integrated free-free emission coefficient
# eps_nu = 6.8e-38 Z^2 n_e n_i T^{-1/2} exp(-h nu / kT) gbar
#          [erg s^-1 cm^-3 Hz^-1], T in Kelvin.
FF_COEFF = 6.8e-38

# Numerical floors
_T_FLOOR_EV = 1.0e-3                # 1 meV; below this a zone emits nothing
_RHO_FLOOR = 1.0e-12                # g/cm^3


# ---------------------------------------------------------------------------
# Gaunt factor
# ---------------------------------------------------------------------------
def gaunt_ff(E_eV: np.ndarray, Te_eV: np.ndarray) -> np.ndarray:
    """
    Thermally-averaged free-free Gaunt factor, non-relativistic Born limit.

        gbar(u) = (sqrt(3)/pi) * exp(u/2) * K_0(u/2),     u = E / kTe

    Accurate to a few percent over the range relevant to ICF coronal and
    hot-spot continuum emission (u ~ 0.01 - 10).  Asymptotes correctly to
    the logarithmic small-u form and decays for u >> 1.

    Parameters
    ----------
    E_eV : array
        Photon energies, eV.  Broadcast against ``Te_eV``.
    Te_eV : array
        Electron temperature, eV.

    Returns
    -------
    array
        Gaunt factor, same broadcast shape.
    """
    u = np.asarray(E_eV, dtype=float) / np.maximum(np.asarray(Te_eV, dtype=float),
                                                   _T_FLOOR_EV)
    u = np.clip(u, 1.0e-8, 7.0e2)
    half = 0.5 * u
    # exp(u/2) * K0(u/2) computed via the exponentially-scaled Bessel to
    # avoid overflow/underflow at large u.  kv(0, x) * exp(x) is stable.
    g = (np.sqrt(3.0) / np.pi) * kv(0, half) * np.exp(half)
    return np.nan_to_num(g, nan=1.0, posinf=1.0, neginf=1.0)


# ---------------------------------------------------------------------------
# Planck function (source function under LTE continuum)
# ---------------------------------------------------------------------------
def planck_E(E_eV: np.ndarray, Te_eV: np.ndarray) -> np.ndarray:
    """
    Planck specific intensity per unit photon energy.

    Returns [erg s^-1 cm^-2 sr^-1 eV^-1].
    """
    E_eV = np.asarray(E_eV, dtype=float)
    Te = np.maximum(np.asarray(Te_eV, dtype=float), _T_FLOOR_EV)
    nu = E_eV / H_PLANCK_EVS                       # Hz
    x = np.clip(E_eV / Te, 1.0e-10, 7.0e2)
    B_nu = (2.0 * H_PLANCK_ERGS * nu**3 / C_LIGHT**2) / np.expm1(x)
    return B_nu / H_PLANCK_EVS                      # per eV


# ---------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------
@dataclass
class EmissivityCube:
    """
    Level 0 emissivity / opacity for one time step.

    Attributes
    ----------
    E_eV : (n_E,)
        Photon energy grid, eV.
    j_E : (n_zone, n_E)
        Spectral emissivity, erg s^-1 cm^-3 eV^-1 sr^-1.
    kappa_E : (n_zone, n_E)
        LTE continuum absorption coefficient, cm^-1.
    source_E : (n_zone, n_E)
        Source function j/kappa = B_E, erg s^-1 cm^-2 sr^-1 eV^-1.
    r_bnd : (n_zone+1,)
        Zone boundary radii, cm.
    time_ns : float
        Simulation time of this snapshot.
    """
    E_eV: np.ndarray
    j_E: np.ndarray
    kappa_E: np.ndarray
    source_E: np.ndarray
    r_bnd: np.ndarray
    time_ns: float


# ---------------------------------------------------------------------------
# Zbar resolution
# ---------------------------------------------------------------------------
def resolve_zbar(data, it: Optional[int] = None) -> np.ndarray:
    """
    Get mean ionization Zbar per zone, with graceful fallbacks.

    Order of preference:
      1. ``data.mean_charge``            (EXODUS 'mean_charge', preferred)
      2. ``data.electron_density / data.ion_density``
      3. constant 1.0 with a loud warning

    Returns the full (n_t, n_zone) array when ``it`` is None, else (n_zone,).
    """
    zbar = getattr(data, "mean_charge", None)

    if zbar is None:
        ne = getattr(data, "electron_density", None)
        ni = getattr(data, "ion_density", None)
        if ne is not None and ni is not None:
            zbar = np.asarray(ne, float) / np.maximum(np.asarray(ni, float), 1.0)
            logger.warning("xray: mean_charge absent -- Zbar from n_e/n_i")
        else:
            logger.warning("xray: no mean_charge and no n_e/n_i -- Zbar set to 1.0. "
                           "Free-free emission will be badly underestimated in "
                           "the ablator/corona. Wire 'mean_charge' into "
                           "data_builder._VARIABLE_MAP.")
            rho = np.asarray(data.mass_density, float)
            zbar = np.ones_like(rho)

    zbar = np.asarray(zbar, float)
    return zbar if it is None else zbar[it]


def resolve_ion_density(data, it: Optional[int] = None) -> np.ndarray:
    """Ion number density [cm^-3], falling back to n_e / Zbar."""
    ni = getattr(data, "ion_density", None)
    if ni is not None:
        ni = np.asarray(ni, float)
        return ni if it is None else ni[it]

    ne = getattr(data, "electron_density", None)
    if ne is None:
        raise ValueError("Neither ion_density nor electron_density available; "
                         "cannot build emissivity.")
    zbar = resolve_zbar(data)
    ni = np.asarray(ne, float) / np.maximum(zbar, 1.0e-3)
    logger.warning("xray: ion_density absent -- using n_e / Zbar")
    return ni if it is None else ni[it]


def resolve_electron_density(data, it: Optional[int] = None) -> np.ndarray:
    """Electron number density [cm^-3], falling back to Zbar * n_i."""
    ne = getattr(data, "electron_density", None)
    if ne is not None:
        ne = np.asarray(ne, float)
        return ne if it is None else ne[it]
    ni = resolve_ion_density(data)
    zbar = resolve_zbar(data)
    ne = ni * zbar
    logger.warning("xray: electron_density absent -- using Zbar * n_i")
    return ne if it is None else ne[it]


def resolve_Te(data, it: Optional[int] = None) -> np.ndarray:
    """
    Electron temperature [eV].  Falls back to ion temperature with a warning:
    in the corona (where most of the sub-5 keV self-emission originates)
    Te and Ti decouple, so this fallback is only acceptable for a smoke test.
    """
    Te = getattr(data, "elec_temperature", None)
    if Te is None:
        logger.warning("xray: elec_temperature absent -- falling back to "
                       "ion_temperature. Coronal emission will be wrong.")
        Te = data.ion_temperature
    Te = np.asarray(Te, float)
    return Te if it is None else Te[it]


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------
def build_emissivity(data,
                     it: int,
                     E_eV: np.ndarray,
                     *,
                     zone_slice: Optional[slice] = None,
                     include_absorption: bool = True) -> EmissivityCube:
    """
    Build the Level 0 emissivity / absorption cube for one Helios time index.

    Parameters
    ----------
    data : ICFRunData
        Populated container from ``helios_postprocess.data_builder``.
    it : int
        Time index into ``data.time``.
    E_eV : (n_E,) array
        Photon energy grid in eV.
    zone_slice : slice, optional
        Restrict to a zone range (e.g. capsule-only for halfraum targets,
        excluding the Cu converter and He fill).
    include_absorption : bool
        If True, compute ``kappa_E`` from Kirchhoff (LTE continuum).
        If False, ``kappa_E`` is zeros and the chord solver runs optically thin.

    Returns
    -------
    EmissivityCube
    """
    E_eV = np.atleast_1d(np.asarray(E_eV, dtype=float))

    Te = resolve_Te(data, it)
    ne = resolve_electron_density(data, it)
    ni = resolve_ion_density(data, it)
    zbar = resolve_zbar(data, it)
    r_bnd = np.asarray(data.zone_boundaries[it], float)

    if zone_slice is not None:
        Te, ne, ni, zbar = Te[zone_slice], ne[zone_slice], ni[zone_slice], zbar[zone_slice]
        start = 0 if zone_slice.start is None else zone_slice.start
        stop = len(r_bnd) - 1 if zone_slice.stop is None else zone_slice.stop
        r_bnd = r_bnd[start:stop + 1]

    Te = np.maximum(Te, _T_FLOOR_EV)
    T_K = Te * EV_TO_K

    # (n_zone, n_E) broadcast
    Ez = E_eV[None, :]
    Tz = Te[:, None]

    g = gaunt_ff(Ez, Tz)
    expo = np.exp(-np.clip(Ez / Tz, 0.0, 7.0e2))

    # 4pi-integrated, per Hz -> per sr, per eV
    j_nu = (FF_COEFF * (zbar[:, None] ** 2) * ne[:, None] * ni[:, None]
            * (T_K[:, None] ** -0.5) * expo * g)
    j_E = j_nu / (4.0 * np.pi) / H_PLANCK_EVS

    j_E = np.nan_to_num(j_E, nan=0.0, posinf=0.0, neginf=0.0)

    B_E = planck_E(Ez, Tz)
    if include_absorption:
        kappa_E = np.where(B_E > 0.0, j_E / np.maximum(B_E, 1.0e-300), 0.0)
        kappa_E = np.nan_to_num(kappa_E, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        kappa_E = np.zeros_like(j_E)

    return EmissivityCube(E_eV=E_eV,
                          j_E=j_E,
                          kappa_E=kappa_E,
                          source_E=B_E,
                          r_bnd=r_bnd,
                          time_ns=float(np.asarray(data.time)[it]))
