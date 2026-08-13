"""
X-ray instrument response: filter transmission, photocathode quantum
efficiency, and the composite channel response R(E).

IMPORTANT ACCURACY CAVEAT
-------------------------
The built-in cold photoabsorption model in this module is an ANCHORED POWER
LAW, not tabulated data.  It reproduces mass attenuation coefficients to
roughly +/- 25% between ~0.5 and ~10 keV away from edges, and worse within a
few hundred eV of an absorption edge.  That is adequate for scoping a channel
and for the Level 0 <-> SPECT3D comparison (both rungs use the SAME response,
so it cancels in the ratio), but it is NOT adequate for absolute photometric
comparison against measured streak data.

Before any quantitative comparison to data, download CXRO/Henke tables and
load them with :func:`load_henke_file`.  Get them from

    https://henke.lbl.gov/optical_constants/atten2.html

as two-column (energy_eV, mu_over_rho_cm2_per_g) ASCII, one file per
material, dropped in ``helios_postprocess/xray/data/``.  The loader is the
production path; the analytic model exists so nothing blocks on a download.

Author: Prof T / helios_postprocess
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent / "data"

# ---------------------------------------------------------------------------
# Anchored power-law cold photoabsorption
#
#   mu/rho(E) = A * (E / 1.5 keV)^(-2.9)   above the listed K edge
#
# A is the mass attenuation coefficient at 1.5 keV in cm^2/g (or, for
# materials whose K edge lies above 1.5 keV, the value just above the edge
# scaled back to 1.5 keV using the same exponent).  ``jump`` is the K-edge
# absorption jump ratio applied below the edge.
# ---------------------------------------------------------------------------
_EXPONENT = -2.9

_ELEMENTS: Dict[str, Dict[str, float]] = {
    #        A@1.5keV   Z      A_atomic   K-edge eV   jump
    "H":  dict(A=7.0,    Z=1,  M=1.008,   edge=0.0,    jump=1.0),
    "Be": dict(A=200.0,  Z=4,  M=9.012,   edge=111.0,  jump=1.0),
    "C":  dict(A=700.0,  Z=6,  M=12.011,  edge=284.0,  jump=1.0),
    "N":  dict(A=1200.0, Z=7,  M=14.007,  edge=410.0,  jump=1.0),
    "O":  dict(A=1900.0, Z=8,  M=15.999,  edge=543.0,  jump=1.0),
    "Al": dict(A=2800.0, Z=13, M=26.982,  edge=1560.0, jump=7.0),
    "Si": dict(A=3400.0, Z=14, M=28.086,  edge=1839.0, jump=7.5),
    "Ti": dict(A=9000.0, Z=22, M=47.867,  edge=4966.0, jump=6.5),
    "Cu": dict(A=17000.0, Z=29, M=63.546, edge=8979.0, jump=6.5),
    "I":  dict(A=12000.0, Z=53, M=126.90, edge=33170.0, jump=1.0),
    "Cs": dict(A=11000.0, Z=55, M=132.91, edge=35985.0, jump=1.0),
}

_WARNED = {"approx": False}


def _warn_approx() -> None:
    if not _WARNED["approx"]:
        warnings.warn(
            "helios_postprocess.xray.response is using the ANCHORED POWER-LAW "
            "photoabsorption model (~25% accurate, worse near edges). Load "
            "CXRO/Henke tables before any absolute comparison to measured "
            "data. See module docstring.",
            RuntimeWarning, stacklevel=3)
        _WARNED["approx"] = True


def mu_over_rho_element(element: str, E_eV: np.ndarray) -> np.ndarray:
    """Approximate cold mass attenuation coefficient [cm^2/g] for an element."""
    if element not in _ELEMENTS:
        raise KeyError(f"No photoabsorption data for element '{element}'. "
                       f"Known: {sorted(_ELEMENTS)}")
    _warn_approx()
    p = _ELEMENTS[element]
    E = np.maximum(np.asarray(E_eV, float), 1.0)
    mu = p["A"] * (E / 1500.0) ** _EXPONENT
    if p["edge"] > 0.0 and p["jump"] > 1.0:
        mu = np.where(E < p["edge"], mu / p["jump"], mu)
    return mu


def mu_over_rho_compound(composition: Dict[str, float],
                         E_eV: np.ndarray) -> np.ndarray:
    """
    Mass attenuation for a compound given as {element: atom_count}.

    Example: polyimide/Kapton C22H10N2O5 -> dict(C=22, H=10, N=2, O=5)
    """
    E_eV = np.asarray(E_eV, float)
    total_mass = sum(n * _ELEMENTS[el]["M"] for el, n in composition.items())
    mu = np.zeros_like(E_eV)
    for el, n in composition.items():
        w = n * _ELEMENTS[el]["M"] / total_mass
        mu += w * mu_over_rho_element(el, E_eV)
    return mu


# ---------------------------------------------------------------------------
# Henke table loader (production path)
# ---------------------------------------------------------------------------
@dataclass
class HenkeTable:
    """Tabulated mu/rho vs energy, log-log interpolated."""
    name: str
    E_eV: np.ndarray
    mu_over_rho: np.ndarray

    def __call__(self, E_eV: np.ndarray) -> np.ndarray:
        E = np.clip(np.asarray(E_eV, float), self.E_eV[0], self.E_eV[-1])
        return np.exp(np.interp(np.log(E), np.log(self.E_eV),
                                np.log(np.maximum(self.mu_over_rho, 1e-30))))


def load_henke_file(path, name: Optional[str] = None) -> HenkeTable:
    """
    Load a two-column CXRO/Henke attenuation file.

    Accepts files with '#' or text header lines; the first two whitespace- or
    comma-separated numeric columns are taken as (energy_eV, mu/rho cm^2/g).
    """
    path = Path(path)
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip().replace(",", " ")
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except (ValueError, IndexError):
            continue
    if not rows:
        raise ValueError(f"No numeric rows parsed from {path}")
    arr = np.array(rows, dtype=float)
    order = np.argsort(arr[:, 0])
    return HenkeTable(name=name or path.stem,
                      E_eV=arr[order, 0], mu_over_rho=arr[order, 1])


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------
@dataclass
class Filter:
    """
    A cold absorbing foil.

    Parameters
    ----------
    name : str
    thickness_um : float
    density : float
        g/cm^3
    composition : dict, optional
        {element: atom_count}.  Ignored if ``table`` is given.
    table : HenkeTable, optional
        Tabulated mu/rho -- takes precedence over the analytic model.
    """
    name: str
    thickness_um: float
    density: float
    composition: Optional[Dict[str, float]] = None
    table: Optional[HenkeTable] = None

    @property
    def areal_density(self) -> float:
        """g/cm^2"""
        return self.density * self.thickness_um * 1.0e-4

    def transmission(self, E_eV: np.ndarray) -> np.ndarray:
        if self.table is not None:
            mu = self.table(E_eV)
        elif self.composition is not None:
            mu = mu_over_rho_compound(self.composition, E_eV)
        else:
            raise ValueError(f"Filter '{self.name}' has neither composition "
                             f"nor a Henke table.")
        return np.exp(-np.clip(mu * self.areal_density, 0.0, 7.0e2))


# Convenience constructors for the foils that actually show up on an XRSC
def be_filter(thickness_um: float) -> Filter:
    return Filter("Be", thickness_um, 1.85, composition={"Be": 1})


def al_filter(thickness_um: float) -> Filter:
    return Filter("Al", thickness_um, 2.70, composition={"Al": 1})


def ch_filter(thickness_um: float) -> Filter:
    return Filter("CH", thickness_um, 1.05, composition={"C": 1, "H": 1})


def kapton_filter(thickness_um: float) -> Filter:
    return Filter("Kapton", thickness_um, 1.42,
                  composition={"C": 22, "H": 10, "N": 2, "O": 5})


def ti_filter(thickness_um: float) -> Filter:
    return Filter("Ti", thickness_um, 4.51, composition={"Ti": 1})


# ---------------------------------------------------------------------------
# Photocathode
# ---------------------------------------------------------------------------
@dataclass
class Photocathode:
    """
    Transmission-mode photocathode quantum efficiency, crude model.

    QE(E) ~ QE_peak * (1 - exp(-mu_rho * rho * t)) * exp(-E / E_roll)

    The first factor is x-ray absorption in the cathode layer; the second is a
    stand-in for the falling escape probability of secondary electrons at high
    photon energy.  Parameters are order-of-magnitude only -- replace with a
    measured curve for quantitative work.
    """
    name: str = "CsI"
    thickness_um: float = 0.2
    density: float = 4.51
    composition: Dict[str, float] = field(default_factory=lambda: {"Cs": 1, "I": 1})
    qe_peak: float = 0.15
    e_rolloff_eV: float = 12000.0

    def qe(self, E_eV: np.ndarray) -> np.ndarray:
        mu = mu_over_rho_compound(self.composition, E_eV)
        absorbed = 1.0 - np.exp(-mu * self.density * self.thickness_um * 1e-4)
        return self.qe_peak * absorbed * np.exp(-np.asarray(E_eV, float)
                                                / self.e_rolloff_eV)


# ---------------------------------------------------------------------------
# Composite channel
# ---------------------------------------------------------------------------
@dataclass
class Channel:
    """
    A single detection channel: a filter stack, an optional photocathode, and
    an optional hard energy window.

    ``response(E)`` returns a dimensionless efficiency to be multiplied into
    the spectral intensity before integrating over photon energy.
    """
    name: str
    filters: Sequence[Filter] = ()
    cathode: Optional[Photocathode] = None
    E_min_eV: float = 0.0
    E_max_eV: float = np.inf
    mirror_reflectivity: Optional[float] = None

    def response(self, E_eV: np.ndarray) -> np.ndarray:
        E_eV = np.asarray(E_eV, float)
        R = np.ones_like(E_eV)
        for f in self.filters:
            R = R * f.transmission(E_eV)
        if self.cathode is not None:
            R = R * self.cathode.qe(E_eV)
        if self.mirror_reflectivity is not None:
            R = R * float(self.mirror_reflectivity)
        R = np.where((E_eV >= self.E_min_eV) & (E_eV <= self.E_max_eV), R, 0.0)
        return R

    def mean_energy(self, E_eV: np.ndarray, spectrum: np.ndarray) -> float:
        """Response-weighted mean photon energy of a given spectrum, eV."""
        w = self.response(E_eV) * np.asarray(spectrum, float)
        denom = np.trapezoid(w, E_eV)
        if denom <= 0:
            return float("nan")
        return float(np.trapezoid(w * E_eV, E_eV) / denom)


_DEFAULT_CATHODE = object()


def channels_from_edges(bands, cathode=_DEFAULT_CATHODE,
                        filters: Sequence[Filter] = ()) -> Dict[str, Channel]:
    """
    Build channels from explicit (E_lo, E_hi) pairs in eV.

    These are IDEALIZED rectangular bands, optionally multiplied by a common
    filter stack and photocathode. They are the right tool for scoping which
    part of the spectrum carries information; they are NOT a model of a real
    filtered channel, whose band is defined by the filter's own transmission
    edge and whose response has long tails on both sides. Move to explicit
    ``Filter`` stacks before designing hardware.

    Example
    -------
    >>> chs = channels_from_edges([(2000, 4000), (4000, 8000), (8000, 20000)])
    """
    # Pass cathode=None explicitly for a pure band study with no detector
    # weighting -- useful above ~10 keV, where the CsI QE model rolls off
    # hard and would otherwise hide the signal you are trying to scope.
    if cathode is _DEFAULT_CATHODE:
        cathode = Photocathode()
    out: Dict[str, Channel] = {}
    for lo, hi in bands:
        name = f"{lo / 1000:.3g}-{hi / 1000:.3g}keV"
        out[name] = Channel(name, filters=tuple(filters), cathode=cathode,
                            E_min_eV=float(lo), E_max_eV=float(hi))
    return out


def default_channels() -> Dict[str, Channel]:
    """
    A conventional three-channel filtered XRSC set for direct-drive
    self-emission, roughly 1-2 keV / 2-4 keV / >4 keV.

    Soft channel sees the ablation front and coronal continuum; hard channel
    is dominated by the stagnation flash.  The ratio of the two is the
    poor-man's Te diagnostic.
    """
    cathode = Photocathode()
    return {
        "soft": Channel("soft", filters=[be_filter(12.0)], cathode=cathode,
                        E_min_eV=800.0, E_max_eV=2000.0),
        "mid": Channel("mid", filters=[be_filter(50.0)], cathode=cathode,
                       E_min_eV=2000.0, E_max_eV=4000.0),
        "hard": Channel("hard", filters=[be_filter(50.0), al_filter(25.0)],
                        cathode=cathode, E_min_eV=4000.0),
    }
