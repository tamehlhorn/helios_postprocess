"""
Per-region target composition from a Helios RHW, for the multi-material
neutron-scatter source.
=========================================================================

The neutron down-scatter a detector sees is not fuel-only: neutrons scatter off
the ablator / foam carbon as well, and the fuel D:T ratio and the foam wetting
vary from layer to layer. Those compositions are *in the RHW* -- each
``Parameters for Region`` block carries the initial density, the mean atomic
weight, the radial extent, and the fusion-reactant mass fractions
(``Fus reac mass frac D/T/He3/B11/H``). The non-fusion remainder is the
structural material -- carbon for CD foam / CH / HDC ablators.

This module reads those blocks so the scatter model can build **per-species**
areal densities (rhoR_D, rhoR_T, rhoR_C) instead of a single hand-set fuel
fraction. Example, from a Vulcan wetted-foam target:

    DT gas       : w_D 0.500, w_T 0.500, w_C 0.000   (rho 0.0006, EOS DT.prp)
    DT/CD foam   : w_D 0.417, w_T 0.417, w_C 0.166   (rho 0.30,   EOS DT_0p25_and_CD_0p05.prp)

Author: Prof T (helios_postprocess)
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

logger = logging.getLogger(__name__)

# atomic masses (u) for converting mass fractions -> number densities
M_U = {"D": 2.014102, "T": 3.016049, "C": 12.011, "H": 1.007825,
       "He3": 3.016029, "B11": 11.009305}

_FUS_KEYS = {                     # RHW "Fus reac mass frac X" -> species
    "Fus reac mass frac D": "D",
    "Fus reac mass frac T": "T",
    "Fus reac mass frac He3": "He3",
    "Fus reac mass frac B11": "B11",
    "Fus reac mass frac H": "H",
}


@dataclass
class RegionComposition:
    """Initial composition of one RHW spatial region."""
    name: str
    index: int
    r_min_cm: float
    r_max_cm: float
    density_gcc: float
    mean_A: float
    mass_frac: Dict[str, float]        # species -> initial mass fraction
    eos_file: str = ""

    @property
    def w_D(self) -> float:
        return self.mass_frac.get("D", 0.0)

    @property
    def w_T(self) -> float:
        return self.mass_frac.get("T", 0.0)

    @property
    def w_C(self) -> float:
        """Carbon (structural remainder) mass fraction: 1 - sum(fusion fracs).

        For CD foam / CH / HDC the non-fusion mass is carbon; trace dopants are
        folded into carbon (negligible for neutron elastic scatter)."""
        fus = sum(v for v in self.mass_frac.values())
        return max(0.0, 1.0 - fus)

    def number_fractions(self) -> Dict[str, float]:
        """Atom-number fractions of D, T, C (normalised)."""
        w = {"D": self.w_D, "T": self.w_T, "C": self.w_C}
        n = {s: (w[s] / M_U["C" if s == "C" else s]) for s in w}
        tot = sum(n.values())
        return {s: (n[s] / tot if tot > 0 else 0.0) for s in n}


def _parse_kv(line: str):
    """Return (key, value) for a ``key = value`` RHW line, else (None, None)."""
    if "=" not in line:
        return None, None
    k, v = line.split("=", 1)
    return k.strip(), v.strip()


def parse_rhw_regions(rhw_path) -> List[RegionComposition]:
    """Parse the ``[Spatial Grid Data]`` region blocks of a Helios RHW.

    Returns a list of :class:`RegionComposition` in file order. Robust to the
    legacy text RHW format (the block layout shown in Xcimer Vulcan/HDD decks).
    Regions with no recognised composition are still returned (all-zero fracs)."""
    text = Path(rhw_path).read_text(errors="replace")
    regions: List[RegionComposition] = []

    # split on region headers, keep the name
    parts = re.split(r"Parameters for Region\s*=\s*", text)
    for idx, block in enumerate(parts[1:]):          # parts[0] is the preamble
        lines = block.splitlines()
        name = lines[0].strip()
        kv: Dict[str, str] = {}
        for ln in lines[1:]:
            # stop if we run into the next major section header
            if ln.strip().startswith("[") and "Region" not in ln:
                pass
            k, v = _parse_kv(ln)
            if k is not None and k not in kv:
                kv[k] = v                            # first occurrence wins

        def _f(key, default=0.0):
            try:
                return float(kv.get(key, default))
            except (TypeError, ValueError):
                return default

        mass_frac = {}
        for rhw_key, sp in _FUS_KEYS.items():
            val = _f(rhw_key, 0.0)
            if val > 0:
                mass_frac[sp] = val

        regions.append(RegionComposition(
            name=name, index=idx,
            r_min_cm=_f("Min. radius"), r_max_cm=_f("Max. radius"),
            density_gcc=_f("Density"), mean_A=_f("Mean atomic weight"),
            mass_frac=mass_frac, eos_file=kv.get("EOS filepath", ""),
        ))

    logger.info("Parsed %d RHW regions from %s", len(regions), rhw_path)
    return regions


def summarize_regions(regions: List[RegionComposition]) -> str:
    """One-line-per-region composition summary."""
    out = ["region                    r[cm]           rho     w_D    w_T    w_C"]
    for r in regions:
        out.append(f"{r.name:24s} {r.r_min_cm:.4f}-{r.r_max_cm:.4f} "
                   f"{r.density_gcc:8.4g}  {r.w_D:.3f}  {r.w_T:.3f}  {r.w_C:.3f}")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Run bridge: zones -> regions -> per-species areal densities
# ---------------------------------------------------------------------------

def assign_zones_to_regions(zone_boundaries: np.ndarray,
                            regions: List[RegionComposition]) -> np.ndarray:
    """Map each Lagrangian zone to an RHW region index by its *initial* radius.

    Composition is fixed per Lagrangian zone (Helios: material moves with the
    zone), so the t=0 geometry is the right key. Returns an ``(n_zones,)`` int
    array of region indices, clamped to the valid range."""
    zb0 = np.asarray(zone_boundaries, dtype=float)[0]
    centers = 0.5 * (zb0[:-1] + zb0[1:])
    idx = np.full(centers.shape, -1, dtype=int)
    for j, r in enumerate(regions):
        m = (centers >= r.r_min_cm) & (centers < r.r_max_cm)
        idx[m] = j
    # clamp: anything inside the innermost / outside the outermost region
    idx[centers < regions[0].r_min_cm] = 0
    idx[centers >= regions[-1].r_max_cm] = len(regions) - 1
    idx[idx < 0] = 0
    return idx


def zone_mass_fractions(zone_region_idx: np.ndarray,
                        regions: List[RegionComposition]):
    """Per-zone (w_D, w_T, w_C) mass fractions from each zone's region."""
    wD = np.array([regions[k].w_D for k in zone_region_idx])
    wT = np.array([regions[k].w_T for k in zone_region_idx])
    wC = np.array([regions[k].w_C for k in zone_region_idx])
    return wD, wT, wC


def _bang_index(data) -> int:
    """Time index of peak DT production (bang time)."""
    rate = getattr(data, "fusion_power", None)
    if rate is None:
        rate = getattr(data, "fusion_rate_DT_nHe4", None)
    if rate is None:
        return int(np.asarray(getattr(data, "time")).size) - 1
    r = np.asarray(rate, dtype=float)
    zm = getattr(data, "zone_mass", None)
    series = (r * np.asarray(zm, dtype=float)).sum(axis=1) if zm is not None \
        else r.sum(axis=1)
    return int(np.argmax(series))


def species_areal_densities(data, regions: List[RegionComposition],
                            rhoR_fuel_gcm2: float,
                            bang_index: Optional[int] = None) -> Optional[Dict]:
    """Split a fuel areal density into per-species D/T/C areal densities.

    The composition *ratios* (D:T:C) are taken from the bang-time radial mass
    integral ``sum(rho * w_species * dr)``; the fuel total (D+T) is anchored to
    the caller's neutron-weighted ``rhoR_fuel_gcm2`` (the value the DT-only path
    already reports), so the fuel DSR is unchanged and carbon is added
    consistently. Returns a dict of areal densities (g/cm^2) and the mass split,
    or ``None`` if geometry/composition is unavailable."""
    zb = getattr(data, "zone_boundaries", None)
    rho = getattr(data, "mass_density", None)
    if zb is None or rho is None or not regions:
        return None
    zb = np.asarray(zb, dtype=float)
    rho = np.asarray(rho, dtype=float)
    t = _bang_index(data) if bang_index is None else int(bang_index)
    t = max(0, min(t, rho.shape[0] - 1))

    zidx = assign_zones_to_regions(zb, regions)
    wD, wT, wC = zone_mass_fractions(zidx, regions)
    dr = np.diff(zb[t])                              # cm, per zone (compressed geom)
    col = rho[t] * dr                                # g/cm^2 mass column per zone
    I_D = float(np.sum(col * wD)); I_T = float(np.sum(col * wT))
    I_C = float(np.sum(col * wC))
    fuel = I_D + I_T
    if fuel <= 0:
        return None
    sD, sT, sC = I_D / fuel, I_T / fuel, I_C / fuel
    return {
        "rhoR_D_gcm2": rhoR_fuel_gcm2 * sD,
        "rhoR_T_gcm2": rhoR_fuel_gcm2 * sT,
        "rhoR_C_gcm2": rhoR_fuel_gcm2 * sC,
        "rhoR_fuel_gcm2": rhoR_fuel_gcm2,
        "mass_split": {"D": sD, "T": sT, "C": sC},
        "bang_index": t,
        "zone_region_idx": zidx,
    }


def fuel_atom_fractions(regions: List[RegionComposition]) -> Dict[str, float]:
    """Whole-target fuel D:T atom fractions (frac_D, frac_T) from the RHW,
    mass-weighted across all DT-bearing regions -- so the scatter model needs no
    hand-set fractions. Defaults to equimolar if no fuel is found."""
    wD = sum(r.w_D * r.density_gcc for r in regions)
    wT = sum(r.w_T * r.density_gcc for r in regions)
    nD = wD / M_U["D"]; nT = wT / M_U["T"]
    if (nD + nT) <= 0:
        return {"frac_D": 0.5, "frac_T": 0.5}
    fD = nD / (nD + nT)
    return {"frac_D": fD, "frac_T": 1.0 - fD}
