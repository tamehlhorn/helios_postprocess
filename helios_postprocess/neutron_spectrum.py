"""
Synthetic Neutron Spectrum + neutronics extraction for Helios ICF
=================================================================

Two layers, aligned with the Xcimer target-group workflow:

1. **Extraction (matches K. Keipper's ``neutronics_output.py``)** — pull the
   transport-ready source data out of a Helios run: per-channel volumetric
   fusion rates, per-zone neutron production, and neutron-weighted
   time-averaged radial profiles ``<rho(r)>`` / ``<T_ion(r)>`` (hydro weighted
   by *when* neutrons are born, source rates by *where + when*), saved to a
   ``neutronics_data.npz`` with the same schema Kyle's tool writes. This is the
   hand-off artifact a neutron-transport / synthetic-diagnostic code (NeSST,
   IRIS) consumes. Built on ``rhino.HeliosSphericalSimulation`` for material
   interfaces / bang time / units when RHINO is importable, with a self-
   contained fallback for machines that lack it.

2. **In-house quick-look observables** — an analytic synthetic *birth* spectrum
   (Brysk-broadened primary peaks) giving the nTOF (spectral) ion temperature,
   plus a forward-model down-scatter ratio. This is a fast surrogate for a
   single-scatter tool (NeSST) or a full-transport code (IRIS); it is labelled
   as such. The DSR here is a calibration forward model (DSR = rhoR / coeff),
   **not** a scattered-transport result — treat it as a placeholder until a
   NeSST/IRIS cross-check is available.

Conventions (matched to Kyle / Helios EXODUS)
---------------------------------------------
- DT neutron energy **14.06 MeV** (Helios ``TimeIntFusionProd_n_1406``),
  DD-n **2.45 MeV**.
- ``FusionRate_*`` is reactions/s **per gram** (CLAUDE.md §5b). Volumetric rate
  (reactions/cm^3/s) = ``FusionRate * rho``; per-zone reactions/s =
  ``FusionRate * rho * cell_volume`` = ``fusion_power * zone_mass`` (identical).
- Six channels carried: DT->nHe4, DD->nHe3, DD->pT, TT->2nHe4, DHe3->pHe4,
  pB11->3He4. Neutron-producing: DT, DD->nHe3, TT.

References: Brysk, Plasma Phys. 15, 611 (1973); Ballabio et al., Nucl. Fusion
38, 1723 (1998); Crilly et al., NeSST (arXiv:2605.20432); Murphy et al., RSI 85,
11D901 (2014).

Author: Prof T (helios_postprocess)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.integrate import trapezoid

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants / calibration
# ---------------------------------------------------------------------------

#: Mean neutron birth energy per reaction (MeV). DT matches Helios' 14.06 label.
E_BIRTH_MEV = {"DT": 14.06, "DD": 2.45}

#: Brysk thermal-broadening coefficient per reaction: FWHM[keV] = C * sqrt(Ti[keV]).
#: DT ~177; DD narrower (~ sqrt(E0_DD/E0_DT) x DT ~ 82.5, Ballabio).
BRYSK_FWHM = {"DT": 177.0, "DD": 82.5}

#: FWHM = FWHM_SIGMA * sigma for a Gaussian.
FWHM_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))   # ~2.35482

#: DSR -> rhoR calibration coefficients (g/cm^2 per unit DSR).
DSR_CALIB = {"NIF": 20.4, "OMEGA": 19.4}

#: EXODUS variable names for the six Helios fusion channels (reactions/s/g).
FUSION_RATE_VARS = {
    "DT_nHe4":  "FusionRate_DT_nHe4",
    "DD_nHe3":  "FusionRate_DD_nHe3",
    "DD_pT":    "FusionRate_DD_pT",
    "TT_nnHe4": "FusionRate_TT_nnHe4",
    "DHe3_pHe4": "FusionRate_DHe3_pHe4",
    "pB11_3He4": "FusionRate_pB11_3He4",
}

#: ICFRunData attribute names for the same channels (loaded by data_builder).
FUSION_RATE_ATTRS = {
    "DT_nHe4":  "fusion_power",
    "DD_nHe3":  "fusion_rate_DD_nHe3",
    "DD_pT":    "fusion_rate_DD_pT",
    "TT_nnHe4": "fusion_rate_TT_nnHe4",
    "DHe3_pHe4": "fusion_rate_DHe3_pHe4",
    "pB11_3He4": "fusion_rate_pB11_3He4",
}


def brysk_fwhm_keV(Ti_keV: np.ndarray, reaction: str = "DT") -> np.ndarray:
    """Primary-peak FWHM (keV) for ion temperature ``Ti_keV`` and ``reaction``."""
    Ti_keV = np.asarray(Ti_keV, dtype=float)
    return BRYSK_FWHM[reaction] * np.sqrt(np.clip(Ti_keV, 0.0, None))


def fwhm_keV_to_Ti_keV(fwhm_keV: float, reaction: str = "DT") -> float:
    """Invert the Brysk relation: spectral (nTOF) ion temperature from FWHM."""
    return (float(fwhm_keV) / BRYSK_FWHM[reaction]) ** 2


# ---------------------------------------------------------------------------
# Geometry / rate helpers (Kyle-convention)
# ---------------------------------------------------------------------------

def spherical_com_and_volume(r_boundaries_cm: np.ndarray
                             ) -> Tuple[np.ndarray, np.ndarray]:
    """Volume-weighted shell centroid radius and cell volume from zone
    boundaries ``(n_t, n_zone+1)``. Matches ``neutronics_output.py``."""
    r_boundaries_cm = np.asarray(r_boundaries_cm, dtype=float)
    r_lo = r_boundaries_cm[:, :-1]
    r_hi = r_boundaries_cm[:, 1:]
    denom = r_hi ** 2 + r_hi * r_lo + r_lo ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        r_com = 0.75 * (r_hi ** 3 + r_hi ** 2 * r_lo
                        + r_hi * r_lo ** 2 + r_lo ** 3) / denom
    r_com = np.where(denom > 0, r_com, 0.5 * (r_lo + r_hi))
    cell_volume = 4.0 / 3.0 * np.pi * (r_hi ** 3 - r_lo ** 3)
    return r_com, cell_volume


def volumetric_rate(fusion_rate_per_g_s: np.ndarray, rho_gcc: np.ndarray) -> np.ndarray:
    """reactions/cm^3/s = (reactions/s/g) * (g/cm^3). Kyle's convention."""
    return np.asarray(fusion_rate_per_g_s, dtype=float) * np.asarray(rho_gcc, dtype=float)


def emission_weights(fusion_power: np.ndarray, zone_mass: np.ndarray) -> np.ndarray:
    """Total reactions/s per zone = fusion_power * zone_mass (== vol_rate * cell_vol)."""
    w = np.asarray(fusion_power, dtype=float) * np.asarray(zone_mass, dtype=float)
    return np.clip(w, 0.0, None)


# ---------------------------------------------------------------------------
# Burn history
# ---------------------------------------------------------------------------

def burn_history(time_ns: np.ndarray, weights: np.ndarray,
                 ion_temperature_eV: np.ndarray) -> Dict:
    """Time-resolved burn diagnostics from per-zone emission ``weights``
    (reactions/s per zone). See module tests for the analytic checks."""
    time_ns = np.asarray(time_ns, dtype=float)
    weights = np.asarray(weights, dtype=float)
    Ti_keV = np.asarray(ion_temperature_eV, dtype=float) / 1000.0

    rate = weights.sum(axis=1)
    time_s = time_ns * 1.0e-9
    total = float(trapezoid(rate, time_s)) if rate.size > 1 else 0.0

    if rate.max() > 0:
        bang_time = float(time_ns[int(np.argmax(rate))])
        burn_fwhm = _fwhm_of_curve(time_ns, rate)
    else:
        bang_time = float("nan")
        burn_fwhm = float("nan")

    with np.errstate(invalid="ignore", divide="ignore"):
        wsum_t = weights.sum(axis=1)
        Ti_emission = np.where(wsum_t > 0, (weights * Ti_keV).sum(axis=1) / wsum_t, 0.0)
        total_w = weights.sum()
        Ti_burn_avg = float((weights * Ti_keV).sum() / total_w) if total_w > 0 else float("nan")
        rhoR_weight = rate / rate.sum() if rate.sum() > 0 else np.zeros_like(rate)

    return {"rate_per_s": rate, "total_reactions": total, "bang_time_ns": bang_time,
            "burn_fwhm_ns": burn_fwhm, "Ti_emission_keV": Ti_emission,
            "Ti_burn_avg_keV": Ti_burn_avg, "rhoR_weight": rhoR_weight}


def _fwhm_of_curve(x: np.ndarray, y: np.ndarray) -> float:
    y = np.asarray(y, dtype=float); x = np.asarray(x, dtype=float)
    ymax = y.max()
    if ymax <= 0:
        return float("nan")
    half = 0.5 * ymax
    above = y >= half
    if not above.any():
        return float("nan")
    idx = np.where(above)[0]; i0, i1 = idx[0], idx[-1]
    x_lo = np.interp(half, [y[i0 - 1], y[i0]], [x[i0 - 1], x[i0]]) if i0 > 0 else x[i0]
    x_hi = np.interp(half, [y[i1 + 1], y[i1]], [x[i1 + 1], x[i1]]) if i1 < len(x) - 1 else x[i1]
    return float(abs(x_hi - x_lo))


# ---------------------------------------------------------------------------
# Neutron-weighted time-averaged profiles (Kyle's method)
# ---------------------------------------------------------------------------

def neutron_weighted_profiles(
    r_com_cm: np.ndarray,
    cell_volume_cm3: np.ndarray,
    dt_s: np.ndarray,
    rho_gcc: np.ndarray,
    T_ion_eV: np.ndarray,
    channel_rates_per_g_s: Dict[str, np.ndarray],
    interface_positions_cm: Optional[list] = None,
    n_r: int = 300,
    r_max_cm: Optional[float] = None,
) -> Dict:
    """Neutron-weighted, time-averaged Eulerian profiles per neutron channel.

    Ports the algorithm in ``neutronics_output.py``: a 14 MeV neutron crosses
    the target (~0.3 cm) in ~10 ps, effectively instantaneous vs the ~0.1 ns
    burn, so it sees the *whole* radial profile at birth. Hydro profiles
    (``rho``, ``T_ion``) are therefore weighted by **when** neutrons are born
    (total channel production per timestep, uniform in radius); source rates are
    weighted by **where + when** (spatially resolved).

    Returns ``{channel: {rho_avg_gcc, T_ion_avg_eV, interface_avg_cm,
    rates_avg}}`` for each neutron-producing channel with nonzero yield.
    """
    r_com_cm = np.asarray(r_com_cm, dtype=float)
    cell_volume_cm3 = np.asarray(cell_volume_cm3, dtype=float)
    dt_s = np.asarray(dt_s, dtype=float)
    rho_gcc = np.asarray(rho_gcc, dtype=float)
    T_ion_eV = np.asarray(T_ion_eV, dtype=float)
    ntimes = rho_gcc.shape[0]
    interface_positions_cm = interface_positions_cm or []

    if r_max_cm is None:
        r_max_cm = float(np.nanmax(r_com_cm[0]))
    r_avg_cm = np.linspace(0.0, r_max_cm, n_r)

    vol_rates = {k: volumetric_rate(v, rho_gcc) for k, v in channel_rates_per_g_s.items()}
    # Neutron-producing channels and their per-zone-per-step production weights.
    neutron_channels = {
        name: vol_rates[name] * cell_volume_cm3 * dt_s[:, None]
        for name in ("DT_nHe4", "DD_nHe3", "TT_nnHe4") if name in vol_rates
    }

    avg_results: Dict[str, Optional[dict]] = {}
    for ch_name, ch_w in neutron_channels.items():
        temporal_w = ch_w.sum(axis=1)
        total_tw = temporal_w.sum()
        if total_tw <= 0:
            avg_results[ch_name] = None
            continue

        rho_acc = np.zeros(n_r); T_acc = np.zeros(n_r); tw_sum = 0.0
        iface_acc = [0.0] * len(interface_positions_cm)
        rate_accs = {k: np.zeros(n_r) for k in vol_rates}
        sw_acc = np.zeros(n_r)

        for t in range(ntimes):
            tw = temporal_w[t]
            if tw < 1e-20 * total_tw:
                continue
            r_lag = r_com_cm[t]
            rho_acc += np.interp(r_avg_cm, r_lag, rho_gcc[t], left=0, right=0) * tw
            T_acc += np.interp(r_avg_cm, r_lag, T_ion_eV[t], left=0, right=0) * tw
            tw_sum += tw
            for k, iarr in enumerate(interface_positions_cm):
                iface_acc[k] += float(iarr[t]) * tw
            sw = np.interp(r_avg_cm, r_lag, ch_w[t], left=0, right=0)
            for rname, rarr in vol_rates.items():
                rate_accs[rname] += np.interp(r_avg_cm, r_lag, rarr[t], left=0, right=0) * sw
            sw_acc += sw

        if tw_sum > 0:
            rho_acc /= tw_sum; T_acc /= tw_sum
            iface_acc = [v / tw_sum for v in iface_acc]
        mask = sw_acc > 0
        for rname in rate_accs:
            rate_accs[rname][mask] /= sw_acc[mask]

        avg_results[ch_name] = {"rho_avg_gcc": rho_acc, "T_ion_avg_eV": T_acc,
                                "interface_avg_cm": iface_acc, "rates_avg": rate_accs}

    return {"r_avg_cm": r_avg_cm, "avg_results": avg_results}


# ---------------------------------------------------------------------------
# Synthetic birth spectrum + nTOF Tion  (in-house quick-look)
# ---------------------------------------------------------------------------

def synthesize_birth_spectrum(
    weights: np.ndarray, ion_temperature_eV: np.ndarray,
    energy_grid_MeV: Optional[np.ndarray] = None,
    reaction: str = "DT", weight_floor_frac: float = 1e-6,
) -> Tuple[np.ndarray, np.ndarray]:
    """Burn-weighted primary birth spectrum N(E) (counts/MeV) for ``reaction``.
    Each cell contributes a Brysk-broadened Gaussian at ``E_BIRTH_MEV[reaction]``
    weighted by its reaction count; integrates to ``weights.sum()``."""
    weights = np.asarray(weights, dtype=float)
    Ti_keV = np.asarray(ion_temperature_eV, dtype=float) / 1000.0
    e0 = E_BIRTH_MEV[reaction]

    if energy_grid_MeV is None:
        half = 4.0 if reaction == "DT" else 1.0
        energy_grid_MeV = np.linspace(e0 - half, e0 + half, 1500)
    energy_grid_MeV = np.asarray(energy_grid_MeV, dtype=float)

    spectrum = np.zeros_like(energy_grid_MeV)
    wmax = weights.max()
    if wmax <= 0:
        return energy_grid_MeV, spectrum
    floor = weight_floor_frac * wmax
    sigma_MeV_all = (brysk_fwhm_keV(Ti_keV, reaction) / FWHM_SIGMA) / 1000.0

    for t in range(weights.shape[0]):
        sel = weights[t] > floor
        if not sel.any():
            continue
        w = weights[t, sel]
        sig = np.clip(sigma_MeV_all[t, sel], 1e-4, None)
        dE = energy_grid_MeV[None, :] - e0
        g = (w[:, None] / (sig[:, None] * np.sqrt(2 * np.pi))
             * np.exp(-0.5 * (dE / sig[:, None]) ** 2))
        spectrum += g.sum(axis=0)
    return energy_grid_MeV, spectrum


def infer_ntof_tion(energy_MeV: np.ndarray, spectrum: np.ndarray,
                    reaction: str = "DT",
                    fit_window_MeV: Optional[Tuple[float, float]] = None) -> Dict:
    """Spectral (nTOF-observable) ion temperature from the primary-peak width."""
    energy_MeV = np.asarray(energy_MeV, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    e0 = E_BIRTH_MEV[reaction]
    if fit_window_MeV is None:
        w = 0.8 if reaction == "DT" else 0.3
        fit_window_MeV = (e0 - w, e0 + w)
    lo, hi = fit_window_MeV
    mask = (energy_MeV >= lo) & (energy_MeV <= hi)
    if not mask.any() or spectrum[mask].sum() <= 0:
        return {"Ti_ntof_keV": float("nan"), "fwhm_keV": float("nan"),
                "peak_energy_MeV": float("nan")}
    E = energy_MeV[mask]; S = np.clip(spectrum[mask], 0.0, None)
    norm = trapezoid(S, E)
    mean_E = trapezoid(E * S, E) / norm
    var_E = trapezoid((E - mean_E) ** 2 * S, E) / norm
    fwhm_keV = np.sqrt(max(var_E, 0.0)) * FWHM_SIGMA * 1000.0
    return {"Ti_ntof_keV": fwhm_keV_to_Ti_keV(fwhm_keV, reaction),
            "fwhm_keV": float(fwhm_keV), "peak_energy_MeV": float(mean_E)}


def synthetic_dsr(rhoR_gcm2: float, calibration: str = "NIF") -> Dict:
    """Forward-model DSR = rhoR / coeff (placeholder for a NeSST/IRIS scattered
    spectrum). NOT a transport result: it merely relabels rhoR in detector
    units, so it carries no information beyond the rhoR already known from the
    hydro. Use only until a single-scatter (NeSST) or MC (IRIS) DSR is wired."""
    if calibration not in DSR_CALIB:
        raise ValueError(f"Unknown calibration {calibration!r}; use 'NIF' or 'OMEGA'.")
    coeff = DSR_CALIB[calibration]
    dsr = float(rhoR_gcm2) / coeff
    rhoR_rt = coeff * dsr
    try:
        from . import neutron_downscatter as nds
        res = nds.downscatter_to_areal_density(dsr, calibration=calibration)
        if isinstance(res, dict) and "rho_R" in res:
            rhoR_rt = float(res["rho_R"])
    except Exception as exc:            # pragma: no cover - defensive
        logger.debug("DSR round-trip via neutron_downscatter unavailable: %s", exc)
    return {"DSR": dsr, "coeff": coeff, "rhoR_roundtrip_gcm2": rhoR_rt,
            "calibration": calibration}


def build_downscattered_spectrum(
    energy_MeV: np.ndarray, primary_spectrum: np.ndarray, dsr: float,
    scatter_window_MeV: Tuple[float, float] = (10.0, 13.5),
    primary_window_MeV: Tuple[float, float] = (13.5, 14.5),
) -> np.ndarray:
    """Add a flat single-scatter shoulder so N(E) round-trips through
    ``neutron_downscatter.calculate_downscatter_ratio`` to ``dsr`` (plotting /
    loop-closure only; shape is not a transport result)."""
    energy_MeV = np.asarray(energy_MeV, dtype=float)
    primary_spectrum = np.asarray(primary_spectrum, dtype=float)
    plo, phi = primary_window_MeV
    pmask = (energy_MeV >= plo) & (energy_MeV <= phi)
    N_primary = trapezoid(primary_spectrum[pmask], energy_MeV[pmask]) if pmask.any() else 0.0
    lo, hi = scatter_window_MeV
    win = (energy_MeV >= lo) & (energy_MeV <= hi)
    out = primary_spectrum.copy()
    if win.any() and N_primary > 0:
        out[win] += (dsr * N_primary) / (hi - lo)
    return out


# ---------------------------------------------------------------------------
# Interfaces / bang time — RHINO backbone with fallback
# ---------------------------------------------------------------------------

def _interfaces_bang_via_rhino(sim_path) -> Optional[Tuple[list, float]]:
    """Material-interface radial histories (cm) + bang time (ns) from RHINO.
    Returns None if RHINO is unavailable or the load fails (caller falls back)."""
    try:
        import rhino as rno                                    # private: wtrickey27/RHINO
    except Exception as exc:                                   # pragma: no cover
        logger.info("RHINO not importable (%s) -- using Helios-native interfaces.", exc)
        return None
    try:
        sim = rno.HeliosSphericalSimulation(str(sim_path))
        ifaces = []
        for iface in sim.material_interfaces:
            r_vals = iface.vals * sim.radial_coords.units
            ifaces.append(r_vals.to("cm").magnitude)
        bang_ns = float(sim.bang_time.to("ns").magnitude)
        logger.info("Interfaces/bang time sourced from RHINO (%d interfaces).", len(ifaces))
        return ifaces, bang_ns
    except Exception as exc:                                   # pragma: no cover
        logger.warning("RHINO load failed (%s) -- using Helios-native interfaces.", exc)
        return None


def _interfaces_bang_native(data, r_com_cm, weights) -> Tuple[list, float]:
    """Fallback: interface radial histories from region_interfaces_indices and
    bang time from the burn-rate peak (no RHINO dependency)."""
    ifaces = []
    ri = getattr(data, "region_interfaces_indices", None)
    zb = getattr(data, "zone_boundaries", None)
    if ri is not None and zb is not None:
        zb = np.asarray(zb, dtype=float)
        ri = np.asarray(ri)
        for col in range(ri.shape[1]):
            node = ri[:, col].astype(int)
            ifaces.append(np.array([zb[t, node[t]] for t in range(zb.shape[0])]))
    rate = np.asarray(weights, dtype=float).sum(axis=1)
    t_ns = np.asarray(getattr(data, "time"), dtype=float)
    bang = float(t_ns[int(np.argmax(rate))]) if rate.max() > 0 else float("nan")
    return ifaces, bang


# ---------------------------------------------------------------------------
# Extraction container + entry point (Kyle-schema .npz)
# ---------------------------------------------------------------------------

@dataclass
class NeutronicsData:
    """Transport-ready neutronics extraction — mirrors ``neutronics_output.py``."""
    stime_ns: np.ndarray
    dt_s: np.ndarray
    r_boundaries_cm: np.ndarray
    r_com_cm: np.ndarray
    cell_volume_cm3: np.ndarray
    rho_gcc: np.ndarray
    T_ion_eV: np.ndarray
    T_elec_eV: Optional[np.ndarray]
    interface_positions_cm: list
    channel_rates_vol: Dict[str, np.ndarray]        # reactions/cm^3/s per channel
    dt_neutron_production: np.ndarray               # reactions/s per zone (DT)
    r_avg_cm: np.ndarray
    avg_results: dict
    total_dt_yield: float = 0.0
    total_dd_yield: float = 0.0
    bang_time_ns: float = 0.0
    interfaces_source: str = "native"
    # in-house quick-look observables (labelled non-transport)
    quicklook: dict = field(default_factory=dict)


def extract_neutronics(data=None, sim_path=None, use_rhino: bool = False,
                       n_r: int = 300, calibration: str = "NIF",
                       save_npz: Optional[str] = None) -> Optional[NeutronicsData]:
    """Extract transport-ready neutronics data from a Helios run.

    Parameters
    ----------
    data : ICFRunData, optional
        Populated run data (preferred). Provides arrays + per-channel rates.
    sim_path : path, optional
        Simulation directory; used for the RHINO interface/bang lookup and, if
        ``data`` is None, could be extended to read the .exo directly.
    use_rhino : bool
        Default **False**: helios_postprocess is the open-access, self-contained
        prototyping tool, so the core path never requires the private RHINO
        package. Set True to *opt in* to sourcing interfaces/bang/units from
        ``rhino.HeliosSphericalSimulation`` as a cross-check (falls back to
        Helios-native indices if RHINO is unavailable regardless).
    save_npz : path, optional
        If set, write ``neutronics_data.npz`` (Kyle schema) here.

    Returns None on a no-burn / missing-field run (graceful, like
    ``analyze_laser_intensity``).
    """
    if data is None:
        raise NotImplementedError("Direct .exo reading not yet wired; pass an ICFRunData.")

    time_ns = getattr(data, "time", None)
    zb = getattr(data, "zone_boundaries", None)
    rho = getattr(data, "mass_density", None)
    Ti = getattr(data, "ion_temperature", None)
    zone_mass = getattr(data, "zone_mass", None)
    if any(v is None for v in (time_ns, zb, rho, Ti, zone_mass)):
        logger.warning("Neutronics: missing core fields -- skipping.")
        return None

    time_ns = np.asarray(time_ns, dtype=float)
    r_com_cm, cell_volume_cm3 = spherical_com_and_volume(zb)
    dt_s = np.gradient(time_ns * 1e-9)             # EXODUS cadence approx (Kyle reads exact dt)

    # Per-channel per-gram rates from ICFRunData (missing channels -> zeros).
    shape = np.asarray(rho).shape
    channel_rates_pg = {}
    for ch, attr in FUSION_RATE_ATTRS.items():
        arr = getattr(data, attr, None)
        channel_rates_pg[ch] = np.asarray(arr, dtype=float) if arr is not None else np.zeros(shape)
    channel_rates_vol = {k: volumetric_rate(v, rho) for k, v in channel_rates_pg.items()}

    weights = emission_weights(channel_rates_pg["DT_nHe4"], zone_mass)   # reactions/s per zone
    if weights.max() <= 0:
        logger.warning("Neutronics: zero DT fusion rate (no-burn) -- skipping.")
        return None
    dt_neutron_production = weights

    # Interfaces / bang: RHINO if available, else native.
    interfaces, bang_ns, src = None, None, "native"
    if use_rhino and sim_path is not None:
        got = _interfaces_bang_via_rhino(sim_path)
        if got is not None:
            interfaces, bang_ns = got
            src = "rhino"
    if interfaces is None:
        interfaces, bang_ns = _interfaces_bang_native(data, r_com_cm, weights)

    # Yields (prefer Helios cumulative products).
    total_dt = float(np.asarray(getattr(data, "dt_neutron_count", [0.0]))[-1]) \
        if getattr(data, "dt_neutron_count", None) is not None else 0.0
    total_dd = float(np.asarray(getattr(data, "dd_neutron_count", [0.0]))[-1]) \
        if getattr(data, "dd_neutron_count", None) is not None else 0.0

    prof = neutron_weighted_profiles(
        r_com_cm, cell_volume_cm3, dt_s, np.asarray(rho, dtype=float),
        np.asarray(Ti, dtype=float), channel_rates_pg,
        interface_positions_cm=interfaces, n_r=n_r)

    # In-house quick-look: DT birth spectrum + nTOF Ti + forward-model DSR.
    energy, primary = synthesize_birth_spectrum(weights, Ti, reaction="DT")
    ntof = infer_ntof_tion(energy, primary, reaction="DT")
    rhoR_emission = _rhoR_from_profile(prof, data)
    dsr = (synthetic_dsr(rhoR_emission, calibration) if np.isfinite(rhoR_emission)
           else {"DSR": float("nan"), "rhoR_roundtrip_gcm2": float("nan"),
                 "coeff": DSR_CALIB[calibration], "calibration": calibration})
    hist = burn_history(time_ns, weights, Ti)
    quicklook = {
        "energy_MeV": energy, "birth_spectrum": primary,
        "ntof_tion_keV": ntof["Ti_ntof_keV"], "ntof_fwhm_keV": ntof["fwhm_keV"],
        "Ti_burn_avg_keV": hist["Ti_burn_avg_keV"], "burn_fwhm_ns": hist["burn_fwhm_ns"],
        "rhoR_emission_gcm2": rhoR_emission, "dsr": dsr["DSR"],
        "dsr_note": "forward-model DSR=rhoR/coeff; placeholder for NeSST/IRIS",
    }

    nd = NeutronicsData(
        stime_ns=time_ns, dt_s=dt_s, r_boundaries_cm=np.asarray(zb, dtype=float),
        r_com_cm=r_com_cm, cell_volume_cm3=cell_volume_cm3,
        rho_gcc=np.asarray(rho, dtype=float), T_ion_eV=np.asarray(Ti, dtype=float),
        T_elec_eV=(np.asarray(getattr(data, "elec_temperature"), dtype=float)
                   if getattr(data, "elec_temperature", None) is not None else None),
        interface_positions_cm=interfaces, channel_rates_vol=channel_rates_vol,
        dt_neutron_production=dt_neutron_production, r_avg_cm=prof["r_avg_cm"],
        avg_results=prof["avg_results"], total_dt_yield=total_dt, total_dd_yield=total_dd,
        bang_time_ns=(bang_ns if bang_ns is not None else hist["bang_time_ns"]),
        interfaces_source=src, quicklook=quicklook)

    if save_npz:
        save_neutronics_npz(nd, save_npz)
    logger.info("Neutronics extracted (interfaces=%s): bang %.3f ns, DT yield %.3e, "
                "nTOF Ti %.2f keV, rhoR_emission %.3f g/cm^2",
                src, nd.bang_time_ns, total_dt, ntof["Ti_ntof_keV"], rhoR_emission)
    return nd


def _rhoR_from_profile(prof: Dict, data) -> float:
    """DT-neutron-weighted fuel areal density = integral of <rho(r)>_DT over the
    fuel region, i.e. 'what a DT neutron sees'. Falls back to the scalar
    areal_density_vs_time (emission-weighted) if the profile is unavailable."""
    dt = (prof.get("avg_results") or {}).get("DT_nHe4")
    if dt is not None:
        r = prof["r_avg_cm"]; rho = dt["rho_avg_gcc"]
        # integrate rho dr from the axis out to the ablator/fuel outer interface
        r_out = None
        ivals = dt.get("interface_avg_cm") or []
        if len(ivals) >= 2:
            r_out = float(ivals[-2])           # fuel/ablator interface (inside grid edge)
        mask = (r <= r_out) if r_out else np.ones_like(r, dtype=bool)
        val = float(trapezoid(rho[mask], r[mask]))
        if np.isfinite(val) and val > 0:
            return val
    hist_rhoR = getattr(data, "areal_density_vs_time", None)
    if hist_rhoR is not None:
        return float(np.nanmax(np.asarray(hist_rhoR, dtype=float)))
    return float(getattr(data, "rhoR_cold_fuel_gcm2", float("nan")))


def save_neutronics_npz(nd: NeutronicsData, path) -> str:
    """Write a ``neutronics_data.npz`` with the same schema as Kyle's tool
    (plus our ``quicklook`` block). Returns the path written."""
    path = str(path)
    if Path(path).is_dir():
        path = str(Path(path) / "neutronics_data.npz")
    cr = nd.channel_rates_vol
    np.savez(
        path,
        stime_ns=nd.stime_ns, dt_s=nd.dt_s, r_boundaries_cm=nd.r_boundaries_cm,
        r_com_cm=nd.r_com_cm, cell_volume_cm3=nd.cell_volume_cm3,
        rho_gcc=nd.rho_gcc, T_ion_eV=nd.T_ion_eV,
        T_elec_eV=(nd.T_elec_eV if nd.T_elec_eV is not None else np.zeros_like(nd.rho_gcc)),
        rate_DT_nHe4=cr.get("DT_nHe4"), rate_DD_nHe3=cr.get("DD_nHe3"),
        rate_DD_pT=cr.get("DD_pT"), rate_TT_nnHe4=cr.get("TT_nnHe4"),
        rate_DHe3_pHe4=cr.get("DHe3_pHe4"), rate_pB11_3He4=cr.get("pB11_3He4"),
        dt_neutron_production=nd.dt_neutron_production,
        r_avg_cm=nd.r_avg_cm,
        rho_avg_gcc=(nd.avg_results.get("DT_nHe4") or {}).get("rho_avg_gcc", np.zeros_like(nd.r_avg_cm)),
        T_ion_avg_eV=(nd.avg_results.get("DT_nHe4") or {}).get("T_ion_avg_eV", np.zeros_like(nd.r_avg_cm)),
        avg_results=np.array(nd.avg_results, dtype=object),
        total_dt_yield=nd.total_dt_yield, total_dd_yield=nd.total_dd_yield,
        bang_time_ns=nd.bang_time_ns, interfaces_source=nd.interfaces_source,
        quicklook=np.array(nd.quicklook, dtype=object),
    )
    return path


# ---------------------------------------------------------------------------
# Back-compat quick-look entry point
# ---------------------------------------------------------------------------

def analyze_neutron_spectrum(data, calibration: str = "NIF") -> Optional[Dict]:
    """Quick-look wrapper (birth spectrum + nTOF Ti + forward-model DSR) that
    now routes through :func:`extract_neutronics` so it uses the Kyle-aligned
    profile-based rhoR. Returns the ``quicklook`` dict (or None on no-burn)."""
    nd = extract_neutronics(data=data, use_rhino=False, calibration=calibration)
    if nd is None:
        return None
    out = dict(nd.quicklook)
    out.update({"bang_time_ns": nd.bang_time_ns, "total_reactions": nd.total_dt_yield,
                "calibration": calibration})
    return out
