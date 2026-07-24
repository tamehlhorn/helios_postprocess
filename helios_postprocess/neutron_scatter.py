"""
Single-scattered neutron spectrum and DSR via NeSST (Step 3).
=============================================================

This is the *transport* down-scatter layer. Where
:func:`helios_postprocess.neutron_spectrum.synthetic_dsr` only relabels the
hydro areal density in detector units (DSR = rhoR / coeff, carrying no new
information), this module feeds our birth spectrum and areal density into
**NeSST** (Crilly, https://github.com/aidancrilly/NeSST) to compute the
*physically singly-scattered* neutron spectrum from ENDF/B-VIII.0 n+D and n+T
elastic (plus (n,2n)) cross sections, and the down-scatter ratio that follows
from it.

helios_postprocess owns the ICF-specific *source term* — the burn-weighted
birth spectrum ``I(E)`` and the DT-neutron-weighted areal density ``rhoR`` — and
hands them to NeSST across a thin, documented interface:

    I(E), rhoR, frac_D, frac_T   -->   NeSST   -->   scattered N(E), DSR, rhoR

NeSST convention (verified against the installed 1.1.5 source)
--------------------------------------------------------------
- Energies are in **eV**; ``E0_DT = 14.03 MeV``.
- ``init_DT_scatter(Eout, Ein)`` pre-computes the D and T scattering matrices on
  the given grids (expensive; cached here by grid signature).
- ``DT_sym_scatter_spec(I_E, frac_D, frac_T)`` returns the single-scatter
  spectrum *shape* for an isotropic areal density (its internal ``rhoL`` is
  unity); the physical magnitude is set by multiplying by the scattering
  amplitude ``A_1S = rhoR_2_A1s(rhoR_kg_m2, frac_D, frac_T)``.
- ``A1s_2_rhoR`` inverts that, giving NeSST's first-principles rhoR<->scatter
  calibration (which we compare against the empirical NIF 20.4 g/cm^2 per DSR).

For a 1D spherical Helios run the areal density is isotropic, so the symmetric
model is the default; the anisotropic path (``DT_asym_scatter_spec`` with an
angular ``rhoL(mu)``) is exposed for later 2D/3D work.

Author: Prof T (helios_postprocess)
"""

from __future__ import annotations

import logging
from typing import Callable, Dict, Optional, Tuple

import numpy as np
from scipy.integrate import trapezoid

logger = logging.getLogger(__name__)

try:
    import NeSST as _nesst
    NESST_AVAILABLE = True
except Exception as _exc:                    # pragma: no cover - environment dependent
    _nesst = None
    NESST_AVAILABLE = False
    _NESST_IMPORT_ERROR = _exc

#: 1 g/cm^2 == 10 kg/m^2 (NeSST areal-density unit).
G_CM2_TO_KG_M2 = 10.0
#: DT birth energy (MeV) for window defaults.
E0_DT_MEV = 14.06
#: Standard nTOF DSR windows (MeV): scattered 10-12, primary 13-15.
DSR_SCATTER_WINDOW_MEV = (10.0, 12.0)
DSR_PRIMARY_WINDOW_MEV = (13.0, 15.0)

# Cache of initialised scatter matrices, keyed by a grid signature, so repeated
# calls on the same energy grid do not re-load ENDF data / rebuild matrices.
_MATRIX_CACHE_KEY: Optional[Tuple] = None


def _require_nesst() -> None:
    if not NESST_AVAILABLE:
        raise ImportError(
            "NeSST is required for the transport down-scatter model "
            "(pip install NeSST). Original import error: "
            f"{_NESST_IMPORT_ERROR!r}"
        )


def _grid_signature(E_eV: np.ndarray) -> Tuple:
    return (E_eV.size, float(E_eV[0]), float(E_eV[-1]))


def _ensure_scatter_matrices(E_eV: np.ndarray, include_n2n: bool = True) -> None:
    """Initialise (or reuse cached) NeSST D/T scatter matrices on grid ``E_eV``.

    NeSST stores the matrices in module-global state, so we only re-init when the
    grid (or the ``include_n2n`` choice) changes. ``Eout == Ein == E_eV``.

    The (n,2n) differential grid is ``(N_in, N_mu, N_out)`` -- its memory scales
    as ``N_E**2`` and dominates (a 900-point grid needs ~8 GB). When
    ``include_n2n`` is False we disable those channels before the matrices are
    built (~10x less memory, finer grids possible) and seed a zero ``n2n_dNdE``
    so NeSST's spectrum functions still run. (n,2n) shifts the standard
    10-12 MeV DSR by ~1.6% relative, so the elastic-only mode is a good
    fast/low-memory approximation for the areal-density diagnostic."""
    global _MATRIX_CACHE_KEY
    _require_nesst()
    from NeSST.core import mat_dict
    sig = _grid_signature(E_eV) + (bool(include_n2n),)
    if sig == _MATRIX_CACHE_KEY:
        return
    logger.info("NeSST: initialising DT scatter matrices on %d-point grid "
                "[%.2f, %.2f] MeV (n2n=%s).", E_eV.size, E_eV[0] / 1e6,
                E_eV[-1] / 1e6, include_n2n)
    mats = [mat_dict["D"], mat_dict["T"]]        # lazy-loads ENDF cross sections
    for m in mats:
        m.l_n2n = bool(include_n2n)
    _nesst.init_DT_scatter(E_eV, E_eV)
    if not include_n2n:
        for m in mats:
            m.n2n_dNdE = np.zeros(E_eV.size)     # seed so *_scatter_spec can read it
    _MATRIX_CACHE_KEY = sig


_CARBON_CACHE_KEY: Optional[Tuple] = None


def _ensure_carbon_matrices(E_eV: np.ndarray) -> "object":
    """Initialise (or reuse) the NeSST C12 scatter matrices on ``E_eV`` and
    return the C12 material object. Carbon (n,2n) has a high threshold and is
    left off (elastic only) to save memory."""
    global _CARBON_CACHE_KEY
    _require_nesst()
    from NeSST.core import mat_dict
    sig = _grid_signature(E_eV)
    cmat = mat_dict.get("C12") if hasattr(mat_dict, "get") else None
    if sig != _CARBON_CACHE_KEY or cmat is None:
        logger.info("NeSST: initialising C12 scatter matrices on %d-point grid.",
                    E_eV.size)
        cmat = mat_dict["C12"]
        if hasattr(cmat, "l_n2n"):
            cmat.l_n2n = False
        _nesst.init_mat_scatter(E_eV, E_eV, "C12")
        if hasattr(cmat, "l_n2n") and not cmat.l_n2n:
            cmat.n2n_dNdE = np.zeros(E_eV.size)
        _CARBON_CACHE_KEY = sig
    return mat_dict["C12"]


def clear_matrix_cache() -> None:
    """Force the next call to re-initialise the NeSST scatter matrices."""
    global _CARBON_CACHE_KEY
    _CARBON_CACHE_KEY = None
    global _MATRIX_CACHE_KEY
    _MATRIX_CACHE_KEY = None


# ---------------------------------------------------------------------------
# Core: scattered spectrum from a primary spectrum + areal density
# ---------------------------------------------------------------------------

def scattered_spectrum(
    energy_MeV: np.ndarray,
    primary_spectrum: np.ndarray,
    rhoR_gcm2: float,
    frac_D: float = 0.5,
    frac_T: float = 0.5,
    n_E: int = 500,
    e_lo_MeV: float = 1.0,
    e_hi_MeV: float = 18.0,
    include_n2n: bool = True,
    rhoL_func: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> Dict:
    """Physically singly-scattered DT neutron spectrum from NeSST.

    Parameters
    ----------
    energy_MeV, primary_spectrum
        The birth spectrum ``I(E)`` (counts/MeV), e.g. from
        :func:`neutron_spectrum.synthesize_birth_spectrum`. Its integral is taken
        as the total primary yield and preserved in the returned spectra.
    rhoR_gcm2
        DT-neutron-weighted fuel areal density (g/cm^2) the neutrons traverse.
    frac_D, frac_T
        Fuel atom fractions (equimolar DT = 0.5/0.5; a D-enhanced 70:30 target is
        0.7/0.3). These weight n+D vs n+T scattering *and* the rhoR<->A1s
        calibration, so pass the run's real loading.
    n_E, e_lo_MeV, e_hi_MeV
        Uniform NeSST energy grid: ``n_E`` points over ``[e_lo, e_hi]`` MeV.
        With ``include_n2n`` the (n,2n) matrix memory scales as ``n_E**2``
        (500 pts ~ 1.3 GB, 900 pts OOMs an 8 GB box); with it off, 1200+ pts fit
        in ~0.3 GB.
    include_n2n
        Include D(n,2n)/T(n,2n) channels (default True). They add ~1.6% to the
        standard DSR but cost most of the memory; set False for a fast,
        fine-grid, elastic-only areal-density estimate.
    rhoL_func
        Optional normalised angular areal-density ``f(mu)``, ``mu in [-1, 1]``,
        ``int f dmu = 1`` -> uses the anisotropic model. Default (None) is the
        isotropic (symmetric) model appropriate for a 1D spherical run.

    Returns
    -------
    dict with keys:
        ``energy_MeV``      output grid (MeV)
        ``primary``         primary spectrum on the grid (counts/MeV)
        ``scattered``       total singly-scattered spectrum (counts/MeV)
        ``full``            primary + scattered (counts/MeV)
        ``components``      {'nD','nT','Dn2n','Tn2n'} scattered pieces (counts/MeV)
        ``A_1S``            NeSST scattering amplitude used
        ``rhoR_gcm2``, ``frac_D``, ``frac_T``, ``yield``  (echo of inputs)
    """
    _require_nesst()
    energy_MeV = np.asarray(energy_MeV, dtype=float)
    primary_spectrum = np.asarray(primary_spectrum, dtype=float)

    yield_tot = float(trapezoid(primary_spectrum, energy_MeV))
    if not np.isfinite(yield_tot) or yield_tot <= 0:
        raise ValueError("primary_spectrum has non-positive integral (no yield).")

    E_eV = np.linspace(e_lo_MeV, e_hi_MeV, int(n_E)) * 1e6
    E_MeV = E_eV / 1e6
    _ensure_scatter_matrices(E_eV, include_n2n=include_n2n)

    # Interpolate the primary onto the NeSST grid and normalise to unit integral
    # (magnitude cancels in the ratio; we restore the true yield at the end).
    I_grid = np.interp(E_MeV, energy_MeV, primary_spectrum, left=0.0, right=0.0)
    I_norm = float(trapezoid(I_grid, E_eV))
    if I_norm <= 0:
        raise ValueError("Primary spectrum does not overlap the NeSST energy grid.")
    I_E = I_grid / I_norm                                   # 1/eV, unit yield

    if rhoL_func is None:
        total_shape, (nD, nT, Dn2n, Tn2n) = _nesst.DT_sym_scatter_spec(
            I_E, frac_D, frac_T)
    else:
        total_shape, (nD, nT, Dn2n, Tn2n) = _nesst.DT_asym_scatter_spec(
            I_E, rhoL_func, frac_D, frac_T)

    A_1S = float(_nesst.rhoR_2_A1s(rhoR_gcm2 * G_CM2_TO_KG_M2, frac_D, frac_T))

    # Restore physical yield and convert 1/eV -> 1/MeV (factor 1e6).
    scale = yield_tot * 1e6
    primary_out = I_E * scale
    scattered_out = A_1S * total_shape * scale
    comp = {"nD": nD, "nT": nT, "Dn2n": Dn2n, "Tn2n": Tn2n}
    comp = {k: A_1S * v * scale for k, v in comp.items()}

    return {
        "energy_MeV": E_MeV,
        "primary": primary_out,
        "scattered": scattered_out,
        "full": primary_out + scattered_out,
        "components": comp,
        "A_1S": A_1S,
        "rhoR_gcm2": float(rhoR_gcm2),
        "frac_D": float(frac_D), "frac_T": float(frac_T),
        "yield": yield_tot,
        # original fine-grid birth spectrum, retained so the primary-peak Ti is
        # read at full resolution (the coarse scatter grid broadens the peak).
        "birth_energy_MeV": energy_MeV,
        "birth_spectrum": primary_spectrum,
    }


# ---------------------------------------------------------------------------
# Down-scatter ratio from a scattered spectrum
# ---------------------------------------------------------------------------

def downscatter_ratio(
    energy_MeV: np.ndarray,
    full_spectrum: np.ndarray,
    scattered_spectrum_arr: Optional[np.ndarray] = None,
    scatter_window_MeV: Tuple[float, float] = DSR_SCATTER_WINDOW_MEV,
    primary_window_MeV: Tuple[float, float] = DSR_PRIMARY_WINDOW_MEV,
) -> Dict:
    """DSR = N(scatter window) / N(primary window).

    If ``scattered_spectrum_arr`` is given, its integral is used for the
    numerator (cleaner, since the primary tail is excluded); otherwise the full
    spectrum is integrated over the scatter window. The denominator is always the
    full spectrum over the primary window (dominated by the primary peak)."""
    energy_MeV = np.asarray(energy_MeV, dtype=float)
    full_spectrum = np.asarray(full_spectrum, dtype=float)
    slo, shi = scatter_window_MeV
    plo, phi = primary_window_MeV
    smask = (energy_MeV >= slo) & (energy_MeV <= shi)
    pmask = (energy_MeV >= plo) & (energy_MeV <= phi)

    num_src = (np.asarray(scattered_spectrum_arr, dtype=float)
               if scattered_spectrum_arr is not None else full_spectrum)
    N_scatter = float(trapezoid(num_src[smask], energy_MeV[smask])) if smask.any() else 0.0
    N_primary = float(trapezoid(full_spectrum[pmask], energy_MeV[pmask])) if pmask.any() else 0.0
    dsr = N_scatter / N_primary if N_primary > 0 else float("nan")
    return {
        "DSR": dsr, "N_scatter": N_scatter, "N_primary": N_primary,
        "scatter_window_MeV": scatter_window_MeV,
        "primary_window_MeV": primary_window_MeV,
    }


def dsr_rhoR_slope(
    frac_D: float = 0.5, frac_T: float = 0.5, Tion_keV: float = 5.0,
    rhoR_ref_gcm2: float = 0.5, n_E: int = 500, include_n2n: bool = True,
    scatter_window_MeV: Tuple[float, float] = DSR_SCATTER_WINDOW_MEV,
    primary_window_MeV: Tuple[float, float] = DSR_PRIMARY_WINDOW_MEV,
) -> Dict:
    """First-principles DSR<->rhoR calibration from NeSST.

    In the single-scatter approximation the scattered spectrum is linear in rhoR,
    so ``DSR = k * rhoR``. This runs the model once at ``rhoR_ref`` with a NeSST
    Ballabio primary at ``Tion_keV`` and returns the slope ``k`` and its inverse
    ``coeff = 1/k`` (g/cm^2 per unit DSR) -- the transport analogue of the
    empirical NIF 20.4. Independent of any Helios run."""
    _require_nesst()
    E_eV = np.linspace(1.0, 18.0, int(n_E)) * 1e6
    mean, var, _ = _nesst.DTprimspecmoments(Tion_keV * 1000.0)
    sig = np.sqrt(var)
    prim = np.exp(-0.5 * ((E_eV - mean) / sig) ** 2)
    res = scattered_spectrum(E_eV / 1e6, prim, rhoR_ref_gcm2, frac_D, frac_T,
                             n_E=n_E, include_n2n=include_n2n)
    dsr = downscatter_ratio(res["energy_MeV"], res["full"], res["scattered"],
                            scatter_window_MeV, primary_window_MeV)["DSR"]
    k = dsr / rhoR_ref_gcm2
    return {"slope_DSR_per_gcm2": k, "coeff_gcm2_per_DSR": (1.0 / k if k else float("nan")),
            "rhoR_ref_gcm2": rhoR_ref_gcm2, "DSR_ref": dsr,
            "Tion_keV": Tion_keV, "frac_D": frac_D, "frac_T": frac_T}


# ---------------------------------------------------------------------------
# High-level entry points
# ---------------------------------------------------------------------------

def scatter_from_extraction(
    nd, frac_D: float = 0.5, frac_T: float = 0.5, n_E: int = 500,
    include_n2n: bool = True,
    scatter_window_MeV: Tuple[float, float] = DSR_SCATTER_WINDOW_MEV,
    primary_window_MeV: Tuple[float, float] = DSR_PRIMARY_WINDOW_MEV,
) -> Dict:
    """Run the NeSST scatter model from a :class:`neutron_spectrum.NeutronicsData`.

    Uses the extraction's quick-look birth spectrum and DT-weighted
    ``rhoR_emission_gcm2``. Returns the :func:`scattered_spectrum` dict augmented
    with a ``dsr`` block and the NeSST-inferred rhoR."""
    ql = getattr(nd, "quicklook", None) or {}
    energy = ql.get("energy_MeV")
    primary = ql.get("birth_spectrum")
    rhoR = ql.get("rhoR_emission_gcm2")
    if energy is None or primary is None or not np.isfinite(rhoR):
        raise ValueError("NeutronicsData lacks a usable birth spectrum / rhoR.")
    return scatter_from_spectrum(energy, primary, rhoR, frac_D, frac_T, n_E,
                                 include_n2n, scatter_window_MeV, primary_window_MeV)


def scatter_from_spectrum(
    energy_MeV: np.ndarray, primary_spectrum: np.ndarray, rhoR_gcm2: float,
    frac_D: float = 0.5, frac_T: float = 0.5, n_E: int = 500,
    include_n2n: bool = True,
    scatter_window_MeV: Tuple[float, float] = DSR_SCATTER_WINDOW_MEV,
    primary_window_MeV: Tuple[float, float] = DSR_PRIMARY_WINDOW_MEV,
) -> Dict:
    """Scattered spectrum + DSR + inferred rhoR for a birth spectrum and rhoR."""
    res = scattered_spectrum(energy_MeV, primary_spectrum, rhoR_gcm2,
                             frac_D, frac_T, n_E=n_E, include_n2n=include_n2n)
    dsr = downscatter_ratio(res["energy_MeV"], res["full"], res["scattered"],
                            scatter_window_MeV, primary_window_MeV)
    # NeSST first-principles inverse: rhoR from the scattering amplitude.
    rhoR_nesst = float(_nesst.A1s_2_rhoR(res["A_1S"], frac_D, frac_T)) / G_CM2_TO_KG_M2
    res["dsr"] = dsr
    res["rhoR_input_gcm2"] = float(rhoR_gcm2)
    res["rhoR_from_A1s_gcm2"] = rhoR_nesst          # should round-trip to input
    logger.info("NeSST scatter: rhoR=%.3f g/cm^2 -> DSR(%.0f-%.0f MeV)=%.3f%%  "
                "(A_1S=%.3e)", rhoR_gcm2, scatter_window_MeV[0],
                scatter_window_MeV[1], 100.0 * dsr["DSR"], res["A_1S"])
    return res


# ---------------------------------------------------------------------------
# Multi-material scatter (D + T + ablator/foam carbon)
# ---------------------------------------------------------------------------

def multi_material_scatter(
    energy_MeV: np.ndarray, primary_spectrum: np.ndarray,
    rhoR_D_gcm2: float, rhoR_T_gcm2: float, rhoR_C_gcm2: float = 0.0,
    n_E: int = 900, e_lo_MeV: float = 1.0, e_hi_MeV: float = 18.0,
    scatter_window_MeV: Tuple[float, float] = DSR_SCATTER_WINDOW_MEV,
    primary_window_MeV: Tuple[float, float] = DSR_PRIMARY_WINDOW_MEV,
) -> Dict:
    """Single-scatter spectrum from *per-species* mass areal densities.

    Unlike :func:`scattered_spectrum` (one fuel fraction, D+T only), this takes
    the D, T and carbon mass areal densities separately -- so it correctly
    handles a layer-varying target (wetted foam, doped ablator) and includes the
    ablator / foam **carbon** down-scatter that a detector actually sees.

    D and T go through NeSST's DT machinery (atom fractions derived from
    ``rhoR_D``/``rhoR_T``); carbon goes through ``init_mat_scatter('C12')``.
    Each species is scaled by its own NeSST scattering amplitude
    ``rhoR_2_A1s(rhoR_species)``. Elastic-only (fine grid, low memory).

    Returns a dict with the primary, the per-species scattered spectra
    (``scattered_D/T/C``), the fuel (D+T) and total scattered spectra, the full
    spectrum, and two DSRs:
      ``dsr_fuel``  -- fuel-only (the classic diagnostic; matches D+T scatter),
      ``dsr_total`` -- including carbon (the full down-scattered signal).
    """
    _require_nesst()
    from NeSST.core import mat_dict
    energy_MeV = np.asarray(energy_MeV, dtype=float)
    primary_spectrum = np.asarray(primary_spectrum, dtype=float)

    yield_tot = float(trapezoid(primary_spectrum, energy_MeV))
    if not np.isfinite(yield_tot) or yield_tot <= 0:
        raise ValueError("primary_spectrum has non-positive integral (no yield).")

    E_eV = np.linspace(e_lo_MeV, e_hi_MeV, int(n_E)) * 1e6
    E_MeV = E_eV / 1e6
    _ensure_scatter_matrices(E_eV, include_n2n=False)
    cmat = _ensure_carbon_matrices(E_eV)

    I_grid = np.interp(E_MeV, energy_MeV, primary_spectrum, left=0.0, right=0.0)
    I_norm = float(trapezoid(I_grid, E_eV))
    if I_norm <= 0:
        raise ValueError("Primary spectrum does not overlap the NeSST energy grid.")
    I_E = I_grid / I_norm
    scale = yield_tot * 1e6                             # 1/eV -> 1/MeV, restore yield

    # --- fuel (D + T): atom fractions from the mass areal densities ---
    rhoR_fuel = rhoR_D_gcm2 + rhoR_T_gcm2
    nD_a = rhoR_D_gcm2 / 2.014102
    nT_a = rhoR_T_gcm2 / 3.016049
    frac_D = nD_a / (nD_a + nT_a) if (nD_a + nT_a) > 0 else 0.5
    frac_T = 1.0 - frac_D
    dt_shape, _ = _nesst.DT_sym_scatter_spec(I_E, frac_D, frac_T)
    A1s_DT = float(_nesst.rhoR_2_A1s(rhoR_fuel * G_CM2_TO_KG_M2, frac_D, frac_T))
    scat_DT = A1s_DT * dt_shape * scale

    # --- carbon (C12) ---
    if rhoR_C_gcm2 > 0:
        c_out = _nesst.mat_scatter_spec(cmat, I_E, lambda x: np.ones_like(x))
        c_shape = c_out[0] if isinstance(c_out, tuple) else c_out
        A1s_C = float(cmat.rhoR_2_A1s(rhoR_C_gcm2 * G_CM2_TO_KG_M2))
        scat_C = A1s_C * np.asarray(c_shape, dtype=float) * scale
    else:
        scat_C = np.zeros_like(E_MeV)
        A1s_C = 0.0

    primary_out = I_E * scale
    scat_fuel = scat_DT
    scat_total = scat_fuel + scat_C
    dsr_fuel = downscatter_ratio(E_MeV, primary_out + scat_fuel, scat_fuel,
                                 scatter_window_MeV, primary_window_MeV)
    dsr_total = downscatter_ratio(E_MeV, primary_out + scat_total, scat_total,
                                  scatter_window_MeV, primary_window_MeV)

    logger.info("multi-material scatter: rhoR D=%.3f T=%.3f C=%.3f g/cm^2 -> "
                "DSR_fuel=%.3f%% DSR_total=%.3f%% (carbon adds %.1f%% of the "
                "down-scatter)", rhoR_D_gcm2, rhoR_T_gcm2, rhoR_C_gcm2,
                100 * dsr_fuel["DSR"], 100 * dsr_total["DSR"],
                100 * (dsr_total["DSR"] - dsr_fuel["DSR"]) / dsr_total["DSR"]
                if dsr_total["DSR"] else 0.0)

    return {
        "energy_MeV": E_MeV,
        "primary": primary_out,
        "scattered_D": None, "scattered_T": None,   # combined below (DT coupled)
        "scattered_fuel": scat_fuel, "scattered_C": scat_C,
        "scattered": scat_total, "full": primary_out + scat_total,
        "dsr": dsr_fuel, "dsr_fuel": dsr_fuel, "dsr_total": dsr_total,
        "rhoR_D_gcm2": rhoR_D_gcm2, "rhoR_T_gcm2": rhoR_T_gcm2,
        "rhoR_C_gcm2": rhoR_C_gcm2, "rhoR_fuel_gcm2": rhoR_fuel,
        "frac_D": frac_D, "frac_T": frac_T,
        "A_1S": A1s_DT, "A_1S_C": A1s_C,
        "birth_energy_MeV": energy_MeV, "birth_spectrum": primary_spectrum,
        "yield": yield_tot,
    }


def scattered_tof(
    scatter_result: Dict, distance_m: float = 3.0, t0_ns: float = 0.0,
    irf_fwhm_ns: float = 0.0, n_t: int = 4000, reaction: str = "DT",
) -> Dict:
    """Forward the NeSST *full* spectrum (primary + down-scatter) through the
    nTOF model so the trace carries the late-time scattered shoulder.

    Returns the :func:`neutron_tof.synthetic_ntof` dict with an added
    ``dsr`` block copied from ``scatter_result``."""
    from . import neutron_tof as tof
    from . import neutron_spectrum as ns
    out = tof.synthetic_ntof(
        energy_MeV=scatter_result["energy_MeV"], spectrum=scatter_result["full"],
        distance_m=distance_m, t0_ns=t0_ns, irf_fwhm_ns=irf_fwhm_ns, n_t=n_t,
        reaction=reaction)
    if out is not None:
        out["dsr"] = scatter_result.get("dsr")
        # Re-read the nTOF ion temperature from the *fine* birth spectrum: the
        # coarse scatter grid broadens the primary peak and would inflate Ti.
        be = scatter_result.get("birth_energy_MeV")
        bs = scatter_result.get("birth_spectrum")
        if be is not None and bs is not None:
            out["primary"]["Ti_ntof_keV"] = ns.infer_ntof_tion(
                be, bs, reaction)["Ti_ntof_keV"]
        # arrival time of the scatter-window edge, for gating the tail
        out["scatter_tail_ns"] = float(tof.energy_to_tof(
            scatter_result["dsr"]["scatter_window_MeV"][0], distance_m, t0_ns))
    return out


# ---------------------------------------------------------------------------
# Figure: scattered spectrum + synthetic nTOF (shared by the pipeline and the
# standalone analyze_neutron_scatter.py)
# ---------------------------------------------------------------------------

def plot_scatter_tof(scat: Dict, tof: Optional[Dict], out_path: str,
                     title: str = "") -> str:
    """Two-panel neutron figure -- (L) singly-scattered spectrum with the nD/nT
    backscatter edges and the DSR window, (R) synthetic nTOF trace -- saved to
    ``out_path``. Returns the path written."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    E = scat["energy_MeV"]
    axL.semilogy(E, scat["full"], "k", lw=1.7, label="full")
    axL.semilogy(E, scat["primary"], "--", color="#1f4e79", lw=1.2, label="primary")
    axL.semilogy(E, scat["components"]["nT"], color="#2e7d32", lw=1.2, label="n+T")
    axL.semilogy(E, scat["components"]["nD"], color="#b26a00", lw=1.2, label="n+D")
    axL.axvspan(10, 12, color="#cccccc", alpha=0.4)
    ymax = float(np.max(scat["full"]))
    if ymax > 0:
        axL.set_ylim(ymax * 1e-5, ymax * 5)
    axL.set_xlim(1, 16)
    axL.set_xlabel("neutron energy (MeV)")
    axL.set_ylabel("dN/dE (arb.)")
    axL.set_title(f"scattered spectrum   DSR={100 * scat['dsr']['DSR']:.2f}%")
    axL.legend(fontsize=8)
    axL.grid(alpha=0.25, which="both")

    if tof:
        t = tof["time_ns"]
        y = tof["signal_ideal"]
        ypos = np.clip(y, np.max(y) * 1e-6, None) if np.max(y) > 0 else y
        axR.semilogy(t, ypos, "k", lw=1.4)
        axR.axvline(tof["primary"]["t_peak_ns"], color="#1f4e79", ls="--", lw=1)
        axR.set_xlabel("time of flight (ns)")
        axR.set_ylabel("signal (arb.)")
        axR.set_title(f"synthetic nTOF @ {tof['distance_m']:.1f} m   "
                      f"T_i={tof['primary']['Ti_ntof_keV']:.1f} keV")
        axR.grid(alpha=0.25, which="both")

    if title:
        fig.suptitle(title, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95] if title else None)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return str(out_path)


# ---------------------------------------------------------------------------
# Pipeline entry point: full neutron block for run_analysis
# ---------------------------------------------------------------------------

def _pub_scalar(pub, keys):
    """First matching published value (handles [value, unc] pairs and scalars)."""
    if not pub:
        return None
    for k in keys:
        if k in pub:
            v = pub[k]
            try:
                return float(v[0]) if isinstance(v, (list, tuple)) else float(v)
            except (TypeError, ValueError):
                return None
    return None


def neutron_report(data, frac_D: float = 0.5, frac_T: float = 0.5,
                   distance_m: float = 3.0, n_E: int = 500,
                   include_n2n: bool = True, published: Optional[Dict] = None,
                   plot_path: Optional[str] = None, plot_title: str = "",
                   rhw_path: Optional[str] = None
                   ) -> Tuple[Optional[Dict], str]:
    """Full neutron post-processing block for a Helios run, for the standard
    pipeline (``run_analysis``). Returns ``(metrics, text)``.

    NeSST-optional: with NeSST it reports the transport DSR (10-12/13-15 MeV),
    the rhoR round-trip, the first-principles coefficient, and the nTOF Ti; if
    NeSST is missing it still reports the birth-spectrum nTOF Ti and the
    quick-look DSR, and says the transport DSR needs NeSST. Never raises on a
    scatter failure -- it degrades to the quick-look. Returns ``(None, note)``
    for a no-burn run.

    If ``published`` is given (a mapping of neutron observables, values as
    ``[value, unc]`` or scalars), a compact ours-vs-published overlay is
    appended to the text block.

    If ``plot_path`` is given and NeSST produced a scattered spectrum, the
    two-panel figure (spectrum + nTOF) is written there and the path is recorded
    in ``metrics['plot_path']``."""
    from . import neutron_spectrum as ns
    nd = ns.extract_neutronics(data=data, use_rhino=False)
    if nd is None:
        return None, "  (no-burn run: no DT neutron spectrum to analyze)\n"

    ql = nd.quicklook or {}
    metrics = {
        "bang_time_ns": nd.bang_time_ns,
        "dt_yield": nd.total_dt_yield,
        "dd_yield": nd.total_dd_yield,
        "Ti_burn_avg_keV": ql.get("Ti_burn_avg_keV"),
        "burn_fwhm_ns": ql.get("burn_fwhm_ns"),
        "rhoR_hydro_gcm2": float(ql.get("rhoR_emission_gcm2", float("nan"))),
        "Ti_ntof_keV": ql.get("ntof_tion_keV"),      # birth-spectrum nTOF Ti
        "DSR": ql.get("dsr"),                        # quick-look placeholder
        "dsr_source": "quick-look (rhoR/coeff)",
        "frac_D": frac_D, "frac_T": frac_T, "distance_m": distance_m,
    }

    scat = tof = None
    if NESST_AVAILABLE:
        try:
            # Multi-material path: if an RHW is available, read per-layer
            # composition, split rhoR into D/T/carbon, and include ablator/foam
            # carbon down-scatter (fuel fractions auto-detected -- no manual
            # frac_D/frac_T needed).
            comp = _composition_scatter(nd, data, rhw_path, n_E) \
                if rhw_path else None
            if comp is not None:
                scat = comp
                metrics.update({
                    "DSR": scat["dsr_fuel"]["DSR"],
                    "DSR_total": scat["dsr_total"]["DSR"],
                    "rhoR_C_gcm2": scat["rhoR_C_gcm2"],
                    "frac_D": scat["frac_D"], "frac_T": scat["frac_T"],
                    "coeff_gcm2_per_DSR": (metrics["rhoR_hydro_gcm2"]
                                           / scat["dsr_fuel"]["DSR"])
                    if scat["dsr_fuel"]["DSR"] else float("nan"),
                    "dsr_source": "NeSST multi-material (D+T+C, ENDF/B-VIII.0)",
                    "composition_source": "RHW",
                })
            else:
                scat = scatter_from_extraction(nd, frac_D=frac_D, frac_T=frac_T,
                                               n_E=n_E, include_n2n=include_n2n)
                dsr = scat["dsr"]["DSR"]
                metrics.update({
                    "DSR": dsr,
                    "rhoR_from_DSR_gcm2": scat["rhoR_from_A1s_gcm2"],
                    "coeff_gcm2_per_DSR": (metrics["rhoR_hydro_gcm2"] / dsr)
                    if dsr else float("nan"),
                    "dsr_source": "NeSST single-scatter D+T (ENDF/B-VIII.0)",
                })
            tof = scattered_tof(scat, distance_m=distance_m)
            if tof:
                metrics["Ti_ntof_keV"] = tof["primary"]["Ti_ntof_keV"]
        except Exception as exc:                      # pragma: no cover
            logger.warning("neutron_report: NeSST scatter failed (%s); "
                           "quick-look DSR only.", exc)
            scat = tof = None

    if plot_path and scat is not None:
        try:
            metrics["plot_path"] = plot_scatter_tof(scat, tof, plot_path,
                                                    title=plot_title)
        except Exception as exc:                      # pragma: no cover
            logger.warning("neutron_report: figure not written (%s).", exc)

    return metrics, _format_neutron_block(metrics, published)


def _composition_scatter(nd, data, rhw_path, n_E) -> Optional[Dict]:
    """Multi-material scatter from an RHW: per-layer composition -> per-species
    areal densities -> D+T+C scatter. Returns the scatter dict, or None if the
    RHW / composition / geometry is unavailable (caller falls back to D+T)."""
    try:
        from . import target_composition as tc
    except Exception:                                 # pragma: no cover
        return None
    from pathlib import Path
    if not (rhw_path and Path(rhw_path).exists()):
        return None
    ql = nd.quicklook or {}
    rhoR_fuel = ql.get("rhoR_emission_gcm2")
    energy = ql.get("energy_MeV")
    primary = ql.get("birth_spectrum")
    if energy is None or primary is None or not np.isfinite(rhoR_fuel):
        return None
    regions = tc.parse_rhw_regions(rhw_path)
    if not regions:
        return None
    sp = tc.species_areal_densities(data, regions, float(rhoR_fuel))
    if sp is None:
        return None
    res = multi_material_scatter(
        energy, primary, sp["rhoR_D_gcm2"], sp["rhoR_T_gcm2"],
        sp["rhoR_C_gcm2"], n_E=max(int(n_E), 800))
    res["composition"] = sp
    res["regions_summary"] = tc.summarize_regions(regions)
    return res


def _format_neutron_block(m: Dict, pub: Optional[Dict] = None) -> str:
    dsr = m.get("DSR")
    dsr_pct = (100.0 * dsr) if dsr is not None else float("nan")
    L = []
    L.append("=" * 72)
    L.append("  NEUTRON DIAGNOSTICS  (birth spectrum -> DSR -> nTOF)")
    L.append("=" * 72)
    L.append(f"  bang time              {m['bang_time_ns']:.3f} ns")
    L.append(f"  DT neutron yield       {m['dt_yield']:.3e}")
    if m.get("dd_yield"):
        L.append(f"  DD neutron yield       {m['dd_yield']:.3e}")
    if m.get("burn_fwhm_ns") is not None:
        L.append(f"  burn width (FWHM)      {m['burn_fwhm_ns']:.3f} ns")
    if m.get("Ti_burn_avg_keV") is not None:
        L.append(f"  T_ion (burn-avg)       {m['Ti_burn_avg_keV']:.2f} keV")
    if m.get("Ti_ntof_keV") is not None:
        L.append(f"  T_ion (nTOF)           {m['Ti_ntof_keV']:.2f} keV")
    L.append("-" * 72)
    L.append(f"  rhoR (hydro, DT-wt)    {m['rhoR_hydro_gcm2']:.3f} g/cm^2")
    L.append(f"  DSR fuel (D+T)         {dsr_pct:.3f} %   [{m['dsr_source']}]")
    if m.get("DSR_total") is not None:
        dsr_tot_pct = 100.0 * m["DSR_total"]
        add = dsr_tot_pct - dsr_pct
        L.append(f"  DSR total (+carbon)    {dsr_tot_pct:.3f} %   "
                 f"(ablator/foam C adds {add:+.3f})")
        if m.get("rhoR_C_gcm2") is not None:
            L.append(f"  rhoR carbon (C)        {m['rhoR_C_gcm2']:.3f} g/cm^2")
    if m.get("rhoR_from_DSR_gcm2") is not None:
        L.append(f"  rhoR back from DSR     {m['rhoR_from_DSR_gcm2']:.3f} g/cm^2"
                 "   (round-trip check)")
    if m.get("coeff_gcm2_per_DSR") is not None:
        L.append(f"  first-principles coeff {m['coeff_gcm2_per_DSR']:.2f}"
                 " g/cm^2 per DSR   (empirical NIF 20.4)")
    src = m.get("composition_source")
    auto = "  (fuel D:T auto-detected from RHW)" if src == "RHW" else ""
    L.append(f"  fuel D:T = {m['frac_D']:.2f}:{m['frac_T']:.2f}"
             f"   nTOF @ {m['distance_m']:.1f} m{auto}")

    rows = [
        (["yield_neutrons", "yield", "DT_yield"], m.get("dt_yield"), "DT yield", "{:.3g}"),
        (["Tion", "T_ion", "T_ion_keV", "T_DT_keV"], m.get("Ti_ntof_keV"), "T_ion (keV)", "{:.3g}"),
        (["DSR", "dsr", "DSR_percent"], dsr_pct, "DSR (%)", "{:.3g}"),
        (["bang_time_ns", "bang_time"], m.get("bang_time_ns"), "bang time (ns)", "{:.3g}"),
    ]
    overlay = []
    for keys, ov, label, fmt in rows:
        pv = _pub_scalar(pub, keys)
        if pv is None or ov is None:
            continue
        overlay.append(f"  {label:20s}{fmt.format(ov):>14s}{fmt.format(pv):>16s}")
    if overlay:
        L.append("-" * 72)
        L.append(f"  {'neutron metric':20s}{'ours':>14s}{'published':>16s}")
        L.extend(overlay)
        L.append("  (published yield is neutron count; T_ion vs apparent DT nTOF; "
                 "DSR is the clean apples-to-apples check.)")
    L.append("=" * 72)
    return "\n".join(L) + "\n"
