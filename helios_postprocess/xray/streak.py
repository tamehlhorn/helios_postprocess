"""
Synthetic x-ray streak camera on capsule self-emission.

Builds the master radiance cube I(t, p, E) from a Helios run, then projects it
two ways:

  * IMAGING streak   S(t, x)  -- slit image vs sweep time, one filtered band
  * SPECTRAL streak  S(t, E)  -- spatially integrated spectrum vs sweep time

Both come from the same cube, so the imaging and spectral instruments are
guaranteed consistent with each other and with the Level 0 / SPECT3D rung
they were built from.

Instrument effects applied, in order:
  1. channel response R(E)                     (response.py)
  2. resampling onto a uniform sweep time base (Helios output is NOT uniform
     in time -- it clusters steps around stagnation)
  3. temporal IRF                              (Gaussian, camera time res)
  4. spatial MTF                               (Gaussian, pinhole/slit + optic)

Observables extracted:
  * limb / emission-edge trajectory R_x(t)
  * x-ray bang time and emission FWHM
  * stagnation flash size

Author: Prof T / helios_postprocess
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter1d

from .chords import ChordResult, integrate, make_impact_grid, spatially_integrate
from .emissivity import build_emissivity
from .response import Channel

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
@dataclass
class StreakConfig:
    """
    Instrument and sampling configuration.

    Parameters
    ----------
    t_start_ns, t_stop_ns : float
        Sweep window.  ``None`` uses the full run.
    n_time : int
        Samples on the uniform sweep time base.
    n_impact : int
        Impact parameters (spatial samples along the slit, object plane).
    E_min_eV, E_max_eV, n_energy : float, float, int
        Log-spaced photon energy grid.
    time_resolution_ps : float
        FWHM of the temporal IRF (photocathode + sweep + phosphor + optics).
    spatial_resolution_um : float
        FWHM of the spatial blur referred to the OBJECT plane.
    magnification : float
        Recorded for bookkeeping; the cube is kept in object-plane units.
    optically_thin : bool
        Level 0a (True) or LTE-continuum formal solution (False).
    r_max_cm : float, optional
        Outer radius of the impact grid.  Defaults to the initial outer
        radius of the included zones.
    zone_slice : slice, optional
        Restrict emission to a zone range.  Essential for halfraum targets:
        the Cu converter and He fill would otherwise dominate the signal.
    rho_floor_g_cc : float, optional
        Treat zones below this density as non-emitting.  Suppresses the
        massless-and-very-hot outer corona zones that can manufacture a
        spurious hard tail and a spurious emission halo outside the limb.
    """
    t_start_ns: Optional[float] = None
    t_stop_ns: Optional[float] = None
    n_time: int = 400
    n_impact: int = 128
    E_min_eV: float = 200.0
    E_max_eV: float = 15000.0
    n_energy: int = 96
    time_resolution_ps: float = 15.0
    spatial_resolution_um: float = 10.0
    magnification: float = 10.0
    optically_thin: bool = True
    r_max_cm: Optional[float] = None
    zone_slice: Optional[slice] = None
    rho_floor_g_cc: Optional[float] = None


def _emission_bounded_rmax(data, sel, cfg, enclose: float = 0.99,
                           margin: float = 2.0) -> float:
    """
    Choose the impact-parameter extent from where the light actually is.

    Using the mesh outer boundary is the obvious choice and the wrong one: a
    tenuous corona can sit at tens of mm while every emitting zone is inside
    a few hundred um, which spends the entire impact grid on vacuum and
    quantizes the extracted flash radius onto grid nodes.

    The bound is the largest radius, over the sweep window, enclosing
    ``enclose`` of a bremsstrahlung-like emission proxy (rho^2 sqrt(Te) dV),
    times ``margin``.
    """
    rho = np.asarray(data.mass_density, float)
    rb = np.asarray(data.zone_boundaries, float)
    Te = np.asarray(getattr(data, "elec_temperature",
                            data.ion_temperature), float)
    vol = (4.0 / 3.0) * np.pi * (rb[:, 1:] ** 3 - rb[:, :-1] ** 3)

    sl = cfg.zone_slice if cfg.zone_slice is not None else slice(None)
    stop = (len(rb[0]) - 1 if cfg.zone_slice is None
            or cfg.zone_slice.stop is None else cfg.zone_slice.stop)

    best = 0.0
    step = max(1, sel.size // 60)
    for it in sel[::step]:
        w = (rho[it, sl] ** 2) * np.sqrt(np.maximum(Te[it, sl], 1e-3)) * vol[it, sl]
        tot = w.sum()
        if tot <= 0:
            continue
        c = np.cumsum(w) / tot
        k = int(np.searchsorted(c, enclose))
        k = min(k, w.size - 1)
        best = max(best, float(rb[it, sl.start or 0:stop + 1][k + 1]))

    hard_cap = float(rb[sel[0], stop])
    if best <= 0.0:
        logger.warning("  [xray] emission-bounded r_max found no emission; "
                       "falling back to the mesh outer boundary. Set "
                       "r_max_cm / --r-max-um explicitly.")
        return hard_cap

    chosen = min(best * margin, hard_cap)
    if chosen >= 0.99 * hard_cap:
        logger.warning(
            f"  [xray] emission-bounded r_max hit the mesh outer boundary "
            f"({hard_cap * 1e4:.0f} um). Either the corona genuinely emits "
            f"out there, or the {enclose:.0%} enclosure tail is being set by "
            f"a tenuous halo. Impact-grid resolution on the capsule will be "
            f"poor -- set --r-max-um explicitly (a few times the largest "
            f"emission-edge radius in the window).")
    return chosen


@dataclass
class RadianceCube:
    """
    Master cube on the native (non-uniform) Helios time base.

    Attributes
    ----------
    t_ns : (n_t,)
    p_cm : (n_p,)
    E_eV : (n_E,)
    I : (n_t, n_p, n_E)
        erg s^-1 cm^-2 sr^-1 eV^-1
    tau : (n_t, n_p, n_E)
        Chord optical depth (zeros if run optically thin).
    I_thin : (n_t, n_p, n_E)
        Optically-thin intensity, always retained.
    meta : dict
    """
    t_ns: np.ndarray
    p_cm: np.ndarray
    E_eV: np.ndarray
    I: np.ndarray
    tau: np.ndarray
    I_thin: np.ndarray
    meta: Dict = field(default_factory=dict)

    def power_spectrum(self) -> np.ndarray:
        """Total emitted spectral power vs time, (n_t, n_E) [erg/s/eV]."""
        w = 2.0 * np.pi * self.p_cm
        return 4.0 * np.pi * np.trapezoid(self.I * w[None, :, None],
                                          self.p_cm, axis=1)


# ---------------------------------------------------------------------------
# Cube construction
# ---------------------------------------------------------------------------
def build_radiance_cube(data, cfg: StreakConfig,
                        verbose: bool = True) -> RadianceCube:
    """
    Run the Level 0 forward model over every Helios time step in the window.

    Note on the impact grid: it is held FIXED in cm across all time steps so
    that the cube is a rectangular array and the streak image has a stationary
    spatial axis, exactly like the real instrument.  Chords with p larger than
    the instantaneous outer radius simply return zero.
    """
    t = np.asarray(data.time, float)

    t0 = t[0] if cfg.t_start_ns is None else cfg.t_start_ns
    t1 = t[-1] if cfg.t_stop_ns is None else cfg.t_stop_ns
    sel = np.flatnonzero((t >= t0) & (t <= t1))
    if sel.size == 0:
        raise ValueError(f"No Helios time steps in [{t0}, {t1}] ns")

    E = np.logspace(np.log10(cfg.E_min_eV), np.log10(cfg.E_max_eV), cfg.n_energy)

    r_max = cfg.r_max_cm
    if r_max is None:
        r_max = _emission_bounded_rmax(data, sel, cfg)
    p = make_impact_grid(r_max, cfg.n_impact)

    n_t, n_p, n_E = sel.size, p.size, E.size
    I = np.zeros((n_t, n_p, n_E))
    tau = np.zeros((n_t, n_p, n_E))
    I_thin = np.zeros((n_t, n_p, n_E))

    for k, it in enumerate(sel):
        cube = build_emissivity(data, int(it), E,
                                zone_slice=cfg.zone_slice,
                                include_absorption=not cfg.optically_thin,
                                rho_floor=cfg.rho_floor_g_cc)
        res: ChordResult = integrate(cube, p, optically_thin=cfg.optically_thin)
        I[k] = res.I
        tau[k] = res.tau
        I_thin[k] = res.I_thin if res.I_thin is not None else res.I
        if verbose and (k % max(1, n_t // 10) == 0):
            logger.info(f"  [xray] {k + 1}/{n_t}  t = {t[it]:.3f} ns")

    # Effective spatial sampling near the object, which is what limits the
    # extracted flash radius and edge trajectory.
    dp_min = float(np.min(np.diff(p))) if p.size > 1 else float("nan")
    dp_at_half = float(np.diff(p)[max(0, p.size // 2 - 1)]) if p.size > 2 else np.nan

    meta = dict(optically_thin=cfg.optically_thin,
                rho_floor_g_cc=cfg.rho_floor_g_cc,
                r_max_cm=r_max,
                dp_min_um=dp_min * 1e4,
                dp_mid_um=dp_at_half * 1e4,
                zone_slice=str(cfg.zone_slice),
                n_helios_steps=int(n_t))

    # The synthetic camera cannot resolve anything faster than the hydro dump
    # cadence.  If the EXODUS output is coarser than the requested IRF the
    # temporal resolution is set by the simulation output interval, not by the
    # instrument, and any FWHM or bang-time claim inherits that limit.
    if cfg.r_max_cm is None:
        logger.info(f"  [xray] impact grid bounded by emission: "
                    f"r_max = {r_max * 1e4:.0f} um, spacing "
                    f"{dp_min * 1e4:.2f}-{dp_at_half * 1e4:.1f} um")

    if n_t > 1:
        dt_ps = float(np.median(np.diff(t[sel]))) * 1.0e3
        meta["median_dump_dt_ps"] = dt_ps
        if dt_ps > cfg.time_resolution_ps:
            logger.warning(
                f"  [xray] EXODUS dump interval {dt_ps:.1f} ps exceeds the "
                f"requested IRF {cfg.time_resolution_ps:.1f} ps -- effective "
                f"time resolution is set by the hydro output cadence, not the "
                f"camera. Re-run Helios with denser output before quoting a "
                f"bang time or emission FWHM.")

    return RadianceCube(t_ns=t[sel], p_cm=p, E_eV=E,
                        I=I, tau=tau, I_thin=I_thin, meta=meta)


# ---------------------------------------------------------------------------
# Instrument projection
# ---------------------------------------------------------------------------
@dataclass
class StreakImage:
    """A synthetic streak record."""
    kind: str                       # 'imaging' or 'spectral'
    t_ns: np.ndarray                # (n_time,) uniform sweep base
    axis: np.ndarray                # (n_axis,) cm (imaging) or eV (spectral)
    axis_label: str
    signal: np.ndarray              # (n_time, n_axis)
    channel: Optional[str] = None
    meta: Dict = field(default_factory=dict)


def _uniform_resample(t_src: np.ndarray, y: np.ndarray,
                      t_dst: np.ndarray) -> np.ndarray:
    """Linear interpolation of (n_t, n_axis) onto a uniform time base."""
    out = np.empty((t_dst.size, y.shape[1]), dtype=float)
    for j in range(y.shape[1]):
        out[:, j] = np.interp(t_dst, t_src, y[:, j])
    return out


def _apply_irf(signal: np.ndarray, t_ns: np.ndarray,
               time_res_ps: float, axis_vals: np.ndarray,
               spatial_res: float, spatial: bool) -> np.ndarray:
    """Gaussian temporal IRF and (optionally) spatial MTF, FWHM inputs."""
    out = signal
    if time_res_ps and time_res_ps > 0 and t_ns.size > 2:
        dt_ns = float(np.mean(np.diff(t_ns)))
        sigma_bins = (time_res_ps * 1e-3 / 2.3548) / dt_ns
        if sigma_bins > 0.3:
            out = gaussian_filter1d(out, sigma_bins, axis=0, mode="nearest")
    if spatial and spatial_res and spatial_res > 0 and axis_vals.size > 2:
        dx_um = float(np.mean(np.diff(axis_vals))) * 1e4
        sigma_bins = (spatial_res / 2.3548) / max(dx_um, 1e-9)
        if sigma_bins > 0.3:
            out = gaussian_filter1d(out, sigma_bins, axis=1, mode="nearest")
    return out


def make_imaging_streak(cube: RadianceCube, channel: Channel,
                        cfg: StreakConfig) -> StreakImage:
    """
    Band-integrate the cube through ``channel`` and sweep it.

    Signal units: erg s^-1 cm^-2 sr^-1 (band-integrated radiance), before
    detector gain.  Absolute normalization is deliberately not attempted --
    the useful comparisons (limb trajectory, bang time, shape) are all
    normalization-independent, and absolute photometry needs the real
    response curves.
    """
    R = channel.response(cube.E_eV)
    band = np.trapezoid(cube.I * R[None, None, :], cube.E_eV, axis=2)  # (n_t, n_p)

    t_uni = np.linspace(cube.t_ns[0], cube.t_ns[-1], cfg.n_time)
    swept = _uniform_resample(cube.t_ns, band, t_uni)
    swept = _apply_irf(swept, t_uni, cfg.time_resolution_ps,
                       cube.p_cm, cfg.spatial_resolution_um, spatial=True)

    return StreakImage(kind="imaging", t_ns=t_uni, axis=cube.p_cm,
                       axis_label="impact parameter (cm)", signal=swept,
                       channel=channel.name,
                       meta=dict(magnification=cfg.magnification, **cube.meta))


def make_spectral_streak(cube: RadianceCube, cfg: StreakConfig,
                         channel: Optional[Channel] = None) -> StreakImage:
    """
    Spatially integrate to the emitted spectrum and sweep it.

    Signal units: erg s^-1 eV^-1 (times response if a channel is supplied).
    This is the streaked-spectrometer view; the spatial information is
    integrated away exactly as a slitless / spatially-integrating
    spectrometer would.
    """
    spec = cube.power_spectrum()                      # (n_t, n_E)
    if channel is not None:
        spec = spec * channel.response(cube.E_eV)[None, :]

    t_uni = np.linspace(cube.t_ns[0], cube.t_ns[-1], cfg.n_time)
    swept = _uniform_resample(cube.t_ns, spec, t_uni)
    swept = _apply_irf(swept, t_uni, cfg.time_resolution_ps,
                       cube.E_eV, 0.0, spatial=False)

    return StreakImage(kind="spectral", t_ns=t_uni, axis=cube.E_eV,
                       axis_label="photon energy (eV)", signal=swept,
                       channel=None if channel is None else channel.name,
                       meta=dict(cube.meta))


# ---------------------------------------------------------------------------
# Observable extraction
# ---------------------------------------------------------------------------
def emission_edge_trajectory(img: StreakImage,
                             threshold_frac: float = 0.5) -> np.ndarray:
    """
    Track the outer emission edge R(t) from an imaging streak.

    At each time the radial profile is normalized to its own maximum and the
    outermost crossing of ``threshold_frac`` is found by linear interpolation.

    PHYSICS CAVEAT: this edge is the ablation front / coronal emission
    boundary, NOT the cold shell radius.  Comparing it directly to a
    Lagrangian shell trajectory from the hydro is an apples-to-oranges error
    that has bitten this diagnostic repeatedly in the literature.  Compare it
    to the same edge extracted from the simulation's own emission, which is
    what this function does.

    Returns (n_time,) in cm, NaN where no crossing exists.
    """
    p = img.axis
    out = np.full(img.t_ns.size, np.nan)
    for i, row in enumerate(img.signal):
        peak = row.max()
        if peak <= 0:
            continue
        f = row / peak
        above = np.flatnonzero(f >= threshold_frac)
        if above.size == 0:
            continue
        k = above[-1]
        if k >= p.size - 1:
            out[i] = p[-1]
            continue
        f1, f2 = f[k], f[k + 1]
        if f1 == f2:
            out[i] = p[k]
        else:
            w = (f1 - threshold_frac) / (f1 - f2)
            out[i] = p[k] + w * (p[k + 1] - p[k])
    return out


@dataclass
class BurnHistoryMetrics:
    """X-ray emission history metrics from a streak record."""
    bang_time_ns: float
    fwhm_ps: float
    rise_10_90_ps: float
    peak_signal: float
    flash_radius_cm: float


def emission_history(img: StreakImage) -> np.ndarray:
    """Spatially (or spectrally) integrated signal vs time, (n_time,)."""
    return np.trapezoid(img.signal, img.axis, axis=1)


def burn_metrics(img: StreakImage) -> BurnHistoryMetrics:
    """
    X-ray bang time (peak of the emission history), FWHM, and 10-90 rise.

    The x-ray bang time is NOT identical to the neutron bang time -- x-ray
    emission is weighted differently in Te and rho -- and the offset between
    them is itself a calibration observable worth reporting.
    """
    h = emission_history(img)
    t = img.t_ns
    if h.max() <= 0:
        return BurnHistoryMetrics(np.nan, np.nan, np.nan, 0.0, np.nan)

    ipk = int(np.argmax(h))
    bang = float(t[ipk])
    hn = h / h.max()

    def _crossing(level, side):
        if side == "left":
            idx = np.flatnonzero(hn[:ipk + 1] <= level)
            if idx.size == 0:
                return np.nan
            k = idx[-1]
            if k >= ipk:
                return t[ipk]
            w = (level - hn[k]) / (hn[k + 1] - hn[k])
            return t[k] + w * (t[k + 1] - t[k])
        idx = np.flatnonzero(hn[ipk:] <= level)
        if idx.size == 0:
            return np.nan
        k = ipk + idx[0]
        if k == 0:
            return t[0]
        w = (hn[k - 1] - level) / (hn[k - 1] - hn[k])
        return t[k - 1] + w * (t[k] - t[k - 1])

    t_half_l, t_half_r = _crossing(0.5, "left"), _crossing(0.5, "right")
    fwhm_ps = (t_half_r - t_half_l) * 1e3 if np.isfinite(t_half_r - t_half_l) else np.nan
    t10, t90 = _crossing(0.1, "left"), _crossing(0.9, "left")
    rise_ps = (t90 - t10) * 1e3 if np.isfinite(t90 - t10) else np.nan

    flash_r = np.nan
    if img.kind == "imaging":
        row = img.signal[ipk]
        if row.max() > 0:
            f = row / row.max()
            above = np.flatnonzero(f >= 0.5)
            if above.size:
                flash_r = float(img.axis[above[-1]])

    return BurnHistoryMetrics(bang_time_ns=bang, fwhm_ps=float(fwhm_ps),
                              rise_10_90_ps=float(rise_ps),
                              peak_signal=float(h.max()),
                              flash_radius_cm=float(flash_r))


# ---------------------------------------------------------------------------
# Optical-depth diagnostic -- the hand-off argument to SPECT3D
# ---------------------------------------------------------------------------
def opacity_report(cube: RadianceCube,
                   tau_warn: float = 0.3) -> Dict[str, float]:
    """
    Where does the optically-thin assumption break?

    Returns the peak chord optical depth, the photon energy and time at which
    it occurs, and the fraction of the (t, p, E) volume above ``tau_warn``.
    Also reports the maximum fractional intensity suppression I/I_thin.

    Only meaningful when the cube was built with ``optically_thin=False``.
    """
    if not np.any(cube.tau > 0):
        return dict(valid=0.0)

    imax = np.unravel_index(int(np.argmax(cube.tau)), cube.tau.shape)
    frac_thick = float(np.mean(cube.tau > tau_warn))

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(cube.I_thin > 0, cube.I / cube.I_thin, np.nan)
    min_ratio = float(np.nanmin(ratio))

    return dict(valid=1.0,
                tau_max=float(cube.tau[imax]),
                tau_max_time_ns=float(cube.t_ns[imax[0]]),
                tau_max_impact_cm=float(cube.p_cm[imax[1]]),
                tau_max_energy_eV=float(cube.E_eV[imax[2]]),
                frac_volume_above_warn=frac_thick,
                min_I_over_I_thin=min_ratio)
