"""
pdd_geometric_coupling.py
=========================

Pre-processor that estimates the time-dependent geometric coupling efficiency
of a polar-direct-drive (PDD) or hot-spot direct-drive (HDD) Helios run and
writes a corrected laser-power table that approximates what a 3D ray trace
(HYDRA / LILAC / Xrage) would compute self-consistently.

Physical motivation
-------------------
Helios is a 1D spherical hydrodynamics code.  Its laser ray trace models the
beam as a single cone of half-angle θ originating from a focal point at
distance d from the target centre.  All rays inside that cone propagate
symmetrically inward and refract according to the local plasma profile.

In real PDD / HDD geometry the beams are *fixed in lab frame* and pointed at a
fixed location.  As the capsule implodes its absorbing surface R(t) shrinks
toward the centre.  At late times R(t) is small compared to the beam footprint
in the focal plane, and a progressively larger fraction of the beam misses the
ablation surface entirely.  HYDRA / LILAC / Xrage capture this in 3D
automatically.  Helios cannot: in 1D the cone always intercepts the spherical
target by construction.

This pre-processor estimates the missing geometric correction f_geom(t) from a
no-burn Helios run, applies it to the original power table, and emits

  • a corrected power CSV (time, P_original, P_corrected, f_geom),
  • a diagnostic PDF, and
  • a brief summary of the integrated effect.

Two formulas are supported.  Both are anchored at t=0 so f_geom(0) = 1; the
correction describes the *evolution* of the coupling efficiency, not a static
offset.

  --formula view_factor   (default, Option 2 from the team memo)
      Solid-angle ratio between the sphere and a cone of half-angle θ:
          Ω(t) = 2π(1 - cos(arcsin(R(t)/d)))
          f_geom(t) = min(1, [1 - sqrt(1 - (R(t)/d)^2)] / [1 - cos(θ)])
      Beam-independent; uses only θ and d.

  --formula gaussian      (Option 1 from the team memo)
      Gaussian-beam disk projection onto a target of radius R(t):
          η(t) = 1 - exp(-2 R(t)^2 / w^2)
          f_geom(t) = η(t) / η(0)
      Depends on spot size σ ≡ w/2 in addition to d.

R(t) is by default taken to be the critical-surface radius r_crit(t) returned
by ``helios_postprocess.laser_intensity.analyze_laser_intensity``.  Alternative
choices (``--r-source``):

  • r_crit   (default) — physically the right "edge" the laser sees
  • r_outer            — outermost zone radius (corona edge)
  • r_ablation         — ablation-front radius (steepest negative ∂ρ/∂r outside
                         the hot spot); falls back to r_crit if not computable

Usage
-----
    python3 pdd_geometric_coupling.py <base_path> [options]

where ``<base_path>`` is the simulation path *without* extension, e.g.
``~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab02_foot25_s016/Olson_PDD_20_fab02_foot25_s016``.

Outputs are written alongside the .exo:

    <base>_geomcorr_power.csv        Corrected power table
    <base>_geomcorr_diagnostic.pdf   Diagnostic figure

This is a *diagnostic* tool.  It does not modify the .rhw automatically — the
corrected CSV is intended to be reviewed and then optionally injected into a
follow-up run as a starting point for empirical re-tuning.  Do not interpret
the corrected pulse as a physically calibrated input until it has been
benchmarked against HYDRA or LILAC absorbed-energy histories for the same
target.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from helios_postprocess import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.laser_intensity import analyze_laser_intensity

logger = logging.getLogger("geom_coupling")


# numpy 2.x renamed trapz → trapezoid; helios_postprocess uses both forms
# across modules.  Use whichever is available.
_TRAPZ = getattr(np, 'trapezoid', None) or np.trapz


# ─────────────────────────────────────────────────────────────────────────────
# Core formulas
# ─────────────────────────────────────────────────────────────────────────────

def view_factor_ratio(R: np.ndarray, d: float, theta: float) -> np.ndarray:
    """
    Option 2 — solid-angle view factor of a sphere of radius R(t) at distance
    d from a focal point, for a beam of cone half-angle θ.

    Three regimes:

      • R/d ≥ 1                — focal point inside the corona; the beam must
                                  traverse the corona to reach its focus, so
                                  the corona fully captures the beam.  f = 1.
      • sin(θ) ≤ R/d < 1       — sphere subtends a cone larger than (or equal
                                  to) the beam cone; the beam is fully
                                  captured by the corona.  f = 1.
      • R/d < sin(θ)           — sphere subtends a cone smaller than the beam
                                  cone; only the central portion of the beam
                                  hits the sphere.  f equals the ratio of
                                  the two solid angles:

                                  f = [1 - cos(arcsin(R/d))] / [1 - cos(θ)]
                                    = [1 - √(1 - (R/d)²)]    / [1 - cos(θ)]

    Returns
    -------
    f : ndarray same shape as R
        Coupling fraction in [0, 1].  NaN where R is NaN or non-positive.
    """
    R = np.asarray(R, dtype=float)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    omega_cone = 1.0 - cos_theta  # > 0 for any θ > 0

    with np.errstate(invalid='ignore'):
        # Clamp R/d to [0, 1] for the partial-cap branch; out-of-range values
        # are handled by the where() selectors below.
        ratio_clamped = np.minimum(np.abs(R) / d, 1.0)
        cos_arcsin    = np.sqrt(1.0 - ratio_clamped ** 2)
        f_partial     = (1.0 - cos_arcsin) / omega_cone

        # Build the result: f = 1 in fully-captured regimes, f_partial otherwise
        ratio = R / d
        fully_captured = (ratio >= sin_theta)        # also catches ratio ≥ 1
        f = np.where(fully_captured, 1.0, f_partial)

    # Mask invalid R
    valid = np.isfinite(R) & (R > 0)
    return np.where(valid, f, np.nan)


def gaussian_disk_ratio(R: np.ndarray, w: float) -> np.ndarray:
    """
    Option 1 — Gaussian-beam fraction intercepted by a disk of radius R(t).

        η(t) = 1 - exp(-2 R(t)^2 / w^2)

    where w is the 1/e^2 beam radius in the focal plane.  This is an
    approximation valid in the limit d ≫ R (the projected disk radius in the
    focal plane is approximately R itself).

    Returned values are not normalized; the caller divides by the t=0 value
    to obtain the time-dependent correction.
    """
    with np.errstate(invalid='ignore', divide='ignore'):
        return 1.0 - np.exp(-2.0 * (R / w) ** 2)


# ─────────────────────────────────────────────────────────────────────────────
# R(t) extraction
# ─────────────────────────────────────────────────────────────────────────────

def get_R_of_t(data, source: str, wavelength_um: float
               ) -> Tuple[np.ndarray, str]:
    """
    Return R(t) and a human-readable label describing the choice.

    ``source`` is one of ``r_crit``, ``r_outer``, ``r_ablation``.

    For r_ablation we use steepest negative density gradient outside the hot
    spot; falls back to r_crit if density is not available.
    """
    t = data.time

    if source == 'r_outer':
        R = data.zone_boundaries[:, -1]
        return R, "outer grid radius"

    # r_crit via the existing laser-intensity analysis
    if source in ('r_crit', 'r_ablation'):
        diag = analyze_laser_intensity(data, wavelength_um=wavelength_um)
        if diag is None:
            logger.warning(
                "analyze_laser_intensity returned None — falling back to "
                "outer-zone radius for R(t)."
            )
            return data.zone_boundaries[:, -1], "outer grid radius (fallback)"

        if source == 'r_crit':
            R = diag['r_crit']
            R = _backfill_nans(R, t, data.zone_boundaries[:, -1])
            return R, f"critical-surface radius (λ={wavelength_um:.3f} µm)"

        # r_ablation: steepest negative ∂ρ/∂r outside the hot spot
        rho = getattr(data, 'mass_density', None)
        if rho is None:
            logger.warning(
                "r_ablation requested but mass_density not available — "
                "using r_crit instead."
            )
            return diag['r_crit'], "critical-surface radius (r_ablation unavailable)"
        R = _ablation_front_radius(rho, data.zone_boundaries)
        return R, "ablation-front radius (steepest -∂ρ/∂r)"

    raise ValueError(f"Unknown --r-source: {source}")


def _backfill_nans(R: np.ndarray, t: np.ndarray, R_fallback: np.ndarray) -> np.ndarray:
    """Fill NaNs in R(t) by linear interpolation, with fallback to R_fallback."""
    R = R.copy()
    bad = ~np.isfinite(R)
    if not bad.any():
        return R
    good = ~bad
    if good.sum() < 2:
        return np.where(bad, R_fallback, R)
    R[bad] = np.interp(t[bad], t[good], R[good])
    return R


def _ablation_front_radius(rho: np.ndarray, zone_bounds: np.ndarray) -> np.ndarray:
    """
    Locate the ablation front as the zone with the steepest negative density
    gradient outside the inner 30% of zones (avoids hot-spot boundary
    artefacts).  Returns the outer boundary of that zone.
    """
    nt, nz = rho.shape
    R = np.full(nt, np.nan)
    inner_cut = nz // 3
    for ti in range(nt):
        # Negative of ∂ρ/∂r along the radial direction
        grad = -np.diff(rho[ti])
        idx_candidates = np.arange(inner_cut, nz - 1)
        if idx_candidates.size == 0:
            continue
        best = idx_candidates[np.argmax(grad[idx_candidates])]
        R[ti] = zone_bounds[ti, best + 1]
    return R


# ─────────────────────────────────────────────────────────────────────────────
# Anchor & normalize
# ─────────────────────────────────────────────────────────────────────────────

def normalize_f_geom(eta: np.ndarray,
                     P: np.ndarray
                     ) -> Tuple[np.ndarray, int]:
    """
    Given absolute coupling fractions ``eta(t)`` and original power ``P(t)``,
    return a normalized correction ``f_geom(t) = eta(t) / eta(anchor)``
    clipped to ``[0, 1]``.

    The anchor is the **maximum** of ``eta`` over the laser-on window.  This is
    the moment of strongest geometric coupling — typically when the corona is
    at peak expansion — and reflects the implicit assumption of the
    uncorrected 1D ray-trace that the laser is fully captured by the spherical
    target.  At later times when the absorbing surface shrinks below this
    maximum, f_geom < 1 and the corrected power is scaled down accordingly.

    Returns
    -------
    f_geom : ndarray
        Normalized correction in ``[0, 1]`` with the anchor at 1.
        Values where ``eta`` is NaN are filled with 1 (no correction).
    anchor_idx : int
        Index of the anchor sample (for diagnostic plotting).
    """
    f_geom = np.ones_like(eta, dtype=float)
    on = P > 0
    if not on.any():
        return f_geom, 0

    valid = on & np.isfinite(eta) & (eta > 0)
    if not valid.any():
        logger.warning("No laser-on sample has finite eta(t) — returning "
                       "uncorrected pulse (f_geom ≡ 1).")
        return f_geom, int(np.argmax(on))

    # Anchor at the MAX eta among laser-on samples — i.e. when the geometric
    # coupling is at its strongest.  Anchoring at the first laser-on sample
    # is unsafe because the corona is just forming there and eta is at a local
    # *minimum*, which would make f_geom = eta/anchor exceed 1 and get clipped
    # at all later times — masking the entire correction.
    eta_masked = np.where(valid, eta, -np.inf)
    anchor_idx = int(np.argmax(eta_masked))
    anchor_eta = float(eta[anchor_idx])

    with np.errstate(invalid='ignore'):
        f_geom = np.where(np.isfinite(eta), eta / anchor_eta, 1.0)
    f_geom = np.clip(f_geom, 0.0, 1.0)
    return f_geom, anchor_idx


# ─────────────────────────────────────────────────────────────────────────────
# Diagnostic plot
# ─────────────────────────────────────────────────────────────────────────────

def make_diagnostic_pdf(out_path: Path,
                        t_ns: np.ndarray,
                        R_cm: np.ndarray,
                        R_label: str,
                        d_cm: float,
                        theta_deg: float,
                        f_geom: np.ndarray,
                        P_orig_TW: np.ndarray,
                        P_corr_TW: np.ndarray,
                        formula: str,
                        base_name: str) -> None:
    """Three-panel diagnostic figure: R(t) / α_c(t), f_geom(t), P(t)."""
    alpha_c = np.degrees(np.arcsin(np.clip(R_cm / d_cm, 0.0, 1.0)))

    with PdfPages(out_path) as pdf:
        fig, axes = plt.subplots(3, 1, figsize=(8.5, 10), sharex=True,
                                  gridspec_kw={'hspace': 0.18})

        # ── Panel 1: R(t) and α_c(t)
        ax = axes[0]
        ax.plot(t_ns, R_cm * 1e4, lw=1.8, color='#1A2841',
                label=f'R(t) — {R_label}')
        ax.set_ylabel('R(t)  [µm]', color='#1A2841')
        ax.tick_params(axis='y', labelcolor='#1A2841')
        ax2 = ax.twinx()
        ax2.plot(t_ns, alpha_c, lw=1.5, color='#B5732C', ls='--',
                 label=fr'$\alpha_c(t) = \arcsin(R/d)$,  d={d_cm:.3f} cm')
        ax2.axhline(theta_deg, color='gray', ls=':', lw=1.0,
                    label=fr'beam cone half-angle $\theta = {theta_deg:.1f}^\circ$')
        ax2.set_ylabel(r'$\alpha_c(t)$  [deg]', color='#B5732C')
        ax2.tick_params(axis='y', labelcolor='#B5732C')
        ax.set_title(f"Geometric coupling correction — {base_name}",
                     fontsize=12, weight='bold')
        ax.grid(alpha=0.3)
        # combined legend
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc='upper right', fontsize=9, framealpha=0.9)

        # ── Panel 2: f_geom(t)
        ax = axes[1]
        ax.plot(t_ns, f_geom, lw=2.0, color='#1E6A7E',
                label=f'f_geom(t)  ({formula})')
        ax.axhline(1.0, color='gray', ls=':', lw=0.8)
        ax.set_ylabel('Geometric coupling f_geom(t)')
        ax.set_ylim(0, 1.1)
        ax.grid(alpha=0.3)
        ax.legend(loc='lower left', fontsize=10)
        # Annotate the late-time decay
        if np.any(np.isfinite(f_geom)):
            f_end = np.nanmin(f_geom)
            t_end_idx = int(np.nanargmin(f_geom))
            ax.annotate(f'min f_geom = {f_end:.3f}\n@ t = {t_ns[t_end_idx]:.2f} ns',
                        xy=(t_ns[t_end_idx], f_end),
                        xytext=(0.6, 0.55), textcoords='axes fraction',
                        fontsize=10, ha='left',
                        arrowprops=dict(arrowstyle='->', color='#666',
                                         connectionstyle='arc3,rad=0.2'))

        # ── Panel 3: P_original and P_corrected
        ax = axes[2]
        ax.plot(t_ns, P_orig_TW, lw=1.8, color='#1A2841',
                label='P original')
        ax.plot(t_ns, P_corr_TW, lw=2.0, color='#C44A30',
                label='P corrected (geom-multiplier applied)')
        ax.set_ylabel('Laser power  [TW]')
        ax.set_xlabel('Time  [ns]')
        ax.grid(alpha=0.3)
        ax.legend(loc='upper left', fontsize=10)

        # Integrated effect annotation
        # (trapezoidal integration in TW·ns = kJ)
        E_orig = _TRAPZ(P_orig_TW * 1e12, t_ns * 1e-9) * 1e-3        # kJ
        E_corr = _TRAPZ(P_corr_TW * 1e12, t_ns * 1e-9) * 1e-3
        dE_frac = (E_corr / E_orig - 1.0) if E_orig > 0 else 0.0
        ax.text(0.98, 0.95,
                f"∫P_orig dt  = {E_orig:7.1f} kJ\n"
                f"∫P_corr dt = {E_corr:7.1f} kJ\n"
                f"Δ           = {100*dE_frac:+6.1f} %",
                transform=ax.transAxes, ha='right', va='top',
                fontsize=10, family='monospace',
                bbox=dict(boxstyle='round', facecolor='#F4F6F8',
                          edgecolor='#C8CFD6'))

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Geometric coupling pre-processor for Helios 1D direct-drive runs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('base_path',
                        help='Simulation base path (no extension).  '
                             'A no-burn run is strongly recommended so the '
                             'R(t) trajectory reflects pure hydrodynamic '
                             'implosion, free of alpha-driven late-time motion.')
    parser.add_argument('--r-source', choices=('r_crit', 'r_outer', 'r_ablation'),
                        default='r_crit',
                        help='Which radius to use for R(t)')
    parser.add_argument('--formula', choices=('view_factor', 'gaussian'),
                        default='view_factor',
                        help='Correction formula.  view_factor (Option 2) is '
                             'beam-independent; gaussian (Option 1) requires '
                             'the beam spot size.')
    parser.add_argument('--d', type=float, default=None,
                        help='Focal distance d from origin [cm] '
                             '(default: from beam 1 in .rhw).')
    parser.add_argument('--theta-deg', type=float, default=None,
                        help='Beam cone half-angle θ [deg] '
                             '(default: from beam 1 in .rhw).')
    parser.add_argument('--w-cm', type=float, default=None,
                        help='Gaussian 1/e² beam radius w [cm] '
                             '(default: laser_spot_size_cm from beam 1 in .rhw, '
                             'which Helios defines as the 1/e² radius).')
    parser.add_argument('--wavelength-um', type=float, default=None,
                        help='Laser wavelength [µm] for n_crit '
                             '(default: from .rhw).')
    parser.add_argument('--out-prefix', type=str, default=None,
                        help='Output filename prefix (default: <base>_geomcorr).')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(format="%(name)s [%(levelname)s]: %(message)s",
                        level=log_level)

    base_path = Path(args.base_path).expanduser().resolve()
    exo_path  = base_path.with_suffix('.exo')
    if not exo_path.exists():
        logger.error(f"EXO not found at {exo_path}")
        return 1

    logger.info(f"Loading {exo_path.name}")
    run = HeliosRun(str(exo_path), verbose=False)

    # Load rhw_config explicitly — build_run_data does not parse it
    # automatically, so we must hand it in.  (See examples/energy_balance_diagnostic.py
    # for the same pattern.)
    rhw_path = exo_path.with_suffix('.rhw')
    rhw_config = None
    if rhw_path.exists():
        try:
            from helios_postprocess.rhw_parser import load_rhw_configuration
            rhw_config = load_rhw_configuration(str(rhw_path))
        except Exception as exc:
            logger.warning(f"Could not parse rhw: {exc}")
    else:
        logger.warning(f"No .rhw found at {rhw_path} — geometry must come from CLI flags.")

    data = build_run_data(run, time_unit='s', rhw_config=rhw_config, verbose=False)

    rhw = getattr(data, 'rhw_config', None) or rhw_config
    if rhw is None and (args.d is None or args.theta_deg is None
                         or args.wavelength_um is None):
        logger.error("rhw_config is None and required geometry flags not "
                     "given.  Pass --d, --theta-deg, --w-cm, --wavelength-um manually.")
        return 1

    # ── Resolve geometric parameters
    # Prefer per-beam list (beam 1 by default) if present; fall back to the
    # top-level "laser_*" attributes; fall back to CLI defaults.
    per_beam = getattr(rhw, 'laser_geometry_per_beam', None) if rhw else None
    beam1    = per_beam[0] if per_beam else None

    def _from_rhw(per_beam_key, top_attr, fallback):
        if beam1 is not None and per_beam_key in beam1:
            return beam1[per_beam_key]
        if rhw is not None and hasattr(rhw, top_attr):
            v = getattr(rhw, top_attr)
            if v is not None and v != 0:
                return v
        return fallback

    d_cm        = args.d            if args.d            is not None else _from_rhw('focus_position',   'laser_focus_position_cm',    0.22)
    theta_deg   = args.theta_deg    if args.theta_deg    is not None else _from_rhw('half_cone_angle',  'laser_half_cone_angle_deg',  1.0)
    spot_cm     = _from_rhw('spot_size', 'laser_spot_size_cm', 0.25)
    # In Helios .rhw, "spot size" is the 1/e² Gaussian beam radius (w).
    # σ = w/2 for the conventional Gaussian-RMS interpretation, but the
    # Option-1 formula η = 1 - exp(-2 R²/w²) uses w directly.
    w_cm        = args.w_cm         if args.w_cm         is not None else spot_cm
    wavelength  = args.wavelength_um if args.wavelength_um is not None else _from_rhw('wavelength', 'laser_wavelength_um', 0.351)
    theta = np.radians(theta_deg)

    logger.info(f"Geometry:  d = {d_cm:.4f} cm,  θ = {theta_deg:.2f}°,  "
                f"σ (spot) = {spot_cm:.4f} cm,  w = {w_cm:.4f} cm,  "
                f"λ = {wavelength:.3f} µm")

    # ── R(t)
    R_cm, R_label = get_R_of_t(data, args.r_source, wavelength_um=wavelength)
    if not np.any(np.isfinite(R_cm)):
        logger.error(f"R(t) ({args.r_source}) is all NaN — cannot proceed.")
        return 1
    # data.time is in nanoseconds regardless of `time_unit` (which describes
    # the unit of the *incoming* raw exo times — see data_builder.py).
    t_ns = data.time
    t_s  = t_ns * 1e-9

    # ── f_geom(t) — absolute coupling fraction eta(t), then normalize
    if args.formula == 'view_factor':
        eta = view_factor_ratio(R_cm, d_cm, theta)
    else:  # gaussian
        eta = gaussian_disk_ratio(R_cm, w_cm)

    # ── Laser power on target (sum over beams)
    P_W = getattr(data, 'laser_power_on_target', None)
    if P_W is None:
        logger.error("data.laser_power_on_target is not present — cannot compute "
                     "corrected power table.")
        return 1
    P_orig_TW = np.asarray(P_W, dtype=float) * 1e-12

    f_geom, anchor_idx = normalize_f_geom(eta, P_orig_TW)
    P_corr_TW = P_orig_TW * f_geom

    logger.info(f"Anchor: t = {t_ns[anchor_idx]:.3f} ns "
                f"(R = {R_cm[anchor_idx]*1e4:.1f} µm, "
                f"η_anchor = {eta[anchor_idx]:.4f})")

    # ── Outputs
    prefix    = args.out_prefix or base_path.name + "_geomcorr"
    csv_path  = base_path.parent / f"{prefix}_power.csv"
    pdf_path  = base_path.parent / f"{prefix}_diagnostic.pdf"

    # CSV
    with open(csv_path, 'w') as fh:
        fh.write("# Geometric-coupling-corrected laser power table\n")
        fh.write(f"# Source run:     {base_path.name}.exo\n")
        fh.write(f"# R(t) source:    {R_label}\n")
        fh.write(f"# Formula:        {args.formula}\n")
        fh.write(f"# Geometry:       d = {d_cm:.4f} cm, "
                 f"θ = {theta_deg:.2f}°, w = {w_cm:.4f} cm, "
                 f"λ = {wavelength:.3f} µm\n")
        fh.write("# Columns:        time[ns], R[cm], f_geom, P_orig[TW], P_corr[TW]\n")
        for i in range(len(t_ns)):
            fh.write(f"{t_ns[i]:11.5f}  {R_cm[i]:11.5e}  "
                     f"{f_geom[i] if np.isfinite(f_geom[i]) else float('nan'):8.4f}  "
                     f"{P_orig_TW[i]:11.4f}  {P_corr_TW[i]:11.4f}\n")
    logger.info(f"Wrote {csv_path}")

    # PDF
    make_diagnostic_pdf(
        pdf_path, t_ns, R_cm, R_label, d_cm, theta_deg,
        f_geom, P_orig_TW, P_corr_TW,
        args.formula, base_path.name,
    )
    logger.info(f"Wrote {pdf_path}")

    # ── Summary
    E_orig_kJ = float(_TRAPZ(P_orig_TW * 1e12, t_s)) * 1e-3
    E_corr_kJ = float(_TRAPZ(P_corr_TW * 1e12, t_s)) * 1e-3
    f_min     = float(np.nanmin(f_geom)) if np.any(np.isfinite(f_geom)) else float('nan')

    # f_geom at peak laser power — the most calibration-relevant single number
    on = P_orig_TW > 0
    if on.any():
        t_peak_idx = int(np.argmax(P_orig_TW))
        f_at_peak  = float(f_geom[t_peak_idx])
        t_peak_ns  = float(t_ns[t_peak_idx])
        # Statistics restricted to laser-on window
        f_geom_on  = f_geom[on]
        f_min_on   = float(np.nanmin(f_geom_on))
        f_mean_on  = float(np.nanmean(f_geom_on))
    else:
        f_at_peak = float('nan'); t_peak_ns = float('nan')
        f_min_on = float('nan');  f_mean_on = float('nan')

    print()
    print("─" * 62)
    print(f"  Geometric coupling pre-processor — {base_path.name}")
    print("─" * 62)
    print(f"  R(t) source:           {R_label}")
    print(f"  Correction formula:    {args.formula}")
    print(f"  d (focal distance):    {d_cm:.4f} cm")
    print(f"  θ (cone half-angle):   {theta_deg:.2f}°")
    if args.formula == 'gaussian':
        print(f"  w (1/e² beam radius):  {w_cm:.4f} cm")
    print(f"  Wavelength λ:          {wavelength:.3f} µm")
    print(f"  Sphere-fills-cone at:  R ≥ d·sin(θ) = {d_cm*np.sin(theta)*1e4:.1f} µm")
    print()
    print(f"  --- During laser-on window ---")
    print(f"  f_geom at peak power:  {f_at_peak:.3f}    (t = {t_peak_ns:.2f} ns)")
    print(f"  mean f_geom (on):      {f_mean_on:.3f}")
    print(f"  min  f_geom (on):      {f_min_on:.3f}")
    print(f"  min  f_geom (any t):   {f_min:.3f}")
    print()
    print(f"  ∫P_original dt:        {E_orig_kJ:8.2f} kJ")
    print(f"  ∫P_corrected dt:       {E_corr_kJ:8.2f} kJ")
    if E_orig_kJ > 0:
        print(f"  Δ energy:              {100*(E_corr_kJ/E_orig_kJ - 1.0):+6.2f} %")

    # Diagnostic warning when correction is essentially nil during laser-on
    if np.isfinite(f_min_on) and f_min_on > 0.95 and args.r_source == 'r_crit':
        print()
        print("  ⚠  Correction is essentially nil during the laser-on window.")
        print("     R(t) stays above d·sin(θ) throughout — the critical surface")
        print("     (corona) fills the beam cone the whole time the laser is on.")
        print("     This is the expected behaviour when R(t) = r_crit and the")
        print("     corona is well-expanded.  The geometric-miss problem really")
        print("     manifests at the *imploding shell*, not the critical surface.")
        print()
        print("     To diagnose the shell-tracking miss instead, try:")
        print(f"       --r-source r_ablation     (steepest -∂ρ/∂r outside HS)")
        print(f"     or override R(0) and the focal distance with --d to test")
        print(f"     a tighter geometry.")
    print("─" * 62)
    print(f"  Outputs:")
    print(f"    {csv_path}")
    print(f"    {pdf_path}")
    print("─" * 62)
    print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
