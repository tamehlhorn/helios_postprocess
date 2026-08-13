#!/usr/bin/env python3
"""
Synthetic x-ray streak camera from a Helios run -- Level 0 forward model.

Produces, for one run:
  page 1  imaging streak (space vs sweep time) in three filtered channels
  page 2  spectral streak (photon energy vs sweep time) + snapshot spectra
  page 3  extracted observables: emission-edge trajectory, emission history,
          channel ratio (Te proxy), and the optical-depth diagnostic

Usage
-----
    python3 examples/xray_streak_synthetic.py <base_path> [--window T0 T1]
                                              [--thin] [--capsule-only]

``base_path`` is the run directory containing the .exo, or the .exo itself.

Author: Prof T / helios_postprocess
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from helios_postprocess.xray.response import _DEFAULT_CATHODE
from helios_postprocess.xray import (StreakConfig, build_radiance_cube,
                                     make_imaging_streak, make_spectral_streak,
                                     default_channels, channels_from_edges,
                                     burn_metrics,
                                     emission_edge_trajectory, emission_history,
                                     opacity_report)

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("xray_streak")


def _find_exo(base_path: Path) -> Path:
    if base_path.suffix == ".exo":
        return base_path
    hits = sorted(base_path.glob("*.exo"))
    if not hits:
        raise FileNotFoundError(f"No .exo found under {base_path}")
    return hits[0]


def _norm(a):
    m = np.nanmax(a)
    return a / m if m > 0 else a


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def make_report(cube, imgs, spec, stem, pdf_path):
    """Write the three-page Level 0 streak report."""
    with PdfPages(pdf_path) as pdf:
        # ---------------- page 1: imaging streaks ----------------
        fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
        for ax, (name, img) in zip(np.atleast_1d(axes), imgs.items()):
            S = _norm(img.signal)
            S = np.where(S > 1e-6, S, 1e-6)
            ax.pcolormesh(img.t_ns, img.axis * 1e4, S.T,
                          norm=LogNorm(vmin=1e-5, vmax=1.0),
                          shading="auto", cmap="inferno")
            ax.plot(img.t_ns, emission_edge_trajectory(img) * 1e4,
                    "c-", lw=1.0, label="50% edge")
            ax.set_title(f"{name} channel")
            ax.set_xlabel("sweep time (ns)")
            ax.legend(loc="upper right", fontsize=8)
        np.atleast_1d(axes)[0].set_ylabel("impact parameter (µm)")
        fig.suptitle(f"{stem} — Level 0 imaging streak (self-emission)")
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ---------------- page 2: spectral streak ----------------
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        S = _norm(spec.signal)
        S = np.where(S > 1e-8, S, 1e-8)
        axes[0].pcolormesh(spec.t_ns, spec.axis * 1e-3, S.T,
                           norm=LogNorm(vmin=1e-6, vmax=1.0),
                           shading="auto", cmap="viridis")
        axes[0].set_yscale("log")
        axes[0].set_xlabel("sweep time (ns)")
        axes[0].set_ylabel("photon energy (keV)")
        axes[0].set_title("spectral streak (spatially integrated)")

        hist = emission_history(spec)
        ipk = int(np.argmax(hist))
        for frac, style in ((0.6, "--"), (1.0, "-"), (1.4, ":")):
            j = int(np.clip(ipk * frac, 0, spec.t_ns.size - 1))
            axes[1].loglog(spec.axis * 1e-3, np.maximum(spec.signal[j], 1e-30),
                           style, label=f"t = {spec.t_ns[j]:.3f} ns")
        axes[1].set_xlabel("photon energy (keV)")
        axes[1].set_ylabel("dP/dE (erg s$^{-1}$ eV$^{-1}$)")
        axes[1].set_title("snapshot spectra")
        axes[1].legend(fontsize=8)
        axes[1].grid(alpha=0.3, which="both")
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # ---------------- page 3: observables ----------------
        fig, axes = plt.subplots(2, 2, figsize=(13, 9))

        ax = axes[0, 0]
        for name, img in imgs.items():
            ax.plot(img.t_ns, emission_edge_trajectory(img) * 1e4, label=name)
        ax.set_xlabel("time (ns)"); ax.set_ylabel("emission edge (µm)")
        ax.set_title("emission-edge trajectory (ablation front, NOT cold shell)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axes[0, 1]
        for name, img in imgs.items():
            ax.semilogy(img.t_ns,
                        np.maximum(_norm(emission_history(img)), 1e-8),
                        label=name)
        ax.set_xlabel("time (ns)"); ax.set_ylabel("normalized emission")
        ax.set_title("x-ray emission history")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        ax = axes[1, 0]
        keys = list(imgs)
        h_soft = emission_history(imgs[keys[0]])
        h_hard = emission_history(imgs[keys[-1]])
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(h_soft > 0, h_hard / h_soft, np.nan)
        ax.semilogy(imgs[keys[0]].t_ns, ratio)
        ax.set_xlabel("time (ns)"); ax.set_ylabel(f"{keys[-1]} / {keys[0]}")
        ax.set_title("channel ratio (Te proxy)")
        ax.grid(alpha=0.3)

        ax = axes[1, 1]
        ax.axis("off")
        lines = ["METRICS (Level 0)", ""]
        for name, img in imgs.items():
            m = burn_metrics(img)
            lines.append(f"{name:>5s}: bang {m.bang_time_ns:.4f} ns   "
                         f"FWHM {m.fwhm_ps:.1f} ps   "
                         f"flash R {m.flash_radius_cm * 1e4:.1f} um")
        ms = burn_metrics(spec)
        lines += ["", f"spectral: bang {ms.bang_time_ns:.4f} ns   "
                      f"FWHM {ms.fwhm_ps:.1f} ps", ""]
        rep = opacity_report(cube)
        if rep.get("valid", 0.0) > 0:
            lines += ["OPTICAL DEPTH DIAGNOSTIC",
                      f"  tau_max = {rep['tau_max']:.3g} at "
                      f"{rep['tau_max_energy_eV']:.0f} eV, "
                      f"t = {rep['tau_max_time_ns']:.3f} ns, "
                      f"p = {rep['tau_max_impact_cm'] * 1e4:.0f} um",
                      f"  volume fraction with tau > 0.3: "
                      f"{rep['frac_volume_above_warn'] * 100:.1f}%",
                      f"  min I/I_thin = {rep['min_I_over_I_thin']:.3g}"]
        else:
            lines += ["OPTICAL DEPTH DIAGNOSTIC", "  run was optically thin"]
        floor = cube.meta.get("rho_floor_g_cc")
        lines += ["", f"rho floor: {'none' if floor is None else f'{floor:g} g/cc'}"]
        lines += ["", "Level 0 = free-free continuum only.",
                  "No free-bound edges, no lines, no non-LTE.",
                  "Filter response is the ANCHORED POWER LAW model -",
                  "load Henke tables before absolute comparison."]
        ax.text(0.0, 1.0, "\n".join(lines), va="top", family="monospace",
                fontsize=9)

        fig.suptitle(f"{stem} — extracted observables")
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base_path")
    ap.add_argument("--window", nargs=2, type=float, default=None,
                    metavar=("T0_NS", "T1_NS"))
    ap.add_argument("--thin", action="store_true",
                    help="optically thin; default includes LTE continuum "
                         "absorption so the tau diagnostic is meaningful")
    ap.add_argument("--capsule-only", action="store_true",
                    help="halfraum targets: exclude external regions")
    ap.add_argument("--zone-stop", type=int, default=None,
                    help="manual capsule outer zone index, overriding region-"
                         "name detection")
    ap.add_argument("--n-time", type=int, default=400)
    ap.add_argument("--n-impact", type=int, default=128)
    ap.add_argument("--n-energy", type=int, default=96)
    ap.add_argument("--time-res-ps", type=float, default=15.0)
    ap.add_argument("--time-unit", choices=["auto", "s", "ns"], default="auto",
                    help="unit of time_whole in the .exo (Helios writes 's')")
    ap.add_argument("--bands", nargs="+", default=None,
                    metavar="E_LO,E_HI",
                    help="channel bands in eV, e.g. --bands 2000,4000 "
                         "4000,8000 8000,20000. Default channels are cut for "
                         "a softer OMEGA-class spectrum; check with "
                         "examples/xray_emission_provenance.py first.")
    ap.add_argument("--e-range", nargs=2, type=float, default=None,
                    metavar=("E_MIN_EV", "E_MAX_EV"),
                    help="photon energy grid range (default 200-15000 eV, "
                         "auto-expanded to cover --bands)")
    ap.add_argument("--no-cathode", action="store_true",
                    help="omit photocathode QE from the channel response; "
                         "use for pure band studies above ~10 keV where the "
                         "QE model rolls off and hides the signal")
    ap.add_argument("--r-max-um", type=float, default=None,
                    help="outer extent of the impact-parameter grid. Default "
                         "is bounded by where the emission is, not by the "
                         "mesh outer boundary.")
    ap.add_argument("--rho-floor", type=float, default=None,
                    metavar="G_PER_CC",
                    help="treat zones below this density as non-emitting; "
                         "suppresses the massless very-hot outer corona. "
                         "ALWAYS run with and without and compare.")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    from helios_postprocess.core import HeliosRun
    from helios_postprocess.data_builder import build_run_data

    base = Path(args.base_path).expanduser()
    exo = _find_exo(base)
    outdir = Path(args.out) if args.out else exo.parent
    outdir.mkdir(parents=True, exist_ok=True)

    log.info(f"[xray] reading {exo}")
    run = HeliosRun(str(exo))
    # Helios EXODUS stores time_whole in SECONDS; build_run_data only converts
    # when told to.  Guard rather than trust the flag, because a silent
    # seconds/ns mix-up would put the whole sweep window in the wrong decade.
    t_unit = args.time_unit
    if t_unit == "auto":
        t_unit = "s" if float(np.max(run.times)) < 1.0e-3 else "ns"
        log.info(f"[xray] time unit detected: {t_unit}")
    data = build_run_data(run, time_unit=t_unit, verbose=False)

    zslice = None
    if args.zone_stop is not None:
        zslice = slice(0, int(args.zone_stop))
        log.info(f"[xray] manual zone slice 0:{args.zone_stop}")
    elif args.capsule_only:
        ncr = getattr(data, "n_capsule_regions", None)
        rii = getattr(data, "region_interfaces_indices", None)
        if ncr is not None and rii is not None:
            stop = int(np.asarray(rii)[0, ncr - 1])
            zslice = slice(0, stop)
            log.info(f"[xray] capsule-only: zones 0:{stop}")
        else:
            log.warning("[xray] --capsule-only requested but region info absent")

    band_edges = None
    if args.bands:
        band_edges = [tuple(float(x) for x in b.split(",")) for b in args.bands]

    e_min, e_max = (args.e_range if args.e_range else (200.0, 15000.0))
    if band_edges:
        lo = min(b[0] for b in band_edges)
        hi = max(b[1] for b in band_edges)
        # A band outside the grid would silently integrate to zero, so widen
        # the grid to cover it rather than returning an empty channel.
        if lo < e_min or hi > e_max:
            new_min, new_max = min(e_min, lo * 0.9), max(e_max, hi * 1.05)
            log.warning(f"[xray] widening energy grid "
                        f"{e_min:.0f}-{e_max:.0f} -> "
                        f"{new_min:.0f}-{new_max:.0f} eV to cover --bands")
            e_min, e_max = new_min, new_max

    t0, t1 = (args.window if args.window else (None, None))
    cfg = StreakConfig(t_start_ns=t0, t_stop_ns=t1,
                       E_min_eV=e_min, E_max_eV=e_max,
                       n_time=args.n_time, n_impact=args.n_impact,
                       n_energy=args.n_energy,
                       time_resolution_ps=args.time_res_ps,
                       optically_thin=args.thin,
                       zone_slice=zslice,
                       rho_floor_g_cc=args.rho_floor,
                       r_max_cm=(args.r_max_um * 1e-4
                                 if args.r_max_um else None))

    log.info("[xray] building radiance cube ...")
    cube = build_radiance_cube(data, cfg)

    if band_edges:
        channels = channels_from_edges(
            band_edges, cathode=None if args.no_cathode else _DEFAULT_CATHODE)
        log.info(f"[xray] channels: {list(channels)}"
                 f"{' (no cathode)' if args.no_cathode else ''}")
    else:
        channels = default_channels()
    imgs = {k: make_imaging_streak(cube, ch, cfg) for k, ch in channels.items()}
    spec = make_spectral_streak(cube, cfg)

    stem = exo.stem
    pdf_path = outdir / f"{stem}_xray_streak_level0.pdf"
    make_report(cube, imgs, spec, stem, pdf_path)
    log.info(f"[xray] wrote {pdf_path}")

    npz_path = outdir / f"{stem}_xray_cube_level0.npz"
    np.savez_compressed(npz_path, t_ns=cube.t_ns, p_cm=cube.p_cm,
                        E_eV=cube.E_eV, I=cube.I, tau=cube.tau,
                        I_thin=cube.I_thin)
    log.info(f"[xray] wrote {npz_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
