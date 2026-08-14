#!/usr/bin/env python3
"""
Summarize a saved Level 0 radiance cube.

Reads the ``*_xray_cube_level0.npz`` written by xray_streak_synthetic.py and
reports the observables numerically, so band choices can be iterated in
seconds instead of re-running the forward model.

Reports per channel:
  * x-ray bang time, emission FWHM, 10-90 rise
  * stagnation flash radius (50% of peak on the radial profile)
  * emission-edge trajectory at a few times

And across channels:
  * bang-time spread -- if the channels disagree, they are looking at
    different things, which on an optically thick target is expected rather
    than a bug
  * escape fraction vs photon energy, power-weighted, from I / I_thin

Usage
-----
    python3 examples/xray_cube_summary.py <cube.npz>
            [--bands E_LO,E_HI ...] [--time-res-ps PS] [--n-time N]

Author: Prof T / helios_postprocess
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from helios_postprocess.xray import (RadianceCube, StreakConfig,
                                     make_imaging_streak, make_spectral_streak,
                                     channels_from_edges, burn_metrics,
                                     emission_edge_trajectory,
                                     encircled_radius, radial_profile,
                                     opacity_report)


def _hdr(t: str) -> None:
    print("\n" + "=" * 74)
    print(f"  {t}")
    print("=" * 74)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--bands", nargs="+", default=["2000,4000", "5000,9000",
                                                   "10000,20000"],
                    metavar="E_LO,E_HI")
    ap.add_argument("--time-res-ps", type=float, default=15.0)
    ap.add_argument("--n-time", type=int, default=600)
    ap.add_argument("--profiles", nargs="*", type=float, default=None,
                    metavar="T_NS",
                    help="dump normalized radial profiles at these times "
                         "(no values = 5 spread across the window)")
    ap.add_argument("--cathode", action="store_true",
                    help="include photocathode QE (default: bare bands)")
    args = ap.parse_args()

    d = np.load(Path(args.npz).expanduser())
    cube = RadianceCube(t_ns=d["t_ns"], p_cm=d["p_cm"], E_eV=d["E_eV"],
                        I=d["I"], tau=d["tau"], I_thin=d["I_thin"], meta={})

    print(f"cube  : {Path(args.npz).name}")
    print(f"shape : {cube.I.shape} (t, p, E)")
    print(f"time  : {cube.t_ns[0]:.4f} - {cube.t_ns[-1]:.4f} ns "
          f"({cube.t_ns.size} hydro steps)")
    print(f"energy: {cube.E_eV[0]:.0f} - {cube.E_eV[-1]:.0f} eV")
    print(f"impact: 0 - {cube.p_cm[-1] * 1e4:.1f} um")

    edges = [tuple(float(x) for x in b.split(",")) for b in args.bands]
    if args.cathode:
        channels = channels_from_edges(edges)
    else:
        channels = channels_from_edges(edges, cathode=None)

    cfg = StreakConfig(n_time=args.n_time, time_resolution_ps=args.time_res_ps)

    # ------------------------------------------------------------ per channel
    _hdr("PER-CHANNEL OBSERVABLES")
    print(f" {'channel':>14s} {'bang (ns)':>11s} {'FWHM (ps)':>11s} "
          f"{'rise (ps)':>11s} {'flash R (um)':>13s}")
    imgs, metrics = {}, {}
    for name, ch in channels.items():
        img = make_imaging_streak(cube, ch, cfg)
        m = burn_metrics(img)
        imgs[name], metrics[name] = img, m
        print(f" {name:>14s} {m.bang_time_ns:11.4f} {m.fwhm_ps:11.1f} "
              f"{m.rise_10_90_ps:11.1f} {m.flash_radius_cm * 1e4:13.1f}")

    spec = make_spectral_streak(cube, cfg)
    ms = burn_metrics(spec)
    print(f" {'spectral':>14s} {ms.bang_time_ns:11.4f} {ms.fwhm_ps:11.1f} "
          f"{ms.rise_10_90_ps:11.1f} {'-':>13s}")

    bangs = np.array([m.bang_time_ns for m in metrics.values()])
    spread = (bangs.max() - bangs.min()) * 1e3
    print(f"\n bang-time spread across channels: {spread:.1f} ps")
    if spread > 2.0 * args.time_res_ps:
        print(" -> The channels do NOT peak together. On an optically thick")
        print("    target that is expected: each band samples a different")
        print("    depth, so 'the' x-ray bang time is band-dependent and must")
        print("    be quoted with its channel.")
    else:
        print(" -> Channels peak together within the instrument resolution.")

    # ------------------------------------------------------ edge trajectories
    _hdr("EMISSION-EDGE TRAJECTORY (um)")
    probe_t = np.linspace(cube.t_ns[0], cube.t_ns[-1], 9)
    print(f" {'t (ns)':>9s}" + "".join(f"{n:>14s}" for n in imgs))
    traj = {n: emission_edge_trajectory(img) for n, img in imgs.items()}
    for tp in probe_t:
        row = f" {tp:9.4f}"
        for n, img in imgs.items():
            j = int(np.argmin(np.abs(img.t_ns - tp)))
            v = traj[n][j]
            row += f"{v * 1e4:14.1f}" if np.isfinite(v) else f"{'--':>14s}"
        print(row)

    names = list(imgs)
    if len(names) >= 2:
        a, b = traj[names[0]], traj[names[-1]]
        # Quote the offset AT BANG TIME, not averaged over the window: early
        # in the sweep the profile is an extended corona rather than a disk,
        # and averaging across that regime change produces a number that
        # describes neither.
        img0 = imgs[names[0]]
        jb = int(np.argmin(np.abs(img0.t_ns - metrics[names[0]].bang_time_ns)))
        if np.isfinite(a[jb]) and np.isfinite(b[jb]):
            print(f"\n at bang time, {names[0]} minus {names[-1]}: "
                  f"{(a[jb] - b[jb]) * 1e4:+.1f} um")
            print(" -> A positive offset means the opaque channel's edge sits")
            print("    OUTSIDE the transparent channel's: the tau = 1 surface")
            print("    standing off the emitting fuel.")

    # -------------------------------------------------------- hydro overlay
    if "r_peak_rho_cm" in d:
        _hdr("EMISSION EDGE vs HYDRO RADII (um)")
        print(" An emission edge is not a material surface. An emissivity")
        print(" contour can sweep outward through stationary material as")
        print(" easily as the material can move, so the edge must be read")
        print(" against Lagrangian radii before it is called a trajectory.\n")
        r_pk = d["r_peak_rho_cm"]
        r_reg = d["r_region_cm"] if "r_region_cm" in d else None
        rnames = ([str(x) for x in d["region_names"]]
                  if "region_names" in d else None)
        hdr = f" {'t (ns)':>9s} {'peak rho':>10s}"
        if r_reg is not None:
            for k in range(r_reg.shape[1]):
                lbl = (rnames[k][:9] if rnames and k < len(rnames)
                       else f"reg{k}")
                hdr += f"{lbl:>11s}"
        hdr += f"{'edge ' + names[0][:6]:>13s}"
        print(hdr)
        for tp in probe_t:
            i = int(np.argmin(np.abs(cube.t_ns - tp)))
            row = f" {tp:9.4f} {r_pk[i] * 1e4:10.1f}"
            if r_reg is not None:
                for k in range(r_reg.shape[1]):
                    row += f"{r_reg[i, k] * 1e4:11.1f}"
            j = int(np.argmin(np.abs(imgs[names[0]].t_ns - tp)))
            v = traj[names[0]][j]
            row += f"{v * 1e4:13.1f}" if np.isfinite(v) else f"{'--':>13s}"
            print(row)

    # ------------------------------------------- threshold vs encircled radius
    _hdr("RADIUS DEFINITION: 50% THRESHOLD vs ENCIRCLED ENERGY (um)")
    print(" A threshold edge takes the outermost 50%-of-peak crossing and is")
    print(" sensitive to shallow shoulders; encircled energy integrates the")
    print(" flux. Large disagreement means the emission is not a compact")
    print(" core, and any quoted radius must name its definition.\n")
    print(f" {'t (ns)':>9s}" + "".join(f"{n[:6]+' thr':>13s}{n[:6]+' R50':>13s}"
                                        f"{n[:6]+' R90':>13s}" for n in imgs))
    enc50 = {n: encircled_radius(img, 0.5) for n, img in imgs.items()}
    enc90 = {n: encircled_radius(img, 0.9) for n, img in imgs.items()}
    for tp in probe_t:
        row = f" {tp:9.4f}"
        for n, img in imgs.items():
            j = int(np.argmin(np.abs(img.t_ns - tp)))
            for arr in (traj[n], enc50[n], enc90[n]):
                v = arr[j]
                row += f"{v * 1e4:13.1f}" if np.isfinite(v) else f"{'--':>13s}"
        print(row)

    n0 = list(imgs)[0]
    good = np.isfinite(traj[n0]) & np.isfinite(enc50[n0])
    if good.any():
        ratio = np.nanmedian(traj[n0][good] / enc50[n0][good])
        print(f"\n median threshold-edge / R50 in {n0}: {ratio:.2f}x")
        if ratio > 2.0:
            print(" -> The threshold edge sits far outside the flux-weighted")
            print("    radius. The profile has a broad low-level halo, so the")
            print("    threshold trajectory is tracking that halo, not the")
            print("    emitting core. Prefer R50/R90 for trajectory work.")

    # ------------------------------------------------------------- profiles
    if args.profiles is not None:
        times = (args.profiles if len(args.profiles) > 0
                 else list(np.linspace(cube.t_ns[0], cube.t_ns[-1], 5)))
        _hdr("NORMALIZED RADIAL PROFILES")
        for tp in times:
            print(f"\n t = {tp:.4f} ns")
            for n, img in imgs.items():
                pax, prof = radial_profile(img, tp)
                keep = pax <= min(pax[-1], 3.0 * np.nanmax(
                    [enc90[n][int(np.argmin(np.abs(img.t_ns - tp)))], 1e-4]))
                idx = np.linspace(0, keep.sum() - 1, 12).astype(int)
                vals = " ".join(f"{prof[keep][i]:6.3f}" for i in idx)
                rads = " ".join(f"{pax[keep][i] * 1e4:6.0f}" for i in idx)
                print(f"   {n:>10s} r(um): {rads}")
                print(f"   {'':>10s} I/Imax: {vals}")

    # ----------------------------------------------------- escape fraction
    _hdr("ESCAPE FRACTION vs PHOTON ENERGY (power-weighted)")
    w = 2.0 * np.pi * cube.p_cm
    P_out = np.trapezoid(cube.I * w[None, :, None], cube.p_cm, axis=1)
    P_thin = np.trapezoid(cube.I_thin * w[None, :, None], cube.p_cm, axis=1)
    ipk = int(np.argmax(np.trapezoid(P_out, cube.E_eV, axis=1)))
    with np.errstate(divide="ignore", invalid="ignore"):
        f_esc = np.where(P_thin[ipk] > 0, P_out[ipk] / P_thin[ipk], np.nan)
    print(f" at peak emission, t = {cube.t_ns[ipk]:.4f} ns\n")
    for probe in (1000.0, 2000.0, 4000.0, 6000.0, 10000.0, 15000.0):
        if cube.E_eV[0] <= probe <= cube.E_eV[-1]:
            print(f"   {probe:7.0f} eV : {np.interp(probe, cube.E_eV, f_esc) * 100:6.2f}%")
    tot = np.trapezoid(P_out[ipk], cube.E_eV) / np.trapezoid(P_thin[ipk],
                                                             cube.E_eV)
    print(f"\n band-integrated: {tot * 100:.2f}%")

    rep = opacity_report(cube)
    if rep.get("valid", 0.0) > 0:
        print(f" tau_max = {rep['tau_max']:.3g} at "
              f"{rep['tau_max_energy_eV']:.0f} eV, "
              f"t = {rep['tau_max_time_ns']:.4f} ns")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
