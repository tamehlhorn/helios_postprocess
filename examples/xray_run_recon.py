#!/usr/bin/env python3
"""
X-ray reconnaissance on a Helios run.

Answers, in one pass and before any expensive forward model is run:

  * Is this a bare capsule or a halfraum? Where does the capsule end?
    (For an indirect-drive target the converter shell outshines the capsule
    by orders of magnitude -- get the zone slice wrong and the synthetic
    streak is a picture of the halfraum.)
  * Are the fields the free-free model actually needs present --
    mean_charge, elec_temperature, ion_density, electron_density?
  * What sweep window makes sense? Reports peak-emission-proxy time,
    stagnation proxy, and how badly the EXODUS time base is clustered.
  * What photon energy range carries the signal, and where is the
    optically-thin assumption at risk?

Usage
-----
    python3 examples/xray_run_recon.py <exo_or_dir> [--capsule-only]

Prints a recommended command line for xray_streak_synthetic.py.

Author: Prof T / helios_postprocess
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

logging.basicConfig(level=logging.WARNING, format="%(message)s")


def _find_exo(base: Path) -> Path:
    if base.suffix == ".exo":
        return base
    hits = sorted(base.glob("*.exo"))
    if not hits:
        raise FileNotFoundError(f"No .exo found under {base}")
    if len(hits) > 1:
        print(f"note: {len(hits)} .exo files present, using {hits[0].name}")
    return hits[0]


def _hdr(title: str) -> None:
    print("\n" + "=" * 72)
    print(f"  {title}")
    print("=" * 72)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("base_path")
    ap.add_argument("--capsule-only", action="store_true",
                    help="evaluate the emission diagnostics on the capsule "
                         "zones only")
    ap.add_argument("--zone-stop", type=int, default=None,
                    help="manual capsule outer zone index, overriding region-"
                         "name detection (use when region names do not decode)")
    args = ap.parse_args()

    from helios_postprocess.core import HeliosRun
    from helios_postprocess.data_builder import build_run_data
    from helios_postprocess.xray import (build_emissivity, make_impact_grid,
                                         integrate_formal, resolve_zbar,
                                         resolve_Te)

    exo = _find_exo(Path(args.base_path).expanduser())
    print(f"file: {exo}")

    run = HeliosRun(str(exo))
    t_unit = "s" if float(np.max(run.times)) < 1.0e-3 else "ns"
    data = build_run_data(run, time_unit=t_unit, verbose=False)
    t = np.asarray(data.time, float)

    # ---------------------------------------------------------------- fields
    _hdr("FIELDS REQUIRED BY THE LEVEL 0 MODEL")
    ds = run.dataset
    checks = [
        ("mass_density", data.mass_density, "required"),
        ("elec_temperature", getattr(data, "elec_temperature", None),
         "required -- coronal emission is set by Te, NOT Ti"),
        ("electron_density", getattr(data, "electron_density", None), "required"),
        ("ion_density", getattr(data, "ion_density", None),
         "needed; falls back to n_e/Zbar"),
        ("mean_charge", getattr(data, "mean_charge", None),
         "j_ff scales as Zbar^2 -- run examples/patch_xray_variables.py"),
    ]
    for name, val, note in checks:
        in_exo = name in ds.variables
        state = "loaded" if val is not None else ("IN .exo, NOT LOADED"
                                                  if in_exo else "ABSENT")
        flag = "  " if val is not None else ("!!" if in_exo else "??")
        print(f" {flag} {name:20s} {state:22s} {note}")

    # ---------------------------------------------------------------- regions
    _hdr("TARGET STRUCTURE")
    print(f" target_class       : {getattr(data, 'target_class', 'unknown')}")
    names = getattr(data, "region_names", None)
    if names and not all(any(c.isalpha() for c in n) for n in names):
        print(" !! region names did not decode as text -- halfraum detection")
        print("    cannot be trusted. Use --zone-stop to set the capsule")
        print("    outer zone index by hand.")
        names = None
    rii = getattr(data, "region_interfaces_indices", None)
    n_cap = getattr(data, "n_capsule_regions", None)
    zslice = None
    if names and rii is not None:
        ri = np.asarray(rii)[0].astype(int)
        prev = 0
        for i, b in enumerate(ri):
            tag = ""
            if n_cap is not None:
                tag = "  <- capsule" if i < n_cap else "  <- EXTERNAL"
            print(f"   {i + 1}. {names[i]:24s} zones {prev:4d}-{int(b) - 1:4d} "
                  f"({int(b) - prev:4d}){tag}")
            prev = int(b)
        if n_cap is not None and n_cap < len(ri):
            stop = int(ri[n_cap - 1])
            zslice = slice(0, stop)
            print(f"\n HALFRAUM DETECTED. Capsule is zones 0:{stop}.")
            print(" Run the streak model with --capsule-only, or the converter")
            print(" shell will dominate every channel.")
    else:
        print(" (no region metadata in this file)")
    if args.zone_stop is not None:
        zslice = slice(0, int(args.zone_stop))
        print(f"\n manual override: capsule zones 0:{args.zone_stop}")
    elif not args.capsule_only:
        zslice = None

    # ---------------------------------------------------------------- timing
    _hdr("TIME BASE")
    dt = np.diff(t)
    print(f" steps              : {t.size}")
    print(f" range              : {t[0]:.4f} - {t[-1]:.4f} ns")
    print(f" dt  min/median/max : {dt.min() * 1e3:.3f} / "
          f"{np.median(dt) * 1e3:.3f} / {dt.max() * 1e3:.3f} ps")
    print(f" clustering ratio   : {dt.max() / dt.min():.1f}x  "
          f"(EXODUS concentrates steps near stagnation; the streak model "
          f"resamples onto a uniform sweep base)")

    # emission proxy: rho^2 sqrt(Te) volume integral (bremsstrahlung-like)
    rho = np.asarray(data.mass_density, float)
    Te = resolve_Te(data)
    rb = np.asarray(data.zone_boundaries, float)
    vol = (4.0 / 3.0) * np.pi * (rb[:, 1:] ** 3 - rb[:, :-1] ** 3)
    sl = zslice if zslice is not None else slice(None)
    proxy = np.sum((rho[:, sl] ** 2) * np.sqrt(np.maximum(Te[:, sl], 1e-3))
                   * vol[:, sl], axis=1)
    i_em = int(np.argmax(proxy))

    # Peak-density radius is a far better convergence proxy than the grid
    # outer boundary, which for a halfraum is a static wall.
    r_peak = np.array([0.5 * (rb[i, :-1] + rb[i, 1:])[sl][
        int(np.argmax(rho[i, sl]))] for i in range(t.size)])
    i_min = int(np.argmin(r_peak))

    print(f"\n peak emission proxy: t = {t[i_em]:.4f} ns  (index {i_em})")
    print(f" min peak-rho radius: t = {t[i_min]:.4f} ns  (index {i_min}), "
          f"r = {r_peak[i_min] * 1e4:.1f} um")

    w0 = max(t[0], t[i_em] - 1.0)
    w1 = min(t[-1], t[i_em] + 0.5)
    print(f" suggested window   : --window {w0:.2f} {w1:.2f}")

    # ------------------------------------------------------- spectrum + tau
    _hdr("SPECTRUM AND OPTICAL DEPTH AT PEAK EMISSION")
    E = np.logspace(np.log10(200.0), np.log10(20000.0), 64)
    cube = build_emissivity(data, i_em, E, zone_slice=zslice,
                            include_absorption=True)
    p = make_impact_grid(cube.r_bnd[-1], 48)
    res = integrate_formal(cube, p)

    w = 2.0 * np.pi * p
    spec = 4.0 * np.pi * np.trapezoid(res.I * w[:, None], p, axis=0)
    if spec.max() > 0:
        cum = np.cumsum(spec * np.gradient(E))
        cum /= cum[-1]
        e10 = float(np.interp(0.10, cum, E))
        e50 = float(np.interp(0.50, cum, E))
        e90 = float(np.interp(0.90, cum, E))
        print(f" emitted power 10/50/90 percentile : "
              f"{e10:.0f} / {e50:.0f} / {e90:.0f} eV")
    tau = res.tau
    print(f" chord tau  max                    : {tau.max():.3g}")
    if tau.max() > 0:
        i, j = np.unravel_index(int(np.argmax(tau)), tau.shape)
        print(f"            at                     : p = {p[i] * 1e4:.0f} um, "
              f"E = {E[j]:.0f} eV")
    print(f" fraction of (p,E) with tau > 0.3  : "
          f"{np.mean(tau > 0.3) * 100:.1f}%")
    if np.mean(tau > 0.3) > 0.05:
        print(" -> self-absorption is NOT negligible at this time. Run the")
        print("    forward model WITHOUT --thin, and expect SPECT3D Level 1")
        print("    to move the answer.")
    else:
        print(" -> optically thin over most of the band at this time.")

    zbar = resolve_zbar(data, i_em)
    print(f"\n Zbar range at peak emission       : "
          f"{np.min(zbar[sl]):.2f} - {np.max(zbar[sl]):.2f}")
    print(f" Te range at peak emission (eV)    : "
          f"{np.min(Te[i_em, sl]):.1f} - {np.max(Te[i_em, sl]):.1f}")

    # ---------------------------------------------------------------- next
    _hdr("SUGGESTED NEXT COMMAND")
    if args.zone_stop is not None:
        flag = f" --zone-stop {args.zone_stop}"
    elif zslice is not None or (n_cap is not None and names
                                and n_cap < len(names)):
        flag = " --capsule-only"
    else:
        flag = ""
    print(f" python3 examples/xray_streak_synthetic.py \\\n"
          f"     {exo} \\\n"
          f"     --window {w0:.2f} {w1:.2f}{flag}")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
