#!/usr/bin/env python3
"""
X-ray emission provenance: which zones actually produce the signal?

A streak camera integrates along a chord through everything, so a synthetic
image is only as trustworthy as the emissivity of the zones that dominate it.
Two failure modes are common and both are invisible in the final image:

  1. **Outer-boundary zones.** The outermost Lagrangian zones of an exploding
     corona carry almost no mass but can reach very high Te. Free-free
     emissivity goes as Zbar^2 n_e n_i exp(-E/kTe), so a nearly massless zone
     at tens of keV contributes essentially nothing at low photon energy but
     can manufacture a hard tail that dominates the high-energy percentiles.
     If the hard channel of a synthetic streak is really a picture of the
     last three zones, the diagnostic is measuring the mesh, not the physics.

  2. **Wrong band.** If the filter/channel definitions sit where the target
     does not emit, the synthetic image is dominated by response tails.

This script answers both by decomposing the band-integrated emitted power
zone by zone and aggregating it by region, and by re-computing the spectrum
with the outermost zones trimmed.

Usage
-----
    python3 examples/xray_emission_provenance.py <exo_or_dir>
            [--time NS] [--trim N] [--band E_LO E_HI] [--zone-stop N]

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
    return hits[0]


def _hdr(t: str) -> None:
    print("\n" + "=" * 74)
    print(f"  {t}")
    print("=" * 74)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("base_path")
    ap.add_argument("--time", type=float, default=None,
                    help="ns; default is the peak-emission proxy")
    ap.add_argument("--trim", type=int, default=5,
                    help="number of outermost zones to trim in the "
                         "boundary-artifact check")
    ap.add_argument("--band", nargs=2, type=float, default=[200.0, 30000.0],
                    metavar=("E_LO", "E_HI"))
    ap.add_argument("--zone-stop", type=int, default=None)
    args = ap.parse_args()

    from helios_postprocess.core import HeliosRun
    from helios_postprocess.data_builder import build_run_data
    from helios_postprocess.xray import build_emissivity, resolve_Te

    exo = _find_exo(Path(args.base_path).expanduser())
    run = HeliosRun(str(exo))
    t_unit = "s" if float(np.max(run.times)) < 1.0e-3 else "ns"
    data = build_run_data(run, time_unit=t_unit, verbose=False)

    t = np.asarray(data.time, float)
    rho = np.asarray(data.mass_density, float)
    rb = np.asarray(data.zone_boundaries, float)
    Te = resolve_Te(data)
    vol = (4.0 / 3.0) * np.pi * (rb[:, 1:] ** 3 - rb[:, :-1] ** 3)

    zslice = slice(0, int(args.zone_stop)) if args.zone_stop else slice(None)

    if args.time is None:
        proxy = np.sum((rho[:, zslice] ** 2)
                       * np.sqrt(np.maximum(Te[:, zslice], 1e-3))
                       * vol[:, zslice], axis=1)
        it = int(np.argmax(proxy))
    else:
        it = int(np.argmin(np.abs(t - args.time)))

    print(f"file : {exo.name}")
    print(f"time : {t[it]:.4f} ns  (index {it})")
    print(f"band : {args.band[0]:.0f} - {args.band[1]:.0f} eV")

    E = np.logspace(np.log10(args.band[0]), np.log10(args.band[1]), 128)
    cube = build_emissivity(data, it, E, zone_slice=zslice,
                            include_absorption=False)

    # Optically-thin emitted power per zone: 4*pi * integral(j_E dE) * V
    P_zone = 4.0 * np.pi * np.trapezoid(cube.j_E, E, axis=1) * vol[it, zslice]
    P_tot = P_zone.sum()

    # ------------------------------------------------------------ by region
    _hdr("EMITTED POWER BY REGION (optically thin)")
    names = getattr(data, "region_names", None)
    rii = getattr(data, "region_interfaces_indices", None)
    zmass = getattr(data, "zone_mass", None)

    if names and rii is not None:
        ri = np.asarray(rii)[it].astype(int)
        print(f" {'region':22s} {'zones':>12s} {'% power':>9s} "
              f"{'% mass':>8s} {'<Te> eV':>10s} {'<rho>':>10s}")
        prev = 0
        m_tot = (zmass[it, zslice].sum() if zmass is not None else np.nan)
        for i, b in enumerate(ri):
            lo, hi = prev, int(b)
            prev = hi
            sub = slice(max(lo, 0), min(hi, P_zone.size))
            if sub.start >= sub.stop:
                continue
            pw = P_zone[sub].sum()
            mm = (zmass[it, sub].sum() / m_tot * 100.0
                  if zmass is not None and m_tot > 0 else np.nan)
            w = P_zone[sub]
            Tw = (np.sum(Te[it, sub] * w) / pw) if pw > 0 else np.nan
            print(f" {names[i][:22]:22s} {lo:5d}-{hi - 1:5d} "
                  f"{pw / P_tot * 100:8.2f}% {mm:7.2f}% "
                  f"{Tw:10.1f} {np.mean(rho[it, sub]):10.3g}")
    else:
        print(" (no region metadata)")

    # -------------------------------------------------- outer-boundary check
    _hdr(f"OUTER-BOUNDARY CHECK (last {args.trim} zones)")
    n = P_zone.size
    k = max(0, n - int(args.trim))
    frac_outer = P_zone[k:].sum() / P_tot * 100.0
    m_all = zmass[it, zslice] if zmass is not None else None
    frac_mass = (m_all[k:].sum() / m_all.sum() * 100.0
                 if m_all is not None and m_all.sum() > 0 else np.nan)
    print(f" zones {k}-{n - 1} carry {frac_outer:.3f}% of the "
          f"band-integrated power and {frac_mass:.3f}% of the mass")
    print(f"\n {'zone':>6s} {'r (um)':>10s} {'rho':>11s} {'Te (eV)':>10s} "
          f"{'Zbar':>7s} {'% power':>9s}")
    zc = 0.5 * (rb[it, :-1] + rb[it, 1:])[zslice]
    zbar = np.asarray(getattr(data, "mean_charge", np.ones_like(rho))[it])[zslice]
    for j in range(max(0, n - 8), n):
        print(f" {j:6d} {zc[j] * 1e4:10.1f} {rho[it, zslice][j]:11.3e} "
              f"{Te[it, zslice][j]:10.1f} {zbar[j]:7.2f} "
              f"{P_zone[j] / P_tot * 100:8.3f}%")

    # ------------------------------------------------ spectrum with/without
    _hdr("SPECTRUM SENSITIVITY TO THE OUTERMOST ZONES")
    spec_full = 4.0 * np.pi * np.sum(cube.j_E * vol[it, zslice][:, None], axis=0)
    spec_trim = 4.0 * np.pi * np.sum(cube.j_E[:k] * vol[it, zslice][:k, None],
                                     axis=0)

    def pct(spec):
        if spec.max() <= 0:
            return (np.nan,) * 3
        c = np.cumsum(spec * np.gradient(E))
        c = c / c[-1]
        return (float(np.interp(0.10, c, E)), float(np.interp(0.50, c, E)),
                float(np.interp(0.90, c, E)))

    f10, f50, f90 = pct(spec_full)
    t10, t50, t90 = pct(spec_trim)
    print(f" full      : 10/50/90 = {f10:7.0f} / {f50:7.0f} / {f90:7.0f} eV   "
          f"P = {np.trapezoid(spec_full, E):.4g} erg/s")
    print(f" trimmed   : 10/50/90 = {t10:7.0f} / {t50:7.0f} / {t90:7.0f} eV   "
          f"P = {np.trapezoid(spec_trim, E):.4g} erg/s")

    shift = abs(f90 - t90) / max(f90, 1.0) * 100.0
    print(f" 90th-percentile shift on trimming: {shift:.2f}%")

    # The artifact signature is BRIGHT AND MASSLESS: a large share of the
    # emitted power carried by a negligible share of the mass. Outer zones
    # contributing in proportion to their mass are simply the ablator, and
    # trimming them would be throwing away real signal.
    massless_and_bright = (np.isfinite(frac_mass) and frac_outer > 5.0
                           and frac_outer > 20.0 * max(frac_mass, 1e-6))
    print()
    if massless_and_bright or shift > 10.0:
        print(" -> OUTER-ZONE ARTIFACT SUSPECTED. These zones are far brighter")
        print("    than their mass warrants, or trimming them moves the hard")
        print("    end of the spectrum appreciably. The hard channel would be")
        print("    partly a picture of the mesh boundary.")
        print(f"    Re-run with --zone-stop {k} and compare.")
    elif frac_outer > 5.0:
        print(" -> Outer zones contribute, but roughly in proportion to their")
        print("    mass, and trimming barely moves the spectrum. That is the")
        print("    ablator emitting, not a boundary artifact. Do not trim.")
    else:
        print(" -> Outer zones are not driving the spectrum. The band is set")
        print("    by the bulk of the capsule.")

    # ------------------------------------------------------- band suggestion
    _hdr("BAND SUGGESTION")
    use_trim = massless_and_bright or shift > 10.0
    lo, mid, hi = (t10, t50, t90) if use_trim else (f10, f50, f90)
    e1, e2, e3, e4 = lo * 0.8, mid * 0.7, mid * 1.6, hi * 1.3
    print(f" 90% of the emitted power lies between {lo:.0f} and {hi:.0f} eV.")
    print(" Three channels spanning that range:")
    print(f"   --bands {e1:.0f},{e2:.0f} {e2:.0f},{e3:.0f} {e3:.0f},{e4:.0f}")
    print("\n The default channels (0.8-2 / 2-4 / >4 keV) are cut for a")
    print(" softer, OMEGA-class direct-drive target. Using them on a harder")
    print(" spectrum puts most of the signal in one channel and makes the")
    print(" hard/soft ratio a poor Te proxy.")
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
