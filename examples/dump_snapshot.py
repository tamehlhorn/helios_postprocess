#!/usr/bin/env python
"""
dump_snapshot.py  --  RUN IN YOUR HELIOS ANALYSIS ENV (helios_postprocess + NeSST).

Parse a Helios run into a portable shell profile + the NeSST single-scatter
reference, and write one .npz for compare_openmc_nesst.py (openmc-env).

Usage
-----
    python dump_snapshot.py <run_dir_or_exo> [--rhw FILE.rhw] \
        [--time-index N] [--n-E 800] [--out n210808_snapshot.npz]
"""
import argparse
import numpy as np

import helios_to_openmc as hto


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--rhw", default=None)
    ap.add_argument("--time-index", type=int, default=None)
    ap.add_argument("--n-E", type=int, default=800)
    ap.add_argument("--out", default="n210808_snapshot.npz")
    args = ap.parse_args(argv)

    snap, nd, data, rhw = hto.from_helios_snapshot(
        args.run, rhw_path=args.rhw, time_index=args.time_index)
    print(f"[dump] {snap.label}: {len(snap.r_cm)} shells, "
          f"R0 = {snap.r_cm[-1]:.4f} cm, peak rho = {snap.rho_gcc.max():.1f} g/cm^3")
    if rhw:
        print(f"[dump] composition from {rhw}")

    from helios_postprocess import neutron_scatter as nsc
    res = nsc._composition_scatter(nd, data, str(rhw), args.n_E) if rhw else None
    extra = {}
    if res is not None:
        extra.update(nesst_mode="multi-material",
                     nesst_rhoR_D=res["rhoR_D_gcm2"], nesst_rhoR_T=res["rhoR_T_gcm2"],
                     nesst_rhoR_C=res["rhoR_C_gcm2"],
                     nesst_DSR_fuel=res["dsr"]["DSR"], nesst_DSR_total=res["dsr_total"]["DSR"],
                     nesst_scattered=res["scattered"])
    else:
        res = nsc.scatter_from_extraction(nd, n_E=args.n_E)
        extra.update(nesst_mode="fuel-only", nesst_rhoR_D=float("nan"),
                     nesst_DSR_fuel=res["dsr"]["DSR"], nesst_scattered=res["scattered"])
    extra.update(nesst_frac_D=res.get("frac_D", 0.5), nesst_frac_T=res.get("frac_T", 0.5),
                 nesst_yield=res.get("yield", float("nan")),
                 nesst_energy_MeV=res["energy_MeV"], nesst_full=res["full"])

    snap.to_npz(args.out, **extra)
    print(f"[dump] wrote {args.out}   (mode: {extra['nesst_mode']})")
    line = f"[dump] NeSST  DSR_fuel = {100*float(extra['nesst_DSR_fuel']):.2f}%"
    if "nesst_DSR_total" in extra:
        line += (f"   DSR_total = {100*float(extra['nesst_DSR_total']):.2f}%"
                 f"   rhoR_D = {float(extra['nesst_rhoR_D']):.3f}"
                 f"  rhoR_C = {float(extra['nesst_rhoR_C']):.3f} g/cm^2")
    print(line)


if __name__ == "__main__":
    main()
