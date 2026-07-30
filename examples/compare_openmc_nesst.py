#!/usr/bin/env python
"""
compare_openmc_nesst.py  --  RUN IN THE openmc-env.

Load the .npz from dump_snapshot.py, run OpenMC (neutron + photon), and:
  * triangulate rhoR_D: OpenMC / NeSST(center-to-edge) / Method-2(traversed)
  * DSR (OpenMC escaping vs NeSST)
  * absolute-yield unfold: test Y_total = Y(13-15) * e^{4*DSR} (Abu-Shawareb)
  * escaping-gamma diagnostic (carbon 4.4 MeV / H 2.2 MeV lines)

Usage:
    python compare_openmc_nesst.py n210808_snapshot.npz \
        [--particles 400000] [--batches 25] [--dr-min 2e-4] [--prefix cmp]
"""
import argparse
import numpy as np
import helios_to_openmc as hto

_trap = getattr(np, "trapezoid", None) or np.trapz


def _norm_primary(e_MeV, y, lo=13.0, hi=15.0):
    e = np.asarray(e_MeV, float); y = np.asarray(y, float)
    m = (e >= lo) & (e <= hi)
    a = _trap(y[m], e[m]) if m.any() else 0.0
    return y / (a if a > 0 else 1.0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("npz")
    ap.add_argument("--particles", type=int, default=400_000)
    ap.add_argument("--batches", type=int, default=25)
    ap.add_argument("--dr-min", type=float, default=2.0e-4)
    ap.add_argument("--prefix", default="openmc_nesst_cmp")
    ap.add_argument("--no-photons", action="store_true")
    args = ap.parse_args(argv)

    z = np.load(args.npz, allow_pickle=True)
    raw = hto.ProfileSnapshot.from_npz(args.npz)
    snap = hto.resample_min_thickness(raw, dr_min_cm=args.dr_min)
    analytic_rhoR_D = hto.analytic_rhoR_species(snap, "D")
    analytic_rhoR_C = hto.analytic_rhoR_species(snap, "C")

    edges = np.concatenate([[0.0], snap.r_cm]); prod = hto._per_shell_production(snap)
    r_mean = float(np.sum(prod * 0.5 * (edges[:-1] + edges[1:])) / prod.sum())
    print(f"[cmp] geometry: {len(raw.r_cm)} raw zones -> {len(snap.r_cm)} shells "
          f"(dr_min={args.dr_min:g} cm); R0={snap.r_cm[-1]:.4f} cm; mean birth r={r_mean:.4e} cm")

    model = hto.build_model(snap, energy="birth", particles=args.particles,
                            batches=args.batches, photon_transport=not args.no_photons)
    sp = model.run()
    r = hto.read_results(sp)

    nesst_rhoR_D = float(z["nesst_rhoR_D"]) if "nesst_rhoR_D" in z.files else float("nan")
    nesst_dsr_fuel = float(z["nesst_DSR_fuel"])
    nesst_dsr_total = float(z["nesst_DSR_total"]) if "nesst_DSR_total" in z.files else float("nan")
    nesst_E = z["nesst_energy_MeV"]; nesst_full = z["nesst_full"]

    g = lambda x: "  n/a " if not np.isfinite(x) else f"{x:6.3f}"
    p = lambda x: "  n/a " if not np.isfinite(x) else f"{100*x:6.2f}%"

    print("\n" + "=" * 70)
    print("  OpenMC  vs  NeSST  vs  analytic Method-2   (same Helios snapshot)")
    print("=" * 70)
    print(f"  snapshot   {snap.label}")
    print(f"  sigma_el_D (from OpenMC data): {r['sigma_el_D_barn']:.3f} b")
    print("-" * 70)
    print(f"  {'rhoR_D (g/cm^2)':24s}{'OpenMC':>10s}{'NeSST':>10s}{'Method-2':>12s}")
    print(f"  {'':24s}{g(r['rhoR_D_gcm2']):>10s}{g(nesst_rhoR_D):>10s}{g(analytic_rhoR_D):>12s}")
    print(f"  {'  (NeSST=center-to-edge; Method-2/OpenMC=traversed)':50s}")
    print("-" * 70)
    print(f"  {'DSR (10-12/13-15)':24s}{p(r['DSR']):>10s}{p(nesst_dsr_total):>10s}")
    print(f"  {'  NeSST fuel-only':24s}{'':>10s}{p(nesst_dsr_fuel):>10s}")
    print("-" * 70)

    # absolute-yield unfold: Y_total = Y(13-15) * e^{4*DSR}  (Abu-Shawareb NIF)
    N1315 = r["N_13_15_per_source"]; dsr = r["DSR"]
    ratio_meas = (1.0 / N1315) if N1315 > 0 else float("nan")   # Y_total(born)/Y(13-15)
    ratio_paper = float(np.exp(4.0 * dsr)) if np.isfinite(dsr) else float("nan")
    k_eff = (np.log(ratio_meas) / dsr) if (dsr > 0 and np.isfinite(ratio_meas)) else float("nan")
    print("  absolute-yield unfold   Y_total = Y(13-15) * e^(4*DSR)")
    print(f"    escaping in 13-15 MeV / source born   {N1315:.4f}")
    print(f"    Y_total / Y(13-15)  OpenMC (=1/N)      {ratio_meas:.4f}")
    print(f"    Y_total / Y(13-15)  paper e^(4*DSR)    {ratio_paper:.4f}")
    print(f"    effective exponent k (our transport)  {k_eff:.2f}    (paper uses 4)")
    print("    [caveat: source is DT+DD+TT; DD/TT (~1%) never enter 13-15, so")
    print("     1/N slightly over-counts total vs a pure-DT unfold.]")
    print("-" * 70)
    print(f"  neutron leakage / source n     {r['neutron_leak_per_source']:.3f}")
    if not args.no_photons:
        print(f"  escaping gammas / source n     {r['gamma_per_source_neutron']:.4f}")
        print(f"    carbon 4.4 MeV line          {r['gamma_line_4p4_MeV']:.3e}")
        print(f"    H 2.2 MeV capture line       {r['gamma_line_2p2_MeV']:.3e}")
    print("=" * 70)
    print("  NeSST rhoR = center-to-edge int_0^R rho dr (Kyle Method-1). Method-2 &")
    print("  OpenMC = emission-weighted TRAVERSED rhoR (Kyle Method-2), which is what")
    print("  Abu-Shawareb ties DSR to ('areal density emission weighted to larger radii').")
    print("=" * 70 + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yo = _norm_primary(r["neutron_MeV"], r["neutron_dNdE"])
    yn = _norm_primary(nesst_E, nesst_full)
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    ax.axvspan(10, 12, color="#cccccc", alpha=0.4)
    ax.semilogy(r["neutron_MeV"], yo, "k", lw=2.0, label="OpenMC escaping (transport)")
    ax.semilogy(nesst_E, yn, "--", color="#1f4e79", lw=1.8, label="NeSST single-scatter")
    ax.set_xlim(1, 16); ax.set_ylim(1e-4, 5)
    ax.set_xlabel("neutron energy (MeV)"); ax.set_ylabel("dN/dE (norm. to 13-15 MeV)")
    ax.set_title(f"{snap.label}: escaping neutron spectrum")
    ax.legend(fontsize=9); ax.grid(alpha=0.25, which="both")
    fig.tight_layout(); fig.savefig(f"{args.prefix}_neutron.png", dpi=160)
    print(f"[cmp] wrote {args.prefix}_neutron.png")

    if not args.no_photons and r["gamma_per_source_neutron"] > 0:
        fig, ax = plt.subplots(figsize=(8.4, 5.0))
        ax.semilogy(r["gamma_MeV"], r["gamma_dNdE"], color="#7030a0", lw=1.6)
        for e0, lab in ((2.223, "H(n,g) 2.2"), (4.44, "12C(n,n'g) 4.4")):
            ax.axvline(e0, color="#c0392b", ls=":", lw=1.2)
            ax.text(e0 + 0.1, ax.get_ylim()[1] * 0.3, lab, rotation=90, va="top",
                    fontsize=8, color="#c0392b")
        ax.set_xlim(0, 11); ax.set_xlabel("photon energy (MeV)")
        ax.set_ylabel("dN/dE (photons / source n / MeV, 4pi)")
        ax.set_title(f"{snap.label}: escaping gamma spectrum "
                     f"({r['gamma_per_source_neutron']:.3f} g / source n)")
        ax.grid(alpha=0.25, which="both")
        fig.tight_layout(); fig.savefig(f"{args.prefix}_gamma.png", dpi=160)
        print(f"[cmp] wrote {args.prefix}_gamma.png")


if __name__ == "__main__":
    main()
