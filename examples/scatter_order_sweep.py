#!/usr/bin/env python
"""
scatter_order_sweep.py  --  RUN IN THE openmc-env.

Isolate the effect of MULTIPLE scattering on the neutron diagnostic signature,
cleanly and in a single run, by partitioning OpenMC's escaping neutrons by
collision count (openmc.CollisionFilter):

  * single-scatter reference = neutrons that leak with <= 1 collision
        (uncollided primary + once-scattered) -- exactly what NeSST's
        single-scatter model computes.
  * full transport           = all leaking neutrons (multiple scattering,
        (n,2n), and the full down-scatter cascade included).

Same code, same ENDF/B-VIII.0 data, same geometry, one run -> the difference
between the two is multiple scattering and nothing else (no Method-1/2 or
code-to-code confound).

At each areal density it reports:
  - DSR(single) vs DSR(full) and their fractional difference
  - fraction of escaping neutrons that are multiply scattered
  - OpenMC collision-count rhoR_D (all collisions) vs the single-traversal
    analytic rhoR_D  -> the multiple-scatter enhancement of the readout
  - escaping spectra, single vs full

Geometry: an idealized UNIFORM DT (optionally + carbon) sphere with a CENTRAL
point source, so every source neutron traverses exactly rhoR = rho * R and the
single-scatter prediction is exact.  This is the textbook control for isolating
multiple scattering; it is not a Helios profile (use --npz to add the real one).

Usage:
    python scatter_order_sweep.py \
        --rhoR 0.1,0.3,0.6,1.0,1.5,2.0 \
        --particles 200000 --batches 20 --R-um 50 --carbon-frac 0.0 \
        --prefix scatter_sweep

    # sanity-check the analysis pipeline with NO OpenMC (mocked tallies):
    python scatter_order_sweep.py --self-test --prefix selftest

    # also run the real profile through the same collision partition:
    python scatter_order_sweep.py --rhoR 0.3,0.6,1.0 --npz n210808_snapshot.npz
"""
from __future__ import annotations
import argparse
import csv
import numpy as np
import helios_to_openmc as hto

_trap = getattr(np, "trapezoid", None) or np.trapz
EDGES = hto.NEUTRON_EDGES_MeV                       # 0..16 MeV, 0.1 MeV bins (161 edges)
NE = len(EDGES) - 1
# scatter cross-section scale for the optical-depth estimate used only in --self-test
_TAU_PER_GCM2 = 0.24


# --------------------------------------------------------------------------
# geometry / source
# --------------------------------------------------------------------------
def make_uniform_snap(rhoR_D, R_cm, carbon_frac=0.0, n_shell=40, label=None):
    """Uniform sphere whose deuterium areal density (central traversal) is rhoR_D.
    Mass fractions: carbon_frac is the C mass fraction; the rest is 50/50 D:T.
    Production S is placed at the center so analytic_rhoR_species() weights the
    central chord (rho*R)."""
    wC = float(carbon_frac)
    wD = 0.5 * (1.0 - wC)
    wT = 0.5 * (1.0 - wC)
    rho = rhoR_D / (wD * R_cm)                      # rho * wD * R = rhoR_D
    r_out = np.linspace(R_cm / n_shell, R_cm, n_shell)
    S = np.zeros(n_shell); S[0] = 1.0               # central emission
    snap = hto.ProfileSnapshot(
        r_cm=r_out, rho_gcc=np.full(n_shell, rho),
        w_D=np.full(n_shell, wD), w_T=np.full(n_shell, wT),
        w_C=np.full(n_shell, wC), S=S, s_per_shell=True,
        label=label or f"uniform_rhoRD{rhoR_D:.2f}")
    return snap, rho


def central_source(Ti_keV=20.0):
    import openmc
    return _mk_source(openmc, openmc.stats.Point((0.0, 0.0, 0.0)), hto._dt_muir())


def _mk_source(openmc, space, e_dist):
    try:
        return openmc.IndependentSource(space=space, angle=openmc.stats.Isotropic(),
                                        energy=e_dist, particle="neutron")
    except AttributeError:
        return openmc.Source(space=space, angle=openmc.stats.Isotropic(),
                             energy=e_dist, particle="neutron")


def build_sweep_tallies(surfs, ncoll_max=20):
    """D-elastic collision count + collision-resolved leakage current.
    Returns (tallies, have_collision_filter)."""
    import openmc
    outer = surfs[-1]
    tallies = openmc.Tallies()
    t_del = openmc.Tally(name="D_elastic")
    t_del.nuclides = ["H2"]; t_del.scores = ["elastic"]
    tallies.append(t_del)

    filters = [openmc.SurfaceFilter(outer.id),
               openmc.ParticleFilter(["neutron"]),
               openmc.EnergyFilter(EDGES * 1e6)]
    have_coll = False
    try:
        filters.append(openmc.CollisionFilter(list(range(0, ncoll_max + 1))))
        have_coll = True
    except Exception as exc:                                  # pragma: no cover
        print(f"[sweep] WARNING: openmc.CollisionFilter unavailable ({exc}); "
              "single-scatter partition disabled -- full transport only.")
    t_n = openmc.Tally(name="neutron_leakage_coll")
    t_n.filters = filters; t_n.scores = ["current"]
    tallies.append(t_n)
    return tallies, have_coll


# --------------------------------------------------------------------------
# one sweep point
# --------------------------------------------------------------------------
def _partition_leakage(mean, have_coll, ncoll_max):
    """Reshape leakage tally mean -> (full_counts[E], single_counts[E], multi[E])."""
    if have_coll:
        ncb = ncoll_max + 1
        arr = np.asarray(mean).reshape(1, 1, NE, ncb)[0, 0]  # [energy, collision]
        full = arr.sum(axis=1)
        single = arr[:, :2].sum(axis=1)                      # 0 + 1 collisions
        multi = arr[:, 2:].sum(axis=1)
        # neutrons past the last collision bin are dropped; warn if non-negligible
        return full, single, multi
    full = np.asarray(mean).reshape(1, 1, NE)[0, 0]
    return full, None, None


def analyze_point(rhoR_D, snap, C_D, leak_mean, have_coll, ncoll_max, sigma_el_D):
    full, single, multi = _partition_leakage(leak_mean, have_coll, ncoll_max)
    dsr_full = hto.dsr_from_counts(EDGES, full)
    dsr_single = hto.dsr_from_counts(EDGES, single) if have_coll else float("nan")
    rhoR_openmc = C_D * hto.M_D / (hto.N_A * sigma_el_D * 1e-24)
    rhoR_analytic = hto.analytic_rhoR_species(snap, "D")
    tot = float(full.sum())
    multi_frac = float(multi.sum() / tot) if (have_coll and tot > 0) else float("nan")
    return {
        "rhoR_D_target": rhoR_D,
        "rhoR_D_analytic": rhoR_analytic,
        "rhoR_D_openmc": rhoR_openmc,
        "rhoR_enhancement": (rhoR_openmc / rhoR_analytic) if rhoR_analytic > 0 else float("nan"),
        "DSR_single": dsr_single,
        "DSR_full": dsr_full,
        "DSR_frac_diff": ((dsr_full - dsr_single) / dsr_single)
        if (have_coll and dsr_single > 0) else float("nan"),
        "multiscatter_frac": multi_frac,
        "_full": full, "_single": single, "_multi": multi,
    }


def run_point(rhoR_D, args):
    import openmc
    snap, rho = make_uniform_snap(rhoR_D, args.R_cm, args.carbon_frac, args.n_shell)
    mats, geom, cells, surfs = hto.build_target(snap)
    tallies, have_coll = build_sweep_tallies(surfs, args.ncoll_max)
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = int(args.particles)
    settings.batches = int(args.batches)
    settings.photon_transport = False
    settings.source = central_source(args.Ti)
    settings.output = {"tallies": False}
    model = openmc.Model(geom, mats, settings, tallies)
    sp_path = model.run()
    sigma = hto.sigma_el_D_barn_at()
    with openmc.StatePoint(sp_path) as sp:
        C_D = float(sp.get_tally(name="D_elastic").mean.ravel()[0])
        leak = sp.get_tally(name="neutron_leakage_coll").mean
    print(f"[sweep] rhoR_D={rhoR_D:.3f}  rho={rho:.1f} g/cc  done "
          f"(sigma_el_D={sigma:.3f} b)")
    return analyze_point(rhoR_D, snap, C_D, leak, have_coll, args.ncoll_max, sigma)


# --------------------------------------------------------------------------
# self-test: mock the OpenMC tallies so the analysis/plot pipeline is exercised
# without OpenMC.  Physics is schematic, only the trend + code path are tested.
# --------------------------------------------------------------------------
def mock_point(rhoR_D, args):
    snap, rho = make_uniform_snap(rhoR_D, args.R_cm, args.carbon_frac, args.n_shell)
    ctr = 0.5 * (EDGES[:-1] + EDGES[1:])
    tau = _TAU_PER_GCM2 * rhoR_D / max(0.5 - 0.5 * args.carbon_frac, 1e-3) * (0.5)  # ~ total tau
    p_single = 1.0 - np.exp(-tau)
    p_multi = max(tau * tau * 0.5 * np.exp(-tau), 0.0)
    ncb = args.ncoll_max + 1
    arr = np.zeros((NE, ncb))
    # 0 collisions: primary peak near 14 MeV
    prim = np.exp(-((ctr - 14.03) / 0.18) ** 2)
    arr[:, 0] = prim / prim.sum()
    # 1 collision: down-scatter shelf 2-13 MeV, scaled by p_single
    shelf = np.where((ctr > 2) & (ctr < 13), np.exp((ctr - 13) / 4.0), 0.0)
    arr[:, 1] = p_single * shelf / (shelf.sum() or 1.0)
    # >=2 collisions: softer, pushed below 10 MeV, scaled by p_multi
    soft = np.where(ctr < 11, np.exp((ctr - 11) / 3.5), 0.0)
    if ncb > 2:
        arr[:, 2] = p_multi * soft / (soft.sum() or 1.0)
    C_D = (hto.N_A / hto.M_D) * hto.SIGMA_EL_D_BARN * 1e-24 * hto.analytic_rhoR_species(snap, "D") \
        * (1.0 + 0.8 * p_multi)                              # readout inflated by re-scatter
    return analyze_point(rhoR_D, snap, C_D, arr.ravel(), True, args.ncoll_max,
                         hto.SIGMA_EL_D_BARN)


# --------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------
def write_csv(rows, path):
    cols = ["rhoR_D_target", "rhoR_D_analytic", "rhoR_D_openmc", "rhoR_enhancement",
            "DSR_single", "DSR_full", "DSR_frac_diff", "multiscatter_frac"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols)
        for r in rows:
            w.writerow([f"{r[c]:.5g}" for c in cols])
    print(f"[sweep] wrote {path}")


def make_plots(rows, prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = [r["rhoR_D_target"] for r in rows]
    dsr_s = [100 * r["DSR_single"] for r in rows]
    dsr_f = [100 * r["DSR_full"] for r in rows]
    enh = [r["rhoR_enhancement"] for r in rows]
    mfrac = [100 * r["multiscatter_frac"] for r in rows]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    a1.plot(x, dsr_s, "o--", color="#1f4e79", lw=2, label="single-scatter (<=1 coll.)")
    a1.plot(x, dsr_f, "s-", color="#111111", lw=2, label="full transport (all coll.)")
    a1.set_xlabel("fuel  rhoR_D  (g/cm$^2$)"); a1.set_ylabel("DSR  (%)")
    a1.set_title("DSR: single vs multiple scattering"); a1.grid(alpha=0.3)
    a1.legend(fontsize=9)
    a2b = a2.twinx()
    a2.plot(x, mfrac, "^-", color="#990000", lw=2, label="multiply-scattered escaping (%)")
    a2b.plot(x, enh, "d--", color="#C08A1E", lw=2, label="rhoR_D readout enhancement")
    a2.set_xlabel("fuel  rhoR_D  (g/cm$^2$)")
    a2.set_ylabel("multiply-scattered fraction (%)", color="#990000")
    a2b.set_ylabel("OpenMC / analytic rhoR_D", color="#8a6410")
    a2.set_title("multiple-scattering signatures vs rhoR"); a2.grid(alpha=0.3)
    l1, la1 = a2.get_legend_handles_labels(); l2, la2 = a2b.get_legend_handles_labels()
    a2.legend(l1 + l2, la1 + la2, fontsize=8, loc="upper left")
    fig.tight_layout(); fig.savefig(f"{prefix}_summary.png", dpi=160)
    print(f"[sweep] wrote {prefix}_summary.png")

    # spectra single vs full at low/mid/high rhoR
    idx = sorted(set([0, len(rows) // 2, len(rows) - 1]))
    ctr = 0.5 * (EDGES[:-1] + EDGES[1:])
    fig, axes = plt.subplots(1, len(idx), figsize=(4.6 * len(idx), 4.2), squeeze=False)
    for ax, i in zip(axes[0], idx):
        r = rows[i]
        def norm(y):
            m = (ctr >= 13) & (ctr <= 15); a = y[m].sum()
            return y / (a if a > 0 else 1.0)
        ax.axvspan(10, 12, color="#cccccc", alpha=0.4)
        if r["_single"] is not None:
            ax.semilogy(ctr, norm(r["_single"]), "--", color="#1f4e79", lw=1.8,
                        label="single (<=1 coll.)")
        ax.semilogy(ctr, norm(r["_full"]), "k", lw=2.0, label="full")
        ax.set_xlim(1, 16); ax.set_ylim(1e-4, 5)
        ax.set_title(f"rhoR_D = {r['rhoR_D_target']:.2f}")
        ax.set_xlabel("neutron energy (MeV)"); ax.grid(alpha=0.25, which="both")
    axes[0][0].set_ylabel("dN/dE (norm. 13-15 MeV)"); axes[0][0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{prefix}_spectra.png", dpi=160)
    print(f"[sweep] wrote {prefix}_spectra.png")


def print_table(rows):
    print("\n" + "=" * 78)
    print("  MULTIPLE-SCATTERING SIGNATURE vs AREAL DENSITY  (single = <=1 collision)")
    print("=" * 78)
    print(f"  {'rhoR_D':>7} {'DSR_s%':>8} {'DSR_f%':>8} {'dDSR%':>8} "
          f"{'multi%':>8} {'rhoR_enh':>9}")
    for r in rows:
        print(f"  {r['rhoR_D_target']:7.2f} {100*r['DSR_single']:8.3f} "
              f"{100*r['DSR_full']:8.3f} {100*r['DSR_frac_diff']:8.1f} "
              f"{100*r['multiscatter_frac']:8.1f} {r['rhoR_enhancement']:9.3f}")
    print("=" * 78)
    print("  dDSR% = (DSR_full - DSR_single)/DSR_single;  rhoR_enh = OpenMC/analytic.")
    print("  Single scatter (NeSST-equivalent) is exact at low rhoR and breaks down")
    print("  as rhoR rises -- that divergence IS the multiple-scattering signature.\n")


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rhoR", default="0.1,0.3,0.6,1.0,1.5,2.0",
                    help="comma-separated fuel rhoR_D targets (g/cm^2)")
    ap.add_argument("--R-um", type=float, default=50.0, help="sphere radius (microns)")
    ap.add_argument("--carbon-frac", type=float, default=0.0, help="C mass fraction")
    ap.add_argument("--n-shell", type=int, default=40)
    ap.add_argument("--particles", type=int, default=200_000)
    ap.add_argument("--batches", type=int, default=20)
    ap.add_argument("--ncoll-max", type=int, default=20)
    ap.add_argument("--Ti", type=float, default=20.0, help="ion temperature (keV)")
    ap.add_argument("--prefix", default="scatter_sweep")
    ap.add_argument("--npz", default=None, help="also run this real snapshot")
    ap.add_argument("--self-test", action="store_true",
                    help="mock OpenMC; verify analysis/plot pipeline only")
    args = ap.parse_args(argv)
    args.R_cm = args.R_um * 1e-4

    rhoRs = [float(x) for x in args.rhoR.split(",") if x.strip()]
    fn = mock_point if args.self_test else run_point
    rows = [fn(v, args) for v in rhoRs]

    if args.npz and not args.self_test:
        rows.append(_run_real_npz(args))

    print_table(rows)
    write_csv(rows, f"{args.prefix}.csv")
    make_plots(rows, args.prefix)


def _run_real_npz(args):
    import openmc
    raw = hto.ProfileSnapshot.from_npz(args.npz)
    snap = hto.resample_min_thickness(raw)
    mats, geom, cells, surfs = hto.build_target(snap)
    tallies, have_coll = build_sweep_tallies(surfs, args.ncoll_max)
    settings = openmc.Settings()
    settings.run_mode = "fixed source"; settings.particles = int(args.particles)
    settings.batches = int(args.batches); settings.photon_transport = False
    settings.source = hto.build_source(snap, energy="birth")   # real S(r) + birth spectrum
    settings.output = {"tallies": False}
    model = openmc.Model(geom, mats, settings, tallies)
    sp_path = model.run()
    sigma = hto.sigma_el_D_barn_at()
    with openmc.StatePoint(sp_path) as sp:
        C_D = float(sp.get_tally(name="D_elastic").mean.ravel()[0])
        leak = sp.get_tally(name="neutron_leakage_coll").mean
    r = analyze_point(float("nan"), snap, C_D, leak, have_coll, args.ncoll_max, sigma)
    r["rhoR_D_target"] = hto.analytic_rhoR_species(snap, "D")   # label real point by its rhoR
    print(f"[sweep] real snapshot {snap.label}: DSR_single={100*r['DSR_single']:.2f}% "
          f"DSR_full={100*r['DSR_full']:.2f}%")
    return r


if __name__ == "__main__":
    main()
