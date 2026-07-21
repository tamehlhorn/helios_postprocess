#!/usr/bin/env python
"""Cross-check helios_postprocess extraction against Kyle's neutronics_data.npz.

RHINO-free: runs our own ``extract_neutronics`` on a Helios run and diffs its
neutron-weighted profiles, per-channel rates, yields, and bang time against the
``neutronics_data.npz`` produced by K. Keipper's ``neutronics_output.py``.

This is the "our independent implementation == his" validation gate that has to
hold before Will imports our routines into RHINO.

Usage
-----
    python examples/compare_to_kyle_npz.py <run_dir_or_exo> <kyle_npz> [--plot out.png]

``<run_dir_or_exo>`` is the Helios run directory (containing the .exo) or the
.exo path itself. ``<kyle_npz>`` is his neutronics_data.npz for the same run.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------
# Diff helpers (pure — unit-tested without a real .exo)
# --------------------------------------------------------------------------

def reldiff_scalar(a: float, b: float) -> float:
    """Signed relative difference (a-b)/|b|; 0 if both ~0, inf if only b~0."""
    a = float(a); b = float(b)
    if abs(b) < 1e-30:
        return 0.0 if abs(a) < 1e-30 else float("inf")
    return (a - b) / abs(b)


def reldiff_profile(r_ours, y_ours, r_theirs, y_theirs, floor_frac=1e-3):
    """Compare two radial profiles on possibly-different grids.

    Interpolates ours onto their grid and reports max/mean |rel diff| over the
    region where their profile exceeds ``floor_frac`` of its peak (ignores the
    noisy near-zero tails). Returns dict(max, mean, n).
    """
    r_ours = np.asarray(r_ours, float); y_ours = np.asarray(y_ours, float)
    r_theirs = np.asarray(r_theirs, float); y_theirs = np.asarray(y_theirs, float)
    y_ours_i = np.interp(r_theirs, r_ours, y_ours, left=0.0, right=0.0)
    peak = np.nanmax(np.abs(y_theirs)) if y_theirs.size else 0.0
    if peak <= 0:
        return {"max": float("nan"), "mean": float("nan"), "n": 0}
    mask = np.abs(y_theirs) > floor_frac * peak
    if not mask.any():
        return {"max": float("nan"), "mean": float("nan"), "n": 0}
    denom = np.where(np.abs(y_theirs[mask]) > 0, np.abs(y_theirs[mask]), np.nan)
    rd = np.abs(y_ours_i[mask] - y_theirs[mask]) / denom
    return {"max": float(np.nanmax(rd)), "mean": float(np.nanmean(rd)), "n": int(mask.sum())}


def compare_npz(ours: dict, theirs) -> list:
    """Build a list of (label, kind, result) comparison rows from two npz-like
    mappings. ``ours`` is a plain dict; ``theirs`` is an npz mapping (or dict).

    Scalars report signed rel diff; profiles report max/mean |rel diff|.
    """
    def g(m, k, default=None):
        try:
            v = m[k]
            return v.item() if hasattr(v, "item") and np.ndim(v) == 0 else v
        except (KeyError, IndexError, TypeError):
            return default

    rows = []
    # Scalars
    for key, label in [("bang_time_ns", "bang time (ns)"),
                       ("total_dt_yield", "DT yield"),
                       ("total_dd_yield", "DD yield")]:
        a, b = g(ours, key), g(theirs, key)
        if a is not None and b is not None:
            rows.append((label, "scalar", reldiff_scalar(a, b)))

    r_ours, r_theirs = g(ours, "r_avg_cm"), g(theirs, "r_avg_cm")
    # DT-weighted profiles (backward-compat top-level keys in Kyle's schema)
    if r_ours is not None and r_theirs is not None:
        for key, label in [("rho_avg_gcc", "<rho(r)>_DT (g/cc)"),
                           ("T_ion_avg_eV", "<T_ion(r)>_DT (eV)"),
                           ("rate_DT_avg", "<DT rate(r)>")]:
            a, b = g(ours, key), g(theirs, key)
            if a is not None and b is not None:
                rows.append((label, "profile", reldiff_profile(r_ours, a, r_theirs, b)))

    # Per-channel avg_results (dict of dicts) if both present
    a_avg = g(ours, "avg_results"); b_avg = g(theirs, "avg_results")
    if hasattr(a_avg, "item"):
        a_avg = a_avg.item()
    if hasattr(b_avg, "item"):
        b_avg = b_avg.item()
    if isinstance(a_avg, dict) and isinstance(b_avg, dict) and r_ours is not None and r_theirs is not None:
        for ch in ("DT_nHe4", "DD_nHe3", "TT_nnHe4"):
            ca, cb = a_avg.get(ch), b_avg.get(ch)
            if isinstance(ca, dict) and isinstance(cb, dict):
                for fk, fl in [("rho_avg_gcc", f"{ch} <rho(r)>"),
                               ("T_ion_avg_eV", f"{ch} <T_ion(r)>")]:
                    if fk in ca and fk in cb:
                        rows.append((fl, "profile",
                                     reldiff_profile(r_ours, ca[fk], r_theirs, cb[fk])))
    return rows


def format_report(rows: list) -> str:
    out = ["", f"{'field':32s} {'kind':8s} result", "-" * 64]
    for label, kind, res in rows:
        if kind == "scalar":
            out.append(f"{label:32s} {kind:8s} {res:+.3%}")
        else:
            out.append(f"{label:32s} {kind:8s} "
                       f"max |Δ|={res['max']:.2%}  mean |Δ|={res['mean']:.2%}  (n={res['n']})")
    out.append("")
    return "\n".join(out)


# --------------------------------------------------------------------------
# Runner (needs the .exo + pipeline; not exercised by the unit tests)
# --------------------------------------------------------------------------

def _find_exo(path) -> str:
    p = Path(path)
    if p.suffix == ".exo":
        return str(p)
    exos = sorted(p.glob("*.exo"))
    if not exos:
        raise FileNotFoundError(f"No .exo found in {p}")
    return str(exos[0])


def our_neutronics(run_path) -> dict:
    """Run our RHINO-free extraction and return an npz-like dict."""
    from helios_postprocess.core import HeliosRun
    from helios_postprocess.data_builder import build_run_data
    from helios_postprocess import neutron_spectrum as ns

    exo = _find_exo(run_path)
    data = build_run_data(HeliosRun(exo))
    nd = ns.extract_neutronics(data=data, use_rhino=False)
    if nd is None:
        raise RuntimeError("extract_neutronics returned None (no-burn / missing fields).")
    dt = (nd.avg_results.get("DT_nHe4") or {})
    return {
        "bang_time_ns": nd.bang_time_ns,
        "total_dt_yield": nd.total_dt_yield,
        "total_dd_yield": nd.total_dd_yield,
        "r_avg_cm": nd.r_avg_cm,
        "rho_avg_gcc": dt.get("rho_avg_gcc"),
        "T_ion_avg_eV": dt.get("T_ion_avg_eV"),
        "rate_DT_avg": (dt.get("rates_avg") or {}).get("DT_nHe4"),
        "avg_results": nd.avg_results,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description="Diff our neutronics extraction vs Kyle's npz.")
    ap.add_argument("run", help="Helios run directory or .exo path")
    ap.add_argument("kyle_npz", help="Kyle's neutronics_data.npz for the same run")
    ap.add_argument("--plot", default=None, help="optional PNG of the <rho(r)>/<T_i(r)> overlays")
    args = ap.parse_args(argv)

    ours = our_neutronics(args.run)
    theirs = np.load(args.kyle_npz, allow_pickle=True)
    rows = compare_npz(ours, theirs)
    print(format_report(rows))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(theirs["r_avg_cm"] * 1e4, theirs["rho_avg_gcc"], label="Kyle")
        ax[0].plot(ours["r_avg_cm"] * 1e4, ours["rho_avg_gcc"], "--", label="ours")
        ax[0].set_xlabel("r (um)"); ax[0].set_ylabel("<rho(r)>_DT (g/cc)"); ax[0].legend()
        ax[1].plot(theirs["r_avg_cm"] * 1e4, theirs["T_ion_avg_eV"], label="Kyle")
        ax[1].plot(ours["r_avg_cm"] * 1e4, ours["T_ion_avg_eV"], "--", label="ours")
        ax[1].set_xlabel("r (um)"); ax[1].set_ylabel("<T_ion(r)>_DT (eV)"); ax[1].legend()
        fig.tight_layout(); fig.savefig(args.plot, dpi=110)
        print(f"[compare] wrote {args.plot}")


if __name__ == "__main__":
    main()
