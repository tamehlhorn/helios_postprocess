#!/usr/bin/env python
"""Neutron post-processing for a Helios run: birth spectrum -> NeSST DSR -> nTOF.

Runs the full neutron chain on one Helios simulation and reports the
diagnostics-relevant metrics:

    Helios .exo
      -> extract_neutronics          (birth spectrum, nTOF Ti, DT-weighted rhoR)
      -> neutron_scatter (NeSST)      (singly-scattered spectrum, DSR, rhoR)
      -> scattered_tof (neutron_tof)  (synthetic nTOF trace w/ down-scatter tail)

Point it at the run directory (or the .exo). If a ``*_published.json`` sits
next to the run, the neutron-relevant published values are shown alongside ours.

Usage
-----
    python examples/analyze_neutron_scatter.py <run_dir_or_exo> \
        [--frac-D 0.5 --frac-T 0.5] [--distance 3.0] \
        [--published FILE.json] [--plot out.png] [--npz out.npz] \
        [--n-E 500] [--no-n2n]

Fuel fractions default to equimolar DT (0.5 / 0.5) -- correct for a standard NIF
layered capsule such as N210808. For a D-enhanced target pass e.g.
``--frac-D 0.7 --frac-T 0.3``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------
# Analysis (takes an already-extracted NeutronicsData; unit-testable)
# --------------------------------------------------------------------------

def analyze(nd, frac_D=0.5, frac_T=0.5, distance_m=3.0, n_E=500,
            include_n2n=True, data=None, rhw_path=None):
    """Run the NeSST scatter + nTOF chain on a NeutronicsData and return a
    flat dict of headline metrics plus the full scatter/tof result objects.

    If ``data`` and ``rhw_path`` are given, the multi-material path is used:
    per-layer composition from the RHW -> per-species areal densities ->
    D+T+carbon scatter (fuel fractions auto-detected)."""
    from helios_postprocess import neutron_scatter as nsc

    scat = None
    if data is not None and rhw_path:
        scat = nsc._composition_scatter(nd, data, str(rhw_path), max(n_E, 800))
    used_composition = scat is not None
    if scat is None:
        scat = nsc.scatter_from_extraction(
            nd, frac_D=frac_D, frac_T=frac_T, n_E=n_E, include_n2n=include_n2n)
    tof = nsc.scattered_tof(scat, distance_m=distance_m, irf_fwhm_ns=0.0)

    ql = nd.quicklook or {}
    dsr = scat["dsr"]["DSR"]
    rhoR_hydro = float(ql.get("rhoR_emission_gcm2", float("nan")))
    metrics = {
        "bang_time_ns": nd.bang_time_ns,
        "dt_yield": nd.total_dt_yield,
        "dd_yield": nd.total_dd_yield,
        "Ti_burn_avg_keV": ql.get("Ti_burn_avg_keV"),
        "Ti_ntof_keV": (tof or {}).get("primary", {}).get("Ti_ntof_keV"),
        "burn_fwhm_ns": ql.get("burn_fwhm_ns"),
        "rhoR_hydro_gcm2": rhoR_hydro,               # DT-weighted, from profile
        "DSR": dsr,                                   # fuel DSR (first-principles)
        "DSR_percent": 100.0 * dsr,
        "coeff_gcm2_per_DSR": (rhoR_hydro / dsr) if dsr else float("nan"),
        "frac_D": scat.get("frac_D", frac_D),
        "frac_T": scat.get("frac_T", frac_T),
        "distance_m": distance_m,
    }
    if used_composition:
        metrics["DSR_total_percent"] = 100.0 * scat["dsr_total"]["DSR"]
        metrics["rhoR_C_gcm2"] = scat["rhoR_C_gcm2"]
        metrics["composition_source"] = "RHW"
    else:
        metrics["rhoR_from_DSR_gcm2"] = scat["rhoR_from_A1s_gcm2"]
    return metrics, scat, tof


# --------------------------------------------------------------------------
# Published-data overlay (defensive: only shows keys that exist)
# --------------------------------------------------------------------------

# published-json key -> (our metric key, human label). Published values are
# [value, uncertainty]; several aliases are tried per row.
_PUB_ROWS = [
    (["yield", "yield_neutrons", "DT_yield"], "dt_yield", "DT yield"),
    (["Tion", "T_ion", "T_ion_keV", "T_hs"], "Ti_ntof_keV", "T_ion (keV)"),
    (["rhoR", "rhoR_cf", "rhoR_total", "rhoR_tot_gcm2"], "rhoR_hydro_gcm2",
     "rhoR (g/cm^2)"),
    (["DSR", "dsr", "DSR_percent"], "DSR_percent", "DSR (%)"),
    (["bang_time_ns", "bang_time"], "bang_time_ns", "bang time (ns)"),
]


def _pub_value(pub, keys):
    for k in keys:
        if k in pub:
            v = pub[k]
            return (float(v[0]) if isinstance(v, (list, tuple)) else float(v))
    return None


def format_report(metrics, pub=None):
    L = []
    L.append("")
    L.append("=" * 60)
    L.append("  NEUTRON POST-PROCESSING SUMMARY")
    L.append("=" * 60)
    L.append(f"  bang time            {metrics['bang_time_ns']:.3f} ns")
    L.append(f"  DT yield             {metrics['dt_yield']:.3e}")
    if metrics.get("dd_yield"):
        L.append(f"  DD yield             {metrics['dd_yield']:.3e}")
    L.append(f"  burn width (FWHM)    {metrics['burn_fwhm_ns']:.3f} ns")
    L.append(f"  T_ion (burn-avg)     {metrics['Ti_burn_avg_keV']:.2f} keV")
    L.append(f"  T_ion (nTOF)         {metrics['Ti_ntof_keV']:.2f} keV")
    L.append("-" * 60)
    L.append(f"  rhoR (hydro, DT-wt)  {metrics['rhoR_hydro_gcm2']:.3f} g/cm^2")
    if metrics.get("composition_source") == "RHW":
        L.append(f"  DSR fuel (D+T)       {metrics['DSR_percent']:.3f} %")
        L.append(f"  DSR total (+carbon)  {metrics['DSR_total_percent']:.3f} %"
                 f"   (rhoR_C = {metrics['rhoR_C_gcm2']:.3f} g/cm^2)")
    else:
        L.append(f"  DSR (NeSST 10-12/13-15 MeV)   {metrics['DSR_percent']:.3f} %")
        if metrics.get("rhoR_from_DSR_gcm2") is not None:
            L.append(f"  rhoR back from DSR   {metrics['rhoR_from_DSR_gcm2']:.3f}"
                     " g/cm^2   (round-trip check)")
    L.append(f"  first-principles coeff  {metrics['coeff_gcm2_per_DSR']:.2f}"
             " g/cm^2 per DSR   (empirical NIF 20.4)")
    L.append(f"  fuel loading D:T = {metrics['frac_D']:.2f}:{metrics['frac_T']:.2f}"
             f"   nTOF @ {metrics['distance_m']:.1f} m")

    if pub:
        L.append("-" * 60)
        L.append(f"  {'metric':22s}{'ours':>14s}{'published':>16s}")
        for keys, mkey, label in _PUB_ROWS:
            pv = _pub_value(pub, keys)
            ov = metrics.get(mkey)
            if pv is None or ov is None:
                continue
            L.append(f"  {label:22s}{ov:>14.3g}{pv:>16.3g}")
        L.append("  (published units vary by key; compare like-for-like -- "
                 "yield may be normalised, rhoR_cf is cold-fuel.)")
    L.append("=" * 60)
    L.append("")
    return "\n".join(L)


# --------------------------------------------------------------------------
# Run loading + CLI
# --------------------------------------------------------------------------

def _find_exo(path) -> str:
    """Resolve the .exo from any of the forms the user might pass:
    an explicit .exo, a base path (``<dir>/<stem>`` without extension, the
    run_analysis.py convention), a run directory, or a parent to search."""
    p = Path(path).expanduser()
    # 1) explicit .exo
    if p.suffix == ".exo" and p.is_file():
        return str(p)
    # 2) base path -> <path>.exo
    cand = p.with_suffix(".exo")
    if cand.is_file():
        return str(cand)
    # 3) directory: prefer <dir>/<dirname>.exo, else any .exo, else recurse
    if p.is_dir():
        named = p / f"{p.name}.exo"
        if named.is_file():
            return str(named)
        exos = sorted(p.glob("*.exo")) or sorted(p.rglob("*.exo"))
        if exos:
            return str(exos[0])
    # 4) base path whose parent holds the .exo (e.g. dir/stem -> dir/*.exo)
    if p.parent.is_dir():
        exos = sorted(p.parent.glob("*.exo"))
        if exos:
            return str(exos[0])
    raise FileNotFoundError(
        f"No .exo found for {path!r}. Pass the .exo file, the run directory, "
        f"or the base path (without extension). Looked at: {cand}, {p}/*.exo.")


def _find_published(run_path, exo) -> Path | None:
    base = Path(exo).with_suffix("")
    cand = Path(f"{base}_published.json")
    if cand.exists():
        return cand
    hits = sorted(Path(run_path).glob("*_published.json")) \
        if Path(run_path).is_dir() else []
    return hits[0] if hits else None


def _find_rhw(run_path, exo) -> Path | None:
    cand = Path(exo).with_suffix(".rhw")
    if cand.exists():
        return cand
    hits = sorted(Path(run_path).glob("*.rhw")) \
        if Path(run_path).is_dir() else []
    return hits[0] if hits else None


def make_figure(scat, tof, out_png, title=""):
    from helios_postprocess import neutron_scatter as nsc
    nsc.plot_scatter_tof(scat, tof, out_png, title=title)
    print(f"[analyze] wrote {out_png}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Neutron DSR / nTOF post-processing for a Helios run.")
    ap.add_argument("run", help="Helios run directory or .exo path")
    ap.add_argument("--frac-D", type=float, default=0.5, help="D atom fraction (default 0.5)")
    ap.add_argument("--frac-T", type=float, default=0.5, help="T atom fraction (default 0.5)")
    ap.add_argument("--distance", type=float, default=3.0, help="nTOF detector distance (m)")
    ap.add_argument("--published", default=None, help="published.json (else auto-detected)")
    ap.add_argument("--plot", default=None, help="output PNG (spectrum + nTOF)")
    ap.add_argument("--npz", default=None, help="also save the extraction npz here")
    ap.add_argument("--n-E", type=int, default=500, help="NeSST energy-grid points")
    ap.add_argument("--no-n2n", action="store_true", help="elastic only (faster, less memory)")
    ap.add_argument("--rhw", default=None,
                    help="RHW for per-layer composition (else auto-detected); "
                         "enables the multi-material D+T+carbon scatter")
    ap.add_argument("--no-composition", action="store_true",
                    help="ignore the RHW and use fuel-only D+T scatter")
    args = ap.parse_args(argv)

    from helios_postprocess.core import HeliosRun
    from helios_postprocess.data_builder import build_run_data
    from helios_postprocess import neutron_spectrum as ns

    exo = _find_exo(args.run)
    print(f"[analyze] loading {exo}")
    data = build_run_data(HeliosRun(exo))
    nd = ns.extract_neutronics(data=data, use_rhino=False,
                               save_npz=args.npz)
    if nd is None:
        raise SystemExit("extract_neutronics returned None (no-burn / missing fields).")

    rhw = None
    if not args.no_composition:
        rhw = Path(args.rhw) if args.rhw else _find_rhw(args.run, exo)
        if rhw and rhw.exists():
            print(f"[analyze] composition from {rhw}")
        else:
            rhw = None
    metrics, scat, tof = analyze(
        nd, frac_D=args.frac_D, frac_T=args.frac_T, distance_m=args.distance,
        n_E=args.n_E, include_n2n=not args.no_n2n, data=data, rhw_path=rhw)

    pub = None
    pub_path = Path(args.published) if args.published else _find_published(args.run, exo)
    if pub_path and Path(pub_path).exists():
        pub = json.loads(Path(pub_path).read_text())
        print(f"[analyze] published overlay: {pub_path}")

    print(format_report(metrics, pub))

    if args.plot:
        make_figure(scat, tof, args.plot, title=Path(exo).stem)


if __name__ == "__main__":
    main()
