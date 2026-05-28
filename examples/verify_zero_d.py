#!/usr/bin/env python3
"""
verify_zero_d.py
================

Compare Helios's per-zone fusion and bremsstrahlung rates against
Bosch-Hale + NRL analytic for a uniform static (hydro-off) test.

Reads:
  - FusionRate_DT_nHe4   : reactions/s per zone (Helios native)
  - EnExchEleToRadTimeIntg : cumulative electron-to-radiation energy in J
                             (differentiate around the sampled timestep to get W)
  - mass_density, ion_temperature, zone_mass at the sampled timestep

Writes a row (or rows) to notebooks/verification_results.csv that the
foam_vs_ice_investigation notebook reads to render its results dataframe.

Usage
-----
    python verify_zero_d.py <base_path> [options]

Common invocations:
    # Pure DT at the middle timestep of a hydro-off run, default CSV
    python verify_zero_d.py /Volumes/.../DT_static_3keV/DT_static_3keV

    # 1.7%C wetted foam composition (matches fab007 foam), explicit label
    python verify_zero_d.py /Volumes/.../foam_static_3keV/foam_static_3keV \\
        --comp 1.7C-CH --label "foam-CH @3keV"

    # Custom composition (e.g. 5% C atomic, no H)
    python verify_zero_d.py /path/to/run --comp custom --fC 0.05 --fH 0.0

    # Sample a specific time rather than midpoint
    python verify_zero_d.py /path/to/run --time-ns 0.5

CSV columns:
    timestamp, label, run_path, composition, T_keV, rho_gcc,
    metric, helios, analytic, ratio, verdict, source

Re-runs on the same (label, metric) pair are appended (history kept);
the notebook deduplicates on (label, metric, T_keV) keeping the most recent.

Environment: same as run_analysis.py — needs `helios_postprocess` on PYTHONPATH
(invoke from repo root, or set PYTHONPATH=.).
"""
from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

# Ensure helios_postprocess is importable when run from anywhere
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from helios_postprocess import HeliosRun

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AMU_G    = 1.66053906660e-24
E_FUS_J  = 17.6e6 * 1.602176634e-19
M_D, M_T = 2.014, 3.016    # amu
M_C, M_H = 12.000, 1.008   # amu

COMP_PRESETS = {
    # atomic fractions (rest = equimolar DT)
    'pure-DT'  : dict(f_C=0.0,    f_H=0.0   ),
    '1.7C-CH'  : dict(f_C=0.0163, f_H=0.0163),  # CH binder, fab007-foam-like
    '1.7C-CH2' : dict(f_C=0.0149, f_H=0.0299),  # CH2 binder variant
    '1.7C-only': dict(f_C=0.017,  f_H=0.0   ),  # carbon only, no H
}

DEFAULT_CSV = REPO_ROOT / 'notebooks' / 'verification_results.csv'

# ---------------------------------------------------------------------------
# Analytic kernels (Bosch-Hale 1992 + NRL brems)
# ---------------------------------------------------------------------------

def sigmav_DT(T_keV):
    """Bosch-Hale 1992 DT reactivity in cm^3/s. Accepts scalar or array."""
    B_G = 34.3827
    mrc2 = 1124656.0
    C1, C2, C3 = 1.17302e-9, 1.51361e-2, 7.51886e-2
    C4, C5     = 4.60643e-3, 1.35000e-2
    C6, C7     = -1.06750e-4, 1.36600e-5
    T = np.asarray(T_keV, dtype=float)
    theta = T / (1.0 - (T*(C2 + T*(C4 + T*C6))) / (1.0 + T*(C3 + T*(C5 + T*C7))))
    xi = (B_G**2 / (4.0*theta))**(1.0/3.0)
    return C1 * theta * np.sqrt(xi / (mrc2 * T**3)) * np.exp(-3.0*xi)


def analytic_rates(rho_gcc: float, T_keV: float,
                   f_C: float = 0.0, f_H: float = 0.0) -> dict:
    """
    Compute analytic Bosch-Hale fusion + NRL brems for given state and
    atomic composition.  Returns rates per gram (for direct comparison
    with Helios per-zone fields divided by total mass).

    f_C, f_H : atomic fractions of carbon, hydrogen.  Remainder is
               equimolar D+T (i.e. n_D = n_T = (1 - f_C - f_H)/2 * n_total).
    """
    f_DT = 1.0 - f_C - f_H
    if f_DT <= 0:
        raise ValueError(f"f_C + f_H = {f_C+f_H:.3f} >= 1, no DT remaining")

    m_avg_amu = f_DT * (M_D + M_T)/2.0 + f_C * M_C + f_H * M_H
    n_total   = rho_gcc / (m_avg_amu * AMU_G)
    n_DT      = f_DT * n_total
    n_D = n_T = n_DT / 2.0
    n_C       = f_C * n_total
    n_H       = f_H * n_total
    n_e       = n_DT + 6.0*n_C + n_H           # fully ionized
    sumZ2nZ   = n_DT + 36.0*n_C + n_H

    R_fus_cc  = n_D * n_T * sigmav_DT(T_keV)
    P_brem_cc = 1.69e-32 * np.sqrt(T_keV*1000.0) * sumZ2nZ * n_e

    return {
        'R_fus_per_g_s'  : float(R_fus_cc / rho_gcc),
        'P_brem_W_per_g' : float(P_brem_cc / rho_gcc),
        'P_brem_TW_per_g': float(P_brem_cc / rho_gcc * 1e-12),
        'n_total': float(n_total),
        'n_DT'   : float(n_DT),
        'n_e'    : float(n_e),
        'sum_Z2nZ': float(sumZ2nZ),
        'm_avg_amu': float(m_avg_amu),
    }


# ---------------------------------------------------------------------------
# Helios reader
# ---------------------------------------------------------------------------

def _pick_time_idx(time_ns, target_ns):
    if target_ns is None:
        return len(time_ns) // 2
    return int(np.argmin(np.abs(time_ns - target_ns)))


def helios_rates(base_path: Path, target_ns: Optional[float] = None,
                 verbose: bool = False) -> dict:
    """Open Helios .exo and return Helios's per-gram fusion + brems rates."""
    exo = base_path.with_suffix('.exo')
    if not exo.exists():
        raise FileNotFoundError(f"EXODUS file not found: {exo}")

    run = HeliosRun(str(exo), verbose=verbose)
    try:
        time_s   = run.get_variable('time_whole')
        time_ns  = np.asarray(time_s) * 1e9
        n_t      = len(time_ns)
        t_idx    = _pick_time_idx(time_ns, target_ns)

        rho_z    = run.get_variable('mass_density',    time_idx=t_idx)
        T_eV_z   = run.get_variable('ion_temperature', time_idx=t_idx)
        mass_z   = run.get_variable('zone_mass',       time_idx=t_idx)
        Rfus_z   = run.get_variable('FusionRate_DT_nHe4', time_idx=t_idx)

        M_total      = float(np.sum(mass_z))
        rho_avg      = float(np.average(rho_z,  weights=mass_z))
        T_eV_avg     = float(np.average(T_eV_z, weights=mass_z))
        T_eV_minmax  = (float(np.min(T_eV_z)), float(np.max(T_eV_z)))
        rho_minmax   = (float(np.min(rho_z)),  float(np.max(rho_z)))
        rho_spread   = (rho_minmax[1] - rho_minmax[0]) / rho_avg  if rho_avg > 0  else float('nan')
        T_spread     = (T_eV_minmax[1] - T_eV_minmax[0]) / T_eV_avg if T_eV_avg>0 else float('nan')

        # Fusion rate: Helios's FusionRate_DT_nHe4 is in reactions/s/g
        # (mass-specific rate per zone), NOT reactions/s/zone. Confirmed
        # by the diagnostic on a 3 keV / n_DT=1e22 static test where each
        # zone reports ~1.09e+27 = per-gram analytic, regardless of n_zones
        # or zone volume. Take the mass-weighted average across zones
        # for the bulk rate.
        R_fus_pgs = float(np.average(Rfus_z, weights=mass_z))

        # Brems via cumulative-rad differentiation. Use central difference
        # when possible, fall back to forward/backward at endpoints (n_t=2
        # static tests have only the t=0 -> t=t_end pair).
        brems_pg_TW  = None
        brems_source = 'unavailable'
        for var, label in (('EnExchEleToRadTimeIntg',     'EnExchEleToRadTimeIntg/dt'),
                           ('TimeIntRadiationLossAtBds',  'TimeIntRadiationLossAtBds/dt')):
            try:
                E_cum = np.asarray(run.get_variable(var))
            except KeyError:
                continue
            if n_t < 2:
                continue
            if 0 < t_idx < n_t - 1:
                dE = E_cum[t_idx + 1] - E_cum[t_idx - 1]
                dt = (time_ns[t_idx + 1] - time_ns[t_idx - 1]) * 1e-9
                scheme = 'central'
            elif t_idx == 0:
                dE = E_cum[1] - E_cum[0]
                dt = (time_ns[1] - time_ns[0]) * 1e-9
                scheme = 'forward'
            else:                                # t_idx == n_t - 1
                dE = E_cum[-1] - E_cum[-2]
                dt = (time_ns[-1] - time_ns[-2]) * 1e-9
                scheme = 'backward'
            if dt <= 0:
                continue
            P_brem_W = dE / dt
            if M_total > 0:
                brems_pg_TW  = float((P_brem_W / M_total) * 1e-12)
                brems_source = f'{label} ({scheme})'
            break
    finally:
        run.close()

    return {
        'time_idx'        : t_idx,
        'time_ns'         : float(time_ns[t_idx]),
        'rho_gcc'         : rho_avg,
        'T_keV'           : T_eV_avg / 1000.0,
        'M_total_g'       : M_total,
        'rho_spread_pct'  : rho_spread * 100,
        'T_spread_pct'    : T_spread   * 100,
        'R_fus_per_g_s'   : R_fus_pgs,
        'P_brem_TW_per_g' : brems_pg_TW,
        'brems_source'    : brems_source,
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def verdict(ratio: float) -> str:
    r = abs(ratio - 1.0)
    if r < 0.02:  return 'EXACT (<2%)'
    if r < 0.10:  return 'within 10% (convention/physics)'
    if r < 0.30:  return 'within 30% (extra physics)'
    return 'INVESTIGATE'


def print_report(h: dict, a: dict, comp_label: str, base_name: str) -> dict:
    """Print a formatted comparison and return a dict of ratios."""
    print()
    print('=' * 78)
    print(f'  Zero-D verification: {base_name}')
    print('=' * 78)
    print(f'  Composition: {comp_label}')
    print(f'  State (mass-averaged):  T = {h["T_keV"]:.4f} keV   rho = {h["rho_gcc"]:.6f} g/cc')
    print(f'  Uniformity:  T spread = {h["T_spread_pct"]:.2f}%   rho spread = {h["rho_spread_pct"]:.2f}%')
    print(f'  Sampled at t = {h["time_ns"]:.4g} ns (index {h["time_idx"]})')
    print(f'  Total mass = {h["M_total_g"]:.4e} g')
    print()
    print(f'  {"Quantity":<22}  {"Helios":>15}  {"Analytic":>15}  {"H/A":>10}  Verdict')
    print(f'  {"-"*22}  {"-"*15}  {"-"*15}  {"-"*10}  {"-"*32}')

    r_fus = h['R_fus_per_g_s'] / a['R_fus_per_g_s'] if a['R_fus_per_g_s'] > 0 else float('nan')
    print(f'  {"Fusion rate (#/g/s)":<22}  {h["R_fus_per_g_s"]:15.4e}  {a["R_fus_per_g_s"]:15.4e}  {r_fus:10.4f}  {verdict(r_fus)}')

    r_br = float('nan')
    if h['P_brem_TW_per_g'] is not None:
        r_br = h['P_brem_TW_per_g'] / a['P_brem_TW_per_g'] if a['P_brem_TW_per_g'] > 0 else float('nan')
        print(f'  {"Brems (TW/g)":<22}  {h["P_brem_TW_per_g"]:15.4f}  {a["P_brem_TW_per_g"]:15.4f}  {r_br:10.4f}  {verdict(r_br)}')
        print(f'  {"":<22}  source: {h["brems_source"]}')
    else:
        print(f'  {"Brems (TW/g)":<22}  unavailable (no electron-to-rad accumulator in this run)')
    print()
    return {'r_fus': r_fus, 'r_brem': r_br}


CSV_COLUMNS = ['timestamp', 'label', 'run_path', 'composition',
               'T_keV', 'rho_gcc', 'metric', 'helios', 'analytic',
               'ratio', 'verdict', 'source']


def append_csv(csv_path: Path, h: dict, a: dict, ratios: dict,
               label: str, run_path: str, comp_label: str) -> int:
    """Append result rows to CSV.  Returns number of rows added."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    common = dict(
        timestamp   = datetime.now().isoformat(timespec='seconds'),
        label       = label,
        run_path    = run_path,
        composition = comp_label,
        T_keV       = round(h['T_keV'], 5),
        rho_gcc     = round(h['rho_gcc'], 8),
    )
    rows.append({
        **common,
        'metric'  : 'R_fus #/g/s',
        'helios'  : h['R_fus_per_g_s'],
        'analytic': a['R_fus_per_g_s'],
        'ratio'   : ratios['r_fus'],
        'verdict' : verdict(ratios['r_fus']),
        'source'  : 'FusionRate_DT_nHe4',
    })
    if h['P_brem_TW_per_g'] is not None:
        rows.append({
            **common,
            'metric'  : 'P_brem TW/g',
            'helios'  : h['P_brem_TW_per_g'],
            'analytic': a['P_brem_TW_per_g'],
            'ratio'   : ratios['r_brem'],
            'verdict' : verdict(ratios['r_brem']),
            'source'  : h['brems_source'],
        })

    file_exists = csv_path.exists()
    with csv_path.open('a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        if not file_exists:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return len(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description='Zero-D fusion + brems verification against analytic',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split('Usage')[1] if 'Usage' in __doc__ else '',
    )
    ap.add_argument('base_path',
                    help='Helios base path (no extension)')
    ap.add_argument('--time-ns', type=float, default=None,
                    help='Sample at this time in ns (default: midpoint)')
    ap.add_argument('--comp', default='pure-DT',
                    choices=list(COMP_PRESETS) + ['custom'],
                    help='Atomic composition preset (default: pure-DT)')
    ap.add_argument('--fC', type=float, default=None,
                    help='Atomic fraction C (with --comp custom)')
    ap.add_argument('--fH', type=float, default=None,
                    help='Atomic fraction H (with --comp custom)')
    ap.add_argument('--label', default=None,
                    help='Row label in CSV (default: base name)')
    ap.add_argument('--csv', default=None,
                    help=f'Output CSV (default: {DEFAULT_CSV})')
    ap.add_argument('--no-csv', action='store_true',
                    help='Print report only; do not append to CSV')
    ap.add_argument('-v', '--verbose', action='store_true')
    args = ap.parse_args()

    if args.comp == 'custom':
        if args.fC is None or args.fH is None:
            ap.error('--comp custom requires --fC and --fH')
        comp = dict(f_C=args.fC, f_H=args.fH)
        comp_label = f'custom(fC={args.fC},fH={args.fH})'
    else:
        comp = COMP_PRESETS[args.comp]
        comp_label = args.comp

    base = Path(args.base_path).expanduser().resolve()
    label = args.label or base.name
    csv_path = Path(args.csv).expanduser().resolve() if args.csv else DEFAULT_CSV

    h = helios_rates(base, target_ns=args.time_ns, verbose=args.verbose)
    a = analytic_rates(h['rho_gcc'], h['T_keV'], **comp)
    ratios = print_report(h, a, comp_label, base.name)

    if args.no_csv:
        print('  (CSV append skipped via --no-csv)')
    else:
        n = append_csv(csv_path, h, a, ratios, label, str(base), comp_label)
        print(f'  Appended {n} row(s) to {csv_path}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
