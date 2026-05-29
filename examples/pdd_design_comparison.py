#!/usr/bin/env python3
"""
pdd_design_comparison.py
========================

Multi-run comparison for the FL=0.012 PDD design study. Reads the CSV
populated by `dump_burn_region_density.py` and produces:

  1. A multi-panel comparison figure with reference lines:
     - V_peak (km/s)               vs LILAC reference (410), fab007-prod (421)
     - Effective coupling (%)       vs fab007-prod (73.1)
     - Peak total ρR (g/cm²)        vs LILAC (1.05), fab007-prod (1.01)
     - Mean foam ρ at stag (g/cc)   vs fab007 foam-burn (12), fab007-ice (30)
     - Adiabat                       vs fab007-prod (1.95), LILAC (3.0)
     - Total yield (MJ, burn ON)    vs LILAC (87), fab007-prod (26), fab007-ice (81)
       — drawn only if any run has burn-ON data

  2. A printed decision-matrix table classifying each run by yield and
     foam-burn productivity.

Usage
-----
    # Default: read the scan CSV in the repo
    python3 pdd_design_comparison.py

    # Custom CSV
    python3 pdd_design_comparison.py --csv path/to/results.csv

    # Save figure to non-default location
    python3 pdd_design_comparison.py --out comparisons/my_design.png

    # Include reference rows for fab007 production foam, fab007-ice, LILAC
    python3 pdd_design_comparison.py --add-refs

    # Order rows explicitly (otherwise sorted by label)
    python3 pdd_design_comparison.py --order baseline c33 c28 c25 c20
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# --- Reference values (production runs + LILAC published) ---------------------
# Each tuple: (label, V_peak, coupling%, total_rhoR, foam_rho_mean,
#              adiabat, yield_MJ, foam_yield_pct, HS_rhoR)
REFERENCES = {
    'fab007-foam (prod)': dict(V_peak_kms=421, coupling_pct=73.1, peak_total_rhoR=1.01,
                                rho_mean_foam_gcc=12, adiabat=1.95, yield_MJ=26.0,
                                foam_yield_pct=10.1, HS_rhoR_max=0.35),
    'fab007-ice':         dict(V_peak_kms=427, coupling_pct=73.1, peak_total_rhoR=1.0,
                                rho_mean_foam_gcc=30, adiabat=float('nan'), yield_MJ=81.4,
                                foam_yield_pct=float('nan'), HS_rhoR_max=0.71),
    'fab02 (over-drive)': dict(V_peak_kms=463, coupling_pct=84.0, peak_total_rhoR=float('nan'),
                                rho_mean_foam_gcc=float('nan'), adiabat=1.05, yield_MJ=59.0,
                                foam_yield_pct=26.5, HS_rhoR_max=float('nan')),
    'LILAC (reference)':  dict(V_peak_kms=410, coupling_pct=68.0, peak_total_rhoR=1.05,
                                rho_mean_foam_gcc=float('nan'), adiabat=3.0, yield_MJ=87.4,
                                foam_yield_pct=float('nan'), HS_rhoR_max=0.85),
}

# Reference horizontal lines for each panel
REF_LINES = {
    'V_peak_kms':        [(421, 'fab007-prod', '#888888'), (410, 'LILAC', '#cc3333')],
    'coupling_pct':      [(73.1, 'fab007-prod', '#888888')],
    'peak_total_rhoR':   [(1.01, 'fab007-prod', '#888888'), (1.05, 'LILAC', '#cc3333')],
    'rho_mean_foam_gcc': [(12, 'fab007 foam-burn', '#888888'), (30, 'fab007-ice', '#cc3333')],
    'adiabat':           [(1.95, 'fab007-prod', '#888888'), (3.0, 'LILAC', '#cc3333')],
    'yield_MJ':          [(26.0, 'fab007-prod', '#888888'), (81.4, 'fab007-ice', '#1a9850'),
                          (87.4, 'LILAC', '#cc3333')],
    'HS_rhoR_max':       [(0.30, 'ignition threshold', '#888888'),
                          (0.85, 'LILAC', '#cc3333')],
    'foam_yield_pct':    [(10.1, 'fab007-prod', '#888888'),
                          (26.5, 'fab02 over-drive', '#1a9850')],
}


def parse_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return float('nan')


NUMERIC_COLS = ['V_peak_kms', 'coupling_pct', 'peak_total_rhoR',
                'peak_inflight_rhoR', 'rho_peak_all_gcc', 'rho_peak_foam_gcc',
                'rho_mean_foam_gcc', 'rhoR_foam', 'foam_mass_total_mg',
                'adiabat', 'HS_rhoR_max', 'yield_MJ', 'foam_yield_pct',
                't_stag_ns', 'bang_time_ns']


def load_csv(csv_path: Path) -> List[dict]:
    rows = []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            # Normalize: every numeric column becomes a float (nan for missing/empty)
            for col in NUMERIC_COLS:
                r[col] = parse_float(r.get(col, ''))
            rows.append(r)
    return rows


def _fmt(x, fmt='{:.1f}'):
    """Format a value cleanly: nan or None → em-dash; numeric → fmt."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return '—'
    try:
        return fmt.format(float(x))
    except (ValueError, TypeError):
        return str(x)


def reorder(rows: List[dict], order: Optional[List[str]]) -> List[dict]:
    if order is None:
        return sorted(rows, key=lambda r: r['label'])
    label_to_row = {r['label']: r for r in rows}
    out = []
    for key in order:
        # Allow partial matches (e.g. 'c20' matches 'wf_fl012_c20')
        match = [r for lbl, r in label_to_row.items() if key in lbl]
        if match:
            out.append(match[0])
    # Append any unmatched rows at the end
    placed = set(id(r) for r in out)
    for r in rows:
        if id(r) not in placed:
            out.append(r)
    return out


def add_reference_rows(rows: List[dict], which: List[str]) -> List[dict]:
    """Append reference rows (fab007-prod, fab007-ice, LILAC) to the list."""
    keymap = {'fab007-foam': 'fab007-foam (prod)',
              'fab007-ice':  'fab007-ice',
              'fab02':       'fab02 (over-drive)',
              'lilac':       'LILAC (reference)'}
    out = list(rows)
    for w in which:
        key = keymap.get(w.lower(), w)
        if key in REFERENCES:
            ref = dict(REFERENCES[key])
            ref['label'] = key
            ref['run_path'] = '(reference)'
            ref['is_reference'] = True
            out.append(ref)
    return out


def classify(row: dict) -> str:
    """Return a one-word decision verdict for a row."""
    def _nan(x):
        return x is None or (isinstance(x, float) and np.isnan(x))
    y = row.get('yield_MJ')
    foam_share = row.get('foam_yield_pct')
    if _nan(y):
        rho_foam = row.get('rho_mean_foam_gcc')
        if _nan(rho_foam):
            return ''
        return ('compression OK' if rho_foam >= 25 else
                'marginal compression' if rho_foam >= 18 else
                'compression deficit')
    # Burn-on results
    if y >= 70:
        verdict = 'EXCELLENT'
    elif y >= 50:
        verdict = 'good'
    elif y >= 30:
        verdict = 'marginal'
    else:
        verdict = 'POOR'
    if not _nan(foam_share) and foam_share >= 25:
        verdict += ' (strong foam burn)'
    return verdict


def make_figure(rows: List[dict], out_path: Path) -> None:
    labels = [r['label'].replace('wf_fl012_', '').replace(' (reference)', '') for r in rows]
    n = len(rows)
    x = np.arange(n)

    # Build panel definitions: (column_name, title, ylabel, ref_lines_key)
    base_panels = [
        ('V_peak_kms',        'Peak implosion velocity',         'V_peak (km/s)',         'V_peak_kms'),
        ('coupling_pct',      'Effective coupling',              'Coupling (%)',          'coupling_pct'),
        ('peak_total_rhoR',   'Peak total ρR',                   'ρR (g/cm²)',            'peak_total_rhoR'),
        ('rho_mean_foam_gcc', 'Mean foam ρ at stagnation',       'ρ_foam (g/cc)',         'rho_mean_foam_gcc'),
        ('adiabat',           'Adiabat (cold-fuel, mass-avg)',   'Adiabat',                'adiabat'),
    ]
    # If any row has yield, add burn-on panels
    has_yield     = any(not np.isnan(r.get('yield_MJ', float('nan'))) for r in rows)
    has_foamshare = any(not np.isnan(r.get('foam_yield_pct', float('nan'))) for r in rows)
    if has_yield:
        base_panels.append(('yield_MJ', 'Total yield (burn ON)', 'Yield (MJ)', 'yield_MJ'))
    if has_foamshare:
        base_panels.append(('foam_yield_pct', 'Foam yield share', 'Foam yield (%)', 'foam_yield_pct'))

    n_panels = len(base_panels)
    ncols = 3
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.0 * nrows), squeeze=False)

    bar_colors = []
    for r in rows:
        if r.get('is_reference'):
            bar_colors.append('#cccccc')
        else:
            bar_colors.append('#1f77b4')

    for i, (col, title, ylbl, ref_key) in enumerate(base_panels):
        ax = axes[i // ncols, i % ncols]
        vals = [parse_float(r.get(col)) for r in rows]
        # Bar chart
        bars = ax.bar(x, [v if not np.isnan(v) else 0 for v in vals],
                      width=0.7, color=bar_colors, edgecolor='#333',
                      linewidth=0.6)
        # Numeric labels
        for j, (b, v) in enumerate(zip(bars, vals)):
            if not np.isnan(v):
                ax.text(b.get_x() + b.get_width()/2, b.get_height(),
                        f'{v:.1f}' if abs(v) >= 10 else f'{v:.2f}',
                        ha='center', va='bottom', fontsize=8.5)
        # Reference lines
        for ref_val, ref_lbl, ref_color in REF_LINES.get(ref_key, []):
            ax.axhline(ref_val, color=ref_color, ls='--', lw=1.2, alpha=0.8,
                       label=f'{ref_lbl}: {ref_val:.2f}'.rstrip('0').rstrip('.'))
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
        ax.set_ylabel(ylbl)
        ax.set_title(title, fontweight='bold', fontsize=10.5)
        ax.grid(alpha=0.25, axis='y')
        ax.legend(fontsize=8, loc='best')

    # Hide unused subplots
    for k in range(n_panels, nrows * ncols):
        axes[k // ncols, k % ncols].set_visible(False)

    fig.suptitle('PDD design study: FL=0.012 scan vs production references',
                 fontweight='bold', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def print_table(rows: List[dict]) -> None:
    def _nan(x):
        return x is None or (isinstance(x, float) and np.isnan(x))
    has_yield = any(not _nan(r.get('yield_MJ')) for r in rows)
    print()
    print('=' * 110)
    print('  PDD design comparison')
    print('=' * 110)
    if has_yield:
        hdr = f"{'label':22s} {'V_peak':>7s} {'coup':>6s} {'totRR':>7s} {'foam_ρ':>7s} {'adia':>5s} {'yield':>6s} {'HS_ρR':>6s} {'foam%':>6s}  verdict"
    else:
        hdr = f"{'label':22s} {'V_peak':>7s} {'coup':>6s} {'totRR':>7s} {'foam_ρ':>7s} {'adia':>5s}  verdict"
    print(hdr)
    print('-' * 110)
    for r in rows:
        lbl = r['label'].replace('wf_fl012_','')
        if has_yield:
            print(f"{lbl:22s} {_fmt(r.get('V_peak_kms')):>7s} {_fmt(r.get('coupling_pct')):>6s} "
                  f"{_fmt(r.get('peak_total_rhoR'), '{:.2f}'):>7s} {_fmt(r.get('rho_mean_foam_gcc')):>7s} "
                  f"{_fmt(r.get('adiabat'), '{:.2f}'):>5s} {_fmt(r.get('yield_MJ'), '{:.1f}'):>6s} "
                  f"{_fmt(r.get('HS_rhoR_max'), '{:.2f}'):>6s} {_fmt(r.get('foam_yield_pct'), '{:.1f}'):>6s}  "
                  f"{classify(r)}")
        else:
            print(f"{lbl:22s} {_fmt(r.get('V_peak_kms')):>7s} {_fmt(r.get('coupling_pct')):>6s} "
                  f"{_fmt(r.get('peak_total_rhoR'), '{:.2f}'):>7s} {_fmt(r.get('rho_mean_foam_gcc')):>7s} "
                  f"{_fmt(r.get('adiabat'), '{:.2f}'):>5s}  {classify(r)}")
    print()


def main():
    ap = argparse.ArgumentParser(description='Plot + tabulate the PDD design study scan.')
    repo_root = Path(__file__).resolve().parent.parent
    ap.add_argument('--csv', type=Path,
                    default=repo_root / 'notebooks' / 'pdd_scan_results.csv',
                    help='Input CSV (default: notebooks/pdd_scan_results.csv)')
    ap.add_argument('--out', type=Path,
                    default=repo_root / 'comparisons' / 'pdd_design_comparison.png',
                    help='Output figure path')
    ap.add_argument('--add-refs', nargs='*', default=None,
                    metavar='REF',
                    help='Append reference rows for comparison. Valid: '
                         'fab007-foam fab007-ice fab02 lilac. With no value, adds all four.')
    ap.add_argument('--order', nargs='+', default=None,
                    help='Order rows by these label substrings (e.g. baseline c33 c28 c25 c20)')
    args = ap.parse_args()

    if not args.csv.exists():
        ap.error(f'CSV not found: {args.csv}')
    rows = load_csv(args.csv)
    if not rows:
        print('No rows in CSV.')
        return

    if args.add_refs is not None:
        which = args.add_refs if args.add_refs else ['fab007-foam', 'fab007-ice', 'fab02', 'lilac']
        rows = add_reference_rows(rows, which)

    rows = reorder(rows, args.order)

    print_table(rows)
    make_figure(rows, args.out)
    print(f'Wrote figure: {args.out}')


if __name__ == '__main__':
    sys.exit(main())
