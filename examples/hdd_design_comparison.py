#!/usr/bin/env python3
"""
hdd_design_comparison.py
========================

Multi-run HDD comparison tool — analog of pdd_design_comparison.py with
Thomas et al. (Vulcan HDD) reference values, RHINO-measured WT_cthomas
metrics, and HDD-relevant reference horizontal lines.

Reads `notebooks/hdd_scan_results.csv` (populated by
`examples/dump_burn_region_density.py`) and produces:

  1. Multi-panel comparison figure with reference horizontal lines for:
     - V_peak (km/s)            vs Thomas (410), RHINO-on-WT (404)
     - V_peak at CR=1.5         vs RHINO-on-WT (404)
     - V_impl RHINO (km/s)      vs Thomas (410), RHINO-on-WT (404)
       (W. Trickey convention: shell = rho>peak/e; v=sqrt(2KE/m_sh);
        turning point pre-stagnation)
     - Adiabat                  vs Thomas (6.0)
     - Adiabat at CR=1.5        vs RHINO-on-WT (4.13)
     - Min shell adiabat RHINO  vs Thomas (6.0), RHINO-on-WT (4.13)
       (W. Trickey convention: at timestep where inner shell surface
        reaches R_0/1.5, min adiabat over rho>peak/e shell zones)
     - Imploded DT mass (mg)    vs Thomas (3.0)
     - Total yield (MJ)         vs Thomas (256), WT_cthomas (196)
     - HS rhoR peak             vs Thomas implied (~3), RHINO-on-WT (6.4)
     - Gain                     vs Thomas (65)

  2. Printed decision-matrix table classifying each row.

Usage
-----
    # Default: read the HDD scan CSV in the repo
    python3 hdd_design_comparison.py

    # Custom CSV
    python3 hdd_design_comparison.py --csv path/to/results.csv

    # Save figure to non-default location
    python3 hdd_design_comparison.py --out comparisons/my_hdd.png

    # Append reference rows for direct visual comparison
    python3 hdd_design_comparison.py --add-refs

    # Order rows explicitly
    python3 hdd_design_comparison.py --order wt_cthomas vi6_baseline vi6_fl012 vi6_c25
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


# --- Reference values for HDD ----------------------------------------------
# Thomas et al. Vulcan HDD hybrid target with-alpha column from
# VI_6_published.json. RHINO numbers are Will Trickey's postprocess on
# the WT_cthomas_baseline .exo specifically -- not generic refs, but
# useful for cross-tool sanity checks.
REFERENCES = {
    'Thomas (Vulcan HDD)': dict(
        V_peak_kms=410, V_peak_kms_cr15=410, V_impl_rhino_kms=410,
        adiabat=6.0, adiabat_cr15=6.0, adiabat_min_rhino=6.0,
        peak_total_rhoR=1.60,    # rhoR_cf from JSON
        imploded_DT_mg=3.0,
        coupling_pct=97.0,
        yield_MJ=256.0,
        gain=65.0,
        HS_rhoR_max=float('nan'),  # not in Thomas table
    ),
    'RHINO on WT_cthomas': dict(
        V_peak_kms=404, V_peak_kms_cr15=404, V_impl_rhino_kms=404,
        adiabat=4.13, adiabat_cr15=4.13, adiabat_min_rhino=4.13,
        peak_total_rhoR=6.4,
        imploded_DT_mg=3.0,
        coupling_pct=85.5,
        yield_MJ=196.0,
        gain=63.2,
        HS_rhoR_max=6.4,
    ),
}

# Reference horizontal lines for each panel
REF_LINES = {
    'V_peak_kms':        [(410, 'Thomas', '#cc3333'),
                          (404, 'RHINO on WT', '#1a9850')],
    'V_peak_kms_cr15':   [(410, 'Thomas', '#cc3333'),
                          (404, 'RHINO on WT', '#1a9850')],
    'V_impl_rhino_kms':  [(410, 'Thomas', '#cc3333'),
                          (404, 'RHINO on WT', '#1a9850')],
    'adiabat':           [(6.0, 'Thomas', '#cc3333'),
                          (4.13, 'RHINO on WT', '#1a9850')],
    'adiabat_cr15':      [(6.0, 'Thomas', '#cc3333'),
                          (4.13, 'RHINO on WT', '#1a9850')],
    'adiabat_min_rhino': [(6.0, 'Thomas', '#cc3333'),
                          (4.13, 'RHINO on WT', '#1a9850')],
    'imploded_DT_mg':    [(3.0, 'Thomas', '#cc3333')],
    'yield_MJ':          [(256.0, 'Thomas', '#cc3333'),
                          (196.0, 'WT_cthomas (Helios)', '#888888')],
    'HS_rhoR_max':       [(6.4, 'RHINO on WT', '#1a9850'),
                          (0.30, 'ignition threshold', '#888888')],
    'coupling_pct':      [(97.0, 'Thomas', '#cc3333')],
}


def parse_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return float('nan')


NUMERIC_COLS = ['V_peak_kms', 'V_peak_kms_cr15', 'V_impl_rhino_kms',
                'coupling_pct', 'peak_total_rhoR', 'peak_inflight_rhoR',
                'rho_peak_all_gcc', 'rho_peak_foam_gcc', 'rho_mean_foam_gcc',
                'rhoR_foam', 'foam_mass_total_mg', 'imploded_DT_mg',
                'adiabat', 'adiabat_cr15', 'adiabat_min_rhino',
                'HS_rhoR_max', 'yield_MJ', 'foam_yield_pct',
                't_stag_ns', 'bang_time_ns']


def load_csv(csv_path: Path) -> List[dict]:
    rows = []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            for col in NUMERIC_COLS:
                r[col] = parse_float(r.get(col, ''))
            rows.append(r)
    return rows


def _fmt(x, fmt='{:.1f}'):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return '—'
    try:
        return fmt.format(float(x))
    except (ValueError, TypeError):
        return str(x)


def reorder(rows: List[dict], order: Optional[List[str]]) -> List[dict]:
    if order is None:
        return rows[:]
    out = []
    placed = set()
    for key in order:
        for r in rows:
            if key in r['label'] and id(r) not in placed:
                out.append(r)
                placed.add(id(r))
                break
    for r in rows:
        if id(r) not in placed:
            out.append(r)
    return out


def add_reference_rows(rows: List[dict], which: List[str]) -> List[dict]:
    keymap = {'thomas':  'Thomas (Vulcan HDD)',
              'rhino':   'RHINO on WT_cthomas',
              'rhino-wt':'RHINO on WT_cthomas'}
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
    """Return a one-word HDD-context verdict."""
    def _nan(x):
        return x is None or (isinstance(x, float) and np.isnan(x))
    y = row.get('yield_MJ')
    impl = row.get('imploded_DT_mg')
    HSrR = row.get('HS_rhoR_max')

    if _nan(y):
        return ''
    # Mass-deficit warning
    if not _nan(impl) and impl < 2.0:
        return 'mass-deficit (<2 mg imploded)'
    # No ignition
    if not _nan(HSrR) and HSrR < 0.3:
        return 'no ignition (HS ρR < 0.3)'
    # By yield
    if y >= 200:
        verdict = 'EXCELLENT (near Thomas)'
    elif y >= 100:
        verdict = 'good'
    elif y >= 50:
        verdict = 'marginal'
    elif y >= 20:
        verdict = 'low yield'
    else:
        verdict = 'POOR'
    return verdict


def make_figure(rows: List[dict], out_path: Path) -> None:
    labels = [r['label'] for r in rows]
    n = len(rows)
    x = np.arange(n)

    panels = [
        ('V_peak_kms',        'Peak velocity (legacy)',       'V_peak (km/s)',     'V_peak_kms'),
        ('V_peak_kms_cr15',   'Peak velocity at CR=1.5',     'V_peak (km/s)',     'V_peak_kms_cr15'),
        ('V_impl_rhino_kms',  'Implosion velocity (RHINO)',  'V (km/s)',          'V_impl_rhino_kms'),
        ('adiabat',           'Adiabat (legacy)',             'α',                 'adiabat'),
        ('adiabat_cr15',      'Adiabat at CR=1.5',           'α',                 'adiabat_cr15'),
        ('adiabat_min_rhino', 'Min shell adiabat (RHINO)',   'α',                 'adiabat_min_rhino'),
        ('coupling_pct',      'Effective coupling',           'Coupling (%)',      'coupling_pct'),
        ('imploded_DT_mg',    'Imploded DT mass',             'mass (mg)',         'imploded_DT_mg'),
        ('HS_rhoR_max',       'Peak HS ρR (T>4.5 keV)',      'ρR (g/cm²)',       'HS_rhoR_max'),
        ('yield_MJ',          'Total yield',                  'Yield (MJ)',        'yield_MJ'),
    ]

    n_panels = len(panels)
    ncols = 3
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(5.4 * ncols, 4.2 * nrows),
                              squeeze=False)

    bar_colors = ['#cccccc' if r.get('is_reference') else '#1f77b4'
                  for r in rows]

    for i, (col, title, ylbl, ref_key) in enumerate(panels):
        ax = axes[i // ncols, i % ncols]
        vals = [parse_float(r.get(col)) for r in rows]
        bars = ax.bar(x, [v if not np.isnan(v) else 0 for v in vals],
                      width=0.7, color=bar_colors, edgecolor='#333', linewidth=0.6)
        for j, (b, v) in enumerate(zip(bars, vals)):
            if not np.isnan(v):
                ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                        f'{v:.1f}' if abs(v) >= 10 else f'{v:.2f}',
                        ha='center', va='bottom', fontsize=8.5)
        for ref_val, ref_lbl, ref_color in REF_LINES.get(ref_key, []):
            ax.axhline(ref_val, color=ref_color, ls='--', lw=1.2, alpha=0.8,
                       label=f'{ref_lbl}: {ref_val:g}')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=8.5)
        ax.set_ylabel(ylbl)
        ax.set_title(title, fontweight='bold', fontsize=10.5)
        ax.grid(alpha=0.25, axis='y')
        ax.legend(fontsize=8, loc='best')

    for k in range(n_panels, nrows * ncols):
        axes[k // ncols, k % ncols].set_visible(False)

    fig.suptitle('HDD design comparison: scan rows vs Thomas / RHINO references',
                 fontweight='bold', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def print_table(rows: List[dict]) -> None:
    def _nan(x):
        return x is None or (isinstance(x, float) and np.isnan(x))
    print()
    print('=' * 130)
    print('  HDD design comparison')
    print('=' * 130)
    hdr = (f"{'label':24s} {'V_peak':>7s} {'V@cr1.5':>8s} {'V_RHINO':>8s} "
           f"{'adi':>5s} {'adi@cr1.5':>9s} {'a_RHINO':>8s} "
           f"{'coup%':>6s} {'impl_DT':>8s} {'HS_ρR':>6s} {'yield':>7s}  verdict")
    print(hdr)
    print('-' * 150)
    for r in rows:
        print(f"{r['label']:24s} "
              f"{_fmt(r.get('V_peak_kms')):>7s} "
              f"{_fmt(r.get('V_peak_kms_cr15')):>8s} "
              f"{_fmt(r.get('V_impl_rhino_kms')):>8s} "
              f"{_fmt(r.get('adiabat'), '{:.2f}'):>5s} "
              f"{_fmt(r.get('adiabat_cr15'), '{:.2f}'):>9s} "
              f"{_fmt(r.get('adiabat_min_rhino'), '{:.2f}'):>8s} "
              f"{_fmt(r.get('coupling_pct')):>6s} "
              f"{_fmt(r.get('imploded_DT_mg'), '{:.2f}'):>8s} "
              f"{_fmt(r.get('HS_rhoR_max'), '{:.2f}'):>6s} "
              f"{_fmt(r.get('yield_MJ')):>7s}  "
              f"{classify(r)}")
    print()


def main():
    ap = argparse.ArgumentParser(description='Plot + tabulate the HDD design study.')
    repo_root = Path(__file__).resolve().parent.parent
    ap.add_argument('--csv', type=Path,
                    default=repo_root / 'notebooks' / 'hdd_scan_results.csv',
                    help='Input CSV (default: notebooks/hdd_scan_results.csv)')
    ap.add_argument('--out', type=Path,
                    default=repo_root / 'comparisons' / 'hdd_design_comparison.png',
                    help='Output figure path')
    ap.add_argument('--add-refs', nargs='*', default=None,
                    metavar='REF',
                    help='Append reference rows. Valid: thomas, rhino. '
                         'With no value, adds both.')
    ap.add_argument('--order', nargs='+', default=None,
                    help='Order rows by these label substrings')
    args = ap.parse_args()

    if not args.csv.exists():
        ap.error(f'CSV not found: {args.csv}')
    rows = load_csv(args.csv)
    if not rows:
        print('No rows in CSV.')
        return

    if args.add_refs is not None:
        which = args.add_refs if args.add_refs else ['thomas', 'rhino']
        rows = add_reference_rows(rows, which)

    rows = reorder(rows, args.order)

    print_table(rows)
    make_figure(rows, args.out)
    print(f'Wrote figure: {args.out}')


if __name__ == '__main__':
    sys.exit(main())
