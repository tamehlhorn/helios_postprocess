"""
scan_summary.py -- Aggregate key calibration metrics across a directory of
Helios runs into one CSV + ranked stdout table.

Built for the foot-pulse / HS rhoR closure sweep: walks every <base>_summary.txt
under a given root, pulls the shock train (foot/ramp/peak), peak rhoR and HS
rhoR, peak ablation pressure, peak velocity, adiabat, and HS pressure at
ignition. Emits scan_summary.csv next to the script's CWD (override with
--out) and prints a ranked table sorted by |total rhoR - LILAC|.

Usage:
    python3 ~/helios_postprocess/examples/scan_summary.py \\
        ~/Sims/Xcimer/Olson_PDD                    # walks recursively
    python3 ~/helios_postprocess/examples/scan_summary.py \\
        ~/Sims/Xcimer/Olson_PDD --out /tmp/scan.csv

LILAC reference values (Olson 2021 PDD anchor) used for distance ranking:
    total rhoR = 1.05 g/cm^2
    HS rhoR    = 0.85 g/cm^2
    foot shock = 7.50 ns
    ramp shock = 10.0 ns
    peak shock = 13.0 ns
Override via --lilac-total-rhor / --lilac-hs-rhor.
"""

import argparse
import csv
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional


# Default LILAC PDD reference values (Olson 2021)
LILAC = {
    'total_rhoR_gcm2':   1.05,
    'hs_rhoR_T_mask':    0.85,
    't_foot_ns':         7.50,
    't_ramp_ns':         10.0,
    't_peak_ns':         13.0,
}


# ── Summary-file parsers ─────────────────────────────────────────────────────

_RE_SCALAR = re.compile(
    r'^\s+(?P<label>[A-Za-z][^\n=]*?)\s{2,}'
    r'(?P<value>[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)'
)

# Scalar key -> list of regex patterns matched against summary lines.
# First match wins. The summary text is the public API here; labels are
# from icf_output.py (_metric helper writes "  <label:36s> <value> <unit>").
# When the summary format changes, update these patterns.
SCALAR_PATTERNS: Dict[str, List[re.Pattern]] = {
    'peak_velocity_kms':    [re.compile(r'Peak implosion velocity\s+(\d+\.\d+)\s*km/s')],
    'adiabat':              [re.compile(r'Mass-avg adiabat \(peak-vel\)\s+(\d+\.\d+)')],
    'base_adiabat':         [re.compile(r'Base adiabat \(at breakout\)\s+(\d+\.\d+)')],
    'hs_pressure_ignition': [re.compile(r'HS pressure at ignition\s+(\d+\.\d+)\s*Gbar')],
    'hs_radius_ignition_um':[re.compile(r'HS radius at ignition\s+(\d+\.\d+)')],
    'peak_total_rhoR':      [re.compile(r'Peak total ρR\s+(\d+\.\d+)')],
    'peak_hs_rhoR_T_mask':  [re.compile(r'Peak HS ρR \(T>4\.5 keV\)\s+(\d+\.\d+)')],
    'hot_spot_pressure':    [re.compile(r'Hot-spot pressure\s+(\d+\.\d+)\s*Gbar')],
    'yield_MJ':             [re.compile(r'Yield\s+(\d+\.\d+)\s*MJ'),
                             re.compile(r'Fusion energy output.*?(\d+\.\d+)\s*MJ')],
    'gain':                 [re.compile(r'Target gain\s+(\d+\.\d+)'),
                             re.compile(r'Fusion Gain\s+(\d+\.\d+)')],
    'stag_time_ns':         [re.compile(r'Stagnation time\s+(\d+\.\d+)\s*ns')],
    'bang_time_ns':         [re.compile(r'Bang time\s+(\d+\.\d+)\s*ns')],
    't_foot_shock_ns':      [],  # filled from SHOCK TRAIN block below
    't_ramp_shock_ns':      [],
    't_peak_shock_ns':      [],
}


# Shock-train table row: leading whitespace, class word, then five whitespace-
# separated columns. Anchored on the class name so we don't accidentally
# match other tables.
_RE_SHOCK_ROW = re.compile(
    r'^\s+(?P<cls>foot|ramp|peak)\s+'
    r'(?P<t>[-+]?\d+\.\d+)\s+'
    r'(?P<r>[-+]?\d+\.\d+)\s+'
    r'(?P<P>[-+]?\d+\.\d+)\s+'
    r'(?P<pr>[-+]?\d+\.\d+)\s+'
    r'(?P<n>\d+)\s*$'
)


def parse_summary(summary_path: Path) -> Dict[str, Optional[float]]:
    """Parse one _summary.txt into a flat metrics dict (None when missing)."""
    metrics: Dict[str, Optional[float]] = {k: None for k in SCALAR_PATTERNS}
    in_shock_block = False

    with open(summary_path) as f:
        for line in f:
            # Track entry/exit of the SHOCK TRAIN block to avoid grabbing
            # similarly-formatted rows from other tables.
            if line.startswith('SHOCK TRAIN'):
                in_shock_block = True
                continue
            if in_shock_block and line.strip() == '':
                in_shock_block = False
                continue

            if in_shock_block:
                m = _RE_SHOCK_ROW.match(line)
                if m:
                    key = f"t_{m.group('cls')}_shock_ns"
                    metrics[key] = float(m.group('t'))

            # Free-text scalar patterns
            for key, patterns in SCALAR_PATTERNS.items():
                if metrics[key] is not None or not patterns:
                    continue
                for pat in patterns:
                    m = pat.search(line)
                    if m:
                        try:
                            metrics[key] = float(m.group(1))
                            break
                        except (IndexError, ValueError):
                            pass

    return metrics


def find_summary_files(root: Path) -> List[Path]:
    """Recursive walk for *_summary.txt files under root."""
    return sorted(root.rglob('*_summary.txt'))


def run_name_from_summary(path: Path) -> str:
    """Strip the _summary.txt suffix; the parent dir matches by convention."""
    return path.stem[:-len('_summary')] if path.stem.endswith('_summary') else path.stem


# ── Ranking ────────────────────────────────────────────────────────────────


def distance_to_lilac(m: Dict[str, Optional[float]]) -> float:
    """Composite distance from LILAC anchor.

    Squared relative errors on the metrics we have a published value for.
    Missing metrics contribute a fixed penalty so runs with no shock train
    don't outrank converged runs.
    """
    score = 0.0
    pairs = [
        ('peak_total_rhoR',     LILAC['total_rhoR_gcm2']),
        ('peak_hs_rhoR_T_mask', LILAC['hs_rhoR_T_mask']),
        ('t_foot_shock_ns',     LILAC['t_foot_ns']),
        ('t_ramp_shock_ns',     LILAC['t_ramp_ns']),
    ]
    for key, ref in pairs:
        v = m.get(key)
        if v is None or not math.isfinite(v) or ref == 0:
            score += 1.0   # missing-data penalty (RMS error of 100%)
        else:
            score += ((v - ref) / ref) ** 2
    return math.sqrt(score)


# ── Main ───────────────────────────────────────────────────────────────────


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('root', type=Path,
                    help='Root directory to scan for *_summary.txt files.')
    ap.add_argument('--out', type=Path, default=Path('scan_summary.csv'),
                    help='Output CSV path (default: ./scan_summary.csv).')
    ap.add_argument('--lilac-total-rhor', type=float, default=None,
                    help='Override LILAC total rhoR target (g/cm^2).')
    ap.add_argument('--lilac-hs-rhor',    type=float, default=None,
                    help='Override LILAC HS rhoR target  (g/cm^2).')
    args = ap.parse_args(argv)

    if args.lilac_total_rhor is not None:
        LILAC['total_rhoR_gcm2'] = args.lilac_total_rhor
    if args.lilac_hs_rhor is not None:
        LILAC['hs_rhoR_T_mask'] = args.lilac_hs_rhor

    summaries = find_summary_files(args.root)
    if not summaries:
        print(f"No *_summary.txt under {args.root}", file=sys.stderr)
        return 1

    rows = []
    for sp in summaries:
        name = run_name_from_summary(sp)
        m = parse_summary(sp)
        rows.append({'run': name, **m,
                     'distance_to_lilac': distance_to_lilac(m)})

    rows.sort(key=lambda r: r['distance_to_lilac'])

    columns = [
        'run',
        't_foot_shock_ns', 't_ramp_shock_ns', 't_peak_shock_ns',
        'peak_total_rhoR', 'peak_hs_rhoR_T_mask',
        'hs_pressure_ignition', 'hs_radius_ignition_um',
        'hot_spot_pressure',
        'peak_velocity_kms', 'adiabat', 'base_adiabat',
        'yield_MJ', 'gain',
        'stag_time_ns', 'bang_time_ns',
        'distance_to_lilac',
    ]

    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, restval='')
        w.writeheader()
        for r in rows:
            w.writerow({k: ('' if r.get(k) is None else r[k]) for k in columns})
    print(f"Wrote: {args.out}  ({len(rows)} runs)")

    # ── Ranked stdout table ────────────────────────────────────────────────
    print()
    print(f"Ranked by |distance to LILAC|  "
          f"(total ρR={LILAC['total_rhoR_gcm2']:.2f}, "
          f"HS ρR={LILAC['hs_rhoR_T_mask']:.2f}, "
          f"foot={LILAC['t_foot_ns']:.1f} ns, ramp={LILAC['t_ramp_ns']:.1f} ns)")
    print()
    hdr = (f"  {'run':<46s}  {'foot':>5s}  {'ramp':>5s}  {'peak':>5s}  "
           f"{'totρR':>6s}  {'hsρR':>6s}  {'V':>5s}  {'α':>5s}  {'dist':>5s}")
    print(hdr)
    print('  ' + '-' * (len(hdr) - 2))
    def _f(v, w=5, p=2):
        if v is None or not isinstance(v, (int, float)) or not math.isfinite(v):
            return f"{'-':>{w}}"
        return f"{v:>{w}.{p}f}"
    for r in rows:
        print(f"  {r['run']:<46s}  "
              f"{_f(r.get('t_foot_shock_ns'))}  "
              f"{_f(r.get('t_ramp_shock_ns'))}  "
              f"{_f(r.get('t_peak_shock_ns'))}  "
              f"{_f(r.get('peak_total_rhoR'))}  "
              f"{_f(r.get('peak_hs_rhoR_T_mask'))}  "
              f"{_f(r.get('peak_velocity_kms'), 5, 0)}  "
              f"{_f(r.get('adiabat'))}  "
              f"{_f(r.get('distance_to_lilac'))}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
