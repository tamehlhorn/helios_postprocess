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
    'peak_velocity_kms': 410.0,
    'adiabat':           3.0,
}

# Composite distance is a weighted-RMS fractional error across metrics that
# have both a sim value and a LILAC reference. Defaults emphasize HS rhoR
# (the stated closure metric); foot/ramp shock timing are primary
# diagnostics; velocity and adiabat are secondary. Peak shock is optional
# -- it contributes to distance ONLY when sim detected it (most Helios
# alpha=1.05 runs don't produce a 3rd shock).
#
# A run with all four foot/ramp shock arrivals AND all four implosion
# scalars matching LILAC perfectly gets dist=0. Each metric off by 100%
# adds its weight to the squared sum. Missing data on a mandatory metric
# adds a full-weight penalty; missing peak shock adds nothing.
DEFAULT_WEIGHTS = {
    'hs_rhoR_T_mask':    2.0,   # primary closure metric
    'total_rhoR_gcm2':   1.0,
    't_foot_ns':         1.0,
    't_ramp_ns':         1.0,
    'peak_velocity_kms': 0.5,
    'adiabat':           0.5,
    't_peak_ns':         0.5,   # optional: contributes only if detected
}

# Adiabat outliers from broken cold-fuel zone normalization (alpha > 20 or
# < 0.5) are treated as missing data rather than fed into the distance.
ADIABAT_SANITY_RANGE = (0.5, 20.0)

# Multiplicative discount applied when ALL three shocks (foot, ramp, peak)
# were detected at the gas/ice interface -- reflects structural similarity
# to LILAC's 3-shock R-T pattern.
THREE_SHOCK_BONUS = 0.85


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


def distance_to_lilac(m: Dict[str, Optional[float]],
                      weights: Dict[str, float] = None) -> float:
    """Composite distance from LILAC anchor.

    Weighted-RMS fractional error across metrics with a LILAC reference.
    Missing metrics on weighted dimensions contribute a full-weight
    penalty (1.0 = 100% off). Peak shock is special: missing peak shock
    adds nothing (most alpha=1.05 runs don't produce a third shock, and
    we don't want to penalize them for it). A multiplicative discount is
    applied for runs that detected all three shocks at the gas/ice
    interface.
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS

    # (sim_key_in_parsed_metrics, LILAC_reference_key, is_optional)
    pairs = [
        ('peak_total_rhoR',      'total_rhoR_gcm2',   False),
        ('peak_hs_rhoR_T_mask',  'hs_rhoR_T_mask',    False),
        ('t_foot_shock_ns',      't_foot_ns',         False),
        ('t_ramp_shock_ns',      't_ramp_ns',         False),
        ('peak_velocity_kms',    'peak_velocity_kms', False),
        ('adiabat',              'adiabat',           False),
        ('t_peak_shock_ns',      't_peak_ns',         True),   # optional
    ]

    score = 0.0
    weight_sum = 0.0
    for sim_key, ref_key, optional in pairs:
        w = weights.get(ref_key, 0.0)
        if w <= 0:
            continue
        ref = LILAC.get(ref_key, 0.0)
        if ref == 0:
            continue
        v = m.get(sim_key)
        # Adiabat sanity gate: outliers from broken zone-normalisation
        # treated as missing rather than letting them dominate the score.
        if ref_key == 'adiabat' and v is not None:
            lo, hi = ADIABAT_SANITY_RANGE
            if not (lo <= v <= hi):
                v = None
        if v is None or not isinstance(v, (int, float)) or not math.isfinite(v):
            if optional:
                continue                # missing optional metric: skip entirely
            score      += w * 1.0       # missing mandatory metric: full penalty
            weight_sum += w
            continue
        score      += w * ((v - ref) / ref) ** 2
        weight_sum += w

    base = math.sqrt(score / weight_sum) if weight_sum > 0 else float('inf')

    # 3-shock structural bonus: runs that detected all three shocks at the
    # gas/ice interface get a small multiplicative discount, reflecting
    # similarity to LILAC's R-T pattern even if individual times are off.
    n_shocks_detected = sum(
        1 for k in ('t_foot_shock_ns', 't_ramp_shock_ns', 't_peak_shock_ns')
        if m.get(k) is not None and math.isfinite(m[k])
    )
    if n_shocks_detected == 3:
        base *= THREE_SHOCK_BONUS
    return base


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
    ap.add_argument('--lilac-velocity',   type=float, default=None,
                    help='Override LILAC peak velocity target (km/s).')
    ap.add_argument('--lilac-adiabat',    type=float, default=None,
                    help='Override LILAC adiabat target.')
    # Per-metric weights (0 disables a dimension entirely).
    ap.add_argument('--w-hs-rhor',  type=float, default=DEFAULT_WEIGHTS['hs_rhoR_T_mask'],
                    help='Weight on HS rhoR error (default 2.0 -- the closure metric).')
    ap.add_argument('--w-tot-rhor', type=float, default=DEFAULT_WEIGHTS['total_rhoR_gcm2'],
                    help='Weight on total rhoR error (default 1.0).')
    ap.add_argument('--w-foot',     type=float, default=DEFAULT_WEIGHTS['t_foot_ns'],
                    help='Weight on foot shock time error (default 1.0).')
    ap.add_argument('--w-ramp',     type=float, default=DEFAULT_WEIGHTS['t_ramp_ns'],
                    help='Weight on ramp shock time error (default 1.0).')
    ap.add_argument('--w-peak',     type=float, default=DEFAULT_WEIGHTS['t_peak_ns'],
                    help='Weight on peak shock time error (default 0.5, optional).')
    ap.add_argument('--w-velocity', type=float, default=DEFAULT_WEIGHTS['peak_velocity_kms'],
                    help='Weight on peak velocity error (default 0.5).')
    ap.add_argument('--w-adiabat',  type=float, default=DEFAULT_WEIGHTS['adiabat'],
                    help='Weight on adiabat error (default 0.5).')
    args = ap.parse_args(argv)

    if args.lilac_total_rhor is not None: LILAC['total_rhoR_gcm2']   = args.lilac_total_rhor
    if args.lilac_hs_rhor    is not None: LILAC['hs_rhoR_T_mask']    = args.lilac_hs_rhor
    if args.lilac_velocity   is not None: LILAC['peak_velocity_kms'] = args.lilac_velocity
    if args.lilac_adiabat    is not None: LILAC['adiabat']           = args.lilac_adiabat

    weights = {
        'hs_rhoR_T_mask':    args.w_hs_rhor,
        'total_rhoR_gcm2':   args.w_tot_rhor,
        't_foot_ns':         args.w_foot,
        't_ramp_ns':         args.w_ramp,
        't_peak_ns':         args.w_peak,
        'peak_velocity_kms': args.w_velocity,
        'adiabat':           args.w_adiabat,
    }

    summaries = find_summary_files(args.root)
    if not summaries:
        print(f"No *_summary.txt under {args.root}", file=sys.stderr)
        return 1

    rows = []
    for sp in summaries:
        name = run_name_from_summary(sp)
        m = parse_summary(sp)
        rows.append({'run': name, **m,
                     'distance_to_lilac': distance_to_lilac(m, weights)})

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
    print(f"LILAC targets: totρR={LILAC['total_rhoR_gcm2']:.2f}, "
          f"HSρR={LILAC['hs_rhoR_T_mask']:.2f}, "
          f"foot={LILAC['t_foot_ns']:.1f}, ramp={LILAC['t_ramp_ns']:.1f}, "
          f"peak={LILAC['t_peak_ns']:.1f} ns, "
          f"V={LILAC['peak_velocity_kms']:.0f} km/s, α={LILAC['adiabat']:.2f}")
    active_w = ", ".join(f"{k.split('_')[0]}={w}"
                         for k, w in weights.items() if w > 0)
    print(f"Weights:       {active_w}    "
          f"3-shock bonus ×{THREE_SHOCK_BONUS}")
    print(f"α sanity range: {ADIABAT_SANITY_RANGE} "
          "(outside → treated as missing)")
    print()
    hdr = (f"  {'run':<46s}  {'foot':>5s}  {'ramp':>5s}  {'peak':>5s}  "
           f"{'totρR':>6s}  {'hsρR':>6s}  {'V':>5s}  {'α':>5s}  "
           f"{'3sh':>3s}  {'dist':>5s}")
    print(hdr)
    print('  ' + '-' * (len(hdr) - 2))
    def _f(v, w=5, p=2):
        if v is None or not isinstance(v, (int, float)) or not math.isfinite(v):
            return f"{'-':>{w}}"
        return f"{v:>{w}.{p}f}"
    def _alpha_str(v):
        if v is None or not isinstance(v, (int, float)) or not math.isfinite(v):
            return f"{'-':>5}"
        lo, hi = ADIABAT_SANITY_RANGE
        if not (lo <= v <= hi):
            return f"{v:>4.1f}*"     # asterisk = flagged broken (excluded from dist)
        return f"{v:>5.2f}"
    for r in rows:
        n_sh = sum(1 for k in ('t_foot_shock_ns', 't_ramp_shock_ns', 't_peak_shock_ns')
                   if r.get(k) is not None
                   and isinstance(r[k], (int, float)) and math.isfinite(r[k]))
        marker = '✓✓✓' if n_sh == 3 else f"{n_sh}/3"
        print(f"  {r['run']:<46s}  "
              f"{_f(r.get('t_foot_shock_ns'))}  "
              f"{_f(r.get('t_ramp_shock_ns'))}  "
              f"{_f(r.get('t_peak_shock_ns'))}  "
              f"{_f(r.get('peak_total_rhoR'))}  "
              f"{_f(r.get('peak_hs_rhoR_T_mask'))}  "
              f"{_f(r.get('peak_velocity_kms'), 5, 0)}  "
              f"{_alpha_str(r.get('adiabat'))}  "
              f"{marker:>3s}  "
              f"{_f(r.get('distance_to_lilac'))}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
