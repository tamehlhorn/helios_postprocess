#!/usr/bin/env python3
"""
make_pdd_scan_rhw.py
====================

Generate RHW files for a PDD geometric-parameter scan starting from a
production-style template (e.g. fab007). Designed for the foam-burn
recovery design study where we tune geometry to compensate for the
PROPACEOS foam EOS+opacity deficits relative to LILAC.

Per-run modifications:
  - Cone half-angle (beam 1 of the Laser Source Data block)
  - Spot size (beam 1)
  - Focus position (beam 1)
  - Flux limiter, all 4 regions (single value)
  - Burn on/off, all 4 regions (single value)
  - Max simulation time (truncate to save Studio time for no-burn runs)

The pulse table, target geometry, EOS/opacity tables, and beams 2-3
are NOT modified. This isolates geometric defocus / FL effects.

Usage examples
--------------
    # Single scan row
    python3 make_pdd_scan_rhw.py \\
        --template fab007.rhw \\
        --row "c37 cone=37 spot=0.18 focus=0.22 fl=0.012 burn=off max_t=1.5e-8" \\
        --outdir .

    # Full 5-row scan in one command (RECOMMENDED for the design study)
    python3 make_pdd_scan_rhw.py \\
        --template fab007.rhw \\
        --row "wf_fl012_baseline cone=37 spot=0.18 focus=0.22 fl=0.012 burn=off max_t=1.5e-8" \\
        --row "wf_fl012_c33      cone=33 spot=0.16 focus=0.20 fl=0.012 burn=off max_t=1.5e-8" \\
        --row "wf_fl012_c28      cone=28 spot=0.16 focus=0.18 fl=0.012 burn=off max_t=1.5e-8" \\
        --row "wf_fl012_c25      cone=25 spot=0.14 focus=0.16 fl=0.012 burn=off max_t=1.5e-8" \\
        --row "wf_fl012_c20      cone=20 spot=0.14 focus=0.15 fl=0.012 burn=off max_t=1.5e-8" \\
        --outdir .
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path
from typing import Optional


def _replace_first_field(text: str, field_prefix: str, new_value_str: str) -> str:
    """Replace the value on the FIRST matching line (preserves indent + alignment)."""
    lines = text.splitlines(keepends=True)
    pat = re.compile(r'^(\s*' + re.escape(field_prefix) + r'\s*=\s*)\S.*$')
    for i, line in enumerate(lines):
        m = pat.match(line.rstrip('\n').rstrip('\r'))
        if m:
            nl = '\r\n' if line.endswith('\r\n') else ('\n' if line.endswith('\n') else '')
            lines[i] = f'{m.group(1)}{new_value_str}{nl}'
            return ''.join(lines)
    raise ValueError(f"Field '{field_prefix}' not found in template")


def _replace_all_fields(text: str, field_prefix: str, new_value_str: str) -> str:
    """Replace the value on EVERY matching line (for multi-region fields)."""
    lines = text.splitlines(keepends=True)
    pat = re.compile(r'^(\s*' + re.escape(field_prefix) + r'\s*=\s*)\S.*$')
    n_replaced = 0
    for i, line in enumerate(lines):
        m = pat.match(line.rstrip('\n').rstrip('\r'))
        if m:
            nl = '\r\n' if line.endswith('\r\n') else ('\n' if line.endswith('\n') else '')
            lines[i] = f'{m.group(1)}{new_value_str}{nl}'
            n_replaced += 1
    if n_replaced == 0:
        raise ValueError(f"Field '{field_prefix}' not found in template")
    return ''.join(lines)


def write_scan_rhw(template_path: Path, out_path: Path, *,
                   cone_deg: float, spot_cm: float, focus_cm: float,
                   fl_prism: float, burn_on: bool, max_sim_time_s: float) -> None:
    """Read template, modify scan fields, write to out_path."""
    text = template_path.read_text()

    # Laser geometry (beam 1 only — first occurrence of each field)
    text = _replace_first_field(text, 'Half cone angle',  f'{cone_deg:.5e}')
    text = _replace_first_field(text, 'Spot size',        f'{spot_cm:.5e}')
    text = _replace_first_field(text, 'Focus position',   f'{focus_cm:.5e}')

    # Flux limiter — all 4 regions (single value per scan row)
    text = _replace_all_fields(text, 'Flux limiter mult.', f'{fl_prism:.4f}')

    # Burn flags — all 4 regions. For burn OFF, set both to 0 everywhere.
    # (CH Skin region in production templates already has these = 0; leaving
    #  it at 0 is correct.)
    if burn_on:
        # Restore production pattern: DT regions = 1, CH Skin = 0.
        # The simplest robust approach is to leave the template's burn flags
        # alone (they're already correct for production runs).
        pass
    else:
        text = _replace_all_fields(text, 'Fusion reactions on',  '0')
        text = _replace_all_fields(text, 'Fusion transport on',  '0')
        text = _replace_all_fields(text, 'Use alpha deposition', '0')
        text = _replace_all_fields(text, 'Use non alpha deposition', '0')

    # Max simulation time (single occurrence)
    text = _replace_first_field(text, 'Max simulation time', f'{max_sim_time_s:.5e}')

    out_path.write_text(text)


def parse_row(spec: str) -> dict:
    """Parse a row spec like 'label cone=37 spot=0.18 focus=0.22 fl=0.012 burn=off max_t=1.5e-8'."""
    tokens = spec.split()
    if not tokens:
        raise ValueError("Empty row spec")
    label = tokens[0]
    params = {}
    for tok in tokens[1:]:
        if '=' not in tok:
            raise ValueError(f"Bad token (need key=value): {tok}")
        k, v = tok.split('=', 1)
        params[k.strip()] = v.strip()
    # Required keys + defaults
    return dict(
        label=label,
        cone_deg=float(params['cone']),
        spot_cm=float(params['spot']),
        focus_cm=float(params['focus']),
        fl_prism=float(params['fl']),
        burn_on=params.get('burn', 'off').lower() in ('on', 'true', '1'),
        max_sim_time_s=float(params.get('max_t', '1.5e-8')),
    )


def main():
    ap = argparse.ArgumentParser(
        description='Generate PDD scan RHW files from a production template.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--template', required=True, type=Path,
                    help='Production-style RHW template (e.g. fab007 .rhw)')
    ap.add_argument('--row', action='append', required=True, dest='rows',
                    metavar='SPEC',
                    help='Scan row spec: "label cone=37 spot=0.18 focus=0.22 fl=0.012 burn=off max_t=1.5e-8". '
                         'Use --row multiple times for multiple rows.')
    ap.add_argument('--outdir', type=Path, default=Path('.'),
                    help='Output directory (default: current)')
    ap.add_argument('-n', '--dry-run', action='store_true',
                    help='Print what would be written, but do not write')
    args = ap.parse_args()

    if not args.template.exists():
        ap.error(f'Template not found: {args.template}')

    args.outdir.mkdir(parents=True, exist_ok=True)
    print(f'Template:   {args.template}')
    print(f'Output dir: {args.outdir.resolve()}')
    print(f'Rows:       {len(args.rows)}')
    print()
    print(f'{"label":24s} {"cone":>5s} {"spot":>6s} {"focus":>6s} {"FL":>7s} {"burn":>5s} {"max_t [s]":>10s}')
    print('-' * 80)

    for spec in args.rows:
        cfg = parse_row(spec)
        out_path = args.outdir / f"{cfg['label']}.rhw"
        print(f"  {cfg['label']:22s} {cfg['cone_deg']:5.1f} {cfg['spot_cm']:6.3f} {cfg['focus_cm']:6.3f} "
              f"{cfg['fl_prism']:7.4f}  {'on' if cfg['burn_on'] else 'off':>5s}  {cfg['max_sim_time_s']:.3e}")
        if not args.dry_run:
            write_scan_rhw(args.template, out_path,
                           cone_deg=cfg['cone_deg'], spot_cm=cfg['spot_cm'],
                           focus_cm=cfg['focus_cm'], fl_prism=cfg['fl_prism'],
                           burn_on=cfg['burn_on'], max_sim_time_s=cfg['max_sim_time_s'])
    if args.dry_run:
        print('\n(dry-run: nothing written)')


if __name__ == '__main__':
    sys.exit(main())
