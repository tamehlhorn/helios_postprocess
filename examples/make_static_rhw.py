#!/usr/bin/env python3
"""
make_static_rhw.py
==================

Generate Helios RHW input files for zero-D static (hydro-off) verification
tests at varying T and/or rho, starting from an existing template RHW.

Edits exactly three fields:
    Temperature            = <T_keV>      (line ~33; ion T)
    Electron Temperature   = <T_keV>      (line ~37; electron T)
    Density                = <rho_gcc>    (line ~49; g/cc)

All other settings (zoning, EOS path, hydro off, fusion-reactions-on,
output cadence, etc.) are preserved from the template.

Usage examples
--------------
    # Single file, change T only
    python3 make_static_rhw.py --template DT_static_3kev.rhw \\
        --T 10 --outdir .

    # Sweep T at fixed rho (default rho from template)
    python3 make_static_rhw.py --template DT_static_3kev.rhw \\
        --T 5 7 10 15 --outdir .

    # Density sweep at fixed T = 10 keV
    python3 make_static_rhw.py --template DT_static_3kev.rhw \\
        --T 10 --rho 0.418 4.18 41.8 \\
        --tag-format "DT_static_10keV_n{n_DT_label}" --outdir .

    # Both sweeps with foam template
    python3 make_static_rhw.py --template DT_static_WF_3kev.rhw \\
        --T 5 7 10 15 --prefix foam_static --outdir .

Output naming: by default, "<prefix>_<T>keV.rhw" (or with --tag-format).
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path
from typing import Iterable


AMU_G = 1.66053906660e-24
M_AVG_DT_AMU = (2.014 + 3.016) / 2.0   # equimolar DT


def rho_for_nDT(n_DT_cm3: float) -> float:
    """Mass density of pure-DT plasma at given total D+T number density."""
    return n_DT_cm3 * M_AVG_DT_AMU * AMU_G


def nDT_for_rho(rho_gcc: float) -> float:
    """Inverse: total DT number density for a given mass density."""
    return rho_gcc / (M_AVG_DT_AMU * AMU_G)


def _replace_field(text: str, field_prefix: str, new_value_str: str) -> str:
    """
    Replace the value on the FIRST line whose stripped form starts with
    `field_prefix` followed by '='.  Preserves the leading whitespace and
    the column-alignment of '='.

    Example field_prefix: 'Temperature', 'Electron Temperature', 'Density'
    """
    lines = text.splitlines(keepends=True)
    pat_str = r'^(\s*' + re.escape(field_prefix) + r'\s*=\s*)\S.*$'
    pat = re.compile(pat_str)
    for i, line in enumerate(lines):
        m = pat.match(line.rstrip('\n').rstrip('\r'))
        if m:
            # Preserve newline character(s)
            nl = ''
            if line.endswith('\r\n'):
                nl = '\r\n'
            elif line.endswith('\n'):
                nl = '\n'
            lines[i] = f'{m.group(1)}{new_value_str}{nl}'
            return ''.join(lines)
    raise ValueError(f"Field '{field_prefix} =' not found in template")


def write_static_rhw(template_path: Path, out_path: Path,
                     T_keV: float, rho_gcc: float,
                     max_sim_time_s: float = 1.0e-9) -> None:
    """Read template, replace the relevant fields, write to out_path.

    Replaces:
        Temperature, Electron Temperature, Density,
        Max simulation time  (defaults to 1 ns for verification snapshots --
        the template inherits 1e-4 s = 100 us, which can cause multi-minute
        Helios runs at higher T even with hydro off).
    """
    text = template_path.read_text()
    # Order matters: 'Electron Temperature' must be matched before 'Temperature'
    # because the latter is a prefix of the former. The regex is anchored to
    # the LINE START (with whitespace) so this is safe -- but to be defensive
    # we replace 'Electron Temperature' first, then 'Temperature'.
    text = _replace_field(text, 'Electron Temperature', f'{T_keV:g}')
    text = _replace_field(text, 'Temperature',          f'{T_keV:g}')
    text = _replace_field(text, 'Density',              f'{rho_gcc:g}')
    text = _replace_field(text, 'Max simulation time',  f'{max_sim_time_s:.5e}')
    out_path.write_text(text)


def main():
    ap = argparse.ArgumentParser(description='Generate zero-D static RHW files',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--template', required=True, type=Path,
                    help='Source RHW template (e.g. DT_static_3kev.rhw)')
    ap.add_argument('--T', type=float, nargs='+', required=True,
                    help='Temperature(s) in keV. Multiple values -> multiple files.')
    ap.add_argument('--rho', type=float, nargs='+', default=None,
                    help='Mass density(s) in g/cc.  Multiple values -> multiple files. '
                         'If omitted, density is held at template value (0.0418 g/cc for the DT 3keV template).')
    ap.add_argument('--prefix', default=None,
                    help='Output filename prefix. Default: inferred from template stem.')
    ap.add_argument('--tag-format', default=None,
                    help='Override default filename format. Available substitutions: '
                         '{T_keV}, {T_int}, {rho_gcc}, {n_DT_label}. Example: '
                         '"DT_static_{T_int}keV_n{n_DT_label}"')
    ap.add_argument('--outdir', type=Path, default=Path('.'),
                    help='Output directory (default: current)')
    ap.add_argument('--max-sim-time', type=float, default=1.0e-9,
                    help='Helios "Max simulation time" in seconds (default: 1e-9 = 1 ns). '
                         'The template default of 1e-4 s is far too long for zero-D rate '
                         'extraction and can produce multi-minute runs at high T.')
    ap.add_argument('-n', '--dry-run', action='store_true',
                    help='Print what would be written, but do not write')
    args = ap.parse_args()

    if not args.template.exists():
        ap.error(f'Template not found: {args.template}')

    # Determine the template's default density (for the case when --rho not given)
    template_text = args.template.read_text()
    m = re.search(r'^\s*Density\s*=\s*(\S+)', template_text, flags=re.M)
    template_rho = float(m.group(1)) if m else 0.0418

    rhos = args.rho if args.rho else [template_rho]
    Ts   = args.T

    # Decide naming
    default_prefix = args.template.stem.replace('_3kev', '').replace('_3keV', '')
    prefix = args.prefix or default_prefix
    if args.tag_format:
        name_fmt = args.tag_format
    elif len(rhos) > 1:
        name_fmt = prefix + '_{T_int}keV_n{n_DT_label}'
    else:
        name_fmt = prefix + '_{T_int}keV'

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f'Template: {args.template}  (default rho = {template_rho:g} g/cc)')
    print(f'Combinations: {len(Ts)} x {len(rhos)} = {len(Ts)*len(rhos)} runs')
    print(f'Output dir: {args.outdir}')
    print()

    for T_keV in Ts:
        for rho_gcc in rhos:
            n_DT = nDT_for_rho(rho_gcc)
            # Build a short n_DT label like 1e22, 1e23, ...
            n_DT_label = f'{n_DT:.0e}'.replace('+0', '').replace('+', '').replace('e0', 'e')
            tag = name_fmt.format(
                T_keV=T_keV, T_int=int(round(T_keV)),
                rho_gcc=rho_gcc, n_DT_label=n_DT_label)
            out_path = args.outdir / f'{tag}.rhw'
            print(f'  {out_path.name}:  T = {T_keV:g} keV   rho = {rho_gcc:g} g/cc   n_DT = {n_DT:.3e} cm^-3')
            if not args.dry_run:
                write_static_rhw(args.template, out_path, T_keV, rho_gcc,
                                 max_sim_time_s=args.max_sim_time)

    if args.dry_run:
        print('\n(dry-run: nothing written)')


if __name__ == '__main__':
    sys.exit(main())
