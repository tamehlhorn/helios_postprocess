"""
apply_geomcorr_to_rhw.py
========================

Inject the time-dependent geometric coupling correction produced by
``pdd_geometric_coupling.py`` into the laser pulse table of an .rhw file.

For each beam in the ``[Laser Source Data]`` block, this script replaces the
power-table row with the corrected values:

    P_new(t) = P_orig(t) × f_geom(t)

where f_geom(t) is interpolated from the ``geomcorr_power.csv`` (column
``f_geom``) onto the rhw's per-beam time grid.  The correction is applied
identically to all beams since it is geometric and depends only on the cone
half-angle and focal distance, both of which are shared across beams in
typical PDD configurations.

Usage
-----
    python3 apply_geomcorr_to_rhw.py <input.rhw> [options]

Options
    --csv  PATH        Path to <base>_geomcorr_power.csv.  Default: auto-
                       detected as <input>_geomcorr_power.csv where <input>
                       is the rhw stem.
    --out  PATH        Output rhw path.  Default: <input>_geomcorr.rhw
    --dry-run          Print what would change without writing the output.
    --verbose          Print per-beam before/after pulse-energy summaries.

After the patched rhw is written, run Helios with:

    /path/to/Helios -b -i <output.rhw> -d <results_dir> -o <run_name> -x

(``-x`` overwrites the results folder if it already exists.)

This is a *diagnostic* tool: the corrected pulse is a test of the geometric
correction hypothesis, not a physically calibrated input.  Use no-burn runs
for the corrected-pulse test so the implosion trajectory cleanly reflects
the modified drive without alpha-heating amplification.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


# Regex for "Parameters for beam = N" or "Parameters for beam: = N" lines.
# Helios .rhw files use the form "Parameters for beam:        = 1" (with a
# trailing colon after "beam"), so the colon is required-to-tolerate.
BEAM_HEADER_RE  = re.compile(r'Parameters for beam\s*:?\s*=\s*(\d+)', re.IGNORECASE)
TABLE_FORMAT_RE = re.compile(r'\[table format', re.IGNORECASE)
TABLE_ROWS_RE   = re.compile(r'#\s*table rows\s*=\s*(\d+)', re.IGNORECASE)


def load_geomcorr_csv(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load the geomcorr CSV and return (t_ns, f_geom) arrays.

    CSV columns (with # comment header):
        time[ns], R[cm], f_geom, P_orig[TW], P_corr[TW]
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            try:
                vals = [float(v) for v in s.split()]
                if len(vals) >= 3:
                    rows.append(vals)
            except ValueError:
                continue
    if not rows:
        raise RuntimeError(f"No numeric rows found in {path}")
    arr = np.array(rows)
    t_ns   = arr[:, 0]
    f_geom = arr[:, 2]
    return t_ns, f_geom


def is_numeric_row(line: str, n_cols: int) -> bool:
    """True if `line` parses as `n_cols` space-separated floats."""
    s = line.strip()
    if not s or s.startswith('#') or s.startswith('['):
        return False
    try:
        vals = [float(v) for v in s.split()]
        return len(vals) == n_cols
    except ValueError:
        return False


def format_power_row(powers: np.ndarray, template_line: str) -> str:
    """
    Format the corrected power row, attempting to preserve the column-width
    of the original line.  Falls back to a sensible default if the original
    formatting can't be inferred.
    """
    # Try to infer column width from the original line: count whitespace runs
    # and use 6 decimal digits.
    orig_vals = template_line.strip().split()
    n_vals = len(orig_vals)
    if n_vals != len(powers):
        # Mismatch — use safe default
        return '  '.join(f'{p:.6e}' for p in powers) + '\n'

    # Compute width based on longest token in the original
    widths = [len(t) for t in orig_vals]
    avg_width = max(max(widths), 14)

    # Try to preserve the leading indent
    indent = ''
    for c in template_line:
        if c in (' ', '\t'):
            indent += c
        else:
            break

    formatted = indent + ' '.join(f'{p:>{avg_width}.6e}' for p in powers)
    if not formatted.endswith('\n'):
        formatted += '\n'
    return formatted


class BeamPulseTable:
    """A located pulse-table block inside the rhw file."""
    def __init__(self, beam_id: int,
                  time_line_idx: int, power_line_idx: int,
                  times_s: np.ndarray, powers_TW: np.ndarray):
        self.beam_id        = beam_id
        self.time_line_idx  = time_line_idx
        self.power_line_idx = power_line_idx
        self.times_s        = times_s
        self.powers_TW      = powers_TW

    @property
    def n_pts(self) -> int:
        return len(self.times_s)


def locate_pulse_tables(lines: List[str]) -> List[BeamPulseTable]:
    """
    Scan an .rhw file (already split into lines) and return one
    ``BeamPulseTable`` per beam whose power table was successfully parsed.

    The state machine follows the format used by helios_postprocess.rhw_parser
    (see _parse_laser_geometry):

        [Laser Source Data]
            ...
            Parameters for beam = N
                ...
                [table format ...]
                # table rows = M
                <comment / header lines>
                t1 t2 ... tM        <- time row (seconds)
                P1 P2 ... PM        <- power row (TW)
                ... possibly more rows (unused) ...
            Parameters for beam = N+1
            ...
        [End Laser Source Data]
    """
    tables: List[BeamPulseTable] = []
    in_laser   = False
    beam_id    = None
    in_table   = False
    n_rows     = 0
    rows_seen  = 0
    time_idx   = None
    time_vals  = None

    for i, raw in enumerate(lines):
        line = raw.rstrip('\n')
        if '[Laser Source Data]' in line:
            in_laser = True; continue
        if '[End Laser Source Data]' in line:
            in_laser = False; break
        if not in_laser:
            continue

        m = BEAM_HEADER_RE.search(line)
        if m:
            beam_id = int(m.group(1))
            in_table = False
            n_rows = 0
            rows_seen = 0
            time_idx = None
            time_vals = None
            continue

        if TABLE_FORMAT_RE.search(line):
            in_table = True
            n_rows = 0
            rows_seen = 0
            time_idx = None
            time_vals = None
            continue

        if not in_table or beam_id is None:
            continue

        # Number-of-rows comment
        m2 = TABLE_ROWS_RE.search(line)
        if m2:
            n_rows = int(m2.group(1))
            continue

        # Numeric data row
        if n_rows > 0 and is_numeric_row(line, n_rows):
            vals = np.array([float(v) for v in line.strip().split()])
            if rows_seen == 0:
                time_idx  = i
                time_vals = vals.copy()
                rows_seen += 1
            elif rows_seen == 1:
                tables.append(BeamPulseTable(
                    beam_id=beam_id,
                    time_line_idx=time_idx,
                    power_line_idx=i,
                    times_s=time_vals,
                    powers_TW=vals,
                ))
                in_table = False
                rows_seen += 1
            continue

    return tables


def apply_correction(table: BeamPulseTable,
                      csv_t_ns: np.ndarray,
                      csv_f_geom: np.ndarray) -> np.ndarray:
    """
    Compute P_new = P_orig × f_geom(t).  f_geom is interpolated onto the
    rhw time grid; samples before the first CSV time get f_geom = 1
    (no correction), and samples beyond the last CSV time are held at
    the CSV's last value (deep post-laser tail).
    """
    t_rhw_ns = table.times_s * 1e9
    f_at_t = np.interp(t_rhw_ns,
                        csv_t_ns,
                        csv_f_geom,
                        left=1.0,
                        right=float(csv_f_geom[-1]),)
    f_at_t = np.clip(f_at_t, 0.0, 1.0)
    return table.powers_TW * f_at_t


def trapz(y, x):
    """Trapezoid integration with numpy 1.x / 2.x compatibility."""
    fn = getattr(np, 'trapezoid', None) or np.trapz
    return fn(y, x)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Inject geomcorr_power.csv into an .rhw laser pulse table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('input_rhw',
                        help='Path to the input .rhw (typically the rhw of the '
                             'no-burn run that produced the geomcorr CSV).')
    parser.add_argument('--csv', default=None,
                        help='Path to <base>_geomcorr_power.csv.  Default: '
                             '<input_stem>_geomcorr_power.csv in same dir.')
    parser.add_argument('--out', default=None,
                        help='Output .rhw path.  Default: <input_stem>_geomcorr.rhw')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print what would change without writing the output.')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    in_rhw = Path(args.input_rhw).expanduser().resolve()
    if not in_rhw.exists():
        print(f"ERROR: input rhw not found: {in_rhw}", file=sys.stderr)
        return 1

    if args.csv:
        csv_path = Path(args.csv).expanduser().resolve()
    else:
        csv_path = in_rhw.with_name(in_rhw.stem + "_geomcorr_power.csv")
    if not csv_path.exists():
        print(f"ERROR: geomcorr CSV not found: {csv_path}", file=sys.stderr)
        print("       Pass --csv to point at it explicitly, or run "
              "pdd_geometric_coupling.py first.", file=sys.stderr)
        return 1

    if args.out:
        out_rhw = Path(args.out).expanduser().resolve()
    else:
        out_rhw = in_rhw.with_name(in_rhw.stem + "_geomcorr.rhw")

    # ── Load CSV
    csv_t_ns, csv_f_geom = load_geomcorr_csv(csv_path)
    print(f"Loaded f_geom(t) from  {csv_path.name}")
    print(f"  Samples:     {len(csv_t_ns)}")
    print(f"  Time range:  {csv_t_ns[0]:.3f} – {csv_t_ns[-1]:.3f} ns")
    print(f"  f_geom min:  {float(np.nanmin(csv_f_geom)):.4f}")
    print(f"  f_geom max:  {float(np.nanmax(csv_f_geom)):.4f}")
    print()

    # ── Load rhw lines
    with open(in_rhw) as fh:
        lines = fh.readlines()

    # ── Locate pulse tables
    tables = locate_pulse_tables(lines)
    if not tables:
        print("ERROR: no laser pulse tables located in the input rhw.",
              file=sys.stderr)
        print("       Dumping bracketed section markers and the first 30 lines",
              file=sys.stderr)
        print("       so we can diagnose the format mismatch:", file=sys.stderr)
        print("", file=sys.stderr)
        print(f"--- All bracketed [...] section markers in {in_rhw.name} ---",
              file=sys.stderr)
        section_re = re.compile(r'^\s*(\[[^\]]+\])')
        for i, raw in enumerate(lines):
            m = section_re.match(raw)
            if m:
                print(f"  line {i+1:5d}: {m.group(1)}", file=sys.stderr)
        print("", file=sys.stderr)
        # If we didn't find a [Laser Source Data] marker, that's the problem
        any_laser = any('[Laser Source Data]' in ln for ln in lines)
        any_beam  = any(BEAM_HEADER_RE.search(ln) for ln in lines)
        any_tab   = any(TABLE_FORMAT_RE.search(ln) for ln in lines)
        any_rows  = any(TABLE_ROWS_RE.search(ln) for ln in lines)
        print("--- Marker-presence summary ---", file=sys.stderr)
        print(f"  '[Laser Source Data]' substring : {'FOUND' if any_laser else 'NOT FOUND'}",
              file=sys.stderr)
        print(f"  'Parameters for beam = N'       : {'FOUND' if any_beam  else 'NOT FOUND'}",
              file=sys.stderr)
        print(f"  '[table format'                 : {'FOUND' if any_tab   else 'NOT FOUND'}",
              file=sys.stderr)
        print(f"  '# table rows = N'              : {'FOUND' if any_rows  else 'NOT FOUND'}",
              file=sys.stderr)
        print("", file=sys.stderr)
        print("Tell Claude which markers were NOT FOUND and paste a 20-line "
              "snippet of the rhw around the laser-pulse section.  The state-"
              "machine markers are configurable.", file=sys.stderr)
        return 1

    print(f"Located {len(tables)} pulse table(s) in {in_rhw.name}:")
    for tab in tables:
        t_ns = tab.times_s * 1e9
        on   = tab.powers_TW > 0
        t_on = t_ns[on]
        E_kJ = float(trapz(tab.powers_TW * 1e12, tab.times_s)) * 1e-3
        print(f"  beam {tab.beam_id}: {tab.n_pts} points, "
              f"{t_on[0]:.2f}–{t_on[-1]:.2f} ns on, "
              f"peak {tab.powers_TW.max():.1f} TW, "
              f"∫P dt = {E_kJ:.1f} kJ")
    print()

    # ── Apply correction beam-by-beam, build replacement-line map
    replacements: dict[int, str] = {}
    total_E_orig = 0.0
    total_E_corr = 0.0
    for tab in tables:
        p_new = apply_correction(tab, csv_t_ns, csv_f_geom)
        new_line = format_power_row(p_new, lines[tab.power_line_idx])
        replacements[tab.power_line_idx] = new_line

        E_orig = float(trapz(tab.powers_TW * 1e12, tab.times_s)) * 1e-3
        E_corr = float(trapz(p_new          * 1e12, tab.times_s)) * 1e-3
        total_E_orig += E_orig
        total_E_corr += E_corr

        on = tab.powers_TW > 0
        f_at_on = (p_new[on] / tab.powers_TW[on])
        if args.verbose:
            print(f"  beam {tab.beam_id}: "
                  f"∫P_orig = {E_orig:7.2f} kJ → "
                  f"∫P_corr = {E_corr:7.2f} kJ  "
                  f"(Δ {100*(E_corr/E_orig - 1):+5.2f} %; "
                  f"mean f_geom on-pulse = {float(f_at_on.mean()):.3f})")

    print(f"Total laser energy: "
          f"{total_E_orig:.2f} kJ → {total_E_corr:.2f} kJ "
          f"(Δ {100*(total_E_corr/total_E_orig - 1):+5.2f} %)")
    print()

    # ── Emit corrected rhw
    new_lines = list(lines)
    for idx, new in replacements.items():
        new_lines[idx] = new

    if args.dry_run:
        print(f"--dry-run: would write {out_rhw}")
        for tab in tables:
            print(f"  beam {tab.beam_id} power row at line {tab.power_line_idx+1}:")
            print(f"    OLD: {lines[tab.power_line_idx].rstrip()}")
            print(f"    NEW: {replacements[tab.power_line_idx].rstrip()}")
        return 0

    out_rhw.parent.mkdir(parents=True, exist_ok=True)
    with open(out_rhw, 'w') as fh:
        fh.writelines(new_lines)

    print(f"Wrote {out_rhw}")
    print()
    print("To run Helios on the corrected pulse:")
    print(f"  /path/to/Helios -b -i {out_rhw} \\")
    print(f"      -d /path/to/results_root -o {in_rhw.stem}_geomcorr -x")
    print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
