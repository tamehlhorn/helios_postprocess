"""
apply_foot_beam_split.py
========================

Restructure a single-beam Olson PDD .rhw into a 2-beam configuration to
test whether geometric defocus on the FOOT PHASE ONLY (with the
ramp+peak phase held at the calibrated fab007 geometry) slows the foot
shock toward LILAC's 7.5 ns gas/ice arrival.

Motivation
----------
The fab007 calibration uses time-averaged geometric defocus (cone 37°,
spot 0.18, focus 0.22) to mimic the 3D beam-target miss that Helios's 1D
ray-trace can't reproduce.  In reality the foot-phase miss is smaller
(capsule still large) and the peak-phase miss is larger (capsule
shrunk).  The fab007 geometry is therefore a compromise — too defocused
in the foot, not defocused enough at peak.

If we apply EXTRA defocus only on the foot beam — same incident foot
power (25 TW), but spread over more solid angle and a larger spot — the
foot shock launches with lower local ablation pressure, propagates more
slowly through the ablator + ice, and arrives later at the gas/ice
interface.  Total foot-phase delivered energy stays the same; only the
intensity-at-corona is reduced.  This is the GEOMETRIC analogue of
simply dropping foot power, and is more physically motivated because it
mimics what 3D beam-miss actually does (same incident power, less
reaches the capsule per unit area).

The script restructures the rhw so:

  Beam 1 = foot only (0 → t_split) at OVERRIDE geometry
  Beam 2 = ramp + peak (t_split → end) at ORIGINAL (calibrated) geometry

Total delivered energy is preserved to better than 0.1% (numerical
resampling error). Each beam's pulse table is resampled to ``n_points``
uniformly across the full original time range, with the inactive phase
filled in with zeros.

Usage
-----
    python3 apply_foot_beam_split.py <input.rhw> \\
        --t-split  5.0 \\
        --foot-cone 42 \\
        --foot-spot 0.22 \\
        [--foot-focus 0.22] \\
        [--out PATH] \\
        [--n-points 200] \\
        [--dry-run] [--verbose]

Times in nanoseconds at the CLI; converted to seconds in the rhw.
If a geometry override (--foot-cone / --foot-spot / --foot-focus) is
omitted, that parameter is held at the original beam-1 value.  Omitting
ALL three overrides produces a pure pulse-split (useful as a sanity
check that splitting alone doesn't perturb the implosion).

After the patched rhw is written, run Helios with:

    /path/to/Helios -b -i <output.rhw> -d <results_dir> -o <run_name> -x

Use no-burn runs (Burn flag OFF) for this study — the foot/ramp/peak
shock arrivals at gas/ice are diagnostic of pulse-shape physics, not
of alpha-coupling, and no-burn runs land in ~30 min on Mac Studio.

Sister tool: ``apply_geomcorr_to_rhw.py`` applies a TIME-DEPENDENT power
multiplier to all beams uniformly.  This tool MAKES IT A 2-BEAM SETUP
so the foot phase can be perturbed independently of the rest of the
pulse — a different lever, complementary in principle.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


# ── rhw markup ──────────────────────────────────────────────────────────
BEAM_HEADER_RE  = re.compile(r'Parameters for beam\s*:?\s*=\s*(\d+)', re.IGNORECASE)
N_BEAMS_RE      = re.compile(r'(Number of laser beams\s*=\s*)(\d+)', re.IGNORECASE)
LASER_START_RE  = re.compile(r'\[Laser Source Data\]', re.IGNORECASE)
LASER_END_RE    = re.compile(r'\[End Laser Source Data\]', re.IGNORECASE)
TABLE_FORMAT_RE = re.compile(r'\[table format', re.IGNORECASE)
TABLE_ROWS_RE   = re.compile(r'#\s*table rows\s*=\s*(\d+)', re.IGNORECASE)
HALF_CONE_RE    = re.compile(r'(Half cone angle\s*=\s*)(\S+)', re.IGNORECASE)
SPOT_SIZE_RE    = re.compile(r'(Spot size\s*=\s*)(\S+)', re.IGNORECASE)
FOCUS_POS_RE    = re.compile(r'(Focus position\s*=\s*)(\S+)', re.IGNORECASE)


# ── helpers ─────────────────────────────────────────────────────────────
def trapz(y, x):
    """Trapezoid integration with numpy 1.x / 2.x compatibility."""
    fn = getattr(np, 'trapezoid', None) or np.trapz
    return fn(y, x)


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


def find_n_beams_line(lines: List[str]) -> Tuple[Optional[int], Optional[int]]:
    """Locate the `Number of laser beams = N` line inside [Laser Source Data].

    Returns (line_idx, n_beams) or (None, None) if not found.
    """
    in_laser = False
    for i, line in enumerate(lines):
        if LASER_START_RE.search(line):
            in_laser = True
            continue
        if LASER_END_RE.search(line):
            in_laser = False
            continue
        if not in_laser:
            continue
        m = N_BEAMS_RE.search(line)
        if m:
            return i, int(m.group(2))
    return None, None


def find_beam_block(lines: List[str], beam_id: int) -> Tuple[Optional[int], Optional[int]]:
    """Return (start_idx, end_idx) bounding beam `beam_id`'s parameter block
    within [Laser Source Data].

    `end_idx` is exclusive: it points to either the next `Parameters for
    beam:` header line, or the `[End Laser Source Data]` line.
    """
    in_laser = False
    in_target = False
    start = None
    for i, line in enumerate(lines):
        if LASER_START_RE.search(line):
            in_laser = True
            continue
        if LASER_END_RE.search(line):
            if in_target:
                return start, i
            return None, None
        if not in_laser:
            continue
        m = BEAM_HEADER_RE.search(line)
        if m:
            this_id = int(m.group(1))
            if this_id == beam_id:
                in_target = True
                start = i
            elif in_target:
                return start, i
    return None, None


def find_all_beam_blocks(lines: List[str]
                         ) -> Tuple[List[Tuple[int, int, int]], Optional[int]]:
    """Enumerate every beam block within [Laser Source Data].

    Returns (blocks, laser_end_idx) where blocks is a list of
    (beam_id, start_idx, end_idx) tuples and laser_end_idx is the index
    of the `[End Laser Source Data]` line (None if not found).
    """
    blocks: List[Tuple[int, int, int]] = []
    in_laser = False
    cur_id = None
    cur_start = None
    laser_end = None
    for i, line in enumerate(lines):
        if LASER_START_RE.search(line):
            in_laser = True
            continue
        if LASER_END_RE.search(line):
            if cur_id is not None:
                blocks.append((cur_id, cur_start, i))
                cur_id = None
            laser_end = i
            break
        if not in_laser:
            continue
        m = BEAM_HEADER_RE.search(line)
        if m:
            if cur_id is not None:
                blocks.append((cur_id, cur_start, i))
            cur_id = int(m.group(1))
            cur_start = i
    return blocks, laser_end


def is_pulse_active(lines: List[str], start: int, end: int) -> bool:
    """True if the beam block has a pulse table with any non-zero power
    entry.  Helios's NIF PDD convention uses 3 beams with beams 2 and 3
    as placeholders (all-zero pulse tables) and `Laser power model is on
    = 0`; only beam 1 is physically driving.  This helper distinguishes
    a real-active beam from a placeholder.
    """
    rows_idx, time_idx, power_idx, times_s, powers_TW = locate_pulse_table(
        lines, start, end)
    if powers_TW is None:
        return False
    return bool(np.any(powers_TW > 0))


def locate_pulse_table(lines: List[str], start: int, end: int
                       ) -> Tuple[Optional[int], Optional[int], Optional[int],
                                  Optional[np.ndarray], Optional[np.ndarray]]:
    """Locate the pulse table within lines[start:end].

    Returns (rows_meta_idx, time_idx, power_idx, times_s, powers_TW) or
    a tuple of Nones if no table was located.

    The rhw pulse table layout is:
        [table format=2]:    Time-dependent laser powers:
            table is 3D  = 0
            # table rows = K   <-- rows_meta_idx
            # table cols = 2
            t1 t2 ... tK       <-- time_idx  (seconds)
            P1 P2 ... PK       <-- power_idx (TW)
    """
    in_table = False
    n_rows = 0
    rows_seen = 0
    rows_meta_idx = None
    time_idx = None
    time_vals = None

    for i in range(start, end):
        line = lines[i]
        if TABLE_FORMAT_RE.search(line):
            in_table = True
            n_rows = 0
            rows_seen = 0
            rows_meta_idx = None
            time_idx = None
            time_vals = None
            continue
        if not in_table:
            continue
        m = TABLE_ROWS_RE.search(line)
        if m:
            n_rows = int(m.group(1))
            rows_meta_idx = i
            continue
        if n_rows > 0 and is_numeric_row(line, n_rows):
            vals = np.array([float(v) for v in line.strip().split()])
            if rows_seen == 0:
                time_idx = i
                time_vals = vals.copy()
                rows_seen += 1
            else:
                return rows_meta_idx, time_idx, i, time_vals, vals
    return None, None, None, None, None


def format_data_row(values: np.ndarray, template_line: str, prec: int = 6) -> str:
    """Format a numeric data row, preserving the leading indent of
    `template_line` and emitting values in fixed-width scientific notation.
    """
    indent = ''
    for c in template_line:
        if c in (' ', '\t'):
            indent += c
        else:
            break
    width = prec + 8   # " 3.290000e+02" = 13 chars for prec=6
    body = ' '.join(f'{v:>{width}.{prec}e}' for v in values)
    return indent + body + '\n'


def format_table_rows_line(template_line: str, n_new: int) -> str:
    """Replace the integer in a `# table rows = N` line with `n_new`."""
    new = re.sub(r'(table rows\s*=\s*)\d+',
                 lambda m: f'{m.group(1)}{n_new}',
                 template_line, count=1, flags=re.IGNORECASE)
    if not new.endswith('\n'):
        new += '\n'
    return new


def replace_value_after_eq(line: str, regex: re.Pattern, new_value: str) -> str:
    """If `regex` matches `line` (capturing the prefix + current value), rewrite
    the value while preserving everything before and after.  Returns the new
    line, or the original if no match.

    Regex must have two groups: group(1) = "...= " prefix; group(2) = current
    numeric value as captured by `\\S+`.
    """
    m = regex.search(line)
    if not m:
        return line
    return line[:m.start(2)] + new_value + line[m.end(2):]


def split_pulse(times_s: np.ndarray, powers_TW: np.ndarray,
                t_split_s: float, n_each: int
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resample the original pulse onto a fine uniform grid spanning the full
    time range, and return (t_grid, P_beam1, P_beam2).

    Beam 1 carries the pulse only for t <= t_split; beam 2 carries it only
    for t >= t_split.  Both beams span the full original time range, with
    the inactive phase filled with zeros so each beam's table covers the
    same temporal extent (cleaner for Helios's table interpolation than
    truncating).

    Energy budget: trapz(P_beam1 + P_beam2) reconstructs the original pulse
    energy within numerical interpolation error (<0.1% at n_each=200).
    """
    t_min, t_max = float(times_s.min()), float(times_s.max())
    t_grid = np.linspace(t_min, t_max, n_each)
    P_orig = np.interp(t_grid, times_s, powers_TW)
    P_b1   = np.where(t_grid <= t_split_s, P_orig, 0.0)
    P_b2   = np.where(t_grid >= t_split_s, P_orig, 0.0)
    return t_grid, P_b1, P_b2


# ── main ────────────────────────────────────────────────────────────────
def main() -> int:
    parser = argparse.ArgumentParser(
        description="Split a single-beam .rhw into 2 beams: foot-only "
                    "beam 1 (optionally with override geometry) + "
                    "ramp+peak beam 2 (original geometry preserved).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('input_rhw',
                        help='Path to the single-beam input .rhw.')
    parser.add_argument('--t-split', type=float, default=5.0,
                        help='Split time between foot beam and ramp+peak beam, '
                             'in nanoseconds.  Default 5.0 ns matches the '
                             'Olson PDD foot/ramp boundary.')
    parser.add_argument('--foot-cone', type=float, default=None,
                        help='Override beam-1 (foot) half-cone angle, in '
                             'degrees.  Default: keep original.')
    parser.add_argument('--foot-spot', type=float, default=None,
                        help='Override beam-1 (foot) spot size, in cm.  '
                             'Default: keep original.')
    parser.add_argument('--foot-focus', type=float, default=None,
                        help='Override beam-1 (foot) focus position, in cm.  '
                             'Default: keep original.')
    parser.add_argument('--out', default=None,
                        help='Output .rhw path.  Default: '
                             '<input>_footsplit_c<cone>_s<spot>.rhw')
    parser.add_argument('--n-points', type=int, default=200,
                        help='Points per beam in the resampled pulse table.  '
                             '200 ≈ 70 ps resolution over a 14 ns pulse.')
    parser.add_argument('--dry-run', action='store_true',
                        help='Report what would change without writing.')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    in_rhw = Path(args.input_rhw).expanduser().resolve()
    if not in_rhw.exists():
        print(f"ERROR: input rhw not found: {in_rhw}", file=sys.stderr)
        return 1

    # ── Default output filename
    if args.out:
        out_rhw = Path(args.out).expanduser().resolve()
    else:
        suffix_bits = ['footsplit']
        if args.foot_cone is not None:
            suffix_bits.append(f'c{int(round(args.foot_cone))}')
        if args.foot_spot is not None:
            suffix_bits.append(f's{int(round(args.foot_spot * 100)):02d}')
        if args.foot_focus is not None:
            suffix_bits.append(f'd{int(round(args.foot_focus * 100)):02d}')
        suffix = '_'.join(suffix_bits)
        out_rhw = in_rhw.with_name(in_rhw.stem + f'_{suffix}.rhw')

    # ── Load rhw
    with open(in_rhw) as fh:
        lines = fh.readlines()

    # ── Locate `Number of laser beams = N`
    n_beams_idx, n_beams = find_n_beams_line(lines)
    if n_beams_idx is None:
        print("ERROR: could not locate `Number of laser beams = N` line "
              "inside [Laser Source Data].", file=sys.stderr)
        return 1

    # ── Enumerate beam blocks and identify the ACTIVE beam
    # PDD targets typically declare 3 beams (NIF PDD convention) with
    # beams 2 and 3 as zero-pulse placeholders.  We split the single
    # active beam, leave the placeholders alone, and append the new
    # ramp+peak beam at the end of the laser data block.
    beam_blocks, laser_end_idx = find_all_beam_blocks(lines)
    if not beam_blocks:
        print("ERROR: no beam blocks found inside [Laser Source Data].",
              file=sys.stderr)
        return 1
    if laser_end_idx is None:
        print("ERROR: could not locate [End Laser Source Data] marker.",
              file=sys.stderr)
        return 1

    active = [(bid, s, e) for (bid, s, e) in beam_blocks
              if is_pulse_active(lines, s, e)]

    if len(active) == 0:
        print("ERROR: input rhw has no beam with a non-zero pulse table.",
              file=sys.stderr)
        return 1
    if len(active) > 1:
        active_ids = ', '.join(str(bid) for bid, _, _ in active)
        print(f"ERROR: input rhw has multiple active beams (IDs {active_ids}). "
              f"This script splits a SINGLE active beam into foot + ramp+peak. "
              f"Re-run on a single-active-beam calibration rhw.",
              file=sys.stderr)
        return 1

    active_id, b1_start, b1_end = active[0]
    new_beam_id = max(bid for bid, _, _ in beam_blocks) + 1

    if n_beams > 1:
        n_inactive = n_beams - 1
        print(f"Input has {n_beams} declared beams; beam {active_id} is the "
              f"only active one (the other {n_inactive} are zero-pulse "
              f"placeholders).  Splitting beam {active_id}; new ramp+peak "
              f"beam will be added as beam {new_beam_id}.")
        print()

    # ── Locate pulse table within beam 1
    rows_idx, time_idx, power_idx, times_s, powers_TW = locate_pulse_table(
        lines, b1_start, b1_end)
    if rows_idx is None:
        print("ERROR: could not locate pulse table inside beam 1 block.",
              file=sys.stderr)
        return 1

    t_min, t_max = float(times_s.min()), float(times_s.max())
    t_split_s = args.t_split * 1e-9
    if not (t_min < t_split_s < t_max):
        print(f"ERROR: --t-split {args.t_split} ns is outside the pulse "
              f"range [{t_min*1e9:.3f}, {t_max*1e9:.3f}] ns.", file=sys.stderr)
        return 1

    # ── Split the pulse
    t_grid, P_beam1, P_beam2 = split_pulse(times_s, powers_TW,
                                            t_split_s, args.n_points)

    # ── Energy book-keeping
    E_orig   = float(trapz(powers_TW * 1e12, times_s)) * 1e-3   # kJ
    E_b1     = float(trapz(P_beam1   * 1e12, t_grid))  * 1e-3
    E_b2     = float(trapz(P_beam2   * 1e12, t_grid))  * 1e-3
    E_sum    = E_b1 + E_b2

    print(f"Input rhw:  {in_rhw}")
    print(f"Output:     {out_rhw}")
    print()
    print(f"Original pulse: {E_orig:.2f} kJ over "
          f"{t_min*1e9:.2f}–{t_max*1e9:.2f} ns "
          f"(peak {powers_TW.max():.1f} TW)")
    print(f"Split at t = {args.t_split:.3f} ns "
          f"(P_orig(t_split) ≈ "
          f"{float(np.interp(t_split_s, times_s, powers_TW)):.2f} TW):")
    print(f"  Beam 1 (foot, t ≤ t_split):        {E_b1:7.2f} kJ "
          f"({100*E_b1/E_orig:5.1f}%)")
    print(f"  Beam 2 (ramp+peak, t ≥ t_split):   {E_b2:7.2f} kJ "
          f"({100*E_b2/E_orig:5.1f}%)")
    print(f"  Sum (vs original):                  {E_sum:7.2f} kJ "
          f"(Δ = {100*(E_sum/E_orig - 1):+5.3f}%)")
    print()

    if args.foot_cone is None and args.foot_spot is None and args.foot_focus is None:
        print("WARNING: no foot-geometry overrides specified.  Output rhw "
              "will be a pure 2-beam pulse split with identical geometry "
              "on both beams — useful as a regression check that the split "
              "itself doesn't perturb the implosion, but the early-shock "
              "study needs at least one of --foot-cone / --foot-spot / "
              "--foot-focus to be specified to actually defocus the foot.")
        print()

    # ── Build modified beam-1 lines (geometry overrides + foot pulse)
    new_lines = list(lines)

    for i in range(b1_start, b1_end):
        orig_line = lines[i]
        edited = orig_line
        if args.foot_cone is not None and HALF_CONE_RE.search(orig_line):
            edited = replace_value_after_eq(
                orig_line, HALF_CONE_RE, f'{args.foot_cone}')
        elif args.foot_spot is not None and SPOT_SIZE_RE.search(orig_line):
            edited = replace_value_after_eq(
                orig_line, SPOT_SIZE_RE, f'{args.foot_spot}')
        elif args.foot_focus is not None and FOCUS_POS_RE.search(orig_line):
            edited = replace_value_after_eq(
                orig_line, FOCUS_POS_RE, f'{args.foot_focus}')
        new_lines[i] = edited

    # Beam 1 pulse table → foot only
    new_lines[rows_idx]  = format_table_rows_line(lines[rows_idx], args.n_points)
    new_lines[time_idx]  = format_data_row(t_grid,  lines[time_idx])
    new_lines[power_idx] = format_data_row(P_beam1, lines[power_idx])

    # ── Build the new ramp+peak beam (clone of the ORIGINAL active-beam
    # block, then override header ID + pulse table).  Geometry is preserved
    # from the original active beam since beam N+1 inherits the calibrated
    # (e.g. fab007) ramp+peak geometry.
    b2_block = list(lines[b1_start:b1_end])

    # Header line: "Parameters for beam = <active_id>" → "= <new_beam_id>"
    header_replaced = False
    for j, line in enumerate(b2_block):
        if BEAM_HEADER_RE.search(line):
            b2_block[j] = re.sub(
                r'(Parameters for beam\s*:?\s*=\s*)\d+',
                f'\\g<1>{new_beam_id}', line, count=1, flags=re.IGNORECASE)
            if not b2_block[j].endswith('\n'):
                b2_block[j] += '\n'
            header_replaced = True
            break
    if not header_replaced:
        print("ERROR: could not locate beam header in the cloned new-beam "
              "block.", file=sys.stderr)
        return 1

    # New beam pulse table → ramp+peak (uses ORIGINAL template lines for
    # indent/formatting, intentionally — the template was cloned from the
    # unmodified active-beam block).
    b2_rows_rel  = rows_idx  - b1_start
    b2_time_rel  = time_idx  - b1_start
    b2_power_rel = power_idx - b1_start
    b2_block[b2_rows_rel]  = format_table_rows_line(lines[rows_idx], args.n_points)
    b2_block[b2_time_rel]  = format_data_row(t_grid,  lines[time_idx])
    b2_block[b2_power_rel] = format_data_row(P_beam2, lines[power_idx])

    # ── Update "Number of laser beams = N" → "= N+1"
    new_n_beams = n_beams + 1
    new_lines[n_beams_idx] = re.sub(
        r'(Number of laser beams\s*=\s*)\d+',
        f'\\g<1>{new_n_beams}', new_lines[n_beams_idx], count=1,
        flags=re.IGNORECASE)
    if not new_lines[n_beams_idx].endswith('\n'):
        new_lines[n_beams_idx] += '\n'

    # ── Splice the new beam block in just BEFORE [End Laser Source Data].
    # Inserting at the end (rather than directly after the active beam)
    # avoids renumbering any existing placeholder beams.
    output_lines = (new_lines[:laser_end_idx]
                    + b2_block
                    + new_lines[laser_end_idx:])

    if args.verbose:
        print(f"Active beam {active_id} block: lines {b1_start+1}–{b1_end} "
              f"({b1_end - b1_start} lines)")
        print(f"Pulse table:  rows-meta line {rows_idx+1}, "
              f"times line {time_idx+1}, powers line {power_idx+1}")
        print(f"Number-of-beams line: {n_beams_idx+1} "
              f"({lines[n_beams_idx].rstrip()!r} → "
              f"{new_lines[n_beams_idx].rstrip()!r})")
        print(f"New ramp+peak beam: ID = {new_beam_id}, "
              f"{len(b2_block)} lines, inserted at line "
              f"{laser_end_idx+1} of the output (just before "
              f"[End Laser Source Data])")
        if args.foot_cone is not None:
            print(f"Foot cone:    overridden to {args.foot_cone}°")
        if args.foot_spot is not None:
            print(f"Foot spot:    overridden to {args.foot_spot} cm")
        if args.foot_focus is not None:
            print(f"Foot focus:   overridden to {args.foot_focus} cm")
        print()

    if args.dry_run:
        print(f"--dry-run: would write {out_rhw}")
        print(f"  Total output lines: {len(output_lines)} "
              f"(input had {len(lines)})")
        return 0

    out_rhw.parent.mkdir(parents=True, exist_ok=True)
    with open(out_rhw, 'w') as fh:
        fh.writelines(output_lines)

    print(f"Wrote {out_rhw}")
    print()
    print("To run Helios on the split-pulse rhw:")
    print(f"  ~/Codes/Prism/Helios_11.0.0/Helios.app/Contents/MacOS/Helios \\")
    print(f"      -b -i {out_rhw} \\")
    print(f"      -d <results_root> -o {out_rhw.stem} -x")
    print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
