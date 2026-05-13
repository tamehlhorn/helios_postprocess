"""
edit_rhw_flux_limiter.py
========================

Edit per-region flux-limiter (`Flux limiter mult.`) values in a Helios .rhw
file.  Reads the input file, locates each `Parameters for Region = <name>`
block, and for regions named via --set rewrites the flux-limiter multiplier
line.  Writes the result to a new output path (refuses to overwrite the input).

Designed as the safer alternative to sed/awk for FL parameter sweeps:
matches by region name (not by numeric value), so it can't accidentally
change unrelated 0.06's elsewhere in the file.

Usage
-----
    python examples/edit_rhw_flux_limiter.py <in.rhw> <out.rhw> \\
        --set "DT ice"=0.06 \\
        --set "Wetted foam"=0.06 \\
        --set "CD shell"=0.06

    # Add --dry-run to preview without writing.

Region names must match the .rhw exactly as they appear after
`Parameters for Region =` (case-sensitive, internal whitespace preserved).
The script prints the regions it found in the file plus the edits it made,
so you can verify the spelling before committing the change.

Numeric format
--------------
New FL values are written with 3 decimal places (e.g. 0.060) to match the
typical .rhw convention.  If your .rhw uses a different precision and you
need bit-exact preservation elsewhere in the file, adjust `_format_value`.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


# Matches:  '  Parameters for Region = DT ice'   (allowing leading whitespace)
_RE_REGION = re.compile(r'^\s*Parameters for Region\s*=\s*(.+?)\s*$',
                         re.IGNORECASE)

# Matches:  '  Flux limiter mult.     = 0.040    (anything)'
# Captures: (prefix-up-to-and-including '= '), (numeric value), (trailing chars)
_RE_FL_MULT = re.compile(
    r'^(\s*Flux limiter mult\.?\s*=\s*)([-+0-9.eE]+)(.*)$',
    re.IGNORECASE,
)


def _format_value(v: float) -> str:
    """Render FL value to 3 decimal places (matches Helios .rhw convention)."""
    return f"{v:.3f}"


def _parse_set_arg(s: str) -> tuple:
    """Parse 'region name=value' into (region_name, value).

    Uses rsplit so region names containing '=' (unlikely but possible)
    don't break parsing — the rightmost '=' separates name from value.
    """
    if '=' not in s:
        raise argparse.ArgumentTypeError(
            f"Expected '<region name>=<value>', got: {s!r}")
    name, val = s.rsplit('=', 1)
    name = name.strip()
    try:
        value = float(val.strip())
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Could not parse FL value as float: {val!r}")
    return (name, value)


def edit_rhw_flux_limiter(in_path: Path,
                           out_path: Path,
                           region_updates: dict,
                           dry_run: bool = False) -> int:
    """
    Rewrite the flux-limiter mult value for each region in region_updates.

    Parameters
    ----------
    in_path : Path
        Input .rhw file.
    out_path : Path
        Output .rhw file (must differ from in_path).
    region_updates : dict[str, float]
        Mapping of region name (exact match) to new flux-limiter value.
    dry_run : bool
        If True, print intended changes but do not write the output file.

    Returns
    -------
    int
        Number of flux-limiter lines actually edited.
    """
    if in_path.resolve() == out_path.resolve():
        raise ValueError("out_path must differ from in_path "
                          "(refusing to overwrite the input).")
    if not in_path.exists():
        raise FileNotFoundError(f"Input .rhw not found: {in_path}")

    lines = in_path.read_text().splitlines(keepends=True)
    target_names = set(region_updates.keys())

    current_region = None
    seen_regions = []           # in order of appearance
    matched_names = set()
    changes = []                # list of dicts for reporting
    new_lines = []

    for i, line in enumerate(lines, start=1):
        m_region = _RE_REGION.match(line)
        if m_region:
            current_region = m_region.group(1)
            seen_regions.append(current_region)
            new_lines.append(line)
            continue

        m_flmult = _RE_FL_MULT.match(line)
        if (m_flmult is not None
                and current_region is not None
                and current_region in target_names):
            old_value = float(m_flmult.group(2))
            new_value = region_updates[current_region]
            tail = m_flmult.group(3)
            # Reassemble the line: prefix + new value + trailing chars.
            new_line = f"{m_flmult.group(1)}{_format_value(new_value)}{tail}"
            # Preserve original line terminator if any.
            if line.endswith('\n') and not new_line.endswith('\n'):
                new_line += '\n'

            changes.append({
                'line_no': i,
                'region': current_region,
                'old': old_value,
                'new': new_value,
            })
            matched_names.add(current_region)
            new_lines.append(new_line)
        else:
            new_lines.append(line)

    # ── Reporting ──
    print(f"Input:  {in_path}")
    print(f"Output: {out_path}")
    print(f"Regions found in file (in order): {seen_regions}")
    print(f"Region updates requested: "
          f"{ {k: _format_value(v) for k, v in region_updates.items()} }")
    print()

    if changes:
        print("Edits:")
        for c in changes:
            print(f"  line {c['line_no']:>5}   region {c['region']!r:<18}  "
                  f"FL {c['old']:.4f} -> {c['new']:.4f}")
    else:
        print("No edits made — no matching region/flux-limiter pair found.")

    missing = target_names - matched_names
    if missing:
        print()
        print(f"WARNING: requested region(s) not edited (name mismatch?): "
              f"{sorted(missing)}")
        print(f"         Available region names in file: {seen_regions}")

    if dry_run:
        print()
        print("[dry-run] Skipping write.")
        return len(changes)

    out_path.write_text(''.join(new_lines))
    print()
    print(f"Wrote {len(new_lines)} lines to {out_path}")
    return len(changes)


def main():
    ap = argparse.ArgumentParser(
        description="Edit per-region flux-limiter values in a Helios .rhw file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n"
               "  python edit_rhw_flux_limiter.py \\\n"
               "      HDD26_DTI40_1ns130_FL04_lrm4_nb.rhw \\\n"
               "      HDD26_DTI40_1ns130_FL06_lrm4_nb.rhw \\\n"
               '      --set "DT ice"=0.06 \\\n'
               '      --set "Wetted foam"=0.06 \\\n'
               '      --set "CD shell"=0.06 \\\n'
               "      --dry-run",
    )
    ap.add_argument('in_rhw',  type=Path, help="Input .rhw file path.")
    ap.add_argument('out_rhw', type=Path,
                    help="Output .rhw file path (must differ from input).")
    ap.add_argument('--set', action='append', dest='sets',
                    type=_parse_set_arg, default=[],
                    metavar='"REGION"=VALUE',
                    help="Region name and new FL value.  May be repeated.")
    ap.add_argument('--dry-run', action='store_true',
                    help="Show intended changes but do not write the output file.")
    args = ap.parse_args()

    if not args.sets:
        ap.error("At least one --set '<REGION>'=<VALUE> is required.")

    # If multiple --set for the same region were passed, the last one wins.
    region_updates = dict(args.sets)

    try:
        n = edit_rhw_flux_limiter(args.in_rhw, args.out_rhw,
                                    region_updates,
                                    dry_run=args.dry_run)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    return 0 if (n > 0 or args.dry_run) else 1


if __name__ == '__main__':
    sys.exit(main())
