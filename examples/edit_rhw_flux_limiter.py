"""
edit_rhw_flux_limiter.py
========================

Edit per-region flux-limiter (`Flux limiter mult.`) values in a Helios .rhw
file.  Auto-detects file format and handles both:

  * Legacy text format (Helios 11.0.0 and earlier): locates each
    ``Parameters for Region = <name>`` block and rewrites the
    ``Flux limiter mult.`` line for named regions.
  * JSON workspace format (Helios 11.1.0+): walks the JSON tree to
    find each ``Spatial region element[N]`` block and updates its
    ``Flux limiter mult`` field.

Writes the result to a new output path (refuses to overwrite the input).

Designed as the safer alternative to sed/awk for FL parameter sweeps:
matches by region name (not by numeric value), so it can't accidentally
change unrelated 0.06's elsewhere in the file.

Usage
-----
Per-region edits (both formats):
    python examples/edit_rhw_flux_limiter.py <in.rhw> <out.rhw> \\
        --set "DT ice"=0.06 \\
        --set "Wetted foam"=0.06 \\
        --set "CD shell"=0.06

Apply the same FL to every region (convenience for FL sweeps):
    python examples/edit_rhw_flux_limiter.py <in.rhw> <out.rhw> \\
        --all 0.04

    # Add --dry-run to preview without writing.

Region names must match the .rhw exactly as they appear after
``Parameters for Region =`` (text format) or as the ``Region name``
value in the JSON. The script prints the regions it found in the file
plus the edits it made, so you can verify the spelling before
committing the change.

Numeric format
--------------
New FL values are written with 3 decimal places (e.g. 0.060) in text
format (Helios convention) and as native JSON numbers in JSON format.
If your .rhw uses a different text precision and you need bit-exact
preservation elsewhere in the file, adjust ``_format_value``.
"""

from __future__ import annotations

import argparse
import json
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


def _discover_text_regions_and_apply_all(lines, value: float) -> dict:
    """Scan a text-format .rhw for region names and return a uniform-value
    region_updates dict mapping each region to ``value``."""
    region_names = []
    for line in lines:
        m = _RE_REGION.match(line)
        if m:
            region_names.append(m.group(1))
    return {name: value for name in region_names}


def _walk_json_for_regions(obj, out):
    """Collect ('Spatial region element[N]', dict) entries into out as
    (region_id, region_dict_reference) tuples. Walks nested structures."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if (isinstance(k, str) and k.startswith('Spatial region element[')
                    and isinstance(v, dict)):
                try:
                    rid = int(k[len('Spatial region element['):-1])
                except ValueError:
                    continue
                out.append((rid, v))
            if isinstance(v, (dict, list)):
                _walk_json_for_regions(v, out)
    elif isinstance(obj, list):
        for v in obj:
            if isinstance(v, (dict, list)):
                _walk_json_for_regions(v, out)


def _edit_json_format(in_path: Path,
                      out_path: Path,
                      text: str,
                      region_updates: dict = None,
                      apply_all: float = None,
                      dry_run: bool = False) -> int:
    """JSON-format flux-limiter editor. Walks the workspace and updates
    ``Flux limiter mult`` per region per the caller's spec.

    Same semantics as the text path: region_updates takes a name->value
    map, apply_all assigns a uniform value to every region. Mutually
    exclusive.
    """
    if region_updates is None and apply_all is None:
        raise ValueError("Must supply either region_updates or apply_all.")
    if apply_all is not None and region_updates is not None:
        raise ValueError("region_updates and apply_all are mutually exclusive.")

    data = json.loads(text)
    regions = []
    _walk_json_for_regions(data, regions)
    regions.sort(key=lambda x: x[0])

    if not regions:
        raise ValueError("No 'Spatial region element[N]' blocks found in JSON.")

    region_names_in_file = [r.get('Region name', f'Region {rid}')
                             for rid, r in regions]

    if apply_all is not None:
        # Uniform value to every region (by name)
        region_updates = {name: apply_all for name in region_names_in_file}

    target_names = set(region_updates.keys())
    changes = []
    matched_names = set()
    for rid, r in regions:
        name = r.get('Region name', f'Region {rid}')
        if name not in target_names:
            continue
        new_val = region_updates[name]
        old_val = float(r.get('Flux limiter mult', 0.0))
        if old_val != new_val:
            r['Flux limiter mult'] = new_val
        changes.append({
            'region_id': rid,
            'region':    name,
            'old':       old_val,
            'new':       new_val,
        })
        matched_names.add(name)

    # Reporting
    print(f"Input:  {in_path}  (JSON format)")
    print(f"Output: {out_path}")
    print(f"Regions found in file (in order): {region_names_in_file}")
    print(f"Region updates requested: "
          f"{ {k: _format_value(v) for k, v in region_updates.items()} }")
    print()
    if changes:
        print("Edits:")
        for c in changes:
            print(f"  region[{c['region_id']}] {c['region']!r:<18}  "
                  f"FL {c['old']:.4f} -> {c['new']:.4f}")
    else:
        print("No edits made — no matching region/flux-limiter pair found.")

    missing = target_names - matched_names
    if missing:
        print()
        print(f"WARNING: requested region(s) not edited (name mismatch?): "
              f"{sorted(missing)}")
        print(f"         Available region names in file: {region_names_in_file}")

    if dry_run:
        print()
        print("[dry-run] Skipping write.")
        return len(changes)

    out_path.write_text(json.dumps(data, indent=2))
    print()
    print(f"Wrote JSON to {out_path}")
    return len(changes)


def edit_rhw_flux_limiter(in_path: Path,
                           out_path: Path,
                           region_updates: dict = None,
                           apply_all: float = None,
                           dry_run: bool = False) -> int:
    """
    Rewrite the flux-limiter mult value for each region.

    Format dispatch: auto-detected from the first non-whitespace
    character ('{' -> JSON, otherwise legacy text).

    Parameters
    ----------
    in_path : Path
        Input .rhw file.
    out_path : Path
        Output .rhw file (must differ from in_path).
    region_updates : dict[str, float] | None
        Mapping of region name (exact match) to new flux-limiter value.
        Mutually exclusive with apply_all.
    apply_all : float | None
        If given, apply this value to every region (FL-sweep convenience).
        Mutually exclusive with region_updates.
    dry_run : bool
        If True, print intended changes but do not write the output file.

    Returns
    -------
    int
        Number of flux-limiter entries actually edited.
    """
    if in_path.resolve() == out_path.resolve():
        raise ValueError("out_path must differ from in_path "
                          "(refusing to overwrite the input).")
    if not in_path.exists():
        raise FileNotFoundError(f"Input .rhw not found: {in_path}")

    text = in_path.read_text()
    first_char = text.lstrip()[:1] if text.lstrip() else ''
    if first_char == '{':
        return _edit_json_format(in_path, out_path, text,
                                  region_updates=region_updates,
                                  apply_all=apply_all,
                                  dry_run=dry_run)

    # Legacy text format: apply_all expands into region_updates after
    # the first pass discovers the region list.
    if region_updates is None and apply_all is None:
        raise ValueError("Must supply either region_updates or apply_all.")
    if apply_all is not None and region_updates is not None:
        raise ValueError("region_updates and apply_all are mutually exclusive.")

    lines = text.splitlines(keepends=True)
    if apply_all is not None:
        region_updates = _discover_text_regions_and_apply_all(lines, apply_all)
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
                    help="Region name and new FL value.  May be repeated. "
                         "Mutually exclusive with --all.")
    ap.add_argument('--all', type=float, default=None, metavar='VALUE',
                    help="Apply this FL value to every region. "
                         "Mutually exclusive with --set. "
                         "Convenience for FL sweeps.")
    ap.add_argument('--dry-run', action='store_true',
                    help="Show intended changes but do not write the output file.")
    args = ap.parse_args()

    if not args.sets and args.all is None:
        ap.error("Supply either --set '<REGION>'=<VALUE> (may repeat) or "
                  "--all <VALUE>.")
    if args.sets and args.all is not None:
        ap.error("--set and --all are mutually exclusive.")

    # If multiple --set for the same region were passed, the last one wins.
    region_updates = dict(args.sets) if args.sets else None

    try:
        n = edit_rhw_flux_limiter(args.in_rhw, args.out_rhw,
                                    region_updates=region_updates,
                                    apply_all=args.all,
                                    dry_run=args.dry_run)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    return 0 if (n > 0 or args.dry_run) else 1


if __name__ == '__main__':
    sys.exit(main())
