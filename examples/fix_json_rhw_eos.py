"""
fix_json_rhw_eos.py
===================

Patch per-region EOS-file pointers (``EOS filepath`` and ``EOS lib
directory``) in a Helios 11.1.0+ JSON workspace .rhw file. Writes the
result to a new output path (refuses to overwrite the input).

First use case (June 28 2026): repair the DT-CH foam region's EOS in
the 11.1.0 PDD calibration runs, which inherited a pure-DT SESAME
pointer instead of the PROPACEOS DT-CH foam table used in 11.0.0
production.  See CLAUDE.md Session Update for June 27 2026 / "Open
questions for Prism" for the diagnosis.

Usage
-----
Edit one region by name:

    python examples/fix_json_rhw_eos.py <in.rhw> <out.rhw> \\
        --region "DT-CH foam" \\
        --filepath DT_CH_foam.prp \\
        --lib-dir /Users/.../PROPACEOS_6.1.0 \\
        [--data-type 0]                       \\
        [--dry-run]

Multiple regions in one invocation:

    python examples/fix_json_rhw_eos.py in.rhw out.rhw \\
        --set "DT-CH foam:DT_CH_foam.prp:/path/to/PROPACEOS_6.1.0" \\
        --set "CH Skin:CH.prp:/path/to/PROPACEOS_6.1.0"

JSON keys touched per region:
    'EOS filepath'        (str)   — the table filename
    'EOS lib directory'   (str)   — the parent directory
    'EOS data type'       (int)   — 0=PROPACEOS, 1=SESAME (optional;
                                    only updated if --data-type given)

Region match is by exact ``Region name`` string (case-sensitive,
internal whitespace preserved), matching how edit_rhw_flux_limiter.py
matches.  The script prints the regions it found in the file plus the
edits made.

This is the JSON sibling of any text-format EOS-edit tooling; for
Helios 11.0.0 text .rhw the EOS pointers are line-edited directly.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional


# JSON keys within "Spatial region element[N]"
_KEY_FILEPATH  = 'EOS filepath'
_KEY_LIB_DIR   = 'EOS lib directory'
_KEY_DATA_TYPE = 'EOS data type'


def _parse_set_arg(s: str) -> tuple:
    """Parse 'REGION:FILEPATH:LIB_DIR[:DATA_TYPE]' into a tuple.

    Returns (region_name, filepath, lib_dir, data_type_or_None).
    DATA_TYPE is optional; when given, must be 0 or 1.
    """
    parts = s.split(':')
    if len(parts) < 3 or len(parts) > 4:
        raise argparse.ArgumentTypeError(
            f"Expected 'REGION:FILEPATH:LIB_DIR[:DATA_TYPE]', got: {s!r}")
    region   = parts[0].strip()
    filepath = parts[1].strip()
    lib_dir  = parts[2].strip()
    data_type: Optional[int] = None
    if len(parts) == 4 and parts[3].strip():
        try:
            data_type = int(parts[3].strip())
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"DATA_TYPE must be 0 (PROPACEOS) or 1 (SESAME), "
                f"got: {parts[3]!r}")
        if data_type not in (0, 1):
            raise argparse.ArgumentTypeError(
                f"DATA_TYPE must be 0 or 1, got: {data_type}")
    return (region, filepath, lib_dir, data_type)


def _walk_json_for_regions(obj, out):
    """Collect ('Spatial region element[N]', dict) into out.

    Mutates `out` with (region_id, region_dict_reference) tuples.
    """
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


def fix_json_rhw_eos(in_path: Path,
                     out_path: Path,
                     region_updates: dict,
                     dry_run: bool = False) -> int:
    """
    Patch per-region EOS pointers in a JSON .rhw and write the result.

    Parameters
    ----------
    in_path : Path
        Input .rhw file (must be JSON workspace format).
    out_path : Path
        Output .rhw file (must differ from input).
    region_updates : dict[str, dict]
        Mapping of region name -> {'filepath': str, 'lib_dir': str,
        'data_type': int|None}. Region name is matched exactly against
        each region's ``Region name`` JSON value.
    dry_run : bool
        If True, print intended changes but do not write the output file.

    Returns
    -------
    int
        Number of regions actually changed (at least one field updated).
    """
    if in_path.resolve() == out_path.resolve():
        raise ValueError("out_path must differ from in_path "
                          "(refusing to overwrite the input).")
    if not in_path.exists():
        raise FileNotFoundError(f"Input .rhw not found: {in_path}")

    text = in_path.read_text()
    first_char = text.lstrip()[:1] if text.lstrip() else ''
    if first_char != '{':
        raise ValueError(
            f"Input is not JSON workspace format (starts with {first_char!r}). "
            "Use the text-format tools for Helios 11.0.0 .rhw files.")

    data = json.loads(text)
    regions = []
    _walk_json_for_regions(data, regions)
    regions.sort(key=lambda x: x[0])

    if not regions:
        raise ValueError("No 'Spatial region element[N]' blocks found in JSON.")

    region_names_in_file = [r.get('Region name', f'Region {rid}')
                             for rid, r in regions]
    target_names = set(region_updates.keys())
    changes = []
    matched_names = set()

    for rid, r in regions:
        name = r.get('Region name', f'Region {rid}')
        if name not in target_names:
            continue
        spec = region_updates[name]
        per_field = {}
        old_fp = str(r.get(_KEY_FILEPATH, ''))
        new_fp = spec['filepath']
        if old_fp != new_fp:
            r[_KEY_FILEPATH] = new_fp
            per_field['filepath'] = (old_fp, new_fp)
        old_lib = str(r.get(_KEY_LIB_DIR, ''))
        new_lib = spec['lib_dir']
        if old_lib != new_lib:
            r[_KEY_LIB_DIR] = new_lib
            per_field['lib_dir'] = (old_lib, new_lib)
        if spec['data_type'] is not None:
            old_dt = int(r.get(_KEY_DATA_TYPE, -1))
            new_dt = int(spec['data_type'])
            if old_dt != new_dt:
                r[_KEY_DATA_TYPE] = new_dt
                per_field['data_type'] = (old_dt, new_dt)
        changes.append({
            'region_id': rid,
            'region':    name,
            'fields':    per_field,
        })
        matched_names.add(name)

    # Reporting
    print(f"Input:  {in_path}")
    print(f"Output: {out_path}")
    print(f"Regions found in file (in order): {region_names_in_file}")
    print(f"Regions requested for update:     {sorted(target_names)}")
    print()
    if changes:
        print("Edits:")
        for c in changes:
            if not c['fields']:
                print(f"  region[{c['region_id']}] {c['region']!r:<14}  "
                      f"(no field changed — already matches)")
                continue
            print(f"  region[{c['region_id']}] {c['region']!r:<14}")
            for field, (old, new) in c['fields'].items():
                print(f"    {field:>10s}  {old!r:>40s}  ->  {new!r}")
    else:
        print("No edits made — no matching region in the file.")

    missing = target_names - matched_names
    if missing:
        print()
        print(f"WARNING: requested region(s) not found in file: {sorted(missing)}")
        print(f"         Available region names: {region_names_in_file}")

    if dry_run:
        print()
        print("[dry-run] Skipping write.")
        return len(changes)

    out_path.write_text(json.dumps(data, indent=2))
    print()
    print(f"Wrote JSON to {out_path}")
    return len(changes)


def main():
    ap = argparse.ArgumentParser(
        description="Patch per-region EOS pointers in a Helios 11.1.0+ JSON .rhw.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Single region (most common):\n"
            "  python examples/fix_json_rhw_eos.py in.rhw out.rhw \\\n"
            "      --region \"DT-CH foam\" \\\n"
            "      --filepath DT_CH_foam.prp \\\n"
            "      --lib-dir /Users/me/Codes/Prism/DataTables/PROPACEOS_6.1.0 \\\n"
            "      --data-type 0\n\n"
            "  # Multiple regions in one shot:\n"
            "  python examples/fix_json_rhw_eos.py in.rhw out.rhw \\\n"
            "      --set 'DT-CH foam:DT_CH_foam.prp:/path/to/PROPACEOS_6.1.0:0' \\\n"
            "      --set 'CH Skin:CH.prp:/path/to/PROPACEOS_6.1.0:0'\n"
        ),
    )
    ap.add_argument('in_rhw',  type=Path, help="Input JSON .rhw file path.")
    ap.add_argument('out_rhw', type=Path,
                    help="Output .rhw file path (must differ from input).")
    # Single-region form
    ap.add_argument('--region',    type=str, default=None, metavar='NAME',
                    help="Region name to update (single-region form).")
    ap.add_argument('--filepath',  type=str, default=None, metavar='FILE',
                    help="EOS table filename for --region.")
    ap.add_argument('--lib-dir',   type=str, default=None, metavar='DIR',
                    dest='lib_dir',
                    help="EOS table parent directory for --region.")
    ap.add_argument('--data-type', type=int, default=None, metavar='N',
                    dest='data_type', choices=[0, 1],
                    help="EOS data type for --region (0=PROPACEOS, 1=SESAME).")
    # Multi-region form
    ap.add_argument('--set', action='append', dest='sets',
                    type=_parse_set_arg, default=[],
                    metavar='"REGION:FILE:DIR[:TYPE]"',
                    help="Multi-region form. May be repeated.")
    ap.add_argument('--dry-run', action='store_true',
                    help="Show intended changes but do not write output.")
    args = ap.parse_args()

    region_updates = {}
    if args.region is not None:
        if not (args.filepath and args.lib_dir):
            ap.error("--region requires --filepath and --lib-dir.")
        region_updates[args.region] = {
            'filepath': args.filepath,
            'lib_dir':  args.lib_dir,
            'data_type': args.data_type,
        }
    for (region, filepath, lib_dir, data_type) in args.sets:
        region_updates[region] = {
            'filepath': filepath,
            'lib_dir':  lib_dir,
            'data_type': data_type,
        }
    if not region_updates:
        ap.error("Supply --region/--filepath/--lib-dir OR --set 'REGION:FILE:DIR[:TYPE]' "
                  "(may repeat).")

    try:
        n = fix_json_rhw_eos(args.in_rhw, args.out_rhw, region_updates,
                              dry_run=args.dry_run)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    return 0 if (n > 0 or args.dry_run) else 1


if __name__ == '__main__':
    sys.exit(main())
