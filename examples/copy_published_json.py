"""
copy_published_json.py -- Propagate a template published JSON to every Helios
run folder under a root directory.

Walks <root> recursively for <base>.exo files; for each one, copies the source
JSON to <base>_published.json next to it. Existing files are preserved unless
--force is passed; --dry-run reports the plan without writing anything.

Useful for batch-reprocessing legacy runs against the same LILAC anchor:
fill in shock-arrival and rhoR reference values in one canonical JSON, then
distribute to every sim folder before running run_analysis.py in a loop.

Usage:
    python3 ~/helios_postprocess/examples/copy_published_json.py \\
        --source ~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab04_foot25_s016_burn/Olson_PDD_20_fab04_foot25_s016_burn_published.json \\
        --root   ~/Sims/Xcimer/Olson_PDD

    # Preview first
    python3 ~/helios_postprocess/examples/copy_published_json.py \\
        --source <path> --root <root> --dry-run

    # Replace stale copies
    python3 ~/helios_postprocess/examples/copy_published_json.py \\
        --source <path> --root <root> --force
"""

import argparse
import json
import shutil
import sys
from pathlib import Path


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--source', type=Path, required=True,
                    help='Source <base>_published.json template to copy.')
    ap.add_argument('--root', type=Path, required=True,
                    help='Root directory to walk recursively for *.exo files.')
    ap.add_argument('--force', action='store_true',
                    help='Overwrite existing _published.json files at targets.')
    ap.add_argument('--dry-run', action='store_true',
                    help='Report the plan without writing anything.')
    args = ap.parse_args(argv)

    # Validate source is a real JSON file we can parse (fail fast).
    if not args.source.is_file():
        print(f"ERROR: source not found: {args.source}", file=sys.stderr)
        return 1
    try:
        with open(args.source) as f:
            json.load(f)
    except json.JSONDecodeError as e:
        print(f"ERROR: source is not valid JSON: {args.source}\n  {e}",
              file=sys.stderr)
        return 1

    if not args.root.is_dir():
        print(f"ERROR: root is not a directory: {args.root}", file=sys.stderr)
        return 1

    source_resolved = args.source.resolve()
    exos = sorted(args.root.rglob('*.exo'))
    if not exos:
        print(f"No *.exo files under {args.root}", file=sys.stderr)
        return 1

    plan = []   # (action, target_path)
    for exo in exos:
        base   = exo.stem                              # filename without .exo
        target = exo.parent / f"{base}_published.json"
        if target.resolve() == source_resolved:
            plan.append(('skip_same', target)); continue
        if target.exists() and not args.force:
            plan.append(('skip_exists', target)); continue
        plan.append(('copy', target))

    # Execute
    n_copy = n_exists = n_same = 0
    for action, target in plan:
        if action == 'copy':
            verb = "WOULD COPY" if args.dry_run else "COPY"
            print(f"{verb}: -> {target}")
            if not args.dry_run:
                shutil.copy2(args.source, target)
            n_copy += 1
        elif action == 'skip_exists':
            n_exists += 1
        elif action == 'skip_same':
            n_same += 1

    print()
    print(f"Source: {args.source}")
    print(f"Found .exo files:    {len(exos)}")
    print(f"  Copied:            {n_copy}{'  (dry-run)' if args.dry_run else ''}")
    print(f"  Skipped (exists):  {n_exists}  (use --force to overwrite)")
    print(f"  Skipped (source):  {n_same}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
