"""
fix_json_rhw_geometry.py
========================

Patch laser-beam geometry fields (``Spot size``, ``Half cone angle``,
optional ``Focus position``) in a Helios 11.1.0+ JSON workspace .rhw file.
Writes the result to a new output path (refuses to overwrite the input).

Designed for the FL transition work (June 2026): the 11.1.0 GUI saves
.rhw as JSON, and we frequently need to re-anchor an inherited workspace
to PDD_20_s016 geometry (cone=20 deg, spot=0.16 cm) before running.

Usage
-----
    python examples/fix_json_rhw_geometry.py <in.rhw> <out.rhw> \\
        --spot 0.16 --cone 20 [--focus 0.22] [--beam N] [--dry-run]

If ``--beam N`` is omitted, the patch applies to ALL beams in the file.
``--focus`` is optional — omit it to leave the focus position unchanged.

Numeric values are written with 4 decimal places (matches the precision
Helios uses in JSON output). Region/structural keys are left untouched.

This is the JSON sibling of edit_rhw_flux_limiter.py's text path; for
the legacy text format use the corresponding Helios 11.0.0 tools.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional


# JSON keys within "Laser beam element[N]"
_KEY_SPOT  = 'Spot size'
_KEY_CONE  = 'Half cone angle'
_KEY_FOCUS = 'Focus position'


def _format_value(v: float) -> str:
    """Render to 4 decimal places (matches Helios 11.1.0 JSON precision)."""
    return f"{v:.4f}"


def _patch_beam(beam: dict,
                spot_cm: Optional[float],
                cone_deg: Optional[float],
                focus_cm: Optional[float]) -> dict:
    """Return a dict describing the changes made to a single beam block.

    The beam dict is mutated in place. The returned dict is for reporting
    and contains keys 'spot', 'cone', 'focus' each mapping to (old, new)
    when changed, omitted when not changed.
    """
    changes = {}
    if spot_cm is not None and _KEY_SPOT in beam:
        old = float(beam[_KEY_SPOT])
        if old != spot_cm:
            beam[_KEY_SPOT] = spot_cm
            changes['spot'] = (old, spot_cm)
    if cone_deg is not None and _KEY_CONE in beam:
        old = float(beam[_KEY_CONE])
        if old != cone_deg:
            beam[_KEY_CONE] = cone_deg
            changes['cone'] = (old, cone_deg)
    if focus_cm is not None and _KEY_FOCUS in beam:
        old = float(beam[_KEY_FOCUS])
        if old != focus_cm:
            beam[_KEY_FOCUS] = focus_cm
            changes['focus'] = (old, focus_cm)
    return changes


def fix_json_rhw_geometry(in_path: Path,
                          out_path: Path,
                          spot_cm: Optional[float] = None,
                          cone_deg: Optional[float] = None,
                          focus_cm: Optional[float] = None,
                          beam_id: Optional[int] = None,
                          dry_run: bool = False) -> int:
    """
    Patch laser-beam geometry in a JSON .rhw and write the result.

    Parameters
    ----------
    in_path : Path
        Input .rhw file (must be JSON workspace format).
    out_path : Path
        Output .rhw file path (must differ from input).
    spot_cm, cone_deg, focus_cm : float | None
        New values; None means "leave this field untouched".
    beam_id : int | None
        If given, patch only ``Laser beam element[beam_id]``; otherwise
        patch every beam in the file.
    dry_run : bool
        If True, print intended changes but do not write the output file.

    Returns
    -------
    int
        Number of beams actually changed (at least one field updated).
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

    # Find all "Laser beam element[N]" blocks. Helios 11.1.0 nests them
    # under "Laser source data" but use a generic walk to be robust to
    # future restructurings.
    targets = []     # list of (beam_id, beam_dict)
    _collect_beam_blocks(data, targets)
    targets.sort(key=lambda x: x[0])

    if not targets:
        raise ValueError("No 'Laser beam element[N]' blocks found in JSON.")

    if beam_id is not None:
        targets = [(bid, b) for (bid, b) in targets if bid == beam_id]
        if not targets:
            raise ValueError(
                f"Beam id {beam_id} not found. "
                f"Available: {[bid for (bid, _) in _collect_all_beams(data)]}")

    print(f"Input:  {in_path}")
    print(f"Output: {out_path}")
    print(f"Beams in file: {[bid for (bid, _) in _collect_all_beams(data)]}")
    print(f"Patching:      {[bid for (bid, _) in targets]}")
    print(f"Updates:       "
          f"spot={spot_cm}, cone={cone_deg}, focus={focus_cm}")
    print()

    n_changed = 0
    for bid, beam in targets:
        changes = _patch_beam(beam, spot_cm, cone_deg, focus_cm)
        if changes:
            n_changed += 1
            print(f"  Beam {bid}:")
            for field, (old, new) in changes.items():
                print(f"    {field:>6s}  {old:>9.4f}  ->  {new:>9.4f}")
        else:
            print(f"  Beam {bid}: no change (values already match or "
                  f"fields not present)")

    if dry_run:
        print()
        print("[dry-run] Skipping write.")
        return n_changed

    out_path.write_text(json.dumps(data, indent=2))
    print()
    print(f"Wrote JSON to {out_path}")
    return n_changed


def _collect_beam_blocks(obj, out):
    """Walk JSON and collect ('Laser beam element[N]', dict) entries.

    Mutates `out` with (beam_id, beam_dict_reference).
    """
    if isinstance(obj, dict):
        for k, v in obj.items():
            if (isinstance(k, str) and k.startswith('Laser beam element[')
                    and isinstance(v, dict)):
                try:
                    bid = int(k[len('Laser beam element['):-1])
                except ValueError:
                    continue
                out.append((bid, v))
            if isinstance(v, (dict, list)):
                _collect_beam_blocks(v, out)
    elif isinstance(obj, list):
        for v in obj:
            if isinstance(v, (dict, list)):
                _collect_beam_blocks(v, out)


def _collect_all_beams(obj):
    """Convenience: full beam inventory for reporting."""
    found = []
    _collect_beam_blocks(obj, found)
    found.sort(key=lambda x: x[0])
    return found


def main():
    ap = argparse.ArgumentParser(
        description="Patch laser-beam geometry in a Helios 11.1.0+ JSON .rhw.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n"
               "  python examples/fix_json_rhw_geometry.py \\\n"
               "      Olson_PDD_20_fab06_foot25_s016_burn_beta.rhw \\\n"
               "      Olson_PDD_20_fab06_foot25_s016_burn_beta_v2.rhw \\\n"
               "      --spot 0.16 --cone 20 --dry-run",
    )
    ap.add_argument('in_rhw',  type=Path, help="Input JSON .rhw file path.")
    ap.add_argument('out_rhw', type=Path,
                    help="Output .rhw file path (must differ from input).")
    ap.add_argument('--spot',  type=float, default=None, metavar='CM',
                    help="New spot size in cm.")
    ap.add_argument('--cone',  type=float, default=None, metavar='DEG',
                    help="New half cone angle in degrees.")
    ap.add_argument('--focus', type=float, default=None, metavar='CM',
                    help="New focus position in cm (optional).")
    ap.add_argument('--beam',  type=int, default=None, metavar='N',
                    help="If given, patch only Laser beam element[N]; "
                         "otherwise patch every beam in the file.")
    ap.add_argument('--dry-run', action='store_true',
                    help="Show intended changes but do not write output.")
    args = ap.parse_args()

    if args.spot is None and args.cone is None and args.focus is None:
        ap.error("At least one of --spot, --cone, --focus is required.")

    try:
        n = fix_json_rhw_geometry(args.in_rhw, args.out_rhw,
                                   spot_cm=args.spot,
                                   cone_deg=args.cone,
                                   focus_cm=args.focus,
                                   beam_id=args.beam,
                                   dry_run=args.dry_run)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    return 0 if (n > 0 or args.dry_run) else 1


if __name__ == '__main__':
    sys.exit(main())
