#!/usr/bin/env python
"""
patch_claude_md_task2_docs.py -- CLAUDE.md documentation update for Task 2.

Applies three targeted edits to CLAUDE.md, all inside the existing
Laser Intensity Diagnostics section:

  Edit 1: Insert "Pipeline integration (April 2026 -- Task 2)" subsection
          AFTER lesson #5 ("Peak coronal intensity grows through the pulse...")
          and BEFORE the next top-level "## " or end of file.

  Edit 2: Remove the "plot_laser_intensity.py" row from the
          "Diagnostic Scripts (repo root)" table and add a retirement note.

  Edit 3: Replace the plot_laser_intensity.py usage example block with
          run_analysis.py example + note that standalone is retired.

Safety posture (stricter than the code patches):

  - Refuses to write if ANY anchor appears zero times or more than once
  - Refuses to write if any edit's sentinel is already present
    (interpreted as "already applied; re-running would duplicate")
  - Prints before/after line counts; no silent side effects
  - All file I/O limited to CLAUDE.md in the current directory (or via
    --repo arg); no shell, no network

Run from the repo root:
    python patch_claude_md_task2_docs.py
    python patch_claude_md_task2_docs.py --dry-run     # show what would change
    python patch_claude_md_task2_docs.py --repo ~/path # explicit repo root
"""
import sys
import argparse
from pathlib import Path


# ===================================================================
# Anchors (short distinctive strings expected to appear EXACTLY ONCE)
# ===================================================================

# Edit 1 anchor: last lesson before we insert the new subsection.
# The lesson block ends with a blank line after "threshold analysis."
# We anchor on the final-sentence line.
ANCHOR_E1 = "   threshold analysis."

# Edit 2 anchor: exact table row for plot_laser_intensity
ANCHOR_E2 = "| `plot_laser_intensity.py` | I(r,t) reconstruction using Methods 1 & 2; 3-page PDF with cross-check |"

# Edit 3 anchor: the specific 3-line example block that uses the retired script.
ANCHOR_E3 = (
    "python3 ~/helios_postprocessor/plot_laser_intensity.py \\\n"
    "  ~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6/VI_6.exo\n"
    "python3 ~/helios_postprocessor/plot_laser_intensity.py \\\n"
    "  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_2021_01a/Olson_PDD_2021_01a.exo \\\n"
    "  --wavelength_um 0.351 --ntimes 6"
)

# ===================================================================
# Replacement text
# ===================================================================

# Edit 1: insert this block AFTER the anchor line. Keeps the existing
# "threshold analysis." paragraph intact, adds a blank line, then the
# new subsection.
EDIT1_OLD = ANCHOR_E1
EDIT1_NEW = """   threshold analysis.

### Pipeline integration (April 2026 -- Task 2)

All logic above lives in the automated pipeline. The standalone
`plot_laser_intensity.py` has been retired; reconstruction runs as part
of `run_analysis.py`.

| Layer | Contribution |
|---|---|
| `data_builder.py` | Loads `laser_attenuation_coeff` + `laser_power_on_target` via `_VARIABLE_MAP`; squeezes beam axis with labeled log line |
| `helios_postprocess/laser_intensity.py` | Module with `clean_attenuation`, `compute_method1`, `compute_method2`, `find_critical_radius_*`, `analyze_laser_intensity` entry point |
| `ICFAnalyzer.analyze_laser_intensity()` | Called after `analyze_drive_phase`; populates scalars and histories on `data`; caches 2D arrays in `_laser_intensity_arrays` for the plotter |
| `ICFPlotter._plot_laser_intensity()` | 3 PDF pages (log10 I(r,t) heatmap with r_crit overlay; P_laser + I histories; Method 1 vs 2 cross-check at peak power) |
| `ICFOutputGenerator` | `LASER INTENSITY` section in `<name>_summary.txt`, between `LASER CONFIGURATION` and `EOS MODELS` |
| `burn_averaged_metrics.py` | `I_at_crit_peak_Wcm2` and `I_grid_outer_peak_Wcm2` keys plumbed through `histories` -> `sim_metrics` -> `compare_with_published` |

**Renamed attribute:** `I_outer_*` -> `I_grid_outer_*` throughout the
pipeline code. The simulation grid extends well into vacuum (e.g., 0.8
cm while the capsule outer radius is 0.2 cm), so `I = P / (4*pi*R_grid^2)`
is much lower than intensity at the capsule / critical surface.
`I at critical surface` remains the primary physical quantity;
grid-outer is reported for audit purposes only. Inside
`laser_intensity.py` the local variable is still named `I_outer` since
it is the incident ray at the grid outer boundary (Method 2
Beer-Lambert's starting point); all attributes, log messages, and
summary output use `I_grid_outer`.

**M1 filter change:** the old standalone script's threshold
(`alpha > 0.01 * alpha_max_t`) cut at ~100 cm^-1 at the critical
surface, masking the entire absorbing layer. Replaced with an absolute
floor `ALPHA_MIN_M1 = 1e-2 cm^-1` in
`helios_postprocess/laser_intensity.py`. Coronal-noise exclusion is
preserved without losing the absorbing layer itself.

**Published-JSON keys:**

    "I_at_crit_peak_Wcm2":     [value, unc]   # primary comparison target
    "I_grid_outer_peak_Wcm2":  [0.0, 0.0]     # leave at 0; geometry-dependent,
                                              # not cross-code comparable

**Graceful degradation:** if a simulation lacks `LaserPwrOnTargetForBeam`
or `laserAttinuationCoeff` (pre-upgrade `.exo`), `analyze_laser_intensity`
logs a warning and `None`-populates attributes. Plotter emits a
placeholder page; summary block is omitted; compare rows are filtered by
the zero-skip rule. No crashes.

**Cross-check:** on `Olson_PDD_26b_burn`, the pipeline's peak-power
r_crit (0.099 cm / 992 um) sits inside the drive-phase formula-method
1-sigma band (0.079 +/- 0.039 cm), consistent with r_crit expanding
outward at peak drive relative to the all-timesteps mean."""

# Sentinel: a distinctive string from the new subsection that must NOT
# already exist in the file (otherwise the edit was already applied).
EDIT1_SENTINEL = "### Pipeline integration (April 2026 -- Task 2)"


# Edit 2: remove the plot_laser_intensity.py row, add retirement caveat.
# We replace the single row with an empty string + append the caveat
# immediately after the table closes. To do this cleanly we anchor on
# the row + following table delimiter structure.
EDIT2_OLD = """| `plot_laser_intensity.py` | I(r,t) reconstruction using Methods 1 & 2; 3-page PDF with cross-check |
| `helios_postprocessor_guide.docx` | User guide for collaborators |"""

EDIT2_NEW = """| `helios_postprocessor_guide.docx` | User guide for collaborators |

*(`plot_laser_intensity.py` retired April 2026 -- fully integrated into pipeline; see `helios_postprocess/laser_intensity.py` and the `analyze_laser_intensity` phase in `ICFAnalyzer`.)*"""

EDIT2_SENTINEL = "`plot_laser_intensity.py` retired April 2026"


# Edit 3: replace the old usage example block.
EDIT3_OLD = ANCHOR_E3
EDIT3_NEW = """# Standalone diagnostic scripts (still supported):
python3 ~/helios_postprocessor/plot_adiabat_shock.py \\
  ~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/VI_6/VI_6.exo
python3 ~/helios_postprocessor/plot_laser_deposition.py \\
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_2021_01a/Olson_PDD_2021_01a.exo

# Laser intensity now part of the main pipeline -- no separate script needed:
python3 ~/helios_postprocessor/examples/run_analysis.py \\
  ~/Sims/Xcimer/Olson_PDD/Olson_PDD_26b_burn/Olson_PDD_26b_burn
# produces <base>_report.pdf with 3 intensity pages,
#          <base>_summary.txt with LASER INTENSITY section,
#          <base>_comparison.pdf with intensity rows (if _published.json
#            includes \"I_at_crit_peak_Wcm2\")"""

EDIT3_SENTINEL = "# Standalone diagnostic scripts (still supported):"


# ===================================================================
# Driver
# ===================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--repo', default='.',
                    help='Repo root (default: current directory)')
    ap.add_argument('--dry-run', action='store_true',
                    help='Report what would change; do not write')
    args = ap.parse_args()

    target = (Path(args.repo).expanduser() / 'CLAUDE.md').resolve()
    print(f'Target: {target}')
    if not target.exists():
        print(f'  [ABORT] CLAUDE.md not found at {target}')
        return 1
    src = target.read_text()
    orig_lines = src.count('\n')
    print(f'  Before: {orig_lines} lines')
    print()

    # ---- Pre-flight: every anchor must appear EXACTLY once
    print('=== Pre-flight: anchor uniqueness ===')
    checks = [
        ('Edit 1 anchor', EDIT1_OLD, 1),
        ('Edit 2 anchor', EDIT2_OLD, 1),
        ('Edit 3 anchor', EDIT3_OLD, 1),
    ]
    fail = False
    for label, anchor, expected in checks:
        n = src.count(anchor)
        status = 'OK' if n == expected else 'FAIL'
        print(f'  [{status}]  {label}: found {n} time(s), expected {expected}')
        if n != expected:
            fail = True
    if fail:
        print()
        print('  [ABORT] At least one anchor has wrong multiplicity.')
        print('          File may have been partially edited already. Inspect')
        print('          with git diff, or run: git checkout -- CLAUDE.md')
        return 2
    print()

    # ---- Pre-flight: no sentinel already present
    print('=== Pre-flight: idempotency check (sentinels) ===')
    sentinels = [
        ('Edit 1 sentinel', EDIT1_SENTINEL),
        ('Edit 2 sentinel', EDIT2_SENTINEL),
        ('Edit 3 sentinel', EDIT3_SENTINEL),
    ]
    already = False
    for label, sentinel in sentinels:
        if sentinel in src:
            print(f'  [ALREADY] {label}: found in file (edit likely applied)')
            already = True
        else:
            print(f'  [FRESH]   {label}: absent; edit will apply')
    if already:
        print()
        print('  [ABORT] At least one edit already applied.')
        print('          Refusing to apply partially; revert with')
        print('          `git checkout -- CLAUDE.md` and rerun, or apply')
        print('          the missing edits manually.')
        return 3
    print()

    # ---- Apply all three edits
    print('=== Applying edits ===')
    new_src = src
    for label, old, new in [
        ('Edit 1: insert pipeline integration subsection', EDIT1_OLD, EDIT1_NEW),
        ('Edit 2: retire plot_laser_intensity in scripts table', EDIT2_OLD, EDIT2_NEW),
        ('Edit 3: replace retired usage examples', EDIT3_OLD, EDIT3_NEW),
    ]:
        before = new_src.count(old)
        if before != 1:
            print(f'  [INTERNAL] {label}: anchor count changed mid-flight '
                  f'({before}). Aborting.')
            return 4
        new_src = new_src.replace(old, new, 1)
        print(f'  [APPLIED] {label}')

    # Sanity: all sentinels must now be present exactly once
    post_fail = False
    for label, sentinel in sentinels:
        n = new_src.count(sentinel)
        if n != 1:
            print(f'  [POST] {label}: {n} occurrences after apply (expected 1)')
            post_fail = True
    if post_fail:
        print('  [ABORT] Post-apply sentinel count wrong. Not writing.')
        return 5

    new_lines = new_src.count('\n')
    delta = new_lines - orig_lines
    print()
    print(f'  After:  {new_lines} lines  (delta {delta:+d})')

    if args.dry_run:
        print()
        print('  [DRY-RUN] No file written.')
        return 0

    target.write_text(new_src)
    print(f'  [WROTE]  {target.name}')
    print()
    print('Next steps:')
    print('  git diff --stat CLAUDE.md    # expect ~+60, -5')
    print('  git diff CLAUDE.md | less    # review the three sections changed')
    print('  git add CLAUDE.md')
    print('  git commit -m "CLAUDE.md: Task 2 pipeline integration; retire plot_laser_intensity.py"')
    print('  git push')
    return 0


if __name__ == '__main__':
    sys.exit(main())
