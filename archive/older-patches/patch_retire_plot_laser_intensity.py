#!/usr/bin/env python
"""
patch_retire_plot_laser_intensity.py -- Task 2 final cleanup

Delete the standalone plot_laser_intensity.py from the repo root now that
its functionality is fully integrated into the pipeline via
helios_postprocess/laser_intensity.py.

Safety checks BEFORE deletion:
  1. Confirm helios_postprocess/laser_intensity.py exists (the successor).
  2. Confirm helios_postprocess/icf_analysis.py defines analyze_laser_intensity.
  3. Confirm plot_laser_intensity.py has been superseded (file hash not
     checked -- just that the package equivalents exist).

Deletion is via Path.unlink(), NOT shell. This is a single-file operation
against a known path with no globbing.

After running, do `git status` to see the deletion; `git add -u` stages
the removal; commit and push.

If the script has been modified locally in ways you want to preserve,
cancel this step (don't run) and archive it manually first.

Run from repo root:
    python patch_retire_plot_laser_intensity.py
"""
from pathlib import Path
import sys

REPO = Path('~/Codes/helios_postprocessor').expanduser()
TARGET = REPO / 'plot_laser_intensity.py'
SUCCESSOR = REPO / 'helios_postprocess' / 'laser_intensity.py'
ANALYZER = REPO / 'helios_postprocess' / 'icf_analysis.py'

print(f'Target (to delete):  {TARGET}')
print(f'Successor (package): {SUCCESSOR}')
print()

# ---- Pre-flight checks ----
if not TARGET.exists():
    print('  [ALREADY GONE] plot_laser_intensity.py is not in the repo root.')
    print('                  Nothing to do.')
    sys.exit(0)

if not SUCCESSOR.exists():
    print('  [ABORT]  helios_postprocess/laser_intensity.py does NOT exist.')
    print('           Refusing to delete the standalone script before its')
    print('           package successor is in place. Run Stage A first.')
    sys.exit(1)

if not ANALYZER.exists():
    print('  [ABORT]  helios_postprocess/icf_analysis.py not found.')
    sys.exit(1)

analyzer_src = ANALYZER.read_text()
if 'def analyze_laser_intensity' not in analyzer_src:
    print('  [ABORT]  icf_analysis.py does not define analyze_laser_intensity.')
    print('           Run Stage B first.')
    sys.exit(1)

print('  [CHECK]  Successor module exists.')
print('  [CHECK]  ICFAnalyzer.analyze_laser_intensity method present.')
print()


# ---- Report file size before deletion, for the record ----
size = TARGET.stat().st_size
print(f'  Deleting: {TARGET.name} ({size} bytes)')
TARGET.unlink()
print('  [DELETED]')


print('\n-----------------------------------------------------------------------')
print('Next steps:')
print('  cd ~/Codes/helios_postprocessor')
print('  git status      # should show "deleted: plot_laser_intensity.py"')
print('  git add -u      # stage the deletion')
print('  git diff --cached --stat')
print('  git commit -m "Retire standalone plot_laser_intensity.py"')
print('               "  (logic now in helios_postprocess/laser_intensity.py"')
print('               "   and integrated PDF pipeline)"')
print('  git push')
print('-----------------------------------------------------------------------')
