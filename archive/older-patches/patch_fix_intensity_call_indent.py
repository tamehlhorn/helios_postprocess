#!/usr/bin/env python
"""
patch_fix_intensity_call_indent.py -- fix IndentationError introduced by
patch_laser_intensity_plotter.py (Stage C.2).

Problem: the plotter-patch Edit 2 inserted `self._plot_laser_intensity(pdf)`
with a hardcoded 8-space indent. The actual call site in create_full_report
uses 12-space indent (nested inside the method's try block). Result:

    line 101:             if self.data.laser_energy_deposited is not None:
    line 102:                 self._plot_laser_deposition(pdf)
    line 103:         self._plot_laser_intensity(pdf)     <-- wrong, 8 spaces
    line 104:             if self.data.laser_power_delivered is not None:

Python parser hits line 103 with dedent-to-8, then line 104 with indent-to-12
and throws IndentationError. The module then fails to import anywhere in the
package.

Fix: replace the 8-space line with a 12-space line to match the surrounding
block's indentation. No other changes needed -- the method definition at
line 1060 is fine.

Idempotent: sentinel checks for both the broken state and the fixed state.

Run from repo root:
    python patch_fix_intensity_call_indent.py
"""
from pathlib import Path

TARGET = Path(
    '~/Codes/helios_postprocessor/helios_postprocess/icf_plotting.py'
).expanduser()
print(f'Target: {TARGET}\n')

# The broken line (8-space indent) sits between two 16-space-indented lines.
# Use the triple-line context to pin the exact location.
BROKEN = (
    '            if self.data.laser_energy_deposited is not None:\n'
    '                self._plot_laser_deposition(pdf)\n'
    '        self._plot_laser_intensity(pdf)\n'
    '            if self.data.laser_power_delivered is not None:\n'
)

# Fixed: 12-space indent on the intensity call, matching the surrounding
# if-statements' indent level. Sibling, not child, of the deposition if-block.
FIXED = (
    '            if self.data.laser_energy_deposited is not None:\n'
    '                self._plot_laser_deposition(pdf)\n'
    '            self._plot_laser_intensity(pdf)\n'
    '            if self.data.laser_power_delivered is not None:\n'
)

# Sentinel for idempotency: the FIXED state has '            self._plot_laser_intensity'
# (12 spaces), while the BROKEN has '        self._plot_laser_intensity' (8 spaces).
# Check for the exact 12-space version as sentinel.
FIXED_SENTINEL = '            self._plot_laser_intensity(pdf)'
BROKEN_SENTINEL = '        self._plot_laser_intensity(pdf)\n'


src = TARGET.read_text()

if FIXED in src:
    print('  [ALREADY FIXED] 12-space indent already present; no change.')
elif BROKEN in src:
    new_src = src.replace(BROKEN, FIXED, 1)
    TARGET.write_text(new_src)
    print('  [FIXED]   Replaced 8-space indent with 12-space indent.')
    print(f'  [WROTE]   {TARGET.name}')
else:
    # Neither exact context matched -- maybe the surrounding code differs.
    # Report both versions to the user so they can pick the right fix by hand.
    print('  [UNEXPECTED] Neither broken nor fixed exact pattern matched.')
    print('                The file may differ from what the diagnosis assumed.')
    print('                Please paste the output of:')
    print('                    grep -n "_plot_laser_intensity\\|_plot_laser_deposition\\|'
          'laser_energy_deposited\\|laser_power_delivered" \\')
    print('                        helios_postprocess/icf_plotting.py | head -15')


print('\n-----------------------------------------------------------------------')
print('Verify after running:')
print('  python3 -c "import helios_postprocess.icf_plotting; print(\\"import OK\\")"')
print('')
print('Expected: "import OK". If still an IndentationError, paste the error and')
print('the output of:')
print('  python3 <<PY')
print('  lines = open(\\"helios_postprocess/icf_plotting.py\\").readlines()')
print('  for i, line in enumerate(lines[95:115], start=96):')
print('      print(f\\"{i:4d}: {line!r}\\")')
print('  PY')
print('-----------------------------------------------------------------------')
