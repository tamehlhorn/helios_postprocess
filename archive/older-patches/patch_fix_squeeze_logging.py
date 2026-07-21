#!/usr/bin/env python
"""
patch_fix_squeeze_logging.py -- Stage A cleanup

Fixes a bug introduced by patch_laser_intensity_data.py (Stage A):

  When 2c added the laser_power_on_target squeeze block, str.replace orphaned
  the `if verbose: logger.info(...)` block that belonged to laser_power_delivered.
  The result: laser_power_delivered's squeeze no longer logs, and the log line
  is now (mis)attached to the laser_power_on_target block with the wrong label.

Final structure after this patch:

    if data.laser_power_delivered is not None and data.laser_power_delivered.ndim > 1:
        data.laser_power_delivered = data.laser_power_delivered.squeeze()
        if verbose:
            logger.info(f"  \u2713 laser_power_delivered      squeezed \u2192 {data.laser_power_delivered.shape}")
    if data.laser_power_on_target is not None and data.laser_power_on_target.ndim > 1:
        # (n_times, n_beam) -> (n_times,) by taking beam 0
        data.laser_power_on_target = data.laser_power_on_target[:, 0]
        if verbose:
            logger.info(f"  \u2713 laser_power_on_target      squeezed \u2192 {data.laser_power_on_target.shape}")

Idempotent: verifies the broken state before applying; skips if already fixed.

Run from repo root:
    python patch_fix_squeeze_logging.py
"""
from pathlib import Path

TARGET = Path('~/Codes/helios_postprocessor/helios_postprocess/data_builder.py').expanduser()
print(f'Target: {TARGET}\n')

src = TARGET.read_text()

# The broken block as it currently exists after Stage A patch
BROKEN = (
    '    if data.laser_power_delivered is not None and data.laser_power_delivered.ndim > 1:\n'
    '        data.laser_power_delivered = data.laser_power_delivered.squeeze()\n'
    '    if data.laser_power_on_target is not None and data.laser_power_on_target.ndim > 1:\n'
    '        # (n_times, n_beam) -> (n_times,) by taking beam 0\n'
    '        data.laser_power_on_target = data.laser_power_on_target[:, 0]\n'
    '        if verbose:\n'
    '            logger.info(f"  \u2713 laser_power_delivered      squeezed \u2192 {data.laser_power_delivered.shape}")\n'
)

FIXED = (
    '    if data.laser_power_delivered is not None and data.laser_power_delivered.ndim > 1:\n'
    '        data.laser_power_delivered = data.laser_power_delivered.squeeze()\n'
    '        if verbose:\n'
    '            logger.info(f"  \u2713 laser_power_delivered      squeezed \u2192 {data.laser_power_delivered.shape}")\n'
    '    if data.laser_power_on_target is not None and data.laser_power_on_target.ndim > 1:\n'
    '        # (n_times, n_beam) -> (n_times,) by taking beam 0\n'
    '        data.laser_power_on_target = data.laser_power_on_target[:, 0]\n'
    '        if verbose:\n'
    '            logger.info(f"  \u2713 laser_power_on_target      squeezed \u2192 {data.laser_power_on_target.shape}")\n'
)

# Sentinel for idempotency: the fixed state has TWO "squeezed" log lines,
# and the on_target one has its OWN label.
FIXED_SENTINEL = 'laser_power_on_target      squeezed'

if FIXED_SENTINEL in src:
    print('  [ALREADY FIXED] on_target squeeze log line is present; no change.')
elif BROKEN in src:
    src2 = src.replace(BROKEN, FIXED, 1)
    TARGET.write_text(src2)
    print('  [FIXED]         Orphaned verbose block re-attached; added on_target log line.')
    print(f'  [WROTE]         {TARGET.name}')
else:
    print('  [UNEXPECTED]    Neither the broken pattern nor the fixed sentinel was found.')
    print('                  The file may be in a state I did not anticipate.')
    print('                  Please paste the output of:')
    print('                    grep -B 2 -A 10 "laser_power_on_target is not None" \\')
    print('                      ~/Codes/helios_postprocessor/helios_postprocess/data_builder.py')

print('\n-----------------------------------------------------------------------')
print('Verify after running:')
print('  grep -B 2 -A 10 "laser_power_delivered is not None" \\')
print('    ~/Codes/helios_postprocessor/helios_postprocess/data_builder.py')
print('')
print('Expected: TWO separate if-blocks, each with its OWN "if verbose:" logger line.')
print('-----------------------------------------------------------------------')
