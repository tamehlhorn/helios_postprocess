#!/usr/bin/env python
"""
patch_laser_intensity_output.py -- Task 2 Stage C.3

Add LASER INTENSITY block to the summary text, inserted immediately after
LASER CONFIGURATION (beam 1) block (which currently ends at line ~203) and
before EOS MODELS (line 205).

The new block reports:
  - Wavelength and derived n_crit
  - Peak I at critical surface           [W/cm^2]
  - I at r_crit at peak laser power      [W/cm^2]
  - Peak coronal I (Method 2, max over r,t)
  - Peak I at grid outer boundary        (clearly labeled as grid, not capsule)
  - r_crit at peak laser power           [um]

Graceful no-op: if any scalar is None (analyze_laser_intensity skipped),
the section is omitted entirely rather than printed with "-" placeholders.

Idempotent via sentinel checks.

Run from repo root:
    python patch_laser_intensity_output.py
"""
from pathlib import Path

TARGET = Path('~/Codes/helios_postprocessor/helios_postprocess/icf_output.py').expanduser()
print(f'Target: {TARGET}\n')


def apply(src, old, new, sentinel, label):
    if sentinel in src and old not in src:
        print(f'  [ALREADY] {label}')
        return src, False
    if old not in src:
        print(f'  [MISS]    {label}: OLD pattern not found')
        return src, False
    print(f'  [OK]      {label}')
    return src.replace(old, new, 1), True


src = TARGET.read_text()
orig = src


# ====================================================================
# Insert LASER INTENSITY block between LASER CONFIGURATION and EOS MODELS.
# The existing structure ends the LASER CONFIGURATION block with:
#     _a('')
#
#         # ---- EOS models ----
#         if getattr(d, 'eos_models', None):
# Anchor on the EOS models opener to place the new block just before it.
# ====================================================================

ANCHOR = """        # ---- EOS models ----
        if getattr(d, 'eos_models', None):"""

NEW_BLOCK = '''        # ---- Laser intensity (from analyze_laser_intensity) ----
        _I_crit_peak   = getattr(d, 'I_at_crit_peak', None)
        _I_crit_atpk   = getattr(d, 'I_at_crit_at_peak_power', None)
        _I_peak_coro   = getattr(d, 'I_peak_coronal', None)
        _I_grid_outer  = getattr(d, 'I_grid_outer_peak', None)
        _t_peak_pwr    = getattr(d, 't_peak_power_ns', None)
        _ncr           = getattr(d, 'ncr_intensity', None)
        _r_crit_hist   = getattr(d, 'r_crit_intensity_history', None)

        if _I_crit_peak is not None:
            _a('LASER INTENSITY')
            _a('-' * width)
            if _ncr is not None and _ncr > 0:
                _a(f"  {'n_crit':<36s} {_ncr:>10.3e} cm^-3")
            _a(f"  {'Peak I at critical surface':<36s} "
               f"{_I_crit_peak:>10.3e} W/cm^2")
            if _I_crit_atpk is not None:
                _a(f"  {'I at r_crit at peak laser power':<36s} "
                   f"{_I_crit_atpk:>10.3e} W/cm^2")
            if _I_peak_coro is not None:
                _a(f"  {'Peak coronal I (max over r, t)':<36s} "
                   f"{_I_peak_coro:>10.3e} W/cm^2")
            if _I_grid_outer is not None:
                _a(f"  {'Peak I at grid outer boundary':<36s} "
                   f"{_I_grid_outer:>10.3e} W/cm^2")
            # r_crit at peak power, for cross-check vs drive-phase value
            if _t_peak_pwr is not None and _r_crit_hist is not None:
                try:
                    import numpy as _np
                    if _np.isfinite(_t_peak_pwr) and hasattr(d, 'time'):
                        _tidx = int(_np.argmin(_np.abs(d.time - _t_peak_pwr)))
                        _rc = float(_r_crit_hist[_tidx])
                        if _np.isfinite(_rc):
                            _a(f"  {'r_crit at peak laser power':<36s} "
                               f"{_rc * 1e4:>10.1f} um "
                               f"(t={_t_peak_pwr:.2f} ns)")
                except Exception:
                    pass
            _a('')

        # ---- EOS models ----
        if getattr(d, 'eos_models', None):'''

src, _ = apply(src, ANCHOR, NEW_BLOCK,
               sentinel="'LASER INTENSITY'",
               label='insert LASER INTENSITY block before EOS MODELS')


if src != orig:
    TARGET.write_text(src)
    print(f'  [WROTE]   {TARGET.name}')


print('\n-----------------------------------------------------------------------')
print('Verify after running:')
print('  grep -B 1 -A 10 "LASER INTENSITY" \\')
print('    ~/Sims/Xcimer/Olson_PDD/Olson_PDD_26b_burn/Olson_PDD_26b_burn_summary.txt')
print('')
print('Expected 6 new lines in the summary, between LASER CONFIGURATION and')
print('EOS MODELS sections, showing:')
print('  Peak I at critical surface, I at r_crit at peak laser power,')
print('  Peak coronal I, Peak I at grid outer boundary, r_crit at peak power.')
print('-----------------------------------------------------------------------')
