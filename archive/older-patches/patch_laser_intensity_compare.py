#!/usr/bin/env python
"""
patch_laser_intensity_compare.py -- Task 2 Stage C.4

Add two laser-intensity keys to the implosion_rows list in
compare_with_published() so that _published.json can include reference
values for cross-checking against LILAC/HYDRA publications.

Keys added:
  'I_at_crit_peak_Wcm2'     -- peak intensity at critical surface [W/cm^2]
  'I_grid_outer_peak_Wcm2'  -- peak intensity at grid outer boundary [W/cm^2]

These require that sim_metrics contains the same keys. We pull them from
the data object inside compare_with_published's callers - but more
sustainably, we add them to whatever dict builds sim_metrics. The
simplest hook is to inject them at the top of compare_with_published
from data attributes if present.

Pragmatic approach: we assume sim_metrics.get(..., 0.0) works -- if the
caller has pre-populated the metrics dict with these keys (from
data.I_at_crit_peak / data.I_grid_outer_peak), they get compared. If not,
the 'pub_val <= 0 and sim_val <= 0' filter skips the row. No harm done.

Idempotent via sentinel checks.

Run from repo root:
    python patch_laser_intensity_compare.py
"""
from pathlib import Path

TARGET = Path(
    '~/Codes/helios_postprocessor/helios_postprocess/burn_averaged_metrics.py'
).expanduser()
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
# EDIT: add two entries to the implosion_rows list.
# Anchor on the HS radius line (last entry before the closing ]), which
# gives us a unique spot with stable whitespace.
# ====================================================================

OLD = """        ('HS radius at ignition (\u03bcm)',    sim_metrics.get('hs_radius_ignition_um', 0.0),
         'hs_radius_ignition_um',    '.1f'),
    ]"""

NEW = """        ('HS radius at ignition (\u03bcm)',    sim_metrics.get('hs_radius_ignition_um', 0.0),
         'hs_radius_ignition_um',    '.1f'),
        # --- Laser intensity diagnostics (Task 2 Stage C) ---
        ('I at r_crit peak (W/cm\u00b2)', sim_metrics.get('I_at_crit_peak_Wcm2', 0.0),
         'I_at_crit_peak_Wcm2',       '.2e'),
        ('I grid outer peak (W/cm\u00b2)', sim_metrics.get('I_grid_outer_peak_Wcm2', 0.0),
         'I_grid_outer_peak_Wcm2',    '.2e'),
    ]"""

src, _ = apply(src, OLD, NEW,
               sentinel='I_at_crit_peak_Wcm2',
               label='implosion_rows: add laser intensity keys')


if src != orig:
    TARGET.write_text(src)
    print(f'  [WROTE]   {TARGET.name}')


print('\n-----------------------------------------------------------------------')
print('Published JSON key names (for a run.published.json):')
print('')
print('  "I_at_crit_peak_Wcm2":     [1.5e+15, 0.0]   # [value, uncertainty]')
print('  "I_grid_outer_peak_Wcm2":  [1.0e+15, 0.0]')
print('')
print('To populate sim_metrics with these, add this in whatever code builds')
print('the metrics dict (usually run_analysis.py or icf_analysis.py):')
print('')
print('    metrics["I_at_crit_peak_Wcm2"]    = getattr(data, "I_at_crit_peak", 0.0) or 0.0')
print('    metrics["I_grid_outer_peak_Wcm2"] = getattr(data, "I_grid_outer_peak", 0.0) or 0.0')
print('')
print('(A follow-up patch will do that wiring automatically once we see')
print('where the metrics dict is constructed. For now the new rows will')
print('simply be skipped when sim_val == 0 and pub_val == 0.)')
print('-----------------------------------------------------------------------')
