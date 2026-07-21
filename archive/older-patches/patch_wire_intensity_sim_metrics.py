#!/usr/bin/env python
"""
patch_wire_intensity_sim_metrics.py -- Task 2 Stage C.5

Activate the dormant comparison rows `I at r_crit peak (W/cm^2)` and
`I grid outer peak (W/cm^2)` in compare_with_published() by populating
`sim_metrics` with the corresponding keys.

Two edits to helios_postprocess/burn_averaged_metrics.py:

  5a. extract_histories_from_run_data(data) -- inject two keys into the
      returned dict, pulled directly from data.I_at_crit_peak and
      data.I_grid_outer_peak (populated by ICFAnalyzer.analyze_laser_intensity).

  5b. calculate_burn_averaged_metrics(histories) -- pass the two keys
      through from histories into the returned sim_metrics dict.

Missing-value policy: getattr(data, key, None) or 0.0 ensures the compare
table's "both zero -> skip" filter continues to work for runs where
analyze_laser_intensity was skipped (e.g., old .exo without the required
variables).

Idempotent via sentinel checks.

Run from repo root:
    python patch_wire_intensity_sim_metrics.py
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
# 5a. extract_histories_from_run_data: inject two new keys at end of dict.
# Anchor includes the closing brace of the dict + blank line + next def
# so the OLD pattern uniquely identifies the tail.
# ====================================================================
print('=== 5a. extract_histories_from_run_data: inject intensity keys ===\n')

OLD = """        'CR_max':                getattr(data, 'comp_ratio', 0.0),
    }


def _check_required(data, attr_names):"""

NEW = """        'CR_max':                getattr(data, 'comp_ratio', 0.0),
        # Laser intensity metrics (from analyze_laser_intensity); 0.0 if skipped
        'I_at_crit_peak_Wcm2':    getattr(data, 'I_at_crit_peak', None) or 0.0,
        'I_grid_outer_peak_Wcm2': getattr(data, 'I_grid_outer_peak', None) or 0.0,
    }


def _check_required(data, attr_names):"""

src, _ = apply(src, OLD, NEW,
               sentinel="'I_at_crit_peak_Wcm2':    getattr(data",
               label="extract_histories: add I_at_crit_peak_Wcm2, I_grid_outer_peak_Wcm2")


# ====================================================================
# 5b. calculate_burn_averaged_metrics: pass keys through from histories.
# Anchor includes the last existing implosion metric + the # Burn-rate
# weighting profile comment to uniquely identify the insertion point.
# ====================================================================
print('\n=== 5b. calculate_burn_averaged_metrics: pass through intensity keys ===\n')

OLD = """        'imploded_DT_mass_mg':   histories.get('imploded_DT_mass_mg', 0.0),
        # Burn-rate weighting profile
        'burn_fraction':    burn_fraction,"""

NEW = """        'imploded_DT_mass_mg':   histories.get('imploded_DT_mass_mg', 0.0),
        # Laser intensity metrics (pass through from histories)
        'I_at_crit_peak_Wcm2':    histories.get('I_at_crit_peak_Wcm2',    0.0),
        'I_grid_outer_peak_Wcm2': histories.get('I_grid_outer_peak_Wcm2', 0.0),
        # Burn-rate weighting profile
        'burn_fraction':    burn_fraction,"""

src, _ = apply(src, OLD, NEW,
               sentinel="'I_at_crit_peak_Wcm2':    histories.get",
               label="calculate_burn_averaged_metrics: pass through intensity keys")


if src != orig:
    TARGET.write_text(src)
    print(f'\n  [WROTE]   {TARGET.name}')
else:
    print('\n  [NO CHANGE] (all edits already applied)')


print('\n-----------------------------------------------------------------------')
print('Next step:')
print('  1. Also update the published JSON for 26b_burn to include the two new keys,')
print('     plus fix adiabat. Suggested additions:')
print()
print('     "adiabat":                 [3.0, 0.0],')
print('     "I_at_crit_peak_Wcm2":     [1.5e+15, 0.0],   # LILAC reference estimate')
print('     "I_grid_outer_peak_Wcm2":  [1.0e+15, 0.0],   # update with real number')
print()
print('  2. Re-run on Mac Studio:')
print('     python3 examples/run_analysis.py \\\\')
print('       ~/Sims/Xcimer/Olson_PDD/Olson_PDD_26b_burn/Olson_PDD_26b_burn')
print()
print('  3. In the comparison table expect two new rows, e.g.:')
print('     I at r_crit peak (W/cm^2)       1.46e+15     1.50e+15     -2.5')
print('     I grid outer peak (W/cm^2)      4.03e+13     1.00e+15    -96.0')
print()
print('  (The grid-outer delta will be nonsensical because published')
print('   "grid outer" is LILAC-geometry-specific; expect to leave that')
print('   at 0.0 or zero the entry out of the JSON once confirmed.)')
print('-----------------------------------------------------------------------')
