#!/usr/bin/env python
"""
patch_rename_I_outer.py -- Task 2 Stage C.1

Rename I_outer_* to I_grid_outer_* to clarify that the reference surface is
the simulation grid outer boundary (which may extend well into vacuum),
NOT the capsule outer surface.

Rationale: in Olson_PDD_26b_burn the grid extends to ~0.8 cm while the
capsule radius is ~0.2 cm, so I = P / (4 pi R_grid^2) yields 4e13 W/cm^2
(seems too low), while the physically meaningful I at critical surface
is ~1.5e15 W/cm^2. The rename prevents future confusion.

Edits (idempotent):
  helios_postprocess/laser_intensity.py:
    - Variable `I_outer` in analyze_laser_intensity() internal logic: unchanged
      (still correct internally as grid-outer-boundary incident intensity
      feeding Method 2 Beer-Lambert; this IS where the ray starts)
    - Output dict keys: 'I_outer' -> 'I_grid_outer',
                        'peak_I_outer' -> 'peak_I_grid_outer'
    - Docstring updated

  helios_postprocess/icf_analysis.py:
    - data attribute I_outer_peak     -> I_grid_outer_peak
    - data attribute I_outer_history  -> I_grid_outer_history
    - Log message "Peak I_outer" -> "Peak I at grid outer boundary"
    - References to result['peak_I_outer'] / result['I_outer']
      updated to new key names

Attributes UNCHANGED (the physically correct ones):
    I_at_crit_peak, I_at_crit_at_peak_power, I_peak_coronal,
    I_at_crit_history, ncr_intensity, t_peak_power_ns,
    r_crit_intensity_history

Run from repo root:
    python patch_rename_I_outer.py
"""
from pathlib import Path

REPO = Path('~/Codes/helios_postprocessor').expanduser()
LI = REPO / 'helios_postprocess' / 'laser_intensity.py'
ICF = REPO / 'helios_postprocess' / 'icf_analysis.py'


def apply(src, old, new, sentinel, label):
    if sentinel in src and old not in src:
        print(f'  [ALREADY] {label}')
        return src, False
    if old not in src:
        print(f'  [MISS]    {label}: OLD pattern not found')
        return src, False
    print(f'  [OK]      {label}')
    return src.replace(old, new, 1), True


# ====================================================================
# EDIT GROUP 1: laser_intensity.py  -- rename output dict keys + docstring
# ====================================================================
print('=== laser_intensity.py: rename output dict keys ===\n')
print(f'Target: {LI}')

src = LI.read_text()
orig = src

# 1a. Rename the docstring entry describing the scalar
OLD = """      I_outer                    (nt,)       Incident outer-boundary intensity"""
NEW = """      I_grid_outer               (nt,)       Incident intensity at grid outer boundary"""
src, _ = apply(src, OLD, NEW,
               sentinel='Incident intensity at grid outer boundary',
               label='docstring: I_outer line')

# 1b. Rename docstring entry for peak_I_outer scalar
OLD = """      peak_I_outer               float       max over t"""
NEW = """      peak_I_grid_outer          float       max over t (at grid outer boundary)"""
src, _ = apply(src, OLD, NEW,
               sentinel='max over t (at grid outer boundary)',
               label='docstring: peak_I_outer line')

# 1c. Rename the two dict keys in the return statement
OLD = """    return dict(
        I1=I1, I2=I2, I_outer=I_outer, alpha_zone=alpha_zone,
        r_crit=r_crit, ncr=ncr,
        I_at_crit_vs_t=I_at_crit_vs_t,
        I_peak_coronal_vs_t=I_peak_coronal_vs_t,
        peak_I_outer=peak_I_outer,
        peak_I_at_crit=peak_I_at_crit,
        peak_I_coronal=peak_I_coronal,
        I_at_crit_at_peak_power=I_at_crit_at_peak_power,
        t_peak_power_ns=t_peak_ns,
        wavelength_um=wavelength_um,
    )"""
NEW = """    return dict(
        I1=I1, I2=I2, I_grid_outer=I_outer, alpha_zone=alpha_zone,
        r_crit=r_crit, ncr=ncr,
        I_at_crit_vs_t=I_at_crit_vs_t,
        I_peak_coronal_vs_t=I_peak_coronal_vs_t,
        peak_I_grid_outer=peak_I_outer,
        peak_I_at_crit=peak_I_at_crit,
        peak_I_coronal=peak_I_coronal,
        I_at_crit_at_peak_power=I_at_crit_at_peak_power,
        t_peak_power_ns=t_peak_ns,
        wavelength_um=wavelength_um,
    )"""
src, _ = apply(src, OLD, NEW,
               sentinel='I_grid_outer=I_outer',
               label='return dict: I_outer -> I_grid_outer, peak_I_outer -> peak_I_grid_outer')

if src != orig:
    LI.write_text(src)
    print(f'  [WROTE]   {LI.name}')


# ====================================================================
# EDIT GROUP 2: icf_analysis.py  -- rename consumer-side attributes + logs
# ====================================================================
print('\n=== icf_analysis.py: rename consumer-side attributes + log labels ===\n')
print(f'Target: {ICF}')

src = ICF.read_text()
orig = src

# 2a. Attribute assignment: I_outer_peak -> I_grid_outer_peak
#     Use the exact spacing from the current code (aligned columns)
OLD = "        self.data.I_outer_peak            = result['peak_I_outer']"
NEW = "        self.data.I_grid_outer_peak       = result['peak_I_grid_outer']"
src, _ = apply(src, OLD, NEW,
               sentinel="self.data.I_grid_outer_peak",
               label="attribute: I_outer_peak -> I_grid_outer_peak")

# 2b. Attribute assignment: I_outer_history
OLD = "        self.data.I_outer_history          = result['I_outer']"
NEW = "        self.data.I_grid_outer_history     = result['I_grid_outer']"
src, _ = apply(src, OLD, NEW,
               sentinel="self.data.I_grid_outer_history",
               label="attribute: I_outer_history -> I_grid_outer_history")

# 2c. Log message line
OLD = '        logger.info(f"  Peak I_outer                    : "'
NEW = '        logger.info(f"  Peak I at grid outer boundary   : "'
src, _ = apply(src, OLD, NEW,
               sentinel="Peak I at grid outer boundary",
               label="log message: Peak I_outer")

# 2d. Log message value reference
OLD = "                    f\"{result['peak_I_outer']:.3e} W/cm^2\")"
NEW = "                    f\"{result['peak_I_grid_outer']:.3e} W/cm^2\")"
src, _ = apply(src, OLD, NEW,
               sentinel="result['peak_I_grid_outer']",
               label="log message: peak_I_outer value ref")

# 2e. None-out list in the graceful-skip branch
OLD = """            for k in ('I_outer_peak', 'I_at_crit_peak', 'I_at_crit_at_peak_power',
                      'I_peak_coronal', 't_peak_power_ns', 'ncr_intensity',
                      'r_crit_intensity_history', 'I_at_crit_history',
                      'I_outer_history'):"""
NEW = """            for k in ('I_grid_outer_peak', 'I_at_crit_peak', 'I_at_crit_at_peak_power',
                      'I_peak_coronal', 't_peak_power_ns', 'ncr_intensity',
                      'r_crit_intensity_history', 'I_at_crit_history',
                      'I_grid_outer_history'):"""
src, _ = apply(src, OLD, NEW,
               sentinel="'I_grid_outer_peak'",
               label="graceful-skip list: I_outer_* -> I_grid_outer_*")

if src != orig:
    ICF.write_text(src)
    print(f'  [WROTE]   {ICF.name}')


print('\n-----------------------------------------------------------------------')
print('Verify after running:')
print('  grep -n "I_outer\\|I_grid_outer" helios_postprocess/laser_intensity.py \\')
print('     helios_postprocess/icf_analysis.py')
print('')
print('Expected: the word "I_outer" should only appear as the LOCAL variable')
print('          name inside analyze_laser_intensity() (that is by design --')
print('          it is the incident ray at the grid outer boundary). All')
print('          public attributes and log labels should read I_grid_outer.')
print('-----------------------------------------------------------------------')
