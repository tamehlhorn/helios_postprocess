#!/usr/bin/env python
"""
patch_laser_intensity_data.py -- Task 2 Stage A (data side)

Edits to helios_postprocess/data_builder.py only:
  2a. ICFRunData  -- add laser_attenuation_coeff, laser_power_on_target fields
  2b. _VARIABLE_MAP -- load LaserPwrOnTargetForBeam and laserAttinuationCoeff
  2c. build_run_data post-load -- squeeze laser_power_on_target beam axis

Does NOT touch icf_analysis.py, icf_plotting.py, or icf_output.py.
Those follow once loading is verified clean.

Idempotent: sentinel checks skip already-applied edits.

Run from repo root:
    python patch_laser_intensity_data.py
"""
from pathlib import Path

TARGET = Path('~/Codes/helios_postprocessor/helios_postprocess/data_builder.py').expanduser()
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


# ---- 2a. ICFRunData: new array fields --------------------------------------
print('=== 2a. ICFRunData: laser_attenuation_coeff, laser_power_on_target ===')
OLD = """        self.electron_density: Optional[np.ndarray] = None        # (n_times, n_zones) cm\u207b\u00b3"""
NEW = """        self.electron_density: Optional[np.ndarray] = None        # (n_times, n_zones) cm\u207b\u00b3
        # Laser intensity reconstruction inputs (raw from EXODUS; cleaned in laser_intensity.py)
        self.laser_attenuation_coeff: Optional[np.ndarray] = None  # (n_times, n_beam, n_bnd) 1/cm
        self.laser_power_on_target: Optional[np.ndarray] = None    # (n_times,) W  -- LaserPwrOnTargetForBeam, beam 0"""
src, _ = apply(src, OLD, NEW,
               sentinel='self.laser_attenuation_coeff: Optional',
               label='ICFRunData laser_attenuation / on_target fields')


# ---- 2b. _VARIABLE_MAP: two new entries ------------------------------------
print('\n=== 2b. _VARIABLE_MAP: LaserPwrOnTargetForBeam, laserAttinuationCoeff ===')
OLD = """    ("electron_density",      ["electron_density", "elec_density", "n_e"],        False),"""
NEW = """    ("electron_density",      ["electron_density", "elec_density", "n_e"],        False),
    ("laser_power_on_target", ["LaserPwrOnTargetForBeam",
                               "laser_power_on_target"],                           False),
    ("laser_attenuation_coeff", ["laserAttinuationCoeff",
                                 "laser_attenuation_coeff"],                       False),"""
src, _ = apply(src, OLD, NEW,
               sentinel='"LaserPwrOnTargetForBeam"',
               label='_VARIABLE_MAP entries')


# ---- 2c. build_run_data: squeeze beam axis on laser_power_on_target --------
print('\n=== 2c. build_run_data: squeeze laser_power_on_target beam axis ===')
OLD = """    if data.laser_power_delivered is not None and data.laser_power_delivered.ndim > 1:
        data.laser_power_delivered = data.laser_power_delivered.squeeze()"""
NEW = """    if data.laser_power_delivered is not None and data.laser_power_delivered.ndim > 1:
        data.laser_power_delivered = data.laser_power_delivered.squeeze()
    if data.laser_power_on_target is not None and data.laser_power_on_target.ndim > 1:
        # (n_times, n_beam) -> (n_times,) by taking beam 0
        data.laser_power_on_target = data.laser_power_on_target[:, 0]"""
src, _ = apply(src, OLD, NEW,
               sentinel='data.laser_power_on_target = data.laser_power_on_target[:, 0]',
               label='squeeze laser_power_on_target')


if src != orig:
    TARGET.write_text(src)
    print(f'\n  [WROTE]   {TARGET.name}')
else:
    print('\n  [NO CHANGE] (all edits already applied)')


print('\n-----------------------------------------------------------------------')
print('Next steps:')
print('  1. Copy laser_intensity.py to helios_postprocess/ in the repo.')
print('  2. Commit both files and push.')
print('  3. Pull on Mac Studio and re-run on Olson_PDD_26b_burn.')
print('  4. In the run log, confirm these two new lines appear:')
print('       laser_attenuation_coeff       <- laserAttinuationCoeff ...')
print('       laser_power_on_target         <- LaserPwrOnTargetForBeam ...')
print('  5. The PDF/summary/CSV will be unchanged (that is Stages B/C/D).')
print('-----------------------------------------------------------------------')
