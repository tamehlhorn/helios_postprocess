#!/usr/bin/env python
"""
patch_laser_intensity_analysis.py -- Task 2 Stage B (analysis wiring)

Edits two files:
  helios_postprocess/icf_analysis.py -- add ICFAnalyzer.analyze_laser_intensity()
  examples/run_analysis.py           -- add analyzer.analyze_laser_intensity() call

The new method:
  - Calls helios_postprocess.laser_intensity.analyze_laser_intensity(self.data)
  - Populates scalars on self.data:
      I_outer_peak, I_at_crit_peak, I_at_crit_at_peak_power, I_peak_coronal,
      t_peak_power_ns, ncr_intensity
  - Populates time-history arrays on self.data:
      r_crit_intensity_history, I_at_crit_history, I_outer_history
  - Caches 2D arrays for the plotter in a private dict on self.data:
      self.data._laser_intensity_arrays = {I1, I2, alpha_zone, zcen, r_crit, ...}
    Stage C plotter reads and drops this; not persisted to CSV/summary.
  - Logs one cross-check line: r_crit at peak power, to compare against
    analyze_drive_phase's "Critical radius (formula)" message.

Graceful no-op: if required inputs are missing (.exo files that pre-date
Stage A, or non-laser drive), logs a warning and returns.

Idempotent via sentinel checks.

Run from repo root:
    python patch_laser_intensity_analysis.py
"""
from pathlib import Path

REPO = Path('~/Codes/helios_postprocessor').expanduser()
ICF_ANALYSIS = REPO / 'helios_postprocess' / 'icf_analysis.py'
RUNNER = REPO / 'examples' / 'run_analysis.py'


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
# EDIT 1: icf_analysis.py -- add import + analyze_laser_intensity method
# ====================================================================
print('=== Edit 1: icf_analysis.py ===\n')
print(f'Target: {ICF_ANALYSIS}')

src = ICF_ANALYSIS.read_text()
orig = src

# 1a. Add the import. Anchor on `def analyze_drive_phase` which is at line 58;
# we inject the import immediately before the first `class ICFAnalyzer` line.
# More robust: inject after any existing `from .` import near the top.
OLD_IMPORT_ANCHOR = '    def analyze_drive_phase(self):'
NEW_METHOD = '''    def analyze_laser_intensity(self):
        """
        Reconstruct laser intensity I(r, t) and summary diagnostics.

        Called after analyze_drive_phase so that the drive-phase logger
        output (critical density, absorbed energy) precedes intensity
        metrics in the run log.

        Stores on self.data:
          Scalars:
            I_outer_peak, I_at_crit_peak, I_at_crit_at_peak_power,
            I_peak_coronal, t_peak_power_ns, ncr_intensity
          Histories (nt,):
            r_crit_intensity_history, I_at_crit_history, I_outer_history
          2D arrays (cached for plotter, not persisted):
            self.data._laser_intensity_arrays dict

        Graceful no-op with warning if required inputs absent.
        """
        from .laser_intensity import analyze_laser_intensity as _ali
        import numpy as _np

        logger.info("Analyzing laser intensity...")

        wavelength_um = getattr(self.data, 'laser_wavelength_um', 0.351)
        if not wavelength_um or wavelength_um <= 0:
            wavelength_um = 0.351
        result = _ali(self.data, wavelength_um=wavelength_um)

        if result is None:
            logger.warning("  Skipping: required inputs missing "
                           "(laser_power_source, laser_attenuation_coeff, "
                           "laser_power_on_target, or zone_boundaries).")
            # Populate with None so downstream code can check cleanly
            for k in ('I_outer_peak', 'I_at_crit_peak', 'I_at_crit_at_peak_power',
                      'I_peak_coronal', 't_peak_power_ns', 'ncr_intensity',
                      'r_crit_intensity_history', 'I_at_crit_history',
                      'I_outer_history'):
                setattr(self.data, k, None)
            self.data._laser_intensity_arrays = None
            return

        # Scalars
        self.data.I_outer_peak            = result['peak_I_outer']
        self.data.I_at_crit_peak          = result['peak_I_at_crit']
        self.data.I_at_crit_at_peak_power = result['I_at_crit_at_peak_power']
        self.data.I_peak_coronal          = result['peak_I_coronal']
        self.data.t_peak_power_ns         = result['t_peak_power_ns']
        self.data.ncr_intensity           = result['ncr']

        # Histories
        self.data.r_crit_intensity_history = result['r_crit']
        self.data.I_at_crit_history        = result['I_at_crit_vs_t']
        self.data.I_outer_history          = result['I_outer']

        # Cached 2D arrays for plotter (not persisted)
        self.data._laser_intensity_arrays = dict(
            I1=result['I1'],
            I2=result['I2'],
            alpha_zone=result['alpha_zone'],
            r_crit=result['r_crit'],
            I_peak_coronal_vs_t=result['I_peak_coronal_vs_t'],
        )

        # Report
        logger.info(f"  Wavelength: {result['wavelength_um']:.3f} um  "
                    f"(n_crit = {result['ncr']:.3e} cm^-3)")
        t_pk = result['t_peak_power_ns']
        if _np.isfinite(t_pk):
            # Find r_crit at peak power for cross-check against drive-phase value
            try:
                t_idx = int(_np.argmin(_np.abs(self.data.time - t_pk)))
                r_crit_at_peak = result['r_crit'][t_idx]
                if _np.isfinite(r_crit_at_peak):
                    logger.info(f"  r_crit at peak power (t={t_pk:.2f} ns): "
                                f"{r_crit_at_peak:.4f} cm  "
                                f"(cross-check vs drive-phase formula method)")
            except Exception:
                pass
        logger.info(f"  Peak I_outer                    : "
                    f"{result['peak_I_outer']:.3e} W/cm^2")
        logger.info(f"  Peak I at critical surface      : "
                    f"{result['peak_I_at_crit']:.3e} W/cm^2")
        logger.info(f"  I at r_crit at peak laser power : "
                    f"{result['I_at_crit_at_peak_power']:.3e} W/cm^2")

    def analyze_drive_phase(self):'''

src, _ = apply(src, OLD_IMPORT_ANCHOR, NEW_METHOD,
               sentinel='def analyze_laser_intensity(self):',
               label='icf_analysis: insert analyze_laser_intensity method')

if src != orig:
    ICF_ANALYSIS.write_text(src)
    print(f'  [WROTE]   {ICF_ANALYSIS.name}')


# ====================================================================
# EDIT 2: run_analysis.py -- add call after analyze_drive_phase
# ====================================================================
print('\n=== Edit 2: examples/run_analysis.py ===\n')
print(f'Target: {RUNNER}')

src = RUNNER.read_text()
orig = src

OLD = '    analyzer.analyze_drive_phase()\n    analyzer.analyze_stagnation_phase()\n'
NEW = ('    analyzer.analyze_drive_phase()\n'
       '    analyzer.analyze_laser_intensity()\n'
       '    analyzer.analyze_stagnation_phase()\n')

src, _ = apply(src, OLD, NEW,
               sentinel='analyzer.analyze_laser_intensity()',
               label='run_analysis: insert analyze_laser_intensity call')

if src != orig:
    RUNNER.write_text(src)
    print(f'  [WROTE]   {RUNNER.name}')


print('\n-----------------------------------------------------------------------')
print('Next steps:')
print('  1. git diff helios_postprocess/icf_analysis.py examples/run_analysis.py')
print('  2. Read the diff. Expect ~80 new lines in icf_analysis.py, 1 line in runner.')
print('  3. git add, commit, push.')
print('  4. Pull on Mac Studio and re-run on Olson_PDD_26b_burn.')
print('  5. In the log, look for a new section after "Drive phase" output:')
print('       Analyzing laser intensity...')
print('         Wavelength: 0.350 um  (n_crit = 9.10e+21 cm^-3)')
print('         r_crit at peak power (t=XX.XX ns): 0.XXXX cm')
print('         Peak I_outer                    : X.XXe+14 W/cm^2')
print('         Peak I at critical surface      : X.XXe+14 W/cm^2')
print('         I at r_crit at peak laser power : X.XXe+14 W/cm^2')
print('  6. Confirm r_crit at peak power agrees with the earlier drive-phase line:')
print('       Critical radius (formula): 0.0787 \u00b1 0.0386 cm')
print('-----------------------------------------------------------------------')
