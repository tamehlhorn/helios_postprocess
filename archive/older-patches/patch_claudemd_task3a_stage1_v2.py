#!/usr/bin/env python
"""
patch_claudemd_task3a_stage1.py -- Document Task 3.A Stage 1 in CLAUDE.md.

Two edits to CLAUDE.md:

  1. Update stale "todo" bullet for I-profile integration to DONE (Task 2
     already delivered this as part of _plot_laser_intensity).

  2. Insert new top-level section "## Implosion & Ablation Diagnostics
     (Task 3.A, April 2026)" between the end of Laser Intensity Diagnostics
     (line ~804) and the "## Dependencies" section header (line ~806).

The new section documents Stage 1 deliverables, reference values on
PDD_26b_burn, the P_abl semantic note (176.69 Gbar vs Lindl's 6.77 Gbar
prediction), and queued Stages 2-10.

Idempotent. Refuses if sentinels already present.

Run from repo root:
    python patch_claudemd_task3a_stage1.py
    python patch_claudemd_task3a_stage1.py --dry-run
"""
import sys
import argparse
from pathlib import Path


# ===================================================================
# Edit 1: mark stale "integrate I-profile" bullet as DONE
# ===================================================================
# Anchor: the exact bullet text. It's marked with "- Integrate..." at
# column 0 (top-level bullet), immediately before "## Dependencies".

EDIT1_OLD = "- Integrate I-profile output into ICFAnalyzer as an optional page in the main PDF."
EDIT1_NEW = "- Integrate I-profile output into ICFAnalyzer as an optional page in the main PDF. (DONE April 2026 -- Task 2 Stage C.2, see `_plot_laser_intensity` in icf_plotting.py.)"
EDIT1_SENTINEL = "(DONE April 2026 -- Task 2 Stage C.2"


# ===================================================================
# Edit 2: insert new top-level section before "## Dependencies"
# ===================================================================
# Anchor: "## Dependencies" header plus the line before it. The line
# before is the bullet we just updated in Edit 1. To keep the anchor
# stable even after Edit 1 modifies that line, anchor ONLY on the
# "## Dependencies" header plus its first two content lines.

EDIT2_OLD = """## Dependencies

**Required**: numpy, scipy, matplotlib, netCDF4
**Optional**: scikit-learn (RANSAC shock fitting in icf_plotting.py -- guarded with `_HAS_SKLEARN`)"""

EDIT2_NEW = """## Implosion & Ablation Diagnostics (Task 3.A, April 2026)

Timing milestones and ablation-physics scalars. Collaborator request
in this priority order:

1. First shock breakout from fuel to gas
2. Second shock breakout (same interface)
3. Subsequent shock breakouts (3+ shock designs)
4. Shock flash at r=0
5. Time of peak velocity
6. Ablation pressure in ablator (CH / foam phase)
7. Ablation pressure in fuel (when ablator fully consumed)

### Stage 1 (shipped April 2026)

Three additions wired through `ICFAnalyzer` and summary:

| Attribute on ICFRunData | Source | Exposure |
|---|---|---|
| `t_peak_velocity_ns` | `time[peak_velocity_index]` during `analyze_implosion_phase` | Augmented "Peak implosion velocity" log line; new "Peak velocity time" row in TIMING summary |
| `ablation_pressure_Gbar` | `(ion + rad pressure)[t, ablation_front_indices[t]] * 1e-8` for each timestep | History on `data`, available to plotters |
| `P_abl_peak_Gbar` | `max(ablation_pressure_Gbar)` over valid timesteps | New "Peak ablation pressure" log line at end of `_track_ablation_front` |

**PDD_26b_burn reference values (Stage 1 validation):**

    t_peak_velocity_ns          12.000 ns
    P_abl_peak_Gbar             176.69 Gbar

### Important semantic note on `P_abl_peak_Gbar`

Stage 1 measures **total pressure at the ablation-front zone**, where
the ablation front is defined as the steepest negative density gradient
outside the hot spot. This tracks the outer boundary of the DENSE shell,
which sits inside the compressed region and is closer to shell pressure
than to the classical Lindl momentum-balance ablation pressure.

Lindl's scaling `P_abl = 57 (I/1e14)^(2/3) (lambda/um)^(-2/3)` Mbar
at PDD_26b's peak critical-surface intensity (I_crit = 1.46e15 W/cm^2,
lambda = 0.35 um) predicts 6.77 Gbar -- about 26x lower than the
176.69 Gbar reported here. The difference is physical, not a bug:
shell pressure > ablation drive pressure by the time the shell is
converging. Stage 4 (material-split) will clarify the interpretation
by separating ablator-phase and fuel-phase averages.

### Stage 2-10 queued

- **Stage 2:** Debug `_compute_shock_breakout` -- currently returns
  zeros for PDD_26b due to `self.data.time * 1e9` double-multiplication
  (`data.time` is already in ns per `data_builder.py:344`).
- **Stage 3:** Extend shock breakout to N-shock detection, producing
  `shock_breakout_times_ns` list + per-shock scalars.
- **Stage 4:** Material-split `P_abl` via
  `material_index[ablation_front_indices[t]]`; adds
  `P_abl_ablator_peak_Gbar`, `P_abl_fuel_peak_Gbar`,
  `t_fuel_ablation_start_ns`, `fuel_mass_ablated_mg`.
- **Stage 5:** Shock flash at r=0 -- peak pressure at innermost zone
  before stagnation -> `t_shock_flash_ns`.
- **Stage 6:** TIMING MILESTONES and ABLATION sections in summary text;
  comparison JSON keys for the new scalars.
- **Stage 7:** New ABLATION plotter page with material-phase shading.
- **Stage 8:** Wire `adiabat_history.py` module into `ICFAnalyzer`.
- **Stage 9:** Integrate `plot_adiabat_shock.py` into pipeline
  (overlaps with Stages 2-5).
- **Stage 10:** Integrate `plot_laser_deposition.py` into pipeline.

### Note on shock breakout bug

The `analyze_first_shock` call in `_compute_shock_breakout` converts
`self.data.time * 1e9`, which would be correct if `data.time` were in
seconds. But `data_builder.py` line 344 already multiplies by 1e9
during load, so `data.time` is in nanoseconds. The double-multiplication
makes every timestep ~1e10, which fails the foot-pulse mask
`(t_ns >= 4.0) & (t_ns <= 6.0)` trivially. Result: NaN breakout values,
which the error handler converts to 0.0. Stage 2 will fix this with a
one-line change.

## Dependencies

**Required**: numpy, scipy, matplotlib, netCDF4
**Optional**: scikit-learn (RANSAC shock fitting in icf_plotting.py -- guarded with `_HAS_SKLEARN`)"""

EDIT2_SENTINEL = "## Implosion & Ablation Diagnostics (Task 3.A, April 2026)"


# ===================================================================
# Driver
# ===================================================================

def apply(src, label, old, new, sentinel):
    if sentinel in src:
        print(f'  [ALREADY] {label}')
        return src, False
    n = src.count(old)
    if n == 0:
        print(f'  [MISS]    {label}: OLD pattern not found')
        return src, False
    if n > 1:
        print(f'  [AMBIG]   {label}: OLD pattern found {n} times; refusing')
        return src, False
    print(f'  [OK]      {label}')
    return src.replace(old, new, 1), True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default='.',
                    help='Repo root (default: current directory)')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    target = (Path(args.repo).expanduser() / 'CLAUDE.md').resolve()
    print(f'Target: {target}')
    if not target.exists():
        print(f'  [ABORT] CLAUDE.md not found at {target}')
        return 1
    src = target.read_text()
    orig_lines = src.count('\n')
    print(f'  Before: {orig_lines} lines\n')

    new_src = src

    print('=== Applying edits ===')
    new_src, ch1 = apply(new_src,
                         'Edit 1: mark stale I-profile bullet DONE',
                         EDIT1_OLD, EDIT1_NEW, EDIT1_SENTINEL)
    new_src, ch2 = apply(new_src,
                         'Edit 2: insert Implosion & Ablation Diagnostics section',
                         EDIT2_OLD, EDIT2_NEW, EDIT2_SENTINEL)

    if new_src == src:
        print('\n  (no changes)')
        return 0

    new_lines = new_src.count('\n')
    print(f'\n  After:  {new_lines} lines  (delta {new_lines - orig_lines:+d})')

    if args.dry_run:
        print('\n  [DRY-RUN] No file written.')
        return 0

    target.write_text(new_src)
    print(f'  [WROTE]   {target.name}')
    print()
    print('Next steps:')
    print('  git diff --stat CLAUDE.md')
    print('  git diff CLAUDE.md | less')
    print('  git add CLAUDE.md')
    print('  git commit -m "CLAUDE.md: document Task 3.A Stage 1; mark Task 2 I-profile todo DONE"')
    print('  git push')
    return 0


if __name__ == '__main__':
    sys.exit(main())
