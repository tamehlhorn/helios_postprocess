#!/usr/bin/env python
"""
patch_task3a_stage1.py -- Task 3.A Stage 1: peak velocity time + ablation pressure.

Three additive edits (no rewrites, no deletions):

  1a. ICFRunData.__init__ in data_builder.py:
      Declare three new attributes so they exist as None before analysis
      populates them. Prevents AttributeError if user reads before analyze.

      + self.t_peak_velocity_ns: Optional[float] = None
      + self.ablation_pressure_Gbar: Optional[np.ndarray] = None
      + self.P_abl_peak_Gbar: Optional[float] = None

  1b. analyze_implosion_phase in icf_analysis.py:
      After peak_velocity_index is computed, save the corresponding time.
      Expose as data.t_peak_velocity_ns (data.time is already in ns).

  1c. _track_ablation_front in icf_analysis.py:
      After ablation_front_indices is populated, compute pressure at the
      ablation front for each timestep as total = (ion + rad) pressure,
      convert J/cm^3 -> Gbar, store history + peak scalar.

  1d. icf_output.py TIMING block:
      Add "Peak velocity time" row after Bang time.

  1e. icf_output.py IMPLOSION block:
      Add "Peak ablation pressure" row (placement adjacent to peak velocity).

Idempotent via sentinel strings. All anchors include before+after context
lines so indentation is locked in.

Run from repo root:
    python patch_task3a_stage1.py
    python patch_task3a_stage1.py --dry-run
"""
import sys
import argparse
from pathlib import Path


def apply(src, label, old, new, sentinel=None):
    """Return (new_src, applied_bool). Refuses on wrong-count or already-applied."""
    if sentinel and sentinel in src:
        print(f'  [ALREADY] {label}')
        return src, False
    n = src.count(old)
    if n == 0:
        print(f'  [MISS]    {label}: OLD pattern not found (0 hits)')
        return src, False
    if n > 1:
        print(f'  [AMBIG]   {label}: OLD pattern found {n} times; refusing')
        return src, False
    print(f'  [OK]      {label}')
    return src.replace(old, new, 1), True


def patch_file(target: Path, edits: list, dry_run: bool):
    """edits: list of (label, old, new, sentinel) tuples."""
    print(f'\n=== {target.name} ===')
    if not target.exists():
        print(f'  [ABORT] file not found: {target}')
        return False
    src = target.read_text()
    orig = src
    any_changed = False
    for label, old, new, sentinel in edits:
        src, changed = apply(src, label, old, new, sentinel)
        any_changed = any_changed or changed
    if not any_changed:
        print('  (no edits applied; file unchanged)')
        return True
    if dry_run:
        print(f'  [DRY-RUN] would write {target.name}')
        return True
    target.write_text(src)
    print(f'  [WROTE]  {target.name}')
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default='.',
                    help='Repo root (default: current directory)')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    repo = Path(args.repo).expanduser().resolve()
    print(f'Repo: {repo}')

    db  = repo / 'helios_postprocess' / 'data_builder.py'
    an  = repo / 'helios_postprocess' / 'icf_analysis.py'
    out = repo / 'helios_postprocess' / 'icf_output.py'

    # ================================================================
    # 1a. data_builder.py: declare three new attributes on ICFRunData
    # ================================================================
    # Anchor: the existing peak_velocity_index line, which should be
    # adjacent to where we add t_peak_velocity_ns. We extend the
    # dataclass declaration by inserting three new attributes after
    # an existing implosion-metric attribute.
    #
    # We'll anchor on a distinctive existing declaration and add the
    # new three lines immediately after it. Use a common pattern: the
    # peak_implosion_velocity attribute declaration.

    DB_OLD = (
        "    peak_implosion_velocity: Optional[float] = None"
    )
    DB_NEW = (
        "    peak_implosion_velocity: Optional[float] = None\n"
        "    t_peak_velocity_ns: Optional[float] = None\n"
        "    ablation_pressure_Gbar: Optional[np.ndarray] = None\n"
        "    P_abl_peak_Gbar: Optional[float] = None"
    )
    DB_SENTINEL = "t_peak_velocity_ns: Optional[float] = None"

    # ================================================================
    # 1b. icf_analysis.py: save t_peak_velocity_ns in analyze_implosion
    # ================================================================
    # Anchor: two lines around the peak_implosion_velocity assignment.
    # Include the logger.info that follows for additional uniqueness.

    AN_OLD_PEAK = (
        "            self.data.peak_velocity_index = int(pre_bang_indices[peak_idx_in_prebang])\n"
        "            self.data.peak_implosion_velocity = float(min_velocities[peak_idx_in_prebang]) * 1e-5\n"
        "            logger.info(f\"Peak implosion velocity (pre-bang): \"\n"
        "                        f\"{abs(self.data.peak_implosion_velocity):.2f} km/s\")"
    )
    AN_NEW_PEAK = (
        "            self.data.peak_velocity_index = int(pre_bang_indices[peak_idx_in_prebang])\n"
        "            self.data.peak_implosion_velocity = float(min_velocities[peak_idx_in_prebang]) * 1e-5\n"
        "            # Time of peak velocity (data.time is already in ns per data_builder)\n"
        "            self.data.t_peak_velocity_ns = float(self.data.time[self.data.peak_velocity_index])\n"
        "            logger.info(f\"Peak implosion velocity (pre-bang): \"\n"
        "                        f\"{abs(self.data.peak_implosion_velocity):.2f} km/s \"\n"
        "                        f\"at t={self.data.t_peak_velocity_ns:.3f} ns\")"
    )
    AN_PEAK_SENTINEL = "self.data.t_peak_velocity_ns = float(self.data.time"

    # ================================================================
    # 1c. icf_analysis.py: populate ablation_pressure_Gbar history
    # ================================================================
    # Anchor: the end of _track_ablation_front, right before the
    # next method definition (_compute_ifar). We add the pressure
    # computation just before the function ends.
    #
    # Use the two distinctive logger.info lines at the end of the
    # method as the anchor, then the blank line + next def.

    AN_OLD_ABL = (
        "            logger.info(f\"Ablation front tracked: {len(valid_radii)}/{n_times} timesteps\")\n"
        "            logger.info(f\"Initial radius: {valid_radii[0]:.4f} cm, \"\n"
        "                        f\"min radius: {np.min(valid_radii):.4f} cm\")\n"
        "        else:\n"
        "            logger.warning(\"No valid ablation front positions found\")\n"
        "    \n"
        "    def _compute_ifar(self):"
    )
    AN_NEW_ABL = (
        "            logger.info(f\"Ablation front tracked: {len(valid_radii)}/{n_times} timesteps\")\n"
        "            logger.info(f\"Initial radius: {valid_radii[0]:.4f} cm, \"\n"
        "                        f\"min radius: {np.min(valid_radii):.4f} cm\")\n"
        "        else:\n"
        "            logger.warning(\"No valid ablation front positions found\")\n"
        "\n"
        "        # ---- Ablation pressure history (Task 3.A Stage 1) ----\n"
        "        # Total pressure (ion + rad, where rad already includes electron+radiation\n"
        "        # per the 3-component convention in data_builder) at ablation-front zone.\n"
        "        # Converts J/cm^3 -> Gbar via the 1e-8 factor used throughout.\n"
        "        if (self.data.ion_pressure is not None\n"
        "                and self.data.rad_pressure is not None\n"
        "                and ablation_front_indices is not None):\n"
        "            total_p = self.data.ion_pressure + self.data.rad_pressure\n"
        "            n_t = len(self.data.time)\n"
        "            p_abl = np.zeros(n_t)\n"
        "            for t in range(n_t):\n"
        "                i_abl = int(ablation_front_indices[t])\n"
        "                if i_abl > 0:\n"
        "                    p_abl[t] = total_p[t, i_abl] * 1e-8  # J/cm^3 -> Gbar\n"
        "            self.data.ablation_pressure_Gbar = p_abl\n"
        "            valid_p = p_abl[p_abl > 0]\n"
        "            if len(valid_p) > 0:\n"
        "                self.data.P_abl_peak_Gbar = float(np.max(valid_p))\n"
        "                logger.info(f\"Peak ablation pressure: \"\n"
        "                            f\"{self.data.P_abl_peak_Gbar:.2f} Gbar\")\n"
        "            else:\n"
        "                self.data.P_abl_peak_Gbar = 0.0\n"
        "                logger.warning(\"Ablation pressure history is all zero\")\n"
        "\n"
        "    def _compute_ifar(self):"
    )
    AN_ABL_SENTINEL = "# ---- Ablation pressure history (Task 3.A Stage 1) ----"

    # ================================================================
    # 1d. icf_output.py: add peak velocity time to TIMING block
    # ================================================================
    # Anchor: TIMING _a() calls. Add new row after Bang time.
    # Full-statement anchor including both existing lines locks indent.

    OUT_OLD_TIMING = (
        "        _a(self._metric('Stagnation time',      d.stag_time,   'ns',   fmt='.3f'))\n"
        "        _a(self._metric('Bang time',            d.bang_time,   'ns',   fmt='.3f'))"
    )
    OUT_NEW_TIMING = (
        "        _a(self._metric('Stagnation time',      d.stag_time,   'ns',   fmt='.3f'))\n"
        "        _a(self._metric('Bang time',            d.bang_time,   'ns',   fmt='.3f'))\n"
        "        _a(self._metric('Peak velocity time',   d.t_peak_velocity_ns, 'ns', fmt='.3f'))"
    )
    OUT_TIMING_SENTINEL = "Peak velocity time"

    # ================================================================
    # 1e. icf_output.py: add peak ablation pressure row
    # ================================================================
    # For Stage 1 we place this adjacent to existing peak velocity line
    # in the IMPLOSION block. We need to find that line first via grep
    # next session if the anchor differs; for now we use a broad anchor
    # on "Peak implosion velocity" style label. If this anchor misses,
    # the patch safely no-ops and we deliver Stage 1.5 for output.
    #
    # NOTE: intentionally skipping 1e in this patch because we haven't
    # greped the exact IMPLOSION label format in icf_output.py. We'll
    # add it in a small follow-up patch once confirmed.

    # ---- Apply ----
    ok = True
    ok &= patch_file(db, [
        ('1a data_builder.py: declare ICFRunData.t_peak_velocity_ns + '
         'ablation_pressure_Gbar + P_abl_peak_Gbar',
         DB_OLD, DB_NEW, DB_SENTINEL),
    ], args.dry_run)

    ok &= patch_file(an, [
        ('1b icf_analysis.py: save t_peak_velocity_ns at peak',
         AN_OLD_PEAK, AN_NEW_PEAK, AN_PEAK_SENTINEL),
        ('1c icf_analysis.py: populate ablation_pressure_Gbar history',
         AN_OLD_ABL, AN_NEW_ABL, AN_ABL_SENTINEL),
    ], args.dry_run)

    ok &= patch_file(out, [
        ('1d icf_output.py: TIMING peak velocity time row',
         OUT_OLD_TIMING, OUT_NEW_TIMING, OUT_TIMING_SENTINEL),
    ], args.dry_run)

    if not ok:
        print('\n  [FAIL] At least one edit did not apply cleanly.')
        return 1

    print('\n-----------------------------------------------------------------------')
    if args.dry_run:
        print('Dry run complete. No files written. Rerun without --dry-run to apply.')
    else:
        print('Stage 1 applied. Next steps:')
        print('  git diff --stat helios_postprocess/*.py')
        print('  git diff helios_postprocess/*.py | less')
        print('')
        print('  # On Mac Studio after commit:')
        print('  python3 examples/run_analysis.py \\\\')
        print('    ~/Sims/Xcimer/Olson_PDD/Olson_PDD_26b_burn/Olson_PDD_26b_burn')
        print('')
        print('  Expected new log lines:')
        print('    Peak implosion velocity (pre-bang): 478.36 km/s at t=12.00 ns')
        print('    Peak ablation pressure: ~300-500 Gbar (approximate)')
        print('')
        print('  Expected new TIMING line in summary:')
        print('    Peak velocity time       12.00 ns')
    print('-----------------------------------------------------------------------')
    return 0


if __name__ == '__main__':
    sys.exit(main())
