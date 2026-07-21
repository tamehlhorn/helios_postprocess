#!/usr/bin/env python3
"""
Patch A2: refine _compute_shock_breakout based on the diagnostic results
from CH_sphere_ca1_f06.

Two changes:

(1) Threshold floor: 1e-3 Gbar (1 Mbar) instead of 1e-4 Gbar (0.1 Mbar).
    The diagnostic showed P_rear has a clear two-stage rise -- a slow
    preheat ramp from 0 to ~0.4 Mbar over 0-0.8 ns, then a sharp jump to
    ~3 Mbar when the actual shock arrives at ~0.84 ns. Old floor of 0.1
    Mbar caught the preheat at 0.28 ns. New floor of 1 Mbar correctly
    catches shock arrival at 0.84 ns. Multiplier also bumped from 10x to
    100x P_rear[0] for the same reason.

(2) Ablation pressure = spatial peak of the field, not static drive-zone
    probe. The static probe (zones N-5..N-1) sees coronal blowoff, not
    the ablation front -- which moves inward through the Lagrangian grid
    as material ablates. The actual ablation pressure is the pressure
    AT the ablation front, wherever it happens to be. Take np.max over
    space at each timestep, then peak over time up to breakout.

Run from /Users/mehlhorn/Codes/helios_postprocess:
    python3 apply_patch_A2.py
"""
from pathlib import Path
import sys

ia_path = Path('helios_postprocess/icf_analysis.py')
if not ia_path.is_file():
    print(f"ERROR: {ia_path} not found. Run from repo root.")
    sys.exit(1)

ia = ia_path.read_text()

# Anchor: the threshold + drive-probe-peak block in the Patch A version.
# This is what's currently on disk.
old = '''            # Threshold: 10x initial rear-probe pressure, with absolute floor of
            # 1e-4 Gbar (0.1 kbar / 100 bar) so that pathologically low or zero
            # initial values don't trigger false positives on numerical noise.
            P_rear_initial = float(P_rear[0])
            P_threshold = max(10.0 * P_rear_initial, 1e-4)   # Gbar

            # Time floor: after foot launch, to skip IC transients
            t_floor = getattr(self.data, "laser_foot_start_ns", None) or 0.1

            hits = np.where((t_ns >= t_floor) & (P_rear > P_threshold))[0]
            if len(hits) == 0:
                logger.info(
                    f"Shock breakout: not detected ({geometry_label}); "
                    f"max rear-probe P = {np.max(P_rear):.4e} Gbar, "
                    f"threshold = {P_threshold:.4e} Gbar"
                )
                return

            i_b = int(hits[0])
            self.data.shock_breakout_index         = i_b
            self.data.shock_breakout_time_ns       = float(t_ns[i_b])
            self.data.shock_breakout_P_gas_Gbar    = float(P_rear[i_b])    # rear face
            self.data.shock_breakout_P_ice_Gbar    = float(P_drive[i_b])   # drive side
            self.data.shock_breakout_pressure_Gbar = float(P_rear[i_b])    # alias

            # Drive-side foot-pulse peak: max drive-probe pressure from t_floor
            # up to (and including) breakout. For capsules this is the foot
            # shock strength in the ice; for solid shells, the ablation drive
            # peak. The latter should approximately match the published
            # "ablation pressure" diagnostic (~110 Mbar for the CH sphere).
            pre_breakout_mask = (t_ns >= t_floor) & (t_ns <= t_ns[i_b])
            if np.any(pre_breakout_mask):
                self.data.shock_foot_pressure_Gbar = float(np.max(P_drive[pre_breakout_mask]))'''

new = '''            # Threshold: 100x initial rear-probe pressure, with absolute floor of
            # 1e-3 Gbar (1 Mbar). The 1 Mbar floor is calibrated to filter
            # out the radiation/electron-conduction preheat ramp (which is
            # typically 0.1-0.5 Mbar) and trigger only on actual shock
            # arrival. Diagnostic on CH sphere CA1 / f=0.06 showed clean
            # 0.84 ns breakout with this threshold (matches published 0.83 ns).
            P_rear_initial = float(P_rear[0])
            P_threshold = max(100.0 * P_rear_initial, 1e-3)   # Gbar

            # Time floor: after foot launch, to skip IC transients
            t_floor = getattr(self.data, "laser_foot_start_ns", None) or 0.1

            hits = np.where((t_ns >= t_floor) & (P_rear > P_threshold))[0]
            if len(hits) == 0:
                logger.info(
                    f"Shock breakout: not detected ({geometry_label}); "
                    f"max rear-probe P = {np.max(P_rear):.4e} Gbar, "
                    f"threshold = {P_threshold:.4e} Gbar"
                )
                return

            i_b = int(hits[0])
            self.data.shock_breakout_index         = i_b
            self.data.shock_breakout_time_ns       = float(t_ns[i_b])
            self.data.shock_breakout_P_gas_Gbar    = float(P_rear[i_b])    # rear face
            self.data.shock_breakout_P_ice_Gbar    = float(P_drive[i_b])   # drive side
            self.data.shock_breakout_pressure_Gbar = float(P_rear[i_b])    # alias

            # Ablation pressure = spatial peak of the pressure field at any
            # timestep up to breakout. The static drive-zone probe doesn't
            # work for solid-shell targets because the outer Lagrangian zones
            # blow off into corona; the actual ablation front travels inward
            # as material ablates. Taking np.max over space captures the
            # ablation front pressure regardless of where it currently sits.
            # For capsules this matches the foot-shock peak in the shell;
            # for solid shells, the steady-state ablation pressure (~110-150
            # Mbar at 1e15 W/cm^2 drive on CH).
            pre_breakout_mask = (t_ns >= t_floor) & (t_ns <= t_ns[i_b])
            if np.any(pre_breakout_mask):
                # Spatial peak at each timestep -> peak over the pre-breakout window
                P_field_max_t = np.max(total_P[pre_breakout_mask, :], axis=1) * 1e-8  # Gbar
                self.data.shock_foot_pressure_Gbar = float(np.max(P_field_max_t))'''

# Idempotent: detect already-applied state
if "100.0 * P_rear_initial" in ia:
    print("Patch A2 already applied -- nothing to do.")
    sys.exit(0)

count = ia.count(old)
if count == 0:
    print("ERROR: anchor not found. Patch A may have been modified since.")
    print("Inspect with:")
    print("  grep -A 6 'Threshold: 10x initial' helios_postprocess/icf_analysis.py")
    sys.exit(1)

assert count == 1, f"Anchor not unique (found {count} matches)"

ia = ia.replace(old, new)
ia_path.write_text(ia)

verify = Path('helios_postprocess/icf_analysis.py').read_text()
assert "100.0 * P_rear_initial" in verify, "Verify failed: new threshold not on disk"
assert "Spatial peak at each timestep" in verify, "Verify failed: new ablation calc not on disk"

print("Patch A2 applied and verified.")
print()
print("Next:")
print("  git diff helios_postprocess/icf_analysis.py | head -60")
print("  git add helios_postprocess/icf_analysis.py")
print("  git commit -m 'Patch A2: tune breakout threshold (1 Mbar floor) and use spatial-peak ablation pressure'")
print("  git push")
