#!/usr/bin/env python3
"""
Patch: redefine shock_foot_pressure_Gbar as spatial peak at breakout instant.

Previously: max over (time, space) of plasma pressure up to and including
            breakout. Captured the early-time transient ablation spike
            (~150 Mbar for CH sphere at t ~ 0.13 ns).

Now:        max over space at the breakout instant t_b. Captures the
            quasi-steady ablation drive at the moment the shock reaches
            the rear face. Matches typical ICF "ablation pressure"
            reporting convention (~107 Mbar for CH sphere at t = 0.83 ns,
            matching Tom's published Helios reference plot).

Also updates the diagnostic name in the summary output:
    'Ice-side foot peak pressure'  ->  'Ablation pressure at breakout'
This is more accurate -- it's the spatial peak at breakout, not the
ice-side probe history peak.

Run from /Users/mehlhorn/Codes/helios_postprocess:
    python3 apply_patch_ablation_at_breakout.py
"""
from pathlib import Path
import sys

ia_path = Path("helios_postprocess/icf_analysis.py")
io_path = Path("helios_postprocess/icf_output.py")

for f in [ia_path, io_path]:
    if not f.is_file():
        print(f"ERROR: {f} not found. Run from repo root.")
        sys.exit(1)

# ============================================================================
# Edit 1: redefine the ablation-pressure calculation in _compute_shock_breakout
# ============================================================================
ia = ia_path.read_text()

old_calc = '''            # Ablation pressure = spatial peak of the pressure field at any
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

new_calc = '''            # Ablation pressure at breakout instant = spatial peak of the
            # plasma pressure field AT the breakout timestep i_b. This is
            # the quasi-steady ablation drive at the moment the shock
            # reaches the rear face of the shell, which matches the
            # typical ICF reporting convention. Time-of-peak ablation
            # (early in the laser pulse during turn-on) is generally
            # higher than this value (~30-40% for square-pulse drives on
            # CH); the current convention reports the steady-state value
            # at a well-defined physical event (shock breakout) instead.
            self.data.shock_foot_pressure_Gbar = float(np.max(total_P[i_b, :]) * 1e-8)'''

if new_calc in ia and old_calc not in ia:
    print("Edit 1: already applied -- skipping")
else:
    count = ia.count(old_calc)
    if count == 0:
        print("ERROR: anchor not found in icf_analysis.py")
        sys.exit(1)
    if count > 1:
        print(f"ERROR: anchor matched {count} times (expected 1)")
        sys.exit(1)
    ia = ia.replace(old_calc, new_calc, 1)
    ia_path.write_text(ia)
    print("Edit 1 applied: shock_foot_pressure_Gbar now spatial peak at breakout instant")

# ============================================================================
# Edit 2: update the summary output label to match the new definition
# ============================================================================
io = io_path.read_text()

old_label = "        if getattr(d, 'shock_foot_pressure_Gbar', 0.0) > 0:\n            _a(self._metric('Ice-side foot peak pressure',  1000.0 * d.shock_foot_pressure_Gbar,                   'Mbar', fmt='.2f'))"

new_label = "        if getattr(d, 'shock_foot_pressure_Gbar', 0.0) > 0:\n            _a(self._metric('Ablation pressure at breakout', 1000.0 * d.shock_foot_pressure_Gbar,                  'Mbar', fmt='.2f'))"

if new_label in io and old_label not in io:
    print("Edit 2: already applied -- skipping")
else:
    count = io.count(old_label)
    if count == 0:
        print("ERROR: anchor not found in icf_output.py")
        print("Expected to find: 'Ice-side foot peak pressure' label.")
        sys.exit(1)
    if count > 1:
        print(f"ERROR: anchor matched {count} times (expected 1)")
        sys.exit(1)
    io = io.replace(old_label, new_label, 1)
    io_path.write_text(io)
    print("Edit 2 applied: summary label updated to 'Ablation pressure at breakout'")

# ============================================================================
# Verify
# ============================================================================
print()
print("Verifying on disk...")
ia_v = ia_path.read_text()
io_v = io_path.read_text()
assert "spatial peak of the\n            # plasma pressure field AT the breakout" in ia_v, "verify failed: new comment missing"
assert "shock_foot_pressure_Gbar = float(np.max(total_P[i_b, :])" in ia_v, "verify failed: new calc missing"
assert "Ablation pressure at breakout" in io_v, "verify failed: new label missing"
assert "Ice-side foot peak pressure" not in io_v, "verify failed: old label still present"
print("Verified.")
print()
print("Next:")
print("  git diff helios_postprocess/ | head -60")
print("  git add helios_postprocess/")
print("  git commit -m 'Ablation pressure: report spatial peak at breakout instant (matches ICF convention)'")
print("  git push")
