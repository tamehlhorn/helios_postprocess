#!/usr/bin/env python3
"""
Adiabat cleanup (combined): finishes the ice-only patch and re-applies
plasma-pressure where it didn't take in the two adiabat methods.

Reads the current on-disk state of icf_analysis.py and converts both
_compute_adiabat() and _compute_adiabat_at_breakout() to:

  - DT ice layer only (ri[:, 0] : ri[:, 1] for capsules with 3+ regions)
  - Plasma pressure (ion + electron, no radiation)
  - No density mask -- mass-weighted average over full ice

Anchors match the actual byte-exact text on disk as of the diagnostic in
this conversation thread.

Run from /Users/mehlhorn/Codes/helios_postprocess:
    python3 apply_patch_adiabat_ice_only_v2.py
"""
from pathlib import Path
import sys

ia_path = Path("helios_postprocess/icf_analysis.py")
if not ia_path.is_file():
    print(f"ERROR: {ia_path} not found. Run from repo root.")
    sys.exit(1)

text = ia_path.read_text()

# ============================================================================
# Edit 1: replace the density-mask + pressure-computation block in
# _compute_adiabat() with ice-only + plasma-pressure
# ============================================================================
old_compute_adiabat = '''            rho_candidate  = self.data.mass_density[eval_idx, fuel_lo:fuel_hi]
            rho_peak       = float(np.max(rho_candidate))
            rho_threshold  = rho_peak / np.e
            shell_mask_loc = rho_candidate > rho_threshold

            if not np.any(shell_mask_loc):
                logger.warning("Density-based shell is empty -- cannot compute adiabat")
                return

            # -- Total pressure in Mbar over the shell mask --
            p_tot = self.data.ion_pressure[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]
            if self.data.rad_pressure is not None:
                p_tot = p_tot + self.data.rad_pressure[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]
            p_Mbar = p_tot * 1e-5                        # J/cm3 -> Mbar

            rho  = rho_candidate[shell_mask_loc]
            mass = self.data.zone_mass[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]

            logger.info(
                f"Adiabat shell (density-based): {int(np.sum(shell_mask_loc))} zones "
                f"[within fuel-candidate zones {fuel_lo}..{fuel_hi - 1}], "
                f"rho_peak={rho_peak:.3f} g/cc, threshold={rho_threshold:.3f} g/cc"
            )'''

new_compute_adiabat = '''            # -- Plasma pressure (ion + electron) in Mbar over the DT ice layer --
            # Mass-weighted average over the full ice region (no density mask).
            # Per Olson 2021 convention: the DT ice payload is the
            # thermonuclear fuel that needs to ignite; its adiabat is the
            # diagnostic of interest regardless of partial-ablation state.
            p_tot  = self.data.plasma_pressure[eval_idx, fuel_lo:fuel_hi]
            p_Mbar = p_tot * 1e-5                        # J/cm3 -> Mbar

            rho  = self.data.mass_density[eval_idx, fuel_lo:fuel_hi]
            mass = self.data.zone_mass[eval_idx, fuel_lo:fuel_hi]

            n_ice = fuel_hi - fuel_lo
            rho_peak_ice = float(np.max(rho)) if n_ice > 0 else 0.0
            logger.info(
                f"Adiabat (DT ice only, plasma pressure): {n_ice} zones "
                f"[{fuel_lo}..{fuel_hi - 1}], rho_peak_ice={rho_peak_ice:.3f} g/cc"
            )'''

n = text.count(old_compute_adiabat)
if n == 0:
    if "Adiabat (DT ice only, plasma pressure)" in text:
        print("Edit 1: already applied -- skipping")
    else:
        print("ERROR: Edit 1 anchor not found")
        sys.exit(1)
elif n > 1:
    print(f"ERROR: Edit 1 anchor matched {n} times")
    sys.exit(1)
else:
    text = text.replace(old_compute_adiabat, new_compute_adiabat, 1)
    print("Edit 1 applied: _compute_adiabat -> ice-only + plasma_pressure, no density mask")

# ============================================================================
# Edit 2: simplify the cold-mask-comment block in _compute_adiabat
# ============================================================================
# After Edit 1 above, the comment block describing the no-longer-relevant
# cold_mask is misleading. Simplify it.
old_coldmask = '''            # Density-based shell selection above already ensures we have
            # compressed cold fuel (rho > rho_peak / e), so no additional
            # cold_mask is needed. Averaging over the full shell.
            self.data.adiabat_mass_averaged_ice = float(np.average(alpha, weights=mass))'''

new_coldmask = '''            # Mass-weighted average over the DT ice layer.
            self.data.adiabat_mass_averaged_ice = float(np.average(alpha, weights=mass))'''

if old_coldmask in text:
    text = text.replace(old_coldmask, new_coldmask, 1)
    print("Edit 2 applied: simplified comment in _compute_adiabat")
elif "Mass-weighted average over the DT ice layer" in text:
    print("Edit 2: already applied -- skipping")
else:
    print("Edit 2: anchor not found (cosmetic only -- proceeding)")

# ============================================================================
# Edit 3: _compute_adiabat_at_breakout -- ice-only zone selection
# ============================================================================
# The fuel_lo/fuel_hi block is at lines ~661-666; on disk it currently uses
# ri[:, -2] for the upper bound (foam-inclusive). Change to ri[:, 1] for ice
# only.
old_breakout_zones = '''            if ri is not None and ri.shape[1] >= 3:
                fuel_lo = int(ri[idx_b, 0])
                fuel_hi = int(ri[idx_b, -2])
            elif ri is not None and ri.shape[1] == 2:
                fuel_lo = int(ri[idx_b, 0])
                fuel_hi = int(ri[idx_b, -1])
            else:
                fuel_lo = 0
                fuel_hi = n_zones'''

new_breakout_zones = '''            # DT ice layer only (Region 2). Same convention as _compute_adiabat.
            if ri is not None and ri.shape[1] >= 3:
                fuel_lo = int(ri[idx_b, 0])               # gas/fuel interface
                fuel_hi = int(ri[idx_b, 1])               # ice/foam interface
            elif ri is not None and ri.shape[1] == 2:
                fuel_lo = int(ri[idx_b, 0])
                fuel_hi = int(ri[idx_b, -1])              # no ablator: fuel to outer edge
            else:
                fuel_lo = 0
                fuel_hi = n_zones'''

n = text.count(old_breakout_zones)
if n == 0:
    if "fuel_hi = int(ri[idx_b, 1])" in text:
        print("Edit 3: already applied -- skipping")
    else:
        print("ERROR: Edit 3 anchor not found")
        sys.exit(1)
elif n > 1:
    print(f"ERROR: Edit 3 anchor matched {n} times")
    sys.exit(1)
else:
    text = text.replace(old_breakout_zones, new_breakout_zones, 1)
    print("Edit 3 applied: _compute_adiabat_at_breakout uses ice-only zones")

# ============================================================================
# Edit 4: _compute_adiabat_at_breakout -- ice-only + plasma-pressure body
# ============================================================================
# The pressure block on disk uses ion+rad and the density mask. Replace
# with plasma_pressure over the full ice layer.
old_breakout_body = '''            rho_candidate  = self.data.mass_density[idx_b, fuel_lo:fuel_hi]
            rho_peak       = float(np.max(rho_candidate))
            rho_threshold  = rho_peak / np.e
            shell_mask_loc = rho_candidate > rho_threshold

            if not np.any(shell_mask_loc):
                return

            p_tot = self.data.ion_pressure[idx_b, fuel_lo:fuel_hi][shell_mask_loc]
            if self.data.rad_pressure is not None:
                p_tot = p_tot + self.data.rad_pressure[idx_b, fuel_lo:fuel_hi][shell_mask_loc]
            p_Mbar = p_tot * 1e-5

            rho  = rho_candidate[shell_mask_loc]
            mass = self.data.zone_mass[idx_b, fuel_lo:fuel_hi][shell_mask_loc]'''

new_breakout_body = '''            # Plasma pressure (ion + electron) over the full DT ice layer.
            p_tot  = self.data.plasma_pressure[idx_b, fuel_lo:fuel_hi]
            p_Mbar = p_tot * 1e-5

            rho  = self.data.mass_density[idx_b, fuel_lo:fuel_hi]
            mass = self.data.zone_mass[idx_b, fuel_lo:fuel_hi]'''

n = text.count(old_breakout_body)
if n == 0:
    if "Plasma pressure (ion + electron) over the full DT ice layer" in text:
        print("Edit 4: already applied -- skipping")
    else:
        print("ERROR: Edit 4 anchor not found")
        sys.exit(1)
elif n > 1:
    print(f"ERROR: Edit 4 anchor matched {n} times")
    sys.exit(1)
else:
    text = text.replace(old_breakout_body, new_breakout_body, 1)
    print("Edit 4 applied: _compute_adiabat_at_breakout body uses plasma_pressure")

# ============================================================================
# Edit 5: log message at end of _compute_adiabat_at_breakout
# ============================================================================
old_breakout_log = '''            logger.info(
                f"Base adiabat at breakout (t={self.data.shock_breakout_time_ns:.3f} ns): "
                f"{self.data.adiabat_at_breakout:.2f}  "
                f"[{int(np.sum(shell_mask_loc))} shell zones, "
                f"rho_peak={rho_peak:.3f} g/cc]"
            )'''

new_breakout_log = '''            n_ice = fuel_hi - fuel_lo
            rho_peak_ice = float(np.max(rho)) if n_ice > 0 else 0.0
            logger.info(
                f"Base adiabat at breakout (t={self.data.shock_breakout_time_ns:.3f} ns): "
                f"{self.data.adiabat_at_breakout:.2f}  "
                f"[{n_ice} ice zones, rho_peak_ice={rho_peak_ice:.3f} g/cc]"
            )'''

if old_breakout_log in text:
    text = text.replace(old_breakout_log, new_breakout_log, 1)
    print("Edit 5 applied: breakout adiabat log message updated")
elif "ice zones, rho_peak_ice" in text:
    print("Edit 5: already applied -- skipping")
else:
    print("Edit 5: anchor not found (cosmetic only -- proceeding)")

ia_path.write_text(text)

# ============================================================================
# Verify on disk
# ============================================================================
verify = ia_path.read_text()
checks = [
    ("Edit 1 disk-verify (ice-only adiabat with plasma_pressure)",
     "Adiabat (DT ice only, plasma pressure)"),
    ("Edit 3 disk-verify (breakout uses ri[:, 1])",
     "fuel_hi = int(ri[idx_b, 1])"),
    ("Edit 4 disk-verify (breakout uses plasma_pressure)",
     "Plasma pressure (ion + electron) over the full DT ice layer"),
]
all_ok = True
for name, marker in checks:
    if marker in verify:
        print(f"  [verified] {name}")
    else:
        print(f"  [MISSING]  {name}")
        all_ok = False

# Confirm the two old call sites no longer use ion+rad
remaining_old = verify.count("self.data.ion_pressure[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]")
remaining_old += verify.count("self.data.ion_pressure[idx_b, fuel_lo:fuel_hi][shell_mask_loc]")
print(f"  Remaining ion+rad call sites in adiabat methods: {remaining_old} (expected 0)")

# Count plasma_pressure references in adiabat methods
n_plasma = verify.count("self.data.plasma_pressure[eval_idx")
n_plasma += verify.count("self.data.plasma_pressure[idx_b")
print(f"  plasma_pressure references in adiabat methods: {n_plasma} (expected 2)")

if not all_ok or remaining_old > 0:
    print()
    print("WARNING: state may be inconsistent. Inspect with:")
    print("  grep -n 'plasma_pressure\\|ion_pressure\\[\\(eval_idx\\|idx_b\\)' helios_postprocess/icf_analysis.py")
    sys.exit(1)

print()
print("Adiabat cleanup complete. Both _compute_adiabat and _compute_adiabat_at_breakout")
print("now use ice-only zones with plasma pressure, no density mask.")
print()
print("Next:")
print("  git diff helios_postprocess/icf_analysis.py | head -80")
print("  git add helios_postprocess/icf_analysis.py")
print("  git commit -m 'Adiabat: ice-only zones with plasma pressure (Olson 2021 convention)'")
print("  git push")
