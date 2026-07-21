#!/usr/bin/env python3
"""
Patch: adiabat is computed over the DT ice layer only (Region 2).

Previous behavior: density-defined cold shell spanning all fuel-candidate
zones (ri[:, 0] to ri[:, -2]) -- for 4-region targets like Olson PDD_26b
this includes ice + DT-CH foam payload. The foam compression pulled the
mass-weighted adiabat down below 1.0.

New behavior (per Tom's direction, April 2026): for capsule targets with
3+ regions, evaluate adiabat strictly over the inner DT ice layer
(ri[:, 0] to ri[:, 1], i.e. zones between the gas/fuel interface and the
ice/foam interface). No density mask -- mass-weighted average over the
full ice layer including any zones that partially ablated.

Rationale: the DT ice is the thermonuclear fuel that has to ignite. Olson
2021's reference adiabat 3.0 ± 0.5 is for the DT ice layer specifically.
Including the foam in the average produces a definitionally different
number that can't be compared cleanly to the published reference.

For 2-region targets (no ablator) the old behavior is preserved -- ice
spans ri[:, 0] to ri[:, -1] which is the same as before.

For single-region targets the existing skip is preserved.

Run from /Users/mehlhorn/Codes/helios_postprocess:
    python3 apply_patch_adiabat_ice_only.py
"""
from pathlib import Path
import sys

ia_path = Path("helios_postprocess/icf_analysis.py")
if not ia_path.is_file():
    print(f"ERROR: {ia_path} not found. Run from repo root.")
    sys.exit(1)

text = ia_path.read_text()

# ============================================================================
# Edit 1: _compute_adiabat (peak-velocity adiabat)
# ============================================================================
old_compute_adiabat = '''            # -- Zone selection: density-based shell (same criterion as IFAR) --
            # Cold fuel = zones where rho > rho_peak / e, restricted to the
            # fuel-candidate range (between gas/fuel interface and the fuel/ablator
            # interface). This is target-geometry-agnostic: if a nominally-fuel
            # region has partially ablated, its low-density zones fall below the
            # density threshold and are correctly excluded. Works for 2-region,
            # 3-region, or 5-region targets without special casing.
            if ri is not None and ri.shape[1] >= 3:
                fuel_lo = int(ri[eval_idx, 0])            # gas/fuel interface
                fuel_hi = int(ri[eval_idx, -2])           # fuel/ablator interface (outermost fuel)
            elif ri is not None and ri.shape[1] == 2:
                fuel_lo = int(ri[eval_idx, 0])
                fuel_hi = int(ri[eval_idx, -1])           # no ablator: fuel to outer edge
            else:
                fuel_lo = 0
                fuel_hi = n_zones'''

new_compute_adiabat = '''            # -- Zone selection: DT ice layer only (Region 2) --
            # Per Olson 2021 convention, adiabat is evaluated over the cold DT
            # ice fuel only -- not including any foam ablator layer. For 3+
            # region targets ri[:, 1] is the ice/foam interface (outer edge of
            # ice); for 2-region targets there is no separate ablator so ice
            # extends to the outer boundary.
            if ri is not None and ri.shape[1] >= 3:
                fuel_lo = int(ri[eval_idx, 0])            # gas/fuel interface
                fuel_hi = int(ri[eval_idx, 1])            # ice/foam interface  (inner edge of foam)
            elif ri is not None and ri.shape[1] == 2:
                fuel_lo = int(ri[eval_idx, 0])
                fuel_hi = int(ri[eval_idx, -1])           # no ablator: fuel to outer edge
            else:
                fuel_lo = 0
                fuel_hi = n_zones'''

n = text.count(old_compute_adiabat)
if n == 0:
    if "fuel_hi = int(ri[eval_idx, 1])" in text and "DT ice layer only" in text:
        print("Edit 1: already applied -- skipping")
    else:
        print("ERROR: Edit 1 anchor not found in icf_analysis.py")
        sys.exit(1)
elif n > 1:
    print(f"ERROR: Edit 1 anchor matched {n} times")
    sys.exit(1)
else:
    text = text.replace(old_compute_adiabat, new_compute_adiabat, 1)
    print("Edit 1 applied: _compute_adiabat now uses ice-only zones")

# ============================================================================
# Edit 2: drop the density-mask shell selection in _compute_adiabat
# ============================================================================
old_density_mask = '''            if fuel_hi <= fuel_lo:
                logger.warning("Fuel-candidate range is empty -- cannot compute adiabat")
                return

            rho_candidate  = self.data.mass_density[eval_idx, fuel_lo:fuel_hi]
            rho_peak       = float(np.max(rho_candidate))
            rho_threshold  = rho_peak / np.e
            shell_mask_loc = rho_candidate > rho_threshold

            if not np.any(shell_mask_loc):
                logger.warning("Density-based shell is empty -- cannot compute adiabat")
                return

            # -- Plasma pressure (ion + electron) in Mbar over the shell mask --
            # Conventional adiabat uses plasma pressure, not plasma+radiation.
            p_tot = self.data.plasma_pressure[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]
            p_Mbar = p_tot * 1e-5                        # J/cm3 -> Mbar

            rho  = rho_candidate[shell_mask_loc]
            mass = self.data.zone_mass[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]

            logger.info(
                f"Adiabat shell (density-based): {int(np.sum(shell_mask_loc))} zones "
                f"[within fuel-candidate zones {fuel_lo}..{fuel_hi - 1}], "
                f"rho_peak={rho_peak:.3f} g/cc, threshold={rho_threshold:.3f} g/cc"
            )'''

new_no_density_mask = '''            if fuel_hi <= fuel_lo:
                logger.warning("Ice region is empty -- cannot compute adiabat")
                return

            # -- Plasma pressure (ion + electron) in Mbar over the ice layer --
            # Mass-weighted average over all ice zones (no density mask).
            # This matches the Olson 2021 convention: the DT ice payload is
            # the thermonuclear fuel that needs to ignite; its adiabat is the
            # diagnostic of interest regardless of partial-ablation state.
            p_tot = self.data.plasma_pressure[eval_idx, fuel_lo:fuel_hi]
            p_Mbar = p_tot * 1e-5                        # J/cm3 -> Mbar

            rho  = self.data.mass_density[eval_idx, fuel_lo:fuel_hi]
            mass = self.data.zone_mass[eval_idx, fuel_lo:fuel_hi]

            n_ice = fuel_hi - fuel_lo
            rho_peak_ice = float(np.max(rho)) if n_ice > 0 else 0.0
            logger.info(
                f"Adiabat (DT ice only): {n_ice} zones [{fuel_lo}..{fuel_hi - 1}], "
                f"rho_peak_ice={rho_peak_ice:.3f} g/cc"
            )'''

n = text.count(old_density_mask)
if n == 0:
    if "Adiabat (DT ice only)" in text:
        print("Edit 2: already applied -- skipping")
    else:
        print("ERROR: Edit 2 anchor not found")
        sys.exit(1)
elif n > 1:
    print(f"ERROR: Edit 2 anchor matched {n} times")
    sys.exit(1)
else:
    text = text.replace(old_density_mask, new_no_density_mask, 1)
    print("Edit 2 applied: density mask removed in _compute_adiabat")

# ============================================================================
# Edit 3: drop the cold_mask comment block (it was already simplified earlier
# in our session but the comment may still reference density-shell logic)
# ============================================================================
# Idempotent: only apply if old comment present
old_coldmask_comment = '''            # Density-based shell selection above already ensures we have
            # compressed cold fuel (rho > rho_peak / e), so no additional
            # cold_mask is needed. Averaging over the full shell.
            self.data.adiabat_mass_averaged_ice = float(np.average(alpha, weights=mass))'''

new_coldmask_comment = '''            # Mass-weighted average over the DT ice layer.
            self.data.adiabat_mass_averaged_ice = float(np.average(alpha, weights=mass))'''

if old_coldmask_comment in text:
    text = text.replace(old_coldmask_comment, new_coldmask_comment, 1)
    print("Edit 3 applied: simplified comment in _compute_adiabat")

# ============================================================================
# Edit 4: _compute_adiabat_at_breakout (base adiabat at shock breakout)
# ============================================================================
old_breakout_zones = '''            # Density-based shell at the breakout timestep (same logic as _compute_adiabat)
            # Single-region targets have no fuel layer -- skip cleanly.
            if ri is not None and ri.shape[1] < 2:
                logger.info("Base adiabat: single-region target, no fuel layer -- skipping")
                return

            if ri is not None and ri.shape[1] >= 3:
                fuel_lo = int(ri[idx_b, 0])
                fuel_hi = int(ri[idx_b, -2])
            elif ri is not None and ri.shape[1] == 2:
                fuel_lo = int(ri[idx_b, 0])
                fuel_hi = int(ri[idx_b, -1])
            else:
                fuel_lo = 0
                fuel_hi = n_zones

            if fuel_hi <= fuel_lo:
                return

            rho_candidate  = self.data.mass_density[idx_b, fuel_lo:fuel_hi]
            rho_peak       = float(np.max(rho_candidate))
            rho_threshold  = rho_peak / np.e
            shell_mask_loc = rho_candidate > rho_threshold

            if not np.any(shell_mask_loc):
                return

            # Plasma pressure (ion + electron) -- conventional for adiabat
            p_tot = self.data.plasma_pressure[idx_b, fuel_lo:fuel_hi][shell_mask_loc]
            p_Mbar = p_tot * 1e-5

            rho  = rho_candidate[shell_mask_loc]
            mass = self.data.zone_mass[idx_b, fuel_lo:fuel_hi][shell_mask_loc]'''

new_breakout_zones = '''            # DT ice layer only (Region 2). Same convention as _compute_adiabat
            # but evaluated at the shock-breakout timestep idx_b.
            if ri is not None and ri.shape[1] < 2:
                logger.info("Base adiabat: single-region target, no fuel layer -- skipping")
                return

            if ri is not None and ri.shape[1] >= 3:
                fuel_lo = int(ri[idx_b, 0])               # gas/fuel interface
                fuel_hi = int(ri[idx_b, 1])               # ice/foam interface
            elif ri is not None and ri.shape[1] == 2:
                fuel_lo = int(ri[idx_b, 0])
                fuel_hi = int(ri[idx_b, -1])              # no ablator: fuel to outer edge
            else:
                fuel_lo = 0
                fuel_hi = n_zones

            if fuel_hi <= fuel_lo:
                return

            # Plasma pressure (ion + electron) over the full ice layer -- no density mask
            p_tot = self.data.plasma_pressure[idx_b, fuel_lo:fuel_hi]
            p_Mbar = p_tot * 1e-5

            rho  = self.data.mass_density[idx_b, fuel_lo:fuel_hi]
            mass = self.data.zone_mass[idx_b, fuel_lo:fuel_hi]'''

n = text.count(old_breakout_zones)
if n == 0:
    if "DT ice layer only (Region 2). Same convention as _compute_adiabat" in text:
        print("Edit 4: already applied -- skipping")
    else:
        print("ERROR: Edit 4 anchor not found")
        sys.exit(1)
elif n > 1:
    print(f"ERROR: Edit 4 anchor matched {n} times")
    sys.exit(1)
else:
    text = text.replace(old_breakout_zones, new_breakout_zones, 1)
    print("Edit 4 applied: _compute_adiabat_at_breakout uses ice-only zones")

# ============================================================================
# Edit 5: _compute_adiabat_at_breakout log line update for ice-only convention
# ============================================================================
old_breakout_log = '''            self.data.adiabat_at_breakout = float(np.average(alpha, weights=mass))

            logger.info(
                f"Base adiabat at breakout (t={self.data.shock_breakout_time_ns:.3f} ns): "
                f"{self.data.adiabat_at_breakout:.2f}  "
                f"[{int(np.sum(shell_mask_loc))} shell zones, "
                f"rho_peak={rho_peak:.3f} g/cc]"
            )'''

new_breakout_log = '''            self.data.adiabat_at_breakout = float(np.average(alpha, weights=mass))

            n_ice = fuel_hi - fuel_lo
            rho_peak_ice = float(np.max(rho)) if n_ice > 0 else 0.0
            logger.info(
                f"Base adiabat at breakout (t={self.data.shock_breakout_time_ns:.3f} ns): "
                f"{self.data.adiabat_at_breakout:.2f}  "
                f"[{n_ice} ice zones, rho_peak_ice={rho_peak_ice:.3f} g/cc]"
            )'''

n = text.count(old_breakout_log)
if n == 0:
    if "ice zones, rho_peak_ice" in text:
        print("Edit 5: already applied -- skipping")
    else:
        print("ERROR: Edit 5 anchor not found  (the log line may differ from expected)")
        # Don't exit -- this is a cosmetic edit; the calculation is already correct
        print("       (This is cosmetic only -- skipping but proceeding)")
elif n > 1:
    print(f"ERROR: Edit 5 anchor matched {n} times")
    sys.exit(1)
else:
    text = text.replace(old_breakout_log, new_breakout_log, 1)
    print("Edit 5 applied: log message updated")

ia_path.write_text(text)

# ============================================================================
# Verify on disk
# ============================================================================
verify = ia_path.read_text()
checks = [
    ("Edit 1 disk-verify", "fuel_hi = int(ri[eval_idx, 1])"),
    ("Edit 2 disk-verify", "Adiabat (DT ice only)"),
    ("Edit 4 disk-verify", "fuel_hi = int(ri[idx_b, 1])"),
]
for name, marker in checks:
    if marker in verify:
        print(f"  [verified] {name}")
    else:
        print(f"  [MISSING]  {name}: {marker!r}")
        sys.exit(1)

# Also confirm the foam-inclusive logic is gone
bad_markers = [
    "ri[eval_idx, -2]",   # old upper bound for ice
    "ri[idx_b, -2]",      # old upper bound for ice
]
for marker in bad_markers:
    if marker in verify:
        print(f"  [WARNING]  Old foam-inclusive marker still present: {marker!r}")
        # Continue rather than abort — there may be other valid uses

print()
print("Patch complete. Adiabat is now ice-layer-only for multi-region targets.")
print()
print("Next:")
print("  git diff helios_postprocess/icf_analysis.py | head -80")
print("  git add helios_postprocess/icf_analysis.py")
print("  git commit -m 'Adiabat: evaluate over DT ice layer only (Olson 2021 convention)'")
print("  git push")
