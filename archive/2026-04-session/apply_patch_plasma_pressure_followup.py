#!/usr/bin/env python3
"""
Plasma-pressure patch — follow-up
=================================
Completes the remaining 4 call sites in icf_analysis.py that the previous
plasma_pressure patch couldn't anchor to (whitespace mismatches).

These are stagnation-era call sites (hot-spot pressure, internal energy,
ignition pressure, neutron-averaged pressure). The CH-sphere
flux-limiter test doesn't exercise any of these, but we want them
plasma-only too for consistency on capsule runs.

Run from /Users/mehlhorn/Codes/helios_postprocess:
    python3 apply_patch_plasma_pressure_followup.py

Idempotent: safe to re-run. Will report which edits are still needed.
"""
from pathlib import Path
import sys

ia = Path("helios_postprocess/icf_analysis.py")
if not ia.is_file():
    print(f"ERROR: {ia} not found. Run from repo root.")
    sys.exit(1)

text = ia.read_text()

# Each entry: (description, old_text, new_text)
# Anchors copied byte-exact from the file via sed.

edits = [
    # ----- Edit 5: hot-spot pressure at stagnation (line ~1271) -----
    # Note the trailing space after '+' and before newline on line 1271
    ("Hot-spot stagnation pressure: use plasma_pressure",
     "                pressure_Gbar = (self.data.ion_pressure[stag_idx] + \n"
     "                           self.data.rad_pressure[stag_idx]) * 1e-8  # J/cm³ → Gbar",
     "                # Plasma pressure (ion + electron) at stagnation -- conventional definition\n"
     "                pressure_Gbar = self.data.plasma_pressure[stag_idx] * 1e-8  # J/cm³ → Gbar"),

    # ----- Edit 6: hot-spot internal energy (line ~1287) -----
    ("Hot-spot internal energy: use plasma_pressure",
     "                if self.data.ion_pressure is not None:\n"
     "                    # E_int = P / (γ-1) × Volume  for ideal gas\n"
     "                    # Or directly: e_specific × mass  if we have SIE.\n"
     "                    # Use  E = (3/2) n k T × V  ≈ (3/2) P V  for each species\n"
     "                    vol = (4.0 / 3.0) * np.pi * (boundaries[1:]**3 - boundaries[:-1]**3)\n"
     "                    p_ion  = self.data.ion_pressure[stag_idx]\n"
     "                    p_elec = self.data.rad_pressure[stag_idx]  # includes elec + rad\n"
     "                    hs_energy_J = np.sum((1.5 * (p_ion[hot_spot_mask] + p_elec[hot_spot_mask]))\n"
     "                                         * vol[hot_spot_mask])",
     "                if self.data.plasma_pressure is not None:\n"
     "                    # E_int = (3/2) P V  for ideal gas  (P = plasma pressure = ion + electron)\n"
     "                    vol = (4.0 / 3.0) * np.pi * (boundaries[1:]**3 - boundaries[:-1]**3)\n"
     "                    p_plasma = self.data.plasma_pressure[stag_idx]\n"
     "                    hs_energy_J = np.sum((1.5 * p_plasma[hot_spot_mask]) * vol[hot_spot_mask])"),

    # ----- Edit 7: ignition hot-spot pressure (line ~1492-1494) -----
    # Indentation: 16 spaces, line break is BEFORE the '+', not after
    ("Ignition hot-spot pressure: use plasma_pressure",
     "                # HS pressure: mass-averaged total pressure in Gbar\n"
     "                p_total = (self.data.ion_pressure[ign_idx]\n"
     "                           + self.data.rad_pressure[ign_idx]) * 1e-8  # Gbar",
     "                # HS pressure: mass-averaged plasma pressure (ion+electron) in Gbar\n"
     "                p_total = self.data.plasma_pressure[ign_idx] * 1e-8  # Gbar"),

    # ----- Edit 8: neutron-averaged pressure (line ~1734) -----
    ("Neutron-averaged pressure: use plasma_pressure",
     "                    # Include electron (rad) pressure for total pressure\n"
     "                    p_total = self.data.ion_pressure[t, z]\n"
     "                    if self.data.rad_pressure is not None:\n"
     "                        p_total += self.data.rad_pressure[t, z]\n"
     "                    pressure_sum += p_total * fusion_rate[t, z] * dt",
     "                    # Plasma pressure (ion + electron) -- conventional\n"
     "                    p_total = self.data.plasma_pressure[t, z]\n"
     "                    pressure_sum += p_total * fusion_rate[t, z] * dt"),
]

n_applied = 0
n_skipped = 0
for desc, old, new in edits:
    # Idempotency: if old anchor is gone but new text is present, already applied
    if old not in text and new in text:
        print(f"  [skip-already-applied] {desc}")
        n_skipped += 1
        continue
    count = text.count(old)
    if count == 0:
        print(f"  [ANCHOR-MISSING] {desc}")
        print(f"    First 100 chars of expected anchor:")
        print(f"      {old[:100]!r}")
        sys.exit(1)
    if count > 1:
        print(f"  [ANCHOR-NOT-UNIQUE] {desc} (found {count} matches)")
        sys.exit(1)
    text = text.replace(old, new, 1)
    n_applied += 1
    print(f"  [ok] {desc}")

ia.write_text(text)

# Verify
verify = ia.read_text()
n_plasma_total = verify.count("self.data.plasma_pressure")
print()
print(f"Total plasma_pressure references in icf_analysis.py: {n_plasma_total}")
print(f"Edits applied this run: {n_applied}")
print(f"Edits already applied:  {n_skipped}")

# All four edits should produce these specific patterns
assert "self.data.plasma_pressure[stag_idx] * 1e-8" in verify, "stag_idx pressure anchor missing"
assert "self.data.plasma_pressure[ign_idx] * 1e-8" in verify, "ign_idx pressure anchor missing"
assert "self.data.plasma_pressure[t, z]" in verify, "neutron-averaged pressure anchor missing"
assert "p_plasma = self.data.plasma_pressure[stag_idx]" in verify, "internal energy anchor missing"

# Make sure no '+ self.data.rad_pressure' remnants in the analysis pipeline
# (other than the legacy data_builder.py setup, which we're not touching)
remaining = verify.count("+ self.data.rad_pressure")
remaining += verify.count(".ion_pressure[stag_idx] +")
remaining += verify.count(".ion_pressure[ign_idx]")
print(f"Remaining ion+rad combinations in icf_analysis.py: {remaining}")
print()

if n_applied > 0:
    print("Next:")
    print("  git diff helios_postprocess/icf_analysis.py | head -80")
    print("  git add helios_postprocess/icf_analysis.py")
    print("  git commit -m 'Plasma pressure follow-up: hot-spot/ignition/neutron-avg call sites'")
    print("  git push")
else:
    print("Nothing to commit -- patch was already complete.")
