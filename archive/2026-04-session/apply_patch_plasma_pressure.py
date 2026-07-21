#!/usr/bin/env python3
"""
Plasma-pressure convention patch
=================================
Switches the ICF postprocessor to the conventional definition:

    plasma_pressure = ion_pressure + electron_pressure

(This excludes radiation pressure. Radiation pressure typically contributes
< 10% in foot-pulse-era physics but ~30-40% in hot coronas, which was
inflating ablation-pressure and shock-breakout diagnostics. Per Tom's
direction, we also use plasma pressure -- not plasma+radiation -- for the
hot-spot pressure diagnostic.)

This patch:

1. Adds `plasma_pressure` (n_times, n_zones, J/cm^3) as a derived quantity on
   ICFRunData, computed from ion_pressure + elec_pressure once at load time.
2. Leaves `rad_pressure` semantically unchanged for backward compat (it
   currently contains "everything that isn't ion" for legacy reasons -- see
   data_builder header comment). Going forward, code should use
   `plasma_pressure` for ion+electron, and `data.rad_pressure_true` (also
   added here) for the actual radiation pressure component if needed.
3. Updates all six call sites that previously did `ion + rad` to use
   `plasma_pressure` instead.

Files modified:
    helios_postprocess/data_builder.py            (add plasma_pressure derivation)
    helios_postprocess/icf_analysis.py            (5 call sites)
    helios_postprocess/burn_averaged_metrics.py   (1 call site)

Run from /Users/mehlhorn/Codes/helios_postprocess:
    python3 apply_patch_plasma_pressure.py

Effect on key WfCDT_01b numbers:
    Before                           After (plasma-only)
    ----------------------------------------------------
    Hot-spot pressure   103.4 Gbar   roughly half (the 3-component
                                     rad_pressure had electron pressure
                                     baked in, so this number is mainly
                                     re-expressing what was always
                                     dominated by electrons)
    Adiabat (peak-vel)  0.33         lower (P ratio includes less)
    Base adiabat        0.21         lower
    Shock breakout P    in Mbar      lower

These are conventional-units changes -- the underlying simulation hydro is
unchanged. Helios still has the correct dynamics.
"""
from pathlib import Path
import sys

repo = Path(".")
db = repo / "helios_postprocess" / "data_builder.py"
ia = repo / "helios_postprocess" / "icf_analysis.py"
bm = repo / "helios_postprocess" / "burn_averaged_metrics.py"

for f in [db, ia, bm]:
    if not f.is_file():
        print(f"ERROR: {f} not found. Run from repo root.")
        sys.exit(1)


def patch_file(path: Path, edits: list[tuple[str, str, str]]) -> int:
    """Apply a list of (description, old, new) edits to a file. Return count."""
    text = path.read_text()
    n_applied = 0
    for desc, old, new in edits:
        if new in text and old not in text:
            print(f"  [skip-already-applied] {desc}")
            continue
        count = text.count(old)
        if count == 0:
            print(f"  [ANCHOR-MISSING] {desc}")
            print(f"    First 80 chars of expected anchor: {old[:80]!r}")
            sys.exit(1)
        if count > 1:
            print(f"  [ANCHOR-NOT-UNIQUE] {desc} (found {count} matches)")
            sys.exit(1)
        text = text.replace(old, new, 1)
        n_applied += 1
        print(f"  [ok] {desc}")
    path.write_text(text)
    return n_applied


# ============================================================================
# data_builder.py: add plasma_pressure as a derived quantity
# ============================================================================
print(f"\n--- {db} ---")

db_edits = [
    # 1. Add the field to ICFRunData.__init__
    ("Add plasma_pressure field to ICFRunData",
     "        self.elec_pressure: Optional[np.ndarray] = None       # (n_times, n_zones) J/cm³",
     "        self.elec_pressure: Optional[np.ndarray] = None       # (n_times, n_zones) J/cm³\n"
     "        self.plasma_pressure: Optional[np.ndarray] = None     # (n_times, n_zones) J/cm³  ion + electron\n"
     "        self.rad_pressure_true: Optional[np.ndarray] = None   # (n_times, n_zones) J/cm³  radiation only"),

    # 2. After the 3-component handling, derive plasma_pressure and rad_pressure_true
    ("Derive plasma_pressure and rad_pressure_true after rad_pressure assembly",
     "    elif actual_rad is not None:\n"
     "        data.rad_pressure = actual_rad\n"
     "        if verbose:\n"
     "            logger.info(\"  ✓ rad_pressure              ← rad_pressure only\")\n"
     "    else:\n"
     "        if data.ion_pressure is not None:\n"
     "            data.rad_pressure = np.zeros_like(data.ion_pressure)\n"
     "            logger.warning(\"  ⚠ No elec/rad pressure found — rad_pressure set to zeros\")\n"
     "        else:\n"
     "            data.rad_pressure = None",
     "    elif actual_rad is not None:\n"
     "        data.rad_pressure = actual_rad\n"
     "        if verbose:\n"
     "            logger.info(\"  ✓ rad_pressure              ← rad_pressure only\")\n"
     "    else:\n"
     "        if data.ion_pressure is not None:\n"
     "            data.rad_pressure = np.zeros_like(data.ion_pressure)\n"
     "            logger.warning(\"  ⚠ No elec/rad pressure found — rad_pressure set to zeros\")\n"
     "        else:\n"
     "            data.rad_pressure = None\n"
     "\n"
     "    # ------------------------------------------------------------------\n"
     "    # Derived: plasma_pressure = ion + electron  (CONVENTIONAL DEFINITION)\n"
     "    # This is what should be used for shock breakout, ablation pressure,\n"
     "    # adiabat, hot-spot pressure -- all the standard ICF diagnostics.\n"
     "    # The legacy `rad_pressure` field still carries elec+rad combined for\n"
     "    # backward compatibility but should be considered deprecated.\n"
     "    # `rad_pressure_true` carries only the actual radiation pressure if\n"
     "    # available, otherwise zeros.\n"
     "    # ------------------------------------------------------------------\n"
     "    if data.ion_pressure is not None and data.elec_pressure is not None:\n"
     "        data.plasma_pressure = data.ion_pressure + data.elec_pressure\n"
     "        if verbose:\n"
     "            logger.info(\"  ✓ plasma_pressure           ← ion_pressure + elec_pressure (conventional)\")\n"
     "    elif data.ion_pressure is not None:\n"
     "        # Fallback: assume Helios treats ion = electron locally so use ion alone\n"
     "        data.plasma_pressure = data.ion_pressure\n"
     "        logger.warning(\"  ⚠ No electron pressure found -- plasma_pressure = ion_pressure only\")\n"
     "\n"
     "    if actual_rad is not None:\n"
     "        data.rad_pressure_true = actual_rad\n"
     "    elif data.ion_pressure is not None:\n"
     "        data.rad_pressure_true = np.zeros_like(data.ion_pressure)"),
]
patch_file(db, db_edits)


# ============================================================================
# icf_analysis.py: 5 call sites
# ============================================================================
print(f"\n--- {ia} ---")

ia_edits = [
    # 1. _track_shock_fronts (line ~432)
    ("_track_shock_fronts: use plasma_pressure",
     "            pressure = (self.data.ion_pressure + self.data.rad_pressure) * 1e-8  # J/cm³ → Gbar",
     "            # Plasma pressure (ion + electron) -- conventional for shock dynamics\n"
     "            pressure = self.data.plasma_pressure * 1e-8  # J/cm³ → Gbar"),

    # 2. _compute_adiabat
    ("_compute_adiabat: use plasma_pressure",
     "            # -- Total pressure in Mbar over the shell mask --\n"
     "            p_tot = self.data.ion_pressure[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]\n"
     "            if self.data.rad_pressure is not None:\n"
     "                p_tot = p_tot + self.data.rad_pressure[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]\n"
     "            p_Mbar = p_tot * 1e-5                        # J/cm3 -> Mbar",
     "            # -- Plasma pressure (ion + electron) in Mbar over the shell mask --\n"
     "            # Conventional adiabat uses plasma pressure, not plasma+radiation.\n"
     "            p_tot = self.data.plasma_pressure[eval_idx, fuel_lo:fuel_hi][shell_mask_loc]\n"
     "            p_Mbar = p_tot * 1e-5                        # J/cm3 -> Mbar"),

    # 3. _compute_adiabat_at_breakout
    ("_compute_adiabat_at_breakout: use plasma_pressure",
     "            p_tot = self.data.ion_pressure[idx_b, fuel_lo:fuel_hi][shell_mask_loc]\n"
     "            if self.data.rad_pressure is not None:\n"
     "                p_tot = p_tot + self.data.rad_pressure[idx_b, fuel_lo:fuel_hi][shell_mask_loc]\n"
     "            p_Mbar = p_tot * 1e-5",
     "            # Plasma pressure (ion + electron) -- conventional for adiabat\n"
     "            p_tot = self.data.plasma_pressure[idx_b, fuel_lo:fuel_hi][shell_mask_loc]\n"
     "            p_Mbar = p_tot * 1e-5"),

    # 4. _compute_shock_breakout
    ("_compute_shock_breakout: use plasma_pressure",
     "            total_P = self.data.ion_pressure + (\n"
     "                self.data.rad_pressure if self.data.rad_pressure is not None else 0.0\n"
     "            )",
     "            # Plasma pressure (ion + electron) -- conventional for shock/ablation\n"
     "            total_P = self.data.plasma_pressure"),

    # 5. _compute_hot_spot_properties: stagnation hot-spot pressure
    ("_compute_hot_spot_properties: use plasma_pressure",
     "                pressure_Gbar = (self.data.ion_pressure[stag_idx] +\n"
     "                           self.data.rad_pressure[stag_idx]) * 1e-8  # J/cm³ → Gbar",
     "                # Plasma pressure (ion + electron) at stagnation -- conventional definition\n"
     "                pressure_Gbar = self.data.plasma_pressure[stag_idx] * 1e-8  # J/cm³ → Gbar"),

    # 6. _compute_hot_spot_properties: internal energy summation
    ("_compute_hot_spot_properties: use plasma_pressure for internal energy",
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
     "                    # E_int = (3/2) P V  for ideal gas (plasma = ion + electron)\n"
     "                    vol = (4.0 / 3.0) * np.pi * (boundaries[1:]**3 - boundaries[:-1]**3)\n"
     "                    p_plasma = self.data.plasma_pressure[stag_idx]\n"
     "                    hs_energy_J = np.sum((1.5 * p_plasma[hot_spot_mask]) * vol[hot_spot_mask])"),

    # 7. ignition HS pressure
    ("Ignition hot-spot pressure: use plasma_pressure",
     "                # HS pressure: mass-averaged total pressure in Gbar\n"
     "                p_total = (self.data.ion_pressure[ign_idx]\n"
     "                           + self.data.rad_pressure[ign_idx]) * 1e-8  # Gbar",
     "                # HS pressure: mass-averaged plasma pressure (ion+electron) in Gbar\n"
     "                p_total = self.data.plasma_pressure[ign_idx] * 1e-8  # Gbar"),

    # 8. Neutron-averaged pressure
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
patch_file(ia, ia_edits)


# ============================================================================
# burn_averaged_metrics.py: 1 call site
# ============================================================================
print(f"\n--- {bm} ---")

bm_edits = [
    ("burn_averaged_metrics: use plasma_pressure",
     "    P_ion = data.ion_pressure                    # (n_times, n_zones) J/cm³\n"
     "    P_rad = data.rad_pressure                    # (n_times, n_zones) J/cm³\n"
     "    P_tot = P_ion + P_rad                        # total pressure",
     "    # Plasma pressure (ion + electron) -- conventional ICF definition\n"
     "    P_tot = data.plasma_pressure                 # (n_times, n_zones) J/cm³  ion + electron"),
]
patch_file(bm, bm_edits)


# ============================================================================
# Verify on disk
# ============================================================================
print("\n--- Verifying on-disk state ---")

assert "plasma_pressure: Optional[np.ndarray]" in db.read_text(), "data_builder field missing"
assert "data.plasma_pressure = data.ion_pressure + data.elec_pressure" in db.read_text(), "data_builder derivation missing"
print("  data_builder.py: plasma_pressure declared and derived")

ia_text = ia.read_text()
n_plasma = ia_text.count("self.data.plasma_pressure")
print(f"  icf_analysis.py: {n_plasma} plasma_pressure references")
assert n_plasma >= 7, f"Expected at least 7, got {n_plasma}"

bm_text = bm.read_text()
assert "data.plasma_pressure" in bm_text, "burn_averaged_metrics not patched"
print("  burn_averaged_metrics.py: plasma_pressure reference present")

print("\nPatch complete.")
print()
print("Next:")
print("  git diff --stat")
print("  git add helios_postprocess/")
print("  git commit -m 'Plasma pressure convention: use ion+electron throughout (excludes radiation)'")
print("  git push")
print()
print("CLAUDE.md update separately -- see message for the bullet to add.")
