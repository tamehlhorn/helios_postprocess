#!/usr/bin/env python3
"""
Patch B: guard extract_histories_from_run_data() for single-region targets.

Run from the repo root (where helios_postprocess/ lives):
    python3 apply_patch_B.py

Safe to re-run — asserts loudly if anchor is missing or already applied.
"""
from pathlib import Path
import sys

bm_path = Path('helios_postprocess/burn_averaged_metrics.py')
if not bm_path.is_file():
    print(f"ERROR: {bm_path} not found. Run from repo root.")
    sys.exit(1)

bm = bm_path.read_text()

old = (
    "    ri = data.region_interfaces_indices          # (n_times, n_interfaces)\n"
    "    T_ion = data.ion_temperature                 # (n_times, n_zones) eV\n"
    "    rho   = data.mass_density                    # (n_times, n_zones) g/cm\u00b3\n"
    "    P_ion = data.ion_pressure                    # (n_times, n_zones) J/cm\u00b3\n"
    "    P_rad = data.rad_pressure                    # (n_times, n_zones) J/cm\u00b3\n"
    "    P_tot = P_ion + P_rad                        # total pressure\n"
    "    zmass = data.zone_mass                       # (n_times, n_zones) g\n"
    "    zbnd  = data.zone_boundaries                 # (n_times, n_nodes) cm\n"
    "    vel   = data.velocity                        # (n_times, n_nodes) cm/s\n"
)

new = (
    "    ri = data.region_interfaces_indices          # (n_times, n_interfaces)\n"
    "\n"
    "    # Hot-spot averaging requires a defined gas/fuel interface (>=2 regions).\n"
    "    # Single-region targets (e.g. CH-only slab/sphere flux-limiter tests) have\n"
    "    # no hot spot to extract histories for -- return empty result.\n"
    "    if ri is None or ri.shape[1] < 2:\n"
    "        logger.info(\"Hot-spot history extraction: single-region target, no hot spot -- skipping\")\n"
    "        return {\n"
    "            'time_ns':       np.asarray(data.time),\n"
    "            'T_hs_keV':      np.full(n_times, np.nan),\n"
    "            'P_hs_Gbar':     np.full(n_times, np.nan),\n"
    "            'rho_hs_gcc':    np.full(n_times, np.nan),\n"
    "            'rhoR_cf':       np.full(n_times, np.nan),\n"
    "            'r_hs_cm':       np.full(n_times, np.nan),\n"
    "            'v_shell_kms':   np.full(n_times, np.nan),\n"
    "            'single_region': True,\n"
    "        }\n"
    "\n"
    "    T_ion = data.ion_temperature                 # (n_times, n_zones) eV\n"
    "    rho   = data.mass_density                    # (n_times, n_zones) g/cm\u00b3\n"
    "    P_ion = data.ion_pressure                    # (n_times, n_zones) J/cm\u00b3\n"
    "    P_rad = data.rad_pressure                    # (n_times, n_zones) J/cm\u00b3\n"
    "    P_tot = P_ion + P_rad                        # total pressure\n"
    "    zmass = data.zone_mass                       # (n_times, n_zones) g\n"
    "    zbnd  = data.zone_boundaries                 # (n_times, n_nodes) cm\n"
    "    vel   = data.velocity                        # (n_times, n_nodes) cm/s\n"
)

count = bm.count(old)
if count == 0:
    # Check for idempotent reapplication
    if "single_region': True" in bm:
        print("Patch B already applied -- nothing to do.")
        sys.exit(0)
    print("ERROR: anchor not found. File may have been modified.")
    sys.exit(1)

assert count == 1, f"Anchor not unique (found {count} matches)"

bm = bm.replace(old, new)
bm_path.write_text(bm)

# Re-read and verify on disk
verify = Path('helios_postprocess/burn_averaged_metrics.py').read_text()
assert "single_region': True" in verify, "Verify failed: edit not on disk"

print("Patch B applied and verified.")
print("Next steps:")
print("  git diff helios_postprocess/burn_averaged_metrics.py | head -40")
print("  git add helios_postprocess/burn_averaged_metrics.py")
print("  git commit -m 'Patch B: guard extract_histories_from_run_data for single-region targets'")
print("  git push")
