#!/usr/bin/env python3
"""
Patch B2: fix the single-region stub dict in extract_histories_from_run_data()
to use the key names that calculate_burn_averaged_metrics() actually expects.

Patch B (previous) used invented key names (T_hs_keV, P_hs_Gbar, etc.).
The real expected keys, per calculate_burn_averaged_metrics() reading the
dict at line 386 of burn_averaged_metrics.py:
    - 'time_ns'
    - 'temperature_keV'
    - 'pressure_Gbar'
    - 'density_gcc'
    - 'areal_density_gcm2'
    - 'radius_um'

This patch replaces the stub with the correct keys. Returns a dict whose
shape matches what the downstream function expects, populated with NaN
arrays so the integrals return 0.0 / NaN cleanly without crashing.

Run from the repo root:
    python3 apply_patch_B2.py
"""
from pathlib import Path
import sys

bm_path = Path('helios_postprocess/burn_averaged_metrics.py')
if not bm_path.is_file():
    print(f"ERROR: {bm_path} not found. Run from repo root.")
    sys.exit(1)

bm = bm_path.read_text()

# Anchor on the Patch B stub we just added. Match as few lines as possible
# to keep the anchor robust to indentation variations.
old = (
    "    if ri is None or ri.shape[1] < 2:\n"
    "        import numpy as _np\n"
    "        logger.info(\"Hot-spot history extraction: single-region target -- skipping\")\n"
    "        return {\n"
    "            'time_ns':       _np.asarray(data.time),\n"
    "            'T_hs_keV':      _np.full(n_times, _np.nan),\n"
    "            'P_hs_Gbar':     _np.full(n_times, _np.nan),\n"
    "            'rho_hs_gcc':    _np.full(n_times, _np.nan),\n"
    "            'rhoR_cf':       _np.full(n_times, _np.nan),\n"
    "            'r_hs_cm':       _np.full(n_times, _np.nan),\n"
    "            'v_shell_kms':   _np.full(n_times, _np.nan),\n"
    "            'single_region': True,\n"
    "        }\n"
)

new = (
    "    if ri is None or ri.shape[1] < 2:\n"
    "        # Single-region target (e.g. CH-only slab/sphere flux-limiter test).\n"
    "        # Return a dict with the keys that calculate_burn_averaged_metrics()\n"
    "        # expects, filled with NaN so the Simpson-weighted integrals\n"
    "        # evaluate to 0 / NaN cleanly without crashing.\n"
    "        logger.info(\"Hot-spot history extraction: single-region target -- returning NaN stub\")\n"
    "        nan_arr = np.full(n_times, np.nan)\n"
    "        zero_arr = np.zeros(n_times)\n"
    "        return {\n"
    "            'time_ns':             np.asarray(data.time),\n"
    "            'temperature_keV':     zero_arr,   # zero -> burn_rate = 0 -> clean zero integral\n"
    "            'pressure_Gbar':       nan_arr,\n"
    "            'density_gcc':         zero_arr,\n"
    "            'areal_density_gcm2':  nan_arr,\n"
    "            'radius_um':           nan_arr,\n"
    "            'CR_max':              0.0,\n"
    "            'energy_output_MJ':    0.0,\n"
    "            'laser_energy_MJ':     getattr(data, 'laser_energy_delivered_MJ', 0.0),\n"
    "            'target_gain':         0.0,\n"
    "            'stagnation_time_ns':  0.0,\n"
    "            'single_region':       True,\n"
    "        }\n"
)

count = bm.count(old)
if count == 0:
    if "temperature_keV" in bm and "single_region" in bm:
        print("Patch B2 already applied -- nothing to do.")
        sys.exit(0)
    print("ERROR: anchor not found. File may have been modified since Patch B.")
    print("Expected to find the Patch B stub with 'T_hs_keV'. Inspect with:")
    print("  grep -A 15 'if ri is None or ri.shape\\[1\\] < 2:' helios_postprocess/burn_averaged_metrics.py")
    sys.exit(1)

assert count == 1, f"Anchor not unique (found {count} matches)"

bm = bm.replace(old, new)
bm_path.write_text(bm)

# Verify on disk
verify = Path('helios_postprocess/burn_averaged_metrics.py').read_text()
assert "'temperature_keV':" in verify, "Verify failed: correct keys not on disk"
assert "T_hs_keV" not in verify, "Verify failed: old stub still present"
print("Patch B2 applied and verified.")
print()
print("Next steps:")
print("  git diff helios_postprocess/burn_averaged_metrics.py | head -40")
print("  git add helios_postprocess/burn_averaged_metrics.py")
print("  git commit -m 'Patch B2: single-region stub dict uses correct keys (fix STEP 5 KeyError)'")
print("  git push")
