#!/usr/bin/env python3
"""
Patch B2 (corrected anchor): replace the broken single-region stub in
extract_histories_from_run_data() with one that uses the keys
calculate_burn_averaged_metrics() actually expects.

Run from /Users/mehlhorn/Codes/helios_postprocess:
    python3 apply_patch_B2_v2.py
"""
from pathlib import Path
import sys

bm_path = Path('helios_postprocess/burn_averaged_metrics.py')
if not bm_path.is_file():
    print(f"ERROR: {bm_path} not found. Run from repo root.")
    sys.exit(1)

bm = bm_path.read_text()

# Anchor matching exactly what's currently on disk (single-quoted, plain np)
old = '''    if ri is None or ri.shape[1] < 2:
        logger.info("Hot-spot history extraction: single-region target, no hot spot -- skipping")
        return {
            'time_ns':       np.asarray(data.time),
            'T_hs_keV':      np.full(n_times, np.nan),
            'P_hs_Gbar':     np.full(n_times, np.nan),
            'rho_hs_gcc':    np.full(n_times, np.nan),
            'rhoR_cf':       np.full(n_times, np.nan),
            'r_hs_cm':       np.full(n_times, np.nan),
            'v_shell_kms':   np.full(n_times, np.nan),
            'single_region': True,
        }'''

new = '''    if ri is None or ri.shape[1] < 2:
        # Single-region target (e.g. CH-only slab/sphere flux-limiter test).
        # Return a dict with the keys calculate_burn_averaged_metrics() expects,
        # filled with NaN/zero so its Simpson-weighted integrals evaluate to
        # 0 / NaN cleanly without crashing.
        logger.info("Hot-spot history extraction: single-region target -- returning NaN stub")
        nan_arr  = np.full(n_times, np.nan)
        zero_arr = np.zeros(n_times)
        return {
            'time_ns':            np.asarray(data.time),
            'temperature_keV':    zero_arr,    # zero -> burn_rate = 0 -> clean zero integral
            'pressure_Gbar':      nan_arr,
            'density_gcc':        zero_arr,
            'areal_density_gcm2': nan_arr,
            'radius_um':          nan_arr,
            'CR_max':             0.0,
            'energy_output_MJ':   0.0,
            'laser_energy_MJ':    getattr(data, 'laser_energy_delivered_MJ', 0.0),
            'target_gain':        0.0,
            'stagnation_time_ns': 0.0,
            'single_region':      True,
        }'''

# Idempotent re-application
if "'temperature_keV':" in bm and "'T_hs_keV'" not in bm:
    print("Patch B2 v2 already applied -- nothing to do.")
    sys.exit(0)

count = bm.count(old)
if count == 0:
    print("ERROR: anchor not found. Current stub may differ from expected.")
    print("Inspect with:")
    print("  grep -A 12 \"ri is None or ri.shape\\[1\\] < 2\" helios_postprocess/burn_averaged_metrics.py")
    sys.exit(1)

assert count == 1, f"Anchor not unique (found {count} matches)"

bm = bm.replace(old, new)
bm_path.write_text(bm)

verify = Path('helios_postprocess/burn_averaged_metrics.py').read_text()
assert "'temperature_keV':" in verify, "Verify failed: 'temperature_keV' not on disk"
assert "'T_hs_keV'" not in verify, "Verify failed: old key 'T_hs_keV' still present"

print("Patch B2 v2 applied and verified.")
print()
print("Next:")
print("  git diff helios_postprocess/burn_averaged_metrics.py | head -50")
print("  git add helios_postprocess/burn_averaged_metrics.py")
print("  git commit -m 'Patch B2 v2: single-region stub uses correct dict keys'")
print("  git push")
