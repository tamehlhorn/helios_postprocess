"""
Patch 3 — Add stagnation_time_ns as comparison metric in burn_averaged_metrics.py

Two changes:
  1. Wire data.stag_time into the metrics dict returned by
     extract_histories_from_run_data() so it flows through to
     compare_with_published().
  2. Add stagnation_time_ns row to the implosion_rows table in
     compare_with_published().

Usage: python3 ~/patch3_stagtime.py
"""
import os, sys

path = os.path.expanduser(
    "~/helios_postprocessor/helios_postprocess/burn_averaged_metrics.py"
)

with open(path) as f:
    src = f.read()

# ── Change 1: wire stag_time into the returned histories/metrics dict ──────
# The stag_time is already used (line 269) to find stag_idx inside
# extract_histories_from_run_data().  We just need to return it so that
# calculate_burn_averaged_metrics / compare_with_published can use it.
# Search for the return block that builds the dict and add the field.

OLD_RETURN = '''    return {
        'time_ns':          time_ns,
        'burn_rate':        burn_rate,'''

NEW_RETURN = '''    return {
        'time_ns':          time_ns,
        'burn_rate':        burn_rate,
        'stagnation_time_ns': getattr(data, 'stag_time', 0.0),'''

change1_ok = OLD_RETURN in src
if not change1_ok:
    print("WARNING: Change 1 search string not found.")
    print("  Looked for the 'return {' block starting with 'time_ns'.")
    print("  Patch 3 change 1 NOT applied — manual wiring needed.")
else:
    src = src.replace(OLD_RETURN, NEW_RETURN, 1)
    print("Change 1 applied: stagnation_time_ns added to extract_histories return dict.")

# ── Change 2: pass stagnation_time_ns through calculate_burn_averaged_metrics ──
# Search for where sim_metrics is built in calculate_burn_averaged_metrics.
# Look for the yield_MJ line which is a reliable anchor.

OLD_METRICS = "    metrics['yield_MJ']"
if OLD_METRICS not in src:
    print("WARNING: Change 2 anchor 'metrics[\"yield_MJ\"]' not found — skipping.")
else:
    # Find the line and insert stagnation_time_ns near the top of the metrics dict
    OLD_STAG_PASSTHRU = "    metrics['yield_MJ']"
    NEW_STAG_PASSTHRU = (
        "    metrics['stagnation_time_ns'] = histories.get('stagnation_time_ns', 0.0)\n"
        "    metrics['yield_MJ']"
    )
    src = src.replace(OLD_STAG_PASSTHRU, NEW_STAG_PASSTHRU, 1)
    print("Change 2 applied: stagnation_time_ns passed through calculate_burn_averaged_metrics.")

# ── Change 3: add row to compare_with_published implosion_rows ─────────────
OLD_ROWS = '''        ('Peak velocity (km/s)',      sim_metrics.get('peak_velocity_kms', 0.0),
         'peak_velocity_kms',         '.1f'),'''

NEW_ROWS = '''        ('Stagnation time (ns)',       sim_metrics.get('stagnation_time_ns', 0.0),
         'stagnation_time_ns',        '.3f'),
        ('Peak velocity (km/s)',      sim_metrics.get('peak_velocity_kms', 0.0),
         'peak_velocity_kms',         '.1f'),'''

change3_ok = OLD_ROWS in src
if not change3_ok:
    print("WARNING: Change 3 search string not found.")
    print("  Looked for 'Peak velocity (km/s)' row in implosion_rows.")
    print("  Patch 3 change 3 NOT applied — manual addition needed.")
else:
    src = src.replace(OLD_ROWS, NEW_ROWS, 1)
    print("Change 3 applied: stagnation_time_ns row added to comparison table.")

with open(path, 'w') as f:
    f.write(src)

print("\npatch3_stagtime.py done.")
print("  Update your _published.json files: rename '_stagnation_time_ns' -> 'stagnation_time_ns'")
