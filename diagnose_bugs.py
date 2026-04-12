"""
Diagnostic script -- run on Mac Studio to extract the relevant code sections
before applying patches. Shows the exact lines we need to fix.

Usage:
    python3 ~/diagnose_bugs.py
"""
import re, os

BASE = os.path.expanduser("~/helios_postprocessor/helios_postprocess")

def show_section(filename, search_term, context=20):
    path = os.path.join(BASE, filename)
    with open(path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if search_term in line:
            start = max(0, i - 3)
            end = min(len(lines), i + context)
            print(f"\n{'='*60}")
            print(f"  {filename}  (found '{search_term}' at line {i+1})")
            print(f"{'='*60}")
            for j, l in enumerate(lines[start:end], start=start+1):
                marker = ">>>" if (j-1 == i) else "   "
                print(f"{marker} {j:4d}  {l}", end="")
            break
    else:
        print(f"\n[NOT FOUND] '{search_term}' in {filename}")

# 1. Hot-spot radius calculation
show_section("icf_analysis.py", "_compute_hot_spot_properties", context=60)

# 2. cr_inflight calculation
show_section("icf_analysis.py", "cr_inflight", context=20)

# 3. stagnation_time in burn_averaged_metrics
show_section("burn_averaged_metrics.py", "stagnation_time", context=15)

# 4. compare_with_published metric list
show_section("burn_averaged_metrics.py", "METRIC_DISPLAY", context=40)

# 5. stag_time field on data
show_section("burn_averaged_metrics.py", "stag_time", context=10)

print("\nDone.")
