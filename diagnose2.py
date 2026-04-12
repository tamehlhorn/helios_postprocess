"""
Second diagnostic -- run on Mac Studio.
Usage: python3 ~/diagnose2.py
"""
import os

BASE = os.path.expanduser("~/helios_postprocessor/helios_postprocess")

def show_def(filename, method_name, context=60):
    path = os.path.join(BASE, filename)
    with open(path) as f:
        lines = f.readlines()
    search = f"def {method_name}"
    for i, line in enumerate(lines):
        if search in line:
            end = min(len(lines), i + context)
            print(f"\n{'='*60}")
            print(f"  {filename} :: {method_name}  (line {i+1})")
            print(f"{'='*60}")
            for j, l in enumerate(lines[i:end], start=i+1):
                print(f"  {j:4d}  {l}", end="")
            break
    else:
        print(f"\n[NOT FOUND] def {method_name} in {filename}")

def show_around(filename, search_term, context=30):
    path = os.path.join(BASE, filename)
    with open(path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if search_term in line:
            start = max(0, i - 3)
            end = min(len(lines), i + context)
            print(f"\n{'='*60}")
            print(f"  {filename}  ('{search_term}' at line {i+1})")
            print(f"{'='*60}")
            for j, l in enumerate(lines[start:end], start=start+1):
                marker = ">>>" if (j-1 == i) else "   "
                print(f"{marker} {j:4d}  {l}", end="")
            return
    print(f"\n[NOT FOUND] '{search_term}' in {filename}")

# 1. The actual hot-spot properties method
show_def("icf_analysis.py", "_compute_hot_spot_properties", context=70)

# 2. The comparison table / metric loop in compare_with_published
show_def("burn_averaged_metrics.py", "compare_with_published", context=80)

# 3. Check what ablation_front_radius looks like on data
show_around("icf_analysis.py", "ablation_front_radius", context=10)

print("\nDone.")
