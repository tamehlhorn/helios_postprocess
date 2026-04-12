"""
Diagnostic 3 — find return structure in burn_averaged_metrics.py
Usage: python3 ~/diagnose3.py
"""
import os

path = os.path.expanduser(
    "~/helios_postprocessor/helios_postprocess/burn_averaged_metrics.py"
)

def show_def(method_name, context=80):
    with open(path) as f:
        lines = f.readlines()
    search = f"def {method_name}"
    for i, line in enumerate(lines):
        if search in line:
            end = min(len(lines), i + context)
            print(f"\n{'='*60}")
            print(f"  {method_name}  (line {i+1})")
            print(f"{'='*60}")
            for j, l in enumerate(lines[i:end], start=i+1):
                print(f"  {j:4d}  {l}", end="")
            return
    print(f"[NOT FOUND] def {method_name}")

show_def("extract_histories_from_run_data", context=120)
show_def("calculate_burn_averaged_metrics", context=80)
