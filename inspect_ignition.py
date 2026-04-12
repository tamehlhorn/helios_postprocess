#!/usr/bin/env python3
"""
Diagnostic: inspect ignition criterion and burn propagation storage
Run on MacBook: python ~/Codes/helios_postprocessor/inspect_ignition.py
Run on Mac Studio: python3 ~/helios_postprocessor/inspect_ignition.py
"""
import subprocess, pathlib, sys

repo = pathlib.Path.home() / ("Codes/helios_postprocessor" if "mehlhorn" == __import__("os").environ.get("USER","") else "helios_postprocessor")
pkg  = repo / "helios_postprocess"

def grep(pattern, filepath, context=3):
    result = subprocess.run(
        ["grep", "-n", "-A", str(context), "-B", str(context), pattern, str(filepath)],
        capture_output=True, text=True
    )
    return result.stdout

print("=" * 70)
print("1. IGNITION CRITERION in icf_analysis.py")
print("   Looking for rhoR threshold and fractional criterion")
print("=" * 70)
for pattern in ["ignition", "0\\.3", "rhoR_hs", "50", "fraction.*rhoR", "rhoR.*fraction"]:
    hit = grep(pattern, pkg / "icf_analysis.py")
    if hit:
        print(f"\n--- pattern: '{pattern}' ---")
        print(hit[:800])

print("\n" + "=" * 70)
print("2. BURN PROPAGATION storage -- what scalars are set on data?")
print("   Looking for ignition_time, hs_pressure_ignition, hs_radius_ignition")
print("=" * 70)
for pattern in ["ignition_time", "hs_pressure_ign", "hs_radius_ign", "P_hs_ign", "radius_ign"]:
    hit = grep(pattern, pkg / "icf_analysis.py")
    if hit:
        print(f"\n--- pattern: '{pattern}' ---")
        print(hit[:600])

print("\n" + "=" * 70)
print("3. ICFRunData dataclass -- ignition-related fields")
print("=" * 70)
for pattern in ["ignition", "hs_radius_ign", "P_hs_ign"]:
    hit = grep(pattern, pkg / "data_builder.py")
    if hit:
        print(f"\n--- data_builder.py '{pattern}' ---")
        print(hit[:600])

print("\n" + "=" * 70)
print("4. compare_with_published -- current metric keys handled")
print("=" * 70)
hit = grep("compare_with_published\|P_hs_ignition\|hs_radius_ignition\|METRIC_LABELS\|metric_labels", pkg / "burn_averaged_metrics.py")
print(hit[:2000])

print("\nDone. Paste output back to Claude.")
