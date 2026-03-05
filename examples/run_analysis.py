#!/usr/bin/env python3
"""
Helios ICF Post-Processing Runner
==================================

Runs the full analysis pipeline on a Helios simulation, producing:
  - PDF report with diagnostic plots
  - Text summary of key metrics
  - CSV time-history data
  - Burn-averaged metrics (printed to console)
  - Optional comparison with published target design data

Usage
-----
    python run_analysis.py <base_path>

where <base_path> is the path WITHOUT extension.  All files are derived:
    <base_path>.exo           — EXODUS simulation output (required)
    <base_path>.rhw           — RHW input file (optional)
    <base_path>_report.pdf    — output: diagnostic plots
    <base_path>_summary.txt   — output: text summary
    <base_path>_history.csv   — output: time histories
    <base_path>_published.json — input: published data for comparison (optional)

Examples
--------
    python run_analysis.py ~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9
    python run_analysis.py ~/Sims/Xcimer/Vulcan_HDD/VI_6/VI_6

Author: Prof T
Date: March 2026
"""

import sys
import json
import logging
from pathlib import Path

import numpy as np

from helios_postprocess import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.icf_analysis import ICFAnalyzer
from helios_postprocess.icf_plotting import ICFPlotter
from helios_postprocess.icf_output import ICFOutputGenerator
from helios_postprocess.burn_averaged_metrics import (
    extract_histories_from_run_data,
    calculate_burn_averaged_metrics,
    compare_with_published,
)


def main(base_path: str):
    """
    Run full post-processing pipeline.

    Parameters
    ----------
    base_path : str
        Path without extension, e.g. '~/Sims/Xcimer/Olson_PDD/Olson_PDD_9/Olson_PDD_9'
    """
    base = Path(base_path).expanduser().resolve()
    name = base.name

    # Derived paths
    exo_path       = base.with_suffix('.exo')
    rhw_path       = base.with_suffix('.rhw')
    report_path    = base.parent / f"{name}_report.pdf"
    summary_path   = base.parent / f"{name}"
    published_path = base.parent / f"{name}_published.json"

    print("=" * 80)
    print(f"  Helios Post-Processing: {name}")
    print("=" * 80)

    if not exo_path.exists():
        print(f"ERROR: EXODUS file not found: {exo_path}")
        sys.exit(1)

    print(f"  EXODUS:    {exo_path}")
    print(f"  RHW:       {rhw_path}  {'(found)' if rhw_path.exists() else '(not found — skipping)'}")
    print(f"  Published: {published_path}  {'(found)' if published_path.exists() else '(not found — skipping comparison)'}")
    print()

    logging.basicConfig(level=logging.INFO, format='%(name)s: %(message)s')

    # ── Load RHW if available ──
    rhw_config = None
    drive_temperature = None
    drive_time = None
    if rhw_path.exists():
        try:
            from helios_postprocess.rhw_parser import load_rhw_configuration
            rhw_config = load_rhw_configuration(rhw_path)
            drive_temperature = rhw_config.drive_temperature
            drive_time = rhw_config.drive_time
            print(f"  RHW config: {rhw_config.drive_type}, Burn {'ON' if rhw_config.burn_enabled else 'OFF'}")
            if drive_time is not None:
                print(f"  Drive temperature: {len(drive_time)} points")
            print()
        except Exception as e:
            print(f"  WARNING: Could not parse RHW file: {e}")
            print()

    # ── Step 1: Load EXODUS data ──
    print("-" * 80)
    print("STEP 1: Loading simulation data")
    print("-" * 80)

    run = HeliosRun(str(exo_path), verbose=True)
    data = build_run_data(
        run, time_unit='s',
        rhw_config=rhw_config,
        drive_temperature=drive_temperature,
        drive_time=drive_time,
    )
    run.close()
    print()

    # ── Step 2: Run analysis pipeline ──
    print("-" * 80)
    print("STEP 2: Running analysis pipeline")
    print("-" * 80)

    analyzer = ICFAnalyzer(data)
    analyzer.analyze_drive_phase()
    analyzer.analyze_stagnation_phase()
    analyzer.analyze_burn_phase()
    analyzer.compute_performance_metrics()
    print()

    # ── Step 3: Generate PDF report ──
    print("-" * 80)
    print("STEP 3: Generating PDF report")
    print("-" * 80)

    plotter = ICFPlotter(data, {})
    plotter.create_full_report(str(report_path))
    print(f"  Report: {report_path}")
    print()

    # ── Step 4: Write summary + CSV ──
    print("-" * 80)
    print("STEP 4: Writing summary and CSV output")
    print("-" * 80)

    output = ICFOutputGenerator(data)
    output.write_all(str(summary_path))
    print(f"  Summary: {summary_path}_summary.txt")
    print(f"  History: {summary_path}_history.csv")
    print()

    # ── Step 5: Burn-averaged and implosion metrics ──
    print("-" * 80)
    print("STEP 5: Burn-averaged and implosion metrics")
    print("-" * 80)

    histories = extract_histories_from_run_data(data)
    metrics = calculate_burn_averaged_metrics(histories)

    print()
    print("  Burn-averaged quantities:")
    print(f"    ⟨T_hs⟩          = {metrics['T_burn_avg']:.2f} keV")
    print(f"    ⟨P_hs⟩          = {metrics['P_burn_avg']:.1f} Gbar")
    print(f"    ⟨ρR_cf⟩         = {metrics['rhoR_burn_avg']:.4f} g/cm²")
    print(f"    CR_max          = {metrics['CR_max']:.1f}")
    print(f"    Yield           = {metrics['yield_MJ']:.3f} MJ")
    print(f"    Gain            = {metrics['target_gain']:.3f}")
    print()
    print("  Implosion metrics:")
    print(f"    Peak velocity   = {metrics['peak_velocity_kms']:.1f} km/s")
    print(f"    Adiabat         = {metrics['adiabat']:.2f}")
    if metrics['fraction_absorbed_pct'] > 0:
        print(f"    Frac absorbed   = {metrics['fraction_absorbed_pct']:.1f}%")
    if metrics['inflight_KE_kJ'] > 0:
        print(f"    In-flight KE    = {metrics['inflight_KE_kJ']:.1f} kJ")
    if metrics['hydro_efficiency_pct'] > 0:
        print(f"    Hydro eff       = {metrics['hydro_efficiency_pct']:.1f}%")
    if metrics['imploded_DT_mass_mg'] > 0:
        print(f"    Imploded DT     = {metrics['imploded_DT_mass_mg']:.2f} mg")
    print()

    # ── Step 6: Compare with published data (if available) ──
    if published_path.exists():
        print("-" * 80)
        print("STEP 6: Comparison with published data")
        print("-" * 80)

        try:
            with open(published_path, 'r') as f:
                pub_raw = json.load(f)

            # Extract laser energy and remove non-metric keys
            pub_laser = pub_raw.pop('laser_energy_MJ', None)
            published_metrics = {}
            for key, val in pub_raw.items():
                if key.startswith('_'):
                    continue   # skip comments
                if isinstance(val, (list, tuple)) and len(val) == 2:
                    published_metrics[key] = tuple(val)

            table = compare_with_published(
                metrics, published_metrics,
                laser_energy_MJ=pub_laser,
            )
            print(table)
            print()

        except Exception as e:
            print(f"  WARNING: Could not load published data: {e}")
            print()

    # ── Done ──
    print("=" * 80)
    print(f"  Analysis complete: {name}")
    print(f"  Outputs in: {base.parent}")
    print("=" * 80)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    main(sys.argv[1])
