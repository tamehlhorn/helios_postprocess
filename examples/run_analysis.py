#!/usr/bin/env python3
"""
Helios ICF Post-Processing Runner
==================================

Runs the full analysis pipeline on a Helios simulation, producing:
  - PDF report with diagnostic plots
  - Comparison PDF with burn-averaged metrics and published-data table (if JSON provided)
  - Text summary of key metrics (includes comparison if JSON provided)
  - CSV time-history data

Usage
-----
    python run_analysis.py <base_path>

where <base_path> is the path WITHOUT extension.  All files are derived:
    <base_path>.exo               — EXODUS simulation output (required)
    <base_path>.rhw               — RHW input file (optional)
    <base_path>_report.pdf        — output: diagnostic plots
    <base_path>_comparison.pdf    — output: comparison table + burn metrics (if published JSON)
    <base_path>_summary.txt       — output: text summary (includes comparison)
    <base_path>_history.csv       — output: time histories
    <base_path>_published.json    — input: published data for comparison (optional)

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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
    """Run full post-processing pipeline."""
    base = Path(base_path).expanduser().resolve()
    name = base.name

    # Derived paths
    exo_path        = base.with_suffix('.exo')
    rhw_path        = base.with_suffix('.rhw')
    report_path     = base.parent / f"{name}_report.pdf"
    comparison_path = base.parent / f"{name}_comparison.pdf"
    summary_path    = base.parent / f"{name}"
    published_path  = base.parent / f"{name}_published.json"

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

    # ── Append burn-averaged metrics to summary text file ──
    summary_txt = f"{summary_path}_summary.txt"
    _append_metrics_to_summary(summary_txt, name, metrics)

    # ── Step 6: Compare with published data (if available) ──
    published_metrics = None
    pub_laser = None

    if published_path.exists():
        print("-" * 80)
        print("STEP 6: Comparison with published data")
        print("-" * 80)

        try:
            with open(published_path, 'r') as f:
                pub_raw = json.load(f)

            pub_laser = pub_raw.pop('laser_energy_MJ', None)
            published_metrics = {}
            for key, val in pub_raw.items():
                if key.startswith('_'):
                    continue
                if isinstance(val, (list, tuple)) and len(val) == 2:
                    published_metrics[key] = tuple(val)

            table = compare_with_published(
                metrics, published_metrics,
                laser_energy_MJ=pub_laser,
            )
            print(table)
            print()

            # Append comparison table to summary text file
            _append_comparison_to_summary(summary_txt, table)

            # Create comparison PDF
            _create_comparison_pdf(
                str(comparison_path), name,
                metrics, published_metrics, pub_laser, histories
            )
            print(f"  Comparison PDF: {comparison_path}")
            print()

        except Exception as e:
            print(f"  WARNING: Could not load published data: {e}")
            import traceback; traceback.print_exc()
            print()

    # ── Done ──
    print("=" * 80)
    print(f"  Analysis complete: {name}")
    print(f"  Outputs in: {base.parent}")
    outputs = [f"    {name}_report.pdf", f"    {name}_summary.txt", f"    {name}_history.csv"]
    if published_metrics is not None:
        outputs.append(f"    {name}_comparison.pdf")
    print("\n".join(outputs))
    print("=" * 80)


# ── Helper: append burn-averaged metrics to summary text ──

def _append_metrics_to_summary(summary_path: str, name: str, metrics: dict):
    """Append burn-averaged and implosion metrics to the summary text file."""
    with open(summary_path, 'a') as f:
        f.write("\n")
        f.write("=" * 72 + "\n")
        f.write("  BURN-AVERAGED & IMPLOSION METRICS\n")
        f.write("=" * 72 + "\n")
        f.write("\n")
        f.write("  Burn-averaged quantities:\n")
        f.write(f"    <T_hs>          = {metrics['T_burn_avg']:.2f} keV\n")
        f.write(f"    <P_hs>          = {metrics['P_burn_avg']:.1f} Gbar\n")
        f.write(f"    <rhoR_cf>       = {metrics['rhoR_burn_avg']:.4f} g/cm2\n")
        f.write(f"    CR_max          = {metrics['CR_max']:.1f}\n")
        f.write(f"    Yield           = {metrics['yield_MJ']:.3f} MJ\n")
        f.write(f"    Gain            = {metrics['target_gain']:.3f}\n")
        f.write("\n")
        f.write("  Implosion metrics:\n")
        f.write(f"    Peak velocity   = {metrics['peak_velocity_kms']:.1f} km/s\n")
        f.write(f"    Adiabat         = {metrics['adiabat']:.2f}\n")
        if metrics['fraction_absorbed_pct'] > 0:
            f.write(f"    Frac absorbed   = {metrics['fraction_absorbed_pct']:.1f}%\n")
        if metrics['inflight_KE_kJ'] > 0:
            f.write(f"    In-flight KE    = {metrics['inflight_KE_kJ']:.1f} kJ\n")
        if metrics['hydro_efficiency_pct'] > 0:
            f.write(f"    Hydro eff       = {metrics['hydro_efficiency_pct']:.1f}%\n")
        if metrics['imploded_DT_mass_mg'] > 0:
            f.write(f"    Imploded DT     = {metrics['imploded_DT_mass_mg']:.2f} mg\n")
        f.write("\n")


def _append_comparison_to_summary(summary_path: str, table: str):
    """Append comparison table to the summary text file."""
    with open(summary_path, 'a') as f:
        f.write("\n")
        f.write(table)
        f.write("\n")


# ── Helper: create comparison PDF ──

def _create_comparison_pdf(output_path: str, name: str,
                           metrics: dict, published_metrics: dict,
                           pub_laser: float, histories: dict):
    """Create a PDF with comparison table page and diagnostic bar chart."""

    with PdfPages(output_path) as pdf:
        # ── Page 1: Comparison table ──
        fig = plt.figure(figsize=(11, 8.5))

        title = f"{name} — Comparison with Published Target Design"
        fig.text(0.5, 0.95, title, fontsize=14, fontweight='bold',
                 ha='center', va='top', family='monospace')

        # Build table data
        laser_sim = metrics.get('laser_energy_MJ', 0.0)
        laser_pub = pub_laser if pub_laser else laser_sim
        sim_gain = metrics['yield_MJ'] / laser_pub if laser_pub > 0 else 0.0

        all_rows = []

        # Implosion section
        imp_defs = [
            ('Peak velocity (km/s)',  metrics.get('peak_velocity_kms', 0.0),  'peak_velocity_kms',     '.1f'),
            ('Adiabat',               metrics.get('adiabat', 0.0),            'adiabat',               '.2f'),
            ('Fraction absorbed (%)', metrics.get('fraction_absorbed_pct',0.0),'fraction_absorbed_pct', '.1f'),
            ('In-flight KE (kJ)',     metrics.get('inflight_KE_kJ', 0.0),     'inflight_KE_kJ',       '.1f'),
            ('Hydro efficiency (%)',  metrics.get('hydro_efficiency_pct',0.0), 'hydro_efficiency_pct', '.1f'),
            ('Imploded DT mass (mg)', metrics.get('imploded_DT_mass_mg',0.0), 'imploded_DT_mass_mg',  '.2f'),
        ]

        for label, sim_val, pub_key, fmt in imp_defs:
            pub_entry = published_metrics.get(pub_key)
            if pub_entry is None:
                continue
            pv, pu = _to_tuple(pub_entry)
            if pv <= 0 and sim_val <= 0:
                continue
            all_rows.append((label, sim_val, pv, pu, fmt))

        # Burn-averaged section
        burn_defs = [
            ('<T_hs> (keV)',       metrics['T_burn_avg'],    'T_hs',    '.1f'),
            ('<P_hs> (Gbar)',      metrics['P_burn_avg'],    'P_hs',    '.0f'),
            ('<rhoR_cf> (g/cm2)',  metrics['rhoR_burn_avg'], 'rhoR_cf', '.2f'),
            ('CR_max',             metrics['CR_max'],         'CR_max',  '.1f'),
            ('Yield (MJ)',         metrics['yield_MJ'],       'yield',   '.1f'),
            ('Fusion Gain',        sim_gain,                   'gain',    '.1f'),
        ]

        for label, sim_val, pub_key, fmt in burn_defs:
            pub_entry = published_metrics.get(pub_key)
            if pub_entry is None:
                continue
            pv, pu = _to_tuple(pub_entry)
            if pv <= 0:
                continue
            all_rows.append((label, sim_val, pv, pu, fmt))

        # Render table
        if all_rows:
            col_labels = ['Metric', 'Simulation', 'Published', 'Delta (%)']
            cell_text = []
            cell_colors = []
            for label, sv, pv, pu, fmt in all_rows:
                sv_str = format(sv, fmt) if sv > 0 else '—'
                pv_str = format(pv, fmt)
                if pu > 0:
                    pv_str += f' +/- {format(pu, fmt)}'
                if pv > 0 and sv > 0:
                    delta = 100 * (sv - pv) / pv
                    d_str = f'{delta:+.1f}%'
                    # Color: green if within 20%, yellow 20-50%, red >50%
                    if abs(delta) < 20:
                        color = '#d4edda'
                    elif abs(delta) < 50:
                        color = '#fff3cd'
                    else:
                        color = '#f8d7da'
                else:
                    d_str = '—'
                    color = '#ffffff'
                cell_text.append([label, sv_str, pv_str, d_str])
                cell_colors.append(['#ffffff', '#ffffff', '#ffffff', color])

            ax = fig.add_axes([0.08, 0.10, 0.84, 0.78])
            ax.axis('off')

            table = ax.table(
                cellText=cell_text,
                colLabels=col_labels,
                cellColours=cell_colors,
                colColours=['#cce5ff'] * 4,
                loc='upper center',
                cellLoc='center',
            )
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1.0, 1.6)

            # Bold header
            for j in range(4):
                table[0, j].set_text_props(fontweight='bold')

            # Laser energy footnote
            fig.text(0.5, 0.06,
                     f"Laser energy:  Simulation = {laser_sim:.3f} MJ,  Published = {laser_pub:.3f} MJ",
                     fontsize=9, ha='center', style='italic')

        pdf.savefig(fig)
        plt.close(fig)

        # ── Page 2: Burn history plots ──
        fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
        fig.suptitle(f'{name} — Burn-Averaged Diagnostics',
                     fontsize=14, fontweight='bold')

        time = histories['time_ns']
        bf = metrics['burn_fraction']
        bf_max = bf.max() if bf.max() > 0 else 1.0
        bf_norm = bf / bf_max

        # Temperature
        ax = axes[0, 0]
        ax.plot(time, histories['temperature_keV'], 'b-', lw=2)
        ax.axhline(metrics['T_burn_avg'], color='b', ls='--',
                    label=f'Burn-avg: {metrics["T_burn_avg"]:.1f} keV')
        pub_T = _pub_val(published_metrics, 'T_hs')
        if pub_T > 0:
            ax.axhline(pub_T, color='r', ls='--',
                        label=f'Published: {pub_T:.1f} keV')
        ax.fill_between(time, 0, np.max(histories['temperature_keV']) * 1.1,
                        bf_norm, alpha=0.15, color='orange')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Temperature (keV)')
        ax.set_title('Hot Spot Temperature')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Pressure
        ax = axes[0, 1]
        ax.plot(time, histories['pressure_Gbar'], 'b-', lw=2)
        ax.axhline(metrics['P_burn_avg'], color='b', ls='--',
                    label=f'Burn-avg: {metrics["P_burn_avg"]:.0f} Gbar')
        pub_P = _pub_val(published_metrics, 'P_hs')
        if pub_P > 0:
            ax.axhline(pub_P, color='r', ls='--',
                        label=f'Published: {pub_P:.0f} Gbar')
        ax.fill_between(time, 0, np.max(histories['pressure_Gbar']) * 1.1,
                        bf_norm, alpha=0.15, color='orange')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Pressure (Gbar)')
        ax.set_title('Hot Spot Pressure')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Areal density
        ax = axes[1, 0]
        ax.plot(time, histories['areal_density_gcm2'], 'b-', lw=2)
        ax.axhline(metrics['rhoR_burn_avg'], color='b', ls='--',
                    label=f'Burn-avg: {metrics["rhoR_burn_avg"]:.2f} g/cm2')
        pub_rhoR = _pub_val(published_metrics, 'rhoR_cf')
        if pub_rhoR > 0:
            ax.axhline(pub_rhoR, color='r', ls='--',
                        label=f'Published: {pub_rhoR:.2f} g/cm2')
        rhoR_max = np.max(histories['areal_density_gcm2'])
        if rhoR_max > 0:
            ax.fill_between(time, 0, rhoR_max * 1.1, bf_norm,
                            alpha=0.15, color='orange')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel(r'$\rho R$ (g/cm$^2$)')
        ax.set_title('Cold Fuel Areal Density')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Burn rate
        ax = axes[1, 1]
        ax.plot(time, bf_norm, 'orange', lw=2)
        ax.fill_between(time, 0, bf_norm, alpha=0.3, color='orange')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Normalized Burn Rate')
        ax.set_title('Fusion Burn Profile')
        ax.grid(True, alpha=0.3)

        plt.tight_layout(rect=[0, 0, 1, 0.94])
        pdf.savefig(fig)
        plt.close(fig)


def _to_tuple(val):
    if isinstance(val, (list, tuple)) and len(val) >= 2:
        return (float(val[0]), float(val[1]))
    elif isinstance(val, (int, float)):
        return (float(val), 0.0)
    return (0.0, 0.0)


def _pub_val(published_metrics, key):
    """Extract published value (first element), or 0."""
    if published_metrics is None:
        return 0.0
    entry = published_metrics.get(key)
    if entry is None:
        return 0.0
    v, _ = _to_tuple(entry)
    return v


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    main(sys.argv[1])
