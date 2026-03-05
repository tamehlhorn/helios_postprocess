"""
Example: Compare Helios Simulation with Published Target Design

Demonstrates the complete workflow for comparing Helios simulation results
with published ICF target design data using the v3.0 pipeline.

Usage:
    python compare_with_published_integrated.py <path_to_exo_file>

Author: Prof T
Date: November 2025 (original), March 2026 (updated for v3.0 pipeline)
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from helios_postprocess import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.icf_analysis import ICFAnalyzer
from helios_postprocess.burn_averaged_metrics import (
    extract_histories_from_run_data,
    calculate_burn_averaged_metrics,
    compare_with_published,
)


def main(exo_file: str = "Vulcan_HDD_NB_111125_AI_2.exo"):
    """
    Main workflow for burn-averaged analysis and comparison.

    Parameters
    ----------
    exo_file : str
        Path to Helios ExodusII output file
    """
    print("=" * 80)
    print("HELIOS SIMULATION COMPARISON WITH PUBLISHED TARGET DESIGN")
    print("=" * 80)
    print()

    # ========================================================================
    # STEP 1: Define Published Target Design Values
    # ========================================================================
    print("STEP 1: Loading Published Target Design")
    print("-" * 80)

    # Published table values (with alpha-heating)
    published_metrics = {
        'T_hs':    (46.7,  4.8),     # ⟨T_hs⟩ in keV (value, uncertainty)
        'P_hs':    (2720,  212),     # ⟨P_hs⟩ in Gbar
        'rhoR_cf': (1.60,  0.46),   # ⟨ρR_cf⟩ in g/cm²
        'CR_max':  (20.1,  9.6),    # CR_max
        'yield':   (256,   0.6),    # Yield in MJ
        'gain':    (65,    0),      # Prompt fusion gain
    }

    laser_energy_MJ = 4.0

    print("Published Design Specifications:")
    print(f"  ⟨T_hs⟩       = {published_metrics['T_hs'][0]:.1f} ± {published_metrics['T_hs'][1]:.1f} keV")
    print(f"  ⟨P_hs⟩       = {published_metrics['P_hs'][0]:.0f} ± {published_metrics['P_hs'][1]:.0f} Gbar")
    print(f"  ⟨ρR_cf⟩      = {published_metrics['rhoR_cf'][0]:.2f} ± {published_metrics['rhoR_cf'][1]:.2f} g/cm²")
    print(f"  CR_max       = {published_metrics['CR_max'][0]:.1f} ± {published_metrics['CR_max'][1]:.1f}")
    print(f"  Yield        = {published_metrics['yield'][0]:.1f} ± {published_metrics['yield'][1]:.1f} MJ")
    print(f"  Fusion Gain  = {published_metrics['gain'][0]:.0f}")
    print(f"  Laser Energy = {laser_energy_MJ:.1f} MJ")
    print()

    # ========================================================================
    # STEP 2: Load and Analyze with Pipeline
    # ========================================================================
    print("STEP 2: Loading and Analyzing Simulation")
    print("-" * 80)

    try:
        run = HeliosRun(exo_file, verbose=True)
        print(f"  Loaded: {exo_file}")
        print(f"  Time range: {run.times[0]*1e9:.2f} - {run.times[-1]*1e9:.2f} ns")
        print(f"  Time steps: {run.n_times}")
    except FileNotFoundError:
        print(f"  ERROR: Could not find file '{exo_file}'")
        print(f"  Usage: python {Path(__file__).name} <path_to_exo_file>")
        return

    data = build_run_data(run, time_unit='s')
    run.close()

    analyzer = ICFAnalyzer(data)
    analyzer.analyze_drive_phase()
    analyzer.analyze_stagnation_phase()
    analyzer.analyze_burn_phase()
    analyzer.compute_performance_metrics()

    print("  Pipeline analysis complete")
    print()

    # ========================================================================
    # STEP 3: Extract Histories and Compute Burn-Averaged Metrics
    # ========================================================================
    print("STEP 3: Computing Burn-Averaged Metrics")
    print("-" * 80)
    print("  <Q> = ∫ Q(t) · Ṙ(t) dt / ∫ Ṙ(t) dt")
    print("  where Ṙ(t) ∝ ρ² · <σv>(T)")
    print()

    histories = extract_histories_from_run_data(data)
    sim_metrics = calculate_burn_averaged_metrics(histories, ion_fraction=0.5)

    print("Simulation Results:")
    print(f"  ⟨T_hs⟩       = {sim_metrics['T_burn_avg']:.1f} keV")
    print(f"  ⟨P_hs⟩       = {sim_metrics['P_burn_avg']:.0f} Gbar")
    print(f"  ⟨ρR_cf⟩      = {sim_metrics['rhoR_burn_avg']:.2f} g/cm²")
    print(f"  CR_max       = {sim_metrics['CR_max']:.1f}")
    print(f"  Yield        = {sim_metrics['yield_MJ']:.1f} MJ")
    print(f"  Fusion Gain  = {sim_metrics['yield_MJ']/laser_energy_MJ:.1f}")
    print()

    # ========================================================================
    # STEP 4: Comparison Table
    # ========================================================================
    print("STEP 4: Comparison with Published Data")
    print("-" * 80)

    comparison = compare_with_published(
        sim_metrics, published_metrics, laser_energy_MJ=laser_energy_MJ
    )
    print(comparison)
    print()

    # ========================================================================
    # STEP 5: Ignition Criteria
    # ========================================================================
    print("STEP 5: Ignition Criteria Assessment")
    print("-" * 80)

    criteria = {
        'Temperature > 5 keV':        sim_metrics['T_burn_avg'] > 5.0,
        'Pressure > 100 Gbar':        sim_metrics['P_burn_avg'] > 100.0,
        'Areal density > 0.3 g/cm²':  sim_metrics['rhoR_burn_avg'] > 0.3,
        'Convergence ratio > 15':      sim_metrics['CR_max'] > 15.0,
    }

    for criterion, satisfied in criteria.items():
        status = "PASS" if satisfied else "FAIL"
        print(f"  {criterion:<35} {status}")

    print()
    if all(criteria.values()):
        print("  ALL IGNITION CRITERIA SATISFIED")
    else:
        print("  SOME IGNITION CRITERIA NOT MET")
    print()

    # ========================================================================
    # STEP 6: Diagnostic Plots
    # ========================================================================
    print("STEP 6: Generating Diagnostic Plots")
    print("-" * 80)

    plot_comparison(histories, sim_metrics, published_metrics, laser_energy_MJ)

    print("  Saved: burn_comparison_diagnostics.png")
    print()

    # ========================================================================
    # Summary
    # ========================================================================
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print()
    print(f"  Burn-averaged T  = {sim_metrics['T_burn_avg']:.1f} keV")
    print(f"  Burn-averaged P  = {sim_metrics['P_burn_avg']:.0f} Gbar")
    print(f"  Fusion yield     = {sim_metrics['yield_MJ']:.1f} MJ")
    print(f"  Energy gain      = {sim_metrics['yield_MJ']/laser_energy_MJ:.1f}x")

    delta_T = 100 * (sim_metrics['T_burn_avg'] - published_metrics['T_hs'][0]) / published_metrics['T_hs'][0]
    delta_P = 100 * (sim_metrics['P_burn_avg'] - published_metrics['P_hs'][0]) / published_metrics['P_hs'][0]
    delta_Y = 100 * (sim_metrics['yield_MJ'] - published_metrics['yield'][0]) / published_metrics['yield'][0]

    print(f"  vs published:  T {delta_T:+.1f}%,  P {delta_P:+.1f}%,  Yield {delta_Y:+.1f}%")
    print("=" * 80)


def plot_comparison(histories, sim_metrics, published_metrics, laser_energy_MJ):
    """Generate diagnostic comparison plots."""

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Helios Simulation vs Published Design',
                 fontsize=16, fontweight='bold')

    time = histories['time_ns']
    burn_frac = sim_metrics['burn_fraction']
    bf_max = burn_frac.max() if burn_frac.max() > 0 else 1.0
    bf_norm = burn_frac / bf_max

    # ── Temperature ──
    ax = axes[0, 0]
    ax.plot(time, histories['temperature_keV'], 'b-', lw=2, label='Simulation')
    ax.axhline(sim_metrics['T_burn_avg'], color='b', ls='--',
               label=f'Burn-avg: {sim_metrics["T_burn_avg"]:.1f} keV')
    pub_T = published_metrics['T_hs'][0]
    ax.axhline(pub_T, color='r', ls='--', label=f'Published: {pub_T:.1f} keV')
    ax.fill_between(time, 0, np.max(histories['temperature_keV']) * 1.1,
                    bf_norm, alpha=0.2, color='orange')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Temperature (keV)')
    ax.set_title('Hot Spot Temperature')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Pressure ──
    ax = axes[0, 1]
    ax.plot(time, histories['pressure_Gbar'], 'b-', lw=2)
    ax.axhline(sim_metrics['P_burn_avg'], color='b', ls='--',
               label=f'Burn-avg: {sim_metrics["P_burn_avg"]:.0f} Gbar')
    pub_P = published_metrics['P_hs'][0]
    ax.axhline(pub_P, color='r', ls='--', label=f'Published: {pub_P:.0f} Gbar')
    ax.fill_between(time, 0, np.max(histories['pressure_Gbar']) * 1.1,
                    bf_norm, alpha=0.2, color='orange')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Pressure (Gbar)')
    ax.set_title('Hot Spot Pressure')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Areal Density ──
    ax = axes[0, 2]
    ax.plot(time, histories['areal_density_gcm2'], 'b-', lw=2)
    ax.axhline(sim_metrics['rhoR_burn_avg'], color='b', ls='--',
               label=f'Burn-avg: {sim_metrics["rhoR_burn_avg"]:.2f} g/cm²')
    pub_rhoR = published_metrics['rhoR_cf'][0]
    ax.axhline(pub_rhoR, color='r', ls='--',
               label=f'Published: {pub_rhoR:.2f} g/cm²')
    rhoR_max = np.max(histories['areal_density_gcm2'])
    if rhoR_max > 0:
        ax.fill_between(time, 0, rhoR_max * 1.1, bf_norm, alpha=0.2, color='orange')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('ρR (g/cm²)')
    ax.set_title('Cold Fuel Areal Density')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Convergence Ratio ──
    ax = axes[1, 0]
    r0 = histories['initial_radius_um']
    r = histories['radius_um']
    valid = r > 0
    CR = np.where(valid, r0 / r, 0)
    ax.plot(time, CR, 'b-', lw=2)
    ax.axhline(sim_metrics['CR_max'], color='b', ls='--',
               label=f'Max CR: {sim_metrics["CR_max"]:.1f}')
    pub_CR = published_metrics['CR_max'][0]
    ax.axhline(pub_CR, color='r', ls='--', label=f'Published: {pub_CR:.1f}')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Convergence Ratio')
    ax.set_title('Convergence Ratio History')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Burn Rate Profile ──
    ax = axes[1, 1]
    ax.plot(time, bf_norm, 'orange', lw=2)
    ax.fill_between(time, 0, bf_norm, alpha=0.3, color='orange')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Normalized Burn Rate')
    ax.set_title('Fusion Burn Profile')
    ax.grid(True, alpha=0.3)

    # ── Normalized Comparison Bar Chart ──
    ax = axes[1, 2]
    names = ['T_hs\n(keV)', 'P_hs\n(Gbar)', 'ρR_cf\n(g/cm²)', 'CR_max', 'Yield\n(MJ)']
    sim_vals = [
        sim_metrics['T_burn_avg'],
        sim_metrics['P_burn_avg'],
        sim_metrics['rhoR_burn_avg'],
        sim_metrics['CR_max'],
        sim_metrics['yield_MJ'],
    ]
    pub_vals = [
        published_metrics['T_hs'][0],
        published_metrics['P_hs'][0],
        published_metrics['rhoR_cf'][0],
        published_metrics['CR_max'][0],
        published_metrics['yield'][0],
    ]
    ratios = [s / p if p > 0 else 0 for s, p in zip(sim_vals, pub_vals)]

    x = np.arange(len(names))
    ax.bar(x, ratios, 0.35, label='Simulation / Published', color='steelblue')
    ax.axhline(1.0, color='r', ls='--', lw=2, label='Published')
    ax.set_ylabel('Ratio to Published')
    ax.set_title('Normalized Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=9)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('burn_comparison_diagnostics.png', dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    exo_file = sys.argv[1] if len(sys.argv) > 1 else "Vulcan_HDD_NB_111125_AI_2.exo"
    main(exo_file)
