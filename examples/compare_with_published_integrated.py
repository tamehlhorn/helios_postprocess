"""
Example: Compare Helios Simulation with Published Target Design

This script demonstrates the complete workflow for comparing your Helios
simulation results with published ICF target design data.

Author: Prof T
Date: November 2025
"""

import sys
from pathlib import Path

# Add package to path if needed
# sys.path.insert(0, str(Path(__file__).parent))

from helios_postprocess import (
    HeliosRun,
    extract_hot_spot_histories,
    calculate_burn_averaged_metrics,
    compare_with_published
)
import matplotlib.pyplot as plt
import numpy as np


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
    
    # From your published table (with alpha-heating values)
    published_metrics = {
        'T_hs': (46.7, 4.8),      # ⟨T_hs⟩ in keV (value, uncertainty)
        'P_hs': (2720, 212),      # ⟨P_hs⟩ in Gbar
        'rhoR_cf': (1.60, 0.46),  # ⟨ρR_cf⟩ in g/cm²
        'CR_max': (20.1, 9.6),    # CR_max
        'yield': (256, 0.6),      # Yield in MJ
        'gain': (65, 0)           # Prompt fusion gain
    }
    
    laser_energy_MJ = 4.0  # From table
    
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
    # STEP 2: Load Helios Simulation
    # ========================================================================
    print("STEP 2: Loading Helios Simulation Data")
    print("-" * 80)
    
    try:
        run = HeliosRun(exo_file, verbose=True)
        print(f"✓ Loaded: {exo_file}")
        print(f"  Time range: {run.times[0]*1e9:.2f} - {run.times[-1]*1e9:.2f} ns")
        print(f"  Number of time steps: {run.n_times}")
        print()
    except FileNotFoundError:
        print(f"✗ ERROR: Could not find file '{exo_file}'")
        print("  Please provide the correct path to your ExodusII file.")
        print()
        print("Usage:")
        print(f"  python {Path(__file__).name} <path_to_exo_file>")
        return
    
    # ========================================================================
    # STEP 3: Extract Hot Spot Time Histories
    # ========================================================================
    print("STEP 3: Extracting Hot Spot Time Histories")
    print("-" * 80)
    
    # Hot spot threshold: 1 keV = 1000 eV
    T_threshold_eV = 1000.0
    
    print(f"Hot spot threshold: {T_threshold_eV} eV (1 keV)")
    print("Extracting:")
    print("  - Hot spot temperature (neutron-averaged at each time)")
    print("  - Hot spot pressure (neutron-averaged)")
    print("  - Hot spot density (neutron-averaged)")
    print("  - Cold fuel areal density")
    print("  - Hot spot radius, volume, mass")
    print()
    
    try:
        histories = extract_hot_spot_histories(
            run,
            T_threshold=T_threshold_eV,
            time_indices=None  # Use all time steps
        )
        
        print()
        print("✓ Extraction successful")
        print(f"  Extracted {len(histories['time_ns'])} time steps")
        print(f"  Time range: {histories['time_ns'][0]:.2f} - {histories['time_ns'][-1]:.2f} ns")
        print(f"  Peak temperature: {np.max(histories['temperature_keV']):.1f} keV")
        print(f"  Peak pressure: {np.max(histories['pressure_Gbar']):.0f} Gbar")
        print(f"  Min radius: {np.min(histories['radius_um'][histories['radius_um']>0]):.1f} μm")
        print()
        
    except Exception as e:
        print(f"✗ ERROR during extraction: {e}")
        print("  This may be due to:")
        print("    - Incompatible variable names in ExodusII file")
        print("    - Missing required fields")
        print("    - Geometry assumptions")
        print("  Please check the extract_hot_spot_histories() function")
        print("  and adjust for your specific ExodusII format.")
        return
    
    # ========================================================================
    # STEP 4: Calculate Burn-Averaged Metrics
    # ========================================================================
    print("STEP 4: Calculating Burn-Averaged Metrics")
    print("-" * 80)
    
    print("Computing burn-weighted temporal averages:")
    print("  <Q> = ∫ Q(t) · Ṙ(t) dt / ∫ Ṙ(t) dt")
    print("  where Ṙ(t) ∝ ρ² · <σv>(T)")
    print()
    
    sim_metrics = calculate_burn_averaged_metrics(
        histories,
        ion_fraction=0.5  # 50-50 DT mix
    )
    
    print("Simulation Results:")
    print(f"  ⟨T_hs⟩       = {sim_metrics['T_burn_avg']:.1f} keV")
    print(f"  ⟨P_hs⟩       = {sim_metrics['P_burn_avg']:.0f} Gbar")
    print(f"  ⟨ρR_cf⟩      = {sim_metrics['rhoR_burn_avg']:.2f} g/cm²")
    print(f"  CR_max       = {sim_metrics['CR_max']:.1f}")
    print(f"  Yield        = {sim_metrics['yield_MJ']:.1f} MJ")
    print(f"  Fusion Gain  = {sim_metrics['yield_MJ']/laser_energy_MJ:.1f}")
    print()
    
    # ========================================================================
    # STEP 5: Generate Comparison Table
    # ========================================================================
    print("STEP 5: Comparison with Published Data")
    print("-" * 80)
    
    comparison = compare_with_published(
        sim_metrics,
        published_metrics,
        laser_energy_MJ=laser_energy_MJ
    )
    print(comparison)
    print()
    
    # ========================================================================
    # STEP 6: Check Ignition Criteria
    # ========================================================================
    print("STEP 6: Ignition Criteria Assessment")
    print("-" * 80)
    
    criteria = {
        'Temperature > 5 keV': sim_metrics['T_burn_avg'] > 5.0,
        'Pressure > 100 Gbar': sim_metrics['P_burn_avg'] > 100.0,
        'Areal density > 0.3 g/cm²': sim_metrics['rhoR_burn_avg'] > 0.3,
        'Convergence ratio > 15': sim_metrics['CR_max'] > 15.0
    }
    
    for criterion, satisfied in criteria.items():
        status = "✓ PASS" if satisfied else "✗ FAIL"
        print(f"  {criterion:<35} {status}")
    
    print()
    if all(criteria.values()):
        print("  ✓✓✓ ALL IGNITION CRITERIA SATISFIED ✓✓✓")
    else:
        print("  ⚠ SOME IGNITION CRITERIA NOT MET")
    print()
    
    # ========================================================================
    # STEP 7: Generate Diagnostic Plots
    # ========================================================================
    print("STEP 7: Generating Diagnostic Plots")
    print("-" * 80)
    
    plot_comparison(histories, sim_metrics, published_metrics, laser_energy_MJ)
    
    print("✓ Plots saved:")
    print("  - burn_comparison_diagnostics.png")
    print()
    
    # ========================================================================
    # STEP 8: Summary
    # ========================================================================
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print()
    print("Results Summary:")
    print(f"  Your simulation achieved:")
    print(f"    - {sim_metrics['T_burn_avg']:.1f} keV burn-averaged temperature")
    print(f"    - {sim_metrics['P_burn_avg']:.0f} Gbar burn-averaged pressure")
    print(f"    - {sim_metrics['yield_MJ']:.1f} MJ fusion yield")
    print(f"    - {sim_metrics['yield_MJ']/laser_energy_MJ:.1f}× energy gain")
    print()
    
    # Calculate percent differences
    delta_T = 100 * (sim_metrics['T_burn_avg'] - published_metrics['T_hs'][0]) / published_metrics['T_hs'][0]
    delta_P = 100 * (sim_metrics['P_burn_avg'] - published_metrics['P_hs'][0]) / published_metrics['P_hs'][0]
    delta_Y = 100 * (sim_metrics['yield_MJ'] - published_metrics['yield'][0]) / published_metrics['yield'][0]
    
    print(f"  Comparison with published design:")
    print(f"    - Temperature: {delta_T:+.1f}%")
    print(f"    - Pressure: {delta_P:+.1f}%")
    print(f"    - Yield: {delta_Y:+.1f}%")
    print()
    print("=" * 80)


def plot_comparison(histories, sim_metrics, published_metrics, laser_energy_MJ):
    """Generate diagnostic comparison plots."""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Helios Simulation vs Published Design', 
                 fontsize=16, fontweight='bold')
    
    time = histories['time_ns']
    burn_fraction = sim_metrics['burn_fraction']
    
    # Temperature
    ax = axes[0, 0]
    ax.plot(time, histories['temperature_keV'], 'b-', linewidth=2, label='Simulation')
    ax.axhline(sim_metrics['T_burn_avg'], color='b', linestyle='--',
               label=f'Burn-avg: {sim_metrics["T_burn_avg"]:.1f} keV')
    pub_T = published_metrics['T_hs'][0]
    ax.axhline(pub_T, color='r', linestyle='--',
               label=f'Published: {pub_T:.1f} keV')
    ax.fill_between(time, 0, np.max(histories['temperature_keV']) * 1.1,
                     burn_fraction / burn_fraction.max(), alpha=0.2, color='orange')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Temperature (keV)')
    ax.set_title('Hot Spot Temperature')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Pressure
    ax = axes[0, 1]
    ax.plot(time, histories['pressure_Gbar'], 'b-', linewidth=2)
    ax.axhline(sim_metrics['P_burn_avg'], color='b', linestyle='--',
               label=f'Burn-avg: {sim_metrics["P_burn_avg"]:.0f} Gbar')
    pub_P = published_metrics['P_hs'][0]
    ax.axhline(pub_P, color='r', linestyle='--',
               label=f'Published: {pub_P:.0f} Gbar')
    ax.fill_between(time, 0, np.max(histories['pressure_Gbar']) * 1.1,
                     burn_fraction / burn_fraction.max(), alpha=0.2, color='orange')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Pressure (Gbar)')
    ax.set_title('Hot Spot Pressure')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Areal Density
    ax = axes[0, 2]
    ax.plot(time, histories['areal_density_gcm2'], 'b-', linewidth=2)
    ax.axhline(sim_metrics['rhoR_burn_avg'], color='b', linestyle='--',
               label=f'Burn-avg: {sim_metrics["rhoR_burn_avg"]:.2f} g/cm²')
    pub_rhoR = published_metrics['rhoR_cf'][0]
    ax.axhline(pub_rhoR, color='r', linestyle='--',
               label=f'Published: {pub_rhoR:.2f} g/cm²')
    ax.fill_between(time, 0, np.max(histories['areal_density_gcm2']) * 1.1,
                     burn_fraction / burn_fraction.max(), alpha=0.2, color='orange')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('ρR (g/cm²)')
    ax.set_title('Cold Fuel Areal Density')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Convergence Ratio
    ax = axes[1, 0]
    CR_history = histories['initial_radius_um'] / histories['radius_um']
    ax.plot(time, CR_history, 'b-', linewidth=2)
    ax.axhline(sim_metrics['CR_max'], color='b', linestyle='--',
               label=f'Max CR: {sim_metrics["CR_max"]:.1f}')
    pub_CR = published_metrics['CR_max'][0]
    ax.axhline(pub_CR, color='r', linestyle='--',
               label=f'Published: {pub_CR:.1f}')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Convergence Ratio')
    ax.set_title('Convergence Ratio History')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Burn Rate Profile
    ax = axes[1, 1]
    ax.plot(time, burn_fraction / burn_fraction.max(), 'orange', linewidth=2)
    ax.fill_between(time, 0, burn_fraction / burn_fraction.max(),
                     alpha=0.3, color='orange')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Normalized Burn Rate')
    ax.set_title('Fusion Burn Profile')
    ax.grid(True, alpha=0.3)
    
    # Comparison Bar Chart
    ax = axes[1, 2]
    metrics_names = ['T_hs\n(keV)', 'P_hs\n(Gbar)', 'ρR_cf\n(g/cm²)', 'CR_max', 'Yield\n(MJ)']
    sim_values = [
        sim_metrics['T_burn_avg'],
        sim_metrics['P_burn_avg'],
        sim_metrics['rhoR_burn_avg'],
        sim_metrics['CR_max'],
        sim_metrics['yield_MJ']
    ]
    pub_values = [
        published_metrics['T_hs'][0],
        published_metrics['P_hs'][0],
        published_metrics['rhoR_cf'][0],
        published_metrics['CR_max'][0],
        published_metrics['yield'][0]
    ]
    
    # Normalize for comparison
    sim_norm = [s/p for s, p in zip(sim_values, pub_values)]
    
    x = np.arange(len(metrics_names))
    width = 0.35
    
    ax.bar(x, sim_norm, width, label='Simulation/Published', color='steelblue')
    ax.axhline(1.0, color='r', linestyle='--', linewidth=2, label='Published')
    ax.set_ylabel('Ratio to Published')
    ax.set_title('Normalized Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(metrics_names, fontsize=9)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('burn_comparison_diagnostics.png', dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    # Get ExodusII file from command line or use default
    if len(sys.argv) > 1:
        exo_file = sys.argv[1]
    else:
        exo_file = "Vulcan_HDD_NB_111125_AI_2.exo"
    
    main(exo_file)
