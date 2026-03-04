"""
ICF Plotting Module
===================
Comprehensive plotting for ICF analysis:
- Time histories
- Contour plots
- Shock tracking
- Performance summaries
- Critical density position tracking

Adapted for HeliosRun / data_builder / icf_analysis pipeline.
Key changes from legacy version:
  - sklearn made optional (only needed for RANSAC shock fitting)
  - Pressures displayed in Gbar (not Mbar) per CLAUDE.md
  - Temperatures displayed in keV (not eV) per ICF convention
  - Hot spot metrics in Gbar / keV consistent with icf_analysis.py
  - Handles laser_energy_deposited as 1-D or 2-D

Author: Prof T / Xcimer ICF Analysis
Date: 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
from numpy.polynomial import Polynomial
import logging

try:
    from sklearn.linear_model import RANSACRegressor, LinearRegression
    _HAS_SKLEARN = True
except ImportError:
    _HAS_SKLEARN = False

logger = logging.getLogger(__name__)


class ICFPlotter:
    """Plotter for ICF simulation results."""
    
    def __init__(self, run_data, config):
        """
        Initialize plotter with run data and configuration.
        
        Parameters:
        -----------
        run_data : ICFRunData
            Container with simulation data and metrics
        config : dict
            Configuration parameters
        """
        self.data = run_data
        self.config = config
        
        # Set plotting style
        plt.style.use('seaborn-v0_8-darkgrid')
        self.default_figsize = (10, 7)
        
    def create_full_report(self, output_path: str):
        """
        Create comprehensive PDF report with all plots.
        
        Parameters:
        -----------
        output_path : str
            Path for output PDF file
        """
        logger.info(f"Creating PDF report: {output_path}")
        
        with PdfPages(output_path) as pdf:
            # Title page
            self._create_title_page(pdf)
            
            # Performance summary
            self._create_summary_page(pdf)
            
            # Time history plots
            self._plot_density_history(pdf)
            self._plot_temperature_history(pdf)
            self._plot_pressure_history(pdf)
            self._plot_velocity_history(pdf)
            self._plot_radius_history(pdf)
            self._plot_critical_density_position(pdf)
            self._plot_ablation_front(pdf)
            self._plot_combined_surface_tracking(pdf)
      
            # Fusion and burn
            # Fusion diagnostics (if available)
            if self.data.fusion_power is not None:
                self._plot_fusion_power(pdf)
                self._plot_neutron_production(pdf)
                self._plot_alpha_heating(pdf)
            
            # Laser/drive
            if self.data.laser_energy_deposited is not None:
                self._plot_laser_deposition(pdf)
            if self.data.laser_power_delivered is not None:
                self._plot_laser_power(pdf)
            
            # Drive temperature (from RHW file)
            if self.data.drive_temperature is not None:
                self._plot_drive_temperature(pdf)
            
            # Contour plots
            self._plot_density_contour(pdf)
            self._plot_pressure_contour(pdf)
            self._plot_velocity_contour(pdf)
            self._plot_temperature_contour(pdf)
            
            # Pressure gradient and shocks
            self._plot_pressure_gradient_contour(pdf)
            self._plot_shock_tracking(pdf)
            
            # Areal density evolution
            self._plot_areal_density_evolution(pdf)
            
        logger.info("PDF report complete")
        
    def _create_title_page(self, pdf):
        """Create title page with key analysis data."""
        fig, ax = plt.subplots(figsize=self.default_figsize)
        ax.axis('off')
        
        # Main title
        title_text = f"ICF Analysis Report\n\n{self.data.filename}"
        ax.text(0.5, 0.85, title_text, 
                ha='center', va='center', fontsize=24, weight='bold')
        
        subtitle = "Xcimer Inertial Confinement Fusion\nPost-Processing Analysis"
        ax.text(0.5, 0.72, subtitle,
                ha='center', va='center', fontsize=14)
        
        # Key analysis data box
        y_start = 0.60
        
        # Configuration data if available
        if self.data.rhw_config is not None:
            config_text = "SIMULATION CONFIGURATION\n"
            config_text += f"Drive Type: {self.data.rhw_config.drive_type}\n"
            config_text += f"Fusion Reactions: {'ENABLED' if self.data.rhw_config.burn_enabled else 'DISABLED'}\n"
            ax.text(0.5, y_start, config_text,
                   ha='center', va='top', fontsize=11,
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
            y_start -= 0.12
        
        # Key metrics
        metrics_text = "KEY RESULTS\n\n"
        metrics_text += f"Bang Time: {self.data.bang_time:.3f} ns\n"
        metrics_text += f"Stagnation Time: {self.data.stag_time:.3f} ns\n"
        metrics_text += f"Maximum Density: {self.data.max_density:.2f} g/cc\n"
        metrics_text += f"Compression Ratio: {self.data.comp_ratio:.1f}\n"
        metrics_text += f"Hot Spot Pressure: {self.data.hot_spot_pressure:.1f} Gbar\n"
        
        if self.data.energy_output > 0:
            metrics_text += f"Fusion Energy: {self.data.energy_output:.3f} MJ\n"
        if self.data.dt_neutron_yield > 0:
            metrics_text += f"DT Neutron Yield: {self.data.dt_neutron_yield:.3e}\n"
        if self.data.target_gain > 0:
            metrics_text += f"Target Gain: {self.data.target_gain:.2f}\n"
        
        ax.text(0.5, y_start, metrics_text,
               ha='center', va='top', fontsize=10, family='monospace',
               bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
        
        # Footer with timestamp
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ax.text(0.5, 0.05, f"Generated: {timestamp}",
               ha='center', va='bottom', fontsize=9, style='italic', color='gray')
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        
    def _create_summary_page(self, pdf):
        """Create performance summary page."""
        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis('off')
        
        # Title
        ax.text(0.5, 0.95, 'Performance Summary', 
                ha='center', fontsize=18, weight='bold', transform=ax.transAxes)
        
        # Build metrics list - add configuration section if available
        metrics = []
        
        # RHW Configuration (if available)
        if self.data.rhw_config is not None:
            config_metrics = [
                f"Drive Type: {self.data.rhw_config.drive_type}",
                f"Fusion Reactions: {'ENABLED' if self.data.rhw_config.burn_enabled else 'DISABLED'}",
            ]
            if self.data.drive_temperature is not None:
                t_ns = self.data.drive_time * 1e9
                config_metrics.append(
                    f"Drive Temperature: {np.min(self.data.drive_temperature):.1f} - "
                    f"{np.max(self.data.drive_temperature):.1f} eV"
                )
                config_metrics.append(f"Drive Profile: {len(self.data.drive_time)} points, "
                                     f"{t_ns[0]:.3f} - {t_ns[-1]:.3f} ns")
            metrics.append(('Configuration', config_metrics))
        
        # Append existing metrics
        metrics.extend([
            ('Timing Metrics', [
                f"Bang Time: {self.data.bang_time:.4f} ns",
                f"Stagnation Time: {self.data.stag_time:.4f} ns",
                f"Burn Width (FWHM): {self.data.burn_width:.4f} ns",
            ]),
            ('Compression and Density', [
                f"Maximum Density: {self.data.max_density:.4f} g/cc",
                f"Compression Ratio: {self.data.comp_ratio:.2f}",
                f"Core Radius (stagnation): {self.data.core_radius:.4f} cm",
            ]),
            ('Hot Spot Properties', [
                f"Hot Spot Radius (stagnation): {self.data.stagnation_hot_spot_radius:.4f} cm",
                f"Hot Spot Pressure: {self.data.hot_spot_pressure:.2f} Gbar",
                f"Hot Spot Areal Density: {self.data.hot_spot_areal_density:.4f} g/cm²",
                f"Hot Spot Internal Energy: {self.data.hot_spot_internal_energy:.2f} kJ",
            ]),
            ('Areal Densities', [
                f"Time-Averaged ρR: {self.data.time_ave_areal_density:.4f} g/cm²",
                f"Bang Time ρR: {self.data.bang_time_areal_density:.4f} g/cm²",
                f"Bang Time Fuel ρR: {self.data.bang_time_fuel_areal_density:.4f} g/cm²",
                f"Bang Time HDC ρR: {self.data.bang_time_HDC_areal_density:.4f} g/cm²",
                f"Bang Time Hot Spot ρR: {self.data.bang_time_hs_areal_density:.4f} g/cm²",
            ]),
            ('Neutron-Averaged Quantities', [
                f"Neutron-Avg Fuel ρR: {self.data.neutron_ave_fuel_areal_density:.4f} g/cm²",
                f"Neutron-Avg Temperature: {self.data.neutron_ave_ion_temperature:.2f} keV",
                f"Neutron-Avg Pressure: {self.data.neutron_ave_pressure:.2f} Gbar",
            ]),
            ('Energy and Performance', [
                f"Fusion Energy Output: {self.data.energy_output:.3f} MJ",
                f"DT Neutron Yield: {self.data.dt_neutron_yield:.3e}",
                f"Laser Energy Delivered: {self.data.laser_energy:.3f} MJ",
                f"Target Gain: {self.data.target_gain:.3f}",
                f"Maximum DT Temperature: {self.data.max_dt_temp:.2f} keV",
            ]),
            ('Mass Fractions', [
                f"Unablated Ablator Mass Fraction: {self.data.unablated_ablatar_mass:.4f}",
                f"Unablated Fuel Mass Fraction: {self.data.unablated_fuel_mass:.4f}",
                f"Stagnated Fuel Mass Fraction: {self.data.stagnated_fuel_mass:.4f}",
            ]),
            ('Additional Metrics', [
                f"Adiabat (mass-avg, ice): {self.data.adiabat_mass_averaged_ice:.4f}",
            ]),
        ])  # Close the metrics.extend() list
        
        y_pos = 0.88
        for section_title, section_metrics in metrics:
            # Section header
            ax.text(0.05, y_pos, section_title, fontsize=12, weight='bold',
                   transform=ax.transAxes)
            y_pos -= 0.03
            
            # Metrics
            for metric in section_metrics:
                ax.text(0.08, y_pos, metric, fontsize=10,
                       transform=ax.transAxes, family='monospace')
                y_pos -= 0.025
            
            y_pos -= 0.01  # Extra space between sections
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        
    def _plot_density_history(self, pdf):
        """Plot density evolution over time."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=self.default_figsize)
        
        # Mass-averaged density
        mass_avg_density = np.zeros(len(self.data.time))
        for t in range(len(self.data.time)):
            mass_avg_density[t] = np.average(
                self.data.mass_density[t],
                weights=self.data.zone_mass[t]
            )
        
        ax1.plot(self.data.time, mass_avg_density, 'b-', linewidth=2)
        ax1.axvline(self.data.bang_time, color='r', linestyle='--', 
                   label='Bang Time', alpha=0.7)
        ax1.axvline(self.data.stag_time, color='g', linestyle='--',
                   label='Stagnation', alpha=0.7)
        ax1.set_xlabel('Time (ns)')
        ax1.set_ylabel('Mass-Averaged Density (g/cc)')
        ax1.set_title('Density Evolution')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Peak density in each region
        n_zones = self.data.mass_density.shape[1]
        peak_density = np.max(self.data.mass_density, axis=1)
        
        ax2.plot(self.data.time, peak_density, 'r-', linewidth=2)
        ax2.axvline(self.data.bang_time, color='r', linestyle='--', alpha=0.7)
        ax2.axvline(self.data.stag_time, color='g', linestyle='--', alpha=0.7)
        ax2.set_xlabel('Time (ns)')
        ax2.set_ylabel('Peak Density (g/cc)')
        ax2.set_title('Peak Density vs Time')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        
    def plot_critical_density_position(data, output_path):
        """
        Plot critical density position evolution over time.
        
        Shows both the formula-based and algorithm-based methods for finding
        the critical density surface (laser absorption surface).
        
        Parameters
        ----------
        data : ICFRunData
            Container with simulation data including critical density positions
        output_path : Path
            Output directory for saving the plot
        """
        if (data.critical_density_radius_formula is None and 
            data.critical_density_radius_algorithm is None):
            return None
        
        fig, ax = plt.subplots(figsize=(10, 7))
        
        # Plot formula-based critical density position
        if data.critical_density_radius_formula is not None:
            valid_mask_formula = data.critical_density_radius_formula > 0
            if np.any(valid_mask_formula):
                ax.plot(data.time[valid_mask_formula], 
                       data.critical_density_radius_formula[valid_mask_formula],
                       'b-', linewidth=2, label='Formula Method (248 nm)',
                       alpha=0.8)
        
        # Plot algorithm-based critical density position
        if data.critical_density_radius_algorithm is not None:
            valid_mask_algo = data.critical_density_radius_algorithm > 0
            if np.any(valid_mask_algo):
                ax.plot(data.time[valid_mask_algo],
                       data.critical_density_radius_algorithm[valid_mask_algo],
                       'r--', linewidth=2, label='Algorithm Method (n_e crossing)',
                       alpha=0.8)
        
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Critical Density Surface Radius (cm)', fontsize=12)
        ax.set_title('Critical Density Position Evolution\n(Laser Absorption Surface)', 
                     fontsize=14, fontweight='bold')
        
        # Add critical density value to legend
        if hasattr(data, 'critical_density_value'):
            n_c_text = f'n_c = {data.critical_density_value:.2e} cm⁻³'
            ax.text(0.98, 0.98, n_c_text,
                   transform=ax.transAxes,
                   fontsize=10,
                   verticalalignment='top',
                   horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        save_path = output_path / 'critical_density_position.png'
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        return fig


    def plot_time_histories_with_critical_density(data, output_path):
        """
        Enhanced time history plot including critical density position.
        
        Creates a multi-panel plot showing:
        - Density evolution with critical density position overlay
        - Temperature evolution
        - Pressure evolution
        - Velocity evolution with critical density position overlay
        
        Parameters
        ----------
        data : ICFRunData
            Container with simulation data
        output_path : Path
            Output directory for saving plots
        """
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
        
        # Plot 1: Density with critical density position
        ax = axes[0]
        if data.mass_density is not None and data.zone_centers is not None:
            # Get mass-weighted average density
            mass_avg_density = np.zeros(len(data.time))
            max_density = np.zeros(len(data.time))
            
            for t in range(len(data.time)):
                mass = data.zone_mass[t]
                density = data.mass_density[t]
                mass_avg_density[t] = np.average(density, weights=mass)
                max_density[t] = np.max(density)
            
            ax.plot(data.time, mass_avg_density, 'b-', linewidth=2, label='Mass-weighted avg')
            ax.plot(data.time, max_density, 'r--', linewidth=2, label='Maximum')
            
            # Overlay critical density position on secondary y-axis
            if data.critical_density_radius_formula is not None:
                ax2 = ax.twinx()
                valid_mask = data.critical_density_radius_formula > 0
                if np.any(valid_mask):
                    ax2.plot(data.time[valid_mask],
                            data.critical_density_radius_formula[valid_mask],
                            'g-', linewidth=2, alpha=0.6, label='Critical density surface')
                    ax2.set_ylabel('Critical Density Radius (cm)', fontsize=12, color='g')
                    ax2.tick_params(axis='y', labelcolor='g')
                    ax2.legend(loc='upper right', fontsize=9)
            
            ax.set_xlabel('Time (ns)', fontsize=12)
            ax.set_ylabel('Density (g/cm³)', fontsize=12)
            ax.set_title('Density Evolution with Critical Density Surface', fontsize=13, fontweight='bold')
            ax.legend(loc='upper left', fontsize=10)
            ax.grid(True, alpha=0.3)
        
        # Plot 2: Temperature (eV → keV)
        ax = axes[1]
        if data.ion_temperature is not None:
            max_temp = np.max(data.ion_temperature, axis=1) / 1e3
            mass_avg_temp = np.zeros(len(data.time))
            
            for t in range(len(data.time)):
                mass = data.zone_mass[t]
                temp = data.ion_temperature[t]
                mass_avg_temp[t] = np.average(temp, weights=mass) / 1e3
            
            ax.plot(data.time, mass_avg_temp, 'b-', linewidth=2, label='Mass-weighted avg')
            ax.plot(data.time, max_temp, 'r--', linewidth=2, label='Maximum')
            
            ax.set_xlabel('Time (ns)', fontsize=12)
            ax.set_ylabel('Temperature (keV)', fontsize=12)
            ax.set_title('Temperature Evolution', fontsize=13, fontweight='bold')
            ax.legend(loc='best', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.set_yscale('log')
        
        # Plot 3: Pressure (J/cm³ → Gbar)
        ax = axes[2]
        if data.ion_pressure is not None:
            total_p = data.ion_pressure
            if data.rad_pressure is not None:
                total_p = data.ion_pressure + data.rad_pressure
            max_pressure = np.max(total_p, axis=1) * 1e-8
            mass_avg_pressure = np.zeros(len(data.time))
            
            for t in range(len(data.time)):
                mass = data.zone_mass[t]
                mass_avg_pressure[t] = np.average(total_p[t], weights=mass) * 1e-8
            
            ax.plot(data.time, mass_avg_pressure, 'b-', linewidth=2, label='Mass-weighted avg')
            ax.plot(data.time, max_pressure, 'r--', linewidth=2, label='Maximum')
            
            ax.set_xlabel('Time (ns)', fontsize=12)
            ax.set_ylabel('Pressure (Gbar)', fontsize=12)
            ax.set_title('Pressure Evolution', fontsize=13, fontweight='bold')
            ax.legend(loc='best', fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.set_yscale('log')
        
        # Plot 4: Velocity with critical density position
        ax = axes[3]
        if data.velocity is not None:
            min_velocity = np.min(data.velocity, axis=1) * (-1e-5)  # cm/s → km/s inward
            max_velocity = np.max(data.velocity, axis=1) * (-1e-5)
            
            ax.plot(data.time, min_velocity, 'b-', linewidth=2, label='Minimum (inward)')
            ax.plot(data.time, max_velocity, 'r--', linewidth=2, label='Maximum (outward)')
            ax.axhline(y=0, color='k', linestyle=':', alpha=0.5)
            
            # Overlay critical density position on secondary y-axis
            if data.critical_density_radius_algorithm is not None:
                ax2 = ax.twinx()
                valid_mask = data.critical_density_radius_algorithm > 0
                if np.any(valid_mask):
                    ax2.plot(data.time[valid_mask],
                            data.critical_density_radius_algorithm[valid_mask],
                            'm-', linewidth=2, alpha=0.6, label='Critical density surface (algo)')
                    ax2.set_ylabel('Critical Density Radius (cm)', fontsize=12, color='m')
                    ax2.tick_params(axis='y', labelcolor='m')
                    ax2.legend(loc='upper right', fontsize=9)
            
            ax.set_xlabel('Time (ns)', fontsize=12)
            ax.set_ylabel('Radial Velocity (km/s)', fontsize=12)
            ax.set_title('Velocity Evolution with Critical Density Surface', fontsize=13, fontweight='bold')
            ax.legend(loc='upper left', fontsize=10)
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        save_path = output_path / 'time_histories_with_critical_density.png'
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        return fig
        
    def _plot_temperature_history(self, pdf):
        """Plot temperature evolution."""
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        # Mass-averaged temperature (raw data in eV, display in keV)
        mass_avg_temp = np.zeros(len(self.data.time))
        peak_temp = np.zeros(len(self.data.time))
        
        for t in range(len(self.data.time)):
            mass_avg_temp[t] = np.average(
                self.data.ion_temperature[t],
                weights=self.data.zone_mass[t]
            )
            peak_temp[t] = np.max(self.data.ion_temperature[t])
        
        ax.plot(self.data.time, mass_avg_temp / 1e3, 'b-', linewidth=2, label='Mass-Averaged')
        ax.plot(self.data.time, peak_temp / 1e3, 'r-', linewidth=2, label='Peak')
        ax.axvline(self.data.bang_time, color='k', linestyle='--', 
                  label='Bang Time', alpha=0.7)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Ion Temperature (keV)')
        ax.set_title('Temperature Evolution')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_pressure_history(self, pdf):
        """Plot pressure evolution."""
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        # Total pressure (ion + radiation)
        total_pressure = self.data.ion_pressure + self.data.rad_pressure
        
        # Mass-averaged and peak
        mass_avg_pressure = np.zeros(len(self.data.time))
        peak_pressure = np.zeros(len(self.data.time))
        
        for t in range(len(self.data.time)):
            mass_avg_pressure[t] = np.average(
                total_pressure[t] * 1e-8,  # J/cm³ → Gbar
                weights=self.data.zone_mass[t]
            )
            peak_pressure[t] = np.max(total_pressure[t]) * 1e-8
        
        ax.plot(self.data.time, mass_avg_pressure, 'b-', linewidth=2, label='Mass-Averaged')
        ax.plot(self.data.time, peak_pressure, 'r-', linewidth=2, label='Peak')
        ax.axvline(self.data.bang_time, color='k', linestyle='--', alpha=0.7)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Pressure (Gbar)')
        ax.set_title('Pressure Evolution')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_velocity_history(self, pdf):
        """Plot velocity evolution."""
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        # Convert to km/s
        velocity_km_s = -self.data.velocity * 1e-7
        
        # Peak implosion velocity
        peak_vel = np.max(np.abs(velocity_km_s), axis=1)
        
        ax.plot(self.data.time, peak_vel, 'b-', linewidth=2)
        ax.axvline(self.data.bang_time, color='r', linestyle='--', 
                  label='Bang Time', alpha=0.7)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Peak Implosion Velocity (km/s)')
        ax.set_title('Implosion Velocity')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_critical_density_position(self, pdf):
        """Plot critical density position evolution over time."""
        if (self.data.critical_density_radius_formula is None and 
            self.data.critical_density_radius_algorithm is None):
            return
    
        fig, ax = plt.subplots(figsize=self.default_figsize)
    
        # Plot formula-based method
        if self.data.critical_density_radius_formula is not None:
            valid_mask = self.data.critical_density_radius_formula > 0
            if np.any(valid_mask):
                ax.plot(self.data.time[valid_mask], 
                    self.data.critical_density_radius_formula[valid_mask],
                   'b-', linewidth=2, label='Formula Method (n_e crossing)',
                   alpha=0.8)
    
        # Plot algorithm-based method
        if self.data.critical_density_radius_algorithm is not None:
            valid_mask = self.data.critical_density_radius_algorithm > 0
            if np.any(valid_mask):
                ax.plot(self.data.time[valid_mask],
                self.data.critical_density_radius_algorithm[valid_mask],
               'r--', linewidth=2, label='Algorithm Method (laser front)',
               alpha=0.8)
                
        # Mark key times
        if self.data.bang_time > 0:
            ax.axvline(self.data.bang_time, color='orange', linestyle='--', 
                  alpha=0.5, label='Bang Time')
        if self.data.stag_time > 0:
            ax.axvline(self.data.stag_time, color='green', linestyle='--',
                  alpha=0.5, label='Stagnation')
    
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Critical Density Surface Radius (cm)', fontsize=12)
        ax.set_title('Critical Density Position Evolution\n(Laser Absorption Surface)', 
                     fontsize=14, fontweight='bold')
    
        # Add critical density value annotation
        if self.data.critical_density_value > 0:
            n_c_text = f'$n_c$ = {self.data.critical_density_value:.2e} cm$^{{-3}}$'
            ax.text(0.98, 0.98, n_c_text,
               transform=ax.transAxes,
               fontsize=10,
               verticalalignment='top',
               horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
    
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    
    def _plot_ablation_front(self, pdf):
        """Plot ablation front position evolution over time."""
        if self.data.ablation_front_radius is None:
            logger.warning("Ablation front data not available for plotting")
            return
        
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        # Plot ablation front radius
        valid_mask = self.data.ablation_front_radius > 0
        if np.any(valid_mask):
            ax.plot(self.data.time[valid_mask], 
                   self.data.ablation_front_radius[valid_mask],
                   'purple', linewidth=2.5, label='Ablation Front',
                   alpha=0.8)
        
        # Mark key times
        if self.data.bang_time > 0:
            ax.axvline(self.data.bang_time, color='orange', linestyle='--', 
                      alpha=0.5, label='Bang Time')
        if self.data.stag_time > 0:
            ax.axvline(self.data.stag_time, color='green', linestyle='--',
                      alpha=0.5, label='Stagnation')
        
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Ablation Front Radius (cm)', fontsize=12)
        ax.set_title('Ablation Front Position vs Time\n(Minimum Density Scale Length)', 
                    fontsize=14, fontweight='bold')
        
        # Add annotation explaining the ablation front
        info_text = 'Ablation front: Location of steepest\ndensity gradient (min scale length)'
        ax.text(0.02, 0.98, info_text,
               transform=ax.transAxes,
               fontsize=9,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='lavender', alpha=0.7))
        
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Set reasonable axis limits if requested
        if self.config.get('ablation_front_xlim'):
            ax.set_xlim(self.config['ablation_front_xlim'])
        if self.config.get('ablation_front_ylim'):
            ax.set_ylim(self.config['ablation_front_ylim'])
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    
    def _plot_combined_surface_tracking(self, pdf):
        """
        Plot ablation front and critical density surfaces on the same axes.
        Shows how different physical surfaces evolve during the implosion.
        
        Note: Plot is limited to stagnation/bang time to avoid algorithm breakdown
        and to maintain appropriate scale for early-time details.
        """
        # Check if we have data to plot
        has_ablation = self.data.ablation_front_radius is not None
        has_crit_formula = self.data.critical_density_radius_formula is not None
        has_crit_algo = self.data.critical_density_radius_algorithm is not None
        
        if not (has_ablation or has_crit_formula or has_crit_algo):
            logger.warning("No surface tracking data available for combined plot")
            return
        
        # Determine time cutoff (stagnation or bang time)
        # This prevents algorithm breakdown and keeps scale appropriate
        time_cutoff = None
        cutoff_label = None
        
        if self.data.stag_time > 0:
            time_cutoff = self.data.stag_time
            cutoff_label = 'Stagnation'
        elif self.data.bang_time > 0:
            time_cutoff = self.data.bang_time
            cutoff_label = 'Bang Time'
        
        if time_cutoff is None:
            logger.warning("No stagnation or bang time available - plotting full time range")
            time_mask = np.ones(len(self.data.time), dtype=bool)
        else:
            time_mask = self.data.time <= time_cutoff
            logger.info(f"Combined surface plot limited to {cutoff_label} = {time_cutoff:.3e} ns")
        
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        # Track what we plotted for legend
        plotted_something = False
        
        # Plot ablation front
        if has_ablation:
            valid_mask = (self.data.ablation_front_radius > 0) & time_mask
            if np.any(valid_mask):
                ax.plot(self.data.time[valid_mask], 
                       self.data.ablation_front_radius[valid_mask],
                       color='purple', linewidth=2.5, 
                       label='Ablation Front (min scale length)',
                       alpha=0.9, linestyle='-')
                plotted_something = True
        
        # Plot critical density - formula method
        if has_crit_formula:
            valid_mask = (self.data.critical_density_radius_formula > 0) & time_mask
            if np.any(valid_mask):
                ax.plot(self.data.time[valid_mask], 
                       self.data.critical_density_radius_formula[valid_mask],
                       color='blue', linewidth=2, 
                       label='Critical Density (n$_e$ crossing)',
                       alpha=0.8, linestyle='-')
                plotted_something = True
        
        # Plot critical density - algorithm method
        if has_crit_algo:
            valid_mask = (self.data.critical_density_radius_algorithm > 0) & time_mask
            if np.any(valid_mask):
                ax.plot(self.data.time[valid_mask],
                       self.data.critical_density_radius_algorithm[valid_mask],
                       color='red', linewidth=2, 
                       label='Critical Density (laser front)',
                       alpha=0.8, linestyle='--')
                plotted_something = True
        
        if not plotted_something:
            logger.warning("No valid surface tracking data to plot")
            plt.close(fig)
            return
        
        # Mark key times (only if within plot range)
        if self.data.bang_time > 0 and self.data.bang_time <= (time_cutoff if time_cutoff else np.inf):
            ax.axvline(self.data.bang_time, color='orange', linestyle=':', 
                      linewidth=1.5, alpha=0.6, label='Bang Time')
        if self.data.stag_time > 0 and self.data.stag_time <= (time_cutoff if time_cutoff else np.inf):
            ax.axvline(self.data.stag_time, color='green', linestyle=':', 
                      linewidth=1.5, alpha=0.6, label='Stagnation')
        
        # Labels and title
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Surface Position (cm)', fontsize=12)
        
        # Title indicates time limit
        title = 'Combined Surface Tracking\nAblation Front & Critical Density Evolution'
        if time_cutoff is not None:
            title += f'\n(Limited to {cutoff_label})'
        ax.set_title(title, fontsize=14, fontweight='bold')
        
        # Add info box explaining the surfaces
        info_text = 'Multiple physical surfaces:\n'
        if has_ablation:
            info_text += '• Ablation front (steepest ρ gradient)\n'
        if has_crit_formula or has_crit_algo:
            info_text += '• Critical density (laser absorption)\n'
        if time_cutoff:
            info_text += f'\nPlot stops at {cutoff_label}\n'
            info_text += 'to show early-time details'
        else:
            info_text += '\nThese are distinct physical surfaces!'
        
        ax.text(0.02, 0.98, info_text,
               transform=ax.transAxes,
               fontsize=9,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        # Legend
        ax.legend(loc='best', fontsize=9, framealpha=0.9)
        ax.grid(True, alpha=0.3)
        
        # Set x-axis limit to cutoff time (or user override)
        if self.config.get('combined_surface_xlim'):
            ax.set_xlim(self.config['combined_surface_xlim'])
        elif time_cutoff is not None:
            # Add 10% padding for visibility
            ax.set_xlim([0, time_cutoff * 1.1])
        
        # Set y-axis limits (let matplotlib auto-scale or use user config)
        if self.config.get('combined_surface_ylim'):
            ax.set_ylim(self.config['combined_surface_ylim'])
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_radius_history(self, pdf):
        """Plot radius evolution for different regions."""
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        # Plot outer boundaries of each region
        if self.data.region_interfaces_indices is not None:
            for region in range(min(3, self.data.region_interfaces_indices.shape[1])):
                radii = []
                times = []
                for t in range(len(self.data.time)):
                    try:
                        idx = int(self.data.region_interfaces_indices[t, region])
                        if idx < len(self.data.zone_boundaries[t]):
                            radii.append(self.data.zone_boundaries[t, idx])
                            times.append(self.data.time[t])
                    except:
                        continue
                
                if radii:
                    labels = ['Hot Spot', 'Fuel', 'Ablator']
                    ax.plot(times, radii, linewidth=2, label=labels[region])
        
        ax.axvline(self.data.bang_time, color='r', linestyle='--', 
                  label='Bang Time', alpha=0.7)
        ax.axvline(self.data.stag_time, color='g', linestyle='--',
                  label='Stagnation', alpha=0.7)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Radius (cm)')
        ax.set_title('Region Boundary Evolution')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_fusion_power(self, pdf):
        """Plot fusion power history."""
        if self.data.fusion_power is None:
            return
            
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        total_fusion_power = np.sum(self.data.fusion_power, axis=1) * 1e-13  # TW
        
        ax.plot(self.data.time, total_fusion_power, 'b-', linewidth=2)
        ax.axvline(self.data.bang_time, color='r', linestyle='--',
                  label='Bang Time', alpha=0.7)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Fusion Power (TW)')
        ax.set_title('Fusion Power vs Time')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_neutron_production(self, pdf):
        """Plot neutron production rate."""
        if self.data.neutron_production_rate is None:
            return
            
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        total_neutron_rate = np.sum(self.data.neutron_production_rate, axis=1)
        
        ax.plot(self.data.time, total_neutron_rate, 'g-', linewidth=2)
        ax.axvline(self.data.bang_time, color='r', linestyle='--',
                  label='Bang Time', alpha=0.7)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Neutron Production Rate (1/s)')
        ax.set_title('Neutron Production vs Time')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_alpha_heating(self, pdf):
        """Plot alpha particle heating."""
        if self.data.alpha_heating_power is None:
            return
            
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        total_alpha_power = np.sum(self.data.alpha_heating_power, axis=1) * 1e-13  # TW
        
        ax.plot(self.data.time, total_alpha_power, 'm-', linewidth=2)
        ax.axvline(self.data.bang_time, color='r', linestyle='--',
                  label='Bang Time', alpha=0.7)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Alpha Heating Power (TW)')
        ax.set_title('Alpha Particle Heating vs Time')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_laser_deposition(self, pdf):
        """Plot laser energy deposition."""
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        led = self.data.laser_energy_deposited
        if led.ndim == 2:
            total_laser_energy = np.sum(led, axis=1) * 1e-6  # J → MJ
        else:
            total_laser_energy = led * 1e-6  # already per-timestep, J → MJ
        
        ax.plot(self.data.time, total_laser_energy, 'orange', linewidth=2)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Cumulative Laser Energy (MJ)')
        ax.set_title('Laser Energy Deposition')
        ax.grid(True, alpha=0.3)
        
        pdf.savefig(fig)
        plt.close(fig)

    def _plot_laser_power(self, pdf):
        """Plot laser power delivered."""
        fig, ax = plt.subplots(figsize=self.default_figsize)

        power_TW = self.data.laser_power_delivered * 1e-12   # W → TW

        ax.plot(self.data.time, power_TW, color='darkorange', linewidth=1.5)

        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Laser Power (TW)')
        ax.set_title('Laser Power Delivered')
        ax.set_xlim(self.data.time[0], self.data.time[-1])
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

        # Annotate peak
        peak_idx = np.argmax(power_TW)
        ax.annotate(f'{power_TW[peak_idx]:.1f} TW',
                    xy=(self.data.time[peak_idx], power_TW[peak_idx]),
                    xytext=(0.7, 0.9), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='->', color='black'),
                    fontsize=11, ha='center')

        pdf.savefig(fig)
        plt.close(fig)
    
    def _plot_drive_temperature(self, pdf):
        """Plot drive temperature from RHW file."""
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        # Convert time to ns for consistency
        time_ns = self.data.drive_time * 1e9
        
        ax.plot(time_ns, self.data.drive_temperature, 'red', linewidth=2, label='Drive Temperature')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Drive Temperature (eV)')
        ax.set_title('Drive Temperature vs Time')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add configuration info as text annotation
        if self.data.rhw_config is not None:
            config_text = f"Drive Type: {self.data.rhw_config.drive_type}\n"
            config_text += f"Burn: {'ON' if self.data.rhw_config.burn_enabled else 'OFF'}"
            ax.text(0.02, 0.98, config_text, 
                   transform=ax.transAxes, fontsize=10,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_density_contour(self, pdf):
        """Create density contour plot."""
        try:
            self._create_contour_plot(
                np.log10(np.maximum(self.data.mass_density, 1e-30)),
                "Log Density",
                "log10(g/cc)",
                "viridis",
                pdf
            )
        except Exception as e:
            logger.warning(f"Could not create density contour: {e}")
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.text(0.5, 0.5, 'Density contour unavailable',
                   ha='center', va='center', transform=ax.transAxes)
            pdf.savefig(fig)
            plt.close(fig)
        
    def _plot_pressure_contour(self, pdf):
        """Create pressure contour plot."""
        try:
            total_pressure = self.data.ion_pressure + self.data.rad_pressure
            self._create_contour_plot(
                np.log10(np.maximum(total_pressure * 1e-8, 1e-30)),
                "Log Pressure",
                "log10(Gbar)",
                "turbo",
                pdf
            )
        except Exception as e:
            logger.warning(f"Could not create pressure contour: {e}")
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.text(0.5, 0.5, 'Pressure contour unavailable',
                   ha='center', va='center', transform=ax.transAxes)
            pdf.savefig(fig)
            plt.close(fig)
        
    def _plot_velocity_contour(self, pdf):
        """Create velocity contour plot."""
        # Velocity is on zone boundaries, need to handle carefully
        try:
            n_times = len(self.data.time)
            # Get minimum number of zones across all times
            min_zones = min(self.data.velocity.shape[1] - 1, 
                          min(len(self.data.zone_boundaries[t]) - 2 for t in range(n_times)))
            
            # Convert to zone centers and km/s
            velocity_centers = np.zeros((n_times, min_zones))
            for t in range(n_times):
                vel = self.data.velocity[t, :min_zones+1]
                # Average velocity at boundaries to get zone center values
                velocity_centers[t] = (vel[:-1] + vel[1:]) / 2 * (-1e-7)  # Convert to km/s
            
            self._create_contour_plot(
                velocity_centers,
                "Velocity",
                "km/s",
                "RdBu_r",
                pdf,
                symmetric=True
            )
        except Exception as e:
            logger.warning(f"Could not create velocity contour: {e}")
            # Create empty plot with message
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.text(0.5, 0.5, 'Velocity contour unavailable',
                   ha='center', va='center', transform=ax.transAxes)
            pdf.savefig(fig)
            plt.close(fig)
        
    def _plot_temperature_contour(self, pdf):
        """Create temperature contour plot."""
        try:
            self._create_contour_plot(
                np.log10(np.maximum(self.data.ion_temperature / 1e3, 1e-30)),
                "Log Temperature",
                "log10(keV)",
                "hot",
                pdf
            )
        except Exception as e:
            logger.warning(f"Could not create temperature contour: {e}")
            fig, ax = plt.subplots(figsize=(12, 8))
            ax.text(0.5, 0.5, 'Temperature contour unavailable',
                   ha='center', va='center', transform=ax.transAxes)
            pdf.savefig(fig)
            plt.close(fig)
        
    def _create_contour_plot(self, data, title, units, cmap, pdf, symmetric=False):
        """Generic contour plot creator."""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        try:
            # Create meshgrid for contour
            n_times, n_zones = data.shape
            time_grid = self.data.time[:n_times]
            
            # Find maximum zone count for regular grid
            max_zones = min(n_zones, 
                          min(len(self.data.zone_boundaries[t]) - 1 for t in range(n_times)))
            
            # Create regular grid for interpolation
            zone_centers = np.zeros((n_times, max_zones))
            data_regular = np.zeros((n_times, max_zones))
            
            for t in range(n_times):
                zones = self.data.zone_boundaries[t]
                centers = (zones[:-1] + zones[1:]) / 2
                
                # Use only zones that exist in all timesteps
                zone_centers[t] = centers[:max_zones]
                data_regular[t] = data[t, :max_zones]
            
            # Create meshgrid
            T, R = np.meshgrid(time_grid, zone_centers[0, :], indexing='ij')
            
            # Interpolate data onto regular radius grid
            Z = np.zeros_like(T)
            for t in range(n_times):
                Z[t, :] = data_regular[t, :]
            
            # Handle NaN and infinite values
            mask = np.isfinite(Z)
            if not np.any(mask):
                logger.warning(f"No finite data for {title} contour plot")
                ax.text(0.5, 0.5, f'No valid data for {title}',
                       ha='center', va='center', transform=ax.transAxes)
                pdf.savefig(fig)
                plt.close(fig)
                return
            
            # Contour plot
            if symmetric:
                vmax = np.nanmax(np.abs(Z[mask]))
                vmin = -vmax
                levels = 50
            else:
                vmin, vmax = np.nanpercentile(Z[mask], [1, 99])
                levels = 50
            
            contour = ax.contourf(T, R, Z,
                                 levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
            
            # Mark important times
            ax.axvline(self.data.bang_time, color='red', linestyle='--',
                      linewidth=2, alpha=0.7, label='Bang Time')
            ax.axvline(self.data.stag_time, color='lime', linestyle='--',
                      linewidth=2, alpha=0.7, label='Stagnation')
            
            ax.set_xlabel('Time (ns)', fontsize=12)
            ax.set_ylabel('Radius (cm)', fontsize=12)
            ax.set_title(f'{title} Evolution', fontsize=14, weight='bold')
            ax.legend(loc='upper right')
            
            cbar = plt.colorbar(contour, ax=ax)
            cbar.set_label(units, fontsize=11)
            
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            
        except Exception as e:
            logger.warning(f"Could not create {title} contour plot: {e}")
            ax.text(0.5, 0.5, f'{title} contour unavailable',
                   ha='center', va='center', transform=ax.transAxes)
            pdf.savefig(fig)
            plt.close(fig)
        
    def _plot_pressure_gradient_contour(self, pdf):
        """Plot pressure gradient contour with shock indicators."""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        try:
            n_times = len(self.data.time)
            total_pressure = self.data.ion_pressure + self.data.rad_pressure
            
            # Find minimum zone count across all times
            min_zones = min(len(self.data.zone_boundaries[t]) - 1 for t in range(n_times))
            
            # Compute pressure gradient
            pressure_gradient = np.zeros((n_times, min_zones))
            zone_centers_grid = np.zeros((n_times, min_zones))
            
            for t in range(n_times):
                pressure = total_pressure[t, :min_zones] * 1e-8  # Gbar
                zones = self.data.zone_boundaries[t]
                zone_centers = (zones[:-1] + zones[1:])[:min_zones] / 2
                
                dp_dr = np.gradient(pressure, zone_centers)
                pressure_gradient[t] = np.abs(dp_dr)
                zone_centers_grid[t] = zone_centers
            
            # Create meshgrid
            T, R = np.meshgrid(self.data.time[:n_times], zone_centers_grid[0, :], indexing='ij')
            Z = np.zeros_like(T)
            
            for t in range(n_times):
                Z[t, :] = pressure_gradient[t, :]
            
            # Handle invalid values
            Z_plot = np.log10(Z + 1e-10)
            Z_plot[~np.isfinite(Z_plot)] = -10
            
            # Contour plot
            contour = ax.contourf(T, R, Z_plot,
                                 levels=50, cmap='Grays')
            
            # Overlay shock tracking points
            if self.data.shock_times:
                ax.scatter(self.data.shock_times, self.data.shock_radii,
                          c='red', s=2, alpha=0.5, label='Detected Shocks')
            
            ax.axvline(self.data.stag_time, color='lime', linestyle='--',
                      linewidth=2, alpha=0.7, label='Stagnation')
            
            ax.set_xlabel('Time (ns)', fontsize=12)
            ax.set_ylabel('Radius (cm)', fontsize=12)
            ax.set_title('Pressure Gradient (Log Scale)', fontsize=14, weight='bold')
            ax.legend()
            
            cbar = plt.colorbar(contour, ax=ax)
            cbar.set_label('log10(|dP/dr|) [Gbar/cm]', fontsize=11)
            
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            
        except Exception as e:
            logger.warning(f"Could not create pressure gradient contour: {e}")
            ax.text(0.5, 0.5, 'Pressure gradient plot unavailable',
                   ha='center', va='center', transform=ax.transAxes)
            pdf.savefig(fig)
            plt.close(fig)
        
    def _plot_shock_tracking(self, pdf):
        """Plot shock tracking analysis with RANSAC line fitting."""
        if not self.data.shock_times:
            logger.info("No shock data to plot")
            return
            
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Convert to numpy arrays
        shock_times = np.array(self.data.shock_times)
        shock_radii = np.array(self.data.shock_radii)
        
        # Plot all shock points
        ax.plot(shock_times, shock_radii, '.', color='lightgray',
               markersize=2, label='All Shock Points', alpha=0.5)
        
        # Filter for RANSAC analysis (focus on main implosion phase)
        time_max = min(self.data.stag_time * 1.2, np.max(shock_times))
        mask = (shock_times <= time_max) & (shock_radii >= 0.02)
        
        if np.sum(mask) > 20 and _HAS_SKLEARN:
            X_filtered = np.column_stack((shock_times[mask], shock_radii[mask]))
            
            # Extract multiple shock fronts using RANSAC
            remaining_points = X_filtered.copy()
            line_models = []
            line_points = []
            max_lines = 3
            min_line_size = 15
            colors = ['tab:orange', 'tab:blue', 'tab:green']
            
            for i in range(max_lines):
                if len(remaining_points) < min_line_size:
                    break
                
                X = remaining_points[:, 0].reshape(-1, 1)
                y = remaining_points[:, 1]
                
                try:
                    ransac = RANSACRegressor(
                        estimator=LinearRegression(),
                        residual_threshold=0.0008,
                        random_state=i
                    )
                    ransac.fit(X, y)
                    inlier_mask = ransac.inlier_mask_
                    
                    if np.sum(inlier_mask) < min_line_size:
                        break
                    
                    line_models.append(ransac.estimator_)
                    line_points.append(remaining_points[inlier_mask])
                    
                    # Plot this shock front
                    x_cluster = remaining_points[inlier_mask, 0]
                    y_cluster = remaining_points[inlier_mask, 1]
                    ax.plot(x_cluster, y_cluster, '.', color=colors[i],
                           markersize=4, label=f'Shock {i+1}')
                    
                    # Plot fit line
                    x_line = np.linspace(np.min(x_cluster), np.max(x_cluster), 100)
                    y_line = ransac.estimator_.predict(x_line.reshape(-1, 1))
                    ax.plot(x_line, y_line, '-', color=colors[i], linewidth=2)
                    
                    # Remove inliers for next iteration
                    remaining_points = remaining_points[~inlier_mask]
                    
                except Exception as e:
                    logger.warning(f"RANSAC iteration {i} failed: {e}")
                    break
        elif np.sum(mask) > 20 and not _HAS_SKLEARN:
            logger.info("sklearn not available — shock tracking plotted without RANSAC fitting")
        
        # Mark stagnation time
        ax.axvline(self.data.stag_time, color='lime', linestyle='--',
                  linewidth=2, alpha=0.7, label='Stagnation')
        
        ax.set_xlim(0, time_max)
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Shock Radius (cm)', fontsize=12)
        ax.set_title('Shock Tracking Analysis', fontsize=14, weight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        
    def _plot_areal_density_evolution(self, pdf):
        """Plot areal density evolution over time."""
        fig, ax = plt.subplots(figsize=self.default_figsize)
        
        # Compute total areal density at each time
        areal_density = np.zeros(len(self.data.time))
        
        for t in range(len(self.data.time)):
            density = self.data.mass_density[t]
            zones = self.data.zone_boundaries[t]
            dr = zones[1:] - zones[:-1]
            areal_density[t] = np.sum(density * dr)
        
        ax.plot(self.data.time, areal_density, 'b-', linewidth=2)
        ax.axvline(self.data.bang_time, color='r', linestyle='--',
                  label='Bang Time', alpha=0.7)
        ax.axvline(self.data.stag_time, color='g', linestyle='--',
                  label='Stagnation', alpha=0.7)
        
        # Mark time-averaged value
        ax.axhline(self.data.time_ave_areal_density, color='purple',
                  linestyle=':', linewidth=2, label='Time-Averaged', alpha=0.7)
        
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Areal Density (g/cm²)')
        ax.set_title('Areal Density Evolution')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        pdf.savefig(fig)
        plt.close(fig)
