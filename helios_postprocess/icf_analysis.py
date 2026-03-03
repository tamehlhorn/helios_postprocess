"""
ICF Analysis Module with Complete Critical Density Tracking
============================================================
Physics analysis functions including BOTH formula-based and laser absorption
methods for tracking the critical density surface.

Adapted for HeliosRun / data_builder pipeline.  Key changes from legacy version:
  - Accepts ICFRunData from helios_postprocess.data_builder
  - Laser energy units: J → MJ  (× 1e-6, not × 1e-13 for ergs)
  - Handles laser_energy_deposited as 1-D (scalar/time) or 2-D (zone-resolved)
  - Neutron-averaged pressure: J/cm³ → Mbar  (× 1e-5, not erg/cm³ path)
  - Fusion yield: prefers TimeIntFusionProd_n_1406 when available
  - rad_pressure is aliased from elec_pressure by data_builder
    (total pressure = ion + elec per CLAUDE.md conventions)

Author: Prof T / Xcimer ICF Analysis
Date: 2025
"""

import numpy as np
from scipy.signal import find_peaks
from sklearn.linear_model import RANSACRegressor
import logging

logger = logging.getLogger(__name__)

# Physical constants
ELECTRON_MASS = 9.109e-31  # kg
PERMITTIVITY_FREE_SPACE = 8.854e-12  # F/m
SPEED_OF_LIGHT = 3.0e8  # m/s
ELEMENTARY_CHARGE = 1.602e-19  # C

class ICFAnalyzer:
    """Main analysis class for ICF simulation data."""
    
    def __init__(self, data, config=None):
        """
        Initialize analyzer with run data.
        
        Parameters
        ----------
        data : ICFRunData
            Container with simulation data
        config : dict, optional
            Configuration options
        """
        self.data = data
        self.config = config or {}
        
    def analyze_drive_phase(self):
        """Analyze drive: laser energy, absorption."""
        logger.info("Analyzing drive phase...")
        
        # Report RHW configuration if available
        if self.data.rhw_config is not None:
            logger.info(f"Drive configuration: {self.data.rhw_config.drive_type}")
            logger.info(f"Fusion reactions: {'ENABLED' if self.data.rhw_config.burn_enabled else 'DISABLED'}")
            
            if self.data.drive_temperature is not None:
                t_drive_ns = self.data.drive_time * 1e9  # Convert to ns
                logger.info(f"Drive temperature profile: {len(self.data.drive_time)} points, "
                          f"{t_drive_ns[0]:.3f} to {t_drive_ns[-1]:.3f} ns")
                logger.info(f"Drive temperature range: {np.min(self.data.drive_temperature):.1f} to "
                          f"{np.max(self.data.drive_temperature):.1f} eV")
        
        # Compute laser energy delivered
        # Helios stores laser energy in Joules (CLAUDE.md); convert J → MJ (× 1e-6).
        # EnLaserDepositedTimeIntg may be 1-D (scalar per time step) or
        # 2-D (zone-resolved) depending on the EXODUS variable.
        if self.data.laser_energy_deposited is not None:
            led = self.data.laser_energy_deposited
            if led.ndim == 2:
                total_laser = np.sum(led, axis=1)          # sum over zones
            else:
                total_laser = led                           # already per-timestep
            self.data.laser_energy = np.max(total_laser) * 1e-6   # J → MJ
            logger.info(f"Laser energy delivered: {self.data.laser_energy:.3f} MJ")
        else:
            self.data.laser_energy = 0.0
            logger.info("No laser energy data available")
        
        # Track critical density position (laser absorption surface)
        self._track_critical_density_position()
        
    def _track_critical_density_position(self):
        """
        Track the position of the critical density surface over time using BOTH methods.
        
        Method 1 (Formula): Calculate n_c from physics and find where n_e = n_c
        Method 2 (Algorithm): Track laser absorption front from notebook algorithm
        
        For a 248 nm KrF laser:
        n_c = (4π² m_e ε₀ c²) / (λ² e²) ≈ 1.85 × 10²² cm⁻³
        """
        logger.info("Computing critical density position...")
        
        # Calculate critical density for 248 nm KrF laser
        laser_wavelength = 248e-9  # meters
        n_critical_formula = (4 * np.pi**2 * ELECTRON_MASS * PERMITTIVITY_FREE_SPACE * 
                             SPEED_OF_LIGHT**2) / \
                            (laser_wavelength**2 * ELEMENTARY_CHARGE**2)
        
        # Convert to cm^-3
        n_critical_cm3 = n_critical_formula * 1e-6
        self.data.critical_density_value = n_critical_cm3
        
        logger.info(f"Critical density (248 nm): {n_critical_cm3:.3e} cm⁻³")
        
        # Initialize radius arrays
        n_times = len(self.data.time)
        critical_radius_formula = np.zeros(n_times)
        critical_radius_algorithm = np.zeros(n_times)
        
        # METHOD 1: Formula-based (find where electron density = n_c)
        method1_success = False
        if self.data.electron_density is not None and self.data.zone_centers is not None:
            logger.info("Using formula method (electron density crossing)...")
            
            for t_idx in range(n_times):
                try:
                    n_e = self.data.electron_density[t_idx]  # cm^-3
                    radii = self.data.zone_centers[t_idx]    # cm
                    
                    if len(n_e) == 0 or len(radii) == 0:
                        continue
                    
                    # Find where n_e crosses n_critical (from outside moving in)
                    if np.max(n_e) > n_critical_cm3:
                        # Find outermost point where n_e > n_c
                        above_critical = n_e > n_critical_cm3
                        
                        if np.any(above_critical):
                            critical_indices = np.where(above_critical)[0]
                            
                            if len(critical_indices) > 0:
                                # Use the outermost (last) index
                                outer_idx = critical_indices[-1]
                                
                                # Linear interpolation for sub-zone accuracy
                                if outer_idx < len(n_e) - 1:
                                    r1, r2 = radii[outer_idx], radii[outer_idx + 1]
                                    n1, n2 = n_e[outer_idx], n_e[outer_idx + 1]
                                    
                                    if abs(n2 - n1) > 1e-10:
                                        alpha = (n_critical_cm3 - n1) / (n2 - n1)
                                        r_critical = r1 + alpha * (r2 - r1)
                                        critical_radius_formula[t_idx] = r_critical
                                    else:
                                        critical_radius_formula[t_idx] = r1
                                else:
                                    critical_radius_formula[t_idx] = radii[outer_idx]
                    
                except Exception as e:
                    logger.debug(f"Formula method failed at t={self.data.time[t_idx]:.3f} ns: {e}")
                    continue
            
            valid_formula = critical_radius_formula[critical_radius_formula > 0]
            if len(valid_formula) > 0:
                method1_success = True
                logger.info(f"Formula method: {len(valid_formula)}/{n_times} timesteps computed")
                logger.info(f"Critical radius (formula): {np.mean(valid_formula):.4f} ± " +
                          f"{np.std(valid_formula):.4f} cm")
        else:
            logger.warning("Formula method unavailable: missing electron_density or zone_centers")
        
        # METHOD 2: Algorithm-based (from notebook - laser absorption front)
        method2_success = False
        if self.data.laser_power_source is not None and self.data.zone_centers is not None:
            logger.info("Using algorithm method (laser absorption front)...")
            
            for t_idx in range(n_times):
                try:
                    laser_pwr = self.data.laser_power_source[t_idx]  # W/cm^3
                    radii = self.data.zone_centers[t_idx]  # cm
                    
                    if len(laser_pwr) == 0 or len(radii) == 0:
                        continue
                    
                    # Find the maximum radius where laser_power_source == 0
                    # This is the absorption surface (critical density surface)
                    # Moving from innermost to outermost zone
                    zone_idx = 0
                    while zone_idx < len(laser_pwr) - 1 and laser_pwr[zone_idx] == 0:
                        zone_idx += 1
                    
                    if zone_idx > 0 and zone_idx < len(radii):
                        # The critical surface is between zone_idx-1 and zone_idx
                        # Use the center of this transition
                        critical_radius_algorithm[t_idx] = radii[zone_idx]
                    
                except Exception as e:
                    logger.debug(f"Algorithm method failed at t={self.data.time[t_idx]:.3f} ns: {e}")
                    continue
            
            valid_algorithm = critical_radius_algorithm[critical_radius_algorithm > 0]
            if len(valid_algorithm) > 0:
                method2_success = True
                logger.info(f"Algorithm method: {len(valid_algorithm)}/{n_times} timesteps computed")
                logger.info(f"Critical radius (algorithm): {np.mean(valid_algorithm):.4f} ± " +
                          f"{np.std(valid_algorithm):.4f} cm")
        else:
            logger.warning("Algorithm method unavailable: missing laser_power_source or zone_centers")
        
        # Store results
        if method1_success:
            self.data.critical_density_radius_formula = critical_radius_formula
        else:
            self.data.critical_density_radius_formula = None
            
        if method2_success:
            self.data.critical_density_radius_algorithm = critical_radius_algorithm
        else:
            self.data.critical_density_radius_algorithm = None
        
        if not (method1_success or method2_success):
            logger.warning("Critical density tracking failed for both methods")
    
    # =========================================================================
    # REST OF THE ANALYSIS METHODS (UNCHANGED)
    # =========================================================================
    # Copy the rest of your existing icf_analysis.py methods here:
    # - analyze_implosion_phase()
    # - analyze_stagnation_phase()
    # - analyze_burn_phase()
    # - compute_performance_metrics()
    # - All helper methods (_track_shock_fronts, etc.)
    
    def analyze_implosion_phase(self):
        """Analyze implosion: velocity, shocks, adiabat."""
        logger.info("Analyzing implosion phase...")
        
        # Track peak implosion velocity
        if self.data.velocity is not None:
            # Find maximum inward velocity (most negative)
            min_velocities = np.min(self.data.velocity, axis=1)
            self.data.peak_implosion_velocity = np.min(min_velocities) * 1e-5  # Convert to km/s
            logger.info(f"Peak implosion velocity: {abs(self.data.peak_implosion_velocity):.2f} km/s")
        
        # Detect and track shock fronts
        self._track_shock_fronts()
        
        # Compute adiabat
        self._compute_adiabat()
        
        # Track ablation front position
        self._track_ablation_front()
        
    def _track_shock_fronts(self):
        """Track shock fronts using pressure gradient detection with RANSAC fitting."""
        if self.data.ion_pressure is None or self.data.zone_boundaries is None:
            logger.warning("Cannot track shocks: missing pressure or position data")
            return
            
        logger.info("Tracking shock fronts...")
        
        # Compute total pressure
        if self.data.rad_pressure is not None:
            pressure = (self.data.ion_pressure + self.data.rad_pressure) * 1e-5  # Convert to Mbar
        else:
            pressure = self.data.ion_pressure * 1e-5
        
        shock_positions = []
        shock_times = []
        shock_velocities = []
        
        threshold = self.config.get('shock_detection_threshold', 50.0)  # Mbar/cm
        
        for t_idx in range(len(self.data.time)):
            try:
                p = pressure[t_idx]
                boundaries = self.data.zone_boundaries[t_idx]
                radii = (boundaries[:-1] + boundaries[1:]) / 2
                
                if len(p) < 5:
                    continue
                
                # Compute pressure gradient
                dp_dr = np.gradient(p, radii)
                
                # Find peaks in pressure gradient (potential shocks)
                peaks, properties = find_peaks(np.abs(dp_dr), height=threshold, distance=5)
                
                if len(peaks) > 0:
                    for peak in peaks:
                        shock_positions.append(radii[peak])
                        shock_times.append(self.data.time[t_idx])
                        
                        # Estimate shock velocity
                        if t_idx > 0 and len(shock_times) > 1:
                            dt = self.data.time[t_idx] - self.data.time[t_idx-1]
                            if dt > 0:
                                dr = radii[peak] - shock_positions[-2]
                                v_shock = dr / (dt * 1e-9)  # cm/s
                                shock_velocities.append(v_shock)
                
            except Exception as e:
                logger.debug(f"Shock tracking failed at t={self.data.time[t_idx]:.3f} ns: {e}")
                continue
        
        self.data.shock_times = shock_times
        self.data.shock_radii = shock_positions
        if len(shock_velocities) > 0:
            self.data.shock_velocities = shock_velocities
        
        if len(shock_times) > 0:
            logger.info(f"Tracked {len(shock_times)} shock front points")
    
    def _compute_adiabat(self):
        """Compute adiabat (entropy parameter) in ice region."""
        if self.data.ion_pressure is None or self.data.mass_density is None:
            logger.warning("Cannot compute adiabat: missing pressure or density data")
            return
        
        logger.info("Computing adiabat...")
        
        # Adiabat α = P / P_Fermi where P_Fermi = (ħ²/5m)(3π²)^(2/3) * (ρ/m_u)^(5/3)
        # Simplified: α ≈ P / (ρ^(5/3))
        
        try:
            # Use early time (before strong compression)
            early_idx = min(10, len(self.data.time) // 10)
            
            pressure = self.data.ion_pressure[early_idx] * 1e-5  # Mbar
            density = self.data.mass_density[early_idx]  # g/cc
            
            # Mass-weighted average
            mass = self.data.zone_mass[early_idx]
            
            # α = P / (ρ^(5/3)) with appropriate units
            adiabat_zones = pressure / (density ** (5/3))
            self.data.adiabat_mass_averaged_ice = np.average(adiabat_zones, weights=mass)
            
            logger.info(f"Mass-averaged adiabat: {self.data.adiabat_mass_averaged_ice:.4f}")
            
        except Exception as e:
            logger.warning(f"Could not compute adiabat: {e}")
            self.data.adiabat_mass_averaged_ice = 0.0
    
    def _track_ablation_front(self):
        """
        Track ablation front position over time.
        
        The ablation front is defined as the location with the minimum (most negative)
        density scale length. This represents the steepest density gradient.
        
        Based on algorithm from original notebook.
        """
        if self.data.scale_length is None:
            logger.warning("Cannot track ablation front: scale length not computed")
            return
        
        logger.info("Tracking ablation front position...")
        
        n_times, n_zones = self.data.scale_length.shape
        
        # Find minimum scale length at each time step
        ablation_front_indices = np.zeros(n_times, dtype=int)
        min_scale_length = np.zeros(n_times)
        
        for t in range(n_times):
            # Initialize with first zone (take absolute value as starting point)
            min_scale_length[t] = -np.abs(self.data.scale_length[t, 0])
            ablation_front_indices[t] = 0
            
            # Find most negative scale length
            for i in range(n_zones):
                scale_len = self.data.scale_length[t, i]
                
                # Look for negative scale lengths that are more negative than current min
                if scale_len > min_scale_length[t] and scale_len < 0:
                    ablation_front_indices[t] = i
                    min_scale_length[t] = scale_len
        
        # Convert indices to radii
        ablation_front_radius = np.zeros(n_times)
        for t in range(n_times):
            idx = int(ablation_front_indices[t])
            if idx < len(self.data.zone_boundaries[t]) - 1:
                # Use zone center
                ablation_front_radius[t] = (self.data.zone_boundaries[t, idx] + 
                                           self.data.zone_boundaries[t, idx + 1]) / 2
        
        # Smooth the ablation front profile
        ablation_front_radius = self._smooth_ablation_front(ablation_front_radius)
        
        # Store results
        self.data.ablation_front_radius = ablation_front_radius
        self.data.min_scale_length = min_scale_length
        
        # Calculate some statistics
        valid_radii = ablation_front_radius[ablation_front_radius > 0]
        if len(valid_radii) > 0:
            logger.info(f"Ablation front tracked: {len(valid_radii)}/{n_times} timesteps")
            logger.info(f"Initial radius: {valid_radii[0]:.4f} cm")
            if len(valid_radii) > 1:
                logger.info(f"Final radius: {valid_radii[-1]:.4f} cm")
        else:
            logger.warning("No valid ablation front positions found")
    
    def _smooth_ablation_front(self, radii: np.ndarray, iterations: int = 10) -> np.ndarray:
        """
        Smooth ablation front radius profile.
        
        Uses iterative smoothing algorithm from original notebook.
        
        Parameters
        ----------
        radii : np.ndarray
            Raw ablation front radii
        iterations : int
            Number of smoothing iterations
            
        Returns
        -------
        np.ndarray
            Smoothed ablation front radii
        """
        if self.data.stag_time == 0:
            logger.warning("Stagnation time not set, using middle of simulation for smoothing")
            stag_idx = len(radii) // 2
        else:
            # Find stagnation index
            stag_idx = np.argmin(np.abs(self.data.time - self.data.stag_time))
        
        smoothed = radii.copy()
        n_times = len(smoothed)
        
        for iteration in range(iterations):
            # Smooth backward from stagnation
            for t in range(stag_idx, 0, -1):
                if smoothed[t-1] < smoothed[t]:
                    smoothed[t-1] = smoothed[t]
            
            # Smooth forward from stagnation
            for t in range(stag_idx, n_times - 2):
                if smoothed[t] > smoothed[t+1]:
                    smoothed[t] = smoothed[t+1]
            
            # Average out plateaus
            for t in range(1, n_times - 1):
                if smoothed[t] == smoothed[t+1] and t > 0:
                    smoothed[t] = (smoothed[t+1] + smoothed[t-1]) / 2
        
        return smoothed
    
    def analyze_stagnation_phase(self):
        """Analyze stagnation: bang time, max compression, hot spot."""
        logger.info("Analyzing stagnation phase...")
        
        # Determine stagnation time (max density)
        self._find_stagnation_time()
        
        # Determine bang time (peak fusion power)
        self._find_bang_time()
        
        # Compute hot spot properties at stagnation
        self._compute_hot_spot_properties()
        
        # Compute areal densities
        self._compute_areal_densities()
    
    def _find_stagnation_time(self):
        """Find stagnation time from maximum density."""
        if self.data.mass_density is None:
            logger.warning("Cannot find stagnation time: missing density data")
            return
        
        # Find time of maximum mass-averaged density
        mass_avg_density = np.zeros(len(self.data.time))
        for t in range(len(self.data.time)):
            mass_avg_density[t] = np.average(
                self.data.mass_density[t],
                weights=self.data.zone_mass[t]
            )
        
        stag_idx = np.argmax(mass_avg_density)
        self.data.stag_time = self.data.time[stag_idx]
        self.data.max_density = mass_avg_density[stag_idx]
        
        logger.info(f"Stagnation time: {self.data.stag_time:.3f} ns")
        logger.info(f"Maximum density: {self.data.max_density:.2f} g/cc")
    
    def _find_bang_time(self):
        """Find bang time from peak fusion power."""
        if self.data.fusion_power is None:
            logger.warning("Could not determine bang time: no fusion data")
            self.data.bang_time = self.data.stag_time
            return
        
        # Find time of peak fusion power
        total_fusion = np.sum(self.data.fusion_power, axis=1)
        bang_idx = np.argmax(total_fusion)
        self.data.bang_time = self.data.time[bang_idx]
        
        logger.info(f"Bang time: {self.data.bang_time:.3f} ns")
    
    def _compute_hot_spot_properties(self):
        """Compute hot spot properties at stagnation."""
        if self.data.stag_time == 0.0:
            logger.warning("Cannot compute hot spot properties: stagnation time not determined")
            return
        
        # Find timestep closest to stagnation
        stag_idx = np.argmin(np.abs(self.data.time - self.data.stag_time))
        
        temp_threshold = self.config.get('hot_spot_temp_threshold', 1000.0)  # eV
        
        try:
            temperatures = self.data.ion_temperature[stag_idx]
            hot_spot_mask = temperatures > temp_threshold
            
            if np.any(hot_spot_mask):
                boundaries = self.data.zone_boundaries[stag_idx]
                radii = (boundaries[:-1] + boundaries[1:]) / 2
                
                # Hot spot radius (outermost hot zone)
                hot_radii = radii[hot_spot_mask]
                self.data.stagnation_hot_spot_radius = np.max(hot_radii)
                
                # Hot spot pressure (mass-averaged)
                pressure = (self.data.ion_pressure[stag_idx] + 
                           self.data.rad_pressure[stag_idx]) * 1e-5  # Mbar
                mass = self.data.zone_mass[stag_idx]
                
                self.data.hot_spot_pressure = np.average(
                    pressure[hot_spot_mask],
                    weights=mass[hot_spot_mask]
                )
                
                # Hot spot areal density
                density = self.data.mass_density[stag_idx]
                dr = boundaries[1:] - boundaries[:-1]
                self.data.hot_spot_areal_density = np.sum(
                    density[hot_spot_mask] * dr[hot_spot_mask]
                )
                
                logger.info(f"Hot spot radius: {self.data.stagnation_hot_spot_radius:.4f} cm")
                logger.info(f"Hot spot pressure: {self.data.hot_spot_pressure:.2f} Mbar")
                
        except Exception as e:
            logger.warning(f"Could not compute hot spot properties: {e}")
    
    def _compute_areal_densities(self):
        """
        Compute areal density (ρR) for different regions vs time.
        
        Areal density at radius r: ∫[r to R_outer] ρ(r') dr'
        This is the column density from r outward to the outer boundary.
        """
        logger.info("Computing areal densities...")
        
        if self.data.zone_boundaries is None or self.data.mass_density is None:
            logger.warning("Cannot compute areal densities: missing boundary or density data")
            return
        
        n_times, n_zones = self.data.mass_density.shape
        
        # Compute areal density at each zone boundary vs time
        areal_density_array = np.zeros((n_times, n_zones + 1))
        
        for t in range(n_times):
            boundaries = self.data.zone_boundaries[t]
            density = self.data.mass_density[t]
            zone_masses = self.data.zone_mass[t]
            
            # For each boundary position, sum mass outward divided by area at that position
            for i in range(len(boundaries)):
                r_i = boundaries[i]
                area_i = 4 * np.pi * r_i**2 if r_i > 0 else 1e-30
                
                # Sum mass from this radius outward
                mass_outward = np.sum(zone_masses[i:])
                areal_density_array[t, i] = mass_outward / area_i
        
        # Store for later use
        self.data.areal_density_vs_time = areal_density_array
        
        # Compute key values at specific times
        # Bang time areal density
        if self.data.bang_time > 0:
            bang_idx = np.argmin(np.abs(self.data.time - self.data.bang_time))
            
            # Total areal density (at center)
            self.data.bang_time_areal_density = areal_density_array[bang_idx, 0]
            
            # Fuel areal density (at fuel boundary)
            if self.data.region_interfaces_indices is not None:
                fuel_idx = int(self.data.region_interfaces_indices[bang_idx, 1])
                self.data.bang_time_fuel_areal_density = areal_density_array[bang_idx, fuel_idx]
                
                # Hot spot areal density (at hot spot boundary)
                hs_idx = int(self.data.region_interfaces_indices[bang_idx, 0])
                self.data.bang_time_hs_areal_density = areal_density_array[bang_idx, hs_idx]
                
                logger.info(f"Bang time areal densities:")
                logger.info(f"  Total (center):  {self.data.bang_time_areal_density:.4f} g/cm²")
                logger.info(f"  Fuel (boundary): {self.data.bang_time_fuel_areal_density:.4f} g/cm²")
                logger.info(f"  Hot spot:        {self.data.bang_time_hs_areal_density:.4f} g/cm²")
        
        # Time-averaged areal density (simple average over burn period)
        if self.data.bang_time > 0 and self.data.burn_width > 0:
            t_start = self.data.bang_time - self.data.burn_width / 2
            t_end = self.data.bang_time + self.data.burn_width / 2
            
            mask = (self.data.time >= t_start) & (self.data.time <= t_end)
            if np.any(mask):
                self.data.time_ave_areal_density = np.mean(areal_density_array[mask, 0])
                logger.info(f"Time-averaged areal density (burn period): "
                          f"{self.data.time_ave_areal_density:.4f} g/cm²")
    
    def analyze_burn_phase(self):
        """Analyze burn: fusion yield, burn width, neutron-weighted quantities."""
        logger.info("Analyzing burn phase...")
        
        # Check RHW configuration
        if self.data.rhw_config is not None:
            if not self.data.rhw_config.burn_enabled:
                logger.warning("⚠️  RHW configuration indicates fusion reactions are DISABLED")
                logger.warning("    Burn analysis may not be meaningful for this simulation")
        
        if self.data.fusion_power is None:
            logger.warning("No fusion data available")
            
            # Additional context if RHW config is available
            if self.data.rhw_config is not None and not self.data.rhw_config.burn_enabled:
                logger.info("This is expected - fusion was disabled in the simulation")
            
            return
        
        # Compute fusion yield
        self._compute_fusion_yield()
        
        # Compute burn width
        self._compute_burn_width()
        
        # Compute neutron-averaged quantities
        self._compute_neutron_averaged_quantities()
        
        # Maximum DT temperature
        if self.data.ion_temperature is not None:
            self.data.max_dt_temp = np.max(self.data.ion_temperature)
            logger.info(f"Maximum DT temperature: {self.data.max_dt_temp:.2f} eV")
    
    def _compute_fusion_yield(self):
        """
        Compute total fusion energy output.

        Preferred method (per CLAUDE.md):
            yield = TimeIntFusionProd_n_1406 × 17.6 MeV × 1.602e-19 MJ/MeV

        Fallback: integrate zone-level FusionRate_DT_nHe4 over time
        (less precise due to time-grid discretisation).
        """
        # ------------------------------------------------------------------
        # Preferred: time-integrated DT neutron count (already in EXODUS)
        # ------------------------------------------------------------------
        if hasattr(self.data, 'dt_neutron_count') and self.data.dt_neutron_count is not None:
            neutron_count = self.data.dt_neutron_count
            # May be zone-resolved (2-D) or scalar (1-D / 0-D)
            if neutron_count.ndim >= 2:
                total_neutrons = np.sum(neutron_count[-1])   # last time step, all zones
            elif neutron_count.ndim == 1:
                total_neutrons = neutron_count[-1]            # last time step
            else:
                total_neutrons = float(neutron_count)

            MeV_per_reaction = 17.6
            MJ_per_MeV = 1.602e-19
            self.data.energy_output = total_neutrons * MeV_per_reaction * MJ_per_MeV
            logger.info(f"Fusion energy output (from neutron count): "
                        f"{self.data.energy_output:.3f} MJ  "
                        f"({total_neutrons:.3e} DT neutrons)")
            return

        # ------------------------------------------------------------------
        # Fallback: integrate fusion rate over time (CLAUDE.md warns this is
        # imprecise due to the time grid — use only when neutron count unavailable)
        # ------------------------------------------------------------------
        if self.data.fusion_power is not None:
            logger.warning("TimeIntFusionProd_n_1406 not available — "
                           "falling back to time-integration of fusion rate "
                           "(less precise, see CLAUDE.md)")
            total_fusion = np.sum(self.data.fusion_power, axis=1)
            fusion_energy = np.trapz(total_fusion, self.data.time * 1e-9)  # ns→s

            # NOTE: The result of this integration has units that depend on
            # what FusionRate_DT_nHe4 actually represents in Helios.
            # If it is reactions/s/cm³ the integral gives total reactions;
            # multiply by 17.6 MeV * 1.602e-13 J/MeV then × 1e-6 for MJ.
            # If Helios stores it as power (W), trapz gives Joules.
            # We assume reaction-rate here (consistent with variable name).
            MeV_per_reaction = 17.6
            J_per_MeV = 1.602e-13
            energy_J = fusion_energy * MeV_per_reaction * J_per_MeV
            self.data.energy_output = energy_J * 1e-6   # J → MJ
            logger.info(f"Fusion energy output (time-integrated): "
                        f"{self.data.energy_output:.3f} MJ")
        else:
            self.data.energy_output = 0.0
            logger.warning("No fusion data available for yield calculation")
    
    def _compute_burn_width(self):
        """Compute burn width (FWHM of fusion power)."""
        total_fusion = np.sum(self.data.fusion_power, axis=1)
        peak_power = np.max(total_fusion)
        
        threshold = self.config.get('burn_width_threshold', 0.5) * peak_power
        
        above_threshold = total_fusion > threshold
        if np.any(above_threshold):
            indices = np.where(above_threshold)[0]
            t_start = self.data.time[indices[0]]
            t_end = self.data.time[indices[-1]]
            self.data.burn_width = t_end - t_start
            
            logger.info(f"Burn width (FWHM): {self.data.burn_width:.3f} ns")
    
    def _compute_neutron_averaged_quantities(self):
        """
        Compute neutron-weighted average quantities using the original notebook algorithm.
        
        These quantities are averaged over time AND space, weighted by:
        - Fusion rate in each zone (neutrons produced per zone)
        - Zone mass (for temperature and pressure)
        
        Normalization: divide by total neutron yield
        
        Key differences from simple time averaging:
        - Weights by WHERE and WHEN fusion occurs
        - Gives conditions during actual fusion reactions
        """
        logger.info("Computing neutron-averaged quantities...")
        
        # Check if we have fusion power data (used as fusion rate)
        if self.data.fusion_power is None:
            logger.warning("No fusion power data available for neutron averaging")
            logger.info("Neutron-averaged quantities will be zero")
            return
        
        # fusion_power is in Watts, but we need it as a rate (neutrons/s per zone)
        # For DT fusion: each reaction produces 1 neutron and 17.6 MeV energy
        # So: N_neutrons/s = Power(W) / (17.6 MeV * 1.602e-13 J/MeV)
        MeV_per_reaction = 17.6
        J_per_MeV = 1.602e-13
        energy_per_reaction = MeV_per_reaction * J_per_MeV
        
        # Fusion rate array: neutrons/s per zone
        fusion_rate = self.data.fusion_power / energy_per_reaction  # [time, space]
        
        # Total fusion rate vs time (sum over all zones)
        fusion_rate_vec = np.sum(fusion_rate, axis=1)  # [time]
        
        # Check if there's any fusion
        if np.max(fusion_rate_vec) == 0:
            logger.warning("Fusion rate is zero at all times")
            logger.info("Neutron-averaged quantities cannot be computed")
            return
        
        # Compute total neutron yield
        time_seconds = self.data.time * 1e-9  # Convert ns to s
        dt_array = np.diff(time_seconds)
        
        # Neutron yield = integral of fusion rate over time
        neutron_yield = 0.0
        for t in range(len(dt_array)):
            neutron_yield += fusion_rate_vec[t] * dt_array[t]
        
        if neutron_yield == 0:
            logger.warning("Total neutron yield is zero")
            return
        
        logger.info(f"Total neutron yield: {neutron_yield:.3e} neutrons")
        
        # 1. Neutron-averaged FUEL areal density
        # Fuel includes first two layers: DT vapor + DT Ice (or DT vapor + wetted foam)
        if self.data.region_interfaces_indices is not None:
            # Outer spatial index of the second zone at beginning time
            fuel_index = int(self.data.region_interfaces_indices[0, 1])
            
            # Need areal density array - compute if not available
            if not hasattr(self.data, 'areal_density_vs_time'):
                logger.warning("Areal density array not computed - will compute simplified version")
                # Simplified: use bang_time value as proxy
                if self.data.bang_time_fuel_areal_density > 0:
                    self.data.neutron_ave_fuel_areal_density = self.data.bang_time_fuel_areal_density
                    logger.info(f"Neutron-averaged fuel areal density (approximated): "
                              f"{self.data.neutron_ave_fuel_areal_density:.4f} g/cm²")
            else:
                # Use time-resolved areal density at fuel boundary
                areal_density_fuel = self.data.areal_density_vs_time[:, fuel_index]
                
                areal_sum = 0.0
                for t in range(len(dt_array)):
                    areal_sum += areal_density_fuel[t] * fusion_rate_vec[t] * dt_array[t]
                
                self.data.neutron_ave_fuel_areal_density = areal_sum / neutron_yield
                logger.info(f"Neutron-averaged fuel areal density: "
                          f"{self.data.neutron_ave_fuel_areal_density:.4f} g/cm²")
        else:
            logger.warning("Region interfaces not identified - cannot compute fuel areal density")
        
        # 2. Neutron-averaged ion temperature
        # Double loop over time and space, weighted by fusion_rate and zone_mass
        if self.data.ion_temperature is not None and self.data.zone_mass is not None:
            ion_temp_sum = 0.0
            
            n_times = len(dt_array)
            for t in range(n_times):
                dt = dt_array[t]
                n_zones = self.data.ion_temperature.shape[1]
                
                # Handle variable zone counts
                n_zones = min(n_zones, fusion_rate.shape[1], self.data.zone_mass.shape[1])
                
                for z in range(n_zones):
                    ion_temp_sum += (self.data.ion_temperature[t, z] * 
                                    fusion_rate[t, z] * 
                                    dt * 
                                    self.data.zone_mass[t, z])
            
            self.data.neutron_ave_ion_temperature = ion_temp_sum / neutron_yield
            logger.info(f"Neutron-averaged ion temperature: "
                       f"{self.data.neutron_ave_ion_temperature:.2f} eV")
        else:
            logger.warning("Ion temperature or zone mass not available")
        
        # 3. Neutron-averaged pressure
        # Same double loop, weighted by fusion_rate and zone_mass
        if self.data.ion_pressure is not None and self.data.zone_mass is not None:
            pressure_sum = 0.0
            
            n_times = len(dt_array)
            for t in range(n_times):
                dt = dt_array[t]
                n_zones = self.data.ion_pressure.shape[1]
                
                # Handle variable zone counts
                n_zones = min(n_zones, fusion_rate.shape[1], self.data.zone_mass.shape[1])
                
                for z in range(n_zones):
                    pressure_sum += (self.data.ion_pressure[t, z] * 
                                    fusion_rate[t, z] * 
                                    dt * 
                                    self.data.zone_mass[t, z])
            
            neutron_ave_pressure = pressure_sum / neutron_yield
            
            # Convert units: Helios pressure is in J/cm³
            # 1 Mbar = 1e5 J/cm³  →  P_Mbar = P_Jcm3 × 1e-5
            self.data.neutron_ave_pressure = neutron_ave_pressure * 1e-5
            
            logger.info(f"Neutron-averaged pressure: "
                       f"{self.data.neutron_ave_pressure:.2f} Mbar")
        else:
            logger.warning("Ion pressure or zone mass not available")
    
    def compute_performance_metrics(self):
        """Compute overall performance metrics."""
        logger.info("Computing performance metrics...")
        
        # Target gain
        if self.data.laser_energy > 0 and self.data.energy_output > 0:
            self.data.target_gain = self.data.energy_output / self.data.laser_energy
            logger.info(f"Target gain: {self.data.target_gain:.3f}")
        
        # Compression ratio
        if self.data.mass_density is not None and self.data.max_density > 0:
            initial_density = np.mean(self.data.mass_density[0])
            if initial_density > 0:
                self.data.comp_ratio = self.data.max_density / initial_density
                logger.info(f"Compression ratio: {self.data.comp_ratio:.2f}")
