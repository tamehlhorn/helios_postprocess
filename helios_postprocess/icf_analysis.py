"""
ICF Analysis Module with Complete Critical Density Tracking
============================================================
Physics analysis functions including BOTH formula-based and laser absorption
methods for tracking the critical density surface.

Adapted for HeliosRun / data_builder pipeline.  Key changes from legacy version:
  - Accepts ICFRunData from helios_postprocess.data_builder
  - Laser energy units: J → MJ  (× 1e-6, not × 1e-13 for ergs)
  - Handles laser_energy_deposited as 1-D (scalar/time) or 2-D (zone-resolved)
  - Pressures stored/reported in Gbar (× 1e-8 from J/cm³) per CLAUDE.md
  - Temperatures stored/reported in keV (÷ 1000 from eV) per ICF convention
  - Neutron-averaged pressure: J/cm³ → Gbar  (× 1e-8)
  - Fusion yield: prefers TimeIntFusionProd_n_1406 when available
  - FusionRate_DT_nHe4 treated as reaction rate (not power) for neutron weighting
  - rad_pressure is aliased from elec_pressure by data_builder
    (total pressure = ion + elec per CLAUDE.md conventions)

Author: Prof T / Xcimer ICF Analysis
Date: 2025
"""

import numpy as np
from scipy.signal import find_peaks
import logging

try:
    from sklearn.linear_model import RANSACRegressor
    _HAS_SKLEARN = True
except ImportError:
    _HAS_SKLEARN = False

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
        
        # Track peak implosion velocity (only before bang time)
        if self.data.velocity is not None:
            if self.data.bang_time > 0:
                pre_bang = self.data.time <= self.data.bang_time
            else:
                pre_bang = np.ones(len(self.data.time), dtype=bool)
            min_velocities = np.min(self.data.velocity[pre_bang], axis=1)
            self.data.peak_implosion_velocity = np.min(min_velocities) * 1e-5  # cm/s → km/s
            logger.info(f"Peak implosion velocity (pre-bang): "
                        f"{abs(self.data.peak_implosion_velocity):.2f} km/s")
        
        # Detect and track shock fronts
        self._track_shock_fronts()
        
        # Compute adiabat
        self._compute_adiabat()
        
        # Track ablation front position
        self._track_ablation_front()

        # Compute in-flight aspect ratio
        self._compute_ifar()
        
    def _track_shock_fronts(self):
        """Track shock fronts using pressure gradient detection with RANSAC fitting."""
        if self.data.ion_pressure is None or self.data.zone_boundaries is None:
            logger.warning("Cannot track shocks: missing pressure or position data")
            return
            
        logger.info("Tracking shock fronts...")
        
        # Compute total pressure in Gbar
        if self.data.rad_pressure is not None:
            pressure = (self.data.ion_pressure + self.data.rad_pressure) * 1e-8  # J/cm³ → Gbar
        else:
            pressure = self.data.ion_pressure * 1e-8
        
        shock_positions = []
        shock_times = []
        shock_velocities = []
        
        threshold = self.config.get('shock_detection_threshold', 0.05)  # Gbar/cm
        
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
        """
        Compute mass-averaged adiabat (entropy parameter) in cold DT fuel.

        alpha = P / P_Fermi   where P_Fermi = 2.17 (rho/rho_0)^(5/3) Mbar
        for equimolar DT ice (rho_0 = 0.205 g/cc, Lindl convention).

        Evaluated at peak implosion velocity in the cold fuel (DT Solid / DT ice)
        region, excluding the hot spot and ablator.  At this time the first shock
        has traversed the fuel and the shell is in-flight, giving a physically
        meaningful measurement of fuel entropy before stagnation reheating.
        """
        if self.data.ion_pressure is None or self.data.mass_density is None:
            logger.warning("Cannot compute adiabat: missing pressure or density data")
            return

        logger.info("Computing adiabat...")

        try:
            # -- Time: peak implosion velocity (most negative shell velocity) --
            # This is when the shell is in-flight after shock transit but before
            # stagnation reheating -- the standard time for adiabat evaluation.
            ri = self.data.region_interfaces_indices
            vel = self.data.velocity
            n_zones = self.data.mass_density.shape[1]

            if vel is not None:
                # Zone-center velocities
                if vel.shape[1] > n_zones:
                    v_zone = 0.5 * (vel[:, :n_zones] + vel[:, 1:n_zones+1])
                else:
                    v_zone = vel[:, :n_zones]

                # Find time of peak inward velocity in fuel zones only
                if ri is not None and ri.shape[1] >= 2:
                    fuel_start = int(ri[0, 0])
                    fuel_end = int(ri[0, -2])
                    if fuel_end > fuel_start:
                        min_v_per_time = np.min(v_zone[:, fuel_start:fuel_end], axis=1)
                    else:
                        min_v_per_time = np.min(v_zone, axis=1)
                else:
                    min_v_per_time = np.min(v_zone, axis=1)

                # Restrict to times before stagnation to avoid post-bounce
                stag_t = self.data.stag_time if self.data.stag_time > 0 else self.data.time[-1]
                pre_stag = self.data.time < stag_t
                if np.any(pre_stag):
                    masked_v = min_v_per_time.copy()
                    masked_v[~pre_stag] = 0.0  # zero out post-stagnation
                    eval_idx = np.argmin(masked_v)
                else:
                    eval_idx = np.argmin(min_v_per_time)
            else:
                # Fallback: mid-implosion
                stag_t = self.data.stag_time if self.data.stag_time > 0 else self.data.time[-1] / 2
                eval_idx = np.argmin(np.abs(self.data.time - stag_t * 0.5))

            # -- Zone range: cold fuel = between hot-spot and ablator boundaries --
            if ri is not None and ri.shape[1] >= 2:
                z_start = int(ri[eval_idx, 0])           # outer edge of hot spot
                z_end   = int(ri[eval_idx, 1])            # outer edge of cold fuel (exclude ablated region)
            else:
                z_start = 0
                z_end   = n_zones // 2

            if z_end <= z_start:
                logger.warning("Cold fuel region is empty -- cannot compute adiabat")
                return

            # -- Total pressure in Mbar --
            p_tot = self.data.ion_pressure[eval_idx, z_start:z_end]
            if self.data.rad_pressure is not None:
                p_tot = p_tot + self.data.rad_pressure[eval_idx, z_start:z_end]
            p_Mbar = p_tot * 1e-5                        # J/cm3 -> Mbar

            rho  = self.data.mass_density[eval_idx, z_start:z_end]   # g/cc
            mass = self.data.zone_mass[eval_idx, z_start:z_end]

            # -- Fermi pressure for equimolar DT --
            rho0 = 0.205                                  # g/cc (DT ice density)
            p_fermi = 2.17 * (rho / rho0) ** (5.0 / 3.0)  # Mbar

            with np.errstate(divide='ignore', invalid='ignore'):
                alpha = p_Mbar / p_fermi
                alpha[~np.isfinite(alpha)] = 0.0

            self.data.adiabat_mass_averaged_ice = np.average(alpha, weights=mass)
            logger.info(f"Mass-averaged adiabat (cold fuel, t={self.data.time[eval_idx]:.2f} ns): "
                        f"{self.data.adiabat_mass_averaged_ice:.2f}")

        except Exception as e:
            logger.warning(f"Could not compute adiabat: {e}")
            self.data.adiabat_mass_averaged_ice = 0.0

    def _track_ablation_front(self):
        """
        Track ablation front position over time.

        The ablation front is the boundary between the hot, expanding corona
        (low density) and the cold, unablated dense shell.

        Algorithm: at each timestep, find the steepest *negative* density
        gradient  (-dρ/dr  maximised) in the outer portion of the target
        (outside the hot-spot region).  The zone outer boundary at that
        location is the ablation front radius.
        """
        if self.data.mass_density is None or self.data.zone_boundaries is None:
            logger.warning("Cannot track ablation front: missing density or boundary data")
            return

        logger.info("Tracking ablation front position...")

        n_times, n_zones = self.data.mass_density.shape
        ri = self.data.region_interfaces_indices

        ablation_front_radius = np.zeros(n_times)
        ablation_front_indices = np.zeros(n_times, dtype=int)

        for t in range(n_times):
            rho = self.data.mass_density[t]
            boundaries = self.data.zone_boundaries[t]
            zone_centers = 0.5 * (boundaries[:-1] + boundaries[1:])

            # Only search outside the hot-spot region
            if ri is not None:
                search_start = int(ri[t, 0])          # outer edge of hot spot
            else:
                search_start = n_zones // 4           # fallback: skip inner quarter

            if search_start >= n_zones - 2:
                continue

            # Compute dρ/dr on the search range
            r_search   = zone_centers[search_start:]
            rho_search = rho[search_start:]
            drho_dr = np.gradient(rho_search, r_search)

            # Ablation front = steepest density drop (most negative dρ/dr)
            idx_in_search = np.argmin(drho_dr)
            ablation_idx = search_start + idx_in_search

            # Require a meaningful negative gradient (not just flat noise)
            if drho_dr[idx_in_search] < 0:
                ablation_front_indices[t] = ablation_idx
                # Use outer boundary of that zone
                ablation_front_radius[t] = boundaries[ablation_idx + 1]

        # Smooth the profile
        ablation_front_radius = self._smooth_ablation_front(ablation_front_radius)

        # Store results
        self.data.ablation_front_radius = ablation_front_radius
        self.data.ablation_front_indices = ablation_front_indices

        valid_radii = ablation_front_radius[ablation_front_radius > 0]
        if len(valid_radii) > 0:
            logger.info(f"Ablation front tracked: {len(valid_radii)}/{n_times} timesteps")
            logger.info(f"Initial radius: {valid_radii[0]:.4f} cm, "
                        f"min radius: {np.min(valid_radii):.4f} cm")
        else:
            logger.warning("No valid ablation front positions found")
    
    def _compute_ifar(self):
        """
        Compute In-Flight Aspect Ratio (IFAR) at peak implosion velocity.

        IFAR = R_shell / Delta_R  where:
          R_inner = zone boundary at hot-spot/fuel interface (ri[t, 0])
          R_outer = ablation front radius (from _track_ablation_front)
          Delta_R = R_outer - R_inner  (shell thickness)
          R_shell = (R_inner + R_outer) / 2  (mid-shell radius)

        Evaluated at peak inward velocity (pre-stagnation), same time as adiabat.
        Requires _track_ablation_front() to have been called first.
        """
        abl_front = self.data.ablation_front_radius
        if abl_front is None:
            logger.warning("Cannot compute IFAR: ablation front not tracked")
            return

        ri = self.data.region_interfaces_indices
        vel = self.data.velocity
        zbnd = self.data.zone_boundaries

        if ri is None or vel is None or zbnd is None:
            logger.warning("Cannot compute IFAR: missing required data")
            return

        try:
            n_zones = self.data.mass_density.shape[1]

            # Zone-center velocities
            if vel.shape[1] > n_zones:
                v_zone = 0.5 * (vel[:, :n_zones] + vel[:, 1:n_zones+1])
            else:
                v_zone = vel[:, :n_zones]

            # Find peak inward velocity in fuel region (pre-stagnation)
            if ri.shape[1] >= 2:
                fuel_start = int(ri[0, 0])
                fuel_end = int(ri[0, -2])
                if fuel_end > fuel_start:
                    min_v_per_time = np.min(v_zone[:, fuel_start:fuel_end], axis=1)
                else:
                    min_v_per_time = np.min(v_zone, axis=1)
            else:
                min_v_per_time = np.min(v_zone, axis=1)

            stag_t = self.data.stag_time if self.data.stag_time > 0 else self.data.time[-1]
            pre_stag = self.data.time < stag_t
            if np.any(pre_stag):
                masked_v = min_v_per_time.copy()
                masked_v[~pre_stag] = 0.0
                t_pv = np.argmin(masked_v)
            else:
                t_pv = np.argmin(min_v_per_time)

            # Shell boundaries at peak velocity
            hs_bnd = int(ri[t_pv, 0])              # hot-spot / fuel interface (node index)
            R_inner = zbnd[t_pv, hs_bnd]            # cm

            R_outer = abl_front[t_pv]                # cm (from ablation front tracker)

            if R_outer <= R_inner or R_inner <= 0:
                logger.warning(f"IFAR: invalid shell boundaries at t={self.data.time[t_pv]:.2f} ns "
                             f"(R_inner={R_inner:.5f}, R_outer={R_outer:.5f} cm)")
                return

            delta_R = R_outer - R_inner
            R_shell = 0.5 * (R_inner + R_outer)

            self.data.ifar = R_shell / delta_R
            logger.info(f"IFAR at peak v_imp (t={self.data.time[t_pv]:.2f} ns): "
                        f"{self.data.ifar:.1f}  "
                        f"(R_shell={R_shell*1e4:.1f} um, dR={delta_R*1e4:.1f} um)")

        except Exception as e:
            logger.warning(f"Could not compute IFAR: {e}")
            self.data.ifar = 0.0

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
        """Analyze stagnation: bang time, max compression, hot spot, implosion."""
        logger.info("Analyzing stagnation phase...")
        
        # Determine stagnation time (min hot-spot volume)
        self._find_stagnation_time()
        
        # Determine bang time (peak fusion power)
        self._find_bang_time()
        
        # Compute hot spot properties at stagnation
        self._compute_hot_spot_properties()
        
        # Compute areal densities
        self._compute_areal_densities()
        
        # Implosion diagnostics (velocity, shocks, adiabat, ablation front)
        # — needs bang_time and stag_time from above
        self.analyze_implosion_phase()
    
    def _find_stagnation_time(self):
        """
        Find stagnation time = minimum hot-spot volume (before alpha heating).

        For a burning capsule, peak density can occur AFTER significant alpha
        deposition.  True stagnation is when the implosion stalls:
          - If region interfaces available: minimise hot-spot outer radius
          - Fallback: minimise minimum zone boundary > 0 (innermost shell)

        Peak density is recorded separately as max_density (may differ in time).
        """
        if self.data.mass_density is None:
            logger.warning("Cannot find stagnation time: missing density data")
            return

        # Peak density (may not coincide with stagnation for igniting targets)
        peak_density_vs_time = np.max(self.data.mass_density, axis=1)
        peak_idx = np.argmax(peak_density_vs_time)
        self.data.max_density = peak_density_vs_time[peak_idx]

        # Stagnation: minimum hot-spot outer radius
        ri = self.data.region_interfaces_indices
        if ri is not None and self.data.zone_boundaries is not None:
            # Hot-spot outer boundary = zone_boundaries[:, hs_bnd]
            hs_bnd = ri[:, 0].astype(int)
            hs_radius = np.array([
                self.data.zone_boundaries[t, hs_bnd[t]]
                for t in range(len(self.data.time))
            ])
            stag_idx = np.argmin(hs_radius)
        else:
            # Fallback: peak density time
            stag_idx = peak_idx

        self.data.stag_time = self.data.time[stag_idx]

        if ri is not None:
            hs_r = self.data.zone_boundaries[stag_idx, int(ri[stag_idx, 0])]
            logger.info(f"Stagnation time: {self.data.stag_time:.3f} ns  "
                        f"(min HS radius = {hs_r:.4f} cm)")
        else:
            logger.info(f"Stagnation time: {self.data.stag_time:.3f} ns  (from peak density)")
        logger.info(f"Peak density: {self.data.max_density:.2f} g/cc  "
                    f"(at t = {self.data.time[peak_idx]:.3f} ns)")
    
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
        """
        Compute hot spot properties at stagnation (minimum HS volume).

        All quantities here use stag_time, which is the moment the
        implosion stalls — BEFORE significant alpha-heating runaway.
        """
        if self.data.stag_time == 0.0:
            logger.warning("Cannot compute hot spot properties: stagnation time not determined")
            return
        
        # Find timestep closest to stagnation
        stag_idx = np.argmin(np.abs(self.data.time - self.data.stag_time))
        logger.info(f"Computing hot spot properties at stagnation "
                    f"(t = {self.data.time[stag_idx]:.3f} ns, index {stag_idx})")
        
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
                
                # Core radius = outer boundary of the outermost hot-spot zone
                hot_zone_indices = np.where(hot_spot_mask)[0]
                self.data.core_radius = boundaries[hot_zone_indices[-1] + 1]
                
                # Hot spot pressure (mass-averaged)
                # CLAUDE.md: report pressure in Gbar.  1 Gbar = 1e8 J/cm³
                pressure_Gbar = (self.data.ion_pressure[stag_idx] + 
                           self.data.rad_pressure[stag_idx]) * 1e-8  # J/cm³ → Gbar
                mass = self.data.zone_mass[stag_idx]
                
                self.data.hot_spot_pressure = np.average(
                    pressure_Gbar[hot_spot_mask],
                    weights=mass[hot_spot_mask]
                )
                
                # Hot spot areal density
                density = self.data.mass_density[stag_idx]
                dr = boundaries[1:] - boundaries[:-1]
                self.data.hot_spot_areal_density = np.sum(
                    density[hot_spot_mask] * dr[hot_spot_mask]
                )
                
                # Hot spot internal energy  (ion + electron)
                # specific_internal_energy × mass, summed over hot-spot zones,
                # then converted J → kJ.
                hs_energy_J = 0.0
                if self.data.ion_pressure is not None:
                    # E_int = P / (γ-1) × Volume  for ideal gas
                    # Or directly: e_specific × mass  if we have SIE.
                    # Use  E = (3/2) n k T × V  ≈ (3/2) P V  for each species
                    vol = (4.0 / 3.0) * np.pi * (boundaries[1:]**3 - boundaries[:-1]**3)
                    p_ion  = self.data.ion_pressure[stag_idx]
                    p_elec = self.data.rad_pressure[stag_idx]  # includes elec + rad
                    hs_energy_J = np.sum((1.5 * (p_ion[hot_spot_mask] + p_elec[hot_spot_mask]))
                                         * vol[hot_spot_mask])
                self.data.hot_spot_internal_energy = hs_energy_J * 1e-3  # J → kJ
                
                logger.info(f"Core radius: {self.data.core_radius:.4f} cm")
                logger.info(f"Hot spot radius: {self.data.stagnation_hot_spot_radius:.4f} cm")
                logger.info(f"Hot spot pressure: {self.data.hot_spot_pressure:.2f} Gbar")
                logger.info(f"Hot spot internal energy: {self.data.hot_spot_internal_energy:.2f} kJ")
                
        except Exception as e:
            logger.warning(f"Could not compute hot spot properties: {e}")
    
    def _compute_areal_densities(self):
        """
        Compute areal density (ρR) for different regions vs time.

        ICF convention: ρR = ∫ ρ(r) dr  (line integral of density).

        We store a cumulative ρR array so that ρR from zone i to zone j
        is simply  cumulative[j+1] - cumulative[i].

        Region mapping (generalised):
          - Hot spot          = innermost region  (zones 0 .. hs_bnd-1)
          - Cold fuel         = middle regions    (zones hs_bnd .. fuel_bnd-1)
          - Ablator           = outermost region  (zones fuel_bnd .. end)
          - "Fuel" (total DT) = hot spot + cold fuel
        """
        logger.info("Computing areal densities...")

        if self.data.zone_boundaries is None or self.data.mass_density is None:
            logger.warning("Cannot compute areal densities: missing boundary or density data")
            return

        n_times, n_zones = self.data.mass_density.shape

        # Cumulative ρR array:  cumulative[t, i] = ∫_0^{r_i} ρ dr
        # Shape (n_times, n_zones+1) — index 0 is at the centre (=0), index n_zones at the outer edge.
        cumulative = np.zeros((n_times, n_zones + 1))
        for t in range(n_times):
            dr = np.diff(self.data.zone_boundaries[t])          # Δr per zone
            rho_dr = self.data.mass_density[t] * dr              # ρ Δr per zone
            cumulative[t, 1:] = np.cumsum(rho_dr)

        self.data.areal_density_vs_time = cumulative             # store for neutron averaging

        # ---- Identify region boundaries ----
        ri = self.data.region_interfaces_indices                  # (n_times, n_regions)
        has_regions = ri is not None and ri.shape[1] >= 2

        # ---- Values at bang time ----
        if self.data.bang_time > 0:
            bt = np.argmin(np.abs(self.data.time - self.data.bang_time))

            # Total ρR (centre to outer edge)
            self.data.bang_time_areal_density = cumulative[bt, -1]

            if has_regions:
                hs_bnd   = int(ri[bt, 0])                        # hot-spot outer node
                fuel_bnd = int(ri[bt, -2])                        # fuel / ablator interface
                n_total  = n_zones                                # outer edge index

                self.data.bang_time_hs_areal_density   = cumulative[bt, hs_bnd]
                self.data.bang_time_fuel_areal_density  = cumulative[bt, fuel_bnd] - cumulative[bt, hs_bnd]
                self.data.bang_time_HDC_areal_density   = cumulative[bt, n_total]  - cumulative[bt, fuel_bnd]

                logger.info("Bang time areal densities:")
                logger.info(f"  Total:          {self.data.bang_time_areal_density:.4f} g/cm²")
                logger.info(f"  Hot spot:       {self.data.bang_time_hs_areal_density:.4f} g/cm²")
                logger.info(f"  Cold fuel:      {self.data.bang_time_fuel_areal_density:.4f} g/cm²")
                logger.info(f"  Ablator (CH):   {self.data.bang_time_HDC_areal_density:.4f} g/cm²")
            else:
                logger.info(f"Total bang-time ρR: {self.data.bang_time_areal_density:.4f} g/cm²")
    
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
        
        # Time-averaged ρR over burn period (needs burn_width from above)
        if (self.data.bang_time > 0 and self.data.burn_width > 0
                and self.data.areal_density_vs_time is not None):
            t_start = self.data.bang_time - self.data.burn_width / 2
            t_end   = self.data.bang_time + self.data.burn_width / 2
            mask = (self.data.time >= t_start) & (self.data.time <= t_end)
            if np.any(mask):
                self.data.time_ave_areal_density = np.mean(
                    self.data.areal_density_vs_time[mask, -1])
                logger.info(f"Time-averaged ρR (burn period): "
                            f"{self.data.time_ave_areal_density:.4f} g/cm²")
        
        # Maximum DT temperature
        if self.data.ion_temperature is not None:
            self.data.max_dt_temp = np.max(self.data.ion_temperature) / 1000.0  # eV → keV
            logger.info(f"Maximum DT temperature: {self.data.max_dt_temp:.2f} keV")

        # Burn propagation: hot-spot ρR (T_ion > 4.5 keV) and ignition time
        self._compute_burn_propagation()
    
    def _compute_burn_propagation(self):
        """
        Compute burn-propagation diagnostics following Olson et al.

        Hot-spot ρR is defined as the areal density of all zones with
        T_ion > 4.5 keV at each timestep.  This differs from the region-
        based hot-spot ρR (which uses region_interfaces_indices).

        Key outputs stored on self.data:
          - hot_spot_rhoR_vs_time  (n_times,)  g/cm²
          - total_rhoR_vs_time    (n_times,)  g/cm²
          - ignition_time         ns  — when HS ρR first crosses 0.3 g/cm²
          - ignition_hs_radius    cm  — HS radius at ignition
          - ignition_hs_pressure  Gbar — mass-avg HS pressure at ignition
          - burn_propagation_time ns  — when HS ρR ≈ total ρR
        """
        if (self.data.ion_temperature is None
                or self.data.mass_density is None
                or self.data.zone_boundaries is None):
            logger.warning("Cannot compute burn propagation: missing data")
            return

        logger.info("Computing burn propagation (T_ion > 4.5 keV hot-spot ρR)...")

        T_THRESHOLD_EV = 4500.0   # 4.5 keV in eV (ion_temperature stored in eV)
        RHO_R_IGNITION = 0.3      # g/cm² — ignition criterion

        n_times, n_zones = self.data.mass_density.shape
        hs_rhoR = np.zeros(n_times)
        total_rhoR = np.zeros(n_times)

        for t in range(n_times):
            dr = np.diff(self.data.zone_boundaries[t])              # Δr per zone
            rho_dr = self.data.mass_density[t] * dr                  # ρΔr per zone
            total_rhoR[t] = np.sum(rho_dr)

            hot_mask = self.data.ion_temperature[t] > T_THRESHOLD_EV
            hs_rhoR[t] = np.sum(rho_dr[hot_mask])

        self.data.hot_spot_rhoR_vs_time = hs_rhoR
        self.data.total_rhoR_vs_time = total_rhoR

        # ---- Ignition time: first crossing of ρR_hs = 0.3 g/cm² ----
        above = hs_rhoR >= RHO_R_IGNITION
        if np.any(above):
            ign_idx = np.argmax(above)      # first index where True

            # Linear interpolation for precise crossing time
            if ign_idx > 0:
                r0, r1 = hs_rhoR[ign_idx - 1], hs_rhoR[ign_idx]
                t0, t1 = self.data.time[ign_idx - 1], self.data.time[ign_idx]
                if r1 != r0:
                    frac = (RHO_R_IGNITION - r0) / (r1 - r0)
                    self.data.ignition_time = t0 + frac * (t1 - t0)
                else:
                    self.data.ignition_time = self.data.time[ign_idx]
            else:
                self.data.ignition_time = self.data.time[ign_idx]

            logger.info(f"Ignition time (ρR_hs ≥ 0.3 g/cm²): "
                        f"{self.data.ignition_time:.3f} ns")

            # ---- Hot-spot properties at ignition ----
            # Use the first timestep at or past ignition
            hot_mask = self.data.ion_temperature[ign_idx] > T_THRESHOLD_EV
            if np.any(hot_mask):
                boundaries = self.data.zone_boundaries[ign_idx]
                # HS radius: outer boundary of outermost hot zone
                hot_indices = np.where(hot_mask)[0]
                self.data.ignition_hs_radius = boundaries[hot_indices[-1] + 1]
                logger.info(f"Hot-spot radius at ignition: "
                            f"{self.data.ignition_hs_radius * 1e4:.1f} μm")

                # HS pressure: mass-averaged total pressure in Gbar
                p_total = (self.data.ion_pressure[ign_idx]
                           + self.data.rad_pressure[ign_idx]) * 1e-8  # Gbar
                mass = self.data.zone_mass[ign_idx]
                total_hs_mass = np.sum(mass[hot_mask])
                if total_hs_mass > 0:
                    self.data.ignition_hs_pressure = (
                        np.sum(p_total[hot_mask] * mass[hot_mask]) / total_hs_mass)
                    logger.info(f"Hot-spot pressure at ignition: "
                                f"{self.data.ignition_hs_pressure:.1f} Gbar")
        else:
            logger.info("No ignition detected (ρR_hs never reached 0.3 g/cm²)")

        # ---- Complete propagation time: when HS ρR ≈ total ρR ----
        # Defined as the first time after ignition when hs_rhoR >= 0.95 * total_rhoR
        if self.data.ignition_time > 0:
            post_ign = self.data.time >= self.data.ignition_time
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.where(total_rhoR > 0, hs_rhoR / total_rhoR, 0)
            propagated = post_ign & (ratio >= 0.95)
            if np.any(propagated):
                prop_idx = np.argmax(propagated)
                self.data.burn_propagation_time = self.data.time[prop_idx]
                logger.info(f"Complete burn propagation: "
                            f"{self.data.burn_propagation_time:.3f} ns")
            else:
                logger.info("Burn did not fully propagate through fuel")

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
            self.data.dt_neutron_yield = total_neutrons
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
            self.data.dt_neutron_yield = fusion_energy   # total reactions
            logger.info(f"Fusion energy output (time-integrated): "
                        f"{self.data.energy_output:.3f} MJ")
        else:
            self.data.energy_output = 0.0
            logger.warning("No fusion data available for yield calculation")
    
    def _compute_burn_width(self):
        """
        Compute burn width (FWHM of fusion rate pulse).

        Uses interpolation for sub-timestep accuracy, with RMS fallback
        if FWHM cannot be determined (mirrors corrected logic in core.py).
        """
        total_fusion = np.sum(self.data.fusion_power, axis=1)
        peak_power = np.max(total_fusion)

        if peak_power == 0:
            logger.warning("Peak fusion rate is zero — cannot compute burn width")
            return

        half_max = peak_power / 2.0
        above_half = total_fusion >= half_max

        if np.sum(above_half) >= 2:
            crossings = np.where(above_half)[0]

            # Interpolate left crossing
            if crossings[0] > 0:
                idx_l = crossings[0] - 1
                t_left = np.interp(
                    half_max,
                    [total_fusion[idx_l], total_fusion[idx_l + 1]],
                    [self.data.time[idx_l], self.data.time[idx_l + 1]]
                )
            else:
                t_left = self.data.time[crossings[0]]

            # Interpolate right crossing
            if crossings[-1] < len(self.data.time) - 1:
                idx_r = crossings[-1]
                t_right = np.interp(
                    half_max,
                    [total_fusion[idx_r + 1], total_fusion[idx_r]],
                    [self.data.time[idx_r + 1], self.data.time[idx_r]]
                )
            else:
                t_right = self.data.time[crossings[-1]]

            fwhm = t_right - t_left

            if fwhm > 0:
                self.data.burn_width = fwhm
                logger.info(f"Burn width (FWHM, interpolated): {fwhm:.4f} ns")
                return

        # ------- RMS fallback -------
        time = self.data.time
        rate = total_fusion
        total_rate = np.sum(rate)
        if total_rate > 0:
            t_mean = np.sum(time * rate) / total_rate
            t_sq_mean = np.sum(time**2 * rate) / total_rate
            rms_width = np.sqrt(t_sq_mean - t_mean**2)
            fwhm_from_rms = 2.355 * rms_width          # Gaussian relation

            self.data.burn_width = fwhm_from_rms
            logger.info(f"Burn width (RMS × 2.355 fallback): {fwhm_from_rms:.4f} ns")
        else:
            logger.warning("Cannot compute burn width — zero total fusion rate")
    
    def _compute_neutron_averaged_quantities(self):
        """
        Compute neutron-weighted average quantities.

        These are averaged over time AND space, weighted by the DT fusion
        reaction rate in each zone.  Gives the conditions "seen" by the
        average neutron produced during the implosion.

        FusionRate_DT_nHe4 from Helios is already a reaction rate
        (reactions/s per zone, volume-integrated), NOT power in Watts.
        Do NOT divide by energy_per_reaction.

        Proper neutron-weighted average of quantity Q:
            <Q>_n  =  Σ_t Σ_z  Q[t,z] · R[t,z] · Δt
                      ─────────────────────────────────
                       Σ_t Σ_z  R[t,z] · Δt

        where R[t,z] is the zone-level DT fusion rate.
        """
        logger.info("Computing neutron-averaged quantities...")

        if self.data.fusion_power is None:
            logger.warning("No fusion rate data available for neutron averaging")
            return

        # fusion_power is FusionRate_DT_nHe4 — already reactions/s per zone
        fusion_rate = self.data.fusion_power                     # (n_times, n_zones)
        time_seconds = self.data.time * 1e-9                     # ns → s
        dt_array = np.diff(time_seconds)                         # (n_times - 1,)

        # Total neutron yield  (double sum: zones × time)
        neutron_yield = 0.0
        for t in range(len(dt_array)):
            neutron_yield += np.sum(fusion_rate[t]) * dt_array[t]

        if neutron_yield == 0:
            logger.warning("Fusion rate is zero — cannot compute neutron-averaged quantities")
            return

        logger.info(f"Total neutron yield: {neutron_yield:.3e} reactions")

        # ---- Neutron-averaged FUEL areal density ----
        # Cold-fuel ρR = cumulative[fuel_bnd] - cumulative[hs_bnd]
        # This is what DT neutrons born in the hot spot traverse.
        if (self.data.region_interfaces_indices is not None and
                self.data.areal_density_vs_time is not None):
            ri = self.data.region_interfaces_indices
            cumul = self.data.areal_density_vs_time

            areal_sum = 0.0
            for t in range(len(dt_array)):
                hs_bnd   = int(ri[t, 0])
                fuel_bnd = int(ri[t, -2])
                fuel_rhoR = cumul[t, fuel_bnd] - cumul[t, hs_bnd]
                areal_sum += fuel_rhoR * np.sum(fusion_rate[t]) * dt_array[t]

            self.data.neutron_ave_fuel_areal_density = areal_sum / neutron_yield
            logger.info(f"Neutron-averaged fuel areal density: "
                        f"{self.data.neutron_ave_fuel_areal_density:.4f} g/cm²")
        else:
            logger.warning("Region interfaces not identified — cannot compute fuel areal density")

        # ---- Neutron-averaged ion temperature ----
        if self.data.ion_temperature is not None:
            n_times = len(dt_array)
            n_zones = min(self.data.ion_temperature.shape[1], fusion_rate.shape[1])

            temp_sum = 0.0
            for t in range(n_times):
                dt = dt_array[t]
                for z in range(n_zones):
                    temp_sum += self.data.ion_temperature[t, z] * fusion_rate[t, z] * dt

            self.data.neutron_ave_ion_temperature = temp_sum / neutron_yield
            # Convert eV → keV for reporting (CLAUDE.md convention)
            self.data.neutron_ave_ion_temperature /= 1000.0
            logger.info(f"Neutron-averaged ion temperature: "
                        f"{self.data.neutron_ave_ion_temperature:.2f} keV")
        else:
            logger.warning("Ion temperature not available for neutron averaging")

        # ---- Neutron-averaged pressure ----
        if self.data.ion_pressure is not None:
            n_times = len(dt_array)
            n_zones = min(self.data.ion_pressure.shape[1], fusion_rate.shape[1])

            pressure_sum = 0.0
            for t in range(n_times):
                dt = dt_array[t]
                for z in range(n_zones):
                    # Include electron (rad) pressure for total pressure
                    p_total = self.data.ion_pressure[t, z]
                    if self.data.rad_pressure is not None:
                        p_total += self.data.rad_pressure[t, z]
                    pressure_sum += p_total * fusion_rate[t, z] * dt

            neutron_ave_pressure = pressure_sum / neutron_yield
            # Convert: Helios J/cm³ → Gbar  (1 Gbar = 1e8 J/cm³)
            self.data.neutron_ave_pressure = neutron_ave_pressure * 1e-8
            logger.info(f"Neutron-averaged pressure: "
                        f"{self.data.neutron_ave_pressure:.2f} Gbar")
        else:
            logger.warning("Ion pressure not available for neutron averaging")
    
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

        # ---- Mass fractions (require region interfaces + ablation front) ----
        ri = self.data.region_interfaces_indices
        if (ri is not None and ri.shape[1] >= 2
                and self.data.zone_mass is not None
                and self.data.zone_boundaries is not None):

            stag_idx = np.argmin(np.abs(self.data.time - self.data.stag_time)) \
                       if self.data.stag_time > 0 else 0

            n_zones = self.data.zone_mass.shape[1]
            fuel_bnd = int(ri[0, -2])      # fuel / ablator interface (Lagrangian)
            hs_bnd   = int(ri[stag_idx, 0])  # hot-spot outer boundary at stagnation

            # Initial masses
            initial_fuel_mass    = np.sum(self.data.zone_mass[0, :fuel_bnd])
            initial_ablator_mass = np.sum(self.data.zone_mass[0, fuel_bnd:])

            # ---- Determine "inside the shell" mask at stagnation ----
            #
            # In a Lagrangian code the zone INDEX is a mass coordinate that
            # never changes.  Zones with index <= ablation-front index are
            # by definition unablated (inside the shell).
            #
            # Previous approach used a *radius* comparison
            #   inside_shell = zone_outer <= abl_r
            # which fails when the smoothed ablation-front radius is
            # inflated by plateau-averaging with large pre-stagnation values,
            # causing every zone to be classified as "inside."
            #
            # Primary method: use the raw (unsmoothed) ablation-front ZONE
            # INDEX.  Fallback: density-based scan outward from hot spot.

            abl_idx = None   # ablation-front zone index at stagnation

            # Method 1 — tracked ablation-front zone index
            if (self.data.ablation_front_indices is not None
                    and self.data.ablation_front_indices[stag_idx] > 0):
                abl_idx = int(self.data.ablation_front_indices[stag_idx])

            # Method 2 — density fallback: outermost zone (outside HS)
            # whose density exceeds 1 g/cc is the last unablated zone.
            if abl_idx is None or abl_idx <= hs_bnd:
                rho_stag = self.data.mass_density[stag_idx]
                abl_idx = hs_bnd  # default: nothing outside HS is unablated
                for z in range(n_zones - 1, hs_bnd, -1):
                    if rho_stag[z] > 1.0:
                        abl_idx = z
                        break

            # Build the mask using Lagrangian zone index
            zone_indices = np.arange(n_zones)
            inside_shell = zone_indices <= abl_idx

            # Also store the ablation-front radius for diagnostics / plotting
            bnd_stag = self.data.zone_boundaries[stag_idx]
            abl_r = bnd_stag[min(abl_idx + 1, len(bnd_stag) - 1)]

            # ---- Unablated ablator fraction ----
            # Ablator zones (fuel_bnd .. end) that remain inside the shell
            ablator_mask = np.zeros(n_zones, dtype=bool)
            ablator_mask[fuel_bnd:] = True
            unablated_abl_mass = np.sum(
                self.data.zone_mass[stag_idx, ablator_mask & inside_shell])
            if initial_ablator_mass > 0:
                self.data.unablated_ablatar_mass = unablated_abl_mass / initial_ablator_mass

            # ---- Unablated fuel fraction ----
            # Fuel zones (0 .. fuel_bnd-1) that are inside the shell
            fuel_mask = np.zeros(n_zones, dtype=bool)
            fuel_mask[:fuel_bnd] = True
            unablated_fuel_mass = np.sum(
                self.data.zone_mass[stag_idx, fuel_mask & inside_shell])
            if initial_fuel_mass > 0:
                self.data.unablated_fuel_mass = unablated_fuel_mass / initial_fuel_mass

            # ---- Stagnated fuel fraction ----
            # Dense fuel assembled around the hot spot at stagnation — the
            # confinement shell.  Defined as fuel zones between the hot-spot
            # boundary (hs_bnd) and the ablation front (abl_idx).
            #
            # NOTE: no temperature filter.  For igniting capsules, alpha
            # heating warms the dense shell above 1 keV well before peak
            # burn; a temperature cut would incorrectly reject the entire
            # confinement shell.  The hot spot is already excluded by the
            # outside_hs mask.
            outside_hs = zone_indices >= hs_bnd

            stagnated_mask = fuel_mask & inside_shell & outside_hs
            stag_fuel_mass = np.sum(self.data.zone_mass[stag_idx, stagnated_mask])
            if initial_fuel_mass > 0:
                self.data.stagnated_fuel_mass = stag_fuel_mass / initial_fuel_mass

            # ---- Diagnostics ----
            n_inside   = int(np.sum(inside_shell))
            n_abl_in   = int(np.sum(ablator_mask & inside_shell))
            n_fuel_in  = int(np.sum(fuel_mask & inside_shell))
            logger.info(f"Ablation front at stagnation: zone {abl_idx}, "
                        f"r = {abl_r:.4f} cm  "
                        f"(HS bnd = zone {hs_bnd}, fuel/abl bnd = zone {fuel_bnd})")
            logger.info(f"Inside shell: {n_inside}/{n_zones} zones  "
                        f"(fuel: {n_fuel_in}/{fuel_bnd}, "
                        f"ablator: {n_abl_in}/{n_zones - fuel_bnd})")
            logger.info(f"Mass fractions — unablated fuel: {self.data.unablated_fuel_mass:.4f}, "
                        f"unablated ablator: {self.data.unablated_ablatar_mass:.4f}, "
                        f"stagnated fuel: {self.data.stagnated_fuel_mass:.4f}")
