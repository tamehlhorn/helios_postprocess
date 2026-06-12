"""
ICF Data Builder — Bridge from HeliosRun to Legacy Pipeline
============================================================
Bulk-loads all time steps from a HeliosRun object into an ICFRunData
container whose attribute names match the legacy ICFAnalyzer, ICFPlotter,
and OutputGenerator interfaces.

Usage
-----
    from helios_postprocess import HeliosRun
    from helios_postprocess.data_builder import build_run_data

    run = HeliosRun('simulation.exo', verbose=True)
    data = build_run_data(run)

    # data is now ready for:
    #   ICFAnalyzer(data).analyze_drive_phase()
    #   ICFPlotter(data, config).create_full_report('report.pdf')
    #   OutputGenerator(data, config).save_metrics_json('metrics.json')

Author: Prof T / Xcimer ICF Analysis
Date: 2025
"""

import numpy as np
import warnings
import logging
from pathlib import Path
from typing import Optional, Dict

logger = logging.getLogger(__name__)


# ============================================================================
# ICFRunData Container
# ============================================================================

class ICFRunData:
    """
    Data container matching the attribute interface expected by ICFAnalyzer,
    ICFPlotter, and OutputGenerator.

    Attributes are grouped into:
      1. 2-D arrays  — shape (n_times, n_zones) or (n_times, n_nodes)
      2. 1-D arrays  — shape (n_times,) or (n_zones,)
      3. Scalar metrics — floats initialised to 0.0 (set by ICFAnalyzer)
      4. Optional arrays — set to None when variable is absent from the EXODUS file
      5. Metadata — filename, rhw_config, etc.
    """

    def __init__(self):
        # ------------------------------------------------------------------
        # 1. Core 2-D arrays (loaded by build_run_data)
        # ------------------------------------------------------------------
        self.time: Optional[np.ndarray] = None               # (n_times,)  in nanoseconds
        self.mass_density: Optional[np.ndarray] = None        # (n_times, n_zones) g/cm³
        self.ion_temperature: Optional[np.ndarray] = None     # (n_times, n_zones) eV
        self.elec_temperature: Optional[np.ndarray] = None    # (n_times, n_zones) eV
        self.ion_pressure: Optional[np.ndarray] = None        # (n_times, n_zones) J/cm³
        self.elec_pressure: Optional[np.ndarray] = None       # (n_times, n_zones) J/cm³
        self.plasma_pressure: Optional[np.ndarray] = None     # (n_times, n_zones) J/cm³  ion + electron
        self.rad_pressure_true: Optional[np.ndarray] = None   # (n_times, n_zones) J/cm³  radiation only
        self.rad_pressure: Optional[np.ndarray] = None        # (n_times, n_zones) J/cm³
        # NOTE: rad_pressure stores the "non-ion" component so that
        # total_pressure = ion_pressure + rad_pressure works everywhere.
        # If the EXODUS file has a true rad_pressure variable, this is
        # elec_pressure + rad_pressure.  Otherwise it's elec_pressure alone.
        self.velocity: Optional[np.ndarray] = None            # (n_times, n_nodes) cm/s — node-centered
        self.zone_boundaries: Optional[np.ndarray] = None     # (n_times, n_nodes) cm — node-centered
        self.zone_mass: Optional[np.ndarray] = None           # (n_times, n_zones) g
        self.zone_centers: Optional[np.ndarray] = None        # (n_times, n_zones) cm — derived

        # ------------------------------------------------------------------
        # 2. Optional physics arrays (None when absent)
        # ------------------------------------------------------------------
        self.fusion_power: Optional[np.ndarray] = None            # (n_times, n_zones) — zone-level DT fusion rate
        self.laser_energy_deposited: Optional[np.ndarray] = None  # (n_times, n_zones) J  (time-integrated)
        self.laser_energy_delivered_cum: Optional[np.ndarray] = None  # (n_times,) J  Helios cumulative integral of LaserEnDeliveredTimeInt
        self.laser_power_source: Optional[np.ndarray] = None      # (n_times, n_zones) W/cm³
        self.laser_power_delivered: Optional[np.ndarray] = None   # (n_times,) W — total laser power on target
        self.laser_power_delivered_per_beam: Optional[np.ndarray] = None  # (n_times, n_beam) W — preserved for multi-beam
        self.electron_density: Optional[np.ndarray] = None        # (n_times, n_zones) cm⁻³
        self.ion_density: Optional[np.ndarray] = None             # (n_times, n_zones) cm⁻³
        self.mean_charge: Optional[np.ndarray] = None              # (n_times, n_zones) Z̄
        # Laser intensity reconstruction inputs (raw from EXODUS; cleaned in laser_intensity.py)
        self.laser_attenuation_coeff: Optional[np.ndarray] = None  # (n_times, n_beam, n_bnd) 1/cm
        self.laser_power_on_target: Optional[np.ndarray] = None           # (n_times,) W — total summed across beams
        self.laser_power_on_target_per_beam: Optional[np.ndarray] = None  # (n_times, n_beam) W — preserved for multi-beam
        self.neutron_production_rate: Optional[np.ndarray] = None # (n_times, n_zones)
        self.alpha_heating_power: Optional[np.ndarray] = None     # (n_times, n_zones)
        self.alpha_heating_ion: Optional[np.ndarray] = None       # (n_times, n_zones) particle → ion
        self.alpha_heating_ele: Optional[np.ndarray] = None       # (n_times, n_zones) particle → elec
        self.laser_wavelength_um: float = 0.0
        self.laser_geometry_per_beam: Optional[list] = None
        self.laser_spot_size_cm: float = 0.0
        self.laser_half_cone_angle_deg: float = 0.0
        self.laser_focus_position_cm: float = 0.0
        self.laser_power_multiplier: float = 1.0
        self.laser_spatial_profile: str = "Uniform"
        self.laser_foot_power_TW: float = 0.0
        self.laser_peak_power_TW: float = 0.0
        self.laser_foot_start_ns: float = 0.0
        self.laser_foot_end_ns: float = 0.0
        self.laser_peak_start_ns: float = 0.0
        self.laser_peak_end_ns: float = 0.0
        self.laser_pulse_duration_ns: float = 0.0
        self.flux_limiter: float = 0.0              # from .rhw (first-region value, legacy)
        self.flux_limiter_enabled: bool = False     # from .rhw (any region enabled)
        self.flux_limiter_per_region: Optional[list] = None  # [{region, enabled, value}, ...]
        # Active-source counts (May 2026). RHW reports max-array dims for
        # beams/IDD even when only one is driven; these scalars expose what's
        # *actually* on. Set by build_run_data + rhw_config propagation.
        self.n_total_beams:   int = 0
        self.n_active_beams:  int = 0
        self.rad_source_rmin_on: bool = False
        self.rad_source_rmax_on: bool = False
        self.n_active_idd_sources: int = 0
        self.alpha_deposition_local: bool = False
        self.alpha_deposition_nonlocal: bool = False
        self.eos_models: list = None                              # [{region, type, file}]
        self.dt_neutron_count: Optional[np.ndarray] = None        # time-integrated DT neutron yield
        self.dd_neutron_count: Optional[np.ndarray] = None        # time-integrated DD neutron yield
        self.dt_neutron_count_zone: Optional[np.ndarray] = None   # (n_times, n_zones) zone-resolved
        # Boundary-tally cumulative quantities — direct EXODUS measurements
        # used for energy-balance closure (replace inferred fractions).
        self.radiation_energy_at_boundary_cum: Optional[np.ndarray] = None  # (n_times,) J  cumulative rad through grid boundaries
        self.particle_energy_escaped_cum: Optional[np.ndarray] = None       # (n_times,) J  cumulative particle (neutron+) escape
        self.scale_length: Optional[np.ndarray] = None            # (n_times, n_zones) cm — derived
        self.region_interfaces_indices: Optional[np.ndarray] = None  # (n_times, n_regions)
        self.material_index: Optional[np.ndarray] = None          # (n_zones,) material ID per zone
        self.region_names: Optional[list] = None                   # list of region name strings
        # ----- Target classification (for capsule vs halfraum-capsule analysis routing) -----
        # `target_class` = "capsule" | "halfraum_capsule"
        #   capsule:           grid contains only the imploding payload (gas + fuel + ablator)
        #   halfraum_capsule:  grid contains capsule + external structure (gas gap + converter shell)
        # `n_capsule_regions` = number of regions belonging to the capsule itself.
        #   For "capsule":           equals the total number of regions
        #   For "halfraum_capsule":  less than total; outermost (n_total - n_capsule) regions
        #                            are external structure and excluded from implosion analyses
        # These fields are auto-populated below from region_names; analyzers use the derived
        # `capsule_outer_idx` and `fuel_ablator_idx` properties to get the correct ri[:,k]
        # column for the capsule outer surface and fuel/ablator interface respectively.
        self.target_class: str = "capsule"
        self.n_capsule_regions: int = 0
        self.areal_density_vs_time: Optional[np.ndarray] = None  # (n_times, n_zones+1)

        # ------------------------------------------------------------------
        # 3. Scalar metrics — set to sensible defaults, populated by ICFAnalyzer
        # ------------------------------------------------------------------
        # Timing
        self.bang_time: float = 0.0              # ns
        self.stag_time: float = 0.0              # ns
        self.burn_width: float = 0.0             # ns  (FWHM)

        # Compression / density
        self.max_density: float = 0.0            # g/cm³
        self.comp_ratio: float = 0.0
        self.cr_inflight: float = 0.0        # R0 / R_ablfront at peak velocity
        self.core_radius: float = 0.0            # cm

        # Hot spot (at stagnation)
        self.stagnation_hot_spot_radius: float = 0.0  # cm
        self.hot_spot_pressure: float = 0.0            # Gbar
        self.hot_spot_areal_density: float = 0.0       # g/cm²
        self.hot_spot_internal_energy: float = 0.0     # kJ

        # Areal densities (ρR = ∫ρ dr)
        self.time_ave_areal_density: float = 0.0       # g/cm²
        self.bang_time_areal_density: float = 0.0      # g/cm² — total
        self.bang_time_fuel_areal_density: float = 0.0 # g/cm² — cold fuel (hs_bnd to fuel_bnd)
        self.bang_time_HDC_areal_density: float = 0.0  # g/cm² — ablator
        self.bang_time_hs_areal_density: float = 0.0   # g/cm² — hot spot

        # Burn propagation (hot-spot ρR defined by T_ion > 4.5 keV)
        self.hot_spot_rhoR_vs_time: Optional[np.ndarray] = None  # (n_times,) g/cm²
        self.total_rhoR_vs_time: Optional[np.ndarray] = None     # (n_times,) g/cm²
        self.ignition_time: float = 0.0            # ns  (when HS ρR crosses 0.3 g/cm²)
        self.ignition_hs_radius: float = 0.0       # cm  (at ignition time)
        self.ignition_hs_pressure: float = 0.0     # Gbar (at ignition time)
        self.ignition_T_ion_onaxis_keV: float = 0.0  # T_ion at r=0 at ignition;
                                                     # apples-to-apples with
                                                     # published-figure central T
        self.ignition_T_ion_hs_avg_keV: float = 0.0  # mass-averaged T_ion over
                                                     # T>4.5 keV zones at ignition;
                                                     # apples-to-apples with
                                                     # T_ion_hs_at_ignition_*_keV
                                                     # per-code published values
        self.ignition_index: int = -1                # timestep index for the
                                                     # radial-profile CSV writer
        self.peak_density_at_ignition: float = 0.0   # g/cm³ — max ρ at ignition;
                                                     # falls back to peak ρ at
                                                     # stagnation for no-burn runs
        self.peak_density_at_ignition_is_stagnation: bool = False  # True if the
                                                     # value above is the no-burn
                                                     # fallback (used by summary
                                                     # writer to label the row)
        self.burn_propagation_time: float = 0.0    # ns  (when HS ρR = total ρR)

        # Neutron-averaged quantities
        self.neutron_ave_fuel_areal_density: float = 0.0  # g/cm²
        self.neutron_ave_ion_temperature: float = 0.0     # keV
        self.neutron_ave_pressure: float = 0.0            # Gbar

        # Energy / performance
        self.energy_output: float = 0.0    # MJ  (fusion yield)
        self.dt_neutron_yield: float = 0.0 # total DT neutrons produced
        self.laser_energy: float = 0.0     # MJ
        # Note: data.laser_energy is currently max(absorbed) in MJ -- legacy
        # naming, see analyze_drive_phase. The two explicit scalars below
        # disambiguate when both sources are available.
        self.laser_energy_delivered_MJ: float = 0.0   # max of LaserEnDeliveredTimeInt
        self.laser_energy_absorbed_MJ:  float = 0.0   # max of EnLaserDepositedTimeIntg
        # Coupling diagnostics from analyze_drive_phase
        self.eff_avg_coupling_pct:  float = 0.0       # 100 * E_abs(end) / E_del(end)
        self.eff_peak_coupling_pct: float = 0.0       # max(E_abs(t)/E_del(t)) after IC
        self.rad_energy: float = 0.0       # MJ  (radiation / IDD deposited)
        self.target_gain: float = 0.0
        self.max_dt_temp: float = 0.0      # keV

        # Mass fractions
        self.unablated_ablatar_mass: float = 0.0
        self.unablated_fuel_mass: float = 0.0
        self.stagnated_fuel_mass: float = 0.0

        # Implosion
        self.peak_implosion_velocity: float = 0.0  # km/s -- legacy: peak inward
                                                   # zone velocity in shell, pre-bang.
                                                   # For high-CR HDD-class targets
                                                   # this is inflated by late-time
                                                   # shock-convergence spikes; prefer
                                                   # peak_implosion_velocity_at_cr15
                                                   # for cross-code comparison.
        self.peak_implosion_velocity_at_cr15: float = 0.0  # km/s -- mass-avg shell
                                                   # velocity at CR=1.5 snapshot.
                                                   # NOTE: empirically reports much
                                                   # earlier-than-peak velocities;
                                                   # prefer implosion_velocity_rhino
                                                   # for cross-tool comparison.
        self.t_peak_velocity_at_cr15_ns: float = 0.0  # ns
        self.implosion_velocity_rhino_kms: float = 0.0  # km/s -- RHINO/W.Trickey
                                                   # convention: shell defined as
                                                   # rho > rho_peak/e per timestep;
                                                   # v_shell = sqrt(2*KE_sh/m_sh);
                                                   # turning point (peak) of
                                                   # v_shell(t) pre-stagnation.
                                                   # Apples-to-apples with RHINO's
                                                   # reported "Implosion velocity".
        self.t_implosion_velocity_rhino_ns: float = 0.0  # ns
        self.adiabat_min_rhino: float = 0.0        # dimensionless -- RHINO/W.Trickey
                                                   # "Min shell adiabat (CR=1.5)".
                                                   # Trigger = first time the
                                                   # gas/cold-fuel Lagrangian
                                                   # boundary radius reaches
                                                   # R_inner(0)/1.5. Adiabat per
                                                   # shell zone (rho>peak/e in DT
                                                   # cold-fuel region) using
                                                   # standard DT Fermi formula,
                                                   # min over those zones.
        self.t_adiabat_min_rhino_ns: float = 0.0   # ns
        self.r_inner_initial_cm: float = 0.0       # cm (audit)
        self.r_inner_at_cr15_cm: float = 0.0       # cm (audit)
        self.t_breakout_rhino_ns: float = 0.0      # ns -- RHINO breakout-time
                                                   # audit value (peak density
                                                   # crosses 1%-threshold inner
                                                   # shell surface at t=0)
        # ── RHINO diagnostics suite (June 2026) ──────────────────
        # Seven additional RHINO-convention metrics from
        # ICFAnalyzer._compute_rhino_diagnostics. See that method for
        # the algorithmic definitions; all match W. Trickey's RHINO
        # postprocessor (private repo wtrickey27/RHINO) conventions.
        self.shell_inner_pos_history_cm: Optional[np.ndarray] = None
        self.shell_outer_pos_history_cm: Optional[np.ndarray] = None
        self.shell_thickness_history_cm: Optional[np.ndarray] = None
        self.shell_mass_history_mg:      Optional[np.ndarray] = None
        self.shell_velocity_history_kms: Optional[np.ndarray] = None
        self.cr_inner_history:           Optional[np.ndarray] = None
        self.t_max_shell_velocity_rhino_ns: float = 0.0
        self.stag_time_rhino_ns: float = 0.0    # RHINO stagnation_time
                                                # (shell-velocity minimum;
                                                # different from this
                                                # pipeline's stag_time
                                                # which is HS radius min)
        self.assembled_mass_rhino_mg: float = 0.0
        self.burn_fraction_rhino: float = 0.0
        self.ablation_pressure_at_cr_3p5_Mbar: float = 0.0  # Vulcan HDD
                                                # design convention
        self.t_at_cr_3p5_ns: float = 0.0       # audit

        # ── Will list extensions (June 4 2026) ───────────────────
        # Items from W. Trickey's variable list (June 2026) that were
        # not previously named scalars. Definitions independent of the
        # five open convention questions ("shell" / "hotspot" /
        # hydro-efficiency denominator / stag-time / bang-time-burn-off);
        # the convention-dependent items wait on Will's response.
        self.cr_outer: float = 0.0             # outer convergence ratio
                                               # R_outer(t=0) / R_outer(t_stag)
                                               # using shell_outer_pos_history_cm
                                               # (rho > rho_peak/e). RHINO equiv:
                                               # derivable from sim.shell_outer.
        self.laser_overlapped_intensity_Wcm2: float = 0.0  # peak power /
                                               # (4π R0²), where R0 is the
                                               # initial capsule outer radius.
                                               # Will's "overlapped intensity".
        self.shell_mass_at_stagnation_mg: float = 0.0     # shell_mass_history
                                               # evaluated at stag_time_rhino_ns
        self.hot_spot_mass_at_stagnation_mg: float = 0.0  # Σ zone_mass[stag,
                                               # 0:ri[stag, 0]] -- mass inside
                                               # gas/cold-fuel boundary at stag
        self.stag_time_areal_density: float = 0.0        # g/cm² -- total ρR
                                               # at stagnation (parallel to
                                               # bang_time_areal_density)
        self.stag_time_fuel_areal_density: float = 0.0   # g/cm² -- cold-fuel ρR
                                               # at stagnation
        self.neutron_ave_electron_temperature: float = 0.0  # keV -- mirror of
                                               # neutron_ave_ion_temperature
                                               # using elec_temperature in the
                                               # same neutron-rate-weighted
                                               # average

        # ── W. Trickey June 2026 shell convention (Phys. Conv. §18b) ──
        # Inner = 1% pre-breakout / (1/e) post-breakout of peak density.
        # Outer = first inflection point in ρ(r) going outward from ρ_peak
        # (d²ρ/dr² = 0). Replaces Lagrangian gas/foam interface because foam
        # ablates inward during shock breakout, making the Lagrangian
        # boundary a poor measure of the dense shell.
        self.shell_inner_will_history_cm:  Optional[np.ndarray] = None  # (n_times,) cm
        self.shell_outer_will_history_cm:  Optional[np.ndarray] = None  # (n_times,) cm
        self.shell_thickness_will_history_cm: Optional[np.ndarray] = None
        self.shell_mass_will_history_mg:   Optional[np.ndarray] = None  # (n_times,) mg
        self.shell_mass_will_at_stagnation_mg: float = 0.0  # at stag_time_rhino_ns
        self.adiabat_mass_avg_will_cr15: float = 0.0  # mass-avg adiabat
                                               # (RHINO partially_ionized
                                               # formula) over Will-shell at
                                               # cr_inner = 1.5
        self.sound_speed_shell_will_cr15_kms: float = 0.0  # mass-avg
                                               # c_s = sqrt(γP/ρ), γ=5/3
                                               # (ideal-gas approx for shell
                                               # conditions), over Will-shell
                                               # at cr_inner = 1.5
        self.t_will_shell_cr15_ns: float = 0.0  # audit: timestep at which
                                               # cr_inner=1.5 was found for
                                               # the will-shell metrics above

        # ── W. Trickey ignition-product timing (Phys. Conv. §18d/e) ──
        # Peak hotspot quantities are evaluated at the time the ignition
        # PRODUCT peaks, not at a fixed kinematic landmark. Two products:
        #   Lawson:      rhoR_hs(t) × T_i_hs(t)    (volume-avg T_i)
        #   Confinement: P_hs(t)    × R_hs(t)      (volume-avg P_hs)
        # The HS region is the Lagrangian gas/cold-fuel interior
        # (zones [0, ri[t, 0])) per existing convention.
        self.t_peak_rhoR_Ti_ns: float = 0.0
        self.peak_rhoR_Ti_gcm2_keV: float = 0.0       # the product value
        self.rhoR_hs_at_peak_rhoR_Ti_gcm2: float = 0.0
        self.T_i_volume_avg_at_peak_rhoR_Ti_keV: float = 0.0
        self.T_i_mass_avg_at_peak_rhoR_Ti_keV: float = 0.0

        self.t_peak_Phs_Rhs_ns: float = 0.0
        self.peak_Phs_Rhs_Gbar_um: float = 0.0        # the product value
        self.R_hs_at_peak_Phs_Rhs_um: float = 0.0
        self.P_hs_volume_avg_at_peak_Phs_Rhs_Gbar: float = 0.0
        self.P_hs_mass_avg_at_peak_Phs_Rhs_Gbar: float = 0.0
        # Adiabat using RHINO's fully_ionized_dt convention (n_e = ρ/m_avg_ion
        # instead of actual electron_density). This is RHINO's default --
        # matches Will Trickey's RHINO native output. Pure DT zones give
        # identical results to partially_ionized; multi-material foam zones
        # diverge because fully_ionized_dt assumes Z̄=1 (under-counts e- in
        # CD foam).
        self.adiabat_min_rhino_fully_ionized: float = 0.0
        self.adiabat_mass_avg_rhino_fully_ionized: float = 0.0
        self.ifar: float = 0.0                      # in-flight aspect ratio at peak v_imp
        self.adiabat_mass_averaged_ice: float = 0.0  # legacy: at peak velocity
                                                     # Lindl convention
                                                     # (potentially shock-inflated
                                                     # on high-CR targets)
        self.adiabat_mass_averaged_ice_cr15: float = 0.0  # Thomas/RHINO convention:
                                                     # at CR=1.5, before late-time
                                                     # shock pre-heating
                                                     # Lindl convention
        self.adiabat_mass_averaged_ice_rhino_formula: float = 0.0  # same selection
                                                     # as adiabat_mass_averaged_ice
                                                     # but using proper degenerate
                                                     # electron gas Fermi pressure
                                                     # (RHINO partially_ionized).
                                                     # Typically ~14x larger than
                                                     # the Lindl value at ICF
                                                     # densities.
        self.adiabat_mass_averaged_ice_cr15_rhino_formula: float = 0.0  # at CR=1.5
                                                     # mass-avg with RHINO Fermi
                                                     # formula. Cross-tool
                                                     # comparable to publication
                                                     # references that use the
                                                     # physics Fermi (Thomas etc).
        self.adiabat_at_breakout_rhino_formula: float = 0.0  # base adiabat at
                                                     # shock breakout, RHINO
                                                     # Fermi formula.
        self.shock_breakout_time_ns: float = 0.0     # ns  (first shock exits DT ice into hot spot)
        self.shock_breakout_pressure_Gbar: float = 0.0  # Gbar (post-shock pressure at breakout)
        self.shock_breakout_mach: float = 0.0        # Mach number at breakout
        self.shock_foot_pressure_Gbar: float = 0.0   # Gbar (peak pressure in DT ice at foot end)

        # Critical density tracking (populated by ICFAnalyzer)
        self.critical_density_value: float = 0.0
        self.critical_density_radius_formula: Optional[np.ndarray] = None
        self.critical_density_radius_algorithm: Optional[np.ndarray] = None

        # Ablation front (populated by ICFAnalyzer)
        self.ablation_front_radius: Optional[np.ndarray] = None
        self.ablation_front_indices: Optional[np.ndarray] = None

        # Shock tracking (populated by ICFAnalyzer)
        self.shock_times: list = []
        self.shock_radii: list = []
        self.shock_velocities: list = []
        self.first_shock: Optional[Dict] = None  # populated by ICFAnalyzer._analyze_first_shock

        # Multi-shock train (populated by ICFAnalyzer._compute_shock_train)
        self.shock_trajectories: list = []
        self.shock_coalescence_events: list = []
        self.shock_breakouts: list = []
        self.shock_breakout_times_ns: np.ndarray = np.array([], dtype=float)
        # Consolidated foot/ramp/peak events + class-keyed scalars
        self.shock_events: list = []
        self.t_foot_shock_ns: float = float('nan')
        self.t_ramp_shock_ns: float = float('nan')
        self.t_peak_shock_ns: float = float('nan')

        # ------------------------------------------------------------------
        # 4. Metadata
        # ------------------------------------------------------------------
        self.filename: str = ""
        self.rhw_config = None          # RHWConfig dataclass (from rhw_parser)
        self.drive_temperature: Optional[np.ndarray] = None  # eV  (from .rhw file)
        self.drive_time: Optional[np.ndarray] = None         # seconds

    # ------------------------------------------------------------------
    # Helper properties for region indexing (used by ICFAnalyzer).
    #
    # For "capsule" targets (current default behaviour), these reduce to
    # the historical ri[:, -1] and ri[:, -2] respectively, so existing
    # call sites that switch to using these properties produce identical
    # results.
    #
    # For "halfraum_capsule" targets where external structure (gas gap +
    # converter shell) lives outside the capsule, these point to the
    # capsule outer surface and the fuel/ablator interface within the
    # capsule rather than the grid edge -- which is what implosion
    # diagnostics actually want.
    # ------------------------------------------------------------------
    @property
    def capsule_outer_idx(self) -> int:
        """Column index in region_interfaces_indices for the capsule outer surface.
        For standard capsules this is the last column (== grid outer); for
        halfraum-capsule targets it is the inner edge of the first external
        region (== outer surface of the ablator).
        """
        if self.region_interfaces_indices is None:
            return -1
        n_total = self.region_interfaces_indices.shape[1]
        if self.n_capsule_regions <= 0 or self.n_capsule_regions > n_total:
            return n_total - 1
        return self.n_capsule_regions - 1

    @property
    def fuel_ablator_idx(self) -> int:
        """Column index in region_interfaces_indices for the fuel/ablator interface.

        The interface index is the OUTER edge of the LAST DT-bearing region.
        Identified by scanning region_names for the first non-DT region (i.e. the
        first region whose name doesn't contain 'DT'). For targets with multiple
        ablator layers (e.g. WT_cthomas with 'CD ablator' + 'CD ablator_2'),
        this correctly returns the inner ablator boundary rather than mistakenly
        pointing to a boundary inside the ablator stack.

        Backwards-compatible with the legacy n_capsule_regions-2 formula for
        standard 4-region targets (Olson_PDD_20, VI_6) where there is exactly
        one ablator region.
        """
        if self.region_interfaces_indices is None:
            return -2
        n_total = self.region_interfaces_indices.shape[1]

        # Primary: scan region_names for the first region whose name does not
        # contain 'dt'. The interface immediately INNER to that region is the
        # fuel/ablator boundary (= outer edge of the last DT region).
        names = getattr(self, 'region_names', None)
        if names is not None and len(names) > 0:
            for i, name in enumerate(names):
                lower = str(name).lower()
                if 'dt' in lower:
                    continue
                # First non-DT region — interface before it
                if i > 0:
                    return i - 1
                return 0   # entirely ablator from index 0 (unusual)
            # All regions DT-bearing — no ablator; return outer-most region
            return max(n_total - 1, 0)

        # Fallback: legacy n_capsule_regions-2 formula
        if self.n_capsule_regions <= 1 or self.n_capsule_regions > n_total:
            return max(n_total - 2, 0)
        return self.n_capsule_regions - 2

    def capsule_outer_node(self, t_idx: int) -> int:
        """Outermost grid node belonging to the capsule at timestep `t_idx`.
        For standard capsules this is the grid outer node (n_zones); for
        halfraum-capsule targets it is the node at the capsule/external boundary.
        Returns n_zones (or 0 if mass_density not loaded) when ri is unavailable.
        """
        if self.mass_density is None:
            return 0
        n_zones = self.mass_density.shape[1]
        if (self.target_class != "halfraum_capsule"
                or self.region_interfaces_indices is None):
            return n_zones
        return int(self.region_interfaces_indices[t_idx, self.capsule_outer_idx])


# ============================================================================
# Variable-name mapping: Helios EXODUS → ICFRunData attribute
# ============================================================================

# Each entry:  (ICFRunData_attr, [list_of_candidate_EXODUS_names], required?)
#
# Candidates are tried in order; the first match wins.
# 'required' means build_run_data will warn (not error) if missing.

_VARIABLE_MAP = [
    # --- Core zone-centered arrays ---
    ("mass_density",      ["mass_density", "dens", "density"],                    True),
    ("ion_temperature",   ["ion_temperature", "temp", "temperature"],             True),
    ("elec_temperature",  ["elec_temperature", "electron_temperature"],           False),
    ("ion_pressure",      ["ion_pressure", "pres", "pressure"],                   True),
    ("elec_pressure",     ["elec_pressure", "electron_pressure"],                 False),
    ("zone_mass",         ["zone_mass", "mass", "cell_mass"],                     True),

    # --- Node-centered arrays (N+1 points) ---
    ("zone_boundaries",   ["zone_boundaries", "coord", "coordx"],                True),
    ("velocity",          ["fluid_velocity", "velocity", "vel"],                  True),

    # --- Optional physics arrays ---
    ("fusion_power",          ["FusionRate_DT_nHe4", "fusion_rate_DT",
                               "fusion_power", "neutron_rate"],                   False),
    ("laser_energy_deposited", ["EnLaserDepositedTimeIntg",
                                "laser_energy_deposited"],                        False),
    # Helios's own cumulative integral of delivered laser power, used as
    # the authoritative E_delivered(t) for coupling diagnostics. Has the
    # same time-stepping as EnLaserDepositedTimeIntg so the
    # E_abs(t) / E_del(t) ratio is self-consistent.
    ("laser_energy_delivered_cum", ["LaserEnDeliveredTimeInt",
                                    "laser_energy_delivered_cum"],                False),
    ("laser_power_source",    ["LaserPwrSrc", "laser_power_source",
                               "laser_power"],                                    False),
    ("laser_power_delivered", ["LaserPwrDeliveredForBeam",
                               "laser_power_delivered"],                           False),
    ("electron_density",      ["electron_density", "elec_density", "n_e"],        False),
    ("ion_density",           ["ion_density", "n_i"],                              False),
    ("mean_charge",           ["mean_charge", "Zbar", "z_bar"],                    False),
    ("laser_power_on_target", ["LaserPwrOnTargetForBeam",
                               "laser_power_on_target"],                           False),
    ("laser_attenuation_coeff", ["laserAttinuationCoeff",
                                 "laser_attenuation_coeff"],                       False),
    ("neutron_production_rate", ["neutron_production_rate", "neutron_rate",
                                 "FusionRate_DD_nHe3"],                           False),
    ("alpha_heating_power",   ["alpha_heating_power", "alpha_power",
                               "alpha_heating"],                                  False),
    ("alpha_heating_ion",     ["pt_particle_heating_ion"],                         False),
    ("alpha_heating_ele",     ["pt_particle_heating_ele"],                         False),
    ("dt_neutron_count",      ["TimeIntFusionProd_n_1406",
                               "dt_neutron_count"],                               False),
    ("dd_neutron_count",      ["TimeIntFusProd_n_0245",
                               "dd_neutron_count"],                               False),
    ("dt_neutron_count_zone", ["TimeIntFusionProd_n_1406_zone"],                  False),
    # Boundary-tally cumulative quantities (J, direct measurements)
    ("radiation_energy_at_boundary_cum",
                              ["TimeIntRadiationLossAtBds"],                       False),
    ("particle_energy_escaped_cum",
                              ["particle_time_int_energy_escaped"],                False),
]


# ============================================================================
# netCDF direct-load helper (for variables with non-standard names)
# ============================================================================

def _try_load_nc(ds, data, attr_name, nc_var_name, verbose=True):
    """Try to load a variable directly from a netCDF4 Dataset."""
    if nc_var_name in ds.variables:
        try:
            arr = np.asarray(ds.variables[nc_var_name][:], dtype=np.float64)
            setattr(data, attr_name, arr)
            if verbose:
                logger.info(f"  ✓ {attr_name:30s} ← '{nc_var_name}' {arr.shape}")
        except Exception as exc:
            logger.warning(f"  ✗ Failed to load '{nc_var_name}' → {attr_name}: {exc}")
    else:
        if verbose:
            logger.debug(f"  – '{nc_var_name}' not in EXODUS file")


# ============================================================================
# Builder
# ============================================================================

def build_run_data(
    run,
    *,
    time_unit: str = "auto",
    compute_derived: bool = True,
    rhw_config=None,
    drive_temperature: Optional[np.ndarray] = None,
    drive_time: Optional[np.ndarray] = None,
    verbose: bool = True,
) -> ICFRunData:
    """
    Bulk-load every time step from a HeliosRun into an ICFRunData container.

    Parameters
    ----------
    run : HeliosRun
        Open HeliosRun instance (from ``helios_postprocess.core``).
    time_unit : str
        Unit of ``run.times``.  ``'auto'`` (default) inspects magnitude:
        max < 1e-3 → seconds (multiplied by 1e9); else assumed ns.
        ``'s'`` forces seconds→ns; ``'ns'`` forces passthrough.
        ``data.time`` on the returned container is always in nanoseconds.
    compute_derived : bool
        If True, compute ``zone_centers`` and ``scale_length`` from loaded data.
    rhw_config : object, optional
        Parsed RHW configuration (from ``rhw_parser``), attached as metadata.
    drive_temperature : np.ndarray, optional
        Drive temperature profile (eV) from RHW file.
    drive_time : np.ndarray, optional
        Time array for drive profile (seconds).
    verbose : bool
        Print progress messages.

    Returns
    -------
    ICFRunData
        Fully populated container ready for ICFAnalyzer / ICFPlotter / OutputGenerator.
    """
    data = ICFRunData()

    # ------------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------------
    data.filename = Path(run.filepath).stem
    data.rhw_config = rhw_config
    if rhw_config is not None:
        data.laser_wavelength_um       = rhw_config.laser_wavelength_um
        data.laser_geometry_per_beam   = getattr(rhw_config, 'laser_geometry_per_beam', None)
        data.laser_spot_size_cm        = rhw_config.laser_spot_size_cm
        data.laser_half_cone_angle_deg = rhw_config.laser_half_cone_angle_deg
        data.laser_focus_position_cm   = rhw_config.laser_focus_position_cm
        data.laser_power_multiplier    = rhw_config.laser_power_multiplier
        data.laser_spatial_profile     = rhw_config.laser_spatial_profile
        data.laser_foot_power_TW       = rhw_config.laser_foot_power_TW
        data.laser_peak_power_TW       = rhw_config.laser_peak_power_TW
        data.laser_foot_start_ns       = rhw_config.laser_foot_start_ns
        data.laser_foot_end_ns         = rhw_config.laser_foot_end_ns
        data.laser_peak_start_ns       = rhw_config.laser_peak_start_ns
        data.laser_peak_end_ns         = rhw_config.laser_peak_end_ns
        # Flux limiter (from .rhw)
        data.flux_limiter            = getattr(rhw_config, 'flux_limiter', 0.0)
        data.flux_limiter_enabled    = getattr(rhw_config, 'flux_limiter_enabled', False)
        data.flux_limiter_per_region = getattr(rhw_config, 'flux_limiter_per_region', None)
        data.laser_pulse_duration_ns   = rhw_config.laser_pulse_duration_ns
        data.eos_models                = rhw_config.eos_models
        data.alpha_deposition_local    = rhw_config.alpha_deposition_local
        data.alpha_deposition_nonlocal = rhw_config.alpha_deposition_nonlocal
        # IDD source flags (rad-source-at-Rmin / Rmax). RHW reports both
        # booleans even when neither is on; we expose the active count.
        data.rad_source_rmin_on    = bool(getattr(rhw_config, 'rad_source_rmin_on', False))
        data.rad_source_rmax_on    = bool(getattr(rhw_config, 'rad_source_rmax_on', False))
        data.n_active_idd_sources  = int(data.rad_source_rmin_on) + int(data.rad_source_rmax_on)
        
    data.drive_time = drive_time
    data.drive_temperature = drive_temperature

    # ------------------------------------------------------------------
    # Time array
    # ------------------------------------------------------------------
    raw_times = np.asarray(run.times, dtype=np.float64)
    if time_unit == "auto":
        # Implosion times are O(1-30) ns; raw seconds put max ~1e-8.
        # Anything below 1e-3 is unambiguously in seconds.
        if raw_times.size > 0 and float(np.max(raw_times)) < 1e-3:
            data.time = raw_times * 1e9
        else:
            data.time = raw_times.copy()
    elif time_unit == "s":
        data.time = raw_times * 1e9          # seconds → nanoseconds
    elif time_unit == "ns":
        data.time = raw_times.copy()
    else:
        raise ValueError(f"Unknown time_unit: {time_unit!r}  (use 'auto', 's', or 'ns')")

    n_times = len(data.time)
    if verbose:
        logger.info(f"Loading {n_times} time steps  "
                     f"({data.time[0]:.4f} – {data.time[-1]:.4f} ns)")

    # Build a set of available variable names once (avoids repeated queries)
    available_vars = set(run.list_variables())

    # ------------------------------------------------------------------
    # Load 2-D arrays via the variable map
    # ------------------------------------------------------------------
    for attr_name, candidates, required in _VARIABLE_MAP:
        loaded = False
        for var_name in candidates:
            if var_name in available_vars:
                try:
                    arr = run.get_variable(var_name)          # all time steps
                    setattr(data, attr_name, np.asarray(arr, dtype=np.float64))
                    loaded = True
                    if verbose:
                        logger.info(f"  ✓ {attr_name:30s} ← '{var_name}' {arr.shape}")
                    break
                except Exception as exc:
                    logger.warning(f"  ✗ Failed to load '{var_name}' → {attr_name}: {exc}")

        if not loaded:
            if required:
                logger.warning(f"  ⚠ REQUIRED variable '{attr_name}' not found "
                               f"(tried: {candidates})")
            elif verbose:
                logger.debug(f"  – Optional '{attr_name}' not available")

    # ------------------------------------------------------------------
    # Post-processing: collapse beam-indexed arrays to 1-D total
    #
    # Single-beam runs:  shape (n_t, 1) -> squeeze to (n_t,).
    # Multi-beam runs:   shape (n_t, k>1) -> sum across beam axis to give
    #                    the total laser power on the target, with the
    #                    per-beam table preserved on `*_per_beam` for
    #                    diagnostics.
    #
    # Prior behaviour silently kept multi-beam arrays as 2-D for
    # `laser_power_delivered` (np.squeeze is a no-op on shape (n_t, k>1))
    # and silently dropped beams 1..k-1 for `laser_power_on_target` by
    # taking [:, 0]. Both broke downstream argmax/peak-power logic.
    # ------------------------------------------------------------------
    def _collapse_beam_axis(arr, name):
        """Return (total_1d, per_beam_or_None). Logs the action."""
        if arr is None or arr.ndim == 1:
            return arr, None
        n_beam = arr.shape[-1]
        if n_beam == 1:
            total = arr.squeeze(axis=-1)
            if verbose:
                logger.info(f"  ✓ {name:<28s} squeezed → {total.shape}")
            return total, None
        per_beam = arr.copy()
        total = arr.sum(axis=-1)
        if verbose:
            logger.info(f"  ✓ {name:<28s} summed across {n_beam} beams → {total.shape}, "
                        f"per-beam preserved → {per_beam.shape}")
        return total, per_beam

    data.laser_power_delivered, data.laser_power_delivered_per_beam = \
        _collapse_beam_axis(data.laser_power_delivered, "laser_power_delivered")
    data.laser_power_on_target, data.laser_power_on_target_per_beam = \
        _collapse_beam_axis(data.laser_power_on_target, "laser_power_on_target")

    # Count beams that actually carry power. Helios's RHW format always
    # reports the max array dim (typically 3 beams) even when PDD targets
    # only drive beam 1. Threshold = 1 W on peak delivered power, well
    # above any numerical noise; sub-MW beams aren't doing physics.
    BEAM_POWER_THRESHOLD_W = 1.0
    pb = data.laser_power_delivered_per_beam
    if pb is not None and pb.ndim == 2:
        peak_per_beam = np.max(pb, axis=0)             # (n_beam,)
        active_mask   = peak_per_beam > BEAM_POWER_THRESHOLD_W
        data.n_total_beams  = int(pb.shape[1])
        data.n_active_beams = int(active_mask.sum())
        if verbose and data.n_total_beams != data.n_active_beams:
            inactive_idx = np.where(~active_mask)[0]
            logger.info(
                f"  ✓ Active beams: {data.n_active_beams}/{data.n_total_beams}  "
                f"(inactive beam idx: {list(inactive_idx)}, peak <= "
                f"{BEAM_POWER_THRESHOLD_W:.0e} W)"
            )
    elif data.laser_power_delivered is not None:
        data.n_total_beams  = 1
        data.n_active_beams = 1
    else:
        data.n_total_beams  = 0
        data.n_active_beams = 0

    # ------------------------------------------------------------------
    # rad_pressure: non-ion pressure component
    # ------------------------------------------------------------------
    # Downstream code computes total_pressure = ion_pressure + rad_pressure
    # everywhere.  We set rad_pressure = (everything that isn't ion):
    #
    #   - If EXODUS has *both* elec_pressure and rad_pressure:
    #       data.rad_pressure = elec_pressure + exodus_rad_pressure
    #   - If EXODUS has only elec_pressure:
    #       data.rad_pressure = elec_pressure      (PDD 8 path)
    #   - If neither:
    #       data.rad_pressure = zeros
    # ------------------------------------------------------------------
    actual_rad = None
    if "rad_pressure" in available_vars:
        try:
            actual_rad = np.asarray(run.get_variable("rad_pressure"), dtype=np.float64)
            if verbose:
                logger.info(f"  ✓ {'actual rad_pressure':30s} ← 'rad_pressure' {actual_rad.shape}")
        except Exception as exc:
            logger.warning(f"  ✗ Failed to load rad_pressure: {exc}")

    if data.elec_pressure is not None and actual_rad is not None:
        data.rad_pressure = data.elec_pressure + actual_rad
        if verbose:
            logger.info("  ✓ rad_pressure              ← elec_pressure + rad_pressure (3-component)")
    elif data.elec_pressure is not None:
        data.rad_pressure = data.elec_pressure
        if verbose:
            logger.info("  ✓ rad_pressure              ← elec_pressure (2-component)")
    elif actual_rad is not None:
        data.rad_pressure = actual_rad
        if verbose:
            logger.info("  ✓ rad_pressure              ← rad_pressure only")
    else:
        if data.ion_pressure is not None:
            data.rad_pressure = np.zeros_like(data.ion_pressure)
            logger.warning("  ⚠ No elec/rad pressure found — rad_pressure set to zeros")
        else:
            data.rad_pressure = None

    # ------------------------------------------------------------------
    # Derived: plasma_pressure = ion + electron  (CONVENTIONAL DEFINITION)
    # This is what should be used for shock breakout, ablation pressure,
    # adiabat, hot-spot pressure -- all the standard ICF diagnostics.
    # The legacy `rad_pressure` field still carries elec+rad combined for
    # backward compatibility but should be considered deprecated.
    # `rad_pressure_true` carries only the actual radiation pressure if
    # available, otherwise zeros.
    # ------------------------------------------------------------------
    if data.ion_pressure is not None and data.elec_pressure is not None:
        data.plasma_pressure = data.ion_pressure + data.elec_pressure
        if verbose:
            logger.info("  ✓ plasma_pressure           ← ion_pressure + elec_pressure (conventional)")
    elif data.ion_pressure is not None:
        # Fallback: assume Helios treats ion = electron locally so use ion alone
        data.plasma_pressure = data.ion_pressure
        logger.warning("  ⚠ No electron pressure found -- plasma_pressure = ion_pressure only")

    if actual_rad is not None:
        data.rad_pressure_true = actual_rad
    elif data.ion_pressure is not None:
        data.rad_pressure_true = np.zeros_like(data.ion_pressure)

    # ------------------------------------------------------------------
    # Region interfaces, material index, region names
    # ------------------------------------------------------------------
    # These EXODUS variables have non-standard names (spaces, long labels)
    # so we load them directly from the underlying netCDF dataset.
    # ------------------------------------------------------------------
    ds = run.dataset if hasattr(run, 'dataset') else None
    if ds is None and hasattr(run, 'ds'):
        ds = run.ds

    if ds is not None:
        # Region interface indices: (n_times, n_regions)
        _try_load_nc(ds, data, "region_interfaces_indices",
                     "Indices at region interfaces", verbose)

        # Material index: (n_zones,) — static, maps each zone to a material
        _try_load_nc(ds, data, "material_index",
                     "Material index", verbose)

        # Region names: (n_regions, 100) char array → list of strings
        if "name_spatial_regn" in ds.variables:
            try:
                raw = ds.variables["name_spatial_regn"][:]
                names = []
                for row in raw:
                    name = b"".join(row).decode("utf-8", errors="ignore").rstrip("\x00").strip()
                    if name:
                        names.append(name)
                data.region_names = names if names else None
                if verbose and data.region_names:
                    logger.info(f"  ✓ {'region_names':30s} ← {data.region_names}")
            except Exception as exc:
                logger.warning(f"  ✗ Could not read region names: {exc}")

        # ----------- Detect target class from region names ------------
        # Halfraum/hohlraum-style targets append external structural regions
        # (low-density gas fill + a high-Z converter shell) outside the capsule.
        # We scan region names from the OUTSIDE inward and tag regions whose
        # names match known external-structure keywords; the first non-external
        # region from the outside marks the capsule outer surface.
        #
        # If no external regions are detected the run is treated as a
        # standard capsule (existing behaviour).
        _EXTERNAL_REGION_KEYWORDS = (
            'pseudo void', 'void',
            'cu shell', 'pb shell', 'au shell', 'u shell', 'ta shell',
            'hohlraum', 'halfraum',
            'he fill', 'helium', 'hohlraum gas',
        )
        if data.region_names:
            n_total = len(data.region_names)
            n_external = 0
            for name in reversed(data.region_names):
                name_lc = name.lower().strip()
                if any(kw in name_lc for kw in _EXTERNAL_REGION_KEYWORDS):
                    n_external += 1
                else:
                    break  # stop at first non-external region from outside
            data.n_capsule_regions = n_total - n_external
            if n_external > 0 and data.n_capsule_regions >= 2:
                data.target_class = "halfraum_capsule"
                if verbose:
                    external_names = data.region_names[data.n_capsule_regions:]
                    capsule_names = data.region_names[:data.n_capsule_regions]
                    logger.info(f"  ✓ {'target_class':30s} ← halfraum_capsule")
                    logger.info(f"    capsule  ({data.n_capsule_regions}): {capsule_names}")
                    logger.info(f"    external ({n_external}): {external_names}")
            else:
                data.target_class = "capsule"
                if verbose:
                    logger.info(f"  ✓ {'target_class':30s} ← capsule "
                                f"({data.n_capsule_regions} regions)")
    else:
        logger.warning("  ⚠ No netCDF dataset handle — region/material info unavailable")

    # ------------------------------------------------------------------
    # Derived arrays
    # ------------------------------------------------------------------
    if compute_derived:
        _compute_zone_centers(data, verbose)
        _compute_scale_length(data, verbose)

    # ------------------------------------------------------------------
    # Quick sanity summary
    # ------------------------------------------------------------------
    if verbose:
        _log_summary(data)

    return data


# ============================================================================
# Derived-quantity helpers
# ============================================================================

def _compute_zone_centers(data: ICFRunData, verbose: bool = True) -> None:
    """
    Compute zone centres from zone boundaries.

    zone_boundaries has shape (n_times, N+1) for N zones.
    zone_centers will have shape (n_times, N).
    """
    if data.zone_boundaries is None:
        return

    try:
        bnd = data.zone_boundaries               # (n_times, n_nodes)
        data.zone_centers = 0.5 * (bnd[:, :-1] + bnd[:, 1:])
        if verbose:
            logger.info(f"  ✓ zone_centers               derived  {data.zone_centers.shape}")
    except Exception as exc:
        logger.warning(f"  ✗ Could not compute zone_centers: {exc}")
        data.zone_centers = None


def _compute_scale_length(data: ICFRunData, verbose: bool = True) -> None:
    """
    Compute density scale length L = ρ / (dρ/dr) at each zone and time step.

    Used by the ablation-front tracker in ICFAnalyzer.
    """
    if data.mass_density is None or data.zone_centers is None:
        return

    try:
        n_times, n_zones = data.mass_density.shape
        scale_length = np.zeros((n_times, n_zones), dtype=np.float64)

        for t in range(n_times):
            rho = data.mass_density[t]
            r   = data.zone_centers[t]

            # Central-difference gradient (numpy handles endpoints)
            drho_dr = np.gradient(rho, r)

            # Avoid divide-by-zero
            with np.errstate(divide="ignore", invalid="ignore"):
                L = rho / drho_dr
                L[~np.isfinite(L)] = 0.0

            scale_length[t] = L

        data.scale_length = scale_length
        if verbose:
            logger.info(f"  ✓ scale_length               derived  {scale_length.shape}")

    except Exception as exc:
        logger.warning(f"  ✗ Could not compute scale_length: {exc}")
        data.scale_length = None


# ============================================================================
# Summary logger
# ============================================================================

def _log_summary(data: ICFRunData) -> None:
    """Log a short summary of what was loaded."""
    n_times = len(data.time) if data.time is not None else 0
    n_zones = data.mass_density.shape[1] if data.mass_density is not None else "?"

    lines = [
        "",
        "═" * 60,
        "  ICFRunData Summary",
        "═" * 60,
        f"  File:            {data.filename}",
        f"  Time steps:      {n_times}",
        f"  Zones:           {n_zones}",
    ]

    if data.time is not None and n_times > 0:
        lines.append(f"  Time range:      {data.time[0]:.4f} – {data.time[-1]:.4f} ns")

    # List loaded 2-D arrays
    loaded = []
    missing = []
    for attr, _, _ in _VARIABLE_MAP:
        val = getattr(data, attr, None)
        if val is not None:
            loaded.append(attr)
        else:
            missing.append(attr)

    lines.append(f"  Loaded arrays:   {len(loaded)}")
    if missing:
        lines.append(f"  Missing arrays:  {', '.join(missing)}")

    # Derived
    derived = []
    if data.zone_centers is not None:
        derived.append("zone_centers")
    if data.scale_length is not None:
        derived.append("scale_length")
    if derived:
        lines.append(f"  Derived arrays:  {', '.join(derived)}")

    # Region / material info
    if data.region_interfaces_indices is not None:
        ri = data.region_interfaces_indices[0].astype(int)
        n_regions = len(ri)
        lines.append(f"  Regions:         {n_regions}")
        # Build zone range description
        prev = 0
        for i, boundary in enumerate(ri):
            name = data.region_names[i] if data.region_names and i < len(data.region_names) else f"Region {i+1}"
            lines.append(f"    {i+1}. {name:20s}  zones {prev:3d}–{int(boundary)-1:3d}  ({int(boundary)-prev:3d} zones)")
            prev = int(boundary)

    if data.material_index is not None:
        unique_mats = np.unique(data.material_index.astype(int))
        lines.append(f"  Materials:       {len(unique_mats)} ({', '.join(str(m) for m in unique_mats)})")

    lines.append("═" * 60)

    for line in lines:
        logger.info(line)
