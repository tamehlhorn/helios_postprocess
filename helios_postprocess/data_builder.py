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
from typing import Optional

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
        self.rad_pressure: Optional[np.ndarray] = None        # (n_times, n_zones) J/cm³ — alias, see note below
        self.velocity: Optional[np.ndarray] = None            # (n_times, n_nodes) cm/s — node-centered
        self.zone_boundaries: Optional[np.ndarray] = None     # (n_times, n_nodes) cm — node-centered
        self.zone_mass: Optional[np.ndarray] = None           # (n_times, n_zones) g
        self.zone_centers: Optional[np.ndarray] = None        # (n_times, n_zones) cm — derived

        # ------------------------------------------------------------------
        # 2. Optional physics arrays (None when absent)
        # ------------------------------------------------------------------
        self.fusion_power: Optional[np.ndarray] = None            # (n_times, n_zones) — zone-level DT fusion rate
        self.laser_energy_deposited: Optional[np.ndarray] = None  # (n_times, n_zones) J  (time-integrated)
        self.laser_power_source: Optional[np.ndarray] = None      # (n_times, n_zones) W/cm³
        self.electron_density: Optional[np.ndarray] = None        # (n_times, n_zones) cm⁻³
        self.neutron_production_rate: Optional[np.ndarray] = None # (n_times, n_zones)
        self.alpha_heating_power: Optional[np.ndarray] = None     # (n_times, n_zones)
        self.dt_neutron_count: Optional[np.ndarray] = None        # time-integrated DT neutron yield
        self.dd_neutron_count: Optional[np.ndarray] = None        # time-integrated DD neutron yield
        self.scale_length: Optional[np.ndarray] = None            # (n_times, n_zones) cm — derived
        self.region_interfaces_indices: Optional[np.ndarray] = None  # (n_times, n_regions)
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
        self.core_radius: float = 0.0            # cm

        # Hot spot (at stagnation)
        self.stagnation_hot_spot_radius: float = 0.0  # cm
        self.hot_spot_pressure: float = 0.0            # Mbar
        self.hot_spot_areal_density: float = 0.0       # g/cm²
        self.hot_spot_internal_energy: float = 0.0     # kJ

        # Areal densities
        self.time_ave_areal_density: float = 0.0       # g/cm²
        self.bang_time_areal_density: float = 0.0      # g/cm²
        self.bang_time_fuel_areal_density: float = 0.0 # g/cm²
        self.bang_time_HDC_areal_density: float = 0.0  # g/cm²
        self.bang_time_hs_areal_density: float = 0.0   # g/cm²

        # Neutron-averaged quantities
        self.neutron_ave_fuel_areal_density: float = 0.0  # g/cm²
        self.neutron_ave_ion_temperature: float = 0.0     # eV
        self.neutron_ave_pressure: float = 0.0            # Mbar

        # Energy / performance
        self.energy_output: float = 0.0    # MJ  (fusion yield)
        self.laser_energy: float = 0.0     # MJ
        self.rad_energy: float = 0.0       # MJ  (radiation / IDD deposited)
        self.target_gain: float = 0.0
        self.max_dt_temp: float = 0.0      # eV

        # Mass fractions
        self.unablated_ablatar_mass: float = 0.0
        self.unablated_fuel_mass: float = 0.0
        self.stagnated_fuel_mass: float = 0.0

        # Implosion
        self.peak_implosion_velocity: float = 0.0  # km/s
        self.adiabat_mass_averaged_ice: float = 0.0

        # Critical density tracking (populated by ICFAnalyzer)
        self.critical_density_value: float = 0.0
        self.critical_density_radius_formula: Optional[np.ndarray] = None
        self.critical_density_radius_algorithm: Optional[np.ndarray] = None

        # Ablation front (populated by ICFAnalyzer)
        self.ablation_front_radius: Optional[np.ndarray] = None
        self.min_scale_length: Optional[np.ndarray] = None

        # Shock tracking (populated by ICFAnalyzer)
        self.shock_times: list = []
        self.shock_radii: list = []
        self.shock_velocities: list = []

        # ------------------------------------------------------------------
        # 4. Metadata
        # ------------------------------------------------------------------
        self.filename: str = ""
        self.rhw_config = None          # RHWConfig dataclass (from rhw_parser)
        self.drive_temperature: Optional[np.ndarray] = None  # eV  (from .rhw file)
        self.drive_time: Optional[np.ndarray] = None         # seconds


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
    ("laser_power_source",    ["LaserPwrSrc", "laser_power_source",
                               "laser_power"],                                    False),
    ("electron_density",      ["electron_density", "elec_density", "n_e"],        False),
    ("neutron_production_rate", ["neutron_production_rate", "neutron_rate",
                                 "FusionRate_DD_nHe3"],                           False),
    ("alpha_heating_power",   ["alpha_heating_power", "alpha_power",
                               "alpha_heating"],                                  False),
    ("dt_neutron_count",      ["TimeIntFusionProd_n_1406",
                               "dt_neutron_count"],                               False),
    ("dd_neutron_count",      ["TimeIntFusProd_n_0245",
                               "dd_neutron_count"],                               False),
]


# ============================================================================
# Builder
# ============================================================================

def build_run_data(
    run,
    *,
    time_unit: str = "ns",
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
        Unit of ``run.times``.  ``'s'`` → converted to ns;
        ``'ns'`` → used as-is (Helios default after netCDF read is seconds,
        but some pipelines have already converted).
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
    data.drive_temperature = drive_temperature
    data.drive_time = drive_time

    # ------------------------------------------------------------------
    # Time array
    # ------------------------------------------------------------------
    raw_times = np.asarray(run.times, dtype=np.float64)
    if time_unit == "s":
        data.time = raw_times * 1e9          # seconds → nanoseconds
    elif time_unit == "ns":
        data.time = raw_times.copy()
    else:
        raise ValueError(f"Unknown time_unit: {time_unit!r}  (use 's' or 'ns')")

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
    # rad_pressure alias
    # ------------------------------------------------------------------
    # The legacy modules reference ``data.rad_pressure`` (ion + rad) for
    # total-pressure calculations.  In Helios, the physically meaningful
    # second component is *electron* pressure (CLAUDE.md: "Total pressure:
    # ion_pressure + elec_pressure").  We map elec_pressure → rad_pressure
    # so existing downstream code works without modification.
    #
    # If a true radiation-pressure variable is ever added to the EXODUS
    # file, a more explicit split should be introduced.
    # ------------------------------------------------------------------
    if data.elec_pressure is not None:
        data.rad_pressure = data.elec_pressure
        if verbose:
            logger.info("  ✓ rad_pressure              ← elec_pressure (alias)")
    else:
        # Fallback: zero array so (ion + rad) == ion-only and nothing crashes
        if data.ion_pressure is not None:
            data.rad_pressure = np.zeros_like(data.ion_pressure)
            logger.warning("  ⚠ elec_pressure not found — rad_pressure set to zeros")
        else:
            data.rad_pressure = None

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

    lines.append("═" * 60)

    for line in lines:
        logger.info(line)
