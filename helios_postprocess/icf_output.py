"""
ICF Output Module
=================
Text and CSV output for ICF simulation analysis results.

Produces two files from a populated ICFRunData object:
  - ``{basename}_summary.txt``  — all scalar metrics, human-readable
  - ``{basename}_history.csv``  — time-history columns in a single CSV

Unit conventions (per CLAUDE.md):
  - Temperatures: keV  (except ion_temperature arrays which stay in eV)
  - Pressures:    Gbar (1 Gbar = 1e8 J/cm³)
  - Velocities:   km/s
  - Areal density: g/cm²
  - Drive temperature: eV

Adapted for HeliosRun / data_builder / icf_analysis pipeline.

Author: Prof T / Xcimer ICF Analysis
Date: 2025
"""

import csv
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import logging

logger = logging.getLogger(__name__)


class ICFOutputGenerator:
    """Generate text summary and CSV time-history files from ICFRunData."""

    def __init__(self, data, config: Optional[dict] = None):
        """
        Parameters
        ----------
        data : ICFRunData
            Populated container (after ICFAnalyzer has run).
        config : dict, optional
            Reserved for future formatting options.
        """
        self.data = data
        self.config = config or {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def write_all(self, base_path: str) -> None:
        """
        Write both summary and time-history files.

        Parameters
        ----------
        base_path : str
            Path *without* extension.  E.g. ``/path/to/Olson_PDD_9``
            produces ``Olson_PDD_9_summary.txt`` and
            ``Olson_PDD_9_history.csv``.
        """
        base = str(base_path)
        # Strip any trailing extension the caller may have left on
        for ext in ('.exo', '.txt', '.csv', '.pdf'):
            if base.lower().endswith(ext):
                base = base[:-len(ext)]
                break

        summary_path = base + '_summary.txt'
        history_path = base + '_history.csv'

        self.write_summary(summary_path)
        self.write_time_histories(history_path)

    # ------------------------------------------------------------------
    # Summary text file
    # ------------------------------------------------------------------

    def write_summary(self, path: str) -> None:
        """Write scalar metrics to a human-readable text file."""
        logger.info(f"Writing summary: {path}")
        d = self.data

        lines = []
        _a = lines.append   # shorthand

        width = 60
        _a('=' * width)
        _a(f'  ICF Analysis Summary — {d.filename}')
        _a(f'  Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
        _a('=' * width)
        _a('')

        # ---- Simulation overview ----
        _a('SIMULATION OVERVIEW')
        _a('-' * width)
        _a(f'  {"File:":<36s} {d.filename}')
        if d.time is not None:
            _a(f'  {"Time steps:":<36s} {len(d.time)}')
            _a(f'  {"Time range:":<36s} {d.time[0]:.4f} – {d.time[-1]:.4f} ns')
        if d.mass_density is not None:
            _a(f'  {"Zones:":<36s} {d.mass_density.shape[1]}')
        if d.region_names:
            _a(f'  {"Regions:":<36s} {len(d.region_names)}')
            ri0 = d.region_interfaces_indices[0].astype(int) if d.region_interfaces_indices is not None else None
            for i, name in enumerate(d.region_names):
                if ri0 is not None and i < len(ri0):
                    z_start = 0 if i == 0 else int(ri0[i - 1])
                    z_end = int(ri0[i]) - 1 if i < len(ri0) else '?'
                    _a(f'    {i+1}. {name:<20s} zones {z_start}–{z_end}')
                else:
                    _a(f'    {i+1}. {name}')
        _a('')
        
        # ---- Alpha burn model ----
        _local    = getattr(d, 'alpha_deposition_local',    False)
        _nonlocal = getattr(d, 'alpha_deposition_nonlocal', False)
        if _local and _nonlocal:
            burn_str = "Full (local + non-local transport)"
        elif _local:
            burn_str = "Local only (instantaneous — use with caution)"
        elif _nonlocal:
            burn_str = "Non-local transport (time/space dependent)"
        else:
            burn_str = "Disabled"
        _a(f"  {'Alpha burn model':<36s} {burn_str}")
        _a('')

        # ---- Timing ----
        _a('TIMING')
        _a('-' * width)
        _a(self._metric('Stagnation time',      d.stag_time,   'ns',   fmt='.3f'))
        _a(self._metric('Bang time',            d.bang_time,   'ns',   fmt='.3f'))
        _a(self._metric('Burn width (FWHM)',    d.burn_width,  'ns',   fmt='.4f'))
        _a('')

        # ---- Compression / density ----
        _a('COMPRESSION')
        _a('-' * width)
        _a(self._metric('Peak density',         d.max_density,  'g/cc', fmt='.2f'))
        _a(self._metric('Convergence ratio CR',  d.comp_ratio,   '',     fmt='.1f'))
        _a(self._metric('Core radius',          d.core_radius,  'cm',   fmt='.4f'))
        _a('')

        # ---- Hot spot at stagnation ----
        _a('HOT SPOT (at stagnation)')
        _a('-' * width)
        _a(self._metric('Hot-spot radius',             d.stagnation_hot_spot_radius, 'cm',   fmt='.4f'))
        _a(self._metric('Hot-spot pressure',           d.hot_spot_pressure,          'Gbar', fmt='.2f'))
        _a(self._metric('Hot-spot areal density',      d.hot_spot_areal_density,     'g/cm²', fmt='.4f'))
        _a(self._metric('Hot-spot internal energy',    d.hot_spot_internal_energy,   'kJ',   fmt='.2f'))
        _a('')

        # ---- Areal densities at bang time ----
        _a('AREAL DENSITIES (at bang time)')
        _a('-' * width)
        _a(self._metric('Total ρR',            d.bang_time_areal_density,      'g/cm²', fmt='.4f'))
        _a(self._metric('Hot-spot ρR',         d.bang_time_hs_areal_density,   'g/cm²', fmt='.4f'))
        _a(self._metric('Cold-fuel ρR',        d.bang_time_fuel_areal_density, 'g/cm²', fmt='.4f'))
        _a(self._metric('Ablator ρR',          d.bang_time_HDC_areal_density,  'g/cm²', fmt='.4f'))
        _a(self._metric('Time-avg ρR (burn)',  d.time_ave_areal_density,       'g/cm²', fmt='.4f'))
        _a('')

        # ---- Neutron-averaged quantities ----
        _a('NEUTRON-AVERAGED QUANTITIES')
        _a('-' * width)
        _a(self._metric('⟨Ti⟩_n',             d.neutron_ave_ion_temperature,   'keV',   fmt='.2f'))
        _a(self._metric('⟨P⟩_n',              d.neutron_ave_pressure,          'Gbar',  fmt='.2f'))
        _a(self._metric('⟨ρR⟩_n (fuel)',      d.neutron_ave_fuel_areal_density,'g/cm²', fmt='.4f'))
        _a('')

        # ---- Energy / performance ----
        _a('ENERGY & PERFORMANCE')
        _a('-' * width)
        _a(self._metric('Laser energy',        d.laser_energy,     'MJ',  fmt='.3f'))
        _a(self._metric('Fusion yield',        d.energy_output,    'MJ',  fmt='.3f'))
        _a(self._metric('DT neutron yield',    d.dt_neutron_yield, '',    fmt='.3e'))
        _a(self._metric('Target gain',         d.target_gain,      '',    fmt='.3f'))
        _a(self._metric('Max DT temperature',  d.max_dt_temp,      'keV', fmt='.2f'))
        _a('')

        # ---- Laser configuration ----
        if d.laser_wavelength_um > 0:
            _a('LASER CONFIGURATION (beam 1)')
            _a('-' * width)
            _a(self._metric('Wavelength',       d.laser_wavelength_um,       'um',  fmt='.3f'))
            _a(f"  {'Focus position':<30} {d.laser_focus_position_cm:>15.4f} cm")
            _a(self._metric('Half-cone angle',  d.laser_half_cone_angle_deg, 'deg', fmt='.2f'))
            _a(f"  {'Spot radius':<30} {d.laser_spot_size_cm:>15.4f} cm  ({d.laser_spatial_profile})")
            _a(self._metric('Power multiplier', d.laser_power_multiplier,    '',    fmt='.4f'))
            # Flux limiter: bypass _metric (0.06 is a valid value; _metric renders 0.0 as '—')
            if getattr(d, 'flux_limiter_enabled', False):
                _a(f"  {'Flux limiter (f)':<36s} {d.flux_limiter:>10.3f}")
            elif getattr(d, 'flux_limiter', 0.0) > 0:
                _a(f"  {'Flux limiter (f)':<36s} {d.flux_limiter:>10.3f}  (flag off)")
            else:
                _a(f"  {'Flux limiter (f)':<36s} {'(not set)':>10s}")
            # Pulse shape
            if d.laser_peak_power_TW > 0:
                _a(f"  {'Pulse shape (beam 1)':<30}")
                if d.laser_foot_power_TW > 0:
                    _a(f"    {'Foot power':<28} {d.laser_foot_power_TW:>10.1f} TW"
                       f"  ({d.laser_foot_start_ns:.2f} – {d.laser_foot_end_ns:.2f} ns)")
                _a(f"    {'Peak power':<28} {d.laser_peak_power_TW:>10.1f} TW"
                   f"  ({d.laser_peak_start_ns:.2f} – {d.laser_peak_end_ns:.2f} ns)")
                _a(f"    {'Pulse duration':<28} {d.laser_pulse_duration_ns:>10.2f} ns")
            _a('')

        # ---- EOS models ----
        if getattr(d, 'eos_models', None):
            _a('EOS MODELS')
            _a('-' * width)
            for m in d.eos_models:
                _a(f"  {m['region']:<25} {m['type']:<10} {m['file']}")
            _a('')

        # ---- Mass fractions ----
        _a('MASS FRACTIONS (at stagnation)')
        _a('-' * width)
        _a(self._metric('Initial DT mass',     d.initial_fuel_mass_mg,    'mg', fmt='.3f'))
        _a(self._metric('Initial ablator mass',d.initial_ablator_mass_mg, 'mg', fmt='.3f'))
        _a(self._metric('Unablated fuel',      d.unablated_fuel_mass,     '', fmt='.4f'))
        _a(self._metric('Unablated ablator',   d.unablated_ablatar_mass,  '', fmt='.4f'))
        _a('')

        # ---- Implosion ----
        _a('IMPLOSION')
        _a('-' * width)
        _a(self._metric('Peak implosion velocity',  abs(d.peak_implosion_velocity), 'km/s', fmt='.2f'))
        _a(self._metric('In-flight CR', getattr(d, 'cr_inflight', None), '', fmt='.1f'))
        _a(self._metric('Mass-avg adiabat (ice)',   d.adiabat_mass_averaged_ice,    '',     fmt='.2f'))
        if getattr(d, 'shock_foot_pressure_Gbar', 0.0) > 0:
            _a(self._metric('Shock foot pressure (DT ice)', d.shock_foot_pressure_Gbar, 'Gbar', fmt='.4f'))
        if getattr(d, 'shock_breakout_time_ns', 0.0) > 0:
            _a(self._metric('Shock breakout time',     d.shock_breakout_time_ns,       'ns',   fmt='.3f'))
            _a(self._metric('Shock breakout pressure', d.shock_breakout_pressure_Gbar, 'Gbar', fmt='.3f'))
            _a(self._metric('Shock breakout Mach',     d.shock_breakout_mach,          '',     fmt='.2f'))
        _a('')
        if (hasattr(d, 'alpha_onset_time_ns') and
                d.alpha_onset_time_ns is not None and
                d.alpha_onset_time_ns > 0 and
                d.alpha_onset_time_ns < d.stag_time * 1e9):
            _a(self._metric('Alpha onset time', d.alpha_onset_time_ns, 'ns', fmt='.3f'))
            _a('  (implosion metrics evaluated pre-onset)')
        
        # ---- Burn propagation (Olson et al. convention) ----
        _a('BURN PROPAGATION (T_ion > 4.5 keV)')
        _a('-' * width)
        _a(self._metric('Ignition time (ρR_hs ≥ 0.3)',  d.ignition_time,          'ns',   fmt='.3f'))
        _a(self._metric('HS radius at ignition',
                         d.ignition_hs_radius * 1e4 if d.ignition_hs_radius > 0 else 0.0,
                         'μm', fmt='.0f'))
        _a(self._metric('HS pressure at ignition',      d.ignition_hs_pressure,   'Gbar', fmt='.1f'))
        _a(self._metric('Complete propagation time',     d.burn_propagation_time,  'ns',   fmt='.3f'))
        _a('')

        _a('=' * width)
        _a('  End of summary')
        _a('=' * width)

        with open(path, 'w') as f:
            f.write('\n'.join(lines) + '\n')

        logger.info(f"Summary written: {path}")

    # ------------------------------------------------------------------
    # CSV time histories
    # ------------------------------------------------------------------

    def write_time_histories(self, path: str) -> None:
        """
        Write a single CSV with one row per timestep.

        Columns (all present even if data is missing — filled with NaN):
          time_ns              — simulation time
          laser_power_W        — total laser power on target
          fusion_power_total   — summed fusion rate over all zones
          hs_radius_cm         — hot-spot outer radius
          total_rhoR_gcm2      — total areal density
          fuel_rhoR_gcm2       — cold-fuel areal density
          hs_rhoR_gcm2         — hot-spot areal density (T_ion > 4.5 keV)
          peak_density_gcc     — peak density (over zones) at each time
          peak_ion_temp_keV    — peak ion temperature at each time
          ablation_front_r_cm  — ablation front radius
          hs_pressure_Gbar     — mass-averaged hot-spot pressure
          dt_neutron_count     — cumulative DT neutron yield
        """
        logger.info(f"Writing time histories: {path}")
        d = self.data

        if d.time is None:
            logger.warning("No time array — cannot write histories")
            return

        n_times = len(d.time)
        nan_col = np.full(n_times, np.nan)

        # ---- Build each column ----

        # Laser power (W)
        if d.laser_power_delivered is not None:
            laser_power = d.laser_power_delivered.copy()
        else:
            laser_power = nan_col.copy()

        # Fusion power (total reaction rate summed over zones)
        if d.fusion_power is not None:
            fusion_total = np.sum(d.fusion_power, axis=1)
        else:
            fusion_total = nan_col.copy()

        # Hot-spot outer radius (cm)
        hs_radius = self._compute_hs_radius_vs_time()

        # Total areal density (g/cm²)
        if d.areal_density_vs_time is not None:
            total_rhoR = d.areal_density_vs_time[:, -1]
        else:
            total_rhoR = nan_col.copy()

        # Cold-fuel areal density (g/cm²)
        fuel_rhoR = self._compute_fuel_rhoR_vs_time()

        # Peak density (g/cc)
        if d.mass_density is not None:
            peak_density = np.max(d.mass_density, axis=1)
        else:
            peak_density = nan_col.copy()

        # Peak ion temperature (keV)
        if d.ion_temperature is not None:
            peak_ion_temp = np.max(d.ion_temperature, axis=1) / 1000.0  # eV → keV
        else:
            peak_ion_temp = nan_col.copy()

        # Ablation front radius (cm)
        if d.ablation_front_radius is not None:
            abl_front = d.ablation_front_radius.copy()
        else:
            abl_front = nan_col.copy()

        # Hot-spot pressure (Gbar, mass-averaged)
        hs_pressure = self._compute_hs_pressure_vs_time()

        # Cumulative DT neutron count
        if d.dt_neutron_count is not None:
            dt_neutrons = d.dt_neutron_count.copy()
            if dt_neutrons.ndim >= 2:
                # Zone-resolved → sum over zones at each time
                dt_neutrons = np.sum(dt_neutrons, axis=1)
        else:
            dt_neutrons = nan_col.copy()

        # Hot-spot ρR (T_ion > 4.5 keV definition)
        if d.hot_spot_rhoR_vs_time is not None:
            hs_rhoR = d.hot_spot_rhoR_vs_time.copy()
        else:
            hs_rhoR = nan_col.copy()

        # ---- Write CSV ----
        columns = [
            ('time_ns',             d.time),
            ('laser_power_W',       laser_power),
            ('fusion_power_total',  fusion_total),
            ('hs_radius_cm',        hs_radius),
            ('total_rhoR_gcm2',     total_rhoR),
            ('fuel_rhoR_gcm2',      fuel_rhoR),
            ('hs_rhoR_gcm2',        hs_rhoR),
            ('peak_density_gcc',    peak_density),
            ('peak_ion_temp_keV',   peak_ion_temp),
            ('ablation_front_r_cm', abl_front),
            ('hs_pressure_Gbar',    hs_pressure),
            ('dt_neutron_count',    dt_neutrons),
        ]

        header = [name for name, _ in columns]
        arrays = [arr for _, arr in columns]

        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for i in range(n_times):
                row = []
                for arr in arrays:
                    val = arr[i]
                    if np.isnan(val):
                        row.append('')
                    else:
                        row.append(f'{val:.6e}')
                writer.writerow(row)

        logger.info(f"Time histories written: {path}  "
                    f"({n_times} rows × {len(header)} columns)")

    # ------------------------------------------------------------------
    # Helpers for time-history computation
    # ------------------------------------------------------------------

    def _compute_hs_radius_vs_time(self) -> np.ndarray:
        """Hot-spot outer radius from region interface indices."""
        d = self.data
        n_times = len(d.time)

        if (d.region_interfaces_indices is None
                or d.zone_boundaries is None):
            return np.full(n_times, np.nan)

        ri = d.region_interfaces_indices
        hs_bnd = ri[:, 0].astype(int)
        hs_radius = np.array([
            d.zone_boundaries[t, hs_bnd[t]]
            for t in range(n_times)
        ])
        return hs_radius

    def _compute_fuel_rhoR_vs_time(self) -> np.ndarray:
        """Cold-fuel areal density vs time (between HS and fuel/ablator boundary)."""
        d = self.data
        n_times = len(d.time)

        if (d.areal_density_vs_time is None
                or d.region_interfaces_indices is None):
            return np.full(n_times, np.nan)

        ri = d.region_interfaces_indices
        cumulative = d.areal_density_vs_time

        fuel_rhoR = np.zeros(n_times)
        for t in range(n_times):
            hs_bnd   = int(ri[t, 0])
            fuel_bnd = int(ri[t, -2])
            fuel_rhoR[t] = cumulative[t, fuel_bnd] - cumulative[t, hs_bnd]

        return fuel_rhoR

    def _compute_hs_pressure_vs_time(self) -> np.ndarray:
        """
        Mass-averaged hot-spot pressure (Gbar) at each timestep.

        Hot spot defined by region_interfaces_indices[:, 0] (the first
        region boundary — zones 0 to hs_bnd-1).
        """
        d = self.data
        n_times = len(d.time)

        if (d.ion_pressure is None or d.rad_pressure is None
                or d.zone_mass is None
                or d.region_interfaces_indices is None):
            return np.full(n_times, np.nan)

        ri = d.region_interfaces_indices
        hs_pressure = np.zeros(n_times)

        for t in range(n_times):
            hs_bnd = int(ri[t, 0])
            if hs_bnd <= 0:
                hs_pressure[t] = np.nan
                continue

            # Total pressure in Gbar for zones 0..hs_bnd-1
            p_total = (d.ion_pressure[t, :hs_bnd]
                       + d.rad_pressure[t, :hs_bnd]) * 1e-8
            mass = d.zone_mass[t, :hs_bnd]
            total_mass = np.sum(mass)

            if total_mass > 0:
                hs_pressure[t] = np.sum(p_total * mass) / total_mass
            else:
                hs_pressure[t] = np.nan

        return hs_pressure

    # ------------------------------------------------------------------
    # Formatting helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _metric(label: str, value, unit: str = '', fmt: str = '.4f') -> str:
        """Format a single metric line for the summary file."""
        if value is None or (isinstance(value, float) and value == 0.0):
            val_str = '—'
        else:
            val_str = f'{value:{fmt}}'
        if unit:
            val_str = f'{val_str} {unit}'
        return f'  {label:<36s} {val_str}'
