"""
Parser for Helios .rhw (rad-hydro workspace) input files.

Extracts configuration information including:
- Drive type (direct drive vs indirect drive)
- Fusion reaction settings
- Time-dependent drive temperature data
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class RHWConfiguration:
    """Configuration data extracted from RHW file."""
    
    is_direct_drive: bool
    burn_enabled: bool
    # Laser ray-trace geometry (beam 1)
    laser_wavelength_um: float = 0.35
    laser_spot_size_cm: float = 0.0
    laser_half_cone_angle_deg: float = 1.0
    laser_focus_position_cm: float = 0.0
    laser_power_multiplier: float = 1.0
    laser_spatial_profile: str = "Uniform"
    laser_peak_power_TW: float = 0.0
    laser_foot_power_TW: float = 0.0
    laser_foot_start_ns: float = 0.0
    laser_foot_end_ns: float = 0.0
    laser_peak_start_ns: float = 0.0
    laser_peak_end_ns: float = 0.0
    laser_pulse_duration_ns: float = 0.0
    laser_geometry_per_beam: list = None      # [{beam_id, wavelength, spot_size, ...}, ...]
    eos_models: list = None  # [{region, type, file}, ...]
    alpha_deposition_local: bool = False     # Use alpha deposition = 1
    alpha_deposition_nonlocal: bool = False  # Use non alpha deposition = 1
    flux_limiter: float = 0.0                # First-region value (legacy; backward compat)
    flux_limiter_enabled: bool = False       # True if ANY region has flux limiter enabled
    flux_limiter_per_region: list = None     # [{region, enabled, value}, ...] or None
    drive_time: Optional[np.ndarray] = None
    drive_temperature: Optional[np.ndarray] = None
    drive_location: str = ""                  # "Rmin" or "Rmax"; empty if no radiation drive
    drive_flux_multiplier: float = 1.0        # scales effective sigma*T^4 (Helios convention)
    # Individual rad-source flags (May 2026 -- previously discarded after the
    # Rmin/Rmax winner was chosen). Used by data.n_active_idd_sources.
    rad_source_rmin_on: bool = False
    rad_source_rmax_on: bool = False
    source_file: Optional[str] = None
    
    @property
    def drive_type(self) -> str:
        """Human-readable drive type."""
        return "Direct Drive" if self.is_direct_drive else "Indirect Drive"


class RHWParser:
    """Parser for Helios .rhw input files."""
    
    def __init__(self, file_path: Path):
        """
        Initialize RHW parser.
        
        Parameters
        ----------
        file_path : Path
            Path to .rhw file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"RHW file not found: {file_path}")
    
    def parse(self) -> RHWConfiguration:
        """
        Parse RHW file and extract configuration.
        
        Returns
        -------
        RHWConfiguration
            Parsed configuration data
        """
        logger.info(f"Parsing RHW file: {self.file_path}")
        
        with open(self.file_path, 'r') as f:
            lines = f.readlines()
        
        # Extract configuration flags
        is_direct_drive = self._parse_direct_drive(lines)
        burn_enabled = self._parse_burn_status(lines)
        laser_params = self._parse_laser_geometry(lines)
        laser_geometry_per_beam = self._parse_laser_geometry_per_beam(lines)
        eos_models = self._parse_eos_models(lines)
        alpha_local, alpha_nonlocal = self._parse_alpha_transport(lines)
        flux_enabled, flux_value = self._parse_flux_limiter(lines)
        flux_per_region = self._parse_flux_limiter_per_region(lines)
        
        # Extract drive temperature data (radiation source from [Rad Source Data] block)
        drive_time, drive_temp, drive_location, drive_flux_mult = \
            self._parse_drive_temperature(lines)

        # Per-position rad-source flags (preserved separately so downstream
        # code can count active IDD sources without re-parsing).
        rmin_on, rmax_on = self._parse_rad_source_flags(lines)

        config = RHWConfiguration(
            is_direct_drive=is_direct_drive,
            burn_enabled=burn_enabled,
            drive_time=drive_time,
            drive_temperature=drive_temp,
            drive_location=drive_location,
            drive_flux_multiplier=drive_flux_mult,
            source_file=str(self.file_path),
            laser_wavelength_um=laser_params.get("wavelength", 0.35),
            laser_spot_size_cm=laser_params.get("spot_size", 0.0),
            laser_half_cone_angle_deg=laser_params.get("half_cone_angle", 1.0),
            laser_focus_position_cm=laser_params.get("focus_position", 0.0),
            laser_power_multiplier=laser_params.get("power_multiplier", 1.0),
            laser_spatial_profile=laser_params.get("spatial_profile", "Uniform"),
            laser_foot_power_TW=laser_params.get("foot_power_TW", 0.0),
            laser_peak_power_TW=laser_params.get("peak_power_TW", 0.0),
            laser_foot_start_ns=laser_params.get("foot_start_ns", 0.0),
            laser_foot_end_ns=laser_params.get("foot_end_ns", 0.0),
            laser_peak_start_ns=laser_params.get("peak_start_ns", 0.0),
            laser_peak_end_ns=laser_params.get("peak_end_ns", 0.0),
            laser_pulse_duration_ns=laser_params.get("pulse_duration_ns", 0.0),
            laser_geometry_per_beam=laser_geometry_per_beam if laser_geometry_per_beam else None,
            alpha_deposition_local=alpha_local,
            alpha_deposition_nonlocal=alpha_nonlocal,
            flux_limiter=flux_value,
            flux_limiter_enabled=flux_enabled,
            flux_limiter_per_region=flux_per_region if flux_per_region else None,
            rad_source_rmin_on=rmin_on,
            rad_source_rmax_on=rmax_on,
        )
        
        logger.info(f"Configuration: {config.drive_type}, "
                   f"Burn {'ON' if burn_enabled else 'OFF'}")
        if drive_time is not None:
            logger.info(f"Drive temperature data ({drive_location}): {len(drive_time)} points, "
                       f"{drive_time[0]*1e9:.3f} – {drive_time[-1]*1e9:.3f} ns, "
                       f"peak {np.max(drive_temp):.1f} eV, "
                       f"flux multiplier {drive_flux_mult:.3f}")
        
        return config
    
    def _parse_alpha_transport(self, lines: list) -> tuple:
        """Parse burn model from fusion transport and alpha deposition flags.
        
        Classification:
          Fusion transport = 0              -> No burn
          Fusion transport = 1, alpha dep=1 -> Local (instantaneous)
          Fusion transport = 1, alpha dep=0 -> Non-local transport
        """
        fusion_transport = False
        alpha_local = False
        for line in lines:
            s = line.strip()
            if s.startswith('Fusion transport on'):
                try:
                    if int(s.split('=')[-1].strip()) == 1:
                        fusion_transport = True
                except: pass
            elif s.startswith('Use alpha deposition'):
                try:
                    if int(s.split('=')[-1].strip()) == 1:
                        alpha_local = True
                except: pass
        if not fusion_transport:
            return False, False   # no burn
        if alpha_local:
            return True, False    # local instantaneous
        return False, True        # non-local transport
    
    def _parse_flux_limiter(self, lines: list) -> tuple:
        """
        Aggregate flux-limiter state across regions (backward-compat helper).

        Use _parse_flux_limiter_per_region for the full per-region distribution.
        This aggregate returns (any_enabled, first_value) — kept so downstream
        code that only reads scalar fields continues to work, but without the
        "varies across regions" warning since the per-region parser exposes
        the variation explicitly.

        Returns
        -------
        (enabled, value) : (bool, float)
            enabled = True if any region has flux limiter enabled.
            value   = first region's flux-limiter multiplier (legacy convention).
            (False, 0.0) if no flux-limiter line is present.
        """
        per_region = self._parse_flux_limiter_per_region(lines)
        if not per_region:
            return (False, 0.0)
        any_enabled = any(r['enabled'] for r in per_region)
        return (any_enabled, per_region[0]['value'])

    def _parse_eos_models(self, lines: list) -> list:
        """Parse EOS model type and file for each spatial region."""
        models = []
        current_region = None
        eos_type = None
        eos_file = None
        for line in lines:
            s = line.strip()
            if s.startswith('Parameters for Region ='):
                if current_region and eos_type is not None:
                    models.append({
                        'region': current_region,
                        'type': 'SESAME' if eos_type == 1 else 'PROPACEOS',
                        'file': eos_file or ''
                    })
                current_region = s.split('=', 1)[-1].strip()
                eos_type = None
                eos_file = None
            elif s.startswith('EOS data type'):
                try: eos_type = int(s.split('=')[-1].strip())
                except: pass
            elif s.startswith('EOS filepath'):
                eos_file = s.split('=', 1)[-1].strip()
        if current_region and eos_type is not None:
            models.append({
                'region': current_region,
                'type': 'SESAME' if eos_type == 1 else 'PROPACEOS',
                'file': eos_file or ''
            })
        return models if models else None

    def _parse_flux_limiter_per_region(self, lines: list) -> list:
        """
        Parse electron thermal flux limiter per spatial region.

        Mirrors _parse_eos_models region-aware pattern. Each material's
        flux-limiter pair sits inside its 'Parameters for Region =' block:

            Parameters for Region = DT vapor
              ...
              Use flux limiter       = 1
              Flux limiter mult.     = 0.06
            Parameters for Region = CD shell
              ...
              Use flux limiter       = 1
              Flux limiter mult.     = 0.02

        Returns
        -------
        list of dict
            One entry per region with 'region' (str), 'enabled' (bool),
            'value' (float). Empty list if no flux-limiter blocks found.
        """
        out = []
        current_region = None
        cur_enabled = None
        cur_value = None

        def _flush():
            nonlocal current_region, cur_enabled, cur_value
            if current_region is not None and cur_value is not None:
                out.append({
                    'region': current_region,
                    'enabled': bool(cur_enabled) if cur_enabled is not None else False,
                    'value': float(cur_value),
                })
            cur_enabled = None
            cur_value = None

        for line in lines:
            s = line.strip()
            lstrip = s.lower()
            if s.startswith('Parameters for Region ='):
                _flush()
                current_region = s.split('=', 1)[-1].strip()
            elif lstrip.startswith('use flux limiter'):
                try:
                    cur_enabled = (int(s.split('=')[-1].strip()) == 1)
                except (ValueError, IndexError):
                    pass
            elif lstrip.startswith('flux limiter mult'):
                try:
                    cur_value = float(s.split('=')[-1].strip())
                except (ValueError, IndexError):
                    pass
        _flush()
        return out

    def _parse_laser_geometry(self, lines: list) -> dict:
        """Parse beam-1 ray-trace parameters from [Laser Source Data] block."""
        import re
        params = {}
        in_laser = False
        beam1_found = False
        kv = re.compile(r'^\s*(.+?)\s*=\s*([^\s#]+)')
        for line in lines:
            if '[Laser Source Data]' in line:
                in_laser = True; continue
            if '[End Laser Source Data]' in line:
                break
            if not in_laser:
                continue
            m = kv.match(line)
            if m and 'Parameters for beam' in m.group(1):
                beam1_found = (m.group(2).strip() == '1'); continue
            if beam1_found and m and 'Parameters for beam' in m.group(1):
                break
            if not beam1_found or not m:
                continue
            key, val = m.group(1).strip(), m.group(2).strip()
            try:
                fval = float(val)
            except ValueError:
                continue
            lkey = key.lower()
            if 'wavelength' in lkey:
                params['wavelength'] = fval
            elif 'spot size' in lkey:
                params['spot_size'] = fval
            elif 'half cone' in lkey:
                params['half_cone_angle'] = fval
            elif 'focus position' in lkey:
                params['focus_position'] = fval
            elif 'power table multiplier' in lkey:
                params['power_multiplier'] = fval
            elif 'spatial profile model' in lkey:
                params['spatial_profile'] = 'Gaussian' if fval == 1.0 else 'Uniform'

        # ── Parse beam-1 power table ────────────────────────────────────────
        # Table format: times on one line, powers on next, both space-separated
        # Times are in seconds; convert to ns. Powers are in TW.
        in_laser = False
        beam1_found = False
        in_table = False
        n_rows = 0
        data_lines = []

        for line in lines:
            if '[Laser Source Data]' in line:
                in_laser = True; continue
            if '[End Laser Source Data]' in line:
                break
            if not in_laser:
                continue
            m2 = kv.match(line)
            if m2 and 'Parameters for beam' in m2.group(1):
                beam1_found = (m2.group(2).strip() == '1'); continue
            if beam1_found and m2 and 'Parameters for beam' in m2.group(1):
                break
            if not beam1_found:
                continue
            # Detect table start
            if '[table format' in line:
                in_table = True; data_lines = []; n_rows = 0; continue
            if not in_table:
                continue
            stripped = line.strip()
            if stripped.startswith('#') and 'table rows' in stripped:
                try: n_rows = int(stripped.split('=')[1].strip())
                except: pass
                continue
            # Skip header lines (start with letters or #)
            if stripped and not stripped[0].replace('-','').replace('.','').replace('e','').replace('E','').replace('+','').replace(' ','').isdigit():
                if not stripped.replace(' ','').replace('	','').replace('e','').replace('E','').replace('+','').replace('-','').replace('.','').isnumeric():
                    continue
            # Collect numeric data lines
            try:
                vals = [float(v) for v in stripped.split()]
                if len(vals) == n_rows and n_rows > 0:
                    data_lines.append(vals)
                    if len(data_lines) == 2:
                        break
            except ValueError:
                continue

        if len(data_lines) == 2 and n_rows > 0:
            import numpy as np
            times_ns = np.array(data_lines[0]) * 1e9   # s -> ns
            powers_TW = np.array(data_lines[1])         # already TW

            # Remove zero-power tail
            nonzero = powers_TW > 0
            if nonzero.any():
                params['pulse_duration_ns'] = float(times_ns[nonzero][-1])

            # Identify foot: first non-zero power level (below 50% of peak)
            peak_power = float(powers_TW.max())
            params['peak_power_TW'] = peak_power
            foot_mask = (powers_TW > 0) & (powers_TW < 0.5 * peak_power)
            if foot_mask.any():
                foot_power = float(powers_TW[foot_mask].mean())
                params['foot_power_TW'] = foot_power
                foot_times = times_ns[foot_mask]
                params['foot_start_ns'] = float(foot_times[0])
                params['foot_end_ns']   = float(foot_times[-1])

            # Peak: power >= 50% of max
            peak_mask = powers_TW >= 0.5 * peak_power
            if peak_mask.any():
                peak_times = times_ns[peak_mask]
                params['peak_start_ns'] = float(peak_times[0])
                params['peak_end_ns']   = float(peak_times[-1])

        return params
    
    def _parse_laser_geometry_per_beam(self, lines: list) -> list:
        """
        Parse ray-trace geometry for all beams in [Laser Source Data].

        Mirrors _parse_laser_geometry but loops over every 'Parameters
        for beam = N' block instead of stopping at beam 1. Returns
        geometry only (not the per-beam pulse table — that data comes
        from EXODUS via laser_power_delivered/on_target).

        Returns
        -------
        list of dict
            One entry per beam, in beam-number order. Each dict has:
            'beam_id' (int), 'wavelength' (um), 'spot_size' (cm),
            'half_cone_angle' (deg), 'focus_position' (cm),
            'power_multiplier' (float), 'spatial_profile' (str).
            Empty list if no laser block is found.
        """
        import re
        kv = re.compile(r'^\s*(.+?)\s*=\s*([^\s#]+)')
        beams = []
        current = None
        current_id = None

        def _flush():
            nonlocal current, current_id
            if current is not None and current_id is not None:
                current['beam_id'] = current_id
                beams.append(current)
            current = None
            current_id = None

        in_laser = False
        for line in lines:
            if '[Laser Source Data]' in line:
                in_laser = True
                continue
            if '[End Laser Source Data]' in line:
                _flush()
                break
            if not in_laser:
                continue
            m = kv.match(line)
            if m and 'Parameters for beam' in m.group(1):
                _flush()
                try:
                    current_id = int(m.group(2).strip())
                    current = {}
                except ValueError:
                    current_id = None
                    current = None
                continue
            if current is None or not m:
                continue
            key, val = m.group(1).strip(), m.group(2).strip()
            try:
                fval = float(val)
            except ValueError:
                continue
            lkey = key.lower()
            if 'wavelength' in lkey:
                current['wavelength'] = fval
            elif 'spot size' in lkey:
                current['spot_size'] = fval
            elif 'half cone' in lkey:
                current['half_cone_angle'] = fval
            elif 'focus position' in lkey:
                current['focus_position'] = fval
            elif 'power table multiplier' in lkey:
                current['power_multiplier'] = fval
            elif 'spatial profile model' in lkey:
                current['spatial_profile'] = 'Gaussian' if fval == 1.0 else 'Uniform'
        _flush()
        return beams
    
    def _parse_direct_drive(self, lines: list) -> bool:
        """
        Determine if direct drive (laser power model) is enabled.
        
        Parameters
        ----------
        lines : list
            Lines from RHW file
        
        Returns
        -------
        bool
            True if direct drive is enabled
        """
        for line in lines:
            if "Laser power model is on" in line:
                # Check if the value is 1 (enabled)
                is_enabled = "1" in line
                logger.debug(f"Direct drive status: {is_enabled} - {line.strip()}")
                return is_enabled
        
        logger.warning("Could not find 'Laser power model' flag, assuming indirect drive")
        return False
    
    def _parse_burn_status(self, lines: list) -> bool:
        """
        Determine if fusion reactions are enabled.
        
        Parameters
        ----------
        lines : list
            Lines from RHW file
        
        Returns
        -------
        bool
            True if fusion reactions are enabled
        """
        for line in lines:
            if "Fusion reactions on" in line:
                is_enabled = "1" in line
                logger.debug(f"Burn status: {is_enabled} - {line.strip()}")
                return is_enabled
        
        logger.warning("Could not find 'Fusion reactions' flag, assuming burn disabled")
        return False
    
    def _parse_rad_source_flags(self, lines: list) -> Tuple[bool, bool]:
        """Return (rmin_on, rmax_on) booleans from the [Rad Source Data] block.

        Separate from `_parse_drive_temperature` so callers can count
        active IDD sources without depending on the drive-temperature
        table being parseable.
        """
        rad_start = self._find_first_line(lines, "[Rad Source Data]")
        if rad_start is None:
            return False, False
        rad_end = self._find_first_line(lines, "[End Rad Source Data]",
                                        start=rad_start) or len(lines)
        block = lines[rad_start:rad_end]
        rmin = bool(self._extract_keyed_int(
            block, "Rad source model at Rmin is on", default=0))
        rmax = bool(self._extract_keyed_int(
            block, "Rad source model at Rmax is on", default=0))
        return rmin, rmax

    def _parse_drive_temperature(self, lines: list) -> Tuple[Optional[np.ndarray],
                                                              Optional[np.ndarray],
                                                              str, float]:
        """
        Parse the [Rad Source Data] block for time-dependent radiation drive.

        The block format (Helios Format ID = 4) is::

            [Rad Source Data]:
               Format ID = 4
               Rad source model at Rmin is on  = 0|1
               Rad source model at Rmax is on  = 0|1
               Flux multiplier                 = <float>
               ...
             [table format=2]:    Time-dependent Drive Temperatures at Rmin:
              # table rows = N_min
              ... (metadata)
              <N_min time values>     (whitespace-separated, may span lines)
              <N_min temperature values>
             [table format=2]:    Time-dependent Drive Temperatures at Rmax:
              # table rows = N_max
              ...
              <N_max time values>
              <N_max temperature values>
             [table format=2]:    Time-dependent Drive and Brightness Temperatures at Rmin:
               ...
            [End Rad Source Data]

        We pick whichever boundary has its "is on" flag set; if both are on,
        Rmax wins (the more common indirect-drive HDD configuration).

        Returns
        -------
        (time_s, temperature_eV, location, flux_multiplier)
            time_s          : (N,) seconds, or None if no active drive table.
            temperature_eV  : (N,) eV, raw table values (NOT pre-multiplied
                              by flux_multiplier — caller applies that).
            location        : "Rmin", "Rmax", or "" if no drive.
            flux_multiplier : float (default 1.0). Scales effective sigma*T^4
                              when computing radiation power flux on the
                              boundary.
        """
        rad_start = self._find_first_line(lines, "[Rad Source Data]")
        if rad_start is None:
            logger.debug("No [Rad Source Data] block found in RHW")
            return None, None, "", 1.0

        rad_end = self._find_first_line(lines, "[End Rad Source Data]",
                                        start=rad_start)
        if rad_end is None:
            rad_end = len(lines)

        block = lines[rad_start:rad_end]

        rmin_on = self._extract_keyed_int(block,
                                          "Rad source model at Rmin is on",
                                          default=0)
        rmax_on = self._extract_keyed_int(block,
                                          "Rad source model at Rmax is on",
                                          default=0)
        flux_mult = self._extract_keyed_float(block,
                                              "Flux multiplier",
                                              default=1.0)

        if not rmin_on and not rmax_on:
            logger.info("Radiation drive disabled at both Rmin and Rmax")
            return None, None, "", flux_mult

        if rmin_on and rmax_on:
            logger.warning("Both Rmin and Rmax radiation drives enabled; "
                           "using Rmax")
            location = "Rmax"
        else:
            location = "Rmax" if rmax_on else "Rmin"

        time_arr, temp_arr = self._parse_drive_temp_table(block, location)
        if time_arr is None or temp_arr is None:
            logger.warning(f"Drive temperature table at {location} is empty "
                           f"or unparseable")
            return None, None, location, flux_mult

        return time_arr, temp_arr, location, flux_mult

    def _parse_drive_temp_table(self, block: list,
                                 location: str) -> Tuple[Optional[np.ndarray],
                                                          Optional[np.ndarray]]:
        """
        Read the Time-dependent Drive Temperatures table at the given boundary.

        The block layout after the marker line is::

            [table format=2]:  Time-dependent Drive Temperatures at <loc>:
              table is 3D  = 0
              # table rows = N
              # table cols = 2
              interp model = 0
              column fmt id     = 0, 0
              column precision  = 3, 3
              <N time tokens, possibly across multiple lines>
              <N temperature tokens, possibly across multiple lines>
            [table format=2]: <next table>            ← terminates this block

        We find the next '[' line as the section terminator, then in two
        passes pull out N from the metadata and concatenate every numeric
        line into a single token stream.
        """
        marker = f"Time-dependent Drive Temperatures at {location}:"
        idx = None
        for i, line in enumerate(block):
            if marker in line and "Brightness" not in line:
                idx = i
                break
        if idx is None:
            return None, None

        # Find end of this table block: next '[' line (next table or [End)
        end_idx = len(block)
        for j in range(idx + 1, len(block)):
            if block[j].lstrip().startswith("["):
                end_idx = j
                break
        section = block[idx + 1:end_idx]

        # Pass 1: extract # table rows
        n_rows = None
        for line in section:
            if "# table rows" in line and "=" in line:
                try:
                    n_rows = int(line.split("=", 1)[1].strip())
                    break
                except (ValueError, IndexError):
                    continue
        if not n_rows or n_rows <= 0:
            return None, None

        # Pass 2: collect numeric data lines.
        # A data line is one whose first whitespace-separated token parses as
        # a float. This handles two RHW format variants:
        #   - format=2: bare numeric blocks (no labels between time and temp)
        #   - format=1: 'Time (sec)' and 'Drive Temperature (eV)' label lines
        #     between the two arrays — these must be skipped.
        data_lines = []
        for line in section:
            stripped = line.strip()
            if not stripped:
                continue
            if "=" in stripped or stripped.startswith("#"):
                continue                # metadata or comment
            try:
                float(stripped.split()[0])
            except (ValueError, IndexError):
                continue                # label line, e.g. "Time (sec)"
            data_lines.append(stripped)
        if not data_lines:
            return None, None

        try:
            values = np.fromstring(" ".join(data_lines), sep=" ")
        except Exception as e:
            logger.error(f"Could not parse drive temperature data at {location}: {e}")
            return None, None

        if len(values) < 2 * n_rows:
            logger.warning(f"Drive temperature table at {location} reports "
                           f"{n_rows} rows but only {len(values)} numeric "
                           f"values were found (need {2 * n_rows})")
            return None, None

        time_arr = values[:n_rows]
        temp_arr = values[n_rows:2 * n_rows]

        # Helios pads the drive table with sentinel times (~1e18 s) to mean
        # "this drive level extends indefinitely after the last real point."
        # Trim those so log output and downstream consumers see only physical times.
        SENTINEL_T_S = 1e-3       # 1 ms — well past any physical simulation duration
        keep = time_arr < SENTINEL_T_S
        if not keep.all():
            n_trim = int((~keep).sum())
            logger.info(f"Drive temperature at {location}: trimmed {n_trim} "
                        f"sentinel time value(s) (t ≥ {SENTINEL_T_S*1e9:.0f} ns)")
            time_arr = time_arr[keep]
            temp_arr = temp_arr[keep]

        return time_arr, temp_arr
    
    @staticmethod
    def _find_first_line(lines: list, marker: str,
                         start: int = 0) -> Optional[int]:
        """Index of the first line containing marker (substring), or None."""
        for i in range(start, len(lines)):
            if marker in lines[i]:
                return i
        return None

    @staticmethod
    def _extract_keyed_int(lines: list, key: str, default: int = 0) -> int:
        """Find first 'key ... = <int>' line; return parsed int or default."""
        for line in lines:
            if key in line and "=" in line:
                try:
                    return int(line.split("=", 1)[1].strip())
                except (ValueError, IndexError):
                    continue
        return default

    @staticmethod
    def _extract_keyed_float(lines: list, key: str,
                              default: float = 0.0) -> float:
        """Find first 'key ... = <float>' line; return parsed float or default."""
        for line in lines:
            if key in line and "=" in line:
                try:
                    return float(line.split("=", 1)[1].strip())
                except (ValueError, IndexError):
                    continue
        return default

    def _find_drive_temp_section(self, lines: list) -> Optional[int]:
        """Find the start of the drive temperature section."""
        # Try both possible capitalizations
        markers = [
            "Drive Temperature (eV)",
            "Drive temperature (eV)"
        ]
        
        for marker in markers:
            for i, line in enumerate(lines):
                if marker in line:
                    logger.debug(f"Found drive temperature marker at line {i}")
                    return i
        
        return None
    
    def _find_drive_temp_end(self, lines: list, start_idx: int) -> Optional[int]:
        """Find the end of the drive temperature section."""
        marker = "Time-dependent Drive and Brightness Temperatures at Rmin"
        
        for i in range(start_idx, len(lines)):
            if marker in lines[i]:
                logger.debug(f"Found drive temperature end marker at line {i}")
                return i
        
        return None


def load_rhw_configuration(file_path: Path) -> RHWConfiguration:
    """
    Convenience function to load RHW configuration.
    
    Parameters
    ----------
    file_path : Path
        Path to .rhw file
    
    Returns
    -------
    RHWConfiguration
        Parsed configuration
    """
    parser = RHWParser(file_path)
    return parser.parse()
