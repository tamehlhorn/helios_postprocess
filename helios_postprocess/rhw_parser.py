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
    drive_time: Optional[np.ndarray] = None
    drive_temperature: Optional[np.ndarray] = None
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
        
        # Extract drive temperature data
        drive_time, drive_temp = self._parse_drive_temperature(lines)
        
        config = RHWConfiguration(
            is_direct_drive=is_direct_drive,
            burn_enabled=burn_enabled,
            drive_time=drive_time,
            drive_temperature=drive_temp,
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
            laser_pulse_duration_ns=laser_params.get("pulse_duration_ns", 0.0)
        )
        
        logger.info(f"Configuration: {config.drive_type}, "
                   f"Burn {'ON' if burn_enabled else 'OFF'}")
        if drive_time is not None:
            logger.info(f"Drive temperature data: {len(drive_time)} points, "
                       f"{drive_time[0]:.3e} to {drive_time[-1]:.3e} s")
        
        return config
    
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
    
    def _parse_drive_temperature(self, lines: list) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Extract time-dependent drive temperature data.
        
        Parameters
        ----------
        lines : list
            Lines from RHW file
        
        Returns
        -------
        tuple
            (time_array, temperature_array) or (None, None) if not found
        """
        # Find the drive temperature section
        temp_start_idx = self._find_drive_temp_section(lines)
        if temp_start_idx is None:
            logger.warning("Could not find drive temperature data section")
            return None, None
        
        # Find the end marker
        temp_end_idx = self._find_drive_temp_end(lines, temp_start_idx)
        if temp_end_idx is None:
            logger.warning("Could not find end of drive temperature section")
            return None, None
        
        # Calculate number of lines in the data block
        n_lines = temp_end_idx - temp_start_idx
        if n_lines <= 0:
            logger.warning("Invalid drive temperature data block")
            return None, None
        
        try:
            # Extract time data (before the temperature section)
            time_text = ''.join(lines[temp_start_idx - n_lines + 1 : temp_start_idx])
            time_array = np.fromstring(time_text, sep=' ')
            
            # Extract temperature data (after the marker line)
            temp_text = ''.join(lines[temp_start_idx + 1 : temp_start_idx + n_lines])
            temp_array = np.fromstring(temp_text, sep=' ')
            
            if len(time_array) != len(temp_array):
                logger.warning(f"Time and temperature arrays have different lengths: "
                             f"{len(time_array)} vs {len(temp_array)}")
                # Try to use the shorter length
                min_len = min(len(time_array), len(temp_array))
                time_array = time_array[:min_len]
                temp_array = temp_array[:min_len]
            
            logger.debug(f"Parsed drive temperature: {len(time_array)} points")
            return time_array, temp_array
            
        except Exception as e:
            logger.error(f"Error parsing drive temperature data: {e}")
            return None, None
    
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
