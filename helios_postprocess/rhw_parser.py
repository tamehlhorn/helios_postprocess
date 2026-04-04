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
        
        # Extract drive temperature data
        drive_time, drive_temp = self._parse_drive_temperature(lines)
        
        config = RHWConfiguration(
            is_direct_drive=is_direct_drive,
            burn_enabled=burn_enabled,
            drive_time=drive_time,
            drive_temperature=drive_temp,
            source_file=str(self.file_path)
        )
        
        logger.info(f"Configuration: {config.drive_type}, "
                   f"Burn {'ON' if burn_enabled else 'OFF'}")
        if drive_time is not None:
            logger.info(f"Drive temperature data: {len(drive_time)} points, "
                       f"{drive_time[0]:.3e} to {drive_time[-1]:.3e} s")
        
        return config
    
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
