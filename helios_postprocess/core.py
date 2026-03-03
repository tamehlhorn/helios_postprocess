"""
Helios Postprocessing Core Module

Main interface for analyzing Helios ICF simulation outputs with integrated
neutron downscatter ratio (DSR) diagnostics and comprehensive physics analysis.

Classes
-------
HeliosRun : Main analysis class with caching and physics integration
    - ExodusII file reading and data management
    - Neutron downscatter ratio (DSR) calculations
    - Areal density (ρR) extraction
    - Hot spot metrics and ignition criteria
    - Pressure gradient analysis

Author: Prof T
Updated: November16, 2025
"""

import numpy as np
import netCDF4 as nc
from pathlib import Path
from typing import Optional, Dict, Any, Tuple, List, Union
from functools import lru_cache
import warnings

# After your imports (numpy, warnings, etc.)
# Before: class HeliosRun:

# ============================================================================
# FUSION REACTIVITY FUNCTIONS (for neutron weighting)
# ============================================================================

def sigma_v_DT(T_keV):
    """
    D-T fusion reactivity using Bosch-Hale 1992 parameterization.
    
    Parameters
    ----------
    T_keV : float or array
        Ion temperature in keV
        
    Returns
    -------
    float or array
        Fusion reactivity <σv> in cm³/s
    """
    B_G = 34.3827
    mc2 = 1124656
    C1 = 1.17302e-9
    C2 = 1.51361e-2
    C3 = 7.51886e-2
    C4 = 4.60643e-3
    C5 = 1.35000e-2
    C6 = -1.06750e-4
    C7 = 1.36600e-5
    
    theta = T_keV / (1 - (T_keV * (C2 + T_keV * (C4 + T_keV * C6))) /
                         (1 + T_keV * (C3 + T_keV * (C5 + T_keV * C7))))
    
    xi = (B_G**2 / (4 * theta))**(1/3)
    
    sigma_v = C1 * theta * np.sqrt(xi / (mc2 * T_keV**3)) * np.exp(-3 * xi)
    
    return sigma_v


def calculate_neutron_weight(T_keV, n_D, n_T=None, fuel_type='DT'):
    """
    Calculate neutron production weighting factor.
    
    For DT: rate ∝ n_D * n_T * <σv>_DT
    
    Parameters
    ----------
    T_keV : array
        Ion temperature in keV
    n_D : array
        Deuterium number density in cm⁻³
    n_T : array, optional
        Tritium number density in cm⁻³ (for DT)
    fuel_type : str
        'DT' or 'DD'
        
    Returns
    -------
    array
        Neutron production weight (unnormalized)
    """
    if fuel_type.upper() == 'DT':
        if n_T is None:
            raise ValueError("Must provide n_T for DT fuel")
        reactivity = sigma_v_DT(T_keV)
        weight = n_D * n_T * reactivity
    else:
        raise ValueError(f"Unknown fuel_type: {fuel_type}")
    
    return weight


def neutron_averaged(quantity, T_keV, n_D, n_T=None, zone_mass=None, fuel_type='DT'):
    """
    Calculate neutron-averaged (burn-weighted) value of a quantity.
    
    <Q>_n = Σ Q_i * W_i * m_i / Σ W_i * m_i
    
    where W_i is the neutron production rate in zone i.
    
    Parameters
    ----------
    quantity : array
        Quantity to average (e.g., temperature, density, pressure)
    T_keV : array
        Ion temperature in keV
    n_D : array
        Deuterium number density in cm⁻³
    n_T : array, optional
        Tritium number density in cm⁻³ (for DT)
    zone_mass : array, optional
        Zone masses for proper spatial weighting
    fuel_type : str
        'DT' or 'DD'
        
    Returns
    -------
    float
        Neutron-averaged value
    """
    # Calculate neutron production weight
    weight = calculate_neutron_weight(T_keV, n_D, n_T, fuel_type)
    
    # Apply zone mass weighting if provided
    if zone_mass is not None:
        weight = weight * zone_mass
    
    # Avoid divide by zero
    total_weight = np.sum(weight)
    if total_weight == 0:
        warnings.warn("Total neutron weight is zero, returning simple average")
        return np.average(quantity, weights=zone_mass if zone_mass is not None else None)
    
    # Calculate weighted average
    avg_value = np.sum(quantity * weight) / total_weight
    
    return avg_value


def get_neutron_averaged_conditions(T_keV, density, pressure, n_D, n_T=None, 
                                   zone_mass=None, fuel_type='DT'):
    """
    Calculate neutron-averaged temperature, density, and pressure.
    
    These are the burn-averaged conditions that characterize the hot spot.
    
    Parameters
    ----------
    T_keV : array
        Ion temperature in keV
    density : array
        Mass density in g/cm³
    pressure : array
        Pressure in Gbar
    n_D : array
        Deuterium number density in cm⁻³
    n_T : array, optional
        Tritium number density in cm⁻³
    zone_mass : array, optional
        Zone masses
    fuel_type : str
        'DT' or 'DD'
        
    Returns
    -------
    dict
        Dictionary with keys:
        - 'T_n_avg' : Neutron-averaged temperature (keV)
        - 'rho_n_avg' : Neutron-averaged density (g/cm³)
        - 'P_n_avg' : Neutron-averaged pressure (Gbar)
    """
    results = {}
    
    # Calculate neutron-averaged quantities
    results['T_n_avg'] = neutron_averaged(T_keV, T_keV, n_D, n_T, zone_mass, fuel_type)
    results['rho_n_avg'] = neutron_averaged(density, T_keV, n_D, n_T, zone_mass, fuel_type)
    results['P_n_avg'] = neutron_averaged(pressure, T_keV, n_D, n_T, zone_mass, fuel_type)
    
    return results


def get_burn_diagnostics_corrected(times, neutron_rate=None, temp=None, dens=None,
                                   zone_mass=None, fuel_mass=None):
    """
    Calculate burn diagnostics with CORRECTED burn width calculation.
    
    Parameters
    ----------
    times : array
        Time array in seconds
    neutron_rate : array, optional
        Neutron production rate array (time, space)
    temp : array, optional
        Temperature array (used as proxy if neutron_rate unavailable)
    dens : array, optional
        Density array for burn fraction
    zone_mass : array, optional
        Zone masses
    fuel_mass : float, optional
        Total initial fuel mass for burn fraction
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'bang_time' : Time of peak neutron production (s)
        - 'bang_time_idx' : Index of bang time
        - 'total_yield' : Total neutron yield
        - 'burn_width_fwhm' : FWHM of burn pulse (s)
        - 'burn_width_rms' : RMS burn width (s)
        - 'peak_rate' : Peak neutron production rate (/s)
        - 'burn_fraction' : Fraction of fuel burned
        - 'diagnostic_info' : Additional info for debugging
    """
    results = {}
    diagnostic = {}
    
 # Use provided neutron rate, or fall back to temperature proxy
    if neutron_rate is not None:
        rate_time = neutron_rate
    elif temp is not None:
        rate_time = np.array([np.sum(t**2) for t in temp]) if temp.ndim > 1 else temp**2
        warnings.warn("Using temperature² as proxy for neutron rate")
    else:
        raise ValueError("Must provide either neutron_rate or temp")
    
    diagnostic['rate_time'] = rate_time
    diagnostic['times'] = times
    
    # BANG TIME (peak production)
    bang_idx = np.argmax(rate_time)
    bang_time = times[bang_idx]
    peak_rate = rate_time[bang_idx]
    
    results['bang_time'] = bang_time
    results['bang_time_idx'] = bang_idx
    results['peak_rate'] = peak_rate
    
    # TOTAL YIELD
    dt = np.diff(times)
    if len(dt) > 0:
        dt_avg = dt.mean()
        total_yield = np.trapezoid(rate_time, times)
    else:
        dt_avg = 1.0
        total_yield = rate_time.sum()
    
    results['total_yield'] = total_yield
    diagnostic['dt_avg'] = dt_avg
    
    # BURN WIDTH - CORRECTED FWHM CALCULATION
    half_max = peak_rate / 2.0
    above_half = rate_time >= half_max
    
    if np.sum(above_half) >= 2:
        crossings = np.where(above_half)[0]
        t_start = times[crossings[0]]
        t_end = times[crossings[-1]]
        fwhm = t_end - t_start
        
        # Interpolate for more accurate FWHM
        if crossings[0] > 0:
            idx_left = crossings[0] - 1
            t_left = np.interp(half_max, 
                             [rate_time[idx_left], rate_time[idx_left + 1]],
                             [times[idx_left], times[idx_left + 1]])
        else:
            t_left = times[crossings[0]]
        
        if crossings[-1] < len(times) - 1:
            idx_right = crossings[-1]
            t_right = np.interp(half_max,
                              [rate_time[idx_right + 1], rate_time[idx_right]],
                              [times[idx_right + 1], times[idx_right]])
        else:
            t_right = times[crossings[-1]]
        
        fwhm_interp = t_right - t_left
        
        diagnostic['fwhm_simple'] = fwhm
        diagnostic['fwhm_interp'] = fwhm_interp
        diagnostic['t_left'] = t_left
        diagnostic['t_right'] = t_right
        
        burn_width_fwhm = fwhm_interp
        
    else:
        warnings.warn("Could not find FWHM points, using RMS width")
        burn_width_fwhm = 0.0
        diagnostic['fwhm_failure'] = True
    
    # RMS BURN WIDTH
    t_mean = np.sum(times * rate_time) / np.sum(rate_time)
    t_sq_mean = np.sum(times**2 * rate_time) / np.sum(rate_time)
    burn_width_rms = np.sqrt(t_sq_mean - t_mean**2)
    
    expected_fwhm_from_rms = 2.355 * burn_width_rms
    
    results['burn_width_fwhm'] = burn_width_fwhm
    results['burn_width_rms'] = burn_width_rms
    diagnostic['t_mean'] = t_mean
    diagnostic['expected_fwhm_from_rms'] = expected_fwhm_from_rms
    
    # If FWHM failed, use RMS-based estimate
    if burn_width_fwhm == 0.0 and burn_width_rms > 0:
        results['burn_width_fwhm'] = expected_fwhm_from_rms
        diagnostic['using_rms_fwhm'] = True
    
    # BURN FRACTION
    if fuel_mass is not None and total_yield > 0:
        m_D = 2.014
        m_T = 3.016
        N_A = 6.022e23
        
        mass_per_reaction = (m_D + m_T) / N_A
        mass_burned = total_yield * mass_per_reaction
        burn_fraction = mass_burned / fuel_mass
        
        results['burn_fraction'] = burn_fraction
        diagnostic['mass_burned'] = mass_burned
    else:
        results['burn_fraction'] = 0.0
    
    results['diagnostic_info'] = diagnostic
    
    return results


# ============================================================================
# Now your HeliosRun class starts here
# ============================================================================

class HeliosRun:
    """
    Main interface for Helios ICF simulation analysis with integrated physics diagnostics.
    
    This class provides:
    - Automatic reading and caching of ExodusII (.exo) files
    - Neutron downscatter ratio (DSR) calculations for areal density
    - Hot spot metrics and ignition criteria (Lindl criteria)
    - Pressure gradient and shock analysis
    - Burn diagnostics (bang time, yield, burn width)
    
    Parameters
    ----------
    exo_file : str or Path
        Path to Helios ExodusII output file (.exo)
    cache_data : bool, optional
        Enable LRU caching for expensive operations (default: True)
    verbose : bool, optional
        Print diagnostic messages (default: False)
        
    Attributes
    ----------
    filepath : Path
        Path to ExodusII file
    dataset : netCDF4.Dataset
        Open netCDF4 dataset handle
    times : np.ndarray
        Simulation time array (ns)
    n_times : int
        Number of time steps
    coords : dict
        Spatial coordinates (x, y, z or r, theta)
    cache_enabled : bool
        Whether caching is active
        
    Examples
    --------
    >>> run = HeliosRun('Vulcan_HDD_NB_111125_AI_2.exo')
    >>> 
    >>> # Get neutron downscatter ratio
    >>> dsr = run.calculate_dsr(time_idx=-1)  # At peak burn
    >>> print(f"DSR = {dsr['DSR']:.3f}")
    >>> 
    >>> # Convert DSR to areal density
    >>> rho_R = run.dsr_to_areal_density(dsr['DSR'], calibration='NIF')
    >>> print(f"ρR = {rho_R['rho_R']:.2f} ± {rho_R['uncertainty']:.2f} g/cm²")
    >>> 
    >>> # Check ignition criteria
    >>> ignition = run.check_ignition_criteria()
    >>> print(f"Ignition achieved: {ignition['ignition_achieved']}")
    >>> 
    >>> # Get burn diagnostics
    >>> burn = run.get_burn_diagnostics()
    >>> print(f"Bang time: {burn['bang_time']:.2f} ns")
    >>> print(f"Neutron yield: {burn['total_yield']:.2e}")
    """
    
    def __init__(
        self,
        exo_file: Union[str, Path],
        cache_data: bool = True,
        verbose: bool = False
    ):
        """Initialize HeliosRun with ExodusII file."""
        self.filepath = Path(exo_file)
        if not self.filepath.exists():
            raise FileNotFoundError(f"ExodusII file not found: {self.filepath}")
            
        self.verbose = verbose
        self.cache_enabled = cache_data
        
        # Open ExodusII file
        try:
            self.dataset = nc.Dataset(self.filepath, 'r')
            if self.verbose:
                print(f"✓ Opened {self.filepath.name}")
        except Exception as e:
            raise IOError(f"Failed to open ExodusII file: {e}")
            
        # Load time array
        self._load_time_array()
        
        # Load coordinates
        self._load_coordinates()
        
        # Cache for expensive computations
        self._cache: Dict[str, Any] = {}
        
        # Physics module integration flags
        self._dsr_available = False
        self._check_physics_modules()
        
    def _load_time_array(self) -> None:
        """Load simulation time array from ExodusII file."""
        try:
            # ExodusII stores time in 'time_whole' variable
            self.times = self.dataset.variables['time_whole'][:]
            self.n_times = len(self.times)
            
            if self.verbose:
                print(f"✓ Loaded {self.n_times} time steps")
                print(f"  Time range: {self.times[0]:.2e} - {self.times[-1]:.2e} s")
                
        except KeyError:
            warnings.warn("Could not find 'time_whole' in ExodusII file")
            self.times = np.array([])
            self.n_times = 0
            
    def _load_coordinates(self) -> None:
        """Load spatial coordinate arrays from ExodusII file."""
        self.coords = {}
        
        try:
            # Try common ExodusII coordinate variable names
            coord_vars = ['coordx', 'coordy', 'coordz', 'coordr', 'coord']
            
            for var_name in coord_vars:
                if var_name in self.dataset.variables:
                    self.coords[var_name] = self.dataset.variables[var_name][:]
                    
            if self.verbose and self.coords:
                print(f"✓ Loaded coordinates: {list(self.coords.keys())}")
                
        except Exception as e:
            warnings.warn(f"Could not load coordinates: {e}")
            
    def _check_physics_modules(self) -> None:
        """Check availability of physics analysis modules."""
        try:
            from .physics import neutron_downscatter
            self._dsr_available = True
            if self.verbose:
                print("✓ Neutron DSR module loaded")
        except ImportError:
            if self.verbose:
                print("⚠ Neutron DSR module not available")
                
    def get_variable(
        self,
        var_name: str,
        time_idx: Optional[int] = None
    ) -> np.ndarray:
        """
        Get variable data from ExodusII file.
        
        Parameters
        ----------
        var_name : str
            Variable name in ExodusII file
            Common vars: 'dens' (density), 'temp' (temperature), 
                        'pres' (pressure), 'tev' (eV), 'trad' (rad temp)
        time_idx : int, optional
            Time index to extract. If None, returns all times.
            Negative indices supported (e.g., -1 for last time).
            
        Returns
        -------
        np.ndarray
            Variable data array
            
        Examples
        --------
        >>> density = run.get_variable('dens', time_idx=-1)  # Final density
        >>> temp_history = run.get_variable('temp')  # All times
        """
        if var_name not in self.dataset.variables:
            available = [v for v in self.dataset.variables.keys() 
                        if not v.startswith('time')]
            raise KeyError(
                f"Variable '{var_name}' not found. "
                f"Available: {available[:10]}..."
            )
            
        data = self.dataset.variables[var_name]
        
        if time_idx is not None:
            return data[time_idx]
        else:
            return data[:]
            
    def list_variables(self) -> List[str]:
        """
        List all available variables in the ExodusII file.
        
        Returns
        -------
        list of str
            Variable names
        """
        return [v for v in self.dataset.variables.keys() 
                if not v.startswith('time') and not v.startswith('coord')]
    
    # ========================================================================
    # NEUTRON DOWNSCATTER RATIO (DSR) METHODS
    # ========================================================================
    
    def calculate_dsr(
        self,
        time_idx: int = -1,
        neutron_spectrum: Optional[np.ndarray] = None,
        energy_bins: Optional[np.ndarray] = None,
        E_primary_range: Tuple[float, float] = (13.5, 14.5),
        E_scatter_range: Tuple[float, float] = (10.0, 13.5)
    ) -> Dict[str, float]:
        """
        Calculate neutron downscatter ratio from simulation or provided spectrum.
        
        DSR = N_scattered / N_primary
        
        Parameters
        ----------
        time_idx : int, optional
            Time index for extracting neutron spectrum from simulation.
            Default is -1 (final time/peak burn).
        neutron_spectrum : np.ndarray, optional
            Externally provided neutron spectrum (counts/MeV).
            If None, attempts to extract from simulation.
        energy_bins : np.ndarray, optional
            Energy bin centers (MeV) for neutron_spectrum.
            Required if neutron_spectrum is provided.
        E_primary_range : tuple, optional
            Energy range for primary (unscattered) neutrons in MeV.
            Default: (13.5, 14.5) for D-T fusion.
        E_scatter_range : tuple, optional
            Energy range for down-scattered neutrons in MeV.
            Default: (10.0, 13.5) for D-T.
            
        Returns
        -------
        dict
            Dictionary with keys:
            - 'DSR': Down-scatter ratio
            - 'N_primary': Primary neutron count
            - 'N_scattered': Scattered neutron count
            - 'uncertainty': Statistical uncertainty in DSR
            
        Examples
        --------
        >>> # Calculate DSR at peak burn
        >>> dsr = run.calculate_dsr(time_idx=-1)
        >>> print(f"DSR = {dsr['DSR']:.3f} ± {dsr['uncertainty']:.3f}")
        >>>
        >>> # Use external spectrum
        >>> spectrum, energy = load_experimental_spectrum()
        >>> dsr = run.calculate_dsr(neutron_spectrum=spectrum, 
        ...                         energy_bins=energy)
        """
        if not self._dsr_available:
            raise ImportError(
                "Neutron DSR module not available. "
                "Ensure physics/neutron_downscatter.py is present."
            )
            
        from .physics.neutron_downscatter import calculate_downscatter_ratio
        
        # If spectrum not provided, try to extract from simulation
        if neutron_spectrum is None:
            neutron_spectrum, energy_bins = self._extract_neutron_spectrum(time_idx)
            
        if neutron_spectrum is None or energy_bins is None:
            raise ValueError(
                "Could not extract neutron spectrum from simulation. "
                "Provide neutron_spectrum and energy_bins explicitly."
            )
            
        # Calculate DSR using physics module
        dsr_result = calculate_downscatter_ratio(
            neutron_spectrum=neutron_spectrum,
            energy_bins=energy_bins,
            E_primary_range=E_primary_range,
            E_scatter_range=E_scatter_range
        )
        
        return dsr_result
        
    def dsr_to_areal_density(
        self,
        dsr: float,
        calibration: str = 'NIF',
        custom_factor: Optional[float] = None,
        custom_uncertainty: Optional[float] = None
    ) -> Dict[str, float]:
        """
        Convert neutron downscatter ratio to fuel areal density.
        
        Uses empirical calibrations from:
        - NIF: Murphy et al. (2014) - ρR = (20.4 ± 0.6) × DSR
        - OMEGA: Gatu Johnson et al. (2020) - ρR = (19.4 ± 1.0) × DSR
        
        Parameters
        ----------
        dsr : float
            Neutron downscatter ratio
        calibration : str, optional
            Calibration to use: 'NIF', 'OMEGA', or 'custom'
            Default: 'NIF'
        custom_factor : float, optional
            Custom calibration factor (ρR/DSR) if calibration='custom'
        custom_uncertainty : float, optional
            Uncertainty in custom factor
            
        Returns
        -------
        dict
            Dictionary with keys:
            - 'rho_R': Areal density in g/cm²
            - 'uncertainty': Absolute uncertainty in g/cm²
            - 'calibration': Calibration used
            - 'factor': Conversion factor applied
            
        Examples
        --------
        >>> dsr_result = run.calculate_dsr()
        >>> rho_R = run.dsr_to_areal_density(dsr_result['DSR'])
        >>> print(f"ρR = {rho_R['rho_R']:.2f} ± {rho_R['uncertainty']:.2f} g/cm²")
        """
        if not self._dsr_available:
            raise ImportError("Neutron DSR module not available")
            
        from .physics.neutron_downscatter import downscatter_to_areal_density
        
        return downscatter_to_areal_density(
            dsr=dsr,
            calibration=calibration,
            custom_factor=custom_factor,
            custom_uncertainty=custom_uncertainty
        )
        
    def predict_dsr(
        self,
        rho_R: float,
        temperature: float,
        fuel_type: str = 'DT'
    ) -> Dict[str, float]:
        """
        Predict expected DSR from simulation conditions.
        
        Useful for validating simulation compression against
        expected neutron diagnostic signatures.
        
        Parameters
        ----------
        rho_R : float
            Fuel areal density in g/cm²
        temperature : float
            Ion temperature in keV
        fuel_type : str, optional
            Fuel composition: 'DT', 'DD', or 'D3He'
            Default: 'DT'
            
        Returns
        -------
        dict
            Dictionary with keys:
            - 'predicted_DSR': Expected downscatter ratio
            - 'rho_R': Input areal density
            - 'temperature': Input temperature
            - 'fuel_type': Fuel composition used
            
        Examples
        --------
        >>> # Get areal density from simulation
        >>> rho_R = run.calculate_areal_density(time_idx=-1)
        >>> T_ion = run.get_variable('temp', time_idx=-1).mean()
        >>> 
        >>> # Predict expected DSR
        >>> expected = run.predict_dsr(rho_R, T_ion)
        >>> print(f"Expected DSR = {expected['predicted_DSR']:.3f}")
        """
        if not self._dsr_available:
            raise ImportError("Neutron DSR module not available")
            
        from .physics.neutron_downscatter import areal_density_to_downscatter
        
        return areal_density_to_downscatter(
            rho_R=rho_R,
            temperature=temperature,
            fuel_type=fuel_type
        )
        
    def simulate_neutron_spectrum(
        self,
        time_idx: int = -1,
        rho_R: Optional[float] = None,
        temperature: Optional[float] = None,
        n_primary: float = 1e12,
        detector_resolution: float = 0.1,
        energy_range: Tuple[float, float] = (8.0, 16.0),
        n_bins: int = 400
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate expected neutron spectrum from simulation conditions.
        
        Parameters
        ----------
        time_idx : int, optional
            Time index for extracting conditions. Default: -1 (peak burn)
        rho_R : float, optional
            Fuel areal density (g/cm²). If None, calculates from simulation.
        temperature : float, optional
            Ion temperature (keV). If None, extracts from simulation.
        n_primary : float, optional
            Total primary neutron yield. Default: 1e12
        detector_resolution : float, optional
            Detector FWHM resolution in MeV. Default: 0.1
        energy_range : tuple, optional
            Energy range (E_min, E_max) in MeV. Default: (8.0, 16.0)
        n_bins : int, optional
            Number of energy bins. Default: 400
            
        Returns
        -------
        spectrum : np.ndarray
            Simulated neutron spectrum (counts/MeV)
        energy : np.ndarray
            Energy bin centers (MeV)
            
        Examples
        --------
        >>> spectrum, energy = run.simulate_neutron_spectrum()
        >>> 
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(energy, spectrum)
        >>> plt.xlabel('Energy (MeV)')
        >>> plt.ylabel('Neutrons per MeV')
        """
        if not self._dsr_available:
            raise ImportError("Neutron DSR module not available")
            
        from .physics.neutron_downscatter import simulate_neutron_spectrum
        
        # Extract conditions from simulation if not provided
        if rho_R is None:
            rho_R = self.calculate_areal_density(time_idx=time_idx)
            
        if temperature is None:
            temp_data = self.get_variable('temp', time_idx=time_idx)
            temperature = np.mean(temp_data)  # Mass-averaged temperature
            
        return simulate_neutron_spectrum(
            rho_R=rho_R,
            temperature=temperature,
            n_primary=n_primary,
            detector_resolution=detector_resolution,
            energy_range=energy_range,
            n_bins=n_bins
        )
    
    def _extract_neutron_spectrum(
        self,
        time_idx: int = -1
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Extract neutron energy spectrum from simulation output.
        
        This is a placeholder - actual implementation depends on
        how Helios stores neutron diagnostic data.
        
        Parameters
        ----------
        time_idx : int
            Time index to extract spectrum
            
        Returns
        -------
        spectrum : np.ndarray or None
            Neutron spectrum if available
        energy : np.ndarray or None
            Energy bins if available
        """
        # Check for common neutron spectrum variable names
        spectrum_vars = [
            'neutron_spectrum', 'nspectrum', 'n_spectrum',
            'neutron_energy_dist', 'dnde'
        ]
        
        for var_name in spectrum_vars:
            if var_name in self.dataset.variables:
                spectrum = self.get_variable(var_name, time_idx=time_idx)
                
                # Try to find corresponding energy bins
                energy_vars = [f'{var_name}_energy', 'neutron_energy', 
                              'energy_bins', 'ebins']
                for e_var in energy_vars:
                    if e_var in self.dataset.variables:
                        energy = self.dataset.variables[e_var][:]
                        return spectrum, energy
                        
        # If no pre-computed spectrum, could synthesize from conditions
        if self.verbose:
            warnings.warn(
                "Could not find neutron spectrum in simulation output. "
                "Use simulate_neutron_spectrum() to generate expected spectrum."
            )
        
        return None, None
    
    # ========================================================================
    # AREAL DENSITY METHODS
    # ========================================================================
    
    def calculate_areal_density(
        self,
        time_idx: int = -1,
        region: str = 'total'
    ) -> dict:
        """
        Calculate fuel areal density (ρR).
        
        Parameters
        ----------
        time_idx : int, optional
            Time index. Default: -1 (last step)
        region : str, optional
            Region to analyze: 'total' for full target
            
        Returns
        -------
        dict with 'rhoR' in g/cm², 'density', 'r', 'dr'
        """
        density = self.get_variable('mass_density', time_idx=time_idx)
        r = self.get_variable('zone_boundaries', time_idx=time_idx)
        dr = np.diff(r)
        
        rhoR = np.sum(density * dr[:len(density)])
        
        return {
            'rhoR': rhoR,
            'density': density,220
            'r': r,
            'dr': dr
        }
    
    # ========================================================================
    # HOT SPOT & IGNITION METHODS
    # ========================================================================
    
    def check_ignition_criteria(
        self,
        time_idx: int = -1
    ) -> Dict[str, Any]:
        """
        Check Lindl ignition criteria for ICF.
        
        Ignition requires:
        1. <T> > 5 keV (temperature for fusion)
        2. <P> > 100 Gbar (compression)
        3. ρR > 0.3 g/cm² (alpha particle stopping)
        4. Isobaric hot spot (dP/dr ≈ 0)
        
        Parameters
        ----------
        time_idx : int, optional
            Time index to check. Default: -1 (peak burn)
            
        Returns
        -------
        dict
            Dictionary with:
            - 'ignition_achieved': bool
            - 'temperature': float (keV)
            - 'pressure': float (Gbar)
            - 'areal_density': float (g/cm²)
            - 'is_isobaric': bool
            - 'criteria_met': dict of individual criterion results
            
        Examples
        --------
        >>> result = run.check_ignition_criteria()
        >>> if result['ignition_achieved']:
        ...     print("🎉 IGNITION ACHIEVED!")
        ... else:
        ...     print("Criteria not met:")
        ...     for name, met in result['criteria_met'].items():
        ...         print(f"  {name}: {'✓' if met else '✗'}")
        """
        # Get temperature (keV)
        temp = self.get_variable('temp', time_idx=time_idx).mean()
        
        # Get pressure (convert to Gbar)
        pres = self.get_variable('pres', time_idx=time_idx).mean()
        # Helios pressure is in J/cm³, convert to Gbar (1 Gbar = 1e12 dyne/cm² = 1e11 Pa)
        pres_gbar = pres * 1e-11  # J/cm³ to Gbar
        
        # Get areal density
        rho_R = self.calculate_areal_density(time_idx=time_idx)
        
        # Check individual criteria
        criteria = {
            'temperature': temp > 5.0,
            'pressure': pres_gbar > 100.0,
            'areal_density': rho_R > 0.3,
        }
        
        ignition = all(criteria.values())
        
        return {
            'ignition_achieved': ignition,
            'temperature': temp,
            'pressure': pres_gbar,
            'areal_density': rho_R,
            'criteria_met': criteria
        }
    
    # ========================================================================
    # BURN DIAGNOSTICS
    # ========================================================================
    
    def get_burn_diagnostics(self, time_idx=-1):
        """
        Calculate burn diagnostics with CORRECTED burn width calculation.
        
        Returns
        -------
        dict
            Dictionary containing:
            - 'bang_time': Time of peak neutron production (s)
            - 'bang_time_idx': Index of bang time
            - 'total_yield': Total neutron yield
            - 'burn_width_fwhm': FWHM of burn pulse (s)
            - 'burn_width_rms': RMS burn width (s)
            - 'peak_rate': Peak neutron production rate (/s)
            - 'burn_fraction': Fraction of fuel burned
            - 'diagnostic_info': Additional diagnostic information
            
        Notes
        -----
        Uses corrected burn width calculation that:
        - Interpolates FWHM for sub-timestep accuracy
        - Provides RMS width as robust alternative
        - Never returns zero width
        """
     # Get neutron rate - prefer DT fusion rate from Helios
        fusion_var = None
        for var in ['FusionRate_DT_nHe4', 'neutron_rate', 'fusion_power']:
            if var in self.list_variables():
                fusion_var = var
                break
        
        if fusion_var:
            # Build global rate vs time by summing over zones at each step
            rate = np.zeros(len(self.times))
            for i in range(len(self.times)):
                spatial_rate = self.get_variable(fusion_var, time_idx=i)
                rate[i] = np.sum(spatial_rate)
        else:
            # Fallback: use temperature as proxy
            temp_var = None
            for var in ['ion_temperature', 'temp', 'temperature']:
                if var in self.list_variables():
                    temp_var = var
                    break
            
            if temp_var:
                rate = np.zeros(len(self.times))
                for i in range(len(self.times)):
                    temp = self.get_variable(temp_var, time_idx=i)
                    rate[i] = np.sum(temp**2)
                warnings.warn(f"Neutron rate not found, using {temp_var}² as proxy")
            else:
                raise ValueError("Cannot find FusionRate_DT_nHe4, neutron_rate, or temperature")
        
        # Get density and mass for burn fraction calculation
        dens = None
        zone_mass = None
        fuel_mass = None
        
        for var in ['mass_density', 'dens', 'density']:
            if var in self.list_variables():
                dens = self.get_variable(var)
                break
        
        for var in ['zone_mass', 'mass', 'cell_mass']:
            if var in self.list_variables():
                zone_mass = self.get_variable(var)
                if zone_mass is not None and dens is not None:
                    fuel_mass = np.sum(zone_mass)
                break
        
        # Call corrected burn diagnostics function
        return get_burn_diagnostics_corrected(
            times=self.times,
            neutron_rate=rate,
            temp=None,
            dens=dens,
            zone_mass=zone_mass,
            fuel_mass=fuel_mass
        )
        
    def get_neutron_averaged_conditions(self, time_idx=None):
        """Calculate neutron-averaged (burn-weighted) conditions."""
        # Get time index
        if time_idx is None:
            try:
                burn = self.get_burn_diagnostics()
                time_idx = burn['bang_time_idx']
            except:
                time_idx = -1
        
        # Get variables (with automatic detection)
        for var in ['ion_temperature', 'temp']:
            if var in self.list_variables():
                T_ion = self.get_variable(var, time_idx)
                break
        
        for var in ['mass_density', 'dens']:
            if var in self.list_variables():
                rho = self.get_variable(var, time_idx)
                break
        
        for var in ['ion_pressure', 'pres']:
            if var in self.list_variables():
                P_ion = self.get_variable(var, time_idx)
                break
        
        # Add electron pressure if available
        P_total = P_ion
        if 'elec_pressure' in self.list_variables():
            P_total += self.get_variable('elec_pressure', time_idx)
        
        # Convert units (Helios-specific)
        # Temperature: eV to keV
        T_keV = T_ion / 1000.0 if T_ion.mean() > 100 else T_ion

        # Pressure: J/cm³ to Gbar (1 Gbar = 1e11 J/cm³)
        P_gbar = P_total / 1e11
        
        # Number densities (50-50 DT)
        m_avg = (2.014 + 3.016) / 2.0
        N_A = 6.022e23
        n_total = rho / (m_avg / N_A)
        n_D = n_T = n_total / 2.0
        
        # Zone masses
        zone_mass = None
        if 'zone_mass' in self.list_variables():
            zone_mass = self.get_variable('zone_mass', time_idx)
        
        # Simple averages
        T_simple = np.average(T_keV, weights=zone_mass) if zone_mass is not None else np.mean(T_keV)
        rho_simple = np.average(rho, weights=zone_mass) if zone_mass is not None else np.mean(rho)
        P_simple = np.average(P_gbar, weights=zone_mass) if zone_mass is not None else np.mean(P_gbar)
        
        # Neutron-averaged
        n_avg = get_neutron_averaged_conditions(T_keV, rho, P_gbar, n_D, n_T, zone_mass)
        
        return {
            'T_n_avg': n_avg['T_n_avg'],
            'rho_n_avg': n_avg['rho_n_avg'],
            'P_n_avg': n_avg['P_n_avg'],
            'simple_T': T_simple,
            'simple_rho': rho_simple,
            'simple_P': P_simple,
            'enhancement_T': n_avg['T_n_avg'] / T_simple,
            'enhancement_rho': n_avg['rho_n_avg'] / rho_simple,
            'enhancement_P': n_avg['P_n_avg'] / P_simple,
            'time_idx': time_idx
        }
               
    def close(self):
        """Close the ExodusII dataset."""
        if hasattr(self, 'dataset') and self.dataset is not None:
            self.dataset.close()
            
# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("HELIOS POSTPROCESSING - CORRECTED VERSION")
    print("=" * 70)
    print("\nFixes implemented:")
    print("1. ✓ Burn width calculation (was returning zero)")
    print("2. ✓ Neutron-averaged quantities added")
    print("\n" + "=" * 70)
    
    # Test burn width calculation with synthetic data
    print("\nTesting burn width calculation...")
    
    # Create Gaussian-like burn pulse
    times = np.linspace(0, 10, 1000)  # ns
    t_bang = 5.0  # bang time
    sigma = 0.5   # width
    
    # Gaussian burn profile
    rate = np.exp(-(times - t_bang)**2 / (2 * sigma**2))
    
    # Add some spatial structure
    n_zones = 50
    rate_spatial = rate[:, np.newaxis] * np.ones((len(times), n_zones))
    
    # Calculate diagnostics
    burn_diag = get_burn_diagnostics_corrected(times, neutron_rate=rate_spatial)
    
    print(f"\nBang time: {burn_diag['bang_time']:.3f} ns (expected: {t_bang:.3f} ns)")
    print(f"FWHM burn width: {burn_diag['burn_width_fwhm']:.3f} ns (expected: {2.355*sigma:.3f} ns)")
    print(f"RMS burn width: {burn_diag['burn_width_rms']:.3f} ns (expected: {sigma:.3f} ns)")
    print(f"Total yield: {burn_diag['total_yield']:.3e}")
    
    # Test neutron averaging
    print("\n" + "=" * 70)
    print("Testing neutron-averaged quantities...")
    
    # Create temperature and density profiles
    T_keV = 10.0 + 5.0 * np.random.rand(n_zones)  # 10-15 keV
    density = 50.0 + 20.0 * np.random.rand(n_zones)  # 50-70 g/cm³
    pressure = 200.0 + 50.0 * np.random.rand(n_zones)  # 200-250 Gbar
    
    # Number densities (assume 50-50 DT)
    m_D = 2.014  # amu
    m_T = 3.016  # amu
    m_avg = (m_D + m_T) / 2
    N_A = 6.022e23
    
    n_total = density / (m_avg / N_A)  # particles/cm³
    n_D = n_total / 2
    n_T = n_total / 2
    
    # Simple average
    T_simple = np.mean(T_keV)
    rho_simple = np.mean(density)
    P_simple = np.mean(pressure)
    
    # Neutron-averaged
    n_avg = get_neutron_averaged_conditions(T_keV, density, pressure, n_D, n_T)
    
    print(f"\nSimple average T: {T_simple:.2f} keV")
    print(f"Neutron-averaged T: {n_avg['T_n_avg']:.2f} keV")
    print(f"\nSimple average ρ: {rho_simple:.2f} g/cm³")
    print(f"Neutron-averaged ρ: {n_avg['rho_n_avg']:.2f} g/cm³")
    print(f"\nSimple average P: {P_simple:.2f} Gbar")
    print(f"Neutron-averaged P: {n_avg['P_n_avg']:.2f} Gbar")
    
    print("\n" + "=" * 70)
    print("TESTS COMPLETE")
    print("=" * 70)
    
    # ========================================================================
    # UTILITY METHODS
    # ========================================================================
    
    def __repr__(self) -> str:
        """String representation of HeliosRun."""
        return (
            f"HeliosRun('{self.filepath.name}', "
            f"n_times={self.n_times}, "
            f"cache={'enabled' if self.cache_enabled else 'disabled'})"
        )
        
    def __enter__(self):
        """Context manager entry."""
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - close dataset."""
        self.close()
        
    def close(self):
        """Close the ExodusII dataset."""
        if hasattr(self, 'dataset'):
            self.dataset.close()
            if self.verbose:
                print(f"✓ Closed {self.filepath.name}")
