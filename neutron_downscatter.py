"""
Neutron Down-Scatter Ratio (DSR) Diagnostics for ICF
====================================================

The down-scatter ratio is a critical diagnostic for measuring fuel areal density.
Neutrons born at 14.1 MeV (D-T) lose energy via elastic scattering in the fuel.

DSR = N_scattered / N_primary

Where:
- N_primary: Neutrons detected at birth energy (14.1 MeV)
- N_scattered: Neutrons detected at lower energies (10-14 MeV)

Physics:
- More areal density → More scattering → Higher DSR
- DSR ∝ ρR^α where α ~ 0.5-1.0
- Typical DSR: 0.01-0.5 for ICF (1-50% down-scattered)

Applications:
- Validate areal density measurements
- Independent check of compression
- Fuel-ablator mix assessment
- Burn profile characterization

Energy ranges (D-T):
- Primary: 13.5-14.5 MeV (unscattered)
- Down-scatter: 10-13.5 MeV (scattered)
- Thermal: < 10 MeV (multiple scatters)

References:
- Murphy et al., Rev. Sci. Instrum. 85, 11D901 (2014) - NIF calibration
- Gatu Johnson et al., Phys. Plasmas 27, 042703 (2020) - OMEGA
- Frenje et al., Phys. Plasmas 17, 056311 (2010) - DSR physics

Author: Prof T
"""

import numpy as np
from typing import Optional, Tuple, Dict, List
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from scipy.signal import find_peaks


def calculate_downscatter_ratio(
    neutron_spectrum: np.ndarray,
    energy_bins: np.ndarray,
    E_primary_range: Tuple[float, float] = (13.5, 14.5),
    E_scatter_range: Tuple[float, float] = (10.0, 13.5)
) -> Dict[str, float]:
    """
    Calculate neutron down-scatter ratio from spectrum
    
    DSR = ∫(E_scatter) N(E) dE / ∫(E_primary) N(E) dE
    
    Parameters
    ----------
    neutron_spectrum : np.ndarray
        Neutron count spectrum, shape (n_bins,)
        Units: counts/MeV or neutrons/MeV
    energy_bins : np.ndarray
        Energy bin centers in MeV, shape (n_bins,)
    E_primary_range : Tuple[float, float], optional
        Energy range for primary peak [E_min, E_max] in MeV
        Default: (13.5, 14.5) for D-T
    E_scatter_range : Tuple[float, float], optional
        Energy range for down-scattered neutrons in MeV
        Default: (10.0, 13.5) for D-T
    
    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'DSR': Down-scatter ratio (dimensionless)
        - 'N_primary': Total primary neutrons
        - 'N_scatter': Total down-scattered neutrons
        - 'E_primary_mean': Mean energy of primary peak [MeV]
        - 'E_scatter_mean': Mean energy of scattered neutrons [MeV]
    
    Examples
    --------
    >>> energy = np.linspace(8, 16, 400)
    >>> spectrum = np.exp(-((energy - 14.1)**2) / (2 * 0.2**2))
    >>> result = calculate_downscatter_ratio(spectrum, energy)
    >>> print(f"DSR = {result['DSR']:.3f}")
    """
    # Input validation
    if len(neutron_spectrum) != len(energy_bins):
        raise ValueError("Spectrum and energy bins must have same length")
    
    if E_primary_range[0] >= E_primary_range[1]:
        raise ValueError("Primary energy range must be [E_min, E_max]")
    
    if E_scatter_range[0] >= E_scatter_range[1]:
        raise ValueError("Scatter energy range must be [E_min, E_max]")
    
    # Select primary peak region
    primary_mask = (energy_bins >= E_primary_range[0]) & (energy_bins <= E_primary_range[1])
    N_primary = trapezoid(neutron_spectrum[primary_mask], energy_bins[primary_mask])
    
    if N_primary == 0:
        raise ValueError("No neutrons detected in primary energy range")
    
    # Select down-scatter region
    scatter_mask = (energy_bins >= E_scatter_range[0]) & (energy_bins <= E_scatter_range[1])
    N_scatter = trapezoid(neutron_spectrum[scatter_mask], energy_bins[scatter_mask])
    
    # Calculate DSR
    DSR = N_scatter / N_primary
    
    # Calculate mean energies
    E_primary_mean = trapezoid(
        energy_bins[primary_mask] * neutron_spectrum[primary_mask],
        energy_bins[primary_mask]
    ) / N_primary
    
    if N_scatter > 0:
        E_scatter_mean = trapezoid(
            energy_bins[scatter_mask] * neutron_spectrum[scatter_mask],
            energy_bins[scatter_mask]
        ) / N_scatter
    else:
        E_scatter_mean = 0.0
    
    return {
        'DSR': DSR,
        'N_primary': N_primary,
        'N_scatter': N_scatter,
        'E_primary_mean': E_primary_mean,
        'E_scatter_mean': E_scatter_mean
    }


def downscatter_to_areal_density(
    DSR: float,
    calibration: str = 'NIF',
    custom_coeff: Optional[float] = None
) -> Dict[str, float]:
    """
    Convert down-scatter ratio to areal density
    
    Uses empirically-determined calibration relationships from major ICF facilities.
    
    Parameters
    ----------
    DSR : float
        Down-scatter ratio (dimensionless)
    calibration : str, optional
        Calibration to use: 'NIF', 'OMEGA', or 'custom'
        Default: 'NIF'
    custom_coeff : float, optional
        Custom calibration coefficient for ρR = coeff × DSR
        Required if calibration='custom'
    
    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'rho_R': Areal density [g/cm²]
        - 'uncertainty': Calibration uncertainty [g/cm²]
        - 'calibration': Calibration source used
    
    Notes
    -----
    Standard calibrations:
    - NIF: ρR = (20.4 ± 0.6) × DSR  (Murphy et al. 2014)
    - OMEGA: ρR = (19.4 ± 1.0) × DSR  (Gatu Johnson et al. 2020)
    
    The calibration depends on:
    - Energy range used for DSR measurement
    - Detector geometry and resolution
    - Target design (shell vs direct drive)
    
    Examples
    --------
    >>> result = downscatter_to_areal_density(DSR=0.15, calibration='NIF')
    >>> print(f"ρR = {result['rho_R']:.2f} ± {result['uncertainty']:.2f} g/cm²")
    """
    # Calibration coefficients and uncertainties
    calibrations = {
        'NIF': {'coeff': 20.4, 'uncertainty': 0.6},
        'OMEGA': {'coeff': 19.4, 'uncertainty': 1.0}
    }
    
    if calibration == 'custom':
        if custom_coeff is None:
            raise ValueError("custom_coeff required when calibration='custom'")
        coeff = custom_coeff
        uncertainty = 0.0  # User must specify their own
        calib_source = f'custom (coeff={custom_coeff:.2f})'
    elif calibration in calibrations:
        coeff = calibrations[calibration]['coeff']
        uncertainty = calibrations[calibration]['uncertainty']
        calib_source = calibration
    else:
        raise ValueError(f"Unknown calibration: {calibration}. Use 'NIF', 'OMEGA', or 'custom'")
    
    # Convert DSR to areal density
    rho_R = coeff * DSR
    
    # Propagate uncertainty (linear relationship)
    rho_R_uncertainty = uncertainty * DSR if uncertainty > 0 else 0.0
    
    return {
        'rho_R': rho_R,
        'uncertainty': rho_R_uncertainty,
        'calibration': calib_source
    }


def areal_density_to_downscatter(
    rho_R: float,
    temperature: float,
    fuel_Z: float = 1.5,
    fuel_A: float = 2.5
) -> Dict[str, float]:
    """
    Predict down-scatter ratio from areal density using physics model
    
    Uses elastic scattering cross-sections to predict DSR for a given
    areal density and fuel temperature.
    
    Parameters
    ----------
    rho_R : float
        Areal density [g/cm²]
    temperature : float
        Ion temperature [keV]
    fuel_Z : float, optional
        Average atomic number of fuel (default: 1.5 for D-T)
    fuel_A : float, optional
        Average atomic mass of fuel (default: 2.5 for D-T)
    
    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'DSR_predicted': Predicted down-scatter ratio
        - 'scatter_fraction_10_12': Fraction scattered into 10-12 MeV
        - 'scatter_fraction_12_13p5': Fraction scattered into 12-13.5 MeV
        - 'mean_free_path': Neutron mean free path [g/cm²]
    
    Notes
    -----
    Model assumptions:
    - Single elastic scattering approximation (valid for DSR < 0.5)
    - Isotropic scattering in center-of-mass frame
    - Uniform fuel composition
    
    For high areal densities (ρR > 1 g/cm²), multiple scattering
    becomes important and this model underestimates DSR.
    
    Examples
    --------
    >>> result = areal_density_to_downscatter(rho_R=0.30, temperature=5.0)
    >>> print(f"Predicted DSR = {result['DSR_predicted']:.3f}")
    """
    # Neutron scattering cross-section (approximate for D-T fuel)
    # σ_elastic ≈ 2 barns for 14 MeV neutrons on D-T
    sigma_elastic_barn = 2.0  # barns
    sigma_elastic_cm2 = sigma_elastic_barn * 1e-24  # cm²
    
    # Number density [particles/cm³]
    # n = ρ / (A × m_u) where m_u = 1.66e-24 g
    m_u = 1.66e-24  # atomic mass unit in grams
    
    # Mean free path λ = 1 / (n × σ)
    # In areal density units: λ_ρR = A × m_u / σ
    lambda_rhoR = (fuel_A * m_u) / sigma_elastic_cm2  # g/cm²
    
    # Fraction of neutrons that scatter once
    # P(scatter) = 1 - exp(-ρR / λ) ≈ ρR / λ for small ρR/λ
    scatter_probability = 1.0 - np.exp(-rho_R / lambda_rhoR)
    
    # Energy loss distribution for elastic scattering
    # For D-T: Maximum energy transfer to deuterium (A=2)
    # E_final / E_initial = 1 - 2/(A+1)² × (1 - cos(θ))
    # For isotropic scattering, this gives a distribution
    
    # Approximate fractions for different energy ranges (from kinematic calculations)
    # These are rough estimates - real calculation requires full kinematics
    fraction_10_12 = 0.30  # ~30% of scattered neutrons fall in 10-12 MeV
    fraction_12_13p5 = 0.45  # ~45% fall in 12-13.5 MeV
    
    # Total DSR (using 10-13.5 MeV range)
    DSR_predicted = scatter_probability * (fraction_10_12 + fraction_12_13p5)
    
    return {
        'DSR_predicted': DSR_predicted,
        'scatter_fraction_10_12': scatter_probability * fraction_10_12,
        'scatter_fraction_12_13p5': scatter_probability * fraction_12_13p5,
        'mean_free_path': lambda_rhoR
    }


def simulate_neutron_spectrum(
    rho_R: float,
    temperature: float,
    n_primary: float = 1e12,
    detector_resolution: float = 0.1,
    energy_range: Tuple[float, float] = (8.0, 16.0),
    n_bins: int = 400
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate neutron energy spectrum including down-scatter
    
    Generates a synthetic neutron spectrum with primary peak and
    down-scattered component based on areal density and temperature.
    
    Parameters
    ----------
    rho_R : float
        Areal density [g/cm²]
    temperature : float
        Ion temperature [keV]
    n_primary : float, optional
        Total number of primary (14.1 MeV) neutrons
        Default: 1e12
    detector_resolution : float, optional
        Detector energy resolution (FWHM) in MeV
        Default: 0.1 MeV
    energy_range : Tuple[float, float], optional
        Energy range [E_min, E_max] in MeV
        Default: (8.0, 16.0)
    n_bins : int, optional
        Number of energy bins
        Default: 400
    
    Returns
    -------
    spectrum : np.ndarray
        Neutron counts in each energy bin
    energy_bins : np.ndarray
        Energy bin centers [MeV]
    
    Notes
    -----
    The spectrum includes:
    1. Primary peak at 14.1 MeV (Gaussian broadened)
    2. Doppler broadening from ion temperature (ΔE ~ sqrt(T))
    3. Down-scattered component (10-13.5 MeV)
    4. Detector resolution broadening
    
    Examples
    --------
    >>> spectrum, energy = simulate_neutron_spectrum(
    ...     rho_R=0.35, temperature=5.0, n_primary=1e12
    ... )
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(energy, spectrum)
    >>> plt.xlabel('Energy [MeV]')
    >>> plt.ylabel('Counts')
    >>> plt.show()
    """
    # Create energy bins
    energy_bins = np.linspace(energy_range[0], energy_range[1], n_bins)
    dE = energy_bins[1] - energy_bins[0]
    
    # Primary neutron energy (D-T fusion)
    E_primary = 14.1  # MeV
    
    # Doppler broadening from ion temperature
    # ΔE/E ≈ sqrt(T/1000) for D-T
    doppler_width = E_primary * np.sqrt(temperature / 1000.0)
    
    # Total broadening (Doppler + detector)
    # Quadrature sum: σ_total² = σ_doppler² + σ_detector²
    sigma_doppler = doppler_width / 2.355  # FWHM to σ
    sigma_detector = detector_resolution / 2.355
    sigma_total = np.sqrt(sigma_doppler**2 + sigma_detector**2)
    
    # Primary peak (Gaussian)
    primary_peak = n_primary * np.exp(-((energy_bins - E_primary)**2) / (2 * sigma_total**2))
    primary_peak /= (sigma_total * np.sqrt(2 * np.pi))  # Normalize
    
    # Down-scattered component
    # Use physics model to get DSR
    dsr_result = areal_density_to_downscatter(rho_R, temperature)
    n_scattered = n_primary * dsr_result['DSR_predicted']
    
    # Scattered neutron distribution (approximate as uniform in 10-13.5 MeV)
    scatter_mask = (energy_bins >= 10.0) & (energy_bins <= 13.5)
    n_scatter_bins = np.sum(scatter_mask)
    
    scattered_component = np.zeros_like(energy_bins)
    if n_scatter_bins > 0:
        # Uniform distribution in scatter range
        scattered_component[scatter_mask] = n_scattered / n_scatter_bins
        
        # Apply detector broadening to scattered component
        # (Convolve with Gaussian - simplified as direct broadening)
        for i, E in enumerate(energy_bins):
            if scatter_mask[i]:
                gaussian = np.exp(-((energy_bins - E)**2) / (2 * sigma_detector**2))
                gaussian /= (sigma_detector * np.sqrt(2 * np.pi))
                scattered_component += scattered_component[i] * gaussian * dE
    
    # Total spectrum
    spectrum = primary_peak + scattered_component
    
    return spectrum, energy_bins


def compute_scattering_cross_section(
    E_neutron: float,
    target_A: float,
    scattering_type: str = 'elastic'
) -> float:
    """
    Compute neutron scattering cross-section
    
    Calculates elastic or inelastic scattering cross-section for
    neutrons on light nuclei (relevant for D-T fuel).
    
    Parameters
    ----------
    E_neutron : float
        Neutron energy [MeV]
    target_A : float
        Target nucleus mass number (1=H, 2=D, 3=T)
    scattering_type : str, optional
        Type of scattering: 'elastic' or 'inelastic'
        Default: 'elastic'
    
    Returns
    -------
    float
        Cross-section [barns]
    
    Notes
    -----
    Uses empirical parameterizations for:
    - n-D elastic scattering
    - n-T elastic scattering
    - (n,2n) reactions
    
    Valid for neutron energies 1-20 MeV.
    
    References
    ----------
    ENDF/B-VIII.0 nuclear data evaluations
    
    Examples
    --------
    >>> sigma = compute_scattering_cross_section(E_neutron=14.1, target_A=2)
    >>> print(f"n-D elastic: {sigma:.2f} barns")
    """
    if scattering_type == 'elastic':
        # Elastic scattering cross-sections (approximate)
        # Based on ENDF/B data for 14 MeV neutrons
        
        if abs(target_A - 2.0) < 0.1:  # Deuterium
            # n-D elastic at 14 MeV ≈ 0.6 barns
            # Weak energy dependence in 10-20 MeV range
            sigma = 0.6 * (14.0 / E_neutron)**0.1
            
        elif abs(target_A - 3.0) < 0.1:  # Tritium
            # n-T elastic at 14 MeV ≈ 1.2 barns
            sigma = 1.2 * (14.0 / E_neutron)**0.1
            
        elif abs(target_A - 1.0) < 0.1:  # Hydrogen
            # n-p elastic at 14 MeV ≈ 0.7 barns
            sigma = 0.7 * (14.0 / E_neutron)**0.05
            
        else:
            # Generic light nucleus (A < 10)
            # Rough estimate: σ ∝ A^(2/3)
            sigma = 2.0 * (target_A / 2.5)**(2/3) * (14.0 / E_neutron)**0.1
    
    elif scattering_type == 'inelastic':
        # (n,2n) reactions have threshold around 3-6 MeV
        # Cross-sections peak around 14 MeV
        
        if E_neutron < 3.0:
            sigma = 0.0  # Below threshold
        elif abs(target_A - 2.0) < 0.1:  # Deuterium
            # D(n,2n) very small
            sigma = 0.05 * (E_neutron / 14.0)**2
        elif abs(target_A - 3.0) < 0.1:  # Tritium
            # T(n,2n) cross-section
            sigma = 0.2 * np.exp(-((E_neutron - 14.0)**2) / 20.0)
        else:
            sigma = 0.1 * (E_neutron / 14.0)**1.5
    
    else:
        raise ValueError(f"Unknown scattering_type: {scattering_type}")
    
    return sigma


def extract_spectrum_peaks(
    neutron_spectrum: np.ndarray,
    energy_bins: np.ndarray,
    prominence: float = 0.1,
    width: int = 5
) -> Dict[str, np.ndarray]:
    """
    Automatically detect primary and scattered peaks in spectrum
    
    Uses peak-finding algorithm to identify the primary (14.1 MeV)
    peak and any significant down-scattered features.
    
    Parameters
    ----------
    neutron_spectrum : np.ndarray
        Neutron count spectrum
    energy_bins : np.ndarray
        Energy bin centers [MeV]
    prominence : float, optional
        Minimum peak prominence (fraction of max)
        Default: 0.1 (10% of maximum)
    width : int, optional
        Minimum peak width in bins
        Default: 5
    
    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary containing:
        - 'peak_energies': Energy of detected peaks [MeV]
        - 'peak_heights': Height of detected peaks
        - 'peak_widths': FWHM of detected peaks [MeV]
        - 'primary_idx': Index of primary peak (highest energy peak)
    
    Examples
    --------
    >>> spectrum, energy = simulate_neutron_spectrum(rho_R=0.3, temperature=5.0)
    >>> peaks = extract_spectrum_peaks(spectrum, energy)
    >>> print(f"Primary at {peaks['peak_energies'][peaks['primary_idx']]:.2f} MeV")
    """
    # Find peaks using scipy
    max_height = np.max(neutron_spectrum)
    min_prominence = prominence * max_height
    
    peak_indices, properties = find_peaks(
        neutron_spectrum,
        prominence=min_prominence,
        width=width
    )
    
    if len(peak_indices) == 0:
        return {
            'peak_energies': np.array([]),
            'peak_heights': np.array([]),
            'peak_widths': np.array([]),
            'primary_idx': None
        }
    
    # Extract peak energies and heights
    peak_energies = energy_bins[peak_indices]
    peak_heights = neutron_spectrum[peak_indices]
    
    # Estimate peak widths (FWHM)
    peak_widths_bins = properties['widths']
    dE = energy_bins[1] - energy_bins[0]
    peak_widths = peak_widths_bins * dE
    
    # Identify primary peak (should be highest energy peak)
    primary_idx = np.argmax(peak_energies)
    
    return {
        'peak_energies': peak_energies,
        'peak_heights': peak_heights,
        'peak_widths': peak_widths,
        'primary_idx': primary_idx
    }


def calculate_downscatter_asymmetry(
    DSR_measurements: List[float],
    detector_angles: List[float]
) -> Dict[str, float]:
    """
    Analyze asymmetry in down-scatter ratio from multiple lines of sight
    
    Multi-angle DSR measurements can reveal fuel asymmetries and
    non-uniform compression.
    
    Parameters
    ----------
    DSR_measurements : List[float]
        DSR values from different detectors
    detector_angles : List[float]
        Detector angles [degrees] relative to laser axis
        (e.g., [0, 45, 90, 135] for 4 detectors)
    
    Returns
    -------
    Dict[str, float]
        Dictionary containing:
        - 'DSR_mean': Average DSR across all angles
        - 'DSR_std': Standard deviation of DSR
        - 'asymmetry_percent': Asymmetry as percentage (100 × std/mean)
        - 'max_min_ratio': Ratio of max to min DSR
        - 'rho_R_mean': Mean areal density [g/cm²]
        - 'rho_R_std': Std dev of areal density [g/cm²]
    
    Notes
    -----
    Asymmetry metrics:
    - < 10%: Excellent symmetry
    - 10-20%: Good symmetry
    - 20-30%: Moderate asymmetry
    - > 30%: Poor symmetry (likely degraded performance)
    
    Examples
    --------
    >>> DSR = [0.15, 0.16, 0.14, 0.15]  # 4 detectors
    >>> angles = [0, 90, 180, 270]
    >>> result = calculate_downscatter_asymmetry(DSR, angles)
    >>> print(f"Asymmetry: {result['asymmetry_percent']:.1f}%")
    """
    DSR_array = np.array(DSR_measurements)
    
    if len(DSR_array) < 2:
        raise ValueError("Need at least 2 DSR measurements for asymmetry analysis")
    
    # Calculate DSR statistics
    DSR_mean = np.mean(DSR_array)
    DSR_std = np.std(DSR_array, ddof=1)  # Sample std dev
    asymmetry_percent = 100.0 * DSR_std / DSR_mean if DSR_mean > 0 else 0.0
    
    max_min_ratio = np.max(DSR_array) / np.min(DSR_array) if np.min(DSR_array) > 0 else np.inf
    
    # Convert to areal density (using NIF calibration)
    rho_R_array = np.array([
        downscatter_to_areal_density(dsr, calibration='NIF')['rho_R']
        for dsr in DSR_array
    ])
    
    rho_R_mean = np.mean(rho_R_array)
    rho_R_std = np.std(rho_R_array, ddof=1)
    
    return {
        'DSR_mean': DSR_mean,
        'DSR_std': DSR_std,
        'asymmetry_percent': asymmetry_percent,
        'max_min_ratio': max_min_ratio,
        'rho_R_mean': rho_R_mean,
        'rho_R_std': rho_R_std
    }
