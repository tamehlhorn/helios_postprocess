"""
helios_postprocess.xray -- synthetic x-ray diagnostics from Helios runs.

Validation ladder:
    Level 0  this package, pure Python, free-free continuum only
    Level 1  SPECT3D continuum / LTE
    Level 2  SPECT3D full collisional-radiative non-LTE

Level 0 and the SPECT3D rungs land on the same ``RadianceCube`` contract
(t, impact parameter, photon energy) so a comparison is a diff, not a
re-implementation.

Typical use::

    from helios_postprocess.core import HeliosRun
    from helios_postprocess.data_builder import build_run_data
    from helios_postprocess.xray import (StreakConfig, build_radiance_cube,
                                         make_imaging_streak,
                                         make_spectral_streak,
                                         default_channels, burn_metrics)

    run = HeliosRun('run.exo')
    data = build_run_data(run)
    cfg = StreakConfig(t_start_ns=10.0, t_stop_ns=15.0, optically_thin=False)
    cube = build_radiance_cube(data, cfg)
    img = make_imaging_streak(cube, default_channels()['mid'], cfg)
    print(burn_metrics(img))
"""

from .emissivity import (EmissivityCube, build_emissivity, gaunt_ff, planck_E,
                         resolve_zbar, resolve_Te)
from .chords import (ChordResult, integrate, integrate_thin, integrate_formal,
                     make_impact_grid, spatially_integrate)
from .response import (Channel, Filter, Photocathode, HenkeTable,
                       load_henke_file, default_channels,
                       be_filter, al_filter, ch_filter, kapton_filter,
                       ti_filter, mu_over_rho_element, mu_over_rho_compound)
from .streak import (StreakConfig, RadianceCube, StreakImage,
                     BurnHistoryMetrics, build_radiance_cube,
                     make_imaging_streak, make_spectral_streak,
                     emission_edge_trajectory, emission_history,
                     burn_metrics, opacity_report)
from .spect3d_io import (export_profiles, load_spect3d_spectrum,
                         load_spect3d_cube, compare_cubes)

__all__ = [
    "EmissivityCube", "build_emissivity", "gaunt_ff", "planck_E",
    "resolve_zbar", "resolve_Te",
    "ChordResult", "integrate", "integrate_thin", "integrate_formal",
    "make_impact_grid", "spatially_integrate",
    "Channel", "Filter", "Photocathode", "HenkeTable", "load_henke_file",
    "default_channels", "be_filter", "al_filter", "ch_filter",
    "kapton_filter", "ti_filter", "mu_over_rho_element", "mu_over_rho_compound",
    "StreakConfig", "RadianceCube", "StreakImage", "BurnHistoryMetrics",
    "build_radiance_cube", "make_imaging_streak", "make_spectral_streak",
    "emission_edge_trajectory", "emission_history", "burn_metrics",
    "opacity_report",
    "export_profiles", "load_spect3d_spectrum", "load_spect3d_cube",
    "compare_cubes",
]

__version__ = "0.1.0"
