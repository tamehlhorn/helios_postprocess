"""
Helios Postprocessing Package
=============================

Analysis pipeline for Helios ICF simulation outputs (EXODUS/netCDF4).

Pipeline
--------
    HeliosRun  →  build_run_data  →  ICFRunData
                                        ↓
                                  ICFAnalyzer(data)
                                    .analyze_drive_phase()
                                    .analyze_stagnation_phase()
                                    .analyze_burn_phase()
                                    .compute_performance_metrics()
                                        ↓
                                  ICFPlotter(data, config)
                                    .create_full_report(pdf_path)
                                        ↓
                                  ICFOutputGenerator(data)
                                    .write_all(base_path)

Additional Physics Modules
--------------------------
    energetics           — Kinetic energy, hydro efficiency, PdV work
    neutron_downscatter  — Down-scatter ratio (DSR) diagnostics
    pressure_gradients   — Shock identification, RT instability assessment
    burn_averaged_metrics — Temporal burn-averaging, published-data comparison

Quick Start
-----------
>>> from helios_postprocess import HeliosRun
>>> from helios_postprocess.data_builder import build_run_data
>>> from helios_postprocess.icf_analysis import ICFAnalyzer
>>> from helios_postprocess.icf_plotting import ICFPlotter
>>> from helios_postprocess.icf_output import ICFOutputGenerator
>>>
>>> run = HeliosRun('simulation.exo', verbose=True)
>>> data = build_run_data(run, time_unit='s')
>>> run.close()
>>>
>>> analyzer = ICFAnalyzer(data)
>>> analyzer.analyze_drive_phase()
>>> analyzer.analyze_stagnation_phase()
>>> analyzer.analyze_burn_phase()
>>> analyzer.compute_performance_metrics()
>>>
>>> ICFPlotter(data, {}).create_full_report('report.pdf')
>>> ICFOutputGenerator(data).write_all('simulation')

Version: 3.0.0 (March 2026)
Author: Prof T / Xcimer ICF Analysis
"""

# Core EXODUS reader
from .core import HeliosRun

# Data bridge
from .data_builder import build_run_data, ICFRunData

# Analysis pipeline (class-based)
from .icf_analysis import ICFAnalyzer
from .icf_plotting import ICFPlotter
from .icf_output import ICFOutputGenerator

# Burn-averaged metrics (functional — works on ICFRunData after analysis)
from .burn_averaged_metrics import (
    extract_histories_from_run_data,
    extract_hot_spot_histories,      # legacy compat wrapper
    calculate_burn_averaged_metrics,
    compare_with_published,
)

# Additional physics modules (functional)
from . import energetics
from . import neutron_downscatter
from . import pressure_gradients
from . import adiabat_history

# RHW parser (optional — only needed if .rhw files are available)
try:
    from .rhw_parser import RHWParser, RHWConfiguration, load_rhw_configuration
except ImportError:
    pass

__version__ = "3.0.0"

__all__ = [
    # Pipeline classes
    "HeliosRun",
    "build_run_data",
    "ICFRunData",
    "ICFAnalyzer",
    "ICFPlotter",
    "ICFOutputGenerator",
    # Burn-averaged metrics
    "extract_histories_from_run_data",
    "extract_hot_spot_histories",
    "calculate_burn_averaged_metrics",
    "compare_with_published",
    # Physics modules
    "energetics",
    "neutron_downscatter",
    "pressure_gradients",
    # RHW support
    "RHWParser",
    "RHWConfiguration",
    "load_rhw_configuration",
]
