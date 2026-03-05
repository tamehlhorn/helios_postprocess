# Archived Legacy Modules

These 3 standalone modules are from v2.0 of helios_postprocessor.
Their physics is fully incorporated into `icf_analysis.py` and they
are **no longer imported** by the active package.

| Module | Superseded By |
|--------|---------------|
| `areal_density.py` | `ICFAnalyzer._compute_areal_densities()` |
| `burn.py` | `ICFAnalyzer.analyze_burn_phase()` |
| `hot_spot.py` | `ICFAnalyzer._compute_hot_spot_properties()` |

Note: `burn_averaged_metrics.py` formerly imported these modules.
It was refactored in v3.0 to use `ICFRunData` instead.

Archived March 2026.
