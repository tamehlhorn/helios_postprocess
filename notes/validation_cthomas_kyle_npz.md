# Validation — extraction vs K. Keipper `neutronics_output.py` (2026-07-20)

**Run:** `cthomas_baseline_auto_zone_alpha_fT30` (Vulcan HDD, C. Thomas baseline).
**Reference:** Kyle Keipper's `neutronics_data.npz` (from `neutronics_output.py`,
RHINO backend), cross-checked against the RHINO analysis report for the same run.

## Self-consistency vs the RHINO report

| quantity | npz | RHINO report | |
|---|---|---|---|
| bang time (ns) | 23.601 | 23.6 | ✓ |
| DT neutron yield | 1.294e18 | 1.29e18 | ✓ |
| DD neutron yield | 9.361e15 | — | |

Grid: 1423 timesteps × 189 zones; all six fusion channels present (DHe3 and
pB11 identically zero). Peak ρ ≈ 100 g/cc, peak T_ion ≈ 16 keV.

## Routine verification — identical-input test

Feeding Kyle's own arrays (ρ, T_i, per-channel rates, geometry, dt) into our
`neutron_spectrum.neutron_weighted_profiles` reproduces his stored DT-weighted
profiles to machine precision:

| profile | max \|Δ\| | mean \|Δ\| |
|---|---|---|
| `<rho(r)>` | 0.00% | 0.00% |
| `<T_ion(r)>` | 0.00% | 0.00% |
| `<DT rate(r)>` | 0.00% | 0.00% |
| `<DD rate(r)>` | 0.00% | 0.00% |
| `<TT rate(r)>` | 0.00% | 0.00% |

**⇒ our port of Kyle's neutron-weighted-averaging algorithm is exact.** This is
the verification gate for importing the routine into RHINO.

## In-house observables (added; Kyle's extraction stops before these)

From the same arrays, our routines produce the nTOF spectral observables:

- nTOF **T_DT = 8.39 keV** (primary peak 14.060 MeV, FWHM 513 keV — Brysk-consistent)
- nTOF **T_DD = 7.61 keV** (peak 2.450 MeV)
- burn-averaged ⟨T_i⟩ = 8.52 keV, burn FWHM ≈ 363 ps
- **T_DT > T_DD** ordering reproduced (DD samples the cooler, denser fuel),
  the same sense as the N210808 data (10.9 vs 8.94 keV).

## Notes

- Bang time from our EXODUS-snapshot argmax is 35 ps later than RHINO's
  fine-grid bang (23.636 vs 23.601) — a definition difference within one output
  cadence, not an error. Reconcile later if identical values are wanted.
- Absolute yield must come from the cumulative product
  `TimeIntFusionProd_n_1406` (1.294e18, which `extract_neutronics` reads); the
  EXODUS rate integral under-counts ~18× (CLAUDE.md convention 8). Confirmed.

## Status

- **VERIFIED:** neutron-weighted averaging routines (identical-input, 0.00%).
- **PENDING (last mile):** end-to-end from the raw `.exo` through our reader
  (`HeliosRun` → `build_run_data` → `extract_neutronics`) vs Kyle's npz. Run
  `examples/compare_to_kyle_npz.py <run_dir> <kyle_npz>` on the run to close it —
  this checks our EXODUS variable mapping / COM formula / geometry against his.
