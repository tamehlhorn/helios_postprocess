# Olson 2021 Fig 7 / Fig 8 — digitized reference data

Reference radial line-outs and time-history ρR curves digitized from
Olson, R.E. et al., 2021 (or related reference comparing LILAC, xRAGE,
HYDRA for the Polar Direct Drive (PDD) target design that
`Olson_PDD_20_*` family is calibrated against).

## What's here

### Figure 7 — radial line-outs at ignition (ρR_hs = 0.3 g/cm²)

3 panels (one per code) showing:
- T_ion(r) in keV (red curve, shares y-axis)
- 0.1×ρ(r) in g/cm³ (black curve, scaled by 0.1 to fit shared axis)
- Annotation: "hot spot radius = 120 µm" with arrow at ρR=0.3

Files:
- `fig7_LILAC_rho_r.csv`  — ρ(r) digitized from LILAC panel
- `fig7_xRAGE_rho_r.csv`  — ρ(r) digitized from xRAGE panel
- `fig7_HYDRA_rho_r.csv`  — ρ(r) digitized from HYDRA panel
- `fig7_LILAC_Tion_r.csv` — T_ion(r) digitized from LILAC panel (optional)
- (similar for xRAGE / HYDRA T_ion if needed)

Format:
    r_um,rho_gcc
    0.0,20
    20.0,18
    ...

### Figure 8 — ρR(t) time histories

3 panels (one per code) showing:
- Total ρR(t) (black curve)
- Hot-spot ρR(t) (red curve, T_ion > 4.5 keV)

Files:
- `fig8_LILAC_rhoR_t.csv`
- `fig8_xRAGE_rhoR_t.csv`
- `fig8_HYDRA_rhoR_t.csv`

Format:
    t_ns,rhoR_total_gcm2,rhoR_HS_gcm2

## Digitization status

Initial values were estimated by visual inspection of the screen-rendered
figures (May 2026). Accuracy: roughly ±15–20% on individual ρ readings,
±25% on integrated mass via `examples/estimate_imploded_mass_from_rhor.py`.

**For production use, re-digitize with WebPlotDigitizer
(<https://apps.automeris.io/wpd/>) and replace these files in place.**
The pipeline reads CSVs by path; format is stable.

## Why these are here

Olson 2021 reports ρR(t) and T_ion(r) directly but does NOT publish
imploded-DT-mass for any of the three codes. The Helios calibration work
needs a reference imploded mass for the "structural mass deficit"
comparison test in `burn_averaged_metrics.compare_with_published`.

`examples/estimate_imploded_mass_from_rhor.py` integrates these ρ(r)
profiles via M = 4π ∫ρ(r) r² dr to derive reference imploded-mass
estimates for each code, which then go into the `*_published.json` files
as `imploded_DT_mass_mg` entries.
