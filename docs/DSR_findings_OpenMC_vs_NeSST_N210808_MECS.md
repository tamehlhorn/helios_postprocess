# DSR, areal density, and yield unfold: OpenMC transport vs NeSST on N210808

Findings note (Kyle Keipper, Will). Written 2026-07-30, updated 2026-07-31
(added the self-consistent first-principles coefficient and the final pipeline
DSR after the traversed-rhoR switch). Run: 090425_KK_opt_2burn (N210808 clean-1D
burn), bang-time snapshot, ENDF/B-VIII.0, 400k particles, neutron + photon
transport.

## Bottom line
Three methods on the same Helios snapshot agree once the areal-density
definition is pinned down. The OpenMC/NeSST "discrepancy" was the Method-1 vs
Method-2 distinction from Kyle's slide, now resolved and wired into the
pipeline: run_analysis STEP 5b now reports DSR 4.52% (fuel) / 5.02% (+carbon)
on the emission-weighted (traversed) rhoR, matching OpenMC (4.3%). A
first-principles reproduction of the NIF yield unfold shows the standard 4*DSR
exponent is design-specific (~5 for this HDC-ablator target).

## 1. Two areal densities
| quantity | OpenMC (transport) | NeSST (as shipped) | analytic Method-2 |
|---|---|---|---|
| rhoR_D (g/cm^2) | 0.372 | 0.567 | 0.297 |
| DSR (10-12/13-15) | 4.31% | 10.07% (fuel 9.51%) | - |

NeSST rhoR_emission = center-to-edge int_0^R rho dr (Kyle Method 1). OpenMC and
the analytic Method-2 (isotropic 4pi angle-averaged chord integral, validated:
central source->rho*R, uniform->0.75x) give the emission-weighted TRAVERSED
rhoR. Both traversed values (0.30-0.37) ~ half the center-to-edge (0.57); DSR
follows (4.3% vs 10%). Abu-Shawareb: "the areal density of the system is
emission weighted to larger radii, which is why...the areal density decreases
as the yield increases" -> DSR tracks traversed rhoR (Method 2), not
center-to-edge. Residual OpenMC(0.372) vs Method-2(0.297) = multiple-scatter +
energy-dependent sigma (transport-only).

Pipeline run (STEP 5b, traversed mode) on the same snapshot: traversed rhoR
D=0.298, T=0.447, C=0.664 g/cm^2 (fuel D+T = 0.745, total = 1.409); DSR 4.52%
fuel / 5.02% total. The pipeline traversed rhoR_D (0.298) reproduces the
analytic Method-2 (0.297) to 1%.

ACTION: NeSST unchanged (pure rhoR->DSR mapper). Changed our code:
_composition_scatter now feeds it the traversed rhoR (new
target_composition.traversed_species_areal_densities), center-to-edge retained,
selectable via rhoR_mode. Lowers run_analysis STEP 5b DSR to emission-weighted.

## 1a. First-principles DSR<->rhoR coefficient (self-consistent)
The empirical NIF calibration is rhoR = 20.4 * DSR (g/cm^2 per unit DSR). Our
first-principles value must use the SAME areal-density convention on both sides.
Using the traversed (Method-2) rhoR that the DSR is now derived from:

| coefficient | value (g/cm^2 per DSR) | definition |
|---|---|---|
| total (D+T+C) | 28.1 | total traversed rhoR / total DSR (1.409 / 0.0502) |
| fuel (D+T) | 16.5 | fuel traversed rhoR / fuel DSR (0.745 / 0.0452) |

These bracket the empirical 20.4, as expected for a clean over-igniting 1-D burn
where rhoR/DSR run high. NOTE (reporting fix): the pipeline previously printed a
single "31.4" that mixed conventions -- center-to-edge fuel rhoR (1.42) over the
traversed fuel DSR (0.0452). That number was neither Method 1 nor Method 2 and
has been corrected to the two self-consistent values above (commit: neutron_scatter
coefficient uses the traversed convention on both sides). The empirical relation
is tied to the total, emission-weighted rhoR, so 28.1 is the value to compare
against 20.4; the fuel-only 16.5 is the internal D+T check.

## 2. First-principles yield unfold  Y_total = Y(13-15) e^(4*DSR)
| | value |
|---|---|
| escaping 13-15 MeV / born | 0.786 |
| Y_total/Y(13-15) OpenMC (1/N) | 1.272 |
| paper e^(4*DSR) | 1.188 |
| effective exponent k | 5.6 (~5.1-5.2 after DD/TT correction) |

k >> 4 because n-12C elastic loses little energy per collision, pushing 14 MeV
neutrons into 12-13 MeV (out of the unscattered window, above the 10-12 DSR
window), inflating scattered-out/DSR. Consequence: standard 4*DSR unfold
under-estimates total yield ~7% (1.19 vs 1.27) for a heavy-ablator target. A
DT-only source run would pin the exponent exactly.

## 3. Gamma signature
0.0070 photons / source neutron escaping, ~540:1 dominated by 12C(n,n'g) 4.4 MeV
(6.8e-3) over H(n,g) 2.2 MeV (1.3e-5). 4.4 MeV line = clean HDC-ablator carbon
signature.

## Caveats
Clean-1D burn (Y~15.6 MJ vs 1.37 measured); rhoR/DSR run high vs the real shot
(published DSR 2.87%). Demonstrates method + corrections, not an experimental
match. One snapshot. Surface-tracking warnings negligible (DSR stable to ~3%
across two geometry rebuilds). sigma_el_D = 0.641 b from ENDF/B-VIII.0 at
14.06 MeV.

## Recommendations
1. Use emission-weighted (Method-2/OpenMC) rhoR vs measured DSR (pipeline now
   default).
2. Don't assume 4*DSR unfold for carbon ablators; derive per design (~5 here).
3. Pin which rhoR the experimental DSR<->rhoR calibration is tied to (the
   emission-weighted total, per the 28.1 vs 20.4 comparison above).
4. Next: degraded/under-igniting N210808 (lower rhoR) to map how the
   Method-1/2 gap, the unfold exponent, and the DSR<->rhoR coefficient scale
   with yield.

Reproduce: examples/dump_snapshot.py (Helios env) -> examples/compare_openmc_nesst.py (openmc-env).
