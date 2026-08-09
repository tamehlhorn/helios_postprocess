# Analytic ρR_D vs. OpenMC transport — how they fit together

**To:** Kyle Keipper
**From:** Tom Mehlhorn
**Re:** Neutron sampling from Helios → OpenMC, and where the Method-2 areal-density formula fits
**Date:** 2026-07-29 (updated 2026-07-31: σ_el,D corrected to the ENDF/B-VIII.0 value)

## Bottom line

Kyle's Method-2 areal density is the right deterministic ρR_D and we keep it. The
adjustment is a change of ROLE, not physics: Method 2 computes an analytic
DIAGNOSTIC (a scalar ρR_D), whereas OpenMC needs from Helios a SOURCE (birth
positions + energy spectrum) plus GEOMETRY (ρ(r) and composition). OpenMC then
performs the same line integral by Monte Carlo and carries it beyond single
scatter. Cleanest arrangement: use Method 2's ingredients to build the OpenMC
source, let OpenMC transport, keep Method-2 ρR_D as the independent benchmark.

## What each side computes

Method 2 (analytic): for a neutron born at r_s along chord direction cosine μ,
r(l,μ,r_s)=√(r_s²+l²+2 l r_s μ), and
  ρR_D(r_s,t) = (1/4π)∫dΩ ∫₀^{l_exit(μ)} ρ_D(r(l,μ,r_s)) dl
(= ½∫₋₁¹dμ ∫ρ_D dl for spherical symmetry), S-weighted over radius and
burn-averaged in time. Exact for the mean deuterium column an isotropically
emitted neutron traverses: no scattering, one species, frozen profile.

OpenMC: samples birth neutrons from S(r,t) + fusion spectrum, transports through
ρ(r). Computes the same ∫ρ dl plus full down-scatter spectrum, (n,2n), secondary
photons, energy-angle correlations, multiple scattering. Method 2's ingredients
(sample r_s ∝ S, isotropic emission) ARE the OpenMC source recipe; only the
∫ρ dl integral is made redundant by transport.

## Formal correspondence + built-in validation

Single-scatter, single-species, static-profile limit: same integral two ways.
Expected D elastic scatters per source neutron C_D = (N_A/M_D)·σ_el,D·ρR_D, so
  ρR_D(OpenMC) = C_D · M_D / (N_A · σ_el,D),  σ_el,D(14.06 MeV) ≈ 0.641 b
(ENDF/B-VIII.0; the same value used in the findings note, not the ~0.8 b round
number in earlier drafts). Run OpenMC on the same snapshot, tally D-elastic,
convert; the two ρR_D should agree. Divergence (OpenMC a few % high at large ρR)
is the multiply-scattered tail — the same residual slope in the first-principles
ρR↔DSR coefficient, and the reason to run transport at high ρR. Agreement
validates the source/geometry wiring; the gap measures the transport correction.

On that coefficient: with the DSR now taken on the traversed (Method-2) rhoR,
the self-consistent first-principles ρR↔DSR value on N210808 is 28.1 g/cm² per
DSR (total, D+T+C) / 16.5 (fuel D+T only), bracketing the empirical NIF 20.4 —
consistent with this being a clean over-igniting 1-D burn. See the findings note
for the numbers and the reporting fix.

## Hazard: do NOT pre-attenuate

Never use ρR_D to attenuate/scale the source before OpenMC — transport attenuates
during the run, so a pre-attenuated source depletes the population twice. Same
failure mode as the single-scatter double-count we already fixed (undepleted
primary + scattered → unphysical >1.3× escaping). Hand OpenMC the UNDEPLETED
birth source (full yield, all channels); ρR_D stays a diagnostic check.

## Assumptions to keep visible

1. Volume Jacobian: if S is per-volume emissivity, radial birth pdf p(r) ∝
   S(r)·4πr². Method-2 weighting reads as per-shell (∫S dr_s) — keep conventions
   consistent or the source biases toward center.
2. Frozen profile: ρ(r) at emission time. Excellent (neutron transit ~1-2 ps vs
   ~100 ps hydro) but an assumption.
3. Burn-averaging is nonlinear: ρR of time-averaged profile ≠ time-average of
   per-snapshot ρR. Run a few emission-weighted snapshots, S(t)-weight tallies.
4. Species: Method 2 as drawn is D-only (right for fuel ρR_D). Down-scatter and
   the OpenMC geometry need all scatterers along the chord — D, T, and carbon.
5. Isotropy: fine for thermal DT.

## Recommended workflow

One set of Helios inputs, two paths, cross-checked:
1. Per snapshot extract r(shells), ρ(r), (w_D,w_T,w_C), emissivity S(r).
2. Analytic (Kyle): Method-2 ρR_D, S-weight over radius, burn-average → ρR_D.
3. Transport (OpenMC): nested-sphere geometry from ρ(r)+composition; source from
   S(r) (radial, isotropic) + fusion spectrum; feed UNDEPLETED source; tally
   D-elastic, leakage spectrum, DSR.
4. Compare: convert OpenMC D-elastic → ρR_D vs analytic. Agreement validates
   both; residual = multiple-scatter correction.
5. Hand escaping source to the chamber OpenMC model (Kyle/Owen) for TBR/damage.

Implemented as a working sketch in `helios_to_openmc.py` (steps 3-4): builds
geometry + S(r)-weighted isotropic source, runs neutron transport, returns ρR_D
from the D-elastic tally for a one-line comparison against Method 2.

## What Helios hands OpenMC

- Source spatial: birth ∝ S(r,t) (with 4πr² Jacobian), isotropic direction.
- Source energy: fusion birth spectrum — DT 14.06 (Muir/Ballabio) + DD 2.45 + TT
  continuum, exact per-channel yields, undepleted.
- Geometry: concentric shells with ρ(r) and (w_D,w_T,w_C) at chosen time(s).
- NOT an input: a pre-computed ρR. OpenMC produces it; Method 2 checks it.
