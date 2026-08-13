# Outer-zone emission contamination in synthetic x-ray streak records

**Status:** open question on the 4 MJ Vulcan run, criterion and mitigation
implemented, not yet evaluated against data.
**Raised by:** recon on `V_80TW_250TW_550TW_450TW_P3-50ps_P2-100ps_tm.exo`,
which reported Te spanning 304 eV to **80 keV** at peak emission, with the
emitted-power percentiles at 3.1 / 6.9 / 15.4 keV.

## The concern

In a Lagrangian mesh the outermost zones of an expanding corona carry
negligible mass but can reach very high electron temperature. Thermal
free-free emissivity is

```
j_E  ∝  Zbar² · n_e · n_i · T^(-1/2) · exp(-E / kTe)
```

so a hot tenuous zone is penalised hard by the `n²` prefactor but rewarded by
the exponential, and the reward grows with photon energy. Comparing a corona
(subscript c) to the bulk (subscript b), the ratio of contributions at photon
energy `E` is

```
R(E)  ≈  (n_c / n_b)² · (V_c / V_b) · √(T_b / T_c) · exp( E/kT_b − E/kT_c )
```

The exponential factor is the whole story: with `kT_b = 3 keV` and
`kT_c = 80 keV`, `R` grows by `exp(E/kT_b)` — a factor of ~123 between the
soft band and 15 keV. A population contributing 0.1% of the soft signal can
therefore contribute over 10% of the hard signal.

## Why a streak camera is the wrong instrument to be careless with

Two consequences, and the second is the one that matters for our purposes.

**Hard-channel photometry.** A filtered hard channel exists precisely to
sample the exponential tail, so it is maximally exposed to this population.

**Emission-edge trajectory.** The corona sits at *large impact parameter* —
outside the capsule limb. `emission_edge_trajectory` takes the outermost 50%
crossing of the normalized radial profile, so a bright coronal halo pushes
the measured edge outward and flattens the apparent convergence. That edge
`R(t)` is the main observable we wanted from this diagnostic. Contaminating
it would not produce an obviously wrong-looking streak image; it would
produce a plausible one with a biased trajectory.

## What is NOT established

**A very hot outer zone is not by itself evidence of contamination.** The
`n²` penalty is severe. Worked through the code with a bulk at 1 g/cm³ and
3 keV against a 80 keV corona filling 215× the volume, at 30 keV:

| corona ρ (g/cm³) | corona share at 500 eV | corona share at 30 keV |
|---|---|---|
| 1e-6 | 1.0e-6 | 2.5% |
| 1e-4 | 1.0e-6 | 2.5% |
| 1e-3 | 1.0e-4 | 72% |
| 1e-2 | 1.0e-2 | 99.6% |

The takeover is sharp and it is controlled by coronal **density**, not by
coronal temperature. At 1e-4 g/cm³ there is no artifact worth discussing; at
1e-3 the hard channel is a picture of the corona. So the 80 keV reading on
the Vulcan run raises the question but does not answer it — the density of
those zones decides, and that has to be measured, not assumed.

An independent reason for caution at 80 keV: thermal free-free from a fluid
code assumes a Maxwellian at the zone temperature and an LTE-ish ionization
balance. Neither is safe in a tenuous 80 keV corona, so even if those zones
*do* dominate the hard tail, Level 0 is not the right model for them — and
neither is SPECT3D's LTE rung.

## Detection

`examples/xray_emission_provenance.py` decides empirically:

- decomposes band-integrated emitted power zone by zone, aggregated by region
- prints power share **and mass share** of the outermost zones
- recomputes the spectrum with the outer zones trimmed and reports the shift
  in the 90th percentile

The verdict is deliberately **mass-aware**. The artifact signature is
*bright and massless* — a large share of the power carried by a negligible
share of the mass. Outer zones emitting roughly in proportion to their mass
are simply the ablator, and trimming them would discard real signal. An
earlier version of the check flagged on power share alone and produced a
false positive on exactly that case.

## Mitigation

`rho_floor` (as `StreakConfig.rho_floor_g_cc`, or `--rho-floor` on the
runner) treats zones below a density threshold as neither emitting nor
absorbing. Masking emission without masking absorption would leave the zone
shaping the emergent spectrum after being declared non-emitting, so both are
masked together.

**The floor is a mask, not physics.** Standing rule: run with and without,
and treat any hard-channel result that depends on the floor as unconverged
rather than as a measurement. The floor value belongs in the reported
metadata; it is carried in `RadianceCube.meta` and printed on the report page
for that reason.

## Connection to the boundary-condition work

How far the corona expands — and therefore its density and emitting volume —
is set by the outer boundary treatment, the same clamped-vs-free choice that
turned out to dominate the absorbed-energy ledger on HDD lrm4 (clamped 50%,
free 82%, published 97%). If the outer zones do prove to carry the hard tail,
this diagnostic inherits that sensitivity, and any hard-channel result would
have to be quoted against a stated boundary condition.

## Resolution on the 4 MJ Vulcan run

Provenance at peak emission (16.80 ns) closed the question: the outermost
five zones carry **0.000% of the band-integrated power** against 0.223% of
the mass, and trimming them shifts the 90th percentile by **0.00%**. Those
zones sit at r = 15-21 mm with rho = 1e-6 to 1e-7 g/cm3 and Te = 1784 eV —
tenuous, but nowhere near hot enough for the exponential to overcome the n^2
penalty at that density. No artifact. No floor needed on this run.

The 80 keV reading that raised the alarm is the DT gas region (power-weighted
<Te> = 39.8 keV), which is the hot spot, not the boundary — dense (22.9
g/cm3) and physical. The alarm was raised on the right criterion and answered
the right way.

## Next step

Run the provenance script on the Vulcan target and read the mass share, not
just the power share. Three outcomes:

1. **Outer zones bright and massless** → set a floor, report both, and treat
   the hard channel as provisional.
2. **Outer zones bright in proportion to mass** → that is the ablator; do not
   trim, and the 15.4 keV percentile is real.
3. **Outer zones negligible** → the question closes and the hard channel is
   clean.
