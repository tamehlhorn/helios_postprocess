# Single vs. multiple scattering in the neutron diagnostic — closeout study

**Purpose.** Answer, with a curve rather than an inference, whether there is a
*definite* difference between the neutron diagnostic signature computed with
multiple scattering (OpenMC full transport) and with single scattering only
(NeSST), and identify where the single-scatter approximation stops being
adequate. This is the capstone that closes the neutron-diagnostics
postprocessing task before moving to x-ray imaging.

Driver: T. Mehlhorn / M. (Xcimer) — "do we see a definite difference between
OpenMC (multiple scattering) and NeSST (single scattering)?"

## The clean way to isolate multiple scattering

Comparing NeSST directly to OpenMC conflates *scattering order* with
code-to-code differences (cross-section handling, geometry sampling). Instead we
partition **one** OpenMC run by collision count with `openmc.CollisionFilter`:

- **single-scatter reference** = neutrons that leak with **≤ 1 collision**
  (uncollided primary + once-scattered) — exactly the population NeSST's
  single-scatter model computes;
- **full transport** = all leaking neutrons (multiple scattering, (n,2n), full
  down-scatter cascade).

Same code, same ENDF/B-VIII.0 data, same geometry, one run — so the difference
between the two curves *is* multiple scattering and nothing else.

**Geometry.** An idealized uniform DT sphere (optionally + carbon ablator mass
fraction) with a **central point source**, so every source neutron traverses
exactly ρR = ρ·R and the single-scatter prediction is analytically exact. Radius
fixed (default 50 µm), density scaled to hit each target fuel ρR_D. This is the
textbook control, not a Helios profile — use `--npz` to additionally run the
real N210808 snapshot through the same partition.

## Quantities reported at each areal density

| column | meaning |
|---|---|
| `DSR_single` | DSR (10–12 / 13–15) from ≤1-collision neutrons |
| `DSR_full` | DSR from all leaking neutrons |
| `DSR_frac_diff` | (DSR_full − DSR_single) / DSR_single — the multiple-scatter correction to the *diagnostic* |
| `multiscatter_frac` | fraction of escaping neutrons with ≥2 collisions |
| `rhoR_D_openmc` | ρR_D from the D-elastic **collision count** (counts all collisions) |
| `rhoR_D_analytic` | single-traversal analytic ρR_D (Method-2, exact here) |
| `rhoR_enhancement` | `rhoR_D_openmc / rhoR_D_analytic` — multiple-scatter inflation of the collision-count readout |

## How to run (openmc-env)

```
cd examples
conda activate openmc      # the osx-64 OpenMC env
# pipeline check with NO OpenMC (mocked tallies) — verifies analysis + plots:
python scatter_order_sweep.py --self-test --prefix selftest

# the real sweep:
python scatter_order_sweep.py --rhoR 0.1,0.3,0.6,1.0,1.5,2.0 \
    --particles 200000 --batches 20 --R-um 50 --prefix scatter_sweep

# carbon ablator variant (adds 12-13 MeV carbon pile-up to the down-scatter):
python scatter_order_sweep.py --rhoR 0.1,0.3,0.6,1.0,1.5,2.0 --carbon-frac 0.3 \
    --prefix scatter_sweep_C

# fold in the real N210808 profile as an extra point:
python scatter_order_sweep.py --rhoR 0.3,0.6,1.0 --npz n210808_snapshot.npz
```

Outputs: `<prefix>.csv`, `<prefix>_summary.png` (DSR single vs full, plus
multiply-scattered fraction and ρR readout enhancement vs ρR), and
`<prefix>_spectra.png` (escaping spectra, single vs full, at low/mid/high ρR).

Runtime is a few minutes per point at 200k×20 on the Studio; drop to
`--particles 50000 --batches 10` for a quick look.

## What we expect, and the closeout criterion

The double-scatter probability scales like (ρR·σ)², so the multiple-scatter
correction should be **negligible at low ρR and grow with ρR**. Concretely:

- At this design's areal density (traversed fuel ρR_D ≈ 0.3, total ≈ 1.4 g/cm²,
  optical depth τ ≈ 0.3) we expect the **DSR** correction to be small — a few
  percent — because doubly-scattered neutrons mostly land *below* the 10–12 MeV
  window (in the sub-10-MeV tail and the backscatter edge), not in the DSR
  numerator. This is why NeSST and OpenMC agree on the DSR for N210808.
- The **collision-count ρR_D readout** is more sensitive: our N210808 comparison
  already showed OpenMC 0.372 vs the single-traversal analytic 0.297 (~25%). The
  sweep will quantify how much of that 25% is multiple scattering (re-scattered
  neutrons registering extra D-elastic events) versus the energy-dependent cross
  section.
- The single-scatter approximation should visibly break — DSR_frac_diff climbing
  into the 10%+ range — as fuel ρR_D approaches and exceeds ~1 g/cm² (strongly
  compressed, igniting designs). That crossover is the deliverable number.

**Task closes** when the sweep confirms (a) DSR_single ≈ DSR_full within Monte
Carlo error at low ρR, (b) a monotonic, (ρR)²-like divergence at high ρR, and
(c) the escaping-spectrum sub-10-MeV tail is systematically enhanced under full
transport. Then: NeSST is validated as the single-scatter model, its domain of
validity (low-to-moderate ρR) is mapped, and OpenMC is the tool of record where
ρR is high — a clean, defensible closeout.

## Results (fill after running)

| ρR_D (g/cm²) | DSR_single % | DSR_full % | ΔDSR % | multi-scatter % | ρR_D enh. |
|---|---|---|---|---|---|
| 0.1 |  |  |  |  |  |
| 0.3 |  |  |  |  |  |
| 0.6 |  |  |  |  |  |
| 1.0 |  |  |  |  |  |
| 1.5 |  |  |  |  |  |
| 2.0 |  |  |  |  |  |
| N210808 (real) |  |  |  |  |  |

_Paste `scatter_sweep.csv` here and drop in `scatter_sweep_summary.png` /
`_spectra.png`._

## Caveats

Idealized uniform sphere with a central source — chosen to make single scatter
exact and isolate the physics, not to reproduce a Helios profile (the `--npz`
run does that). `CollisionFilter` counts all collision types; at these energies
these are dominantly elastic scatters, so collision count ≈ scatter order.
(n,2n) secondaries are included in "full" and contribute to the low-energy
floor. Photon transport is off for the sweep (neutron signature only); turn it
back on in the main pipeline for the gamma diagnostic.

Reproduce: `examples/scatter_order_sweep.py` (openmc-env).
