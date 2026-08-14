# Outer-shell decompression before peak compression — 4 MJ Vulcan

**Status:** observed in the Level 0 synthetic streak, full time resolution,
single run. Not yet checked against a second target or against SPECT3D.
**Run:** `V_80TW_250TW_550TW_450TW_P3-50ps_P2-100ps_tm.exo`, window
15.80–17.00 ns, 421 hydro steps, 256 impact parameters, bands 2–4 / 5–9 /
10–20 keV, no photocathode weighting.

## The observation

Turnaround times, computed over every step rather than a probe grid:

| marker | t_min (ns) | r_min (µm) |
|---|---|---|
| DT gas outer (hot-spot boundary) | 16.8300 | 44.2 |
| DT/CD foam outer | 16.6975 | 133.6 |
| edge 2–4 keV | 16.6234 | 100.8 |
| edge 5–9 keV | 16.5913 | 91.8 |
| edge 10–20 keV | 16.5833 | 79.3 |

The **outer surface of the fuel shell turns around 132 ps before the hot-spot
boundary does.** While the core is still converging to 44 µm, the outer foam
is already expanding: 133.6 µm at 16.70 ns, 224.5 at 16.85, 622.4 at 17.00.

Shell thickness (foam outer minus DT gas outer) tells the same story more
directly:

| t (ns) | thickness (µm) |
|---|---|
| 16.55 | 101.5 |
| 16.70 | 86.0 |
| 16.85 | 178.7 |

The layer reaches its minimum thickness ~130 ps before peak compression and
then roughly doubles during the final deceleration. The shell is coming apart
on the outside while the inside is still doing work.

## Why this is worth having

Every other handle on the compression deficit in this project is an
integrated scalar compared code-to-code — rho_R_cf floors at 0.67 (PDD) and
0.74 (HDD lrm4), hot-spot rho_R far below the 0.3 g/cm2 ignition threshold.
Those say *that* compression is short, not *when* or *how* it is lost.

A turnaround-time split is a time-resolved, single-channel, in-principle
measurable signature of the same physics. It says the loss happens during
final deceleration, on the outer surface, before peak compression — which is
the behaviour expected from a shell that is too thick or riding too high an
adiabat to hold together through deceleration.

## What is NOT established

**This is one run, one code, one dimension.** No mix, no low modes, no
comparison target. The next check is whether the split appears in the
PDD_20_s016 anchor and in HDD lrm4, and whether its magnitude tracks the
adiabat when the prepulse is varied.

**The emission edge is not a clean proxy for the material interface.** The
edges turn around at 16.583–16.623 ns, i.e. **75–115 ps ahead of the foam
interface at 16.6975**, and the lead is channel-dependent (harder bands turn
earlier). So the edge is neither pure material motion nor a pure emissivity
contour: the geometry is bounded by the material, but the 50%-of-peak contour
sweeps outward through the outer foam as it heats, ahead of the material.

An earlier reading of the coarse probe grid concluded "material motion, not
an emissivity contour sweeping." At full time resolution that is too clean.
The correct statement: **the edge tracks the interface in magnitude but leads
it in phase**, and the lead is set by heating rather than by hydrodynamics.
Any velocity inferred by differentiating the emission-edge trajectory
inherits that error.

**Consequence for measurement.** The shell-decompression signature above is
stated from the Lagrangian interfaces, which a real streak camera does not
see. To claim it as an observable we need the mapping from the *measured*
edge turnaround to the *material* turnaround — which depends on the
emissivity threshold, hence on opacity and ionization, hence on exactly the
physics Level 0 approximates most crudely. That mapping is the strongest
current argument for the SPECT3D rung.

## Method notes

- Peak-rho radius (16.8875 ns, 51.6 µm) is `argmax` over zones and jumps when
  the density maximum moves between layers. Not a Lagrangian marker; use
  region interfaces.
- Threshold edge and encircled R50 sit at the uniform-disk ratio of 0.71 from
  ~16.2 ns onward, so the profile is a sharp-edged disk and the threshold is
  well posed there. Before ~16.2 ns the profile is an extended corona and the
  two disagree by up to 4x — the sweep window spans two regimes and
  window-averaged statistics describe neither.
- Bang time is 16.8177 ns (2–4 keV), 12 ps before the hot-spot boundary
  reaches minimum radius. Channel spread is 10 ps, consistent with the 66%
  band-integrated escape fraction: all three bands see essentially the same
  volume at stagnation.
