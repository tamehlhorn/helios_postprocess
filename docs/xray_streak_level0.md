# Synthetic x-ray streak camera — Level 0 design note

Opens the x-ray diagnostics workstream. Structure deliberately mirrors the
neutron path: a transparent analytic forward model first, validated against
the heavyweight code, so we know *where* the heavyweight code is doing work.

## Validation ladder

| Rung | Engine | Physics | Purpose |
|---|---|---|---|
| **Level 0** | this package | thermal free-free continuum, LTE continuum absorption | fast, inspectable, no black box |
| **Level 1** | SPECT3D | continuum + LTE opacities, self-absorption | isolates the opacity contribution |
| **Level 2** | SPECT3D | full collisional-radiative non-LTE, lines | isolates NLTE + line contribution |

The neutron-side analogue is the single-vs-multiple-scatter DSR sweep: the
deliverable is not "SPECT3D gives the right answer" but a **quantitative map
of where Level 0 fails**, so we know when the cheap model is admissible.

## What the diagnostic measures

An XRSC on self-emission records `I(x, t)`: one spatial axis through a slit
against sweep time, band-integrated over the filter/photocathode response.
From a 1D run the extractable observables are:

- **emission-edge trajectory `R(t)`** — the ablation front / coronal emission
  boundary during the drive. **Not the cold shell.** Comparing this edge to a
  Lagrangian shell trajectory is the standard error with this diagnostic; the
  extraction routine compares like to like by pulling the same 50% edge from
  the simulated emission.
- **x-ray bang time and emission FWHM** — the x-ray bang time is not the
  neutron bang time (different Te and ρ weighting), and the offset between
  them is itself a calibration observable.
- **stagnation flash size** — hot-spot emission radius, pairs directly with
  the HS ρR deficit that is the current binding residual.
- **channel ratio** — hard/soft band ratio as a coarse Te(t) proxy when a
  multi-channel filtered configuration is used.

A 1D spherical run is perfectly symmetric, so there is no low-mode or mix
information in the synthetic image. That is a scope limit, not a defect.

## Physics in Level 0

Free-free emissivity, 4π-integrated form (Rybicki & Lightman 5.14b) reduced
to per-steradian, per-eV:

```
j_E = 6.8e-38 · Zbar² · n_e · n_i · T_K^(-1/2) · exp(-E/kTe) · ḡ_ff / (4π h_eV)
```

with the thermally-averaged Born/Elwert Gaunt factor
`ḡ(u) = (√3/π)·e^(u/2)·K₀(u/2)`, `u = E/kTe`. Absorption follows from
Kirchhoff, `κ_E = j_E / B_E`, so the source function is exactly the Planck
function and the "with absorption" mode is an LTE continuum calculation.

**Absent by design** — these are what the SPECT3D rungs supply, and what the
Level 0 ↔ Level 1 comparison quantifies:

- free-bound (radiative recombination) continuum and its edges
- bound-bound line emission and absorption
- non-LTE / collisional-radiative ionization balance
- real opacity tables rather than pure continuum

## Geometry

1D spherical symmetry makes the chord integral analytic — no ray tracer. The
path length of a chord of impact parameter `p` through the shell bounded by
`r_k`, `r_{k+1}` is `√(r_{k+1}²−p²) − √(max(r_k²−p²,0))`, appearing twice.

Two solvers: an optically-thin sum, and an ordered short-characteristic
formal solution (far side inward to the tangent zone, then back out) that
also accumulates chord optical depth. **The τ map is the hand-off argument to
SPECT3D**: it identifies the (t, p, E) region where the thin assumption fails
before any SPECT3D run is made.

## Module layout

```
helios_postprocess/xray/
    emissivity.py    Helios state -> j_E(r), κ_E(r), source function
    chords.py        spherical chord integration, thin + formal solvers
    response.py      filter transmission, photocathode QE, channel response
    streak.py        radiance cube, sweep/IRF/MTF, observable extraction
    spect3d_io.py    profile export + SPECT3D reader onto the same contract
examples/
    xray_run_recon.py            pre-flight: fields, regions, window, tau
    xray_emission_provenance.py  which zones make the light; band selection
    xray_streak_synthetic.py     runner, 3-page PDF + .npz cube
    patch_xray_variables.py    wires mean_charge / ion_density into data_builder
tests/
    test_xray_level0.py        14 tests, no .exo required
```

The master object is `RadianceCube`: `I(t, p, E)` in
erg s⁻¹ cm⁻² sr⁻¹ eV⁻¹. Both the imaging streak (integrate over E with the
channel response) and the spectral streak (integrate over p) are projections
of the same cube, so the two instruments are guaranteed self-consistent.
SPECT3D output is read onto the identical contract, so the cross-validation
is a diff rather than a re-implementation.

## Known weak points

1. **Filter/photocathode response.** The built-in cold photoabsorption model
   is an anchored power law (`μ/ρ ∝ A·(E/1.5 keV)^−2.9` with K-edge jumps),
   good to roughly ±25% away from edges. It **cancels in the Level 0 ↔
   SPECT3D ratio** because both rungs use the same response, but it is not
   adequate for absolute comparison to measured data. Replace with CXRO/Henke
   tables (`load_henke_file`) before any photometric claim.
2. **Zbar.** Free-free scales as Zbar², so `mean_charge` matters more here
   than anywhere else in the pipeline. It is not currently in
   `_VARIABLE_MAP`; run `examples/patch_xray_variables.py` first. Without it
   the module falls back to `n_e/n_i` and warns.
3. **Te vs Ti.** Coronal emission is set by Te. Falling back to Ti (which the
   module will do, loudly) makes the drive-phase signal wrong.
4. **Non-uniform Helios time base.** EXODUS clusters timesteps around
   stagnation. The cube is built on the native steps and resampled onto a
   uniform sweep base before the IRF is applied — a real streak camera samples
   uniformly, and applying an IRF on a non-uniform base would be wrong.
5. **Halfraum targets.** For HDD the Cu converter and He fill dominate the
   emission. Use `--capsule-only` / `zone_slice`, or the streak is a picture
   of the halfraum, not the capsule. Detection keys off region names; if they
   fail to decode, set the boundary by hand with `--zone-stop N`.
6. **Outer-boundary zones.** The outermost Lagrangian zones carry little
   mass but can reach very high Te, and their contribution grows as
   `exp(E/kT_bulk)` with photon energy, so they can take over the hard tail.
   Whether they do is controlled by coronal *density*, not temperature — the
   `n²` penalty is severe and a very hot outer zone is not by itself evidence
   of contamination. Run `xray_emission_provenance.py` to decide empirically;
   mitigate with `--rho-floor`, and always report with and without it. Full
   write-up in `docs/claude_xray_outer_zone_emission.md`.
7. **Channel placement.** `default_channels()` is cut for a softer,
   OMEGA-class direct-drive spectrum. On a harder target most of the signal
   lands in one channel and the hard/soft ratio stops working as a Te proxy.
   Use `--bands`/`channels_from_edges` guided by the provenance percentiles.
8. **Hydro dump cadence.** The synthetic camera cannot resolve anything
   faster than the EXODUS output interval. If the dump interval exceeds the
   requested IRF the model warns, and any bang-time or FWHM number inherits
   the hydro cadence as its real resolution limit.

## Open items for the SPECT3D side

- Does SPECT3D on the Studio import the Helios output natively, or does it
  want the exported profile (`spect3d_io.export_profiles` writes
  `r_in, r_out, ρ, Te, Ti, Zbar, n_e` per snapshot)?
- Can it be driven headless/batch? If GUI-only, the time sweep has to be
  exported snapshot-by-snapshot and the loop restructured.
- Confirm the export column layout once against a real file before trusting
  `load_spect3d_cube` on a batch — the reader takes explicit column indices
  because the headers vary by version.
