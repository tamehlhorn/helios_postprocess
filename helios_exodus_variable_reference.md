# Helios EXODUS Variable Reference

A reference for all EXODUS (`.exo`, netCDF-based) variables consumed by the
`helios_postprocess` pipeline. Information assembled from
`helios_postprocess/data_builder.py` (`_VARIABLE_MAP` and `_try_load_nc`
calls), verified against run logs from `Olson_PDD_26b_burn` (initial) and
`VI_6` / `HDD26_DTI40_1ns130_FL04_lrm4_nb` (May 2026 updates).

**Source of truth:** for the current authoritative mapping, see
`helios_postprocess/data_builder.py` (`_VARIABLE_MAP` and direct netCDF
loads). This document can become stale if the map changes; re-generate
from source if in doubt.

**Conventions**:
- Shapes are `(n_times, n_zones)` for zone-centered arrays,
  `(n_times, n_zones+1)` for node-centered arrays,
  `(n_times,)` for scalar-per-timestep, `(n_zones,)` for static-per-zone.
- `n_times` is the number of EXODUS output timesteps
  (typically 162–296 for current runs; varies by output cadence).
- `n_zones` is the number of Lagrangian zones
  (350 for Olson PDD, 410 for Xcimer HDD designs).
- Helios stores most physical quantities in CGS + ergs. The pipeline
  converts on load or during analysis (see Unit Conversions section).

---

## 1. Core zone-centered arrays (required)

These must be present or the pipeline emits warnings and downstream
analysis may fail.

| Pipeline attribute | EXODUS candidates (first match wins) | Shape | Unit (EXODUS) | Description |
|---|---|---|---|---|
| `mass_density` | `mass_density`, `dens`, `density` | (n_t, n_z) | g/cm³ | Zone-averaged mass density |
| `ion_temperature` | `ion_temperature`, `temp`, `temperature` | (n_t, n_z) | eV | Zone-averaged ion temperature |
| `ion_pressure` | `ion_pressure`, `pres`, `pressure` | (n_t, n_z) | J/cm³ | Zone-averaged ion pressure. Multiply by `1e-8` → Gbar |
| `zone_mass` | `zone_mass`, `mass`, `cell_mass` | (n_t, n_z) | g | Lagrangian mass per zone (constant in time for each zone) |

---

## 2. Node-centered arrays (required)

Evaluated at zone boundaries — hence one more point than there are zones.

| Pipeline attribute | EXODUS candidates | Shape | Unit | Description |
|---|---|---|---|---|
| `zone_boundaries` | `zone_boundaries`, `coord`, `coordx` | (n_t, n_z+1) | cm | Radial positions of Lagrangian zone boundaries |
| `velocity` | `fluid_velocity`, `velocity`, `vel` | (n_t, n_z+1) | cm/s | Radial fluid velocity at zone boundaries. Negative = inward (imploding). Multiply by `1e-5` → km/s |

---

## 3. Temperature and pressure — secondary components (optional)

These are used to reconstruct the 3-component total pressure
(`ion + elec + rad`). Pipeline behavior when missing is documented in the
notes below.

| Pipeline attribute | EXODUS candidates | Shape | Unit | Description |
|---|---|---|---|---|
| `elec_temperature` | `elec_temperature`, `electron_temperature` | (n_t, n_z) | eV | Zone-averaged electron temperature |
| `elec_pressure` | `elec_pressure`, `electron_pressure` | (n_t, n_z) | J/cm³ | Zone-averaged electron pressure |
| `rad_pressure_true` | `rad_pressure` (EXODUS native) | (n_t, n_z) | J/cm³ | **Pure radiation pressure only** (no electron contribution). Loaded directly by `data_builder` and stored on `data.rad_pressure_true`. Use this for `U_rad = 3·P_rad·V` calculations in the energy-balance diagnostic. |

**Pipeline-computed `data.rad_pressure` (non-ion pressure sum):**

| EXODUS available | `data.rad_pressure` equals | Convention name |
|---|---|---|
| `elec_pressure` AND `rad_pressure` | `elec_pressure + rad_pressure` | 3-component |
| `elec_pressure` only | `elec_pressure` | 2-component (PDD_8 path) |
| `rad_pressure` only | `rad_pressure` | rad-only |
| neither | zeros | fallback |

Downstream, `total_pressure = ion_pressure + data.rad_pressure`. Used
for hot-spot pressure, ablation pressure, and shock breakout detection.

**Important distinction for energy accounting:** `data.rad_pressure` is
typically the 3-component mixed quantity; for the energy-ledger
`U_rad = 3·P_rad·V`, use `data.rad_pressure_true` (pure radiation) to
avoid double-counting electron pressure (which is already in `U_plasma`).

---

## 4. Laser diagnostics (optional)

Loaded when the simulation includes a laser drive. Direct-drive runs
populate these; indirect-drive runs may not.

| Pipeline attribute | EXODUS candidates | Shape | Unit | Description |
|---|---|---|---|---|
| `laser_energy_deposited` | `EnLaserDepositedTimeIntg` | (n_t,) | J | Cumulative laser energy deposited up to time t. Final value ≈ absorbed energy |
| `laser_power_source` | `LaserPwrSrc` | (n_t, n_z) | W/cm³ | Laser power deposited per unit volume per zone. Spatial profile of absorption |
| `laser_power_delivered` | `LaserPwrDeliveredForBeam` | (n_t, n_beam) → squeezed to (n_t,) | W | Laser power crossing the grid outer boundary into the plasma, per beam. Helios stores with trailing beam axis; `data_builder` squeezes for n_beam=1 |
| `laser_power_on_target` | `LaserPwrOnTargetForBeam` | (n_t, n_beam) → squeezed to (n_t,) | W | Laser power incident at the grid outer boundary (before absorption). Used for Beer-Lambert reconstruction in `laser_intensity.py` |
| `laser_attenuation_coeff` | `laserAttinuationCoeff` | (n_t, n_beam, n_z+2) | cm⁻¹ | Inverse-Bremsstrahlung absorption coefficient α per zone. Has TWO ghost zones: leading ghost at r=0 (index 0), plus opaque-sentinel ceiling (`1e30`) inside the critical surface. Critical to understanding Method 2 reconstruction — see "laserAttinuationCoeff conventions" section below |
| `electron_density` | `electron_density`, `elec_density`, `n_e` | (n_t, n_z) | cm⁻³ | Free-electron number density. Used for critical-surface location (primary method): `ne >= ncr = 1.115e21 / λ_μm²` |

### `laserAttinuationCoeff` conventions (CRITICAL)

Mis-indexing this array silently breaks diagnostics. Key points:

- Shape is `(n_t, n_beam, n_z+2)` — two extra zones relative to
  standard zone-centered arrays.
- **Index 0**: leading ghost at r=0 (ignore).
- **Opaque region**: zones inside the critical surface contain
  `1e30` (sentinel, not a physical value). Treat as "opaque — no
  light reaches here."
- **Coronal region**: zones outside the turning point contain
  physically small α (~1e-7 cm⁻¹ or less). Numerical noise
  dominates below ~1e-2 cm⁻¹.
- The absorbing layer is the narrow band where α is physically
  meaningful (typically ~0.1 to ~100 cm⁻¹). This is where Method 1
  (`I = P_src/α`) is valid.
- Method-1 filter in `helios_postprocess/laser_intensity.py`:
  `α > ALPHA_MIN_M1 = 1e-2 cm⁻¹`. Below that, treated as noise.

---

## 5. Fusion and neutron diagnostics (optional)

Populated when `BURN` is enabled in the input deck.

| Pipeline attribute | EXODUS candidates | Shape | Unit | Description |
|---|---|---|---|---|
| `fusion_power` | `FusionRate_DT_nHe4`, `fusion_rate_DT`, `fusion_power`, `neutron_rate` | (n_t, n_z) | **reactions/s/g** (mass-specific per zone) | DT fusion reaction rate, normalized per gram of plasma in each zone. CORRECTED May 27 2026 via zero-D verification (`examples/verify_zero_d.py`): a uniform-T uniform-ρ static sphere reports a per-zone value that matches the Bosch-Hale per-gram analytic to 0.04%, definitively NOT a per-zone-per-second or per-cm³-per-second quantity. To get total reactions/s: `np.sum(fusion_power[t] * zone_mass[t])`. For per-gram bulk rate: `np.average(fusion_power[t], weights=zone_mass[t])`. See CLAUDE.md Physics Convention §5b. |
| `neutron_production_rate` | `neutron_production_rate`, `neutron_rate`, `FusionRate_DD_nHe3` | (n_t, n_z) | reactions/cm³/s | DD fusion reaction rate density |
| `alpha_heating_power` | `alpha_heating_power`, `alpha_power`, `alpha_heating` | (n_t, n_z) | W/cm³ | Alpha-particle heating power density. NOT LOADED in PDD_26b_burn — pipeline falls back to `alpha_heating_ion + alpha_heating_ele` |
| `alpha_heating_ion` | `pt_particle_heating_ion` | (n_t, n_z) | W/cm³ | Alpha-particle heating of ions |
| `alpha_heating_ele` | `pt_particle_heating_ele` | (n_t, n_z) | W/cm³ | Alpha-particle heating of electrons |
| `dt_neutron_count` | `TimeIntFusionProd_n_1406` | (n_t,) | count | Cumulative DT neutron count (time-integrated, full volume) |
| `dd_neutron_count` | `TimeIntFusProd_n_0245` | (n_t,) | count | Cumulative DD neutron count |
| `dt_neutron_count_zone` | `TimeIntFusionProd_n_1406_zone` | (n_t, n_z) | count | Per-zone cumulative DT neutron count |

**Yield and gain calculation**: The pipeline uses
`dt_neutron_count[-1] * 14.06 MeV/neutron` (converted to MJ) as the
authoritative yield. This is Helios's internal time-integral —
strictly superior to re-integrating from the EXODUS-sampled
`fusion_power`, because `fusion_power` is only sampled at the output
timesteps (not every simulation step).

---

## 6. Geometry and region metadata (optional but commonly present)

Loaded via `_try_load_nc` because the EXODUS names contain spaces or
are non-standard identifiers.

| Pipeline attribute | EXODUS variable | Shape | Unit | Description |
|---|---|---|---|---|
| `region_interfaces_indices` | `Indices at region interfaces` | (n_t, n_regions) | int | Zone indices at which one region ends and the next begins. For a 4-region target (Vapor/Solid/Foam/CH), shape is (n_t, 4). `ri[t, 0]` = inner edge of solid DT = outer edge of vapor, etc. **Note:** `ri[0, k]` is time-invariant for Lagrangian zones — zones don't change region membership |
| `material_index` | `Material index` | (n_z,) | int | Per-zone material ID. Static in time (Lagrangian). In PDD_26b: three unique values (1, 2, 3). Mapping (to verify): 1 = DT, 2 = DT-CH foam, 3 = CH skin |
| `region_names` | `name_spatial_regn` | list of str | — | Human-readable region names. Example for PDD_26b: `['DT Vapor', 'DT Solid', 'DT-CH foam', 'CH Skin']` |

Notes on regions vs materials:
- A 4-region design (Vapor, Solid, Foam, CH) can have only 3 materials
  because DT Vapor and DT Solid share chemistry (both `material_index = 1`).
- For "is the ablation front in fuel yet?" questions, use
  `material_index[ablation_front_indices[t]]` — the classification is
  cleaner than region-based.

---

## 7. Time axis

| Pipeline attribute | Source | Shape | Unit | Description |
|---|---|---|---|---|
| `time` | `time_whole` (netCDF variable), converted on load | (n_t,) | ns | Simulation time. Stored in EXODUS as SECONDS; `data_builder` multiplies by 1e9 on load when `time_unit='s'` is passed. **All downstream code assumes `data.time` is in ns.** Do NOT multiply by 1e9 again — a known bug in `_compute_shock_breakout` did this and produced zero shock times; was fixed in Stage 2 |

---

## 8. Derived arrays (computed by `data_builder`, not from EXODUS)

| Pipeline attribute | Formula | Shape | Unit | Description |
|---|---|---|---|---|
| `zone_centers` | `(zone_boundaries[:, :-1] + zone_boundaries[:, 1:]) / 2` | (n_t, n_z) | cm | Radial position of zone centers |
| `scale_length` | `\|ρ / (dρ/dr)\|` via `np.gradient` | (n_t, n_z) | cm | Density scale length. Useful for ablation-front identification and for adiabat evaluation |

---

## 9. Unit conversions in the pipeline

Whenever you read a value from an `ICFRunData` attribute, it is typically
in EXODUS/CGS units. The analysis code and summary outputs report in
these derived units:

| Quantity | EXODUS (internal) | Reported/derived unit | Conversion factor |
|---|---|---|---|
| Time | s (via EXODUS) | ns (on data.time) | × 1e9 |
| Temperature | eV | keV (scalars) | × 1e-3 |
| Pressure | J/cm³ | Gbar (scalars) | × 1e-8 |
| Velocity | cm/s | km/s (scalars) | × 1e-5 |
| Length | cm | µm (when small; cm for shell radius etc.) | × 1e4 for µm |
| Density | g/cm³ | g/cm³ or g/cc (unchanged) | — |
| Areal density (ρR) | ∫ρ dr (g/cm²) | g/cm² (unchanged) | — |

---

## 10. What's NOT loaded (notable absences)

Helios output includes many more variables than the pipeline currently
consumes. If a new diagnostic is needed, check the EXODUS file first
before computing indirectly. Common variables often present but not
loaded:

- `zone_temperature` — radiation temperature (if separate from ion/elec)
- `heat_flux_*` — conductive heat fluxes, useful for flux-limiter
  validation (currently only the parameter value is read from `.rhw`)
- `opacity` / `absorption_coeff` by frequency group (radiation transport)
- Particle-specific densities: `n_D`, `n_T`, `n_He4` (could be useful
  for isotope tracking in mixed-ablator designs)

To inspect a specific `.exo` file's full variable list:

    import netCDF4
    ds = netCDF4.Dataset('path/to/run.exo')
    for v in sorted(ds.variables):
        print(f"{v:40s} {ds.variables[v].shape}")

For radiation/escape/boundary-tally discovery specifically:

    candidates = [n for n in ds.variables
                  if any(k in n.lower() for k in
                         ['rad', 'escape', 'bound', 'flux', 'lost',
                          'timeint', 'enrad', 'rmax'])]
    for n in sorted(candidates):
        v = ds.variables[n]
        print(f'  {n:55s}  shape={v.shape}  units={getattr(v,"units","?")}')

---

## 11. ICFRunData non-EXODUS attributes (from .rhw parsing)

These come from the companion `.rhw` input file (via `rhw_parser.py`),
NOT from the EXODUS output. Listed here because they appear alongside
the EXODUS-derived attributes in `ICFRunData`.

| Attribute | Source in .rhw | Unit | Description |
|---|---|---|---|
| `laser_wavelength_um` | laser block | µm | Laser wavelength. Default 0.351 for Nd-glass 3ω; 0.248 for KrF |
| `laser_spot_size_cm` | laser block (beam 1) | cm | Spot size (1/e radius for Gaussian, else uniform radius). For multi-beam configs, this is the beam-1 value; full per-beam data on `laser_geometry_per_beam` |
| `laser_half_cone_angle_deg` | laser block | deg | Half-angle of cone |
| `laser_focus_position_cm` | laser block | cm | Axial focus offset from capsule center |
| `laser_power_multiplier` | laser block | — | Scalar multiplier applied to pulse shape |
| `laser_spatial_profile` | laser block | str | Profile type: "uniform", "gaussian", etc. |
| `laser_foot_power_TW` | pulse shape | TW | Power during foot-pulse phase |
| `laser_peak_power_TW` | pulse shape | TW | Peak power in drive phase |
| `laser_foot_start_ns`, `laser_foot_end_ns`, `laser_peak_start_ns`, `laser_peak_end_ns` | pulse shape | ns | Pulse timing breakpoints |
| `flux_limiter` | heat flux block | — | Flux-limiter coefficient f (typically 0.06 or 0.08). For multi-region configs, this is the first region's value; see `flux_limiter_per_region` |
| `flux_limiter_enabled` | heat flux block | bool | Whether flux limiter is active |
| `flux_limiter_per_region` | heat flux blocks per region | list of dicts | Per-region FL values: `[{region, enabled, value}, ...]`. Surfaces variation across regions (e.g., HDD26: gas=0.06, ice/foam/CD=0.04). Reported in summary output as "Flux limiter (f) per region" table |
| `laser_geometry_per_beam` | laser blocks per beam | list of dicts | Per-beam ray-trace geometry: `[{beam_id, wavelength, spot_size, half_cone_angle, focus_position, power_multiplier, spatial_profile}, ...]`. Surfaces multi-beam configurations (e.g., HDD26 focal zooming with 3 beams). Reported in summary output as "RAY-TRACE GEOMETRY (per beam)" table |
| `drive_temperature` | Rad Source Data block | (n_points,) eV | Indirect-drive radiation-temperature profile (prepulse). Time array also stored. Sentinel times (t ≥ 1 ms ≈ 1e18 s) are filtered out by parser as "extend last value indefinitely" markers |

---

## 12. Boundary-tally and global cumulative quantities (May 2026 addition)

These are scalar-per-timestep `(n_times,)` quantities that measure energy
**leaving the grid** or **summed over the whole grid**. They're the
right-hand-side terms in the global energy-balance closure equation and
are consumed by `examples/energy_balance_diagnostic.py`.

### Currently loaded (wired into `data_builder`)

| Pipeline attribute | EXODUS variable | Shape | Unit | Description |
|---|---|---|---|---|
| `radiation_energy_at_boundary_cum` | `TimeIntRadiationLossAtBds` | (n_t,) | J cumulative | Total radiation energy that has crossed the grid boundaries (Rmin + Rmax combined) up to time t. Used directly in the energy-ledger closure as the rad-escape channel. Confirmed J despite empty units field by closure agreement with computed channels. |
| `particle_energy_escaped_cum` | `particle_time_int_energy_escaped` | (n_t,) | J cumulative | Total kinetic energy of escaped charged particles up to time t. **Does NOT include neutrons** — verified by VI_6 reporting ~0 kJ despite 28.8 MJ DT fusion. Semantically ambiguous across Helios versions; the energy-balance diagnostic uses it only when it captures > 30% of fusion energy, otherwise falls back to nominal 14.1/17.6 DT neutron fraction. |

### Energy-balance closure equation

After Section 12 wiring, the diagnostic computes:

```
E_absorbed + E_alpha_deposited
  = (KE_inward + KE_outward + U_plasma + U_rad)   ← in-plasma channels
  +  E_rad_boundary                                ← direct measurement
  +  residual_gap                                  ← closure residual
```

where:
- `E_alpha_deposited = E_fusion_total × (3.5/17.6)` (≈ 0.199) — nominal
  alpha fraction. Used since `particle_time_int_energy_escaped` doesn't
  track neutron escape consistently.
- `E_neutron_escaped = E_fusion_total × (14.1/17.6)` (≈ 0.801) —
  treated as a separate escape channel.
- `U_plasma = Σ (3/2)(P_ion + P_elec) V_zone` — ideal-gas approximation.
  Slightly overstates in dense ignited regions where EOS-internal
  ionization energy matters; residual_gap is typically small negative
  there (~ −1% of absorbed for VI_6 stagnation).
- `U_rad = Σ 3·P_rad_true·V_zone` — uses pure radiation pressure
  (`data.rad_pressure_true`), NOT the 3-component `data.rad_pressure`
  that mixes in electron pressure.
- KE inward/outward split by sign of zone-centered velocity
  (Helios convention: v < 0 = inward).

`residual_gap` should be < 1% of absorbed for a well-closed ledger.
Observed values across runs tested so far:

| Run | Snapshot | Residual (% of absorbed) |
|---|---|---|
| VI_6 | peak v | +0.6% |
| VI_6 | stagnation | −0.8% |
| HDD26 lrm1_burn | peak v | +0.9% |
| HDD26 lrm1_burn | stagnation | +0.9% |
| HDD26 lrm4_nb | peak v | +0.9% |
| HDD26 lrm4_nb | stagnation | +1.5% |

### Available but not yet wired — wishlist for further accounting

These EXODUS variables have been identified by inspection but are not yet
mapped to `ICFRunData` attributes. Wire as needed when the corresponding
accounting becomes important.

| EXODUS variable | Shape | Likely units | Use case |
|---|---|---|---|
| `FreqIntgRadEnLossRmax` | (n_t,) | J/s? | Frequency-integrated radiation power loss at outer boundary only. Lets us split the combined `TimeIntRadiationLossAtBds` into outer-vs-inner contributions. |
| `FreqIntgRadEnLossRmin` | (n_t,) | J/s? | Same, inner boundary. For most direct-drive runs this should be ~0 (no inner sink); useful as a sanity check. |
| `EnTotRadiation` | (n_t,) | J | Total radiation energy resident in the grid. Cross-check on our `U_rad = 3·P_rad·V` ideal-blackbody calculation. |
| `EnRadSinkTimeIntg` | (n_t,) | J cumulative | Cumulative radiation sink energy in zones. Useful for tracking absorbed-rad / re-emitted-rad balance. |
| `EnExchEleToRadTimeIntg` | (n_t,) | J cumulative | Cumulative electron-to-radiation energy exchange. Indicates how much plasma thermal energy converts to in-grid radiation (vs. radiation that just escapes). |
| `EnJouleHeatingTimeIntg` | (n_t,) | J cumulative | Cumulative Joule heating. Source term if magnetic effects are non-negligible (usually small in ICF). |
| `EnMagTimeIntg` | (n_t,) | J cumulative | Cumulative magnetic-field energy. Same — usually negligible but should close if non-zero. |
| `HeatFluxAtRegionBdEle` | (n_t, n_regions) | J/(s·cm²)? | Electron heat flux at the inter-region boundaries. Enables per-region energy ledger (corona vs ablator vs fuel). |
| `HeatFluxAtRegionBdIon` | (n_t, n_regions) | J/(s·cm²)? | Ion heat flux at inter-region boundaries (same use case). |
| `RadNetCoolingRateRegion` | (n_t, n_regions) | J/(s·cm²) | Net radiation cooling rate per region. Useful for understanding where the trapped coronal thermal energy is converted to radiation. |
| `TimeIntFusionProd_He4_0352` | (n_t,) | reactions | Cumulative DT alpha count (3.52 MeV alphas). Direct cross-check on our DT-neutron-derived alpha-energy estimate. |
| `TimeIntFusionProd_He4_0367` | (n_t,) | reactions | Cumulative DD alpha count (3.67 MeV branch). |
| `TimeIntFusionProd_T_0101_zone`, `TimeIntFusionProd_p_*_zone`, etc. | (n_t, n_z) | reactions | Per-zone fusion-product counts — enable spatial fusion-burn diagnostics. |
| `LaserEnDeliveredTimeInt` | (n_t,) | J/geofac | Cumulative laser energy delivered at outer boundary (vs. absorbed). Allows direct frac-absorbed time-series tracking. |
| `rad_energy_density` | (n_t, n_z+2) | J/cm³ | Per-zone radiation energy density (node-centered with two ghosts). Direct alternative to `U_rad = 3·P_rad·V` — bypasses the blackbody-equilibrium assumption. |
| `radiation_temperature` | (n_t, n_z) | eV | Per-zone radiation temperature. Useful for radiation-front tracking and for cross-checking blackbody assumption against actual rad field. |

---

## Version

- Generated: April 2026, end of Task 3.A Stage 1 session.
- May 2026 updates: added Section 12 (boundary-tally and global cumulative
  quantities); added `rad_pressure_true`, `flux_limiter_per_region`,
  `laser_geometry_per_beam`, `drive_temperature` to existing sections;
  noted EXODUS particle-escape variable does NOT include neutrons.
- Authoritative source: `helios_postprocess/data_builder.py` at the
  corresponding session date.
- Verified against: `Olson_PDD_26b_burn` (original), `VI_6`,
  `HDD26_DTI40_1ns130_FL04_lrm4_nb` (May 2026 updates).
- Re-generate if `_VARIABLE_MAP` or `_try_load_nc` calls are modified.
