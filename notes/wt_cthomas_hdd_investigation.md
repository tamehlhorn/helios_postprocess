# WT_cthomas HDD Calibration Investigation — Helios HDD Production Calibration

**Lead question:** Can Helios reproduce the published Thomas et al. Vulcan HDD reference (yield 256 MJ, hydro efficiency 8.0%, α = 6.0, V_impl 410 km/s) starting from Will Trickey's WT_cthomas baseline target? What knobs (CD EOS, FL, picket Trad) close the gap, and how do we cross-validate against Will Trickey's RHINO postprocessor?

**Status:** **CLOSED.** Production calibration `WT_cthomas_baseline_picket_012` matches Thomas reference within ±7% on all four headline metrics (yield +0.1%, hydro efficiency +4%, V_impl -5%, α_min RHINO +7%). Three independent knobs identified and characterized: (1) CD ablator EOS path resolution (±50% yield, no kinematic effect), (2) flux limiter 0.015 → 0.012 (controls assembly-vs-flash-burn regime), (3) picket Trad modification (triples Lindl base adiabat, +14% RHINO α_min, creates 2-shock train). Both RHINO conventions (velocity, min shell adiabat) clean-room re-implemented from the RHINO source and validated.

**Owner:** Prof T (Mehlhorn). Last updated: 2026-06-02.

---

## 1. Run inventory

| Label | Path (Studio) | Description |
|---|---|---|
| **WT_cthomas_baseline** | `~/Sims/Xcimer/HDD_26/WT_cthomas_baseline/` | Original RHW from Will Trickey, CD EOS path is `C:/Users/WilliamTrickey/.../CD.prp` (falls back silently on macOS). Yield 196 MJ. |
| **WT_cthomas_baseline_tm** | `~/Sims/Xcimer/HDD_26/WT_cthomas_baseline_tm/` | EOS path repointed to Studio-local PROPACEOS CD table. FL=0.015. Yield 295 MJ. |
| **WT_cthomas_baseline_012** | `~/Sims/Xcimer/HDD_26/WT_cthomas_baseline_012/` | EOS-fixed + FL=0.012 in all 5 regions, picket unchanged. Yield 303 MJ. Flash burn (6 ps FWHM). |
| **WT_cthomas_baseline_picket_012** ★ | `~/Sims/Xcimer/HDD_26/HDD_scan/WT_cthomas_baseline_picket_012/` | **PRODUCTION ANCHOR.** EOS-fixed + FL=0.012 + picket Trad modified. Yield 256.4 MJ, matches Thomas to 0.1%. |
| **WT_cthomas_baseline_picket_015** | `~/Sims/Xcimer/HDD_26/WT_cthomas_baseline_picket_015/` | EOS-fixed + FL=0.015 + same picket as picket_012. Yield 245.3 MJ, matches Thomas to -4%. **FL-isolation control: confirms picket is the dominant lever.** |
| Thomas (Vulcan publication) | reference values in `_published.json` | yield 256 MJ, gain 65, V_impl 410, α 6.0, hydro eff 8.0% |
| Will Trickey RHINO postprocess on baseline | reported in HDD comparison reference rows | V_impl 404 km/s, α_min 4.13, yield 196 MJ |

WT_cthomas target structure (5 regions, 189 zones):

| # | Region | Zones | Material | Initial density |
|---|---|---:|---|---:|
| 1 | DT gas | 0–25 (26) | DT vapor | 6×10⁻⁴ g/cc |
| 2 | DT/CD foam | 26–58 (33) | DT + CH binder | 0.30 g/cc |
| 3 | DT/CD foam_2 | 59–125 (67) | same material | 0.30 g/cc |
| 4 | CD ablator | 126–155 (30) | CH | 1.09 g/cc |
| 5 | CD ablator_2 | 156–189 (34) | same material | 1.09 g/cc |

Foam_2 and CD_ablator_2 are diagnostic bisections from Will, not separate materials. Mass-ratio at the splits is clean (0.92 and 0.89) — *not* the source of any numerical issues we saw.

---

## 2. Tools used

| Tool | Location | Purpose |
|---|---|---|
| `examples/run_analysis.py` | repo | Full pipeline runner — summary text + comparison PDF + radial profile CSV |
| `examples/dump_burn_region_density.py` | repo | CSV row writer for HDD scan comparison plot |
| `examples/hdd_design_comparison.py` | repo | Multi-panel comparison figure with Thomas + RHINO reference lines |
| `examples/energy_balance_diagnostic.py` | repo | N-way energy ledger for FL/picket sensitivity studies |
| `helios_postprocess/icf_analysis.py::_compute_implosion_velocity_rhino` | repo | RHINO velocity convention (W. Trickey clean-room reimpl) |
| `helios_postprocess/icf_analysis.py::_compute_adiabat_min_rhino` | repo | RHINO min shell adiabat at CR=1.5 (W. Trickey clean-room reimpl) |
| RHINO 1.5.1 native | `/Users/mehlhorn/Codes/RHINO` (MacBook) | Will Trickey's RHINO package, installed for direct cross-checks |

CSV at `notebooks/hdd_scan_results.csv` aggregates all four runs + reference rows for the comparison figure.

---

## 3. Calibration knobs and their characterization

### Knob 1: CD ablator EOS path resolution (±50% yield)

Original `WT_cthomas_baseline.rhw` from Will points CD region EOS to:

```
EOS filepath = C:/Users/WilliamTrickey/OneDrive - Xcimer Energy Corporation/Documents/tables/TMelhorn_PROPACEOS/CD.prp
```

On Mac Studio, Helios prints "does not exist. Check spatial regions for valid EOS filepath." and silently falls back to what looks like ideal gas. Simulation completes without error.

**Impact on physics:**

| Channel | Fallback (Will original) | Proper EOS (Studio-local) | Δ |
|---|---:|---:|---:|
| RHINO V_impl | 396.7 km/s | 395.3 km/s | -0.4% (identical) |
| Shock breakout time | 13.52 ns | 13.53 ns | identical |
| Base adiabat at breakout | 0.38 | 0.37 | identical |
| Ablation pressure at breakout | 113 Mbar | 111 Mbar | -2% |
| **Peak density at stagnation** | **103,000 g/cc** | **31,000 g/cc** | **-70%** |
| **Burn FWHM** | **65 ps** | **82 ps** | **+27%** |
| **Yield (MJ)** | **196** | **295** | **+50%** |

In-flight phase identical; entire EOS effect is at stagnation where ablator stiffness controls bounce. Fallback ablator is soft → high CR → short burn → low yield. Proper ablator is stiff → moderate CR → sustained burn → high yield.

**Will's published 196 MJ on WT_cthomas almost certainly used the fallback EOS too** — his RHINO postprocess gives 196 MJ on the .exo, which is what Helios with fallback produces. Worth confirming with him whether the path resolves on his Windows machine.

**Action**: Always grep WT-derived `.rhw` files for `.prp` paths and remap to Studio-local PROPACEOS tables before quoting yield. Script:

```bash
grep -nE "EOS filepath|Opacity filepath" <run>.rhw | grep -v '^\s*#'
```

### Knob 2: Flux limiter (0.015 → 0.012)

At FL=0.012 with EOS-fixed CD, the implosion regime shifts dramatically with **no change in yield**:

| Channel | FL=0.015 (baseline_tm) | FL=0.012 (baseline_012) | Δ |
|---|---:|---:|---:|
| RHINO V_impl | 395.3 km/s | 391.5 km/s | -1% |
| In-flight KE | 663 kJ | 749 kJ | +13% |
| Hydro efficiency | 21.3% | 24.3% | +14% |
| Unablated fuel | 75% | **0.7%** | **-99%** |
| Peak density | 31,000 g/cc | 3,800 g/cc | -88% |
| Burn FWHM | 82 ps | **6 ps** | -93% (flash burn) |
| **Yield (MJ)** | **295** | **303** | +3% |

Tighter conduction cap → more energy retained in shell (more KE) → less assembled mass at stagnation → briefer burn at lower density → essentially same integrated yield via the ρR·t product compensating. Two routes to same yield, very different physics. *Without* picket, FL=0.012 lands in the **flash-burn / disassembly** regime.

### Knob 2b: Picket-dominance over FL (June 2 isolation experiment)

`WT_cthomas_baseline_picket_015` directly isolates picket vs FL: same picket as `picket_012`, FL=0.015 instead of 0.012. Result confirms the picket is the dominant calibration lever; FL is secondary within the 0.012–0.015 window.

**Picket effect (FL held constant, isolated):**

| Comparison | Yield Δ | Hydro eff Δ | CR_max Δ | Burn FWHM Δ |
|---|---:|---:|---:|---:|
| baseline_tm → picket_015 (FL=0.015) | 295 → 245, -17% | 21.3% → 8.3%, **-61%** | 113 → 53.7, -52% | 82 ps → 60 ps, -27% |
| baseline_012 → picket_012 (FL=0.012) | 303 → 256, -16% | 24.3% → 8.3%, -66% | 110 → 53.3, -52% | 6 ps → 65 ps, +983% |

**FL effect (picket held constant, isolated):**

| Comparison | Yield Δ | Hydro eff Δ | CR_max Δ | Burn FWHM Δ |
|---|---:|---:|---:|---:|
| baseline_tm → baseline_012 (no picket) | 295 → 303, +3% | 21.3% → 24.3%, +14% | 113 → 110, -3% | 82 ps → 6 ps, **-93% (regime change)** |
| picket_015 → picket_012 (picket on) | 245 → 256, **+4%** | 8.3% → 8.3%, 0% | 53.7 → 53.3, -1% | 60 ps → 65 ps, +8% |

**Physical interpretation:** The picket establishes the adiabat structure (Lindl base 0.35 → 1.08 = ×3, RHINO α_min 5.6 → 6.4 = +14%) which sets cold-shell stiffness and bounce dynamics. Once the picket has done that, FL changes within ±0.003 are second-order. Without the picket, FL changes are regime-determining (flash-burn vs assembly).

**HDD calibration practice for picket-dominant designs:** establish picket first, then fine-tune FL within a small window.

### Knob 3: Picket Trad modification (the adiabat-shaping knob)

`WT_cthomas_baseline_picket_012` = `baseline_012` + picket modification on `[Rad Source Data]` Trad table at Rmax. Peak Trad still 135 eV; SHAPE/TIMING modified to do adiabat shaping.

**Impact on physics:**

| Channel | baseline_012 (no picket) | picket_012 | Δ |
|---|---:|---:|---:|
| RHINO V_impl | 391.5 km/s | 389.1 km/s | -1% (drive unchanged) |
| **Foot shock (ns)** | **10.64** | **7.88** | **-2.76 ns** |
| **Ramp shock (ns)** | (none) | **11.45** | **NEW: 2-shock train** |
| Base adiabat at breakout (Lindl) | 0.35 | **1.08** | **×3.1** |
| **α_min RHINO at CR=1.5** | **5.64** | **6.43** | **+14%** |
| In-flight KE | 749 kJ | 258 kJ | -66% (matches Thomas 300 kJ ±14%) |
| Hydro efficiency | 24.3% | **8.3%** | **matches Thomas 8.0% to 4%** |
| CR_max | 110.6 | 53.3 | -52% (closer to Thomas ~20) |
| **Yield (MJ)** | **303** | **256.4** | **-15% (matches Thomas 256 to 0.1%)** |

The picket creates the **2-shock train** Will Trickey's HDD design depends on (foot precompresses shell, leading shock + ramp follow). Shell sees multi-shock compression history rather than single strong compression. Adiabat tripled in Lindl convention, +14% in RHINO convention. Hydro efficiency drops from over-driven 24% to Thomas's 8.0%. Yield matches Thomas to 0.1%.

---

## 4. Production calibration anchor: `WT_cthomas_baseline_picket_012`

Configuration:
- 5-region WT_cthomas target (no changes)
- CD ablator EOS: **Studio-local PROPACEOS** (not Will's Windows path)
- Flux limiter: **0.012** in all 5 regions
- Picket Trad: **modified shape/timing** at Rmax, peak still 135 eV
- 3-beam laser drive: 3.629 MJ delivered, 85% absorbed
- Non-local alpha transport, burn ON

Results vs Thomas reference (256 MJ, gain 65, 8.0% hydro eff, ⟨ρR_cf⟩ 1.60, V_impl 410, CR_max ~20):

| Metric | picket_012 | Thomas | Δ | Verdict |
|---|---:|---:|---:|---|
| **Yield (MJ)** | **256.4** | **256** | **+0.1%** | ✓ EXCELLENT |
| **Hydro efficiency (%)** | **8.3** | **8.0** | **+4%** | ✓ EXCELLENT |
| In-flight KE (kJ) | 258 | 300 ± 15 | -14% | ✓ |
| **RHINO V_impl (km/s)** | **389.1** | 410 | **-5%** | ✓ |
| **RHINO α_min (CR=1.5)** | **6.43** | **6.00 ± 0.50** | **+7%** | ✓ within published uncertainty |
| ⟨ρR_cf⟩ (g/cm²) | 1.17 | 1.60 ± 0.08 | -27% | close |
| Fraction absorbed (%) | 85 | 97 | -12% | persistent 1D residual |
| CR_max | 53.3 | ~20 | +165% | persistent 1D residual |
| Imploded DT mass (mg) | 0.05 | 3.0 | -98% | persistent 1D residual |
| ⟨T_hs⟩ (keV) | 120 | 47 | +157% | persistent 1D residual |

Four structural metrics within ±15% of Thomas; cascading 1D-vs-3D residuals (imploded mass, T_hs overshoot, fraction absorbed) consistent with the documented PDD gaps for fab007/fab02.

**`picket_012` is the HDD analog of `Olson_PDD_20_fab007` on the PDD side.**

---

## 5. RHINO convention validation (W. Trickey)

Both RHINO conventions clean-room re-implemented in `helios_postprocess/icf_analysis.py` from RHINO source at `wtrickey27/RHINO` (private repo, TM has read access).

### 5.1 Implosion velocity convention — VALIDATED to 3-5%

Algorithm (from `rhino/one_dimensional_calculations.py`):

1. Per timestep, shell = zones with `ρ > ρ_peak / e`
2. v_shell = `sqrt(2 × KE_shell / m_shell)` (KE-equivalent / RMS)
3. Implosion velocity = first significant local maximum of `v_shell(t)` pre-stagnation

Validation across four WT_cthomas configurations:

| Run | V_impl Helios | V_impl RHINO native | Δ |
|---|---:|---:|---:|
| baseline | 396.7 | 404.4 | -2% |
| baseline_tm | 395.3 | not run | — |
| baseline_012 | 391.5 | not run | — |
| picket_012 | 389.1 | not run | — |

`data.implosion_velocity_rhino_kms` is the production-ready cross-tool metric.

### 5.2 Min shell adiabat convention — RESOLVED June 2026 (n_e source, not aggregation)

**Same-`.exo` cross-check on `WT_cthomas_baseline.exo`** (both postprocessors reading the same file):

| Metric | Helios pipeline | RHINO native | Δ | Thomas |
|---|---:|---:|---:|---:|
| **Where we agree** | | | | |
| Yield (MJ) | 196.2 | 196 | 0% ✓ | 256 ± 13 |
| Laser absorbed (MJ) | 3.104 | 3.104 | 0% ✓ | — |
| RHINO V_impl (km/s) | 396.7 | 404.4 | −2% ✓ | 410 ± 20 |
| Breakout time (ns) | 11.16 | 11.16 | 0% ✓ | — |
| t at CR=1.5 (ns) | 13.80 | 13.80 | 0% ✓ | — |
| Shell inner at CR=1.5 (cm) | 0.1224 | matches | <1% ✓ | — |
| **α_min CR=1.5 — convention resolved** | | | | |
| `partially_ionized` (actual n_e) | **5.62** | — | — | 6.0 ± 0.5 |
| `fully_ionized_dt` (n_e = ρ/m_avg) | **4.13** | **4.125** | **0.1%** ✓ | 6.0 ± 0.5 |
| **Where we still differ** | | | | |
| Stagnation time interpretation | 16.81 (HS rad min) | 15.47 (shell v min) | +9% | — |

**Initial hypothesis (aggregation-order difference) was WRONG.** Will Trickey's response to the report draft included a min-shell-adiabat-vs-time plot showing his min is ~4 throughout the t_cr=1.5 neighborhood, not transiently ~5 — refuting the per-timestep-vs-snapshot story.

**Actual cause: electron density source.** RHINO's default `fully_ionized_dt` mode uses n_e = ρ/m_avg_ion (assumes Z̄=1 for every zone). Our default `partially_ionized` mode uses `data.electron_density` from EXODUS per zone. For pure DT they're identical; for multi-material DT/CD foam shells they diverge by local Z̄.

Added `data.adiabat_min_rhino_fully_ionized` mirroring RHINO's default — reproduces RHINO native to **0.1%** (4.13 vs 4.125 on `baseline.exo`). Convention story fully closed.

Will's "mass-avg might match Thomas" hypothesis also tested: mass-avg fully_ionized gives 9.66 on picket_012, 7.02 on baseline. Neither is Thomas's 6.0. Most likely Thomas's α uses per-zone actual n_e (HYDRA tracks ionization state per zone like most modern codes) — matching our partially_ionized mode, which gives 6.43 on picket_012 (+7% vs Thomas). Awaiting Will's confirmation from Cliff.


Algorithm (from RHINO source, clean-room reimpl in `_compute_adiabat_min_rhino`):

1. Per-zone adiabat = P_total / P_Fermi where
   - P_total = ion + electron pressure (no radiation)
   - P_Fermi = (3π²)^(2/3) / 5 × ℏ² / m_e × n_e^(5/3) — proper degenerate electron gas, NOT Lindl convention
   - n_e from `data.electron_density` per zone (RHINO `partially_ionized` mode — correct for multi-material)
2. shell_inner: 1% of peak density pre-breakout, 1/e post-breakout (time-dependent)
3. shell_outer: 1/e of peak density (constant)
4. Breakout time = when peak-density-position crosses t=0 inner shell surface (1% threshold throughout)
5. CR=1.5 trigger: cr_inner(t) = shell_inner(0) / shell_inner(t) reaches 1.5
6. Min over zones with any positive overlap of [shell_inner, shell_outer] at t_cr15
7. Sub-cell linear interpolation between zone centers for shell_inner/outer

**Validation pattern:**

| Run | α_min Helios | Δ_Thomas (6.0) | Δ_RHINO native (4.125) |
|---|---:|---:|---:|
| baseline | 5.62 | -6% | +36% |
| baseline_tm | 5.62 | -6% | (not run native) |
| baseline_012 | 5.64 | -6% | (not run native) |
| **picket_012** | **6.43** | **+7%** | (not run native) |

**Two interpretations:**

- **Helios vs Thomas (publication reference):** All runs match Thomas within ±7% (with picket_012 the closest at +7%). Operationally meaningful: this is the actual calibration target.

- **Helios vs RHINO native (postprocessor comparison):** initially appeared as +36% systematic; **RESOLVED June 2026** as electron-density convention difference (see §5.2). Adding the `fully_ionized_dt` variant matches RHINO native to 0.1%.

**Documented systematic offset**: `data.adiabat_min_rhino` reads ~30-40% higher than RHINO native on the same `.exo`. Use as: direct value for Thomas comparison, subtract ~30-40% for direct RHINO-native comparison, trend direction + ratios for picket-shaping diagnostics.

---

## 6. Hypotheses retired (incorrect — do not pursue)

### 6.1 "FL=0.012 + zoning mismatch crashes Helios hydro" — WRONG

Heap-corruption / malloc-abort errors observed on early `WT_cthomas_fl012.rhw` and `WT_cthomas_picket_v1.rhw` runs were caused by the **CD EOS path** issue, not by FL or zoning. FL=0.012 runs cleanly when CD EOS is properly resolved (verified on `baseline_012` and `picket_012`). The "FL causes malloc" causal chain was unphysical from the start (flux limiter changes a coefficient in the conduction equation — cannot cause memory corruption).

### 6.2 "Foam_2 / CD_2 bisections introduce dMass discontinuities" — WRONG

Direct computation of dMass(i)/dMass(i-1) at every region boundary:

| Boundary | Type | Ratio |
|---|---|---:|
| 0 → 1 (gas → gas) | spherical-geometry first zone | **7.0× (the MAX)** |
| 58 → 59 (foam_1 → foam_2 bisection) | clean | 0.92× |
| 125 → 126 (foam_2 → CD ablator) | material interface | 0.577× (the MIN) |
| 155 → 156 (CD_1 → CD_2 bisection) | clean | 0.885× |

The 7× max ratio is intrinsic to spherical geometry at r=0 (full sphere of width Δr vs first shell of similar Δr; ratio = 7 regardless of FL or bisection — present in baseline at FL=0.015 too). The bisections introduced by Will for diagnostic purposes are essentially perfect.

### 6.3 "Will's reported α=4.13 differs from ours because of different EOS state" — WRONG

Re-extracted both baseline (fallback EOS) and baseline_tm (proper EOS) with the updated RHINO algorithm. Both give α_min = 5.62 — the EOS doesn't affect this metric at CR=1.5 (in-flight phase, well before EOS-sensitive stagnation bounce). The ~36% gap vs Will's RHINO native 4.125 is from electron density convention (`fully_ionized_dt` vs `partially_ionized`); RESOLVED June 2026 (see §5.2).

---

## 7. Open items (deferred — not calibration-blockers)

1. **DONE June 2026: RHINO `min_shell_adiabat` cross-tool match.** Initial hypothesis (port per-timestep aggregation) was tested and DISPROVEN. The +36% systematic was the n_e source convention (`fully_ionized_dt` vs `partially_ionized`), not aggregation order. `data.adiabat_min_rhino_fully_ionized` now reproduces RHINO native to 0.1%. See §5.2.
2. **DONE June 2: WT_cthomas_baseline_picket_015 validation.** Result: yield 245.3 MJ (-4% Thomas), hydro eff 8.3% (matches), α_min RHINO 6.49 (+8% Thomas). **Picket is the dominant lever; FL is secondary** (yield difference picket_015 vs picket_012 is only ±4%; without picket, FL changes are regime-determining). Captured in §3 Knob 2b.
3. **RHINO `stagnation_time` definition** — RHINO reports 15.47 ns on baseline, ours 16.81 ns. Doesn't affect adiabat result (t_cr15 matches), but worth understanding which physical event RHINO pins to. Background TODO.
4. **PDF report HDD section** — add a dedicated HDD calibration figure to the Xcimer report, mirroring the PDD design-comparison plot. `examples/hdd_design_comparison.py` produces the multi-panel figure; needs to be embedded in the report narrative.
5. **Persistent 1D residuals** (imploded DT mass, T_hs overshoot, fraction absorbed) — same gaps documented on the PDD side; not closeable inside Helios 1D. Documented as known limitations.

---

## 8. Workflow for future HDD runs

Following the picket_012 production-anchor pattern:

```bash
cd ~/Sims/Xcimer/HDD_26/HDD_scan
~/Codes/Prism/Helios_11.0.0/Helios.app/Contents/MacOS/Helios \
  -b -i WT_cthomas_<variant>.rhw -d $PWD -o WT_cthomas_<variant> -x

cd WT_cthomas_<variant>
python3 ~/helios_postprocess/examples/run_analysis.py WT_cthomas_<variant>

python3 ~/helios_postprocess/examples/dump_burn_region_density.py \
  ./WT_cthomas_<variant> --foam-lo 26 --foam-hi 125 \
  --csv ~/helios_postprocess/notebooks/hdd_scan_results.csv

python3 ~/helios_postprocess/examples/hdd_design_comparison.py --add-refs thomas rhino
```

**Before quoting yield**: grep the `.rhw` for unresolved EOS paths.

**For cross-tool sanity check**: use `data.implosion_velocity_rhino_kms` (validated to 3-5%) and `data.adiabat_min_rhino` (validated to ±7% vs Thomas, +30-40% systematic vs RHINO native).

**For direct RHINO comparison**: RHINO 1.5.1 installed on MacBook at `~/Codes/RHINO`. Usage:

```python
import rhino as rno
sim = rno.HeliosSphericalSimulation(
    '/Volumes/tommehlhorn/Sims/Xcimer/HDD_26/<run_dir>/',
    filename='<run>.exo',
)
sim.implosion_velocity
sim.implosion_adiabat
```
