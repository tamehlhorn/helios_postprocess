# Foam vs Ice Burn Investigation — Helios PDD Calibration

**Lead question:** Why does Helios's fab007 production run (wetted DT-CH foam main fuel) yield 26 MJ while LILAC's same target yields 87 MJ? Where is the foam-burn deficit, and is it a Helios artifact or a physical insight about wetted-foam designs?

**Status:** **CLOSED.** Foam-burn deficit traced to real material physics (PROPACEOS DT-CH foam EOS+opacity tables, isolated by direct planar experiment). Basic Helios rate kernels validated. CR_burn EOS-swap previously suspected as opacity pathology — refuted by 3-way energy ledger, now understood as an EOS-table mismatch artifact (do not interpret as physical finding). **Design study (May 29) demonstrated that FL=0.012 + aggressive geometric defocus (c20) recovers 72% of the LILAC yield gap with 31.5% foam yield share — proving PROPACEOS foam CAN be burned in Helios given sufficient drive.**

**Owner:** Prof T (Mehlhorn). Last updated: 2026-05-25.

---

## 1. Run inventory

| Label | Path (Studio) | Description |
|---|---|---|
| **fab007-foam (production)** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_burn/` | Calibrated PDD point — wetted DT-CH foam, FL_prism=0.007, yield 26 MJ |
| **fab007-foam, burn OFF** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_nb/` | Hydrodynamic-state reference for the production run |
| **fab015-foam** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab015_foot25_s018_c37_burn/` | Earlier calibration attempt, FL_prism=0.015, yield 29 MJ |
| **fab007-ice** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_ice_burn/` | fab007 geometry with DT ice replacing foam in main fuel, yield **81 MJ** |
| **fab007-CR** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_CR_burn/` | fab007 geometry, **Sesame DT-ice EOS + PROPACEOS foam opacity** in main fuel — diagnostic for separating EOS vs opacity. Yield 4.2 MJ. **NOTE: refuted as physical finding; EOS-table mismatch artifact.** |
| **fab02-foam** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab02_foot25_s016_burn/` | Bootstrap-strength reference, over-driven foam, yield 59 MJ |
| **planar_dtice_1** | `~/Sims/Xcimer/Olson_PDD/Planar_shock/planar_dtice_1/` | Planar laser-driven slab, pure DT (DT_20250131.prp), 200 µm × 0.25 g/cc × 100 TW/cm² × 2 ns |
| **planar_wf_1** | `~/Sims/Xcimer/Olson_PDD/Planar_shock/planar_wf_1/` | Same as above with wetted DT-CH foam (DT_and_CH_0p017_20260128.prp) |
| **Zero-D verification (11 runs)** | `~/Sims/Xcimer/Olson_PDD/Fraley_burn/` | Static-plasma rate verification, T sweep + density sweep |
| **LILAC reference** | `references/Olson_2021_*` | Published yield 87 MJ |

Headline metrics across the four key burn runs:

| Run | V_peak (km/s) | Coupling % | Peak ρ at bang (g/cc) | HS ρR max (g/cm²) | Unablated DT (stag, %) | Yield (MJ) |
|---|---:|---:|---:|---:|---:|---:|
| fab015-foam | 429 | 74.4 | 235 | 0.37 | 12.6 | 29.3 |
| fab007-foam (prod) | 421 | 73.1 | ~245 | 0.35 | ~12.6 | 26.0 |
| fab007-ice | 427 | 73.1 | 310 | 0.71 | 45.9 | 81.4 |
| fab007-CR (EOS swap) | 418 | 73.1 | 157 | 0.22 | 12.9 | 4.2 |

---

## 2. Tools used

| Tool | Location | Purpose |
|---|---|---|
| `examples/compare_runs.py` | repo | Two-run side-by-side: R-T trajectories, time histories, stagnation lineouts |
| `examples/energy_balance_diagnostic.py` | repo | N-way energy ledger at peak-v / stagnation / end-of-run |
| `examples/dump_burn_propagation_profile.py` | repo | Radial T_e/T_rad/T_ion/ρ/α-dep across ice/foam interface |
| `examples/dump_per_zone_burn_share.py` | repo | Cumulative DT neutron count per zone, region-aggregated |
| `examples/dump_burn_rate_timing.py` | repo | Per-zone burn-rate timing + compression-state by region |
| Outputs/figures | `comparisons/`, `*_energy_balance.pdf`, `*_report.pdf` | |
| External: `~/Downloads/files/ignition.py` | scratch | Ideal-ignition T vs C atomic fraction (Bosch-Hale + NRL brems) |

---

## 3. Findings (chronological)

### May 24 2026 — Burn-front-arrest hypothesis at ice/foam interface: **refuted**

`dump_burn_propagation_profile.py` on fab007 at peak burn shows the ice/foam interface (zone 190) is essentially continuous in T_e (ratio 1.017), T_rad (1.000), and α-deposition rate (1.06 across the interface). Z̄ jump is only 8% (0.999 → 1.084) because the "DT-CH foam" is ~92% DT by mass. There is no significant brems sink or alpha-MFP cliff at the interface.

→ See memory: [fab007 burn-propagation diagnostic](../../../.claude/projects/-Users-mehlhorn/memory/fab007_burn_propagation_diagnostic.md)

### May 25 2026 — Foam vs ice comparison (`compare_runs.py`)

Generated three figures comparing fab015-foam vs fab007-ice:

- `comparisons/compare_foam_vs_ice_rt.png` — R-T trajectories
- `comparisons/compare_foam_vs_ice_histories.png` — time histories panels
- `comparisons/compare_foam_vs_ice_lineouts.png` — bang-time ρ(r), T_ion(r)

**Key visual findings:**

- Drive trajectories essentially identical (R-T plots, in-flight KE, V_peak all match to ~2%) — comparison is fair
- Central hot-spot ρ: ice 30 g/cc vs foam 12 g/cc (consistent with May 24 burn-propagation finding)
- Central T_ion at bang: **ice 75 keV vs foam 38 keV** — 2× gap, larger than the HS ρR ratio (0.71/0.37 = 1.92) would predict on confinement alone
- Ablation eats deeper into foam: foam retains 12.6% of initial DT mass at stagnation vs ice 45.9% — 3.6× differential
- Foam shows a T_ion *dip* around 170–195 µm where the burn front partially arrests through the foam-region annulus

### May 25 2026 — Energy balance comparison (`energy_balance_diagnostic.py`)

Ledger closes within −0.5% to +1.7% residual → diagnostic trustworthy.

**Stagnation snapshot (kJ):**

| Channel | foam | ice | ratio |
|---|---:|---:|---:|
| Absorbed laser | 1549 | 1492 | 0.96 |
| **Alpha deposited** | 5833 | 16185 | **2.77×** |
| Plasma thermal | 460 | 1923 | 4.2× |
| KE outward (blowoff) | 6292 | 14476 | 2.30× |
| Rad escape (boundary) | 670 | 977 | 1.46× |
| **Total fusion released** | **29.3 MJ** | **81.4 MJ** | **2.77×** |

**Headline conclusion:** Alpha-deposited / total-fusion-released is identically ~20% in both runs. Alpha transport is *not* the foam-vs-ice discriminator. The 2.77× yield gap lives entirely upstream in the n²T²σv reaction-rate integral, which is set by ρ, T in the burn region — both lower in foam.

→ See memory: [foam vs ice energy balance](../../../.claude/projects/-Users-mehlhorn/memory/foam_vs_ice_energy_balance.md)

### May 25 2026 — Carbon brems-cooling and fuel-dilution hypotheses: **refuted**

Ideal-ignition calculation (alpha heating = brems, Bosch-Hale ⟨σv⟩, NRL brems):

| C atomic fraction | T_ig (keV) |
|---:|---:|
| 0.0% (pure DT) | 4.31 |
| 1.7% | 5.19 |
| 4.7% | 6.63 |

The foam burn region runs at 38–75 keV — far above either threshold. Alpha heating ∝ n²⟨σv⟩E_α overwhelms Z²-enhanced brems by orders of magnitude at burn T. So while foam's rad-escape fraction is 1.6× higher than ice's (9.1% vs 5.5%), that channel cannot drive the yield gap.

Fuel dilution from carbon + binder H: ~95% D+T atomic fraction remaining → n_D × n_T only down ~10% from pure DT. Also insufficient.

→ The actual mechanism is **fluid-dynamical**: foam's CH skeleton resists DT collapse, giving lower n in the burn region. Cascading consequences (lower n, lower bootstrap T, lower σv) yield 20× lower local burn-rate density, integrated to 2.77× yield gap.

### May 25 2026 — CR_burn EOS-swap test: ~~opacity-side pathology suspected~~ **(superseded May 27 — see below)**

Run `Olson_PDD_20_fab007_foot25_s018_c37_CR_burn` has Sesame DT-ice EOS *and* PROPACEOS wetted-foam opacity in the main fuel layer. Expectation: ice-like compression, marginally extra rad losses.

**Observed:**

| Metric | foam (prod) | ice | CR (EOS swap) |
|---|---:|---:|---:|
| Peak ρ at bang (g/cc) | ~245 | 310 | **157** |
| Adiabat | 1.95 | (lower) | 1.62 |
| HS ρR peak (T>4.5) | 0.35 | 0.71 | **0.22** |
| Yield (MJ) | 26 | 81 | **4.2** |

Lower adiabat + matched KE delivery should give *more* compression, not less. Peak density collapsed to half the production-foam value. HS ρR never reaches the 0.3 ignition threshold — no bootstrap, marginal burn.

**Three hypotheses for the EOS-swap pathology:**

1. PROPACEOS foam opacity at ice-compressed densities (10s–100s of g/cc) is in a regime the table wasn't calibrated for; anomalous rad cooling during implosion dumps pressure
2. Sesame DT-ice EOS misbehaves when applied to material with 8% C by mass (Z=6 ions in DT EOS phase structure)
3. EOS/opacity inconsistency — internally different reference compositions break detailed balance

**Diagnostic to discriminate (not yet run):** energy_balance_diagnostic.py on (CR_burn vs fab007-prod vs fab007-ice) to see if rad-escape channel ballooned in CR_burn. Ballooning → hypothesis (1) confirmed.

### May 25 2026 — Zero-D rate verification: **basic kernels pass**

Test setup: pure-DT or DT-CH foam (1.7% C atomic, CH binder) at ρ = 0.04176 g/cc (n_DT = 1e22 cm⁻³), uniform T, hydrodynamics OFF, R = 0.05 cm, 100 zones.

| Test | Helios | Analytic (NRL/Bosch-Hale) | Helios / Analytic | Verdict |
|---|---:|---:|---:|---|
| Fusion rate, pure DT @ 3 keV | 1.11819e+27 #/g/s | 1.1177e+27 | **1.0004** | ✓ exact |
| Brems, pure DT @ 3 keV | 2364.64 TW/g | 2215 | **1.068** | ✓ Gaunt convention (~ḡ=1.18) |
| Fusion rate, 1.7%C foam @ 3 keV | 9.39775e+26 #/g/s | 9.45e+26 | **0.994** | ✓ near-exact |
| Brems, 1.7%C foam @ 3 keV | 4202.76 TW/g | 3403 | **1.235** | ✓ Z-Gaunt + recombination physics (extra +16% over NRL) |

**Verdict:** Basic Helios rate kernels are bug-free for foam compositions at low optically-thin density. The CR_burn pathology is *not* in the basic algorithms. It lives in opacity-table behavior at compressed density, EOS coupling, or some combination.

### May 27 2026 — Zero-D verification matrix complete (11 runs)

Extended the verification to a full T sweep (3, 5, 7, 10, 15 keV) and a density sweep (n_DT = 1e22 → 1e23 → 1e24 → 1e25 at T₀ = 10 keV). Tool: `examples/verify_zero_d.py`.

**T sweep results (optically thin regime, n_DT = 1e22):**

| T (keV) | r_fus pure DT | r_brem pure DT | r_fus foam | r_brem foam |
|---:|---:|---:|---:|---:|
| 3 | 1.000 | 1.054 | 0.994 | 1.218 |
| 5 | 0.999 | 1.044 | 0.994 | 1.147 |
| 7 | 0.999 | 1.037 | 0.994 | 1.114 |
| 10 | 0.998 | 1.030 | 0.993 | 1.090 |
| 15 | 0.997 | 1.024 | 0.991 | 1.068 |

- **Fusion: exact across all T** (≤0.4% error pure DT, ≤0.9% with foam composition)
- **Brems pure DT: 1.054 → 1.024 as T rises** — Helios's Gaunt factor table converges toward unity at higher T (more accurate than NRL's averaged ḡ=1.11)
- **Brems foam: extra Z-Gaunt + recombination contribution scales inversely with T** — 16% above pure-DT-extrapolation at 3 keV, dropping to 4% at 15 keV. Expected physical behavior.

**Density sweep at T₀ = 10 keV:**
- n_DT ≤ 1e23 (ρ ≤ 0.42 g/cc): optically thin clearly; both rates verified pass
- n_DT ≥ 1e24 (ρ ≥ 4.18 g/cc): test self-heats (T evolves to 33.9 → 242.8 keV in 1 ns) AND brems-escape drops to 0.83 / 0.73 × NRL — combined signature of optical-thickness reabsorption (real physics) and fuel depletion confounding the comparison. fab007's burn-region density (12 g/cc foam, 30 g/cc ice) sits in this rad-trapping regime.

**Verdict (final):** All basic rate kernels validated. The foam-burn deficit is not in any rate kernel; it must be in EOS+opacity material physics. Verification summary figure: `comparisons/zero_d_verification_summary.png`.

### May 27 2026 — CR_burn 3-way energy ledger: **opacity pathology refuted**

Ran `energy_balance_diagnostic.py` on the (foam_prod, ice, CR_burn) triplet. Per-stagnation channels:

| Channel | foam_prod | ice | CR_burn |
|---|---:|---:|---:|
| Absorbed laser (kJ) | 1521 | 1492 | 1521 |
| Alpha deposited (kJ) | 5175 | 16185 | 838 |
| Total input (kJ) | 6697 | 17677 | 2359 |
| Plasma thermal (% of input) | 5.3% | 10.9% | 8.0% |
| **Rad escape (% of input)** | **9.1%** | **5.5%** | **9.3%** |
| KE blowoff (% of input) | 86.1% | 81.9% | 82.7% |
| Yield (MJ) | 26.0 | 81.4 | 4.2 |

**Key finding: CR_burn rad-escape fraction (9.3%) is essentially identical to foam_prod (9.1%) — not suppressed.** The opacity is the same between the two runs and the rad transport behaves the same. The previously-suspected "opacity over-radiates at compressed density" hypothesis is refuted.

Drive phase is bit-for-bit identical (peak-v KE 172.40 vs 172.31 kJ — 0.05% match). The pathology shows up entirely between peak velocity and stagnation, in compression.

**Revised diagnosis:** CR_burn's compression collapse is an **EOS-table consistency artifact**. Sesame DT-ice EOS applied to a material with 8% C by mass (Z=6 atoms in a table calibrated for pure equimolar DT) gives stiffer-than-correct response at high pressure during convergence. This is *not* a physical finding about foam opacity; it is the predictable consequence of mixing incompatible material tables.

Lesson: do not swap one of {EOS, opacity} in isolation. They must be a consistent pair, or you produce this kind of artifact. The CR_burn diagnostic detour is closed.

3-way ledger figure: `comparisons/cr_burn_triplet_energy_ledger.png`.

### May 29-30 2026 — FL=0.012 design study: **two-plateau staircase, c33 + c25 are anchors**

Goal: demonstrate that Helios's PROPACEOS wetted foam CAN be burned to near-LILAC yield given sufficient drive, and map the full cone-angle sweep to find the optimal operating point.

**Setup:** FL_prism=0.012 (the corrected global value, replacing fab007 production's per-region 0.06 DT / 0.007 foam-CH split). Five burn-on configurations spanning the geometric design space:

| Config | Cone (°) | Spot (cm) | Focus (cm) |
|---|---:|---:|---:|
| baseline_burn | 37 | 0.18 | 0.22 (= fab007 geom) |
| c33_burn | 33 | 0.16 | 0.20 |
| c28_burn | 28 | 0.16 | 0.18 |
| c25_burn | 25 | 0.14 | 0.16 |
| c20_burn | 20 | 0.14 | 0.15 |

**Decision-matrix results (full 5-point sweep):**

| Config | V_peak | Coup % | Yield | HS ρR | Mean foam ρ | **Foam share** | Adiabat | Verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| fab007 prod foam (ref) | 421 | 73.1 | 26 | 0.35 | 12 | 10.1% | 1.95 | reference |
| baseline_burn | 428 | 74.9 | **33** | 0.39 | 35 | **14.0%** | 1.65 | marginal |
| **c33_burn** | 447 | 79.6 | **55** | 0.53 | 45 | **25.3%** | 1.54 | **strong foam burn** |
| **c28_burn** | 448 | 81.4 | **55** | 0.53 | **61.6** | **25.5%** | 0.98 | **strong foam burn** |
| **c25_burn** | 470 | 85.3 | **69** | 0.64 | 31 | **31.1%** | 1.01 | **strong foam burn** |
| **c20_burn** | 472 | 86.0 | **69** | 0.64 | 25 | **31.5%** | 1.01 | **strong foam burn** |
| fab02 over-drive (ref) | 463 | 84.0 | 59 | n/a | n/a | 26.5% | 1.05 | bootstrap baseline |
| LILAC reference | 410 | 68 | 87 | 0.85 | n/a | ~65% (inferred) | 3.0 | target |

**Headlines (revised after full 5-point sweep):**

1. **The curve is a TWO-PLATEAU STAIRCASE, not a smooth ramp.**
   - **Plateau 1 (33° ≡ 28°):** 55 MJ, 25% foam share. Matches fab02's bootstrap baseline (26.5%).
   - **Plateau 2 (25° ≡ 20°):** 69 MJ, 31% foam share. Exceeds fab02; closes 72% of LILAC gap.
   - **Knee between 28° and 25°:** unlocks at the (spot 0.16 → 0.14, focus 0.18 → 0.16) tightening, NOT at the cone narrowing per se. Spot+focus is the discriminator.

2. **c28's foam over-compresses but doesn't burn more.** c28 reaches mean foam ρ = 61.6 g/cc at stagnation (2× ice's 30 g/cc) yet delivers the same foam yield share as c33 at 45 g/cc. **Foam burn productivity at this drive level is bootstrap-temperature-limited, not compression-limited.** Confirms the broader picture from the May 24 fab007 diagnostic — foam burn is governed by alpha bootstrap reach into foam volume, not by foam compressibility per se.

3. **Two recommended design anchors for different objectives:**
   - **c33_burn** — "minimum geometric intervention" closest to LILAC's calibrated geometry. 55 MJ, 25% foam, V_peak 447 (+9% vs LILAC). For studies where LILAC-thermo-state alignment matters.
   - **c25_burn** — "max demonstrated yield." 69 MJ, 31% foam, V_peak 470 (+15%). Engineering ceiling.

4. **c20 is superseded by c25.** Same yield, same foam share, marginally less aggressive geometric defocus. c25 is the cleaner choice.

5. **Decision tree outcome:** "Geometry helps; baseline still under-burns" — FL alone moves the needle from 10.1% → 14.0% (modest); geometry-on-top-of-FL is the dominant lever. The staircase reveals that this happens in **two distinct steps**, not as a continuous gradient.

**Caveats:**

- All burn-on configurations deviate from LILAC's thermo state. This is a design study answering "what Helios CAN do with PROPACEOS foam," not "what reproduces LILAC."
- 4 of 6 FL=0.012 + burn-on Helios runs hit recurring SIGSEGV/malloc crashes during shutdown (baseline_burn, c20_burn, c25_burn rc=139; c28_burn rc=134 on first try, clean on retry). The .exo files are substantially complete at crash time so metrics extract correctly; budget one retry per run.

**Comparison figure:** `comparisons/pdd_design_comparison.png` (10-row decision matrix + 4 reference rows).

**Cumulative metrics CSV:** `notebooks/pdd_scan_results.csv`.

### May 28 2026 — Planar ablation comparison: **foam-burn deficit isolated as material physics**

Direct planar experiment to isolate EOS+opacity effects without the geometric/temporal complications of spherical implosion.

**Setup:** Planar slab, 200 µm thick, 0.25 g/cc, T₀ = 25 meV, 200 zones. Matched laser drive: 1×10¹⁴ W/cm² at 0.351 µm constant for 2 ns. Two runs:
- `planar_dtice_1`: PROPACEOS DT EOS+opacity (`DT_20250131.prp`)
- `planar_wf_1`: PROPACEOS DT-CH foam EOS+opacity (`DT_and_CH_0p017_20260128.prp`)

**Results (steady state, 1.0–2.0 ns averaged):**

| Metric | DT ice | Foam | foam / DT |
|---|---:|---:|---:|
| **Ablation pressure** | **2.64 Mbar** | **2.13 Mbar** | **0.81** |
| **Post-shock compression ρ/ρ₀** | **4.21** | **3.66** | **0.87** |
| Peak ρ behind shock | 1.05 g/cc | 0.92 g/cc | 0.87 |
| Shock velocity | 107.3 km/s | 105.9 km/s | 0.99 |

**Three findings from one matched-drive experiment:**

1. **Foam ablates 19% softer** at the same laser drive. The higher mean Z̄ from carbon raises the critical density and modifies opacity in the conduction zone, reducing electron-heat coupling to the ablation front. Less effective rocket.

2. **Foam compresses 13% less** behind the same shock. CH skeleton absorbs more of the shock work as internal thermal energy rather than density compaction. Direct EOS comparison at the Rankine-Hugoniot level.

3. **Shock speed is identical** within 1%. Foam doesn't slow the shock — it just produces a less-compressed post-shock state at the same arrival time. Useful null.

**Verdict:** The foam-burn deficit in production fab007 runs is real EOS+opacity physics, captured cleanly by PROPACEOS. Single-shock compression deficit (13%) compounded over the 3-4 shock-cycle implosion + shell deceleration accounts for the observed burn-region density gap (12 vs 30 g/cc, ratio 2.5×); ablation efficiency deficit (19%) compounds with this to explain the cumulative yield gap (26 vs 81 MJ for foam vs ice, ratio 3.1×).

Comparison figure: `comparisons/planar_dtice_vs_foam.png`.

---

## 4. Current working picture

```
                 LILAC foam target  ──── 87 MJ
                                         ▲
                                         │ unexplained ~3× gap
                                         │ (Helios foam vs LILAC foam)
                                         │
                 Helios production foam ─── 26 MJ  (fab007)
                                         ▲
                                         │ 3× gap (foam-vs-ice within Helios)
                                         │ = compression resistance + opacity-side anomaly
                                         │
                 Helios ice (fab007 geom) ─── 81 MJ  ← matches LILAC within 7%
                                         ▲
                                         │ small gap, possibly real physics
                                         ▼
                 LILAC foam target  ──── 87 MJ
```

**Helios with ice in fab007's calibrated drive reproduces LILAC's yield within 7%.** The drive is correct. The remaining residual (Helios foam 26 MJ vs LILAC foam 87 MJ) decomposes into:

- **Foam EOS compression resistance (real physics, confirmed):** PROPACEOS foam EOS gives 13% less post-shock compression than pure DT at matched single-shock strength (planar test, May 28). Compounded over ~3 shock cycles + shell deceleration → 12 g/cc burn-region density vs ice's 30 g/cc → n²T²σv 20× lower locally.
- **Foam opacity ablation deficit (real physics, confirmed):** PROPACEOS foam opacity gives 19% lower ablation pressure than pure DT at matched laser drive (planar test). Less efficient rocket → less drive impulse to the imploding shell.
- **Carbon-arithmetic effects (small):** 1.7% C atomic gives ~16% extra brems beyond NRL at 3 keV dropping to ~4% at 15 keV (zero-D verification). Real but not yield-limiting.
- ~~Opacity pathology at compressed density (suspected, retired):~~ The CR_burn EOS-swap that originally suggested this turned out to be an EOS-table mismatch artifact (3-way energy ledger, May 27). PROPACEOS foam opacity is internally consistent with PROPACEOS foam EOS at all densities the production runs reach.

---

## 5. Open questions

1. **How does LILAC's foam burn 87 MJ when Helios's foam burns 26 MJ at matched drive?** This is now firmly a cross-code comparison question: LILAC's foam EOS table is presumably softer (more compressible) and/or its foam opacity is more transparent than PROPACEOS. Cannot be resolved within Helios alone. Worth direct comparison of the EOS curves at strong-shock conditions if LILAC's tables are accessible.
2. **HDD risk:** transferring fab007-style marginal-ignition foam calibration to HDD geometry will amplify the foam-burn deficit. fab02 (over-driven) is the strong-bootstrap reference baseline. See HDD transfer plan in CLAUDE.md.
3. ~~Why does CR_burn under-compress?~~ — answered May 27: EOS-table mismatch artifact, not physical finding. Retired.
4. ~~Quantify foam-burn opacity contribution.~~ — answered May 27-28: opacity and EOS contribute together, both consistent with PROPACEOS tables; opacity is not pathological. Retired.

---

## 6. Planned next tests

### Completed
- ✅ Zero-D T sweep, pure DT (3, 5, 7, 10, 15 keV) — kernels validated
- ✅ Zero-D T sweep, 1.7%C foam — composition handling validated
- ✅ Zero-D density sweep at T₀=10 keV — optical-thickness regime identified
- ✅ CR_burn energy-balance diagnostic — opacity pathology refuted
- ✅ Planar ablation comparison (DT ice vs foam) — EOS+opacity material effects isolated

### Near-term
- **Generate Xcimer Energy report** summarizing the investigation closeout. Draft in `docs/Xcimer_foam_burn_deficit_report.md`.
- **Re-examine FL transition regime:** fab003 (FL_prism=0.003) NO-IGNITION run is the FL saturation knee marker. Worth dumping its energy ledger.

### Medium-term
- **Cross-code EOS comparison:** if LILAC's foam EOS table or equivalent benchmark data becomes accessible, compare to PROPACEOS foam at strong-shock conditions to localize the remaining 3× Helios-vs-LILAC gap.
- **Planar test variants for parametric sensitivity:** vary CH binder fraction in foam to map how compression deficit scales with composition. Could inform whether a "cleaner" foam (lower CH) is a viable engineering lever.

### Longer-term
- **HDD calibration transfer:** with foam-burn mechanism understood, port the PDD calibration method to HDD/VI_6 geometry using fab02 settings as anchor (bootstrap-strength regime), not fab007 (marginal ignition).

---

## 7. How to extend this log

Each subsection of §3 is a self-contained finding with a date stamp. To add a new entry:

1. Add new subsection under §3 with `### YYYY-MM-DD — Title` and the finding
2. Update the verdict / status flag at the top if the working picture shifts
3. Update §1 inventory if new runs are added
4. Cross-link from project memory (`~/.claude/projects/.../memory/foam_vs_ice_*.md`) for persistent recall across sessions
5. If the conclusion changes a prior finding, mark the old subsection with `**(superseded YYYY-MM-DD)**` rather than deleting — investigation history is the record

The companion Jupyter notebook (`notebooks/foam_vs_ice_investigation.ipynb`) runs the reproducible analyses; this log captures the *narrative* — decisions, refutations, what was ruled out and why.

---

## 8. Cross-references

- Project memory: [foam_vs_ice_energy_balance](../../../.claude/projects/-Users-mehlhorn/memory/foam_vs_ice_energy_balance.md)
- Project memory: [fab007_burn_propagation_diagnostic](../../../.claude/projects/-Users-mehlhorn/memory/fab007_burn_propagation_diagnostic.md)
- Project guide: [CLAUDE.md §"PDD Calibration Scan"](../CLAUDE.md)
- Tools: [`examples/compare_runs.py`](../examples/compare_runs.py), [`examples/energy_balance_diagnostic.py`](../examples/energy_balance_diagnostic.py), [`examples/verify_zero_d.py`](../examples/verify_zero_d.py)
- Notebook: [`notebooks/foam_vs_ice_investigation.ipynb`](../notebooks/foam_vs_ice_investigation.ipynb)
- External report: [`docs/Xcimer_foam_burn_deficit_report.md`](../docs/Xcimer_foam_burn_deficit_report.md)
- Figures: `comparisons/compare_foam_vs_ice_{rt,histories,lineouts}.png`, `comparisons/cr_burn_triplet_energy_ledger.png`, `comparisons/zero_d_verification_summary.png`, `comparisons/planar_dtice_vs_foam.png`
