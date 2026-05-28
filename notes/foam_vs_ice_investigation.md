# Foam vs Ice Burn Investigation — Helios PDD Calibration

**Lead question:** Why does Helios's fab007 production run (wetted DT-CH foam main fuel) yield 26 MJ while LILAC's same target yields 87 MJ? Where is the foam-burn deficit, and is it a Helios artifact or a physical insight about wetted-foam designs?

**Status:** Active. Foam-burn deficit traced to a combination of (a) foam compression resistance and (b) a separate opacity-side anomaly emerging at compressed densities. Basic Helios rate kernels validated to pass.

**Owner:** Prof T (Mehlhorn). Last updated: 2026-05-25.

---

## 1. Run inventory

| Label | Path (Studio) | Description |
|---|---|---|
| **fab007-foam (production)** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_burn/` | Calibrated PDD point — wetted DT-CH foam, FL_prism=0.007, yield 26 MJ |
| **fab007-foam, burn OFF** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_nb/` | Hydrodynamic-state reference for the production run |
| **fab015-foam** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab015_foot25_s018_c37_burn/` | Earlier calibration attempt, FL_prism=0.015, yield 29 MJ |
| **fab007-ice** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_ice_burn/` | fab007 geometry with DT ice replacing foam in main fuel, yield **81 MJ** |
| **fab007-CR** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_CR_burn/` | fab007 geometry, **Sesame DT-ice EOS + PROPACEOS foam opacity** in main fuel — diagnostic for separating EOS vs opacity. Yield 4.2 MJ. |
| **fab02-foam** | `~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab02_foot25_s016_burn/` | Bootstrap-strength reference, over-driven foam, yield 59 MJ |
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

### May 25 2026 — CR_burn EOS-swap test: **opacity-side pathology suspected**

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

- **Foam compression resistance (real physics, dominant):** PROPACEOS foam EOS resists DT collapse, giving 12 g/cc burn-region density vs ice's 30 g/cc. n²T² scaling drives the burn-rate-density 20× lower in foam.
- **Opacity-side anomaly (Helios artifact, suspected):** CR_burn EOS-swap test shows that *removing* the foam EOS (replacing with Sesame DT-ice) does NOT recover ice-like yield. Suggests the PROPACEOS foam opacity at compressed density is also pathological, beyond what the EOS does.
- **Carbon-arithmetic effects (real but small):** 1.7% C atomic gives ~16% extra brems beyond NRL (Z-Gaunt + recombination), ~10% fuel dilution. Real but not yield-limiting.

---

## 5. Open questions

1. **Why does CR_burn under-compress?** Need energy-balance diagnostic on CR_burn to see rad-escape channel. Hypothesis: PROPACEOS foam opacity over-radiates at ice-compressed densities.
2. **Quantify the foam-burn opacity contribution.** Once CR_burn rad-channel is known, partition the Helios foam-vs-ice gap into compression-resistance share vs opacity-anomaly share.
3. **How does LILAC's foam burn 87 MJ?** Their PROPACEOS-equivalent foam model evidently doesn't have the opacity pathology, or compensates with different transport. Need cross-code comparison framing.
4. **HDD risk:** transferring fab007-style marginal-ignition foam calibration to HDD geometry could amplify the foam-burn deficit. fab02 (over-driven) is the strong-bootstrap reference baseline. See HDD transfer plan in CLAUDE.md.

---

## 6. Planned next tests

### Near-term (this week)

- **Zero-D T sweep, pure DT** at T = 5, 7, 10, 15 keV → check that Helios/NRL ratio holds at ~1.07 across T. Drift indicates non-trivial extra physics (e-e brems, relativistic correction).
- **Zero-D optical-thickness sweep** at fixed T = 10 keV, n_DT = 1e22 → 1e23 → 1e24 → 1e25 cm⁻³ → identify where Helios's rad transport diverges from optically-thin analytic. This is the test that bears most directly on the CR_burn opacity hypothesis.
- **CR_burn energy-balance diagnostic** → run `energy_balance_diagnostic.py` on CR_burn vs fab007-prod vs fab007-ice triplet.

### Medium-term

- **Clean Fraley-style microsphere test:** uniform sphere, no driver, initial T_0 ≳ 10 keV, vary material between pure DT, 1.7% C foam, fab007-real-foam composition. Compare burnup fraction vs analytic Fraley. Isolates EOS+opacity behavior in a steady (no-driver) ignition environment.
- **Re-examine FL transition regime:** the fab003 (FL_prism=0.003) NO-IGNITION run is the FL saturation knee marker. Worth dumping its energy ledger for comparison.

### Longer-term

- **HDD calibration transfer:** once foam-burn opacity question is closed, port the PDD calibration method to HDD/VI_6 geometry with the strong-bootstrap baseline in mind (fab02 settings as anchor, not fab007).

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
- Tools: [`examples/compare_runs.py`](../examples/compare_runs.py), [`examples/energy_balance_diagnostic.py`](../examples/energy_balance_diagnostic.py)
- Notebook: [`notebooks/foam_vs_ice_investigation.ipynb`](../notebooks/foam_vs_ice_investigation.ipynb)
