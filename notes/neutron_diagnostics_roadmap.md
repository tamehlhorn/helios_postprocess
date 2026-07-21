# Feature roadmap from Gopalaswamy et al., PoP 32, 092702 (2025)

"Interpreting scattered neutron spectra on OMEGA cryogenic implosions with
measured backgrounds and higher-order scattering." (IRIS Monte-Carlo
postprocessor; LLE/Imperial/LANL/MIT.) Read against your two goals:

- **A. Correct the *emergent* (within-target) neutron spectrum for diagnostics** — TOF included.
- **B. Provide a fuller physics *source* for downstream transport** — tritium breeding, first-wall damage.

The paper is, in effect, the physics spec for both. Its central finding is that
the scattered spectrum is a rich, model-dependent object: an **analytic
single-scatter model with ice-block attenuation reproduces IRIS second-order to
within ~10% at OMEGA ρR (~150 mg/cm²)** — which is exactly the fidelity a
homegrown module can reach, and validates the "keep homegrown, align" choice.
It also draws the line where you *must* hand off to Monte Carlo (IRIS/NeSST):
high ρR, higher-than-second-order scatter, and angle-resolved images.

## The emergent-spectrum stack (what "within-target physics" actually means)

The measured spectrum (Fig. 2) = primaries + all scattering channels:

| Layer | Physics | Key relation | Analytic? |
|---|---|---|---|
| Primaries | DT (14.03 MeV), DD (2.45 MeV), **TT continuum** | Ballabio widths+shift; TT shape from ENDF / Appelbe | yes |
| Flow shift | bulk-fluid velocity shifts the primary peak | per-LOS mean-energy shift | yes |
| Elastic down-scatter | nD, nT off compressed fuel | E(θ)=E₀[(cosθ+√(A²−sin²θ))/(A+1)]² (Eq. B5) | yes (single-scatter) |
| Backscatter edge | nT edge at **3.5 MeV** (=E₀·[(A−1)/(A+1)]², A≈3) → ρR | edge amplitude ∝ ρR | yes |
| Edge broadening | shell **temperature** broadens the edge | σ=2√2·((A−1)√A/(A+1)²)·√(E₀⟨T_s⟩) (Eq. 9) | yes |
| Inelastic (n,2n) | D(n,2n), T(n,2n) — **neutron multiplication** | ENDF cross-sections | yes (rates); MC (spectrum) |
| Attenuation | differential forward/back through ρL | ice-block L(r,θ)=√(R²−r²sin²θ)−r·cosθ (Eq. B1) | yes |
| Higher-order scatter | 2nd order suppresses edge ~10–20% at 150 mg/cm² | fraction ~(σρR/m_i)ⁿ; ρR bias ≈ e^(−ρR/14 g·cm⁻²) | **MC only** |

## Where our module is vs where this points

We have: DT (now 14.06) + DD birth peaks, Brysk widths, emission weighting,
Kyle-aligned neutron-weighted `<ρ(r)>`/`<T_i(r)>` profiles, and a **placeholder**
DSR = ρR/20.4. The paper shows the placeholder is exactly the thing to replace:
real ρR inference comes from a **scattered-spectrum model**, not a fixed
calibration constant, and the observable is subtly *not* ρR (see ρL below).

## Recommended features — Goal A (diagnostics / TOF)

Ordered by value, all analytic (single-scatter regime = OMEGA/PDD-class ρR):

1. **TOF forward model** (Eqs. 2–6). τ=ct/d, E(τ)=m_nc²(1/√(1−(1/τ)²)−1),
   dN/dτ=dN/dE·dE/dτ, then `dN/dτ_sim = IRF ⊗ (dN/dE·dE/dτ) + S_bg`. Small,
   unambiguous, testable; turns any synthetic spectrum into a synthetic nTOF
   trace directly comparable to Forrest/Knauer P7/H10 records. **Do first.**
2. **TT continuum + DD primary spectral shapes** (ENDF/Appelbe). Completes the
   3–10 MeV region that the TT primary and backgrounds dominate (Eq. 7).
3. **Single-scatter scattered spectrum**: fold each primary through nD/nT elastic
   σ(E) and the Eq. B5 kinematics, weighted by the neutron-weighted `<ρ(r)>`
   profile we already produce → the down-scatter shoulder + **nT backscatter
   edge**. This replaces the placeholder DSR with a physical scattered spectrum;
   ρR then comes from the edge amplitude and the 10–12 MeV DSR band.
4. **Ice-block attenuation** (Appendix B): differential forward/back attenuation
   (back attenuated ~6× more; effective ρR +50% back / −60% forward, Fig. 15).
   Needed to get the edge amplitude — hence ρR — right.
5. **Backscatter-edge thermal broadening** (Eq. 9): the edge width reports shell
   `⟨T_s⟩` — a genuinely new diagnostic beyond ρR.
6. **ρL vs ρR bookkeeping** (Appendix D): a scattered-neutron detector measures a
   neutron-path-length-averaged ρL, not geometric ρR (1D: ρR≈(0.97±0.02)ρL).
   Report both, with the two-region (η,ξ) mapping, so our numbers mean the same
   thing Will's/LLE's do.
7. **Backgrounds + detector response** (Eqs. 5,8; NLO, gating): only needed if we
   want to compare to *raw* scintillator traces rather than unfolded spectra.

## Recommended features — Goal B (transport source: breeding, first-wall)

The transport source is the **same emergent spectrum, but exported full
energy–angle–time-resolved over 4π**, not just along diagnostic LOS. Additions:

1. **Full-4π emergent source export**: the primaries + within-target-scattered
   spectrum as `S(E, Ω, t)` at the target boundary — the source term a chamber
   /blanket transport code (breeding, dpa) ingests. Our neutron-weighted profiles
   + single-scatter engine already produce the ingredients; this is the export
   format (extend the `.npz`, and/or an MCNP/OpenMC SDEF-style source writer).
2. **(n,2n) neutron multiplication** (D(n,2n), T(n,2n)): matters for total
   neutron economy into a blanket — track the multiplied population and its soft
   spectrum, not just DT/DD/TT births. Channels already loaded in `data_builder`.
3. **Time resolution**: keep the burn-history time axis on the source (pulsed
   first-wall loading / activation timing), which our `burn_history` already has.
4. **Spectrum shape to threshold reactions**: breeding cares about both the 14 MeV
   peak (⁷Li(n,n'α)T, threshold ~2.5 MeV) and the down-scattered tail (⁶Li(n,α)T);
   first-wall dpa cares about the full fluence spectrum. So the *within-target
   down-scatter we compute for Goal A is also what shifts the breeding/damage
   spectrum* — the two goals share one engine.

**Caveat for Goal B specifically:** IFE-class fuel ρR (≫ OMEGA's 150 mg/cm²) puts
you out of the single-scatter regime — the emergent source is already heavily
multiply-scattered before it leaves the target. There, our analytic engine gives
a first-order source; the honest path is to either (a) add an approximate
second-order correction, or (b) hand the neutron-weighted profiles to IRIS/NeSST
and let MC produce the source. Flag this explicitly wherever a high-ρR source is
exported.

## The analytic / Monte-Carlo boundary (from the paper's own validation)

- **Analytic single-scatter + ice-block ≈ IRIS 2nd-order above the edge** (Fig. 5)
  → our homegrown engine is legitimate for OMEGA/PDD-class ρR diagnostics.
- **Below the nT edge (≤3.5–5 MeV), and at high ρR, you need MC** — rescattered
  neutrons that the analytic method can't transport (it doesn't conserve flux).
  Using a 1st-order model on a真 2nd-order spectrum under-infers ρR by
  ≈ e^(−ρR/14 g·cm⁻²) (~10% at 150 mg/cm², worse at IFE ρR).
- **NeSST is the drop-in for the single-scatter scattered spectrum** (same regime,
  community-validated, DT/DD/TT + elastic/(n,2n)); **IRIS is the full-MC source**
  for high ρR / angle-resolved / >2nd order. This paper is the bridge that tells
  us our homegrown engine and NeSST should agree, and where IRIS is required.

## Suggested sequencing

1. TOF forward model (Eqs. 2–6) — small, high-value, unblocks synthetic nTOF.
2. Single-scatter scattered spectrum (Eq. B5 + nD/nT σ(E)) on the neutron-weighted
   `<ρ(r)>` profile → real DSR + nT edge; retire the placeholder DSR.
3. Ice-block attenuation (Appendix B) + edge thermal broadening (Eq. 9).
4. ρL/ρR bookkeeping (Appendix D) + TT/DD spectral shapes.
5. 4π source export (+ (n,2n) multiplication) for the transport hand-off.
6. Cross-validate the whole engine against NeSST on one run; document where IRIS
   is required (high-ρR IFE source, >2nd-order, images).

Every step reuses the Kyle-aligned neutron-weighted profiles we just built, and
each is independently testable against the paper's figures (edge at 3.5 MeV,
ρR-bias ≈ e^(−ρR/14), forward/back attenuation ratio ~6×).

## Supplement — Forrest et al., RSI 83, 10D919 (2012)

The foundational OMEGA nTOF scattered-spectrum paper (ref. 27 of Gopalaswamy;
same IRIS-on-LILAC method). It **confirms** the stack above and adds a few
concrete, implementable specifics:

**Grounds our DSR in cross-section physics (biggest new item).** The
down-scatter yield is a cross-section-weighted line integral over the fuel
ion-density profile,
`Y_n'/Y_n = ∫∫ (σ_d·n_d + σ_t·n_t) dr dt` (Eq. 4), and areal density follows as
`⟨ρR⟩ = 5·(m_p/(σ_d+σ_t))·(Y_n'/Y_n)` mg/cm² (Eq. 5). **This is where the
"~20 g/cm² per DSR" constant comes from** — it's `5·m_p/(σ_d+σ_t)` at the
down-scatter energy. Two upgrades for us: (a) replace the hardcoded 20.4 with the
ENDF-derived `σ_d+σ_t` value; (b) better, evaluate Eq. 4 **directly on our
Kyle-aligned neutron-weighted `<n_d(r)>`/`<n_t(r)>` profiles** — turning the
placeholder DSR into a physics DSR built from the profiles we already produce.
Still fully analytic. This retires the "circular DSR" objection: it's no longer
ρR relabeled, it's the cross-section-weighted profile integral.

**Explicit kinematic edges (test anchors):** nT edge **3.53 MeV**, nD edge
**1.56 MeV**; scattered fraction `E_n'/E_n = 1 − 4A/(1+A)²·cos²θ` (Eq. 6, same
kinematics as Gopalaswamy B5). Elastic ranges: nD 1.56–14.06, nT 3.53–14.06 MeV.

**Practical TT shape:** below 6 MeV the TT spectrum is well approximated by a
decaying exponential `A·e^(−t/τ)`, amplitude set by reactivity (∝ Ti⁴) — a
lighter first cut than the full ENDF/Appelbe shape.

**Extra inelastic channel:** deuteron breakup **D(n,2n)p** matters below 2 MeV —
a neutron-multiplying channel populating the low-energy tail (relevant to the
Goal-B breeding/first-wall source).

**Detector realism (details the TOF-forward-model feature):** gated PMT to
suppress the DT primary; oxygenated-xylene low-afterglow scintillator (×10⁵ at
100 ns); mid-beam collimator + gold-lined housing; MCNP-modeled chamber-scatter
background; energy-dependent IRF; light output ∝ neutron energy with a
finite-thickness interaction-probability correction (~5%); ~20% calibration
drift from scintillator oxygen depletion. Error budget **~15%** total (DT yield
5%, light-output model 5%, statistics ~2%, nD/nT σ at 180° ~10%). These are the
knobs — and the uncertainty budget — of the synthetic-nTOF forward model.

**Cross-check geometry:** MRS measures *forward*-scatter (10–12 MeV recoil
protons/deuterons); nTOF measures *back*-scatter (1–6 MeV) — different LOS, the
asymmetry/coverage point again. `χ_1D = ρR^0.8·(Ti/4.4 keV)^1.8` (Eq. 1) is the
Ti-based Lawson form to report alongside ρR.

Net: **no change to the roadmap** — this paper sharpens step 2 (the single-scatter
DSR is the cross-section-weighted integral over our profiles, Eqs. 4–5, not a
fixed constant) and step 1 (concrete detector-response terms), and adds D(n,2n)p
to the Goal-B channel list.

Sources: Gopalaswamy et al., *Phys. Plasmas* **32**, 092702 (2025),
doi:10.1063/5.0288729; Forrest et al., *Rev. Sci. Instrum.* **83**, 10D919
(2012), doi:10.1063/1.4742926.
