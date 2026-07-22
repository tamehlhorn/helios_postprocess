# Nuclear data — n+D / n+T elastic (for down-scatter modeling)

`endf_nD_nT_elastic.npz` — elastic scattering cross sections and angular
distributions used by the down-scatter model (Step 3). Self-contained: loaded
with numpy alone, no ENDF parser or NeSST needed at runtime.

## Provenance

- **Library:** ENDF/B-VIII.0 (n+¹H₂ MAT 128, n+¹H₃ MAT 131).
- **Evaluation:** Young / Hale / Chadwick R-matrix (CH97/CH99, DIST-FEB18).
- **Source of the raw files:** bundled with NeSST (Crilly, arXiv:2605.20432),
  `n-001_H_002.endf`, `n-001_H_003.endf`.
- **Extracted with:** the `endf` package (ENDF-6 parser) — MF3/MT2 cross section,
  MF4/MT2 angular distribution (Legendre, CM frame LCT=2).

### Version note (VII.1 vs VIII.0)
The Gopalaswamy/IRIS work cites **ENDF/B-VII.1** (Chadwick et al., Nucl. Data
Sheets 112, 2887, 2011). NeSST ships **VIII.0**. For the n+D and n+T *elastic*
channels these are the *same* underlying Hale R-matrix evaluation, so the
difference at ICF energies is sub-percent and does not affect the down-scatter.
If an exact byte-match to an IRIS run is ever required, the VII.1 files can be
substituted; the table format here is version-agnostic.

## Contents (npz keys)

| key | shape | meaning |
|---|---|---|
| `nD_E_eV`, `nD_sigma_barn` | (n,) | elastic σ(E), n+D |
| `nD_AWR`, `nD_LCT` | scalar | mass ratio (m_target/m_n), frame (2=CM) |
| `nD_ang_E_eV` | (341,) | energy grid for angular data |
| `nD_ang_al` | (341, 11) | Legendre coeffs a_l(E), l=1..; a_0≡1 |
| `nT_*` | — | same set for n+T (130 E-pts, max l=8) |
| `_library`, `_source`, `_reaction`, `_note` | str | provenance |

## Key values

| | σ_el @14.06 MeV | σ_el @3.53 MeV | AWR | backscatter E′=((A−1)/(A+1))²·14.06 |
|---|---|---|---|---|
| **n+D** | 0.641 b | 1.964 b | 1.9968 | 1.56 MeV (nD edge) |
| **n+T** | 0.936 b | 2.447 b | 2.9896 | 3.50 MeV (nT edge) |

n+T elastic rises ~2.6× from the 14 MeV DT birth energy to its 3.5 MeV
backscatter edge — the amplitude that encodes ρR (Gopalaswamy Fig. 4).
