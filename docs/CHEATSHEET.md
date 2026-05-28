# helios_postprocess — Cheat Sheet

Quick reference for the commands you run most often. Detailed docs are in
[`CLAUDE.md`](../CLAUDE.md), [`README.md`](../README.md), and
[`SETUP_GUIDE.md`](../SETUP_GUIDE.md).

---

## Running Helios (the simulator)

Studio install: `~/Codes/Prism/Helios_11.0.0/`. The binary lives inside the
.app bundle: `Helios.app/Contents/MacOS/Helios`.

**Single batch-mode run:**
```bash
cd ~/Codes/Prism/Helios_11.0.0
./Helios.app/Contents/MacOS/Helios -b -i <input>.rhw -d <run_dir> -o <run_name> -x
```

**Flags:**
| Flag | Meaning |
|---|---|
| `-b` | Batch mode (no GUI) |
| `-i <file>.rhw` | RHW input deck |
| `-d <dir>` | Output directory (created if missing) |
| `-o <name>` | Output basename — produces `<name>.exo`, `.log`, `.rhc`, `.hyw` |
| `-x` | Exit on completion |

**For M-series Macs** (per Prism `README_Mac.txt`), runs can be faster under explicit x86_64 emulation:
```bash
arch -x86_64 ./Helios.app/Contents/MacOS/Helios -b -i ... -d ... -o ... -x
```

**One-time tweak to stop macOS throttling background batch jobs:**
```bash
defaults write NSGlobalDomain NSAppSleepDisabled -bool yes
```

**Example: full Olson PDD run:**
```bash
cd ~/Codes/Prism/Helios_11.0.0
./Helios.app/Contents/MacOS/Helios -b \
  -i ~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_burn.rhw \
  -d ~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab007_foot25_s018_c37_burn \
  -o Olson_PDD_20_fab007_foot25_s018_c37_burn -x
```

---

## Postprocessing — one-run analyses

Run from the repo root with `PYTHONPATH=.` (or `pip install -e .` once and
the import works anywhere). MacBook uses `python`; Studio uses `python3`.

| Task | Command |
|---|---|
| Full PDF report + summary + history CSV | `python examples/run_analysis.py <base_path>` |
| Append per-run row to the run-scan summary | `python examples/scan_summary.py <parent_dir> --out scan_summary.csv` |
| Burn-propagation profile (R / T_e / T_rad / α-dep across interfaces) | `python examples/dump_burn_propagation_profile.py <base_path>` |
| Per-zone burn share + region aggregation | `python examples/dump_per_zone_burn_share.py <base_path>` |
| Burn-rate timing + compression state | `python examples/dump_burn_rate_timing.py <base_path>` |
| Flux-limiter regime check | `python examples/check_flux_limiter.py <base_path>` |

`<base_path>` is the path WITHOUT extension. The script derives the
`.exo`, `.rhw`, `_published.json`, output PDF/CSV/TXT paths from it.

---

## Postprocessing — multi-run comparisons

| Task | Command |
|---|---|
| Two-run side-by-side (R-T, histories, lineouts) | `python examples/compare_runs.py <baseA> <baseB> --labels foam ice` |
| N-way energy ledger (peak-v / stagnation / end-of-run snapshots + PDF) | `python examples/energy_balance_diagnostic.py <base1> <base2> [<base3> ...]` |
| Absorbed-energy comparison | `python examples/compare_absorbed_energy.py <baseA> <baseB>` |
| Published-data comparison (per-run, JSON next to .exo) | runs automatically as part of `run_analysis.py` if `<base>_published.json` exists |

---

## Zero-D verification (basic rate kernels)

| Task | Command |
|---|---|
| Generate RHW input from template at T, ρ | `python examples/make_static_rhw.py --template <tmpl>.rhw --T 5 7 10 15 --outdir .` |
| Generate RHW with density sweep | `python examples/make_static_rhw.py --template <tmpl>.rhw --T 10 --rho 0.418 4.18 41.8 --outdir .` |
| Verify one Helios run vs Bosch-Hale + NRL | `python examples/verify_zero_d.py <base_path> --comp pure-DT --label <name>` |
| Run the full 11-test verification matrix end-to-end (Helios + verify + push CSV) | `bash examples/run_zero_d_matrix.sh --commit` |
| Resume a single failed test | `bash examples/run_zero_d_matrix.sh --rhw <stem> --label <lbl> --comp <preset>` |
| Re-extract only (skip Helios reruns) | `bash examples/run_zero_d_matrix.sh --skip-helios` |

Composition presets: `pure-DT`, `1.7C-CH` (DT + 1.7% C atomic + 1.7% H from CH binder).

---

## Standard paths

### Mac Studio (`tommehlhorn`, run machine, `python3`)
- Helios install: `~/Codes/Prism/Helios_11.0.0/`
- helios_postprocess: `~/helios_postprocess/`
- Olson PDD sims: `~/Sims/Xcimer/Olson_PDD/<run_name>/<run_name>.exo`
- VI_6 / HDD sims: `~/Sims/Xcimer/Xcimer_Sims/D_Montgomery/<run_name>/<run_name>.exo`
- Zero-D verification runs: `~/Sims/Xcimer/Olson_PDD/Fraley_burn/<test>/<test>.exo`

### MacBook (dev machine, `python` via Anaconda)
- helios_postprocess: `~/Codes/helios_postprocess/`
- Local sample sims: `~/Sims/Helios_Sims/Xcimer_Sims/Olson_PDD/Olson_PDD_9/Olson_PDD_9`
- Studio share (read-only via mount): `/Volumes/tommehlhorn/` ⇒ Studio's `~/`

---

## Workflow loop

Edit on MacBook → push → pull on Studio → run.

```bash
# MacBook (after editing code)
git add <files>
git commit -m "..."
git push

# Studio (to pick up code changes)
cd ~/helios_postprocess
git pull

# Studio (after running sims or extracting CSV updates)
git add notebooks/verification_results.csv   # or other output artifacts
git commit -m "..."
git push

# MacBook (to pull data back for the notebook)
cd ~/Codes/helios_postprocess
git pull
```

---

## Notebooks & narrative logs

| File | Purpose |
|---|---|
| [`notebooks/foam_vs_ice_investigation.ipynb`](../notebooks/foam_vs_ice_investigation.ipynb) | Reproducible foam-vs-ice analysis + verification table |
| [`notes/foam_vs_ice_investigation.md`](../notes/foam_vs_ice_investigation.md) | Narrative log: findings, chronology, current picture |
| [`notebooks/verification_results.csv`](../notebooks/verification_results.csv) | Zero-D Helios-vs-analytic results, auto-appended by `verify_zero_d.py` |

Open the notebook with `jupyter lab notebooks/foam_vs_ice_investigation.ipynb`.

---

## Common gotchas

- **`zsh: parse error near '\n'`** — you pasted a multi-line command that began with `<placeholder>` text. Replace placeholders before pasting.
- **`zsh: command not found: #`** — your shell isn't treating `#` as a comment in interactive mode. Strip comment lines from pasted blocks, or wrap in a script file.
- **`ModuleNotFoundError: No module named 'helios_postprocess'`** — invoke with `PYTHONPATH=. python3 examples/...` from the repo root, or `pip install -e .` once.
- **`FileNotFoundError: EXODUS file not found`** — the .exo doesn't exist yet; the run hasn't been simulated. Run Helios first (see top of this doc), then the postprocessor.
- **Helios binary "not found"** — it's nested in the .app bundle. Use `Helios.app/Contents/MacOS/Helios`, not bare `Helios`.

---

## See also

- [`CLAUDE.md`](../CLAUDE.md) — project guide, physics conventions, calibration history
- [`SETUP_GUIDE.md`](../SETUP_GUIDE.md) — install / setup
- [`helios_exodus_variable_reference.md`](../helios_exodus_variable_reference.md) — what each EXODUS variable means + units
