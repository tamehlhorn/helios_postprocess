#!/usr/bin/env python3
"""
dump_burn_region_density.py
===========================

Extract foam-region density and other no-burn implosion metrics at
stagnation for one or more PDD scan runs. Designed for the foam-burn
recovery design study where the headline metric is "did the foam
compress to ice-equivalent burn-region density?"

Per-run extracts (at the stagnation timestep, defined as minimum
hot-spot outer radius):

  - V_peak                       (km/s)
  - peak in-flight ρR             (g/cm²)
  - peak total ρR                 (g/cm²)
  - peak ρ at stagnation          (g/cc)  - over all zones
  - peak ρ in foam region         (g/cc)  - zones 191-320 by default
  - mean ρ in foam region (mass-weighted, at stagnation)
  - peak ρR_foam                  (g/cm²) - cumulative ρR within foam zones
  - mass-averaged adiabat (cold fuel, at peak v) - if available
  - effective coupling %         - from absorbed/delivered

Foam-zone range is configurable (defaults match the Olson_PDD_20 design:
zones 191-320 are the DT-CH foam main fuel layer).

Usage
-----
    python3 dump_burn_region_density.py \\
        ~/Sims/Xcimer/Olson_PDD/PDD_scan/wf_fl012_baseline/wf_fl012_baseline/wf_fl012_baseline \\
        ~/Sims/Xcimer/Olson_PDD/PDD_scan/wf_fl012_c33/wf_fl012_c33/wf_fl012_c33 \\
        --csv ~/helios_postprocess/notebooks/pdd_scan_results.csv

Each base_path is the simulation base WITHOUT extension. Appends one
row per run to the CSV (deduplicates by run label).
"""
from __future__ import annotations
import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np


def extract_metrics(base_path: Path,
                    foam_zone_lo: int = 191,
                    foam_zone_hi: int = 320,
                    verbose: bool = False) -> dict:
    """Load a run and extract no-burn implosion metrics at stagnation."""
    from helios_postprocess import HeliosRun
    from helios_postprocess.data_builder import build_run_data
    from helios_postprocess.icf_analysis import ICFAnalyzer

    exo_path = base_path.with_suffix('.exo')
    rhw_path = base_path.with_suffix('.rhw')
    if not exo_path.exists():
        raise FileNotFoundError(f"EXODUS file not found: {exo_path}")

    # Optional RHW for drive/coupling diagnostics
    rhw_config = drive_T = drive_t = None
    if rhw_path.exists():
        try:
            from helios_postprocess.rhw_parser import load_rhw_configuration
            rhw_config = load_rhw_configuration(rhw_path)
            drive_T = rhw_config.drive_temperature
            drive_t = rhw_config.drive_time
        except Exception:
            pass

    run = HeliosRun(str(exo_path), verbose=verbose)
    data = build_run_data(run, time_unit='s',
                          rhw_config=rhw_config,
                          drive_temperature=drive_T,
                          drive_time=drive_t)
    run.close()

    # Light analysis -- just enough to identify stagnation and ablation front
    analyzer = ICFAnalyzer(data)
    analyzer.analyze_drive_phase()
    analyzer.analyze_stagnation_phase()
    analyzer.compute_performance_metrics()

    time_ns = data.time
    rho     = data.mass_density                  # (n_t, n_z) g/cc
    zbnd    = data.zone_boundaries               # (n_t, n_z+1) cm
    zmass   = data.zone_mass                     # (n_t, n_z) g
    n_zones = rho.shape[1]

    # ---- Stagnation timestep ----
    stag_idx = None
    if getattr(data, 'stag_time', 0.0) > 0:
        stag_idx = int(np.argmin(np.abs(time_ns - data.stag_time)))
    else:
        # Fallback: minimum hot-spot outer radius across the run
        if data.region_interfaces_indices is not None and data.region_interfaces_indices.shape[1] >= 1:
            ri = data.region_interfaces_indices
            hs_r = np.array([zbnd[t, int(ri[t, 0])] for t in range(len(time_ns))])
            stag_idx = int(np.argmin(hs_r))
        else:
            stag_idx = int(np.argmax(np.max(rho, axis=1)))

    t_stag_ns = float(time_ns[stag_idx])

    # ---- Foam-region selection ----
    foam_lo = max(0, foam_zone_lo)
    foam_hi = min(n_zones, foam_zone_hi + 1)
    foam_slice = slice(foam_lo, foam_hi)

    # ---- Per-stag metrics ----
    rho_stag      = rho[stag_idx]
    rho_peak_all  = float(np.max(rho_stag))
    rho_peak_foam = float(np.max(rho_stag[foam_slice]))
    foam_mass     = zmass[stag_idx, foam_slice]
    if np.sum(foam_mass) > 0:
        rho_mean_foam = float(np.sum(rho_stag[foam_slice] * foam_mass) / np.sum(foam_mass))
    else:
        rho_mean_foam = float('nan')

    # Total foam mass at t=0 vs at stagnation (Lagrangian -- should be conserved)
    foam_mass_total_g = float(np.sum(zmass[stag_idx, foam_slice]))

    # ---- Cumulative areal density in foam zones at stagnation ----
    # ρR_foam = integral from foam_lo node to foam_hi node of rho dr
    zwidths = np.diff(zbnd[stag_idx])
    rhoR_foam = float(np.sum(rho_stag[foam_slice] * zwidths[foam_slice]))

    # ---- Peak total ρR over time ----
    if data.areal_density_vs_time is not None:
        peak_total_rhoR = float(np.nanmax(data.areal_density_vs_time[:, -1]))
    elif data.areal_density is not None:
        peak_total_rhoR = float(np.nanmax(data.areal_density))
    else:
        peak_total_rhoR = float('nan')

    # ---- V_peak (in-flight, from ICFAnalyzer) ----
    v_peak_kms = abs(float(getattr(data, 'peak_implosion_velocity', float('nan'))))

    # ---- Adiabat (mass-averaged ice, at peak v) ----
    adiabat = float(getattr(data, 'adiabat_mass_averaged_ice', float('nan')))

    # ---- Effective coupling ----
    coupling_pct = float(getattr(data, 'eff_avg_coupling_pct', float('nan')))

    # ---- Peak in-flight ρR (pre-stagnation max of total ρR) ----
    if data.areal_density_vs_time is not None and getattr(data, 'peak_velocity_index', None):
        pv_idx = int(data.peak_velocity_index)
        cum = data.areal_density_vs_time[:, -1]
        peak_inflight_rhoR = float(np.nanmax(cum[:pv_idx + 1]))
    else:
        peak_inflight_rhoR = float('nan')

    # ---- Burn-on extras (auto-detected; NaN if burn was OFF) ----
    yield_MJ        = float('nan')
    bang_time_ns    = float('nan')
    HS_rhoR_max     = float('nan')
    foam_yield_pct  = float('nan')
    ignition_flag   = ''

    bt = getattr(data, 'bang_time', 0.0) or 0.0
    if bt > 0:
        bang_time_ns = float(bt)
    if hasattr(data, 'hot_spot_rhoR_vs_time') and data.hot_spot_rhoR_vs_time is not None:
        v = data.hot_spot_rhoR_vs_time
        if np.any(~np.isnan(v)):
            HS_rhoR_max = float(np.nanmax(v))
            ignition_flag = 'YES' if HS_rhoR_max >= 0.3 else 'no'

    # Total fusion yield -- prefer Helios's energy_output if populated, else from neutron count
    energy_output = getattr(data, 'energy_output', 0.0) or 0.0
    if energy_output > 0:
        yield_MJ = float(energy_output)
    elif hasattr(data, 'dt_neutron_count') and data.dt_neutron_count is not None:
        nc = np.asarray(data.dt_neutron_count)
        if nc.ndim >= 2:
            n_total = nc[-1].sum()
        else:
            n_total = nc[-1]
        if n_total > 0:
            yield_MJ = float(n_total * 17.6e6 * 1.602176634e-19 * 1e-6)

    # Per-zone fusion product (DT alphas) for foam yield share
    try:
        if hasattr(data, 'dt_neutron_count_per_zone') and data.dt_neutron_count_per_zone is not None:
            # If pipeline already exposes per-zone cumulative DT neutron count
            nz = np.asarray(data.dt_neutron_count_per_zone)[-1]
            if nz.sum() > 0:
                foam_yield_pct = 100.0 * float(nz[foam_slice].sum() / nz.sum())
    except Exception:
        pass

    def _r(x, n):
        return round(float(x), n) if not (x is None or (isinstance(x, float) and np.isnan(x))) else ''

    return dict(
        label                 = base_path.name,
        run_path              = str(base_path),
        timestamp             = datetime.now().isoformat(timespec='seconds'),
        t_stag_ns             = _r(t_stag_ns, 4),
        bang_time_ns          = _r(bang_time_ns, 4),
        V_peak_kms            = _r(v_peak_kms, 1),
        peak_inflight_rhoR    = _r(peak_inflight_rhoR, 4),
        peak_total_rhoR       = _r(peak_total_rhoR, 4),
        rho_peak_all_gcc      = _r(rho_peak_all, 2),
        rho_peak_foam_gcc     = _r(rho_peak_foam, 2),
        rho_mean_foam_gcc     = _r(rho_mean_foam, 2),
        rhoR_foam             = _r(rhoR_foam, 4),
        foam_mass_total_mg    = _r(foam_mass_total_g * 1e3, 4),
        adiabat               = _r(adiabat, 2),
        coupling_pct          = _r(coupling_pct, 1),
        HS_rhoR_max           = _r(HS_rhoR_max, 4),
        yield_MJ              = _r(yield_MJ, 3),
        foam_yield_pct        = _r(foam_yield_pct, 1),
        ignition              = ignition_flag,
        foam_zone_range       = f'{foam_lo}-{foam_hi - 1}',
    )


CSV_COLUMNS = [
    'timestamp', 'label', 'run_path',
    't_stag_ns', 'bang_time_ns', 'V_peak_kms',
    'peak_inflight_rhoR', 'peak_total_rhoR',
    'rho_peak_all_gcc', 'rho_peak_foam_gcc', 'rho_mean_foam_gcc',
    'rhoR_foam', 'foam_mass_total_mg', 'adiabat', 'coupling_pct',
    'HS_rhoR_max', 'yield_MJ', 'foam_yield_pct', 'ignition',
    'foam_zone_range',
]


def append_csv(csv_path: Path, row: dict) -> None:
    """Append/update a row in the CSV, keyed by label."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    existing = []
    if csv_path.exists():
        with csv_path.open() as f:
            reader = csv.DictReader(f)
            for r in reader:
                if r['label'] != row['label']:
                    existing.append(r)
    existing.append({c: row.get(c, '') for c in CSV_COLUMNS})
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for r in existing:
            writer.writerow(r)


def print_row(row: dict) -> None:
    print(f"\n  {row['label']}:")
    print(f"    t_stag = {row['t_stag_ns']} ns" + (f"    bang = {row['bang_time_ns']} ns" if row.get('bang_time_ns') not in ('', None) else ''))
    print(f"    V_peak                      = {row['V_peak_kms']:>8} km/s")
    print(f"    peak in-flight rhoR (g/cm2) = {row['peak_inflight_rhoR']:>8}")
    print(f"    peak total rhoR (g/cm2)     = {row['peak_total_rhoR']:>8}")
    print(f"    peak rho (all, g/cc)        = {row['rho_peak_all_gcc']:>8}")
    print(f"    PEAK rho in FOAM (g/cc)     = {row['rho_peak_foam_gcc']:>8}     <-- design metric")
    print(f"    mean rho in FOAM (g/cc)     = {row['rho_mean_foam_gcc']:>8}     <-- design metric")
    print(f"    foam rhoR contribution      = {row['rhoR_foam']:>8} g/cm2")
    print(f"    foam mass (total)           = {row['foam_mass_total_mg']:>8} mg")
    print(f"    adiabat (mass-avg ice)      = {row['adiabat']:>8}")
    print(f"    effective coupling %        = {row['coupling_pct']:>8}")
    # Burn-on extras (only if present)
    if row.get('yield_MJ') not in ('', None):
        print(f"    --- burn ON ---")
        print(f"    Total yield (MJ)            = {row['yield_MJ']:>8}     <-- headline")
        print(f"    HS rhoR peak (g/cm2)        = {row['HS_rhoR_max']:>8}")
        print(f"    Foam yield share (%)        = {row['foam_yield_pct']:>8}")
        print(f"    Ignition                    = {row['ignition']:>8}")


def main():
    ap = argparse.ArgumentParser(description='Extract foam-region density metrics from PDD scan runs.')
    ap.add_argument('base_paths', nargs='+',
                    help='One or more run base paths (without extension)')
    ap.add_argument('--csv', type=Path,
                    help='Append results to this CSV (default: print only)')
    ap.add_argument('--foam-lo', type=int, default=191,
                    help='Lowest foam zone index (default: 191 for Olson_PDD_20)')
    ap.add_argument('--foam-hi', type=int, default=320,
                    help='Highest foam zone index inclusive (default: 320)')
    ap.add_argument('-v', '--verbose', action='store_true')
    args = ap.parse_args()

    rows = []
    for bp_str in args.base_paths:
        bp = Path(bp_str).expanduser().resolve()
        try:
            row = extract_metrics(bp, foam_zone_lo=args.foam_lo,
                                  foam_zone_hi=args.foam_hi, verbose=args.verbose)
        except FileNotFoundError as e:
            print(f"\n  SKIP {bp.name}: {e}")
            continue
        except Exception as e:
            print(f"\n  ERROR on {bp.name}: {e}")
            continue
        rows.append(row)
        print_row(row)

    if args.csv and rows:
        for row in rows:
            append_csv(args.csv, row)
        print(f"\nWrote/updated {len(rows)} row(s) to {args.csv}")


if __name__ == '__main__':
    sys.exit(main())
