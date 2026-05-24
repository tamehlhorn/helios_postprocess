"""
dump_burn_rate_timing.py
========================

Per-zone burn-rate timing and compression-state diagnostic, aggregated
by region.  Discriminates between mechanism (1)/(3) "foam compression
deficit / inertia" and mechanism (2) "burn duration too short to reach
foam bulk" for the PDD priority 0a investigation.

Motivation (May 24 2026)
------------------------
`dump_per_zone_burn_share.py` showed that fab007 burns DT_ice at
36,786 kJ/mg-DT but DT-CH_foam at only 642 kJ/mg-DT — a 57× specific-
yield gap.  The earlier `dump_burn_propagation_profile.py` showed the
ice/foam interface is continuous in T_e, T_rad, α-deposition.  So the
arrest is BODY-LEVEL not INTERFACE-LEVEL: the bulk foam DT mass never
reaches burn conditions.

Three remaining candidate mechanisms:

  (1) Foam compression deficit.  ρ_foam ≪ ρ_ice at burn time; n²σv
      collapses.  Diagnostic: foam ρ_max/ρ_0 (compression ratio at peak
      compression) compared to ice ρ_max/ρ_0.  An EOS problem with
      PROPACEOS for DT-CH foam would show up as anomalously low foam
      compression ratio.
  (2) Burn duration too short.  Foam zones DO get to high T_ion, but
      only after ice has burnt out and disassembly starts.  Diagnostic:
      timeline of first-zone-in-region crossing 4.5 keV vs bang time.
  (3) Foam mass inertia / heat capacity.  Same heating power, 7× more
      mass to warm; foam doesn't reach burn T in the available time.
      Diagnostic: per-zone t_peak_burn vs t_cross_4.5keV — large lag in
      foam zones implies they were heating slowly.

(1) and (3) are physically nearly the same; (2) is independently
distinguishable.

Usage
-----
    python3 dump_burn_rate_timing.py <base_path> [options]

Produces alongside the .exo:

    <base>_burn_rate_timing.csv     -- per-zone timing + compression
    <base>_burn_rate_timing.txt     -- console summary mirrored to disk

Options
-------
    --t-bootstrap-keV FLOAT
        T_ion threshold for "alpha bootstrap catches" in keV.  Default
        4.5, matching the Olson burn-propagation convention used
        elsewhere in this codebase.
    --rate-frac FLOAT
        Fraction of peak global rate at which a zone is considered to
        be "burning meaningfully" for the rise/fall timing.  Default 0.1.
    --out-base PATH
    --verbose
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from helios_postprocess import HeliosRun


# EXODUS variable name resolution
ZBND_NAMES        = ['zone_boundaries', 'coord', 'coordx']
ZBAR_NAMES        = ['mean_charge', 'Zbar']
RHO_NAMES         = ['mass_density', 'dens', 'density']
TION_NAMES        = ['ion_temperature', 'temp', 'temperature']
PI_NAMES          = ['ion_pressure', 'pres', 'pressure']
PE_NAMES          = ['elec_pressure', 'electron_pressure']
ZONE_MASS_NAMES   = ['zone_mass', 'mass', 'cell_mass']
FUSION_NAMES      = ['FusionRate_DT_nHe4', 'fusion_rate_DT', 'fusion_power']
RI_NAME           = 'Indices at region interfaces'

# Canonical Olson PDD region zoning (last-zone-of-region convention)
OLSON_PDD_CANONICAL = (150, 190, 320)


def first_available(run: HeliosRun, candidates: List[str]) -> Optional[str]:
    avail = set(run.dataset.variables.keys())
    for name in candidates:
        if name in avail:
            return name
    return None


def load_var(run: HeliosRun, candidates: List[str], required: bool = False,
             label: str = '') -> Tuple[Optional[str], Optional[np.ndarray]]:
    name = first_available(run, candidates)
    if name is None:
        if required:
            raise KeyError(
                f"Required variable '{label or candidates[0]}' not in EXODUS. "
                f"Tried: {candidates}")
        return None, None
    return name, np.asarray(run.get_variable(name))


def load_region_zones(run: HeliosRun) -> Tuple[Tuple[int, ...], str]:
    if RI_NAME in run.dataset.variables:
        ri_var = np.asarray(run.dataset.variables[RI_NAME][:])
        ri_raw = tuple(int(x) for x in ri_var[0, :-1])
        ri_zone = tuple(x - 1 for x in ri_raw)
        return ri_zone, (f"EXODUS '{RI_NAME}' at t=0: raw={list(ri_raw)} "
                         f"→ last-zone-of-region {list(ri_zone)}")
    zbar_name, zbar = load_var(run, ZBAR_NAMES)
    if zbar is not None:
        dz = np.diff(zbar[0])
        jumps = np.where(np.abs(dz) > 0.1)[0]
        coalesced = []
        last = -10
        for j in jumps:
            if j - last > 1:
                coalesced.append(int(j))
            last = j
        if coalesced:
            return tuple(coalesced), f"Z̄-jump autodetect: {coalesced}"
    return OLSON_PDD_CANONICAL, f"fallback Olson canonical {list(OLSON_PDD_CANONICAL)}"


def label_region(zone_idx: int, ri_zone: Tuple[int, ...],
                 region_names: List[str]) -> str:
    for k, last_zone in enumerate(ri_zone):
        if zone_idx <= last_zone:
            return region_names[k] if k < len(region_names) else f"R{k+1}"
    last_k = len(ri_zone)
    return region_names[last_k] if last_k < len(region_names) else f"R{last_k+1}"


def first_crossing_idx(arr_t: np.ndarray, threshold: float) -> int:
    """Return the first index where arr_t crosses `threshold` (rising), or
    -1 if it never crosses."""
    above = arr_t >= threshold
    if not above.any():
        return -1
    return int(np.argmax(above))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=("Per-zone burn-rate timing and compression-state "
                     "diagnostic.  Discriminates compression-deficit vs "
                     "burn-duration mechanisms for PDD priority 0a."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('base_path',
                        help='Run base path WITHOUT .exo extension.')
    parser.add_argument('--t-bootstrap-keV', type=float, default=4.5,
                        help='T_ion threshold for alpha-bootstrap onset.')
    parser.add_argument('--rate-frac', type=float, default=0.1,
                        help='Fraction of global peak rate at which a zone '
                             'is considered burning meaningfully.')
    parser.add_argument('--out-base', default=None,
                        help='Output base (without extension).')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    base = Path(args.base_path).expanduser().resolve()
    exo = base.with_suffix('.exo')
    if not exo.exists():
        print(f"ERROR: EXODUS file not found: {exo}", file=sys.stderr)
        return 1

    out_base = (Path(args.out_base).expanduser().resolve()
                if args.out_base else base)
    csv_path = out_base.parent / (out_base.name + '_burn_rate_timing.csv')
    txt_path = out_base.parent / (out_base.name + '_burn_rate_timing.txt')

    print(f"Input EXODUS: {exo}")
    print(f"Output CSV:   {csv_path}")
    print(f"Output TXT:   {txt_path}")
    print()

    run = HeliosRun(str(exo))

    zbnd_name,  zbnd  = load_var(run, ZBND_NAMES, required=True)
    rho_name,   rho   = load_var(run, RHO_NAMES,  required=True)
    Tion_name,  T_ion = load_var(run, TION_NAMES, required=True)
    Pi_name,    P_i   = load_var(run, PI_NAMES,   required=True)
    Pe_name,    P_e   = load_var(run, PE_NAMES,   required=False)
    zmass_name, zmass = load_var(run, ZONE_MASS_NAMES, required=True)
    fusion_name, fusion = load_var(run, FUSION_NAMES, required=True)

    if P_e is None:
        P_e = np.zeros_like(P_i)

    t_s   = np.asarray(run.dataset.variables['time_whole'][:])
    t_ns  = t_s * 1e9
    n_t, n_z = rho.shape

    # ── Region structure
    ri_zone, ri_src = load_region_zones(run)
    region_names: List[str]
    if len(ri_zone) == 3:
        region_names = ['DT_vapor', 'DT_ice', 'DT-CH_foam', 'CH_skin']
    else:
        region_names = [f"R{k+1}" for k in range(len(ri_zone) + 1)]
    n_regions = len(ri_zone) + 1
    region_zone_range: List[Tuple[int, int]] = []
    last_zone = -1
    for k in range(n_regions):
        z_lo = last_zone + 1
        z_hi = ri_zone[k] if k < len(ri_zone) else n_z - 1
        region_zone_range.append((z_lo, z_hi))
        last_zone = z_hi

    # ── Volume-integrated burn rate (global) for bang time
    V_zone = (4.0 / 3.0) * np.pi * (zbnd[:, 1:n_z + 1] ** 3
                                    - zbnd[:, :n_z] ** 3)        # (n_t, n_z) cm³
    rate_zone = fusion * V_zone                                  # reactions/s, per zone
    rate_global = rate_zone.sum(axis=1)                          # reactions/s
    bang_idx = int(np.argmax(rate_global))
    bang_t_ns = float(t_ns[bang_idx])
    peak_global_rate = float(rate_global[bang_idx])

    rate_threshold = args.rate_frac * peak_global_rate / n_z  # rate/zone threshold

    # ── Per-zone scalars
    t_peak_burn       = np.full(n_z, np.nan)
    peak_local_rate   = np.zeros(n_z)
    t_peak_rho        = np.full(n_z, np.nan)
    peak_rho          = np.zeros(n_z)
    init_rho          = rho[0, :].copy()
    t_peak_Tion       = np.full(n_z, np.nan)
    peak_Tion         = np.zeros(n_z)
    t_cross_bootstrap = np.full(n_z, np.nan)
    P_total_at_peak_rho = np.zeros(n_z)
    P_total_at_bang     = (P_i[bang_idx] + P_e[bang_idx])

    for z in range(n_z):
        rate_z = rate_zone[:, z]
        i_peak = int(np.argmax(rate_z))
        peak_local_rate[z] = float(rate_z[i_peak])
        if peak_local_rate[z] > 0:
            t_peak_burn[z] = float(t_ns[i_peak])

        i_peak_rho = int(np.argmax(rho[:, z]))
        peak_rho[z] = float(rho[i_peak_rho, z])
        t_peak_rho[z] = float(t_ns[i_peak_rho])
        P_total_at_peak_rho[z] = float(P_i[i_peak_rho, z] + P_e[i_peak_rho, z])

        i_peak_T = int(np.argmax(T_ion[:, z]))
        peak_Tion[z] = float(T_ion[i_peak_T, z]) / 1000.0  # eV → keV
        t_peak_Tion[z] = float(t_ns[i_peak_T])

        ix = first_crossing_idx(T_ion[:, z] / 1000.0, args.t_bootstrap_keV)
        if ix >= 0:
            t_cross_bootstrap[z] = float(t_ns[ix])

    compression_ratio = np.where(init_rho > 0, peak_rho / init_rho, 0.0)

    # ── Region aggregates
    def _agg(arr, mask, weights=None):
        """Return (count, min, median, max, mean) over `arr[mask]`.
        `weights` is for the mass-weighted mean; if None, uses unweighted."""
        vals = arr[mask]
        valid = vals[np.isfinite(vals)]
        if valid.size == 0:
            return 0, np.nan, np.nan, np.nan, np.nan
        if weights is not None:
            w = weights[mask]
            wv = w[np.isfinite(vals)]
            mean = float(np.sum(valid * wv) / np.sum(wv)) if np.sum(wv) > 0 else np.nan
        else:
            mean = float(np.mean(valid))
        return (int(valid.size), float(np.min(valid)),
                float(np.median(valid)), float(np.max(valid)), mean)

    # ── Console summary
    summary: List[str] = []
    def _emit(s: str = "") -> None:
        summary.append(s)
        print(s)

    _emit("=" * 84)
    _emit("BURN-RATE TIMING + COMPRESSION-STATE DIAGNOSTIC  (priority 0a, mech (1)/(3) vs (2))")
    _emit("=" * 84)
    _emit(f"Run:           {base.name}")
    _emit(f"EXODUS:        {exo}")
    _emit(f"Timesteps:     n_t = {n_t}, n_zones = {n_z}")
    _emit(f"Bang time:     t = {bang_t_ns:.4f} ns  (idx {bang_idx}), "
          f"global peak rate {peak_global_rate:.3e} reactions/s")
    _emit(f"T_bootstrap:   {args.t_bootstrap_keV:.2f} keV (T_ion threshold)")
    _emit(f"Region map:    {ri_src}")
    _emit("")

    # ── Region compression table (EOS-direct diagnostic)
    _emit("Compression state (sensitive to EOS — PROPACEOS foam vs DT cryo):")
    _emit(f"  {'region':>13} {'zones':>10}  {'rho0_gcc':>10} "
          f"{'rho_max_gcc':>12} {'<CR>':>7} {'CR_med':>7} {'CR_max':>7} "
          f"{'t@CRmax_ns':>11}")
    _emit("  " + "-" * 84)
    for k in range(n_regions):
        z_lo, z_hi = region_zone_range[k]
        mask = np.zeros(n_z, dtype=bool); mask[z_lo:z_hi + 1] = True
        rho0_mean = float(np.mean(init_rho[mask]))
        rho_max_max = float(np.max(peak_rho[mask]))
        _, cr_min, cr_med, cr_max, cr_mean = _agg(compression_ratio, mask,
                                                  weights=zmass[0])
        # Time at which the highest-compression zone in this region hit its
        # peak ρ — informative for the EOS-vs-timing question
        z_at_crmax = z_lo + int(np.argmax(peak_rho[z_lo:z_hi + 1]))
        t_at_crmax = t_peak_rho[z_at_crmax]
        _emit(f"  {region_names[k]:>13} {z_lo:>4d}..{z_hi:<4d}  "
              f"{rho0_mean:>10.4f} {rho_max_max:>12.3f} "
              f"{cr_mean:>7.2f} {cr_med:>7.2f} {cr_max:>7.2f} "
              f"{t_at_crmax:>11.4f}")
    _emit("")
    _emit("  EOS-anomaly signal: if DT-CH_foam's CR_max is much smaller than")
    _emit("  DT_ice's CR_max (say <30%), PROPACEOS foam table is the suspect.")
    _emit("  If foam CR_max is comparable to ice but offset in TIME (peak ρ")
    _emit("  reached far AFTER ice's), the compression is fine but bootstrap")
    _emit("  arrived late → mechanism (2).")
    _emit("")

    # ── Region timing table
    _emit("Burn-rate timing per region:")
    _emit(f"  {'region':>13} {'n_burned':>10} {'n_xT':>6} "
          f"{'t_first_xT_ns':>14} {'t_med_xT_ns':>13} "
          f"{'t_first_burn_ns':>16} {'t_med_burn_ns':>15}")
    _emit("  " + "-" * 84)
    for k in range(n_regions):
        z_lo, z_hi = region_zone_range[k]
        mask = np.zeros(n_z, dtype=bool); mask[z_lo:z_hi + 1] = True
        n_burned = int(np.sum(peak_local_rate[mask] > 0))
        crossings = t_cross_bootstrap[mask]
        crossings_valid = crossings[np.isfinite(crossings)]
        n_xT = crossings_valid.size
        t_first_xT = float(np.min(crossings_valid)) if n_xT > 0 else float('nan')
        t_med_xT   = float(np.median(crossings_valid)) if n_xT > 0 else float('nan')
        peaks = t_peak_burn[mask]
        peaks_valid = peaks[np.isfinite(peaks)]
        t_first_burn = float(np.min(peaks_valid)) if peaks_valid.size > 0 else float('nan')
        t_med_burn   = float(np.median(peaks_valid)) if peaks_valid.size > 0 else float('nan')
        _emit(f"  {region_names[k]:>13} {n_burned:>5d}/{z_hi-z_lo+1:<4d} "
              f"{n_xT:>5d} {t_first_xT:>14.4f} {t_med_xT:>13.4f} "
              f"{t_first_burn:>16.4f} {t_med_burn:>15.4f}")
    _emit("")
    _emit("  Reading the timing: ice should reach T_bootstrap WELL BEFORE bang time,")
    _emit(f"  foam shortly after.  If t_first_xT_foam ≳ bang_time ({bang_t_ns:.3f} ns)")
    _emit("  the foam never gets hot in time → mechanism (2).")
    _emit("  If t_first_xT_foam < bang_time but n_xT_foam ≪ n_zones_foam, only the")
    _emit("  innermost foam shell crosses while the rest doesn't → mechanism (1)/(3).")
    _emit("")

    # ── Snapshot at bang time per region (compression and pressure)
    _emit("Snapshot at bang time (mass-weighted avg by zone-mass at bang):")
    _emit(f"  {'region':>13} {'rho_avg_gcc':>13} {'rho_max_gcc':>13} "
          f"{'Tion_avg_keV':>14} {'Tion_max_keV':>14} {'P_avg_Gbar':>12}")
    _emit("  " + "-" * 84)
    rho_bang = rho[bang_idx]
    Tion_bang_keV = T_ion[bang_idx] / 1000.0
    P_total_bang = (P_i[bang_idx] + P_e[bang_idx])   # J/cm³
    mass_bang = zmass[bang_idx]
    for k in range(n_regions):
        z_lo, z_hi = region_zone_range[k]
        sl = slice(z_lo, z_hi + 1)
        m = mass_bang[sl]
        m_sum = float(np.sum(m)) if np.sum(m) > 0 else 1.0
        rho_avg = float(np.sum(rho_bang[sl] * m) / m_sum)
        rho_max = float(np.max(rho_bang[sl]))
        Tion_avg = float(np.sum(Tion_bang_keV[sl] * m) / m_sum)
        Tion_max = float(np.max(Tion_bang_keV[sl]))
        P_avg_Gbar = float(np.sum(P_total_bang[sl] * m) / m_sum) * 1e-8
        _emit(f"  {region_names[k]:>13} {rho_avg:>13.3f} {rho_max:>13.3f} "
              f"{Tion_avg:>14.3f} {Tion_max:>14.3f} {P_avg_Gbar:>12.3f}")
    _emit("")

    # ── Top-10 lowest-compression burning zones in foam
    foam_lo, foam_hi = region_zone_range[2] if n_regions >= 3 else (0, -1)
    if foam_hi >= foam_lo:
        foam_burning_mask = (peak_local_rate > 0)
        foam_burning_mask[:foam_lo] = False
        foam_burning_mask[foam_hi + 1:] = False
        if foam_burning_mask.any():
            idx_sorted = np.argsort(compression_ratio[foam_burning_mask])
            foam_idx = np.where(foam_burning_mask)[0]
            _emit("Lowest-CR foam zones that DID fire (any rate > 0):")
            _emit(f"  {'rank':>4} {'zone':>5} {'rho0':>7} {'rho_max':>9} "
                  f"{'CR':>6} {'t_peakRho':>10} {'t_xTbs':>9} "
                  f"{'t_peakBurn':>11} {'peak_Tion_keV':>14}")
            _emit("  " + "-" * 76)
            for rank, kk in enumerate(idx_sorted[:5], 1):
                z = int(foam_idx[kk])
                _emit(f"  {rank:>4d} {z:>5d} {init_rho[z]:>7.4f} "
                      f"{peak_rho[z]:>9.3f} {compression_ratio[z]:>6.1f} "
                      f"{t_peak_rho[z]:>10.4f} "
                      f"{t_cross_bootstrap[z] if np.isfinite(t_cross_bootstrap[z]) else float('nan'):>9.4f} "
                      f"{t_peak_burn[z]:>11.4f} {peak_Tion[z]:>14.3f}")
            _emit("")
    _emit("=" * 84)

    # ── CSV
    with open(csv_path, 'w') as f:
        f.write(f"# Burn-rate timing and compression diagnostic for {base.name}\n")
        f.write(f"# bang_t_ns = {bang_t_ns:.6f}  (idx {bang_idx})\n")
        f.write(f"# T_bootstrap_keV = {args.t_bootstrap_keV}\n")
        f.write(f"# region interfaces (last-zone): {list(ri_zone)}\n")
        f.write(f"# region names: {region_names}\n")
        f.write("zone_idx,region,rho0_gcc,peak_rho_gcc,compression_ratio,"
                "t_peak_rho_ns,P_at_peak_rho_Gbar,"
                "peak_Tion_keV,t_peak_Tion_ns,t_cross_bootstrap_ns,"
                "peak_local_rate_per_s,t_peak_burn_ns,rho_at_bang_gcc,"
                "Tion_at_bang_keV\n")
        for z in range(n_z):
            reg = label_region(z, ri_zone, region_names)
            f.write(f"{z},{reg},"
                    f"{init_rho[z]:.5e},{peak_rho[z]:.5e},{compression_ratio[z]:.5e},"
                    f"{t_peak_rho[z]:.5e},{P_total_at_peak_rho[z]*1e-8:.5e},"
                    f"{peak_Tion[z]:.5e},{t_peak_Tion[z]:.5e},"
                    f"{t_cross_bootstrap[z]:.5e},"
                    f"{peak_local_rate[z]:.5e},{t_peak_burn[z]:.5e},"
                    f"{rho[bang_idx, z]:.5e},{T_ion[bang_idx, z]/1000.0:.5e}\n")

    with open(txt_path, 'w') as f:
        f.write("\n".join(summary) + "\n")

    print(f"Wrote {csv_path}  ({n_z} zones)")
    print(f"Wrote {txt_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
