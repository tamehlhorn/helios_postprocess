"""
dump_per_zone_burn_share.py
===========================

Per-zone cumulative DT neutron count at end-of-run, aggregated by region.
Answers: of the total yield, which zones actually burned?

Motivation (priority 0a follow-up, May 24 2026)
-----------------------------------------------
The first burn-propagation diagnostic (`dump_burn_propagation_profile.py`)
showed the ice/foam interface in fab007 is essentially continuous at peak
burn: T_e ratio 1.017, T_rad 1.000, α-deposition rate ratio 1.062, Z̄
jump only 8% (not the factor-of-3 the original priority-0a framing
assumed).  Both leading hypotheses (brems sink at Z̄ jump, non-local
α-transport range-out) are refuted by that one-time-slice diagnostic.

But that diagnostic is a point-in-time snapshot.  It doesn't directly
say whether the foam ACTUALLY produced yield commensurate with its
mass.  The DT-CH wetted foam is mostly DT (0.222 g/cc DT in a 0.020
g/cc CH matrix, ~92% DT by mass), so per-zone yield share is the right
quantitative cut: if foam zones contribute proportionally to their
DT mass, "burn propagates into foam" is true and the residual is in
compression/timing.  If foam zones contribute proportionally LESS than
ice zones, "burn front stalls" survives but in a refined form (it
stalls in TIME, not at the material interface).

EXODUS gives us this directly: `TimeIntFusionProd_n_1406_zone` is the
per-zone cumulative DT reaction count, shape (n_t, n_zones).  The
last timestep value is the final answer for that run.

Each DT reaction releases 17.6 MeV = 2.819e-12 J total energy
(3.5 MeV α + 14.1 MeV neutron).  We dump per-zone yield in kJ, the
fractional contribution of each zone to the total, and a region-
aggregated breakdown.

Usage
-----
    python3 dump_per_zone_burn_share.py <base_path> [options]

Produces alongside the .exo:

    <base>_per_zone_burn_share.csv     -- zone-by-zone yield share
    <base>_per_zone_burn_share.txt     -- console summary mirrored to disk

Options
-------
    --time-idx INT
        Use this timestep index instead of the last (-1).  Negative
        supported.  Default -1 (end of simulation).
    --out-base PATH
        Override output base.  Default: alongside the .exo.
    --verbose

Standalone — uses HeliosRun, doesn't require ICFAnalyzer to have run.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from helios_postprocess import HeliosRun


# EXODUS variable name resolution.  Same order/fallbacks as
# dump_burn_propagation_profile.py.
ZBND_NAMES         = ['zone_boundaries', 'coord', 'coordx']
ZBAR_NAMES         = ['mean_charge', 'Zbar']
ZONE_MASS_NAMES    = ['zone_mass', 'mass', 'cell_mass']
DT_COUNT_ZONE_NAME = 'TimeIntFusionProd_n_1406_zone'
DT_COUNT_TOTAL_NAME = 'TimeIntFusionProd_n_1406'
RI_NAME            = 'Indices at region interfaces'
MATERIAL_NAME      = 'Material index'

# DT fusion reaction energy split
Q_DT_J         = 17.6e6 * 1.602e-19   # ≈ 2.819e-12 J / reaction
Q_ALPHA_FRAC   = 3.5 / 17.6
Q_NEUTRON_FRAC = 14.1 / 17.6

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


def load_region_zones(run: HeliosRun, n_z: int
                      ) -> Tuple[Tuple[int, ...], str]:
    """Resolve last-zone-of-region indices (zone-index convention)."""
    if RI_NAME in run.dataset.variables:
        ri_var = np.asarray(run.dataset.variables[RI_NAME][:])
        ri_raw = tuple(int(x) for x in ri_var[0, :-1])
        ri_zone = tuple(x - 1 for x in ri_raw)
        ri_src = (f"EXODUS '{RI_NAME}' at t=0: raw={list(ri_raw)} "
                  f"(node/slice-end) → last-zone-of-region {list(ri_zone)}")
        return ri_zone, ri_src
    # Z̄-jump autodetect fallback
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
    return OLSON_PDD_CANONICAL, (
        f"fallback Olson_PDD canonical {list(OLSON_PDD_CANONICAL)}")


def label_region(zone_idx: int, ri_zone: Tuple[int, ...],
                 region_names: List[str]) -> str:
    for k, last_zone in enumerate(ri_zone):
        if zone_idx <= last_zone:
            return region_names[k] if k < len(region_names) else f"R{k+1}"
    last_k = len(ri_zone)
    return region_names[last_k] if last_k < len(region_names) else f"R{last_k+1}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=("Per-zone DT yield share at end-of-run, region-aggregated. "
                     "Answers 'which zones actually burned' for the PDD priority "
                     "0a investigation."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('base_path',
                        help='Run base path WITHOUT .exo extension.')
    parser.add_argument('--time-idx', type=int, default=-1,
                        help='Timestep index (negative supported).  Default -1.')
    parser.add_argument('--out-base', default=None,
                        help='Output base path (without extension).')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    base = Path(args.base_path).expanduser().resolve()
    exo = base.with_suffix('.exo')
    if not exo.exists():
        print(f"ERROR: EXODUS file not found: {exo}", file=sys.stderr)
        return 1

    out_base = (Path(args.out_base).expanduser().resolve()
                if args.out_base else base)
    csv_path = out_base.parent / (out_base.name + '_per_zone_burn_share.csv')
    txt_path = out_base.parent / (out_base.name + '_per_zone_burn_share.txt')

    print(f"Input EXODUS: {exo}")
    print(f"Output CSV:   {csv_path}")
    print(f"Output TXT:   {txt_path}")
    print()

    run = HeliosRun(str(exo))

    zbnd_name, zbnd = load_var(run, ZBND_NAMES, required=True,
                               label='zone_boundaries')
    zmass_name, zone_mass = load_var(run, ZONE_MASS_NAMES, required=True,
                                     label='zone_mass')

    if DT_COUNT_ZONE_NAME not in run.dataset.variables:
        print(f"ERROR: '{DT_COUNT_ZONE_NAME}' not in EXODUS.  This run "
              f"probably has burn OFF.  Burn-share is undefined.",
              file=sys.stderr)
        return 1
    dt_count_zone = np.asarray(run.get_variable(DT_COUNT_ZONE_NAME))  # (n_t, n_z)

    dt_count_total = None
    if DT_COUNT_TOTAL_NAME in run.dataset.variables:
        dt_count_total = np.asarray(run.get_variable(DT_COUNT_TOTAL_NAME))

    n_t, n_z = dt_count_zone.shape
    t_idx = args.time_idx if args.time_idx >= 0 else n_t + args.time_idx
    if not (0 <= t_idx < n_t):
        print(f"ERROR: time-idx {args.time_idx} out of range for n_t={n_t}",
              file=sys.stderr)
        return 1

    t_s = np.asarray(run.dataset.variables['time_whole'][:])
    t_ns = t_s * 1e9
    t_at_idx = float(t_ns[t_idx])

    # Per-zone cumulative count and energy at the chosen timestep
    count_zone   = dt_count_zone[t_idx]                          # (n_z,)
    yield_zone_J = count_zone * Q_DT_J
    yield_total_zone_J = float(yield_zone_J.sum())

    yield_total_global_J = (float(dt_count_total[t_idx]) * Q_DT_J
                            if dt_count_total is not None else None)

    if yield_total_zone_J <= 0:
        print(f"WARNING: zero total yield at t_idx={t_idx} "
              f"(t={t_at_idx:.4f} ns).  Either burn never caught or you "
              f"selected a pre-burn timestep.")

    # Region structure
    ri_zone, ri_src = load_region_zones(run, n_z)
    # Region names — read from EXODUS material assignment if possible, else
    # use Olson PDD canonical names for 4-region targets.
    region_names: List[str]
    if len(ri_zone) == 3:
        region_names = ['DT_vapor', 'DT_ice', 'DT-CH_foam', 'CH_skin']
    else:
        region_names = [f"R{k+1}" for k in range(len(ri_zone) + 1)]

    # Per-zone initial radius (cm) and radius at chosen timestep
    r0_um   = 0.5 * (zbnd[0, :-1] + zbnd[0, 1:])  * 1e4
    rT_um   = 0.5 * (zbnd[t_idx, :-1] + zbnd[t_idx, 1:]) * 1e4
    zmass_g = zone_mass[t_idx]                                   # (n_z,) g

    # Cumulative running sum (fraction of total, from inner to outer)
    cum_frac = (np.cumsum(yield_zone_J) / yield_total_zone_J
                if yield_total_zone_J > 0 else np.zeros(n_z))

    # ── Aggregate by region
    n_regions = len(ri_zone) + 1
    region_yield_J = np.zeros(n_regions)
    region_mass_g  = np.zeros(n_regions)
    region_zone_range: List[Tuple[int, int]] = []

    last_zone = -1
    for k in range(n_regions):
        z_lo = last_zone + 1
        z_hi = ri_zone[k] if k < len(ri_zone) else n_z - 1
        region_yield_J[k] = yield_zone_J[z_lo:z_hi + 1].sum()
        region_mass_g[k]  = zmass_g[z_lo:z_hi + 1].sum()
        region_zone_range.append((z_lo, z_hi))
        last_zone = z_hi

    # ── Console summary
    summary: List[str] = []
    def _emit(s: str = "") -> None:
        summary.append(s)
        print(s)

    _emit("=" * 78)
    _emit("PER-ZONE BURN SHARE  (priority 0a quantitative follow-up)")
    _emit("=" * 78)
    _emit(f"Run:         {base.name}")
    _emit(f"EXODUS:      {exo}")
    _emit(f"Timesteps:   n_t = {n_t}, n_zones = {n_z}")
    _emit(f"Chosen idx:  {t_idx}  →  t = {t_at_idx:.4f} ns")
    _emit("")
    _emit(f"Region interfaces: {ri_src}")
    _emit("")
    _emit("Totals:")
    _emit(f"  Σ per-zone DT reactions     = {count_zone.sum():.4e}")
    _emit(f"  Σ per-zone DT yield         = "
          f"{yield_total_zone_J:.4e} J  "
          f"= {yield_total_zone_J*1e-6:.4f} MJ "
          f"(total energy released, n + α)")
    _emit(f"  Σ per-zone α-deposited      = "
          f"{yield_total_zone_J*Q_ALPHA_FRAC*1e-6:.4f} MJ  "
          f"(nominal 3.5/17.6 split)")
    _emit(f"  Σ per-zone neutron-escape   = "
          f"{yield_total_zone_J*Q_NEUTRON_FRAC*1e-6:.4f} MJ")
    if yield_total_global_J is not None:
        _emit(f"  Global '{DT_COUNT_TOTAL_NAME}' yield = "
              f"{yield_total_global_J*1e-6:.4f} MJ "
              f"(cross-check, should match Σ to floating-point)")
        if yield_total_zone_J > 0:
            _emit(f"  Closure                     = "
                  f"{100*yield_total_zone_J/yield_total_global_J:.3f}% of global")
    _emit("")

    # Region table
    _emit("Region-aggregated burn share:")
    _emit(f"  {'region':>12}  {'zones':>10}  {'mass_mg':>9} "
          f"{'yield_kJ':>11} {'yield_MJ':>10} {'%_total':>9} "
          f"{'yield/mass kJ/mg':>17}")
    _emit("  " + "-" * 86)
    total_y = yield_total_zone_J
    for k in range(n_regions):
        name = region_names[k]
        z_lo, z_hi = region_zone_range[k]
        m_mg = region_mass_g[k] * 1e3
        y_J  = region_yield_J[k]
        y_kJ = y_J * 1e-3
        y_MJ = y_J * 1e-6
        pct  = 100 * y_J / total_y if total_y > 0 else 0.0
        ymm  = (y_kJ / m_mg) if m_mg > 0 else 0.0
        _emit(f"  {name:>12}  {z_lo:>4d}..{z_hi:<4d}  "
              f"{m_mg:>9.4f} {y_kJ:>11.3f} {y_MJ:>10.4f} "
              f"{pct:>8.3f}% {ymm:>17.3f}")
    _emit("")
    _emit("Interpretation:")
    _emit("  - Yield/mass ratio across DT-bearing regions tells us where the")
    _emit("    burn substrate actually fired.")
    _emit("  - DT_ice + DT-CH_foam are both ~92%+ DT by mass; if their")
    _emit("    yield/mass ratios are within ~2-3× of each other, burn IS")
    _emit("    propagating into the foam.")
    _emit("  - If DT-CH_foam yield/mass ≪ DT_ice yield/mass, the foam DT mass")
    _emit("    didn't fire — even though the interface looks continuous in 3T.")
    _emit("    That would point at TIMING (foam doesn't get hot until burn is")
    _emit("    almost over) or COMPRESSION (foam ρ too low, n²σv too small).")
    _emit("")

    # Top-10 yielding zones
    top_idx = np.argsort(yield_zone_J)[::-1][:10]
    _emit("Top-10 highest-yield zones:")
    _emit(f"  {'rank':>4} {'zone':>5} {'region':>13} "
          f"{'r0_um':>8} {'r_um':>9} "
          f"{'yield_kJ':>10} {'%_total':>9} {'%_cum':>8}")
    _emit("  " + "-" * 76)
    for rank, k in enumerate(top_idx, 1):
        reg = label_region(int(k), ri_zone, region_names)
        pct = 100 * yield_zone_J[k] / total_y if total_y > 0 else 0.0
        _emit(f"  {rank:>4d} {int(k):>5d} {reg:>13} "
              f"{r0_um[k]:>8.2f} {rT_um[k]:>9.2f} "
              f"{yield_zone_J[k]*1e-3:>10.3f} {pct:>8.3f}% "
              f"{100*cum_frac[k]:>7.2f}%")
    _emit("")

    # Where cumulative crosses 50% and 90% (yield radius)
    if total_y > 0:
        i50 = int(np.searchsorted(cum_frac, 0.50))
        i90 = int(np.searchsorted(cum_frac, 0.90))
        _emit(f"Yield-cumulative radius:")
        _emit(f"  50% of yield from zones 0..{i50}  (r0 ≤ {r0_um[i50]:.2f} µm, "
              f"r_at_t ≤ {rT_um[i50]:.2f} µm)   region: "
              f"{label_region(i50, ri_zone, region_names)}")
        _emit(f"  90% of yield from zones 0..{i90}  (r0 ≤ {r0_um[i90]:.2f} µm, "
              f"r_at_t ≤ {rT_um[i90]:.2f} µm)   region: "
              f"{label_region(i90, ri_zone, region_names)}")
        _emit("")
    _emit("=" * 78)

    # ── Write CSV
    with open(csv_path, 'w') as f:
        f.write(f"# Per-zone burn share for {base.name}\n")
        f.write(f"# t_idx = {t_idx}  t_ns = {t_at_idx:.6f}\n")
        f.write(f"# region interfaces (last-zone): {list(ri_zone)}\n")
        f.write(f"# region names: {region_names}\n")
        f.write(f"# Q_DT = 17.6 MeV = {Q_DT_J:.6e} J per reaction\n")
        f.write(f"# Σ yield (zone) = {yield_total_zone_J:.6e} J\n")
        if yield_total_global_J is not None:
            f.write(f"# Σ yield (global) = {yield_total_global_J:.6e} J\n")
        f.write("zone_idx,region,r0_um,r_at_t_um,zone_mass_mg,"
                "dt_count_cumulative,yield_J,yield_kJ,frac_pct,cum_frac_pct\n")
        for k in range(n_z):
            reg = label_region(k, ri_zone, region_names)
            pct = 100 * yield_zone_J[k] / total_y if total_y > 0 else 0.0
            f.write(f"{k},{reg},{r0_um[k]:.5e},{rT_um[k]:.5e},"
                    f"{zmass_g[k]*1e3:.5e},{count_zone[k]:.5e},"
                    f"{yield_zone_J[k]:.5e},{yield_zone_J[k]*1e-3:.5e},"
                    f"{pct:.5f},{100*cum_frac[k]:.5f}\n")

    with open(txt_path, 'w') as f:
        f.write("\n".join(summary) + "\n")

    print(f"Wrote {csv_path}  ({n_z} zones)")
    print(f"Wrote {txt_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
