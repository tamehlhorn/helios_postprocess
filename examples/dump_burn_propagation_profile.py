"""
dump_burn_propagation_profile.py
================================

Zone-by-zone radial dump of T_e, T_rad, T_ion, ρ and α-deposition rate
at peak burn time, with special focus on the DT-ice / DT-CH-foam
material interface.

Motivation (priority 0a of the PDD closeout, May 23 2026)
---------------------------------------------------------
The fab007 production calibration (`Olson_PDD_20_fab007_foot25_s018_c37_burn`)
ignites the DT-ice layer but the burn front does **not** propagate across
the Z̄ discontinuity into the DT-CH wetted foam.  The ice→foam-swap
experiment collapses burn entirely, confirming Helios's burn substrate is
the pure-DT ice ONLY.  Three candidate mechanisms (ranked by hypothesis
quality):

    (1) Bremsstrahlung losses across the Z̄ jump — radiative sink.
    (2) Non-local α-transport mean-free-path scaling into Z̄≈3 mixed
        material — alphas deposit before reaching the foam.
    (3) PROPACEOS foam EOS/opacity vs LILAC's foam EOS.
    (4) 1D Lagrangian mix model at the interface.

Diagnostic strategy (priority 0a part (i)): pull T_e and T_rad and the
α-deposition rate radial profiles at peak burn.

    * Sharp T_e drop (2-5×) across the ice/foam boundary → mechanism (1).
    * T_e continuous but α-deposition-rate falls off sharply → mech (2).
    * Both continuous → look elsewhere (EOS, mix, timing).

Usage
-----
    python3 dump_burn_propagation_profile.py <base_path> [options]

`base_path` is the run's EXODUS base path WITHOUT extension (the same
convention as `run_analysis.py`):

    ~/Sims/Xcimer/Olson_PDD/Olson_PDD_20/Olson_PDD_20_fab007_foot25_s018_c37_burn

This produces (alongside the .exo file):

    <base>_burn_propagation_profile.csv     -- zone-by-zone radial profile
    <base>_burn_propagation_profile.txt     -- console summary mirrored to disk

and prints the same summary to stdout.  Both files are co-located with
the .exo so the diagnostic stays with the run.

Outputs are at peak burn time by default (argmax of the global
DT fusion rate).  Pass `--time-mode ignition` to use the moment ρR_hs
first crosses 0.3 g/cm² instead, or `--t-ns <T>` for an explicit time.

Options
-------
    --time-mode {peak-burn,ignition,t-ns}
        Which timestep to dump.  Default peak-burn.
    --t-ns FLOAT
        Explicit time in nanoseconds; overrides --time-mode.
    --region-zones I1,I2,I3
        Region/material interface ZONE indices (1-D, time-invariant
        Lagrangian).  If omitted, read from the EXODUS
        `Indices at region interfaces` variable; if THAT is also absent
        fall back to Olson_PDD canonical zoning (150,190,320) and warn.
    --ice-foam-half-window N
        Number of zones on each side of the ice/foam interface to print
        in the console summary table (default 8).
    --out-base PATH
        Override the output base path (default: alongside the .exo).
    --verbose

Standalone — uses HeliosRun for raw .exo access; does NOT require the
full ICFAnalyzer pipeline to have run on this base.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from helios_postprocess import HeliosRun


# ── EXODUS variable name candidates ─────────────────────────────────────
# First match wins.  Names taken from `helios_exodus_variable_reference.md`.
RHO_NAMES    = ['mass_density', 'dens', 'density']
TION_NAMES   = ['ion_temperature', 'temp', 'temperature']
TE_NAMES     = ['elec_temperature', 'electron_temperature']
TRAD_NAMES   = ['radiation_temperature']
ZBND_NAMES   = ['zone_boundaries', 'coord', 'coordx']
ZBAR_NAMES   = ['mean_charge', 'Zbar']
FUSION_NAMES = ['FusionRate_DT_nHe4', 'fusion_rate_DT', 'fusion_power']
ALPHA_TOT_NAMES = ['alpha_heating_power', 'alpha_power', 'alpha_heating']
ALPHA_ION_NAMES = ['pt_particle_heating_ion']
ALPHA_ELE_NAMES = ['pt_particle_heating_ele']
RI_NAME         = 'Indices at region interfaces'

# Canonical Olson_PDD zone-region boundaries (zones 0..N inclusive end-exclusive)
# Region 1 (DT vapor):     0   .. 150
# Region 2 (DT solid ice): 151 .. 190
# Region 3 (DT-CH foam):   191 .. 320
# Region 4 (CH skin):      321 .. 350
OLSON_PDD_CANONICAL = (150, 190, 320)


# ── helpers ─────────────────────────────────────────────────────────────
def first_available(run: HeliosRun, candidates: List[str]) -> Optional[str]:
    """Return the first variable name in `candidates` present in the EXODUS
    file, or None if none are present."""
    avail = set(run.dataset.variables.keys())
    for name in candidates:
        if name in avail:
            return name
    return None


def load_var(run: HeliosRun, candidates: List[str], required: bool = False,
             label: str = '') -> Tuple[Optional[str], Optional[np.ndarray]]:
    """Resolve a candidate list to a variable and return (name, data)."""
    name = first_available(run, candidates)
    if name is None:
        if required:
            raise KeyError(
                f"Required variable '{label or candidates[0]}' not in EXODUS. "
                f"Tried: {candidates}")
        return None, None
    data = run.get_variable(name)
    return name, np.asarray(data)


def autodetect_region_zones(zbar0: np.ndarray, n_zones: int) -> Tuple[int, ...]:
    """Detect region/material interface zone indices from Z̄ at t=0.

    Returns the zone index BEFORE each Z̄ jump (i.e. last zone of the
    region on the inside of each interface).  These are time-invariant
    in Helios's Lagrangian scheme.
    """
    # Smooth tiny floating-point variations within a single material
    # region (Z̄ should be exactly constant within a pure region; with
    # mixtures e.g. DT-CH foam it's a smooth Zbar set by composition).
    # We flag a "jump" wherever |dZ̄/dz| > 0.1 (well above mixture
    # smoothness, well below the ~factor-of-3 inter-region jumps).
    dz = np.diff(zbar0)
    jumps = np.where(np.abs(dz) > 0.1)[0]
    if jumps.size == 0:
        return ()
    # Coalesce adjacent indices (a jump can occupy 1-2 zones depending on
    # how the mesh handles the boundary)
    coalesced = []
    last = -10
    for j in jumps:
        if j - last > 1:
            coalesced.append(int(j))
        last = j
    return tuple(coalesced)


def find_peak_burn_index(run: HeliosRun, zbnd: np.ndarray
                         ) -> Tuple[int, float, float, str]:
    """Find the bang time index = argmax of volume-integrated fusion rate.

    Returns (idx, t_ns, peak_rate, var_name_used).
    """
    fusion_name = first_available(run, FUSION_NAMES)
    if fusion_name is None:
        raise KeyError(f"No fusion-rate variable found (tried {FUSION_NAMES}).")

    fusion_zone = np.asarray(run.get_variable(fusion_name))   # (n_t, n_z)
    n_t, n_z = fusion_zone.shape

    # Volume-weight to get global reactions/s, since FusionRate_DT_nHe4
    # is reactions/cm³/s per zone.  Volume changes in time (Lagrangian
    # zone compresses), so this is a per-timestep volume integral.
    V_zone = (4.0 / 3.0) * np.pi * (zbnd[:, 1:n_z + 1] ** 3
                                    - zbnd[:, :n_z] ** 3)
    fusion_global = (fusion_zone * V_zone).sum(axis=1)
    idx = int(np.argmax(fusion_global))
    # `time_whole` (seconds) → ns
    t_s = np.asarray(run.dataset.variables['time_whole'][:])
    t_ns = float(t_s[idx] * 1e9)
    return idx, t_ns, float(fusion_global[idx]), fusion_name


def find_time_index(times_ns: np.ndarray, t_ns: float) -> int:
    """Return the closest time index to `t_ns`."""
    return int(np.argmin(np.abs(times_ns - t_ns)))


def find_ignition_index(run: HeliosRun, zbnd: np.ndarray, rho: np.ndarray,
                        T_ion: np.ndarray, ri_zone: Optional[Tuple[int, ...]]
                        ) -> Optional[int]:
    """Replicate the ICFAnalyzer ignition criterion: first time ρR_hs ≥ 0.3
    g/cm² where the hot-spot boundary is region-interface index 0 (i.e. the
    inner-most material interface).  Returns None if the file lacks the
    needed structure or never crosses the threshold.
    """
    if ri_zone is None or len(ri_zone) == 0:
        return None
    hs_bnd_zone = ri_zone[0]  # last zone of the hot-spot region
    n_t, n_z = rho.shape
    # ρR from r=0 inward boundary to the hot-spot outer boundary
    # (= node index hs_bnd_zone + 1 since region interface k = "last
    # zone of region k", node hs_bnd_zone+1 is the inner edge of next region)
    dr = np.diff(zbnd, axis=1)                                  # (n_t, n_z)
    rhoR_cum = np.cumsum(rho * dr, axis=1)                      # ρR(r) per t
    # ρR enclosed inside the HS = cumulative up to the HS-outer zone
    rhoR_hs = rhoR_cum[:, hs_bnd_zone]
    crossings = np.where(rhoR_hs >= 0.3)[0]
    return int(crossings[0]) if crossings.size else None


def label_region(zone_idx: int, ri_zone: Tuple[int, ...],
                 region_names: List[str]) -> str:
    """Look up the region name for a given zone index."""
    for k, last_zone in enumerate(ri_zone):
        if zone_idx <= last_zone:
            return region_names[k] if k < len(region_names) else f"R{k+1}"
    last_k = len(ri_zone)
    return region_names[last_k] if last_k < len(region_names) else f"R{last_k+1}"


# ── main ────────────────────────────────────────────────────────────────
def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Dump zone-by-zone (r, ρ, T_ion, T_e, T_rad, α-deposition-rate) "
            "at peak burn time, with focus on the ice/foam material interface. "
            "Diagnostic for PDD closeout priority 0a (burn-propagation arrest)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('base_path',
                        help='Run base path WITHOUT .exo extension.')
    parser.add_argument('--time-mode',
                        choices=['peak-burn', 'ignition', 't-ns'],
                        default='peak-burn',
                        help='Which timestep to dump.')
    parser.add_argument('--t-ns', type=float, default=None,
                        help='Explicit time in ns (forces --time-mode t-ns).')
    parser.add_argument('--region-zones', default=None,
                        help='Comma-separated zone indices for region '
                             'interfaces (e.g. "150,190,320").  Overrides '
                             'auto-detection.')
    parser.add_argument('--ice-foam-half-window', type=int, default=8,
                        help='Zones on each side of the ice/foam boundary '
                             'to print in the console summary.')
    parser.add_argument('--out-base', default=None,
                        help='Output base path (without extension).  '
                             'Default: alongside the input .exo.')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    base = Path(args.base_path).expanduser().resolve()
    exo = base.with_suffix('.exo')
    if not exo.exists():
        print(f"ERROR: EXODUS file not found: {exo}", file=sys.stderr)
        return 1

    out_base = (Path(args.out_base).expanduser().resolve()
                if args.out_base else base)
    csv_path = out_base.parent / (out_base.name + '_burn_propagation_profile.csv')
    txt_path = out_base.parent / (out_base.name + '_burn_propagation_profile.txt')

    print(f"Input EXODUS: {exo}")
    print(f"Output CSV:   {csv_path}")
    print(f"Output TXT:   {txt_path}")
    print()

    run = HeliosRun(str(exo))

    # ── Required geometry + state arrays
    zbnd_name, zbnd = load_var(run, ZBND_NAMES, required=True, label='zone_boundaries')
    rho_name,  rho  = load_var(run, RHO_NAMES,  required=True, label='mass_density')
    Tion_name, T_ion = load_var(run, TION_NAMES, required=True, label='ion_temperature')

    Te_name,   T_e   = load_var(run, TE_NAMES)
    Trad_name, T_rad = load_var(run, TRAD_NAMES)
    Zbar_name, Zbar  = load_var(run, ZBAR_NAMES)

    # ── α-deposition rate: prefer the directly-loaded variable, else
    # sum ion + electron heating, else fail with a useful message
    alpha_tot_name, alpha_tot = load_var(run, ALPHA_TOT_NAMES)
    alpha_ion_name, alpha_ion = load_var(run, ALPHA_ION_NAMES)
    alpha_ele_name, alpha_ele = load_var(run, ALPHA_ELE_NAMES)

    if alpha_tot is None:
        if alpha_ion is not None and alpha_ele is not None:
            alpha_tot = alpha_ion + alpha_ele
            alpha_src = (f"{alpha_ion_name} + {alpha_ele_name}  "
                         f"(no direct alpha_heating_power in this run)")
        elif alpha_ion is not None:
            alpha_tot = alpha_ion
            alpha_src = (f"{alpha_ion_name} only (ele not present — "
                         f"undercount of ~30-50% expected)")
        else:
            alpha_tot = np.zeros_like(rho)
            alpha_src = ("(NO alpha-deposition variables found — "
                         "burn flag OFF or alpha-transport disabled?)")
    else:
        alpha_src = alpha_tot_name

    # ── Times
    t_s   = np.asarray(run.dataset.variables['time_whole'][:])
    t_ns  = t_s * 1e9
    n_t, n_z = rho.shape

    # ── Region interface zone indices (time-invariant)
    if args.region_zones is not None:
        ri_zone = tuple(int(x) for x in args.region_zones.split(','))
        ri_src = f"--region-zones {args.region_zones} (CLI override)"
    elif RI_NAME in run.dataset.variables:
        ri_var = np.asarray(run.dataset.variables[RI_NAME][:])
        # Helios stores `ri[t, k]` as the NODE index of the inter-region
        # interface, which equals the slice-end (first-zone-of-NEXT-region).
        # E.g. Olson_PDD_9 stores 151 for the vapor/solid boundary because
        # zones 0..150 are vapor and zone 151 is the first ice zone.  This
        # script's `ri_zone` convention is "last zone of region k", so
        # subtract 1.  Drop the final entry (= n_zones, the outer
        # simulation boundary, not a true material interface).
        ri_raw = tuple(int(x) for x in ri_var[0, :-1])
        ri_zone = tuple(x - 1 for x in ri_raw)
        ri_src = (f"EXODUS '{RI_NAME}' at t=0: raw={list(ri_raw)} "
                  f"(node/slice-end) → last-zone-of-region {list(ri_zone)}")
    elif Zbar is not None:
        ri_zone = autodetect_region_zones(Zbar[0], n_z)
        ri_src = f"auto-detected from Z̄(t=0) jumps: {list(ri_zone)}"
        if ri_zone != OLSON_PDD_CANONICAL:
            print(f"  NOTE: auto-detected zones {list(ri_zone)} differ "
                  f"from Olson_PDD canonical {list(OLSON_PDD_CANONICAL)}; "
                  f"verify against run's target file.")
    else:
        ri_zone = OLSON_PDD_CANONICAL
        ri_src = (f"fallback (no Z̄, no '{RI_NAME}' variable): assuming "
                  f"Olson_PDD canonical {list(OLSON_PDD_CANONICAL)}")
        print("  WARNING: " + ri_src)

    # ── Region names for a 4-region Olson target.  If the target has a
    # different region count we still print, but with generic labels.
    if len(ri_zone) == 3:
        region_names = ['DT_vapor', 'DT_ice', 'DT-CH_foam', 'CH_skin']
    else:
        region_names = [f"R{k+1}" for k in range(len(ri_zone) + 1)]

    # ── Choose timestep
    if args.t_ns is not None:
        t_idx = find_time_index(t_ns, args.t_ns)
        t_mode_used = f"t-ns={args.t_ns} (closest EXODUS sample)"
    elif args.time_mode == 'peak-burn':
        t_idx, t_ns_peak, peak_rate, fusion_name = find_peak_burn_index(run, zbnd)
        t_mode_used = (f"peak-burn  (argmax volume-integrated {fusion_name}, "
                       f"peak rate {peak_rate:.3e} reactions/s)")
    elif args.time_mode == 'ignition':
        i_ign = find_ignition_index(run, zbnd, rho, T_ion, ri_zone)
        if i_ign is None:
            print("ERROR: --time-mode ignition requested but ρR_hs never "
                  "crosses 0.3 g/cm² (no ignition detected).", file=sys.stderr)
            return 1
        t_idx = i_ign
        t_mode_used = f"ignition (ρR_hs first crosses 0.3 g/cm²)"
    else:
        print(f"ERROR: --time-mode {args.time_mode} requires --t-ns to be "
              f"specified.", file=sys.stderr)
        return 1

    t_at_idx = float(t_ns[t_idx])

    # ── Build profile arrays at this timestep
    zb   = zbnd[t_idx]                                  # (n_z+1,) cm
    zc   = 0.5 * (zb[:-1] + zb[1:])                     # zone centers, cm
    dr   = np.diff(zb)                                  # zone widths, cm
    rho_t   = rho[t_idx]                                # g/cm³
    Ti_keV  = T_ion[t_idx] / 1000.0
    Te_keV  = (T_e[t_idx] / 1000.0
               if T_e is not None else np.full(n_z, np.nan))
    Trad_eV = (T_rad[t_idx]
               if T_rad is not None else np.full(n_z, np.nan))
    Zbar_t  = Zbar[t_idx] if Zbar is not None else np.full(n_z, np.nan)
    alpha_t = alpha_tot[t_idx]
    alpha_i = alpha_ion[t_idx] if alpha_ion is not None else np.full(n_z, np.nan)
    alpha_e = alpha_ele[t_idx] if alpha_ele is not None else np.full(n_z, np.nan)

    # ── Console summary
    summary = []
    def _emit(s: str = "") -> None:
        summary.append(s)
        print(s)

    _emit("=" * 76)
    _emit("BURN-PROPAGATION DIAGNOSTIC (priority 0a)")
    _emit("=" * 76)
    _emit(f"Run:         {base.name}")
    _emit(f"EXODUS:      {exo}")
    _emit(f"Timesteps:   n_t = {n_t}, n_zones = {n_z}")
    _emit(f"Time chosen: t = {t_at_idx:.4f} ns  (idx {t_idx})   mode: {t_mode_used}")
    _emit("")
    _emit("Variables found:")
    _emit(f"  zone_boundaries     -> {zbnd_name}")
    _emit(f"  mass_density        -> {rho_name}")
    _emit(f"  ion_temperature     -> {Tion_name}")
    _emit(f"  elec_temperature    -> {Te_name or '(missing)'}")
    _emit(f"  radiation_temperature -> {Trad_name or '(missing)'}")
    _emit(f"  mean_charge         -> {Zbar_name or '(missing)'}")
    _emit(f"  alpha_deposition    -> {alpha_src}")
    _emit("")
    _emit(f"Region interfaces (last zone of each region):")
    _emit(f"  {ri_src}")
    if len(ri_zone) >= 1:
        for k, lz in enumerate(ri_zone):
            r_um = zb[lz + 1] * 1e4
            _emit(f"    after zone {lz:3d} (r = {r_um:7.2f} µm): "
                  f"{region_names[k]} → {region_names[k+1] if k+1 < len(region_names) else '...'}")
    _emit("")

    # Identify ice/foam interface = end of region 2 (DT ice) in a 4-region target
    ice_foam_zone = ri_zone[1] if len(ri_zone) >= 2 else None
    if ice_foam_zone is not None:
        ifz = ice_foam_zone
        win = args.ice_foam_half_window
        z0 = max(0, ifz - win)
        z1 = min(n_z - 1, ifz + win)
        _emit(f"Ice/foam interface at zone {ifz} (radius "
              f"{zb[ifz+1]*1e4:.2f} µm at this timestep)")
        _emit(f"Printing zones [{z0}..{z1}] (window ±{win}):")
        _emit("")
        _emit(f"  {'zone':>4} {'r_um':>8} {'rho_gcc':>8} {'Tion_keV':>9} "
              f"{'Te_keV':>8} {'Trad_eV':>9} {'Zbar':>6} {'alpha_Wpcc':>12} "
              f"{'region':>12}")
        _emit("  " + "-" * 80)
        for k in range(z0, z1 + 1):
            reg = label_region(k, ri_zone, region_names)
            marker = '  <-- ICE/FOAM' if k == ifz else (
                     '  (foam side)'  if k == ifz + 1 else '')
            _emit(f"  {k:>4d} {zc[k]*1e4:>8.2f} {rho_t[k]:>8.3f} "
                  f"{Ti_keV[k]:>9.3f} {Te_keV[k]:>8.3f} {Trad_eV[k]:>9.2f} "
                  f"{Zbar_t[k]:>6.3f} {alpha_t[k]:>12.3e} {reg:>12s}{marker}")

        # Quick ratio diagnostic across the boundary
        if ifz + 1 < n_z and ifz - 0 >= 0:
            def _ratio(a, b):
                if not np.isfinite(a) or not np.isfinite(b) or b == 0:
                    return float('nan')
                return float(a / b)
            r_Te    = _ratio(Te_keV[ifz],   Te_keV[ifz + 1])
            r_Trad  = _ratio(Trad_eV[ifz],  Trad_eV[ifz + 1])
            r_Tion  = _ratio(Ti_keV[ifz],   Ti_keV[ifz + 1])
            r_rho   = _ratio(rho_t[ifz],    rho_t[ifz + 1])
            r_alpha = _ratio(alpha_t[ifz],  alpha_t[ifz + 1])
            r_Zbar  = _ratio(Zbar_t[ifz + 1], Zbar_t[ifz])  # foam/ice
            _emit("")
            _emit(f"Ratios across ice/foam boundary (ice / foam, except Z̄):")
            _emit(f"  T_ion        : {r_Tion:.3f}")
            _emit(f"  T_e          : {r_Te:.3f}    "
                  f"  <-- mechanism (1) trigger if ≫ 1")
            _emit(f"  T_rad        : {r_Trad:.3f}")
            _emit(f"  ρ            : {r_rho:.3f}")
            _emit(f"  α-dep rate   : {r_alpha:.3f}    "
                  f"  <-- mechanism (2) trigger if ≫ 1")
            _emit(f"  Z̄ (foam/ice): {r_Zbar:.3f}    "
                  f"  <-- magnitude of the discontinuity")
            _emit("")
            _emit("Interpretation cheat-sheet:")
            _emit("  T_e ratio  ≥ 2-5 and T_rad continuous  →  bremsstrahlung sink at Z̄ jump")
            _emit("  T_e continuous, α-dep ratio ≥ 2-3      →  non-local α-transport range-out")
            _emit("  both continuous                        →  not (1) or (2); look at EOS / mix / timing")
    else:
        _emit("(No ice/foam interface identified — single-region or "
              "non-Olson target; CSV still written.)")

    _emit("")
    _emit("=" * 76)

    # ── Write CSV
    with open(csv_path, 'w') as f:
        f.write(f"# Burn-propagation diagnostic profile\n")
        f.write(f"# Run: {base.name}\n")
        f.write(f"# t_ns = {t_at_idx:.6f}  (idx {t_idx} / {n_t})\n")
        f.write(f"# time-mode: {t_mode_used}\n")
        f.write(f"# region interfaces (last-zone): {list(ri_zone)}\n")
        f.write(f"# region names: {region_names}\n")
        f.write(f"# alpha source: {alpha_src}\n")
        f.write("zone_idx,r_left_um,r_center_um,r_right_um,dr_um,"
                "rho_gcc,T_ion_keV,T_e_keV,T_rad_eV,Zbar,"
                "alpha_dep_total_Wpcc,alpha_dep_ion_Wpcc,alpha_dep_ele_Wpcc,"
                "region\n")
        for k in range(n_z):
            reg = label_region(k, ri_zone, region_names)
            f.write(f"{k},"
                    f"{zb[k]*1e4:.5f},{zc[k]*1e4:.5f},{zb[k+1]*1e4:.5f},"
                    f"{dr[k]*1e4:.5f},"
                    f"{rho_t[k]:.5e},{Ti_keV[k]:.5e},{Te_keV[k]:.5e},"
                    f"{Trad_eV[k]:.5e},{Zbar_t[k]:.5e},"
                    f"{alpha_t[k]:.5e},{alpha_i[k]:.5e},{alpha_e[k]:.5e},"
                    f"{reg}\n")

    # ── Mirror summary to .txt for paste-back
    with open(txt_path, 'w') as f:
        f.write("\n".join(summary) + "\n")

    print(f"Wrote {csv_path}  ({n_z} zones)")
    print(f"Wrote {txt_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
