"""
check_flux_limiter.py -- Probe whether the file-reported flux limiter f
is the value the simulation is actually using, and confirm the direction
of the Prism ×4 convention factor.

What we can measure cleanly
---------------------------
Helios doesn't output the electron heat flux q_e directly, so we cannot
read the saturated flux off the EXODUS arrays. What we CAN compute:

  q_FS = α n_e kT_e v_th_e          (free-streaming limit, α≈0.6)
  q_SH = κ_SH(T_e) |dT_e/dr|        (Spitzer-Harm classical estimate)

The ratio q_SH/q_FS tells us *if* the limiter is saturating:
  q_SH/q_FS >> 1   →  limiter is fully active (q_actual = f × q_FS)
  q_SH/q_FS ≲ 1    →  limiter inactive (q_actual = q_SH, no cap)

For the *value* of f, the cleanest hardware-independent indicator is
the absolute corona temperature. At the critical surface in steady
state with the limiter saturating, an energy balance gives:

    T_e(r_crit) ∝ f^(-2/3)

so a run at f=0.015 should have T_crit roughly (0.06/0.015)^(2/3) ≈
2.5× higher than a run at f=0.06. If you pair two runs whose only
RHW difference is f, comparing their T_crit values resolves the
convention:

  T_crit ratio matches (f_file ratio)^(-2/3)
      →  file values are the actual physical f (no conversion)

  T_crit ratio matches (4·f_file)^(-2/3)
      →  Prism ×4 hypothesis is correct (standard = 4 × file)

  T_crit nearly identical between the runs
      →  both files are in the "essentially unlimited" regime; FL
         knob isn't engaging at all (suggests f_standard ≫ 0.1 for
         both, regardless of convention)

The q_e profile from P_outside(r)/(4π r²) is also plotted as an
*upper bound* — it's the heat flux that would flow inward if 100% of
the absorbed laser power became electron conduction with no corona
heating, radiation, or outward flow. Use only for "is the limiter
saturating" sanity, not for quantitative f.

Usage
-----
    python3 ~/helios_postprocess/examples/check_flux_limiter.py \\
        <base1> [<base2> ...] \\
        [--out check_flux_limiter.pdf]
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from helios_postprocess.core import HeliosRun
from helios_postprocess.data_builder import build_run_data


# Physical constants (SI)
EV          = 1.602e-19      # J per eV
M_E         = 9.109e-31      # kg, electron mass
ALPHA_FS    = 0.6            # Goncharov / ICF convention for q_FS prefactor


def free_streaming_flux_W_per_cm2(n_e_cm3: np.ndarray, T_e_eV: np.ndarray) -> np.ndarray:
    """q_FS = α n_e kT_e v_th_e, returned in W/cm²."""
    n_e_m3 = n_e_cm3 * 1e6                                  # cm^-3 -> m^-3
    kT_J   = T_e_eV * EV                                    # J
    v_th   = np.sqrt(np.maximum(kT_J, 0.0) / M_E)           # m/s
    q_W_m2 = ALPHA_FS * n_e_m3 * kT_J * v_th                # W/m²
    return q_W_m2 * 1e-4                                    # -> W/cm²


def analyze_run(base: Path) -> dict:
    """Return a diagnostic dict for one run."""
    run = HeliosRun(str(base) + ".exo")
    # Load the RHW so data.flux_limiter_per_region / data.flux_limiter get
    # populated. run_analysis.py does this in two steps; mirror that.
    rhw_config = None
    rhw_path = Path(str(base) + ".rhw")
    if rhw_path.exists():
        try:
            from helios_postprocess.rhw_parser import load_rhw_configuration
            rhw_config = load_rhw_configuration(rhw_path)
        except Exception as e:
            print(f"  WARNING: RHW parse failed for {base.name}: {e}",
                  file=sys.stderr)
    data = build_run_data(run, time_unit='s', rhw_config=rhw_config)
    run.close()

    # ---- time and laser deposition geometry ----
    t_raw = np.asarray(data.time)
    t_ns  = t_raw * 1e9 if float(np.max(t_raw)) < 1e-3 else t_raw

    P_total = np.asarray(data.laser_power_delivered)
    if P_total.ndim == 2:
        P_total = P_total.sum(axis=1)
    i_peak = int(np.argmax(P_total))

    # zone centers / boundaries at peak power
    zb       = np.asarray(data.zone_boundaries[i_peak])
    zone_c   = 0.5 * (zb[:-1] + zb[1:])
    n_zones  = len(zone_c)

    # ---- profile quantities ----
    T_e_eV = np.asarray(data.elec_temperature[i_peak])
    n_e    = np.asarray(data.electron_density[i_peak])
    rho    = np.asarray(data.mass_density[i_peak])

    # ---- Cumulative absorbed power outside each radius ----
    # LaserPwrSrc is W/cm³ per zone. The laser deposits OUTSIDE the
    # conduction zone (at / near the critical surface). The heat flowing
    # INWARD through a radius r in the conduction zone equals the total
    # power absorbed OUTSIDE that radius (steady-state spherical
    # conservation with no inner sources/sinks):
    #     q_e(r) = [ P_absorbed at r' > r ] / (4 π r²)
    P_src = np.asarray(data.laser_power_source[i_peak])
    dV    = (4.0 / 3.0) * np.pi * (zb[1:]**3 - zb[:-1]**3)   # cm³ per zone
    P_per_zone = P_src * dV                                  # W per zone
    P_inside   = np.cumsum(P_per_zone)                       # W absorbed inside r
    P_outside  = P_inside[-1] - P_inside                     # W absorbed outside r

    # ---- implied inward heat flux ----
    safe_r = np.maximum(zone_c, 1e-8)
    q_e    = P_outside / (4.0 * np.pi * safe_r**2)           # W/cm²

    # ---- free-streaming limit ----
    q_FS   = free_streaming_flux_W_per_cm2(n_e, T_e_eV)      # W/cm²

    # ---- Spitzer-Harm classical heat flux from local gradient ----
    # κ_SH = 9.4e-13 × T_e^(5/2) / (Z·lnΛ)  [W/(cm·eV)]  (Atzeni 4.27)
    # dT_e/dr in eV/cm from np.gradient on zone centres.
    Z, ln_Lambda = 1.0, 10.0
    kappa_SH    = 9.4e-13 * np.maximum(T_e_eV, 0.0)**2.5 / (Z * ln_Lambda)
    dT_dr       = np.gradient(T_e_eV, zone_c)                # eV/cm
    q_SH        = kappa_SH * np.abs(dT_dr)                   # W/cm²

    # ---- conduction zone: between ablation front and critical surface ----
    ri = data.region_interfaces_indices
    hs_outer_node = int(ri[i_peak, 0]) if ri is not None else 0
    abl_idx       = int(data.ablation_front_indices[i_peak]) \
                    if data.ablation_front_indices is not None \
                    else max(hs_outer_node, n_zones // 2)
    abl_idx = max(abl_idx, hs_outer_node + 1)
    abl_idx = min(abl_idx, n_zones - 1)

    # Critical surface from electron density crossing
    # n_crit at 0.351 µm = 9.10e21 cm^-3
    n_crit = 9.10e21
    crossings = np.where(n_e >= n_crit)[0]
    if crossings.size > 0:
        crit_idx = int(np.max(crossings))    # outermost zone with ne >= n_crit
    else:
        crit_idx = max(abl_idx, n_zones - 5)

    # Conduction zone is [abl_idx + 1, crit_idx]. Bail out if degenerate.
    cz_lo = abl_idx + 1
    cz_hi = max(cz_lo + 2, crit_idx)
    cz_hi = min(cz_hi, n_zones - 1)
    cz_slice = slice(cz_lo, cz_hi + 1)

    # Effective flux limiter from the median of the ratio in the
    # conduction zone (robust against single-zone spikes / log-T noise).
    ratio = np.where(q_FS > 0, q_e / q_FS, np.nan)
    cz_ratio = ratio[cz_slice]
    if cz_ratio.size > 0 and np.any(np.isfinite(cz_ratio)):
        f_implied = float(np.nanmedian(cz_ratio))
        f_max     = float(np.nanmax(cz_ratio))
    else:
        f_implied = float('nan')
        f_max     = float('nan')

    # Pull the RHW-reported f (first non-zero region)
    fl_per = getattr(data, 'flux_limiter_per_region', None) or []
    f_file = None
    for entry in fl_per:
        v = entry.get('value', 0.0)
        if v > 0:
            f_file = float(v); break
    if f_file is None:
        f_file = float(getattr(data, 'flux_limiter', 0.0)) or None

    # ---- saturation ratio q_SH / q_FS in conduction zone ----
    sat_ratio = np.where(q_FS > 0, q_SH / q_FS, np.nan)
    cz_sat = sat_ratio[cz_slice]
    sat_median = float(np.nanmedian(cz_sat)) if cz_sat.size else float('nan')

    return {
        'name':           base.name,
        't_ns_peak':      float(t_ns[i_peak]),
        'P_peak_TW':      float(P_total[i_peak]) * 1e-12,
        'zone_c':         zone_c,
        'T_e_eV':         T_e_eV,
        'n_e':            n_e,
        'rho':            rho,
        'q_e':            q_e,
        'q_FS':           q_FS,
        'q_SH':           q_SH,
        'ratio':          ratio,
        'sat_ratio':      sat_ratio,
        'sat_median_cz':  sat_median,
        'cz_lo':          cz_lo,
        'cz_hi':          cz_hi,
        'r_abl':          float(zone_c[abl_idx]),
        'r_crit':         float(zone_c[crit_idx]),
        'T_crit_keV':     float(T_e_eV[crit_idx]) * 1e-3,
        'T_abl_keV':      float(T_e_eV[abl_idx])  * 1e-3,
        'f_file':         f_file,
        'f_implied':      f_implied,
        'f_max':          f_max,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('bases', nargs='+', type=Path,
                    help='Run base paths (without .exo extension).')
    ap.add_argument('--out', type=Path,
                    default=Path('check_flux_limiter.pdf'))
    args = ap.parse_args(argv)

    results = []
    for b in args.bases:
        try:
            results.append(analyze_run(b))
        except Exception as e:
            print(f"ERROR loading {b}: {e}", file=sys.stderr)
            return 1

    # --- stdout summary ---
    print()
    print(f"  {'run':<46s}  {'f_file':>7s}  {'T_crit':>7s}  {'T_abl':>7s}  "
          f"{'q_SH/q_FS':>10s}  {'q_e/q_FS_max':>13s}")
    print('  ' + '-' * 105)
    for r in results:
        f_file_str = (f"{r['f_file']:>7.4f}"
                      if (r['f_file'] is not None and r['f_file'] > 0)
                      else "   —   ")
        print(f"  {r['name']:<46s}  "
              f"{f_file_str}  "
              f"{r['T_crit_keV']:>7.2f}  "
              f"{r['T_abl_keV']:>7.3f}  "
              f"{r['sat_median_cz']:>10.2e}  "
              f"{r['f_implied']:>13.2e}")
    print()
    print("Columns:")
    print("  T_crit, T_abl       : Te at critical surface / ablation front (keV)")
    print("  q_SH/q_FS           : median in conduction zone -- >>1 means")
    print("                        limiter IS saturating; ≲1 means inactive")
    print("  q_e/q_FS_max        : upper-bound implied ratio assuming all")
    print("                        absorbed power flows inward (saturation check only)")
    print()

    # --- Pairwise convention verdict (needs ≥2 runs) ---
    if len(results) >= 2:
        valid = [r for r in results
                 if r['f_file'] is not None and r['f_file'] > 0
                 and math.isfinite(r['T_crit_keV']) and r['T_crit_keV'] > 0]
        if len(valid) >= 2:
            # Pick the two extreme f_file values
            valid.sort(key=lambda r: r['f_file'])
            lo, hi = valid[0], valid[-1]
            f_ratio   = hi['f_file']      / lo['f_file']
            T_ratio   = lo['T_crit_keV']  / hi['T_crit_keV']  # T scales as f^(-2/3)
            T_ratio_expected_direct  = f_ratio ** (2.0 / 3.0)        # if file f IS actual
            T_ratio_expected_x4      = (f_ratio) ** (2.0 / 3.0)      # ×4 gives same ratio!
            print(f"Pairwise convention check:")
            print(f"  {lo['name']:<46s}  f_file={lo['f_file']:.4f}  T_crit={lo['T_crit_keV']:.2f} keV")
            print(f"  {hi['name']:<46s}  f_file={hi['f_file']:.4f}  T_crit={hi['T_crit_keV']:.2f} keV")
            print(f"  f ratio (hi/lo)               = {f_ratio:.2f}")
            print(f"  T_crit ratio (lo/hi)          = {T_ratio:.2f}")
            print(f"  expected ratio f^(-2/3)       = {T_ratio_expected_direct:.2f}")
            print()
            # The T_crit ratio depends on the RATIO of f values, which is
            # invariant under the ×4 convention. So this pair test can't
            # distinguish direction; it only confirms saturation.
            # The ABSOLUTE T_crit value is the convention indicator:
            #   ~30 keV  -> f_actual ≈ 0.015 (heavy limiting)
            #   ~12 keV  -> f_actual ≈ 0.06
            #   ~5 keV   -> f_actual ≈ 0.24 (light limiting)
            #   ~2 keV   -> essentially unlimited
            print("Absolute T_crit indicator (for each run):")
            for r in valid:
                T = r['T_crit_keV']
                if   T >= 20:  inferred = "≈ 0.015 (heavy limiting)"
                elif T >= 8:   inferred = "≈ 0.06  (standard ICF)"
                elif T >= 3:   inferred = "≈ 0.24  (light limiting)"
                else:          inferred = "≪ 1     (essentially unlimited)"
                print(f"  {r['name']:<46s}  T_crit={T:5.2f} keV  →  f_actual {inferred}")
            print()

    # --- figure: q_e / q_FS profile, T_e profile, conduction-zone band ---
    fig, (axR, axT) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True,
        gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.08},
    )
    colors = plt.cm.tab10(np.linspace(0, 1, max(10, len(results))))

    for i, r in enumerate(results):
        c = colors[i % len(colors)]
        r_um = r['zone_c'] * 1e4
        # Plot the SATURATION ratio q_SH/q_FS (the physically meaningful one).
        axR.plot(r_um, r['sat_ratio'], color=c, lw=1.6, ls='-',
                 label=f"{r['name']}  q_SH/q_FS")
        # And the q_e_max/q_FS upper bound, dashed.
        axR.plot(r_um, r['ratio'], color=c, lw=1.0, ls='--', alpha=0.5)
        # Conduction-zone shade
        axR.axvspan(r['zone_c'][r['cz_lo']] * 1e4,
                    r['zone_c'][r['cz_hi']] * 1e4,
                    color=c, alpha=0.05)
        # f_file horizontal reference (where saturated q would land)
        if r['f_file'] is not None and r['f_file'] > 0:
            axR.axhline(r['f_file'], color=c, lw=0.8, ls=':',
                        alpha=0.6)
        axT.plot(r_um, r['T_e_eV'] * 1e-3, color=c, lw=1.4)

    axR.set_ylabel(r'$q_{SH}/q_{FS}$ (solid),  $q_e^{max}/q_{FS}$ (dashed)',
                   fontsize=11)
    axR.set_yscale('log')
    axR.set_ylim(1e-3, 1e5)
    axR.axhline(1.0, color='k', lw=0.5, alpha=0.4)
    axR.grid(True, alpha=0.3, which='both')
    axR.legend(fontsize=8, loc='upper left')
    axR.set_title('Flux-limiter saturation indicator (solid q_SH/q_FS) '
                  'and upper bound (dashed). Dotted = f_file.',
                  fontsize=11)

    axT.set_xlabel('Radius (µm)', fontsize=12)
    axT.set_ylabel(r'$T_e$ (keV)', fontsize=12)
    axT.set_yscale('log')
    axT.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    plt.close(fig)
    print(f"\nSaved: {args.out}")
    print("\nReading the plot:")
    print("  - Each curve = q_e(r)/q_FS(r) at peak laser power")
    print("  - Shaded band = conduction zone (between ablation front and r_crit)")
    print("  - Dotted horizontal line = f reported in RHW file")
    print("  - If curve plateaus AT the dotted line in the shaded band -> file f is actual")
    print("  - If plateau is 4x ABOVE -> Prism ×4 convention is correct")
    print("  - If plateau is 4x BELOW -> convention is inverted")
    return 0


if __name__ == '__main__':
    sys.exit(main())
