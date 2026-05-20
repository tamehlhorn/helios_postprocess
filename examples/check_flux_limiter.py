"""
check_flux_limiter.py -- Diagnose the *effective* electron flux limiter
implied by a Helios run, and compare it against the value reported in the
RHW file.

Built to test the Prism team's hypothesis that the file-reported `f` is
off by a factor of 4 from the ICF community standard. The direction of
that factor (×4 or ÷4) is what this script confirms.

Method
------
Helios doesn't output the electron heat flux directly, but in steady
state inside the critical surface the inward-flowing heat is just the
absorbed laser energy:

    q_e(r) = [ integrated absorbed power inside r ] / (4 π r²)

In the conduction zone (between the ablation front and the critical
surface) this q_e is exactly what the flux limiter caps:

    q_e ≈ f_actual × q_FS

where q_FS is the free-streaming thermal flux (Goncharov / Atzeni
convention):

    q_FS = α × n_e × kT_e × v_th_e         α ≈ 0.6
    v_th_e = sqrt(kT_e / m_e)

So the ratio q_e / q_FS in the conduction zone gives **f_actual**
directly. Compare that to the f written in the RHW file:

    if  f_actual ≈ f_file            : no convention conversion needed
    if  f_actual ≈ 4 × f_file        : Prism ×4 hypothesis is correct
    if  f_actual ≈ f_file / 4        : convention is inverted

Usage
-----
    python3 ~/helios_postprocess/examples/check_flux_limiter.py \\
        <base1> [<base2> ...] \\
        [--out check_flux_limiter.pdf]

Each positional argument is a run base path (without .exo extension).
Comparing two runs that differ ONLY in flux limiter is the cleanest
test, but a single run can also be diagnosed.
"""

import argparse
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
    run  = HeliosRun(str(base) + ".exo")
    data = build_run_data(run)

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

    # ---- Cumulative absorbed power inside each radius ----
    # LaserPwrSrc is W/cm³ per zone. Multiply by zone volume (spherical
    # shell) and cumulate from r=0 outward.
    P_src = np.asarray(data.laser_power_source[i_peak])
    dV    = (4.0 / 3.0) * np.pi * (zb[1:]**3 - zb[:-1]**3)   # cm³ per zone
    P_per_zone = P_src * dV                                  # W per zone
    P_inside   = np.cumsum(P_per_zone)                       # W absorbed inside r

    # ---- implied inward heat flux ----
    # q_e(r) = P_absorbed(<r) / (4 π r²)
    safe_r = np.maximum(zone_c, 1e-8)
    q_e    = P_inside / (4.0 * np.pi * safe_r**2)            # W/cm²

    # ---- free-streaming limit ----
    q_FS   = free_streaming_flux_W_per_cm2(n_e, T_e_eV)      # W/cm²

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

    return {
        'name':       base.name,
        't_ns_peak':  float(t_ns[i_peak]),
        'P_peak_TW':  float(P_total[i_peak]) * 1e-12,
        'zone_c':     zone_c,
        'T_e_eV':     T_e_eV,
        'n_e':        n_e,
        'rho':        rho,
        'q_e':        q_e,
        'q_FS':       q_FS,
        'ratio':      ratio,
        'cz_lo':      cz_lo,
        'cz_hi':      cz_hi,
        'r_abl':      float(zone_c[abl_idx]),
        'r_crit':     float(zone_c[crit_idx]),
        'f_file':     f_file,
        'f_implied':  f_implied,
        'f_max':      f_max,
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
    print(f"{'run':<46s}  {'t_pk':>5s}  {'P_pk':>5s}  "
          f"{'r_abl':>6s}  {'r_crit':>6s}  "
          f"{'f_file':>7s}  {'f_implied':>10s}  {'ratio':>6s}  verdict")
    print('-' * 120)
    for r in results:
        if r['f_file'] is None or r['f_file'] <= 0:
            ratio_str = '   —  '
            verdict   = '(no f in RHW)'
        else:
            ratio = r['f_implied'] / r['f_file']
            ratio_str = f"{ratio:>6.2f}"
            if 0.7 <= ratio <= 1.5:
                verdict = 'NO CONVERSION   (f_file already actual)'
            elif 3.0 <= ratio <= 5.0:
                verdict = 'PRISM ×4 CORRECT (standard = 4·file)'
            elif 0.20 <= ratio <= 0.33:
                verdict = 'CONVENTION INVERTED (standard = file/4)'
            else:
                verdict = '(no clean match — inspect plot)'
        print(f"  {r['name']:<44s}  "
              f"{r['t_ns_peak']:>5.2f}  "
              f"{r['P_peak_TW']:>5.0f}  "
              f"{r['r_abl']*1e4:>6.1f}  "
              f"{r['r_crit']*1e4:>6.1f}  "
              f"{(r['f_file'] if r['f_file'] is not None else float('nan')):>7.4f}  "
              f"{r['f_implied']:>10.4f}  "
              f"{ratio_str}  {verdict}")

    # --- figure: q_e / q_FS profile, T_e profile, conduction-zone band ---
    fig, (axR, axT) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True,
        gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.08},
    )
    colors = plt.cm.tab10(np.linspace(0, 1, max(10, len(results))))

    for i, r in enumerate(results):
        c = colors[i % len(colors)]
        r_um = r['zone_c'] * 1e4
        axR.plot(r_um, r['ratio'], color=c, lw=1.6,
                 label=(f"{r['name']}  (f_file={r['f_file']:.4f},  "
                        f"f_impl={r['f_implied']:.4f})"
                        if r['f_file'] else
                        f"{r['name']}  (f_impl={r['f_implied']:.4f})"))
        # Conduction-zone bar
        axR.axvspan(r['zone_c'][r['cz_lo']] * 1e4,
                    r['zone_c'][r['cz_hi']] * 1e4,
                    color=c, alpha=0.05)
        # f_file horizontal reference
        if r['f_file'] is not None and r['f_file'] > 0:
            axR.axhline(r['f_file'], color=c, lw=0.8, ls=':',
                        alpha=0.6)
        axT.plot(r_um, r['T_e_eV'] * 1e-3, color=c, lw=1.4)

    axR.set_ylabel(r'$q_e / q_{FS}$  (effective f)', fontsize=12)
    axR.set_yscale('log')
    axR.set_ylim(1e-3, 10)
    axR.grid(True, alpha=0.3, which='both')
    axR.legend(fontsize=8, loc='upper left')
    axR.set_title('Effective flux limiter from q_e / q_FS profile  '
                  '(shaded = conduction zone, dotted = f from RHW)',
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
