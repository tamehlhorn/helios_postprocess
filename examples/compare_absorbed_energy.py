"""
compare_absorbed_energy.py -- Overlay cumulative absorbed laser energy vs time
for two or more Helios runs.

Built to diagnose where two runs with identical RHW power tables but
different beam geometry (cone, spot, focus) diverge in *effective* coupling.
The interesting curves:

  Top    : E_absorbed(t)  (MJ)            -- raw absorbed energy
  Middle : E_delivered(t)  vs t (MJ)      -- power-table integral
  Bottom : E_absorbed / E_delivered (%)   -- instantaneous coupling fraction

Usage:
    python3 ~/helios_postprocess/examples/compare_absorbed_energy.py \\
        ~/Sims/Xcimer/Olson_PDD/_Working_archive/PDD_TM_4nb/PDD_TM_4nb \\
        ~/Sims/Xcimer/Olson_PDD/Olson_PDD_20_fab02_foot25_s016_burn/Olson_PDD_20_fab02_foot25_s016_burn

Each positional argument is the run's base path (without extension); the
script reads <base>.exo.

Output: absorbed_energy_compare.pdf in CWD, or override with --out.
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from helios_postprocess.core import HeliosRun
from helios_postprocess.data_builder import build_run_data


def load_run(base: Path):
    """Return (label, t_ns, E_absorbed_MJ, E_delivered_MJ, P_delivered_W)."""
    run  = HeliosRun(str(base) + ".exo")
    data = build_run_data(run)
    # Normalise time to ns (build_run_data may leave it in seconds depending
    # on default time_unit branch).
    t_raw = np.asarray(data.time)
    t_ns  = t_raw * 1e9 if float(np.max(t_raw)) < 1e-3 else t_raw

    # Absorbed energy: EnLaserDepositedTimeIntg (J, cumulative).
    # May be (n_times,) or (n_times, n_zones) depending on the loader.
    e_abs = data.laser_energy_deposited
    if e_abs is None:
        raise RuntimeError(f"{base}: laser_energy_deposited not loaded")
    if e_abs.ndim == 2:
        e_abs = e_abs.sum(axis=1)
    E_abs_MJ = np.asarray(e_abs) * 1e-6           # J -> MJ

    # Delivered: prefer Helios's own cumulative integral
    # (LaserEnDeliveredTimeInt, loaded as laser_energy_delivered_cum). It
    # shares the EnLaserDepositedTimeIntg time grid so the absorbed/
    # delivered ratio is self-consistent at every step. Fallback to a
    # python-side trapezoid integral of laser_power_delivered if the new
    # variable isn't in this EXODUS file.
    E_del_cum = getattr(data, 'laser_energy_delivered_cum', None)
    if E_del_cum is not None:
        E_del = np.asarray(E_del_cum)
        if E_del.ndim == 2:
            E_del = E_del.sum(axis=1)
        E_del_MJ = E_del * 1e-6                                # J -> MJ
        P_del_arr = np.asarray(data.laser_power_delivered) \
                    if data.laser_power_delivered is not None \
                    else np.zeros_like(t_ns)
    else:
        P_del = data.laser_power_delivered
        if P_del is None:
            P_del_arr = np.zeros_like(t_ns)
            E_del_MJ  = np.zeros_like(t_ns)
        else:
            P_del_arr = np.asarray(P_del)
            if P_del_arr.ndim == 2:
                P_del_arr = P_del_arr.sum(axis=1)
            # cumtrapz: W*ns = 1e-9 J; * 1e-6 -> MJ; total factor 1e-15
            E_del = np.zeros_like(t_ns)
            if len(t_ns) > 1:
                E_del[1:] = np.cumsum(0.5 * (P_del_arr[1:] + P_del_arr[:-1])
                                      * (t_ns[1:] - t_ns[:-1])) * 1e-9
            E_del_MJ = E_del * 1e-6

    return base.name, t_ns, E_abs_MJ, E_del_MJ, P_del_arr


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('bases', nargs='+', type=Path,
                    help='Run base paths (without .exo extension).')
    ap.add_argument('--out', type=Path, default=Path('absorbed_energy_compare.pdf'),
                    help='Output PDF path.')
    ap.add_argument('--xlim', type=float, nargs=2, metavar=('T_LO', 'T_HI'),
                    help='Manual x-axis limits in ns.')
    args = ap.parse_args(argv)

    if len(args.bases) < 2:
        print("Need at least 2 run base paths to compare.", file=sys.stderr)
        return 1

    runs = []
    for base in args.bases:
        try:
            runs.append(load_run(base))
        except Exception as e:
            print(f"ERROR loading {base}: {e}", file=sys.stderr)
            return 1

    # Print divergence summary
    print()
    print(f"{'run':<46s}  {'t_end':>6s}  {'E_del':>7s}  {'E_abs':>7s}  "
          f"{'eff_avg':>8s}  {'eff_peak':>9s}")
    print('-' * 90)
    for label, t_ns, E_abs, E_del, _ in runs:
        E_abs_final = float(E_abs[-1])
        E_del_final = float(E_del[-1])
        eff_avg = 100.0 * E_abs_final / E_del_final if E_del_final > 0 else 0.0
        # Peak instantaneous coupling fraction, after a 0.5 ns warmup
        mask = t_ns > 0.5
        if np.any(mask) and np.any(E_del[mask] > 1e-3):
            inst = np.where(E_del > 1e-3, E_abs / np.maximum(E_del, 1e-12), 0.0)
            eff_peak = 100.0 * float(np.max(inst[mask]))
        else:
            eff_peak = 0.0
        print(f"{label:<46s}  {t_ns[-1]:>6.2f}  {E_del_final:>7.3f}  "
              f"{E_abs_final:>7.3f}  {eff_avg:>7.1f}%  {eff_peak:>8.1f}%")

    # --- Figure ---
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(9, 9), sharex=True,
        gridspec_kw={'height_ratios': [2, 2, 1.5], 'hspace': 0.08},
    )
    colors = plt.cm.tab10(np.linspace(0, 1, max(10, len(runs))))

    for i, (label, t_ns, E_abs, E_del, _) in enumerate(runs):
        c = colors[i % len(colors)]
        ax1.plot(t_ns, E_abs, color=c, lw=1.8, label=label)
        ax2.plot(t_ns, E_del, color=c, lw=1.8, ls='--')
        # Avoid divide-by-zero before laser turn-on
        eff = np.where(E_del > 1e-3, 100.0 * E_abs / np.maximum(E_del, 1e-12), np.nan)
        ax3.plot(t_ns, eff, color=c, lw=1.5)

    ax1.set_ylabel('E_absorbed (MJ)', fontsize=12)
    ax1.set_title('Absorbed laser energy: cumulative comparison', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8, loc='upper left')

    ax2.set_ylabel('E_delivered (MJ)\n(dashed)', fontsize=11)
    ax2.grid(True, alpha=0.3)

    ax3.set_ylabel('absorbed /\ndelivered (%)', fontsize=11)
    ax3.set_xlabel('Time (ns)', fontsize=12)
    ax3.set_ylim(0, 110)
    ax3.axhline(100.0, color='k', lw=0.6, alpha=0.4)
    ax3.grid(True, alpha=0.3)

    # Sensible x-axis: from t=0 to max time across runs, or user override
    if args.xlim is not None:
        for ax in (ax1, ax2, ax3):
            ax.set_xlim(args.xlim[0], args.xlim[1])
    else:
        t_max = max(float(r[1][-1]) for r in runs)
        for ax in (ax1, ax2, ax3):
            ax.set_xlim(0, t_max)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    plt.close(fig)
    print(f"\nSaved: {args.out}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
