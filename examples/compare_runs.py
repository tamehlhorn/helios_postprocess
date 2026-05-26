#!/usr/bin/env python3
"""
Compare two Helios runs side-by-side: R-T trajectories, time histories,
stagnation lineouts.  Loads each run via the standard `build_run_data`
+ `ICFAnalyzer` pipeline so all derived attributes (ablation-front
indices, hot-spot \u03c1R history, alpha-onset time, etc.) are populated.

Usage
-----
    python compare_runs.py <base_path_A> <base_path_B> [--labels A B] \
                           [--outdir comparisons]

Each `base_path` is the simulation base path WITHOUT extension (same
convention as run_analysis.py).  Outputs three PNG figures under
`<outdir>/`:

    compare_<A>_vs_<B>_rt.png         -- R-T side-by-side
    compare_<A>_vs_<B>_histories.png  -- ablation / implosion / burn / \u03b1
    compare_<A>_vs_<B>_lineouts.png   -- bang-time \u03c1(r) and T_ion(r)
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from helios_postprocess import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.icf_analysis import ICFAnalyzer

# Symbols (escape form so the source stays ASCII -- avoids paste/encoding issues)
RHO    = '\u03c1'   # rho
MU     = '\u00b5'   # micro
SUP2   = '\u00b2'   # superscript 2
ALPHA  = '\u03b1'   # alpha

# Region palette: copied from icf_plotting._plot_zone_trajectories so the
# R-T figure matches the per-run report.
REGION_COLORS = ['#00CC00', '#0066FF', '#00BBDD', '#FF2200',
                 '#FF8800', '#AA00FF', '#888888', '#FFCC00']

# Per-run overlay colors (used in histories + lineouts).
RUN_COLORS = ('#1f77b4', '#d62728')      # A blue, B red
RUN_STYLES = ('-', '--')                 # A solid, B dashed


# -- Data loading -----------------------------------------------------------

def load_run(base_path: Path, verbose: bool = False):
    """Load EXODUS + RHW for one run; return populated ICFRunData."""
    exo_path = base_path.with_suffix('.exo')
    rhw_path = base_path.with_suffix('.rhw')
    if not exo_path.exists():
        raise FileNotFoundError(f"EXODUS file not found: {exo_path}")

    rhw_config = drive_T = drive_t = None
    if rhw_path.exists():
        try:
            from helios_postprocess.rhw_parser import load_rhw_configuration
            rhw_config = load_rhw_configuration(rhw_path)
            drive_T = rhw_config.drive_temperature
            drive_t = rhw_config.drive_time
        except Exception as e:
            print(f"  WARNING: Could not parse RHW for {base_path.name}: {e}")

    run = HeliosRun(str(exo_path), verbose=verbose)
    data = build_run_data(
        run, time_unit='s',
        rhw_config=rhw_config,
        drive_temperature=drive_T,
        drive_time=drive_t,
    )
    run.close()

    analyzer = ICFAnalyzer(data)
    analyzer.analyze_drive_phase()
    analyzer.analyze_laser_intensity()
    analyzer.analyze_stagnation_phase()
    analyzer.analyze_burn_phase()
    analyzer.compute_performance_metrics()
    return data


# -- Derived per-time series (not on data directly) -------------------------

def unablated_fuel_fraction_vs_time(data) -> np.ndarray:
    """
    Generalize ICFAnalyzer's scalar 'unablated_fuel_mass' to a (n_times,)
    time series.  Uses the same zone-index recipe as
    compute_performance_metrics (line ~2278 onward).
    """
    ri = data.region_interfaces_indices
    zmass = data.zone_mass
    n_times, n_zones = zmass.shape
    fuel_bnd = int(ri[0, getattr(data, 'fuel_ablator_idx', -2)])
    initial_fuel_mass = float(np.sum(zmass[0, :fuel_bnd]))
    if initial_fuel_mass <= 0:
        return np.full(n_times, np.nan)

    abl_indices = data.ablation_front_indices
    rho = data.mass_density
    out = np.zeros(n_times)
    zone_idx_arr = np.arange(n_zones)

    for t in range(n_times):
        hs_bnd = int(ri[t, 0])
        abl_idx = None
        if abl_indices is not None and abl_indices[t] > 0:
            abl_idx = int(abl_indices[t])
        if abl_idx is None or abl_idx <= hs_bnd:
            # density fallback (matches analyzer)
            abl_idx = hs_bnd
            for z in range(n_zones - 1, hs_bnd, -1):
                if rho[t, z] > 1.0:
                    abl_idx = z
                    break
        inside = zone_idx_arr <= abl_idx
        fuel_mask = zone_idx_arr < fuel_bnd
        out[t] = np.sum(zmass[t, fuel_mask & inside]) / initial_fuel_mass
    return out


def convergence_ratio_vs_time(data) -> np.ndarray:
    """CR(t) = R0 / r_hs(t).  R0 from t=0 hot-spot boundary node."""
    ri = data.region_interfaces_indices
    zbnd = data.zone_boundaries
    n_times = zbnd.shape[0]
    hs_node_0 = int(ri[0, 0])
    R0 = zbnd[0, hs_node_0]
    out = np.full(n_times, np.nan)
    for t in range(n_times):
        hs_node = int(ri[t, 0])
        if hs_node <= 0:
            continue
        r_hs = zbnd[t, hs_node]
        if r_hs > 0:
            out[t] = R0 / r_hs
    return out


def peak_inward_shell_velocity_kms_vs_time(data) -> np.ndarray:
    """
    |min(v_zone[t, shell])| in km/s, where shell zones are those with
    index >= hot-spot boundary at that timestep.  Mirrors the
    peak-velocity logic in analyze_implosion_phase.
    """
    ri = data.region_interfaces_indices
    vel = data.velocity
    n_zones = data.zone_mass.shape[1]
    if vel.shape[1] > n_zones:
        v_zone = 0.5 * (vel[:, :n_zones] + vel[:, 1:n_zones + 1])
    else:
        v_zone = vel[:, :n_zones]
    n_times = v_zone.shape[0]
    out = np.zeros(n_times)
    for t in range(n_times):
        hs_bnd = int(ri[t, 0])
        if hs_bnd >= n_zones:
            continue
        shell = v_zone[t, hs_bnd:]
        if shell.size == 0:
            continue
        out[t] = -min(shell.min(), 0.0) * 1e-5     # cm/s -> km/s
    return out


def fusion_power_TW_vs_time(data) -> np.ndarray:
    """
    data.fusion_power is FusionRate_DT_nHe4 (reactions/s per zone).
    Total fusion power = sum_z(rate) * 17.6 MeV * 1.602e-13 J/MeV.
    """
    if data.fusion_power is None:
        return np.zeros(len(data.time))
    rate = np.sum(data.fusion_power, axis=1)
    return rate * 17.6 * 1.602e-13 * 1e-12          # W -> TW


def alpha_heating_sum_vs_time(data) -> Optional[np.ndarray]:
    """Sum over zones of alpha_heating_ion + alpha_heating_ele.

    Native units follow Helios's pt_particle_heating_{ion,ele} -- used here
    in raw form to match how ICFAnalyzer.analyze_implosion_phase consumes
    the same quantity (line ~396).  Magnitudes between runs are directly
    comparable; absolute power scale not asserted.
    """
    ion = getattr(data, 'alpha_heating_ion', None)
    ele = getattr(data, 'alpha_heating_ele', None)
    if ion is None or ele is None:
        return None
    return (ion + ele).sum(axis=1)


# -- Plot helpers -----------------------------------------------------------

def _zone_region_assignment(data) -> Tuple[np.ndarray, list, list]:
    """Lagrangian zone -> region index (using t=0 interfaces)."""
    ri0 = data.region_interfaces_indices[0].astype(int)
    n_zones = data.zone_mass.shape[1]
    zone_region = np.zeros(n_zones, dtype=int)
    prev = 0
    for reg_idx, bnd in enumerate(ri0):
        zone_region[prev:int(bnd)] = reg_idx
        prev = int(bnd)
    region_names = data.region_names or [f'Region {i + 1}' for i in range(len(ri0))]
    colors = REGION_COLORS[:len(ri0)]
    return zone_region, list(region_names), colors


def _plot_rt_panel(ax, data, title: str, ylim_um: float):
    """Render one R-T trajectory panel (mimics icf_plotting._plot_zone_trajectories)."""
    time = data.time
    boundaries = data.zone_boundaries
    n_times, n_nodes = boundaries.shape
    n_zones = n_nodes - 1

    zone_region, region_names, colors = _zone_region_assignment(data)

    target_lines = 150
    stride = max(1, n_zones // target_lines)
    ri0 = data.region_interfaces_indices[0].astype(int)
    boundary_nodes = {0, n_nodes - 1}
    boundary_nodes.update(int(b) for b in ri0)

    bnd_um = boundaries * 1e4
    legend_added = set()
    for z in range(n_zones):
        node = z + 1
        if node not in boundary_nodes and z % stride != 0:
            continue
        reg = zone_region[z]
        lbl = region_names[reg] if reg not in legend_added else None
        legend_added.add(reg)
        ax.plot(time, bnd_um[:, node], color=colors[reg],
                linewidth=0.4, alpha=0.85, label=lbl)
    ax.plot(time, bnd_um[:, 0], color=colors[0], linewidth=0.4, alpha=0.85)

    if data.stag_time > 0:
        ax.axvline(data.stag_time, color='black', linestyle='--',
                   linewidth=2, alpha=0.9,
                   label=f'Stagnation ({data.stag_time:.2f} ns)', zorder=10)
    if data.bang_time > 0:
        ax.axvline(data.bang_time, color='darkred', linestyle='--',
                   linewidth=2, alpha=0.9,
                   label=f'Bang time ({data.bang_time:.2f} ns)', zorder=10)

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(f'Radius ({MU}m)')
    ax.set_title(title, fontweight='bold')
    ax.set_xlim(time[0], time[-1])
    ax.set_ylim(0, ylim_um)
    ax.grid(True, alpha=0.15)
    ax.legend(fontsize=8, loc='upper right')


def figure_rt(dA, dB, labelA, labelB, out_path: Path):
    """R-T side-by-side, run A left, run B right, shared y range."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)
    # Use a generous shared y limit large enough for both runs' early corona.
    max_um = max(
        np.nanpercentile(dA.zone_boundaries[:, -1] * 1e4, 99),
        np.nanpercentile(dB.zone_boundaries[:, -1] * 1e4, 99),
    )
    ylim = float(min(np.ceil(max_um / 500.0) * 500.0, 5000.0))

    _plot_rt_panel(axes[0], dA, labelA, ylim)
    _plot_rt_panel(axes[1], dB, labelB, ylim)
    fig.suptitle(f'R-T trajectories: {labelA} vs {labelB}',
                 fontweight='bold', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)


def figure_histories(dA, dB, labelA, labelB, out_path: Path):
    """Four time-history panels overlaid."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    runs = [(dA, labelA, RUN_COLORS[0], RUN_STYLES[0]),
            (dB, labelB, RUN_COLORS[1], RUN_STYLES[1])]

    # --- Ablation panel: unablated DT fuel fraction + ablation-front radius
    ax = axes[0, 0]
    ax2 = ax.twinx()
    for d, lab, col, ls in runs:
        frac = unablated_fuel_fraction_vs_time(d)
        ax.plot(d.time, frac, color=col, linestyle=ls, linewidth=2,
                label=f'{lab}: unablated DT')
        if d.ablation_front_radius is not None:
            afr = d.ablation_front_radius * 1e4   # cm -> um
            ax2.plot(d.time, afr, color=col, linestyle=':',
                     linewidth=1.5, alpha=0.7,
                     label=f'{lab}: abl-front r')
    ax.set_ylabel('Unablated DT fraction')
    ax.set_ylim(0, 1.05)
    ax2.set_ylabel(f'Ablation-front radius ({MU}m)', color='gray')
    ax2.tick_params(axis='y', colors='gray')
    ax.set_title('Ablation')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper left')
    ax2.legend(fontsize=7, loc='upper right')

    # --- Implosion panel: CR(t) + peak inward shell velocity
    ax = axes[0, 1]
    ax2 = ax.twinx()
    for d, lab, col, ls in runs:
        cr = convergence_ratio_vs_time(d)
        ax.plot(d.time, cr, color=col, linestyle=ls, linewidth=2,
                label=f'{lab}: CR')
        vk = peak_inward_shell_velocity_kms_vs_time(d)
        ax2.plot(d.time, vk, color=col, linestyle=':',
                 linewidth=1.5, alpha=0.7,
                 label=f'{lab}: |v_inward|')
    ax.set_ylabel('Convergence ratio (R$_0$/r$_{hs}$)')
    ax.set_yscale('log')
    ax2.set_ylabel('Peak inward shell |v| (km/s)', color='gray')
    ax2.tick_params(axis='y', colors='gray')
    ax.set_title('Implosion')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=8, loc='upper left')
    ax2.legend(fontsize=7, loc='upper right')

    # --- Burn panel: total & HS rhoR + total fusion power
    ax = axes[1, 0]
    ax2 = ax.twinx()
    for d, lab, col, ls in runs:
        if d.total_rhoR_vs_time is not None:
            ax.plot(d.time, d.total_rhoR_vs_time, color=col, linestyle=ls,
                    linewidth=2, label=f'{lab}: total {RHO}R')
        if d.hot_spot_rhoR_vs_time is not None:
            ax.plot(d.time, d.hot_spot_rhoR_vs_time, color=col, linestyle=ls,
                    linewidth=1.2, alpha=0.6,
                    label=f'{lab}: HS {RHO}R')
        fp = fusion_power_TW_vs_time(d)
        if np.any(fp > 0):
            ax2.semilogy(d.time, np.clip(fp, 1e-6, None), color=col,
                         linestyle=':', linewidth=1.5, alpha=0.7,
                         label=f'{lab}: fusion P')
    ax.set_ylabel(f'{RHO}R (g/cm{SUP2})')
    ax2.set_ylabel('Fusion power (TW)', color='gray')
    ax2.tick_params(axis='y', colors='gray')
    ax.set_xlabel('Time (ns)')
    ax.set_title('Burn')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper left')
    ax2.legend(fontsize=7, loc='upper right')

    # --- Alpha-heating panel
    ax = axes[1, 1]
    for d, lab, col, ls in runs:
        ah = alpha_heating_sum_vs_time(d)
        if ah is not None and np.any(ah > 0):
            ax.semilogy(d.time, np.clip(ah, 1e-3, None),
                        color=col, linestyle=ls, linewidth=2,
                        label=f'{lab}: {ALPHA} heat')
        onset = getattr(d, 'alpha_onset_time_ns', None)
        if onset is not None and onset > 0:
            ax.axvline(onset, color=col, linestyle=':', alpha=0.6,
                       label=f'{lab}: {ALPHA} onset {onset:.2f} ns')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(f'{ALPHA} heating (sum over zones)')
    ax.set_title(f'{ALPHA} heating')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=8, loc='upper left')

    fig.suptitle(f'Time histories: {labelA} vs {labelB}',
                 fontweight='bold', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)


def _rhoR_03_radius(data, t_idx: int) -> float:
    """Radius (um) where cumulative \u03c1R from r=0 crosses 0.3 g/cm\u00b2,
    evaluated at timestep t_idx.  Returns NaN if it never crosses."""
    cumul = data.areal_density_vs_time
    if cumul is None:
        return float('nan')
    zbnd = data.zone_boundaries
    crho = cumul[t_idx]                 # (n_zones+1,)
    nodes = zbnd[t_idx]                 # (n_zones+1,) cm
    over = np.where(crho >= 0.3)[0]
    if len(over) == 0:
        return float('nan')
    return float(nodes[over[0]] * 1e4)  # cm -> um


def figure_lineouts(dA, dB, labelA, labelB, out_path: Path):
    """Bang-time \u03c1(r), T_ion(r) overlaid."""
    fig, axes = plt.subplots(2, 1, figsize=(10, 8.5), sharex=True)
    runs = [(dA, labelA, RUN_COLORS[0], RUN_STYLES[0]),
            (dB, labelB, RUN_COLORS[1], RUN_STYLES[1])]

    for d, lab, col, ls in runs:
        if d.bang_time <= 0:
            continue
        bang_idx = int(np.argmin(np.abs(d.time - d.bang_time)))
        zbnd = d.zone_boundaries[bang_idx]
        zc_cm = 0.5 * (zbnd[:-1] + zbnd[1:])
        zc_um = zc_cm * 1e4

        # rho
        axes[0].plot(zc_um, d.mass_density[bang_idx], color=col,
                     linestyle=ls, linewidth=2,
                     label=f'{lab} (t={d.bang_time:.2f} ns)')
        # T_ion in keV
        T_keV = d.ion_temperature[bang_idx] / 1000.0
        axes[1].plot(zc_um, T_keV, color=col,
                     linestyle=ls, linewidth=2,
                     label=f'{lab} (t={d.bang_time:.2f} ns)')

        # Markers per run: HS radius (solid thin) and rhoR=0.3 (dot-dash thin)
        hs_node = int(d.region_interfaces_indices[bang_idx, 0])
        r_hs_um = float(d.zone_boundaries[bang_idx, hs_node] * 1e4)
        r03_um = _rhoR_03_radius(d, bang_idx)
        for ax in axes:
            ax.axvline(r_hs_um, color=col, linestyle=ls, alpha=0.35, linewidth=1)
            if not np.isnan(r03_um):
                ax.axvline(r03_um, color=col, linestyle='-.', alpha=0.35,
                           linewidth=1)

    axes[0].set_ylabel(f'{RHO} (g/cm{SUP2})'.replace('cm\u00b2', 'cc'))
    axes[0].set_yscale('log')
    axes[0].set_ylim(1e-3, None)
    axes[0].set_title(f'Stagnation lineouts (bang time): {labelA} vs {labelB}',
                      fontweight='bold')
    axes[0].grid(True, alpha=0.3, which='both')
    axes[0].legend(fontsize=9)

    axes[1].set_xlabel(f'Radius ({MU}m)')
    axes[1].set_ylabel('T$_{ion}$ (keV)')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=9)

    # Restrict x range to interesting near-core region
    axes[1].set_xlim(0, 500)

    # Caption line explaining marker conventions
    fig.text(0.5, 0.005,
             f'Vertical markers (per run): solid/dashed = HS outer radius; '
             f'dash-dot = radius where cumulative {RHO}R = 0.3 g/cm{SUP2}',
             ha='center', fontsize=8, style='italic')

    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)


# -- Main -------------------------------------------------------------------

def _sanity_snapshot(d, label: str) -> str:
    """One-line summary used for end-of-run sanity check."""
    bang_idx = (int(np.argmin(np.abs(d.time - d.bang_time)))
                if d.bang_time > 0 else 0)
    rho_peak = (float(np.max(d.mass_density[bang_idx]))
                if d.mass_density is not None else float('nan'))
    hs_rhoR_peak = (float(np.nanmax(d.hot_spot_rhoR_vs_time))
                    if d.hot_spot_rhoR_vs_time is not None else float('nan'))
    unablated_final = float(getattr(d, 'unablated_fuel_mass', float('nan')))
    return (f"  [{label}] bang={d.bang_time:.2f} ns  "
            f"peak {RHO}@bang={rho_peak:.1f} g/cc  "
            f"max HS {RHO}R={hs_rhoR_peak:.3f} g/cm{SUP2}  "
            f"unablated DT (stag scalar)={unablated_final*100:.1f}%")


def main():
    parser = argparse.ArgumentParser(
        description='Compare two Helios runs side-by-side.')
    parser.add_argument('base_a', help='Base path of run A (no extension)')
    parser.add_argument('base_b', help='Base path of run B (no extension)')
    parser.add_argument('--labels', nargs=2, metavar=('LBL_A', 'LBL_B'),
                        default=None,
                        help='Short labels for the two runs (default: dir names)')
    parser.add_argument('--outdir', default='comparisons',
                        help='Output directory (relative paths resolved from CWD)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose HeliosRun loader')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(name)s: %(message)s')

    base_a = Path(args.base_a).expanduser().resolve()
    base_b = Path(args.base_b).expanduser().resolve()
    if args.labels:
        labelA, labelB = args.labels
    else:
        labelA, labelB = base_a.name, base_b.name

    outdir = Path(args.outdir).expanduser()
    if not outdir.is_absolute():
        outdir = (Path.cwd() / outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print('=' * 78)
    print(f'  Comparing:  A = {labelA}')
    print(f'              B = {labelB}')
    print(f'  Output dir: {outdir}')
    print('=' * 78)

    print(f'Loading A: {base_a}')
    dA = load_run(base_a, verbose=args.verbose)
    print(f'Loading B: {base_b}')
    dB = load_run(base_b, verbose=args.verbose)

    print()
    print('Sanity snapshot:')
    print(_sanity_snapshot(dA, labelA))
    print(_sanity_snapshot(dB, labelB))
    print()

    rt_png        = outdir / f'compare_{labelA}_vs_{labelB}_rt.png'
    histories_png = outdir / f'compare_{labelA}_vs_{labelB}_histories.png'
    lineouts_png  = outdir / f'compare_{labelA}_vs_{labelB}_lineouts.png'

    print(f'Writing {rt_png.name} ...')
    figure_rt(dA, dB, labelA, labelB, rt_png)
    print(f'Writing {histories_png.name} ...')
    figure_histories(dA, dB, labelA, labelB, histories_png)
    print(f'Writing {lineouts_png.name} ...')
    figure_lineouts(dA, dB, labelA, labelB, lineouts_png)

    print()
    print('Outputs:')
    for p in (rt_png, histories_png, lineouts_png):
        print(f'  {p}')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
