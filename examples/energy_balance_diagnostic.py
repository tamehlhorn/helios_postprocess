"""
energy_balance_diagnostic.py
============================

Side-by-side energy-flow comparison between two Helios runs.

Tracks where absorbed laser energy ends up across the available channels:

    KE_inward      : kinetic energy of zones moving inward (the implosion)
    KE_outward     : kinetic energy of zones moving outward (ablation blowoff)
    U_plasma       : plasma thermal energy, ideal-gas approx (3/2)(P_i + P_e) V
    U_rad          : radiation thermal energy, 3 P_rad V
    E_fusion       : cumulative fusion energy released (DT reactions × 17.6 MeV)

The four "where" channels should approximately sum to the cumulative
absorbed laser energy.  The remainder is reported as `gap` and is
dominated by:

    1. Radiation leaving the outer boundary (not tracked in zone-resident
       quantities, since zones-outside-grid carry energy that exits).
    2. Difference between the ideal-gas plasma thermal model used here and
       the true EOS-internal energy (ionization, electronic excitation,
       molecular binding).  For high-T plasmas this is small but nonzero.

Usage
-----
    python energy_balance_diagnostic.py <run1_basepath> <run2_basepath>

Example
-------
    python energy_balance_diagnostic.py \
        ~/Sims/Xcimer/Vulcan/VI_6/VI_6 \
        ~/Sims/Xcimer/HDD_26/HDD26_DTI40_1ns130_FL02_lrm1_burn/HDD26_DTI40_1ns130_FL02_lrm1_burn

Produces
--------
    <run1>_vs_<run2>_energy_balance.pdf  (3-page PDF in cwd)
    Console table with snapshots at peak-velocity, stagnation, end-of-run.

Notes
-----
This script is intentionally standalone — it uses the existing
helios_postprocess library to load runs but does not modify or wire
into the main analysis pipeline.  Run it independently on any pair of
existing .exo files.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from helios_postprocess import HeliosRun
from helios_postprocess.data_builder import build_run_data

logger = logging.getLogger("energy_balance")
logging.basicConfig(level=logging.INFO,
                    format='%(name)s: %(message)s')

# DT fusion reaction energy
Q_DT_J = 17.6e6 * 1.602e-19   # 17.6 MeV per reaction in J  (~2.82e-12)
# Q_DT split: 3.5 MeV alpha (deposited) + 14.1 MeV neutron (escapes 1D target)
Q_ALPHA_FRAC   = 3.5 / 17.6   # ~0.199
Q_NEUTRON_FRAC = 14.1 / 17.6  # ~0.801


# ── Energy ledger ────────────────────────────────────────────────────────

def compute_energy_ledger(data) -> dict:
    """
    Compute time-resolved energy ledger from an ICFRunData.

    All output quantities in Joules unless noted.  Each is shape (n_times,).
    """
    rho   = data.mass_density          # (T, Z) g/cm^3
    mass  = data.zone_mass             # (T, Z) g
    vfac  = data.velocity              # (T, Z+1) cm/s
    rb    = data.zone_boundaries       # (T, Z+1) cm
    P_i   = data.ion_pressure          # (T, Z) erg/cm^3
    P_e   = data.elec_pressure         # (T, Z) erg/cm^3
    t_ns  = data.time                  # (T,) ns

    n_t, n_z = rho.shape

    # Zone-centered velocity
    if vfac.shape[1] > n_z:
        v_zone = 0.5 * (vfac[:, :n_z] + vfac[:, 1:n_z + 1])
    else:
        v_zone = vfac[:, :n_z]

    # Spherical shell volumes (T, Z), cm^3
    V_zone = (4.0 / 3.0) * np.pi * (rb[:, 1:n_z + 1] ** 3 - rb[:, :n_z] ** 3)

    # Per-zone KE [erg], split by sign of velocity
    ke_per_zone = 0.5 * mass * v_zone ** 2
    inward  = v_zone < 0
    outward = v_zone > 0

    ke_inward  = np.where(inward,  ke_per_zone, 0.0).sum(axis=1) * 1e-7   # erg -> J
    ke_outward = np.where(outward, ke_per_zone, 0.0).sum(axis=1) * 1e-7

    # Plasma thermal (ideal gas): U = (3/2)(P_i + P_e) V  
    # Helios pressures are J/cm³, volumes cm³ → product is J directly.
    u_plasma = ((3.0 / 2.0) * (P_i + P_e) * V_zone).sum(axis=1)

    # Radiation: U_rad = 3 P_rad V (P_rad = a T^4 / 3, so u_rad = a T^4 = 3 P_rad)
    # Pure radiation pressure (data.rad_pressure is the 3-component
    # elec+rad sum; we need rad alone for u_rad = 3 P_rad V).
    P_rad = getattr(data, 'rad_pressure_true', None)
    if P_rad is not None and not np.all(P_rad == 0):
        u_rad = (3.0 * P_rad * V_zone).sum(axis=1)
    else:
        u_rad = np.zeros(n_t)
        logger.info("  rad_pressure_true missing or all-zero — U_rad set to 0")

    # Cumulative absorbed laser energy (already in J per data_builder convention)
    led = data.laser_energy_deposited
    if led.ndim > 1:
        # Sum any non-time axes (per-beam axis if present)
        e_absorbed = led.sum(axis=tuple(range(1, led.ndim)))
    else:
        e_absorbed = led.copy()
    e_absorbed = e_absorbed.astype(float)

    # Cumulative fusion energy released: DT neutron count × 17.6 MeV
    # Cumulative fusion energy released: DT neutron count × 17.6 MeV
    ndt = getattr(data, 'dt_neutron_count', None)
    if ndt is not None:
        e_fusion_cum = np.asarray(ndt, dtype=float) * Q_DT_J
    else:
        e_fusion_cum = np.zeros(n_t)

    # Direct boundary tallies — preferred over inferred fractions when present.
    rad_b = getattr(data, 'radiation_energy_at_boundary_cum', None)
    e_rad_boundary = (np.asarray(rad_b, dtype=float)
                      if rad_b is not None else np.zeros(n_t))

    part_esc = getattr(data, 'particle_energy_escaped_cum', None)
    e_particle_esc = (np.asarray(part_esc, dtype=float)
                      if part_esc is not None else np.zeros(n_t))

    # Use direct particle-escape measurement when nonzero; otherwise fall
    # back to nominal DT 3.5/14.1 MeV split.
    if e_particle_esc.max() > 1.0:
        e_neutron = e_particle_esc
        e_alpha   = e_fusion_cum - e_neutron
    else:
        e_alpha   = e_fusion_cum * Q_ALPHA_FRAC
        e_neutron = e_fusion_cum * Q_NEUTRON_FRAC

    # Closure with direct tallies:
    #   E_absorbed + E_alpha_deposited (sources kept in plasma)
    #     = Σ in-plasma channels + E_rad_boundary (escape) + residual_gap
    # residual_gap ≈ 0 means the ledger is closed by the measured tallies;
    # non-zero residual indicates a still-untracked channel or a model error
    # (e.g. ideal-gas U_plasma missing EOS-internal energy).
    sum_channels = ke_inward + ke_outward + u_plasma + u_rad
    residual_gap = (e_absorbed + e_alpha) - sum_channels - e_rad_boundary

    # Closure: absorbed + alpha-deposited = Σ channels + rad escape
    # Neutrons escape entirely — track separately.
    sum_channels = ke_inward + ke_outward + u_plasma + u_rad
    gap_rad      = (e_absorbed + e_alpha) - sum_channels  # actual rad through boundary

    return dict(
        time=t_ns,
        ke_inward=ke_inward,
        ke_outward=ke_outward,
        u_plasma=u_plasma,
        u_rad=u_rad,
        e_absorbed=e_absorbed,
        e_fusion_cum=e_fusion_cum,
        e_alpha=e_alpha,
        e_neutron=e_neutron,
        e_rad_boundary=e_rad_boundary,
        sum_channels=sum_channels,
        residual_gap=residual_gap,
        gap=residual_gap,   # legacy alias so existing print/plot code still works
    )


def snapshot(ledger: dict, t_target_ns: float) -> dict:
    """Energy channels at the timestep nearest t_target_ns."""
    idx = int(np.argmin(np.abs(ledger['time'] - t_target_ns)))
    return {
        'idx':           idx,
        'time':          float(ledger['time'][idx]),
        'e_absorbed':    float(ledger['e_absorbed'][idx]),
        'ke_inward':     float(ledger['ke_inward'][idx]),
        'ke_outward':    float(ledger['ke_outward'][idx]),
        'u_plasma':      float(ledger['u_plasma'][idx]),
        'u_rad':         float(ledger['u_rad'][idx]),
        'e_fusion':      float(ledger['e_fusion_cum'][idx]),
        'sum_channels':  float(ledger['sum_channels'][idx]),
        'e_alpha':        float(ledger['e_alpha'][idx]),
        'e_neutron':      float(ledger['e_neutron'][idx]),
        'e_rad_boundary': float(ledger['e_rad_boundary'][idx]),
        'gap':           float(ledger['gap'][idx]),
    }


# ── Plotting ─────────────────────────────────────────────────────────────

def plot_ledger_timeseries(ax, ledger: dict, run_label: str):
    """Cumulative energy channels vs time."""
    t = ledger['time']
    kJ = 1.0 / 1e3

    ax.plot(t, ledger['e_absorbed']  * kJ, 'k-',          lw=2.5, label='Absorbed laser (cum.)')
    ax.plot(t, ledger['ke_inward']   * kJ, color='tab:blue',   lw=1.5, label='KE inward')
    ax.plot(t, ledger['ke_outward']  * kJ, color='tab:red',    lw=1.5, label='KE outward (blowoff)')
    ax.plot(t, ledger['u_plasma']    * kJ, color='tab:green',  lw=1.5, label='Plasma thermal')
    ax.plot(t, ledger['u_rad']       * kJ, color='tab:orange', lw=1.5, label='Radiation')
    ax.plot(t, ledger['sum_channels'] * kJ, 'k--', lw=1, alpha=0.6, label='Σ channels')

    if ledger['e_fusion_cum'].max() > 1e3:   # > 1 kJ
        ax.plot(t, ledger['e_fusion_cum'] * kJ, color='tab:purple', lw=1.5, label='Fusion released')

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Energy (kJ)')
    ax.set_title(f'{run_label}: energy ledger')
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3)


def plot_snapshot_bars(ax, snap1: dict, snap2: dict,
                       label1: str, label2: str, title: str):
    """Side-by-side bar chart of energy-channel fractions of absorbed."""
    channels = ['KE inward', 'KE outward', 'U plasma', 'U rad', 'Fusion', 'Gap']

    def fractions(s):
        e = s['e_absorbed']
        if e <= 0:
            return [0] * 6
        return [100 * s['ke_inward']  / e,
                100 * s['ke_outward'] / e,
                100 * s['u_plasma']   / e,
                100 * s['u_rad']      / e,
                100 * s['e_fusion']   / e,
                100 * s['gap']        / e]

    f1, f2 = fractions(snap1), fractions(snap2)
    x = np.arange(len(channels))
    w = 0.38

    ax.bar(x - w/2, f1, w, label=label1, color='tab:blue', alpha=0.75)
    ax.bar(x + w/2, f2, w, label=label2, color='tab:red',  alpha=0.75)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(channels, rotation=18, fontsize=9)
    ax.set_ylabel('% of absorbed laser')
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')


# ── Console output ───────────────────────────────────────────────────────

def print_side_by_side(snap1: dict, snap2: dict,
                       label1: str, label2: str, title: str):
    """Side-by-side ASCII comparison table for one snapshot moment."""
    print()
    print('=' * 86)
    print(f'  {title}')
    print('=' * 86)

    def fmt(v):  # kJ
        return f"{v/1e3:>9.2f}"

    e1, e2 = snap1['e_absorbed'], snap2['e_absorbed']

    def row(name, v1, v2):
        if e1 > 0:
            s1 = f"{fmt(v1)} kJ ({100 * v1 / e1:>5.1f}%)"
        else:
            s1 = f"{fmt(v1)} kJ (  ---)"
        if e2 > 0:
            s2 = f"{fmt(v2)} kJ ({100 * v2 / e2:>5.1f}%)"
        else:
            s2 = f"{fmt(v2)} kJ (  ---)"
        print(f"  {name:<27s} {s1:>26s}    {s2:>26s}")

    print(f"  {'':<27s} {label1:>26s}    {label2:>26s}")
    print(f"  {'-' * 81}")
    print(f"  {'Time (ns)':<27s} {snap1['time']:>26.3f}    {snap2['time']:>26.3f}")
    print(f"  {'Absorbed laser (cum.)':<27s} {fmt(e1):>16s} kJ          {fmt(e2):>16s} kJ")
    row('KE inward',          snap1['ke_inward'],   snap2['ke_inward'])
    row('KE outward (blowoff)', snap1['ke_outward'], snap2['ke_outward'])
    row('Plasma thermal',     snap1['u_plasma'],    snap2['u_plasma'])
    row('Radiation',          snap1['u_rad'],       snap2['u_rad'])
    row('Fusion released',    snap1['e_fusion'],    snap2['e_fusion'])
    row('  ↳ Alpha (deposited)', snap1['e_alpha'],    snap2['e_alpha'])
    row('  ↳ Neutron (escaped)', snap1['e_neutron'],  snap2['e_neutron'])
    row('Rad escape (boundary tally)',
                              snap1['e_rad_boundary'], snap2['e_rad_boundary'])
    print(f"  {'-' * 81}")
    row('Σ channels (in-plasma)', snap1['sum_channels'], snap2['sum_channels'])
    row('Residual gap (closure)', snap1['gap'],      snap2['gap'])
    print(f"  {'-' * 81}")
    row('Σ channels',         snap1['sum_channels'], snap2['sum_channels'])
    row('Gap (rad escape, etc)', snap1['gap'],      snap2['gap'])


# ── Loading ──────────────────────────────────────────────────────────────

def load_run(base_path: str):
    """Load an ICFRunData via HeliosRun + ICFRunBuilder."""
    bp = Path(base_path).expanduser()
    exo_path = bp if bp.suffix == '.exo' else bp.with_suffix('.exo')

    logger.info(f"Loading {exo_path.name}")
    run = HeliosRun(str(exo_path), verbose=False)

    rhw_path = exo_path.with_suffix('.rhw')
    rhw_config = None
    if rhw_path.exists():
        try:
            from helios_postprocess.rhw_parser import load_rhw_configuration
            rhw_config = load_rhw_configuration(str(rhw_path))
        except Exception as e:
            logger.warning(f"  Could not parse rhw: {e}")

    data = build_run_data(run, time_unit='s', rhw_config=rhw_config, verbose=False)
    return data


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 3:
        print(__doc__)
        print("\nNeeds exactly two run base paths.\n")
        sys.exit(1)

    base1, base2 = sys.argv[1], sys.argv[2]
    name1 = Path(base1).name
    name2 = Path(base2).name

    data1 = load_run(base1)
    data2 = load_run(base2)

    logger.info(f"Computing energy ledger for {name1}")
    led1 = compute_energy_ledger(data1)
    logger.info(f"Computing energy ledger for {name2}")
    led2 = compute_energy_ledger(data2)

    # Snapshot times — fall back to ledger-derived heuristics if not on data
    def _peak_v_time(data, led):
        v = getattr(data, 'peak_velocity_time', None)
        if v is not None and v > 0:
            return float(v)
        return float(led['time'][int(np.argmax(led['ke_inward']))])

    def _stag_time(data, led):
        v = getattr(data, 'stag_time', None)
        if v is not None and v > 0:
            return float(v)
        return float(led['time'][-1])

    t_pv1, t_pv2     = _peak_v_time(data1, led1), _peak_v_time(data2, led2)
    t_stag1, t_stag2 = _stag_time(data1, led1),   _stag_time(data2, led2)
    t_end1, t_end2   = float(led1['time'][-1]),   float(led2['time'][-1])

    # ── Console ──
    print()
    print('=' * 86)
    print('  Energy-balance diagnostic')
    print(f'  Run 1: {name1}')
    print(f'  Run 2: {name2}')
    print('=' * 86)

    snap1_pv  = snapshot(led1, t_pv1)
    snap2_pv  = snapshot(led2, t_pv2)
    print_side_by_side(snap1_pv,   snap2_pv,
                       name1, name2, 'AT PEAK VELOCITY')

    snap1_stag = snapshot(led1, t_stag1)
    snap2_stag = snapshot(led2, t_stag2)
    print_side_by_side(snap1_stag, snap2_stag,
                       name1, name2, 'AT STAGNATION')

    snap1_end = snapshot(led1, t_end1)
    snap2_end = snapshot(led2, t_end2)
    print_side_by_side(snap1_end,  snap2_end,
                       name1, name2, 'AT END OF SIMULATION')

    # ── PDF ──
    out_path = Path.cwd() / f"{name1}_vs_{name2}_energy_balance.pdf"
    logger.info(f"Writing PDF: {out_path}")

    with PdfPages(out_path) as pdf:
        # Page 1: time-series for both runs
        fig, axes = plt.subplots(2, 1, figsize=(10, 11))
        plot_ledger_timeseries(axes[0], led1, name1)
        plot_ledger_timeseries(axes[1], led2, name2)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2: bar charts at peak v and stagnation
        fig, axes = plt.subplots(2, 1, figsize=(10, 11))
        plot_snapshot_bars(axes[0], snap1_pv, snap2_pv,
                           name1, name2,
                           f'At peak velocity '
                           f'(run1: {snap1_pv["time"]:.2f} ns, '
                           f'run2: {snap2_pv["time"]:.2f} ns)')
        plot_snapshot_bars(axes[1], snap1_stag, snap2_stag,
                           name1, name2,
                           f'At stagnation '
                           f'(run1: {snap1_stag["time"]:.2f} ns, '
                           f'run2: {snap2_stag["time"]:.2f} ns)')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 3: bar chart at end of run
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        plot_snapshot_bars(ax, snap1_end, snap2_end,
                           name1, name2,
                           f'At end of simulation '
                           f'(run1: {snap1_end["time"]:.2f} ns, '
                           f'run2: {snap2_end["time"]:.2f} ns)')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    logger.info(f"Done.  PDF at {out_path}")


if __name__ == "__main__":
    main()
