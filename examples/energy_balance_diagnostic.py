"""
energy_balance_diagnostic.py
============================

Energy-flow comparison across two or more Helios runs (N-way).

Tracks where absorbed laser energy ends up across the available channels:

    KE_inward      : kinetic energy of zones moving inward (the implosion)
    KE_outward     : kinetic energy of zones moving outward (ablation blowoff)
    U_plasma       : plasma thermal energy, ideal-gas approx (3/2)(P_i + P_e) V
    U_rad          : radiation thermal energy, 3 P_rad V
    E_fusion       : cumulative fusion energy released (DT reactions × 17.6 MeV)
                     split into alpha (deposited in plasma) and neutron (escapes)
    E_rad_boundary : cumulative radiation through the outer grid boundary
                     (direct EXODUS tally when available)

Closure condition:

    E_absorbed + E_alpha_deposited
        = (KE_in + KE_out + U_plasma + U_rad) + E_rad_boundary + residual_gap

A residual_gap near zero means the ledger closes with the tracked channels.
Non-zero residual indicates either an untracked channel or an EOS-internal-energy
difference (ideal-gas U_plasma vs true material internal energy).

Usage
-----
    python energy_balance_diagnostic.py [--denom-mode MODE] <run1> <run2> [<run3> ...]

    --denom-mode {auto,absorbed,total_input}
        Normalization denominator for channel fractions.  Default is `auto`:
        peak-velocity snapshot uses `absorbed` (familiar "% of absorbed laser",
        cleanly readable as hydro efficiency etc.); stagnation and end-of-run
        snapshots use `total_input` = E_absorbed + E_alpha_deposited so
        percentages don't blow up once fusion ignites.

Accepts 2 or more run base paths.

Example
-------
    python energy_balance_diagnostic.py \\
        ~/Sims/Xcimer/Olson_PDD_20/s016_baseline/Olson_PDD_20_..._s016 \\
        ~/Sims/Xcimer/Olson_PDD_20/s016_FL04/Olson_PDD_20_..._s016_FL04 \\
        ~/Sims/Xcimer/HDD_26/HDD26_DTI40_1ns130_FL04_lrm4_nb/HDD26_..._lrm4_nb

Produces
--------
    <stem>_energy_balance.pdf in cwd
       — page 1: cumulative energy timeseries, one subplot per run
       — page 2: grouped bar charts at peak-velocity and stagnation
       — page 3: grouped bar chart at end of simulation
    Console table with snapshots at peak-velocity, stagnation, end-of-run,
    showing all runs side by side per snapshot.

Notes
-----
Standalone — uses the existing helios_postprocess library to load runs but
does not modify or wire into the main analysis pipeline.

For 2 runs the output filename is preserved as <name1>_vs_<name2>_energy_balance.pdf
for backward compatibility with prior outputs.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Iterable

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

# Per-run colors for bar charts and timeseries-comparison overlays.
RUN_COLORS = ['tab:blue', 'tab:red', 'tab:green',
              'tab:orange', 'tab:purple', 'tab:brown',
              'tab:pink', 'tab:olive']


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

    # Use direct particle-escape measurement only when it captures a
    # sensible fraction of fusion energy. For DT we expect neutron escape
    # ≈ 0.8 × fusion; if the EXODUS tally is much smaller it's likely
    # not populated for this run (varies by Helios version / alpha-
    # transport configuration), so fall back to the nominal fractional
    # split rather than concluding "all fusion stayed as alpha."
    if e_fusion_cum[-1] > 1.0:
        direct_frac = float(e_particle_esc[-1] / e_fusion_cum[-1])
    else:
        direct_frac = 0.0

    if direct_frac > 0.3:
        # Direct tally consistent with DT neutron escape (~0.8 nominal)
        e_neutron = e_particle_esc
        e_alpha   = e_fusion_cum - e_neutron
    else:
        # Direct tally not populated or partial — use nominal DT 14.1/17.6 split
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
        'idx':            idx,
        'time':           float(ledger['time'][idx]),
        'e_absorbed':     float(ledger['e_absorbed'][idx]),
        'ke_inward':      float(ledger['ke_inward'][idx]),
        'ke_outward':     float(ledger['ke_outward'][idx]),
        'u_plasma':       float(ledger['u_plasma'][idx]),
        'u_rad':          float(ledger['u_rad'][idx]),
        'e_fusion':       float(ledger['e_fusion_cum'][idx]),
        'sum_channels':   float(ledger['sum_channels'][idx]),
        'e_alpha':        float(ledger['e_alpha'][idx]),
        'e_neutron':      float(ledger['e_neutron'][idx]),
        'e_rad_boundary': float(ledger['e_rad_boundary'][idx]),
        'gap':            float(ledger['gap'][idx]),
    }


# ── Plotting ─────────────────────────────────────────────────────────────

def plot_ledger_timeseries(ax, ledger: dict, run_label: str):
    """Cumulative energy channels vs time."""
    t = ledger['time']
    kJ = 1.0 / 1e3

    ax.plot(t, ledger['e_absorbed']   * kJ, 'k-',              lw=2.5, label='Absorbed laser (cum.)')
    ax.plot(t, ledger['ke_inward']    * kJ, color='tab:blue',   lw=1.5, label='KE inward')
    ax.plot(t, ledger['ke_outward']   * kJ, color='tab:red',    lw=1.5, label='KE outward (blowoff)')
    ax.plot(t, ledger['u_plasma']     * kJ, color='tab:green',  lw=1.5, label='Plasma thermal')
    ax.plot(t, ledger['u_rad']        * kJ, color='tab:orange', lw=1.5, label='Radiation')
    ax.plot(t, ledger['sum_channels'] * kJ, 'k--', lw=1, alpha=0.6, label='Σ channels')

    if ledger['e_rad_boundary'].max() > 1e3:   # > 1 kJ
        ax.plot(t, ledger['e_rad_boundary'] * kJ, color='tab:brown', lw=1.5,
                label='Rad escape (bdry)')
    if ledger['e_fusion_cum'].max() > 1e3:     # > 1 kJ
        ax.plot(t, ledger['e_fusion_cum'] * kJ, color='tab:purple', lw=1.5,
                label='Fusion released')

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Energy (kJ)')
    ax.set_title(f'{run_label}: energy ledger')
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3)


def plot_snapshot_bars(ax, snapshots: list, labels: list, title: str,
                        denom_mode: str = 'absorbed'):
    """
    Grouped bar chart of energy-channel fractions.

    denom_mode='absorbed'     -> fractions of E_absorbed (laser only).
                                  Includes Fusion as a channel for context.
                                  Best for pre-burn snapshots.
    denom_mode='total_input'  -> fractions of (E_absorbed + E_alpha).
                                  Rad escape is a channel; Fusion is omitted
                                  (alpha is in the denominator, neutrons leave).
                                  Best for post-ignition snapshots.
    """
    n = len(snapshots)

    if denom_mode == 'total_input':
        channels = ['KE inward', 'KE outward', 'U plasma', 'U rad',
                    'Rad escape', 'Gap']

        def denom(s):
            return s['e_absorbed'] + s['e_alpha']

        def fractions(s):
            d = denom(s)
            if d <= 0:
                return [0.0] * len(channels)
            return [
                100 * s['ke_inward']      / d,
                100 * s['ke_outward']     / d,
                100 * s['u_plasma']       / d,
                100 * s['u_rad']          / d,
                100 * s['e_rad_boundary'] / d,
                100 * s['gap']            / d,
            ]

        y_label = '% of total in-plasma input  (E_abs + E_alpha)'
    else:  # 'absorbed'
        channels = ['KE inward', 'KE outward', 'U plasma', 'U rad',
                    'Fusion', 'Rad escape', 'Gap']

        def fractions(s):
            e = s['e_absorbed']
            if e <= 0:
                return [0.0] * len(channels)
            return [
                100 * s['ke_inward']      / e,
                100 * s['ke_outward']     / e,
                100 * s['u_plasma']       / e,
                100 * s['u_rad']          / e,
                100 * s['e_fusion']       / e,
                100 * s['e_rad_boundary'] / e,
                100 * s['gap']            / e,
            ]

        y_label = '% of absorbed laser'

    x = np.arange(len(channels))
    total_w = 0.82
    w = total_w / max(n, 1)

    for i, (s, label) in enumerate(zip(snapshots, labels)):
        offset = (i - (n - 1) / 2.0) * w
        ax.bar(x + offset, fractions(s), w, label=label,
               color=RUN_COLORS[i % len(RUN_COLORS)], alpha=0.78)

    ax.axhline(0, color='k', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(channels, rotation=18, fontsize=9)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(fontsize=8, loc='best')
    ax.grid(True, alpha=0.3, axis='y')


# ── Console output ───────────────────────────────────────────────────────

def print_snapshot_table(snapshots: list, labels: list, title: str,
                          denom_mode: str = 'absorbed',
                          col_w: int = 26, name_w: int = 28):
    """
    N-column comparison table for one snapshot moment.

    denom_mode='absorbed':     percentages relative to E_absorbed (laser only).
                                Existing layout: Fusion shown as a row, alpha and
                                neutron split, rad-escape listed below channels.

    denom_mode='total_input':  percentages relative to (E_absorbed + E_alpha).
                                Layout switches to SOURCES / CHANNELS /
                                INFORMATIONAL sections.  Sum tracked includes
                                rad escape as a channel; percentages sum to
                                ~100% (closure).

    The 'absorbed' layout is best pre-ignition.  The 'total_input' layout is
    best post-ignition where fusion-released >> absorbed-laser would otherwise
    make 'absorbed' percentages blow past 100%.
    """
    n = len(snapshots)
    line_w = name_w + (col_w + 1) * n + 2

    print()
    print('=' * line_w)
    print(f'  {title}')
    if denom_mode == 'total_input':
        print(f'  Denominator: total in-plasma input  (E_abs + E_alpha)')
    else:
        print(f'  Denominator: absorbed laser  (E_abs)')
    print('=' * line_w)

    def fmt_cell(v: float, d: float) -> str:
        if d > 0:
            return f"{v/1e3:>9.2f} kJ ({100 * v / d:>5.1f}%)"
        else:
            return f"{v/1e3:>9.2f} kJ (  ---)"

    def fmt_plain(v: float) -> str:
        return f"{v/1e3:>9.2f} kJ"

    # Header row
    header = f"  {'':<{name_w}s}"
    for label in labels:
        header += f" {label:>{col_w}s}"
    print(header)
    print('  ' + '-' * (line_w - 2))

    # Time row (no %, no kJ)
    time_row = f"  {'Time (ns)':<{name_w}s}"
    for s in snapshots:
        time_row += f" {s['time']:>{col_w}.3f}"
    print(time_row)

    # Per-run denominator for the % columns
    if denom_mode == 'total_input':
        denoms = [s['e_absorbed'] + s['e_alpha'] for s in snapshots]
    else:
        denoms = [s['e_absorbed'] for s in snapshots]

    def row(name: str, key: str):
        line = f"  {name:<{name_w}s}"
        for s, d in zip(snapshots, denoms):
            line += f" {fmt_cell(s[key], d):>{col_w}s}"
        print(line)

    def row_plain(name: str, key: str):
        """Absolute kJ only, no percentage."""
        line = f"  {name:<{name_w}s}"
        for s in snapshots:
            cell = fmt_plain(s[key])
            line += f" {cell:>{col_w}s}"
        print(line)

    def row_value(name: str, vals: list):
        """Absolute kJ from a precomputed list of values, no percentage."""
        line = f"  {name:<{name_w}s}"
        for v in vals:
            cell = fmt_plain(v)
            line += f" {cell:>{col_w}s}"
        print(line)

    if denom_mode == 'total_input':
        # ── SOURCES ──
        print(f"  {'-- SOURCES --':<{name_w}s}")
        row_plain('Absorbed laser (cum.)',  'e_absorbed')
        row_plain('Alpha deposited (cum.)', 'e_alpha')
        print('  ' + '-' * (line_w - 2))
        row_value('TOTAL INPUT (denom)',    denoms)
        print()

        # ── CHANNELS ──
        print(f"  {'-- CHANNELS --':<{name_w}s}")
        row('KE inward',                    'ke_inward')
        row('KE outward (blowoff)',         'ke_outward')
        row('Plasma thermal',               'u_plasma')
        row('Radiation',                    'u_rad')
        row('Rad escape (boundary tally)',  'e_rad_boundary')
        print('  ' + '-' * (line_w - 2))

        # Sum tracked = in-plasma channels + rad escape (closure denominator side)
        sum_tracked = [s['sum_channels'] + s['e_rad_boundary'] for s in snapshots]
        line = f"  {'Sum tracked':<{name_w}s}"
        for v, d in zip(sum_tracked, denoms):
            line += f" {fmt_cell(v, d):>{col_w}s}"
        print(line)
        row('Residual gap (closure)',       'gap')
        print()

        # ── INFORMATIONAL ──
        print(f"  {'-- INFORMATIONAL (not in budget) --':<{name_w}s}")
        row_plain('Neutron escaped (cum.)', 'e_neutron')
        row_plain('Total fusion released',  'e_fusion')
    else:
        # 'absorbed' mode — existing layout
        row_plain('Absorbed laser (cum.)',  'e_absorbed')
        row('KE inward',                    'ke_inward')
        row('KE outward (blowoff)',         'ke_outward')
        row('Plasma thermal',               'u_plasma')
        row('Radiation',                    'u_rad')
        row('Fusion released',              'e_fusion')
        row('  ↳ Alpha (deposited)',        'e_alpha')
        row('  ↳ Neutron (escaped)',        'e_neutron')
        row('Rad escape (boundary tally)',  'e_rad_boundary')
        print('  ' + '-' * (line_w - 2))
        row('Σ channels (in-plasma)',       'sum_channels')
        row('Residual gap (closure)',       'gap')


# Legacy alias — keeps any external callers working.
def print_side_by_side(snap1: dict, snap2: dict,
                       label1: str, label2: str, title: str):
    """Backwards-compatible 2-run wrapper."""
    print_snapshot_table([snap1, snap2], [label1, label2], title)


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


# ── Snapshot-time helpers ────────────────────────────────────────────────

def _peak_v_time(data, led) -> float:
    v = getattr(data, 'peak_velocity_time', None)
    if v is not None and v > 0:
        return float(v)
    return float(led['time'][int(np.argmax(led['ke_inward']))])


def _stag_time(data, led) -> float:
    v = getattr(data, 'stag_time', None)
    if v is not None and v > 0:
        return float(v)
    return float(led['time'][-1])


# ── Output-filename helper ───────────────────────────────────────────────

def _build_output_stem(names: list) -> str:
    """
    Compose a sensible PDF stem from N run names.

    For 2–3 runs, chain with '_vs_' (preserves prior 2-way filename).
    For 4+ runs, use '<first>_and_<N-1>_others' to avoid pathological lengths.
    """
    if len(names) <= 3:
        return '_vs_'.join(names) + '_energy_balance'
    return f'{names[0]}_and_{len(names) - 1}_others_energy_balance'


def _parse_args(argv):
    """Lightweight arg parsing — split out --denom-mode from base paths."""
    denom_mode = 'auto'
    paths = []
    i = 1
    while i < len(argv):
        tok = argv[i]
        if tok == '--denom-mode':
            if i + 1 >= len(argv):
                raise ValueError("--denom-mode requires an argument "
                                  "(one of: auto, absorbed, total_input)")
            denom_mode = argv[i + 1]
            i += 2
        elif tok.startswith('--denom-mode='):
            denom_mode = tok.split('=', 1)[1]
            i += 1
        elif tok in ('-h', '--help'):
            print(__doc__)
            sys.exit(0)
        else:
            paths.append(tok)
            i += 1

    if denom_mode not in ('auto', 'absorbed', 'total_input'):
        raise ValueError(f"--denom-mode must be one of "
                          f"auto/absorbed/total_input (got {denom_mode!r})")
    return paths, denom_mode


def _resolve_mode(global_mode: str, snapshot_kind: str) -> str:
    """
    Map a global denom_mode + snapshot kind to the actual mode for that table.

    'auto' -> 'absorbed' for peak_v, 'total_input' for stag/end.
    Explicit modes pass through unchanged.
    """
    if global_mode in ('absorbed', 'total_input'):
        return global_mode
    # 'auto'
    if snapshot_kind == 'peak_v':
        return 'absorbed'
    return 'total_input'


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    try:
        base_paths, denom_mode = _parse_args(sys.argv)
    except ValueError as e:
        print(__doc__)
        print(f"\n{e}\n")
        sys.exit(1)

    if len(base_paths) < 2:
        print(__doc__)
        print("\nNeeds at least two run base paths.\n")
        sys.exit(1)

    names = [Path(bp).name for bp in base_paths]
    n_runs = len(base_paths)

    # ── Load and compute ledgers ──
    datas = []
    ledgers = []
    for bp, name in zip(base_paths, names):
        d = load_run(bp)
        datas.append(d)
        logger.info(f"Computing energy ledger for {name}")
        ledgers.append(compute_energy_ledger(d))

    # ── Snapshot times (per run; peak-v / stagnation land at different t in each) ──
    t_pv   = [_peak_v_time(d, l) for d, l in zip(datas, ledgers)]
    t_stag = [_stag_time(d, l)   for d, l in zip(datas, ledgers)]
    t_end  = [float(l['time'][-1]) for l in ledgers]

    snaps_pv   = [snapshot(l, t) for l, t in zip(ledgers, t_pv)]
    snaps_stag = [snapshot(l, t) for l, t in zip(ledgers, t_stag)]
    snaps_end  = [snapshot(l, t) for l, t in zip(ledgers, t_end)]

    mode_pv   = _resolve_mode(denom_mode, 'peak_v')
    mode_stag = _resolve_mode(denom_mode, 'stag')
    mode_end  = _resolve_mode(denom_mode, 'end')

    # ── Console ──
    print()
    print('=' * 86)
    print('  Energy-balance diagnostic')
    for i, name in enumerate(names, start=1):
        print(f'  Run {i}: {name}')
    print(f'  Denominator mode: {denom_mode}'
          + ('  (peak-v: absorbed, stag/end: total_input)'
             if denom_mode == 'auto' else ''))
    print('=' * 86)

    print_snapshot_table(snaps_pv,   names, 'AT PEAK VELOCITY',     denom_mode=mode_pv)
    print_snapshot_table(snaps_stag, names, 'AT STAGNATION',        denom_mode=mode_stag)
    print_snapshot_table(snaps_end,  names, 'AT END OF SIMULATION', denom_mode=mode_end)

    # ── PDF ──
    out_path = Path.cwd() / f"{_build_output_stem(names)}.pdf"
    logger.info(f"Writing PDF: {out_path}")

    with PdfPages(out_path) as pdf:
        # Page 1: time-series, one subplot per run, stacked vertically.
        fig, axes = plt.subplots(n_runs, 1, figsize=(10, 5.5 * n_runs))
        if n_runs == 1:
            axes = [axes]
        for ax, led, name in zip(axes, ledgers, names):
            plot_ledger_timeseries(ax, led, name)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 2: grouped bar charts at peak-v and stagnation.
        fig, axes = plt.subplots(2, 1, figsize=(10 + 0.4 * max(n_runs - 2, 0), 11))
        pv_title = ('At peak velocity ('
                    + ', '.join(f'{name}: {t:.2f} ns'
                                for name, t in zip(names, t_pv))
                    + ')')
        stag_title = ('At stagnation ('
                      + ', '.join(f'{name}: {t:.2f} ns'
                                  for name, t in zip(names, t_stag))
                      + ')')
        plot_snapshot_bars(axes[0], snaps_pv,   names, pv_title,
                            denom_mode=mode_pv)
        plot_snapshot_bars(axes[1], snaps_stag, names, stag_title,
                            denom_mode=mode_stag)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Page 3: grouped bar chart at end of simulation.
        fig, ax = plt.subplots(1, 1,
                                figsize=(10 + 0.4 * max(n_runs - 2, 0), 6))
        end_title = ('At end of simulation ('
                     + ', '.join(f'{name}: {t:.2f} ns'
                                 for name, t in zip(names, t_end))
                     + ')')
        plot_snapshot_bars(ax, snaps_end, names, end_title,
                            denom_mode=mode_end)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    logger.info(f"Done.  PDF at {out_path}")


if __name__ == "__main__":
    main()
