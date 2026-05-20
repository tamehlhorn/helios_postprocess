"""
plot_shock_trajectories.py -- Overlay multi-shock trajectories on R-T axes.

Built to compare Helios's shock train against the LILAC/xRAGE/HYDRA R-T
panels in Olson 2021, which show 3 distinct shocks launched through the
cold fuel. Per-shock breakouts at the gas/ice interface are annotated.

Usage:
    python  ~/Codes/helios_postprocess/examples/plot_shock_trajectories.py <base_path>
    python3 ~/helios_postprocess/examples/plot_shock_trajectories.py    <base_path>

<base_path> is the path WITHOUT extension (the script reads <base>.exo).

Output: <base>_shock_trajectories.pdf next to the .exo file.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from helios_postprocess.core import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.icf_analysis import ICFAnalyzer


def main(base):
    print(f"Loading {base}.exo ...")
    run  = HeliosRun(base + ".exo")
    data = build_run_data(run)

    analyzer = ICFAnalyzer(data)
    analyzer.analyze_drive_phase()
    analyzer.analyze_stagnation_phase()   # runs analyze_implosion_phase -> _compute_shock_train

    # data.time nominally in ns but build_run_data default leaves it in
    # seconds for some EXODUS files. Auto-detect by magnitude.
    _raw_t      = np.asarray(data.time)
    t_ns        = _raw_t * 1e9 if float(np.max(_raw_t)) < 1e-3 else _raw_t
    ri          = data.region_interfaces_indices
    trajectories = data.shock_trajectories
    breakouts    = data.shock_breakouts
    coalescences = data.shock_coalescence_events
    events       = getattr(data, 'shock_events', []) or []

    # Map each trajectory id -> shock class via the merged-id list on each
    # consolidated event. Anything not in a consolidated event is uncategorised.
    traj_class = {}
    for ev in events:
        for tr_id in ev.get('merged_trajectory_ids', []):
            traj_class[tr_id] = ev['class']

    print(f"Detected {len(trajectories)} trajectories, "
          f"{len(coalescences)} coalescence events, "
          f"{len(breakouts)} raw breakouts -> {len(events)} consolidated events")
    for k, tr in enumerate(trajectories):
        n = tr['indices'].size
        cls = traj_class.get(k, '-')
        print(f"  traj #{k:>2d}: n={n:>3d}, "
              f"t=[{tr['time_ns'][0]:6.3f}, {tr['time_ns'][-1]:6.3f}] ns, "
              f"r=[{tr['radius'].min()*1e4:7.1f}, {tr['radius'].max()*1e4:7.1f}] um, "
              f"P_ratio_max={tr['P_ratio'].max():5.2f}, "
              f"class={cls:<5s} end={tr['reason_ended']}")
    print("Consolidated breakouts:")
    for e in events:
        print(f"  {e['class']:<5s} (#{e['order']}) at t = {e['time_ns']:6.3f} ns, "
              f"r = {e['radius']*1e4:6.1f} um, "
              f"P_post_max = {e['P_post_Gbar_max']*1000.0:7.2f} Mbar, "
              f"P_ratio_max = {e['P_ratio_max']:5.2f}  "
              f"(merged {e['n_merged']} raw)")

    # ---- Figure: R-T overlay + inward-velocity panel ------------------------
    fig, (ax, axv) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True,
        gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.08},
    )

    # Region boundaries vs time (Lagrangian envelope) on R-T panel
    cap_idx = data.capsule_outer_idx if ri is not None else 0
    n_nodes = data.zone_boundaries.shape[1]
    if ri is not None:
        r_gas_ice = np.array([
            data.zone_boundaries[t, int(ri[t, 0])] * 1e4 for t in range(len(t_ns))
        ])  # um
        ax.plot(t_ns, r_gas_ice, color="k", lw=1.2, ls="--",
                label="gas/ice interface")
        # outer ablator boundary -- clamp because Helios stores ri as
        # slice-end indices (one past the last zone), which can equal n_zones
        # while zone_boundaries has only n_zones+1 elements indexed 0..n_zones.
        r_outer = np.array([
            data.zone_boundaries[t, min(int(ri[t, cap_idx]), n_nodes - 1)] * 1e4
            for t in range(len(t_ns))
        ])
        ax.plot(t_ns, r_outer, color="0.5", lw=1.0, ls=":",
                label="ablator outer boundary")

    # Trajectories: r(t) on top, v_inward = -dr/dt on bottom.
    # Color by shock class (foot/ramp/peak) for the classified ones;
    # uncategorised trajectories drawn semi-transparent gray.
    CLASS_COLORS = {
        'foot': 'tab:blue',
        'ramp': 'tab:orange',
        'peak': 'tab:red',
    }
    UNCLASSIFIED_COLOR = '0.55'
    seen_labels = set()
    traj_velocities = {}
    for k, tr in enumerate(trajectories):
        if tr['indices'].size < 3:
            continue
        cls = traj_class.get(k)
        if cls in CLASS_COLORS:
            c     = CLASS_COLORS[cls]
            lw    = 1.8
            alpha = 0.95
            label = f"{cls} shock" if cls not in seen_labels else None
            seen_labels.add(cls)
        else:
            c     = UNCLASSIFIED_COLOR
            lw    = 1.0
            alpha = 0.55
            label = 'other / pile-up' if 'other' not in seen_labels else None
            seen_labels.add('other')
        ax.plot(tr['time_ns'], tr['radius'] * 1e4,
                color=c, lw=lw, alpha=alpha, label=label)
        # cm/ns -> km/s via x1e4
        v_in = -np.gradient(tr['radius'], tr['time_ns']) * 1e4
        traj_velocities[k] = v_in
        axv.plot(tr['time_ns'], v_in, color=c, lw=lw, alpha=alpha)

    # Coalescence markers (R-T panel only)
    for ev_co in coalescences:
        ax.plot(ev_co['time_ns'], ev_co['radius'] * 1e4,
                marker="x", color="black", ms=8, mew=1.5)

    # One marker per consolidated event (foot/ramp/peak), not per raw breakout.
    for e in events:
        c_e = CLASS_COLORS.get(e['class'], 'firebrick')
        ax.plot(e['time_ns'], e['radius'] * 1e4,
                marker="o", mfc="white", mec=c_e, ms=10, mew=2.0,
                zorder=5)
        ax.annotate(
            f"{e['class']}\nt={e['time_ns']:.2f} ns",
            xy=(e['time_ns'], e['radius'] * 1e4),
            xytext=(8, 8), textcoords="offset points",
            fontsize=9, color=c_e, fontweight='bold',
        )
        # Mirror marker on the velocity panel.
        tr_b = trajectories[e['trajectory_id']]
        v_arr = traj_velocities.get(e['trajectory_id'])
        if v_arr is not None:
            k_b_arr = np.where(tr_b['indices'] == e['t_idx'])[0]
            if k_b_arr.size > 0:
                axv.plot(e['time_ns'], float(v_arr[int(k_b_arr[0])]),
                         marker="o", mfc="white", mec=c_e,
                         ms=10, mew=2.0, zorder=5)

    # Timing marker (stagnation) on both panels
    _stag_raw = float(data.stag_time) if data.stag_time > 0 else 0.0
    t_stag = _stag_raw * 1e9 if 0 < _stag_raw < 1e-3 else (_stag_raw or float(t_ns[-1]))
    if t_stag > 0:
        for _a in (ax, axv):
            _a.axvline(t_stag, color="darkred", lw=1, ls=":", alpha=0.6)
        ax.text(t_stag + 0.05, 0.97, "stagnation", fontsize=8,
                color="darkred", rotation=90, va="top",
                transform=ax.get_xaxis_transform())

    name = os.path.basename(base)
    ax.set_ylabel("Radius (µm)", fontsize=12)
    ax.set_title(f"{name} -- Helios shock train (cf. LILAC 3-shock R-T)",
                 fontsize=12)
    if ri is not None:
        y_top = float(data.zone_boundaries[0, min(int(ri[0, cap_idx]), n_nodes - 1)]) * 1e4 * 1.1
        ax.set_ylim(bottom=0, top=y_top)
    else:
        ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper right", ncol=2)

    # Headline summary box: one line per foot/ramp/peak slot. Missing slots
    # show "—" so the comparison against LILAC's 3-shock pattern is obvious.
    summary_lines = ["Shock arrivals at gas/ice:"]
    for cls in ("foot", "ramp", "peak"):
        ev_cls = next((e for e in events if e['class'] == cls), None)
        if ev_cls is not None:
            summary_lines.append(
                f"  {cls:<5s} {ev_cls['time_ns']:5.2f} ns  "
                f"P={ev_cls['P_post_Gbar_max']*1000:5.1f} Mbar"
            )
        else:
            summary_lines.append(f"  {cls:<5s}   —")
    ax.text(
        0.02, 0.02, "\n".join(summary_lines),
        transform=ax.transAxes, fontsize=9, family="monospace",
        va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="0.6", alpha=0.9),
    )

    axv.set_xlabel("Time (ns)", fontsize=12)
    axv.set_ylabel("Inward shock\nspeed (km/s)", fontsize=11)
    axv.axhline(0.0, color="k", lw=0.6, alpha=0.4)
    axv.grid(True, alpha=0.3)
    axv.set_xlim(left=max(0.0, float(t_ns[0])),
                 right=min(float(t_ns[-1]), t_stag + 0.3))

    fig.tight_layout()

    out = base + "_shock_trajectories.pdf"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    # CSV of consolidated events, ready for batch comparison across runs.
    csv_path = base + "_shock_events.csv"
    with open(csv_path, "w") as f:
        f.write("order,class,time_ns,radius_um,P_post_Gbar_max,"
                "P_ratio_max,n_merged,trajectory_id\n")
        for e in events:
            f.write(f"{e['order']},{e['class']},{e['time_ns']:.4f},"
                    f"{e['radius']*1e4:.2f},{e['P_post_Gbar_max']:.6f},"
                    f"{e['P_ratio_max']:.3f},{e['n_merged']},"
                    f"{e['trajectory_id']}\n")
    print(f"Saved: {csv_path}")

    # One-line summary that batch sweeps can grep for.
    def _fmt(t):
        return f"{t:5.2f}" if t is not None else "    -"
    foot = next((e['time_ns'] for e in events if e['class'] == 'foot'), None)
    ramp = next((e['time_ns'] for e in events if e['class'] == 'ramp'), None)
    peak = next((e['time_ns'] for e in events if e['class'] == 'peak'), None)
    print(
        f"[shock_summary] {os.path.basename(base)}  "
        f"foot={_fmt(foot)}  ramp={_fmt(ramp)}  peak={_fmt(peak)}  (ns)"
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: plot_shock_trajectories.py <base_path>")
        sys.exit(1)
    main(sys.argv[1])
