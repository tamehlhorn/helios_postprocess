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

    t_ns        = data.time
    ri          = data.region_interfaces_indices
    trajectories = data.shock_trajectories
    breakouts    = data.shock_breakouts
    coalescences = data.shock_coalescence_events

    print(f"Detected {len(trajectories)} trajectories, "
          f"{len(coalescences)} coalescence events, "
          f"{len(breakouts)} gas/ice breakouts")
    for b in breakouts:
        print(f"  trajectory #{b['trajectory_id']:>2d}  "
              f"breakout at t = {b['time_ns']:6.3f} ns, "
              f"r = {b['radius']*1e4:6.1f} um, "
              f"P_post = {b['P_post_Gbar']:6.1f} Gbar, "
              f"P_ratio = {b['P_ratio']:5.2f}")

    # ---- Figure: R-T overlay --------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))

    # Region boundaries vs time (Lagrangian envelope)
    if ri is not None:
        r_gas_ice = np.array([
            data.zone_boundaries[t, int(ri[t, 0])] * 1e4 for t in range(len(t_ns))
        ])  # um
        ax.plot(t_ns, r_gas_ice, color="k", lw=1.2, ls="--",
                label="gas/ice interface")
        # outer ablator boundary
        cap_idx = data.capsule_outer_idx
        r_outer = np.array([
            data.zone_boundaries[t, int(ri[t, cap_idx])] * 1e4 for t in range(len(t_ns))
        ])
        ax.plot(t_ns, r_outer, color="0.5", lw=1.0, ls=":",
                label="ablator outer boundary")

    # Trajectories
    colors = plt.cm.tab10(np.linspace(0, 1, max(10, len(trajectories))))
    for k, tr in enumerate(trajectories):
        # Only plot trajectories of meaningful length
        if tr['indices'].size < 3:
            continue
        ax.plot(tr['time_ns'], tr['radius'] * 1e4,
                color=colors[k % len(colors)], lw=1.4, alpha=0.9,
                label=f"shock #{k}")

    # Coalescence markers
    for ev in coalescences:
        ax.plot(ev['time_ns'], ev['radius'] * 1e4,
                marker="x", color="black", ms=8, mew=1.5)

    # Breakout annotations
    for b in breakouts:
        ax.plot(b['time_ns'], b['radius'] * 1e4,
                marker="o", mfc="white", mec="firebrick", ms=8, mew=1.5)
        ax.annotate(
            f"BO #{b['trajectory_id']}\nt={b['time_ns']:.2f} ns",
            xy=(b['time_ns'], b['radius'] * 1e4),
            xytext=(8, 8), textcoords="offset points",
            fontsize=8, color="firebrick",
        )

    # Timing markers
    t_stag = float(data.stag_time) if data.stag_time > 0 else float(t_ns[-1])
    if t_stag > 0:
        ax.axvline(t_stag, color="darkred", lw=1, ls=":", alpha=0.6)
        ax.text(t_stag + 0.05, 0.97, "stagnation", fontsize=8,
                color="darkred", rotation=90, va="top",
                transform=ax.get_xaxis_transform())

    name = os.path.basename(base)
    ax.set_xlabel("Time (ns)", fontsize=12)
    ax.set_ylabel("Radius (µm)", fontsize=12)
    ax.set_title(f"{name} -- Helios shock train (cf. LILAC 3-shock R-T)",
                 fontsize=12)
    ax.set_xlim(left=max(0.0, float(t_ns[0])),
                right=min(float(t_ns[-1]), t_stag + 0.3))
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    fig.tight_layout()

    out = base + "_shock_trajectories.pdf"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: plot_shock_trajectories.py <base_path>")
        sys.exit(1)
    main(sys.argv[1])
