"""
adiabat_history.py -- Plot mass-averaged adiabat vs time for one or more Helios PDD runs.

Usage (Mac Studio):
    python3 ~/helios_postprocessor/examples/adiabat_history.py \
        ~/Sims/Xcimer/Olson_PDD/PDD_TM_9f6nb/PDD_TM_9f6nb \
        ~/Sims/Xcimer/Olson_PDD/PDD_TM_7f6nb/PDD_TM_7f6nb

Usage (MacBook):
    python ~/Codes/helios_postprocessor/examples/adiabat_history.py <base_path> [<base_path> ...]

Each positional argument is a base path WITHOUT extension.
The script also marks the peak-velocity timestep for each run as a dashed vertical line,
which is the key diagnostic for whether the foot/adiabat trend is a timestep-selection
artifact or a genuine physical difference in shock heating.

Output: adiabat_history.pdf in the current working directory.

Physics notes:
  - Adiabat alpha = P / P_Fermi
  - P_Fermi = 2.17 * (rho / rho_0)^(5/3)  Mbar  [Lindl convention]
  - rho_0 = 0.205 g/cc  (equimolar DT ice)
  - Evaluated in cold unablated fuel only: zone indices ri[t,0] to ri[t,1]
  - Pressures on ICFRunData are in J/cm^3; multiply by 1e-8 for Gbar, 1e-7 for Mbar
  - Mass-averaged over cold fuel zones at each timestep
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


# ---------------------------------------------------------------------------
# Physics constants
# ---------------------------------------------------------------------------
RHO0_DT = 0.205          # g/cc, equimolar DT ice (Lindl convention)
P_FERMI_COEFF_MBAR = 2.17  # coefficient in Mbar when rho in g/cc


def compute_adiabat_history(base_path):
    """
    Load a Helios run and compute mass-averaged adiabat of cold unablated fuel
    at every EXODUS timestep.

    Returns
    -------
    time_ns       : 1D array, simulation time in ns
    adiabat_t     : 1D array, mass-averaged adiabat at each timestep
    pv_idx        : int, index of peak implosion velocity timestep
    stag_idx      : int, index of stagnation timestep
    run_label     : str, short name derived from base_path
    """
    exo_path = base_path + ".exo"
    if not os.path.exists(exo_path):
        raise FileNotFoundError(f"EXODUS file not found: {exo_path}")

    run = HeliosRun(exo_path)
    data = build_run_data(run)

    analyzer = ICFAnalyzer(data)
    analyzer.analyze_drive_phase()
    analyzer.analyze_stagnation_phase()   # populates peak_velocity_index, stag_time_index
    analyzer.analyze_burn_phase()

    ri   = data.region_interfaces_indices  # (n_times, n_regions+1), NODE indices
    rho  = data.density                    # (n_times, n_zones), g/cc
    # Pressures in J/cm^3 on ICFRunData; convert to Mbar (x1e-7) for P_Fermi comparison
    P_total = (data.ion_pressure + data.rad_pressure) * 1e-7   # Mbar

    n_times = rho.shape[0]
    adiabat_t = np.full(n_times, np.nan)

    for t in range(n_times):
        hs_node = int(ri[t, 0])   # hot-spot boundary (node index = zone index here)
        cf_node = int(ri[t, 1])   # cold fuel / ablated fuel boundary
        if cf_node <= hs_node:
            continue
        z0, z1 = hs_node, cf_node
        rho_z = rho[t, z0:z1]
        P_z   = P_total[t, z0:z1]
        m_z   = data.zone_mass[t, z0:z1]   # grams
        total_mass = m_z.sum()
        if total_mass == 0 or len(rho_z) == 0:
            continue
        P_Fermi = P_FERMI_COEFF_MBAR * (rho_z / RHO0_DT) ** (5.0 / 3.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            alpha_z = np.where(P_Fermi > 0, P_z / P_Fermi, 0.0)
        adiabat_t[t] = np.sum(alpha_z * m_z) / total_mass

    run_label = os.path.basename(base_path)
    return data.time_ns, adiabat_t, data.peak_velocity_index, data.stag_time_index, run_label


def plot_adiabat_histories(base_paths, output_pdf="adiabat_history.pdf",
                           t_min=8.0, t_max=15.5, alpha_max=8.0):
    """
    Compute and overplot adiabat histories for a list of base paths.
    Dashed vertical lines mark the peak-velocity timestep for each run.
    """
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, bp in enumerate(base_paths):
        color = colors[i % len(colors)]
        try:
            t, alpha, pv_idx, stag_idx, label = compute_adiabat_history(bp)
        except Exception as e:
            print(f"WARNING: could not process {bp}: {e}")
            continue

        # Report scalar values to stdout for cross-checking
        pv_adiabat = alpha[pv_idx] if not np.isnan(alpha[pv_idx]) else float("nan")
        stag_adiabat = alpha[stag_idx] if stag_idx < len(alpha) and not np.isnan(alpha[stag_idx]) else float("nan")
        print(f"{label}")
        print(f"  peak-v idx={pv_idx}  t={t[pv_idx]:.3f} ns  adiabat={pv_adiabat:.3f}")
        print(f"  stag  idx={stag_idx}  t={t[stag_idx]:.3f} ns  adiabat={stag_adiabat:.3f}")

        ax.plot(t, alpha, color=color, label=label)
        ax.axvline(t[pv_idx], color=color, linestyle="--", linewidth=1.0,
                   label=f"{label} peak-v @ {t[pv_idx]:.3f} ns  (α={pv_adiabat:.2f})")

    ax.axhline(3.0, color="black", linestyle=":", linewidth=1.2, label="Target α = 3.0")
    ax.set_xlabel("Time (ns)", fontsize=12)
    ax.set_ylabel("Mass-avg adiabat α (cold DT fuel)", fontsize=12)
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(0, alpha_max)
    ax.set_title("Adiabat history — cold unablated DT fuel [ri(t,0) → ri(t,1)]", fontsize=11)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_pdf)
    print(f"\nSaved: {output_pdf}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    base_paths = [os.path.expanduser(p) for p in sys.argv[1:]]
    plot_adiabat_histories(base_paths)
