"""
Cross-check two independent methods for locating the critical surface r_crit(t):

  Method A (direct):  innermost r where elec_density(r, t) >= n_crit
  Method B (fallback): innermost r where laserAttinuationCoeff(r, t) < 1e20
                       i.e. first zone outside the Helios opaque-sentinel region,
                       stepped inward one zone (matches plot_laser_intensity.py logic)

If the fallback has been valid, |r_A - r_B| should be <= ~1 zone width at all
laser-on timesteps. Any systematic offset is diagnostic.

Usage:
    python3 rcrit_crosscheck.py <path/to/file.exo> [--wavelength_um 0.351]
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset


def compute_rcrit_from_ne(ne, zbnd_zones, n_crit):
    """
    ne: (nt, nzone) electron density [cm^-3] at zone centers
    zbnd_zones: (nt, nzone) zone center radii [cm] (NOT boundaries)
    Returns: (nt,) r_crit [cm], or np.nan where the plasma is everywhere subcritical
    """
    nt, nz = ne.shape
    r_crit = np.full(nt, np.nan)
    for t in range(nt):
        # Work from outer boundary inward; first zone at/above n_crit wins.
        # Equivalently: the outermost overcritical zone when scanning inward.
        overcrit = ne[t] >= n_crit
        if not overcrit.any():
            continue
        # Find the OUTERMOST index where ne is still >= n_crit,
        # i.e. the last True when scanning from inside outward.
        # That zone's outer boundary is r_crit.
        idx = np.where(overcrit)[0].max()
        r_crit[t] = zbnd_zones[t, idx]
    return r_crit


def compute_rcrit_from_att(att, zbnd_zones, sentinel=1e20):
    """
    att: (nt, nzone+1) attenuation at zone boundaries, leading ghost already dropped
         so shape matches zbnd (nt, nzone+1). If it's (nt, nbeam, nzone+2), caller
         should pre-slice.
    zbnd_zones: (nt, nzone) zone center radii
    Returns: (nt,) r_crit via opaque-sentinel method
    """
    nt = att.shape[0]
    r_crit = np.full(nt, np.nan)
    # att at zone centers (average adjacent boundaries)
    att_c = 0.5 * (att[:, :-1] + att[:, 1:])
    for t in range(nt):
        opaque = att_c[t] >= sentinel
        if not opaque.any():
            # laser penetrates fully -- no critical surface reached
            continue
        # Outermost opaque zone, step one zone outward for the surface
        i_op = np.where(opaque)[0].max()
        i_crit = min(i_op + 1, att_c.shape[1] - 1)
        r_crit[t] = zbnd_zones[t, i_crit]
    return r_crit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exo_path")
    ap.add_argument("--wavelength_um", type=float, default=0.351)
    ap.add_argument("--outfile", default=None,
                    help="PDF output path (default: <exo>_rcrit_crosscheck.pdf)")
    args = ap.parse_args()

    n_crit = 1.115e21 / args.wavelength_um**2
    print(f"[crosscheck] opening {args.exo_path}")
    print(f"[crosscheck] n_crit(lambda={args.wavelength_um} um) = {n_crit:.3e} cm^-3")

    nc = Dataset(args.exo_path, "r")

    # Pull arrays
    t_ns = nc.variables["time_whole"][:] * 1e9
    zbnd = nc.variables["zone_boundaries"][:]  # (nt, nzone+1)
    zc = 0.5 * (zbnd[:, :-1] + zbnd[:, 1:])    # (nt, nzone)
    ne = nc.variables["elec_density"][:]       # (nt, nzone)
    att = nc.variables["laserAttinuationCoeff"][:]  # (nt, nbeam, nzone+2)
    pwr = nc.variables["LaserPwrOnTargetForBeam"][:]  # (nt, nbeam)
    nc.close()

    # Trim att: drop leading ghost, take beam 0
    if att.ndim == 3:
        att = att[:, 0, :]
    if att.shape[-1] == zbnd.shape[-1] + 1:
        att = att[:, 1:]
    assert att.shape == zbnd.shape, f"att {att.shape} vs zbnd {zbnd.shape}"

    # Compute both r_crit trajectories
    r_ne = compute_rcrit_from_ne(ne, zc, n_crit)
    r_att = compute_rcrit_from_att(att, zbnd, sentinel=1e20)

    # Laser-on mask: P > 1% of peak
    p_peak = pwr[:, 0].max()
    laser_on = pwr[:, 0] > 0.01 * p_peak

    both_valid = np.isfinite(r_ne) & np.isfinite(r_att) & laser_on
    if not both_valid.any():
        print("[crosscheck] no overlapping laser-on timesteps with both methods valid")
        sys.exit(1)

    # Difference in zone widths (use zone width at r_crit location)
    diff_cm = r_ne[both_valid] - r_att[both_valid]
    # Use zone width at current r_crit position (median zone near r_att)
    dz_at_rcrit = np.empty(int(both_valid.sum()))
    dz_all = np.diff(zbnd, axis=1)
    for k, t in enumerate(np.where(both_valid)[0]):
        # Find zone nearest to r_att[t]
        j = int(np.argmin(np.abs(zc[t] - r_att[t])))
        dz_at_rcrit[k] = dz_all[t, j]
    diff_zones = diff_cm / dz_at_rcrit

    # Summary
    print()
    print(f"  timesteps with laser on:     {int(laser_on.sum())} / {len(t_ns)}")
    print(f"  timesteps both methods valid: {int(both_valid.sum())}")
    print()
    print(f"  r_crit_ne  - r_crit_att (cm):  mean={diff_cm.mean():+.3e}  "
          f"median={np.median(diff_cm):+.3e}  max|d|={np.max(np.abs(diff_cm)):.3e}")
    print(f"  in zone widths:                mean={diff_zones.mean():+.2f}  "
          f"median={np.median(diff_zones):+.2f}  max|d|={np.max(np.abs(diff_zones)):.2f}")
    print()
    print(f"  verdict: ", end="")
    if np.max(np.abs(diff_zones)) <= 1.5:
        print("AGREE (max |dr| <= 1.5 zones) -- sentinel fallback validated")
    elif np.median(np.abs(diff_zones)) <= 1.5:
        print("MOSTLY AGREE (median <= 1.5 zones, some outliers)")
    else:
        print("DISAGREE -- sentinel fallback biased; investigate")

    # Plot
    outfile = args.outfile or args.exo_path.replace(".exo", "_rcrit_crosscheck.pdf")
    with PdfPages(outfile) as pdf:
        fig, axes = plt.subplots(2, 1, figsize=(8.5, 9), sharex=True)

        ax = axes[0]
        ax.plot(t_ns, r_ne * 1e4, label="Method A: ne >= n_crit", lw=1.5)
        ax.plot(t_ns, r_att * 1e4, label="Method B: att sentinel", lw=1.5, ls="--")
        ax.fill_between(t_ns, 0, 1, where=laser_on, alpha=0.1, color="red",
                        transform=ax.get_xaxis_transform(), label="laser on")
        ax.set_ylabel("r_crit (um)")
        ax.set_title(f"r_crit(t) cross-check: {args.exo_path.split('/')[-1]}")
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3)

        ax = axes[1]
        t_valid = t_ns[both_valid]
        ax.plot(t_valid, diff_zones, "o-", ms=3, lw=0.8)
        ax.axhline(0, color="k", lw=0.5)
        ax.axhline(1, color="gray", lw=0.5, ls=":")
        ax.axhline(-1, color="gray", lw=0.5, ls=":")
        ax.set_xlabel("time (ns)")
        ax.set_ylabel("(r_ne - r_att) / dz_local")
        ax.set_title("difference in zone widths at local r_crit")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        pdf.savefig(fig, dpi=100)
        plt.close(fig)

    print(f"\n  plot: {outfile}")


if __name__ == "__main__":
    main()
