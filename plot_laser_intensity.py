#!/usr/bin/env python3
"""
plot_laser_intensity.py

Reconstruct laser intensity I(r, t) [W/cm^2] from Helios EXODUS output.

Two independent methods, both plotted for cross-check:

    Method 1 (exact, local):
        I(r) = LaserPwrSrc(r) / laserAttinuationCoeff(r)
    Derived from P_src = alpha * I by definition. Undefined where alpha = 0.

    Method 2 (Beer-Lambert, spherical):
        I(r) = I_outer * (R_out / r)^2 * exp(-tau(r))
    where
        I_outer   = LaserPwrOnTargetForBeam / (4 pi R_out^2)
        tau(r)    = integral from r to R_out of alpha(r') dr'
    The (R_out/r)^2 term is geometric convergence for the 1D spherical beam.

In regions where alpha > 0 and the 1D spherical assumption holds, Method 1
and Method 2 should agree to within zone-discretization error. Disagreement
is diagnostic: either alpha is being evaluated inconsistently at boundaries
vs centers, or there is absorbed power from a mechanism not captured by the
IB attenuation coefficient (e.g. resonant absorption).

Companion to plot_laser_deposition.py. Output: <base>_laser_intensity.pdf.

Usage:
    python3 plot_laser_intensity.py <path to .exo>
    python3 plot_laser_intensity.py <path to .exo> --wavelength_um 0.351 --ntimes 6

Project convention:
    Mac Studio: python3
    MacBook   : python (Anaconda)
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from netCDF4 import Dataset


# ---------------------------------------------------------------------------
# Constants and defaults
# ---------------------------------------------------------------------------

LAMBDA_UM_DEFAULT = 0.351   # 3w Nd:glass (NIF baseline)
N_CR_COEFF = 1.115e21       # cm^-3 * um^2 ; ncr = N_CR_COEFF / lambda_um^2


# ---------------------------------------------------------------------------
# EXODUS ingestion
# ---------------------------------------------------------------------------

def read_exodus(exo_path: str) -> dict:
    """
    Read the variables needed for intensity reconstruction.

    Returns a dict with:
        time_ns          (nt,)
        zbnd             (nt, nzone+1)     zone boundaries [cm]
        zcen             (nt, nzone)       zone centers    [cm]
        laser_pwr_src    (nt, nzone)       [W/cm^3]
        alpha_zone       (nt, nzone)       IB absorption coeff at zone centers [1/cm]
        P_on_target      (nt,)             LaserPwrOnTargetForBeam, beam 0 [W]
        P_delivered      (nt,)             LaserPwrDeliveredForBeam, beam 0 [W]
        ne               (nt, nzone) or None   electron number density [1/cm^3]
    """
    if not os.path.exists(exo_path):
        raise FileNotFoundError(exo_path)

    with Dataset(exo_path, "r") as ds:
        time_s = np.array(ds.variables["time_whole"][:])
        zbnd = np.array(ds.variables["zone_boundaries"][:])  # (nt, nzone+1)

        pwr_src = np.array(ds.variables["LaserPwrSrc"][:])  # (nt, nzone)

        # laserAttinuationCoeff: (nt, nbeam, nbnd). nbnd is usually nzone+1 but
        # some Helios builds have nzone+2 (extra ghost). Take beam 0.
        atten = np.array(ds.variables["laserAttinuationCoeff"][:])  # (nt, nbeam, nbnd)
        if atten.ndim == 3:
            atten_bnd = atten[:, 0, :]
        else:
            atten_bnd = atten  # fallback if beam axis already collapsed

        P_on_target = np.array(ds.variables["LaserPwrOnTargetForBeam"][:])  # (nt, nbeam)
        P_delivered = np.array(ds.variables["LaserPwrDeliveredForBeam"][:]) # (nt, nbeam)
        if P_on_target.ndim == 2:
            P_on_target = P_on_target[:, 0]
            P_delivered = P_delivered[:, 0]

        try:
            ne = np.array(ds.variables["NumElecDensity"][:])
        except KeyError:
            ne = None

    # Derived
    nt, nzone_p1 = zbnd.shape
    nzone = nzone_p1 - 1
    zcen = 0.5 * (zbnd[:, :-1] + zbnd[:, 1:])

    # Trim attenuation to physical boundaries (nzone+1) and handle Helios sentinels.
    # Helios layout observed on VI_6: att has nzone+2 points = 1 leading ghost (at r=0)
    # + 351 physical boundaries. Opaque zones (past critical surface or outside the
    # laser-reachable region looking inward) are flagged with alpha = 1e30, not a real
    # value. These must be clipped before use in exp(-tau) or the whole sum overflows.
    nbnd = atten_bnd.shape[1]
    if nbnd == nzone + 2:
        atten_bnd = atten_bnd[:, 1:]            # drop leading ghost at r=0
    elif nbnd == nzone + 1:
        pass
    else:
        # Not a shape we understand; truncate/pad to nzone+1 and warn
        print(f"[warn] laserAttinuationCoeff has {nbnd} points along last axis; "
              f"expected {nzone+1}. Truncating.", file=sys.stderr)
        atten_bnd = atten_bnd[:, :nzone + 1]

    # Clip the opaque-sentinel values. 1e30 means "laser cannot propagate here"
    # (inside the critical surface, typically). Cap at a large but finite value so
    # exp(-tau) evaluates to zero cleanly rather than overflowing to -inf.
    ALPHA_MAX = 1.0e6   # [1/cm] -- saturates tau on any physical path length
    SENTINEL_THRESHOLD = 1.0e20
    opaque_mask = atten_bnd > SENTINEL_THRESHOLD
    atten_bnd = np.where(opaque_mask, ALPHA_MAX, atten_bnd)

    alpha_zone = 0.5 * (atten_bnd[:, :-1] + atten_bnd[:, 1:])

    return dict(
        time_ns=time_s * 1e9,
        zbnd=zbnd,
        zcen=zcen,
        laser_pwr_src=pwr_src,
        alpha_zone=alpha_zone,
        P_on_target=P_on_target,
        P_delivered=P_delivered,
        ne=ne,
    )


# ---------------------------------------------------------------------------
# Intensity reconstruction
# ---------------------------------------------------------------------------

def intensity_method1(pwr_src: np.ndarray, alpha_zone: np.ndarray) -> np.ndarray:
    """
    I = P_src / alpha, masked to NaN in zones where the ratio is not
    physically meaningful:
      - alpha <= 0 (no absorption mechanism)
      - alpha at the opaque-sentinel ceiling (laser cannot reach)
      - alpha < 1% of per-timestep maximum (tenuous outer corona where
        both alpha and P_src are numerical noise)

    In the outer corona where alpha is orders of magnitude below peak,
    the ratio P_src/alpha becomes noise-dominated and can return
    unphysical values (10^18+ W/cm^2). Filtering on alpha relative to
    its per-timestep maximum keeps the absorbing-layer signal while
    excluding the noisy far corona.
    """
    OPAQUE_FLOOR = 0.9e6
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha_active = alpha_zone.copy()
        alpha_active[alpha_active >= OPAQUE_FLOOR] = 0.0
        alpha_max_t = np.nanmax(alpha_active, axis=1, keepdims=True)
        alpha_thresh = np.where(alpha_max_t > 0, 0.01 * alpha_max_t, 0.0)

        valid = (
            (alpha_zone > 0.0)
            & (alpha_zone < OPAQUE_FLOOR)
            & (alpha_zone > alpha_thresh)
        )
        I = np.where(valid, pwr_src / alpha_zone, np.nan)
    return I


def intensity_method2(zbnd: np.ndarray,
                      alpha_zone: np.ndarray,
                      I_outer: np.ndarray) -> np.ndarray:
    """
    Beer-Lambert intensity in 1D spherical geometry.

    I(r, t) = I_outer(t) * (R_out(t) / r)^2 * exp(-tau(r, t))
    tau(r, t) = integral_r^{R_out} alpha(r') dr'

    Integrate inward from the outer boundary using the zone-centered alpha.
    tau is evaluated at zone centers.
    """
    nt, nz = alpha_zone.shape
    dr = np.diff(zbnd, axis=1)              # (nt, nz)
    zcen = 0.5 * (zbnd[:, :-1] + zbnd[:, 1:])
    R_out = zbnd[:, -1]                     # (nt,)

    alpha_dr = alpha_zone * dr              # (nt, nz)

    # cum_out[:, i] = sum_{j >= i} alpha_j * dr_j  (inclusive)
    cum_out = np.cumsum(alpha_dr[:, ::-1], axis=1)[:, ::-1]

    # tau at outer edge of zone i = sum_{j > i} alpha_j dr_j
    #                             = cum_out[:, i+1] for i < nz-1, else 0
    tau_outer_edge = np.zeros_like(cum_out)
    tau_outer_edge[:, :-1] = cum_out[:, 1:]
    # tau at zone center = tau at outer edge + half the current zone
    tau_center = tau_outer_edge + 0.5 * alpha_dr

    # Geometric convergence; guard r -> 0
    r_safe = np.where(zcen > 1e-12, zcen, 1e-12)
    geom = (R_out[:, None] / r_safe) ** 2

    I = I_outer[:, None] * geom * np.exp(-tau_center)
    return I


def find_critical_radius(ne: np.ndarray,
                         zcen: np.ndarray,
                         lam_um: float) -> tuple[np.ndarray, float]:
    """
    Return (r_crit_vs_t [cm], ncr [1/cm^3]).

    The critical surface is the outermost zone where ne >= ncr -- i.e. the
    first surface the incoming ray encounters where it cannot propagate.
    """
    ncr = N_CR_COEFF / (lam_um ** 2)
    nt = ne.shape[0]
    r_crit = np.full(nt, np.nan)
    for t in range(nt):
        idx = np.where(ne[t] >= ncr)[0]
        if idx.size > 0:
            i = int(idx.max())
            r_crit[t] = zcen[t, i]
    return r_crit, ncr


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def pick_timesteps(P_delivered: np.ndarray, ntimes: int) -> np.ndarray:
    """Pick ntimes indices spread over the laser-active window."""
    if P_delivered.max() <= 0:
        return np.linspace(0, len(P_delivered) - 1, ntimes).astype(int)
    active = np.where(P_delivered > 0.01 * P_delivered.max())[0]
    if active.size == 0:
        return np.linspace(0, len(P_delivered) - 1, ntimes).astype(int)
    if active.size <= ntimes:
        return active
    return np.linspace(active[0], active[-1], ntimes).astype(int)


def plot_report(data: dict,
                I1: np.ndarray,
                I2: np.ndarray,
                I_outer: np.ndarray,
                r_crit: Optional[np.ndarray],
                ncr: Optional[float],
                lam_um: float,
                ntimes: int,
                pdf_path: str) -> None:
    time_ns = data["time_ns"]
    zcen = data["zcen"]
    P_del = data["P_delivered"]

    tidx = pick_timesteps(P_del, ntimes)

    # ---- Page 1: I(r) at selected times, both methods
    # Compute a sensible y-range from the finite, physically meaningful data
    # (exclude the exp(-tau) underflow region inside the critical surface).
    finite_vals = np.concatenate([
        I1[I1 > 1.0].ravel(),
        I2[I2 > 1.0].ravel(),
    ])
    if finite_vals.size > 0:
        ylim_hi = 10 ** np.ceil(np.log10(np.nanmax(finite_vals)) + 0.3)
        ylim_lo = max(1.0, ylim_hi / 1e8)      # 8 decades of dynamic range
    else:
        ylim_lo, ylim_hi = 1e10, 1e17

    rows = 2
    cols = int(np.ceil(len(tidx) / rows))
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 7), sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, ti in zip(axes, tidx):
        r = zcen[ti]
        ax.semilogy(r, I1[ti], lw=1.4, label=r"$I = P_{\rm src}/\alpha$")
        ax.semilogy(r, I2[ti], lw=1.4, ls="--",
                    label=r"$I_{\rm out}(R_{\rm out}/r)^2 e^{-\tau}$")
        if r_crit is not None and np.isfinite(r_crit[ti]):
            ax.axvline(r_crit[ti], color="red", ls=":",
                       label=f"$r_{{\\rm crit}}$ = {r_crit[ti]*1e4:.1f} um")
        ax.set_xlabel("r (cm)")
        ax.set_ylabel("I (W/cm$^2$)")
        ax.set_title(f"t = {time_ns[ti]:.2f} ns,  "
                     f"P = {P_del[ti]/1e12:.1f} TW")
        ax.set_ylim(ylim_lo, ylim_hi)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")
    # hide unused axes
    for ax in axes[len(tidx):]:
        ax.set_visible(False)
    fig.suptitle(f"Laser intensity reconstruction   (lambda = {lam_um} um)")
    fig.tight_layout()

    # ---- Page 2: time history of I at critical surface and peak coronal I
    fig2, ax2 = plt.subplots(figsize=(10, 5.5))
    ax2.semilogy(time_ns, I_outer, lw=1.6, ls="--",
                 label=r"$I_{\rm outer}$ (incident)")

    # Peak I (Method 2) at each time -- "peak coronal intensity" history
    I_peak_t = np.full(len(time_ns), np.nan)
    for t in range(len(time_ns)):
        row = I2[t]
        mask = np.isfinite(row) & (row > 0)
        if mask.any():
            I_peak_t[t] = np.nanmax(row[mask])
    ax2.semilogy(time_ns, I_peak_t, lw=1.8, color="C2",
                 label="peak coronal I (Method 2)")

    if r_crit is not None:
        I_crit_t = np.full(len(time_ns), np.nan)
        for t in range(len(time_ns)):
            if np.isfinite(r_crit[t]):
                # Find the zone index matching r_crit and read I2 there directly.
                # Avoid np.interp across the exp(-tau) cliff at r_crit.
                idx = int(np.argmin(np.abs(zcen[t] - r_crit[t])))
                # Step outward by one zone if we landed in the opaque region
                while idx < I2.shape[1] - 1 and not (np.isfinite(I2[t, idx])
                                                      and I2[t, idx] > 0):
                    idx += 1
                if np.isfinite(I2[t, idx]) and I2[t, idx] > 0:
                    I_crit_t[t] = I2[t, idx]
        ax2.semilogy(time_ns, I_crit_t, lw=1.8, color="C1",
                     label=r"I at $r_{\rm crit}$ (first transparent zone)")

    ax2.set_xlabel("time (ns)")
    ax2.set_ylabel("I (W/cm$^2$)")
    ax2.set_title("Laser intensity history")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend()
    # Clip y-axis to physically meaningful range (I_outer and I_peak_t).
    # Avoid letting any rogue zone near the critical surface drag the scale
    # to ~1e-300 from exp(-tau) underflow.
    y_candidates = []
    for arr in (I_outer, I_peak_t):
        a = arr[np.isfinite(arr) & (arr > 0)]
        if a.size > 0:
            y_candidates.append(a)
    if y_candidates:
        all_y = np.concatenate(y_candidates)
        y_hi = 10 ** np.ceil(np.log10(np.nanmax(all_y)) + 0.3)
        y_lo = max(1.0, y_hi / 1e8)
        ax2.set_ylim(y_lo, y_hi)
    fig2.tight_layout()

    # ---- Page 3: method cross-check, I1/I2 at a representative time
    t_mid = tidx[len(tidx) // 2]
    r = zcen[t_mid]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = I1[t_mid] / I2[t_mid]
    fig3, ax3 = plt.subplots(figsize=(10, 5.5))
    ax3.plot(r, ratio, lw=1.5)
    ax3.axhline(1.0, color="k", ls=":", alpha=0.6)
    ax3.set_xlabel("r (cm)")
    ax3.set_ylabel(r"$I_{\rm M1} / I_{\rm M2}$")
    ax3.set_title(f"Method cross-check at t = {time_ns[t_mid]:.2f} ns "
                  f"(expect ~1 where alpha > 0)")
    ax3.set_ylim(0, 2)
    ax3.grid(True, alpha=0.3)
    fig3.tight_layout()

    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig)
        pdf.savefig(fig2)
        pdf.savefig(fig3)
    plt.close("all")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Reconstruct laser intensity I(r,t) from EXODUS output.")
    p.add_argument("exo_path", help="Path to Helios EXODUS .exo file")
    p.add_argument("--wavelength_um", type=float, default=LAMBDA_UM_DEFAULT,
                   help=f"Laser wavelength in microns (default {LAMBDA_UM_DEFAULT})")
    p.add_argument("--ntimes", type=int, default=6,
                   help="Number of timesteps to plot (default 6)")
    p.add_argument("-o", "--output", default=None,
                   help="Output PDF path (default: <exo base>_laser_intensity.pdf)")
    args = p.parse_args(argv)

    if args.output is None:
        base = os.path.splitext(args.exo_path)[0]
        args.output = base + "_laser_intensity.pdf"

    print(f"[read]  {args.exo_path}")
    data = read_exodus(args.exo_path)

    # I_outer from LaserPwrOnTargetForBeam / (4 pi R_out^2)
    R_out = data["zbnd"][:, -1]
    with np.errstate(divide="ignore", invalid="ignore"):
        I_outer = np.where(R_out > 0,
                           data["P_on_target"] / (4.0 * np.pi * R_out ** 2),
                           0.0)

    print("[calc]  Method 1: I = P_src / alpha")
    I1 = intensity_method1(data["laser_pwr_src"], data["alpha_zone"])

    print("[calc]  Method 2: Beer-Lambert from outer boundary")
    I2 = intensity_method2(data["zbnd"], data["alpha_zone"], I_outer)

    r_crit, ncr = None, None
    if data["ne"] is not None:
        r_crit, ncr = find_critical_radius(data["ne"], data["zcen"],
                                           args.wavelength_um)
        print(f"[crit]  from NumElecDensity, lambda = {args.wavelength_um} um  "
              f"->  ncr = {ncr:.3e} cm^-3")
    else:
        # Fallback strategies for locating the critical surface when
        # NumElecDensity is not in the EXODUS file:
        #   (a) outermost opaque zone + 1 (works when 1e30 sentinels tag a layer
        #       inside the critical surface)
        #   (b) if no sentinels found, use the innermost zone with significant
        #       laser deposition -- P_src peaks at the absorption layer which
        #       sits just outside r_crit.
        alpha = data["alpha_zone"]
        pwr_src = data["laser_pwr_src"]
        zcen_arr = data["zcen"]
        nt, nz = alpha.shape
        r_crit = np.full(nt, np.nan)
        OPAQUE_FLOOR = 9.0e5
        for t in range(nt):
            opaque = np.where(alpha[t] > OPAQUE_FLOOR)[0]
            if opaque.size > 0:
                idx = int(opaque.max()) + 1    # first zone OUTSIDE opaque region
                if idx < nz:
                    r_crit[t] = zcen_arr[t, idx]
                    continue
            # Strategy (b): use innermost zone with significant deposition
            p = pwr_src[t]
            if p.max() > 0:
                active = np.where(p > 0.01 * p.max())[0]
                if active.size > 0:
                    r_crit[t] = zcen_arr[t, int(active.min())]
        ncr = N_CR_COEFF / (args.wavelength_um ** 2)
        n_found = np.isfinite(r_crit).sum()
        print(f"[crit]  fallback from alpha sentinel / P_src peak "
              f"({n_found}/{nt} timesteps resolved)")
        print(f"[crit]  ncr(lambda={args.wavelength_um} um) = {ncr:.3e} cm^-3")

    print(f"[plot]  {args.output}")
    plot_report(data, I1, I2, I_outer, r_crit, ncr,
                args.wavelength_um, args.ntimes, args.output)

    # --- console summary
    print("")
    print("=== summary ===")
    print(f"  peak I_outer         : {np.nanmax(I_outer):.3e} W/cm^2")
    print(f"  peak I (Method 2)    : {np.nanmax(I2):.3e} W/cm^2")
    if r_crit is not None:
        # max over t of interpolated I2 at r_crit
        I_crit_peak = 0.0
        for t in range(len(data["time_ns"])):
            if np.isfinite(r_crit[t]):
                row = I2[t]
                mask = np.isfinite(row) & (row > 0)
                if mask.any():
                    v = float(np.interp(r_crit[t],
                                        data["zcen"][t, mask], row[mask]))
                    if v > I_crit_peak:
                        I_crit_peak = v
        print(f"  peak I at r_crit     : {I_crit_peak:.3e} W/cm^2")
    print(f"  output PDF           : {args.output}")


if __name__ == "__main__":
    main()
