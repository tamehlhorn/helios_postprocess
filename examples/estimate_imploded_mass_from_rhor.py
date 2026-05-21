"""
estimate_imploded_mass_from_rhor.py -- Estimate imploded DT mass by integrating
a digitized ρ(r) profile.

Built to derive a LILAC / xRAGE / HYDRA imploded-DT-mass reference value
from a digitized Olson 2021 Fig 7/8 (rho vs radius at stagnation or peak
compression). No published imploded-DT-mass exists for these codes;
integrating the published ρ(r) is the cleanest way to estimate one.

The integral:
                  r_outer
    M_DT  =  4π ∫        ρ(r) r² dr        (spherical mass shell)
                 r_inner

For a stagnation / peak-compression profile, integrate from the gas/ice
interface (or wherever ρ first jumps) out to the cold-fuel/ablator
boundary (or wherever ρ drops back below some threshold).

Input format (CSV with header):
    r_um, rho_gcc
    100.0, 0.10
    105.0, 0.50
    ...

Or with two columns named r_cm, rho_gcc; the script auto-detects units.

Usage:
    python3 ~/helios_postprocess/examples/estimate_imploded_mass_from_rhor.py \\
        olson2021_fig7_LILAC.csv \\
        olson2021_fig7_xRAGE.csv \\
        olson2021_fig7_HYDRA.csv

    # Integrate only between specified radii (cm or µm — auto-detected from
    # input column units):
    python3 ~/.../estimate_imploded_mass_from_rhor.py file.csv \\
        --r-inner-um 50 --r-outer-um 250

    # Specify a density threshold for auto-bounds:
    python3 ~/.../estimate_imploded_mass_from_rhor.py file.csv \\
        --rho-threshold 1.0   # g/cc; integrate where ρ > 1 g/cc
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


def load_rhor(path: Path) -> Tuple[np.ndarray, np.ndarray, str]:
    """Read CSV with two columns; return (r_cm, rho_gcc, label).

    Auto-detects radius units from header: 'r_um', 'radius_um', 'r [um]'
    -> µm; 'r_cm', 'r [cm]' -> cm. Defaults to cm if ambiguous and the
    max value is < 1 (capsule-scale radii are typically 1e-2 - 1e-1 cm).
    """
    rows = []
    header = None
    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            row = [c.strip() for c in row if c.strip()]
            if not row:
                continue
            try:
                vals = [float(x) for x in row]
                rows.append(vals)
            except ValueError:
                if header is None:
                    header = [c.lower() for c in row]
                continue
    if not rows:
        raise RuntimeError(f"{path}: no numeric data found")
    arr = np.asarray(rows)
    if arr.shape[1] < 2:
        raise RuntimeError(f"{path}: need >=2 columns (r, rho)")
    r_raw   = arr[:, 0]
    rho_gcc = arr[:, 1]

    # Detect units
    unit = "cm"
    if header is not None:
        joined = " ".join(header)
        if "um" in joined or "µm" in joined or "micron" in joined:
            unit = "um"
    elif float(np.max(r_raw)) > 1.0:
        # >1 in raw column almost certainly means µm
        unit = "um"
    r_cm = r_raw * (1e-4 if unit == "um" else 1.0)

    return r_cm, rho_gcc, path.stem


def integrate_mass(
    r_cm: np.ndarray,
    rho_gcc: np.ndarray,
    r_inner_cm: Optional[float] = None,
    r_outer_cm: Optional[float] = None,
    rho_threshold: Optional[float] = None,
) -> Tuple[float, float, float]:
    """Return (M_grams, r_inner_used_cm, r_outer_used_cm).

    Trapezoidal integration of  M = 4π ∫ ρ(r) r² dr  on the provided
    sample points. ρ in g/cm³, r in cm; result is in grams.

    Window selection priority:
        1. Explicit r_inner_cm / r_outer_cm if given.
        2. If rho_threshold is set, use the contiguous range where
           ρ >= threshold (largest such window).
        3. Otherwise use the full profile.
    """
    order = np.argsort(r_cm)
    r_cm  = r_cm[order]
    rho   = rho_gcc[order]

    if r_inner_cm is None and r_outer_cm is None and rho_threshold is not None:
        mask = rho >= rho_threshold
        if not mask.any():
            raise RuntimeError(
                f"no points with ρ >= {rho_threshold} g/cc; "
                f"max ρ in profile = {rho.max():.3f}"
            )
        idxs = np.where(mask)[0]
        i_lo, i_hi = int(idxs[0]), int(idxs[-1])
        r_inner_cm = float(r_cm[i_lo])
        r_outer_cm = float(r_cm[i_hi])

    if r_inner_cm is None: r_inner_cm = float(r_cm[0])
    if r_outer_cm is None: r_outer_cm = float(r_cm[-1])

    sel = (r_cm >= r_inner_cm) & (r_cm <= r_outer_cm)
    r_sel   = r_cm[sel]
    rho_sel = rho[sel]
    if r_sel.size < 2:
        raise RuntimeError("not enough sample points inside the integration window")

    integrand = rho_sel * r_sel ** 2
    M_g = 4.0 * np.pi * float(np.trapz(integrand, r_sel))
    return M_g, r_inner_cm, r_outer_cm


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('files', nargs='+', type=Path,
                    help='CSV files with (r, rho) columns. r in cm or µm '
                         '(auto-detected from header).')
    ap.add_argument('--r-inner-um', type=float, default=None,
                    help='Inner integration radius (µm). Default: auto.')
    ap.add_argument('--r-outer-um', type=float, default=None,
                    help='Outer integration radius (µm). Default: auto.')
    ap.add_argument('--rho-threshold', type=float, default=None,
                    help='Auto-bound to where ρ ≥ threshold (g/cc). '
                         'Default: full profile.')
    args = ap.parse_args(argv)

    r_inner_cm = args.r_inner_um * 1e-4 if args.r_inner_um is not None else None
    r_outer_cm = args.r_outer_um * 1e-4 if args.r_outer_um is not None else None

    rows: List[Tuple[str, float, float, float, float]] = []
    print()
    print(f"{'profile':<40s}  {'r_in (µm)':>10s}  {'r_out (µm)':>11s}  "
          f"{'M (mg)':>9s}  {'<ρ> (g/cc)':>11s}")
    print('-' * 90)
    for path in args.files:
        if not path.is_file():
            print(f"ERROR: {path} not found", file=sys.stderr)
            return 1
        try:
            r_cm, rho_gcc, label = load_rhor(path)
            M_g, rin, rout = integrate_mass(
                r_cm, rho_gcc,
                r_inner_cm=r_inner_cm,
                r_outer_cm=r_outer_cm,
                rho_threshold=args.rho_threshold,
            )
        except Exception as e:
            print(f"ERROR processing {path}: {e}", file=sys.stderr)
            return 1
        M_mg     = M_g * 1e3
        V_cm3    = (4.0 / 3.0) * np.pi * (rout ** 3 - rin ** 3)
        rho_avg  = M_g / V_cm3 if V_cm3 > 0 else 0.0
        rows.append((label, rin * 1e4, rout * 1e4, M_mg, rho_avg))
        print(f"{label:<40s}  {rin*1e4:>10.1f}  {rout*1e4:>11.1f}  "
              f"{M_mg:>9.3f}  {rho_avg:>11.3f}")

    # Comparison to Helios PDD_20 if we have something to compare to
    HELIOS_PDD20_MG = 0.60
    if rows:
        print()
        print(f"Helios PDD_20 imploded DT (validated):  {HELIOS_PDD20_MG:.2f} mg")
        for label, _, _, M_mg, _ in rows:
            if M_mg > 0:
                delta = 100 * (HELIOS_PDD20_MG - M_mg) / M_mg
                tag = "MATCH" if abs(delta) < 20 else "STRUCTURAL"
                print(f"  vs {label:<36s} = {HELIOS_PDD20_MG/M_mg:5.2f} × ref  "
                      f"(Helios is {delta:+6.1f}% off)  [{tag}]")
    return 0


if __name__ == '__main__':
    sys.exit(main())
