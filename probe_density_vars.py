"""
Probe an EXODUS file for any variable that could give us electron number density,
either directly or via ne = Zbar * ni reconstruction.

Usage:
    python3 probe_density_vars.py <path/to/file.exo>

Prints, for every candidate variable:
  - name
  - shape
  - dtype
  - min / max / median over non-NaN finite values (to sanity-check units)
  - the netCDF attributes Helios writes (if any), which usually include a
    units string and a human description
"""

import sys
import numpy as np
from netCDF4 import Dataset

# Substrings that are plausibly related to electron / ion / charge density.
# Cast a wide net; we'll eyeball the shapes and magnitudes afterward.
KEYWORDS = [
    "elec", "electron",
    "ne_", "_ne", "nelec", "n_e",
    "ion_num", "ion_den", "ionnum", "ionden",
    "numion", "numelec",
    "zbar", "z_bar", "avgion", "avg_ion", "avgcharge", "charge",
    "density", "dens", "nbar",
    "mass_frac", "massfrac", "mfrac",
    "rho",  # mass density -- used to back out ni if atomic mass available
    "atmass", "atomic", "aweight", "awght", "matid", "region",
]


def summarize(var_name: str, var) -> str:
    try:
        arr = var[:]
    except Exception as exc:
        return f"  (could not read: {exc})"

    arr = np.asarray(arr)
    finite = np.isfinite(arr)
    if not finite.any():
        return f"  shape={arr.shape} dtype={arr.dtype}  (no finite values)"

    vals = arr[finite]
    lo, hi, med = float(vals.min()), float(vals.max()), float(np.median(vals))
    lines = [f"  shape={arr.shape} dtype={arr.dtype}"]
    lines.append(f"  min={lo:.3e}  max={hi:.3e}  median={med:.3e}")

    # netCDF attributes (units, description, long_name, etc.)
    attrs = {}
    for a in var.ncattrs():
        try:
            attrs[a] = var.getncattr(a)
        except Exception:
            pass
    if attrs:
        attr_str = ", ".join(f"{k}={v!r}" for k, v in attrs.items())
        lines.append(f"  attrs: {attr_str}")
    return "\n".join(lines)


def main(path: str) -> None:
    print(f"[probe] opening {path}")
    nc = Dataset(path, "r")
    all_names = sorted(nc.variables.keys())
    print(f"[probe] {len(all_names)} total variables in file\n")

    # Find matches (case-insensitive substring match).
    matches = []
    lowered = {n.lower(): n for n in all_names}
    for kw in KEYWORDS:
        for low, orig in lowered.items():
            if kw in low and orig not in matches:
                matches.append(orig)

    # Also surface anything whose units *attribute* mentions cm^-3 or /cc,
    # since that's the dead giveaway for a number density.
    for name in all_names:
        var = nc.variables[name]
        for a in var.ncattrs():
            try:
                val = str(var.getncattr(a)).lower()
            except Exception:
                continue
            if any(tok in val for tok in ("cm^-3", "cm-3", "cm**-3", "1/cc", "/cc", "particles/cm")):
                if name not in matches:
                    matches.append(name)
                    break

    matches.sort()
    print(f"[probe] {len(matches)} candidate variables matched keywords:\n")

    for name in matches:
        var = nc.variables[name]
        print(f"{name}")
        print(summarize(name, var))
        print()

    nc.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python3 probe_density_vars.py <path/to/file.exo>")
        sys.exit(1)
    main(sys.argv[1])
