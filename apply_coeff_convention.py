#!/usr/bin/env python
"""Make the printed 'first-principles coeff' self-consistent.

Problem
-------
In the multi-material (RHW) path, neutron_report divided the *center-to-edge*
hydro column (metrics['rhoR_hydro_gcm2'], Method 1) by a DSR that is now derived
from the *traversed* (Method-2) areal density.  That mixes conventions -- e.g.
1.42 / 0.0452 = 31.4 -- so the number is neither the center-to-edge nor the
traversed coefficient.

Fix
---
Use the SAME areal-density convention on both sides.  The DSR now comes from the
traversed rhoR, so the coefficient divides the traversed rhoR by that DSR.  The
primary reported number is the total (D+T+C) traversed rhoR over the total DSR --
the quantity the empirical NIF 20.4 g/cm^2-per-DSR actually calibrates -- and a
fuel-only (D+T) self-consistent coefficient is reported alongside it.

The single-scatter (no-RHW) path is left untouched: there both the DSR and
rhoR_hydro are center-to-edge, so it is already self-consistent.

Usage:
    python apply_coeff_convention.py /path/to/helios_postprocess

Idempotent: re-running is a no-op once applied.
"""
import sys
from pathlib import Path

repo = Path(sys.argv[1] if len(sys.argv) > 1 else ".").expanduser()
pkg = repo / "helios_postprocess"
ns = pkg / "neutron_scatter.py"
if not ns.exists():
    sys.exit(f"could not find helios_postprocess/neutron_scatter.py under {repo} "
             "-- pass the repo root")

s = ns.read_text()

# ---- 1. multi-material metrics block ------------------------------------
OLD_BLOCK = '''                scat = comp
                metrics.update({
                    "DSR": scat["dsr_fuel"]["DSR"],
                    "DSR_total": scat["dsr_total"]["DSR"],
                    "rhoR_C_gcm2": scat["rhoR_C_gcm2"],
                    "frac_D": scat["frac_D"], "frac_T": scat["frac_T"],
                    "coeff_gcm2_per_DSR": (metrics["rhoR_hydro_gcm2"]
                                           / scat["dsr_fuel"]["DSR"])
                    if scat["dsr_fuel"]["DSR"] else float("nan"),
                    "dsr_source": "NeSST multi-material (D+T+C, ENDF/B-VIII.0)",
                    "composition_source": "RHW",
                })'''

NEW_BLOCK = '''                scat = comp
                # Self-consistent DSR<->rhoR coefficient: use the SAME areal-
                # density convention on both sides.  The DSR is derived from the
                # traversed (Method-2) rhoR, so the coefficient divides the
                # traversed rhoR -- not the center-to-edge hydro column -- by that
                # DSR.  Total (D+T+C) traversed rhoR / total DSR is the quantity
                # the empirical NIF 20.4 g/cm^2-per-DSR calibrates.
                _cmp = (scat.get("composition_traversed")
                        or scat.get("composition") or {})
                _rD = float(_cmp.get("rhoR_D_gcm2", float("nan")))
                _rT = float(_cmp.get("rhoR_T_gcm2", float("nan")))
                _rC = float(_cmp.get("rhoR_C_gcm2", 0.0))
                _rhoR_tr_fuel = _rD + _rT
                _rhoR_tr_tot = _rhoR_tr_fuel + _rC
                _dsr_fuel = scat["dsr_fuel"]["DSR"]
                _dsr_tot = scat["dsr_total"]["DSR"]
                metrics.update({
                    "DSR": _dsr_fuel,
                    "DSR_total": _dsr_tot,
                    "rhoR_C_gcm2": scat["rhoR_C_gcm2"],
                    "frac_D": scat["frac_D"], "frac_T": scat["frac_T"],
                    "rhoR_traversed_fuel_gcm2": _rhoR_tr_fuel,
                    "rhoR_traversed_total_gcm2": _rhoR_tr_tot,
                    "coeff_gcm2_per_DSR": (_rhoR_tr_tot / _dsr_tot)
                    if _dsr_tot else float("nan"),
                    "coeff_fuel_gcm2_per_DSR": (_rhoR_tr_fuel / _dsr_fuel)
                    if _dsr_fuel else float("nan"),
                    "coeff_convention":
                        "traversed Method-2; total (D+T+C) rhoR / total DSR",
                    "dsr_source": "NeSST multi-material (D+T+C, ENDF/B-VIII.0)",
                    "composition_source": "RHW",
                })'''

# ---- 2. formatter print line --------------------------------------------
OLD_PRINT = '''    if m.get("coeff_gcm2_per_DSR") is not None:
        L.append(f"  first-principles coeff {m['coeff_gcm2_per_DSR']:.2f}"
                 " g/cm^2 per DSR   (empirical NIF 20.4)")'''

NEW_PRINT = '''    if m.get("coeff_gcm2_per_DSR") is not None:
        L.append(f"  first-principles coeff {m['coeff_gcm2_per_DSR']:.2f}"
                 " g/cm^2 per DSR   (empirical NIF 20.4)")
        conv = m.get("coeff_convention")
        if conv:
            L.append(f"                         [{conv}]")
        if m.get("coeff_fuel_gcm2_per_DSR") is not None:
            L.append(f"  first-principles coeff {m['coeff_fuel_gcm2_per_DSR']:.2f}"
                     " g/cm^2 per DSR_fuel   (fuel D+T only, self-consistent)")'''

already = ('rhoR_traversed_total_gcm2' in s and 'coeff_convention' in s)
if already:
    print("neutron_scatter.py: coefficient convention already patched -- skipped")
    sys.exit(0)

if OLD_BLOCK not in s:
    sys.exit("ERROR: could not locate the multi-material metrics block verbatim.\n"
             "Your neutron_report may have drifted -- patch not applied. Send me\n"
             "the current lines around 'coeff_gcm2_per_DSR' and I'll re-cut it.")
if OLD_PRINT not in s:
    sys.exit("ERROR: could not locate the formatter print line verbatim.\n"
             "Patch not applied; send me the current '_format_neutron_block' coeff line.")

s = s.replace(OLD_BLOCK, NEW_BLOCK, 1)
s = s.replace(OLD_PRINT, NEW_PRINT, 1)
ns.write_text(s)
print("neutron_scatter.py: coefficient made self-consistent (traversed convention)")
print("  - primary  coeff = total (D+T+C) traversed rhoR / total DSR   (vs NIF 20.4)")
print("  - fuel-only coeff = (D+T) traversed rhoR / DSR_fuel")
print("done. verify: git diff helios_postprocess/neutron_scatter.py")
