#!/usr/bin/env python
"""
patch_flux_limiter_round2.py — fix/complete Task 1

Against live file text verified via repr() output:
  1c. rhw_parser.parse()   — call _parse_flux_limiter + pass to RHWConfiguration(...)
  1e. data_builder.py      — wire flux_limiter from rhw_config (redo with 8-space indent)
  1f. icf_output.py        — add Flux limiter line after Power multiplier, before pulse shape

Idempotent: each edit checks for "already applied" sentinels before writing.

Run on MacBook from repo root: python patch_flux_limiter_round2.py
"""
from pathlib import Path

ROOT = Path('~/Codes/helios_postprocessor/helios_postprocess').expanduser()
print(f'Repo: {ROOT}\n')

def read(p): return Path(p).read_text()
def write(p, s): Path(p).write_text(s); print(f'  [WROTE]  {p.name}')

def apply(src, old, new, sentinel, label):
    """Apply old -> new. If sentinel already present in src, skip as already applied."""
    if sentinel in src and old not in src:
        print(f'  [ALREADY] {label}')
        return src, False
    if old not in src:
        print(f'  [MISS]    {label}: OLD pattern not found')
        return src, False
    print(f'  [OK]      {label}')
    return src.replace(old, new, 1), True


# ── 1c. rhw_parser.parse(): call and pass flux_limiter ───────────────────────
print('=== 1c. rhw_parser.parse() ===')
p = ROOT / 'rhw_parser.py'
src = read(p); orig = src

# 1c-part-1: add _parse_flux_limiter call after _parse_alpha_transport
OLD = '        alpha_local, alpha_nonlocal = self._parse_alpha_transport(lines)\n        \n        # Extract drive temperature data\n'
NEW = '        alpha_local, alpha_nonlocal = self._parse_alpha_transport(lines)\n        flux_enabled, flux_value = self._parse_flux_limiter(lines)\n        \n        # Extract drive temperature data\n'
src, _ = apply(src, OLD, NEW,
               sentinel='flux_enabled, flux_value = self._parse_flux_limiter(lines)',
               label='parse(): call _parse_flux_limiter')

# 1c-part-2: add flux_limiter/flux_limiter_enabled to RHWConfiguration(...)
OLD = '            alpha_deposition_local=alpha_local,\n            alpha_deposition_nonlocal=alpha_nonlocal,\n        )\n'
NEW = '            alpha_deposition_local=alpha_local,\n            alpha_deposition_nonlocal=alpha_nonlocal,\n            flux_limiter=flux_value,\n            flux_limiter_enabled=flux_enabled,\n        )\n'
src, _ = apply(src, OLD, NEW,
               sentinel='flux_limiter=flux_value,',
               label='parse(): pass flux_limiter to RHWConfiguration')

if src != orig:
    write(p, src)


# ── 1e. data_builder.py: wire flux_limiter from rhw_config ───────────────────
print('\n=== 1e. data_builder.py ===')
p = ROOT / 'data_builder.py'
src = read(p); orig = src

# Anchor on laser_peak_end_ns (last per-beam field in the block, stable)
OLD = '        data.laser_peak_end_ns         = rhw_config.laser_peak_end_ns\n'
NEW = ('        data.laser_peak_end_ns         = rhw_config.laser_peak_end_ns\n'
       '        # Flux limiter (from .rhw)\n'
       '        data.flux_limiter         = getattr(rhw_config, \'flux_limiter\', 0.0)\n'
       '        data.flux_limiter_enabled = getattr(rhw_config, \'flux_limiter_enabled\', False)\n')
src, _ = apply(src, OLD, NEW,
               sentinel='data.flux_limiter_enabled = getattr(rhw_config,',
               label='build_run_data: flux_limiter wiring')

if src != orig:
    write(p, src)


# ── 1f. icf_output.py: flux limiter line in LASER CONFIGURATION ──────────────
print('\n=== 1f. icf_output.py ===')
p = ROOT / 'icf_output.py'
src = read(p); orig = src

# Anchor: Power multiplier line + # Pulse shape line (line 193 + 194 exactly).
# Flux limiter goes between them — logically part of static config, not pulse shape.
OLD = ("            _a(self._metric('Power multiplier', d.laser_power_multiplier,    '',    fmt='.4f'))\n"
       "            # Pulse shape\n")
NEW = ("            _a(self._metric('Power multiplier', d.laser_power_multiplier,    '',    fmt='.4f'))\n"
       "            # Flux limiter: bypass _metric (0.06 is a valid value; _metric renders 0.0 as '—')\n"
       "            if getattr(d, 'flux_limiter_enabled', False):\n"
       "                _a(f\"  {'Flux limiter (f)':<36s} {d.flux_limiter:>10.3f}\")\n"
       "            elif getattr(d, 'flux_limiter', 0.0) > 0:\n"
       "                _a(f\"  {'Flux limiter (f)':<36s} {d.flux_limiter:>10.3f}  (flag off)\")\n"
       "            else:\n"
       "                _a(f\"  {'Flux limiter (f)':<36s} {'(not set)':>10s}\")\n"
       "            # Pulse shape\n")
src, _ = apply(src, OLD, NEW,
               sentinel="'Flux limiter (f)':<36s",
               label='LASER CONFIGURATION: flux limiter line')

if src != orig:
    write(p, src)


print('\nDone. Expected sequence in Olson_PDD_26b_summary.txt after re-run:')
print('  LASER CONFIGURATION (beam 1)')
print('  ------------------------------------------------------------')
print('    Wavelength                           0.351 um')
print('    Focus position                       0.2000 cm')
print('    Half-cone angle                      23.00 deg')
print('    Spot radius                          0.0800 cm  (Gaussian)')
print('    Power multiplier                     1.4000')
print('    Flux limiter (f)                          0.060')
print('    Pulse shape (beam 1)')
print('      Foot power                             25.0 TW  (0.00 – 2.00 ns)')
print('      ...')