#!/usr/bin/env python
"""
patch_flux_limiter.py — Task 1 (partial)

Edits:
  1a. rhw_parser.py        — add flux_limiter fields to RHWConfiguration
  1b. rhw_parser.py        — add _parse_flux_limiter() method
  1d. data_builder.py      — add flux_limiter scalars to ICFRunData
  1e. data_builder.py      — wire from rhw_config in build_run_data
  1f. icf_output.py        — add Flux limiter line to LASER CONFIGURATION block

Not yet edited:
  1c. rhw_parser.parse()   — needs parse()'s body (see README at end).

Run: python patch_flux_limiter.py
"""
from pathlib import Path

ROOT = Path('~/Codes/helios_postprocessor/helios_postprocess').expanduser()
print(f'Repo: {ROOT}\n')

def apply(src, old, new, label):
    if old not in src:
        if new.split('\n')[1].strip() in src:
            print(f'  [ALREADY] {label}')
        else:
            print(f'  [MISS]   {label}: OLD pattern not found')
        return src, False
    out = src.replace(old, new, 1)
    print(f'  [OK]     {label}')
    return out, True

def read(p): return Path(p).read_text()
def write(p, s): Path(p).write_text(s); print(f'  [WROTE]  {p.name}')

# ── 1a. rhw_parser.py: RHWConfiguration fields ───────────────────────────────
print('=== 1a. rhw_parser.py: RHWConfiguration fields ===')
p = ROOT / 'rhw_parser.py'
src = read(p); orig = src
OLD = """    alpha_deposition_local: bool = False     # Use alpha deposition = 1
    alpha_deposition_nonlocal: bool = False  # Use non alpha deposition = 1"""
NEW = """    alpha_deposition_local: bool = False     # Use alpha deposition = 1
    alpha_deposition_nonlocal: bool = False  # Use non alpha deposition = 1
    flux_limiter: float = 0.0                # Flux limiter mult. (from .rhw; 0 if off)
    flux_limiter_enabled: bool = False       # Use flux limiter = 1"""
src, _ = apply(src, OLD, NEW, 'RHWConfiguration.flux_limiter fields')

# ── 1b. rhw_parser.py: _parse_flux_limiter method ────────────────────────────
print('\n=== 1b. rhw_parser.py: _parse_flux_limiter method ===')
OLD = "    def _parse_eos_models(self, lines: list) -> list:"
NEW = '''    def _parse_flux_limiter(self, lines: list) -> tuple:
        """
        Parse electron thermal flux limiter from the .rhw file.

        The .rhw has four copies of the flux-limiter block — one per region
        (DT vapor, DT ice, CH skin, foam). Pattern:

            Use flux limiter       = 1
            Flux limiter mult.     = 0.06

        Normally identical across regions. Take the first (enabled, value)
        pair seen. Warn once if a later block has a different value.

        Returns
        -------
        (enabled, value) : (bool, float)
            (False, 0.0) if no flux-limiter line is present.
        """
        import warnings as _w
        enabled_first = None
        value_first = None
        warned = False
        for line in lines:
            lstrip = line.strip().lower()
            if lstrip.startswith('use flux limiter'):
                try:
                    v = int(line.split('=')[-1].strip())
                except (ValueError, IndexError):
                    continue
                if enabled_first is None:
                    enabled_first = (v == 1)
            elif lstrip.startswith('flux limiter mult'):
                try:
                    f = float(line.split('=')[-1].strip())
                except (ValueError, IndexError):
                    continue
                if value_first is None:
                    value_first = f
                elif abs(f - value_first) > 1e-9 and not warned:
                    _w.warn(
                        f'RHW flux limiter varies across regions '
                        f'(first={value_first}, later={f}); reporting first.',
                        stacklevel=2)
                    warned = True
        return (bool(enabled_first) if enabled_first is not None else False,
                float(value_first) if value_first is not None else 0.0)

    def _parse_eos_models(self, lines: list) -> list:'''
src, _ = apply(src, OLD, NEW, '_parse_flux_limiter method')

if src != orig:
    write(p, src)

# ── 1d. data_builder.py: ICFRunData fields ───────────────────────────────────
print('\n=== 1d. data_builder.py: ICFRunData flux_limiter fields ===')
p = ROOT / 'data_builder.py'
src = read(p); orig = src
OLD = "        self.laser_pulse_duration_ns: float = 0.0"
NEW = """        self.laser_pulse_duration_ns: float = 0.0
        self.flux_limiter: float = 0.0              # from .rhw (Flux limiter mult.)
        self.flux_limiter_enabled: bool = False     # from .rhw (Use flux limiter)"""
src, _ = apply(src, OLD, NEW, 'ICFRunData.flux_limiter fields')

# ── 1e. data_builder.py: wire from rhw_config ────────────────────────────────
print('\n=== 1e. data_builder.py: wire flux_limiter from rhw_config ===')
OLD = "            data.laser_power_multiplier    = rhw_config.laser_power_multiplier"
NEW = """            data.laser_power_multiplier    = rhw_config.laser_power_multiplier
            # Flux limiter (from .rhw)
            data.flux_limiter         = getattr(rhw_config, 'flux_limiter', 0.0)
            data.flux_limiter_enabled = getattr(rhw_config, 'flux_limiter_enabled', False)"""
src, _ = apply(src, OLD, NEW, 'build_run_data flux_limiter wiring')

if src != orig:
    write(p, src)

# ── 1f. icf_output.py: LASER CONFIGURATION block ─────────────────────────────
print('\n=== 1f. icf_output.py: LASER CONFIGURATION flux limiter line ===')
p = ROOT / 'icf_output.py'
src = read(p); orig = src
OLD = """            _a(self._metric('Power multiplier', d.laser_power_multiplier,    '',    fmt='.4f'))
            _a('')"""
NEW = """            _a(self._metric('Power multiplier', d.laser_power_multiplier,    '',    fmt='.4f'))
            # Flux limiter: bypass _metric (0.06 is a valid value; _metric renders 0.0 as '—')
            if getattr(d, 'flux_limiter_enabled', False):
                _a(f"  {'Flux limiter (f)':<36s} {d.flux_limiter:>10.3f}")
            elif getattr(d, 'flux_limiter', 0.0) > 0:
                _a(f"  {'Flux limiter (f)':<36s} {d.flux_limiter:>10.3f}  (flag off)")
            else:
                _a(f"  {'Flux limiter (f)':<36s} {'(not set)':>10s}")
            _a('')"""
src, _ = apply(src, OLD, NEW, 'LASER CONFIGURATION: flux limiter line')

if src != orig:
    write(p, src)

print('\nDone.')
print('\nRemaining (needs parse() body to write):')
print('  1c. rhw_parser.parse() — call self._parse_flux_limiter(lines)')
print('\nRun:')
print('  sed -n "68,155p" ~/Codes/helios_postprocessor/helios_postprocess/rhw_parser.py')
print('and paste so the parse() integration patch can be written.')