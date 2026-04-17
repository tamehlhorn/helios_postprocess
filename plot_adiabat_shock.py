"""
Standalone diagnostic plots: adiabat history and first shock in DT ice.
Run against any simulation with a DT ice layer (ri[:,0] to ri[:,1]).

Usage:
    python3 plot_adiabat_shock.py <base_path> [--fuel_col 0] [--foam_col 1]

Example:
    python3 plot_adiabat_shock.py \
        ~/Sims/Xcimer/Olson_PDD/Olson_PDD_2021_01a/Olson_PDD_2021_01a
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from helios_postprocess.core import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.icf_analysis import ICFAnalyzer
from helios_postprocess.adiabat_history import compute_adiabat_history
from helios_postprocess.pressure_gradients import analyze_first_shock

# ── CLI ──────────────────────────────────────────────────────────────────────
base = sys.argv[1] if len(sys.argv) > 1 else None
if base is None:
    print("Usage: python3 plot_adiabat_shock.py <base_path>")
    sys.exit(1)

fuel_col = 0   # ri column for inner DT ice boundary (hot-spot edge)
foam_col = 1   # ri column for outer DT ice boundary

# Optional published reference values (edit as needed)
published = {
    'adiabat':       7.4,    # set to None to skip
    'peak_velocity': 370.0,  # km/s
    'label':         'Olson 2021 (2D HYDRA)',
}

# ── Load and analyze ─────────────────────────────────────────────────────────
print(f"Loading {base}.exo ...")
run  = HeliosRun(base + '.exo')
data = build_run_data(run)
analyzer = ICFAnalyzer(data)
analyzer.analyze_drive_phase()
analyzer.analyze_stagnation_phase()
analyzer.analyze_burn_phase()
analyzer.compute_performance_metrics()

t_ns = data.time * 1e9
ri   = data.region_interfaces_indices
name = os.path.basename(base)

print(f"ri[0,:] = {ri[0,:]}")
print(f"Peak velocity index: {data.peak_velocity_index}, t={t_ns[data.peak_velocity_index]:.3f} ns")
print(f"Stagnation time: {data.stag_time*1e9:.3f} ns")
print(f"Adiabat at peak v: {data.adiabat_mass_averaged_ice:.3f}")

# ── Adiabat history ──────────────────────────────────────────────────────────
adiabat_result = compute_adiabat_history(
    ion_pressure=data.ion_pressure,
    rad_pressure=data.rad_pressure,
    mass_density=data.mass_density,
    zone_mass=data.zone_mass,
    region_interfaces_indices=ri,
    time=t_ns,
    fuel_region_col=fuel_col,
    fuel_outer_col=foam_col,
)

# ── First shock in DT ice ─────────────────────────────────────────────────────
ice_inner_node = int(ri[0, fuel_col])
ice_outer_node = int(ri[0, foam_col])
ice_inner_zone = ice_inner_node
ice_outer_zone = ice_outer_node - 1
total_pressure = data.ion_pressure + data.rad_pressure

shock_result = analyze_first_shock(
    pressure=total_pressure,
    zone_boundaries=data.zone_boundaries,
    density=data.mass_density,
    time=t_ns,
    ablator_outer_zone=ice_outer_zone,
    fuel_inner_zone=ice_inner_node,
    search_inner_zone=ice_inner_zone,
    dP_dr_threshold=1e9,
)

print(f"\nShock breakout time:     {shock_result['breakout_time_ns']}")
print(f"Shock breakout pressure: {shock_result['breakout_pressure_Gbar']}")

# ── Timing markers ───────────────────────────────────────────────────────────
t_foot_end   = data.laser_foot_end_ns   if hasattr(data, 'laser_foot_end_ns')   else 5.0
t_peak_start = data.laser_peak_start_ns if hasattr(data, 'laser_peak_start_ns') else 9.0
t_peak_v     = t_ns[data.peak_velocity_index]
# stag_time and bang_time are set by ICFAnalyzer on the data object
t_stag = data.stag_time * 1e9 if hasattr(data, "stag_time") and data.stag_time > 0 else t_ns[data.peak_velocity_index] + 1.0
t_bang = data.bang_time * 1e9 if hasattr(data, "bang_time") and data.bang_time > 0 else t_stag + 0.1

# ── FIGURE 1: Adiabat history ────────────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(9, 5))
ax2 = ax1.twinx()

# Mask to drive/implosion window only
mask = (t_ns >= 4.0) & (t_ns <= t_stag + 0.3)
t_plot = t_ns[mask]
a_plot = adiabat_result['adiabat'][mask]
amin_plot = adiabat_result['adiabat_min'][mask]
p_plot = adiabat_result['P_fuel_mean_Gbar'][mask]

ax1.plot(t_plot, a_plot,    color='royalblue',  lw=2,   label=r'$\langle\alpha\rangle$ mass-avg (DT ice)')
ax1.plot(t_plot, amin_plot, color='steelblue',  lw=1.5, ls='--', label=r'$\alpha_{min}$ (densest zones)')
ax2.plot(t_plot, p_plot,    color='firebrick',  lw=1.5, ls=':',  label=r'$\langle P \rangle$ DT ice (Gbar)', alpha=0.8)

# Reference line
if published.get('adiabat') is not None:
    ax1.axhline(published['adiabat'], color='darkorange', lw=1.5, ls='-.',
                label=f"Published α={published['adiabat']} ({published['label']})")

# Timing markers
for t_mark, label, col in [
    (t_foot_end,   'Foot end',    'gray'),
    (t_peak_start, 'Peak start',  'gray'),
    (t_peak_v,     'Peak $v_{imp}$', 'navy'),
    (t_stag,       'Stagnation',  'darkred'),
]:
    ax1.axvline(t_mark, color=col, lw=1, ls=':', alpha=0.7)
    ax1.text(t_mark+0.05, 0.97, label, fontsize=7, color=col,
             rotation=90, va='top', transform=ax1.get_xaxis_transform())

ax1.set_xlabel('Time (ns)', fontsize=12)
ax1.set_ylabel(r'Fuel adiabat $\alpha$', fontsize=12, color='royalblue')
ax2.set_ylabel(r'Mean pressure in DT ice (Gbar)', fontsize=11, color='firebrick')
ax2.set_yscale('log')
ax2.set_ylim(1e-4, None)
ax1.tick_params(axis='y', labelcolor='royalblue')
ax2.tick_params(axis='y', labelcolor='firebrick')
ax1.set_xlim(4.0, t_stag + 0.3)
ax1.set_ylim(bottom=0)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper left')
ax1.set_title(f'{name} — DT Ice Adiabat History', fontsize=13)
ax1.grid(True, alpha=0.3)

fig1.tight_layout()
out1 = base + '_adiabat_history.pdf'
fig1.savefig(out1, dpi=150)
print(f"\nSaved: {out1}")

# ── FIGURE 2: First shock in DT ice ─────────────────────────────────────────
# Filter to valid shock detections in drive window
shock_mask = (~np.isnan(shock_result['shock_radius'])) & (t_ns >= 4.5) & (t_ns <= t_stag)
t_sh  = t_ns[shock_mask]
r_sh  = shock_result['shock_radius'][shock_mask] * 1e4   # cm -> um
p_sh  = shock_result['shock_pressure_Gbar'][shock_mask]
m_sh  = shock_result['mach_number'][shock_mask]
v_sh  = shock_result['shock_velocity'][shock_mask] * 100  # cm/ns -> km/s

fig2, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True,
                           gridspec_kw={'hspace': 0.08})

# Panel 1: shock radius
ax = axes[0]
ax.scatter(t_sh, r_sh, s=8, color='steelblue', alpha=0.6, label='Shock radius')
# Overlay DT ice boundaries
r_ice_inner = data.zone_boundaries[:, ice_inner_node] * 1e4
r_ice_outer = data.zone_boundaries[:, ice_outer_node] * 1e4
ax.plot(t_ns, r_ice_inner, 'k--', lw=1, label='DT ice inner boundary')
ax.plot(t_ns, r_ice_outer, 'k:',  lw=1, label='DT ice outer boundary')
ax.set_ylabel('Radius (µm)', fontsize=11)
ax.set_ylim(bottom=0)
ax.legend(fontsize=8, loc='upper right')
ax.set_title(f'{name} — First Shock in DT Ice', fontsize=13)
ax.grid(True, alpha=0.3)

# Panel 2: shock pressure
ax = axes[1]
ax.scatter(t_sh, p_sh, s=8, color='firebrick', alpha=0.6, label='Post-shock P')
ax.set_ylabel('Post-shock pressure (Gbar)', fontsize=11)
ax.set_yscale('log')
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, alpha=0.3, which='both')

# Panel 3: Mach number
ax = axes[2]
valid_mach = ~np.isnan(m_sh)
ax.scatter(t_sh[valid_mach], m_sh[valid_mach], s=8, color='seagreen', alpha=0.6, label='Mach number')
ax.axhline(1.0, color='gray', lw=1, ls='--')
ax.set_ylabel('Mach number', fontsize=11)
ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylim(0.9, None)
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, alpha=0.3)

# Timing markers on all panels
for ax in axes:
    for t_mark, col in [
        (t_foot_end,   'gray'),
        (t_peak_start, 'gray'),
        (t_peak_v,     'navy'),
    ]:
        ax.axvline(t_mark, color=col, lw=1, ls=':', alpha=0.6)

xlim_end = min(t_stag + 0.1, t_ns[-1])
for ax in axes:
    ax.set_xlim(4.5, xlim_end)

fig2.tight_layout()
out2 = base + '_shock_history.pdf'
fig2.savefig(out2, dpi=150)
print(f"Saved: {out2}")

plt.close('all')
print("\nDone.")
