"""
Plot LaserEnTimeIntg spatial deposition profile at several timesteps.
Shows cumulative absorbed energy per zone vs radius, revealing where
in the corona laser energy is deposited.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc

base = sys.argv[1] if len(sys.argv) > 1 else \
    '/Users/tommehlhorn/Sims/Xcimer/Olson_PDD/Olson_PDD_2021_01a/Olson_PDD_2021_01a'

exo = base + '.exo'
print(f'Loading {exo}')
ds = nc.Dataset(exo)

t_ns   = np.array(ds.variables['time_whole'][:]) * 1e9
zbnd   = np.array(ds.variables['zone_boundaries'][:])          # (425, 441) cm
zmas   = np.array(ds.variables['zone_mass'][:])                # (425, 440) g
EnLas  = np.array(ds.variables['LaserEnTimeIntg'][:]) * 1e-6   # (425, 440) -> MJ
LasPwr = np.array(ds.variables['LaserPwrSrc'][:])              # (425, 440) W/cm3
rho    = np.array(ds.variables['mass_density'][:])             # (425, 440) g/cc
ri     = np.array(ds.variables['Indices at region interfaces'][:])  # (425, 4)
dep_cum = np.array(ds.variables['EnLaserDepositedTimeIntg'][:]) * 1e-6  # MJ
deliv   = np.array(ds.variables['LaserEnDeliveredTimeInt'][:]) * 1e-6   # MJ
ds.close()

# Zone centers (cm) and radii (um)
zc_cm = 0.5 * (zbnd[:, :-1] + zbnd[:, 1:])   # (425, 440)
zc_um = zc_cm * 1e4                           # um

# Target name
import os
name = os.path.basename(base)
R0_cm = zbnd[0, -1]  # initial outer radius (cm)

# Pick timesteps to plot
# Early foot: ~2 ns, ~4 ns
# Mid ramp: ~6 ns, ~8 ns
# Peak drive: ~10 ns, ~12 ns
target_times = [2.0, 4.0, 5.5, 7.0, 9.0, 11.0, 12.5]
colors = ['#3B8BD4', '#3B8BD4', '#5DCAA5', '#5DCAA5', '#EF9F27', '#EF9F27', '#D85A30']
linestyles = ['--', '-', '--', '-', '--', '-', '-']
alphas     = [0.5, 0.9, 0.5, 0.9, 0.5, 0.9, 1.0]

# Find nearest index for each target time
tidxs = [np.argmin(np.abs(t_ns - tt)) for tt in target_times]
actual_times = [t_ns[i] for i in tidxs]

# ── Figure layout: 2 rows ────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32)
ax1 = fig.add_subplot(gs[0, :])   # top: incremental deposition per zone
ax2 = fig.add_subplot(gs[1, 0])   # bottom-left: power density profile
ax3 = fig.add_subplot(gs[1, 1])   # bottom-right: cumulative absorbed energy history

fig.suptitle(f'{name} — Laser deposition spatial profiles', fontsize=13, y=0.98)

# ── Panel 1: Cumulative energy deposited per zone at several times ───────────
# Normalize by zone volume to get J/cm3 → energy density
for k, (idx, tt, col, ls, al) in enumerate(zip(
        tidxs, actual_times, colors, linestyles, alphas)):

    r_um = zc_um[idx, :]
    # Energy in MJ per zone. Normalize by initial zone mass for specific deposition (MJ/g)
    spec_dep = EnLas[idx, :] / np.where(zmas[idx, :] > 0, zmas[idx, :], 1e-30)

    # Outer ablator region only (zones where R < R0 and rho still significant)
    outer = r_um < R0_cm * 1e4 * 1.05

    lbl = f'{tt:.1f} ns'
    ax1.plot(r_um[outer], spec_dep[outer] * 1e6,  # MJ/g → J/g
             color=col, ls=ls, alpha=al, lw=1.3, label=lbl)

# Mark initial region boundaries
for col_idx, node_col in enumerate(range(ri.shape[1])):
    node = min(int(ri[0, col_idx]), zbnd.shape[1]-1)
    r_bnd_um = zbnd[0, node] * 1e4
    lbl_names = ['HS / vapor', 'DT ice', 'foam', 'CH skin']
    ax1.axvline(r_bnd_um, color='gray', lw=0.7, ls=':', alpha=0.5)
    ax1.text(r_bnd_um + 5, ax1.get_ylim()[1] if ax1.get_ylim()[1] > 0 else 1,
             lbl_names[col_idx], fontsize=7, color='gray', rotation=90, va='top')

ax1.set_xlabel('Radius (µm)', fontsize=11)
ax1.set_ylabel('Specific energy deposited (J/g)', fontsize=11)
ax1.set_title('Cumulative laser energy per zone (normalized by zone mass)', fontsize=11)
ax1.legend(fontsize=9, ncol=4, loc='upper left')
ax1.set_xlim(left=0)
ax1.grid(True, alpha=0.25)

# ── Panel 2: Instantaneous power density at several times ───────────────────
for k, (idx, tt, col, ls, al) in enumerate(zip(
        tidxs, actual_times, colors, linestyles, alphas)):

    r_um = zc_um[idx, :]
    pwr  = LasPwr[idx, :]  # W/cm3

    # Limit to relevant region
    outer = (r_um < R0_cm * 1e4 * 1.05) & (pwr > 0)
    if not np.any(outer):
        outer = r_um < R0_cm * 1e4 * 1.05

    lbl = f'{tt:.1f} ns'
    ax2.semilogy(r_um[outer], np.maximum(pwr[outer], 1e-10) * 1e-12,
                 color=col, ls=ls, alpha=al, lw=1.3, label=lbl)

# Mark region boundaries at reference time
for col_idx in range(ri.shape[1]):
    node = min(int(ri[0, col_idx]), zbnd.shape[1]-1)
    r_bnd_um = zbnd[0, node] * 1e4
    ax2.axvline(r_bnd_um, color='gray', lw=0.7, ls=':', alpha=0.5)

ax2.set_xlabel('Radius (µm)', fontsize=11)
ax2.set_ylabel('Laser power density (TW/cm³)', fontsize=11)
ax2.set_title('Instantaneous deposition power density', fontsize=11)
ax2.legend(fontsize=8, ncol=2)
ax2.set_xlim(left=0)
ax2.grid(True, alpha=0.25, which='both')

# ── Panel 3: Cumulative absorption history ───────────────────────────────────
# Delivered, deposited, and absorption fraction vs time
laser_on = deliv > 1e-6  # mask where laser is actually on
t_plot = t_ns[laser_on]
d_plot = deliv[laser_on]
p_plot = dep_cum[laser_on]
abs_pct = np.where(d_plot > 1e-6, 100 * p_plot / d_plot, 0)

ax3b = ax3.twinx()
ax3.plot(t_plot, d_plot, color='#3B8BD4', lw=1.5, label='Delivered (MJ)')
ax3.plot(t_plot, p_plot, color='#EF9F27', lw=1.5, label='Deposited (MJ)')
ax3b.plot(t_plot, abs_pct, color='#D85A30', lw=1.2, ls='--', label='Absorption %')

ax3.set_xlabel('Time (ns)', fontsize=11)
ax3.set_ylabel('Cumulative energy (MJ)', fontsize=11)
ax3b.set_ylabel('Absorption fraction (%)', fontsize=11, color='#D85A30')
ax3b.tick_params(axis='y', labelcolor='#D85A30')
ax3.set_title('Cumulative delivered vs deposited energy', fontsize=11)

# Mark vertical lines at the same timesteps
for tt, col in zip(actual_times, colors):
    ax3.axvline(tt, color=col, lw=0.7, ls=':', alpha=0.5)

lines1, labs1 = ax3.get_legend_handles_labels()
lines2, labs2 = ax3b.get_legend_handles_labels()
ax3.legend(lines1 + lines2, labs1 + labs2, fontsize=9, loc='upper left')
ax3.grid(True, alpha=0.25)

# ── Fix panel 1 y-axis after data is plotted ────────────────────────────────
# Re-annotate region boundaries with correct ylim
ymax1 = ax1.get_ylim()[1]
for col_idx in range(ri.shape[1]):
    node = min(int(ri[0, col_idx]), zbnd.shape[1]-1)
    r_bnd_um = zbnd[0, node] * 1e4
    lbl_names = ['HS bnd', 'ice bnd', 'foam bnd', 'outer']
    ax1.text(r_bnd_um + 5, ymax1 * 0.95,
             lbl_names[col_idx], fontsize=7, color='gray', rotation=90, va='top')

# Save
outpath = base + '_laser_deposition.pdf'
fig.savefig(outpath, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {outpath}')

# Also print a summary table of deposition peak location at each timestep
print('\nDeposition peak radius vs time:')
print(f'  {"Time(ns)":>8}  {"R_peak(um)":>11}  {"R_target(um)":>13}  {"R_peak/R_target":>16}  {"Abs%":>6}')
for idx, tt in zip(tidxs, actual_times):
    pwr = LasPwr[idx, :]
    r_um = zc_um[idx, :]
    r_outer_um = zbnd[idx, -1] * 1e4
    if pwr.max() > 0:
        peak_zone = np.argmax(pwr)
        r_peak = r_um[peak_zone]
        ratio = r_peak / r_outer_um if r_outer_um > 0 else 0
        dv = deliv[idx]; dp = dep_cum[idx]
        pct = 100*dp/dv if dv > 1e-6 else 0
        print(f'  {tt:8.2f}  {r_peak:11.1f}  {r_outer_um:13.1f}  {ratio:16.3f}  {pct:6.1f}')
