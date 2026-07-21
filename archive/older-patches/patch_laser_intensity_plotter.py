#!/usr/bin/env python
"""
patch_laser_intensity_plotter.py -- Task 2 Stage C.2

Add _plot_laser_intensity() method to ICFPlotter in icf_plotting.py,
plus register call in create_full_report. Three pages:
  1. Log10 I(r, t) 2D heatmap (Method 2) with r_crit overlay
  2. Time histories: P_laser, I_grid_outer, I_at_crit, I_peak_coronal
  3. Method 1 vs Method 2 cross-check: radial profiles at peak power

Inserted between _plot_laser_deposition (line 1040) and _plot_laser_power
(line 1059), so the PDF reads: deposition -> intensity -> power delivered.

Graceful no-op: if self.data._laser_intensity_arrays is None (Stage B
populated it with None because inputs were missing), emit a placeholder
page and return.

Idempotent via sentinel checks.

Run from repo root:
    python patch_laser_intensity_plotter.py
"""
from pathlib import Path

TARGET = Path('~/Codes/helios_postprocessor/helios_postprocess/icf_plotting.py').expanduser()
print(f'Target: {TARGET}\n')


def apply(src, old, new, sentinel, label):
    if sentinel in src and old not in src:
        print(f'  [ALREADY] {label}')
        return src, False
    if old not in src:
        print(f'  [MISS]    {label}: OLD pattern not found')
        return src, False
    print(f'  [OK]      {label}')
    return src.replace(old, new, 1), True


src = TARGET.read_text()
orig = src


# ====================================================================
# EDIT 1: Insert _plot_laser_intensity method
# Anchor: _plot_laser_power (which comes right after _plot_laser_deposition)
# ====================================================================
print('=== Edit 1: insert _plot_laser_intensity method ===\n')

ANCHOR = '    def _plot_laser_power(self, pdf):'

NEW_METHOD = '''    def _plot_laser_intensity(self, pdf):
        """
        Three pages of laser intensity diagnostics.

        Requires self.data._laser_intensity_arrays populated by
        ICFAnalyzer.analyze_laser_intensity(). If absent, emits a
        placeholder page and returns.
        """
        import numpy as _np
        import matplotlib.pyplot as _plt

        arrays = getattr(self.data, '_laser_intensity_arrays', None)
        if arrays is None:
            fig, ax = _plt.subplots(figsize=(8.5, 11))
            ax.text(0.5, 0.5,
                    "Laser intensity: required inputs missing\\n"
                    "(see run log for details)",
                    ha='center', va='center',
                    transform=ax.transAxes, fontsize=14)
            ax.set_axis_off()
            pdf.savefig(fig); _plt.close(fig)
            return

        t_ns = self.data.time  # already in ns
        zbnd = self.data.zone_boundaries
        zcen = 0.5 * (zbnd[:, :-1] + zbnd[:, 1:])
        r_um = zcen * 1e4  # cm -> um for plot axes

        I1 = arrays['I1']
        I2 = arrays['I2']
        alpha = arrays['alpha_zone']
        r_crit = arrays['r_crit']
        I_peak_coronal_vs_t = arrays['I_peak_coronal_vs_t']

        P_laser = self.data.laser_power_on_target
        I_grid_outer = self.data.I_grid_outer_history
        I_at_crit = self.data.I_at_crit_history
        t_peak_power = self.data.t_peak_power_ns

        # ---- PAGE 1: log10 I2(r, t) heatmap with r_crit overlay -------------
        fig, ax = _plt.subplots(figsize=(11, 8.5))
        with _np.errstate(invalid='ignore', divide='ignore'):
            I2_pos = _np.where((I2 > 0) & _np.isfinite(I2), I2, _np.nan)
            log_I = _np.log10(I2_pos)

        # Color range: 10^11 W/cm^2 to max, clipped at 10^17
        vmin = 11.0
        try:
            vmax = min(17.0, _np.nanmax(log_I))
            if not _np.isfinite(vmax) or vmax <= vmin:
                vmax = 16.0
        except Exception:
            vmax = 16.0

        # Use pcolormesh with time on x, radius on y (one mesh per frame)
        # For performance with 1131 timesteps, reduce to a uniform mesh
        T_mesh, R_mesh = _np.meshgrid(t_ns, _np.arange(log_I.shape[1]))
        # But we want physical radius as y-axis; each row is a timestep so
        # r varies per timestep. Use imshow on median-radius grid for speed.
        # Simple approach: imshow on (time, zone_index) with r_crit overlaid
        # in zone-index units too.
        im = ax.imshow(log_I.T,
                       origin='lower', aspect='auto',
                       extent=[t_ns[0], t_ns[-1], 0, log_I.shape[1]],
                       cmap='magma', vmin=vmin, vmax=vmax,
                       interpolation='nearest')
        cb = _plt.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(r'log$_{10}$ I [W/cm$^2$]')

        # Overlay r_crit in zone-index space
        r_crit_zone = _np.full_like(r_crit, _np.nan, dtype=float)
        for tt in range(len(r_crit)):
            if _np.isfinite(r_crit[tt]):
                r_crit_zone[tt] = _np.argmin(_np.abs(zcen[tt] - r_crit[tt]))
        ax.plot(t_ns, r_crit_zone, color='cyan', lw=1.2, label='r_crit')

        if _np.isfinite(t_peak_power):
            ax.axvline(t_peak_power, color='white', ls='--', lw=0.8,
                       alpha=0.7, label='t(peak P)')

        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('Zone index (radial)')
        ax.set_title(f'Laser intensity I(r, t)  --  Method 2 (Beer-Lambert)\\n'
                     f'File: {getattr(self.data, "filename", "")}')
        ax.legend(loc='upper right', fontsize=9)
        pdf.savefig(fig); _plt.close(fig)

        # ---- PAGE 2: time histories -----------------------------------------
        fig, axes = _plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)

        ax = axes[0]
        if P_laser is not None:
            ax.plot(t_ns, P_laser * 1e-12, color='black', lw=1.2,
                    label='Laser power on target')
        ax.set_ylabel('Laser power [TW]')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_title('Laser power and intensity vs time')

        ax = axes[1]
        with _np.errstate(invalid='ignore', divide='ignore'):
            if I_grid_outer is not None:
                ax.semilogy(t_ns, I_grid_outer, color='steelblue', lw=1.0,
                            label='I at grid outer')
            if I_at_crit is not None:
                ax.semilogy(t_ns, I_at_crit, color='crimson', lw=1.2,
                            label='I at r_crit')
            if I_peak_coronal_vs_t is not None:
                ax.semilogy(t_ns, I_peak_coronal_vs_t,
                            color='darkorange', lw=0.8, ls='--',
                            label='I peak coronal')
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel(r'Intensity [W/cm$^2$]')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3, which='both')
        if _np.isfinite(t_peak_power):
            for a in axes:
                a.axvline(t_peak_power, color='gray', ls=':', lw=0.8, alpha=0.6)

        pdf.savefig(fig); _plt.close(fig)

        # ---- PAGE 3: Method 1 vs Method 2 cross-check at peak power ---------
        fig, axes = _plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)

        # Pick the peak-power timestep
        t_idx = 0
        if _np.isfinite(t_peak_power) and t_ns is not None:
            t_idx = int(_np.argmin(_np.abs(t_ns - t_peak_power)))
        elif P_laser is not None:
            t_idx = int(_np.argmax(P_laser))

        r_row = r_um[t_idx]
        with _np.errstate(invalid='ignore', divide='ignore'):
            I1_row = I1[t_idx]
            I2_row = I2[t_idx]
            alpha_row = alpha[t_idx]

        ax = axes[0]
        ax.semilogy(r_row, I2_row, color='steelblue', lw=1.3,
                    label='Method 2 (Beer-Lambert)')
        m1_finite = _np.isfinite(I1_row) & (I1_row > 0)
        if m1_finite.any():
            ax.semilogy(r_row[m1_finite], I1_row[m1_finite],
                        color='crimson', lw=0.0, marker='.', ms=3,
                        alpha=0.6, label='Method 1 (P_src / alpha)')
        if _np.isfinite(r_crit[t_idx]):
            ax.axvline(r_crit[t_idx] * 1e4, color='black', ls='--', lw=0.8,
                       label=f'r_crit = {r_crit[t_idx]*1e4:.1f} um')
        ax.set_ylabel(r'Intensity [W/cm$^2$]')
        ax.set_title(f'Radial profile at t = {t_ns[t_idx]:.2f} ns '
                     f'(peak laser power)')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3, which='both')

        ax = axes[1]
        ax.semilogy(r_row, _np.maximum(alpha_row, 1e-4),
                    color='black', lw=1.0)
        ax.set_xlabel(r'Radius [um]')
        ax.set_ylabel(r'Attenuation coefficient $\\alpha$ [1/cm]')
        ax.grid(True, alpha=0.3, which='both')
        if _np.isfinite(r_crit[t_idx]):
            ax.axvline(r_crit[t_idx] * 1e4, color='black', ls='--', lw=0.8)

        pdf.savefig(fig); _plt.close(fig)

    def _plot_laser_power(self, pdf):'''

src, _ = apply(src, ANCHOR, NEW_METHOD,
               sentinel='def _plot_laser_intensity(self, pdf):',
               label='insert _plot_laser_intensity method')


# ====================================================================
# EDIT 2: register the call in create_full_report
# Anchor on the existing call to _plot_laser_deposition, since that's the
# natural predecessor. We insert the intensity call right after.
# ====================================================================
print('\n=== Edit 2: register in create_full_report ===\n')

# The call pattern in create_full_report: "self._plot_laser_deposition(pdf)"
# We need to see what it looks like locally. Most likely form:
OLD = 'self._plot_laser_deposition(pdf)'
NEW = ('self._plot_laser_deposition(pdf)\n'
       '        self._plot_laser_intensity(pdf)')

src, _ = apply(src, OLD, NEW,
               sentinel='self._plot_laser_intensity(pdf)',
               label='create_full_report: register _plot_laser_intensity call')


if src != orig:
    TARGET.write_text(src)
    print(f'\n  [WROTE]   {TARGET.name}')


print('\n-----------------------------------------------------------------------')
print('Verify after running:')
print('  grep -n "_plot_laser_intensity\\|_plot_laser_deposition\\|_plot_laser_power" \\')
print('    helios_postprocess/icf_plotting.py | head -10')
print('')
print('Expected: def _plot_laser_deposition, def _plot_laser_intensity,')
print('          def _plot_laser_power (in that order), plus the call in')
print('          create_full_report.')
print('-----------------------------------------------------------------------')
