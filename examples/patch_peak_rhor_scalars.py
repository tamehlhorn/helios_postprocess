from pathlib import Path

# ── Edit 1: _compute_burn_propagation in icf_analysis.py ────────────────────
ANALYSIS = Path("helios_postprocess/icf_analysis.py")

A_OLD = """        self.data.hot_spot_rhoR_vs_time = hs_rhoR
        self.data.total_rhoR_vs_time = total_rhoR

        # ---- Ignition time: first crossing of ρR_hs = 0.3 g/cm² ----"""

A_NEW = """        self.data.hot_spot_rhoR_vs_time = hs_rhoR
        self.data.total_rhoR_vs_time = total_rhoR

        # ---- Peak HS ρR (T_ion > 4.5 keV) and peak total ρR with timing ----
        # Olson 2021 "complete propagation" anchor: peak of HS ρR(t); also
        # peak total ρR(t). Both used directly in cross-code comparison.
        if hs_rhoR.size > 0 and np.any(hs_rhoR > 0):
            pk_hs = int(np.argmax(hs_rhoR))
            self.data.peak_hs_rhoR_T_mask = float(hs_rhoR[pk_hs])
            self.data.peak_hs_rhoR_time_ns = float(self.data.time[pk_hs])
            logger.info(
                f"Peak HS ρR (T_ion>4.5 keV): "
                f"{self.data.peak_hs_rhoR_T_mask:.3f} g/cm² "
                f"at t={self.data.peak_hs_rhoR_time_ns:.3f} ns"
            )
        if total_rhoR.size > 0 and np.any(total_rhoR > 0):
            pk_tot = int(np.argmax(total_rhoR))
            self.data.peak_total_rhoR = float(total_rhoR[pk_tot])
            self.data.peak_total_rhoR_time_ns = float(self.data.time[pk_tot])
            logger.info(
                f"Peak total ρR: "
                f"{self.data.peak_total_rhoR:.3f} g/cm² "
                f"at t={self.data.peak_total_rhoR_time_ns:.3f} ns"
            )

        # ---- Ignition time: first crossing of ρR_hs = 0.3 g/cm² ----"""

text = ANALYSIS.read_text()
if A_OLD not in text:
    raise SystemExit(f"OLD block not found in {ANALYSIS}")
ANALYSIS.write_text(text.replace(A_OLD, A_NEW))
print(f"Patched {ANALYSIS}")

# ── Edit 2: BURN PROPAGATION summary block in icf_output.py ─────────────────
OUTPUT = Path("helios_postprocess/icf_output.py")

B_OLD = """        _a(self._metric('HS pressure at ignition',      d.ignition_hs_pressure,   'Gbar', fmt='.1f'))
        _a(self._metric('Complete propagation time',     d.burn_propagation_time,  'ns',   fmt='.3f'))
        _a('')"""

B_NEW = """        _a(self._metric('HS pressure at ignition',      d.ignition_hs_pressure,   'Gbar', fmt='.1f'))
        _a(self._metric('Complete propagation time',     d.burn_propagation_time,  'ns',   fmt='.3f'))
        # Peak ρR scalars (Olson 2021 ρR-vs-time anchors)
        _a(self._metric('Peak total ρR',                 getattr(d, 'peak_total_rhoR', 0.0),         'g/cm²', fmt='.3f'))
        _a(self._metric('  ... at time',                 getattr(d, 'peak_total_rhoR_time_ns', 0.0), 'ns',    fmt='.3f'))
        _a(self._metric('Peak HS ρR (T>4.5 keV)',        getattr(d, 'peak_hs_rhoR_T_mask', 0.0),     'g/cm²', fmt='.3f'))
        _a(self._metric('  ... at time',                 getattr(d, 'peak_hs_rhoR_time_ns', 0.0),    'ns',    fmt='.3f'))
        _a('')"""

text = OUTPUT.read_text()
if B_OLD not in text:
    raise SystemExit(f"OLD block not found in {OUTPUT}")
OUTPUT.write_text(text.replace(B_OLD, B_NEW))
print(f"Patched {OUTPUT}")