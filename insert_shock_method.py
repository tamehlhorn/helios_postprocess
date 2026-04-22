"""Insert _compute_shock_breakout method into icf_analysis.py before _track_ablation_front."""
import os
import sys

# Run from repo root
path = 'helios_postprocess/icf_analysis.py'
if not os.path.exists(path):
    print(f"ERROR: {path} not found. Run from repo root.")
    sys.exit(1)

with open(path) as f:
    t = f.read()

# Check not already present
if '_compute_shock_breakout' in t and 'def _compute_shock_breakout' in t:
    print("Method already present -- nothing to do.")
    sys.exit(0)

# Find insertion point
marker = '    def _track_ablation_front(self):'
if marker not in t:
    print(f"ERROR: insertion marker not found: {marker!r}")
    sys.exit(1)

new_method = '''    def _compute_shock_breakout(self):
        """
        Compute first shock breakout properties in DT ice layer.
        Stores foot-pulse peak pressure, breakout time/pressure/Mach on ICFRunData.
        Skipped silently if no distinct DT ice layer or data missing.
        """
        try:
            from .pressure_gradients import analyze_first_shock
            ri = self.data.region_interfaces_indices
            if ri is None or ri.shape[1] < 2:
                return
            ice_inner_node = int(ri[0, 0])
            ice_outer_node = int(ri[0, 1])
            if ice_outer_node - ice_inner_node < 5:
                logger.info("Skipping shock breakout: no distinct DT ice layer")
                return
            if self.data.ion_pressure is None or self.data.mass_density is None:
                return
            total_pressure = self.data.ion_pressure + self.data.rad_pressure
            t_ns = self.data.time * 1e9
            result = analyze_first_shock(
                pressure=total_pressure,
                zone_boundaries=self.data.zone_boundaries,
                density=self.data.mass_density,
                time=t_ns,
                ablator_outer_zone=ice_outer_node - 1,
                fuel_inner_zone=ice_inner_node,
                search_inner_zone=ice_inner_node,
                dP_dr_threshold=1e9,
            )
            bt = result["breakout_time_ns"]
            self.data.shock_breakout_time_ns = float(bt) if not np.isnan(bt) else 0.0
            bp = result["breakout_pressure_Gbar"]
            self.data.shock_breakout_pressure_Gbar = float(bp) if not np.isnan(bp) else 0.0
            bm = result["breakout_mach"]
            self.data.shock_breakout_mach = float(bm) if not np.isnan(bm) else 0.0
            # Foot-pulse peak: max shock pressure from 4 ns to foot_end+1 ns
            t_foot_end = getattr(self.data, "laser_foot_end_ns", None) or 6.0
            foot_mask = (t_ns >= 4.0) & (t_ns <= float(t_foot_end) + 1.0)
            p_foot = result["shock_pressure_Gbar"][foot_mask]
            valid = p_foot[~np.isnan(p_foot)]
            self.data.shock_foot_pressure_Gbar = float(np.max(valid)) if len(valid) > 0 else 0.0
            logger.info(
                f"Shock foot P={self.data.shock_foot_pressure_Gbar:.4f} Gbar, "
                f"breakout t={self.data.shock_breakout_time_ns:.3f} ns "
                f"P={self.data.shock_breakout_pressure_Gbar:.3f} Gbar "
                f"M={self.data.shock_breakout_mach:.2f}"
            )
        except Exception as e:
            logger.warning(f"Could not compute shock breakout: {e}")

'''

result = t.replace(marker, new_method + marker)
assert 'def _compute_shock_breakout' in result, "Insertion failed"
print(f"Method inserted. File: {len(result.splitlines())} lines (was {len(t.splitlines())})")

with open(path, 'w') as f:
    f.write(result)

print("Done.")
