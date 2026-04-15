"""Test analyze_first_shock on DT ice layer of PDD_9."""
import numpy as np
from helios_postprocess.core import HeliosRun
from helios_postprocess.data_builder import build_run_data
from helios_postprocess.pressure_gradients import analyze_first_shock

base = "/Users/mehlhorn/Sims/Helios_Sims/Xcimer_Sims/Olson_PDD/Olson_PDD_9/Olson_PDD_9"
run = HeliosRun(base + ".exo")
data = build_run_data(run)

total_pressure = data.ion_pressure + data.rad_pressure
ri = data.region_interfaces_indices

# PDD_9 region nodes: [151, 191, 321, 351]
# DT ice = zones 151-190 (sets cold fuel adiabat)
ice_inner_node = int(ri[0, 0])      # 151: ice/hot-spot boundary = breakout
ice_outer_node = int(ri[0, 1])      # 191: foam/ice boundary = shock entry
ice_inner_zone = ice_inner_node     # first zone of DT ice
ice_outer_zone = ice_outer_node - 1 # last zone of DT ice = 190

print(f"DT ice zones: {ice_inner_zone} to {ice_outer_zone}")
print(f"Shock entry at node {ice_outer_node}, breakout at node {ice_inner_node}")

result = analyze_first_shock(
    pressure=total_pressure,
    zone_boundaries=data.zone_boundaries,
    density=data.mass_density,
    time=data.time * 1e9,
    ablator_outer_zone=ice_outer_zone,
    fuel_inner_zone=ice_inner_node,
    search_inner_zone=ice_inner_zone,
)

print(f"\nBreakout time:      {result['breakout_time_ns']:.4f} ns")
print(f"Breakout pressure:  {result['breakout_pressure_Gbar']:.2f} Gbar")
print(f"Breakout Mach:      {result['breakout_mach']:.2f}")

print("\nShock history (first 10 valid timesteps):")
print(f"{'Time(ns)':>10} {'R_shock(um)':>12} {'v_shock(km/s)':>14} {'P_post(Gbar)':>13} {'Mach':>6} {'rho_jump':>9}")
count = 0
for i in range(len(data.time)):
    if not np.isnan(result['shock_radius'][i]):
        r_um = result['shock_radius'][i] * 1e4
        v = result['shock_velocity'][i]
        v_kms = v * 100 if not np.isnan(v) else float('nan')
        print(f"{data.time[i]*1e9:10.4f} {r_um:12.1f} {v_kms:14.1f} "
              f"{result['shock_pressure_Gbar'][i]:13.2f} "
              f"{result['mach_number'][i]:6.2f} "
              f"{result['density_jump'][i]:9.3f}")
        count += 1
        if count >= 10:
            break
