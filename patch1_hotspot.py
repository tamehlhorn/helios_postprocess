"""
Patch 1 — Fix _compute_hot_spot_properties() in icf_analysis.py

Replaces temperature-mask approach with region-interface approach.
The temperature mask (T > 1 keV) is unreliable for igniting capsules:
alpha heating warms the dense shell above 1 keV, causing the mask to
extend far outside the true hot spot.  Physics Convention #13.

Usage: python3 ~/patch1_hotspot.py
"""
import os, sys

path = os.path.expanduser(
    "~/helios_postprocessor/helios_postprocess/icf_analysis.py"
)

with open(path) as f:
    src = f.read()

OLD = '''        temp_threshold = self.config.get('hot_spot_temp_threshold', 1000.0)  # eV
        
        try:
            temperatures = self.data.ion_temperature[stag_idx]
            hot_spot_mask = temperatures > temp_threshold
            
            if np.any(hot_spot_mask):
                boundaries = self.data.zone_boundaries[stag_idx]
                radii = (boundaries[:-1] + boundaries[1:]) / 2
                
                # Hot spot radius (outermost hot zone)
                hot_radii = radii[hot_spot_mask]
                self.data.stagnation_hot_spot_radius = np.max(hot_radii)
                
                # Core radius = outer boundary of the outermost hot-spot zone
                hot_zone_indices = np.where(hot_spot_mask)[0]
                self.data.core_radius = boundaries[hot_zone_indices[-1] + 1]
                
                # Hot spot pressure (mass-averaged)'''

NEW = '''        try:
            ri = self.data.region_interfaces_indices
            boundaries = self.data.zone_boundaries[stag_idx]

            # Use region interface to define hot-spot zones (Convention #13).
            # Temperature mask is unreliable for igniting capsules — alpha
            # heating warms the dense shell above 1 keV, giving unphysically
            # large radii (e.g. 0.15 cm > R_hs_true for 26b_burn).
            if ri is not None:
                hs_bnd = int(ri[stag_idx, 0])          # node index of HS boundary
                n_zones = self.data.mass_density.shape[1]
                hot_spot_mask = np.zeros(n_zones, dtype=bool)
                hot_spot_mask[:hs_bnd] = True
                # HS outer radius = zone_boundaries at the boundary node
                self.data.core_radius = float(boundaries[hs_bnd])
                self.data.stagnation_hot_spot_radius = self.data.core_radius
            else:
                # Fallback only when region data absent
                temp_threshold = self.config.get('hot_spot_temp_threshold', 1000.0)
                temperatures = self.data.ion_temperature[stag_idx]
                hot_spot_mask = temperatures > temp_threshold
                if np.any(hot_spot_mask):
                    hot_zone_indices = np.where(hot_spot_mask)[0]
                    self.data.core_radius = float(boundaries[hot_zone_indices[-1] + 1])
                    self.data.stagnation_hot_spot_radius = self.data.core_radius

            if np.any(hot_spot_mask):
                # Hot spot pressure (mass-averaged)'''

if OLD not in src:
    print("ERROR: search string not found — file may have changed.")
    print("Looking for:\n", OLD[:200])
    sys.exit(1)

new_src = src.replace(OLD, NEW, 1)

with open(path, 'w') as f:
    f.write(new_src)

print("patch1_hotspot.py applied successfully.")
print("  _compute_hot_spot_properties() now uses region interface (ri[:,0])")
print("  instead of temperature mask.")
