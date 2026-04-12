"""
Patch 2 — Fix cr_inflight in compute_performance_metrics() in icf_analysis.py

Currently uses hot-spot boundary node as Rf (gives ~6.3-6.8).
Should use ablation front zone radius (gives ~2.16 for PDD_9).
Physics Convention #12.

Usage: python3 ~/patch2_crinflight.py
"""
import os, sys

path = os.path.expanduser(
    "~/helios_postprocessor/helios_postprocess/icf_analysis.py"
)

with open(path) as f:
    src = f.read()

OLD = '''            R0 = self.data.zone_boundaries[0, int(self.data.region_interfaces_indices[0, 0])]
            Rf = self.data.zone_boundaries[pv_idx, int(self.data.region_interfaces_indices[pv_idx, 0])]
            if Rf > 0:
                self.data.cr_inflight = R0 / Rf
                logger.info(f"In-flight CR = R0/Rf: {R0:.4f} cm / {Rf:.4f} cm = {self.data.cr_inflight:.2f}")'''

NEW = '''            R0 = self.data.zone_boundaries[0, int(self.data.region_interfaces_indices[0, 0])]
            # Use ablation front radius at peak velocity (Convention #12).
            # The HS boundary radius gives cr_inflight ~ 6-7 (wrong);
            # the ablation front gives ~ 2.16 for PDD_9 (correct).
            abl_indices = getattr(self.data, 'ablation_front_indices', None)
            if abl_indices is not None and int(abl_indices[pv_idx]) > 0:
                abl_idx = int(abl_indices[pv_idx])
                Rf = 0.5 * (self.data.zone_boundaries[pv_idx, abl_idx] +
                            self.data.zone_boundaries[pv_idx, abl_idx + 1])
            else:
                # Fallback if ablation front not yet computed
                Rf = self.data.zone_boundaries[pv_idx, int(self.data.region_interfaces_indices[pv_idx, 0])]
            if Rf > 0:
                self.data.cr_inflight = R0 / Rf
                logger.info(f"In-flight CR = R0/Rf: {R0:.4f} cm / {Rf:.4f} cm = {self.data.cr_inflight:.2f}")'''

if OLD not in src:
    print("ERROR: search string not found — file may have changed.")
    print("Looking for:\n", OLD[:200])
    sys.exit(1)

new_src = src.replace(OLD, NEW, 1)

with open(path, 'w') as f:
    f.write(new_src)

print("patch2_crinflight.py applied successfully.")
print("  cr_inflight now uses ablation_front_indices[pv_idx] for Rf.")
