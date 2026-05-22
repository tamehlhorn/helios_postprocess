"""Populate Olson-2021-derived reference values in PDD_20 published JSON files.

Two keys are filled here:

  imploded_DT_mass_mg          -- integration of Fig 7 ρ(r) profiles
                                  (mean LILAC/xRAGE/HYDRA inside r=160 µm)

  T_ion_onaxis_ignition_keV    -- on-axis T_ion read off Fig 7 at the
                                  ignition (ρR_hs = 0.3 g/cm²) timestep,
                                  averaged across the three codes
"""
import json, glob, os

# Mean of LILAC/xRAGE/HYDRA integrated at r <= 160 µm, with uncertainty ~ spread
REF_MASS_MG  = 1.58   # mean(1.62, 1.63, 1.49)
REF_MASS_UNC = 0.10   # ± std-equivalent

# On-axis T_ion at ignition (Olson 2021 Fig 7 central values)
# LILAC ~15, xRAGE ~12, HYDRA ~14 keV at r=0
REF_T_ON_AXIS_KEV = 14.0
REF_T_ON_AXIS_UNC = 2.0

REFS = {
    "imploded_DT_mass_mg":        [REF_MASS_MG,        REF_MASS_UNC],
    "T_ion_onaxis_ignition_keV":  [REF_T_ON_AXIS_KEV,  REF_T_ON_AXIS_UNC],
}

ROOTS = [
    "/Users/tommehlhorn/Sims/Xcimer/Olson_PDD",   # Studio mount path
]

updated = 0
touched_files = 0
for root in ROOTS:
    if not os.path.isdir(root): continue
    for path in glob.glob(f"{root}/**/*_published.json", recursive=True):
        try:
            with open(path) as f: data = json.load(f)
        except Exception as e:
            print(f"SKIP (parse error): {path}: {e}")
            continue
        file_dirty = False
        for key, val in REFS.items():
            existing = data.get(key, None)
            if existing == [0.0, 0.0] or existing is None:
                data[key] = val
                print(f"  UPDATED {key:<32s} -> {val} in {os.path.basename(path)}")
                updated += 1
                file_dirty = True
            else:
                print(f"  SKIP {key:<32s} (already set in {os.path.basename(path)} = {existing})")
        if file_dirty:
            with open(path, "w") as f:
                json.dump(data, f, indent=4)
            touched_files += 1

print(f"\n{updated} key-value pairs updated across {touched_files} files.")
