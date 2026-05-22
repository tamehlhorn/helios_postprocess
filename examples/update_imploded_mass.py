"""Populate imploded_DT_mass_mg in PDD_20 published JSON files based on Fig 7 integration."""
import json, glob, os

# Mean of LILAC/xRAGE/HYDRA integrated at r <= 160 µm, with uncertainty ~ spread
REF_MASS_MG  = 1.58   # mean(1.62, 1.63, 1.49)
REF_MASS_UNC = 0.10   # ± std-equivalent

ROOTS = [
    "/Users/tommehlhorn/Sims/Xcimer/Olson_PDD",   # Studio mount path
]

updated = 0
for root in ROOTS:
    if not os.path.isdir(root): continue
    for path in glob.glob(f"{root}/**/*_published.json", recursive=True):
        try:
            with open(path) as f: data = json.load(f)
        except Exception as e:
            print(f"SKIP (parse error): {path}: {e}")
            continue
        existing = data.get("imploded_DT_mass_mg", None)
        if existing == [0.0, 0.0] or existing is None:
            data["imploded_DT_mass_mg"] = [REF_MASS_MG, REF_MASS_UNC]
            with open(path, "w") as f:
                json.dump(data, f, indent=4)
            print(f"  UPDATED: {path}")
            updated += 1
        else:
            print(f"  SKIP (already set): {path} = {existing}")

print(f"\n{updated} files updated.")
