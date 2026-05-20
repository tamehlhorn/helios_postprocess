import json
from pathlib import Path

BAM = Path("helios_postprocess/burn_averaged_metrics.py")
text = BAM.read_text()

# ── Edit 1: pull the two attrs off `data` into the histories dict ───────────
H_OLD = """        'P_hs_ignition_Gbar':    getattr(data, 'ignition_hs_pressure', 0.0),
        'hs_radius_ignition_um': getattr(data, 'ignition_hs_radius',   0.0) * 1e4,"""
H_NEW = """        'P_hs_ignition_Gbar':    getattr(data, 'ignition_hs_pressure', 0.0),
        'hs_radius_ignition_um': getattr(data, 'ignition_hs_radius',   0.0) * 1e4,
        'peak_total_rhoR':       getattr(data, 'peak_total_rhoR',       0.0),
        'peak_hs_rhoR_T_mask':   getattr(data, 'peak_hs_rhoR_T_mask',   0.0),"""
if H_OLD not in text:
    raise SystemExit("Edit 1 anchor not found")
text = text.replace(H_OLD, H_NEW)

# ── Edit 2: copy them from histories into sim_metrics ───────────────────────
M_OLD = """        'P_hs_ignition_Gbar':    histories.get('P_hs_ignition_Gbar',    0.0),
        'hs_radius_ignition_um': histories.get('hs_radius_ignition_um', 0.0),"""
M_NEW = """        'P_hs_ignition_Gbar':    histories.get('P_hs_ignition_Gbar',    0.0),
        'hs_radius_ignition_um': histories.get('hs_radius_ignition_um', 0.0),
        'peak_total_rhoR':       histories.get('peak_total_rhoR',       0.0),
        'peak_hs_rhoR_T_mask':   histories.get('peak_hs_rhoR_T_mask',   0.0),"""
if M_OLD not in text:
    raise SystemExit("Edit 2 anchor not found")
text = text.replace(M_OLD, M_NEW)

# ── Edit 3: add two rows to burn_rows so they render in the table ───────────
R_OLD = """        ('CR_max',           sim_metrics['CR_max'],         'CR_max',  '.1f'),
        ('Yield (MJ)',       sim_metrics['yield_MJ'],       'yield',   '.1f'),"""
R_NEW = """        ('CR_max',           sim_metrics['CR_max'],         'CR_max',  '.1f'),
        ('Peak total ρR (g/cm²)',     sim_metrics.get('peak_total_rhoR', 0.0),
         'peak_total_rhoR_gcm2',      '.2f'),
        ('Peak HS ρR T>4.5 (g/cm²)',  sim_metrics.get('peak_hs_rhoR_T_mask', 0.0),
         'peak_hs_rhoR_T_mask_gcm2',  '.2f'),
        ('Yield (MJ)',       sim_metrics['yield_MJ'],       'yield',   '.1f'),"""
if R_OLD not in text:
    raise SystemExit("Edit 3 anchor not found")
text = text.replace(R_OLD, R_NEW)

BAM.write_text(text)
print(f"Patched {BAM}")

# ── Edit 4: add the cluster keys to the published JSON ──────────────────────
# Adjust this path to the run's _published.json (and the 26c master).
PUB = Path("Olson_PDD_26c_burn_published.json")  # repo-local copy if present
if PUB.exists():
    d = json.loads(PUB.read_text())
    d["peak_total_rhoR_gcm2"]     = [1.05, 0.05]
    d["peak_hs_rhoR_T_mask_gcm2"] = [0.85, 0.10]
    PUB.write_text(json.dumps(d, indent=4))
    print(f"Patched {PUB}")
else:
    print(f"(skipped {PUB} — add the two cluster keys by hand to each "
          f"_published.json the runs actually load)")