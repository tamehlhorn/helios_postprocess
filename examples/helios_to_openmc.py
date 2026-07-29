#!/usr/bin/env python
"""
helios_to_openmc.py  --  Helios -> OpenMC source + geometry builder, the NeSST
cross-check readers, the gamma diagnostic, and (new) a geometry conditioner +
an analytic Method-2 rhoR anchor for triangulating discrepancies.

Two-environment design: `import openmc` is LAZY (inside the functions that need
it), so this module imports cleanly in either your Helios env or the openmc-env.

    dump_snapshot.py   (helios env)  : Helios .exo -> ProfileSnapshot + NeSST ref -> .npz
    compare_openmc_nesst.py (openmc-env): .npz -> OpenMC run -> compare vs NeSST + gammas
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np

_trap = getattr(np, "trapezoid", None) or np.trapz

NEUTRON_EDGES_MeV = np.round(np.arange(0.0, 16.0 + 1e-9, 0.1), 6)
GAMMA_EDGES_MeV = np.round(np.arange(0.0, 12.0 + 1e-9, 0.05), 6)
DSR_SCATTER_MeV = (10.0, 12.0)
DSR_PRIMARY_MeV = (13.0, 15.0)

# n-D elastic at 14.06 MeV ~ 0.61 b (ENDF/B-VIII.0).  read_results prefers the
# value read straight from the OpenMC data (sigma_el_D_barn_at); this is fallback.
SIGMA_EL_D_BARN = 0.61
M_D = 2.0141
N_A = 6.02214076e23


# ============================================================ data snapshot ==
@dataclass
class ProfileSnapshot:
    r_cm: np.ndarray            # shell OUTER radii, strictly increasing (cm)
    rho_gcc: np.ndarray         # mass density per shell (g/cm^3)
    w_D: np.ndarray
    w_T: np.ndarray
    w_C: np.ndarray
    S: np.ndarray               # neutron production weight per shell
    s_per_shell: bool = False   # True: S is per-shell production; False: per-volume
    label: str = "helios_snapshot"
    birth_energy_MeV: np.ndarray | None = None
    birth_spectrum: np.ndarray | None = None

    def to_npz(self, path, **extra):
        d = dict(r_cm=self.r_cm, rho_gcc=self.rho_gcc, w_D=self.w_D, w_T=self.w_T,
                 w_C=self.w_C, S=self.S, s_per_shell=self.s_per_shell,
                 label=self.label)
        if self.birth_energy_MeV is not None:
            d["birth_energy_MeV"] = self.birth_energy_MeV
            d["birth_spectrum"] = self.birth_spectrum
        d.update(extra)
        np.savez(path, **d)

    @staticmethod
    def from_npz(path):
        z = np.load(path, allow_pickle=True)
        be = z["birth_energy_MeV"] if "birth_energy_MeV" in z.files else None
        bs = z["birth_spectrum"] if "birth_spectrum" in z.files else None
        return ProfileSnapshot(
            r_cm=z["r_cm"], rho_gcc=z["rho_gcc"], w_D=z["w_D"], w_T=z["w_T"],
            w_C=z["w_C"], S=z["S"], s_per_shell=bool(z["s_per_shell"]),
            label=str(z["label"]), birth_energy_MeV=be, birth_spectrum=bs)


# ------------------------------------- shell helpers (edges / per-shell prod) -
def _edges(snap):
    return np.concatenate([[0.0], np.asarray(snap.r_cm, dtype=float)])


def _per_shell_production(snap):
    """Return S as per-shell production (neutrons born in the shell)."""
    S = np.clip(np.asarray(snap.S, dtype=float), 0.0, None)
    if snap.s_per_shell:
        return S
    e = _edges(snap)
    vshell = (4.0 / 3.0) * np.pi * (e[1:] ** 3 - e[:-1] ** 3)
    return S * vshell


# ==================================================== geometry conditioning ==
def resample_min_thickness(snap, dr_min_cm=2.0e-4):
    """Merge adjacent Lagrangian zones so no OpenMC shell is thinner than
    dr_min_cm, preserving areal density (∫rho dr), composition (mass-weighted),
    and total production (summed).

    Why: raw Helios zoning splits the compressed core into sub-micron shells;
    near-coincident spheres make OpenMC's surface tracking fail ("Newton-Raphson
    iterations exceeded" / lost particles) and under-count scattering in exactly
    the region that carries the rhoR.  Merging to ~µm shells fixes tracking
    without changing the physics.
    """
    e = _edges(snap)
    dr = np.diff(e)
    rho = np.asarray(snap.rho_gcc, float)
    col = rho * dr                                   # g/cm^2 per zone
    prod = _per_shell_production(snap)
    wD, wT, wC = (np.asarray(snap.w_D, float), np.asarray(snap.w_T, float),
                  np.asarray(snap.w_C, float))

    r_out_new, rho_new, wD_n, wT_n, wC_n, S_n = [], [], [], [], [], []
    i, n = 0, len(dr)
    lo_edge = e[0]
    while i < n:
        j = i
        acc_dr = 0.0
        while j < n and acc_dr < dr_min_cm:
            acc_dr += dr[j]
            j += 1
        sl = slice(i, j)
        thick = e[j] - lo_edge
        c = col[sl].sum()
        p = prod[sl].sum()
        r_out_new.append(e[j])
        rho_new.append(c / thick if thick > 0 else 0.0)   # preserves ∫rho dr
        cw = col[sl].sum()
        if cw > 0:
            wD_n.append(float((col[sl] * wD[sl]).sum() / cw))
            wT_n.append(float((col[sl] * wT[sl]).sum() / cw))
            wC_n.append(float((col[sl] * wC[sl]).sum() / cw))
        else:
            wD_n.append(0.0); wT_n.append(0.0); wC_n.append(1.0)
        S_n.append(p)
        lo_edge = e[j]
        i = j

    out = ProfileSnapshot(
        r_cm=np.array(r_out_new), rho_gcc=np.array(rho_new),
        w_D=np.array(wD_n), w_T=np.array(wT_n), w_C=np.array(wC_n),
        S=np.array(S_n), s_per_shell=True, label=snap.label,
        birth_energy_MeV=snap.birth_energy_MeV, birth_spectrum=snap.birth_spectrum)
    return out


# ================================================ analytic Method-2 anchor ==
def analytic_rhoR_species(snap, species="D"):
    """Emission-weighted, outward-radial areal density of a species, computed
    straight from the snapshot (no transport).  This is the radial form of
    Kyle's Method 2 -- a transport-free anchor: if this matches NeSST but OpenMC
    does not, the discrepancy is in the OpenMC geometry/source, not the profile.
    """
    e = _edges(snap)
    dr = np.diff(e)
    w = {"D": snap.w_D, "T": snap.w_T, "C": snap.w_C}[species]
    col = np.asarray(snap.rho_gcc, float) * np.asarray(w, float) * dr   # g/cm^2 per shell
    cum_out = np.cumsum(col[::-1])[::-1]             # rhoR from each shell outward
    Sp = _per_shell_production(snap)
    if Sp.sum() <= 0:
        Sp = np.ones_like(Sp)
    return float(np.sum(Sp * cum_out) / np.sum(Sp))


# ============================================== OpenMC build (openmc-env) ==
def _dt_muir():
    import openmc
    st = openmc.stats
    for make in (lambda: st.muir(e0=14.06e6, m_rat=5.0, kt=20000.0),
                 lambda: st.Muir(e0=14.06e6, m_rat=5.0, kt=20000.0),
                 lambda: st.Normal(14.06e6, 0.1e6)):
        try:
            return make()
        except Exception:
            continue
    return st.Normal(14.06e6, 0.1e6)


def build_target(snap):
    import openmc
    mats, cells, surfs = [], [], []
    prev = None
    n = len(snap.r_cm)
    for i in range(n):
        m = openmc.Material(name=f"{snap.label}_shell{i}")
        wD, wT, wC = float(snap.w_D[i]), float(snap.w_T[i]), float(snap.w_C[i])
        if wD > 0.0:
            m.add_nuclide("H2", wD, "wo")
        if wT > 0.0:
            m.add_nuclide("H3", wT, "wo")
        if wC > 0.0:
            m.add_element("C", wC, "wo")
        if wD <= 0 and wT <= 0 and wC <= 0:
            m.add_element("C", 1.0, "wo")
        m.set_density("g/cm3", max(float(snap.rho_gcc[i]), 1e-6))
        mats.append(m)
        bt = "vacuum" if i == n - 1 else "transmission"
        s = openmc.Sphere(r=float(snap.r_cm[i]), boundary_type=bt)
        surfs.append(s)
        region = -s if prev is None else (+prev & -s)
        cells.append(openmc.Cell(name=f"shell{i}", fill=m, region=region))
        prev = s
    return openmc.Materials(mats), openmc.Geometry(cells), cells, surfs


def build_source(snap, energy="birth"):
    import openmc
    edges = _edges(snap)
    prod = _per_shell_production(snap)
    if prod.sum() <= 0.0:
        prod = np.diff(edges)
    r_dist = openmc.stats.Tabular(edges, prod, interpolation="histogram")
    space = openmc.stats.SphericalIndependent(
        r=r_dist, cos_theta=openmc.stats.Uniform(-1.0, 1.0),
        phi=openmc.stats.Uniform(0.0, 2.0 * np.pi), origin=(0.0, 0.0, 0.0))
    if energy == "birth" and snap.birth_energy_MeV is not None:
        E = np.asarray(snap.birth_energy_MeV, float) * 1e6
        p = np.clip(np.asarray(snap.birth_spectrum, float), 0.0, None)
        e_dist = openmc.stats.Tabular(E, p, interpolation="linear-linear")
    elif isinstance(energy, tuple):
        E, p = energy
        e_dist = openmc.stats.Tabular(np.asarray(E, float), np.asarray(p, float),
                                      interpolation="linear-linear")
    else:
        e_dist = _dt_muir()
    try:
        return openmc.IndependentSource(space=space, angle=openmc.stats.Isotropic(),
                                        energy=e_dist, particle="neutron")
    except AttributeError:
        return openmc.Source(space=space, angle=openmc.stats.Isotropic(),
                             energy=e_dist, particle="neutron")


def build_tallies(surfs):
    import openmc
    outer = surfs[-1]
    tallies = openmc.Tallies()
    t_del = openmc.Tally(name="D_elastic")
    t_del.nuclides = ["H2"]; t_del.scores = ["elastic"]
    tallies.append(t_del)
    t_n = openmc.Tally(name="neutron_leakage")
    t_n.filters = [openmc.SurfaceFilter(outer.id), openmc.ParticleFilter(["neutron"]),
                   openmc.EnergyFilter(NEUTRON_EDGES_MeV * 1e6)]
    t_n.scores = ["current"]; tallies.append(t_n)
    t_g = openmc.Tally(name="gamma_leakage")
    t_g.filters = [openmc.SurfaceFilter(outer.id), openmc.ParticleFilter(["photon"]),
                   openmc.EnergyFilter(GAMMA_EDGES_MeV * 1e6)]
    t_g.scores = ["current"]; tallies.append(t_g)
    return tallies


def build_model(snap, energy="birth", particles=200_000, batches=20,
                photon_transport=True):
    import openmc
    mats, geom, cells, surfs = build_target(snap)
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = int(particles)
    settings.batches = int(batches)
    settings.photon_transport = bool(photon_transport)
    settings.source = build_source(snap, energy)
    settings.output = {"tallies": False}
    return openmc.Model(geom, mats, settings, build_tallies(surfs))


# =================================================== readers / pure cores ==
def sigma_el_D_barn_at(E_eV=14.06e6):
    """n-D elastic cross section (barn) at E_eV, read from the configured OpenMC
    data library.  Falls back to SIGMA_EL_D_BARN if anything is unavailable."""
    try:
        import openmc.data
        lib = openmc.data.DataLibrary.from_xml()
        rec = lib.get_by_material("H2", data_type="neutron")
        h2 = openmc.data.IncidentNeutron.from_hdf5(rec["path"])
        temp = sorted(h2[2].xs.keys())[0]            # e.g. '294K'
        return float(h2[2].xs[temp](E_eV))
    except Exception:
        return SIGMA_EL_D_BARN


def _counts_window(edges_MeV, counts, lo, hi):
    ctr = 0.5 * (edges_MeV[:-1] + edges_MeV[1:])
    m = (ctr >= lo) & (ctr <= hi)
    return float(np.asarray(counts)[m].sum())


def dsr_from_counts(edges_MeV, counts, scatter=DSR_SCATTER_MeV, primary=DSR_PRIMARY_MeV):
    ns_ = _counts_window(edges_MeV, counts, *scatter)
    np_ = _counts_window(edges_MeV, counts, *primary)
    return (ns_ / np_) if np_ > 0 else float("nan")


def _dNdE(edges_MeV, counts):
    ctr = 0.5 * (edges_MeV[:-1] + edges_MeV[1:])
    return ctr, np.asarray(counts, dtype=float) / np.diff(edges_MeV)


def read_results(sp_path, sigma_el_D_barn=None):
    import openmc
    if sigma_el_D_barn is None:
        sigma_el_D_barn = sigma_el_D_barn_at()
    with openmc.StatePoint(sp_path) as sp:
        C_D = float(sp.get_tally(name="D_elastic").mean.ravel()[0])
        nleak = sp.get_tally(name="neutron_leakage").mean.ravel()
        gleak = sp.get_tally(name="gamma_leakage").mean.ravel()
    rhoR_D = C_D * M_D / (N_A * sigma_el_D_barn * 1e-24)
    n_ctr, n_dNdE = _dNdE(NEUTRON_EDGES_MeV, nleak)
    g_ctr, g_dNdE = _dNdE(GAMMA_EDGES_MeV, gleak)
    return {
        "C_D": C_D, "sigma_el_D_barn": sigma_el_D_barn, "rhoR_D_gcm2": rhoR_D,
        "DSR": dsr_from_counts(NEUTRON_EDGES_MeV, nleak),
        "neutron_MeV": n_ctr, "neutron_dNdE": n_dNdE,
        "neutron_leak_per_source": float(np.sum(nleak)),
        "gamma_MeV": g_ctr, "gamma_dNdE": g_dNdE,
        "gamma_per_source_neutron": float(np.sum(gleak)),
        "gamma_line_4p4_MeV": _counts_window(GAMMA_EDGES_MeV, gleak, 4.2, 4.6),
        "gamma_line_2p2_MeV": _counts_window(GAMMA_EDGES_MeV, gleak, 2.0, 2.4),
    }


# ================================================= Helios parse (helios env) ==
def _find_exo(path) -> str:
    p = Path(path).expanduser()
    if p.suffix == ".exo" and p.is_file():
        return str(p)
    cand = p.with_suffix(".exo")
    if cand.is_file():
        return str(cand)
    if p.is_dir():
        named = p / f"{p.name}.exo"
        if named.is_file():
            return str(named)
        exos = sorted(p.glob("*.exo")) or sorted(p.rglob("*.exo"))
        if exos:
            return str(exos[0])
    if p.parent.is_dir():
        exos = sorted(p.parent.glob("*.exo"))
        if exos:
            return str(exos[0])
    raise FileNotFoundError(f"No .exo found for {path!r}")


def _find_rhw(run, exo):
    cand = Path(exo).with_suffix(".rhw")
    if cand.exists():
        return str(cand)
    hits = sorted(Path(run).glob("*.rhw")) if Path(run).is_dir() else []
    return str(hits[0]) if hits else None


def from_helios_snapshot(run, rhw_path=None, time_index=None):
    from helios_postprocess.core import HeliosRun
    from helios_postprocess.data_builder import build_run_data
    from helios_postprocess import neutron_spectrum as ns
    from helios_postprocess import target_composition as tc

    exo = _find_exo(run)
    data = build_run_data(HeliosRun(exo))
    nd = ns.extract_neutronics(data=data, use_rhino=False)
    if nd is None:
        raise SystemExit("extract_neutronics returned None (no-burn / missing fields).")

    rhw = rhw_path or _find_rhw(run, exo)
    regions = tc.parse_rhw_regions(rhw) if rhw else []

    zb = np.asarray(data.zone_boundaries, dtype=float)
    rho = np.asarray(data.mass_density, dtype=float)
    t = tc._bang_index(data) if time_index is None else int(time_index)
    t = max(0, min(t, rho.shape[0] - 1))

    r_out = zb[t][1:].astype(float)
    rho_t = rho[t].astype(float)

    if regions:
        zidx = tc.assign_zones_to_regions(zb, regions)
        wD, wT, wC = tc.zone_mass_fractions(zidx, regions)
    else:
        wD = np.full_like(rho_t, 0.5); wT = np.full_like(rho_t, 0.5)
        wC = np.zeros_like(rho_t)

    fp = getattr(data, "fusion_power", None)
    zm = getattr(data, "zone_mass", None)
    if fp is not None:
        fp_t = np.asarray(fp, dtype=float)
        fp_t = fp_t[t] if fp_t.ndim == 2 else fp_t
        if zm is not None:
            zm_a = np.asarray(zm, dtype=float)
            zm_t = zm_a[t] if zm_a.ndim == 2 else zm_a
            S = fp_t * zm_t
        else:
            S = fp_t
    else:
        S = np.ones_like(rho_t)
    S = np.clip(np.asarray(S, dtype=float), 0.0, None)

    keep = (rho_t > 0.0) & (r_out > 0.0)
    r_out, rho_t = r_out[keep], rho_t[keep]
    wD, wT, wC, S = wD[keep], wT[keep], wC[keep], S[keep]
    incr = np.concatenate([[True], np.diff(r_out) > 0.0])
    r_out, rho_t = r_out[incr], rho_t[incr]
    wD, wT, wC, S = wD[incr], wT[incr], wC[incr], S[incr]

    be = bs = None
    ql = getattr(nd, "quicklook", None) or {}
    if ql.get("energy_MeV") is not None and ql.get("birth_spectrum") is not None:
        be = np.asarray(ql["energy_MeV"], dtype=float)
        bs = np.asarray(ql["birth_spectrum"], dtype=float)

    snap = ProfileSnapshot(r_out, rho_t, wD, wT, wC, S, s_per_shell=True,
                           label=Path(exo).stem, birth_energy_MeV=be,
                           birth_spectrum=bs)
    return snap, nd, data, rhw
