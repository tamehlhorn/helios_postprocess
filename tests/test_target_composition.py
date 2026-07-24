"""Tests for target_composition (RHW parse + per-species areal densities) and
the multi-material carbon scatter path."""
import numpy as np
import pytest
from types import SimpleNamespace

from helios_postprocess import target_composition as tc

nesst = pytest.importorskip("NeSST")

_RHW = """[Spatial Grid Data]:
  Parameters for Region =  DT gas
   Min. radius            = 0
   Max. radius            = 0.1875
   Density                = 0.0006
   Mean atomic weight     = 2.51
   EOS filepath           = DT.prp
   Fus reac mass frac D   = 0.5
   Fus reac mass frac T   = 0.5
   Fus reac mass frac H   = 0
  Parameters for Region =  DT/CD foam
   Min. radius            = 0.1875
   Max. radius            = 0.1965
   Density                = 0.3
   Mean atomic weight     = 3.65
   EOS filepath           = DT_0p25_and_CD_0p05.prp
   Fus reac mass frac D   = 0.417
   Fus reac mass frac T   = 0.417
   Fus reac mass frac H   = 0
"""


@pytest.fixture
def rhw(tmp_path):
    p = tmp_path / "vulcan.rhw"
    p.write_text(_RHW)
    return str(p)


def test_parse_regions(rhw):
    regs = tc.parse_rhw_regions(rhw)
    assert len(regs) == 2
    assert regs[0].w_D == 0.5 and regs[0].w_T == 0.5 and regs[0].w_C == 0.0
    # foam carbon = 1 - 0.417 - 0.417
    assert regs[1].w_C == pytest.approx(0.166, abs=1e-3)
    assert regs[1].eos_file == "DT_0p25_and_CD_0p05.prp"


def test_fuel_atom_fractions_from_mass(rhw):
    # equal-mass D:T -> 60:40 by number
    fr = tc.fuel_atom_fractions(tc.parse_rhw_regions(rhw))
    assert fr["frac_D"] == pytest.approx(0.6, abs=0.02)


def _burn_run():
    nt, nz = 30, 60
    r = np.concatenate([np.linspace(0, 0.1875, 46),
                        np.linspace(0.1875, 0.1965, 16)[1:]])
    time = np.linspace(10, 14, nt)
    shape = np.exp(-0.5 * ((time - 12) / 0.2) ** 2)
    zb = np.array([r * (1 - 0.9 * shape[k]) for k in range(nt)])
    rho = np.zeros((nt, nz)); rho[:, :45] = 20 * shape[:, None] + 0.1
    rho[:, 45:] = 5 * shape[:, None] + 0.3
    Ti = np.full((nt, nz), 300.0); Ti[:, :45] = 6000 * shape[:, None] + 300
    dr = np.diff(zb, axis=1)
    zm = rho * (4 * np.pi * (zb[:, :-1] + dr / 2) ** 2 * dr)
    fp = np.zeros((nt, nz)); fp[:, :45] = shape[:, None] * 1e26
    return SimpleNamespace(
        time=time, zone_boundaries=zb, mass_density=rho, ion_temperature=Ti,
        elec_temperature=Ti.copy(), zone_mass=zm, fusion_power=fp,
        fusion_rate_DD_nHe3=fp * 0.01,
        region_interfaces_indices=np.tile(np.array([45, 60, 61]), (nt, 1)),
        dt_neutron_count=np.cumsum(shape) * 1e17,
        areal_density_vs_time=np.full(nt, 1.0))


def test_species_areal_densities(rhw):
    regs = tc.parse_rhw_regions(rhw)
    sp = tc.species_areal_densities(_burn_run(), regs, rhoR_fuel_gcm2=0.9)
    assert sp is not None
    # fuel D+T anchored to the input
    assert sp["rhoR_D_gcm2"] + sp["rhoR_T_gcm2"] == pytest.approx(0.9, rel=1e-6)
    assert sp["rhoR_C_gcm2"] >= 0
    assert set(sp["mass_split"]) == {"D", "T", "C"}


def test_multi_material_report(rhw):
    from helios_postprocess import neutron_scatter as nsc
    m, txt = nsc.neutron_report(_burn_run(), rhw_path=rhw, n_E=700)
    assert m["composition_source"] == "RHW"
    assert "multi-material" in m["dsr_source"]
    assert m["frac_D"] == pytest.approx(0.6, abs=0.05)   # auto-detected
    assert m["DSR_total"] >= m["DSR"]                    # carbon only adds
    assert "DSR total (+carbon)" in txt
    assert "auto-detected" in txt


def test_plot_scatter_tof_handles_both_shapes(tmp_path):
    """Regression: the figure helper must accept the multi-material dict
    (scattered_fuel/scattered_C) as well as the single-material one
    (components) -- previously it KeyError'd on 'components'."""
    from helios_postprocess import neutron_scatter as nsc
    from helios_postprocess import neutron_spectrum as ns
    E, S = ns.synthesize_birth_spectrum(np.array([[1.0]]), np.array([[5000.0]]),
            energy_grid_MeV=np.linspace(1.0, 18.0, 1600), reaction="DT")
    mm = nsc.multi_material_scatter(E, S, 0.36, 0.54, 0.15, n_E=700)   # no "components"
    out = tmp_path / "mm.png"
    nsc.plot_scatter_tof(mm, nsc.scattered_tof(mm), str(out))
    assert out.exists() and out.stat().st_size > 0
    sm = nsc.scatter_from_spectrum(E, S, 0.9, n_E=700, include_n2n=False)  # has "components"
    out2 = tmp_path / "sm.png"
    nsc.plot_scatter_tof(sm, nsc.scattered_tof(sm), str(out2))
    assert out2.exists() and out2.stat().st_size > 0


def test_multi_material_fuel_matches_dt_only(rhw):
    from helios_postprocess import neutron_scatter as nsc
    from helios_postprocess import neutron_spectrum as ns
    E, S = ns.synthesize_birth_spectrum(np.array([[1.0]]), np.array([[5000.0]]),
            energy_grid_MeV=np.linspace(1.0, 18.0, 1600), reaction="DT")
    mm = nsc.multi_material_scatter(E, S, 0.36, 0.54, 0.0, n_E=700)
    dt = nsc.scatter_from_spectrum(E, S, 0.9, frac_D=mm["frac_D"],
                                   frac_T=mm["frac_T"], n_E=700, include_n2n=False)
    # zero carbon -> fuel DSR identical to the D+T-only path
    assert mm["dsr_fuel"]["DSR"] == pytest.approx(dt["dsr"]["DSR"], rel=1e-6)
    assert mm["dsr_total"]["DSR"] == pytest.approx(mm["dsr_fuel"]["DSR"], rel=1e-9)
