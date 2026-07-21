"""Analytic unit tests for helios_postprocess.neutron_spectrum.

Physics validated with constructed inputs -- no .exo needed. Covers the
in-house quick-look primitives (Brysk width, birth spectrum, DSR) and the
Kyle-aligned extraction layer (volumetric rates, spherical geometry,
neutron-weighted profiles, npz schema, RHINO-absent fallback).
"""
import numpy as np
import pytest

from helios_postprocess import neutron_spectrum as ns
from helios_postprocess import neutron_downscatter as nds


# ----------------------- Brysk / birth spectrum -----------------------

def test_brysk_relation_roundtrip():
    for Ti in (2.0, 5.0, 10.0, 20.0):
        fwhm = ns.brysk_fwhm_keV(Ti, "DT")
        assert ns.fwhm_keV_to_Ti_keV(fwhm, "DT") == pytest.approx(Ti, rel=1e-9)
    assert ns.brysk_fwhm_keV(10.0, "DT") == pytest.approx(177.0 * np.sqrt(10.0), rel=1e-9)


def test_dt_birth_energy_matches_helios_label():
    assert ns.E_BIRTH_MEV["DT"] == pytest.approx(14.06)
    assert ns.E_BIRTH_MEV["DD"] == pytest.approx(2.45)


def test_single_temperature_spectrum_infers_that_temperature():
    Ti = 10.0
    E, S = ns.synthesize_birth_spectrum(np.array([[1.0]]), np.array([[Ti * 1000.0]]), reaction="DT")
    assert np.trapezoid(S, E) == pytest.approx(1.0, rel=1e-3)
    out = ns.infer_ntof_tion(E, S, reaction="DT")
    assert out["Ti_ntof_keV"] == pytest.approx(Ti, rel=0.03)
    assert out["peak_energy_MeV"] == pytest.approx(ns.E_BIRTH_MEV["DT"], abs=0.01)


def test_burn_weighted_temperature_two_zone():
    weights = np.array([[3.0, 1.0]]); Ti_eV = np.array([[5000.0, 9000.0]])
    hist = ns.burn_history(np.array([0.0]), weights, Ti_eV)
    assert hist["Ti_burn_avg_keV"] == pytest.approx(6.0, rel=1e-9)
    E, S = ns.synthesize_birth_spectrum(weights, Ti_eV, reaction="DT")
    assert ns.infer_ntof_tion(E, S, reaction="DT")["Ti_ntof_keV"] == pytest.approx(6.0, rel=0.04)


def test_dsr_roundtrip_forward_model():
    res = ns.synthetic_dsr(1.02, "NIF")
    assert res["DSR"] == pytest.approx(1.02 / 20.4, rel=1e-9)
    assert res["rhoR_roundtrip_gcm2"] == pytest.approx(1.02, rel=1e-6)


def test_downscattered_spectrum_recovers_dsr():
    E, primary = ns.synthesize_birth_spectrum(
        np.array([[1.0]]), np.array([[4000.0]]),
        energy_grid_MeV=np.linspace(9.0, 16.0, 3000), reaction="DT")
    full = ns.build_downscattered_spectrum(E, primary, 0.06)
    assert nds.calculate_downscatter_ratio(full, E)["DSR"] == pytest.approx(0.06, rel=0.05)


def test_burn_history_bang_and_fwhm():
    t = np.linspace(10.0, 14.0, 401); t0, sig = 12.0, 0.10
    weights = np.exp(-0.5 * ((t - t0) / sig) ** 2)[:, None]
    hist = ns.burn_history(t, weights, np.full((t.size, 1), 5000.0))
    assert hist["bang_time_ns"] == pytest.approx(t0, abs=0.02)
    assert hist["burn_fwhm_ns"] == pytest.approx(2.3548 * sig, rel=0.02)


# ----------------------- Kyle-aligned extraction -----------------------

def test_volumetric_rate_convention():
    rate_pg = np.array([[2.0, 4.0]]); rho = np.array([[3.0, 0.5]])
    assert np.allclose(ns.volumetric_rate(rate_pg, rho), [[6.0, 2.0]])


def test_spherical_com_and_volume():
    zb = np.array([[0.0, 1.0, 2.0]])
    com, vol = ns.spherical_com_and_volume(zb)
    # outer shell volume 4/3 pi (2^3 - 1^3) = 4/3 pi 7
    assert vol[0, 1] == pytest.approx(4.0 / 3.0 * np.pi * 7.0, rel=1e-9)
    # COM radius lies between the boundaries and is volume-weighted outward
    assert 1.0 < com[0, 1] < 2.0


def test_neutron_weighted_profile_when_weighting():
    # 3 zones, fixed geometry; rho uniform per timestep (1 then 3); DT production
    # temporal weights 3:9 -> <rho> = (3*1 + 9*3)/12 = 2.5 across the covered grid.
    r_com = np.array([[0.05, 0.15, 0.25], [0.05, 0.15, 0.25]])
    cell_vol = np.ones((2, 3))
    dt_s = np.array([1.0, 1.0])
    rho = np.array([[1.0, 1.0, 1.0], [3.0, 3.0, 3.0]])
    Tion = np.full((2, 3), 4000.0)
    # DT rate_pg=1 both steps -> vol_rate = rate*rho -> temporal_w = sum(vol*cellvol*dt)
    # t0: 3*(1*1)=3 ; t1: 3*(1*3)=9
    rates = {"DT_nHe4": np.ones((2, 3))}
    prof = ns.neutron_weighted_profiles(r_com, cell_vol, dt_s, rho, Tion, rates, n_r=101)
    dt = prof["avg_results"]["DT_nHe4"]
    r = prof["r_avg_cm"]
    idx = int(np.argmin(np.abs(r - 0.15)))
    assert dt["rho_avg_gcc"][idx] == pytest.approx(2.5, rel=1e-6)
    assert dt["T_ion_avg_eV"][idx] == pytest.approx(4000.0, rel=1e-6)


class _BurnRun:
    """Minimal ICFRunData-like object for a small igniting run."""
    def __init__(self, nt=40, nz=20):
        t = np.linspace(10.0, 14.0, nt)
        zb = np.tile(np.linspace(0.0, 0.03, nz + 1), (nt, 1))     # fixed geometry
        _, cell_vol = ns.spherical_com_and_volume(zb)
        rho = np.full((nt, nz), 0.05)
        rho[:, :6] = 5.0                                          # dense hot spot
        Ti = np.full((nt, nz), 400.0)
        Ti[:, :6] = 9000.0                                       # 9 keV hot spot
        shape = np.exp(-0.5 * ((t - 12.0) / 0.15) ** 2)
        fp = np.zeros((nt, nz)); fp[:, :6] = shape[:, None] * 1e26
        self.time = t
        self.zone_boundaries = zb
        self.mass_density = rho
        self.ion_temperature = Ti
        self.elec_temperature = Ti.copy()
        self.zone_mass = rho * cell_vol
        self.fusion_power = fp
        self.fusion_rate_DD_nHe3 = fp * 0.01
        # region interfaces: hot-spot(6), fuel/ablator(12), grid-edge stored as
        # a count (nz+1) -- one past the last valid boundary index, as real
        # Helios runs do; exercises the clip in _interfaces_bang_native.
        ri = np.tile(np.array([6, 12, nz + 1]), (nt, 1))
        self.region_interfaces_indices = ri
        self.dt_neutron_count = np.cumsum(shape) * 1e17
        self.areal_density_vs_time = np.full(nt, 1.0)


def test_extract_neutronics_end_to_end_native():
    nd = ns.extract_neutronics(data=_BurnRun(), use_rhino=True, sim_path=None)
    assert nd is not None
    assert nd.interfaces_source == "native"          # no RHINO available here
    assert nd.avg_results["DT_nHe4"] is not None
    assert nd.quicklook["ntof_tion_keV"] == pytest.approx(9.0, rel=0.06)
    assert np.isfinite(nd.quicklook["rhoR_emission_gcm2"])
    assert nd.quicklook["dsr"] == pytest.approx(
        nd.quicklook["rhoR_emission_gcm2"] / 20.4, rel=1e-6)
    # volumetric-rate channel carried (Kyle schema)
    assert "DT_nHe4" in nd.channel_rates_vol and "TT_nnHe4" in nd.channel_rates_vol


def test_resolve_dt_source():
    t = np.array([0.0, 1.0, 3.0, 6.0])                 # ns, nonuniform
    grad = np.gradient(t * 1e-9)

    class _D:
        timestep_size_s = np.array([0.1, 0.2, 0.3, 0.4]) * 1e-9

    class _N:
        timestep_size_s = None

    assert np.allclose(ns._resolve_dt_s(_D(), t, "gradient"), grad)
    assert np.allclose(ns._resolve_dt_s(_D(), t, "exodus"),
                       np.array([0.1, 0.2, 0.3, 0.4]) * 1e-9)
    assert np.allclose(ns._resolve_dt_s(_N(), t, "exodus"), grad)     # graceful fallback
    with pytest.raises(ValueError):
        ns._resolve_dt_s(_D(), t, "bogus")


def test_extract_exodus_dt_source_runs():
    run = _BurnRun()
    run.timestep_size_s = np.gradient(run.time * 1e-9) * 0.5          # valid shape
    nd = ns.extract_neutronics(data=run, use_rhino=False, dt_source="exodus")
    assert nd is not None and nd.avg_results["DT_nHe4"] is not None


def test_extract_neutronics_npz_roundtrip(tmp_path):
    out = ns.extract_neutronics(data=_BurnRun(), use_rhino=False,
                                save_npz=str(tmp_path))
    assert out is not None
    d = np.load(tmp_path / "neutronics_data.npz", allow_pickle=True)
    for key in ("stime_ns", "r_com_cm", "cell_volume_cm3", "rate_DT_nHe4",
                "r_avg_cm", "rho_avg_gcc", "avg_results", "total_dt_yield",
                "bang_time_ns", "quicklook"):
        assert key in d.files
    assert d["avg_results"].item()["DT_nHe4"] is not None
    assert float(d["total_dt_yield"]) > 0


def test_extract_graceful_on_noburn():
    class _NoBurn:
        time = np.array([0.0, 1.0, 2.0])
        zone_boundaries = np.tile(np.linspace(0, 0.03, 5), (3, 1))
        mass_density = np.ones((3, 4))
        ion_temperature = np.full((3, 4), 100.0)
        zone_mass = np.ones((3, 4))
        fusion_power = np.zeros((3, 4))
    assert ns.extract_neutronics(data=_NoBurn(), use_rhino=False) is None
    assert ns.analyze_neutron_spectrum(_NoBurn()) is None
