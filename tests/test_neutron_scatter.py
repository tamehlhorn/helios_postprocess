"""Tests for helios_postprocess.neutron_scatter (NeSST down-scatter adapter).

NeSST-backed tests are guarded with importorskip so the suite still runs where
NeSST is absent; they use a small elastic-only grid (fast, low memory) since the
physics under test -- DSR scaling, rhoR round-trip, calibration -- does not need
the (n,2n) channels or a fine grid. One test exercises the default n2n-on path.
"""
import numpy as np
import pytest

from helios_postprocess import neutron_spectrum as ns
from helios_postprocess import neutron_scatter as nsc

nesst = pytest.importorskip("NeSST")

# Shared small/fast grid config so the NeSST matrix cache is reused across tests.
N_E = 400
KW = dict(n_E=N_E, include_n2n=False)


@pytest.fixture(scope="module")
def birth():
    """A Ti=5 keV DT birth spectrum on a fine grid (unit yield)."""
    E, S = ns.synthesize_birth_spectrum(
        np.array([[1.0]]), np.array([[5000.0]]),
        energy_grid_MeV=np.linspace(1.0, 18.0, 1400), reaction="DT")
    return E, S


def test_dsr_increases_and_is_roughly_linear_in_rhoR(birth):
    E, S = birth
    rhoRs = [0.3, 0.6, 0.9, 1.2]
    # raw single-scatter (no primary depletion) is near-linear in rhoR
    dsrs = [nsc.scatter_from_spectrum(E, S, r, deplete_primary=False, **KW)["dsr"]["DSR"]
            for r in rhoRs]
    assert all(b > a for a, b in zip(dsrs, dsrs[1:]))          # strictly increasing
    slopes = [d / r for d, r in zip(dsrs, rhoRs)]
    assert (max(slopes) - min(slopes)) / np.mean(slopes) < 0.10


def test_rhoR_roundtrips_through_A1s(birth):
    E, S = birth
    for rhoR in (0.3, 0.9, 1.3):
        res = nsc.scatter_from_spectrum(E, S, rhoR, **KW)
        assert res["rhoR_from_A1s_gcm2"] == pytest.approx(rhoR, rel=1e-6)


def test_dsr_in_physical_range_and_matches_nif_ballpark(birth):
    E, S = birth
    res = nsc.scatter_from_spectrum(E, S, 0.9, **KW)
    dsr = res["dsr"]["DSR"]
    # NIF empirical rhoR = 20.4 * DSR -> DSR ~ 0.044 at 0.9 g/cm^2; first-
    # principles single-scatter sits a bit lower (higher coeff). Bracket wide.
    assert 0.02 < dsr < 0.07
    coeff = 0.9 / dsr
    assert 18.0 < coeff < 28.0            # cf. empirical NIF 20.4


def test_windows_are_respected(birth):
    E, S = birth
    res = nsc.scatter_from_spectrum(E, S, 0.9, **KW)
    d = res["dsr"]
    assert d["scatter_window_MeV"] == nsc.DSR_SCATTER_WINDOW_MEV
    assert d["primary_window_MeV"] == nsc.DSR_PRIMARY_WINDOW_MEV
    # numerator uses the scattered array only (no primary leakage): N_scatter
    # equals the scattered integral over the window
    E_out = res["energy_MeV"]; sc = res["scattered"]
    m = (E_out >= 10.0) & (E_out <= 12.0)
    assert d["N_scatter"] == pytest.approx(np.trapezoid(sc[m], E_out[m]), rel=1e-9)


def test_fuel_loading_roundtrips_for_D_enhanced(birth):
    E, S = birth
    res = nsc.scatter_from_spectrum(E, S, 0.9, frac_D=0.7, frac_T=0.3, **KW)
    assert res["rhoR_from_A1s_gcm2"] == pytest.approx(0.9, rel=1e-6)
    assert res["frac_D"] == 0.7 and res["frac_T"] == 0.3


def test_calibration_slope_first_principles(birth):
    cal = nsc.dsr_rhoR_slope(**KW)
    assert cal["slope_DSR_per_gcm2"] > 0
    assert 18.0 < cal["coeff_gcm2_per_DSR"] < 28.0
    # slope and DSR_ref are consistent
    assert cal["DSR_ref"] == pytest.approx(
        cal["slope_DSR_per_gcm2"] * cal["rhoR_ref_gcm2"], rel=1e-9)


def test_scattered_tof_recovers_Ti_and_shows_tail(birth):
    E, S = birth
    res = nsc.scatter_from_spectrum(E, S, 0.9, **KW)
    out = nsc.scattered_tof(res, distance_m=3.0, irf_fwhm_ns=0.5)
    # Ti read from the fine birth spectrum, not the coarse scatter grid
    assert out["primary"]["Ti_ntof_keV"] == pytest.approx(5.0, rel=0.05)
    # DSR carried through
    assert out["dsr"]["DSR"] == pytest.approx(res["dsr"]["DSR"], rel=1e-9)
    # a down-scatter tail exists after the primary, before the 10 MeV edge time
    t = out["time_ns"]; y = out["signal_detector"]
    tp = out["primary"]["t_peak_ns"]
    tail = (t > tp + 3.0) & (t < out["scatter_tail_ns"] + 5.0)
    assert y[tail].max() > 1e-4 * y.max()


def test_full_equals_primary_plus_scattered(birth):
    E, S = birth
    Y_in = np.trapezoid(S, E)
    # undepleted: primary integral == input yield (double-counts scatter)
    res = nsc.scatter_from_spectrum(E, S, 0.9, deplete_primary=False, **KW)
    assert np.allclose(res["full"], res["primary"] + res["scattered"])
    assert np.trapezoid(res["primary"], res["energy_MeV"]) == pytest.approx(Y_in, rel=0.02)
    # depleted, elastic-only: primary attenuated; escaping ~ conserved (~1)
    res2 = nsc.scatter_from_spectrum(E, S, 0.9, **KW)
    assert np.trapezoid(res2["primary"], res2["energy_MeV"]) < Y_in
    assert 0.95 < res2["escaping_fraction"] < 1.06


def test_default_path_includes_n2n(birth):
    """The default (n2n on, n_E=500) runs and gives a comparable DSR."""
    E, S = birth
    res = nsc.scatter_from_spectrum(E, S, 0.9)          # defaults: n2n on
    assert 0.02 < res["dsr"]["DSR"] < 0.07


def test_graceful_error_without_nesst(birth, monkeypatch):
    E, S = birth
    monkeypatch.setattr(nsc, "NESST_AVAILABLE", False)
    monkeypatch.setattr(nsc, "_NESST_IMPORT_ERROR", ImportError("stub"), raising=False)
    with pytest.raises(ImportError, match="NeSST is required"):
        nsc.scatter_from_spectrum(E, S, 0.9, **KW)


def test_rejects_zero_yield_spectrum():
    E = np.linspace(1.0, 18.0, 200)
    with pytest.raises(ValueError, match="non-positive integral"):
        nsc.scattered_spectrum(E, np.zeros_like(E), 0.9, **KW)


# ----------------------- pipeline helper: neutron_report -----------------------

from test_neutron_spectrum import _BurnRun          # noqa: E402


def test_neutron_report_with_nesst():
    m, txt = nsc.neutron_report(_BurnRun(), n_E=N_E, include_n2n=False)
    assert m is not None
    assert "NeSST" in m["dsr_source"]
    assert m["DSR"] > 0 and np.isfinite(m["rhoR_from_DSR_gcm2"])
    assert "NEUTRON DIAGNOSTICS" in txt and "DSR fuel (D+T)" in txt
    assert "round-trip" in txt


def test_neutron_report_published_overlay():
    pub = {"yield_neutrons": [4.86e17, 0], "Tion": [9.0, 0.5],
           "DSR": [2.87, 0.24], "bang_time_ns": [12.0, 0.1]}
    m, txt = nsc.neutron_report(_BurnRun(), n_E=N_E, include_n2n=False, published=pub)
    assert "published" in txt and "DSR (%)" in txt
    # overlay shows both ours and the published DSR value
    assert "2.87" in txt


def test_neutron_report_graceful_without_nesst(monkeypatch):
    monkeypatch.setattr(nsc, "NESST_AVAILABLE", False)
    m, txt = nsc.neutron_report(_BurnRun(), n_E=N_E, include_n2n=False)
    assert m is not None
    assert "quick-look" in m["dsr_source"]
    assert m["DSR"] is not None                      # placeholder DSR still present
    assert "round-trip" not in txt                   # no transport round-trip line


def test_apply_absolute_scale_normalises_to_yield(birth):
    E, S = birth
    # undepleted so the birth primary integrates exactly to the yield
    res = nsc.scatter_from_spectrum(E, S, 0.9, deplete_primary=False, **KW)
    tof = nsc.scattered_tof(res, distance_m=3.0)
    Y = 3.7e17
    dsr_before = res["dsr"]["DSR"]
    f = nsc.apply_absolute_scale(res, tof, Y)
    assert f > 0 and res["absolute"] is True and tof["absolute"] is True
    # birth primary integrates to the yield; DSR (a ratio) is unchanged
    assert np.trapezoid(res["primary"], res["energy_MeV"]) == pytest.approx(Y, rel=1e-6)
    assert res["dsr"]["DSR"] == pytest.approx(dsr_before, rel=1e-9)
    # solid-angle fraction scales linearly
    res2 = nsc.scatter_from_spectrum(E, S, 0.9, deplete_primary=False, **KW)
    nsc.apply_absolute_scale(res2, None, Y, solid_angle_frac=0.5)
    assert np.trapezoid(res2["primary"], res2["energy_MeV"]) == pytest.approx(0.5 * Y, rel=1e-6)


def test_neutron_report_absolute_default(birth):
    # the pipeline scales to the run's DT yield by default
    m, txt = nsc.neutron_report(_BurnRun(), n_E=N_E, include_n2n=False)
    assert m.get("absolute_scale") is True
    assert m.get("dt_yield", 0) > 0


def test_neutron_report_noburn():
    class _NoBurn:
        time = np.array([0.0, 1.0, 2.0])
        zone_boundaries = np.tile(np.linspace(0, 0.03, 5), (3, 1))
        mass_density = np.ones((3, 4))
        ion_temperature = np.full((3, 4), 100.0)
        zone_mass = np.ones((3, 4))
        fusion_power = np.zeros((3, 4))
    m, txt = nsc.neutron_report(_NoBurn())
    assert m is None and "no-burn" in txt


# ----------------------- TT primary spectrum (R-matrix) -----------------------

def test_tt_spectrum_shape_and_endpoint(birth):
    E = np.linspace(0.1, 12.0, 2000)
    tt = nsc.tt_spectrum(E, 10.0, model="Brune")
    assert np.trapezoid(tt, E) == pytest.approx(1.0, rel=1e-6)   # normalised
    # broad continuum peaking a few MeV, endpoint below the 10 MeV DSR window
    assert 2.0 < E[np.argmax(tt)] < 6.0
    assert E[tt > 0.01 * tt.max()].max() < 10.0


def test_tt_dt_neutron_ratio_two_per_reaction():
    from types import SimpleNamespace
    nt, nz = 5, 4
    zm = np.ones((nt, nz)); t = np.linspace(0, 4, nt)
    dtr = np.ones((nt, nz)); ttr = dtr * 0.01           # TT rate 1% of DT
    run = SimpleNamespace(time=t, zone_mass=zm, fusion_power=dtr,
                          fusion_rate_TT_nnHe4=ttr)
    # 2 neutrons per TT reaction -> ratio = 2 * 0.01
    assert nsc.tt_dt_neutron_ratio(run) == pytest.approx(0.02, rel=1e-6)


def test_neutron_report_include_tt_leaves_dsr_unchanged():
    from types import SimpleNamespace
    import test_neutron_spectrum as tns
    base = tns._BurnRun()
    # give the burn run a TT channel (~1% of DT rate)
    base.fusion_rate_TT_nnHe4 = np.asarray(base.fusion_power) * 0.005
    m0, _ = nsc.neutron_report(base, n_E=N_E, include_n2n=False)
    m1, txt = nsc.neutron_report(base, n_E=N_E, include_n2n=False, include_tt=True)
    assert m1.get("tt_dt_ratio", 0) == pytest.approx(0.01, rel=1e-6)   # 2 x 0.005
    # DSR (10-12/13-15 MeV) is above the TT endpoint -> unchanged
    assert m1["DSR"] == pytest.approx(m0["DSR"], rel=1e-6)


# ----------------------- depletion / escaping budget + birth source -----------------------

def test_depletion_gives_physical_escaping_count(birth):
    E, S = birth
    r_dep = nsc.scatter_from_spectrum(E, S, 1.0, **KW)                       # depleted (default)
    r_raw = nsc.scatter_from_spectrum(E, S, 1.0, deplete_primary=False, **KW)
    Eo = r_dep["energy_MeV"]
    # depleted primary is attenuated below the raw primary
    assert np.trapezoid(r_dep["primary"], Eo) < np.trapezoid(r_raw["primary"], Eo)
    # the un-depleted model over-counts (double counts the scattered neutrons)
    assert r_raw["escaping_fraction"] > r_dep["escaping_fraction"] + 0.1
    # (n,2n) multiplication pushes the depleted escaping count above the yield
    r_n2n = nsc.scatter_from_spectrum(E, S, 1.0, n_E=N_E, include_n2n=True)
    assert r_n2n["escaping_fraction"] > 1.0


def test_birth_source_channel_yields():
    import test_neutron_spectrum as tns
    run = tns._BurnRun()
    run.fusion_rate_TT_nnHe4 = np.asarray(run.fusion_power) * 0.004      # TT ~ DT*0.008 neutrons
    nd = ns.extract_neutronics(data=run, use_rhino=False)
    src = nsc.birth_source(nd, run, n_E=1500)
    assert src["Y_DT"] > 0
    assert src["Y_TT"] == pytest.approx(src["tt_dt_ratio"] * src["Y_DT"], rel=1e-6)
    assert src["total_neutrons"] == pytest.approx(src["Y_DT"] + src["Y_DD"] + src["Y_TT"], rel=1e-9)
    # the source spectrum integrates to the total born neutrons
    assert np.trapezoid(src["spectrum"], src["energy_MeV"]) == pytest.approx(
        src["total_neutrons"], rel=0.05)
