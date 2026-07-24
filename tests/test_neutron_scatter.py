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
    dsrs = [nsc.scatter_from_spectrum(E, S, r, **KW)["dsr"]["DSR"] for r in rhoRs]
    # strictly increasing
    assert all(b > a for a, b in zip(dsrs, dsrs[1:]))
    # near-linear: DSR/rhoR (the slope) varies by < 10% across the range
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
    res = nsc.scatter_from_spectrum(E, S, 0.9, **KW)
    assert np.allclose(res["full"], res["primary"] + res["scattered"])
    # yield preserved: primary integral == input yield
    Y_in = np.trapezoid(S, E)
    Y_out = np.trapezoid(res["primary"], res["energy_MeV"])
    assert Y_out == pytest.approx(Y_in, rel=0.02)


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
