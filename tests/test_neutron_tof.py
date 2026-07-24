"""Analytic unit tests for helios_postprocess.neutron_tof."""
import numpy as np
import pytest

from helios_postprocess import neutron_tof as tof
from helios_postprocess import neutron_spectrum as ns


def test_neutron_velocity_and_tof():
    assert tof.neutron_velocity(14.06) == pytest.approx(5.129, rel=1e-3)
    assert tof.neutron_velocity(2.45) == pytest.approx(2.161, rel=1e-3)
    assert float(tof.energy_to_tof(14.06, 3.0)) == pytest.approx(58.5, rel=2e-3)
    assert float(tof.energy_to_tof(14.06, 10.0)) == pytest.approx(195.0, rel=2e-3)


def test_energy_tof_roundtrip():
    E = np.array([2.45, 3.53, 10.0, 14.06])
    t = tof.energy_to_tof(E, 3.0, t0_ns=1.5)
    assert np.allclose(tof.tof_to_energy(t, 3.0, t0_ns=1.5), E, rtol=1e-9)


def test_spectrum_to_tof_conserves_counts():
    E = np.linspace(11.0, 17.0, 4000)
    dNdE = np.exp(-0.5 * ((E - 14.06) / 0.2) ** 2)
    N_E = np.trapezoid(dNdE, E)
    t, dNdt = tof.spectrum_to_tof(E, dNdE, 3.0)
    assert np.trapezoid(dNdt, t) == pytest.approx(N_E, rel=1e-3)
    # arrival time ordering: higher E arrives earlier
    assert t.min() > 0 and np.all(np.diff(t) >= 0)


def test_tion_recovered_from_tof_width():
    Ti = 8.0
    E, S = ns.synthesize_birth_spectrum(np.array([[1.0]]), np.array([[Ti * 1000.0]]),
                                        reaction="DT")
    t, dNdt = tof.spectrum_to_tof(E, S, 3.0)
    m = tof.primary_peak_metrics(t, dNdt, 3.0, "DT")
    assert m["t_peak_ns"] == pytest.approx(58.5, rel=3e-3)
    # single clean Gaussian: the linear TOF-width estimate is good to ~6%
    assert m["Ti_tof_width_keV"] == pytest.approx(Ti, rel=0.07)


def test_apply_irf_preserves_area_and_broadens():
    t = np.linspace(50, 67, 4000)
    y = np.exp(-0.5 * ((t - 58.5) / 0.3) ** 2)
    yc = tof.apply_irf(t, y, irf_fwhm_ns=1.0)
    assert np.trapezoid(yc, t) == pytest.approx(np.trapezoid(y, t), rel=1e-3)
    assert tof._fwhm(t, yc) > tof._fwhm(t, y)          # broadened


def test_synthetic_ntof_from_spectrum():
    Ti = 9.0
    E, S = ns.synthesize_birth_spectrum(np.array([[1.0]]), np.array([[Ti * 1000.0]]),
                                        reaction="DT")
    out = tof.synthetic_ntof(energy_MeV=E, spectrum=S, distance_m=3.0, irf_fwhm_ns=0.5)
    assert out["primary"]["t_peak_ns"] == pytest.approx(58.5, rel=3e-3)
    assert out["primary"]["Ti_ntof_keV"] == pytest.approx(Ti, rel=0.06)
    assert out["signal_detector"].shape == out["time_ns"].shape
    assert np.isfinite(out["signal_ideal"]).all()
