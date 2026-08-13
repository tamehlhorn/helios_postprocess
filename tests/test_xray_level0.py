"""
Level 0 x-ray forward-model tests.

Runs entirely on an analytic mock implosion -- no .exo required -- so the
physics kernels can be verified on a MacBook without simulation data.

    python -m pytest tests/test_xray_level0.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

from helios_postprocess.xray import (build_emissivity, gaunt_ff, planck_E,
                                     integrate_thin, integrate_formal,
                                     make_impact_grid, StreakConfig,
                                     build_radiance_cube, make_imaging_streak,
                                     make_spectral_streak, default_channels,
                                     burn_metrics, emission_edge_trajectory,
                                     mu_over_rho_element)


# ---------------------------------------------------------------------------
# Mock run
# ---------------------------------------------------------------------------
class MockRunData:
    """
    Analytic 1D spherical implosion: a shell of fixed mass converging on a
    hot core, with a hot low-density corona outside.  Not a real Helios run;
    just enough structure to exercise every code path and give an emission
    history with a well-defined peak.
    """

    def __init__(self, n_t=40, n_z=200):
        t = np.linspace(0.0, 2.0, n_t)                      # ns
        R0, Rmin = 0.05, 0.005                              # cm
        # smooth converge-and-bounce trajectory
        R = Rmin + (R0 - Rmin) * np.cos(np.pi * t / 2.0 / t[-1] * 2.0) ** 2
        R = np.clip(R, Rmin, R0)

        self.time = t
        zb = np.empty((n_t, n_z + 1))
        rho = np.empty((n_t, n_z))
        Te = np.empty((n_t, n_z))
        Ti = np.empty((n_t, n_z))
        ne = np.empty((n_t, n_z))
        ni = np.empty((n_t, n_z))
        zbar = np.empty((n_t, n_z))

        for i in range(n_t):
            outer = max(R[i] * 3.0, 0.02)
            zb[i] = np.linspace(0.0, outer, n_z + 1)
            rc = 0.5 * (zb[i, :-1] + zb[i, 1:])

            conv = R0 / max(R[i], Rmin)
            core = np.exp(-(rc / (0.4 * R[i])) ** 2)
            shell = np.exp(-((rc - R[i]) / (0.15 * R[i])) ** 2)

            rho[i] = 1e-4 + 0.25 * conv ** 2 * rho_shape(shell) + 0.02 * core
            Te[i] = 50.0 + 3000.0 * core * conv ** 1.5 + 800.0 * np.exp(
                -((rc - 2.0 * R[i]) / R[i]) ** 2)
            Ti[i] = Te[i] * 1.05
            zbar[i] = np.clip(1.0 + 5.0 * (Te[i] / 1000.0) ** 0.4, 1.0, 6.0)
            ni[i] = rho[i] / (2.5 * 1.6726e-24)
            ne[i] = ni[i] * zbar[i]

        self.zone_boundaries = zb
        self.mass_density = rho
        self.elec_temperature = Te
        self.ion_temperature = Ti
        self.electron_density = ne
        self.ion_density = ni
        self.mean_charge = zbar


def rho_shape(x):
    return x


@pytest.fixture(scope="module")
def mock():
    return MockRunData()


# ---------------------------------------------------------------------------
# Kernels
# ---------------------------------------------------------------------------
def test_gaunt_limits():
    """Gaunt factor is O(1-10), rises logarithmically as u -> 0, falls at u >> 1."""
    E = np.array([10.0, 100.0, 1000.0, 10000.0])
    g = gaunt_ff(E, 1000.0)
    assert np.all(np.isfinite(g))
    assert np.all(g > 0)
    assert g[0] > g[-1]                    # decreasing with u
    assert 1.0 < g[1] < 10.0


def test_planck_wien_and_rj():
    """Planck function has the right limiting behaviour."""
    Te = 1000.0
    E = np.array([10.0, 1000.0, 20000.0])
    B = planck_E(E, Te)
    assert np.all(B > 0)
    # Wien tail must be steeply suppressed
    assert B[2] / B[1] < 1e-4


def test_kirchhoff_source_is_planck(mock):
    """j/kappa must return the Planck function to machine precision."""
    E = np.logspace(2, 4, 24)
    cube = build_emissivity(mock, 10, E, include_absorption=True)
    good = cube.kappa_E > 0
    ratio = cube.j_E[good] / cube.kappa_E[good]
    expected = np.broadcast_to(cube.source_E, cube.j_E.shape)[good]
    assert np.allclose(ratio, expected, rtol=1e-8)


def test_emissivity_scales_with_ne_ni_zbar(mock):
    """Free-free emissivity must scale as Zbar^2 n_e n_i."""
    E = np.array([1000.0])
    c1 = build_emissivity(mock, 10, E)
    mock2 = MockRunData()
    mock2.electron_density = mock2.electron_density * 2.0
    c2 = build_emissivity(mock2, 10, E)
    good = c1.j_E[:, 0] > 0
    assert np.allclose(c2.j_E[good, 0] / c1.j_E[good, 0], 2.0, rtol=1e-10)


# ---------------------------------------------------------------------------
# Chord geometry
# ---------------------------------------------------------------------------
def test_chord_length_uniform_sphere():
    """
    Uniform emissivity in a sphere of radius a must give the analytic
    optically-thin chord integral I(p) = 2 j sqrt(a^2 - p^2).
    """
    a = 1.0
    n_z = 400
    r_bnd = np.linspace(0.0, a, n_z + 1)
    j = np.ones((n_z, 1))

    class C:
        pass
    c = C()
    c.r_bnd, c.j_E = r_bnd, j
    c.kappa_E = np.zeros_like(j)
    c.source_E = np.zeros_like(j)
    c.E_eV = np.array([1000.0])

    p = np.array([0.0, 0.25, 0.5, 0.75, 0.95])
    res = integrate_thin(c, p)
    exact = 2.0 * np.sqrt(a ** 2 - p ** 2)
    assert np.allclose(res.I[:, 0], exact, rtol=2e-3)


def test_formal_matches_thin_at_low_tau(mock):
    """Wherever the chord is optically thin the formal solution must reduce
    to the thin sum.  (The mock is deliberately dense, so soft-photon chords
    are thick -- the test selects the thin region rather than assuming it.)"""
    E = np.logspace(3, 4.3, 16)
    cube = build_emissivity(mock, 20, E, include_absorption=True)
    p = make_impact_grid(cube.r_bnd[-1], 32)
    thin = integrate_thin(cube, p)
    formal = integrate_formal(cube, p)
    m = (formal.tau < 1e-3) & (thin.I > 0)
    assert m.sum() > 50, "no optically-thin samples to compare"
    assert np.allclose(formal.I[m], thin.I[m], rtol=1e-3)


def test_formal_never_exceeds_planck(mock):
    """Emergent intensity cannot exceed the peak source function on the chord."""
    E = np.logspace(2, 3.5, 16)
    cube = build_emissivity(mock, 20, E, include_absorption=True)
    p = make_impact_grid(cube.r_bnd[-1], 24)
    formal = integrate_formal(cube, p)
    Smax = cube.source_E.max(axis=0)
    assert np.all(formal.I <= Smax[None, :] * (1.0 + 1e-9))


def test_absorption_reduces_intensity(mock):
    """Formal solution must be <= thin solution everywhere."""
    E = np.logspace(2, 3, 10)
    cube = build_emissivity(mock, 20, E, include_absorption=True)
    p = make_impact_grid(cube.r_bnd[-1], 24)
    formal = integrate_formal(cube, p)
    assert np.all(formal.I <= formal.I_thin * (1.0 + 1e-6))


# ---------------------------------------------------------------------------
# Streak assembly
# ---------------------------------------------------------------------------
def test_streak_end_to_end(mock):
    cfg = StreakConfig(n_time=120, n_impact=48, n_energy=32,
                       optically_thin=False)
    cube = build_radiance_cube(mock, cfg, verbose=False)
    assert cube.I.shape == (mock.time.size, 48, 32)

    ch = default_channels()["mid"]
    img = make_imaging_streak(cube, ch, cfg)
    assert img.signal.shape == (120, 48)
    assert np.all(np.isfinite(img.signal))

    spec = make_spectral_streak(cube, cfg)
    assert spec.signal.shape == (120, 32)

    m = burn_metrics(img)
    assert np.isfinite(m.bang_time_ns)
    assert cube.t_ns[0] <= m.bang_time_ns <= cube.t_ns[-1]

    edge = emission_edge_trajectory(img)
    assert np.isfinite(edge).sum() > 0.5 * edge.size


def test_edge_trajectory_converges(mock):
    """The emission edge must move inward as the mock shell converges."""
    cfg = StreakConfig(n_time=200, n_impact=64, n_energy=24)
    cube = build_radiance_cube(mock, cfg, verbose=False)
    img = make_imaging_streak(cube, default_channels()["soft"], cfg)
    edge = emission_edge_trajectory(img)
    finite = np.isfinite(edge)
    first = np.nanmean(edge[finite][:10])
    minimum = np.nanmin(edge[finite])
    assert minimum < first


def test_time_base_is_uniform(mock):
    """Streak time base must be uniform even though Helios output is not."""
    cfg = StreakConfig(n_time=100, n_impact=16, n_energy=16)
    cube = build_radiance_cube(mock, cfg, verbose=False)
    img = make_spectral_streak(cube, cfg)
    dt = np.diff(img.t_ns)
    assert np.allclose(dt, dt[0], rtol=1e-10)


# ---------------------------------------------------------------------------
# Response
# ---------------------------------------------------------------------------
def test_filter_transmission_monotonic():
    """Transmission must increase with photon energy away from edges."""
    from helios_postprocess.xray import be_filter
    E = np.array([500.0, 1000.0, 2000.0, 5000.0, 10000.0])
    T = be_filter(25.0).transmission(E)
    assert np.all(np.diff(T) > 0)
    assert np.all((T >= 0) & (T <= 1))


def test_mu_decreases_with_energy():
    E = np.array([1000.0, 3000.0, 9000.0])
    mu = mu_over_rho_element("Be", E)
    assert np.all(np.diff(mu) < 0)


def test_channel_response_bounded():
    for ch in default_channels().values():
        E = np.logspace(2, 4.2, 60)
        R = ch.response(E)
        assert np.all(R >= 0.0)
        assert np.all(R <= 1.0)
