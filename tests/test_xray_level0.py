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


# ---------------------------------------------------------------------------
# Density floor (outer-corona mask)
# ---------------------------------------------------------------------------
def test_rho_floor_masks_emission_and_absorption(mock):
    """A masked zone must neither emit nor absorb."""
    E = np.logspace(2, 4, 20)
    rho = mock.mass_density[20]
    floor = float(np.median(rho))
    c = build_emissivity(mock, 20, E, include_absorption=True, rho_floor=floor)
    masked = rho < floor
    assert masked.any() and (~masked).any()
    assert np.all(c.j_E[masked] == 0.0)
    assert np.all(c.kappa_E[masked] == 0.0)
    assert np.any(c.j_E[~masked] > 0.0)


class CoronaMock:
    """
    Minimal two-population profile: a dense 3 keV bulk plus a tenuous 80 keV
    outer corona, i.e. the configuration the recon found on the 4 MJ Vulcan
    run.  Built explicitly rather than reusing MockRunData, whose outer zones
    are cool -- the artifact under test requires the outer zones to be BOTH
    low density AND very hot.
    """

    def __init__(self, n_bulk=100, n_cor=10, rho_corona=1.0e-3):
        n_z = n_bulk + n_cor
        self.time = np.array([0.0, 1.0])
        rb = np.linspace(0.0, 0.05, n_bulk + 1)
        rb = np.concatenate([rb, np.linspace(0.05, 0.30, n_cor + 1)[1:]])
        self.zone_boundaries = np.tile(rb, (2, 1))

        rho = np.concatenate([np.full(n_bulk, 1.0),
                              np.full(n_cor, float(rho_corona))])
        Te = np.concatenate([np.full(n_bulk, 3.0e3), np.full(n_cor, 8.0e4)])
        zbar = np.concatenate([np.full(n_bulk, 3.5), np.full(n_cor, 3.5)])
        ni = rho / (2.5 * 1.6726e-24)

        self.mass_density = np.tile(rho, (2, 1))
        self.elec_temperature = np.tile(Te, (2, 1))
        self.ion_temperature = np.tile(Te, (2, 1))
        self.ion_density = np.tile(ni, (2, 1))
        self.electron_density = np.tile(ni * zbar, (2, 1))
        self.mean_charge = np.tile(zbar, (2, 1))


def test_hot_tenuous_corona_dominates_the_hard_tail():
    """
    The physics the density floor exists to control: a corona carrying a
    negligible share of the emission at low photon energy can take over the
    hard tail, because the cool bulk is exponentially suppressed there while
    the corona is not.
    """
    m = CoronaMock()
    E = np.array([500.0, 30000.0])
    c = build_emissivity(m, 0, E, include_absorption=False)
    vol = (4.0 / 3.0) * np.pi * (m.zone_boundaries[0, 1:] ** 3
                                 - m.zone_boundaries[0, :-1] ** 3)
    power = c.j_E * vol[:, None]
    cor_soft = power[100:, 0].sum() / power[:, 0].sum()
    cor_hard = power[100:, 1].sum() / power[:, 1].sum()
    assert cor_soft < 1e-3, f"corona already visible at 500 eV ({cor_soft:.3g})"
    assert cor_hard > 0.5, f"corona should own the 30 keV tail ({cor_hard:.3f})"


def test_corona_takeover_is_density_controlled():
    """
    Whether the corona takes over the hard tail is set by (n_c/n_b)^2 against
    the exponential advantage exp(E/kT_b - E/kT_c).  At low enough coronal
    density the n^2 penalty wins and there is no artifact at all -- so the
    presence of a very hot outer zone is NOT by itself evidence of
    contamination.  This is why the provenance script decides empirically
    rather than from the Te profile alone.
    """
    E = np.array([30000.0])

    def hard_fraction(rho_c):
        m = CoronaMock(rho_corona=rho_c)
        c = build_emissivity(m, 0, E, include_absorption=False)
        vol = (4.0 / 3.0) * np.pi * (m.zone_boundaries[0, 1:] ** 3
                                     - m.zone_boundaries[0, :-1] ** 3)
        power = c.j_E * vol[:, None]
        return power[100:, 0].sum() / power[:, 0].sum()

    assert hard_fraction(1.0e-6) < 0.01
    assert hard_fraction(1.0e-4) < 0.10
    assert hard_fraction(1.0e-3) > 0.50


def test_rho_floor_reduces_hard_tail_more_than_soft():
    """Masking the tenuous corona must suppress the hard end far more than
    the soft end -- that asymmetry is the whole point of the floor."""
    m = CoronaMock()
    E = np.array([500.0, 30000.0])
    vol = (4.0 / 3.0) * np.pi * (m.zone_boundaries[0, 1:] ** 3
                                 - m.zone_boundaries[0, :-1] ** 3)
    full = build_emissivity(m, 0, E, include_absorption=False)
    cut = build_emissivity(m, 0, E, include_absorption=False, rho_floor=1.0e-2)
    p_full = (full.j_E * vol[:, None]).sum(axis=0)
    p_cut = (cut.j_E * vol[:, None]).sum(axis=0)
    keep_soft = p_cut[0] / p_full[0]
    keep_hard = p_cut[1] / p_full[1]
    assert keep_soft > 0.95
    assert keep_hard < 0.5
    assert keep_hard < keep_soft


def test_rho_floor_propagates_through_config(mock):
    from helios_postprocess.xray import StreakConfig, build_radiance_cube
    cfg = StreakConfig(n_time=40, n_impact=16, n_energy=12,
                       rho_floor_g_cc=float(np.median(mock.mass_density[20])))
    cube = build_radiance_cube(mock, cfg, verbose=False)
    assert cube.meta["rho_floor_g_cc"] == cfg.rho_floor_g_cc
    assert np.all(np.isfinite(cube.I))
