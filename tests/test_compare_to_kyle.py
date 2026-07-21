"""Unit tests for the diff math in examples/compare_to_kyle_npz.py.

Exercises the pure comparison helpers (no .exo / pipeline needed).
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

_MOD = Path(__file__).resolve().parent.parent / "examples" / "compare_to_kyle_npz.py"
_spec = importlib.util.spec_from_file_location("compare_to_kyle_npz", _MOD)
cmp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cmp)


def test_reldiff_scalar():
    assert cmp.reldiff_scalar(10.0, 8.0) == pytest.approx(0.25)
    assert cmp.reldiff_scalar(0.0, 0.0) == 0.0
    assert cmp.reldiff_scalar(5.0, 0.0) == float("inf")
    assert cmp.reldiff_scalar(8.0, 10.0) == pytest.approx(-0.2)


def test_reldiff_profile_identical_and_scaled():
    r = np.linspace(0, 0.03, 200)
    y = np.exp(-((r - 0.01) / 0.003) ** 2)
    same = cmp.reldiff_profile(r, y, r, y)
    assert same["max"] == pytest.approx(0.0, abs=1e-9)
    scaled = cmp.reldiff_profile(r, 1.1 * y, r, y)
    assert scaled["max"] == pytest.approx(0.10, rel=1e-6)
    assert scaled["mean"] == pytest.approx(0.10, rel=1e-6)


def test_reldiff_profile_regrids():
    # ours on a finer grid, same underlying function -> ~0 after interp
    ra = np.linspace(0, 0.03, 500); rb = np.linspace(0, 0.03, 137)
    fa = 2.0 + ra; fb = 2.0 + rb
    res = cmp.reldiff_profile(ra, fa, rb, fb)
    assert res["max"] < 1e-3


def test_compare_npz_rows():
    r = np.linspace(0, 0.03, 100)
    rho = np.exp(-((r - 0.01) / 0.004) ** 2)
    Ti = 9000.0 * np.exp(-((r - 0.008) / 0.004) ** 2)
    ours = {
        "bang_time_ns": 23.6, "total_dt_yield": 1.29e18, "total_dd_yield": 5.0e15,
        "r_avg_cm": r, "rho_avg_gcc": 1.05 * rho, "T_ion_avg_eV": Ti,
        "rate_DT_avg": rho,
        "avg_results": {"DT_nHe4": {"rho_avg_gcc": 1.05 * rho, "T_ion_avg_eV": Ti}},
    }
    theirs = {
        "bang_time_ns": 23.6, "total_dt_yield": 1.29e18, "total_dd_yield": 5.0e15,
        "r_avg_cm": r, "rho_avg_gcc": rho, "T_ion_avg_eV": Ti, "rate_DT_avg": rho,
        "avg_results": {"DT_nHe4": {"rho_avg_gcc": rho, "T_ion_avg_eV": Ti}},
    }
    rows = cmp.compare_npz(ours, theirs)
    labels = {r[0]: r for r in rows}
    assert labels["bang time (ns)"][2] == pytest.approx(0.0, abs=1e-9)
    assert labels["DT yield"][2] == pytest.approx(0.0, abs=1e-9)
    # rho profile 5% high in ours
    assert labels["<rho(r)>_DT (g/cc)"][2]["max"] == pytest.approx(0.05, rel=1e-6)
    # T_ion identical
    assert labels["<T_ion(r)>_DT (eV)"][2]["max"] == pytest.approx(0.0, abs=1e-9)
    # report renders without error
    assert "bang time" in cmp.format_report(rows)
