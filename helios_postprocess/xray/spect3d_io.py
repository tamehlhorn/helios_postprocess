"""
SPECT3D <-> helios_postprocess adapter.

This module holds BOTH halves of the hand-off:

  * ``export_profiles``   -- write Helios hydro snapshots in a plain,
    documented ASCII form that SPECT3D can ingest if the native Helios
    coupling is unavailable or inconvenient to script.

  * ``load_spect3d_*``    -- read SPECT3D output back onto the SAME
    ``RadianceCube`` contract that the Level 0 model produces, so the
    comparison is a diff rather than a re-implementation.

STATUS: the reader is written against the generic SPECT3D ASCII export
(2-column spectra, 2-D image dumps, and a simple space-energy table).  The
exact column headers vary with SPECT3D version and with what is checked in
the output dialog, so ``load_spect3d_cube`` takes explicit column-index
arguments and prints what it parsed.  Confirm against one real export before
trusting a batch.

Author: Prof T / helios_postprocess
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .streak import RadianceCube

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Export: Helios -> SPECT3D
# ---------------------------------------------------------------------------
def export_profiles(data, times_ns: Sequence[float], outdir,
                    *, zone_slice: Optional[slice] = None,
                    prefix: str = "helios_profile") -> List[Path]:
    """
    Write one ASCII profile per requested time.

    Columns (SI-free, CGS + eV, the units SPECT3D expects for a 1D spherical
    Lagrangian dump):

        r_inner_cm  r_outer_cm  rho_g_cm3  Te_eV  Ti_eV  Zbar  n_e_cm3

    A ``times.txt`` index file is also written so the batch driver knows the
    mapping from file to simulation time.

    Parameters
    ----------
    data : ICFRunData
    times_ns : sequence of float
        Requested times; the NEAREST available Helios step is used and the
        actual time is written into the header (Helios output is not uniform
        in time, so requested != actual in general).
    outdir : path
    zone_slice : slice, optional
        Restrict to the capsule for halfraum targets.
    """
    from .emissivity import (resolve_Te, resolve_zbar,
                             resolve_electron_density)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t = np.asarray(data.time, float)
    written: List[Path] = []
    index_lines = ["# file  requested_ns  actual_ns  time_index"]

    for k, treq in enumerate(times_ns):
        it = int(np.argmin(np.abs(t - float(treq))))

        rb = np.asarray(data.zone_boundaries[it], float)
        rho = np.asarray(data.mass_density[it], float)
        Te = resolve_Te(data, it)
        Ti = np.asarray(data.ion_temperature[it], float)
        zbar = resolve_zbar(data, it)
        ne = resolve_electron_density(data, it)

        if zone_slice is not None:
            rho, Te, Ti, zbar, ne = (rho[zone_slice], Te[zone_slice],
                                     Ti[zone_slice], zbar[zone_slice],
                                     ne[zone_slice])
            start = 0 if zone_slice.start is None else zone_slice.start
            stop = len(rb) - 1 if zone_slice.stop is None else zone_slice.stop
            rb = rb[start:stop + 1]

        fname = outdir / f"{prefix}_{k:03d}.dat"
        header = (f"# Helios profile export for SPECT3D\n"
                  f"# requested_time_ns = {float(treq):.6f}\n"
                  f"# actual_time_ns    = {t[it]:.6f}\n"
                  f"# time_index        = {it}\n"
                  f"# n_zones           = {rho.size}\n"
                  f"# columns: r_in_cm r_out_cm rho_g_cm3 Te_eV Ti_eV Zbar ne_cm3\n")
        block = np.column_stack([rb[:-1], rb[1:], rho, Te, Ti, zbar, ne])
        with open(fname, "w") as fh:
            fh.write(header)
            np.savetxt(fh, block, fmt="%.8e")
        written.append(fname)
        index_lines.append(f"{fname.name}  {float(treq):.6f}  {t[it]:.6f}  {it}")

    (outdir / "times.txt").write_text("\n".join(index_lines) + "\n")
    logger.info(f"[xray] exported {len(written)} profiles to {outdir}")
    return written


# ---------------------------------------------------------------------------
# Import: SPECT3D -> RadianceCube
# ---------------------------------------------------------------------------
def _read_numeric(path: Path, comment_chars: str = "#!;") -> np.ndarray:
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.strip().replace(",", " ")
        if not s or s[0] in comment_chars:
            continue
        parts = s.split()
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"No numeric rows parsed from {path}")
    width = min(len(r) for r in rows)
    return np.array([r[:width] for r in rows], dtype=float)


def load_spect3d_spectrum(path, energy_col: int = 0, value_col: int = 1,
                          energy_unit: str = "eV") -> Tuple[np.ndarray, np.ndarray]:
    """
    Read a single SPECT3D spectrum export.

    Returns (E_eV, value).  ``energy_unit`` may be 'eV' or 'keV'.
    """
    arr = _read_numeric(Path(path))
    E = arr[:, energy_col]
    if energy_unit.lower() == "kev":
        E = E * 1.0e3
    v = arr[:, value_col]
    order = np.argsort(E)
    return E[order], v[order]


def load_spect3d_cube(files: Sequence,
                      times_ns: Sequence[float],
                      *,
                      p_cm: Optional[np.ndarray] = None,
                      space_col: int = 0,
                      energy_col: int = 1,
                      value_col: int = 2,
                      energy_unit: str = "eV",
                      space_unit: str = "cm") -> RadianceCube:
    """
    Assemble a ``RadianceCube`` from per-time SPECT3D space-energy exports.

    Each file is expected to hold a long-format table with one row per
    (spatial position, photon energy) pair.  The routine pivots it into a
    (n_p, n_E) grid.  If the export is spectrum-only (no spatial column),
    pass ``space_col=None`` and a single-element ``p_cm``.

    The result carries zeros in ``tau`` and a copy of ``I`` in ``I_thin``:
    SPECT3D's own opacity treatment is inside ``I``, and the thin reference
    for the comparison comes from the Level 0 cube, not from here.
    """
    times_ns = np.asarray(times_ns, float)
    if len(files) != times_ns.size:
        raise ValueError(f"{len(files)} files vs {times_ns.size} times")

    planes, E_ref, p_ref = [], None, None

    for path in files:
        arr = _read_numeric(Path(path))

        if space_col is None:
            E = arr[:, energy_col]
            v = arr[:, value_col]
            if energy_unit.lower() == "kev":
                E = E * 1.0e3
            order = np.argsort(E)
            E, v = E[order], v[order]
            p_vals = np.array([0.0]) if p_cm is None else np.asarray(p_cm, float)
            plane = v[None, :]
        else:
            x = arr[:, space_col]
            E = arr[:, energy_col]
            v = arr[:, value_col]
            if energy_unit.lower() == "kev":
                E = E * 1.0e3
            if space_unit.lower() in ("um", "micron", "microns"):
                x = x * 1.0e-4
            elif space_unit.lower() == "mm":
                x = x * 0.1
            p_vals = np.unique(x)
            E_vals = np.unique(E)
            plane = np.full((p_vals.size, E_vals.size), np.nan)
            ip = np.searchsorted(p_vals, x)
            ie = np.searchsorted(E_vals, E)
            plane[ip, ie] = v
            plane = np.nan_to_num(plane, nan=0.0)
            E = E_vals

        if E_ref is None:
            E_ref, p_ref = E, p_vals
        else:
            if plane.shape[1] != E_ref.size:
                plane = np.array([np.interp(E_ref, E, row) for row in plane])
            if plane.shape[0] != p_ref.size:
                plane = np.array([np.interp(p_ref, p_vals, plane[:, j])
                                  for j in range(plane.shape[1])]).T
        planes.append(plane)

    I = np.stack(planes, axis=0)
    logger.info(f"[xray] SPECT3D cube: {I.shape} (t, p, E); "
                f"E {E_ref[0]:.1f}-{E_ref[-1]:.1f} eV")
    return RadianceCube(t_ns=times_ns, p_cm=np.asarray(p_ref, float),
                        E_eV=np.asarray(E_ref, float),
                        I=I, tau=np.zeros_like(I), I_thin=I.copy(),
                        meta=dict(source="spect3d", n_files=len(files)))


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------
def compare_cubes(level0: RadianceCube, spect3d: RadianceCube,
                  channel=None) -> Dict[str, float]:
    """
    Quantitative Level 0 vs SPECT3D comparison on the common grid.

    Reports the band-integrated power ratio, the spectral-shape difference,
    and the shift in x-ray bang time -- i.e. exactly the three ways the extra
    physics can change a conclusion drawn from the streak record.
    """
    from .streak import burn_metrics, StreakImage

    t = spect3d.t_ns
    E = spect3d.E_eV

    p0 = level0.power_spectrum()
    ps = spect3d.power_spectrum()

    p0_i = np.array([np.interp(E, level0.E_eV, row) for row in p0])
    p0_i = np.array([np.interp(t, level0.t_ns, p0_i[:, j])
                     for j in range(E.size)]).T

    R = None if channel is None else channel.response(E)
    w = np.ones_like(E) if R is None else R

    band0 = np.trapezoid(p0_i * w[None, :], E, axis=1)
    bands = np.trapezoid(ps * w[None, :], E, axis=1)

    img0 = StreakImage("spectral", t, E, "eV", p0_i * w[None, :])
    imgs = StreakImage("spectral", t, E, "eV", ps * w[None, :])
    m0, ms = burn_metrics(img0), burn_metrics(imgs)

    total0 = np.trapezoid(band0, t)
    totals = np.trapezoid(bands, t)

    return dict(
        power_ratio_spect3d_over_level0=float(totals / total0) if total0 else np.nan,
        peak_ratio=float(bands.max() / band0.max()) if band0.max() else np.nan,
        bang_time_level0_ns=m0.bang_time_ns,
        bang_time_spect3d_ns=ms.bang_time_ns,
        bang_time_shift_ps=float((ms.bang_time_ns - m0.bang_time_ns) * 1e3),
        fwhm_level0_ps=m0.fwhm_ps,
        fwhm_spect3d_ps=ms.fwhm_ps,
    )
