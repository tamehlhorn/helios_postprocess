"""
qhat.py
=======
Fraction of alpha-particle energy deposited in an ICF hot spot.

Python port of Qhat.pro (IDL). Follows Atzeni & Meyer-ter-Vehn,
"The Physics of Inertial Fusion" (Oxford, 2004):

  - Eq. 4.7  : piecewise deposition fraction f(tau)
  - Eq. 4.8  : tau = R_h / L_alpha  (equivalently rhoR / rho*lambda)
  - Eq. 10.136: ion-electron Coulomb log (present in source, unused in result)

The alpha range is evaluated via the Krokhin-Rozanov fit
(rho * lambda_alpha in g/cm^2).

Inputs are broadcastable (scalar or np.ndarray). Returns the deposited
energy fraction in [0, 1].
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
AMU_G = 1.6605e-24          # atomic mass unit [g]
A_DT = 2.5                   # mean atomic mass number of 50/50 D-T


# ---------------------------------------------------------------------------
# Krokhin-Rozanov alpha rho-range fit
# ---------------------------------------------------------------------------
def rho_lambda_alpha(rho: ArrayLike, Te: ArrayLike) -> np.ndarray:
    """
    Krokhin-Rozanov fit for alpha-particle rho-range in DT [g/cm^2].

    Parameters
    ----------
    rho : array-like
        Mass density [g/cm^3].
    Te : array-like
        Electron temperature [keV].

    Returns
    -------
    np.ndarray
        rho * lambda_alpha [g/cm^2].
    """
    rho = np.asarray(rho, dtype=np.float64)
    Te = np.asarray(Te, dtype=np.float64)
    return (
        0.03 * Te
        * (1.0 - 0.24 * np.log10(1.0 + Te))
        * (1.0 + 0.37 * np.log10((1.0 + rho) / (1.0 + 0.01 * Te**2)))
    )


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------
def qhat(
    rho: ArrayLike,
    rhoR: ArrayLike,
    Te: ArrayLike,
) -> np.ndarray:
    """
    Fraction of alpha-particle energy deposited in the hot spot.

    Parameters
    ----------
    rho : array-like
        Hot-spot mass density [g/cm^3].
    rhoR : array-like
        Hot-spot areal density [g/cm^2].
    Te : array-like
        Electron temperature [keV].

    Returns
    -------
    np.ndarray
        Deposited fraction, in [0, 1]. Shape is the broadcast of the inputs.

    Notes
    -----
    The original IDL source also computes electron density

        eden = rho / (amu * A_DT)          # [cm^-3], Z = 1

    and the ion-electron Coulomb logarithm

        logLam_ie = 7.1 - 0.5*log(eden/1e21) + log(Te)     # Atzeni Eq. 10.136

    but neither enters the final fraction. They are omitted here; uncomment
    the commented block below if you want them returned for diagnostics.
    """
    rho = np.asarray(rho, dtype=np.float64)
    rhoR = np.asarray(rhoR, dtype=np.float64)
    Te = np.asarray(Te, dtype=np.float64)

    # --- diagnostic quantities from the original source (unused in result) ---
    # eden = rho / (AMU_G * A_DT)
    # logLam_ie = 7.1 - 0.5 * np.log(eden / 1.0e21) + np.log(Te)
    # Rh = rhoR / rho                       # hot spot radius [cm]
    # L_alpha = 0.107 * Te**1.5 / (rho * logLam_ie)   # alt. range [cm]

    # tau = R_h / L_alpha, computed via the rho-range fit (Eq. 4.8)
    rho_lambda = rho_lambda_alpha(rho, Te)
    tau = rhoR / rho_lambda

    # Piecewise deposition fraction (Atzeni & Meyer-ter-Vehn Eq. 4.7)
    small_tau = 1.5 * tau - 0.8 * tau**2
    large_tau = 1.0 - 0.25 / tau + 1.0 / (160.0 * tau**3)

    fraction = np.where(tau <= 0.5, small_tau, large_tau)

    return fraction


# ---------------------------------------------------------------------------
# Demo / sanity check
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    # Single-point sanity check: marginally igniting hot spot
    rho_c, rhoR_c, Te_c = 100.0, 0.30, 10.0
    f_c = qhat(rho_c, rhoR_c, Te_c)
    print(f"Hot spot: rho = {rho_c:.1f} g/cm^3, "
          f"rhoR = {rhoR_c:.2f} g/cm^2, Te = {Te_c:.1f} keV")
    print(f"  rho*lambda_alpha = {rho_lambda_alpha(rho_c, Te_c):.4f} g/cm^2")
    print(f"  tau              = {rhoR_c / rho_lambda_alpha(rho_c, Te_c):.3f}")
    print(f"  Q-hat fraction   = {float(f_c):.4f}")
    print()

    # Scan over rhoR at fixed rho, Te — broadcasting demo
    rhoR_scan = np.array([0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 1.00])
    f_scan = qhat(rho=100.0, rhoR=rhoR_scan, Te=10.0)
    print("rhoR scan at rho=100 g/cm^3, Te=10 keV:")
    print(f"  {'rhoR':>8s}  {'Q-hat':>8s}")
    for rr, ff in zip(rhoR_scan, f_scan):
        print(f"  {rr:8.3f}  {ff:8.4f}")
    print()

    # Temperature scan at fixed rhoR -- ignition-relevant regime
    Te_scan = np.array([3.0, 5.0, 7.0, 10.0, 15.0, 20.0])
    f_Te = qhat(rho=100.0, rhoR=0.30, Te=Te_scan)
    print("Te scan at rho=100 g/cm^3, rhoR=0.30 g/cm^2:")
    print(f"  {'Te':>8s}  {'Q-hat':>8s}")
    for tt, ff in zip(Te_scan, f_Te):
        print(f"  {tt:8.2f}  {ff:8.4f}")
