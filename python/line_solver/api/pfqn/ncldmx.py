"""
Normalizing constant for mixed open/closed networks with limited load dependence.

The closed-conditional normalizing constant of a mixed limited load-dependent (LLD)
network equals a purely closed load-dependent normalizing constant in which every
queueing station i carries the Bruell-Balbo-Afshari effective capacity rate
mu_i^eff(n) = 1/EC_i(n), where EC is returned by pfqn_ldmx_ec and folds the open
classes into the closed subnetwork. The open classes contribute the separable
prefactor lGopen = sum_i log E_i(0), which reduces to -sum_i log(1-rho_i) in the
load-independent limit. Mean closed metrics follow from the standard normalizing
constant ratios, e.g. X_r = G(N-e_r)/G(N), matching the exact pfqn_mvaldmx solver.

References:
    MATLAB: matlab/src/api/pfqn/pfqn_ncldmx.m
"""

import numpy as np
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

from .mvaldmx import pfqn_ldmx_ec
from .ncld import pfqn_ncld


@dataclass
class PfqnNcldmxResult:
    """Result of a mixed limited load-dependent normalizing constant computation."""
    G: float
    lG: float
    lGopen: float
    method: str = "default"


def pfqn_ncldmx(lam: np.ndarray, D: np.ndarray, N: np.ndarray,
                Z: Optional[np.ndarray] = None,
                mu: Optional[np.ndarray] = None,
                S: Optional[np.ndarray] = None,
                options: Optional[Dict[str, Any]] = None) -> PfqnNcldmxResult:
    """
    Normalizing constant for mixed open/closed networks with limited load dependence.

    Args:
        lam: Arrival rate vector (R,) - 0 on closed classes
        D: Service demand matrix (M x R)
        N: Population vector (R,) - inf for open classes
        Z: Think time vector (R,), optional
        mu: Load-dependent rate matrix (M x >= sum(N_closed)), optional
        S: Number of servers per station (M,), kept for signature parity
        options: Solver options forwarded to pfqn_ncld

    Returns:
        PfqnNcldmxResult with the closed-conditional constant (G, lG), the open-class
        normalizing prefactor lGopen and the method used for the closed-conditional solve.
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    lam = np.asarray(lam, dtype=float).flatten()
    M, R = D.shape
    if Z is None:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float).flatten()

    openClasses = np.where(np.isinf(N))[0]
    closedClasses = np.array([r for r in range(R) if r not in openClasses], dtype=int)
    for r in closedClasses:
        if lam[r] != 0 and N[r] > 0:
            raise ValueError("pfqn_ncldmx: Arrival rate cannot be specified on closed classes.")

    Kc = int(np.sum(N[closedClasses])) if closedClasses.size > 0 else 0

    if mu is None:
        mu = np.ones((M, max(1, Kc)))
    mu = np.atleast_2d(np.asarray(mu, dtype=float))

    # pad mu to at least max(1,Kc) columns, then append one extra column as in
    # pfqn_mvaldmx so pfqn_ldmx_ec returns EC with Nt = max(1,Kc)+1 columns.
    min_cols = max(1, Kc)
    pad_cols = max(mu.shape[1], min_cols) + 1
    mup = np.empty((M, pad_cols))
    ncol = mu.shape[1]
    for j in range(pad_cols):
        src = min(j, ncol - 1)
        mup[:, j] = mu[:, src]

    EC, E, Eprime, _ = pfqn_ldmx_ec(lam, D, mup)
    lGopen = float(np.sum(np.log(E[:, 0])))

    if Kc == 0:
        return PfqnNcldmxResult(G=1.0, lG=0.0, lGopen=lGopen, method="exact")

    Dc = D[:, closedClasses]
    Nc = N[closedClasses]
    Zc = Z[closedClasses]
    muEff = 1.0 / EC[:, :Kc]

    cc = pfqn_ncld(Dc, Nc, Zc, muEff, options)
    return PfqnNcldmxResult(G=cc.G, lG=cc.lG, lGopen=lGopen, method=cc.method)
