"""
Reduction Heuristic (RD) method for load-dependent networks.

This module implements the Reduction Heuristic method for computing
normalizing constants in load-dependent queueing networks.

References:
    Casale, G., "Approximate Mean Value Analysis for Load-Dependent
    Queueing Networks", IEEE TPDS, 2017.
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass

from .mva import pfqn_mva
from .ncld import pfqn_gldsingle
from .nc import pfqn_nc


@dataclass
class RdOptions:
    """Options for RD algorithm."""
    tol: float = 1e-6
    method: str = 'default'


@dataclass
class RdResult:
    """Result from Reduced Decomposition algorithm."""
    lGN: float        # Logarithm of normalizing constant
    Cgamma: float     # Gamma correction factor


def pfqn_rd(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    mu: Optional[np.ndarray] = None,
    options: Optional[RdOptions] = None
) -> Tuple[float, float]:
    """
    Reduction Heuristic (RD) method for load-dependent networks.

    Computes the logarithm of the normalizing constant using the reduction
    heuristic method, which handles load-dependent service rates.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (1 x R)
        Z: Think time vector (1 x R)
        mu: Load-dependent rate matrix (M x sum(N)). If None, assumes
            load-independent rates (all 1.0).
        options: Solver options

    Returns:
        Tuple of (lGN, Cgamma) where:
            lGN: Logarithm of normalizing constant
            Cgamma: Gamma correction factor
    """
    if options is None:
        options = RdOptions()

    L = np.atleast_2d(L).astype(float)
    N = np.atleast_1d(N).astype(float)
    Z = np.atleast_1d(Z).astype(float)

    M = L.shape[0]
    total_pop = int(np.sum(N))

    # If mu not provided, create load-independent rates (all 1.0)
    if mu is None:
        mu = np.ones((M, max(total_pop, 1)))
    else:
        mu = np.atleast_2d(mu).astype(float)

    M, R = L.shape
    lambda_vec = np.zeros(R)

    if total_pop < 0:
        return -np.inf, 0.0

    if total_pop == 0:
        return 0.0, 1.0

    # Normalize L by mu for LI (load-independent) stations
    for ist in range(M):
        if np.all(mu[ist, :min(total_pop, mu.shape[1])] == mu[ist, 0]):
            # LI station
            L[ist, :] = L[ist, :] / mu[ist, 0]
            mu[ist, :] = 1

    # Initialize gamma and related arrays
    gamma = np.ones((M, total_pop))
    mu_work = mu[:, :total_pop].copy()
    mu_work[np.isnan(mu_work)] = np.inf

    s = np.zeros(M, dtype=int)
    for ist in range(M):
        if np.isfinite(mu_work[ist, -1]):
            diffs = np.abs(mu_work[ist, :] - mu_work[ist, -1])
            matches = np.where(diffs < options.tol)[0]
            if len(matches) > 0:
                s[ist] = matches[0] + 1  # 1-based index
            else:
                s[ist] = total_pop
        else:
            s[ist] = total_pop

    y = L.copy()

    for ist in range(M):
        if np.isinf(mu_work[ist, s[ist] - 1]):
            finite_indices = np.where(np.isfinite(mu_work[ist, :]))[0]
            if len(finite_indices) > 0:
                s[ist] = finite_indices[-1] + 1
        y[ist, :] = y[ist, :] / mu_work[ist, s[ist] - 1]

    for ist in range(M):
        gamma[ist, :] = mu_work[ist, :] / mu_work[ist, s[ist] - 1]

    # Compute beta - match MATLAB's computation exactly
    # MATLAB computes beta directly without special-casing gamma=1
    beta = np.ones((M, total_pop))
    for ist in range(M):
        with np.errstate(divide='ignore', invalid='ignore'):
            beta[ist, 0] = gamma[ist, 0] / (1 - gamma[ist, 0])

            for j in range(1, total_pop):
                beta[ist, j] = (1 - gamma[ist, j - 1]) * (gamma[ist, j] / (1 - gamma[ist, j]))

    beta[np.isnan(beta)] = np.inf

    if np.all(beta == np.inf):
        # Fall back to standard NC
        _, lGN = pfqn_nc(L, N, Z, method='default')
        return lGN, 1.0

    # Compute correction factor
    Cgamma = 0.0
    sld = s[s > 1]
    vmax = min(int(np.sum(sld - 1)), total_pop)

    # Debug: compare with MATLAB/Kotlin
    import os
    _debug = os.environ.get('DEBUG_RD')

    # Compute MVA for y values
    mva_result = pfqn_mva(y, N, np.zeros_like(N))
    if hasattr(mva_result, 'X'):
        Y = mva_result.X
    else:
        Y = mva_result[0] if isinstance(mva_result, tuple) else np.zeros(R)

    # rhoN = y * Y' in MATLAB: per-station utilization vector (M x 1)
    Y_flat = Y.flatten()
    rhoN = y @ Y_flat  # (M, R) @ (R,) = (M,)

    # Compute lEN values
    lEN = np.zeros(vmax + 2)
    lEN[0] = 0  # ln(1) = 0

    for vtot in range(1, vmax + 1):
        # Reshape rhoN to column (M, 1) so pfqn_gldsingle sees M stations, 1 class
        lEN[vtot] = np.real(pfqn_gldsingle(rhoN.reshape(-1, 1), vtot, beta).lG)

    # Compute Cgamma
    for vtot in range(vmax + 1):
        EN = np.exp(lEN[vtot])
        Cgamma += ((total_pop - max(0, vtot - 1)) / total_pop) * EN

    # Compute final normalizing constant
    _, lGN = pfqn_nc(y, N, Z, method='default')

    if _debug:
        print(f'DEBUG RD: M={M}, R={R}, total_pop={total_pop}')
        print(f'DEBUG RD: L=\n{L}')
        print(f'DEBUG RD: Z={Z}')
        print(f'DEBUG RD: mu=\n{mu}')
        print(f'DEBUG RD: s={s}')
        print(f'DEBUG RD: y=\n{y}')
        print(f'DEBUG RD: gamma=\n{gamma}')
        print(f'DEBUG RD: beta=\n{beta}')
        print(f'DEBUG RD: vmax={vmax}')
        print(f'DEBUG RD: Y (throughputs)={Y_flat}')
        print(f'DEBUG RD: rhoN={rhoN}')
        print(f'DEBUG RD: lEN={lEN[:vmax+1]}')
        print(f'DEBUG RD: lGN_nc={lGN}')
        print(f'DEBUG RD: Cgamma={Cgamma}')

    lGN = lGN + np.log(Cgamma) if Cgamma > 0 else lGN

    return lGN, Cgamma
