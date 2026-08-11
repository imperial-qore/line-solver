"""
Limited Joint Dependence (LJD) lattice indexing helpers for population-vector lookup tables.

These are lattice-indexing helpers for tabulated dependence functions. The
tables themselves are no longer a model field: a tabulated dependence is
expressed as a class-dependence handle beta_{i,r}(n) that indexes the table
internally (see fes_beta_handle).

Key functions:
    ljd_linearize: Convert population vector to linear index
    fes_beta_handle: Wrap a per-class throughput table as beta_{i,r}(n)

References:
    Original MATLAB: matlab/src/api/pfqn/ljd_linearize.m
"""

import numpy as np
from typing import Union


def ljd_linearize(nvec: np.ndarray, cutoffs: np.ndarray) -> int:
    """
    Convert per-class population vector to linearized index.

    Maps a multi-dimensional population vector to a single linear index
    for efficient lookups in tabulated scaling tables.

    Index formula: idx = 1 + n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...

    Args:
        nvec: Per-class populations [n1, n2, ..., nK]
        cutoffs: Per-class cutoffs [N1, N2, ..., NK]

    Returns:
        1-based linearized index

    References:
        Original MATLAB: matlab/src/api/pfqn/ljd_linearize.m
    """
    nvec = np.asarray(nvec, dtype=int).flatten()
    cutoffs = np.asarray(cutoffs, dtype=int).flatten()

    K = len(nvec)
    idx = 1  # 1-indexed for MATLAB compatibility
    multiplier = 1

    for k in range(K):
        nk = min(nvec[k], cutoffs[k])  # Clamp to cutoff
        idx += nk * multiplier
        multiplier *= (cutoffs[k] + 1)

    return idx


def infradius_h(x: np.ndarray, L: np.ndarray, N: np.ndarray,
                alpha: np.ndarray) -> np.ndarray:
    """
    Helper function for infinite radius computation with logistic transformation.

    Used in normalizing constant computation via integration methods.

    Args:
        x: Logistic transformation parameters
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        alpha: Load-dependent rate matrix

    Returns:
        Evaluated function value for integration

    References:
        Original MATLAB: matlab/src/api/pfqn/infradius_h.m
    """
    # Import here to avoid circular dependency
    from .ncld import pfqn_gld

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    x = np.atleast_2d(np.asarray(x, dtype=float))

    M = L.shape[0]
    Nt = int(np.sum(N))
    beta = N / Nt if Nt > 0 else np.zeros_like(N)

    y = np.zeros(x.shape[0])

    for i in range(x.shape[0]):
        xi = x[i, :]

        # Logistic transformation
        t = np.exp(xi) / (1 + np.exp(xi))
        tb = np.sum(beta * t)

        # Evaluate h function
        z = np.sum(L * np.tile(np.exp(2 * np.pi * 1j * (t - tb)), (M, 1)), axis=1)
        # Reshape z to column vector (M x 1) to match MATLAB convention
        # (M stations, 1 class) - np.atleast_2d would make it (1, M) which is wrong
        gld_result = pfqn_gld(z.reshape(-1, 1), Nt, alpha)
        # Extract scalar G value from PfqnNcResult
        gld_value = gld_result.G if hasattr(gld_result, 'G') else gld_result

        # Jacobian of transformation
        jacobian = np.prod(np.exp(xi) / (1 + np.exp(xi)) ** 2)

        y[i] = np.real(gld_value * jacobian)

    return y


def infradius_hnorm(x: np.ndarray, L: np.ndarray, N: np.ndarray,
                    alpha: np.ndarray) -> np.ndarray:
    """
    Helper function for infinite radius computation with normal CDF (probit) transformation.

    Uses normcdf/normpdf transformation instead of logistic (used in infradius_h).

    Args:
        x: Normal CDF transformation parameters
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        alpha: Load-dependent rate matrix

    Returns:
        Evaluated function value for integration

    References:
        Original MATLAB: matlab/src/api/pfqn/infradius_hnorm.m
    """
    from scipy.stats import norm
    from .ncld import pfqn_gld

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    x = np.atleast_2d(np.asarray(x, dtype=float))

    M = L.shape[0]
    Nt = int(np.sum(N))
    beta = N / Nt if Nt > 0 else np.zeros_like(N)

    y = np.zeros(x.shape[0])

    for i in range(x.shape[0]):
        xi = x[i, :]

        # Probit (normal CDF) transformation
        t = norm.cdf(xi)
        tb = np.sum(beta * t)

        # Evaluate h function
        z = np.sum(L * np.tile(np.exp(2 * np.pi * 1j * (t - tb)), (M, 1)), axis=1)
        # Reshape z to column vector (M x 1) to match MATLAB convention
        gld_result = pfqn_gld(z.reshape(-1, 1), Nt, alpha)
        gld_value = gld_result.G if hasattr(gld_result, 'G') else gld_result

        # Jacobian of normal CDF transformation
        jacobian = np.prod(norm.pdf(xi))

        y[i] = np.real(gld_value * jacobian)

    return y



def fes_beta_handle(scaling_table, cutoffs):
    """
    Wrap a flow-equivalent-server (FES) throughput table as a per-class
    class-dependence function beta_{i,r}(n).

    scaling_table is a sequence of length K in which scaling_table[r] is the
    linearized vector of class-r throughputs X_r(n) of the aggregated
    subnetwork, indexed by ljd_linearize(min(n, cutoffs), cutoffs). cutoffs is
    the per-class population vector the table was tabulated on.

    The returned callable takes the per-class population vector n at the station
    and returns the length-K array [X_1(n), ..., X_K(n)], i.e. Sauer's
    chain-dependent service rates mu_{r,i}(n) (Sauer 1983, "Computational
    Algorithms for State-Dependent Queueing Networks", eq. (40)). The population
    is clamped to cutoffs, so the rate saturates beyond the tabulated range
    exactly as the underlying table intends.

    This is the single class-dependence mechanism used across the solvers: the
    exact convolution (pfqn_conv) reads mu_{r,i}(n) from it, and AMVA-QD reads
    the same handle through pfqn_cdfun.

    References:
        Original MATLAB: matlab/src/api/fes/fes_beta_handle.m
    """
    tables = [None if t is None else np.asarray(t, dtype=float).flatten()
              for t in scaling_table]
    cut = np.asarray(cutoffs, dtype=int).flatten()

    def _beta(n):
        K = len(tables)
        v = np.ones(K)
        nn = np.round(np.asarray(n, dtype=float).flatten()).astype(int)
        if nn.size < cut.size:
            nn = np.concatenate([nn, np.zeros(cut.size - nn.size, dtype=int)])
        elif nn.size > cut.size:
            nn = nn[:cut.size]
        n_clamped = np.maximum(0, np.minimum(nn, cut))
        idx = ljd_linearize(n_clamped, cut)
        tot = float(np.sum(nn))
        for r in range(K):
            tbl = tables[r]
            if tbl is not None and tbl.size > 0 and 1 <= idx <= tbl.size and nn[r] > 0:
                # see _kb/03-api-layer.md for rationale
                v[r] = tbl[idx - 1] * tot / float(nn[r])  # ljd_linearize is 1-based
        return v

    return _beta


__all__ = [
    'ljd_linearize',
    'fes_beta_handle',
    'infradius_h',
    'infradius_hnorm',
]
