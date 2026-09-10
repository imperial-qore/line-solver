"""
Mean busy period of order n for a subnetwork of a product-form queueing network.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_busyp.m``.

References:
    H. Daduna, "Busy Periods for Subnetworks in Stochastic Networks: Mean Value
    Analysis", Journal of the ACM 35(3):668-674, 1988 (Theorem 1 for closed
    Gordon-Newell networks and Theorem 3 for open Jackson networks).
"""

import numpy as np


def _lse(v):
    """log-sum-exp, stable when every entry is -inf."""
    v = np.asarray(v, dtype=float).ravel()
    if v.size == 0:
        return -np.inf
    m = v.max()
    if not np.isfinite(m):
        return m
    return m + np.log(np.exp(v - m).sum())


def _rates(mu, idx, K):
    """Rates of the selected nodes for populations 1..K.

    A rate table shorter than K keeps its last rate, the saturated-server
    convention; a node whose rate does not saturate (an infinite server) must be
    supplied as a callable ``mu(j, kvec)`` instead.
    """
    idx = np.asarray(idx, dtype=int)
    if callable(mu):
        out = np.zeros((idx.size, K))
        kvec = np.arange(1, K + 1)
        for t, j in enumerate(idx):
            out[t, :] = np.asarray(mu(int(j), kvec), dtype=float)
        return out
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    rows = mu[idx, :]
    if rows.shape[1] < K:
        pad = np.repeat(rows[:, -1:], K - rows.shape[1], axis=1)
        rows = np.hstack([rows, pad])
    return rows[:, :K]


def _lgvec(alpha, mu, K):
    """Log normalizing constants of orders 0..K of a set of nodes.

    ``lg[m]`` is the log of the sum over the compositions n_1+...+n_L = m of the
    product over the nodes of prod_{k=1}^{n_i} alpha_i/mu_i(k), the G(m,I) and
    H(m,I) of the paper. The nodes are convolved one at a time in the log
    domain, which avoids the overflow the paper handles with the ratio
    recursions of Corollaries 2 and 4.
    """
    alpha = np.asarray(alpha, dtype=float).ravel()
    lg = np.full(K + 1, -np.inf)
    lg[0] = 0.0
    with np.errstate(divide='ignore'):
        lalpha = np.log(alpha)
        lmu = np.log(np.asarray(mu, dtype=float))
    for i in range(alpha.size):
        li = np.concatenate([[0.0], np.cumsum(lalpha[i] - lmu[i, :K])])
        lgnew = np.full(K + 1, -np.inf)
        for m in range(K + 1):
            lgnew[m] = _lse(lg[m::-1] + li[:m + 1])
        lg = lgnew
    return lg


def _trunc(alpha, mu, subnet, nmax, tol):
    """Truncation order of the open-network sum sum_{m>=n} G(m,I)."""
    K = max(nmax + 8, 16)
    rows = _rates(mu, subnet, K)
    rho = np.asarray(alpha, dtype=float).ravel() / rows[:, -1]
    if rho.max() >= 1:
        raise ValueError('The subnetwork is not stable, its busy period is infinite.')
    while True:
        lg = _lgvec(alpha, _rates(mu, subnet, K), K)
        # decay rate read off the last two orders, the exact ratio for a
        # saturated single-server subnetwork and an upper estimate otherwise
        r = np.exp(lg[K] - lg[K - 1])
        if not r < 1:
            r = rho.max()
        ltail = lg[K] + np.log(r) - np.log1p(-r)
        if ltail - _lse(lg[nmax:K + 1]) < np.log(tol):
            return K, lg
        K = 2 * K
        if K > 1e6:
            raise ValueError('The open busy period sum did not converge, the subnetwork is nearly saturated.')


def pfqn_busyp(alpha, mu, P, N, subnet, n, gamma=None, tol=1e-12):
    """Mean busy period of order ``n`` for the subnetwork ``subnet``.

    The busy period of order n is the time from the instant a job entering the
    subnetwork finds n-1 jobs in it up to the next instant when fewer than n
    jobs remain in it.

    Parameters
    ----------
    alpha : array (J,)
        Relative arrival rates, the solution of x*P = x for a closed network and
        of x = gamma + x*P for an open one.
    mu : array (J, K) or callable
        Load-dependent service rates, ``mu[j, k-1]`` with k jobs at node j, or a
        callable ``mu(j, kvec)`` when the rates do not saturate.
    P : array (J, J)
        Routing matrix.
    N : int or float
        Population; ``numpy.inf`` for an open network.
    subnet : sequence of int
        Zero-based indexes of the nodes forming the subnetwork.
    n : int or sequence of int
        Busy period order(s), 1 <= n <= N.
    gamma : array (J,), optional
        External arrival rates; required for an open network.
    tol : float
        Relative tolerance of the open-network tail truncation.

    Returns
    -------
    b : np.ndarray
        Mean busy period duration(s), same shape as ``n``.
    lG : np.ndarray
        Log normalizing constants of the subnetwork.
    lH : np.ndarray
        Log normalizing constants of the complement, empty for an open network.
    """
    alpha = np.asarray(alpha, dtype=float).ravel()
    J = alpha.size
    P = np.asarray(P, dtype=float)
    nvec = np.atleast_1d(np.asarray(n, dtype=int))
    if N is None:
        N = np.inf
    is_closed = np.isfinite(N)

    subnet = np.unique(np.asarray(subnet, dtype=int))
    if subnet.size == 0:
        raise ValueError('The subnetwork must be non-empty.')
    if is_closed and subnet.size >= J:
        # a closed network needs jobs outside the subnetwork to start a busy period
        raise ValueError('In a closed network the subnetwork must be a proper subset of the nodes.')
    if subnet.min() < 0 or subnet.max() >= J:
        raise ValueError('The subnetwork indexes are out of range.')
    compl = np.setdiff1d(np.arange(J), subnet)

    if is_closed:
        N = int(N)
        if N < 1:
            raise ValueError('The population of a closed network must be a positive integer.')
        if nvec.min() < 1 or nvec.max() > N:
            raise ValueError('The busy period order must be an integer in 1..N.')
    else:
        if gamma is None:
            raise ValueError('An open network requires the external arrival rates gamma.')
        if nvec.min() < 1:
            raise ValueError('The busy period order must be a positive integer.')

    # A(I) for a closed network, C(I) for an open one: both are the total rate
    # at which jobs enter the subnetwork from outside it, which is what starts a
    # busy period. The closed network has no external stream.
    inflow = float(alpha[compl] @ P[np.ix_(compl, subnet)] @ np.ones(subnet.size))
    if gamma is not None:
        gamma = np.asarray(gamma, dtype=float).ravel()
        inflow += float(gamma[subnet].sum())
    if inflow <= 0:
        raise ValueError('No job ever enters the subnetwork, its busy period is undefined.')

    b = np.zeros(nvec.shape)
    if is_closed:
        lG = _lgvec(alpha[subnet], _rates(mu, subnet, N), N)
        lH = _lgvec(alpha[compl], _rates(mu, compl, N), N)
        for t, nt in enumerate(nvec):
            # Theorem 1: sum_{m=n}^{N} G(m,I) H(N-m,I) over G(n-1,I) H(N-n,I) A(I)
            num = _lse(lG[nt:N + 1] + lH[N - nt::-1])
            b[t] = np.exp(num - lG[nt - 1] - lH[N - nt] - np.log(inflow))
    else:
        K, lG = _trunc(alpha[subnet], mu, subnet, int(nvec.max()), tol)
        lH = np.zeros(0)
        for t, nt in enumerate(nvec):
            # Theorem 3: sum_{m=n}^{Inf} G(m,I) over G(n-1,I) C(I)
            num = _lse(lG[nt:K + 1])
            b[t] = np.exp(num - lG[nt - 1] - np.log(inflow))
    if np.isscalar(n) or np.asarray(n).ndim == 0:
        return float(b[0]), lG, lH
    return b, lG, lH
