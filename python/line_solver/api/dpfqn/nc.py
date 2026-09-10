"""
Normalizing constants of discrete-time closed cycles of Bernoulli servers.

Native Python port of matlab/src/api/dpfqn/dpfqn_nc.m and dpfqn_ncld.m.

Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
Springer, 2001, theorem 3.2, corollary 3.4, propositions 3.5 and 3.17-3.19,
corollary 3.20.
"""

from typing import Tuple

import numpy as np

__all__ = ['dpfqn_nc', 'dpfqn_ncld']

_RESCALE = 1e250


def dpfqn_nc(p, N: int) -> Tuple[float, float, np.ndarray]:
    """Buzen-style recursions for the discrete-time closed cycle.

    Args:
        p: per-slot service completion probabilities p_j in (0,1), one per node
        N: number of customers cycling in the J nodes

    Returns:
        (lG, G, G1) with lG the log of the time-stationary normalizing constant
        G(N,J), G that constant, and G1 the vector of arrival constants with
        ``G1[k] = G_1(k,J)`` for k = 0..N.

    With q_j = 1 - p_j the queue length vector has the product form of
    corollary 3.4::

        pi(n_1,...,n_J) = prod_j (q_j/p_j)^n_j (1/q_j)^{1{n_j>0}} / G(N,J)

    The extra factor on the busy nodes is what separates it from the
    continuous-time Gordon-Newell form: a homogeneous cycle is uniform on the
    state space in continuous time and is not here. G obeys proposition 3.18::

        G(k,j) = G(k,j-1) + (q_j/p_j) G(k-1,j) + G(k-1,j-1)

    with G(0,j) = 1 and G(k,0) = 0 for k >= 1. Unlike the continuous-time
    convolution algorithm this recursion is not invariant to the numbering of
    the nodes. The arrival constants obey proposition 3.19::

        G1(k,J) = G(k-1,J-1) + (q_J/p_J) G1(k-1,J),   k >= 3

    with G1(1,J) = 1 and G1(2,J) = q_1/p_1 + sum_{j>=2} 1/p_j. By lemma 7.3 the
    arrival constant is the same at every node, so one family suffices, and

        throughput per slot   X   = G1(N,J) / G(N,J)   (equal at every node)
        utilization           U_j = X / p_j
        tail probability      P(X_j >= k)
                              = (q_j/p_j)^k (1/q_j) G1(N-k+1,J) / G(N,J)

    The last identity is corollary 3.20(a) with its index corrected: as printed
    there the right-hand side evaluates to P(X_j >= k+1). Corollary 3.20(c),
    which transfers the tail from node 1 to node j, holds for k >= 1 only; at
    k = 0 both tails are 1 while the stated ratio is q_1/q_j.

    G and G1 come back in a common scale, unity unless the recursion had to be
    rescaled to keep (q/p)^N inside double range. Ratios such as G1[k]/G are
    exact either way, and lG is always the true log constant.
    """
    p = np.asarray(p, dtype=float).ravel()
    if p.size == 0 or not np.all(np.isfinite(p)) or np.any(p <= 0) or np.any(p >= 1):
        raise ValueError('service probabilities must be real and in the open interval (0,1)')
    if int(N) != N or N < 0:
        raise ValueError('N must be a non-negative integer')
    N = int(N)
    J = p.size
    x = (1.0 - p) / p

    # Proposition 3.18. The recursion is linear and homogeneous in the whole
    # table, so rescaling every entry at once preserves it; G1 is rescaled in
    # step so that the two families stay in one common scale.
    Gt = np.zeros((N + 1, J + 1))
    Gt[0, :] = 1.0
    G1 = np.zeros(N + 1)
    lscale = 0.0
    for k in range(1, N + 1):
        for j in range(1, J + 1):
            Gt[k, j] = Gt[k, j - 1] + x[j - 1] * Gt[k - 1, j] + Gt[k - 1, j - 1]
        if k == 1:
            G1[1] = Gt[0, 0]
        elif k == 2:
            G1[2] = (x[0] + np.sum(1.0 / p[1:])) * Gt[0, 0]
        else:
            G1[k] = Gt[k - 1, J - 1] + x[J - 1] * G1[k - 1]
        mx = Gt[k, :].max()
        if mx > _RESCALE:
            Gt /= mx
            G1 /= mx
            lscale += float(np.log(mx))
    G = float(Gt[N, J])
    lG = float(np.log(G)) + lscale
    return lG, G, G1


def dpfqn_ncld(P, N: int) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Convolution over the discrete-time product form with load dependence.

    Args:
        P: (J, N) array of service probabilities, ``P[j, n-1] = p_j(n)``
        N: number of customers cycling in the J nodes

    Returns:
        (lG, G, W, Gc, Wa) with lG the log normalizing constant, G the vector
        of constants G(k) for k = 0..N, W the (J, N+1) time-stationary node
        weights, Gc the (J, N+1) complement constants of the cycle with node j
        removed, and Wa the (J, N+1) arrival node weights. All four tables are
        in one common scale.

    Theorem 3.2 keeps the product form when the service probability of node j
    depends on its own queue length::

        pi(n_1,...,n_J) = prod_j w_j(n_j) / G(N,J)
        w_j(n) = prod_{h=1}^{n-1} q_j(h) / prod_{h=1}^{n} p_j(h),  w_j(0) = 1

    Note where the state dependence sits: the missing q_j in the numerator is
    tied to the actual queue length, not to the node being non-empty, so the
    tidy (1/q_j)^{1{n>0}} of the state independent case cannot be factored out.
    The marginal law follows from the complement constants::

        P(X_j = n) = W[j, n] * Gc[j, N-n] / G[N]

    and is exact whatever the common scale. Deconvolution is never used, so a
    node with a near-unit service probability does not spoil the accuracy of
    the other marginals.

    Wa holds the arrival weights of proposition 3.5, v_j(n) = w_j(n) q_j(n),
    the law seen by a customer joining node j with himself not counted::

        Parr(X_j = n) = Wa[j, n] * Gc[j, N-1-n] / sum_m Wa[j, m] Gc[j, N-1-m]

    on n = 0..N-1. It is not the time-stationary law, which is the
    discrete-time departure from the continuous-time arrival theorem.
    """
    if int(N) != N or N < 0:
        raise ValueError('N must be a non-negative integer')
    N = int(N)
    P = np.asarray(P, dtype=float)
    if P.ndim != 2 or P.size == 0:
        raise ValueError('P must be a non-empty (J, N) matrix of service probabilities')
    J, cols = P.shape
    if cols < N:
        raise ValueError('P must supply p_j(n) for every n = 1..%d, but has %d columns' % (N, cols))
    if N == 0:
        ones = np.ones((J, 1))
        return 0.0, np.ones(1), ones, ones, ones
    Pn = P[:, :N]
    if not np.all(np.isfinite(Pn)) or np.any(Pn <= 0) or np.any(Pn > 1):
        raise ValueError('service probabilities must be real and in the interval (0,1]')
    if N >= 2 and np.any(Pn[:, :N - 1] >= 1):
        # p_j(n)=1 is admissible only at the last reachable population, where
        # the missing q_j(n) never multiplies any weight.
        raise ValueError('service probabilities below the population bound must be strictly less than 1')

    # Per-node log weights of theorem 3.2, then a per-node shift so that every
    # weight vector peaks at one. The shifts cancel in every ratio below.
    lw = np.zeros((J, N + 1))
    for j in range(J):
        cq = 0.0
        cp = 0.0
        for n in range(1, N + 1):
            cp += float(np.log(Pn[j, n - 1]))
            lw[j, n] = cq - cp
            cq += float(np.log(1.0 - Pn[j, n - 1])) if Pn[j, n - 1] < 1.0 else -np.inf
    shift = lw.max(axis=1)
    W = np.exp(lw - shift[:, None])

    Wa = W.copy()
    for j in range(J):
        Wa[j, 1:] = W[j, 1:] * (1.0 - Pn[j, :N])

    # Prefix/suffix convolution: pre[j] covers nodes 0..j-1 and suf[j] covers
    # nodes j..J-1, so the complement of node j is pre[j] * suf[j+1].
    pre = np.zeros((J + 1, N + 1))
    pre[0, 0] = 1.0
    for j in range(J):
        pre[j + 1, :] = _dt_conv(pre[j, :], W[j, :], N)
    suf = np.zeros((J + 1, N + 1))
    suf[J, 0] = 1.0
    for j in range(J - 1, -1, -1):
        suf[j, :] = _dt_conv(suf[j + 1, :], W[j, :], N)
    G = pre[J, :]
    Gc = np.zeros((J, N + 1))
    for j in range(J):
        Gc[j, :] = _dt_conv(pre[j, :], suf[j + 1, :], N)

    lG = float(np.log(G[N])) + float(shift.sum())
    return lG, G, W, Gc, Wa


def _dt_conv(a: np.ndarray, b: np.ndarray, N: int) -> np.ndarray:
    """Convolution of two population tables truncated at N."""
    c = np.zeros(N + 1)
    for k in range(N + 1):
        c[k] = float(np.dot(a[:k + 1], b[k::-1]))
    return c
