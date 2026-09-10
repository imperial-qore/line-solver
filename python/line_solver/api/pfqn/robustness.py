"""Operational sensitivity of throughput to homogeneous-service-time violations.

R. Suri, "Robustness of Queuing Network Formulas", JACM 30(3):564-594, 1983
(eq. 3.11, Lemma 3.1, problem (P1)).
"""

from typing import Optional

import numpy as np

from .rgf import pfqn_rgf

__all__ = ['pfqn_hst']


def pfqn_hst(L, N, Z: float = 0.0, ist: Optional[int] = None) -> dict:
    """Robustness certificate for a single-class closed product-form solution:
    how far the predicted throughput can move when the homogeneous-service-time
    (HST) assumption fails at one station.

    HST states that the mean service time at station i does not depend on the
    queue length there. Suri perturbs it to S_i(n) = S_i (1 + a_n), one relative
    deviation per queue-length level, and shows (eq. 3.11) that to first order::

        [(1/X0) dX0/da_n] = c_n = P(n_i >= n+1)/u_i - P(n_i >= n),

    with u_i = L_i X0 and P(n_i >= n) = L_i^n G(N-n)/G(N). The naive certificate
    ``|dX0/X0| <= (sum_n |c_n|) d`` follows from ``|a_n| <= d`` alone, and by Lemma 3.1
    that total equals Q_i(N) - Q_i(N-1).

    That bound is loose because an operationally consistent perturbation must
    leave the observed mean service time unchanged, sum_n p_n a_n = 0. The
    constrained problem (P1)::

        max |sum_n c_n a_n|  s.t.  |a_n| <= d,  sum_n p_n a_n = 0

    is a one-constraint linear program, solved here exactly: its optimum sets
    a_n = +/-d according to whether c_n/p_n exceeds a threshold, with at most
    one fractional coordinate.

    Args:
        L: Service demand vector (M,) of the queueing stations.
        N: Population (nonnegative integer scalar).
        Z: Think time (scalar, default 0).
        ist: 0-based station index the perturbation applies to (default: the
            bottleneck, argmax L).

    Returns:
        Dict with keys station, X, U, Q, Pgeq, p, c, total, worst, astar.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = L.size
    N = np.asarray(N, dtype=float).ravel()
    if N.size > 1:
        raise ValueError('pfqn_hst is a single-class method.')
    n = float(N[0])
    if n < 1 or abs(n - round(n)) > 0:
        raise ValueError('pfqn_hst requires an integer population of at least one job.')
    n = int(round(n))
    if ist is None:
        ist = int(np.argmax(L))
    if ist < 0 or ist >= M:
        raise ValueError(f'Station index {ist} is out of range (the model has {M} '
                         f'queueing stations).')
    if L[ist] <= 0:
        raise ValueError(f'Station {ist} has zero demand, so its queue-length '
                         f'marginals are degenerate.')

    _, _, lg = pfqn_rgf(L, n, Z)          # lg[k] = log G(k)

    X = float(np.exp(lg[n - 1] - lg[n]))
    y = L[ist]
    k = np.arange(n + 1)
    Pgeq = np.exp(k * np.log(y) + lg[::-1] - lg[n])
    p = Pgeq - np.append(Pgeq[1:], 0.0)
    U = y * X
    Q = float(Pgeq[1:].sum())

    # eq. (3.11)
    c = np.zeros(n)
    for i in range(1, n + 1):
        Pn1 = Pgeq[i + 1] if i + 1 <= n else 0.0
        c[i - 1] = Pn1 / U - Pgeq[i]
    total = float(np.abs(c).sum())

    # (P1): one equality constraint plus a box. At the optimum
    # a_n = sign(c_n - lambda p_n), so sorting by c_n/p_n and sweeping the split
    # point enumerates every candidate lambda; the constraint fixes the single
    # fractional coordinate at the split.
    pp = p[1:]
    order = np.argsort(-(c / np.maximum(pp, np.finfo(float).tiny)))
    best = 0.0
    astar = np.zeros(n)
    for split in range(n + 1):
        a = -np.ones(n)
        a[order[:split]] = 1.0
        for piv in range(n):
            if pp[piv] <= 0:
                continue
            aa = a.copy()
            rest = float(pp @ aa) - pp[piv] * aa[piv]
            v = -rest / pp[piv]
            if v < -1 or v > 1:
                continue
            aa[piv] = v
            obj = float(c @ aa)
            if abs(obj) > abs(best):
                best = obj
                astar = aa
    if best < 0:
        best = -best
        astar = -astar                    # the feasible set is symmetric

    return {'station': ist, 'X': X, 'U': float(U), 'Q': Q, 'Pgeq': Pgeq,
            'p': p, 'c': c, 'total': total, 'worst': float(best),
            'astar': astar}
