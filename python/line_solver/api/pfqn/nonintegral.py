"""Nonintegral degrees of multiprogramming (Dowdy and Gordon 1984).

L. W. Dowdy, K. D. Gordon, "Algorithms for Nonintegral Degrees of
Multiprogramming in Closed Queuing Networks", Performance Evaluation
4(1):19-28, 1984.
"""

from typing import Tuple

import numpy as np
from scipy.special import gammaln

__all__ = ['pfqn_dnc', 'pfqn_nintmva']


def _dnc_eval(n: float, A: np.ndarray, u: np.ndarray, node: np.ndarray,
              j: np.ndarray) -> float:
    """Partial-fraction series evaluated at a real population n.

    Analytic for n > -1; below that the binomial factor changes sign and the
    log-domain evaluation would lose it, so it is not extended there.
    """
    if n <= -1:
        return float('nan')
    return float(np.sum(A * np.exp(gammaln(n + j) - gammaln(j)
                                   - gammaln(n + 1) + n * np.log(u[node]))))


def pfqn_dnc(L, N) -> Tuple[float, float, float]:
    """Normalizing constant and throughput at a REAL-VALUED population, by
    partial-fraction inversion of the network generating function.

    With distinct loads x_1..x_G of multiplicities m_1..m_G the generating
    function prod_g (1-x_g u)^{-m_g} expands as::

        G(n) = sum_g sum_{j=1..m_g} A_gj C(n+j-1,j-1) x_g^n,

    every term of which is analytic in n, so evaluating at a real n
    interpolates the integral normalizing constants exactly and gives a smooth
    throughput curve X(N) = G(N-1)/G(N) through the integral points. For
    all-distinct loads A_g = prod_{l!=g} x_g/(x_g - x_l) is used directly; with
    repeated loads the coefficients are recovered from G(0..M-1).

    Only the queueing part admits this continuation: the delay sequence Z^n/n!
    is entire and has no partial-fraction expansion, so a think time is not
    accepted here. Use pfqn_nintmva for nonintegral populations with a delay.

    Args:
        L: Service demand vector (M,) of the queueing stations.
        N: Population (real nonnegative scalar; may be fractional).

    Returns:
        Tuple (X, G, lG).
    """
    L = np.asarray(L, dtype=float).ravel()
    N = np.asarray(N, dtype=float).ravel()
    if N.size > 1:
        raise ValueError('pfqn_dnc is a single-class method, but the population '
                         'vector has more than one entry.')
    n = float(N[0])
    L = L[L > 0]
    if L.size == 0:
        raise ValueError('pfqn_dnc requires at least one station with positive demand.')
    if n < 0:
        raise ValueError('pfqn_dnc requires a nonnegative population.')

    M = L.size
    xmax = L.max()
    y = L / xmax

    # Distinct loads merged under a relative tolerance, so numerically
    # coincident loads go to the multiplicity branch rather than to a
    # near-singular partial-fraction denominator.
    ys = np.sort(y)[::-1]
    u = [ys[0]]
    mult = [1]
    for i in range(1, M):
        if ys[i] > u[-1] * (1 - 1e-9):
            mult[-1] += 1
        else:
            u.append(ys[i])
            mult.append(1)
    u = np.array(u)
    mult = np.array(mult)
    Gd = u.size

    if np.all(mult == 1):
        A = np.ones(Gd)
        for g in range(Gd):
            idx = np.arange(Gd) != g
            A[g] = np.prod(u[g] / (u[g] - u[idx]))
        j = np.ones(Gd)
        node = np.arange(Gd)
    else:
        gint = np.array([1.0])
        for i in range(M):
            gint = np.convolve(gint, y[i] ** np.arange(M))[:M]
        node = np.zeros(M, dtype=int)
        j = np.zeros(M)
        c = 0
        for g in range(Gd):
            for jj in range(1, mult[g] + 1):
                node[c] = g
                j[c] = jj
                c += 1
        F = np.zeros((M, M))
        for k in range(M):
            F[k, :] = np.exp(gammaln(k + j) - gammaln(j) - gammaln(k + 1)
                             + k * np.log(u[node]))
        A = np.linalg.solve(F, gint)

    GN = _dnc_eval(n, A, u, node, j)
    GN1 = _dnc_eval(n - 1, A, u, node, j)

    lG = np.log(GN) + n * np.log(xmax)
    G = np.exp(lG)
    if n <= 0 or GN <= 0 or np.isnan(GN1):
        X = float('nan')
    else:
        X = (GN1 / GN) / xmax
    return float(X), float(G), float(lG)


def pfqn_nintmva(L, N, Z: float = 0.0):
    """Exact MVA recursion started from the FRACTIONAL base n0 = N - floor(N),
    giving mean performance measures at a real-valued population ("aMVA").

    The recursion is the standard Reiser-Lavenberg one, stepped in unit
    increments from n = n0 (where the arrival-theorem term is taken as 0, the
    network below the base being empty) up to n = N. At integer N the base is 0
    and the recursion is bit-identical to exact MVA; at fractional N it
    interpolates smoothly through the integral points.

    Unlike pfqn_dnc this accepts a think time. It is single-class: for
    fractional multiclass populations use pfqn_bs, which accepts them directly.

    Args:
        L: Service demand vector (M,) of the queueing stations.
        N: Population (real nonnegative scalar; may be fractional).
        Z: Think time (scalar, default 0).

    Returns:
        Tuple (X, Q, U, R).
    """
    L = np.asarray(L, dtype=float).ravel()
    N = np.asarray(N, dtype=float).ravel()
    if N.size > 1:
        raise ValueError('pfqn_nintmva is a single-class method, but the '
                         'population vector has more than one entry. Use '
                         'pfqn_bs for fractional multiclass populations.')
    n_target = float(N[0])
    Z = float(np.sum(np.asarray(Z, dtype=float)))
    if n_target < 0:
        raise ValueError('pfqn_nintmva requires a nonnegative population.')

    M = L.size
    Q = np.zeros(M)
    R = np.zeros(M)
    X = 0.0
    if n_target == 0:
        return 0.0, Q, np.zeros(M), R

    n = n_target - np.floor(n_target)
    if n == 0:
        n = 1.0
    while n <= n_target + 1e-12:
        R = L * (1 + Q)
        X = n / (Z + R.sum())
        Q = X * R
        n += 1
    return float(X), Q, X * L, R
