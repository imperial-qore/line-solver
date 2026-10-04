"""
Effective arrival processes of a network under the Robust Queueing calculus.

Reference: C. Bandi, D. Bertsimas, N. Youssef (2015), "Robust Queueing Theory",
Operations Research 63(3), 676-700, Theorems 4-7 and 10.

MATLAB: matlab/src/api/npfqn/npfqn_traffic_rqt.m
"""

from typing import Tuple

import numpy as np

__all__ = ['npfqn_traffic_rqt']


def npfqn_traffic_rqt(lambda0, Gamma0, alpha0, F) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Effective arrival process perceived at each node of a single-class open
    queueing network under the Robust Queueing Theory calculus.

    The characterization composes three operators: passage through a queue with
    adversarial servers leaves the uncertainty set unchanged (robust Burke,
    Theorem 4), superposition merges sets by Theorem 5, and thinning by a
    fraction f scales the rate by f and the variability by f^(-1/alpha)
    (Theorem 6). The resulting equations are::

        lambda_j = lambda0_j + sum_i lambda_i f_ij
        Gamma_j  = (1/lambda_j) [ 1{a0_j=ab_j} (lambda0_j Gamma0_j)^(p_j)
                    + sum_i 1{ab_i=ab_j} (lambda_i Gamma_i)^(p_i) f_ij ]^(1/p_j)

    with p_j = ab_j/(ab_j-1) and ab_j the minimum tail coefficient among the
    streams feeding j: the heaviest tail upstream dominates. Both are solved
    exactly rather than iteratively. The rate equations are the usual traffic
    equations, and in the variables z_j = (lambda_j Gamma_j)^(p_j) the
    variability equations are linear as well, so each is one linear system; ab is
    obtained by propagating the minimum to a fixed point.

    Args:
        lambda0: (J,) external arrival rate at each node, 0 where there is none
        Gamma0: (J,) variability parameter of each external arrival process,
            which for a renewal stream is the interarrival standard deviation
        alpha0: (J,) tail coefficient in (1,2] of each external arrival process
        F: (J,J) routing probability matrix, F[i,j] = fraction of the jobs
            leaving node i that go to node j (row sums <= 1)

    Returns:
        Tuple of (lambda, Gamma, alpha), each (J,)
    """
    lambda0 = np.asarray(lambda0, dtype=np.float64).ravel()
    Gamma0 = np.asarray(Gamma0, dtype=np.float64).ravel()
    alpha0 = np.asarray(alpha0, dtype=np.float64).ravel()
    F = np.asarray(F, dtype=np.float64)
    J = lambda0.size

    # traffic equations
    lam = np.linalg.solve(np.eye(J) - F.T, lambda0)
    lam[np.abs(lam) < np.finfo(float).eps] = 0.0

    # effective tail coefficient: the minimum propagated along the routing graph
    alpha = np.where(lambda0 > 0, alpha0, np.inf)
    for _ in range(J):
        prev = alpha.copy()
        for j in range(J):
            for i in range(J):
                if F[i, j] > 0 and lam[i] > 0:
                    alpha[j] = min(alpha[j], alpha[i])
        if np.array_equal(prev, alpha):
            break
    alpha[~np.isfinite(alpha)] = 2.0  # an unreachable node keeps the light-tailed default

    # variability equations, linear in z_j = (lambda_j Gamma_j)^(p_j)
    p = alpha / (alpha - 1)
    z0 = np.zeros(J)
    for j in range(J):
        if lambda0[j] > 0 and abs(alpha0[j] - alpha[j]) < 1e-12:
            z0[j] = (lambda0[j] * Gamma0[j]) ** p[j]
    A = np.zeros((J, J))
    for i in range(J):
        for j in range(J):
            if F[i, j] > 0 and abs(alpha[i] - alpha[j]) < 1e-12:
                A[i, j] = F[i, j]
    z = np.linalg.solve(np.eye(J) - A.T, z0)
    z[z < 0] = 0.0

    Gamma = np.zeros(J)
    for j in range(J):
        if lam[j] > 0:
            Gamma[j] = z[j] ** (1.0 / p[j]) / lam[j]
    return lam, Gamma, alpha
