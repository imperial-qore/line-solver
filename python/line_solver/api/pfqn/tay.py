"""Tay's arrival-instant approximate mean value analysis.

Presented as eqs. 4.8.2-1..3 of the Schweitzer-Serazzi-Broglia survey of
bottleneck analysis and MVA generalisations, where it is benchmarked against
exact, Linearizer and Bard-Schweitzer on Tay's Example 4.
"""

from typing import Optional, Tuple

import numpy as np

__all__ = ['pfqn_tay']


def pfqn_tay(L, N, Z=None, tol: float = 1e-6, maxiter: int = 1000,
             QN0: Optional[np.ndarray] = None):
    """Approximate MVA whose arrival-instant queue lengths come from the
    THROUGHPUT ELASTICITIES rather than from a population-shift heuristic.

    Let E_mkc = (D_mk/X_c) dX_c/dD_mk be the elasticity of the class-c
    throughput with respect to the class-k demand at station m. Tay shows that
    the elasticities satisfy the R linear equations

        E_mkj sum_t B_tj Q_jt (1+Q_jt) =
            -[(delta_jk + Q_jm) B_mk Q_km
              + sum_{c!=j} E_mkc sum_t B_tc Q_jt Q_ct]

    with B_ir = 1/(1 + D_ir X_r/N_r), and that the arrival-instant queue length
    is then simply Q_km^(r) = Q_km + E_mkr, which closes the MVA recursion
    R_rm = D_rm (1 + sum_k Q_km^(r)). One R x R solve per (station, class) pair
    per iteration.

    Delay stations enter through Z only. They are "AS" servers in the survey's
    notation (d_t = 0), so they contribute Z_j X_j to the denominator of the
    elasticity equations but nothing to its numerator.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).
        Z: Think time vector (R,). Default: zeros.
        tol: Convergence tolerance on the queue lengths.
        maxiter: Maximum number of iterations.
        QN0: Initial guess for the queue lengths (M x R).

    Returns:
        Tuple (XN, QN, UN, RN, it, QNarr). QNarr[m,k,r] is the class-k queue
        length at station m as seen by an arriving class-r job: the auxiliary
        quantity the method is tabulated on, NOT the queue length of the model
        re-solved at N - e_r (the same object only for an exact solution).
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float)
        if Z.ndim > 1 and Z.shape[0] > 1:
            Z = Z.sum(axis=0)          # several delay stations aggregate
        Z = Z.ravel()
    if Z.size < R:
        Z = np.full(R, float(Z.sum()))

    XN = np.zeros(R)
    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    RN = np.zeros((M, R))
    Qarr = np.zeros((M, R, R))
    it = 0

    # Empty classes contribute no jobs anywhere and make the elasticity system
    # singular (their denominator is identically zero); solve without them and
    # re-expand, as in pfqn_bs.
    act = np.flatnonzero(N > 0)
    if act.size == 0:
        return XN, QN, UN, RN, it, Qarr
    if act.size < R:
        Xa, Qa, Ua, Ra, it, Qarra = pfqn_tay(L[:, act], N[act], Z[act], tol, maxiter)
        XN[act] = Xa
        QN[:, act] = Qa
        UN[:, act] = Ua
        RN[:, act] = Ra
        Qarr[np.ix_(np.arange(M), act, act)] = Qarra
        return XN, QN, UN, RN, it, Qarr

    if QN0 is None:
        QN = np.tile(N, (M, 1)) / M
    else:
        QN = np.array(QN0, dtype=float)
    XN = N / (Z + L.sum(axis=0) * (1 + QN.sum(axis=0)))

    for it in range(1, maxiter + 1):
        QN_1 = QN.copy()

        B = 1.0 / (1.0 + L * (XN / N)[np.newaxis, :])

        # Denominators of the elasticity system; the delay term Z_j X_j is the
        # AS-server contribution (d_t = 0 leaves B = 1).
        den = (B * QN * (1 + QN)).sum(axis=0) + Z * XN

        # Off-diagonal coupling C[j,c] = sum_t B_tc Q_jt Q_ct
        C = np.einsum('tc,tj,tc->jc', B, QN, QN)

        for m in range(M):
            for k in range(R):
                A = C / den[:, np.newaxis]
                np.fill_diagonal(A, 1.0)
                b = -(np.eye(R)[:, k] + QN[m, :]) * B[m, k] * QN[m, k] / den
                Qarr[m, k, :] = QN[m, k] + np.linalg.solve(A, b)

        for r in range(R):
            RN[:, r] = L[:, r] * (1 + Qarr[:, :, r].sum(axis=1))
        XN = N / (Z + RN.sum(axis=0))
        QN = RN * XN[np.newaxis, :]

        if np.max(np.abs(QN - QN_1)) < tol:
            break

    UN = L * XN[np.newaxis, :]
    return XN, QN, UN, RN, it, Qarr
