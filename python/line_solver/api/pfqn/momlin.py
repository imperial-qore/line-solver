"""
Moment linearizer: approximate mean queue lengths and their second moments
(variance / covariance) for large closed product-form queueing networks.

Means come from the Schweitzer-Bard AMVA fixed point. Second moments use the
exact product-form identity ``Cov[n_ir, n_js] = D_js dQ_ir/dD_js``, with the
demand derivatives obtained by analytically linearizing the AMVA fixed point
(a "moment linearizer" in the sense of Strelen and Akyildiz, developed per
class). Both moments carry the AMVA approximation error; for exact moments on
tractable models use :func:`pfqn_sens`.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_momlin.m`` and the JAR ``Pfqn_momlin``.
"""

import numpy as np


class PfqnMomlin:
    """Result container for :func:`pfqn_momlin`.

    Attributes
    ----------
    Q : np.ndarray (M x R)
        Mean queue length.
    X : np.ndarray (1 x R)
        Throughput per class.
    U, R : np.ndarray (M x R)
        Utilization and residence time.
    QVar : np.ndarray (M x R)
        Queue-length variance.
    QCov : np.ndarray (M x R x M x R)
        Queue-length covariance, ``QCov[i, r, j, s] = Cov[n_ir, n_js]``.
    dQ : np.ndarray (M x R x M x R)
        Demand derivatives, ``dQ[i, r, j, s] = dQ_ir/dD_js``.
    """

    def __init__(self, Q, X, U, R, QVar, QCov, dQ):
        self.Q = Q
        self.X = X
        self.U = U
        self.R = R
        self.QVar = QVar
        self.QCov = QCov
        self.dQ = dQ


def pfqn_momlin(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                tol: float = 1e-8, maxiter: int = 1000) -> PfqnMomlin:
    """Approximate first and second queue-length moments of a closed
    product-form network, scalable to large populations and many classes.

    Parameters
    ----------
    L : (M, R) array
        Service demand matrix.
    N : (R,) array
        Closed population per class.
    Z : (R,) array, optional
        Think time per class (default zeros).
    tol : float
        Convergence tolerance on the queue-length fixed point.
    maxiter : int
        Maximum iterations.

    Returns
    -------
    PfqnMomlin
    """
    L = np.asarray(L, dtype=np.float64)
    N = np.ceil(np.asarray(N, dtype=np.float64).flatten()).astype(int)
    R = len(N)
    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)
    M = L.shape[0]
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()
    if np.any(np.isinf(N)):
        raise ValueError("pfqn_momlin supports closed classes only.")

    # Schweitzer coefficients c[r, s] = (N_s - delta_rs)/N_s
    c = np.ones((R, R))
    for r in range(R):
        for s in range(R):
            if N[s] > 0:
                c[r, s] = (N[s] - (1.0 if r == s else 0.0)) / N[s]
            else:
                c[r, s] = 0.0

    # ---- Schweitzer-Bard AMVA fixed point for the means -----------------
    Q = np.zeros((M, R))
    for r in range(R):
        if N[r] > 0:
            Q[:, r] = N[r] / M
    X = np.zeros(R)
    Rmat = np.zeros((M, R))
    for _ in range(maxiter):
        Qold = Q.copy()
        for r in range(R):
            if N[r] == 0:
                X[r] = 0.0
                Rmat[:, r] = 0.0
                continue
            for i in range(M):
                Rmat[i, r] = L[i, r] * (1.0 + c[r, :] @ Q[i, :])
            X[r] = N[r] / (Z[r] + Rmat[:, r].sum())
            Q[:, r] = X[r] * Rmat[:, r]
        if np.max(np.abs(Q - Qold)) < tol:
            break
    U = np.zeros((M, R))
    for r in range(R):
        U[:, r] = X[r] * L[:, r]

    # ---- analytic linearization of the fixed point ----------------------
    # dQ[:, :, j, s0] = dQ_ir / dD_{j, s0}
    dQ = np.zeros((M, R, M, R))
    for j in range(M):
        for s0 in range(R):
            if N[s0] == 0:
                continue   # empty class: derivative 0
            dq = np.zeros((M, R))
            for _ in range(maxiter):
                dqold = dq.copy()
                for r in range(R):
                    if N[r] == 0:
                        continue
                    dR = np.zeros(M)
                    for i in range(M):
                        dDir = 1.0 if (i == j and r == s0) else 0.0
                        dR[i] = dDir * (1.0 + c[r, :] @ Q[i, :]) + L[i, r] * (c[r, :] @ dq[i, :])
                    dXr = -(X[r] ** 2 / N[r]) * dR.sum()   # dZ/dtheta = 0
                    dq[:, r] = dXr * Rmat[:, r] + X[r] * dR
                if np.max(np.abs(dq - dqold)) < tol:
                    break
            dQ[:, :, j, s0] = dq

    # ---- second moments via the covariance identity ---------------------
    QCov = np.zeros((M, R, M, R))
    for i in range(M):
        for r in range(R):
            for j in range(M):
                for s in range(R):
                    QCov[i, r, j, s] = L[j, s] * dQ[i, r, j, s]
    QVar = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            QVar[i, r] = QCov[i, r, i, r]

    return PfqnMomlin(Q, X.reshape(1, R), U, Rmat, QVar, QCov, dQ)
