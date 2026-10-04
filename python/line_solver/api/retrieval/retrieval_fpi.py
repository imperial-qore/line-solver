"""Fixed-point (FPI) approximation of delayed-hit cache metrics.

Native-Python port of retrieval_fpi.m / Retrieval_fpi.java.
"""
import numpy as np


def _reldiff(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.max(np.abs(b)) if b.size else 0.0
    if denom == 0:
        denom = 1.0
    return np.max(np.abs(a - b)) / denom if a.size else 0.0


def retrieval_fpi(m, lambda_, eta, gamma, max_iter=1000, tol=1e-6):
    """Return (pmiss (n,), phit (h,n), pdh (r+1,n))."""
    m = np.asarray(m, dtype=float).ravel()
    lambda_ = np.asarray(lambda_, dtype=float).ravel()
    eta = np.asarray(eta, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    n = len(lambda_)
    h = len(m)
    r = eta.shape[1] - 1
    eta0 = eta[:, 0]
    etaPS = eta[:, 1:]

    phi = np.full((r + 1, n), 1.0 / ((h + 1) * (r + 2)))
    pij = np.full((h, n), 1.0 / (h + 1))
    pi0 = np.full(n, 1.0 / ((h + 1) * (r + 2)))

    for _ in range(max_iter):
        F = np.ones((r, n))
        for s in range(r):
            phis = phi[s + 1, :]
            F[s, :] = 1 + (phis.sum() - phis)
        D = 1 + lambda_ * eta0
        for s in range(r):
            D = D + lambda_ * etaPS[:, s] * F[s, :]
        theta = gamma / D[:, None]
        oneminus = (1 - pij.sum(axis=0))
        xi = np.zeros(h)
        for j in range(h):
            xi[j] = m[j] / np.sum(theta[:, j] * oneminus)
        txi = theta * xi[None, :]
        pij_new = (txi / (1 + txi.sum(axis=1))[:, None]).T
        pi0_new = (1 - pij_new.sum(axis=0)) / D
        phi_new = np.zeros((r + 1, n))
        phi_new[0, :] = lambda_ * eta0 * pi0_new
        for s in range(r):
            phi_new[s + 1, :] = lambda_ * etaPS[:, s] * F[s, :] * pi0_new
        delta = max(_reldiff(pi0_new, pi0), _reldiff(pij_new, pij), _reldiff(phi_new, phi))
        pij, pi0, phi = pij_new, pi0_new, phi_new
        if not np.isfinite(delta):
            break
        if delta < tol:
            break
    return pi0, pij, phi
