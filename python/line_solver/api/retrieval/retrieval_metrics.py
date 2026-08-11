"""Exact miss/hit/delayed-hit metrics of a delayed-hit (list-based) cache.

Native-Python port of retrieval_metrics.m / Retrieval_metrics.java.
"""
import numpy as np

from .retrieval_nc import retrieval_nc


def retrieval_metrics(m, lambda_, eta, gamma):
    """Return (pmiss (n,), phit (h,n), pdh (r+1,n))."""
    m = np.asarray(m, dtype=float).ravel()
    lambda_ = np.asarray(lambda_, dtype=float).ravel()
    eta = np.asarray(eta, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    n = len(lambda_)
    h = len(m)
    r = eta.shape[1] - 1
    v0 = np.zeros(r)

    E = retrieval_nc(v0, m, lambda_, eta, gamma)
    pmiss = np.zeros(n)
    phit = np.zeros((h, n))
    pdh = np.zeros((r + 1, n))

    for i in range(n):
        keep = [k for k in range(n) if k != i]
        lambda_i = lambda_[keep]
        eta_i = eta[keep, :]
        gamma_i = gamma[keep, :]

        Ei = retrieval_nc(v0, m, lambda_i, eta_i, gamma_i)
        pmiss[i] = Ei / E
        pdh[0, i] = lambda_[i] * eta[i, 0] * Ei / E
        for s in range(1, r + 1):
            vs = np.zeros(r)
            vs[s - 1] = 1
            Eis = retrieval_nc(vs, m, lambda_i, eta_i, gamma_i)
            pdh[s, i] = lambda_[i] * eta[i, s] * Eis / E
        for j in range(h):
            if m[j] > 0:
                mp = m.copy()
                mp[j] -= 1
                Eij = retrieval_nc(v0, mp, lambda_i, eta_i, gamma_i)
                phit[j, i] = m[j] * gamma[i, j] * Eij / E
    return pmiss, phit, pdh
