"""Exact normalizing constant of a delayed-hit (list-based) cache.

Native-Python port of matlab/src/api/retrieval/retrieval_nc.m (and
java/.../jline/api/retrieval/Retrieval_nc.java). Computes E(v,m) by the exact
recurrence (paper eq. mainrec); the plain constant is E(m) = E(0,m).
"""
import numpy as np


def retrieval_nc(v, m, lambda_, eta, gamma):
    """E(v,m).

    Parameters
    ----------
    v : array_like, moment-order vector for the PS stations (zeros(r) for the plain constant)
    m : array_like, cache list capacities (length h)
    lambda_ : array_like, per-item arrival rates (length n)
    eta : ndarray (n, r+1), col 0 = IS aggregate, cols 1..r = PS stations
    gamma : ndarray (n, h), access factors
    """
    v = np.asarray(v, dtype=float).ravel().copy()
    m = np.asarray(m, dtype=float).ravel().copy()
    lambda_ = np.asarray(lambda_, dtype=float).ravel()
    eta = np.asarray(eta, dtype=float)
    gamma = np.asarray(gamma, dtype=float)
    return _sub(v, m, lambda_, eta, gamma, len(lambda_))


def _sub(v, m, lambda_, eta, gamma, k):
    r = len(v)
    h = len(m)
    if m.sum() > k or m.min(initial=0.0) < 0:
        return 0.0
    if k == 0:
        return 1.0
    ki = k - 1
    E = (1.0 + lambda_[ki] * eta[ki, 0]) * _sub(v, m, lambda_, eta, gamma, k - 1)
    for s in range(1, r + 1):
        vp = v.copy()
        vp[s - 1] += 1
        E += lambda_[ki] * eta[ki, s] * (v[s - 1] + 1) * _sub(vp, m, lambda_, eta, gamma, k - 1)
    for j in range(h):
        if m[j] > 0:
            mp = m.copy()
            mp[j] -= 1
            E += gamma[ki, j] * m[j] * _sub(v, mp, lambda_, eta, gamma, k - 1)
    return E
