"""Ray (WKB) asymptotic expansion of the list-based cache normalizing constant.

Native-Python port of matlab/src/api/retrieval/retrieval_rayint.m (and
jar/.../jline/api/retrieval/Retrieval_rayint.java).

Approximates the constant that ``cache_erec`` computes exactly, in the SAME
normalization, so the two are interchangeable::

    E(m,n) = E(m,n-1) + sum_j gamma_{n,j} m_j E(m-1_j,n-1),   E(0,0)=1

Writing ``E = prod_j m_j! * Et``, the relaxation ``Et ~ H exp(phi/eps)`` with
``n = y/eps`` and ``m_j = x_j/eps`` gives the eikonal
``e^{phi_y} = 1 + sum_j gamma_j(y) e^{-phi_j}``, whose rays carry the constants
``xi_j = e^{-phi_j}``.  With ``S(v) = 1 + sum_j gamma_j(v) xi_j``::

    x_j    = int_0^y gamma_j(v) xi_j / S(v) dv          (the saddle conditions)
    phi    = int_0^y log S(v) dv - sum_j x_j log xi_j
    H      = (2 pi)^{-h/2} sqrt(S(y)/S(0)) / sqrt(prod_j xi_j * det A)
    A_{ik} = d x_i / d xi_k

and ``E ~ prod_j m_j! * eps^{h/2} H exp(phi/eps)``.

DISCRETE (``gamma`` an ``n x h`` array).  The ray integrals are the sums they
discretize and the expansion collapses to the Laplace form::

    Et ~ (2 pi)^{-h/2} exp(sum_k log D_k - sum_j m_j log xi_j) / sqrt(det Sigma)

with ``D_k = 1 + sum_j gamma_{k,j} xi_j``, ``sum_k gamma_{k,j} xi_j / D_k = m_j``
and ``Sigma = A diag(xi)`` the Hessian in ``log xi``.  This is the more accurate
of the two forms; the ``sqrt(S(y)/S(0))`` factor is exactly the Euler-Maclaurin
term relating ``sum_k`` to ``int dv`` and is already accounted for.

CONTINUUM (``gamma`` a callable on ``v in [0,1]``).  Composite Simpson quadrature
on the profile itself, the form written in the note.  Costs roughly a factor two
in accuracy but does not need the ``n`` rows.

ACCURACY.  The relative error is ``O(1/n)`` at fixed occupancy but is governed by
the smallest occupancy rather than by ``n``, tracking
``0.14 * (1/min_j m_j + 1/(n - sum_j m_j))``, so a per cent needs every ``m_j``
and ``n - sum_j m_j`` above about 15 and a part in a thousand needs them above
about 150.  Returned as ``out.relerrEst``.  Lists with ``m_j = 0`` contribute
nothing and are dropped before the saddle is solved.

This is the no-fetch (q=0) case, i.e. the same quantity as ``cache_erec``.  The
delayed-hit extension carrying the fetch coordinates is NOT implemented: its
eikonal is known but its amplitude has not been derived.
"""
import warnings
from types import SimpleNamespace

import numpy as np
from scipy.special import gammaln


def retrieval_rayint(gamma, m, n=None, nquad=4097):
    """Ray expansion of the cache normalizing constant.

    Parameters
    ----------
    gamma : ndarray (n, h) or callable
        Access factors ``gamma[k, j]``, or a profile ``v -> (len(v), h)`` on
        ``v in [0, 1]``.
    m : array_like
        Cache list capacities, length h.
    n : int, optional
        Number of items.  Required, and used only, with a callable ``gamma``.
    nquad : int, optional
        Composite Simpson nodes for the continuum form (default 4097, forced odd).

    Returns
    -------
    E : float
        Normalizing constant, same normalization as ``cache_erec``.
    logE : float
        Its natural logarithm, safe for large n.
    out : SimpleNamespace
        ``xi``, ``phi``, ``logdetSigma``, ``S0``, ``Sy``, ``method``,
        ``relerrEst``, ``iter``.
    """
    m = np.asarray(m, dtype=float).ravel()
    if np.any(m < 0) or np.any(np.abs(m - np.rint(m)) > 0):
        raise ValueError("retrieval_rayint: list capacities must be non-negative integers.")

    isprofile = callable(gamma)
    if isprofile:
        if n is None:
            raise ValueError("retrieval_rayint: the number of items n must be given when the "
                             "access factors are a callable.")
        nquad = max(5, 2 * (int(nquad) // 2) + 1)
        h = np.asarray(gamma(np.zeros(1)), dtype=float).reshape(1, -1).shape[1]
    else:
        gamma = np.asarray(gamma, dtype=float)
        if gamma.ndim != 2 or gamma.size == 0:
            raise ValueError("retrieval_rayint: the access factors must be a non-empty (n, h) array.")
        n, h = gamma.shape
    if m.size != h:
        raise ValueError("retrieval_rayint: the capacity vector must have one entry per cache list "
                         "(%d given, %d expected)." % (m.size, h))

    out = SimpleNamespace(xi=np.zeros(h), phi=np.nan, logdetSigma=np.nan,
                          S0=np.nan, Sy=np.nan, method='', relerrEst=np.nan, iter=0)

    msum = float(m.sum())
    if msum > n:
        out.method = 'boundary'
        return 0.0, -np.inf, out
    if msum == 0:
        out.method = 'boundary'
        return 1.0, 0.0, out
    if msum == n:
        raise ValueError("retrieval_rayint: the expansion requires sum(m) < n; at sum(m) = n the "
                         "saddle point escapes to infinity. Use cache_erec for a full cache.")

    keep = m > 0
    mk = m[keep]
    hk = int(mk.size)

    if isprofile:
        v = np.linspace(0.0, 1.0, nquad)
        w = np.ones(nquad)
        w[1:-1:2] = 4.0
        w[2:-2:2] = 2.0
        w /= 3.0 * (nquad - 1)
        G = np.asarray(gamma(v), dtype=float)
        if G.shape[0] != nquad:
            raise ValueError("retrieval_rayint: the access-factor callable must return one row per "
                             "evaluation point.")
        G = G[:, keep]
        tgt = mk / n
    else:
        G = gamma[:, keep]
        w = np.ones(G.shape[0])
        tgt = mk

    xi, iters = _saddle(G, w, tgt)
    S = 1.0 + G @ xi
    Sig = _hessian(G, w, xi)
    logdet = _logdet(Sig)

    if isprofile:
        phi = float(w @ np.log(S) - (mk / n) @ np.log(xi))
        logEt = (-0.5 * hk * np.log(n) - 0.5 * hk * np.log(2 * np.pi) + n * phi
                 - 0.5 * logdet + 0.5 * np.log(S[-1] / S[0]))
        out.method = 'rayint'
    else:
        phi = float(np.sum(np.log(S)) - mk @ np.log(xi))
        logEt = -0.5 * hk * np.log(2 * np.pi) + phi - 0.5 * logdet
        out.method = 'saddle'

    logE = float(logEt + gammaln(m + 1.0).sum())
    with np.errstate(over='ignore'):
        E = float(np.exp(logE))

    out.xi[keep] = xi
    out.phi = phi
    out.logdetSigma = float(logdet)
    out.S0 = float(S[0])
    out.Sy = float(S[-1])
    out.iter = int(iters)
    out.relerrEst = float(0.14 * (1.0 / mk.min() + 1.0 / (n - msum)))
    if min(mk.min(), n - msum) < 2:
        warnings.warn("retrieval_rayint: the smallest occupancy is %d, so the expansion is only "
                      "qualitative here (estimated relative error %.0f%%); cache_erec is exact."
                      % (int(min(mk.min(), n - msum)), 100 * out.relerrEst), RuntimeWarning)
    return E, logE, out


def _saddle(G, w, tgt):
    """Newton on theta = log xi for sum_k w_k G_kj xi_j / (1 + sum_l G_kl xi_l) = tgt_j.

    The objective sum_k w_k log(1 + sum_l G_kl e^{theta_l}) - tgt @ theta is strictly
    convex, so the root is unique and damped Newton converges globally.
    """
    tgt = np.asarray(tgt, dtype=float).ravel()
    gb = w @ G
    slack = max(1.0 - tgt.sum() / w.sum(), 1e-9)
    th = np.log(np.maximum(tgt, 1e-12) / np.maximum(gb * slack, 1e-12))
    tmax = max(1.0, float(np.max(np.abs(tgt))))
    it = 0
    for it in range(1, 201):
        xi = np.exp(th)
        a = G * xi
        S = 1.0 + a.sum(axis=1)
        g = (w[:, None] * a / S[:, None]).sum(axis=0) - tgt
        if np.max(np.abs(g)) <= 1e-12 * tmax:
            break
        H = _hessian(G, w, xi)
        d = -np.linalg.solve(H, g)
        if not np.all(np.isfinite(d)):
            raise ValueError("retrieval_rayint: the saddle-point Newton step is not finite; check "
                             "that the access factors are positive.")
        step = 1.0
        while np.max(np.abs(step * d)) > 2.0:
            step /= 2.0
        th = th + step * d
    return np.exp(th), it


def _hessian(G, w, xi):
    """Hessian in theta = log xi; equals A diag(xi) with A_{ik} = d x_i / d xi_k."""
    a = G * xi
    S = 1.0 + a.sum(axis=1)
    aS = a / S[:, None]
    return np.diag((w[:, None] * aS).sum(axis=0)) - aS.T @ (w[:, None] * aS)


def _logdet(H):
    try:
        L = np.linalg.cholesky(0.5 * (H + H.T))
    except np.linalg.LinAlgError:
        raise ValueError("retrieval_rayint: the saddle-point Hessian is not positive definite; "
                         "the ray map is singular here.")
    return 2.0 * float(np.sum(np.log(np.diag(L))))
