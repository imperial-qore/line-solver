"""
Exact passage-time law along an overtake-free path of a closed tree-like
product-form network.

Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions in
Large Markov Chains", 2002, Sec. 7.1, Theorems 1 and 2, after P. G. Harrison,
J. Appl. Prob. 27, 1990 and H. Duduna, Adv. Appl. Prob. 14, 1982. The
underlying sojourn-time result for overtake-free paths is F. Kelly and
P. Pollett, Adv. Appl. Prob. 15, 1983.

Twin of the MATLAB matlab/src/api/pfqn/pfqn_cyclet_ofree.m. Node indices are
0-based here and 1-based there.
"""

from math import factorial
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.special import gammainc

from ..lti import laplace_invert_cdf, laplace_invert_pdf


def _buzen(y: Sequence[complex], n: int) -> np.ndarray:
    """
    Buzen's convolution: g[k] = G at population k for the node set y, k = 0..n.
    This is the k(y,a,b) recursion of Sec. 7.1 with the node index rolled up,
    k(y,a,b) = k(y,a-1,b) + y_a k(y,a,b-1), k(y,a,0) = 1, k(y,0,b>0) = 0.
    """
    y = np.asarray(y)
    dtype = complex if np.iscomplexobj(y) else float
    g = np.zeros(n + 1, dtype=dtype)
    g[0] = 1.0
    for yi in y:
        for k in range(1, n + 1):
            g[k] = g[k] + yi * g[k - 1]
    return g


def _series_mul(a: np.ndarray, b: np.ndarray, K: int) -> np.ndarray:
    c = np.zeros(K + 1)
    for i in range(K + 1):
        if a[i] == 0.0:
            continue
        for j in range(K + 1 - i):
            c[i + j] += a[i] * b[j]
    return c


def _thm2(v, mu, N, z, tset, Gn1, x):
    """
    Theorem 2 in closed form. The density is a finite sum of terms
    t^k exp(-mu_j t), so its integral is an incomplete gamma and the CDF comes
    out in closed form too rather than by quadrature.
    """
    M = len(v)
    m = len(z)
    off = [i for i in range(M) if i not in z]
    mup = mu[list(z)]
    vp = v[list(z)]

    Gm = _buzen(x[off], N - 1)              # network minus the path
    coef = np.zeros((m, N))                 # coef[j,k] multiplies t^k exp(-mu_j t)
    for j in range(m):
        den = 1.0
        for i in range(m):
            if i != j:
                den *= (mup[i] - mup[j])
        if den == 0.0:
            raise ValueError("Theorem 2 needs distinct service rates on the "
                             "path; two coincide. Use method='lt'.")
        idx = [i for i in range(m) if i != j]
        w = (vp[idx] - vp[j]) / (mup[idx] - mup[j])
        K = _buzen(w, N - 1)                # K^m(j,l), l = 0..N-1
        for c in range(N):
            Gmc = Gm[N - 1 - c]             # G_m(N-c-1)
            if Gmc == 0.0:
                continue
            for i in range(c + 1):
                coef[j, c - i] += Gmc * K[i] / den

    pref = float(np.prod(mup)) / Gn1
    t = np.asarray(tset, dtype=float)
    f = np.zeros(t.shape)
    F = np.zeros(t.shape)
    for j in range(m):
        for k in np.flatnonzero(coef[j, :] != 0.0):
            cjk = coef[j, k]
            f += pref * cjk * (vp[j] ** k) * (t ** k) / factorial(int(k)) * np.exp(-mup[j] * t)
            # int_0^t s^k exp(-mu s) ds = k!/mu^(k+1) * P(k+1, mu t)
            F += pref * cjk * (vp[j] ** k) / (mup[j] ** (k + 1)) * gammainc(k + 1, mup[j] * t)
    return np.maximum(f, 0.0), np.clip(F, 0.0, 1.0)


def _lst(v, mu, N, z, s, Gn1):
    """L(s|z) = prod_j mu_j/(s+mu_j) * G(y(s), N-1) / G(x, N-1)."""
    M = len(v)
    y = (v / mu).astype(complex)
    for j in z:
        y[j] = y[j] * (mu[j] / (s + mu[j]))
    L = _buzen(y, N - 1)[-1] / Gn1
    for j in z:
        L = L * (mu[j] / (s + mu[j]))
    return L


def _moments(v, mu, N, z, Gn1, x, nmom):
    """
    The same Buzen convolution run in the ring of truncated power series in s.
    Every operation in the recursion is an addition or a multiplication, so the
    series ring carries it unchanged, and E[T^q] = (-1)^q q! [s^q] L(s).
    """
    K = nmom
    M = len(v)
    Y = np.zeros((M, K + 1))
    for i in range(M):
        if i in z:
            Y[i, :] = x[i] * ((-1.0 / mu[i]) ** np.arange(K + 1))
        else:
            Y[i, 0] = x[i]

    G = np.zeros((N, K + 1))
    G[0, 0] = 1.0
    for i in range(M):
        for n in range(1, N):
            G[n, :] = G[n, :] + _series_mul(Y[i, :], G[n - 1, :], K)

    L = G[N - 1, :] / Gn1
    for j in z:
        e = (-1.0 / mu[j]) ** np.arange(K + 1)
        L = _series_mul(L, e, K)
    return np.array([((-1) ** q) * factorial(q) * L[q] for q in range(1, nmom + 1)])


def pfqn_cyclet_ofree(v, mu, N, path, tset, method: str = 'auto',
                      nmom: int = 3, pathprob=None, lti_method: str = 'euler',
                      tol: float = 1e-8):
    """
    Exact passage-time density, CDF and moments along an OVERTAKE-FREE PATH of a
    closed single-chain tree-like product-form network with population N.

    v, mu are per-node visit ratios and service rates, path is the node list
    z = (z_1, ..., z_m) with z_1 the root, tset the time grid. path may instead
    be a sequence of paths, in which case pathprob weights them and the outputs
    are the mixture; that is how a cycle time is assembled when the root
    branches.

    THE ONE FACT THAT MAKES ALL THREE ROUTES WORK. Conditional on the path,

        T | z  =  sum_{j in z} Erlang(u_{z_j} + 1, mu_{z_j})

    with u distributed as the network's equilibrium population vector AT N-1
    (the arrival theorem). Hence the transform of Theorem 1 collapses to

        L(s|z) = prod_{j in z} mu_j/(s+mu_j) * G(y(s), N-1) / G(x, N-1)

    where x_i = v_i/mu_i and y_i(s) = x_i mu_i/(s+mu_i) on the path, x_i off it.
    One Buzen convolution per value of s.

    method: 'auto' (default) uses 'exact' when the path rates are separated and
    'lt' otherwise; 'exact' is Theorem 2 in closed form and REQUIRES DISTINCT
    RATES on the path, since its partial fractions divide by
    prod_{i!=j}(mu_i - mu_j); 'lt' inverts the transform above through api/lti.

    MOMENTS ARE NEVER TAKEN FROM THE DENSITY. They come from running the same
    Buzen convolution in the ring of truncated power series in s, so they are
    exact to machine precision, are unaffected by the time grid, and stay valid
    when the rates coincide and Theorem 2 does not apply.

    NOTE ON THE PAPER. The inner sum of Theorem 2 reads (v_j t)^(c-i)/(c-i)! and
    that is CORRECT as printed, however odd the visit ratio looks against a
    time: substituting the service rate instead returns negative densities.
    Verified against a direct mixture-of-Erlangs oracle to 1e-15, and at the
    paper's own N = 18 example against the transform route to 1e-11.

    Returns (f, F, mom, out).
    """
    v = np.asarray(v, dtype=float).ravel()
    mu = np.asarray(mu, dtype=float).ravel()
    M = v.size
    if mu.size != M:
        raise ValueError("v and mu must name the same number of nodes.")
    if np.any(mu <= 0):
        raise ValueError("Every service rate must be positive.")
    if N < 1 or int(N) != N:
        raise ValueError("The population N must be a positive integer.")
    N = int(N)

    first = path[0] if len(path) else None
    if isinstance(first, (list, tuple, np.ndarray)):
        paths = [list(np.asarray(p, dtype=int).ravel()) for p in path]
    else:
        paths = [list(np.asarray(path, dtype=int).ravel())]
    if pathprob is None:
        pathprob = np.ones(len(paths)) / len(paths) if len(paths) > 1 else np.ones(1)
    pathprob = np.asarray(pathprob, dtype=float).ravel()
    if pathprob.size != len(paths):
        raise ValueError("pathprob must carry one probability per path.")

    t = np.atleast_1d(np.asarray(tset, dtype=float)).ravel()
    f = np.zeros(t.shape)
    F = np.zeros(t.shape)
    mom = np.zeros(nmom)
    out: List[Dict] = []

    x = v / mu
    Gn1 = float(_buzen(x, N - 1)[-1].real)
    if Gn1 <= 0:
        raise ValueError("The network normalizing constant at population N-1 "
                         "vanished; check v and mu.")

    for ip, z in enumerate(paths):
        if len(z) == 0:
            raise ValueError("An overtake-free path must contain at least the "
                             "root node.")
        if min(z) < 0 or max(z) >= M or len(set(z)) != len(z):
            raise ValueError("A path must be a set of distinct node indices "
                             "within the network.")

        m = method
        if m == 'auto':
            mp = mu[list(z)]
            if len(mp) == 1:
                m = 'exact'
            else:
                d = np.abs(mp[:, None] - mp[None, :])
                np.fill_diagonal(d, np.inf)
                m = 'exact' if d.min() > tol * mp.max() else 'lt'

        if m == 'exact':
            fi, Fi = _thm2(v, mu, N, z, t, Gn1, x)
        elif m == 'lt':
            Lfun = lambda s, zz=z: _lst(v, mu, N, zz, s, Gn1)
            fi = laplace_invert_pdf(Lfun, t, lti_method)
            Fi = laplace_invert_cdf(Lfun, t, lti_method)
        else:
            raise ValueError("Unknown method: %s. Supported: auto, exact, lt."
                             % method)

        momi = _moments(v, mu, N, z, Gn1, x, nmom)
        f = f + pathprob[ip] * fi
        F = F + pathprob[ip] * Fi
        mom = mom + pathprob[ip] * momi
        out.append({'method': m, 'lG': float(np.log(Gn1)), 'path': list(z)})

    return f, F, mom, out
