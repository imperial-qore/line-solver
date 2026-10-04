"""Multivariate normal rectangle probability for the moment-closure methods.

Twin of `matlab/src/solvers/FLD/fluid_mvn_rectangle.m` and of
`jline.solvers.fluid.moments.MvnRectangle`.
"""

import numpy as np
from scipy.special import ndtr, ndtri

__all__ = ['mvn_rectangle']

# first 100 primes, listed rather than sieved so that the MATLAB and Java twins
# generate the identical lattice
_PRIMES = np.array([
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541],
    dtype=float)


def _chol_psd(C, dtol):
    """Cholesky factor of a symmetric positive SEMI-definite matrix.

    A vanishing pivot leaves a zero row/column, which the caller reads as a
    deterministic coordinate rather than as a failure.
    """
    d = C.shape[0]
    L = np.zeros((d, d))
    for i in range(d):
        v = C[i, i] - float(L[i, :i] @ L[i, :i])
        if v > dtol:
            L[i, i] = np.sqrt(v)
            for j in range(i + 1, d):
                L[j, i] = (C[j, i] - float(L[j, :i] @ L[i, :i])) / L[i, i]
        else:
            L[i, i] = 0.0
            L[i + 1:, i] = 0.0
    return L


def mvn_rectangle(m, C, a, b, npoints=4096):
    """Rectangle probability P(a <= Y <= b) for Y ~ Normal(m, C).

    The integral has no closed form beyond one dimension, so it is evaluated by
    the separation-of-variables transformation of Genz (1992): the Cholesky
    factor of C turns the rectangle into an iterated integral over the unit cube
    whose integrand is a product of normal-CDF differences, and the first
    coordinate is integrated exactly. The remaining cube is integrated with a
    DETERMINISTIC Richtmyer lattice rule, frac(k*sqrt(p_j)) over the first
    primes, averaged with its antithetic reflection. Determinism is required
    here, not merely convenient: the MATLAB, Java and Python twins must return
    the same number, and a randomized rule would make them agree only in
    distribution.

    C may be SINGULAR, which is the common case: a closed population fixes the
    sum of the station coordinates, so the covariance of a station holding a
    whole class is rank deficient. A coordinate whose CONDITIONAL variance
    vanishes is not integrated; it is a hard constraint, contributing 1 when the
    conditional mean falls inside its interval and 0 otherwise.

    Args:
        m: mean vector, length d
        C: covariance matrix, d-by-d, symmetric positive semi-definite
        a: lower corner, length d, -inf allowed
        b: upper corner, length d, +inf allowed
        npoints: lattice points per antithetic pair

    Returns:
        (p, logp) with p in [0,1] and logp = log(p), -inf when p is zero
    """
    m = np.asarray(m, dtype=float).ravel()
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    C = np.atleast_2d(np.asarray(C, dtype=float))
    d = m.size
    if d == 0:
        return 1.0, 0.0
    al = a - m
    bu = b - m
    if np.any(bu <= al):
        return 0.0, -np.inf

    # scale-relative tolerances: dtol decides which coordinate carries noise,
    # ctol whether a deterministic coordinate satisfies its constraint
    scale = max(1.0, float(np.max(np.abs(np.diag(C)))))
    dtol = 1e-12 * scale
    ctol = 1e-6 * np.sqrt(scale)

    L = _chol_psd(C, dtol)
    isInt = np.diag(L) > 0
    intIdx = np.where(isInt)[0]
    nInt = intIdx.size
    lastInt = intIdx[-1] if nInt > 0 else -1
    nw = max(0, nInt - 1)

    if nw == 0:
        W = np.zeros((1, 0))
    else:
        if nw > _PRIMES.size:
            raise ValueError(
                "The lattice rule carries generators for at most %d integration dimensions, "
                "but this rectangle has %d. Aggregate classes before evaluating the cell."
                % (_PRIMES.size, nw))
        k = np.arange(1, npoints + 1, dtype=float).reshape(-1, 1)
        Wbase = np.mod(k * np.sqrt(_PRIMES[:nw]).reshape(1, -1), 1.0)
        W = np.vstack([Wbase, 1.0 - Wbase])  # antithetic reflection

    Np = W.shape[0]
    y = np.zeros((Np, d))
    f = np.ones(Np)
    kw = 0
    for i in range(d):
        s = y[:, :i] @ L[i, :i] if i > 0 else np.zeros(Np)
        if isInt[i]:
            lo = (al[i] - s) / L[i, i]
            hi = (bu[i] - s) / L[i, i]
            dd = ndtr(lo)
            ee = ndtr(hi)
            f = f * np.maximum(0.0, ee - dd)
            if i != lastInt:
                u = dd + W[:, kw] * (ee - dd)
                kw += 1
                # the inverse CDF is evaluated strictly inside the unit interval
                y[:, i] = ndtri(np.clip(u, 1e-15, 1.0 - 1e-15))
        else:
            # zero conditional variance: the coordinate is pinned at s, so the
            # cell is either met or not
            f[(s < al[i] - ctol) | (s > bu[i] + ctol)] = 0.0

    p = float(np.mean(f))
    p = min(max(p, 0.0), 1.0)
    return p, (float(np.log(p)) if p > 0 else -np.inf)
