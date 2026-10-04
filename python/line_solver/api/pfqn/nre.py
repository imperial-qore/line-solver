"""
Normalizing constant via the saddle-tilted Edgeworth (NRE) approximation.

Port of matlab/src/api/pfqn/pfqn_nre.m. Evaluates the Norlund-Rice integral
form of the limited-load-dependent normalizing constant (Casale-Harrison-Ong,
Perform. Eval. 152, 2021, Thm. 6) by steepest descent instead of a Laplace
approximation on the untilted contour, as done by pfqn_nrl and pfqn_nrp. Two
corrections are applied over those methods:

  1. The integrand is invariant under t -> t + c*1, since h is homogeneous of
     degree sum(N) in the class variables and that degree cancels against
     exp(-1i*N*t). The redundant direction is quotiented out, so the integral
     is (R-1)-dimensional, not R-dimensional.
  2. The contour radii are tilted per class to the saddle point, i.e. to the X
     solving X_r*dlog(h)/dX_r = N_r, so that the origin is a stationary point
     of the phase.

A second-order Edgeworth term built from the third and fourth cumulants of the
tilted distribution is then added, giving a relative error of O(1/sum(N)^2)
instead of O(1) for a heuristic Gaussian fit.

All integrand evaluations are at real positive demands, so unlike pfqn_nrl and
pfqn_nrp this method needs no complex arithmetic and runs entirely in the log
domain through pfqn_lldsingle, whose cost is linear rather than quadratic in
the population wherever the rates settle.

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_nre.m
"""

from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np

__all__ = ['pfqn_nre', 'pfqn_nre_full', 'PfqnNreResult']


@dataclass
class PfqnNreResult:
    """The reference's ``[lG,G,lGs,vsad]``.

    ``lG - lGs`` is the Edgeworth correction, so a caller wanting the plain
    saddlepoint estimate reads ``lGs`` rather than re-deriving it. ``vsad`` is
    None on the shortcut arms that never solve a saddle point, matching the
    reference's empty ``[]``.
    """
    lG: float
    G: float
    lGs: float
    vsad: Optional[np.ndarray] = None

# finite-difference step, results are flat over [5e-3,5e-2]
_HSTEP = 2e-2
# beyond 8 classes the fourth-cumulant tensor is no longer affordable
_MAX_DIM = 7


def pfqn_nre(L: np.ndarray, N: np.ndarray, Z: Optional[np.ndarray] = None,
             alpha: Optional[np.ndarray] = None,
             options: Optional[Dict[str, Any]] = None) -> float:
    """
    Logarithm of the normalizing constant of a limited load-dependent model.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) or matrix (D x R) - optional
        alpha: Load-dependent rate matrix (M x sum(N))

    Returns:
        lG: Logarithm of the normalizing constant
    """
    return pfqn_nre_full(L, N, Z, alpha, options).lG


def pfqn_nre_full(L: np.ndarray, N: np.ndarray, Z: Optional[np.ndarray] = None,
                  alpha: Optional[np.ndarray] = None,
                  options: Optional[Dict[str, Any]] = None,
                  vfix: Optional[np.ndarray] = None) -> PfqnNreResult:
    """
    The full form of the reference's four outputs, named alike in the JAR and
    the C++ port.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) or matrix (D x R) - optional
        alpha: Load-dependent rate matrix (M x sum(N))
        options: Solver options
        vfix: Tilt to use instead of solving the saddle-point equation, None for
            the standard estimator. Supplying the tilt obtained at a nearby
            population makes numerator and denominator of a ratio share one
            expansion point, which is the Tierney-Kadane arrangement.

    Returns:
        The constant, the saddlepoint term alone and the tilt actually used
    """
    from .ncld import pfqn_gld, pfqn_lldsingle

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()

    if np.sum(N) < 0:
        return PfqnNreResult(lG=-np.inf, G=0.0, lGs=-np.inf)
    if np.sum(N) == 0:
        return PfqnNreResult(lG=0.0, G=1.0, lGs=0.0)

    Nt = int(round(float(np.sum(N))))
    if alpha is None:
        alpha = np.ones((L.shape[0], Nt))
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    if alpha.shape[1] < Nt:
        raise ValueError("pfqn_nre: the load-dependent rate matrix must have "
                         "at least sum(N) columns.")
    # trim so that every rate used downstream is positive
    alpha = alpha[:, :Nt].copy()

    if Z is not None:
        Z = np.atleast_2d(np.asarray(Z, dtype=float))
        if Z.size > 0 and np.sum(Z) > 0:
            L = np.vstack([L, np.sum(Z, axis=0).reshape(1, -1)])
            # the delay is an infinite server station
            alpha = np.vstack([alpha, np.arange(1, Nt + 1, dtype=float).reshape(1, -1)])

    M, R = L.shape
    if M == 1:
        lg_one = float(pfqn_gld(L, N, alpha, options).lG)
        return PfqnNreResult(lG=lg_one, G=float(np.exp(lg_one)), lGs=lg_one)

    # scale demands in [0,1] per class, the residual factor is exact by homogeneity
    Lmax = np.max(L, axis=0)
    Lmax[Lmax <= 0] = 1.0
    L = L / Lmax[np.newaxis, :]
    lGscale = float(np.dot(N, np.log(Lmax)))

    if R == 1:
        # coefficient extraction is the identity in a single class
        lg_one = float(pfqn_lldsingle(L, np.array([Nt], dtype=float),
                                      alpha, options).lG) + lGscale
        return PfqnNreResult(lG=lg_one, G=float(np.exp(lg_one)), lGs=lg_one)

    d = R - 1  # dimension of the quotient torus
    if d > _MAX_DIM:
        raise ValueError("pfqn_nre: pfqn_nre is limited to 8 classes, "
                         "use nrl or clw beyond that.")

    Nd = N[:d]
    cache = {}
    vbase = np.zeros(d)

    def cgfat(v):
        """log of the single-class LLD constant at the class tilt X=[exp(v),1]"""
        X = np.ones(R)
        X[:d] = np.exp(v)
        return float(pfqn_lldsingle(L.dot(X).reshape(-1, 1),
                                    np.array([Nt], dtype=float), alpha, options).lG)

    def cgf(off):
        """cumulant generating function at vbase+off*hstep, memoised on the stencil"""
        key = tuple(off)
        y = cache.get(key)
        if y is None:
            y = cgfat(vbase + np.asarray(off, dtype=float) * _HSTEP)
            cache[key] = y
        return y

    def unitoff(a, sgn):
        off = [0] * d
        off[a] = sgn
        return off

    def offsum(*offs):
        acc = [0] * d
        for off in offs:
            for r in range(d):
                acc[r] += off[r]
        return acc

    def second_diff(a, b):
        return (cgf(offsum(unitoff(a, 1), unitoff(b, 1)))
                - cgf(offsum(unitoff(a, 1), unitoff(b, -1)))
                - cgf(offsum(unitoff(a, -1), unitoff(b, 1)))
                + cgf(offsum(unitoff(a, -1), unitoff(b, -1)))) / (4 * _HSTEP ** 2)

    # ---- saddle point: minimise the convex F(v) = K(v) - Nd*v ----
    # A tilt supplied by the caller is used as given, so that a ratio of two
    # constants can be expanded about one common point rather than two.
    converged = False
    tilt_given = vfix is not None and np.size(vfix) > 0
    if tilt_given:
        vfix = np.asarray(vfix, dtype=float).flatten()
        if vfix.size < d:
            raise ValueError("pfqn_nre: the supplied tilt must have one entry "
                             "per quotient dimension (R-1).")
        vbase = vfix[:d].copy()
        converged = True
    for _ in range(0 if tilt_given else 100):
        cache.clear()
        grad = np.zeros(d)
        hess = np.zeros((d, d))
        for a in range(d):
            grad[a] = (cgf(unitoff(a, 1)) - cgf(unitoff(a, -1))) / (2 * _HSTEP) - Nd[a]
        for a in range(d):
            for b in range(d):
                hess[a, b] = second_diff(a, b)
        step = -np.linalg.solve(hess, grad)
        F0 = cgf([0] * d) - float(np.dot(Nd, vbase))
        tau = 1.0
        while tau > 1e-10:
            vtry = vbase + tau * step
            if cgfat(vtry) - float(np.dot(Nd, vtry)) <= F0:
                break
            tau = tau / 2
        vbase = vbase + tau * step
        # Newton converges to the root of the differenced gradient, whose own
        # O(hstep^2) bias puts any absolute gradient target out of reach
        if np.linalg.norm(tau * step) < 1e-10:
            converged = True
            break
    if not converged:
        print("pfqn_nre warning: the saddle point search did not converge, "
              "the estimate may be inaccurate.")

    # ---- cumulants of the tilted distribution at the saddle ----
    cache.clear()
    K0 = cgf([0] * d)
    Sigma = np.zeros((d, d))
    for a in range(d):
        for b in range(d):
            Sigma[a, b] = second_diff(a, b)
    Sigma = (Sigma + Sigma.T) / 2
    if np.min(np.linalg.eigvalsh(Sigma)) <= 0:
        raise ValueError("pfqn_nre: the tilted covariance is singular, "
                         "a class has no demand at any station.")

    k3 = np.zeros((d, d, d))
    for a in range(d):
        for b in range(d):
            for c in range(d):
                acc = 0.0
                for s in range(8):
                    sg = [1 - 2 * ((s >> j) & 1) for j in range(3)]
                    off = offsum(unitoff(a, sg[0]), unitoff(b, sg[1]), unitoff(c, sg[2]))
                    acc += sg[0] * sg[1] * sg[2] * cgf(off)
                k3[a, b, c] = acc / (8 * _HSTEP ** 3)

    k4 = np.zeros((d, d, d, d))
    for a in range(d):
        for b in range(d):
            for c in range(d):
                for e in range(d):
                    acc = 0.0
                    for s in range(16):
                        sg = [1 - 2 * ((s >> j) & 1) for j in range(4)]
                        off = offsum(unitoff(a, sg[0]), unitoff(b, sg[1]),
                                     unitoff(c, sg[2]), unitoff(e, sg[3]))
                        acc += sg[0] * sg[1] * sg[2] * sg[3] * cgf(off)
                    k4[a, b, c, e] = acc / (16 * _HSTEP ** 4)

    # ---- second-order Edgeworth factor, see the module header ----
    S = np.linalg.inv(Sigma)
    rho4 = float(np.einsum('abce,ab,ce->', k4, S, S))
    u = np.einsum('ab,abc->c', S, k3)
    rhoA = float(u.dot(S).dot(u))
    rhoB = float(np.einsum('ijk,ia,jb,kc,abc->', k3, S, S, S, k3))
    corr = 1 + rho4 / 8 - (3 * rhoA + 2 * rhoB) / 24
    if corr <= 0:
        print("pfqn_nre warning: the Edgeworth correction is non-positive, "
              "falling back on the saddlepoint term.")
        corr = 1.0

    sign, logdet = np.linalg.slogdet(Sigma)
    lGs = (K0 - float(np.dot(Nd, vbase)) - (d / 2) * np.log(2 * np.pi)
           - 0.5 * logdet + lGscale)
    lG = lGs + np.log(corr)
    return PfqnNreResult(lG=float(lG), G=float(np.exp(lG)), lGs=float(lGs),
                         vsad=np.asarray(vbase, dtype=float).copy())
