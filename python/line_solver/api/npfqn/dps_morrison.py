"""
Two-term heavy-usage asymptotic approximation for a closed queueing network with one
infinite-server (think) station and one discriminatory processor-sharing (DPS) station.

Python port of matlab/src/api/npfqn/npfqn_dps_morrison.m, cross-checked against
jar/src/main/java/jline/api/npfqn/Npfqn_dps_morrison.java and
cpp/include/line/api/npfqn/npfqn_dps_morrison.h (identical term for term, including the sigma
solve and the W_m recursion).

Reference: J.A. Morrison, "Asymptotic analysis of a large closed queueing network with
discriminatory processor sharing", Queueing Systems 9 (1991) 191-214.

The network is NOT product-form, so nothing here computes a normalizing constant: the method
expands the GENERATING FUNCTION of the balance equations. The substitution P(n) = <w,n> f(n)
clears the DPS denominator and turns the balance recursion into a linear PDE with affine
coefficients (eq. 2.5); rescaling z = 1 - xi/sqrt(N) and expanding in powers of N^(-1/2) leaves a
degenerate leading operator whose kernel is the functions of the similarity variable eta, and the
solvability condition along its characteristic gives an ODE for the amplitude (eq. 2.20).
RESULT 1 (eq. 4.11) and RESULT 2 (eq. 4.17) are the two-term approximations returned here.

Scaling. Morrison writes K_j = N b_j and lambda_j = N r_j g_j with usage
rho = sum_j b_j/g_j = 1 - a/sqrt(N). N is bookkeeping only and the approximation is invariant to
it, so this routine fixes N = 1: b = N_pop, g = Z/S, r = 1/Z, a = 1 - rho. Accuracy is governed by
the PHYSICAL regime -- large populations with rho near 1. rho > 1 is admissible, being the
saturated regime of appendix A.
"""

from dataclasses import dataclass, field
from typing import List

import numpy as np
from scipy.special import erfcx

__all__ = ['npfqn_dps_morrison', 'NpfqnDpsMorrisonResult']


@dataclass
class NpfqnDpsMorrisonResult:
    """Mean queue lengths, sojourn times and throughputs, with Morrison's constants."""
    Q: np.ndarray = field(default=None)       # mean class-k jobs at the DPS station
    R: np.ndarray = field(default=None)       # mean class-k sojourn time per DPS visit
    X: np.ndarray = field(default=None)       # per-class throughput
    Qlead: np.ndarray = field(default=None)   # leading-order (one-term) queue lengths
    Rlead: np.ndarray = field(default=None)   # leading-order (one-term) sojourn times
    sigma: np.ndarray = field(default=None)   # the vector sigma of eq. (4.12)
    W: np.ndarray = field(default=None)       # W_0..W_4 of eq. (3.23)
    rho: float = 0.0
    a: float = 0.0
    cB: float = 0.0
    cC: float = 0.0
    cD: float = 0.0
    cH: float = 0.0
    cI: float = 0.0
    cJ: float = 0.0
    cK: float = 0.0
    cL: float = 0.0
    cM: float = 0.0
    cQ: float = 0.0
    delta: float = 0.0
    cR: float = 0.0
    cS: float = 0.0
    cU: float = 0.0
    cA: float = 0.0
    cV: float = 0.0


def _wm(cB, cC, cD, y):
    """W_m of eq. (3.23), m = 0..4.

    Substituting z = sigma s with sigma = sqrt(D/(BC)) normalizes the Gaussian to
    W_m(y) = (B/D)^2 sigma^(m+1) I_m(yh), yh = y sqrt(B/(CD)),
    I_m(yh) = int_0^inf s^m exp(-s^2/2 - yh s) ds,
    so I_0 = sqrt(pi/2) erfcx(yh/sqrt(2)) -- the SCALED complementary error function, which is what
    keeps large yh from overflowing -- with I_1 = 1 - yh I_0 and I_m = (m-1) I_{m-2} - yh I_{m-1}.
    That recursion subtracts nearly equal terms once yh is large, so a loss of positivity (the I_m
    are integrals of positive integrands) triggers a quadrature fallback.
    """
    mmax = 4
    sig = np.sqrt(cD / (cB * cC))
    yh = y * np.sqrt(cB / (cC * cD))

    Iv = np.zeros(mmax + 1)
    Iv[0] = np.sqrt(np.pi / 2) * erfcx(yh / np.sqrt(2))
    if not np.isfinite(Iv[0]):
        raise ValueError(
            "The usage is so far above saturation (rho = %g) that the Morrison expansion "
            "overflows. This model is outside the moderately-heavy regime the approximation is "
            "derived for; use SolverFLD, SolverMVA or SolverCTMC." % (1 - y))
    Iv[1] = 1 - yh * Iv[0]
    for m in range(2, mmax + 1):
        Iv[m] = (m - 1) * Iv[m - 2] - yh * Iv[m - 1]
    if np.any(Iv <= 0):
        from scipy.integrate import quad
        for m in range(mmax + 1):
            Iv[m] = quad(lambda s, m=m: s ** m * np.exp(-s * s / 2 - yh * s), 0.0, np.inf,
                         limit=400)[0]
    return (cB / cD) ** 2 * sig ** np.arange(1, mmax + 2) * Iv


def npfqn_dps_morrison(N, Z, S, w) -> NpfqnDpsMorrisonResult:
    """Evaluate Morrison's two-term approximation.

    Parameters
    ----------
    N : array_like
        Per-class populations, finite and positive.
    Z : array_like
        Per-class mean think times, finite and positive.
    S : array_like
        Per-class mean DPS service times, finite and positive.
    w : array_like
        Per-class DPS weights, finite and positive.

    Returns
    -------
    NpfqnDpsMorrisonResult
        Mean queue lengths, sojourn times and throughputs, with the intermediate constants.
    """
    b = np.asarray(N, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float).ravel()
    S = np.asarray(S, dtype=float).ravel()
    w = np.asarray(w, dtype=float).ravel()
    p = b.size
    if Z.size != p or S.size != p or w.size != p:
        raise ValueError("N, Z, S and w must have the same number of classes.")
    if np.any(~np.isfinite(b)) or np.any(b <= 0):
        raise ValueError("The Morrison approximation requires finite positive class populations "
                         "(closed classes only).")
    if np.any(~np.isfinite(Z)) or np.any(Z <= 0) or np.any(~np.isfinite(S)) or np.any(S <= 0):
        raise ValueError("Think times Z and DPS service times S must be finite and positive.")
    if np.any(~np.isfinite(w)) or np.any(w <= 0):
        raise ValueError("DPS weights must be finite and positive.")

    # Morrison's parameters at the bookkeeping scale N = 1
    r = 1.0 / Z
    g = Z / S
    rho = float(np.sum(b / g))
    a = 1.0 - rho

    # constants, eqs. (2.18), (2.19), (3.11), (3.13)-(3.15)
    cB = float(np.sum(b / (r * g ** 2 * w)))
    cC = float(np.sum(b / (g ** 2 * w)))
    cD = float(np.sum(b / (r * g ** 2)))
    cH = float(np.sum(b / (r ** 2 * g ** 3 * w)))
    cI = float(np.sum(b / (r ** 2 * g ** 3 * w ** 2)))
    cJ = float(np.sum(b / (r * g ** 3 * w ** 2)))
    cK = float(np.sum(b / (g ** 3 * w ** 2)))
    cL = float(np.sum(b / (r * g ** 3 * w)))
    cM = float(np.sum(b / (r ** 2 * g ** 3)))
    cQ = float(np.sum(b / g ** 2))

    # sigma: eq. (4.12) with the normalization (4.13). The p equations have rank p-1 (Morrison
    # p.197), so the last one -- implied by the others -- is REPLACED by (4.13), giving a square
    # nonsingular system. All four codebases use this same scheme so their sigma agree.
    A = np.zeros((p, p))
    rhs = np.zeros(p)
    for i in range(p):
        A[i, i] += rho
        for j in range(p):
            den = r[i] * g[i] * w[i] + r[j] * g[j] * w[j]
            A[i, i] -= w[j] * b[j] * r[j] / den
            A[i, j] -= w[j] * b[i] * r[i] / den
        rhs[i] = rho * (b[i] / g[i]) * (cD / (cB * w[i]) - 1.0)
    A[p - 1, :] = 1.0 / (r * g)
    rhs[p - 1] = 0.0
    sigma = np.linalg.solve(A, rhs)

    # alpha from eq. (4.9), then delta of eq. (3.15)
    alpha = sigma - (b / g) * (cD / (cB * w) - 1.0)
    delta = float(np.sum(alpha / g))

    # eqs. (3.19)-(3.21)
    cR = 3.0 * (cB * cL - cD * cJ) / (cB * cD)
    cS = (2 * cB * (cD * cH - cB * cM) - cD * (cD * cI - cB * cH)) / (2 * cB ** 2 * cD ** 2)
    cU = (cQ - cC * cD / cB - delta) / rho - cD * cR / cB + (a ** 2 - cC * cD / cB) * cS
    cA = cS * cC ** 2 + cR * cC - cK
    cV = cR + 2 * cS * cC

    W = _wm(cB, cC, cD, a)
    W0, W1, W2, W3, W4 = W

    # RESULT 1 (4.11) and RESULT 2 (4.17), at sqrt(N) = 1. NOTE the numerator bracket carries
    # cU*W2: eq. (4.10) of the paper misprints it as U*W1, but (4.7), (4.11), (A6) and (B2) all
    # agree on U*W2, and it is what the derivation from (4.4)-(4.9) gives.
    eps = cB / cD
    num = W1 - eps * (cA / 3 * W4 + a / 2 * cV * W3 + cU * W2)
    den = W0 - eps * (cA / 3 * W3 + a / 2 * cV * W2 + cU * W1 + cS)
    if den == 0 or not np.isfinite(den):
        raise ValueError("The Morrison expansion is degenerate for this model (vanishing "
                         "denominator); the usage is too far from the moderately-heavy regime.")

    Qlead = b * W1 / (g * w * W0)
    Q = b * num / (g * w * den) - b * W2 / (g ** 2 * w ** 2 * W0) - sigma / rho
    Rlead = W1 / (r * g * w * W0)
    R = (num / (r * g * w * den)
         + ((W1 / W0) ** 2 - W2 / W0) / (r * g ** 2 * w ** 2)
         - sigma / (rho * r * b))
    X = r * (b - Q)

    return NpfqnDpsMorrisonResult(
        Q=Q, R=R, X=X, Qlead=Qlead, Rlead=Rlead, sigma=sigma, W=W,
        rho=rho, a=a, cB=cB, cC=cC, cD=cD, cH=cH, cI=cI, cJ=cJ, cK=cK, cL=cL, cM=cM, cQ=cQ,
        delta=delta, cR=cR, cS=cS, cU=cU, cA=cA, cV=cV)
