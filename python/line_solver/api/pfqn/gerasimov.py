"""Gerasimov's residue (closed-form) normalizing constant, generalized to R classes.

A. I. Gerasimov, "On Normalizing Constants in Multiclass Queueing Networks",
Operations Research 43(4):704-711, 1995. Ported at parity from
matlab/src/api/pfqn/pfqn_gerasimov.m.
"""

import math
from typing import List, Tuple

import numpy as np
from scipy.special import gammaln

__all__ = ['pfqn_gerasimov']

# largest log-scale value whose exponential is still a finite double
_LOG_DBL_MAX = math.log(np.finfo(float).max)


def _binom(n: int, k: int) -> float:
    """Binomial coefficient, exact while the value stays representable."""
    if k < 0 or n < 0 or k > n:
        return 0.0
    k = min(k, n - k)
    b = 1.0
    for i in range(1, k + 1):
        b = b * (n - k + i) / i
    if b < 2.0 ** 53:
        b = float(round(b))
    return b


def _compositions(n: int, k: int) -> List[Tuple[int, ...]]:
    """All k-tuples of nonnegative integers summing to n."""
    if k == 0:
        return [()] if n == 0 else []
    if k == 1:
        return [(n,)]
    out = []
    for a in range(n + 1):
        for sub in _compositions(n - a, k - 1):
            out.append((a,) + sub)
    return out


def _merge(F: np.ndarray, m: np.ndarray, c: float, tol: float):
    """Merge proportional affine forms: f_k = lambda f_j is one pole of order
    m_j+m_k, not two nearby simple ones, and lambda^-m_k moves into the scalar."""
    nf = F.shape[0]
    if nf <= 1:
        return F, m, c
    keep = np.ones(nf, dtype=bool)
    for j in range(nf):
        if not keep[j]:
            continue
        pj = int(np.argmax(np.abs(F[j, :])))
        if F[j, pj] == 0.0:
            continue
        for k in range(j + 1, nf):
            if not keep[k]:
                continue
            lam = F[k, pj] / F[j, pj]
            if lam == 0.0:
                continue
            scale = max(np.max(np.abs(F[k, :])), np.max(np.abs(F[j, :])))
            if np.max(np.abs(F[k, :] - lam * F[j, :])) <= tol * scale:
                c = c * lam ** (-float(m[k]))
                m[j] = m[j] + m[k]
                keep[k] = False
    return F[keep, :], m[keep], c


def _step(terms, r: int, Nr: int, Zr: float, tol: float, maxterms: int):
    """One residue elimination: integrate out u_r and return the surviving sum of
    products of affine powers, each form narrowed from r+1 to r columns."""
    out = []
    for (c0, F0, m0) in terms:
        F, m, c = _merge(F0.copy(), m0.copy(), c0, tol)
        A = F[:, :r]           # affine part in u_1..u_(r-1), column 0 = constant
        B = -F[:, r]           # f_j = A_j - B_j u_r
        scale = np.max(np.abs(F), axis=1)
        is_mono = np.all(np.abs(A) <= tol * scale[:, None], axis=1)
        if np.any(is_mono & (np.abs(B) <= tol * scale)):
            raise ValueError('pfqn_gerasimov met an identically zero factor, which '
                             'cannot happen after merging proportional ones.')
        # A factor -B u_r carries no finite pole: it only shifts the exponent.
        shift = 0
        if np.any(is_mono):
            c = c * float(np.prod((-B[is_mono]) ** (-m[is_mono].astype(float))))
            shift = int(np.sum(m[is_mono]))
            keep = ~is_mono
            A, B, m = A[keep, :], B[keep], m[keep]
        S = np.flatnonzero(B != 0.0)    # factors carrying a pole in u_r
        P = np.flatnonzero(B == 0.0)    # factors free of u_r, carried through
        Ntot = Nr + shift
        plist = range(Ntot + 1) if Zr > 0 else (0,)
        for p in plist:
            # Poisson weight Z_r^p/p! through logs: the naive ratio overflows
            # for p >~ 171, reachable when the eliminated class has think time.
            cz = (c * math.exp(p * math.log(Zr) - math.lgamma(p + 1.0))
                  if p > 0 else c)
            Neff = Ntot - p
            if S.size == 0:
                if Neff == 0:
                    out.append((cz, A.copy(), m.copy()))
                continue
            for j in S:
                oth = S[S != j]
                no = oth.size
                if no > 0:
                    Cjl = (A[oth, :] * B[j] - np.outer(B[oth], A[j, :])) / B[j]
                else:
                    Cjl = np.zeros((0, r))
                for k in range(int(m[j])):
                    for nk in _compositions(k, no):
                        coef = cz * (-B[j]) ** (-k) \
                            * _binom(Neff + int(m[j]) - k - 1, Neff) * B[j] ** Neff
                        for li in range(no):
                            coef *= _binom(int(m[oth[li]]) + nk[li] - 1, nk[li]) \
                                * B[oth[li]] ** nk[li]
                        if coef == 0.0:
                            continue
                        Fn = np.vstack((A[j, :][None, :], Cjl, A[P, :]))
                        if no > 0:
                            moth = m[oth] + np.asarray(nk, dtype=int)
                        else:
                            moth = np.zeros(0, dtype=int)
                        mn = np.concatenate((np.array([Neff + int(m[j]) - k], dtype=int),
                                             moth, m[P]))
                        out.append((coef, Fn, mn))
        if len(out) > maxterms:
            raise ValueError('pfqn_gerasimov exceeded maxterms (%d) while eliminating '
                             'class %d; the residue expansion of this model is too '
                             'large. Use pfqn_ca or pfqn_nc.' % (maxterms, r))
    return out


def _base(terms, N1: int, Z1: float, tol: float) -> float:
    """Last class: a univariate coefficient extraction. Summing the residues here too
    would repeat the step above, but convolving the series of each factor returns the
    same number without expanding the multiple poles."""
    G = 0.0
    for (c0, F0, m0) in terms:
        F, m, c = _merge(F0.copy(), m0.copy(), c0, tol)
        A = F[:, 0].copy()
        B = -F[:, 1].copy()
        scale = np.max(np.abs(F), axis=1)
        is_mono = np.abs(A) <= tol * scale
        if np.any(is_mono & (np.abs(B) <= tol * scale)):
            raise ValueError('pfqn_gerasimov met an identically zero factor at the '
                             'innermost coefficient extraction.')
        shift = 0
        if np.any(is_mono):
            c = c * float(np.prod((-B[is_mono]) ** (-m[is_mono].astype(float))))
            shift = int(np.sum(m[is_mono]))
            keep = ~is_mono
            A, B, m = A[keep], B[keep], m[keep]
        Ntot = N1 + shift
        s = np.zeros(Ntot + 1)
        s[0] = 1.0
        if Z1 > 0:
            _n = np.arange(Ntot + 1, dtype=float)
            pois = np.exp(_n * np.log(Z1) - gammaln(_n + 1.0))
            s = np.convolve(s, pois)[:Ntot + 1]
        for j in range(A.size):
            c = c * A[j] ** (-float(m[j]))
            if B[j] == 0.0:
                continue
            ratio = B[j] / A[j]
            seq = np.array([_binom(int(m[j]) + n - 1, n) * ratio ** n for n in range(Ntot + 1)])
            s = np.convolve(s, seq)[:Ntot + 1]
        G += c * s[Ntot]
    return G


def pfqn_gerasimov(L, N, Z=None, tol: float = 1e-12,
                   maxterms: int = 200000) -> Tuple[float, float]:
    """Exact normalizing constant of a closed multiclass product-form network by
    ITERATED RESIDUES of its rational generating function, one class at a time.

    Gerasimov (1995) evaluates

        G(N_1,...,N_R) = (2 pi i)^-R int_G1 ... int_GR
                           prod_s z_s^(N_s-1) prod_i (1 - sum_s x_is/z_s)^-1

    by residues, and gives the resulting CLOSED FORM only for R = 1 (Thm 1-2) and
    R = 2 (Thm 3 for simple poles, Thm 4 for multiple ones), stating that "for three
    or more classes of customers, the normalizing constants can be found by numerical
    methods". This routine implements the residue elimination itself, so the closed
    form is produced for ANY R; at R = 2 it reproduces Thm 3/4 term by term.

    Written as a coefficient of the u_s = 1/z_s series,

        G(N) = [prod_s u_s^(N_s)] exp(sum_s Z_s u_s) prod_i (1 - sum_s x_is u_s)^-1,

    every factor is AFFINE in u, so singling out u_r gives f = A - B u_r with A affine
    in the surviving variables. Partial fractions in u_r map a sum of products of
    affine powers into another one with one variable fewer, and R-1 such steps leave a
    univariate coefficient extraction. At R = 2 the single step returns one term per
    station i, with outer factor x_i2^(N_2+M-1) / prod_{k!=i}(x_i2-x_k2), a pole of
    order N_2+1 at x_i1 and simple poles at the paper's
    z_1ik = (x_k1 x_i2 - x_i1 x_k2)/(x_i2 - x_k2): exactly Thm 3, with the multiple
    poles of Thm 4 (his xi_i < M) handled by the same step. Tied x_i2, vanishing x_i2
    and identical station rows, all outside the paper's hypotheses, are ordinary cases
    here.

    Cost. Let M be the number of stations and order the populations
    N_(1) <= ... <= N_(R). The first elimination turns the single input term into M,
    and every later one multiplies the count by C(S+M-1,M-1) + M-1, where S is the
    total population already eliminated: a pole of order S+1 has to be differentiated
    against the M-1 remaining ones. The innermost extraction then convolves M series
    of length N_(1). Hence R = 1 costs O(M N), Buzen's own cost; R = 2 costs
    O(M^2 N_(1)^2), INDEPENDENT OF N_(2); and R >= 3 costs the same times
    prod_{r=3}^{R} C(N_(r)+M-1, M-1). The R = 2 line is the reason to reach for this
    method: a population removed by residues enters only as a pole ORDER, i.e. through
    binomial coefficients, so it costs nothing at all. On a 4-station two-class model
    at N = [6, 20000] this returns lG in 0.4 ms where pfqn_ca needs 1.8 s, to the same
    1.3e-16. For R >= 3 the term count is polynomial in the populations of degree
    (M-1)(R-2) and exponential in R, which is why the paper stops at two classes and
    why maxterms exists.

    Conditioning. The sum is alternating, exactly as the paper writes it, and two
    decisions keep it usable: near-coincident poles are merged under a RELATIVE
    tolerance, so they are one multiple pole rather than two nearly cancelling simple
    ones, and the class left for the innermost extraction is the one with the SMALLEST
    population, because that population is the degree the final, sign-indefinite series
    is carried to. Measured on 372 random models against pfqn_ca: median 2.0e-16, p90
    4.6e-15, p99 1.1e-11, worst 1.5e-09. On an ill-conditioned demand matrix pfqn_ca or
    pfqn_nc are still the safer routes to the same number.

    Args:
        L: service demand matrix (M x R), L[i,r] = demand of class r at station i.
        N: population vector (R,), nonnegative integers.
        Z: think time vector (R,), default zeros. A delay contributes the entire
           factor exp(sum_s Z_s u_s), handled exactly by convolving its Poisson
           coefficients into each elimination.
        tol: relative tolerance for declaring two affine forms proportional, hence
             one pole rather than two.
        maxterms: cap on the number of residue terms carried between eliminations.
             Exceeding it is an error, not a truncation: a truncated residue sum is
             not a bound or an approximation of G, it is a wrong number.

    Returns:
        (G, lG) the normalizing constant and its logarithm.
    """
    L = np.asarray(L, dtype=float)
    if L.ndim == 1:
        L = L.reshape(-1, 1) if np.asarray(N).size == 1 else L.reshape(1, -1)
    R = L.shape[1]
    N = np.asarray(N, dtype=float).reshape(-1)
    Z = np.zeros(R) if Z is None else np.asarray(Z, dtype=float).reshape(-1)
    if N.size != R:
        raise ValueError('pfqn_gerasimov requires len(N) to match the number of columns of L.')
    if Z.size != R:
        raise ValueError('pfqn_gerasimov requires len(Z) to match the number of columns of L.')
    if np.any(L < 0) or np.any(Z < 0) or np.any(N < 0):
        raise ValueError('pfqn_gerasimov requires nonnegative L, N and Z.')
    if np.any(N != np.round(N)):
        raise ValueError('pfqn_gerasimov requires integer populations.')
    N = np.round(N).astype(int)

    # A class with no jobs is eliminated by evaluating the generating function at
    # u_r = 0, i.e. by deleting its column outright.
    keepr = N > 0
    L, N, Z = L[:, keepr], N[keepr], Z[keepr]
    R = N.size
    if R == 0:
        return 1.0, 0.0
    # A station with no demand at all contributes the factor 1.
    L = L[np.any(L > 0, axis=1), :]

    # Per-class scaling. The residue coefficients carry x_ir^(N_r+M-1), which in
    # double overflows well before G itself does: at x = 4 and N_r = 400 the factor
    # alone is 1e240 while G is finite. Dividing column r by c_r divides G by exactly
    # c_r^N_r (substitute u_r -> u_r/c_r in the generating function), so the scaling
    # is exact and is undone in the log domain at the end.
    cs = np.maximum(L.max(axis=0) if L.shape[0] > 0 else np.zeros(R), Z)
    cs[cs <= 0] = 1.0
    L = L / cs[None, :]
    Z = Z / cs
    lGscale = float(np.sum(N * np.log(cs)))

    # Class order, which decides both the cost and the accuracy.
    #  - The class left for the innermost extraction sets the CONDITIONING. Its
    #    population is the degree the final series is carried to, and the poles of the
    #    reduced problem have arbitrary sign, so that series cancels; the populations
    #    eliminated by residues enter only as pole ORDERS, through binomial
    #    coefficients, and cancel nothing. Basing on N = 100 rather than on N = 6 in
    #    one 4-station model cost 39 nats of lG. The SMALLEST population goes to base.
    #  - Eliminating class r leaves a pole of order N_r+1 that every LATER elimination
    #    has to differentiate, so the rest are eliminated smallest-first to keep the
    #    multiplicities low for as long as possible.
    # Eliminations run from index R down to 2, so indices 2..R hold the remaining
    # populations in DECREASING order and index 1 holds the smallest.
    asc = np.argsort(N, kind='stable')
    ordr = np.concatenate((asc[:1], asc[1:][::-1]))
    L, N, Z = L[:, ordr], N[ordr], Z[ordr]

    M = L.shape[0]
    terms = [(1.0, np.hstack((np.ones((M, 1)), -L)), np.ones(M, dtype=int))]
    # class s (1-based) sits in column s of F = [1, -L]; eliminate R, R-1, ..., 2
    for r in range(R, 1, -1):
        terms = _step(terms, r, int(N[r - 1]), float(Z[r - 1]), tol, maxterms)
        if not terms:
            return 0.0, -np.inf
    Gs = _base(terms, int(N[0]), float(Z[0]), tol)
    if Gs <= 0:
        return 0.0, -np.inf
    lG = float(np.log(Gs) + lGscale)
    # Undoing the scaling can leave the double range; lG is then the only usable
    # form, so report the linear-scale G as infinite rather than overflowing exp.
    G = float('inf') if lG > _LOG_DBL_MAX else float(np.exp(lG))
    return G, lG
