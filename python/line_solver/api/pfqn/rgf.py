"""Recursion by Generating Functions (RGF) for product-form normalizing constants.

Single class: Coury and Harrison (1997), "Asymptotic properties of queuing
networks", IEE Proc.-Comput. Digit. Tech. 144(5):247-254, Property 1.

Multiclass: Harrison and Coury (2002), "On the asymptotic behaviour of closed
multiclass queueing networks", Perf. Eval. 47:131-138, Thm 1, as algorithmised
by Harrison and Lee (2004), "A new recursive algorithm for computing generating
functions in closed multi-class queueing networks", IEEE MASCOTS, eqs. (4)-(5).
Neither carries an infinite server: their generating function is the rational
prod_i (1 - rho_i z)^-m_i, and an entire numerator breaks the residues-sum-to-
zero identity the recursion rests on. Think times are therefore carried by the
truncation of Bertozzi and McKenna (1993), "Multidimensional residues,
generating functions, and their application to queueing networks", SIAM Review
35(2):239-268, eqs. (3.19)-(3.21).
"""

from typing import Optional, Tuple

import numpy as np
from scipy.special import gammaln

__all__ = ['pfqn_rgf', 'pfqn_rgfmc']


def _logconv(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Log-domain linear convolution truncated at the common length."""
    n = u.size
    c = np.full(n, -np.inf)
    for k in range(n):
        t = u[:k + 1] + v[k::-1]
        vm = t.max()
        if np.isinf(vm):
            c[k] = vm
        else:
            c[k] = vm + np.log(np.exp(t - vm).sum())
    return c


def pfqn_rgf(L, N, Z: float = 0.0) -> Tuple[float, float, np.ndarray]:
    """Exact normalizing constant of a single-class closed product-form network
    by convolving per-node generating-function sequences.

    A GROUP of m stations sharing the same demand p collapses into the single
    negative-binomial sequence r(k) = C(k+m-1,k) p^k, so the whole group costs
    one sequence rather than m convolution passes; the delay contributes the
    Poisson sequence Z^k/k!. Cost O(G N^2) against Buzen's O(M N), so RGF is
    the cheaper route on heavily replicated models with moderate populations
    (G N < M). The recursion runs entirely in the log domain.

    Args:
        L: Service demand vector (M,) of the queueing stations.
        N: Population (nonnegative integer scalar).
        Z: Think time (scalar, default 0).

    Returns:
        Tuple (G, lG, lg) with lg the vector of log g(0), ..., log g(N).
    """
    L = np.asarray(L, dtype=float).ravel()
    N = np.asarray(N, dtype=float).ravel()
    if N.size > 1:
        raise ValueError('pfqn_rgf is a single-class method, but the population '
                         'vector has more than one entry.')
    n = float(N[0])
    Z = float(np.sum(np.asarray(Z, dtype=float)))
    if n < 0 or abs(n - round(n)) > 0:
        raise ValueError('pfqn_rgf requires a nonnegative integer population.')
    n = int(round(n))
    if np.any(L < 0) or Z < 0:
        raise ValueError('pfqn_rgf requires nonnegative demands and think time.')

    kk = np.arange(n + 1, dtype=float)
    lg = np.full(n + 1, -np.inf)
    lg[0] = 0.0

    if Z > 0:
        lg = _logconv(lg, kk * np.log(Z) - gammaln(kk + 1))

    Lp = L[L > 0]
    if Lp.size:
        p, mult = np.unique(Lp, return_counts=True)
        for i in range(p.size):
            m = int(mult[i])
            if m == 1:
                lr = kk * np.log(p[i])
            else:
                lr = (gammaln(kk + m) - gammaln(kk + 1) - gammaln(m)
                      + kk * np.log(p[i]))
            lg = _logconv(lg, lr)

    lG = float(lg[-1])
    return float(np.exp(lG)), lG, lg


# --------------------------------------------------------------------------
# Multiclass RGF: iterated residues with think times.
# --------------------------------------------------------------------------

def _slogsum(lv, sv, cond=None):
    """Sum sv[i]*exp(lv[i]) in the log domain; return (log|sum|, sign).

    When ``cond`` is a one-element list it collects the worst CANCELLATION
    RATIO log(sum|t|) - log|sum t| seen. The elimination is exact in exact
    arithmetic, so that ratio is the only thing between a correct lG and a
    confidently wrong one, and pfqn_rgfmc refuses on it rather than guessing.
    """
    lv = np.asarray(lv, dtype=float)
    sv = np.asarray(sv, dtype=float)
    k = np.isfinite(lv) & (sv != 0.0)
    if not k.any():
        return -np.inf, 0.0
    lv = lv[k]
    sv = sv[k]
    mx = lv.max()
    e = np.exp(lv - mx)
    tot = float(np.sum(sv * e))
    if tot == 0.0:
        return -np.inf, 0.0
    lsum = mx + np.log(abs(tot))
    if cond is not None and lv.size > 1:
        labs = mx + np.log(float(np.sum(e)))
        if labs - lsum > cond[0]:
            cond[0] = labs - lsum
    return lsum, (1.0 if tot > 0 else -1.0)


def _slogconv(lu, su, lv, sv, cond=None):
    """Signed log-domain linear convolution, truncated at the common length."""
    n = lu.size
    lc = np.full(n, -np.inf)
    sc = np.zeros(n)
    for k in range(n):
        lc[k], sc[k] = _slogsum(lu[:k + 1] + lv[k::-1], su[:k + 1] * sv[k::-1],
                                cond)
    return lc, sc


def _lbinom(n, r):
    """log C(n,r) for real n >= r >= 0, never a factorial quotient."""
    return gammaln(n + 1.0) - gammaln(r + 1.0) - gammaln(n - r + 1.0)


def _compositions(total, parts):
    """Nonnegative integer vectors of length `parts` summing to `total`."""
    if parts == 0:
        if total == 0:
            yield ()
        return
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, parts - 1):
            yield (first,) + rest


class _Term(object):
    """coef * prod_t ( F[t,0] + sum_{p>=1} F[t,p] z_p ) ** (-m[t])."""

    __slots__ = ('lc', 'sc', 'F', 'm')

    def __init__(self, lc, sc, F, m):
        self.lc = lc
        self.sc = sc
        self.F = np.asarray(F, dtype=float)
        self.m = np.asarray(m, dtype=int)


def _merge(F, m, tol):
    """Fuse proportional affine forms into one factor of summed multiplicity.

    Harrison-Coury Thm 1 assumes rho_iq != rho_lq. When two forms are
    PROPORTIONAL they name the same pole, C_jl vanishes and the partial
    fraction is undefined; merging is what makes the degenerate case the paper
    leaves open in its Conclusion go through with nothing else changed.
    """
    T = F.shape[0]
    used = np.zeros(T, dtype=bool)
    outF = []
    outm = []
    lc = 0.0
    sc = 1.0
    for i in range(T):
        if used[i]:
            continue
        used[i] = True
        Fi = F[i].copy()
        mi = int(m[i])
        pi = int(np.argmax(np.abs(Fi)))
        if Fi[pi] == 0.0:
            raise ValueError('pfqn_rgfmc met an identically zero factor.')
        for j in range(i + 1, T):
            if used[j]:
                continue
            Fj = F[j]
            nj = float(np.max(np.abs(Fj)))
            if nj <= 0.0:
                continue
            r = Fj[pi] / Fi[pi]
            if r == 0.0:
                continue
            if np.all(np.abs(Fj - r * Fi) <= tol * nj):
                lc -= m[j] * np.log(abs(r))
                sc *= 1.0 if r > 0 else (-1.0) ** int(m[j])
                mi += int(m[j])
                used[j] = True
        outF.append(Fi)
        outm.append(mi)
    return np.array(outF, dtype=float), np.array(outm, dtype=int), lc, sc


def _rgf_step(terms, col, kr, Zr, tol, maxterms):
    """Eliminate the class in column `col` of F, leaving columns 0..col-1.

    This is Harrison-Coury Thm 1 written as a partial fraction: the poles of
    z_col sit at A_j/B_j with order m_j, and the residue there re-expresses the
    network as one with a class fewer. The delay rides along as the Poisson
    weight of the truncated exp(Z z), which is EXACT here because coefficients
    above the target degree cannot reach [z_col^Ntot].
    """
    out = []
    for t in terms:
        F, m, dlc, dsc = _merge(t.F, t.m, tol)
        lc = t.lc + dlc
        sc = t.sc * dsc
        A = F[:, :col]
        B = -F[:, col]
        scale = np.max(np.abs(F), axis=1)
        scale[scale == 0.0] = 1.0
        isMono = np.all(np.abs(A) <= tol * scale[:, None], axis=1)
        shift = 0
        if isMono.any():
            lc -= float(np.sum(m[isMono] * np.log(np.abs(B[isMono]))))
            sc *= float(np.prod(np.where(-B[isMono] > 0, 1.0, -1.0)
                                ** m[isMono]))
            shift = int(np.sum(m[isMono]))
            keep = ~isMono
            A = A[keep]
            B = B[keep]
            m = m[keep]
        S = np.flatnonzero(B != 0.0)
        P = np.flatnonzero(B == 0.0)
        Ntot = kr + shift
        plist = range(Ntot + 1) if Zr > 0 else (0,)
        for p in plist:
            lcz = lc + (p * np.log(Zr) - gammaln(p + 1.0) if p > 0 else 0.0)
            n = Ntot - p
            if S.size == 0:
                if n == 0:
                    out.append(_Term(lcz, sc, A[P], m[P]))
                continue
            for j in S:
                oth = S[S != j]
                no = oth.size
                Bj = B[j]
                Aj = A[j]
                Cjl = ((A[oth] * Bj - B[oth][:, None] * Aj[None, :]) / Bj
                       if no else np.zeros((0, col)))
                for k in range(int(m[j])):
                    lbase = (lcz - k * np.log(abs(Bj))
                             + _lbinom(n + m[j] - k - 1, n)
                             + n * np.log(abs(Bj)))
                    sbase = (sc * ((-1.0 if -Bj < 0 else 1.0) ** k)
                             * ((-1.0 if Bj < 0 else 1.0) ** n))
                    for jl in _compositions(k, no):
                        lt = lbase
                        st = sbase
                        for a in range(no):
                            if jl[a]:
                                lt += (_lbinom(m[oth[a]] + jl[a] - 1, jl[a])
                                       + jl[a] * np.log(abs(B[oth[a]])))
                                st *= ((-1.0 if B[oth[a]] < 0 else 1.0)
                                       ** int(jl[a]))
                        newF = ([Aj] + [Cjl[a] for a in range(no)]
                                + [A[q] for q in P])
                        newm = ([n + int(m[j]) - k]
                                + [int(m[oth[a]]) + int(jl[a])
                                   for a in range(no)]
                                + [int(m[q]) for q in P])
                        out.append(_Term(lt, st, np.array(newF, dtype=float),
                                         np.array(newm, dtype=int)))
        if len(out) > maxterms:
            raise ValueError('pfqn_rgfmc exceeded maxterms=%d. The residue '
                             'term count grows as C(S+M-1,M-1) per further '
                             'elimination; use method "ca".' % maxterms)
    return out


def _rgf_base_kernel(loads, mults, N, Z, cond=None):
    """[z^N] exp(Z z) prod_t (1 - p_t z)^-m_t, signed and in the log domain.

    Coury-Harrison Property 1 / Harrison-Lee sec. 3.3.1, with the loads allowed
    to be negative: an eliminated class leaves pole differences of either sign.
    """
    kk = np.arange(N + 1, dtype=float)
    lg = np.full(N + 1, -np.inf)
    sg = np.zeros(N + 1)
    lg[0] = 0.0
    sg[0] = 1.0
    if Z > 0:
        lg, sg = _slogconv(lg, sg, kk * np.log(Z) - gammaln(kk + 1.0),
                           np.ones(N + 1), cond)
    for a in range(loads.size):
        p = loads[a]
        mm = int(mults[a])
        if p == 0:
            continue
        if mm == 1:
            lr = kk * np.log(abs(p))
        else:
            lr = (gammaln(kk + mm) - gammaln(kk + 1.0) - gammaln(float(mm))
                  + kk * np.log(abs(p)))
        sr = np.ones(N + 1) if p > 0 else (-1.0) ** kk
        lg, sg = _slogconv(lg, sg, lr, sr, cond)
    return lg[N], sg[N]


def _rgf_base(terms, k1, Z1, tol, cache, cond):
    """Single-class base case, memoised per Harrison-Lee sec. 3.4."""
    lv = []
    sv = []
    for t in terms:
        F, m, dlc, dsc = _merge(t.F, t.m, tol)
        lc = t.lc + dlc
        sc = t.sc * dsc
        shift = 0
        loads = []
        mults = []
        for a in range(F.shape[0]):
            a0, a1 = F[a, 0], F[a, 1]
            ma = int(m[a])
            sca = max(abs(a0), abs(a1))
            if sca == 0.0:
                raise ValueError('pfqn_rgfmc met an identically zero factor.')
            if abs(a0) <= tol * sca:
                lc -= ma * np.log(abs(a1))
                sc *= (1.0 if a1 > 0 else -1.0) ** ma
                shift += ma
            else:
                lc -= ma * np.log(abs(a0))
                sc *= (1.0 if a0 > 0 else -1.0) ** ma
                if abs(a1) > tol * sca:
                    loads.append(-a1 / a0)
                    mults.append(ma)
        Ntot = k1 + shift
        key = (Ntot, round(Z1, 15),
               tuple(sorted(zip([round(x, 15) for x in loads], mults))))
        hit = cache.get(key)
        if hit is None:
            hit = _rgf_base_kernel(np.array(loads, dtype=float),
                                   np.array(mults, dtype=int), Ntot, Z1, cond)
            cache[key] = hit
        lg, sg = hit
        if sg != 0.0:
            lv.append(lc + lg)
            sv.append(sc * sg)
    return _slogsum(lv, sv, cond)


def pfqn_rgfmc(L, N, Z=None, tol: float = 1e-12, maxterms: int = 1000000,
               maxcancel: float = 15.0) -> Tuple[float, float]:
    """Exact multiclass normalizing constant by recursion on generating functions.

    Eliminates one class at a time by residues (Harrison-Coury 2002 Thm 1,
    algorithmised as Harrison-Lee 2004 eqs. 4-5) until a single class is left,
    then finishes with the multiplicity-aware convolution of Coury-Harrison
    1997 Property 1, memoised by load vector as Harrison-Lee sec. 3.4.

    THINK TIMES ARE NOT IN EITHER PAPER. Their generating function is the
    rational prod_i (1 - rho_i z)^-m_i, and an infinite server multiplies it by
    the entire exp(sum_r Z_r z_r), which destroys the residues-sum-to-zero
    identity (Bertozzi-McKenna 1993 fact IV, p. 246) the recursion rests on.
    The delay is therefore carried by their truncation, eqs. (3.19)-(3.21):
    only the first k_r+1 Taylor coefficients of exp(Z_r z_r) can reach
    [z_r^k_r], so replacing it by that polynomial is EXACT, not an
    approximation, and restores a rational integrand. The cost is that the
    eliminated class's population re-enters the term count, which is precisely
    the population-insensitivity Harrison-Lee sec. 4 advertises; the class kept
    for the base case pays nothing, so the base is chosen to be a populous one
    only when conditioning allows.

    The recursion is exact in exact arithmetic but is an ALTERNATING sum over
    residues, so near-coincident loads over the eliminated class destroy
    significance. `maxcancel` bounds the nats of cancellation tolerated and the
    routine REFUSES beyond it rather than returning a confidently wrong lG.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,), nonnegative integers.
        Z: Think time vector (R,), default zeros.
        tol: relative tolerance for calling two affine forms proportional.
        maxterms: cap on residue terms carried between eliminations.
        maxcancel: nats of cancellation tolerated before refusing.

    Returns:
        Tuple (G, lG).
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    R = L.shape[1]
    Z = np.zeros(R) if Z is None else np.asarray(Z, dtype=float).ravel()
    if N.size != R or Z.size != R:
        raise ValueError('pfqn_rgfmc requires numel(N) and numel(Z) to match '
                         'the number of columns of L.')
    if np.any(L < 0) or np.any(Z < 0) or np.any(N < 0):
        raise ValueError('pfqn_rgfmc requires nonnegative L, N and Z.')
    if np.any(np.abs(N - np.round(N)) > 0):
        raise ValueError('pfqn_rgfmc requires integer populations.')
    keep = N > 0
    L = L[:, keep]
    N = N[keep]
    Z = Z[keep]
    R = int(keep.sum())
    if R == 0:
        return 1.0, 0.0
    L = L[np.any(L > 0, axis=1)]
    M = L.shape[0]
    if M == 0:
        lG = float(np.sum(N * np.log(Z) - gammaln(N + 1.0)))
        return float(np.exp(lG)), lG
    if R == 1:
        return pfqn_rgf(L[:, 0], float(N[0]), float(Z[0]))[:2]
    # Base class = smallest population: that population is the degree the
    # sign-indefinite base series is carried to, so it drives the cancellation.
    order = np.argsort(N, kind='stable')
    L = L[:, order]
    N = N[order]
    Z = Z[order]
    cs = np.maximum(np.max(L, axis=0), Z)
    cs[cs <= 0] = 1.0
    L = L / cs
    Z = Z / cs
    lGscale = float(np.sum(N * np.log(cs)))
    terms = [_Term(0.0, 1.0, np.hstack([np.ones((M, 1)), -L]),
                   np.ones(M, dtype=int))]
    cond = [0.0]
    # F column p carries class p-1, so eliminating class p-1 means column p.
    for col in range(R, 1, -1):
        terms = _rgf_step(terms, col, int(round(N[col - 1])),
                          float(Z[col - 1]), tol, maxterms)
        if not terms:
            return 0.0, -np.inf
    lG, sG = _rgf_base(terms, int(round(N[0])), float(Z[0]), tol, {}, cond)
    if sG == 0.0:
        return 0.0, -np.inf
    if sG < 0 or cond[0] > maxcancel:
        raise ValueError(
            'pfqn_rgfmc: the residue sum cancelled %.1f nats, past the %.1f '
            'allowed, so lG carries no significant digits. The eliminated '
            'classes have near-coincident loads over the stations. Use method '
            '"ca" for the exact convolution.' % (cond[0], maxcancel))
    lG += lGscale
    return float(np.exp(lG)), lG
