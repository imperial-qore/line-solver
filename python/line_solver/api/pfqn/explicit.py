"""Explicit closed-form normalizing constant of a multiclass closed network.

Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic
Methods", ACM SIGMETRICS 2017, Eqs. (15) and (16). Twin of MATLAB
pfqn_explicit.m.
"""

from typing import Optional, Tuple

import numpy as np
from scipy.special import gammaln

from .utils import multichoose

__all__ = ['pfqn_explicit']


def _factln(x):
    return gammaln(np.asarray(x, dtype=float) + 1.0)


def _nchoosekln(n: float, m: float) -> float:
    return gammaln(1.0 + n) - gammaln(1.0 + n - m) - gammaln(1.0 + m)


def _signed_logsumexp(lterm: np.ndarray, sterm: np.ndarray) -> Tuple[float, float, float]:
    """Signed log-sum-exp of S = sum_i sterm(i)*exp(lterm(i)).

    Returns log|S|, sign(S) and the decimal digits lost to cancellation.
    """
    lterm = np.asarray(lterm, dtype=float)
    sterm = np.asarray(sterm, dtype=float)
    keep = np.isfinite(lterm) & (sterm != 0)
    if not np.any(keep):
        return -np.inf, 0.0, 0.0
    lterm = lterm[keep]
    sterm = sterm[keep]
    a = lterm.max()
    s = float(np.sum(sterm * np.exp(lterm - a)))
    sgn = float(np.sign(s))
    if s == 0:
        return -np.inf, 0.0, np.inf
    lS = a + np.log(abs(s))
    # max(exp(lterm-a)) is 1, so -log10|s| is the shortfall of the sum against
    # its largest term. Every one of the n terms carries a rounding error of
    # order eps*max_term, so the digits actually lost are that shortfall PLUS
    # log10(n); dropping the count understates the loss and lets a wrong answer
    # past the guard.
    lossDigits = max(0.0, float(np.log10(lterm.size / abs(s))))
    return lS, sgn, lossDigits


def _gdistinct(th: np.ndarray, Nt: float, K: int) -> Tuple[float, float, float]:
    """Eq. (14): single-class constant at pairwise distinct demands th.

    A zero demand contributes nothing, which also realizes the 0/0 = 0
    convention of Eq. (15) when the zero is repeated.
    """
    lin = np.full(K, -np.inf)
    sgv = np.zeros(K)
    for k in range(K):
        if th[k] <= 0:
            continue
        d = th[k] - np.delete(th, k)
        lin[k] = (Nt + K - 1) * np.log(th[k]) - np.sum(np.log(np.abs(d)))
        sgv[k] = np.prod(np.sign(d))
    return _signed_logsumexp(lin, sgv)


def _grepeated(th: np.ndarray, Nt: float, K: int, tol: float) -> Tuple[float, float, float]:
    """Eq. (16): single-class constant at demands th of arbitrary multiplicity.

    Demands within tol of each other, relatively to the largest one, are merged
    into one distinct value carrying their count.
    """
    ths = np.sort(np.asarray(th, dtype=float).ravel())
    scale = ths[-1]
    if scale <= 0:
        scale = 1.0
    gid = np.cumsum(np.concatenate(([1], (np.diff(ths) > tol * scale).astype(int))))
    Kp = int(gid[-1])
    thd = np.zeros(Kp)
    m = np.zeros(Kp, dtype=int)
    for j in range(Kp):
        sel = gid == (j + 1)
        thd[j] = ths[sel].mean()  # the centroid represents a cluster of near-ties
        m[j] = int(np.sum(sel))
    lin = []
    sgv = []
    for j in range(Kp):
        if thd[j] <= 0:
            # the exponent Nt+K-m_j is at least Nt>=1, so a zero cluster
            # contributes nothing
            continue
        louter = (Nt + K - m[j]) * np.log(thd[j])
        souter = (-1.0) ** (m[j] - 1)
        rs = multichoose(Kp, int(m[j]) - 1)  # every K'-vector r>=0 with sum(r)=m_j-1
        for i in range(rs.shape[0]):
            r = rs[i, :]
            lval = louter + _nchoosekln(Nt + r[j], r[j])
            sval = souter * ((-1.0) ** r[j])
            for k in range(Kp):
                if k == j:
                    continue
                lval = lval + _nchoosekln(m[k] + r[k] - 1, r[k])
                if r[k] > 0:
                    if thd[k] <= 0:
                        lval = -np.inf  # theta_k^r_k vanishes, 0^0=1 is the r_k=0 case
                        break
                    lval = lval + r[k] * np.log(thd[k])
                dd = thd[j] - thd[k]
                lval = lval - (m[k] + r[k]) * np.log(abs(dd))
                sval = sval * np.sign(dd) ** (m[k] + r[k])
            lin.append(lval)
            sgv.append(sval)
    return _signed_logsumexp(np.array(lin), np.array(sgv))


def pfqn_explicit(L, N, tol: Optional[float] = None, method: str = 'auto',
                  maxloss: float = np.inf) -> Tuple[float, float, str, float]:
    """Explicit closed-form normalizing constant of a multiclass closed network.

    Evaluates the two explicit expressions of Casale, "Accelerating Performance
    Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS 2017,
    Eqs. (15) and (16). Both instantiate the divided-difference form of
    Corollary 3.2,

        G(N) = sum_{0<=t<=N} (-1)^(|N|-|t|)/(N_1!...N_R!) prod_r C(N_r,t_r) g_t(|N|)

    by substituting a closed form for the single-class constant g_t(|N|) at the
    induced demands theta_k(t) = sum_r t_r L(k,r). Eq. (15) is Gordon's partial
    fraction and needs the induced demands PAIRWISE DISTINCT; Eq. (16) is the
    general partial-fraction expansion over the distinct values and their
    multiplicities, and reduces to Eq. (15) when every multiplicity is one. The
    choice is automatic: Eq. (16) is used as soon as two induced demands are
    closer than tol relative to the largest one at that t.

    SINGLE CLASS. At R=1 the multiclass constant IS the single-class constant at
    demands L, so the outer sum is skipped: g_t(N) = t^N g_1(N) and
    sum_t (-1)^(N-t) t^N/(t!(N-t)!) = S(N,N) = 1. Running the difference anyway
    would add N alternating terms, and their cancellation, to a closed form that
    carries none of them. What is left is O(K^2) work at any population.

    Only single-server load-independent queues are admissible: infinite servers
    need the integral form of Corollary 3.4 and load-dependent rates need the
    load-dependent generalization of the outer sum.

    NUMERICS. Both expressions alternate in sign with terms far larger than the
    result, so they are evaluated as signed log-sum-exps: this removes the
    floating-point RANGE problem but not the cancellation, which is what makes
    multiprecision arithmetic necessary on all but small models.

    Args:
        L: Service demand matrix (KxR) of single-server load-independent queues.
        N: Population vector (1xR).
        tol: Relative tolerance declaring two induced demands redundant
            (default: machine epsilon).
        method: 'auto' (default), 'distinct' to force Eq. (15), 'repeated' to
            force Eq. (16).
        maxloss: Cancellation budget in decimal digits. Finite values turn the
            warnings into a silent REFUSAL (lG=nan) once the budget is exceeded,
            for callers that hold a fallback; default inf keeps the warnings.

    Returns:
        (lG, G, method, lossDigits): the logarithm of the normalizing constant,
        the constant, the expression actually used ('distinct' or 'repeated'),
        and the decimal digits lost to cancellation.
    """
    N = np.atleast_1d(np.asarray(N, dtype=float)).ravel()
    R = N.size
    lossDigits = 0.0
    if tol is None:
        tol = np.finfo(float).eps  # the tolerance is relative to max(theta)
    if method is None or method == '':
        method = 'auto'
    if method not in ('auto', 'distinct', 'repeated'):
        raise ValueError("Unrecognized method, use 'auto', 'distinct' (Eq. 15) or 'repeated' (Eq. 16).")
    if np.sum(N) < 0:
        return -np.inf, 0.0, 'distinct', 0.0
    if np.sum(N) == 0:
        return 0.0, 1.0, 'distinct', 0.0
    L = np.asarray(L, dtype=float)
    if L.ndim == 1:
        L = L.reshape(-1, 1)
    if L.size == 0:
        return -np.inf, 0.0, 'distinct', 0.0
    if L.shape[1] != R:
        raise ValueError('the demand matrix must have one column per class of N.')
    if np.any(L < 0):
        raise ValueError('the demand matrix must be nonnegative.')
    K = L.shape[0]
    Nt = float(np.sum(N))

    # ---- redundancy scan: are the induced demands pairwise distinct at every t? ----
    isRedundant = False
    if R == 1:
        # the induced demands at t are t*L, so both the tie structure and the
        # relative tolerance are those of L itself, at every t at once
        th = np.sort(L[:, 0])
        scale = th[-1]
        isRedundant = bool(scale > 0 and np.any(np.diff(th) <= tol * scale))
    else:
        for t in _population_lattice(N):
            if t.sum() > 0:
                th = np.sort(L @ t)
                scale = th[-1]
                # scale==0 leaves every induced demand at zero, so g_t(|N|)=0 at
                # sum(N)>0 and the term takes no part in the sum
                if scale > 0 and np.any(np.diff(th) <= tol * scale):
                    isRedundant = True
                    break
    if method == 'auto':
        method = 'repeated' if isRedundant else 'distinct'
    elif method == 'distinct' and isRedundant:
        raise ValueError('Eq. (15) requires pairwise distinct induced demands, but two of them '
                         "agree to within tol. Use 'auto' or 'repeated'.")

    if R == 1:
        # ---- single class: the divided difference is the identity, evaluate g ----
        if method == 'distinct':
            lG, sgn, lossDigits = _gdistinct(L[:, 0], Nt, K)
        else:
            lG, sgn, lossDigits = _grepeated(L[:, 0], Nt, K, tol)
    else:
        # ---- outer divided-difference sum over 0 <= t <= N ----
        lterm = []
        sterm = []
        innerLoss = 0.0
        for t in _population_lattice(N):
            if t.sum() <= 0:
                continue
            th = L @ t
            if th.max() <= 0:
                continue
            if method == 'distinct':
                lg, sg, dl = _gdistinct(th, Nt, K)
            else:
                lg, sg, dl = _grepeated(th, Nt, K, tol)
            innerLoss = max(innerLoss, dl)
            if sg != 0:
                lterm.append(lg - np.sum(_factln(t)) - np.sum(_factln(N - t)))
                sterm.append(sg * ((-1.0) ** (Nt - t.sum())))
        lG, sgn, lossDigits = _signed_logsumexp(np.array(lterm), np.array(sterm))
        lossDigits = max(lossDigits, innerLoss)

    # A caller that named a cancellation budget has a fallback and wants a
    # verdict, not a warning: refuse quietly. lossDigits is inf when the sum
    # vanished identically, which is a total loss rather than a legitimate G=0.
    if np.isfinite(maxloss) and (sgn < 0 or lossDigits > maxloss):
        return np.nan, np.nan, method, lossDigits
    if sgn == 0:
        return -np.inf, 0.0, method, lossDigits
    if sgn < 0:
        from ..io.logging import line_warning
        line_warning('pfqn_explicit',
                     'The explicit expression returned a negative value, double precision is '
                     'exhausted by cancellation (%.1f digits lost). Multiprecision arithmetic '
                     'is required.\n' % lossDigits)
        return np.nan, np.nan, method, lossDigits
    G = float(np.exp(lG))
    if lossDigits > 15:
        from ..io.logging import line_warning
        line_warning('pfqn_explicit',
                     'Cancellation has consumed about %.1f decimal digits, more than double '
                     'precision carries. The result is unreliable, multiprecision arithmetic '
                     'is required.\n' % lossDigits)
    return float(lG), G, method, float(lossDigits)


def _population_lattice(N: np.ndarray):
    """Every integer vector 0 <= t <= N, in lexicographic order."""
    N = np.asarray(N, dtype=int)
    t = np.zeros(N.size, dtype=int)
    while True:
        yield t.astype(float)
        r = N.size - 1
        while r >= 0 and t[r] == N[r]:
            t[r] = 0
            r -= 1
        if r < 0:
            return
        t[r] += 1
