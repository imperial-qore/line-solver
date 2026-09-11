"""Explicit closed-form normalizing constant of a multiclass limited load-dependent network.

Casale, Harrison and Ong, "Facilitating Load-Dependent Queueing Analysis
Through Factorization", Perform. Eval. 2021, Theorem 1, carried over the
divided-difference form of Casale, "Accelerating Performance Inference over
Closed Systems by Asymptotic Methods", ACM SIGMETRICS 2017, Corollary 3.2.
Twin of MATLAB pfqn_explicit_ld.m.
"""

from typing import Optional, Tuple

import numpy as np

from .explicit import (_factln, _gdistinct, _grepeated, _population_lattice,
                       _signed_logsumexp)

__all__ = ['pfqn_explicit_ld']


def _hlld(th, M, Nt, alphaS, vcap, lcum, lbr, sbr, method, tol):
    """Theorem 1, Eq. (8): the single-class LLD constant at induced demands th.

    The finite sum over 0 <= v < s of the fixed-rate constant at the scaled
    demands th/alpha(s), one population level lower for every job held back by
    v. Returns log|h|, sign(h) and the decimal digits lost to cancellation.
    """
    sigma = th / alphaS
    with np.errstate(divide='ignore'):
        lth = np.where(th > 0, np.log(np.abs(th)), -np.inf)
    lterm = []
    sterm = []
    lossDigits = 0.0
    for v in _population_lattice(vcap):
        v = v.astype(int)
        nv = int(v.sum())
        if nv > Nt:
            continue
        lval = 0.0
        sval = 1.0
        for i in range(M):
            vi = v[i]
            if vi > 0:
                if not np.isfinite(lth[i]):
                    # theta_k = 0 kills every v_k>0, and 0^0=1 keeps v_k=0
                    sval = 0.0
                    break
                # kept inside the guard because 0*(-inf) is nan, not 0
                lval += vi * lth[i]
            lval += -lcum[i][vi] + lbr[i][vi]
            sval *= sbr[i][vi]
        if sval == 0 or not np.isfinite(lval):
            continue
        if Nt - nv == 0:
            # g_sigma(0) = 1 by definition. Reading it off the partial fraction
            # instead would spend digits on an alternating sum whose value is
            # known exactly.
            lg, sg, dl = 0.0, 1.0, 0.0
        elif method == 'distinct':
            lg, sg, dl = _gdistinct(sigma, float(Nt - nv), M)
        else:
            lg, sg, dl = _grepeated(sigma, float(Nt - nv), M, tol)
        lossDigits = max(lossDigits, dl)
        if sg == 0:
            continue
        lterm.append(lval + lg)
        sterm.append(sval * sg)
    lh, sh, dl = _signed_logsumexp(np.array(lterm), np.array(sterm))
    return lh, sh, max(lossDigits, dl)


def pfqn_explicit_ld(L, N, mu=None, tol: Optional[float] = None, method: str = 'auto',
                     maxloss: float = np.inf) -> Tuple[float, float, str, float]:
    r"""Explicit closed-form normalizing constant of a multiclass LLD network.

    Load-dependent counterpart of pfqn_explicit. It evaluates the same
    divided-difference form of Casale (SIGMETRICS 2017), Corollary 3.2,

        G(N) = sum_{0<=t<=N} (-1)^(\|N\|-\|t\|)/(N_1!...N_R!) prod_r C(N_r,t_r) h_t(\|N\|)

    but substitutes for the single-class constant h_t(\|N\|) the LIMITED
    LOAD-DEPENDENT closed form of Casale, Harrison and Ong (Perform. Eval.
    2021), Theorem 1, Eq. (8),

        h_theta(N) = sum_{0<=v<s} g_sigma(N-\|v\|) prod_k phi_k(v_k)
        phi_k(v_k) = theta_k^v_k / prod_{t=1..v_k} alpha_k(t) * (1 - alpha_k(v_k)/alpha_k(s_k))

    at the induced demands theta_k(t) = sum_r t_r L(k,r). Here alpha_k(.) =
    mu(k,.) is the load-dependent scaling of station k, s_k the population past
    which it stays constant, sigma_k = theta_k/alpha_k(s_k) the SCALED demands,
    and g_sigma the FIXED-RATE single-class constant at those scaled demands,
    which is exactly what pfqn_explicit evaluates in closed form (Eqs. 15 and
    16). The result is therefore explicit throughout, with no recursion over
    population.

    Two conventions of Theorem 1 are not those of the equilibrium distribution
    and are easy to get wrong. alpha_k(0) is taken as ZERO inside the bracket of
    phi_k, so that phi_k(0) = 1, even though the state probabilities use
    alpha_k(0) = 1; and g_sigma(n) = 0 for n < 0, which caps the outer sum at
    \|v\| <= \|N\|. With alpha_k(n) = min(n,s_k) the expression collapses to
    Gordon's multi-server formula, Oper. Res. 38(5), 1990, Eq. (29), but unlike
    that one it needs neither a multi-server shape nor distinct scaled demands.

    LIMITED LOAD DEPENDENCE. Theorem 1 holds for any s_k with
    alpha_k(n) = alpha_k(s_k) for all n >= s_k, and a LARGER s_k is always
    admissible, so s_k is detected here as the smallest index whose value the
    tail of mu(k,:) repeats to within tol. A station whose rates never settle
    (an infinite server, mu(k,n) = n) gets s_k = \|N\|, which is still exact:
    populations above \|N\| do not occur, so redefining alpha_k there changes
    nothing. It is merely expensive, since the inner sum costs prod_k s_k terms,
    capped by \|v\| <= \|N\|. Think time is not admissible: a delay would have to
    enter g_sigma, whose closed form covers queues only.

    NUMERICS. Both sums alternate in sign with terms far larger than the result,
    so they are evaluated as signed log-sum-exps. phi_k is sign-definite when
    alpha_k increases, as a multi-server station does, and changes sign where
    alpha_k decreases, so a decreasing rate function costs digits in the inner
    sum too.

    SINGLE CLASS. At R=1 the divided difference is the identity, since
    h_theta(N) is homogeneous of degree N in theta exactly as in the fixed-rate
    case, so the outer sum is skipped and Theorem 1 is evaluated once at
    theta = L.

    Args:
        L: Service demand matrix (MxR).
        N: Population vector (1xR).
        mu: Load-dependent rate matrix (Mx sum(N)), alpha_i(j) = mu[i,j-1];
            default all ones.
        tol: Relative tolerance declaring two scaled demands redundant, and the
            rate tail constant (default: machine epsilon).
        method: 'auto' (default), 'distinct' to force Eq. (15), 'repeated' to
            force Eq. (16).
        maxloss: Cancellation budget in decimal digits. Finite values turn the
            warnings into a silent REFUSAL (lG=nan) once the budget is exceeded,
            for callers that hold a fallback; default inf keeps the warnings.

    Returns:
        (lG, G, method, lossDigits): the logarithm of the normalizing constant,
        the constant, the expression used for g_sigma ('distinct' or
        'repeated'), and the decimal digits lost to cancellation.
    """
    N = np.atleast_1d(np.asarray(N, dtype=float)).ravel()
    R = N.size
    lossDigits = 0.0
    if tol is None:
        tol = np.finfo(float).eps  # the tolerance is relative to max(sigma)
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
    M = L.shape[0]
    Nt = int(np.sum(N))
    if mu is None:
        mu = np.ones((M, Nt))
    mu = np.asarray(mu, dtype=float)
    if mu.ndim == 1:
        mu = mu.reshape(1, -1)
    if mu.shape[0] != M:
        raise ValueError('the load-dependent rate matrix must have one row per station of L.')
    if mu.shape[1] < Nt:
        raise ValueError('the load-dependent rate matrix must have at least sum(N) columns.')
    mu = mu[:, :Nt]
    if np.any(mu <= 0):
        raise ValueError('the load-dependent rates must be strictly positive.')

    # ---- s_k: the smallest index whose value the tail of the rate row repeats ----
    # Any larger s_k also satisfies alpha_k(n)=alpha_k(s_k) for n>=s_k, so a
    # missed tie only adds terms; a false tie would be a wrong answer, hence the
    # strict tol.
    s = np.ones(M, dtype=int)
    for i in range(M):
        s[i] = Nt
        tail = mu[i, Nt - 1]
        n = Nt
        while n > 1 and abs(mu[i, n - 2] - tail) <= tol * max(abs(tail), 1.0):
            n -= 1
            s[i] = n
    alphaS = np.array([mu[i, s[i] - 1] for i in range(M)])

    # ---- per-station phi tables, in the log domain, indexed by v_k = 0..s_k-1 ----
    lcum = []  # sum_{t=1..v} log alpha_k(t)
    lbr = []   # log|1 - alpha_k(v)/alpha_k(s_k)|, with alpha_k(0) := 0
    sbr = []
    for i in range(M):
        lcum.append(np.concatenate(([0.0], np.cumsum(np.log(mu[i, :s[i] - 1])))))
        br = 1.0 - np.concatenate(([0.0], mu[i, :s[i] - 1])) / alphaS[i]
        lb = np.full(s[i], -np.inf)
        lb[br != 0] = np.log(np.abs(br[br != 0]))
        lbr.append(lb)
        sbr.append(np.sign(br))
    # g_sigma vanishes below zero population, Eq. (8) caps |v| <= |N|
    vcap = np.minimum(s - 1, Nt).astype(float)

    # ---- redundancy scan: are the SCALED induced demands pairwise distinct? ----
    # The scan MUST form sigma exactly as _hlld does, (L @ t)/alphaS and not
    # (L/alphaS) @ t: the two orderings differ in the last ulp, so an exact tie
    # can clear an eps-relative gap under one and not the other, and Eq. (15)
    # would then divide by that ulp. Measured on L=[[0,1.3],[0.9,0.7]],
    # mu=min(n,2), N=[2,3]: at t=[2,3] both scaled demands are 1.95, the
    # scaled-first ordering reports a 4.4e-16 gap and misses the tie, the
    # demand-first ordering reports 2.2e-16 and catches it.
    isRedundant = False
    if R == 1:
        # the scaled demands at t are t*sigma, so both the tie structure and the
        # relative tolerance are those of sigma itself, at every t at once
        th = np.sort(L[:, 0] / alphaS)
        scale = th[-1]
        isRedundant = bool(scale > 0 and np.any(np.diff(th) <= tol * scale))
    else:
        for t in _population_lattice(N):
            if t.sum() > 0:
                th = np.sort((L @ t) / alphaS)
                scale = th[-1]
                # scale==0 leaves every scaled demand at zero, so the term takes
                # no part in the sum
                if scale > 0 and np.any(np.diff(th) <= tol * scale):
                    isRedundant = True
                    break
    if method == 'auto':
        method = 'repeated' if isRedundant else 'distinct'
    elif method == 'distinct' and isRedundant:
        raise ValueError('Eq. (15) requires pairwise distinct scaled demands, but two of them '
                         "agree to within tol. Use 'auto' or 'repeated'.")

    if R == 1:
        # ---- single class: the divided difference is the identity ----
        lG, sgn, lossDigits = _hlld(L[:, 0], M, Nt, alphaS, vcap, lcum, lbr, sbr, method, tol)
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
            lh, sh, dl = _hlld(th, M, Nt, alphaS, vcap, lcum, lbr, sbr, method, tol)
            innerLoss = max(innerLoss, dl)
            if sh == 0:
                continue
            lterm.append(lh - np.sum(_factln(t)) - np.sum(_factln(N - t)))
            sterm.append(sh * ((-1.0) ** (Nt - t.sum())))
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
        line_warning('pfqn_explicit_ld',
                     'The explicit expression returned a negative value, double precision is '
                     'exhausted by cancellation (%.1f digits lost). Multiprecision arithmetic '
                     'is required.\n' % lossDigits)
        return np.nan, np.nan, method, lossDigits
    G = float(np.exp(lG))
    if lossDigits > 15:
        from ..io.logging import line_warning
        line_warning('pfqn_explicit_ld',
                     'Cancellation has consumed about %.1f decimal digits, more than double '
                     'precision carries. The result is unreliable, multiprecision arithmetic '
                     'is required.\n' % lossDigits)
    return float(lG), G, method, float(lossDigits)
