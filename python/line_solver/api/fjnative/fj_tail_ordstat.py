"""
Tail latency of a k-of-n (quorum) fork-join request.

Port of matlab/src/api/fj/fj_tail_ordstat.m, mirrored by
jline.api.fj.FJ_tail_ordstat.
"""

import numpy as np
from scipy.optimize import brentq
from scipy.special import betainc

from .fj_tail_forktail import fj_tail_forktail, ge_fit


def _ge_cdf(x, alpha, beta):
    """The generalized-exponential CDF (1-exp(-x/beta))**alpha, per branch.

    Clamped into [0,1] so that a rounding excursion cannot leave the beta or
    binomial routines out of domain.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        u = np.exp(alpha * np.log1p(-np.exp(-x / beta)))
    u = np.where(np.isfinite(u), u, 0.0)
    return np.clip(u, 0.0, 1.0)


def _poissbin_upper(u, kreq):
    """P(at least kreq successes) for independent Bernoulli trials of success
    probabilities u, by the convolution recurrence over the trials.

    Every term is non-negative, so no cancellation is introduced.
    """
    u = np.atleast_1d(np.asarray(u, dtype=float)).ravel()
    n = u.size
    pmf = np.zeros(n + 1)
    pmf[0] = 1.0
    for i in range(n):
        pmf[1:i + 2] = pmf[1:i + 2] * (1.0 - u[i]) + pmf[0:i + 1] * u[i]
        pmf[0] *= (1.0 - u[i])
    return float(pmf[kreq:].sum())


def fj_tail_ordstat(ET, VT, K=None, p=99, kreq=None):
    """
    Predict the p-th percentile of a k-of-n fork-join request response time.

    The request forks into N parallel tasks and joins on the KREQ-th of them.
    KREQ = N is the ordinary AND-join and reproduces fj_tail_forktail exactly;
    KREQ = 1 is the first completion.

    Each branch is the same black box ForkTail uses: its task response time is
    fitted by a generalized exponential law F_i(x) = (1-exp(-x/beta_i))**alpha_i
    matched on the branch mean and variance. The request completes once KREQ of
    the N branches have, so its law is the KREQ-th ORDER STATISTIC of
    independent, not identically distributed branch times,

        F_X(x) = P(at least KREQ of the N branches are done by x),

    the upper tail of a Poisson-binomial with success probabilities F_i(x). It
    is evaluated by the convolution recurrence, which adds no cancellation, and
    inverted by bisection. With homogeneous branches the recurrence collapses to
    the regularized incomplete beta function I_{F(x)}(KREQ, N-KREQ+1), used
    directly. At KREQ = N both routes reduce term by term to prod_i F_i(x), the
    product of the branch CDFs that fj_tail_forktail inverts.

    BRANCH INDEPENDENCE is assumed, as in ForkTail: the branches of one request
    are positively correlated through their shared arrival instant, so the true
    quorum percentile is somewhat larger than this one. The same heavy-traffic
    caveat applies.

    Args:
        ET: mean task response time, a scalar (homogeneous branches) or one
            entry per branch
        VT: variance of the task response time, same shape as ET
        K: number of branches; used only when ET is a scalar
        p: percentile, a fraction in (0,1) or a percentage in (0,100)
        kreq: the join fires on the kreq-th branch (default: every branch)

    Returns:
        (xp, alpha, beta): the predicted percentile and the fitted parameters
    """
    p = float(p)
    if p > 1.0:
        p = p / 100.0
    if p <= 0.0 or p >= 1.0:
        raise ValueError("The percentile must lie strictly between 0 and 1 (or 0 and 100).")

    ET = np.atleast_1d(np.asarray(ET, dtype=float)).ravel()
    VT = np.atleast_1d(np.asarray(VT, dtype=float)).ravel()
    if ET.size != VT.size:
        raise ValueError("ET and VT must have the same number of entries.")
    if np.any(ET <= 0) or np.any(VT <= 0):
        raise ValueError("The task response time mean and variance must be positive.")

    nbranch = int(ET.size)
    if nbranch == 1:
        nsib = max(1, int(round(float(np.atleast_1d(np.asarray(K if K is not None else 1)).ravel()[0]))))
    else:
        nsib = nbranch
    if kreq is None:
        kreq = nsib
    kreq = int(round(kreq))
    if kreq < 1 or kreq > nsib:
        raise ValueError("The quorum must satisfy 1 <= kreq <= %d. Got %d." % (nsib, kreq))

    fits = [ge_fit(ET[i], VT[i]) for i in range(nbranch)]
    alpha = np.array([f[0] for f in fits])
    beta = np.array([f[1] for f in fits])

    # The AND-join percentile is an upper bound for every quorum, and the single
    # branch percentile a lower one, so the two bracket the root without a search.
    xhi = fj_tail_forktail(ET if nbranch > 1 else ET[0],
                           VT if nbranch > 1 else VT[0], nsib, p)[0]
    if kreq == nsib:
        return float(xhi), alpha, beta
    xlo = float(np.min(-beta * np.log(1.0 - p ** (1.0 / alpha))))

    if nbranch == 1:
        # homogeneous: the count done by x is Binomial(nsib, F(x)), so the
        # quorum CDF is the regularized incomplete beta of its upper tail
        def cdfk(x):
            return float(betainc(kreq, nsib - kreq + 1, _ge_cdf(x, alpha, beta)[0]))
    else:
        def cdfk(x):
            return _poissbin_upper(_ge_cdf(x, alpha, beta), kreq)

    def residual(x):
        return cdfk(x) - p

    # A single branch may already meet the percentile; the quorum is met earlier.
    while residual(xlo) > 0 and xlo > np.finfo(float).tiny:
        xlo /= 2.0
    xp = brentq(residual, xlo, xhi, xtol=1e-14, rtol=1e-14)
    return float(xp), alpha, beta
