"""
Mean of the k-th smallest of n independent EXPONENTIAL branch completion times, i.e. the
instant a k-of-n (quorum) join fires. k = n is the ordinary AND-join, the maximum, and
k = 1 the minimum.

With lambda_i = 1/ri(i) and m = n-k stragglers allowed,

    E[X_(k)] = sum_{j=m+1..n} (-1)^(j-m-1) C(j-1,m) e_j,
    e_j      = sum_{|S|=j} 1 / sum_{i in S} lambda_i

the inclusion-exclusion identity for the order statistics of independent exponentials. At
m = 0 it collapses to sum_j (-1)^(j-1) e_j, the classical expression for the maximum, TERM
BY TERM: a full join therefore evaluates exactly as it did before this module existed.

The sum has 2^n terms and its signs alternate, so it is evaluated exactly only while the
branch count is small. Beyond MAXEXACT branches a genuine quorum (k < n) is evaluated by
quorum_moments instead, whose Poisson-binomial recurrence adds no cancellation; a full
join keeps the exact path at every n so that no existing result moves.

Port of matlab/src/api/fj/fj_ordstat_exp.m, mirrored by jline.api.fj.FJ_ordstat_exp.

Reference: A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems", ACM
Computing Surveys 47(2), Article 17, 2014, Sec. 3 (Eq. 18-19).
"""

from itertools import combinations
from math import comb, isinf, isnan
from typing import Sequence

# Branch count above which a genuine quorum leaves the exact alternating sum.
MAXEXACT = 15


def fj_ordstat_exp(ri: Sequence[float], k: int) -> float:
    """
    Mean instant a k-of-n join over exponential branches fires.

    Args:
        ri: per-branch mean completion times
        k: quorum, 1 <= k <= len(ri)

    Returns:
        The mean of the k-th order statistic.
    """
    r = [float(x) for x in ri if not (isnan(float(x)) or isinf(float(x)))]
    n = len(r)
    if n == 0:
        return 0.0
    if k < 1 or k > n:
        raise ValueError('fj_ordstat_exp: k must satisfy 1 <= k <= n. Got k=%d, n=%d.' % (k, n))
    # A branch of zero mean completes instantly: it never delays the join and it counts
    # toward the quorum at once. Removing it here keeps the reciprocal below finite, which
    # an exact arithmetic requires and IEEE only tolerates.
    nzero = sum(1 for x in r if x <= 0.0)
    if nzero > 0:
        k -= nzero
        if k <= 0:
            return 0.0
        r = [x for x in r if x > 0.0]
        n = len(r)
    if n == 1:
        return r[0]

    if k < n and n > MAXEXACT:
        from ...lib.thirdparty.fj import quorum_moments
        # Branch times are taken as exponential, so the variance is the square of the mean.
        return quorum_moments(r, [x * x for x in r], k)[0]

    lambdai = [1.0 / x for x in r]
    nstrag = n - k
    total = 0.0
    for j in range(nstrag + 1, n + 1):
        ej = 0.0
        for subset in combinations(lambdai, j):
            ej += 1.0 / sum(subset)
        total += ((-1.0) ** (j - nstrag - 1)) * comb(j - 1, nstrag) * ej
    return total
