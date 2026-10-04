"""
Completion time of a k-of-n (quorum) join over independent, not necessarily identically
distributed branches.

Each branch is summarised by its first two moments and expanded into a discrete step CDF
by a two-moment fit. The CDF of the k-th order statistic is then assembled by the
inclusion-exclusion identity

    F_(k)(t) = sum_{i=k..n} (-1)^(i-k) * C(i-1, k-1) * e_i(F_1(t), ..., F_n(t))

where e_i is the i-th elementary symmetric polynomial of the branch CDFs. For k = n this
collapses to the product of the branch CDFs, i.e. the ordinary AND-join, and for k = 1 to
1 - prod(1 - F_i), i.e. the minimum.

Follows the formulation of Omari, Franks, Woodside and Pan, as implemented in LQNS 6.x
(randomvar.cc). Mirrors jline.api.fj.FJ_quorum in the Java runtime and
fj_quorum_moments.m in the MATLAB runtime.
"""

from math import sqrt
from typing import List, Sequence, Tuple


class StepCDF:
    """
    A discrete step function on a finite, increasing time grid.

    Holds parallel lists t[0..m-1] of abscissae and a[0..m-1] of accumulated values, with
    the value implicitly 0 below t[0] and constant at a[m-1] above t[m-1]. When the values
    are a probability distribution, a is its CDF; intermediate results of the
    inclusion-exclusion sum are not distributions and may leave [0, 1].
    """

    __slots__ = ('t', 'a')

    def __init__(self, t: Sequence[float], a: Sequence[float]):
        if len(t) != len(a):
            raise ValueError("t and a must have equal length, got %d and %d." % (len(t), len(a)))
        self.t = list(t)
        self.a = list(a)

    def is_empty(self) -> bool:
        return len(self.t) == 0

    def at(self, x: float) -> float:
        """Value of the step function at time x."""
        v = 0.0
        for ti, ai in zip(self.t, self.a):
            if ti > x:
                break
            v = ai
        return v

    def times(self, c: float) -> 'StepCDF':
        """Scale the values, leaving the time grid untouched."""
        return StepCDF(list(self.t), [ai * c for ai in self.a])

    def mean(self) -> float:
        m = 0.0
        prev = 0.0
        for ti, ai in zip(self.t, self.a):
            m += (ai - prev) * ti
            prev = ai
        return m

    def variance(self) -> float:
        m = self.mean()
        v = 0.0
        prev = 0.0
        for ti, ai in zip(self.t, self.a):
            v += (ai - prev) * (ti - m) ** 2
            prev = ai
        return max(v, 0.0)


def three_point_fit(mean: float, variance: float) -> StepCDF:
    """
    Two-moment fit of a branch completion time to a three-point discrete distribution.

    Matches the mean and variance exactly. A zero standard deviation degenerates to a
    single deterministic point, and a zero mean yields the empty step function.
    """
    if mean < 0.0 or variance < 0.0:
        raise ValueError("mean and variance must be non-negative, got mean=%r, variance=%r."
                         % (mean, variance))
    if mean == 0.0:
        return StepCDF([], [])
    sd = sqrt(variance)
    if sd == 0.0:
        return StepCDF([mean], [1.0])

    t1 = mean - sd if mean > sd else 0.0
    t2 = mean
    t3 = mean + 2.0 * variance / mean if sd >= mean else mean + 2.0 * sd

    delta = t1 * t1 * (t3 - t2) + t2 * t2 * (t1 - t3) + t3 * t3 * (t2 - t1)
    if delta == 0.0:
        # The three abscissae are not distinct, so the fit is not determined.
        return StepCDF([mean], [1.0])
    temp = variance + mean * mean

    a1 = (temp * (t3 - t2) + t2 * t2 * (mean - t3) + t3 * t3 * (t2 - mean)) / delta
    a3 = (t1 * t1 * (mean - t2) + t2 * t2 * (t1 - mean) + temp * (t2 - t1)) / delta

    return StepCDF([t1, t2, t3], [a1, 1.0 - a3, 1.0])


def _combine(x: StepCDF, y: StepCDF, is_add: bool) -> StepCDF:
    """Pointwise sum or product of two step functions on the union of their grids."""
    nx, ny = len(x.t), len(y.t)
    ts: List[float] = []
    vs: List[float] = []
    i = j = 0
    x_prev = y_prev = 0.0
    while i < nx or j < ny:
        # Advance whichever grid holds the next time, stepping both when they coincide so
        # a shared abscissa is emitted once.
        if j >= ny or (i < nx and x.t[i] < y.t[j]):
            time = x.t[i]
            x_prev = x.a[i]
            i += 1
        elif i >= nx or y.t[j] < x.t[i]:
            time = y.t[j]
            y_prev = y.a[j]
            j += 1
        else:
            time = x.t[i]
            x_prev = x.a[i]
            y_prev = y.a[j]
            i += 1
            j += 1
        ts.append(time)
        vs.append(x_prev + y_prev if is_add else x_prev * y_prev)
    return StepCDF(ts, vs)


def add(x: StepCDF, y: StepCDF) -> StepCDF:
    """Pointwise sum of two step functions."""
    return _combine(x, y, True)


def multiply(x: StepCDF, y: StepCDF) -> StepCDF:
    """
    Pointwise product of two step functions.

    An empty operand annihilates the product, matching the convention that an absent
    branch distribution contributes no mass.
    """
    if x.is_empty() or y.is_empty():
        return StepCDF([], [])
    return _combine(x, y, False)


#: Largest branch count accepted by :func:`quorum_kofn`. Evaluation is cubic in the branch
#: count, so this bounds the work at roughly 1e8 elementary operations. An LQN AND-join with
#: more branches is degenerate; the limit turns a silent hang into an immediate failure.
MAX_BRANCHES = 512


def quorum_kofn(branches: Sequence[StepCDF], k: int) -> StepCDF:
    """
    CDF of the k-th smallest of n independent branch completion times.

    Evaluated pointwise on the union of the branch grids. At each time the number of
    completed branches is Poisson-binomial, so its distribution is built by the recurrence
    q_j <- q_{j-1} * F_i + q_j * (1 - F_i) over branches i, and the k-th order statistic is
    the upper tail sum_{j>=k} q_j.

    This is the same quantity as the inclusion-exclusion identity
    sum_{i=k..n} (-1)^(i-k) C(i-1,k-1) e_i(F_1..F_n) used by LQNS, but every term here is a
    probability in [0, 1] and none is subtracted, so it does not suffer the catastrophic
    cancellation the alternating binomial sum incurs as n grows.

    Args:
        branches: the branch CDFs, one per branch
        k: the quorum count, in [1, len(branches)]

    Raises:
        ValueError: if k is out of range or there are more than MAX_BRANCHES branches
    """
    n = len(branches)
    if n == 0:
        return StepCDF([], [])
    if k < 1 or k > n:
        raise ValueError("k must satisfy 1 <= k <= n. Got k=%d, n=%d." % (k, n))
    if n > MAX_BRANCHES:
        raise ValueError("quorum join has %d branches, above the supported maximum of %d; "
                         "evaluation is cubic in the branch count." % (n, MAX_BRANCHES))

    # A branch with no mass completes instantly, so that it neither delays the join nor
    # suppresses the completion counts below.
    b = [StepCDF([0.0], [1.0]) if br.is_empty() else br for br in branches]

    grid = sorted(set(ti for br in b for ti in br.t))
    m = len(grid)

    # Tabulate each branch along the merged grid by a single monotone walk, so the
    # evaluation below is linear rather than quadratic in the grid size.
    fv = []
    for br in b:
        p = 0
        cur = 0.0
        row = [0.0] * m
        for g in range(m):
            while p < len(br.t) and br.t[p] <= grid[g]:
                cur = br.a[p]
                p += 1
            row[g] = cur
        fv.append(row)

    out: List[float] = []
    for g in range(m):
        q = [0.0] * (n + 1)
        q[0] = 1.0
        for i in range(n):
            f = fv[i][g]
            for j in range(min(i + 1, n), 0, -1):
                q[j] = q[j - 1] * f + q[j] * (1.0 - f)
            q[0] = q[0] * (1.0 - f)
        out.append(sum(q[k:]))
    return StepCDF(grid, out)


def quorum_moments(branch_means: Sequence[float], branch_variances: Sequence[float],
                   k: int) -> Tuple[float, float]:
    """
    Mean and variance of the completion time of a k-of-n join whose branches are given by
    their first two moments.

    This is the entry point used by the layered solver: each branch mean and variance
    comes from the residence time and variance accumulated along that branch of the
    activity graph.

    Returns:
        a (mean, variance) pair for the join completion time
    """
    if len(branch_means) != len(branch_variances):
        raise ValueError("branch_means and branch_variances must have equal length, got %d and %d."
                         % (len(branch_means), len(branch_variances)))
    n = len(branch_means)
    if n == 0:
        return 0.0, 0.0
    branches = [three_point_fit(m, v) for m, v in zip(branch_means, branch_variances)]
    join = quorum_kofn(branches, k)
    return join.mean(), join.variance()
