"""
Exact analysis of loss networks by the Manjunath-Sikdar transform.

Native Python port of matlab/src/api/lossn/lossn_manjunath.m. Calls on route r arrive
Poisson at rate nu_r with unit mean holding time, so nu_r is the offered load,
and a call is admitted only while every constraint holds,

    sum_r A[j,r] n_r <= C[j],   j = 1 ... J.

The admissible set is coordinate convex, so by Kelly's truncation theorem the
stationary distribution is the truncated product form
p(n) = nu^n / n! / g(C) and every metric is a ratio of normalization constants

    g(C)   = sum_{A n <= C} prod_r nu_r**n_r / n_r!
    E[n_r] = nu_r g(C - A e_r) / g(C)
    Loss_r = 1 - g(C - A e_r) / g(C)

because a class r call is blocked exactly when the state cannot absorb one more
unit of its own requirement vector.

WHY IT IS A COEFFICIENT COMPUTATION AND NOT A QUADRATURE. Writing each indicator
as a contour integral turns g(C) into a J-fold integral over the unit circle
whose integrand factorizes into the per-route z-transforms. Inside the circle the
only pole in z_j sits at the origin with order C_j+1, so each integration is a
residue, i.e. a Taylor coefficient. The routine therefore never evaluates an
integral: it builds the generating function as a multivariate power series
truncated at degree C_j in z_j, one shift-and-accumulate convolution per route,
and discharges each '<=' constraint by summing the coefficients of degrees
0 ... C_j along that dimension. Truncation is exact because A is nonnegative, so
a monomial above degree C_j can never contribute to an extracted coefficient.

THE ELIMINATION ORDER IS THE MEMORY BOUND. Contour integrations are interleaved
with the product rather than deferred: variable z_j is created when the first
route with A[j,r] != 0 is multiplied in and integrated out immediately after the
last one. Peak memory is therefore the product of (C_j+1) over the
SIMULTANEOUSLY LIVE links, an induced width of the route-link incidence, not over
all J links. That product is bounded by `max_live_states` and a region above it is
refused by name rather than allowed to exhaust the machine: the algorithm is
exact but not unconditionally cheap, and `lossn_mci` answers the same question at
any size.

Unlike `lossn_erlangfp` this is exact rather than a reduced-load approximation,
and unlike `lossn_mci` it carries no sampling error, which is what matters for
rare blocking: a loss probability of 1e-4 recovered from a sampled throughput is
dominated by the estimator variance.

Reference: D. Manjunath and B. Sikdar, Integral Expressions for the Numerical
Evaluation of Product Form Expressions Over Irregular Multidimensional Integer
Spaces.
"""

import math
from typing import Optional, Tuple

import numpy as np
from scipy import special

# Cap on the product of (C_j+1) over the simultaneously live links, i.e. on the
# number of series coefficients held at once. 2**26 coefficients is half a
# gigabyte at float64, which no region a FiniteCapacityRegion can express reaches
# by accident. Raise it deliberately.
DEFAULT_MAX_LIVE_STATES = 1 << 26


def _integralize(A: np.ndarray, C: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integer view of (A, C): rows that constrain nothing dropped, every remaining
    row divided by the gcd of its entries together with its right-hand side.

    Both steps are exact. The second is what makes the cost tractable when the
    class sizes share a factor: a row (4, 8) <= 20 becomes (1, 2) <= 5 and its
    dimension shrinks from 21 coefficients to 6. A fractional or negative entry
    is refused rather than rounded, since the residue argument counts whole units
    of capacity.
    """
    if np.any(A < 0) or np.any(np.abs(A - np.round(A)) > 1e-9):
        raise ValueError(
            "lossn_manjunath: A must contain nonnegative integers -- the residue argument "
            "counts whole units of capacity. Use lossn_mci, which compares in real "
            "arithmetic, or lossn_erlangfp")
    if np.any(C < 0) or np.any(np.abs(C - np.round(C)) > 1e-9):
        raise ValueError(
            "lossn_manjunath: C must contain nonnegative integers -- the residue argument "
            "counts whole units of capacity. Use lossn_mci, which compares in real "
            "arithmetic, or lossn_erlangfp")

    Ai = np.round(A).astype(np.int64)
    Ci = np.round(C).astype(np.int64)

    # A row of zeros bounds nothing and is dropped, not carried as a
    # one-coefficient dimension: keeping it would leave `first`/`last` undefined
    # for that link.
    keep = np.any(Ai != 0, axis=1)
    Ai = Ai[keep, :]
    Ci = Ci[keep]

    for j in range(len(Ci)):
        g = int(Ci[j])
        for r in np.nonzero(Ai[j, :])[0]:
            g = math.gcd(g, int(Ai[j, r]))
        if g > 1:
            Ai[j, :] //= g
            Ci[j] //= g
    return Ai, Ci


def _series(f, A: np.ndarray, C: np.ndarray, nmax_full: np.ndarray,
            max_live_states: int, peak: list) -> float:
    """
    Coefficient-domain evaluation of the J-fold contour integral at the
    right-hand side `C`, which is the full rule for g(C) and the rule shifted by
    one class requirement for g(C - A e_r).

    The series is an ndarray whose shape is the live grid: axis j has extent
    C[j]+1 while link j is live and 1 while it is not, so a shift along a dead
    axis is always zero and the translated add below is a plain slice.
    """
    R = len(f)
    J = len(C)

    # The elimination order: link j is created at its first route and summed out
    # after its last, so only an induced width of links is ever live.
    first = np.zeros(J, dtype=np.int64)
    last = np.zeros(J, dtype=np.int64)
    for j in range(J):
        idx = np.nonzero(A[j, :])[0]
        if idx.size == 0:
            raise RuntimeError(
                "lossn_manjunath: a constraint row with no nonzero entry reached the series; "
                "the rule was not reduced")
        first[j] = idx[0]
        last[j] = idx[-1]

    dims = [1] * J
    ser = np.ones([1] * J, dtype=np.float64)

    for r in range(R):
        # 1. Create the links whose first route is this one, keeping the existing
        # content at degree zero in the new variable.
        for j in np.nonzero(first == r)[0]:
            newdim = int(C[j]) + 1
            live = int(np.prod(dims)) // dims[j] * newdim
            if live > max_live_states:
                raise RuntimeError(
                    "lossn_manjunath: the exact transform would hold more than %d series "
                    "coefficients at once. Peak memory is the product of (C_j+1) over "
                    "the links live at the same time, so a wide constraint row with a "
                    "large capacity is what costs; raise max_live_states deliberately, "
                    "or use lossn_mci, which is unbiased at any size" % max_live_states)
            newshape = list(dims)
            newshape[j] = newdim
            grown = np.zeros(newshape, dtype=np.float64)
            sl = [slice(None)] * J
            sl[j] = slice(0, 1)
            grown[tuple(sl)] = ser
            ser = grown
            dims = newshape
            peak[0] = max(peak[0], ser.size)

        # 2. Multiply in route r.
        s = A[:, r]
        if not np.any(s != 0):
            # Bounded by no link, so its z-transform is a constant: the whole
            # truncated sequence sums into the series. For a route absent from
            # every row this is the factor exp(nu_r) the caller folded into
            # lGfree, hence a multiply by 1.
            ser = ser * float(np.sum(f[r][:nmax_full[r] + 1]))
        else:
            # The degree of route r is capped by every row it appears in,
            # evaluated at THIS right-hand side: the shifted series for
            # g(C - A e_r) admits strictly fewer calls than g(C).
            nmax = min(int(nmax_full[r]), len(f[r]) - 1)
            for j in np.nonzero(s > 0)[0]:
                nmax = min(nmax, int(C[j]) // int(s[j]))
            nxt = np.zeros_like(ser)
            for n in range(nmax + 1):
                c = float(f[r][n])
                if c == 0.0:
                    continue
                if n == 0:
                    nxt += c * ser
                    continue
                # Shift by n requirement vectors, dropping the coefficients the
                # shift would push past the capacity. Those monomials can never
                # contribute to an extracted coefficient, which is exactly why
                # the truncation is exact and not an approximation.
                shift = s * n
                if np.any(shift >= np.asarray(dims)):
                    break            # the shift only grows with n
                src = tuple(slice(0, dims[j] - int(shift[j])) for j in range(J))
                dst = tuple(slice(int(shift[j]), dims[j]) for j in range(J))
                nxt[dst] += c * ser[src]
            ser = nxt

        # 3. Integrate out the links whose last route was this one. The
        # multiplier (z^{C+1}-1)/(z-1) of a '<=' constraint turns the residue
        # into the partial sum of the coefficients of degrees 0 ... C_j, which is
        # the sum along that dimension.
        for j in np.nonzero(last == r)[0]:
            ser = np.sum(ser, axis=int(j), keepdims=True)
            dims[j] = 1

    if ser.size != 1:
        raise RuntimeError(
            "lossn_manjunath: a link was never integrated out; the elimination order is "
            "inconsistent with the constraint rows")
    return float(ser.reshape(-1)[0])


def lossn_manjunath(nu: np.ndarray, A: np.ndarray, C: np.ndarray,
             max_live_states: Optional[int] = None
             ) -> Tuple[np.ndarray, np.ndarray, float, int]:
    """
    Exact normalization constant, carried load and blocking of a loss network.

    Args:
        nu: Offered load of route (class) r (R,), nonnegative.
        A: Capacity requirement of link j for route r (J, R), nonnegative
            integers.
        C: Available capacity of link j (J,), nonnegative integers.
        max_live_states: Cap on the simultaneously live series coefficients;
            defaults to DEFAULT_MAX_LIVE_STATES.

    Returns:
        Tuple (qlen, loss, lG, niter) where:
            - qlen: Mean carried load E[n_r] per route (R,).
            - loss: Blocking probability per route (R,).
            - lG: Log of the EXACT normalization constant g(C).
            - niter: Always 1; the transform is direct.
    """
    nu = np.asarray(nu, dtype=np.float64).ravel()
    A = np.atleast_2d(np.asarray(A, dtype=np.float64))
    C = np.asarray(C, dtype=np.float64).ravel()
    R = len(nu)
    if A.shape != (len(C), R):
        raise ValueError("lossn_manjunath: A must be %dx%d (J x R), got %s"
                         % (len(C), R, A.shape))
    if np.any(nu < 0):
        raise ValueError("lossn_manjunath: nu must be nonnegative")

    if max_live_states is None:
        max_live_states = DEFAULT_MAX_LIVE_STATES

    qlen = np.zeros(R)
    loss = np.zeros(R)

    Ai, Ci = _integralize(A, C)
    J = len(Ci)

    # A route absent from every remaining row never blocks: its marginal is an
    # untruncated Poisson, so it carries its full offered load and factors
    # exp(nu_r) out of g(C).
    if J == 0:
        free = np.ones(R, dtype=bool)
    else:
        free = ~np.any(Ai != 0, axis=0)
    lGfree = float(np.sum(nu[free]))
    qlen[free] = nu[free]

    if J == 0 or np.all(free):
        return qlen, loss, lGfree, 1

    # Per-route truncation, and the terms f_r(n) = nu_r**n / n!. Each sequence is
    # built in the LOG domain and shifted by its own maximum before
    # exponentiating, so its largest entry is exactly 1; the scale cancels in
    # every ratio below and is added back in log g. Shifting before the
    # exponential rather than after is what keeps a heavy route in range:
    # nu**n/n! peaks near n = nu at roughly exp(nu)/sqrt(2 pi nu), so forming the
    # terms first and dividing by their maximum afterwards overflows to inf for a
    # load above about 700, and inf is not finite, so a rescaling step guarded on
    # finiteness would then decline to run and the inf would reach g, returning
    # nan for every metric. Measured at nu = C = 900.
    nmax_full = np.zeros(R, dtype=np.int64)
    f = [np.ones(1) for _ in range(R)]
    logscale = 0.0
    for r in range(R):
        if free[r]:
            continue
        pos = Ai[:, r] > 0
        v = int(np.min(Ci[pos] // Ai[pos, r]))
        nmax_full[r] = v
        n = np.arange(v + 1)
        if nu[r] == 0.0:
            # No offered load: only the empty term survives, and log(0) has no
            # shift to take.
            fr = np.zeros(v + 1)
            fr[0] = 1.0
        else:
            lfr = n * math.log(nu[r]) - special.gammaln(n + 1)
            m = float(np.max(lfr))
            fr = np.exp(lfr - m)
            logscale += m
        f[r] = fr

    peak = [0]
    G = _series(f, Ai, Ci, nmax_full, max_live_states, peak)
    if G <= 0.0:
        raise RuntimeError(
            "lossn_manjunath: the admissible set is empty -- no state satisfies A n <= C, "
            "so the loss network has no stationary distribution")
    lG = math.log(G) + logscale + lGfree

    for r in range(R):
        if free[r]:
            continue
        Cr = Ci - Ai[:, r]
        if np.any(Cr < 0):
            # A single class r call already exceeds a capacity, so the route is
            # blocked in every state, including the empty one.
            loss[r] = 1.0
            qlen[r] = 0.0
            continue
        Gr = _series(f, Ai, Cr, nmax_full, max_live_states, peak)
        # The scale cancels here, which is what lets the terms be normalized.
        ratio = Gr / G
        qlen[r] = nu[r] * ratio
        loss[r] = 1.0 - ratio

    return qlen, loss, lG, 1


__all__ = ['lossn_manjunath']
