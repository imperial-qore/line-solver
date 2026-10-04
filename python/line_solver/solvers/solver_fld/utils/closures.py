"""
Moment closures for the second-order fluid methods of SolverFLD.

The default fluid methods close the moment hierarchy at first order: the drift
of the mean uses ``min(E[X], c)`` in place of ``E[min(X, c)]``, so no second
moment ever enters and the mean is biased wherever a rate function bends. This
module supplies the closures that reinstate it.

Three distinct non-linear rate terms occur in the drift, and each needs its own
closure:

``min(n, c)``
    the station capacity. Closed by the min-normal closure of Guenther,
    Stefanek and Bradley (EPEW/UKPEW 2012, LNCS 7587:32-47, eq. 4).

``w_j*X_j / sum_m w_m*X_m``
    the PS/DPS capacity share, a RATIO of populations. Closed at second order
    by the delta method; the correction is exactly capacity-conserving.

``w_r / sum_{j backlogged} w_j``
    the GPS capacity share, a function of the backlog INDICATOR rather than of
    the populations. Closed by enumerating the ``2^K`` backlog patterns.

Port of ``fluid_min_closure.m``, ``fluid_capacity_closure.m``,
``fluid_lld_scaling.m``, ``fluid_share_closure.m`` and ``fluid_gps_share.m``.
"""

import numpy as np
from math import erfc, sqrt, pi, exp, isinf

from ....constants import GlobalConstants


def lld_scaling(lldrow, n):
    """Limited load-dependent rate scaling at a CONTINUOUS population.

    ``sn.lldscaling[i]`` tabulates the station rate multiplier at integer
    populations 1..lldlimit and the discrete solvers read it as
    ``lldscaling[i, min(n, lldlimit)]``. The fluid state is continuous, so the
    table is read by linear interpolation between consecutive entries, clamped
    to the first entry below n=1 and to the last above the table end, matching
    the clamping the CTMC already applies.

    Returns ``(a, da)``: the scaling and its derivative, zero on the clamps.
    """
    n = np.atleast_1d(np.asarray(n, dtype=float))
    if lldrow is None or len(lldrow) == 0:
        return np.ones_like(n), np.zeros_like(n)

    lldrow = np.asarray(lldrow, dtype=float).ravel()
    L = len(lldrow)
    a = np.zeros_like(n)
    da = np.zeros_like(n)

    lo = n <= 1
    a[lo] = lldrow[0]
    hi = n >= L
    a[hi] = lldrow[L - 1]

    mid = ~lo & ~hi
    if np.any(mid):
        k = np.floor(n[mid]).astype(int)  # 1 <= k <= L-1
        frac = n[mid] - k
        a[mid] = lldrow[k - 1] * (1 - frac) + lldrow[k] * frac
        da[mid] = lldrow[k] - lldrow[k - 1]
    return a, da


def _min_closure_scalar(n, c, s2, vc, cov_nc):
    """The one-point case of `min_closure`, without the array machinery.

    THE SAME EXPRESSIONS IN THE SAME ORDER as the vector path below. It exists
    because the fluid ODE evaluates this ONCE PER STATION PER RIGHT-HAND SIDE on
    SCALARS, and the vector path spends that call on five `np.full`s, two masks
    and a list comprehension rather than on the closure: measured on
    mqn_singleserver_ps, whose fixed point needs 20 closure passes of 200
    windows each, `min_closure` and its caller were 37% of the whole solve,
    which ran ~80x slower than the JAR on the identical algorithm and could not
    finish inside the parity row's budget at all.

    NOT BIT-IDENTICAL TO THE VECTOR PATH, and that is a property of `np.exp`
    rather than of this arithmetic: numpy evaluates `exp` over an array with its
    SIMD kernel and libm evaluates the scalar, and the two disagree in the last
    ulp on ~5% of arguments. The vector path is not self-consistent either for
    the same reason -- its result depends on the length of the array it is
    handed. The degenerate branch and `dh` on both branches ARE exact; `h`
    agrees to 5.3e-14 relative over a 40k sweep, and that worst case is
    cancellation at a value of order 1e-89. See
    `python/tests/test_fld_min_closure_scalar.py`.
    """
    th2 = s2 - 2 * cov_nc + vc
    if th2 < 0:
        th2 = 0.0  # beyond Cauchy-Schwarz is inadmissible
    # the kink band, and why it is a band rather than a strict test, is the
    # cross-codebase requirement documented on the vector path below
    if not (th2 > 0) or isinf(c):
        return (min(n, c),
                1.0 if (c - n) > GlobalConstants.FineTol * max(1.0, n) else 0.0,
                0.0)
    th = sqrt(th2)
    d = (n - c) / th
    Phid = 0.5 * erfc(-d / sqrt(2))
    phid = exp(-0.5 * d ** 2) / sqrt(2 * pi)
    h = n * (1 - Phid) + c * Phid - th * phid
    dh = 1 - Phid
    # min() is piecewise linear, so its second derivative is carried entirely by
    # the atom at X = Y; smoothing over the normal marginal turns that atom into
    # the density. It stays zero on the degenerate branch, which smooths nothing.
    d2h = -phid / th
    # A POPULATION IS NONNEGATIVE AND THE NORMAL MARGINAL IS NOT. For X >= 0
    # pathwise min(X,c) >= 0, and min() being concave Jensen puts
    # E[min(X,c)] <= min(E[X],c), so the value belongs to [0, min(n,c)]. The
    # normal marginal has no such support, and the mass it places below zero drags
    # the expectation out of that range once the mean falls to about one standard
    # deviation: n = 0, c = 1, th = 0.664 returns -0.019, a station that CREATES
    # work. Project onto the admissible range and take the derivatives of the
    # bound that binds, so the Jacobian still matches h.
    hi = min(n, c)
    at_lo = h < 0
    at_hi = h > hi
    h = min(max(h, 0.0), hi)
    if at_lo:
        dh = 0.0
    elif at_hi:
        dh = 1.0 if n < c else 0.0
    if at_lo or at_hi:
        d2h = 0.0
    return h, dh, d2h


def min_closure(n, c, s2=0.0, vc=0.0, cov_nc=0.0, want_d2=False):
    """Min-normal moment closure of ``E[min(X, Y)]`` for jointly normal X, Y.

    Implemented in the general two-population form of Guenther-Stefanek-Bradley
    eq. (4)::

        E[min(X,Y)] = E[X]*Phi((E[Y]-E[X])/th) + E[Y]*Phi((E[X]-E[Y])/th)
                      - th*phi((E[Y]-E[X])/th)
        th = sqrt(Var[X] - 2*Cov[X,Y] + Var[Y])

    The derivative with respect to ``E[X]`` is ``P(X < Y)``. SolverFLD only ever
    needs the specialisation ``Y = c`` deterministic, but the general arguments
    are kept so this IS the published closure rather than one instance of it.

    With ``th == 0`` the expressions collapse to ``min(n, c)`` and to the
    indicator ``1{n < c}``, so the first-order closure is recovered exactly and
    callers share a single code path. The derivative at the kink is taken as 0,
    the right derivative of ``min()``.
    """
    if (np.ndim(n) == 0 and np.ndim(c) == 0 and np.ndim(s2) == 0
            and np.ndim(vc) == 0 and np.ndim(cov_nc) == 0):
        h, dh, d2h = _min_closure_scalar(float(n), float(c), float(s2),
                                         float(vc), float(cov_nc))
        if want_d2:
            return np.array([h]), np.array([dh]), np.array([d2h])
        return np.array([h]), np.array([dh])

    n = np.atleast_1d(np.asarray(n, dtype=float))
    c = np.atleast_1d(np.asarray(c, dtype=float))
    s2 = np.atleast_1d(np.asarray(s2, dtype=float))
    vc = np.atleast_1d(np.asarray(vc, dtype=float))
    cov_nc = np.atleast_1d(np.asarray(cov_nc, dtype=float))
    sz = max(n.size, c.size, s2.size, vc.size, cov_nc.size)
    if n.size == 1:
        n = np.full(sz, n[0])
    if c.size == 1:
        c = np.full(sz, c[0])
    if s2.size == 1:
        s2 = np.full(sz, s2[0])
    if vc.size == 1:
        vc = np.full(sz, vc[0])
    if cov_nc.size == 1:
        cov_nc = np.full(sz, cov_nc[0])

    th2 = s2 - 2 * cov_nc + vc
    th2 = np.where(th2 < 0, 0.0, th2)  # beyond Cauchy-Schwarz is inadmissible

    h = np.zeros(sz)
    dh = np.zeros(sz)
    d2h = np.zeros(sz)

    # THE INDICATOR CARRIES A BAND, and it is a cross-codebase requirement rather
    # than a modelling choice. A saturated fluid fixed point sits exactly AT n = c,
    # and each engine's ODE stops on its own residual: MATLAB lands at 1.0004 and
    # the C++ port at 1 - 1.8e-13 on the same model, so a strict ``n < c`` reads
    # saturated in one and unsaturated in the other. That flips a whole Jacobian
    # row between zero and unit, and with it the hyperbolicity verdict that decides
    # whether SolverFLD answers with 'minnormal' or falls back to the first-order
    # method. A population within FineTol of the server count IS at the kink.
    deg = ~(th2 > 0) | np.isinf(c)
    h[deg] = np.minimum(n[deg], c[deg])
    dh[deg] = ((c[deg] - n[deg]) >
               GlobalConstants.FineTol * np.maximum(1.0, n[deg])).astype(float)

    sm = ~deg
    if np.any(sm):
        th = np.sqrt(th2[sm])
        d = (n[sm] - c[sm]) / th
        Phid = 0.5 * np.array([erfc(-di / sqrt(2)) for di in np.atleast_1d(d)])
        phid = np.exp(-0.5 * d ** 2) / sqrt(2 * pi)
        hsm = n[sm] * (1 - Phid) + c[sm] * Phid - th * phid
        dhsm = 1 - Phid
        d2hsm = -phid / th
        # the projection onto [0, min(n,c)], and the derivatives of the bound that
        # binds; see the scalar kernel above for why the normal marginal needs it
        hi = np.minimum(n[sm], c[sm])
        at_lo = hsm < 0
        at_hi = hsm > hi
        hsm = np.minimum(np.maximum(hsm, 0.0), hi)
        dhsm = np.where(at_lo, 0.0, dhsm)
        dhsm = np.where(at_hi, (n[sm] < c[sm]).astype(float), dhsm)
        d2hsm = np.where(at_lo | at_hi, 0.0, d2hsm)
        h[sm] = hsm
        dh[sm] = dhsm
        d2h[sm] = d2hsm
    if want_d2:
        return h, dh, d2h
    return h, dh


def capacity_closure(n, c, s2=0.0, lldrow=None, is_inf=False, want_d2=False):
    """Moment closure of the station capacity term ``psi(X)`` and its derivative.

    Every scheduling branch scales the coordinates of a station by
    ``psi(n_i)/n_i``, where ``psi`` is how much work the station clears::

        psi(n) = min(n, c) * alpha(n)   at a queueing station
        psi(n) = n         * alpha(n)   at an infinite server

    Returns ``E[psi(X)]`` for ``X ~ Normal(n, s2)`` with its derivative.
    ``s2 = 0`` gives ``psi(n)`` itself.

    With a tabulated ``alpha`` the scaling is piecewise linear on the integer
    lattice, so ``psi`` is piecewise QUADRATIC with breakpoints at the integers
    and at ``c``, and the expectation is integrated segment by segment against
    the normal using the truncated moments M0, M1, M2.

    The integration must be EXACT, not a fixed quadrature rule: Gauss-Hermite
    applied to a piecewise-linear integrand does not smooth its kinks, it
    relocates them, and the resulting estimate is itself piecewise linear, so
    its second derivative vanishes and the 1/N refinement silently returns a
    null correction.
    """
    n = float(np.asarray(n).ravel()[0])
    has_lld = lldrow is not None and len(lldrow) > 0

    if not has_lld:
        if is_inf:
            return (n, 1.0, 0.0) if want_d2 else (n, 1.0)
        # the scalar kernel rather than min_closure: this is the ODE's hot path
        # and the two agree to the bit, so the round trip through two 1-element
        # arrays buys nothing
        h, dh, d2h = _min_closure_scalar(n, float(c), float(s2), 0.0, 0.0)
        # A population cannot be negative, but the normal marginal puts mass
        # below zero, and there min(X,c) = X < 0. Once n is small against the
        # standard deviation that drives E[min(X,c)] itself negative, giving
        # negative service rates and mass destroyed by the non-negativity clamp
        # of the ODE solver. A floor at zero is the minimal repair. Subtracting
        # the whole negative tail instead is more self-consistent but departs
        # from the published closure everywhere rather than only where it is
        # unusable, and measurably lost accuracy in the mid-load range.
        if h < 0:
            return (0.0, 0.0, 0.0) if want_d2 else (0.0, 0.0)
        return (h, dh, d2h) if want_d2 else (h, dh)

    lldrow = np.asarray(lldrow, dtype=float).ravel()
    L = len(lldrow)
    if not (s2 > 0):
        ph, pdh, pd2h = _psi(n, c, lldrow, L, is_inf)
        return (ph, pdh, pd2h) if want_d2 else (ph, pdh)

    s = sqrt(s2)
    bps = list(range(0, L + 1))
    if (not is_inf) and np.isfinite(c) and c > L:
        bps.append(float(c))
    bps = sorted(set(bps))

    h = 0.0
    dh = 0.0
    d2h = 0.0
    # psi is extended by zero below the first breakpoint, so psi' jumps there too
    # and that atom belongs in psi'' exactly like the interior ones
    b_prev = 0.0
    c_prev = 0.0
    for k in range(len(bps)):
        p = float(bps[k])
        q = float(bps[k + 1]) if k < len(bps) - 1 else np.inf
        A, B, C = _segment(p, q, c, lldrow, L, is_inf)
        M0, M1, M2 = _moments(p, q, n, s)
        h += A * M0 + B * M1 + C * M2
        dh += B * M0 + 2 * C * M1
        d2h += 2 * C * M0
        # psi' jumps across this breakpoint, so psi'' carries an atom there; the
        # segment sum above sees only the quadratic part and would miss it
        jump = (B + 2 * C * p) - (b_prev + 2 * c_prev * p)
        d2h += jump * exp(-0.5 * ((p - n) / s) ** 2) / (s * sqrt(2 * pi))
        b_prev, c_prev = B, C
    return (h, dh, d2h) if want_d2 else (h, dh)


def _segment(p, q, c, lldrow, L, is_inf):
    """psi(u) = A + B*u + C*u^2 on [p, q], from base(u)*alpha(u)."""
    if p >= L:
        a0, a1 = lldrow[L - 1], 0.0
    elif p < 1:
        a0, a1 = lldrow[0], 0.0
    else:
        k = int(np.floor(p))
        a1 = lldrow[k] - lldrow[k - 1]
        a0 = lldrow[k - 1] - a1 * k
    # base(u) is u below the saturation point and c above it; the breakpoints
    # guarantee the segment lies entirely on one side
    if is_inf or (not np.isfinite(c)) or q <= c:
        return 0.0, a0, a1
    return c * a0, c * a1, 0.0


def _moments(p, q, n, s):
    """Truncated moments E[X^j * 1{p < X < q}] for X ~ Normal(n, s^2)."""
    zp = (p - n) / s
    Pp = 0.5 * erfc(-zp / sqrt(2))
    pp = exp(-0.5 * zp * zp) / sqrt(2 * pi)
    if np.isinf(q):
        Pq, pq = 1.0, 0.0
    else:
        zq = (q - n) / s
        Pq = 0.5 * erfc(-zq / sqrt(2))
        pq = exp(-0.5 * zq * zq) / sqrt(2 * pi)
    M0 = Pq - Pp
    M1 = n * M0 + s * (pp - pq)
    if np.isinf(q):
        M2 = (n * n + s * s) * M0 + s * ((p + n) * pp)
    else:
        M2 = (n * n + s * s) * M0 + s * ((p + n) * pp - (q + n) * pq)
    return M0, M1, M2


def _psi(u, c, lldrow, L, is_inf):
    """psi and its derivative at a continuous population, zero below 0."""
    a, da = lld_scaling(lldrow, u)
    a = float(a[0])
    da = float(da[0])
    if is_inf:
        base, dbase = u, 1.0
    else:
        base = min(u, c)
        dbase = 1.0 if u < c else 0.0
    p = base * a
    dp = dbase * a + base * da
    # base and alpha are both piecewise linear, so psi is piecewise quadratic and
    # psi'' = 2*base'*alpha' away from the breakpoints; the atoms AT them are not
    # representable without a marginal, and this first-order branch smooths none.
    d2p = 2.0 * dbase * da
    if u <= 0:
        return 0.0, 0.0, 0.0
    return p, dp, d2p


def _expansion_weight(ratio):
    """How much of the second-order correction the series admits, and d tau/d ratio.

    Every second-order term of the share closure is a term of the series for
    ``E[1/v]``, whose successive terms are in the ratio ``Var(v)/v^2``, so the
    truncation is meaningful below 1 and the terms GROW above it. Nothing in the
    algebra notices: at a near-empty station the corrections come back larger than
    the quantity they correct, and the drift that follows is not integrable.

    One on [0, 1], zero from 4 up, and the C^1 smoothstep between. Both ends
    matter. The lower one has to be EXACTLY one on the whole convergent region, so
    every model already inside it is bit-identical; the upper one has to be reached
    with a vanishing derivative, because the drift is integrated and a kink in it
    is what collapses the step size. The thresholds are the series, not a tuning:
    at ratio 1 successive terms stop shrinking, and at ratio 4 the standard
    deviation of v is twice its mean, where a non-negative v has essentially no
    mass near the point being expanded about.
    """
    lo, hi = 1.0, 4.0
    if ratio <= lo:
        return 1.0, 0.0
    if ratio >= hi:
        return 0.0, 0.0
    t = (ratio - lo) / (hi - lo)
    return 1.0 - t * t * (3.0 - 2.0 * t), -6.0 * t * (1.0 - t) / (hi - lo)


def project_rate(r, xb, capped, tot):
    """Project a jointly closed per-coordinate service share onto its own set.

    That set is ``r >= 0``, ``r <= xb`` where the bound applies, and
    ``sum(r) == tot``.

    THE JOINT CLOSURE IS AN EXPANSION AND CAN LEAVE IT. ``r = s*psi +
    psi'*Cov(S,N)`` adds a term that sums to ZERO over the coordinates, so it
    moves mass between them and its entries can push one past either bound; the
    first-order share ``x_j/n_i*psi`` cannot, being x_j scaled by
    ``psi/n_i <= 1``. Either breach ends the same way, because the integrator
    holds every coordinate non-negative: ``r_j > x_j`` drains coordinate j faster
    than it holds, the state goes negative and the clamp INJECTS mass.

    THE UPPER BOUND HOLDS ONLY WITHOUT LOAD DEPENDENCE, which is what CAPPED
    selects: r is an expected NUMBER in service so ``r_j <= x_j``, but
    ``psi(n) = min(n,c)*alpha(n)`` folds the load-dependent scaling into the same
    variable, and with alpha > 1 the first-order share itself exceeds x_j.

    Clip, then move the residual onto the coordinates that still have slack in
    proportion to it, so ``sum(r) == tot`` survives and the station still clears
    what its capacity closure says it clears. A NO-OP whenever the expansion
    stayed inside the set, which is why models already inside it are
    bit-identical.
    """
    r = np.asarray(r, dtype=float).ravel().copy()
    xb = np.asarray(xb, dtype=float).ravel()
    zt = GlobalConstants.Zero
    if (r >= -zt).all() and ((not capped) or (r <= xb + zt).all()):
        return r
    r = np.maximum(r, 0.0)
    if capped:
        r = np.minimum(r, xb)
    for _ in range(len(r) + 1):
        d = tot - float(r.sum())
        if abs(d) <= zt:
            break
        if d > 0:
            slack = (xb - r) if capped else np.ones(len(r))
        else:
            slack = r
        tsl = float(slack.sum())
        if tsl <= zt:
            break
        r = np.maximum(r + d * slack / tsl, 0.0)
        if capped:
            r = np.minimum(r, xb)
    return r


_SHARE_SCALAR_MAX = 8


def _share_closure_scalar(x, wv, C, n, want_cov=False):
    """The no-Jacobian case of `share_closure` at a SMALL station, in plain floats.

    THE SAME EXPRESSIONS IN THE SAME ASSOCIATION as the array path below,
    including the two matrix-vector products, which are spelled out so that
    ``cuv`` sums over the ROW of C and ``cvv`` over the COLUMN exactly as
    ``C @ wv`` and ``wv @ C @ wv`` do.

    NOT BIT-IDENTICAL, for the same class of reason as `_min_closure_scalar`:
    every sum here is sequential, while numpy's `add.reduce` and its matmul
    accumulate in an unrolled/pairwise ORDER that depends on the length. Over a
    60k sweep covering all five branches the shares agree to 2.8e-14 ABSOLUTE on
    quantities that sum to one (1.2e-11 relative, reached only where a large
    covariance makes the delta-method correction cancel), and no coordinate ever
    changed sign -- which is the part the branches below actually read. That is
    eleven orders below the 1e-4 the ODE is integrated to. See
    `python/tests/test_fld_min_closure_scalar.py`.

    Why it exists: the drift calls this once per sharing station per right-hand
    side, on an array carrying ONE COORDINATE PER CLASS -- length 2 on
    mqn_singleserver_ps -- where the eighteen numpy calls below cost several
    microseconds each and the arithmetic itself costs tens of nanoseconds. It was
    29% of that model's fluid solve after `min_closure` was fixed.
    """
    v = 0.0
    u = [0.0] * n
    for j in range(n):
        u[j] = wv[j] * x[j]
        v += u[j]
    if v <= 0:
        return (np.zeros(n), np.zeros(n)) if want_cov else np.zeros(n)

    s = [u[j] / v for j in range(n)]
    if C is None:
        return (np.array(s), np.zeros(n)) if want_cov else np.array(s)

    nz = False
    for j in range(n):
        for m in range(n):
            if C[j][m]:
                nz = True
                break
        if nz:
            break
    if not nz:
        return (np.array(s), np.zeros(n)) if want_cov else np.array(s)

    v2 = v ** 2
    v3 = v ** 3
    cvv = 0.0
    for j in range(n):
        acc = 0.0
        for m in range(n):
            acc += wv[m] * C[m][j]
        cvv += acc * wv[j]
    # how far into the series this point sits; see `_expansion_weight`
    tau, _ = _expansion_weight(cvv / v2)
    cn = np.zeros(n)
    if want_cov and tau > 0.0:
        # Cov(S_j, N) at the means, from grad(S_j)'*C*1
        cvn = 0.0
        an = [0.0] * n
        for j in range(n):
            acc = 0.0
            for m in range(n):
                acc += C[j][m]
            an[j] = acc
            cvn += wv[j] * acc
        for j in range(n):
            cn[j] = tau * (wv[j] * an[j] / v - u[j] * (cvn / v2))
    if tau <= 0.0:
        return (np.array(s), cn) if want_cov else np.array(s)
    for j in range(n):
        acc = 0.0
        for m in range(n):
            acc += C[j][m] * wv[m]
        s[j] = s[j] + tau * (-(wv[j] * acc) / v2 + (u[j] * cvv) / v3)

    if all(sj >= -GlobalConstants.Zero for sj in s):
        sc = np.array([0.0 if sj < 0 else sj for sj in s])
        return (sc, cn) if want_cov else sc

    # the expansion has left its region of validity for at least one coordinate;
    # clip and renormalise so the shares still sum to one
    T = 0.0
    any_act = False
    for sj in s:
        if sj > 0:
            T += sj
            any_act = True
    if not any_act:
        sc = np.array([u[j] / v for j in range(n)])
        return (sc, cn) if want_cov else sc
    sc = np.array([s[j] / T if s[j] > 0 else 0.0 for j in range(n)])
    return (sc, cn) if want_cov else sc


def _share_out(s, ds, cn, dcn, want_jac, want_cov):
    """The return shape of `share_closure`: s, then ds, cn, dcn as requested."""
    if want_jac and want_cov:
        return s, ds, cn, dcn
    if want_jac:
        return s, ds
    if want_cov:
        return s, cn
    return s


def share_closure(x, wv, C=None, want_jac=False, want_cov=False):
    """Second-order closure of the capacity share of a sharing discipline.

    A DPS station gives coordinate ``j`` the fraction
    ``S_j = w_j*X_j / sum_m w_m*X_m`` of its capacity, and PS is the same with
    unit weights. The first-order closure evaluates that ratio at the mean,
    which is not ``E[S_j]``: the map is a ratio, so Jensen biases it towards the
    coordinates carrying the LARGER weight. With ``u_j = w_j*X_j`` and
    ``v = sum_m u_m`` the delta method gives::

        E[u_j/v] = mu_j/v - Cov(u_j,v)/v^2 + mu_j*Var(v)/v^3 + O(sigma^3)

    INVARIANT to check on any change here: the correction is exactly
    capacity-conserving, because ``sum_j Cov(u_j,v) = Var(v)`` makes the two
    correction terms cancel in the sum, so ``sum_j S_j == 1`` as a
    work-conserving discipline requires.

    The expansion is local and fails when ``v`` is small against its own
    standard deviation. There a raw share can come out negative; it is clipped
    at zero and the survivors renormalised, which preserves conservation.
    """
    # The reductions below are spelled as METHODS rather than as np.sum/np.any/
    # np.all throughout: same reduction on the same buffer, but without the
    # np.* dispatch wrapper, which on the arrays this sees -- one coordinate per
    # class at ONE station, so length 2 here -- costs several times the
    # arithmetic. This is the fluid ODE's second hot spot after min_closure; see
    # the note on _min_closure_scalar for the measurement that found both.
    x = np.asarray(x, dtype=float).ravel()
    wv = np.asarray(wv, dtype=float).ravel()
    n = len(x)
    if not want_jac and n <= _SHARE_SCALAR_MAX:
        return _share_closure_scalar(
            x.tolist(), wv.tolist(),
            None if C is None else np.asarray(C, dtype=float).tolist(), n, want_cov)
    u = wv * x
    v = float(u.sum())
    cn = np.zeros(n) if want_cov else None
    dcn = np.zeros((n, n)) if (want_cov and want_jac) else None
    if v <= 0:
        return _share_out(np.zeros(n), np.zeros((n, n)) if want_jac else None,
                          cn, dcn, want_jac, want_cov)

    s = u / v
    ds = (np.diag(wv) / v - np.outer(u / v ** 2, wv)) if want_jac else None

    if C is None:
        return _share_out(s, ds, cn, dcn, want_jac, want_cov)
    C = np.asarray(C, dtype=float)
    if not C.any():
        return _share_out(s, ds, cn, dcn, want_jac, want_cov)
    cuv = wv * (C @ wv)      # Cov(u_j, v)
    cvv = float(wv @ C @ wv)  # Var(v)

    # How far into the series this point sits, and how much of the second-order
    # correction that leaves admissible. tau depends on x through v alone, since C
    # is held fixed here exactly as the Jacobian is.
    tau, dtau_dratio = _expansion_weight(cvv / v ** 2)
    if tau <= 0 and not want_jac:
        return _share_out(s, ds, cn, dcn, want_jac, want_cov)
    dtau = (dtau_dratio * (-2.0 * cvv / v ** 3)) * wv   # d tau / d x_m

    if want_cov:
        # Cov(S_j, N) at the means, from grad(S_j)'*C*1: S_j = w_j X_j / V, so
        # dS_j/dX_a = w_j*delta_aj/V - u_j*w_a/V^2 and the two pieces contract
        # against C*1 = Cov(X,N) and w'*C*1 = Cov(V,N).
        an = C @ np.ones(n)
        cvn = float(wv @ an)
        cn0 = wv * an / v - u * (cvn / v ** 2)
        cn = tau * cn0
        if want_jac:
            dcn = (tau * (-np.outer(wv * an, wv / v ** 2)
                          - (cvn / v ** 2) * np.diag(wv)
                          + (2.0 * cvn / v ** 3) * np.outer(u, wv))
                   + np.outer(cn0, dtau))

    scorr = -cuv / v ** 2 + (u * cvv) / v ** 3
    s = s + tau * scorr
    if want_jac:
        ds = (ds + tau * ((2.0 / v ** 3) * np.outer(cuv, wv)
                          + (cvv / v ** 3) * np.diag(wv)
                          - (3.0 * cvv / v ** 4) * np.outer(u, wv))
              + np.outer(scorr, dtau))

    # A COORDINATE CARRYING NO MASS MUST NOT DECIDE THE CLIP. With u_j = 0 the
    # plug-in share is zero and the correction leaves s_j = -Cov(u_j,v)/v^2, a
    # quantity of the order of rounding whose SIGN is not meaningful. Letting it
    # select the branch below zeroes that coordinate's whole Jacobian row, and a
    # zero row is an exact zero eigenvalue: the Lyapunov step then reads the fixed
    # point as non-hyperbolic and SolverFLD silently drops from 'minnormal' to the
    # first-order method. Clip only a share that is negative BEYOND the numerical
    # zero.
    if (s >= -GlobalConstants.Zero).all():
        s = np.where(s < 0, 0.0, s)
        return _share_out(s, ds, cn, dcn, want_jac, want_cov)

    # the expansion has left its region of validity for at least one
    # coordinate; clip and renormalise so the shares still sum to one
    act = s > 0
    if not act.any():
        s = u / v
        if want_jac:
            ds = np.diag(wv) / v - np.outer(u / v ** 2, wv)
        return _share_out(s, ds, cn, dcn, want_jac, want_cov)
    T = float(s[act].sum())
    snew = np.zeros(n)
    snew[act] = s[act] / T
    if want_jac:
        dT = np.sum(ds[act, :], axis=0)
        dsnew = np.zeros((n, n))
        dsnew[act, :] = ds[act, :] / T - np.outer(s[act] / T ** 2, dT)
        ds = dsnew
    return _share_out(snew, ds, cn, dcn, want_jac, want_cov)


def gps_share(xk, wk, vk, want_jac=False):
    """Expected capacity share of a GPS station under a normal marginal.

    GPS divides the server by WEIGHT among the classes that are BACKLOGGED, then
    equally among that class's own jobs. The share therefore depends on the
    backlog INDICATOR vector, not on the populations, and that is what makes GPS
    unreachable for a first-order closure: with continuous ``x_k > 0`` every
    class is always backlogged, the indicator is identically one, and the share
    collapses to the constant ``w_r/sum_j w_j`` regardless of load. For GPS the
    second moment is not a correction, it is the entire mechanism.

    The closure is an EXACT enumeration rather than an expansion: the share is
    piecewise CONSTANT over the ``2^K`` backlog patterns, so

        E[S_r] = sum_{A ni r} P(backlog set = A) * w_r / sum_{j in A} w_j

    carries no truncation error once the pattern probabilities are given. Those
    come from the marginals, ``P(N_k >= 1) = Phi((x_k - 1/2)/sigma_k)`` with the
    continuity correction for an integer population, multiplied as if the
    backlogs were independent. That independence is the one approximation here
    and it is not innocuous: in a closed network the station coordinates are
    NEGATIVELY correlated through population conservation.

    The empty pattern contributes zero share, so ``sum_r E[S_r]`` is
    ``1 - P(all classes empty)`` rather than 1. That is deliberate: GPS is
    single-server, so the backlog indicator plays the role ``min(n, c)`` plays
    at a PS station and no separate capacity term is applied by the caller.

    Unlike the DPS ratio closure this is NOT perturbative in sigma: as
    ``sigma_k -> 0`` the probability tends to a step at ``x_k = 1/2`` and its
    derivative ``phi(.)/sigma_k`` diverges, so the Jacobian stiffens at low
    variance.
    """
    xk = np.asarray(xk, dtype=float).ravel()
    wk = np.asarray(wk, dtype=float).ravel()
    vk = np.asarray(vk, dtype=float).ravel()
    K = len(xk)
    s = np.zeros(K)
    ds = np.zeros((K, K))

    # the enumeration is 2^K, so refuse rather than crawl; K here is the number
    # of classes AT ONE STATION, small in every practical model
    if K > 12:
        raise ValueError(
            "GPS closes its capacity share by enumerating the 2^K backlog patterns of a "
            "station, and this station carries %d classes. Above 12 the enumeration is no "
            "longer tractable; use a DPS station with method='closing' instead." % K)

    sw = float(np.sum(wk))
    if sw <= 0:
        return (s, ds) if want_jac else s
    wk = wk / sw

    # P(class k backlogged) = P(N_k >= 1) for an INTEGER population, so the
    # normal approximation needs the continuity correction P(N_k > 1/2);
    # thresholding at 1 instead makes the share vanish for every class with
    # x_k < 1 at sigma = 0, which stalls the server completely and is an
    # ABSORBING state for the ODE. The sigma = 0 fallback is the fluid limit
    # x_k > 0, matching the mean-field statement that positive mass is
    # backlogged.
    p = np.zeros(K)
    dp = np.zeros(K)
    for k in range(K):
        if vk[k] > 0:
            sd = sqrt(vk[k])
            z = (xk[k] - 0.5) / sd
            p[k] = 0.5 * erfc(-z / sqrt(2))
            dp[k] = exp(-0.5 * z * z) / (sqrt(2 * pi) * sd)
        else:
            p[k] = 1.0 if xk[k] > 0 else 0.0

    dsdp = np.zeros((K, K))
    for mask in range(1, 1 << K):
        A = np.array([(mask >> k) & 1 for k in range(K)], dtype=bool)
        W = float(np.sum(wk[A]))
        if W <= 0:
            continue  # every backlogged class in this pattern carries zero weight
        shareA = wk[A] / W
        q = p.copy()
        q[~A] = 1 - p[~A]  # per-class factor of P(A)
        s[A] += float(np.prod(q)) * shareA
        if want_jac:
            idxA = np.where(A)[0]
            for m in range(K):
                qm = q.copy()
                qm[m] = 1.0  # product over j != m
                dsdp[idxA, m] += (2 * int(A[m]) - 1) * float(np.prod(qm)) * shareA

    if want_jac:
        ds = dsdp * dp[None, :]  # chain rule through p_m = Phi((x_m-1/2)/sigma_m)
        return s, ds
    return s
