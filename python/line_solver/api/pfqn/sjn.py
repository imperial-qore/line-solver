"""Shortest-job-next stations in mean value analysis (Kant 1992).

Closed queueing networks in which a subset of the single-server stations schedules
non-preemptively by shortest job next (SJN/SJF), the job size being known on arrival.

The station is carried by the conditional waiting time ``W(x,n)`` of a tagged customer whose
service requirement is ``x``, obtained from the arrival theorem as the sum of the residual life of
the job in service, the work of the queued jobs that will be served before the tagged one, and the
work of the jobs that overtake it while it waits::

    W(x,n) = [ (1+CV^2) s U(n-1)/2 + X(n-1) phi(x,n-1) ] / [ 1 - X(n-1) theta(x) ]
    theta(x) = int_0^x t f(t) dt,   phi(x,n) = int_0^x W(t,n) t f(t) dt
    R(n)     = s + int_0^inf W(x,n) f(x) dx

:func:`pfqn_mvasjn` steps the whole population lattice; :func:`pfqn_amvasjn` replaces it with a
Schweitzer fixed point on the size-resolved queue length. Ported at parity from MATLAB
``pfqn_mvasjn.m`` / ``pfqn_amvasjn.m`` and their ``private/sjn_*.m`` helpers.

Reference: K. Kant, "MVA approximations for SJN scheduling", Performance Evaluation 15(1):41-61,
1992.
"""

import numpy as np
from scipy.special import gammainc, gammaincc, gammaln

__all__ = ['pfqn_mvasjn', 'pfqn_amvasjn', 'SjnOptions', 'SjnStarvationError']


class SjnOptions:
    """Options shared by the two shortest-job-next solvers."""

    def __init__(self, ns=32, lfactor=8.0, prio=None, tol=1e-8, iter_max=1000, umax=0.999):
        #: number of grid subdivisions of the job size axis, even
        self.ns = int(ns)
        #: grid extent, in units of the largest mean service time at the station
        self.lfactor = float(lfactor)
        #: priority levels, one per class, lower is higher priority; None pools the classes
        self.prio = None if prio is None else np.asarray(prio, dtype=int).ravel()
        #: convergence tolerance of the fixed point
        self.tol = float(tol)
        #: iteration cap of the fixed point
        self.iter_max = int(iter_max)
        #: utilization cap at an SJN station, strictly below one
        self.umax = float(umax)

    def validate(self, R):
        if self.ns % 2 != 0:
            raise ValueError('ns must be even, composite Simpson integrates over panels of two '
                             'subdivisions')
        if not (0 < self.umax < 1):
            raise ValueError('umax must lie strictly between zero and one')
        if self.prio is not None:
            if self.prio.size != R:
                raise ValueError('prio must have one priority level per class')
            if np.unique(self.prio).size != R:
                raise ValueError('prio must assign distinct levels, ties across classes are not '
                                 'covered by the SJN priority equations')


class _Fit:
    """Erlang mixture fitted to a mean and a squared coefficient of variation."""

    __slots__ = ('w', 'k', 'mu')

    def __init__(self, w, k, mu):
        self.w = np.asarray(w, dtype=float)
        self.k = np.asarray(k, dtype=int)
        self.mu = np.asarray(mu, dtype=float)

    @property
    def empty(self):
        return self.w.size == 0


def _fit(s, cv2):
    """Two-moment fit: branching Erlang below CV^2 = 1, balanced-means H2 above it.

    The mixture form is what makes theta(x) and the tail integrals closed form.
    """
    if s <= 0:
        return _Fit([], [], [])
    if cv2 < 0:
        raise ValueError('negative squared coefficient of variation')
    if abs(cv2 - 1.0) < 1e-8:
        return _Fit([1.0], [1], [1.0 / s])
    if cv2 < 1:
        k = int(np.ceil(1.0 / cv2))
        p = (k * cv2 - np.sqrt(k * (1 + cv2) - k * k * cv2)) / (1 + cv2)
        mu = (k - p) / s
        return _Fit([p, 1 - p], [k - 1, k], [mu, mu])
    p = 0.5 * (1 + np.sqrt((cv2 - 1) / (cv2 + 1)))
    return _Fit([p, 1 - p], [1, 1], [2 * p / s, 2 * (1 - p) / s])


def _pdf(f, x):
    y = np.zeros_like(x, dtype=float)
    xs = np.maximum(x, np.finfo(float).tiny)
    for j in range(f.w.size):
        k = f.k[j]
        mu = f.mu[j]
        y += f.w[j] * np.exp(k * np.log(mu) + (k - 1) * np.log(xs) - mu * x - gammaln(k))
    return y


def _theta(f, x):
    """The primitive int_0^x t f(t) dt, in closed form."""
    y = np.zeros_like(x, dtype=float)
    for j in range(f.w.size):
        k = f.k[j]
        mu = f.mu[j]
        y += f.w[j] * (k / mu) * gammainc(k + 1, mu * x)
    return y


def _ccdf(f, x):
    y = 0.0
    for j in range(f.w.size):
        y += f.w[j] * gammaincc(f.k[j], f.mu[j] * x)
    return float(y)


def _tailmom(f, lx, c, order):
    """int_lx^inf t^order exp(-c (t-lx)) f(t) dt, in logarithms so exp(c lx) cannot overflow."""
    y = 0.0
    for j in range(f.w.size):
        k = f.k[j]
        mu = f.mu[j]
        rate = mu + c
        g = gammaincc(k + order, rate * lx)
        if g <= 0:
            continue
        lg = c * lx + k * np.log(mu / rate) + np.log(g)
        if order == 1:
            lg += np.log(k / rate)
        y += f.w[j] * np.exp(lg)
    return float(y)


def _simpson(y, dx):
    """Composite Simpson over an even number of subdivisions."""
    n = y.size
    return dx / 3 * (y[0] + y[n - 1] + 4 * np.sum(y[1:n - 1:2]) + 2 * np.sum(y[2:n - 2:2]))


def _cumsimpson(y, dx):
    """Cumulative Simpson: full panels at the odd nodes, half panel at the even ones.

    The primitive must be available at every grid node, the profile being needed again at the next
    population step, which quadrature at arbitrary abscissae could not provide.
    """
    n = y.size
    out = np.zeros(n)
    for i in range(2, n, 2):
        out[i] = out[i - 2] + dx / 3 * (y[i - 2] + 4 * y[i - 1] + y[i])
    for i in range(1, n, 2):
        if i + 1 < n:
            out[i] = out[i - 1] + dx / 12 * (5 * y[i - 1] + 8 * y[i] - y[i + 1])
        else:
            out[i] = out[i - 1] + dx / 12 * (-y[i - 2] + 8 * y[i - 1] + 5 * y[i])
    return out


class _Grid:
    """Job size grid of one SJN station and the population-independent integrals over it."""

    __slots__ = ('lx', 'dx', 'x', 'f', 'theta', 'tail0', 'tail1', 'fit')


def _setup(S, scv, ns, lfactor):
    """Build the grid of one station.

    It spans ``[0, lfactor * max_r s_r]`` because the conditional waiting time has flattened out
    well before that point, its remainder being carried by the analytic tail of the recursion.
    """
    R = S.size
    smax = float(np.max(S))
    if smax <= 0:
        raise ValueError('the station has zero service demand in every class')
    g = _Grid()
    g.lx = lfactor * smax
    g.x = np.linspace(0.0, g.lx, ns + 1)
    g.dx = g.x[1] - g.x[0]
    g.f = np.zeros((ns + 1, R))
    g.theta = np.zeros((ns + 1, R))
    g.tail0 = np.zeros(R)
    g.tail1 = np.zeros(R)
    g.fit = [None] * R
    for r in range(R):
        g.fit[r] = _fit(float(S[r]), float(scv[r]))
        if g.fit[r].empty:
            continue
        g.f[:, r] = _pdf(g.fit[r], g.x)
        g.theta[:, r] = _theta(g.fit[r], g.x)
        g.tail0[r] = _ccdf(g.fit[r], g.lx)
        g.tail1[r] = S[r] - g.theta[-1, r]
    return g


def _singular_message(m):
    return ('the SJN recursion at station %d has no solution: the work brought by jobs no longer '
            'than the tagged one saturates the server, at which point long jobs starve and the '
            'arrival theorem no longer holds. Reduce the load at that station or model it with '
            'SolverCTMC or SolverLDES.' % m)


def _station(m, r, g, S, scv, V, st, beta, useprio, prio):
    """One evaluation of the conditional waiting time equation at an SJN station.

    ``lam_k W_k(x) f_k(x)`` is the density, in the job size, of the queued class-k customers, so
    deflating it by ``beta_k`` turns the same equation into either the exact recursion (beta = 1,
    the state already being the one at n - e_r) or the Schweitzer closure (beta_r = (N_r-1)/N_r,
    the state being the one at N).
    """
    R = S.size
    lamb = beta * st['lam']
    ub = beta * st['U']
    qb = beta * st['Q']
    rl = float(np.sum((1 + scv) * S * ub) / 2)
    if useprio:
        hi = prio < prio[r]
        base = rl + float(np.sum(S[hi] * (qb[hi] - ub[hi])))
        uhi = float(np.sum(ub[hi]))
        num = base + lamb[r] * st['phi'][:, r]
        den = 1 - uhi - lamb[r] * g.theta[:, r]
        numinf = base + lamb[r] * st['phiinf'][r]
        deninf = 1 - uhi - lamb[r] * S[r]
    else:
        num = rl + st['phi'] @ lamb
        den = 1 - g.theta @ lamb
        numinf = rl + float(lamb @ st['phiinf'])
        deninf = 1 - float(np.sum(lamb * S))
    if np.any(den <= 0) or deninf <= 0:
        raise ValueError(_singular_message(m))
    W = num / den
    winf = numinf / deninf
    # eq. (14) of the reference, generalised by differentiating the recursion at the grid edge
    if useprio:
        slope = g.lx * lamb[r] * g.f[-1, r] * (st['W'][-1, r] + W[-1]) / den[-1]
    else:
        slope = g.lx * float(np.sum(lamb * g.f[-1, :] * (st['W'][-1, :] + W[-1]))) / den[-1]
    a = winf
    b = winf - W[-1]
    if b <= 0:
        b = 0.0
        c = 0.0
    elif slope < 0:
        raise ValueError('the conditional waiting time at SJN station %d decreases in the job '
                         'size, which the discipline forbids: the recursion has become '
                         'numerically unstable.' % m)
    else:
        c = slope / b
    phi = _cumsimpson(W * g.x * g.f[:, r], g.dx)
    phiinf = phi[-1] + a * g.tail1[r] - b * _tailmom(g.fit[r], g.lx, c, 1)
    wbar = _simpson(W * g.f[:, r], g.dx) + a * g.tail0[r] - b * _tailmom(g.fit[r], g.lx, c, 0)
    return V[r] * (S[r] + wbar), W, phi, phiinf, np.array([a, b, c])


def _thru(C, N, Z):
    """Throughputs implied by the residence times, keeping Little's law exact."""
    den = Z + np.sum(C, axis=0)
    X = np.zeros(N.size)
    act = N > 0
    X[act] = N[act] / den[act]
    return X


def _cap(C, L, N, Z, sjnset, umax):
    """Enforce U <= umax at every SJN station by inflating its waiting time.

    The response time equation is an open-system one and has no solution once the fraction of the
    server taken by jobs no longer than x reaches one. A closed network never reaches it in
    reality, but the approximation can, because it underestimates the residence time at a congested
    SJN station and the throughput then exceeds the station capacity. What is imposed is the
    utilization law ``sum_r X_r L_mr <= umax``, an exact property of the network.

    The constraint acts on the waiting time, i.e. on the excess ``C - L``, and never on the
    throughput directly, so that ``X (Z + sum_m C) = N`` still holds exactly and no jobs are lost.
    The same factor scales the station's conditional waiting time profile.
    """
    nsjn = len(sjnset)
    kappa = np.ones(nsjn)
    bound = False
    X = _thru(C, N, Z)
    if nsjn == 0:
        return C, X, kappa, bound
    if umax >= 1:
        raise ValueError('the utilization cap must be strictly below one, the response time '
                         'equation is singular at one')

    def rho_at(m, wq, kap):
        saved = C[m, :].copy()
        C[m, :] = L[m, :] + kap * wq
        rho = float(np.sum(_thru(C, N, Z) * L[m, :]))
        C[m, :] = saved
        return rho

    for _ in range(20):
        viol = False
        for q in range(nsjn):
            m = sjnset[q]
            if float(np.sum(X * L[m, :])) <= umax:
                continue
            viol = True
            bound = True
            wq = C[m, :] - L[m, :]
            hi = 2.0
            while rho_at(m, wq, hi) > umax:
                hi *= 2
                if hi > 1e12:
                    raise ValueError('station %d cannot be brought under the utilization cap by '
                                     'any waiting time: its service demands alone saturate it at '
                                     'this population.' % m)
            lo = 1.0
            for _b in range(200):
                mid = (lo + hi) / 2
                if rho_at(m, wq, mid) > umax:
                    lo = mid
                else:
                    hi = mid
            kappa[q] *= hi
            C[m, :] = L[m, :] + hi * wq
            X = _thru(C, N, Z)
        if not viol:
            return C, X, kappa, bound
    raise ValueError('the utilization cap did not settle across the SJN stations')


def _args(L, N, Z, scv, sjnset, V, options):
    L = np.atleast_2d(np.asarray(L, dtype=float))
    if L.shape[0] == 1 and L.shape[1] > 1 and np.ndim(N) == 0:
        L = L.T
    M, R = L.shape
    N = np.atleast_1d(np.asarray(N, dtype=float)).ravel()
    if N.size != R:
        raise ValueError('demand matrix and population vector have different number of classes')
    N = np.round(N)
    if np.any(N < 0):
        raise ValueError('negative class populations')
    Z = np.zeros(R) if Z is None else np.atleast_1d(np.asarray(Z, dtype=float)).ravel()
    if Z.size != R:
        raise ValueError('the think times and the demand matrix disagree on the class count')
    scv = np.ones((M, R)) if scv is None else np.asarray(scv, dtype=float).reshape(M, R)
    V = np.ones((M, R)) if V is None else np.asarray(V, dtype=float).reshape(M, R)
    S = np.zeros((M, R))
    nz = V > 0
    S[nz] = L[nz] / V[nz]
    sjn = np.zeros(0, dtype=int) if sjnset is None else np.asarray(sjnset, dtype=int).ravel()
    if sjn.size and (np.any(sjn < 0) or np.any(sjn >= M)):
        raise ValueError('sjnset contains a station index outside the demand matrix')
    if np.unique(sjn).size != sjn.size:
        raise ValueError('sjnset repeats a station index')
    options = SjnOptions() if options is None else options
    options.validate(R)
    return M, R, L, N, Z, scv, sjn, V, S, options


class SjnStarvationError(RuntimeError):
    """The SJN waiting time equation has no solution at this population."""


def _warn_capped(umax):
    from ..io.logging import line_warning
    line_warning('pfqn_sjn', 'the utilization cap of %g was binding at an SJN station: the station '
                 'is in the starvation regime, where long jobs are held back and the arrival '
                 'theorem is badly violated. The results are stable but their accuracy is not '
                 'warranted, use SolverCTMC or SolverLDES there.' % umax)


def pfqn_mvasjn(L, N, Z=None, scv=None, sjnset=None, V=None, options=None):
    """Mean value analysis with shortest-job-next stations, over the population lattice.

    The recursion is explicit: ``W(.,n)`` needs only ``phi(.,n-1)``, so it is carried alongside the
    population recursion of exact MVA. This costs ``prod(N+1)`` steps; :func:`pfqn_amvasjn` is the
    fixed-point counterpart.

    The service time density is not an input: only its mean and squared coefficient of variation
    are, and the density is reconstructed by the two-moment Erlang-mixture fit the reference
    prescribes. The x-integrals run on a fixed grid by composite Simpson, ``W(.,n)`` being needed
    at the next population so that quadrature rules sampling at arbitrary abscissae cannot be used;
    beyond the grid the profile is closed by the analytic tail ``W = a - b exp(-c (x - Lx))``.

    :param L: service demand matrix (M x R) of the queueing stations
    :param N: population vector (1 x R)
    :param Z: think time vector (1 x R)
    :param scv: squared coefficients of variation of the service times (M x R)
    :param sjnset: zero-based indices of the stations scheduling by SJN
    :param V: visit ratios (M x R), so that the per-visit service time is L/V
    :param options: :class:`SjnOptions`
    :return: (X, Q, U, C, profiles, iter)
    """
    M, R, L, N, Z, scv, sjn, V, S, options = _args(L, N, Z, scv, sjnset, V, options)
    useprio = options.prio is not None
    ns = options.ns
    ngrid = ns + 1
    nsjn = sjn.size
    G = [_setup(S[sjn[q], :], scv[sjn[q], :], ns, options.lfactor) for q in range(nsjn)]

    stride = np.ones(R, dtype=int)
    npop = 1
    for r in range(R):
        stride[r] = npop
        npop *= int(N[r]) + 1
    Xp = np.zeros((npop, R))
    Qp = np.zeros((npop, M, R))
    Up = np.zeros((npop, M, R))
    Cp = np.zeros((npop, M, R))
    Wp = [np.zeros((npop, ngrid, R)) for _ in range(nsjn)]
    Pp = [np.zeros((npop, ngrid, R)) for _ in range(nsjn)]
    Ip = [np.zeros((npop, R)) for _ in range(nsjn)]
    Tp = [np.zeros((npop, R, 3)) for _ in range(nsjn)]

    beta1 = np.ones(R)
    sjnlist = list(sjn)
    for idx in range(1, npop):
        n = np.array([(idx // stride[r]) % (int(N[r]) + 1) for r in range(R)], dtype=float)
        Call = np.zeros((M, R))
        for r in range(R):
            if n[r] == 0:
                continue
            iprev = idx - stride[r]
            for m in range(M):
                if m not in sjnlist:
                    Call[m, r] = L[m, r] * (1 + float(np.sum(Qp[iprev, m, :])))
                    continue
                q = sjnlist.index(m)
                # the population step already supplies the neighbouring profile, no deflation
                st = {'lam': Xp[iprev, :] * V[m, :], 'U': Up[iprev, m, :], 'Q': Qp[iprev, m, :],
                      'W': Wp[q][iprev], 'phi': Pp[q][iprev], 'phiinf': Ip[q][iprev]}
                cmr, wprof, phiprof, phiinf, tail = _station(m, r, G[q], S[m, :], scv[m, :],
                                                             V[m, :], st, beta1, useprio,
                                                             options.prio)
                Call[m, r] = cmr
                Wp[q][idx][:, r] = wprof
                Pp[q][idx][:, r] = phiprof
                Ip[q][idx][r] = phiinf
                Tp[q][idx][r, :] = tail
        Call, Xn, kappa, bound = _cap(Call, L, n, Z, sjn, options.umax)
        if bound:
            # the cap has invalidated the profile the next population step reads back
            raise SjnStarvationError(
                'the utilization cap of %g was binding at an SJN station at population %s: the '
                'station is in the starvation regime, where the conditional waiting time equation '
                'has no solution and the population lattice no valid continuation. Use the '
                'Schweitzer fixed point (pfqn_amvasjn, method \'amva\'), SolverCTMC or SolverLDES.'
                % (options.umax, np.array2string(n)))
        Xp[idx, :] = Xn
        Cp[idx] = Call
        Qp[idx] = Xn[None, :] * Call
        Up[idx] = Xn[None, :] * L

    last = npop - 1
    profiles = [{'station': int(sjn[q]), 'x': G[q].x, 'W': Wp[q][last], 'tail': Tp[q][last]}
                for q in range(nsjn)]
    return Xp[last, :], Qp[last], Up[last], Cp[last], profiles, 1


def pfqn_amvasjn(L, N, Z=None, scv=None, sjnset=None, V=None, options=None):
    """Mean value analysis with shortest-job-next stations, through a Schweitzer fixed point.

    :func:`pfqn_mvasjn` carries the conditional waiting time profile over the whole population
    lattice, which costs ``prod(N+1)`` steps. The closure used here rests on the observation that
    ``lam_k W_k(x,n) f_k(x) dx`` is the mean number of queued class-k customers whose service
    requirement lies in ``(x, x+dx)``, that is, the queue length resolved by job size. Schweitzer's
    assumption is applied to that density rather than to its integral: removing one customer of
    class r scales the class-r size-resolved queue length by ``(N_r-1)/N_r``. Integrating over x
    recovers the usual rule for the aggregate queue lengths, so the closure is the exact analogue
    of the one applied at the ordinary stations.

    What is given up is the population dependence of the SHAPE of ``W(x)``: the closure lets its
    level scale but keeps its shape fixed, whereas the true profile stiffens with the load because
    the denominator sharpens. The error therefore concentrates at high utilization, where the SJN
    approximation is already at its weakest.

    The iteration is started from the product-form Schweitzer solution, not from a light-load
    guess: the latter puts the deflated utilization above one, where the equation has no solution.

    :return: (X, Q, U, C, profiles, iter)
    """
    from .mva import pfqn_bs
    M, R, L, N, Z, scv, sjn, V, S, options = _args(L, N, Z, scv, sjnset, V, options)
    useprio = options.prio is not None
    ns = options.ns
    ngrid = ns + 1
    nsjn = sjn.size
    G = [_setup(S[sjn[q], :], scv[sjn[q], :], ns, options.lfactor) for q in range(nsjn)]

    # start from the product-form Schweitzer solution: a light-load guess would put the deflated
    # utilization above one and the SJN denominator has no solution there
    Xb, Qb, Ub, Cb = pfqn_bs(L, N, Z, options.tol, options.iter_max)[:4]
    X = np.asarray(Xb, dtype=float).ravel().copy()
    Q = np.asarray(Qb, dtype=float).reshape(M, R).copy()
    U = np.asarray(Ub, dtype=float).reshape(M, R).copy()
    C = np.asarray(Cb, dtype=float).reshape(M, R).copy()
    W = [np.zeros((ngrid, R)) for _ in range(nsjn)]
    P = [np.zeros((ngrid, R)) for _ in range(nsjn)]
    Iinf = [np.zeros(R) for _ in range(nsjn)]
    T = [np.zeros((R, 3)) for _ in range(nsjn)]

    sjnlist = list(sjn)
    it = 0
    capped = False
    converged = False
    delta = np.inf
    while not converged and it < options.iter_max:
        it += 1
        Cit = np.zeros((M, R))
        Wit = [w.copy() for w in W]
        Pit = [p.copy() for p in P]
        Iit = [i.copy() for i in Iinf]
        Tit = [t.copy() for t in T]
        for r in range(R):
            if N[r] == 0:
                continue
            beta = np.ones(R)
            beta[r] = (N[r] - 1) / N[r]
            for m in range(M):
                if m not in sjnlist:
                    Cit[m, r] = L[m, r] * (1 + float(np.sum(beta * Q[m, :])))
                    continue
                q = sjnlist.index(m)
                st = {'lam': X * V[m, :], 'U': U[m, :], 'Q': Q[m, :],
                      'W': W[q], 'phi': P[q], 'phiinf': Iinf[q]}
                cmr, wprof, phiprof, phiinf, tail = _station(m, r, G[q], S[m, :], scv[m, :],
                                                             V[m, :], st, beta, useprio,
                                                             options.prio)
                Cit[m, r] = cmr
                Wit[q][:, r] = wprof
                Pit[q][:, r] = phiprof
                Iit[q][r] = phiinf
                Tit[q][r, :] = tail
        Cit, Xit, kappa, bound = _cap(Cit, L, N, Z, sjn, options.umax)
        if bound:
            capped = True
            for q in range(nsjn):
                Wit[q] *= kappa[q]
                Pit[q] *= kappa[q]
                Iit[q] *= kappa[q]
                Tit[q][:, 0:2] *= kappa[q]
        Qit = Xit[None, :] * Cit
        delta = float(np.max(np.abs(Qit - Q)))
        for q in range(nsjn):
            delta = max(delta, float(np.max(np.abs(Wit[q] - W[q]))))
        X = Xit
        C = Cit
        Q = Qit
        U = Xit[None, :] * L
        W, P, Iinf, T = Wit, Pit, Iit, Tit
        converged = delta < options.tol
    if not converged:
        from ..io.logging import line_warning
        line_warning('pfqn_amvasjn', 'the SJN fixed point did not converge in %d iterations, '
                     'residual %g' % (options.iter_max, delta))
    if capped:
        _warn_capped(options.umax)

    profiles = [{'station': int(sjn[q]), 'x': G[q].x, 'W': W[q], 'tail': T[q]}
                for q in range(nsjn)]
    return X, Q, U, C, profiles, it
