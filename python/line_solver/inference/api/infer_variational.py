"""Variational inference for Markovian queueing networks.

Port of the MATLAB `infer_variational`, following I. Perez, G. Casale,
"Variational Inference for Markovian Queueing Networks", Advances in Applied
Probability 53(3), 2021.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
This code is released under the 3-Clause BSD License.
"""

import numpy as np
from scipy.special import digamma, gammaln


class VariationalSpec(object):
    """Inference problem handed to :func:`infer_variational`.

    Station-class pairs are flattened column-major, so that pair (m, r) sits
    at index r*M + m.

    Attributes:
        arcs: (T, 3) integer array of transitions [i, j, c]; i == 0 marks an
            external source and j == 0 a sink. Station and class indices are
            one-based, matching the MATLAB and Java specifications.
        x0: (M, R) initial queue lengths.
        sched: (M,) discipline codes, 0 = infinite server, 1 = shared server
            (PS/FCFS), 2 = external.
        nservers: (M,) number of servers.
        routeprob: (T,) routing probability of each transition.
        arcparam: (T,) index in 1..P of the rate governing the transition,
            0 when the rate is known.
        arcrate: (T,) known rate for transitions with arcparam == 0.
        alpha0, beta0: (P,) Gamma prior shape and rate.
        obsTimes: (K,) observation epochs.
        obsData: (K, M*R) observed queue lengths, NaN where not observed.
        obsRange: (M, R) support size of the uniform contamination.
        epsilon: probability that a reading is faulty.
        capacity: (M, R) upper bound on the queue length, infinite by default.
            In a closed network this is the chain population, and clamping the
            load there keeps the expanded state space from crediting a station
            with more jobs than the network holds.
    """

    def __init__(self, arcs, x0, sched, nservers, routeprob, arcparam, arcrate,
                 alpha0, beta0, obsTimes, obsData, obsRange, epsilon, capacity=None):
        self.arcs = np.asarray(arcs, dtype=int).reshape(-1, 3)
        self.x0 = np.asarray(x0, dtype=float)
        if self.x0.ndim == 1:
            self.x0 = self.x0.reshape(-1, 1)
        self.sched = np.asarray(sched, dtype=int).flatten()
        self.nservers = np.asarray(nservers, dtype=float).flatten()
        self.routeprob = np.asarray(routeprob, dtype=float).flatten()
        self.arcparam = np.asarray(arcparam, dtype=int).flatten()
        self.arcrate = np.asarray(arcrate, dtype=float).flatten()
        self.alpha0 = np.atleast_1d(np.asarray(alpha0, dtype=float)).flatten()
        self.beta0 = np.atleast_1d(np.asarray(beta0, dtype=float)).flatten()
        self.obsTimes = np.asarray(obsTimes, dtype=float).flatten()
        self.obsData = np.asarray(obsData, dtype=float).reshape(len(self.obsTimes), -1)
        self.obsRange = np.asarray(obsRange, dtype=float).reshape(-1, order='F')
        self.epsilon = float(epsilon)
        if capacity is None:
            self.capacity = np.full(self.x0.size, np.inf)
        else:
            self.capacity = np.asarray(capacity, dtype=float).reshape(-1, order='F')


class VariationalOptions(object):
    """Options of :func:`infer_variational`."""

    def __init__(self, **kwargs):
        self.verbose = 0
        self.iter_max = 20
        self.tol = 1e-3
        self.nsamples = 200
        self.delta = 1e-3
        self.floor = 1e-4
        self.rate_max = None
        self.rate_cap_factor = 10.0
        self.unifmax = 30.0
        self.unif_tol = 1e-12
        self.unif_max_terms = 2000
        self.tmax = None
        self.dt = None
        self.ngrid = None
        self.ymax = None
        for k, v in kwargs.items():
            setattr(self, k, v)


class VariationalResult(object):
    """Outcome of :func:`infer_variational`."""

    def __init__(self):
        self.alpha = None
        self.beta = None
        self.rates = None
        self.meanServiceTime = None
        self.bound = None
        self.alphaTrace = None
        self.betaTrace = None
        self.Y = None
        self.nu = None
        self.tgrid = None
        self.qlen = None
        self.iter = 0
        self.converged = False
        self.tailmass = 0.0


def _ups(xic, xis, nservers, sched, cap=np.inf, capstat=np.inf):
    """Load factor Upsilon of a transition leaving a station-class pair."""
    if sched == 2:
        return np.ones_like(np.asarray(xic, dtype=float))
    xic = np.minimum(cap, np.maximum(0.0, xic))
    if sched == 0:
        return xic
    xis = np.minimum(capstat, np.maximum(0.0, xis))
    u = np.zeros(np.broadcast(xic, xis).shape)
    pos = np.broadcast_to(xis, u.shape) > 0
    xicb = np.broadcast_to(xic, u.shape)
    xisb = np.broadcast_to(xis, u.shape)
    u[pos] = xicb[pos] / xisb[pos] * np.minimum(nservers, xisb[pos])
    return u


def _prime(k):
    """k-th prime, k >= 1."""
    n = 0
    c = 1
    p = 2
    while n < k:
        c += 1
        isp = True
        d = 2
        while d * d <= c:
            if c % d == 0:
                isp = False
                break
            d += 1
        if isp:
            n += 1
            p = c
    return p


def _radical_inverse(i, base):
    """Van der Corput radical inverse of i in the given base."""
    r = 0.0
    f = 1.0 / base
    while i > 0:
        r += f * (i % base)
        i //= base
        f /= base
    return r


def _lattice(S, e):
    """Halton points of transition e, sorted ascending with their order."""
    base = _prime(e + 1)
    u = np.array([_radical_inverse(s + 1, base) for s in range(S)])
    perm = np.argsort(u, kind='stable')
    return u[perm], perm


def _sample(q, S, e):
    """Inverse-c.d.f. samples of a marginal on a Halton lattice."""
    G, ny = q.shape
    us, perm = _lattice(S, e)
    ys = np.zeros((G, S))
    for g in range(G):
        c = np.cumsum(q[g, :])
        if c[ny-1] > 0:
            c = c / c[ny-1]
        c[ny-1] = 1.0
        idx = np.searchsorted(c, us, side='left')
        idx = np.minimum(idx, ny - 1)
        ys[g, perm] = idx
    return ys


def _sample_all(Y, S):
    narcs, G, ny = Y.shape
    ys = np.zeros((narcs, G, S))
    for e in range(narcs):
        ys[e] = _sample(Y[e], S, e)
    return ys


def _obs_weight(obsRow, obsRange, epsilon, Aall, sgnE, yvec, opt):
    """Multiplicative jump carried by an observation in the backward pass."""
    ny = yvec.size
    S = Aall.shape[1]
    acc = np.zeros((ny, S))
    for k in np.flatnonzero(~np.isnan(obsRow)):
        x = Aall[k, :][None, :] + sgnE[k] * yvec[:, None]
        hit = (x == obsRow[k])
        feas = (x >= 0) & (x <= obsRange[k])
        p = hit * (1.0 - epsilon) + (~hit & feas) * (epsilon / max(1.0, obsRange[k]))
        acc += np.log(opt.floor + p)
    return np.exp(np.mean(acc, axis=1))


def _damp(sl, Ye, opt):
    """Slack multiplier of the rate cap; unity when the cap is inactive."""
    z = sl / np.maximum(opt.floor, Ye)
    damp = (1.0 + z) / np.exp(z)
    damp[sl == 0] = 1.0
    damp[~np.isfinite(damp) | (damp < 0)] = 0.0
    return damp


def _rescale(v):
    m = np.max(v)
    if m > 0 and np.isfinite(m):
        return v / m
    return v


def _back_uniformize(v0, pd, pu, lt, opt):
    """One uniformization step of the backward sub-generator."""
    ny = v0.size
    w = np.exp(-lt)
    v = w * v0
    u = v0.copy()
    cum = w
    n = 1
    while (1.0 - cum) > opt.unif_tol and n < opt.unif_max_terms:
        un = u * (1.0 - pd)
        un[:ny-1] += u[1:] * pu[:ny-1]
        u = un
        w = w * lt / n
        v = v + w * u
        cum += w
        n += 1
    return v


def _backward(ge, he, sl, Ye, obsIdx, obsw, dt, opt):
    """Backward pass for the Lagrange multipliers, Eq. (14).

    The equation is linear in r, so on a grid cell with frozen coefficients
    it is the action of a matrix exponential. The generator has non-positive
    row sums by Jensen, so uniformization evaluates it without the stiffness
    that an explicit rule suffers when exp(E log Xi) falls orders of
    magnitude below E[Xi].
    """
    G, ny = ge.shape
    r = np.zeros((G, ny))
    v = np.ones(ny)
    for q in np.flatnonzero(obsIdx == G - 1):
        v = v * np.maximum(0.0, obsw[q, :])
    v = _rescale(v)
    r[G-1, :] = v
    for g in range(G-2, -1, -1):
        gv = np.maximum(0.0, ge[g, :])
        hv = np.maximum(0.0, he[g, :] * _damp(sl[g, :], Ye[g, :], opt))
        hv = np.minimum(hv, gv)
        Lam = np.max(gv)
        if Lam > 0:
            ncell = max(1, int(np.ceil(Lam * dt / opt.unifmax)))
            h = dt / ncell
            pd = gv / Lam
            pu = hv / Lam
            for _ in range(ncell):
                v = _back_uniformize(v, pd, pu, Lam * h, opt)
            v = _rescale(v)
        for q in np.flatnonzero(obsIdx == g):
            v = v * np.maximum(0.0, obsw[q, :])
            v = _rescale(v)
        r[g, :] = v
    return r


def _uniformize(v0, p, lt, opt):
    """One uniformization step of the pure-birth chain."""
    ny = v0.size
    w = np.exp(-lt)
    v = w * v0
    u = v0.copy()
    cum = w
    n = 1
    while (1.0 - cum) > opt.unif_tol and n < opt.unif_max_terms:
        un = u * (1.0 - p)
        un[1:] += u[:ny-1] * p[:ny-1]
        u = un
        w = w * lt / n
        v = v + w * u
        cum += w
        n += 1
    return v


def _forward(nue, dt, opt):
    """Forward master equation of an inhomogeneous pure-birth process."""
    G, ny = nue.shape
    q = np.zeros((G, ny))
    v = np.zeros(ny)
    v[0] = 1.0
    q[0, :] = v
    for g in range(G-1):
        rates = np.maximum(0.0, nue[g, :])
        Lam = np.max(rates)
        if Lam <= 0:
            q[g+1, :] = v
            continue
        ncell = max(1, int(np.ceil(Lam * dt / opt.unifmax)))
        h = dt / ncell
        p = rates / Lam
        for _ in range(ncell):
            v = _uniformize(v, p, Lam * h, opt)
        v = np.maximum(0.0, v)
        s = np.sum(v)
        if s > 0:
            v = v / s
        q[g+1, :] = v
    return q


def _trapz(f, dt):
    f = np.asarray(f, dtype=float).flatten()
    if f.size < 2:
        return 0.0
    return dt * (np.sum(f) - 0.5 * f[0] - 0.5 * f[-1])


def _kl_gamma(a, b, a0, b0):
    """KL(Gamma(a,b) || Gamma(a0,b0)) with rate parameterisation."""
    return ((a - a0) * digamma(a) - gammaln(a) + gammaln(a0)
            + a0 * (np.log(b) - np.log(b0)) + a * (b0 - b) / b)


def _setup(spec, options):
    """Validate the specification and fill in the default options."""
    M, R = spec.x0.shape
    narcs = spec.arcs.shape[0]
    if np.any((spec.arcs[:, 0] == 0) & (spec.arcs[:, 1] == 0)):
        raise ValueError('A transition cannot be external at both ends.')
    if spec.sched.size != M or spec.nservers.size != M:
        raise ValueError('sched and nservers must have one entry per station.')
    if spec.routeprob.size != narcs or spec.arcparam.size != narcs or spec.arcrate.size != narcs:
        raise ValueError('routeprob, arcparam and arcrate must have one entry per transition.')
    if spec.alpha0.size != spec.beta0.size:
        raise ValueError('alpha0 and beta0 must have the same length.')
    for e in range(narcs):
        if spec.arcparam[e] == 0 and not spec.arcrate[e] > 0:
            raise ValueError('Transition %d has no parameter and no positive known rate.' % e)
    if spec.obsData.shape[0] != spec.obsTimes.size:
        raise ValueError('obsData must have one row per observation epoch.')
    if spec.obsData.shape[1] != M * R:
        raise ValueError('obsData must have M*R columns.')
    if spec.capacity.size != M * R:
        raise ValueError('capacity must have M*R entries.')

    opt = options if options is not None else VariationalOptions()

    if opt.tmax is None:
        if spec.obsTimes.size == 0:
            raise ValueError('options.tmax is required when there are no observations.')
        opt.tmax = float(np.max(spec.obsTimes))
    if opt.tmax <= 0:
        raise ValueError('options.tmax must be positive.')
    if opt.ngrid is None and opt.dt is None:
        opt.ngrid = 201
    if opt.ngrid is None:
        opt.ngrid = int(round(opt.tmax / opt.dt)) + 1
    opt.ngrid = max(2, int(round(opt.ngrid)))
    opt.dt = opt.tmax / (opt.ngrid - 1)

    xbar = spec.x0.reshape(-1, order='F').copy()
    for k in range(M * R):
        col = spec.obsData[:, k]
        col = col[~np.isnan(col)]
        if col.size > 0:
            xbar[k] = np.mean(col)
    xbars = np.zeros(M)
    for m in range(M):
        xbars[m] = np.sum(xbar[m + np.arange(R) * M])
    opt.xbar = xbar
    opt.xbars = xbars

    if opt.ymax is None:
        fmax = 0.0
        for e in range(narcs):
            if spec.arcparam[e] > 0:
                p = spec.arcparam[e] - 1
                lam = spec.routeprob[e] * spec.alpha0[p] / spec.beta0[p]
            else:
                lam = spec.routeprob[e] * spec.arcrate[e]
            i = spec.arcs[e, 0]
            if i > 0:
                u = float(_ups(xbar[(spec.arcs[e, 2]-1)*M + i - 1], xbars[i-1],
                               spec.nservers[i-1], spec.sched[i-1],
                               spec.capacity[(spec.arcs[e, 2]-1)*M + i - 1]))
            else:
                u = 1.0
            fmax = max(fmax, lam * u * opt.tmax)
        opt.ymax = max(20, int(np.ceil(2*fmax + 5*np.sqrt(max(1.0, fmax)))))
    opt.ymax = max(2, int(round(opt.ymax)))

    if opt.rate_max is None:
        opt.rate_max = opt.rate_cap_factor * opt.ymax / opt.tmax
    return opt


def _rate_moments(spec, opt, e, Ys, Y, sgnClass, sgnStat, x0v, x0s, capStat,
                  arcSched, arcServers, alpha, beta, obsIdx, yvec, want_obs):
    """Conditional rate moments of one transition, and its observation jumps.

    Returns E[Xi | Y^eta = y] and exp(E[log Xi | Y^eta = y]) on the time grid,
    both taken under Q with the transition's own contribution removed, plus the
    multiplicative jump each observation carries in the backward pass.
    """
    narcs, G, ny = Y.shape
    S = Ys.shape[2]
    M, R = spec.x0.shape
    src = spec.arcs[:, 0]
    cls = spec.arcs[:, 2]
    p = spec.arcparam[e]
    if p == 0:
        lam = spec.routeprob[e] * spec.arcrate[e]
        loglam = np.log(lam)
    else:
        lam = spec.routeprob[e] * alpha[p-1] / beta[p-1]
        loglam = np.log(spec.routeprob[e]) + digamma(alpha[p-1]) - np.log(beta[p-1])
    kclass = (cls[e]-1)*M + src[e] - 1 if src[e] > 0 else -1
    ge = np.zeros((G, ny))
    he = np.zeros((G, ny))
    obsw = np.ones((spec.obsTimes.size, ny))
    sgnC = sgnClass.copy()
    sgnC[e, :] = 0.0
    sgnS = sgnStat.copy()
    sgnS[e, :] = 0.0
    for g in range(G):
        Yg = Ys[:, g, :]
        Aall = x0v[:, None] + sgnC.T.dot(Yg)
        if kclass >= 0:
            Asta = x0s[:, None] + sgnS.T.dot(Yg)
            xic = Aall[kclass, :][None, :] + sgnClass[e, kclass] * yvec[:, None]
            xis = Asta[src[e]-1, :][None, :] + sgnStat[e, src[e]-1] * yvec[:, None]
            ups = _ups(xic, xis, arcServers[e], arcSched[e],
                       spec.capacity[kclass], capStat[src[e]-1])
        else:
            ups = np.ones((ny, S))
        # E[Xi] and exp(E[log Xi]) of the SAME rate Xi = delta + lam*Ups.
        # Writing the second as exp(E[log lam]) exp(E[log(Ups + delta/E[lam])])
        # keeps the two consistent wherever Ups is deterministic, which is what
        # stops the backward equation from developing a gradient away from the
        # empty-station boundary.
        ge[g, :] = opt.delta + lam * np.mean(ups, axis=1)
        he[g, :] = np.exp(loglam + np.mean(np.log(ups + opt.delta/lam), axis=1))
        if want_obs:
            for q in np.flatnonzero(obsIdx == g):
                obsw[q, :] = _obs_weight(spec.obsData[q, :], spec.obsRange,
                                         spec.epsilon, Aall, sgnClass[e, :], yvec, opt)
    return ge, he, obsw


def infer_variational(spec, options=None):
    """Variational inference for Markovian queueing networks.

    The network trajectory is reparameterised by the transition counts
    Y^eta, eta=(i,j,c), so that the station marginals decouple. The
    variational family is a product of inhomogeneous pure-birth processes,
    one per transition, times a product of Gamma densities over the unknown
    service rates. The state space is expanded by adding DELTA to every
    feasible rate, so that queue lengths may go negative and the
    approximating measure stays mutually absolutely continuous with the
    target; the original model is recovered as DELTA -> 0.

    Expectations over the other transitions are taken on a deterministic
    Halton lattice mapped through the inverse marginal c.d.f., so the
    estimator carries no random-number stream and is reproducible across the
    MATLAB, Java, Python and C++ implementations.

    Args:
        spec: a :class:`VariationalSpec`.
        options: a :class:`VariationalOptions`, or None for the defaults.

    Returns:
        A :class:`VariationalResult`.
    """
    opt = _setup(spec, options)

    M, R = spec.x0.shape
    narcs = spec.arcs.shape[0]
    P = spec.alpha0.size
    src = spec.arcs[:, 0]
    dst = spec.arcs[:, 1]
    cls = spec.arcs[:, 2]

    sgnClass = np.zeros((narcs, M*R))
    sgnStat = np.zeros((narcs, M))
    for e in range(narcs):
        if dst[e] > 0:
            k = (cls[e]-1)*M + dst[e] - 1
            sgnClass[e, k] += 1.0
            sgnStat[e, dst[e]-1] += 1.0
        if src[e] > 0:
            k = (cls[e]-1)*M + src[e] - 1
            sgnClass[e, k] -= 1.0
            sgnStat[e, src[e]-1] -= 1.0

    G = opt.ngrid
    dt = opt.dt
    ymax = opt.ymax
    S = opt.nsamples
    ny = ymax + 1
    yvec = np.arange(ny, dtype=float)
    tgrid = np.arange(G) * dt

    x0v = spec.x0.reshape(-1, order='F')
    x0s = np.sum(spec.x0, axis=1)

    K = spec.obsTimes.size
    obsIdx = np.minimum(G-1, np.maximum(0, np.round(spec.obsTimes/dt).astype(int)))

    capStat = np.zeros(M)
    for m in range(M):
        capStat[m] = np.sum(spec.capacity[m + np.arange(R)*M])

    arcSched = np.zeros(narcs, dtype=int)
    arcServers = np.ones(narcs)
    for e in range(narcs):
        if src[e] > 0:
            arcSched[e] = spec.sched[src[e]-1]
            arcServers[e] = spec.nservers[src[e]-1]
        else:
            arcSched[e] = 2

    alpha = spec.alpha0.copy()
    beta = spec.beta0.copy()

    def rate_mean(e):
        p = spec.arcparam[e]
        if p == 0:
            return spec.routeprob[e] * spec.arcrate[e]
        return spec.routeprob[e] * alpha[p-1] / beta[p-1]

    def rate_log_mean(e):
        p = spec.arcparam[e]
        if p == 0:
            return np.log(spec.routeprob[e] * spec.arcrate[e])
        return np.log(spec.routeprob[e]) + digamma(alpha[p-1]) - np.log(beta[p-1])

    Y = np.zeros((narcs, G, ny))
    nu = np.zeros((narcs, G, ny))
    slack = np.zeros((narcs, G, ny))
    gexp = np.zeros((narcs, G, ny))
    hexp = np.ones((narcs, G, ny))

    for e in range(narcs):
        lam = rate_mean(e)
        if src[e] > 0:
            u0 = float(_ups(opt.xbar[(cls[e]-1)*M + src[e]-1], opt.xbars[src[e]-1],
                            arcServers[e], arcSched[e]))
        else:
            u0 = 1.0
        nue = np.zeros((G, ny))
        nue[:, :ymax] = max(opt.delta, lam * u0)
        nu[e] = nue
        Y[e] = _forward(nue, dt, opt)

    bound = np.zeros(opt.iter_max)
    alphaTrace = np.zeros((P, opt.iter_max))
    betaTrace = np.zeros((P, opt.iter_max))
    converged = False
    it = 0

    for it in range(1, opt.iter_max+1):
        for e in range(narcs):
            Ys = _sample_all(Y, S)
            ge, he, obsw = _rate_moments(spec, opt, e, Ys, Y, sgnClass, sgnStat,
                                         x0v, x0s, capStat, arcSched, arcServers,
                                         alpha, beta, obsIdx, yvec, True)
            Ye = Y[e]
            r = _backward(ge, he, slack[e], Ye, obsIdx, obsw, dt, opt)

            # Eq. (15). A vanishing multiplier marks a count the future
            # observations rule out; the rate there is zero, which is what
            # keeps the forward pass from placing mass on it.
            nue = np.zeros((G, ny))
            den = r[:, :ymax]
            num = r[:, 1:ny]
            ratio = np.zeros((G, ymax))
            pos = den > 0
            ratio[pos] = num[pos] / den[pos]
            nue[:, :ymax] = he[:, :ymax] * ratio
            nue[~np.isfinite(nue) | (nue < 0)] = 0.0
            sl = np.zeros((G, ny))
            if np.isfinite(opt.rate_max):
                over = nue > opt.rate_max
                if np.any(over):
                    sl[over] = np.maximum(opt.floor, Ye[over]) * np.log(nue[over]/opt.rate_max)
                    nue[over] = opt.rate_max
            nu[e] = nue
            slack[e] = sl
            Y[e] = _forward(nue, dt, opt)

        # conjugate Gamma updates: the shape gains the expected number of
        # firings, the rate the expected exposure time of the station-class
        # pair that the parameter governs
        Ys = _sample_all(Y, S)
        firings = np.zeros(P)
        exposure = np.zeros(P)
        seen = np.zeros((P, M*R))
        for e in range(narcs):
            p = spec.arcparam[e]
            if p == 0:
                continue
            # expected number of firings of the transition over the horizon;
            # taken from the marginal itself, which is exact, rather than by
            # quadrature of the intensity, which a near-deterministic
            # marginal makes inaccurate
            firings[p-1] += float(Y[e][G-1, :].dot(yvec) - Y[e][0, :].dot(yvec))
            kclass = (cls[e]-1)*M + src[e] - 1
            if seen[p-1, kclass] == 0:
                seen[p-1, kclass] = 1
                ue = np.zeros(G)
                for g in range(G):
                    Yg = Ys[:, g, :]
                    Aall = x0v[:, None] + sgnClass.T.dot(Yg)
                    Asta = x0s[:, None] + sgnStat.T.dot(Yg)
                    ue[g] = np.mean(_ups(Aall[kclass, :], Asta[src[e]-1, :],
                                         spec.nservers[src[e]-1], spec.sched[src[e]-1],
                                         spec.capacity[kclass], capStat[src[e]-1]))
                exposure[p-1] += _trapz(ue, dt)
        alpha = spec.alpha0 + firings
        beta = spec.beta0 + exposure
        # the bound is evaluated at the state the iteration ended in, so the
        # rate moments are recomputed against the updated marginals rather
        # than reused from the sweep that produced them
        for e in range(narcs):
            gexp[e], hexp[e], _ = _rate_moments(spec, opt, e, Ys, Y, sgnClass, sgnStat,
                                                x0v, x0s, capStat, arcSched, arcServers,
                                                alpha, beta, obsIdx, yvec, False)
        alphaTrace[:, it-1] = alpha
        betaTrace[:, it-1] = beta

        bound[it-1] = _bound(spec, opt, alpha, beta, Y, nu, gexp, hexp, Ys,
                             sgnClass, x0v, obsIdx, dt)
        if opt.verbose > 0:
            print('infer_variational: iteration %d, lower bound %.6f, max rate %.3f'
                  % (it, bound[it-1], np.max(nu)))
        # The rate update solves a stationarity condition rather than
        # maximising the bound in a block, so the bound need not ascend;
        # convergence is judged on the bound AND on the rate posteriors.
        if it > 1:
            crit = abs(bound[it-1] - bound[it-2]) / max(1.0, abs(bound[it-2]))
            if P > 0:
                prev = alphaTrace[:, it-2] / betaTrace[:, it-2]
                crit = max(crit, float(np.max(np.abs(alpha/beta - prev)
                                              / np.maximum(1e-12, prev))))
            if crit <= opt.tol:
                converged = True
                break

    tailmass = float(np.max(Y[:, :, ny-1])) if narcs > 0 else 0.0
    if tailmass > 1e-6:
        import warnings
        warnings.warn('Transition-count truncation ymax=%d carries mass %.3e, '
                      'increase options.ymax.' % (ymax, tailmass))

    qlen = np.tile(x0v, (G, 1))
    for e in range(narcs):
        my = Y[e].dot(yvec)
        qlen = qlen + np.outer(my, sgnClass[e, :])

    out = VariationalResult()
    out.alpha = alpha
    out.beta = beta
    out.rates = alpha / beta
    out.meanServiceTime = beta / alpha
    out.bound = bound[:it]
    out.alphaTrace = alphaTrace[:, :it]
    out.betaTrace = betaTrace[:, :it]
    out.Y = Y
    out.nu = nu
    out.tgrid = tgrid
    out.qlen = qlen
    out.iter = it
    out.converged = converged
    out.tailmass = tailmass
    return out


def _bound(spec, opt, alpha, beta, Y, nu, gexp, hexp, Ys, sgnClass, x0v, obsIdx, dt):
    """Evidence lower bound: path term, observation term and the divergence
    of the rate posteriors from their priors."""
    narcs, G, ny = Y.shape
    S = Ys.shape[2]
    b = 0.0
    for e in range(narcs):
        Ye = Y[e]
        nue = nu[e]
        ge = gexp[e]
        he = hexp[e]
        term = nue - ge
        pos = nue > 0
        tmp = np.zeros((G, ny))
        tmp[pos] = nue[pos] * np.log(nue[pos] / np.maximum(opt.floor, he[pos]))
        term = term - tmp
        b += _trapz(np.sum(Ye * term, axis=1), dt)
    for k in range(spec.obsTimes.size):
        g = obsIdx[k]
        Yg = Ys[:, g, :]
        Xall = x0v[:, None] + sgnClass.T.dot(Yg)
        acc = np.zeros(S)
        row = spec.obsData[k, :]
        for j in np.flatnonzero(~np.isnan(row)):
            x = Xall[j, :]
            hit = (x == row[j])
            feas = (x >= 0) & (x <= spec.obsRange[j])
            p = hit * (1.0 - spec.epsilon) + (~hit & feas) * (spec.epsilon / max(1.0, spec.obsRange[j]))
            acc += np.log(opt.floor + p)
        b += float(np.mean(acc))
    for p in range(alpha.size):
        b -= _kl_gamma(alpha[p], beta[p], spec.alpha0[p], spec.beta0[p])
    return b
