"""
Sojourn time distribution of the M/G/1 processor-sharing queue.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``qsys_mg1_ps.m`` and the JAR twin ``Qsys_mg1_ps.java``.

Jobs arrive in a Poisson stream of rate ``lam`` at a single egalitarian
processor-sharing server whose service requirement has Laplace-Stieltjes
transform ``bhat(tau)`` and mean ``m1``. Writing ``V(x)`` for the sojourn time
of a tagged job of service requirement ``x`` and ``rho = lam*m1 < 1``, Ott
(1984) and Yashkov (1983) express the conditional transform as::

    E[exp(-s V(x))] = (1-rho) / D(s,x)

where ``D(s,x)`` is the inverse Laplace transform, evaluated at ``x``, of::

    f(tau;s) = ( (1-rho)*tau**2 - (1-rho)*lam*(1-bhat(tau))*tau
                 + s*rho*tau - s*lam*(1-bhat(tau)) )
               / ( tau**2 * (tau - s - lam*(1-bhat(tau))) )

The transform is exact but implicit, since ``f`` must be inverted in ``tau``.
For phase-type service ``f(tau;s)`` is a proper rational function of ``tau``,
the double pole at ``tau = 0`` cancels, and ``D(s,x)`` follows in closed form as
a finite sum of residues (or, for repeated poles, from the matrix exponential of
the companion realization). This makes the M/PH/1-PS queue, hence every service
law that can be fitted with a phase-type distribution, exactly solvable. For a
service transform supplied as a callable, ``f`` is inverted in ``tau``
numerically on a Bromwich contour placed to the right of the dominant
singularity ``tau*(s)``, the unique root of ``tau = s + lam*(1-bhat(tau))`` in
the right half plane, which the fixed-point iteration of that equation reaches
at geometric rate ``rho``.

The conditional sojourn time is atomic on the lattice ``t = (k+1)*x``,
``k = 0,1,2,...``: processor sharing gives every job in the system the same
amount of work, so if the ``k`` jobs present on arrival all outlive the tagged
job and no arrival intervenes, the sojourn is exactly ``(k+1)*x``. For
exponential service the masses are
``A_k = (1-rho)*rho**k*exp(-k*mu*x)*exp(-lam*(k+1)*x)``, the ``k = 0`` term
being the probability of finding the system empty and sharing it with nobody,
which is the only one that stays exact for general service. The ``k = 0`` atom
is removed before inverting in ``s``; the remaining atoms make ``cdfCond`` jump
and leave no density, so ``pdfCond`` is returned as NaN on the lattice.

References:
    T. J. Ott, "The sojourn-time distribution in the M/G/1 queue with processor
    sharing", J. Appl. Prob. 21(2):360-378, 1984.
    S. F. Yashkov, "A derivation of response time distribution for an M/G/1
    processor-sharing queue", Probl. Contr. Inform. Theory 12:133-148, 1983.
    Q. Zhen, C. Knessl, "Asymptotic expansions for the sojourn time
    distribution in the M/G/1-PS queue", Math. Meth. Oper. Res. 74, 2011,
    equations (2.2)-(2.5).
"""

import numpy as np

from scipy.linalg import expm

__all__ = ["qsys_mg1_ps"]


def qsys_mg1_ps(lam, svc, svcparam, x=None, s=None, t=None, nterms=41, pdf=None):
    """Sojourn time distribution of the M/G/1-PS queue.

    Parameters
    ----------
    lam : float
        Poisson arrival rate, positive.
    svc : array_like (n,) or callable
        Phase-type initial probability vector ``alpha``, or a callable
        ``bhat(tau)`` returning the service LST, which must accept complex
        arguments.
    svcparam : array_like (n,n) or float
        Phase-type subgenerator ``T`` when ``svc`` is a vector, or the mean
        service time ``m1`` when ``svc`` is a callable.
    x : array_like, optional
        Service requirements to condition on.
    s : array_like, optional
        Transform arguments at which to tabulate the LST.
    t : array_like, optional
        Times at which to evaluate the sojourn time distribution.
    nterms : int
        Function evaluations per numerical Laplace inversion, odd.
    pdf : callable, optional
        Service density, needed to remove the conditioning when ``svc`` is a
        callable. Filled in automatically on the phase-type path.

    Returns
    -------
    dict
        With keys ``rho``, ``m1``, ``m2``, ``lstCond``, ``lstExcess``,
        ``lstUncond``, ``dominantRoot``, ``x``, ``s``, ``t``, ``lstCondVal``,
        ``lstUncondVal``, ``atomCond``, ``atomUncond``, ``meanCond``,
        ``m2Cond``, ``varCond``, ``meanUncond``, ``m2Uncond``, ``varUncond``,
        ``pdfCond``, ``cdfCond``, ``pdfUncond``, ``cdfUncond``.
    """
    lam = float(lam)
    if not np.isfinite(lam) or lam <= 0:
        raise ValueError("lambda must be a finite positive scalar")
    if nterms % 2 == 0 or nterms < 11:
        raise ValueError("nterms must be an odd integer of at least 11")

    is_ph = not callable(svc)
    if is_ph:
        alpha = np.asarray(svc, dtype=float).flatten()
        T = np.asarray(svcparam, dtype=float)
        n = alpha.size
        if T.shape != (n, n):
            raise ValueError("T must be %d x %d to match alpha" % (n, n))
        if abs(alpha.sum() - 1.0) > 1e-8 or np.any(alpha < -1e-12):
            raise ValueError("alpha must be a probability vector")
        exitrate = -T @ np.ones(n)
        if np.any(exitrate < -1e-10) or np.any(np.diag(T) >= 0):
            raise ValueError("T must be a proper phase-type subgenerator")
        negTinv = np.linalg.inv(-T)
        m1 = float(alpha @ negTinv @ np.ones(n))
        m2 = float(2.0 * alpha @ negTinv @ negTinv @ np.ones(n))

        def bhat(tau):
            return complex(alpha @ np.linalg.solve(tau * np.eye(n) - T, exitrate))

        def bpdf(y):
            return float(alpha @ expm(T * y) @ exitrate)

        # Faddeev-LeVerrier gives det(tau*I-T) and the adjugate in one sweep, so
        # bhat(tau) = nb(tau)/db(tau) as polynomials of degree n-1 and n
        db = np.zeros(n + 1)
        db[0] = 1.0
        nb = np.zeros(n)
        Mk = np.eye(n)
        for k in range(1, n + 1):
            nb[k - 1] = alpha @ Mk @ exitrate
            TM = T @ Mk
            db[k] = -np.trace(TM) / k
            Mk = TM + db[k] * np.eye(n)
    else:
        bhat = svc
        m1 = float(svcparam)
        if not np.isfinite(m1) or m1 <= 0:
            raise ValueError("the mean service time must be a finite positive scalar")
        m2 = float("nan")
        bpdf = pdf
        db = None
        nb = None

    rho = lam * m1
    if rho >= 1:
        raise ValueError("system is unstable: utilization %.6f >= 1" % rho)

    def rootfun(sv):
        return _root(sv, lam, bhat, rho)

    if is_ph:
        def denom(sv, xv):
            return _denom_ph(sv, xv, lam, rho, db, nb)
        ymax = max(40.0 * m1, 40.0 / np.min(-np.real(np.linalg.eigvals(T))))
    else:
        shift = 0.25 * (1.0 - rho) / m1

        def denom(sv, xv):
            return _denom_gen(sv, xv, lam, rho, bhat, rootfun, shift, nterms)
        ymax = 60.0 * m1

    def lst_cond(sv, xv):
        return _lst_cond(sv, xv, rho, denom, False)

    def lst_excess(sv, xv):
        return _lst_cond(sv, xv, rho, denom, True)

    has_pdf = bpdf is not None
    # fixed panelled Gauss-Legendre rule for removing the conditioning: the
    # nodes do not move with s, so the service density is sampled only once
    yq, wq = _quad_nodes(ymax)
    bq = np.array([bpdf(float(yy)) for yy in yq]) if has_pdf else None

    def lst_uncond(sv):
        if not has_pdf:
            raise ValueError("the service density is required to remove the "
                             "conditioning, pass it as the pdf argument")
        return complex(np.sum(wq * bq * np.atleast_1d(lst_cond(sv, yq))))

    x = np.atleast_1d(np.asarray([] if x is None else x, dtype=float)).flatten()
    s = np.atleast_1d(np.asarray([] if s is None else s, dtype=float)).flatten()
    t = np.atleast_1d(np.asarray([] if t is None else t, dtype=float)).flatten()

    mean_uncond = m1 / (1.0 - rho)
    res = {
        "rho": rho,
        "m1": m1,
        "m2": m2,
        "lstCond": lst_cond,
        "lstExcess": lst_excess,
        "lstUncond": lst_uncond,
        "dominantRoot": rootfun,
        "meanUncond": mean_uncond,
        "atomUncond": (1.0 - rho) * float(np.real(bhat(lam))),
        "x": x,
        "s": s,
        "t": t,
    }

    lst_cond_val = np.zeros((x.size, s.size))
    for i in range(x.size):
        for j in range(s.size):
            lst_cond_val[i, j] = np.real(lst_cond(s[j], x[i]))
    res["lstCondVal"] = lst_cond_val
    lst_uncond_val = np.full(s.size, np.nan)
    if has_pdf:
        for j in range(s.size):
            lst_uncond_val[j] = np.real(lst_uncond(s[j]))
    res["lstUncondVal"] = lst_uncond_val

    res["atomCond"] = (1.0 - rho) * np.exp(-lam * x)
    res["meanCond"] = x / (1.0 - rho)
    m2_cond = np.zeros(x.size)
    for i in range(x.size):
        m2_cond[i] = _second_moment(lambda u, xi=x[i]: lst_cond(u, xi),
                                    res["meanCond"][i], mean_uncond)
    res["m2Cond"] = m2_cond
    res["varCond"] = m2_cond - res["meanCond"] ** 2
    res["m2Uncond"] = float("nan")
    res["varUncond"] = float("nan")
    if has_pdf:
        # integrate the conditional second moment, which is far better
        # conditioned than differentiating the quadrature over the density
        m2q = np.array([_second_moment(lambda u, yy=float(yy): lst_cond(u, yy),
                                       float(yy) / (1.0 - rho), mean_uncond) for yy in yq])
        res["m2Uncond"] = float(np.sum(wq * bq * m2q))
        res["varUncond"] = res["m2Uncond"] - mean_uncond ** 2

    if t.size and not is_ph:
        raise ValueError("qsys_mg1_ps: the sojourn time distribution needs "
                         "phase-type service, since inverting a numerically "
                         "inverted transform is unstable in double precision; "
                         "with a transform handle only the LST and its moments "
                         "are available")
    pdf_cond = np.zeros((x.size, t.size))
    cdf_cond = np.zeros((x.size, t.size))
    for i in range(x.size):
        atom = float(res["atomCond"][i])
        xi = float(x[i])
        # V(x) >= x with an atom at x, so invert the excess V(x)-x net of its atom
        def gpdf(u, xi=xi, atom=atom):
            return lst_excess(u, xi) - atom

        def gcdf(u, xi=xi, atom=atom):
            return (lst_excess(u, xi) - atom) / u

        for j in range(t.size):
            if t[j] < xi:
                continue
            if t[j] == xi:
                cdf_cond[i, j] = atom
                continue
            cdf_cond[i, j] = np.real(_ilt(gcdf, t[j] - xi, nterms)) + atom
            # V(x) is atomic on the lattice (k+1)*x, where no density exists
            ratio = t[j] / xi
            if abs(ratio - round(ratio)) < 1e-9:
                pdf_cond[i, j] = np.nan
            else:
                pdf_cond[i, j] = np.real(_ilt(gpdf, t[j] - xi, nterms))
    res["pdfCond"] = pdf_cond
    res["cdfCond"] = cdf_cond

    pdf_uncond = np.full(t.size, np.nan)
    cdf_uncond = np.full(t.size, np.nan)
    if has_pdf:
        for j in range(t.size):
            if t[j] <= 0:
                pdf_uncond[j] = 0.0
                cdf_uncond[j] = 0.0
                continue
            pdf_uncond[j] = np.real(_ilt(lst_uncond, t[j], nterms))
            cdf_uncond[j] = np.real(_ilt(lambda u: lst_uncond(u) / u, t[j], nterms))
    res["pdfUncond"] = pdf_uncond
    res["cdfUncond"] = cdf_uncond
    return res


def _lst_cond(s, x, rho, denom, excess):
    """E[exp(-s V(x))] = (1-rho)/D(s,x), or the transform of V(x)-x. D carries
    its dominant exponential separately so that neither factor overflows."""
    xa = np.asarray(x, dtype=float)
    sh = s if excess else 0.0
    if xa.ndim == 0:
        if xa == 0:
            return 1.0 + 0j
        val, scale = denom(s, float(xa))
        return (1.0 - rho) * np.exp((sh - scale) * float(xa)) / val
    val, scale = denom(s, xa)
    out = (1.0 - rho) * np.exp((sh - scale) * xa) / val
    return np.where(xa == 0, 1.0 + 0j, out)


def _root(s, lam, bhat, rho):
    """Unique root of tau = s + lam*(1-bhat(tau)) in the right half plane."""
    tau = complex(s)
    maxit = max(200, int(np.ceil(3.0 * np.log(1e-15) / np.log(max(rho, 1e-3)))))
    for _ in range(maxit):
        taunew = s + lam * (1.0 - bhat(tau))
        if abs(taunew - tau) <= 1e-14 * max(1.0, abs(taunew)):
            return taunew
        tau = taunew
    raise RuntimeError("qsys_mg1_ps: the dominant root iteration did not converge")


def _denom_ph(s, x, lam, rho, db, nb):
    """D(s,x) = exp(scale*x)*val for phase-type service, exactly, and for a
    whole vector of service requirements at once."""
    n = nb.size
    dbc = np.asarray(db, dtype=complex)
    nbc = np.asarray(nb, dtype=complex)
    dm = dbc - np.concatenate(([0.0], nbc))
    P = np.convolve(np.array([1.0, -(s + lam)], dtype=complex), dbc)
    P = P + np.concatenate(([0.0, 0.0], lam * nbc))
    A = np.convolve(np.array([1.0 - rho, s * rho, 0.0], dtype=complex), dbc) \
        - np.concatenate(([0.0], np.convolve(np.array([1.0 - rho, s], dtype=complex), lam * dm)))
    if np.linalg.norm(A[-2:]) > 1e-6 * max(1.0, np.linalg.norm(A)):
        raise RuntimeError("qsys_mg1_ps: the double pole at the origin did not cancel")
    Ahat = A[:n + 1]
    r = np.roots(P)
    scale = float(np.max(r.real))
    xv = np.atleast_1d(np.asarray(x, dtype=float))
    sep = np.abs(r[:, None] - r[None, :]) + np.diag(np.full(r.size, np.inf))
    if np.all(sep.min(axis=1) > 1e-7 * max(1.0, np.max(np.abs(r)))):
        coef = np.polyval(Ahat, r) / np.polyval(np.polyder(P), r)
        val = np.exp(np.outer(xv, r - scale)) @ coef
    else:
        # repeated poles: use the companion realization of Ahat/P instead
        Pn = P / P[0]
        Ac = np.zeros((n + 1, n + 1), dtype=complex)
        Ac[0, :] = -Pn[1:]
        Ac[1:, :n] = np.eye(n)
        e1 = np.zeros(n + 1, dtype=complex)
        e1[0] = 1.0
        cv = Ahat / P[0]
        val = np.array([cv @ expm((Ac - scale * np.eye(n + 1)) * xx) @ e1 for xx in xv])
    return (complex(val[0]) if np.isscalar(x) or np.ndim(x) == 0 else val), scale


def _denom_gen(s, x, lam, rho, bhat, rootfun, shift, nterms):
    """D(s,x) for a general service transform, by inverting f(tau;s) in tau."""
    scale = float(np.real(rootfun(s))) + shift
    return _ilt(lambda u: _f(u + scale, s, lam, rho, bhat), x, nterms), scale


def _f(tau, s, lam, rho, bhat):
    bh = bhat(tau)
    num = (1.0 - rho) * tau ** 2 - (1.0 - rho) * lam * (1.0 - bh) * tau \
        + s * rho * tau - s * lam * (1.0 - bh)
    return num / (tau ** 2 * (tau - s - lam * (1.0 - bh)))


def _ilt(fun, t, nterms):
    """Abate-Whitt Euler inversion, symmetrized so that complex-valued time
    functions are handled as well as real-valued ones."""
    ne = (nterms - 1) // 2
    eta = np.zeros(2 * ne + 1)
    eta[0] = 0.5
    eta[1:ne + 1] = 1.0
    eta[2 * ne] = 2.0 ** (-ne)
    from math import lgamma, log
    for k in range(1, ne):
        eta[2 * ne - k] = eta[2 * ne - k + 1] + np.exp(
            lgamma(ne + 1) - ne * log(2.0) - lgamma(k + 1) - lgamma(ne - k + 1))
    k = np.arange(2 * ne + 1)
    beta = ne * np.log(10.0) / 3.0 + 1j * np.pi * k
    eta = 10.0 ** (ne / 3.0) * (1 - (k % 2) * 2) * eta
    g = 0.0 + 0j
    for j in range(k.size):
        bj = beta[j] / t
        g += 0.5 * eta[j] * (fun(bj) + fun(np.conj(bj)))
    return g / t


def _second_moment(lst, meanref, meanscale):
    """Second moment from the transform curvature at the origin. The stencil is
    one-sided so that the transform is never sampled at negative arguments,
    where it need not converge, the step is scaled by the conditional mean but
    capped by the unconditional one so that it stays finite as x -> 0, and
    Richardson extrapolation over h and h/2 removes the leading truncation."""
    if meanref <= 0:
        return 0.0
    h = min(1e-2 / meanref, 1.0 / meanscale)
    d1 = _d2_forward(lst, h)
    d2 = _d2_forward(lst, 0.5 * h)
    return float((16.0 * d2 - d1) / 15.0)


def _d2_forward(lst, h):
    f = np.array([np.real(lst(j * h)) for j in range(6)])
    return (45 * f[0] - 154 * f[1] + 214 * f[2] - 156 * f[3]
            + 61 * f[4] - 10 * f[5]) / (12 * h ** 2)


def _quad_nodes(ymax, npanel=8, ng=32):
    """Panelled Gauss-Legendre rule on [0,ymax], with the panels growing
    geometrically so that both ends of an exponentially decaying density are
    resolved. The rule is fixed, so it is identical in every codebase."""
    xg, wg = np.polynomial.legendre.leggauss(ng)
    edges = np.concatenate(([0.0], ymax * 2.0 ** np.arange(-npanel, 1)))
    y = np.zeros((edges.size - 1) * ng)
    w = np.zeros_like(y)
    for k in range(edges.size - 1):
        a, b = edges[k], edges[k + 1]
        y[k * ng:(k + 1) * ng] = 0.5 * (a + b) + 0.5 * (b - a) * xg
        w[k * ng:(k + 1) * ng] = 0.5 * (b - a) * wg
    return y, w
