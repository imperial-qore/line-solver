"""Quantity of work in a load-dependent processor sharing station with blocking.

Stationary distribution of the workload of the single-stage generalized
processor sharing model with Poisson arrivals and blocking, per J.W. Cohen,
"The multiple phase service network with generalized processor sharing",
Acta Informatica 12, 245-284 (1979), Sect. 9, eqs. (9.1)-(9.3).
"""

import numpy as np

__all__ = ['qsys_ldps_workload']

DEFAULT_NGRID = 2001


def qsys_ldps_workload(lambd, B, alpha, N, t=None, ngrid=None):
    """Stationary distribution of the quantity of work in a single-stage
    load-dependent processor sharing station with Poisson arrivals and blocking.

    Assumed model (Cohen 1979, Sect. 9; the model of Sect. 7 with one stage):
    a single service stage fed by a Poisson arrival stream of rate lambd; a
    blocking capacity N, so that a request arriving when N requests are already
    present is lost and leaves no trace on the state; generalized processor
    sharing, so that when x requests are present each accrues service at rate
    f(x) and the stage completes work at total rate x*f(x); and required service
    times i.i.d. with absolutely continuous distribution B of finite mean beta.
    The station is parametrized by the LINE load-dependent total rate scaling
    alpha(x)=x*f(x), the argument of setLoadDependence at a PS station.

    With psi the total amount of service still to be given to the requests
    present, Cohen eqs. (9.1)-(9.3) give

        Pr{psi < y} = sum_{h=0}^{N} p_h Psi^{h*}(y)
        p_h = (rho^h/h!) phi(h) / sum_k (rho^k/k!) phi(k),   rho = lambd*beta
        phi(h) = 1/prod_{k=1}^{h} f(k),   phi(0)=1
        Psi(y) = int_0^y (1-B(v))/beta dv

    with Psi^{h*} the h-fold convolution of Psi and Psi^{0*} degenerate at zero.
    Substituting f(k)=alpha(k)/k the factorial cancels, leaving p_h proportional
    to rho^h/prod_{k=1}^{h} alpha(k), the familiar load-dependent birth-death
    form. Psi is the equilibrium (residual life) distribution of B, so psi is a
    mixture of h-fold convolutions of residual service times with an atom p_0 at
    zero.

    This is the model of Cohen (1979) Sect. 9 only. It is not the weighted
    GPS/DPS discipline of SchedStrategy.GPS, whose per-class weights this
    formula does not represent.

    Args:
        lambd: rate of the Poisson arrival stream (finite, positive)
        B: required service time Distribution (continuous, finite positive mean)
        alpha: rate scaling alpha(n)=n*f(n) for n=1..N (finite, positive)
        N: blocking capacity (finite positive integer)
        t: optional grid at which the CDF is returned. Default: an automatically
           sized grid covering the bulk of the distribution.
        ngrid: optional number of points of the internal uniform quadrature grid
           on which the convolutions are formed. Default 2001. Accuracy is
           second order in the step for an absolutely continuous B, the case
           Cohen assumes, and falls back to first order when B has an atom so
           that 1-B is discontinuous; raise ngrid for those.

    Returns:
        (F, t, p) where F[j] = Pr{psi <= t[j]}, t is the grid F is reported on,
        and p[h] = Pr{x = h} for h = 0..N. Note F[0] = p[0] = Pr{psi = 0} when
        t[0] = 0, since the workload has an atom at zero.
    """
    # gating: reject anything outside the assumed model
    if not np.isscalar(lambd) or not np.isfinite(lambd) or lambd <= 0:
        raise ValueError('lambd must be a finite positive scalar, the rate of '
                         'the Poisson arrival stream.')
    if not hasattr(B, 'getMean') or not hasattr(B, 'evalCDF'):
        raise ValueError('B must be a Distribution giving the required service time.')
    if B.isDisabled() or B.isImmediate():
        raise ValueError('B must be an active service time distribution, not '
                         'Disabled or Immediate.')
    if not B.isContinuous():
        raise ValueError('B must be a continuous distribution: Cohen (1979) '
                         'Sect. 9 assumes an absolutely continuous required '
                         'service time.')
    beta = B.getMean()
    if not np.isfinite(beta) or beta <= 0:
        raise ValueError('B must have a finite positive mean.')
    if int(N) != N or N < 1:
        raise ValueError('N must be a finite positive integer, the blocking '
                         'capacity of the service stage.')
    N = int(N)
    alpha = np.asarray(alpha, dtype=float).ravel()
    if alpha.size < N:
        raise ValueError('alpha must supply the rate scaling for n=1..%d, but '
                         'only %d entries were given.' % (N, alpha.size))
    alpha = alpha[:N]
    if not np.all(np.isfinite(alpha)) or np.any(alpha <= 0):
        raise ValueError('alpha(n) must be finite and strictly positive for '
                         'n=1..N, since every request in a busy stage is served '
                         'at a positive rate.')
    if ngrid is None:
        ngrid = DEFAULT_NGRID
    if int(ngrid) != ngrid or ngrid < 2:
        raise ValueError('ngrid must be an integer of at least 2.')
    ngrid = int(ngrid)

    # stationary number in system, eq. (9.1). Accumulated in logs so that large
    # rho or large N do not overflow before normalization.
    rho = lambd * beta
    logw = np.concatenate(([0.0], np.cumsum(np.log(rho) - np.log(alpha))))
    logw = logw - np.max(logw)
    w = np.exp(logw)
    p = w / np.sum(w)

    # see _kb/03-api-layer.md for rationale
    m1e = beta * (1 + B.getSCV()) / 2
    if not np.isfinite(m1e) or m1e <= 0:
        raise ValueError('B must have a finite second moment: the equilibrium '
                         'residual service time is otherwise undefined.')
    user_grid = t is not None and np.size(t) > 0
    if user_grid:
        t = np.asarray(t, dtype=float).ravel()
        if not np.all(np.isfinite(t)) or np.any(t < 0):
            raise ValueError('t must be a vector of finite non-negative times.')
        tmax = np.max(t)
        if tmax <= 0:
            tmax = m1e
    else:
        nz = np.nonzero(p > 1e-12)[0]
        hmax = int(nz[-1]) if nz.size > 0 else 1
        if hmax < 1:
            hmax = 1
        tmax = m1e * (hmax + 8 * np.sqrt(hmax))
        tmax = max(tmax, 8 * m1e)

    # Internal uniform grid: the convolutions below need a constant step. The
    # result is interpolated back onto t when the caller supplied one.
    tg = np.linspace(0.0, tmax, ngrid)
    dt = tg[1] - tg[0]

    # see _kb/03-api-layer.md for rationale
    e = np.array([(1 - B.evalCDF(float(v))) / beta for v in tg])

    # workload distribution, eq. (9.2). The h=0 term is degenerate at zero,
    # contributing the atom p_0 over the whole non-negative grid.
    Fg = np.full(ngrid, p[0])
    dens = e.copy()
    for h in range(1, N + 1):
        if h > 1:
            dens = _conv_trap(dens, e, dt)
        Fg = Fg + p[h] * _cumtrapz(tg, dens)

    if user_grid:
        return np.interp(t, tg, Fg), t, p
    return Fg, tg, p


def _conv_trap(f, g, dt):
    """Convolution of two densities sampled on a uniform grid, using the
    trapezoidal rule rather than the rectangle rule implied by a bare
    np.convolve. The correction matters here because the equilibrium density
    does not vanish at the origin: e(0)=1/beta. Writing t_i=i*dt,

        (f*g)(t_i) = int_0^{t_i} f(v) g(t_i-v) dv
                   ~ dt*[ sum_{j=0}^{i} f_j g_{i-j} - (f_0 g_i + f_i g_0)/2 ]

    i.e. the raw convolution less half of each endpoint. Without the correction
    each convolution over-counts by dt*f_0*g_i, which accumulates over h and
    drives the mixture CDF above one.
    """
    n = f.size
    c = np.convolve(f, g)[:n] * dt
    return c - dt * (f[0] * g[:n] + f[:n] * g[0]) / 2


def _cumtrapz(x, y):
    """Cumulative trapezoidal integral of y over the grid x, starting at zero."""
    c = np.zeros(x.size)
    c[1:] = np.cumsum(np.diff(x) * (y[1:] + y[:-1]) / 2)
    return c
