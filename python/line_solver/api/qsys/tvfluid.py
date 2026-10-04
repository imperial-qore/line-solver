"""
Time-varying many-server fluid queue and network (Liu and Whitt).

Native Python twin of matlab/src/api/qsys/qsys_gtmtst_fluid.m and
matlab/src/api/npfqn/npfqn_gtmtst_fluid.m: the Gt/Mt/st+GI fluid queue of
Y. Liu and W. Whitt, Queueing Systems 71 (2012), 405-444, with the algorithms
and the network fixed point of Y. Liu and W. Whitt, INFORMS Journal on
Computing 26(1) (2014), 59-73.
"""

from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np


def _grid_eval(f, t: np.ndarray) -> np.ndarray:
    """Evaluate a handle on a grid, accepting an array-aware or a scalar one."""
    if not callable(f):
        return np.full(t.shape, float(f))
    try:
        y = np.asarray(f(t), dtype=float)
        if y.shape == t.shape:
            return y
        if y.size == 1:
            return np.full(t.shape, float(y))
    except Exception:
        pass
    return np.array([float(f(float(ti))) for ti in t])


def qsys_gtmtst_fluid(lambdaFun, sFun, muFun, patienceCcdf: Callable[[Any], Any],
                      T: float, dt: float = None, B0: float = 0.0, w0: float = 0.0,
                      sPrimeFun=None, patiencePdf: Optional[Callable[[Any], Any]] = None,
                      lambdaPast=None) -> Dict[str, Any]:
    """
    The Gt/Mt/st+GI many-server fluid queue.

    Time-varying arrival rate ``lambda(t)``, time-varying staffing ``s(t)``,
    exponential service at the time-varying rate ``mu(t)``, general patience with
    complementary cdf ``F^c``, and unlimited waiting room.

    THE MODEL ALTERNATES BETWEEN TWO REGIMES and the whole algorithm is the
    bookkeeping of that alternation:

    * UNDERLOADED: the queue is empty and every arrival enters service at once,
      so the system is the infinite-server fluid model and ``B`` obeys
      ``B'(t) = lambda(t) - mu(t)B(t)`` (eq. 18 of the reference, in its Mt form).
      It ends when ``B`` reaches ``s`` while ``lambda`` exceeds the rate
      ``Gamma(t) = s'(t) + s(t)mu(t)`` at which capacity frees up (eq. 15).
    * OVERLOADED: every server is busy, ``B(t) = s(t)``, fluid enters service at
      exactly ``Gamma(t)``, and the queue is described by its BOUNDARY WAITING
      TIME ``w(t)``, the age of the oldest fluid still waiting. Content of age
      ``x`` is what arrived ``x`` ago and has not yet abandoned,
      ``q(t,x) = lambda(t-x)F^c(x)``, and the boundary moves by the delay
      differential equation (eq. 21)

          w'(t) = 1 - Gamma(t) / [lambda(t-w(t)) F^c(w(t))].

      It ends when ``w`` returns to 0 with ``lambda`` no longer above ``Gamma``
      (eq. 14).

    WHY w AND NOT Q. The queue content is a functional of ``w``, but not the
    other way round: two systems with the same ``Q`` and different age profiles
    abandon at different rates. Tracking the boundary keeps the age profile
    exact, which is what makes a general patience law admissible at all.

    Args:
        lambdaFun: arrival rate lambda(t)
        sFun: staffing s(t), a positive function or a constant
        muFun: service rate mu(t), a function or a constant
        patienceCcdf: F^c(x) = P(patience > x)
        T: horizon; the model is solved on [0,T]
        dt: grid step, default T/2000
        B0: fluid in service at time 0
        w0: boundary waiting time at time 0, 0 for an empty queue
        sPrimeFun: s'(t); differentiated numerically from sFun when absent
        patiencePdf: the patience density, for the abandonment rate; differenced
            from the ccdf when absent
        lambdaPast: the arrival rate before time 0, needed only when the queue
            starts non-empty; defaults to lambdaFun evaluated at negative times

    Returns:
        Dict on the grid: ``times``, ``regime`` (1 overloaded, 0 underloaded),
        ``B`` (fluid in service), ``Q`` (fluid in queue), ``X = B+Q``, ``w``
        (boundary waiting time), ``v`` (potential waiting time), ``sigma``
        (service completion rate), ``alpha`` (abandonment rate), ``utilization``
        (B/s), ``arrivalRate``, ``staffing``, ``capacityRate`` (Gamma).

    References:
        Y. Liu, W. Whitt (2012). The Gt/GI/st+GI many-server fluid queue.
        Queueing Systems 71, 405-444; Y. Liu, W. Whitt (2014). Algorithms for
        time-varying networks of many-server fluid queues. INFORMS Journal on
        Computing 26(1), 59-73.
    """
    if T <= 0:
        raise ValueError('The horizon T must be positive.')
    if dt is None:
        dt = T / 2000.0
    n = int(round(T / dt)) + 1
    t = np.linspace(0.0, T, n)
    dt = float(t[1] - t[0])

    lam = _grid_eval(lambdaFun, t)
    s = _grid_eval(sFun, t)
    mu = _grid_eval(muFun, t)
    if np.any(s <= 0):
        raise ValueError('The staffing function must be positive.')
    if np.any(mu <= 0):
        raise ValueError('The service rate must be positive.')
    if sPrimeFun is not None:
        sp = _grid_eval(sPrimeFun, t)
    else:
        sp = np.gradient(s, dt)

    lam_at = lambdaFun if callable(lambdaFun) else (lambda u, _c=float(lambdaFun): _c)
    past = lambdaPast if lambdaPast is not None else lam_at

    def lam_of(u: float) -> float:
        """The arrival rate at a possibly negative time."""
        return float(past(u)) if u < 0 else float(lam_at(u))

    fc = lambda x: float(np.asarray(patienceCcdf(x)))
    if patiencePdf is not None:
        fpdf = lambda x: float(np.asarray(patiencePdf(x)))
    else:
        h = 1e-6
        fpdf = lambda x: max(0.0, (fc(max(0.0, x - h)) - fc(x + h)) / (2 * h))

    B = np.zeros(n)
    Q = np.zeros(n)
    w = np.zeros(n)
    alpha = np.zeros(n)
    regime = np.zeros(n, dtype=int)
    B[0] = B0
    w[0] = w0
    gamma = sp + s * mu                       # Gamma(t), eq. (13)

    def queue_from_w(i: int, wi: float) -> float:
        """Q(t) = int_0^w lambda(t-x)F^c(x)dx, the content that has not abandoned."""
        if wi <= 0:
            return 0.0
        m = max(8, int(np.ceil(wi / dt)) + 1)
        x = np.linspace(0.0, wi, m + 1 if m % 2 == 0 else m + 2)
        vals = np.array([lam_of(t[i] - xx) * fc(xx) for xx in x])
        wgt = np.ones(x.size)
        wgt[1:-1:2] = 4.0
        wgt[2:-1:2] = 2.0
        return float((x[-1] - x[0]) / (3.0 * (x.size - 1)) * np.sum(wgt * vals))

    def abandon_from_w(i: int, wi: float) -> float:
        """alpha(t) = int_0^w lambda(t-x) f(x) dx, the fluid whose patience expires."""
        if wi <= 0:
            return 0.0
        m = max(8, int(np.ceil(wi / dt)) + 1)
        x = np.linspace(0.0, wi, m + 1 if m % 2 == 0 else m + 2)
        vals = np.array([lam_of(t[i] - xx) * fpdf(xx) for xx in x])
        wgt = np.ones(x.size)
        wgt[1:-1:2] = 4.0
        wgt[2:-1:2] = 2.0
        return float((x[-1] - x[0]) / (3.0 * (x.size - 1)) * np.sum(wgt * vals))

    # Initial regime: overloaded when the queue is already occupied, or when the
    # servers are full and the arrival rate beats the rate capacity frees up.
    over = w[0] > 0 or (B[0] >= s[0] - 1e-12 and lam[0] > gamma[0])
    regime[0] = 1 if over else 0
    if over:
        B[0] = s[0]
    Q[0] = queue_from_w(0, w[0]) if over else 0.0
    alpha[0] = abandon_from_w(0, w[0]) if over else 0.0

    for i in range(n - 1):
        if regime[i] == 0:
            # Underloaded: B' = lambda - mu B, by RK4 on the grid step.
            def f(tt, bb):
                return float(np.interp(tt, t, lam)) - float(np.interp(tt, t, mu)) * bb
            k1 = f(t[i], B[i])
            k2 = f(t[i] + dt / 2, B[i] + dt * k1 / 2)
            k3 = f(t[i] + dt / 2, B[i] + dt * k2 / 2)
            k4 = f(t[i] + dt, B[i] + dt * k3)
            Bnext = B[i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            wnext = 0.0
            if Bnext >= s[i + 1] and lam[i + 1] > gamma[i + 1]:
                # The servers just filled and the input outruns the freed
                # capacity: eq. (15), the underloaded interval ends here.
                Bnext = s[i + 1]
                regime[i + 1] = 1
            else:
                regime[i + 1] = 0
                Bnext = min(Bnext, s[i + 1])
        else:
            # Overloaded: B = s and the boundary moves by eq. (21).
            def g(tt, ww):
                den = lam_of(tt - ww) * fc(ww)
                if den <= 0:
                    # No fluid of that age survives, so the boundary can only
                    # advance with the clock.
                    return 1.0
                return 1.0 - float(np.interp(tt, t, gamma)) / den
            k1 = g(t[i], w[i])
            k2 = g(t[i] + dt / 2, max(0.0, w[i] + dt * k1 / 2))
            k3 = g(t[i] + dt / 2, max(0.0, w[i] + dt * k2 / 2))
            k4 = g(t[i] + dt, max(0.0, w[i] + dt * k3))
            wnext = w[i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            Bnext = s[i + 1]
            if wnext <= 0 and lam[i + 1] <= gamma[i + 1]:
                # The queue has drained and the input no longer outruns the
                # freed capacity: eq. (14), the overloaded interval ends here.
                wnext = 0.0
                regime[i + 1] = 0
            else:
                wnext = max(wnext, 0.0)
                regime[i + 1] = 1
        B[i + 1] = Bnext
        w[i + 1] = wnext
        Q[i + 1] = queue_from_w(i + 1, wnext) if regime[i + 1] == 1 else 0.0
        alpha[i + 1] = abandon_from_w(i + 1, wnext) if regime[i + 1] == 1 else 0.0

    sigma = mu * B                                  # service completion rate, eq. (3)
    # The potential waiting time of an arrival at t is the u-t at which the
    # boundary reaches it, i.e. the solution of u - w(u) = t. That map is
    # non-decreasing, so one interpolation inverts it.
    entry = t - w
    v = np.maximum(0.0, np.interp(t, entry, t, left=t[0], right=t[-1]) - t)

    return {
        'times': t,
        'regime': regime,
        'B': B,
        'Q': Q,
        'X': B + Q,
        'w': w,
        'v': v,
        'sigma': sigma,
        'alpha': alpha,
        'utilization': B / s,
        'arrivalRate': lam,
        'staffing': s,
        'capacityRate': gamma,
    }


def npfqn_gtmtst_fluid(lambdaFuns: Sequence, sFuns: Sequence, muFuns: Sequence,
                       patienceCcdfs: Sequence, P, T: float, dt: float = None,
                       B0: Optional[Sequence[float]] = None,
                       w0: Optional[Sequence[float]] = None,
                       tol: float = 1e-6, maxIter: int = 100) -> Dict[str, Any]:
    """
    A time-varying open network of many-server fluid queues with abandonment.

    Each queue is the Gt/Mt/st+GI fluid queue of :func:`qsys_gtmtst_fluid`; the
    departure flow of queue i is routed to queue j with the (possibly
    time-varying) proportion ``P[i][j]``, whatever is left leaving the network.

    THE NETWORK IS A FIXED POINT. The total arrival rate of queue j is
    ``lambda_j(t) = lambda_j^0(t) + sum_i sigma_i(t) P_ij(t)`` with
    ``sigma_i = mu_i B_i`` the service completion rate (eqs. 23-24), and
    ``sigma_i`` itself depends on ``lambda_i``. The iteration starts from the
    external rates alone and adds one more traversal of the network per round,
    so the nth iterate is the fluid that has made n transitions; the map is a
    monotone contraction, so the rates increase to the fixed point.

    Args:
        lambdaFuns: external arrival rate of each queue
        sFuns: staffing of each queue
        muFuns: service rate of each queue
        patienceCcdfs: patience ccdf of each queue
        P: routing proportions, either an m x m array or a callable P(t)
            returning one
        T: horizon
        dt: grid step
        B0: initial fluid in service at each queue
        w0: initial boundary waiting time at each queue
        tol: sup-norm tolerance on the arrival-rate iteration
        maxIter: cap on the iterations

    Returns:
        Dict with ``times``, ``queues`` (the per-queue dicts of
        :func:`qsys_gtmtst_fluid`), ``arrivalRates`` (the converged total rates,
        one row per queue), ``iterations`` and ``residual``.

    References:
        Y. Liu, W. Whitt (2014). Algorithms for time-varying networks of
        many-server fluid queues. INFORMS Journal on Computing 26(1), 59-73.
    """
    m = len(lambdaFuns)
    if not (len(sFuns) == len(muFuns) == len(patienceCcdfs) == m):
        raise ValueError('Every queue needs an arrival rate, a staffing, a service rate and a '
                         'patience law.')
    if T <= 0:
        raise ValueError('The horizon T must be positive.')
    if dt is None:
        dt = T / 2000.0
    n = int(round(T / dt)) + 1
    t = np.linspace(0.0, T, n)
    B0 = list(B0) if B0 is not None else [0.0] * m
    w0 = list(w0) if w0 is not None else [0.0] * m

    ext = np.array([_grid_eval(f, t) for f in lambdaFuns])
    if callable(P):
        Pgrid = np.array([np.asarray(P(float(ti)), dtype=float) for ti in t])   # n x m x m
    else:
        Pfixed = np.asarray(P, dtype=float)
        if Pfixed.shape != (m, m):
            raise ValueError('The routing matrix must be m x m.')
        Pgrid = np.repeat(Pfixed[None, :, :], n, axis=0)
    if np.any(Pgrid < -1e-12) or np.any(Pgrid.sum(axis=2) > 1 + 1e-9):
        raise ValueError('The routing matrix must be substochastic.')

    lam = ext.copy()
    queues: List[Dict[str, Any]] = []
    residual = np.inf
    it = 0
    for it in range(1, maxIter + 1):
        queues = []
        sigma = np.zeros((m, n))
        for i in range(m):
            lam_i = lam[i]
            fun = lambda u, _v=lam_i: float(np.interp(u, t, _v, left=_v[0], right=_v[-1]))
            res = qsys_gtmtst_fluid(fun, sFuns[i], muFuns[i], patienceCcdfs[i], T, dt,
                                    B0[i], w0[i])
            queues.append(res)
            sigma[i] = res['sigma']
        # lambda_j = lambda_j^0 + sum_i sigma_i P_ij, eqs. (23)-(24).
        newlam = ext.copy()
        for j in range(m):
            newlam[j] += np.einsum('kn,nk->n', sigma, Pgrid[:, :, j])
        residual = float(np.max(np.abs(newlam - lam)))
        lam = newlam
        if residual < tol:
            break

    return {
        'times': t,
        'queues': queues,
        'arrivalRates': lam,
        'iterations': it,
        'residual': residual,
    }
