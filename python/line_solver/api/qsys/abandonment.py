"""
Multiserver queues with customer abandonment.

Native Python twin of matlab/src/api/qsys/qsys_mgisrgi_whitt.m and
qsys_erlanga.m: the engineering solution of the call-center model
M/GI/s/r+GI of W. Whitt (2005), Management Science 51(2), 221-235, and the
Erlang A model M/M/s/r+M it contains as an exact special case.
"""

from typing import Any, Callable, Dict, Optional, Sequence, Union

import numpy as np

Patience = Union[float, Callable[[float], float], Dict[str, Callable[[float], float]]]


def _resolve_patience(patience: Patience):
    """
    Resolve the three accepted forms of the patience argument into a hazard
    handle or a ccdf handle, plus the exponential flag and its rate.
    """
    if isinstance(patience, dict):
        if 'ccdf' in patience:
            return None, patience['ccdf'], False, float('nan')
        if 'hazard' in patience:
            return patience['hazard'], None, False, float('nan')
        raise ValueError("patience dict must carry key 'ccdf' or 'hazard'")
    if callable(patience):
        return patience, None, False, float('nan')
    theta = float(patience)
    if theta < 0:
        raise ValueError('the patience rate theta must be non-negative')
    return (lambda t, _th=theta: _th), None, True, theta


def _rates_step(j: int, lambda_val: float, delta_prev: float, hazard_fun, ccdf_fun):
    """
    One step of eqs. (3.3)-(3.4) (hazard form) or (3.5)-(3.6) (ccdf form).

    DIVERGENCE from the printed eqs. (3.5)-(3.6): they read
    ``delta_j = int_{(j-1)/lambda}^{j/lambda} h(t) dt`` and
    ``Delta_k = -log F^c(k/lambda)``, which are cumulative hazards, i.e.
    dimensionless, while delta and Delta are rates everywhere else in the paper.
    They are the AVERAGE hazard over an interval of length 1/lambda, so the
    factor lambda is missing. Restoring it makes the ccdf form reduce to the
    exact Erlang A rates for exponential patience, which the paper states this
    approximation does (eq. 7.12); the literal form gives theta/lambda instead
    of theta and is wrong by that factor.
    """
    if ccdf_fun is None:
        delta_j = float(hazard_fun(j / lambda_val))
        delta_tot = delta_prev + delta_j
    else:
        g = float(ccdf_fun(j / lambda_val))
        if g <= 0:
            raise ValueError(
                'the patience ccdf vanishes at t = %g, so every customer has abandoned by '
                'then; supply a hazard handle instead' % (j / lambda_val))
        delta_tot = -lambda_val * np.log(g)
        delta_j = delta_tot - delta_prev
    if delta_j < 0:
        raise ValueError('the patience law produced a negative abandonment rate')
    return delta_j, delta_tot


def _kernel(k: int, smu: float, dlt: np.ndarray, delta: np.ndarray):
    """
    Eqs. (7.10)-(7.11): with k waiting, the total departure rate before the jth
    departure epoch is s*mu + Delta_k - Delta_{j-1}, of which delta_j is the
    share belonging to the customer of interest.
    """
    j = np.arange(1, k + 1)
    rate = smu + dlt[k] - dlt[j - 1]
    return delta[j - 1] / rate, 1.0 / rate


def _ratio(num: float, den: float) -> float:
    """A conditional moment is 0/0 when the conditioning event cannot happen."""
    return 0.0 if den <= 0 else num / den


def _transform(z: complex, w_arr: np.ndarray, sigma: np.ndarray,
               kernels, served: bool) -> complex:
    """
    Eqs. (7.22)-(7.23) when served, eqs. (7.32)-(7.33) otherwise. Both fold the
    same per-position kernel: the wait is a sum of exponentials with rates
    1/m_k(j), truncated at the departure epoch that serves or loses the customer.

    ``kernels[k-1]`` carries the rates, the abandonment shares phi_k(j) and the
    survival products prod_{l<j}(1-phi_k(l)), all precomputed, because the
    inversion evaluates this at dozens of nodes per time point.
    """
    val = 0.0 + 0.0j
    for k in range(1, len(w_arr) + 1):
        rate, phik, surv = kernels[k - 1]
        factor = rate / (rate + z)
        if served:
            val += w_arr[k - 1] * sigma[k - 1] * np.prod(factor)
        else:
            val += w_arr[k - 1] * np.sum(surv * phik * np.cumprod(factor))
    return val


def qsys_mgisrgi_whitt(lambda_val: float, mu: float, s: int, r: float,
                       patience: Patience,
                       wPoints: Optional[Sequence[float]] = None,
                       maxQueue: int = 100000, tol: float = 1e-14,
                       invMethod: str = 'euler', invN: int = 41) -> Dict[str, Any]:
    """
    Engineering solution of the M/GI/s/r+GI queue.

    Poisson arrivals at rate ``lambda_val``, iid general service times of mean
    ``1/mu``, ``s`` servers, ``r`` extra waiting spaces and iid patience times
    with a general distribution.

    The general patience law is replaced by state-dependent Markovian
    abandonment, a customer jth from the end of the queue abandoning at rate
    ``delta_j = h(j/lambda)`` for the patience hazard ``h`` (eq. 3.3), because
    such a customer has been waiting for about ``j/lambda``; the general service
    law is replaced by an exponential of the same mean (Section 5). What is left
    is a birth-and-death process, solved exactly. Only the hazard NEAR THE ORIGIN
    matters, not the mean or the tail of the patience law.

    Args:
        lambda_val: arrival rate
        mu: service rate of one server, the reciprocal of the mean service time
        s: number of servers
        r: extra waiting spaces, ``float('inf')`` for an unbounded queue
        patience: scalar rate (exponential patience, then the answer is exact
            and the model is Erlang A), a callable hazard ``h(t)``, or a dict
            ``{'ccdf': G}`` using the integrated form of eq. (3.6)
        wPoints: times at which to return the waiting-time cdfs
        maxQueue: truncation level used when ``r`` is infinite
        tol: relative tail tolerance for that truncation
        invMethod: Laplace inversion method for the cdfs
        invN: number of inversion nodes

    Returns:
        Dict with the steady-state distribution ``queueLengthDist``, the
        probabilities ``probLoss``/``probNoWait``/``probServed``/``probAbandon``,
        the moments ``meanNumber``/``varNumber``/``meanQueueLength``/
        ``varQueueLength``/``meanWaitServed``/``varWaitServed``/
        ``meanWaitAbandon``/``varWaitAbandon``/``meanWait``/``secondMomentWait``,
        the rates ``utilization``/``throughput``/``abandonRate``, the
        abandonment rates ``abandonRates``/``totalAbandonRates``, and, when
        ``wPoints`` is given, ``cdfWaitServed``/``cdfWaitAbandon``/``cdfWait``.

    References:
        W. Whitt (2005). Engineering solution of a basic call-center model.
        Management Science 51(2), 221-235.
    """
    if lambda_val <= 0:
        raise ValueError('The arrival rate lambda must be positive.')
    if mu <= 0:
        raise ValueError('The service rate mu must be positive.')
    s = int(round(s))
    if s < 1:
        raise ValueError('The number of servers s must be at least 1.')
    if r < 0:
        raise ValueError('The number of extra waiting spaces r must be non-negative.')

    hazard_fun, ccdf_fun, is_exponential, theta = _resolve_patience(patience)

    finite_r = np.isfinite(r)
    rr = int(round(r)) if finite_r else int(maxQueue)

    # The birth-death recursion of eqs. (7.4)-(7.7), unnormalized with x_s = 1.
    x_up = np.zeros(rr + 1)
    x_up[0] = 1.0
    dlt = np.zeros(rr + 1)
    delta = np.zeros(max(rr, 1))
    smu = s * mu
    k_used = rr
    peak = 1.0
    for k in range(rr):
        j = k + 1
        delta[j - 1], dlt[j] = _rates_step(j, lambda_val, dlt[j - 1], hazard_fun, ccdf_fun)
        x_up[k + 1] = lambda_val * x_up[k] / (smu + dlt[j])
        peak = max(peak, x_up[k + 1])
        if (not finite_r) and x_up[k + 1] < tol * peak and k >= 1:
            k_used = j
            break
    if not finite_r:
        if k_used == rr and rr > 0:
            raise ValueError(
                'the queue-length tail is still %g of its peak at the truncation level %d; '
                'with r = inf the patience law must make the chain ergodic (raise maxQueue if '
                'the model is genuinely that large)' % (x_up[rr] / peak, rr))
        x_up = x_up[:k_used + 1]
        dlt = dlt[:k_used + 1]
        delta = delta[:k_used]
        rr = k_used

    # The downward leg, eq. (7.5), over the states where not all servers are busy.
    x_down = np.zeros(s)
    xk = 1.0
    for k in range(s, 0, -1):
        xk = k * mu * xk / lambda_val
        x_down[k - 1] = xk

    x = np.concatenate([x_down, x_up])
    p = x / np.sum(x)
    prob_loss = p[-1] if finite_r else 0.0
    pa = p / (1.0 - prob_loss)          # eq. (7.8), the state seen by an ENTERING customer

    k_all = np.arange(s + rr + 1)
    q_all = np.maximum(0, k_all - s)
    mean_number = float(np.sum(k_all * p))
    var_number = float(np.sum(((k_all - mean_number) ** 2) * p))
    mean_queue = float(np.sum(q_all * p))
    var_queue = float(np.sum(((q_all - mean_queue) ** 2) * p))
    utilization = float(np.sum(np.minimum(k_all, s) * p) / s)

    prob_no_wait = float(np.sum(pa[:s]))        # eq. (7.9), states 0..s-1

    sigma = np.zeros(rr)
    m_sum = np.zeros(rr)
    v_sum = np.zeros(rr)
    ewa1 = np.zeros(rr)
    ewa2 = np.zeros(rr)
    kernels = []
    for k in range(1, rr + 1):
        phik, mk = _kernel(k, smu, dlt, delta)
        surv = np.concatenate([[1.0], np.cumprod(1.0 - phik)[:-1]])
        sigma[k - 1] = float(np.prod(1.0 - phik))
        m_sum[k - 1] = float(np.sum(mk))
        v_sum[k - 1] = float(np.sum(mk ** 2))
        # Eqs. (7.28)-(7.29): abandoning at the jth departure epoch costs the sum
        # of the first j interdeparture times, whose moments accumulate.
        cum_m = np.cumsum(mk)
        cum_v = np.cumsum(mk ** 2)
        w = surv * phik
        ewa1[k - 1] = float(np.sum(w * cum_m))
        ewa2[k - 1] = float(np.sum(w * (cum_v + cum_m ** 2)))
        kernels.append((1.0 / mk, phik, surv))

    w_arr = pa[s:s + rr]                # w_arr[k] = pa_{s+k}, the arrival joins position k+1
    prob_served = prob_no_wait + float(np.sum(w_arr * sigma))
    prob_abandon = 1.0 - prob_served
    ews1 = float(np.sum(w_arr * sigma * m_sum))                     # eq. (7.16)
    ews2 = float(np.sum(w_arr * sigma * (v_sum + m_sum ** 2)))      # eq. (7.17)
    ewa1_tot = float(np.sum(w_arr * ewa1))                          # eq. (7.26)
    ewa2_tot = float(np.sum(w_arr * ewa2))                          # eq. (7.27)

    mean_wait_served = _ratio(ews1, prob_served)
    mean_wait_abandon = _ratio(ewa1_tot, prob_abandon)
    result: Dict[str, Any] = {
        'queueLengthDist': p,
        'probLoss': float(prob_loss),
        'probNoWait': prob_no_wait,
        'probServed': prob_served,
        'probAbandon': prob_abandon,
        'meanNumber': mean_number,
        'varNumber': var_number,
        'meanQueueLength': mean_queue,
        'varQueueLength': var_queue,
        'utilization': utilization,
        'throughput': lambda_val * (1.0 - prob_loss) * prob_served,
        'abandonRate': lambda_val * (1.0 - prob_loss) * prob_abandon,
        'meanWaitServed': mean_wait_served,
        'varWaitServed': max(0.0, _ratio(ews2, prob_served) - mean_wait_served ** 2),
        'meanWaitAbandon': mean_wait_abandon,
        'varWaitAbandon': max(0.0, _ratio(ewa2_tot, prob_abandon) - mean_wait_abandon ** 2),
        'meanWait': ews1 + ewa1_tot,
        'secondMomentWait': ews2 + ewa2_tot,
        'abandonRates': delta,
        'totalAbandonRates': dlt,
        'numWaitingSpaces': rr,
        'isExponentialPatience': is_exponential,
        'patienceRate': theta,
    }

    if wPoints is not None and len(wPoints) > 0:
        from ..lti import laplace_invert

        t = np.atleast_1d(np.asarray(wPoints, dtype=float))
        fs = np.zeros(t.size)
        fa = np.zeros(t.size)
        for i in range(t.size):
            fs[i] = laplace_invert(
                lambda z: _transform(z, w_arr, sigma, kernels, True) / z,
                float(t[i]), invMethod, invN)
            fa[i] = laplace_invert(
                lambda z: _transform(z, w_arr, sigma, kernels, False) / z,
                float(t[i]), invMethod, invN)
        fs = np.clip(fs, 0.0, max(prob_served - prob_no_wait, 0.0))
        fa = np.clip(fa, 0.0, max(prob_abandon, 0.0))
        result['waitPoints'] = t
        result['cdfWaitServed'] = (prob_no_wait + fs) / max(prob_served, np.finfo(float).tiny)
        result['cdfWaitAbandon'] = fa / max(prob_abandon, np.finfo(float).tiny)
        result['cdfWait'] = prob_no_wait + fs + fa

    return result


def qsys_erlanga(lambda_val: float, mu: float, theta: float, s: int,
                 r: float = float('inf'), **kwargs) -> Dict[str, Any]:
    """
    Exact analysis of the Erlang A model M/M/s/r+M.

    Poisson arrivals at rate ``lambda_val``, exponential service of rate ``mu``
    at each of ``s`` servers and exponential patience of rate ``theta``. The
    number in system is the birth-and-death process with death rate
    ``min(k,s)*mu + (k-s)^+ * theta``, so every measure is exact: this is the
    case in which the state-dependent Markovian approximation of
    :func:`qsys_mgisrgi_whitt` reproduces the model rather than approximating it
    (eq. 7.12 of the reference). ``theta = 0`` recovers M/M/s/r, and then a
    finite ``r`` is required whenever ``lambda_val >= s*mu``.

    Args:
        lambda_val: arrival rate
        mu: service rate of one server
        theta: abandonment rate of a waiting customer
        s: number of servers
        r: extra waiting spaces, infinite by default
        **kwargs: passed through to :func:`qsys_mgisrgi_whitt`

    Returns:
        The dict returned by :func:`qsys_mgisrgi_whitt`.

    References:
        W. Whitt (2005). Engineering solution of a basic call-center model.
        Management Science 51(2), 221-235, Section 7 and eq. (7.12). The model
        itself is due to C. Palm (1937, 1957).
    """
    if theta <= 0 and not np.isfinite(r) and lambda_val >= s * mu:
        raise ValueError('without abandonment (theta = 0) and with an infinite waiting room the '
                         'queue is unstable at lambda >= s*mu; give a finite r or a positive theta')
    return qsys_mgisrgi_whitt(lambda_val, mu, s, r, theta, **kwargs)
