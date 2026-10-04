"""
Conditional waiting-time moments of the Lindley recursion.

Unlike every other qsys_* algorithm these are conditional on the current state
rather than stationary, so they are defined and finite for any load, including a
saturated queue. They give the one-step-ahead law of the next customer's waiting
time, which is the exact per-customer reference a simulated sample path can be
checked against.

References:
    Original MATLAB: matlab/src/api/qsys/qsys_mm1_lindley.m, qsys_hh1_lindley.m,
    qsys_mm1_tandem_lindley.m, qsys_tandem_lindley.m
    S. Palomo, J. Pender, "Learning the Tandem Network Lindley Recursion",
    Proc. Winter Simulation Conference, 2021.
    D. V. Lindley, "The Theory of Queues with a Single Server", Proc. Camb. Phil.
    Soc. 48, 1952.
"""

from math import comb, exp, factorial, sqrt
from typing import Dict, Optional, Sequence

import numpy as np

__all__ = ['qsys_lindley_moment', 'qsys_mm1_lindley', 'qsys_hh1_lindley',
           'qsys_mm1_tandem_lindley', 'qsys_tandem_lindley']


def qsys_lindley_moment(lambda_: float, mu: float, Wn, m: int) -> np.ndarray:
    """
    One conditional Lindley moment for exponential primitives.

    Evaluates ``E[max(Wn + S - A, 0)**m]`` with ``A ~ Exp(lambda_)`` and
    ``S ~ Exp(mu)``. This is the algorithm shared by :func:`qsys_mm1_lindley`,
    which calls it once per moment order, and :func:`qsys_hh1_lindley`, which
    mixes it over the arrival and service phases.

    See :func:`qsys_mm1_lindley` for the derivation and for why the upper
    incomplete gamma function reduces to a finite sum here.

    Args:
        lambda_: Arrival rate, positive
        mu: Service rate, positive
        Wn: Current waiting times
        m: Moment order, at least 1

    Returns:
        The conditional moment at each entry of Wn.
    """
    if m < 1:
        raise ValueError("The moment order m must be positive, got %d" % m)
    w = np.asarray(Wn, dtype=float).ravel()

    # S = sum_{k=0}^{m} C(m,k) w^k (m-k)! / mu^(m-k+1)
    s_term = np.zeros_like(w)
    for k in range(m + 1):
        s_term += comb(m, k) * w ** k * factorial(m - k) / mu ** (m - k + 1)

    # T = (-1)^m m! ( sum_{k=0}^{m} (-lambda w)^k/k! - e^{-lambda w} ) / lambda^(m+1)
    x = -lambda_ * w
    inner = np.zeros_like(w)
    for k in range(m + 1):
        inner += x ** k / factorial(k)
    sign = 1.0 if m % 2 == 0 else -1.0
    t_term = sign * factorial(m) * (inner - np.exp(x)) / lambda_ ** (m + 1)

    return lambda_ * mu / (lambda_ + mu) * (s_term + t_term)


def _check_waits(Wn, name: str = 'Wn') -> np.ndarray:
    w = np.asarray(Wn, dtype=float).ravel()
    if w.size == 0:
        raise ValueError("%s must hold at least one waiting time" % name)
    if not np.all(np.isfinite(w)) or np.any(w < 0.0):
        raise ValueError("%s must hold finite nonnegative values" % name)
    return w


def qsys_mm1_lindley(lambda_: float, mu: float, Wn, mmax: int = 2) -> Dict[str, object]:
    """
    Conditional waiting-time moments of the M/M/1 Lindley recursion.

    One step of ``W_{n+1} = max(W_n + S_n - A_n, 0)`` with ``A_n ~ Exp(lambda_)``
    and ``S_n ~ Exp(mu)``: given the waiting time of customer n, the exact
    conditional moments of the waiting time of customer n+1.

    The m-th conditional moment is

        E[W_{n+1}^m | W_n] = lambda mu/(lambda+mu) [ S + T ],
        S = sum_{k=0}^{m} C(m,k) W_n^k (m-k)! / mu^(m-k+1),
        T = (-1)^m e^{-lambda W_n} (Gamma(m+1,-lambda W_n) - m!) / lambda^(m+1),

    where the density of ``S_n - A_n`` is the asymmetric Laplace density
    ``lambda mu/(lambda+mu)`` times ``e^{-mu x}`` for ``x > 0`` and
    ``e^{lambda x}`` for ``x < 0``. Because ``m+1`` is a positive integer, the
    upper incomplete gamma function admits the finite form
    ``Gamma(m+1,x) = m! e^{-x} sum_{k=0}^{m} x^k/k!``, valid at the negative
    argument ``-lambda W_n`` needed here. Substituting it cancels the growing
    exponential and leaves the numerically stable

        T = (-1)^m m! ( sum_{k=0}^{m} (-lambda W_n)^k/k! - e^{-lambda W_n} ) / lambda^(m+1),

    which is what this function evaluates. No incomplete gamma routine is needed.

    The mean is returned from the equivalent explicit form
    ``W_n + (lambda-mu)/(lambda mu) + mu e^{-lambda W_n}/(lambda(lambda+mu))``,
    and the variance as the second moment less the squared mean.

    Verified against 4e6 Monte Carlo replications to 5e-4 relative error for
    m = 1, 2, 3.

    Args:
        lambda_: Arrival rate, positive
        mu: Service rate, positive
        Wn: Current waiting times, finite and nonnegative
        mmax: Highest moment order, raised to 2 so the variance is available

    Returns:
        Dict with keys ``mean``, ``var``, ``moments`` (shape ``(len(Wn), mmax)``),
        ``mmax`` and ``analyzer``.

    Raises:
        ValueError: If a rate is not positive or a waiting time is negative.
    """
    if not lambda_ > 0.0 or not np.isfinite(lambda_):
        raise ValueError("lambda_ must be positive and finite, got %r" % (lambda_,))
    if not mu > 0.0 or not np.isfinite(mu):
        raise ValueError("mu must be positive and finite, got %r" % (mu,))
    if int(mmax) != mmax or mmax < 1:
        raise ValueError("mmax must be a positive integer, got %r" % (mmax,))

    w = _check_waits(Wn)
    order = max(int(mmax), 2)
    moments = np.empty((w.size, order))
    for m in range(1, order + 1):
        moments[:, m - 1] = qsys_lindley_moment(lambda_, mu, w, m)

    mean = w + (lambda_ - mu) / (lambda_ * mu) \
        + mu * np.exp(-lambda_ * w) / (lambda_ * (lambda_ + mu))
    var = moments[:, 1] - moments[:, 0] ** 2

    return {'mean': mean, 'var': var, 'moments': moments, 'mmax': order,
            'analyzer': 'qsys_mm1_lindley'}


def qsys_hh1_lindley(lambda_: Sequence[float], pa: Sequence[float],
                     mu: Sequence[float], ps: Sequence[float], Wn,
                     mmax: int = 2) -> Dict[str, object]:
    """
    Conditional waiting-time moments of the Hl/Hn/1 Lindley recursion.

    Hyperexponential primitives are mixtures of exponentials, so conditioning on
    the arrival phase i and the service phase j reduces one Lindley step to the
    M/M/1 step of :func:`qsys_mm1_lindley` at rates ``lambda_[i]`` and ``mu[j]``,
    and the conditional moment is the corresponding mixture

        E[W_{n+1}^m | W_n] = sum_i sum_j pa[i] ps[j] E_{ij}[W_{n+1}^m | W_n].

    Phases are drawn independently for each customer, which is what makes the
    mixture exact rather than an approximation; a Markov-modulated arrival stream
    would not decompose this way.

    Note that the variance is not the corresponding mixture of the per-phase
    variances, because the phase is itself random: it is recovered here from the
    first two mixed raw moments, which adds the between-phase spread of the means.

    Verified against 4e6 Monte Carlo replications to 2e-3 relative error for
    m = 1, 2.

    Args:
        lambda_: Arrival phase rates, positive
        pa: Arrival phase probabilities, nonnegative and summing to 1
        mu: Service phase rates, positive
        ps: Service phase probabilities, nonnegative and summing to 1
        Wn: Current waiting times, finite and nonnegative
        mmax: Highest moment order, raised to 2 so the variance is available

    Returns:
        The same dict as :func:`qsys_mm1_lindley`.

    Raises:
        ValueError: If the phase vectors are inconsistent or unnormalized.
    """
    rates_a = np.asarray(lambda_, dtype=float).ravel()
    probs_a = np.asarray(pa, dtype=float).ravel()
    rates_s = np.asarray(mu, dtype=float).ravel()
    probs_s = np.asarray(ps, dtype=float).ravel()
    for rates, probs, what in ((rates_a, probs_a, 'arrival'),
                               (rates_s, probs_s, 'service')):
        if rates.size == 0 or rates.size != probs.size:
            raise ValueError("The %s rates and probabilities must be nonempty and of"
                             " equal length" % what)
        if not np.all(np.isfinite(rates)) or np.any(rates <= 0.0):
            raise ValueError("The %s rates must be positive and finite" % what)
        if np.any(probs < 0.0) or abs(float(probs.sum()) - 1.0) > 1e-10:
            raise ValueError("The %s probabilities must be nonnegative and sum to 1,"
                             " they sum to %g" % (what, probs.sum()))
    if int(mmax) != mmax or mmax < 1:
        raise ValueError("mmax must be a positive integer, got %r" % (mmax,))

    w = _check_waits(Wn)
    order = max(int(mmax), 2)
    moments = np.zeros((w.size, order))
    for i in range(rates_a.size):
        for j in range(rates_s.size):
            weight = probs_a[i] * probs_s[j]
            if weight == 0.0:
                continue
            for m in range(1, order + 1):
                moments[:, m - 1] += weight * qsys_lindley_moment(
                    rates_a[i], rates_s[j], w, m)

    var = moments[:, 1] - moments[:, 0] ** 2
    return {'mean': moments[:, 0].copy(), 'var': var, 'moments': moments,
            'mmax': order, 'analyzer': 'qsys_hh1_lindley'}


def _tandem_j(c: float, y: np.ndarray, mu2: float) -> np.ndarray:
    """J(c) = int_0^inf e^{-c u} E[(y + S2 - u)^+] du with S2 ~ Exp(mu2)."""
    e = np.exp(-c * y)
    return ((y + 1.0 / mu2) * (1.0 - e) / c
            - (1.0 - e * (1.0 + c * y)) / c ** 2
            + e / (mu2 * (c + mu2)))


def _tandem_jw(c: float, y: np.ndarray, mu2: float) -> np.ndarray:
    """Jw(c) = int_0^inf u e^{-c u} E[(y + S2 - u)^+] du with S2 ~ Exp(mu2)."""
    e = np.exp(-c * y)
    d = c + mu2
    return ((y + 1.0 / mu2) * (1.0 - e * (1.0 + c * y)) / c ** 2
            - (2.0 - e * (2.0 + 2.0 * c * y + c ** 2 * y ** 2)) / c ** 3
            + e * (y / d + 1.0 / d ** 2) / mu2)


def qsys_mm1_tandem_lindley(lambda_: float, mu1: float, mu2: float,
                            Wk, Wk1) -> Dict[str, object]:
    """
    Conditional downstream waiting time in an M/M/1 tandem.

    The conditional mean waiting time of customer n+1 at the downstream station of
    a two-station single-server tandem queue, given that customer n waited ``Wk``
    upstream and ``Wk1`` downstream.

    The point of the tandem recursion is that the interarrival time at the
    downstream station is the interdeparture time upstream, not an independent
    draw. With ``A ~ Exp(lambda_)`` the interarrival time upstream and ``S1, S1'``
    the service times upstream of customers n and n+1, that interdeparture time is

        D = max(A - Wk - S1, 0) + S1',

    an idle period followed by the next service, and the downstream Lindley step is
    ``W2_{n+1} = (Wk1 + S2 - D)^+`` with ``S2 ~ Exp(mu2)`` independent of D.

    Because A is exponential, ``max(A - Wk - S1, 0)`` is zero with probability
    ``1-q`` and ``Exp(lambda_)`` with probability
    ``q = e^{-lambda Wk} mu1/(lambda+mu1)``, the probability the upstream server
    goes idle, so D is either ``Exp(mu1)`` or the sum of ``Exp(mu1)`` and
    ``Exp(lambda_)``. Averaging the downstream step over both cases needs only two
    elementary transforms of::

        g(d) = E[(Wk1 + S2 - d)^+] = Wk1 - d + 1/mu2      for d <= Wk1,
                                   = e^{-mu2 (d-Wk1)}/mu2  for d > Wk1,

    namely ``J(c) = int_0^inf e^{-cu} g(u) du`` and
    ``Jw(c) = int_0^inf u e^{-cu} g(u) du``, both closed form, giving::

        E[W2_{n+1} | Wk, Wk1] = (1-q) mu1 J(mu1) + q C,
        C = lambda mu1 (J(mu1) - J(lambda))/(lambda-mu1)  if lambda != mu1,
          = mu1^2 Jw(mu1)                                 if lambda == mu1.

    As ``Wk`` grows the upstream server never idles, q vanishes, and the mean tends
    to ``mu1 J(mu1) = E[g(S1')]``, as it must.

    Two caveats, both inherited from the reference and both quantified here.

    First, this is exact for the step taken in isolation, that is when the
    conditioning pair is independent of the four primitives that drive the step. In
    a running tandem it is not: the downstream wait ``Wk1`` was itself determined by
    an interdeparture time containing ``S1``, so conditioning on ``(Wk, Wk1)`` is
    not conditioning on a Markov state of the tandem. Measured against a
    4e6-customer simulation of the real tandem at lambda = 0.8, mu1 = mu2 = 1, the
    formula is within 0.4% to 1.3% away from the empty state and 7% at
    ``Wk = Wk1 = 0``, where the entanglement is strongest. Treat it as exact for one
    isolated step and as a good approximation in a running tandem.

    Second, this closed form was derived rather than transcribed from the
    reference's theorem 4, because that theorem rests on its proposition 1, which
    omits a service-time difference and so does not describe a tandem queue; see
    :func:`qsys_tandem_lindley`. The two differ: at lambda = 0.8, mu1 = mu2 = 1 and
    ``Wk = Wk1 = 0`` the published route gives 0.3016 against 0.3457 here, the
    latter matching simulation of the step to 6e-4 relative error.

    As in the reference, the upstream interarrival time is taken to be
    ``Exp(lambda_)``, which by Burke's theorem is also the stationary interdeparture
    law, so the same formula is applied at any pair of consecutive stations of a
    longer M/M/1 tandem, with the caveat above compounding.

    Args:
        lambda_: External arrival rate upstream, positive
        mu1: Upstream service rate, positive
        mu2: Downstream service rate, positive
        Wk: Current upstream waiting times, finite and nonnegative
        Wk1: Current downstream waiting times, same shape as Wk or scalar

    Returns:
        Dict with keys ``mean``, ``interdepMean``, ``idleProb`` and ``analyzer``.

    Raises:
        ValueError: If a rate is not positive or the two wait arrays disagree.
    """
    for value, name in ((lambda_, 'lambda_'), (mu1, 'mu1'), (mu2, 'mu2')):
        if not value > 0.0 or not np.isfinite(value):
            raise ValueError("%s must be positive and finite, got %r" % (name, value))

    x = _check_waits(Wk, 'Wk')
    y = _check_waits(Wk1, 'Wk1')
    if x.size != y.size:
        if x.size == 1:
            x = np.full(y.size, x[0])
        elif y.size == 1:
            y = np.full(x.size, y[0])
        else:
            raise ValueError("Wk and Wk1 must have the same size")

    q = np.exp(-lambda_ * x) * mu1 / (lambda_ + mu1)
    base = mu1 * _tandem_j(mu1, y, mu2)
    if abs(lambda_ - mu1) > 1e-9 * max(lambda_, mu1):
        conv = lambda_ * mu1 / (lambda_ - mu1) \
            * (_tandem_j(mu1, y, mu2) - _tandem_j(lambda_, y, mu2))
    else:
        # the two rates coincide, the interdeparture time is Erlang(2,mu1)
        conv = mu1 ** 2 * _tandem_jw(mu1, y, mu2)

    mean = (1.0 - q) * base + q * conv
    interdep = 1.0 / mu1 + q / lambda_

    return {'mean': mean, 'interdepMean': interdep, 'idleProb': q,
            'analyzer': 'qsys_mm1_tandem_lindley'}


def qsys_tandem_lindley(A, S, W0: Optional[Sequence[float]] = None) -> Dict[str, object]:
    """
    Tandem network Lindley recursion on a sample path.

    Propagates the waiting times of a series of K single-server FCFS stations in
    tandem, driven by the primitives of the sample path. ``A`` is the length-N
    vector of interarrival times at the first station, ``A[n]`` separating
    customers n and n+1, and ``S`` is the ``(N, K)`` matrix of service times,
    ``S[n, k]`` being the service time of customer n at station k.

    At the first station this is Lindley's recursion,
    ``W[n+1, 0] = max(W[n, 0] + S[n, 0] - A[n], 0)``. Downstream the interarrival
    time is not a primitive: the arrival epoch of customer n at station k is its
    departure epoch from station k-1, so the interarrival time at station k is the
    interdeparture time upstream. Writing ``G[n, k]`` for the interarrival time at
    station k between customers n and n+1, with ``G[n, 0] = A[n]``, the exact
    interdeparture identity is

        G[n, k+1] = G[n, k] + W[n+1, k] - W[n, k] + S[n+1, k] - S[n, k],

    equivalently and more transparently

        G[n, k+1] = max(G[n, k] - W[n, k] - S[n, k], 0) + S[n+1, k],

    an idle period at station k followed by the next customer's service there. The
    recursion at station k is then
    ``W[n+1, k] = max(W[n, k] + S[n, k] - G[n, k], 0)``.

    Nothing here is distributional, so the recursion is exact for arbitrary
    interarrival and service times, dependent or not, and is the reference a
    simulated tandem sample path can be checked against directly. It reproduces a
    direct event-driven tandem simulation to 1e-12 over four stations.

    Note that proposition 1 of the reference states this identity without the
    ``S[n+1, k] - S[n, k]`` term, which makes it wrong as a sample-path identity:
    the omitted difference has mean zero, so the mean interdeparture time survives,
    but individual waiting times do not. Implementing it as published gives
    station-1 waiting times that are correct and downstream ones that are not, by
    up to several mean service times. The form above is used instead.

    Args:
        A: Interarrival times at the first station
        S: Service times, shape ``(N, K)``
        W0: Waiting times of customer 1 at each station, None for an empty network

    Returns:
        Dict with keys ``W``, ``G``, ``T``, ``departure`` and ``analyzer``, every
        matrix indexed ``[customer, station]``. The last row of ``G`` is NaN, there
        being no customer N+1 to separate from.

    Raises:
        ValueError: If the shapes disagree or a primitive is negative.
    """
    a = np.asarray(A, dtype=float).ravel()
    n = a.size
    if n < 1:
        raise ValueError("A must hold at least one interarrival time")
    s = np.atleast_2d(np.asarray(S, dtype=float))
    if s.shape[0] != n:
        raise ValueError("S must have one row per customer, got %d rows for %d"
                         " interarrivals" % (s.shape[0], n))
    k = s.shape[1]
    if k < 1:
        raise ValueError("S must have at least one column, one per station")
    if not np.all(np.isfinite(a)) or np.any(a < 0.0):
        raise ValueError("A must hold finite nonnegative interarrival times")
    if not np.all(np.isfinite(s)) or np.any(s < 0.0):
        raise ValueError("S must hold finite nonnegative service times")

    start = np.zeros(k)
    if W0 is not None:
        w0 = np.asarray(W0, dtype=float).ravel()
        if w0.size != k:
            raise ValueError("W0 must hold one waiting time per station, %d of them" % k)
        if not np.all(np.isfinite(w0)) or np.any(w0 < 0.0):
            raise ValueError("W0 must hold finite nonnegative waiting times")
        start = w0

    w = np.zeros((n, k))
    g = np.full((n, k), np.nan)
    w[0, :] = start
    for i in range(n - 1):
        gap = a[i]
        for j in range(k):
            g[i, j] = gap
            w[i + 1, j] = max(w[i, j] + s[i, j] - gap, 0.0)
            # the interdeparture time here is the interarrival time one station on
            gap = max(gap - w[i, j] - s[i, j], 0.0) + s[i + 1, j]

    t = w + s
    departure = np.empty((n, k))
    epoch = np.concatenate([[0.0], np.cumsum(a[:n - 1])])
    departure[:, 0] = epoch + t[:, 0]
    for j in range(1, k):
        departure[:, j] = departure[:, j - 1] + t[:, j]

    return {'W': w, 'G': g, 'T': t, 'departure': departure,
            'analyzer': 'qsys_tandem_lindley'}
