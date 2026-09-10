"""
Tail bounds for a tandem of two single-server FCFS stations.

Polynomial-exponential upper bounds on the end-to-end waiting and sojourn time
tails of a GI/Hn/1 -> ./Hn/1 tandem. The bound mixes two exponentials with a
polynomial of degree one, which is what lets it follow the concave bend of the
tail on a linear-log scale where a purely exponential bound cannot, and it is
exact in the M/M/1 -> ./M/1 case.

References:
    Original MATLAB: matlab/src/api/qsys/qsys_tandem_ub_ciucu.m
    F. Ciucu, S. Mehri, "On the Distribution of Sojourn Times in Tandem Queues",
    Proc. ACM Meas. Anal. Comput. Syst. 9(2), Article 27, 2025 (ACM SIGMETRICS
    2025). Registered in .citations() as 'tandemub'.
"""

from typing import Callable, Dict, Optional, Sequence

import numpy as np

__all__ = ['qsys_tandem_ub_ciucu']


def qsys_tandem_ub_ciucu(x, lst: Callable[[float], float], p, mu,
                         dlst: Optional[Callable[[float], float]] = None) -> Dict:
    """
    Tail bounds for a GI/Hn/1 -> ./Hn/1 tandem of two FCFS single servers.

    Both stations serve the same hyperexponential law ``Y, Z ~ sum_i p_i
    Exp(mu_i)``, a scalar ``p = 1`` giving exponential service, and the arrivals
    are renewal with a light-tailed interarrival time supplied through its
    Laplace-Stieltjes transform ``E[e^{-s X}]``.

    With ``theta`` the positive root of ``E[e^{theta (Y-X)}] = 1`` and ``alpha =
    E[X e^{-theta X}]``, the test function

        gamma(u,v) = 1{0<=u<=v} [1 - A e^{-theta u} - (B + C u + D v) e^{-theta v}]

    satisfies the integral inequality of Theorem 1(b) of the reference once the
    five sufficient conditions of its Lemma 4 fix A, B, C and D,

        A = 1,  C = theta sum_i p_i/(mu_i-theta) / sum_i p_i mu_i/(mu_i-theta)^2,
        D = max(-C E[U e^{theta V}]/E[V e^{theta V}], 0),  U = Y-X, V = Z-X,
        B = C (1/mu_1 - alpha E[e^{theta Z}])   if D = 0,
          = (C+D)/(mu_1-theta) - theta/mu_1     if D > 0,

    with ``mu_1`` the smallest service rate. Corollary 2 then turns gamma into

        P(S > x) <= sum_i p_i { e^{-mu_i x}
                      + mu_i/(mu_i-theta) (A+B) (e^{-theta x} - e^{-mu_i x})
                      + mu_i/(mu_i-theta)^2 (C+D) (((mu_i-theta)x-1) e^{-theta x}
                                                   + e^{-mu_i x}) }

    and the corresponding closed form for W when the service is exponential.
    ``E[V e^{theta V}]`` is positive at any stable load, so D is always well
    defined: ``h(s) = E[e^{s(Z-X)}]`` is convex with ``h(0) = h(theta) = 1``,
    hence ``h'(theta) > 0``.

    In the M/M/1 -> ./M/1 case the five inequalities hold as equalities, so gamma
    is the exact joint distribution and both bounds are exact, ``P(S > x) =
    (1 + theta x) e^{-theta x}``. Away from it the bound stays sharp: against an
    exact CTMC reference for the Erlang(2)/M/1 -> ./M/1 tandem it is within 2% at
    P(S>x) = 1e-2 and within 0.6% at 5e-10, with the correct asymptotic slope
    ``theta^2/(mu(1-alpha mu))``. Accuracy degrades with service variability, to
    about a factor of two at CV(Y) = 2.

    Args:
        x: Thresholds at which the tails are bounded, nonnegative
        lst: Interarrival transform, ``s -> E[e^{-s X}]`` for s >= 0
        p: Service phase probabilities, nonnegative and summing to one
        mu: Service phase rates, positive
        dlst: ``s -> E[X e^{-s X}]``, minus the derivative of lst; None to obtain
            it by a Richardson-extrapolated central difference, which costs four
            extra transform evaluations and loses roughly four digits

    Returns:
        Dict with keys 'S' (bound on P(S>x), capped at one), 'W' (bound on
        P(W>x), NaN unless the service is exponential), 'theta', 'alpha', and the
        coefficients 'A', 'B', 'C', 'D'.

    Example:
        >>> from math import exp
        >>> r = qsys_tandem_ub_ciucu([5, 10], lambda s: exp(-s * 4 / 3), 1.0, 1.0,
        ...                          lambda s: (4 / 3) * exp(-s * 4 / 3))
        >>> float(round(r['S'][0], 6))
        0.493699
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    p = np.atleast_1d(np.asarray(p, dtype=float))
    mu = np.atleast_1d(np.asarray(mu, dtype=float))
    if p.size != mu.size:
        raise ValueError("p and mu must have the same number of phases")
    if np.any(x < 0):
        raise ValueError("The thresholds x must be nonnegative")
    if np.any(p < 0) or abs(p.sum() - 1.0) > 1e-10:
        raise ValueError("The phase probabilities p must be nonnegative and sum to one")
    if np.any(mu <= 0):
        raise ValueError("The service rates mu must be positive")

    def mgf_y(t):
        return float(np.sum(p * mu / (mu - t)))          # E[e^{t Y}], t < min(mu)

    def residual(t):
        return mgf_y(t) * float(lst(t)) - 1.0

    mu1 = float(mu.min())
    # Stability: E[X] > E[Y] is what makes E[e^{t(Y-X)}] - 1 cross zero on (0,mu1).
    hi = mu1 * (1.0 - 1e-12)
    if residual(hi) <= 0.0:
        raise ValueError("No positive root of E[e^{theta(Y-X)}]=1 below min(mu): the tandem "
                         "is unstable or the service is not the lighter tail")
    lo = mu1 * 1e-12
    while residual(lo) >= 0.0 and lo > mu1 * 1e-16:
        lo = lo / 10.0                                   # walk below the root at zero
    if residual(lo) >= 0.0:
        # E[e^{t(Y-X)}]-1 is convex and vanishes at t=0, so it stays positive on the
        # whole of (0,mu1) exactly when its slope E[Y]-E[X] there is nonnegative.
        raise ValueError("The tandem is unstable, E[X] <= E[Y]: theta = 0 is the only "
                         "root of E[e^{theta(Y-X)}]=1")
    theta = _bisect(residual, lo, hi)
    alpha = float(dlst(theta)) if dlst is not None else _numerical_dlst(lst, theta)

    eexp_z = mgf_y(theta)                                # E[e^{theta Z}]
    ez_exp = float(np.sum(p * mu / (mu - theta) ** 2))   # E[Z e^{theta Z}]
    e_y = float(np.sum(p / mu))
    A = 1.0
    C = theta * float(np.sum(p / (mu - theta))) / ez_exp
    eu_ev = e_y - alpha * eexp_z                         # E[U e^{theta V}]
    ev_ev = ez_exp / eexp_z - alpha * eexp_z             # E[V e^{theta V}] > 0
    D = -C * eu_ev / ev_ev
    if not D > 0.0:
        D = 0.0
    if D > 0.0:
        B = (C + D) / (mu1 - theta) - theta / mu1
    else:
        B = C * (1.0 / mu1 - alpha * eexp_z)

    S = np.zeros_like(x)
    for i in range(p.size):
        m = mu[i]
        S = S + p[i] * (np.exp(-m * x)
                        + m / (m - theta) * (A + B) * (np.exp(-theta * x) - np.exp(-m * x))
                        + m / (m - theta) ** 2 * (C + D)
                        * (((m - theta) * x - 1.0) * np.exp(-theta * x) + np.exp(-m * x)))
    S = np.minimum(S, 1.0)

    if p.size == 1:
        beta = float(lst(mu1))                           # E[e^{-mu X}]
        if D == 0.0:
            W = ((1.0 - 2.0 * theta ** 2 / (mu1 * (mu1 + theta))
                  + theta * (mu1 - theta) / (mu1 + theta) * x) * np.exp(-theta * x)
                 + beta * (theta * mu1 * alpha / (2.0 * (mu1 - theta))
                           - theta / (2.0 * mu1)) * np.exp(-mu1 * x))
        else:
            W = ((1.0 - 2.0 * theta / mu1
                  + 2.0 * theta ** 2 * (2.0 - alpha * mu1)
                  / ((mu1 + theta) ** 2 * (1.0 - alpha * mu1))
                  + theta ** 2 * (mu1 - theta)
                  / (mu1 * (mu1 + theta) * (1.0 - alpha * mu1)) * x) * np.exp(-theta * x))
        W = np.minimum(W, 1.0)
    else:
        W = np.full_like(x, np.nan)                      # the W form is Exp-service only

    return {'S': S, 'W': W, 'theta': theta, 'alpha': alpha,
            'A': A, 'B': B, 'C': C, 'D': D, 'analyzer': 'qsys_tandem_ub_ciucu'}


def _bisect(f, lo: float, hi: float) -> float:
    """Bisection on a bracket with a sign change; deterministic and dependency free."""
    a, b = lo, hi
    for _ in range(200):
        if (b - a) <= 1e-14 * max(1.0, b):
            break
        m = 0.5 * (a + b)
        if f(m) > 0.0:
            b = m
        else:
            a = m
    return 0.5 * (a + b)


def _numerical_dlst(lst, s: float) -> float:
    """Richardson-extrapolated central difference of -lst at s, i.e. E[X e^{-s X}]."""
    h = 1e-3 * (1.0 + s)
    if h > s:
        h = s / 2.0
    d1 = (float(lst(s - h)) - float(lst(s + h))) / (2.0 * h)
    d2 = (float(lst(s - h / 2.0)) - float(lst(s + h / 2.0))) / h
    return (4.0 * d2 - d1) / 3.0
