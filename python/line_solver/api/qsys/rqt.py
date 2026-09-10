"""
Robust Queueing Theory (RQT) worst-case analysis of a G/G/k FCFS queue.

The arrival and service processes are not described by distributions but by
polyhedral uncertainty sets whose shape follows the (generalized) central limit
theorem, so performance analysis becomes a worst-case optimization.

Reference: C. Bandi, D. Bertsimas, N. Youssef (2015), "Robust Queueing Theory",
Operations Research 63(3), 676-700.

MATLAB: matlab/src/api/qsys/qsys_gigk_rqt.m, qsys_gig1_rqt.m,
        qsys_gigk_rqt_gamma.m
"""

from typing import Tuple

import numpy as np

__all__ = ['qsys_gigk_rqt', 'qsys_gig1_rqt', 'qsys_gigk_rqt_gamma']

# Table 1 adaptation regimes: (theta0, theta1, theta2)
_RQT_THETA = {
    'pareto': (-0.05, 1.09, 1.11),
    'normal': (-0.02, 1.03, 1.04),
    'independent': (-0.06, 1.07, 1.07),
    'default': (-0.06, 1.07, 1.07),
}


def qsys_gigk_rqt(lambda_val: float, mu: float, Gamma_a: float, Gamma_s: float,
                  k: int = 1, alpha_a: float = 2.0,
                  alpha_s: float = 2.0) -> Tuple[float, float, float]:
    """
    Robust Queueing Theory worst-case system time of a G/G/k FCFS queue.

    The uncertainty sets are

        U^a = {T : (sum_{i=k+1}^n T_i - (n-k)/lambda)/(n-k)^(1/alpha_a) >= -Gamma_a}
        U^s = {X : (sum_{i=k}^n X_i - (n-k+1)/mu)/(n-k+1)^(1/alpha_s) <= Gamma_s}

    with alpha=2 the finite-variance regime and alpha in (1,2) the heavy-tailed
    one. The returned W is the closed-form bound of Theorem 3 (Theorem 8 when the
    two tail coefficients differ, with alphabar = min(alpha_a,alpha_s)),

        W <= (ab-1)/ab^(ab/(ab-1)) lambda^(1/(ab-1))
             (Gamma_a+Gamma_s/k^(1/ab))^(ab/(ab-1)) / (1-rho)^(1/(ab-1)) + k/lambda,

    which for k=1 reduces to Theorem 2 and, at alphabar=2, to the Kingman-like
    form (lambda/4)(Gamma_a+Gamma_s)^2/(1-rho) + 1/lambda. Sworst is the exact
    worst case over the uncertainty sets, eq. (45), the supremum over the integer
    x >= 1 of

        x/mu + Gamma_s x^(1/alpha_s) - k(x-1)/lambda + Gamma_a (k(x-1))^(1/alpha_a).

    The arrival deviation ADDS to the worst case, since the adversary shortens
    the interarrival times: the sign printed in eq. (12) is easily misread as a
    subtraction of the whole arrival bracket.

    W is a SYSTEM time (waiting plus service), and its additive term is k/lambda
    rather than the mean service time 1/mu.

    Args:
        lambda_val: Arrival rate
        mu: Service rate of each server
        Gamma_a: Variability parameter of the arrival uncertainty set
        Gamma_s: Variability parameter of the service uncertainty set
        k: Number of servers
        alpha_a: Arrival tail coefficient in (1,2]
        alpha_s: Service tail coefficient in (1,2]

    Returns:
        Tuple of (W, rhohat, Sworst)
    """
    if alpha_a <= 1 or alpha_a > 2 or alpha_s <= 1 or alpha_s > 2:
        raise RuntimeError("RQT tail coefficients must lie in (1,2].")

    rho = lambda_val / (k * mu)
    if lambda_val <= 0:
        return 1.0 / mu, 0.0, 1.0 / mu
    if rho >= 1:
        return np.inf, 1.0, np.inf

    # Theorem 8 collapses to Theorem 3 when the two tails agree
    ab = min(alpha_a, alpha_s)
    beta = Gamma_a + Gamma_s / k ** (1.0 / ab)
    if beta <= 0:
        # a nonpositive effective variability leaves only the deterministic term
        W = k / lambda_val
    else:
        W = ((ab - 1) / ab ** (ab / (ab - 1)) * lambda_val ** (1.0 / (ab - 1))
             * beta ** (ab / (ab - 1)) / (1 - rho) ** (1.0 / (ab - 1)) + k / lambda_val)
    rhohat = W * lambda_val / (1 + W * lambda_val)

    def obj(x):
        y = max(x - 1.0, 0.0)
        return (x / mu + Gamma_s * x ** (1.0 / alpha_s) - k * y / lambda_val
                + Gamma_a * (k * y) ** (1.0 / alpha_a))

    # the continuous maximizer of the bounding problem, eq. (16), sizes the scan
    xstar = (lambda_val * beta / (ab * (1 - rho))) ** (ab / (ab - 1)) if beta > 0 else 1.0
    xhi = max(4.0, np.ceil(4.0 * xstar))
    xs = np.unique(np.round(np.logspace(0, np.log10(xhi), 400)))
    gs = np.array([obj(float(x)) for x in xs])
    imax = int(np.argmax(gs))
    Sworst = float(gs[imax])
    # refine on the continuous relaxation, then round back onto the integer lattice
    lo = float(xs[max(0, imax - 1)])
    hi = float(xs[min(len(xs) - 1, imax + 1)])
    if hi > lo:
        xc = _golden_max(obj, lo, hi, 1e-8)
        for x in (np.floor(xc), np.ceil(xc)):
            if x >= 1:
                Sworst = max(Sworst, obj(float(x)))
    return W, rhohat, Sworst


def qsys_gig1_rqt(lambda_val: float, mu: float, Gamma_a: float, Gamma_s: float,
                  alpha_a: float = 2.0, alpha_s: float = 2.0) -> Tuple[float, float, float]:
    """
    Robust Queueing Theory worst-case system time of a G/G/1 FCFS queue, the
    single-server case of qsys_gigk_rqt (Theorem 2 and eq. 12).

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        Gamma_a: Variability parameter of the arrival uncertainty set
        Gamma_s: Variability parameter of the service uncertainty set
        alpha_a: Arrival tail coefficient in (1,2]
        alpha_s: Service tail coefficient in (1,2]

    Returns:
        Tuple of (W, rhohat, Sworst)
    """
    return qsys_gigk_rqt(lambda_val, mu, Gamma_a, Gamma_s, 1, alpha_a, alpha_s)


def qsys_gigk_rqt_gamma(rho: float, mu: float, Gamma_a: float, sigma_s: float,
                        k: int = 1, alpha_a: float = 2.0,
                        regime: str = 'independent') -> float:
    """
    Service variability parameter of the RQT framework, from the first two
    moments, by the adaptation of Section 7.1::

        Gamma_s = (2 (theta0 + theta1 sigma_s^2/k + theta2 Gamma_a^2 rho^2 k))^((a-1)/a)
                  - Gamma_a k^((a-1)/a)

    where (theta0,theta1,theta2) are regressed so that the worst-case system time
    of Theorem 3 approximates the MEAN system time of the corresponding
    stochastic queue. The arrival side needs no adaptation: Gamma_a = sigma_a for
    an external renewal stream. Since the last term cancels Gamma_a at alpha=2,
    the adaptation acts on the sum Gamma_a + Gamma_s/k^(1/alpha) that Theorem 3
    reads.

    THE FACTOR 2 IS NOT IN THE PRINTED FORMULA and is restored here. Section 7.1
    states that the form is motivated by Kingman's bound, which the alpha=2 bound
    of Theorem 3 reproduces when (Gamma_a+Gamma_s)^2 = 2(sigma_a^2+sigma_s^2);
    the published thetas are all near unity, i.e. corrections to that bound
    rather than a substitute for its factor 2. Dropping the factor puts M/M/1
    about 40% BELOW its exact mean system time at rho=0.9, contradicting the
    errors of at most 9.5% that Tables 2-3 report; restoring it gives +4.7%.

    CAUTION: the form is not dimensionally homogeneous, since theta0 is an
    additive constant on a scale of variances, so it is only valid in the time
    unit the regression was run in. It is evaluated here in units of the mean
    service time, 1/mu = 1, and converted back.

    Args:
        rho: Traffic intensity lambda/(k*mu)
        mu: Service rate of each server, which sets the time unit
        Gamma_a: Variability parameter of the arrival uncertainty set
        sigma_s: Standard deviation of the service time
        k: Number of servers
        alpha_a: Effective arrival tail coefficient in (1,2]
        regime: Adaptation regime of Table 1, 'independent', 'normal' or 'pareto'

    Returns:
        The service variability parameter Gamma_s
    """
    key = (regime or 'independent').lower()
    if key not in _RQT_THETA:
        raise RuntimeError("Unknown RQT adaptation regime: %s" % regime)
    t0, t1, t2 = _RQT_THETA[key]

    # evaluate in units of the mean service time, then convert back
    ga = Gamma_a * mu
    ss = sigma_s * mu
    e = (alpha_a - 1) / alpha_a
    b = 2 * (t0 + t1 * ss ** 2 / k + t2 * ga ** 2 * rho ** 2 * k)
    gs = max(b, 0.0) ** e - ga * k ** e
    return gs / mu


def _golden_max(f, a: float, b: float, tol: float) -> float:
    """Golden-section maximization of a unimodal (in practice) function."""
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = f(c)
    fd = f(d)
    for _ in range(500):
        if abs(b - a) <= tol:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = f(d)
    return (a + b) / 2.0
