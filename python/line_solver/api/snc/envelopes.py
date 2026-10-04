"""Stochastic network calculus: MGF arrival and service envelopes.

Native port of matlab/src/api/snc/snc_env_*.m, snc_srv_*.m and the three
min-plus operations, cross-checked against jline.api.snc.

An envelope is the pair ``(sigma(theta), rho(theta))`` in the exponential form

    E[exp(theta*A(s,t))] <= exp(theta*(rho*(t-s) + sigma)),   theta > 0,

and its service counterpart with the sign of theta reversed. Every function
here returns that pair as a tuple with ``theta`` LAST in the signature, so a
callable envelope is ``lambda theta: snc_env_poisson(lam, theta)``; the bound
functions in :mod:`line_solver.api.snc.bounds` take arrival and service as such
callables.

Time is slotted with unit slot length, which is what makes the geometric sum
over the start of the backlogged period converge to
``1/(1-exp(-theta*(rhoS-rhoA)))``; the continuous-time formulation would give
``1/(theta*(rhoS-rhoA))`` instead.

INFEASIBILITY IS SIGNALLED BY ``inf``, never by an exception: the compound
Poisson envelope diverges at ``theta >= mu`` and a composition can go unstable,
and the theta search discards those points.

Reference: M. Fidler, A. Rizk, "A Guide to the Stochastic Network Calculus",
IEEE Communications Surveys and Tutorials 17(1), 92-105, 2015.
"""

import math

import numpy as np

__all__ = [
    'snc_env_poisson', 'snc_env_cpoisson', 'snc_env_tokenbucket', 'snc_env_map',
    'snc_srv_rate', 'snc_srv_exp',
    'snc_leftover', 'snc_conv', 'snc_output',
]


def snc_env_poisson(lam, theta):
    """MGF arrival envelope of a Poisson flow with unit-size jobs.

    ``log E[exp(theta*A(0,t))] = lam*t*(exp(theta)-1)`` exactly, so the envelope
    is tight with a zero burst term.

    :param lam: arrival rate, jobs per slot
    :param theta: Chernoff parameter, theta > 0
    :return: ``(sigma, rho)``, with sigma = 0
    """
    if lam < 0:
        raise ValueError("snc_env_poisson: lam must be nonnegative, got %g." % lam)
    if theta <= 0:
        raise ValueError("snc_env_poisson: theta must be positive, got %g." % theta)
    return 0.0, lam * (math.exp(theta) - 1.0) / theta


def snc_env_cpoisson(lam, mu, theta):
    """MGF arrival envelope of a compound Poisson flow with Exp(mu) job sizes.

    Jobs arrive Poisson at rate ``lam`` and each carries an Exp(mu) amount of
    work, giving ``rho(theta) = lam/(mu-theta)`` for ``0 < theta < mu`` and a
    zero burst. Fed to a constant-rate server of rate mu this is the network
    calculus model of the M/M/1 queue in units of WORK.

    :param lam: job arrival rate, jobs per slot
    :param mu: rate of the Exp job size, so the mean work per job is 1/mu
    :param theta: Chernoff parameter, theta > 0
    :return: ``(sigma, rho)``, with rho = inf when theta >= mu
    """
    if lam < 0 or mu <= 0:
        raise ValueError("snc_env_cpoisson: lam must be nonnegative and mu positive, "
                         "got %g, %g." % (lam, mu))
    if theta <= 0:
        raise ValueError("snc_env_cpoisson: theta must be positive, got %g." % theta)
    if theta >= mu:
        return 0.0, float('inf')  # the job-size MGF diverges, no envelope here
    return 0.0, lam / (mu - theta)


def snc_env_tokenbucket(b, r, theta=None):
    """Deterministic token-bucket arrival envelope.

    ``A(s,t) <= b + r*(t-s)`` with probability one, so the envelope is constant
    in theta. ``theta`` is accepted and ignored, for signature compatibility
    with the stochastic envelopes.

    :param b: bucket depth, units of work
    :param r: token rate, work per slot
    :return: ``(b, r)``
    """
    if b < 0 or r < 0:
        raise ValueError("snc_env_tokenbucket: b and r must be nonnegative, "
                         "got %g, %g." % (b, r))
    if theta is not None and theta <= 0:
        raise ValueError("snc_env_tokenbucket: theta must be positive, got %g." % theta)
    return float(b), float(r)


def snc_env_map(D0, D1, theta):
    """MGF arrival envelope of a MAP/MMPP flow with unit-size jobs.

    With ``A(theta) = D0 + D1*exp(theta)``, lstar its eigenvalue of maximal real
    part and v > 0 the corresponding right Perron eigenvector,

        rho(theta) = lstar/theta,  sigma(theta) = log(max(v)/min(v))/theta.

    The burst term is what the modulating chain contributes: 0 for a one-phase
    MAP (which reproduces :func:`snc_env_poisson` exactly) and positive for an
    MMPP.

    :param D0: hidden-transition generator block
    :param D1: arrival-transition block
    :param theta: Chernoff parameter, theta > 0
    :return: ``(sigma, rho)``
    """
    if theta <= 0:
        raise ValueError("snc_env_map: theta must be positive, got %g." % theta)
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    if D0.shape != D1.shape or D0.shape[0] != D0.shape[1]:
        raise ValueError("snc_env_map: D0 and D1 must be square and of equal size.")
    A = D0 + D1 * math.exp(theta)
    w, V = np.linalg.eig(A)
    imax = int(np.argmax(np.real(w)))
    lstar = float(np.real(w[imax]))
    v = np.real(V[:, imax])
    if v.max() < 0:
        v = -v  # the eigenvector sign is arbitrary, take the positive representative
    if v.min() <= 0:
        raise ValueError("snc_env_map: MAP is not irreducible: the Perron eigenvector "
                         "is not positive.")
    return math.log(v.max() / v.min()) / theta, lstar / theta


def snc_srv_rate(C, theta=None):
    """MGF service envelope of a constant-rate work-conserving server.

    ``S(s,t) = C*(t-s)``, so ``(sigma, rho) = (0, C)`` for every theta.

    :param C: server capacity, work per slot
    :return: ``(0, C)``
    """
    if C <= 0:
        raise ValueError("snc_srv_rate: C must be positive, got %g." % C)
    if theta is not None and theta <= 0:
        raise ValueError("snc_srv_rate: theta must be positive, got %g." % theta)
    return 0.0, float(C)


def snc_srv_exp(mu, theta):
    """MGF service envelope of an exponential server, in JOB units.

    A single server with Exp(mu) service times completes jobs at the epochs of a
    Poisson process of rate mu while busy, so ``rho(theta) =
    mu*(1-exp(-theta))/theta`` and the burst is zero.

    THIS IS THE SERVICE ELEMENT TO USE WHENEVER THE WORK UNIT IS THE JOB.
    Pairing :func:`snc_srv_rate` with a job-counting arrival envelope would model
    an M/D/1 and UNDERSTATE the delay of an exponential server rather than bound
    it. Job units also compose across hops, which service-time work units do not.

    :param mu: service rate, jobs per slot
    :param theta: Chernoff parameter, theta > 0
    :return: ``(0, rho)``
    """
    if mu <= 0:
        raise ValueError("snc_srv_exp: mu must be positive, got %g." % mu)
    if theta <= 0:
        raise ValueError("snc_srv_exp: theta must be positive, got %g." % theta)
    return 0.0, mu * (1.0 - math.exp(-theta)) / theta


def snc_leftover(sigmaS, rhoS, sigmaX, rhoX):
    """Leftover service envelope under blind (arbitrary) multiplexing.

    ``rho = rhoS - rhoX``, ``sigma = sigmaS + sigmaX``. Exact for independent
    flow and cross traffic; the dependent case needs a Hoelder split that this
    elementary version does not implement. A nonpositive rho means the cross
    traffic can exhaust the server, which the bound functions report as a
    violation probability of 1.
    """
    return sigmaS + sigmaX, rhoS - rhoX


def snc_conv(sigma1, rho1, sigma2, rho2, theta, delta=None):
    """Min-plus convolution of two service envelopes (tandem concatenation).

    ``rho = min(rho1,rho2)`` and ``sigma = sigma1 + sigma2 -
    log(1-exp(-theta*|rho1-rho2|))/theta``: the end-to-end burst grows
    additively rather than the per-station delay bounds being summed, which is
    the pay-bursts-only-once result. Equal rates make the series diverge, so the
    slower server is shifted down by ``delta`` (default ``1e-2*min(rho1,rho2)``).
    """
    if theta <= 0:
        raise ValueError("snc_conv: theta must be positive, got %g." % theta)
    if delta is None:
        delta = 1e-2 * min(rho1, rho2)
    if delta <= 0:
        raise ValueError("snc_conv: delta must be positive, got %g." % delta)
    if not math.isfinite(rho1) or not math.isfinite(rho2):
        return float('inf'), min(rho1, rho2)
    gap = abs(rho1 - rho2)
    if gap <= delta:
        gap = delta  # equal rates: shift the slower server down to close the series
        rho = min(rho1, rho2) - delta
    else:
        rho = min(rho1, rho2)
    if rho <= 0:
        return float('inf'), rho
    return sigma1 + sigma2 - math.log(1.0 - math.exp(-theta * gap)) / theta, rho


def snc_output(sigmaA, rhoA, sigmaS, rhoS, theta):
    """Output (departure) arrival envelope of a flow leaving a server.

    The rate is conserved and the server adds burstiness:
    ``sigma = sigmaA + sigmaS - log(1-exp(-theta*(rhoS-rhoA)))/theta``. This is
    what carries a flow across a feed-forward network one hop at a time; for a
    tandem traversed by the same flow, :func:`snc_conv` is tighter.
    """
    if theta <= 0:
        raise ValueError("snc_output: theta must be positive, got %g." % theta)
    if not math.isfinite(rhoA) or not math.isfinite(rhoS) or rhoS <= rhoA:
        return float('inf'), rhoA  # unstable station, no exponential-form envelope
    return sigmaA + sigmaS - math.log(1.0 - math.exp(-theta * (rhoS - rhoA))) / theta, rhoA
