"""
Discrete-time (slotted) single-queue formulas.

Time advances in slots of unit length. These are the only discrete-time queue
solvers in the API; every other "geometric" in the package is matrix-geometric
or a geometric series.

References:
    Original MATLAB: matlab/src/api/qsys/qsys_geogeo1.m, qsys_geoxgeo1.m
    H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
    Springer 2001, section 2.1 (corollary 2.7) and chapter 6.
"""

from typing import Callable, Dict

__all__ = ['qsys_geogeo1', 'qsys_geoxgeo1', 'qsys_geoxgeo1_moments']

# Observation epochs. Compared by name across codebases, matching the Java
# GeoGeo1Convention enum constants.
LAS_DA = 'LAS_DA'
EAS = 'EAS'


def _check_convention(convention: str) -> str:
    if not isinstance(convention, str):
        raise ValueError("convention must be the string 'LAS_DA' or 'EAS'")
    normalized = convention.upper()
    if normalized not in (LAS_DA, EAS):
        raise ValueError("convention must be 'LAS_DA' or 'EAS', got %r" % (convention,))
    return normalized


def qsys_geogeo1(a: float, s: float, convention: str = LAS_DA) -> Dict[str, object]:
    """
    Analyze a discrete-time Geo/Geo/1 queue.

    In each slot an arrival occurs with probability ``a`` and, if the server is
    engaged, a service completion occurs with probability ``s``, independently
    of everything else. The system content at slot boundaries is a discrete
    birth-death chain with a geometric stationary distribution.

    The two conventions are not two systems. They are one system, Daduna's
    LA-rule (events at the end of their slot) with the D/A-rule (departure
    resolved before arrival), observed at two instants. With
    ``X(t+1) = X(t) - D(t) + A(t)``, ``LAS_DA`` is the law of ``X``, taken
    after both events, and ``EAS`` is the law of ``Y(t) = X(t) - D(t)``, taken
    after the departure and before the arrival. They are one departure apart,
    so the mean contents differ by exactly ``a`` and the sojourn times by one
    slot; the queueing delay is the same under both.

    Args:
        a: Per-slot arrival probability, 0 < a < s
        s: Per-slot service completion probability, 0 < s <= 1
        convention: Observation epoch, 'LAS_DA' (default) or 'EAS'

    Returns:
        dict with keys:
            - convention: epoch used
            - arrivalProb, serviceProb: the inputs
            - utilization: a/s
            - throughput: a
            - emptyProb: probability the system is empty at the epoch
            - ratio: geometric decay ratio r = a(1-s)/(s(1-a))
            - meanQueueLength, meanWaitingQueue
            - meanSojournTime, meanWaitingTime, meanServiceTime
            - pmf: callable n -> stationary probability of n jobs
            - analyzer: 'qsys_geogeo1'

    Example:
        >>> r = qsys_geogeo1(0.2, 0.5)
        >>> print(f"{r['meanSojournTime']:.4f}")
        2.6667
    """
    convention = _check_convention(convention)
    if not (0.0 < a <= 1.0):
        raise ValueError("arrival probability a=%r must lie in (0,1]" % (a,))
    if not (0.0 < s <= 1.0):
        raise ValueError("service probability s=%r must lie in (0,1]" % (s,))
    if not a < s:
        raise ValueError("load a/s=%r must be strictly less than 1" % (a / s,))

    rho = a / s
    ratio = a * (1.0 - s) / (s * (1.0 - a))

    # The queueing delay does not depend on the observation epoch; only the
    # accounting of the slot in which service takes place does.
    mean_waiting_time = a * (1.0 - s) / (s * (s - a))
    mean_waiting_queue = a * mean_waiting_time

    if convention == LAS_DA:
        empty_prob = 1.0 - rho
        mean_queue_length = a * (1.0 - a) / (s - a)
        mean_sojourn_time = (1.0 - a) / (s - a)
        mean_service_time = 1.0 / s

        def pmf(n: int) -> float:
            _check_queue_length(n)
            if n == 0:
                return empty_prob
            return empty_prob * (rho / (1.0 - a)) * ratio ** (n - 1)
    else:
        empty_prob = 1.0 - ratio
        mean_queue_length = a * (1.0 - s) / (s - a)
        mean_sojourn_time = (1.0 - s) / (s - a)
        mean_service_time = (1.0 - s) / s

        def pmf(n: int) -> float:
            _check_queue_length(n)
            return (1.0 - ratio) * ratio ** n

    return {
        'convention': convention,
        'arrivalProb': a,
        'serviceProb': s,
        'utilization': rho,
        'throughput': a,
        'emptyProb': empty_prob,
        'ratio': ratio,
        'meanQueueLength': mean_queue_length,
        'meanWaitingQueue': mean_waiting_queue,
        'meanSojournTime': mean_sojourn_time,
        'meanWaitingTime': mean_waiting_time,
        'meanServiceTime': mean_service_time,
        'pmf': pmf,
        'analyzer': 'qsys_geogeo1',
    }


def _check_queue_length(n: int) -> None:
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError("queue length must be a non-negative integer, got %r" % (n,))


def qsys_geoxgeo1(a: float, beta: float, s: float,
                  convention: str = LAS_DA) -> Dict[str, object]:
    """
    Analyze a discrete-time Geo^X/Geo/1 queue with geometric batch sizes.

    In each slot a batch arrives with probability ``a``; the batch size is
    geometric on {1,2,...} with parameter ``beta``, so ``E[X] = 1/beta`` and
    ``E[X(X-1)] = 2(1-beta)/beta**2``. Stability requires
    ``lambda = a*E[X] < s``. At ``beta == 1`` the batch is always a single job
    and the result equals :func:`qsys_geogeo1`.

    Args:
        a: Per-slot probability that a batch arrives, 0 < a <= 1
        beta: Batch-size geometric parameter, 0 < beta <= 1
        s: Per-slot service completion probability, 0 < s <= 1
        convention: Observation epoch, 'LAS_DA' (default) or 'EAS'

    Returns:
        dict, see :func:`qsys_geoxgeo1_moments`

    Example:
        >>> r = qsys_geoxgeo1(0.1, 0.5, 0.9)
        >>> print(f"{r['arrivalRate']:.4f}")
        0.2000
    """
    if not (0.0 < beta <= 1.0):
        raise ValueError("batch geometric parameter beta=%r must lie in (0,1]" % (beta,))
    batch_mean = 1.0 / beta
    batch_second_factorial = 2.0 * (1.0 - beta) / (beta * beta)
    return qsys_geoxgeo1_moments(a, batch_mean, batch_second_factorial, s, convention)


def qsys_geoxgeo1_moments(a: float, batch_mean: float, batch_second_factorial: float,
                          s: float, convention: str = LAS_DA) -> Dict[str, object]:
    """
    Analyze a discrete-time Geo^X/Geo/1 queue for an arbitrary batch law.

    The batch enters the solution only through its first two factorial moments,
    so specifying those is fully general. With ``A(z) = 1-a+a*X(z)`` the pgf of
    the number of jobs arriving in one slot, the slot-boundary content obeys
    ``X(t+1) = X(t) - D(t) + A(t)`` with the departure resolved first, giving::

        P(z) = p0 s (z-1) A(z) / ( z - A(z)(s+(1-s)z) ),  p0 = 1 - lambda/s
        E[N] = lambda + ( a E[X(X-1)]/2 + lambda(1-s) ) / (s - lambda)

    Args:
        a: Per-slot probability that a batch arrives, 0 < a <= 1
        batch_mean: E[X], at least 1 since an arriving batch carries a job
        batch_second_factorial: E[X(X-1)], non-negative and at least
            ``batch_mean**2 - batch_mean``
        s: Per-slot service completion probability, 0 < s <= 1
        convention: Observation epoch, 'LAS_DA' (default) or 'EAS'

    Returns:
        dict with keys:
            - convention, batchArrivalProb, batchMean,
              batchSecondFactorialMoment, serviceProb
            - arrivalRate, throughput: lambda = a*E[X]
            - utilization: lambda/s
            - boundaryEmptyProb: 1 - lambda/s, the empty probability AT THE
              SLOT BOUNDARY under both conventions
            - meanQueueLength, meanWaitingQueue
            - meanSojournTime, meanWaitingTime, meanServiceTime
            - pgf: callable (z, A_of_z) -> P(z), defined for 0 < z <= 1
            - analyzer: 'qsys_geoxgeo1'

    No pmf is returned: for a general batch law the stationary distribution has
    no elementary closed form, so only the generating function is exact.
    """
    convention = _check_convention(convention)
    if not (0.0 < a <= 1.0):
        raise ValueError("batch arrival probability a=%r must lie in (0,1]" % (a,))
    if not (0.0 < s <= 1.0):
        raise ValueError("service probability s=%r must lie in (0,1]" % (s,))
    if not batch_mean >= 1.0:
        raise ValueError("mean batch size E[X]=%r must be at least 1: a batch that "
                         "arrives carries at least one job" % (batch_mean,))
    if batch_second_factorial < 0.0:
        raise ValueError("E[X(X-1)]=%r must be non-negative" % (batch_second_factorial,))
    # E[X^2] = E[X(X-1)] + E[X] >= E[X]^2, so a smaller second factorial moment
    # describes no random variable at all.
    min_second_factorial = batch_mean * batch_mean - batch_mean
    if batch_second_factorial < min_second_factorial - 1e-9 * max(1.0, min_second_factorial):
        raise ValueError("E[X(X-1)]=%r is below E[X]^2-E[X]=%r, so the batch moments "
                         "describe no random variable"
                         % (batch_second_factorial, min_second_factorial))

    lambda_val = a * batch_mean
    if not lambda_val < s:
        raise ValueError("load lambda/s=%r must be strictly less than 1" % (lambda_val / s,))

    rho = lambda_val / s
    boundary_empty_prob = 1.0 - rho

    # P'(1) from the generating function.
    mean_at_boundary = (lambda_val
                        + (a * batch_second_factorial / 2.0 + lambda_val * (1.0 - s))
                        / (s - lambda_val))

    mean_sojourn_at_boundary = mean_at_boundary / lambda_val
    mean_waiting_time = mean_sojourn_at_boundary - 1.0 / s
    mean_waiting_queue = lambda_val * mean_waiting_time

    if convention == LAS_DA:
        mean_queue_length = mean_at_boundary
        mean_sojourn_time = mean_sojourn_at_boundary
        mean_service_time = 1.0 / s
    else:
        # One departure earlier: the epoch drops exactly the departures of the
        # slot, whose rate is lambda, hence one slot of sojourn.
        mean_queue_length = mean_at_boundary - lambda_val
        mean_sojourn_time = mean_sojourn_at_boundary - 1.0
        mean_service_time = (1.0 - s) / s

    def pgf(z: float, batch_pgf_at_z: float) -> float:
        if not (0.0 < z <= 1.0):
            raise ValueError("pgf argument z=%r must lie in (0,1]" % (z,))
        if z == 1.0:
            return 1.0
        denom = z - batch_pgf_at_z * (s + (1.0 - s) * z)
        boundary = boundary_empty_prob * s * (z - 1.0) * batch_pgf_at_z / denom
        if convention == LAS_DA:
            return boundary
        return boundary * (s / z + 1.0 - s) + boundary_empty_prob * s * (1.0 - 1.0 / z)

    return {
        'convention': convention,
        'batchArrivalProb': a,
        'batchMean': batch_mean,
        'batchSecondFactorialMoment': batch_second_factorial,
        'serviceProb': s,
        'arrivalRate': lambda_val,
        'throughput': lambda_val,
        'utilization': rho,
        'boundaryEmptyProb': boundary_empty_prob,
        'meanQueueLength': mean_queue_length,
        'meanWaitingQueue': mean_waiting_queue,
        'meanSojournTime': mean_sojourn_time,
        'meanWaitingTime': mean_waiting_time,
        'meanServiceTime': mean_service_time,
        'pgf': pgf,
        'analyzer': 'qsys_geoxgeo1',
    }
