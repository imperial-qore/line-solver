"""
Fixed-sample-size confidence intervals for steady-state quantiles.

``fquest`` takes a single sample path, ``firquest`` takes independent
replications. Both are automated: the user supplies a dataset, the quantile
probability and the nominal coverage, and the procedure chooses the batch size and
batch count itself and warns when the dataset is too small.

References:
    Original MATLAB: matlab/src/api/sim/sim_fquest.m, sim_firquest.m
    A. Lolos, C. Alexopoulos, D. Goldsman, K. D. Dingec, A. C. Mokashi,
    J. R. Wilson, "A Fixed-Sample-Size Method for Estimating Steady-State
    Quantiles", Proc. Winter Simulation Conference, 2023.
    Same authors, "A Fixed-Sample-Size Procedure for Estimating Steady-State
    Quantiles Based on Independent Replications", Proc. Winter Simulation
    Conference, 2025.
    R. Willink, "A Confidence Interval and Test for the Mean of an Asymmetric
    Distribution", Commun. Statist. Theory Methods 34, 2005.
"""

from math import ceil, exp, sqrt
from typing import Dict, List, Optional, Sequence

import numpy as np

from .dist import tinv
from .sts import DEFAULT_WEIGHT, sts_quantile_areas
from .tests import shapirowilk, vonneumann

__all__ = ['fquest', 'firquest', 'quest_options', 'firquest_batchcounts']


def quest_options(**overrides) -> Dict[str, object]:
    """
    Validate and complete the option set of the QUEST procedures.

    The defaults are the ones the articles report after their own
    experimentation: ``b0 = 50`` gives the warmup randomness test enough power,
    32 batches suffice to estimate the variance parameter while fewer than 10
    make the interval unreliable, and the decaying warmup significance keeps the
    batch size from growing so far that truncation eats a short sample. With
    these values the fourth warmup iteration runs at
    ``beta*exp(-0.2*3**2.3) = 0.025``.

    Recognized keys are ``b0``, ``m0``, ``s``, ``beta``, ``eta``, ``theta``,
    ``weight`` and ``force``.

    Args:
        **overrides: Options to replace

    Returns:
        The completed option dict.

    Raises:
        ValueError: If a key is unknown or a value is inadmissible.
    """
    opt = {'b0': 50, 'm0': 500, 's': [32, 24, 16, 10], 'beta': 0.30,
           'eta': 0.2, 'theta': 2.3, 'weight': DEFAULT_WEIGHT, 'force': True}
    for key, value in overrides.items():
        if key not in opt:
            raise ValueError("Unknown option %r" % (key,))
        opt[key] = value

    if int(opt['b0']) != opt['b0'] or opt['b0'] < 3:
        raise ValueError("b0 must be an integer >= 3, got %r" % (opt['b0'],))
    if int(opt['m0']) != opt['m0'] or opt['m0'] < 1:
        raise ValueError("m0 must be a positive integer, got %r" % (opt['m0'],))
    s = [int(v) for v in opt['s']]
    if not s or any(v < 1 for v in s):
        raise ValueError("s must be a nonempty sequence of positive integers")
    if any(s[i] >= s[i - 1] for i in range(1, len(s))):
        raise ValueError("s must be strictly decreasing")
    opt['s'] = s
    if not 0.0 < opt['beta'] < 1.0:
        raise ValueError("beta must lie in (0,1), got %r" % (opt['beta'],))
    if opt['eta'] < 0.0:
        raise ValueError("eta must be nonnegative, got %r" % (opt['eta'],))
    if opt['theta'] <= 0.0:
        raise ValueError("theta must be positive, got %r" % (opt['theta'],))
    if opt['weight'] == 0.0 or not np.isfinite(opt['weight']):
        raise ValueError("weight must be nonzero and finite, got %r" % (opt['weight'],))
    opt['force'] = bool(opt['force'])
    return opt


def firquest_batchcounts(R: int) -> List[int]:
    """
    Article default batch counts per replication, as a function of R.

    Chosen so that ``R*b`` pooled statistics remain enough to test while every
    replication still contributes at least one batch.

    Args:
        R: Number of replications, at least 2

    Returns:
        The descending batch counts.
    """
    if R < 2:
        raise ValueError("At least 2 replications are required, got %d" % R)
    if R == 2:
        return [14, 11, 8, 5]
    if R == 3:
        return [10, 8, 6, 4]
    if R == 4:
        return [6, 5, 4, 3]
    if R < 10:
        return [5, 4, 3, 2]
    if R < 17:
        return [4, 3, 2, 1]
    if R < 23:
        return [3, 2, 1]
    if R < 33:
        return [2, 1]
    return [1]


def _willink(zeta: float, gamma: float) -> float:
    if abs(gamma) <= 0.001:
        return zeta
    arg = 1.0 + 6.0 * gamma * (zeta - gamma)
    # the cube root is taken on the reals, the argument may turn negative for a
    # strongly skewed and small batch sample
    root = arg ** (1.0 / 3.0) if arg >= 0.0 else -((-arg) ** (1.0 / 3.0))
    return (root - 1.0) / (2.0 * gamma)


def _heuristic_ci(bqe: np.ndarray, centre: float, ap: float, np_est: float,
                  nstar: int, alpha: float, use_autocorr: bool):
    """
    Fallback interval used when a stage test fails.

    Three intervals are formed and the smallest interval containing all of them
    is returned: two symmetric ones of half-width
    ``max(t_{1-alpha/2,K} sqrt(Ap/n), t_{1-alpha/2,K-1} sqrt(Np/n))`` about the
    full-sample quantile and about the average batch quantile, and Willink's
    skewness-adjusted asymmetric interval. Taking the wider of the two variance
    components is deliberately conservative, since neither can be trusted once a
    stage test has failed.

    ``use_autocorr`` applies the residual-correlation factor
    ``max(sqrt((1+phi1)/(1-phi1)), 1)``; pass True for ``fquest``, where the batch
    quantiles come from one sample path, and False for ``firquest``, where they
    come from independent replications and the article drops it.
    """
    K = int(bqe.size)
    if K < 3:
        raise ValueError("The heuristic interval needs at least 3 batch quantiles,"
                         " got %d" % K)

    half = max(tinv(1.0 - alpha / 2.0, K) * sqrt(ap / nstar),
               tinv(1.0 - alpha / 2.0, K - 1) * sqrt(np_est / nstar))

    mean = float(bqe.mean())
    s2 = float(np.sum((bqe - mean) ** 2) / (K - 1))
    s2tilde = float(np.sum((bqe - centre) ** 2) / (K - 1))

    lower = min(centre - half, mean - half)
    upper = max(centre + half, mean + half)
    if s2 <= 0.0:
        return lower, upper

    skew = (K / ((K - 1.0) * (K - 2.0))) * float(np.sum(((bqe - mean) / sqrt(s2)) ** 3))
    gamma = skew / (6.0 * sqrt(K))

    varphi = 1.0
    if use_autocorr:
        phi1 = float(np.sum((bqe[:-1] - mean) * (bqe[1:] - mean)) / ((K - 1) * s2))
        if abs(phi1) < 1.0:
            varphi = max(sqrt((1.0 + phi1) / (1.0 - phi1)), 1.0)

    tq = tinv(1.0 - alpha / 2.0, K - 1)
    scale = varphi * sqrt(s2tilde / K)
    g1 = _willink(tq, gamma) * scale
    g2 = _willink(-tq, gamma) * scale

    lower = min(lower, centre - g1, centre - g2)
    upper = max(upper, centre - g1, centre - g2)
    return lower, upper


def _warmup(path: np.ndarray, b: int, m0: int, p: float, opt) -> (int, bool):
    """Grow the batch size until the signed areas pass the randomness test."""
    n = path.size
    m = m0
    if n < b * m:
        m = n // b
    if m < 1:
        raise ValueError("The sample path is too short for the initial batch count"
                         " b0 = %d" % b)
    ell = 1
    at_max = False
    passed = False
    while True:
        stats = sts_quantile_areas(path[:b * m], b, m, p, opt['weight'])
        sig = opt['beta'] * exp(-opt['eta'] * (ell - 1) ** opt['theta'])
        if not vonneumann(stats['areas'], sig)['reject']:
            passed = True
            break
        if at_max:
            break
        ell += 1
        m_next = int(round(m * sqrt(2.0)))
        if n < b * m_next and m_next != n // b:
            m = n // b
        else:
            m = m_next
            if n < b * m:
                m = n // b
                at_max = True
        if m < 1:
            break
    return m, passed


def _stage_tests(pool, s: Sequence[int], beta: float):
    """
    Run the four stage tests in order, stepping the batch count down through s.

    Returns the final pooled statistics, the batch count reached, and whether all
    four tests passed. ``pool`` is called with a batch count and returns the
    pooled statistics dict, or None when the batch size would fall below 1.
    """
    v = 0
    b = s[v]
    stats = pool(b)
    ok = True
    for stage in range(1, 5):
        while True:
            if stats is None:
                ok = False
                break
            sample = stats['areas'] if stage <= 2 else stats['bqe']
            reject = (vonneumann(sample, beta)['reject'] if stage % 2 == 1
                      else shapirowilk(sample, beta)['reject'])
            if not reject:
                break
            v += 1
            if v >= len(s):
                ok = False
                break
            b = s[v]
            stats = pool(b)
        if not ok:
            break
    return stats, b, ok


def fquest(y: Sequence[float], p: float, alpha: float = 0.05,
           options: Optional[Dict[str, object]] = None) -> Dict[str, object]:
    """
    FQUEST: fixed-sample-size confidence interval for a steady-state quantile.

    Returns a point estimate and an interval for the p-quantile of the
    steady-state marginal distribution of the output process whose single sample
    path is ``y``. The path has arbitrary fixed length; no sequential control of
    the run length is needed.

    The procedure has four blocks. **Warmup** grows the batch size by sqrt(2)
    until the signed STS areas pass von Neumann's randomness test at the decaying
    significance ``beta*exp(-eta*(l-1)**theta)``, which means the areas are
    approximately independent so any initialization bias sits in the first batch.
    **Truncation** deletes that first batch, and is the entire warmup treatment;
    there is no separate transient detector. **Batch-count selection** steps b
    down through ``s`` until four tests pass in order, von Neumann and
    Shapiro-Wilk on the signed areas then on the batched quantile estimators,
    checking the asymptotic properties the interval rests on. **Delivery** returns
    ``ytilde_p(n*) +- t_{1-alpha/2,2b-1} sqrt(V_p/n*)``, or the conservative
    fallback when a test failed.

    Coverage was measured on the article's own test bed, the waiting-time process
    of an M/M/1 queue with lambda = 0.8, mu = 1 started with 113 jobs in system,
    over 500 independent replications at N = 200000, giving a standard error near
    1%: 95.2% at p = 0.5, 96.2% at p = 0.9 and 95.6% at p = 0.99 against a nominal
    95%. The delivered half-width exceeds the empirically needed one by factors of
    1.16, 1.24 and 1.99 respectively, so the interval is conservative and
    increasingly so into the tail, consistent with the half-widths the article
    reports. On i.i.d. Exp(1) data, where ``sigma_p^2 = p(1-p)/f(y_p)^2`` is
    exact, both A_p and N_p are unbiased to within 7%.

    Two properties are worth knowing. Coverage is not monotone in the sample size
    over this range: a larger sample passes the stage tests more often and so
    reaches the conservative fallback less. And a substantial fraction of runs
    takes that fallback at all, 22% to 65% here and rising with p, so a delivered
    interval may well be the heuristic one; ``heuristic`` in the result says
    which. Fewer than about 100 replications cannot resolve a two-point
    difference in coverage, so do not read a small experiment as a defect.

    Applicability is a condition on the output process, not on the model that
    produced it. The theory needs geometric moment contraction (Wu 2005), which
    holds for ARMA series, a broad class of short-range-dependent linear and
    nonlinear processes, many Markov chains, and was proved for M/M/1 and
    non-heavy-tailed G/G/1 waiting times by Dingec et al. (2022); a density that is
    positive and differentiable at the quantile of interest; short-range dependence
    and an FCLT for the indicator process. M/M/1 is only the validation bed, chosen
    because its exact quantiles are known.

    Two practical exclusions follow. **Do not use this on integer-valued output
    such as a queue length**: the marginal has no density, the density-regularity
    condition fails, and the batched quantile has no Bahadur representation. Use it
    on continuous output, that is response, waiting and sojourn times. And
    heavy-tailed service, which can break geometric moment contraction and induce
    long-range dependence, is outside the theory.

    Args:
        y: The sample path
        p: Quantile probability in (0,1)
        alpha: Significance level in (0,1)
        options: Procedure constants, see :func:`quest_options`

    Returns:
        Dict with keys ``estimate``, ``lower``, ``upper``, ``halfwidth``, ``b``,
        ``m``, ``n``, ``R``, ``truncated``, ``Ap``, ``Np``, ``Vp``, ``heuristic``,
        ``warnings`` and ``analyzer``.

    Raises:
        ValueError: If the arguments are out of range or the path is too short.
    """
    opt = quest_options(**(options or {}))
    if not 0.0 < p < 1.0:
        raise ValueError("p must lie in (0,1), got %r" % (p,))
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0,1), got %r" % (alpha,))
    if opt['s'][-1] < 3:
        raise ValueError("min(s) = %d, but the stage tests need at least 3 batches"
                         % opt['s'][-1])

    path = np.asarray(y, dtype=float).ravel()
    N = path.size
    if not np.all(np.isfinite(path)):
        raise ValueError("The sample path must be finite")
    if N < opt['s'][-1] * 2:
        raise ValueError("The sample path holds %d observations, at least %d are needed"
                         % (N, opt['s'][-1] * 2))

    warnings: List[str] = []
    m, passed = _warmup(path, opt['b0'], opt['m0'], p, opt)
    if not passed:
        warnings.append("the warmup randomness test could not be passed at the largest"
                        " admissible batch size, the sample path is too short")

    truncated = m
    tail = path[truncated:]
    nstar_len = tail.size

    def pool(batch_count: int):
        size = nstar_len // batch_count
        if size < 1:
            return None
        return sts_quantile_areas(tail[nstar_len - batch_count * size:],
                                  batch_count, size, p, opt['weight'])

    stats, b, ok = _stage_tests(pool, opt['s'], opt['beta'])
    if stats is None:
        raise ValueError("The sample path is too short to form %d batches" % opt['s'][-1])

    n = stats['n']
    estimate = stats['quantile']
    if ok:
        half = tinv(1.0 - alpha / 2.0, 2 * b - 1) * sqrt(stats['Vp'] / n)
        lower = estimate - half
        upper = estimate + half
        heuristic = False
    else:
        warnings.append("a randomness or normality test failed at b = %d, the delivered"
                        " interval is heuristic" % opt['s'][-1])
        heuristic = True
        if opt['force']:
            lower, upper = _heuristic_ci(stats['bqe'], estimate, stats['Ap'],
                                         stats['Np'], n, alpha, True)
            half = (upper - lower) / 2.0
        else:
            lower = upper = half = float('nan')

    return {'estimate': estimate, 'lower': lower, 'upper': upper, 'halfwidth': half,
            'b': b, 'm': stats['m'], 'n': n, 'R': 1, 'truncated': truncated,
            'Ap': stats['Ap'], 'Np': stats['Np'], 'Vp': stats['Vp'],
            'heuristic': heuristic, 'warnings': warnings, 'analyzer': 'fquest'}


def firquest(y, p: float, alpha: float = 0.05,
             options: Optional[Dict[str, object]] = None) -> Dict[str, object]:
    """
    FIRQUEST: fixed-sample-size quantile interval from independent replications.

    The replicated counterpart of :func:`fquest`. ``y`` is an ``(n, R)`` array
    whose columns are the replicate sample paths, or a sequence of R equal-length
    vectors.

    It differs from :func:`fquest` in four places. The warmup randomness test runs
    independently on each replicate path, and the batch size it settles on may
    differ between replications. Truncation removes the largest of those batch
    sizes from the front of every replication, which is more aggressive on
    purpose: an untruncated transient common to all replications biases every
    replicate estimate the same way, and averaging cannot remove it. The four
    stage tests act on the ``R*b`` signed areas and ``R*b`` replicate batched
    quantile estimators pooled across replications. The delivered interval is
    ``ytilde_p(N*) +- t_{1-alpha/2,2Rb-1} sqrt(Vtilde_p/N*)`` with ``N* = R*b*m``,
    and the heuristic fallback drops the residual-autocorrelation correction since
    the pooled batch quantiles come from independent paths.

    Two defaults differ from FQUEST: ``b0 = 25`` rather than 50, and ``s`` is
    chosen from R by :func:`firquest_batchcounts` rather than fixed, because the
    stage tests act on the pooled statistics and so need fewer batches per
    replication. Supplying either explicitly overrides that choice.

    Independent replications shorten the correlation the estimator has to fight,
    and they parallelize, but they reintroduce initialization bias in every path,
    so a short run length per replication is worse here than in :func:`fquest`.
    The article reports slight undercoverage at p = 0.99 when the total sample is
    under 500000, down to 90.8%.

    Args:
        y: The replicate sample paths, as an (n, R) array or a sequence of R
            equal-length vectors
        p: Quantile probability in (0,1)
        alpha: Significance level in (0,1)
        options: Procedure constants, see :func:`quest_options`

    Returns:
        The same dict as :func:`fquest`, with ``R`` the replication count and
        ``b``, ``m`` per replication, so that ``n == R*b*m``.

    Raises:
        ValueError: If the arguments are out of range or the paths are too short.
    """
    paths = np.asarray(y, dtype=float)
    if paths.ndim == 1:
        raise ValueError("At least 2 replications are required, use fquest for a"
                         " single path")
    if paths.ndim != 2:
        raise ValueError("y must be two dimensional, got %d dimensions" % paths.ndim)
    if paths.shape[0] < paths.shape[1] and paths.shape[0] < 2:
        raise ValueError("At least 2 replications are required, use fquest for a"
                         " single path")
    n, R = paths.shape
    if R < 2:
        raise ValueError("At least 2 replications are required, use fquest for a"
                         " single path")
    if not np.all(np.isfinite(paths)):
        raise ValueError("The sample paths must be finite")

    supplied = dict(options or {})
    supplied.setdefault('b0', 25)
    supplied.setdefault('s', firquest_batchcounts(R))
    opt = quest_options(**supplied)

    if not 0.0 < p < 1.0:
        raise ValueError("p must lie in (0,1), got %r" % (p,))
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0,1), got %r" % (alpha,))
    if R * opt['s'][-1] < 3:
        raise ValueError("R*min(s) = %d pooled batches is below the 3 the stage tests"
                         " need" % (R * opt['s'][-1]))

    warnings: List[str] = []
    m_max = 0
    failed = False
    for r in range(R):
        m, passed = _warmup(paths[:, r], opt['b0'], opt['m0'], p, opt)
        failed = failed or not passed
        m_max = max(m_max, m)
    if failed:
        warnings.append("the warmup randomness test could not be passed in every"
                        " replication, the replicate paths are too short")

    truncated = m_max
    if truncated >= n:
        raise ValueError("The warmup batch size %d exhausts the replication length %d"
                         % (truncated, n))
    tail = paths[truncated:, :]
    nstar_len = tail.shape[0]

    def pool(batch_count: int):
        size = nstar_len // batch_count
        if size < 1:
            return None
        kept = tail[nstar_len - batch_count * size:, :]
        areas = np.empty((batch_count, R))
        bqe = np.empty((batch_count, R))
        for r in range(R):
            stats = sts_quantile_areas(kept[:, r], batch_count, size, p, opt['weight'])
            areas[:, r] = stats['areas']
            bqe[:, r] = stats['bqe']
        total = R * batch_count * size
        flat = kept.ravel(order='F')
        idx = ceil(total * p) - 1
        quantile = float(np.partition(flat, idx)[idx])
        K = R * batch_count
        ap = float(np.sum(areas.ravel(order='F') ** 2) / K)
        np_est = float(size * np.sum((bqe.ravel(order='F') - quantile) ** 2) / (K - 1))
        vp = (K * ap + (K - 1) * np_est) / (2 * K - 1)
        return {'areas': areas.ravel(order='F'), 'bqe': bqe.ravel(order='F'),
                'quantile': quantile, 'Ap': ap, 'Np': np_est, 'Vp': vp,
                'b': batch_count, 'm': size, 'n': total}

    pooled, b, ok = _stage_tests(pool, opt['s'], opt['beta'])
    if pooled is None:
        raise ValueError("The replicate paths are too short to form %d batches each"
                         % opt['s'][-1])

    nstar = pooled['n']
    estimate = pooled['quantile']
    if ok:
        half = tinv(1.0 - alpha / 2.0, 2 * R * b - 1) * sqrt(pooled['Vp'] / nstar)
        lower = estimate - half
        upper = estimate + half
        heuristic = False
    else:
        warnings.append("a randomness or normality test failed at b = %d per"
                        " replication, the delivered interval is heuristic"
                        % opt['s'][-1])
        heuristic = True
        if opt['force']:
            lower, upper = _heuristic_ci(pooled['bqe'], estimate, pooled['Ap'],
                                         pooled['Np'], nstar, alpha, False)
            half = (upper - lower) / 2.0
        else:
            lower = upper = half = float('nan')

    return {'estimate': estimate, 'lower': lower, 'upper': upper, 'halfwidth': half,
            'b': b, 'm': pooled['m'], 'n': nstar, 'R': R, 'truncated': truncated,
            'Ap': pooled['Ap'], 'Np': pooled['Np'], 'Vp': pooled['Vp'],
            'heuristic': heuristic, 'warnings': warnings, 'analyzer': 'firquest'}
