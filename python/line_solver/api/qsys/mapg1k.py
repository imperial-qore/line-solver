"""
Exact MAP/G/1/K and MMAP[K]/G/1/K finite-buffer analysis.

Key functions:
    qsys_mapg1k: exact MAP/G/1/K with tail drop, arbitrary service law
    qsys_mmapg1k: exact per-class loss/throughput for MMAP[K]/G/1/K
    qsys_mapg1k_perflow: per-flow Palm-Khinchin approximation for N MAPs

References:
    Original MATLAB: matlab/src/api/qsys/qsys_mapg1k.m,
                     matlab/src/api/qsys/qsys_mmapg1k.m,
                     matlab/src/api/qsys/qsys_mapg1k_perflow.m
    Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied System
        Innovation 2026, 9, 112.
    Niu, Z.; Cooper, R.B. Transform-Free Analysis of M/G/1/K and Related
        Queues. Mathematics of Operations Research 1993, 18, 486-510.
"""

import warnings
from typing import Any, Dict, List, Sequence

import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln

from ..mc.dtmc import dtmc_solve
from ..mam.map_analysis import map_lambda, map_prob

__all__ = ['qsys_mapg1k', 'qsys_mmapg1k', 'qsys_mapg1k_perflow']


def _svc_mean(svc: Dict[str, Any]) -> float:
    """Mean service time of the descriptor."""
    stype = str(svc['type']).lower()
    if stype == 'gamma':
        return float(svc['alpha']) * float(svc['theta'])
    if stype == 'det':
        return float(svc['d'])
    if stype == 'ph':
        alpha = np.asarray(svc['alpha'], dtype=float).reshape(-1)
        T = np.asarray(svc['T'], dtype=float)
        return float(-alpha @ np.linalg.solve(T, np.ones(T.shape[0])))
    if stype == 'density':
        tmax = float(svc.get('tmax', np.inf))
        return _dquad(lambda x: x, svc['pdf'], tmax)
    raise ValueError("qsys_mapg1k: unsupported service type '%s'." % svc['type'])


def _finite(y):
    y = np.asarray(y, dtype=float)
    return np.where(np.isfinite(y), y, 0.0)


def _dquad(w, pdf, tmax: float) -> float:
    """E[w(S)] for a service law given by a density, under x = exp(u).

    An integrable density may diverge at the origin, which caps adaptive
    quadrature on [0,tmax]. The Jacobian exp(u) turns x^(alpha-1)dx into
    exp(alpha*u)du, which decays smoothly as u -> -Inf for any alpha > 0, so
    the singularity disappears rather than being resolved. The transformed
    integrand tends to 0 at both ends, but in floating point those limits are
    reached as 0*Inf, so the NaN produced there is replaced by its limit.
    """
    ulim = np.log(tmax) if np.isfinite(tmax) else np.inf

    def integrand(u):
        x = np.exp(u)
        return float(_finite(np.asarray(w(x), dtype=float)
                             * np.asarray(pdf(x), dtype=float) * x))

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        val, _ = quad(integrand, -np.inf, ulim, epsabs=1e-300, epsrel=1e-13,
                      limit=400)
    return float(val)


def _cquad(pdf, theta: float, nn, tmax: float) -> np.ndarray:
    nn = np.atleast_1d(np.asarray(nn))
    out = np.zeros(nn.size)
    for i, n in enumerate(nn):
        n = float(n)
        out[i] = _dquad(
            lambda x, n=n: np.exp(-theta * x + n * np.log(theta * x) - gammaln(n + 1)),
            pdf, tmax)
    return out


class _PhBlock:
    """c_n = theta^n * alpha * (theta*I-T)^{-(n+1)} * t, with cached powers."""

    def __init__(self, alpha: np.ndarray, Minv: np.ndarray, tv: np.ndarray,
                 theta: float):
        self.alpha = alpha
        self.Minv = Minv
        self.tv = tv
        self.theta = theta

    def __call__(self, nn) -> np.ndarray:
        nn = np.atleast_1d(np.asarray(nn))
        out = np.zeros(nn.size)
        for i, n in enumerate(nn):
            n = int(n)
            P = np.linalg.matrix_power(self.Minv, n + 1)
            out[i] = self.theta ** n * float(self.alpha @ (P @ self.tv))
        return out


def _guess(m: float, cap: int) -> int:
    """Initial uniformization order: mean plus a generous deviation allowance."""
    return int(min(cap, max(32, np.ceil(m + 10 * np.sqrt(max(m, 1.0)) + 32))))


def _grow(fn, tol: float, cap: int, n0: int) -> np.ndarray:
    """Build c_0..c_N in blocks, stopping when the series sums to 1 within tol
    or when a whole block adds nothing in floating point, i.e. the
    representable series is exhausted. The second criterion terminates paths
    whose terms are known only to quadrature accuracy, where the first can
    never be met.
    """
    cn = np.asarray(fn(np.arange(0, n0 + 1)), dtype=float).reshape(-1)
    while cn.size < cap:
        if abs(1.0 - cn.sum()) <= tol:
            break
        n = cn.size
        add = np.asarray(fn(np.arange(n, min(cap - 1, n + 63) + 1)),
                         dtype=float).reshape(-1)
        if add.size == 0:
            break
        cn = np.concatenate([cn, add])
        if add.sum() <= np.finfo(float).eps * cn.sum():
            break
    return cn


def _service_series(svc: Dict[str, Any], theta: float, tol: float, nmax_cap: int):
    """c_n = E[exp(-theta*S)*(theta*S)^n/n!] for n = 0..nmax, and the mean S.

    sum_{n>=0} c_n = E[exp(-theta*S)*exp(theta*S)] = 1 exactly, which both sets
    the truncation order and certifies it.
    """
    if 'type' not in svc:
        raise ValueError("qsys_mapg1k: service descriptor must have a 'type' field.")
    smean = _svc_mean(svc)
    stype = str(svc['type']).lower()
    if stype == 'gamma':
        al = float(svc['alpha'])
        th = float(svc['theta'])

        def fn(nn):
            nn = np.asarray(nn, dtype=float)
            return np.exp(nn * np.log(theta * th) - gammaln(nn + 1)
                          + gammaln(al + nn) - gammaln(al)
                          - (al + nn) * np.log1p(th * theta))
    elif stype == 'det':
        d = float(svc['d'])

        def fn(nn):
            nn = np.asarray(nn, dtype=float)
            return np.exp(-theta * d + nn * np.log(theta * d) - gammaln(nn + 1))
    elif stype == 'ph':
        alpha = np.asarray(svc['alpha'], dtype=float).reshape(-1)
        T = np.asarray(svc['T'], dtype=float)
        tv = -T @ np.ones(T.shape[0])
        Minv = np.linalg.inv(theta * np.eye(T.shape[0]) - T)
        fn = _PhBlock(alpha, Minv, tv, theta)
    elif stype == 'density':
        tmax = float(svc.get('tmax', np.inf))

        def fn(nn):
            return _cquad(svc['pdf'], theta, nn, tmax)
    else:
        raise ValueError("qsys_mapg1k: unsupported service type '%s'." % svc['type'])

    cn = _grow(fn, tol, nmax_cap, _guess(theta * smean, nmax_cap))
    if abs(1.0 - cn.sum()) > 1e-6:
        warnings.warn('qsys_mapg1k: uniformization series truncated at n=%d with '
                      'residual %g; increase nmax.' % (cn.size - 1, abs(1.0 - cn.sum())))
    return cn, smean


def qsys_mapg1k(D0, D1, svc: Dict[str, Any], K: int, tol: float = 1e-12,
                nmax: int = 200000) -> Dict[str, Any]:
    """
    Exact analysis of a MAP/G/1/K queue with tail drop.

    Markovian arrivals, arbitrary service time distribution F, and a finite
    buffer of K packets (the position held by the packet in transmission
    included). Unlike qsys_mapg1 the service time is NOT fitted to a
    phase-type distribution: F enters exactly, through the functionals A_m and
    Q_m evaluated by uniformization of the arrival MAP.

    Args:
        D0, D1: MAP parameter matrices (M x M), D0 + D1 an irreducible generator
        svc: service time descriptor dict with key 'type':
            'gamma'  : keys alpha (shape), theta (scale)
            'det'    : key d (constant service time)
            'ph'     : keys alpha (1 x p), T (p x p subgenerator)
            'density': key pdf (callable), optional key tmax
        K: buffer size in packets, K >= 1
        tol: uniformization truncation tolerance
        nmax: cap on the uniformization order

    Returns:
        Dict with p0, pK, lossProbability, throughput, lambda, meanServiceTime,
        utilization, rho, nmax, sigma, pKvec, p0vec, plevel, meanQueueLength.

    Method:
        The chain embedded at departure epochs is used, in the state (n,j):
        n = 0..K-1 packets left behind by a departure, j = MAP phase. With A_m
        the matrix of "m arrivals during a service, phase i -> j",
            n >= 1: n' = n-1+min(m, K-n), overflow sum_{m>=K-n} A_m
            n == 0: the phase first jumps by (-D0)^{-1}*D1 (the idle period
                    ends at an arrival), the service then proceeds as from n=1.
        Its stationary law sigma gives, by Markov renewal reward, the cycle
        mean, p0 and pK, where Q_m is the expected time within a service with
        exactly m arrivals so far. Time-stationary p0 and pK follow, so no
        PASTA assumption is needed on the MAP side.

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mapg1k.m
    """
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    M = D0.shape[0]
    if D0.shape[1] != M or D1.shape != (M, M):
        raise ValueError('qsys_mapg1k: D0 and D1 must be square matrices of equal size.')
    if K < 1 or K != round(K):
        raise ValueError('qsys_mapg1k: buffer size K must be a positive integer.')
    K = int(K)
    beta = -np.diag(D0)
    if np.any(beta <= 0):
        raise ValueError('qsys_mapg1k: D0 must have strictly negative diagonal entries.')

    # Uniformization constant: theta >= max_i beta_i keeps I+D0/theta substochastic
    theta = float(np.max(beta))

    # c_n = E[exp(-theta*S)*(theta*S)^n/n!], summing to f(0) = 1
    cn, smean = _service_series(svc, theta, tol, nmax)
    nmax_used = cn.size - 1

    # d_n weights the residual service time: see MATLAB qsys_mapg1k
    tailc = np.cumsum(cn[::-1])[::-1]
    dn = np.concatenate([tailc[1:], [0.0]]) / theta

    # A_m and Q_m for m = 0..K-1, plus B0 = sum_m A_m = E[exp((D0+D1)*S)]
    mmax = max(K - 1, 0)
    A = np.zeros((mmax + 1, M, M))
    Q = np.zeros((mmax + 1, M, M))
    Sn = np.zeros((mmax + 1, M, M))
    Sn[0] = np.eye(M)
    B0 = np.zeros((M, M))
    Qtot = np.zeros((M, M))
    Pn = np.eye(M)
    Pt0 = np.eye(M) + D0 / theta
    Pt1 = D1 / theta
    PD = np.eye(M) + (D0 + D1) / theta
    for n in range(nmax_used + 1):
        for m in range(min(n, mmax) + 1):
            A[m] += Sn[m] * cn[n]
            Q[m] += Sn[m] * dn[n]
        B0 += Pn * cn[n]
        Qtot += Pn * dn[n]     # sum_m Q_m = int_0^inf exp(D*x)*(1-F(x))dx
        if n < nmax_used:
            Snew = np.zeros((mmax + 1, M, M))
            for m in range(min(n + 1, mmax) + 1):
                acc = np.zeros((M, M))
                if m <= n:
                    acc = acc + Sn[m] @ Pt0
                if m >= 1 and m - 1 <= n:
                    acc = acc + Sn[m - 1] @ Pt1
                Snew[m] = acc
            Sn = Snew
            Pn = Pn @ PD

    # 1-D, as the per-class block at the bottom of this file already builds it:
    # with a COLUMN e every `float(row @ ... @ e)` below is a size-1 ARRAY, and
    # numpy 2 no longer coerces one to a scalar, so idle_time and both time_l
    # accumulators raised TypeError instead of computing.
    e = np.ones(M)
    negD0inv = np.linalg.inv(-D0)
    Psi = negD0inv @ D1           # phase at the arrival that ends an idle period
    idle = negD0inv @ e           # expected idle time from each phase

    # Embedded chain at departure epochs, state (n,j) -> index n*M+j
    P = np.zeros((K * M, K * M))
    last = slice((K - 1) * M, K * M)
    for n in range(1, K):
        rows = slice(n * M, (n + 1) * M)
        Bacc = B0.copy()
        for m in range(0, K - n):
            col = slice((n - 1 + m) * M, (n + m) * M)
            P[rows, col] += A[m]
            Bacc -= A[m]
        # Bacc = sum_{m>=K-n} A_m: every further arrival overflows the buffer
        P[rows, last] += Bacc
    rows = slice(0, M)
    Bacc = B0.copy()
    for m in range(0, K - 1):
        col = slice(m * M, (m + 1) * M)
        P[rows, col] += Psi @ A[m]
        Bacc -= A[m]
    P[rows, last] += Psi @ Bacc

    rowdev = float(np.max(np.abs(P.sum(axis=1) - 1.0)))
    if rowdev > 1e-8:
        raise ValueError('qsys_mapg1k: embedded chain rows deviate from 1 by %.2e. '
                         'The uniformization series for A_m has not converged; '
                         'raise nmax.' % rowdev)

    sigma = np.asarray(dtmc_solve(P), dtype=float).reshape(-1)
    sigma0 = sigma[0:M]

    # Markov renewal reward over the interval between successive departures
    idle_time = float(sigma0 @ idle)
    ecyc = smean + idle_time
    T = 1.0 / ecyc
    p0 = idle_time / ecyc

    # Qcum[r] = sum_{m=0}^{r} Q_m
    Qcum = np.zeros((mmax + 1, M, M))
    acc = np.zeros((M, M))
    for m in range(mmax + 1):
        acc = acc + Q[m]
        Qcum[m] = acc
    time_kvec = np.zeros(M)
    for n in range(1, K):
        r = K - n - 1
        time_kvec += sigma[n * M:(n + 1) * M] @ (Qtot - Qcum[r])
    if K >= 2:
        time_kvec += sigma0 @ Psi @ (Qtot - Qcum[K - 2])
    else:
        time_kvec += sigma0 @ Psi @ Qtot
    pkvec = time_kvec / ecyc
    pK = float(pkvec.sum())

    # Time-stationary level law from the same renewal-reward decomposition
    time_l = np.zeros(K + 1)
    time_l[0] = idle_time
    for n in range(1, K):
        sn_row = sigma[n * M:(n + 1) * M]
        for l in range(n, K):
            time_l[l] += float(sn_row @ Q[l - n] @ e)
    s0psi = sigma0 @ Psi
    for l in range(1, K):
        time_l[l] += float(s0psi @ Q[l - 1] @ e)
    time_l[K] = float(time_kvec.sum())
    plevel = time_l / ecyc
    massdev = abs(plevel.sum() - 1.0)
    if massdev > 1e-8:
        raise ValueError('qsys_mapg1k: level distribution has mass %.12f. The Q_m '
                         'series has not converged; raise nmax.' % plevel.sum())
    meanq = float(np.arange(K + 1) @ plevel)
    # Time at level 0 is the idle period alone, whose phase law is (-D0)^{-1}
    p0vec = (sigma0 @ negD0inv) / ecyc

    lam = float(map_lambda(D0, D1))

    return {
        'p0': p0,
        'pK': pK,
        'throughput': T,
        'lossProbability': 1.0 - T / lam,
        'lambda': lam,
        'meanServiceTime': smean,
        'utilization': 1.0 - p0,
        'rho': lam * smean,
        'nmax': nmax_used,
        'sigma': sigma,
        'pKvec': pkvec,
        'p0vec': p0vec,
        'plevel': plevel,
        'meanQueueLength': meanq,
        'analyzer': 'qsys_mapg1k',
    }


def qsys_mmapg1k(D0, D1c: Sequence, svc: Dict[str, Any], K: int,
                 tol: float = 1e-12, nmax: int = 200000) -> Dict[str, Any]:
    """
    Exact per-class throughput and loss ratio of an MMAP[K]/G/1/K queue.

    Two classes of equal arrival rate but different interarrival variability or
    autocorrelation receive different loss ratios. Aggregate-only finite-buffer
    analyses cannot express it: they return a single blocking probability p and
    set T_k = lambda_k*(1-p), making the loss ratio identical by construction.

    Args:
        D0: M x M hidden transition matrix of the arrival MMAP
        D1c: sequence of R matrices, D1c[k] = M x M arrival matrix of class k
        svc: service time descriptor, see qsys_mapg1k
        K: buffer size in packets, K >= 1

    Returns:
        Dict with throughput, lossRatio, lambda (all per class), the aggregate
        quantities, and the level/phase quantities of the driving MAP model.

    Method:
        The aggregate MAP {D0, sum_k D1c[k]} drives qsys_mapg1k, whose embedded
        chain returns the joint law of buffer level and MAP phase. A class-k
        arrival leaves phase i at rate (D1c[k]*e)_i, so
            lambda_k = pi*D1c[k]*e,  L_k = (pKvec*D1c[k]*e)/lambda_k.
        This is exact: no independence between classes is assumed and no PASTA
        argument is used, the phase resolution of pKvec doing the work.

        Assumes a single server and a service law that is iid and independent
        of class.

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mmapg1k.m
    """
    if isinstance(D1c, np.ndarray) and D1c.ndim == 2:
        raise ValueError('qsys_mmapg1k: D1c must be a sequence of per-class D1 matrices.')
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    M = D0.shape[0]
    R = len(D1c)
    D1 = np.zeros((M, M))
    Dk = []
    for k in range(R):
        Dki = np.atleast_2d(np.asarray(D1c[k], dtype=float))
        if Dki.shape != (M, M):
            raise ValueError('qsys_mmapg1k: D1c[%d] must be %dx%d.' % (k, M, M))
        Dk.append(Dki)
        D1 = D1 + Dki

    r = qsys_mapg1k(D0, D1, svc, K, tol=tol, nmax=nmax)

    e = np.ones(M)
    pit = np.asarray(map_prob(D0, D1), dtype=float).reshape(-1)

    lam = np.zeros(R)
    Lk = np.zeros(R)
    Tk = np.zeros(R)
    for k in range(R):
        lam[k] = float(pit @ Dk[k] @ e)
        if lam[k] > 0:
            Lk[k] = float(r['pKvec'] @ Dk[k] @ e) / lam[k]
        else:
            Lk[k] = 0.0
        Tk[k] = lam[k] * (1.0 - Lk[k])

    return {
        'throughput': Tk,
        'lossRatio': Lk,
        'lambda': lam,
        'lambdaAggregate': float(lam.sum()),
        'throughputAggregate': float(Tk.sum()),
        'lossAggregate': r['lossProbability'],
        'p0': r['p0'],
        'pK': r['pK'],
        'pKvec': r['pKvec'],
        'plevel': r['plevel'],
        'meanQueueLength': r['meanQueueLength'],
        'meanServiceTime': r['meanServiceTime'],
        'utilization': r['utilization'],
        'rho': r['rho'],
        'analyzer': 'qsys_mmapg1k',
    }


def qsys_mapg1k_perflow(MAPS: Sequence, svc: Dict[str, Any], K: int,
                        tol: float = 1e-12, nmax: int = 200000) -> Dict[str, Any]:
    """
    Per-flow throughput and loss ratio of a FIFO buffer fed by N flows.

    Flow n is described by its own MAP, so two flows may share an arrival rate
    and still differ in the shape and autocorrelation of their interarrival
    times. The buffer holds K packets including the one in transmission.

    Args:
        MAPS: sequence of N pairs (D0n, D1n); the orders M_n may differ
        svc: service time descriptor, see qsys_mapg1k
        K: buffer size in packets, K >= 1

    Returns:
        Dict with throughput, lossRatio, lambda, p0 and pK per flow, plus the
        aggregate quantities and rho.

    Method:
        The exact model of N flows would need prod_n M_n * (K+1) states.
        Instead one model per flow is solved: flow n is kept exactly as MAP_n
        while the other N-1 flows are replaced by a single Poisson stream of
        rate lambda - lambda_n, justified by the Palm-Khinchin limiting theorem
        on the superposition of many point processes. The superposition yields
            D0 = D0n - lambdaBar_n*I,  D1 = D1n + lambdaBar_n*I,
        which is passed to qsys_mapg1k. The sweep is O(N*(K*M)^3) against the
        O(M^(3N)*K^3) of the exact joint model.

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mapg1k_perflow.m
        Chydzinski, A. Applied System Innovation 2026, 9, 112, Theorem 1.
    """
    N = len(MAPS)
    if N < 1:
        raise ValueError('qsys_mapg1k_perflow: at least one flow is required.')

    lam = np.zeros(N)
    for n in range(N):
        if len(MAPS[n]) < 2:
            raise ValueError('qsys_mapg1k_perflow: MAPS[%d] must be a (D0,D1) pair.' % n)
        lam[n] = float(map_lambda(np.asarray(MAPS[n][0], dtype=float),
                                  np.asarray(MAPS[n][1], dtype=float)))
    lam_tot = float(lam.sum())

    Tn = np.zeros(N)
    Ln = np.zeros(N)
    p0 = np.zeros(N)
    pK = np.zeros(N)
    smean = np.nan
    for n in range(N):
        D0n = np.atleast_2d(np.asarray(MAPS[n][0], dtype=float))
        D1n = np.atleast_2d(np.asarray(MAPS[n][1], dtype=float))
        lam_bar = lam_tot - lam[n]
        # Superposition of MAP_n with a Poisson background of rate lam_bar
        Mn = D0n.shape[0]
        D0 = D0n - lam_bar * np.eye(Mn)
        D1 = D1n + lam_bar * np.eye(Mn)
        r = qsys_mapg1k(D0, D1, svc, K, tol=tol, nmax=nmax)
        smean = r['meanServiceTime']
        p0[n] = r['p0']
        pK[n] = r['pK']
        # Throughput of flow n: the aggregate departure rate (1-p0)/S less the
        # background throughput lam_bar*(1-pK), the background loss ratio being
        # pK by PASTA since the background is Poisson.
        Tn[n] = (1.0 - r['p0']) / smean + r['pK'] * lam_bar - lam_bar
        Ln[n] = 1.0 - Tn[n] / lam[n]

    return {
        'throughput': Tn,
        'lossRatio': Ln,
        'lambda': lam,
        'lambdaAggregate': lam_tot,
        'throughputAggregate': float(Tn.sum()),
        'lossAggregate': float((Ln * lam).sum() / lam_tot),
        'p0': p0,
        'pK': pK,
        'meanServiceTime': smean,
        'rho': lam_tot * smean,
        'analyzer': 'qsys_mapg1k_perflow',
    }
