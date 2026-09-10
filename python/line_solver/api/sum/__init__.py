"""
Summation method (SUM) and closing method for queueing networks.

Native Python implementations of the summation method for closed
queueing networks, including the extended SUM (ESUM) node functions for
non-product-form networks with generally distributed service times, and
of the closing method for open and mixed non-product-form networks.

Key functions:
    sum_closed: summation method for closed networks (SUM/ESUM)
    sum_closing: closing method for open and mixed networks, solved by SUM

References:
    Original MATLAB: matlab/src/api/sum/sum_*.m
    G. Bolch, S. Greiner, H. de Meer, K.S. Trivedi, Queueing Networks and
    Markov Chains, 2nd ed., Wiley, 2006, Secs. 9.2, 10.1.4.4, and 10.1.5.
"""

import math

import numpy as np

__all__ = ['sum_closed', 'sum_closing']


def sum_closed(L, N, Z=None, mi=None, scv=None, tol=1e-6, maxiter=10000):
    """
    Summation method (SUM) for closed queueing networks, including the
    extended SUM (ESUM) node functions for non-product-form networks.

    The method expresses the mean queue length of each station as a
    function of its throughput, Ki = fi(lambdai), and solves the
    population constraint sum_i Ki = K. Single-class models are solved by
    bisection on the system throughput (Bolch et al., Sec. 9.2.1);
    multiclass models by Gauss-Seidel sweeps of per-class bisections on
    the population constraints, a robust alternative to the successive
    substitution of Sec. 9.2.2.

    Node functions:
    - Product-form stations (scv=1, or insensitive disciplines PS/LCFS-PR,
      for which the caller must pass scv=1): Eq. (9.15)/(9.19).
    - FCFS stations with general service (scv!=1): ESUM corrections,
      Eq. (10.88) for -/G/1 and Eq. (10.89) for -/G/m, with
      ai=(1+scv_i)/2 and Erlang-C waiting probability P_mi.
    - Infinite-server stations (mi=inf) and think times Z: Ki=lambdai*Li.

    Args:
        L: MxR service demand matrix, L[i,r] = e[i,r]/mu[i,r]
        N: 1xR population vector
        Z: 1xR think times (aggregated as a delay term)
        mi: Mx1 number of servers (np.inf for infinite-server stations)
        scv: MxR squared coefficient of variation of service times
        tol: convergence tolerance
        maxiter: maximum number of iterations

    Returns:
        (XN, QN, UN, RN, it): 1xR throughputs, MxR mean queue lengths,
        MxR utilizations (per-server for queueing stations, X*L for IS),
        MxR residence times QN/XN, iteration count.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    if L.shape[0] == 1 and L.shape[1] > 1 and np.ndim(N) == 0:
        L = L.T
    M, R = L.shape
    N = np.atleast_1d(np.asarray(N, dtype=float)).ravel()
    Z = np.zeros(R) if Z is None else np.atleast_1d(np.asarray(Z, dtype=float)).ravel()
    mi = np.ones(M) if mi is None else np.atleast_1d(np.asarray(mi, dtype=float)).ravel()
    scv = np.ones((M, R)) if scv is None else np.atleast_2d(np.asarray(scv, dtype=float))
    if scv.shape != (M, R):
        scv = scv.reshape(M, R)
    K = float(np.sum(N[np.isfinite(N)]))

    XN = np.zeros(R)
    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    RN = np.zeros((M, R))
    it = 0

    if K == 0:
        return XN, QN, UN, RN, it

    if R == 1:
        # single class: bisection on the throughput (Sec. 9.2.1)
        lambda_l = 0.0
        lambda_u = np.inf
        for i in range(M):
            if L[i, 0] > 0:
                if np.isinf(mi[i]):
                    lambda_u = min(lambda_u, K / L[i, 0])
                else:
                    lambda_u = min(lambda_u, mi[i] / L[i, 0])
        if Z[0] > 0:
            lambda_u = min(lambda_u, K / Z[0])
        if np.isinf(lambda_u):
            raise ValueError('sum_closed: all service demands are zero.')
        lam = lambda_u
        for it in range(1, maxiter + 1):
            lam = (lambda_l + lambda_u) / 2
            g = lam * Z[0] + np.sum(_sum_node_qlen(L, np.array([lam]), mi, scv, K)[:, 0])
            if abs(g - K) <= tol or (lambda_u - lambda_l) <= tol * lambda_u:
                break
            if g > K:
                lambda_u = lam
            else:
                lambda_l = lam
        XN[0] = lam
    else:
        # multiclass: Gauss-Seidel sweeps of per-class bisections
        for it in range(1, maxiter + 1):
            delta = 0.0
            for r in range(R):
                if N[r] == 0:
                    continue
                ub = np.inf
                for i in range(M):
                    if L[i, r] > 0:
                        if np.isinf(mi[i]):
                            ub = min(ub, K / L[i, r])
                        else:
                            rem = mi[i] - (XN @ L[i, :] - XN[r] * L[i, r])
                            ub = min(ub, max(rem, 0.0) / L[i, r])
                if Z[r] > 0:
                    ub = min(ub, N[r] / Z[r])
                if np.isinf(ub):
                    raise ValueError('sum_closed: all service demands are zero.')
                lambda_old = XN[r]
                lambda_l = 0.0
                lambda_u = ub
                while (lambda_u - lambda_l) > tol * max(ub, 1.0) / 1e3:
                    lam = (lambda_l + lambda_u) / 2
                    XN[r] = lam
                    Qir = _sum_node_qlen(L, XN, mi, scv, K)
                    g = lam * Z[r] + np.sum(Qir[:, r])
                    if g > N[r]:
                        lambda_u = lam
                    else:
                        lambda_l = lam
                XN[r] = (lambda_l + lambda_u) / 2
                delta = max(delta, abs(XN[r] - lambda_old))
            if delta <= tol:
                break

    QN = _sum_node_qlen(L, XN, mi, scv, K)
    for i in range(M):
        for r in range(R):
            if np.isinf(mi[i]):
                UN[i, r] = XN[r] * L[i, r]
            else:
                UN[i, r] = XN[r] * L[i, r] / mi[i]
            if XN[r] > 0:
                RN[i, r] = QN[i, r] / XN[r]
    return XN, QN, UN, RN, it


def sum_closing(lambda0, scva, L, mi=None, scv=None, N=None, Z=None,
                Kclosed=5000, tol=1e-6, maxiter=10000):
    """
    Closing method for open and mixed non-product-form queueing networks
    (Bolch et al., Sec. 10.1.5), solved with the summation method.

    The external world of each open class is replaced by an additional
    -/G/1 station with service rate mu_inf,r = Ropen*lambda0[r], where
    Ropen is the number of open classes, service SCV equal to the
    interarrival time SCV of the open class, and unit visit ratio. The
    resulting closed network is then solved by sum_closed with a large
    population Kclosed for the open classes. Closed classes are passed
    through unchanged, which makes the method applicable to mixed
    networks.

    Args:
        lambda0: 1xR external arrival rates (0 for closed classes)
        scva: 1xR interarrival time SCVs of the open classes (1 if Poisson)
        L: MxR service demand matrix of the original network, with visit
           ratios of open classes normalized per external arrival
        mi: Mx1 number of servers (np.inf for infinite-server stations)
        scv: MxR service time SCVs (pass 1 for insensitive stations)
        N: 1xR populations, np.inf for open classes
        Z: 1xR think times
        Kclosed: closing population for the open classes (default: 5000)
        tol: convergence tolerance
        maxiter: maximum number of iterations

    Returns:
        (XN, QN, UN, RN, TN, it): 1xR throughputs (open classes approach
        lambda0 from below as Kclosed grows), MxR original-station mean
        queue lengths, utilizations and residence times, 1xR mean
        response times TN=sum(QN)/XN, iteration count.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    lambda0 = np.atleast_1d(np.asarray(lambda0, dtype=float)).ravel()
    if L.shape[0] == 1 and L.shape[1] > 1 and lambda0.size == 1:
        L = L.T
    M, R = L.shape
    scva = np.ones(R) if scva is None else np.atleast_1d(np.asarray(scva, dtype=float)).ravel()
    mi = np.ones(M) if mi is None else np.atleast_1d(np.asarray(mi, dtype=float)).ravel()
    scv = np.ones((M, R)) if scv is None else np.atleast_2d(np.asarray(scv, dtype=float))
    if scv.shape != (M, R):
        scv = scv.reshape(M, R)
    N = np.inf * np.ones(R) if N is None else np.atleast_1d(np.asarray(N, dtype=float)).ravel()
    Z = np.zeros(R) if Z is None else np.atleast_1d(np.asarray(Z, dtype=float)).ravel()

    open_classes = np.where(lambda0 > 0)[0]
    Ropen = len(open_classes)
    if Ropen == 0:
        raise ValueError('sum_closing: no open class, use sum_closed for closed networks.')

    # augment with the closing -/G/1 station, visited by open classes only
    Laug = np.vstack([L, np.zeros((1, R))])
    scvaug = np.vstack([scv, np.ones((1, R))])
    Naug = N.copy()
    for r in open_classes:
        Laug[M, r] = 1 / (Ropen * lambda0[r])
        scvaug[M, r] = scva[r]
        Naug[r] = Kclosed
    miaug = np.concatenate([mi, [1.0]])

    XN, QNa, UNa, RNa, it = sum_closed(Laug, Naug, Z, miaug, scvaug, tol, maxiter)

    QN = QNa[:M, :]
    UN = UNa[:M, :]
    RN = RNa[:M, :]
    TN = np.zeros(R)
    for r in range(R):
        if XN[r] > 0:
            TN[r] = np.sum(QN[:, r]) / XN[r]
    return XN, QN, UN, RN, TN, it


def _sum_node_qlen(L, XN, mi, scv, K):
    """Per-station per-class mean queue lengths Ki_r = fir(lambda_r)."""
    M, R = L.shape
    XN = np.asarray(XN, dtype=float).ravel()
    Qir = np.zeros((M, R))
    for i in range(M):
        if np.isinf(mi[i]):
            Qir[i, :] = XN * L[i, :]  # Type 3, Eq. (9.15)
            continue
        m = mi[i]
        Uir = XN * L[i, :]  # class offered loads lambda_r*e_ir/mu_ir
        Ui = np.sum(Uir)
        if Ui == 0:
            continue
        # per-server utilization; the correction factors keep the node
        # functions finite at rho=1 (Ki(1)<=K)
        rho = min(Ui / m, 1.0)
        ci2 = float(np.sum(Uir * scv[i, :]) / Ui)  # demand-weighted node SCV
        ai = (1 + ci2) / 2
        if K <= m:
            # never more than m jobs at a m-server node: no queueing
            Qir[i, :] = Uir
            continue
        if m == 1:
            if ci2 == 1 or K <= 1:
                # Type 1,2,4 with mi=1, Eq. (9.15)/(9.19)
                Qir[i, :] = Uir / (1 - (K - 1) / K * rho)
            else:
                # -/G/1 FCFS, Eq. (10.88)
                den = 1 - (K - 1 - ai) / (K - 1) * rho
                Qir[i, :] = Uir * (1 + rho * ai / den)
        else:
            Pm = _sum_erlangc(int(m), rho)
            if ci2 == 1:
                # Type 1 with mi>1, Eq. (9.15)/(9.19)
                den = 1 - (K - m - 1) / (K - m) * rho
                Qir[i, :] = Uir + (Uir / m) * Pm / den
            else:
                # -/G/m FCFS, Eq. (10.89)
                den = 1 - (K - m - ai) / (K - m) * rho
                Qir[i, :] = Uir + (Uir / m) * ai * Pm / den
    return Qir


def _sum_erlangc(m, rho):
    """Erlang-C probability of waiting for an M/M/m queue (Eq. 6.28)."""
    if rho >= 1:
        return 1.0
    a = m * rho
    s = 0.0
    term = 1.0  # a^k/k!
    for k in range(m):
        if k > 0:
            term *= a / k
        s += term
    last = term * a / m / (1 - rho)  # a^m/(m!(1-rho))
    return last / (s + last)
