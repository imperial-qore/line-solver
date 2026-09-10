"""
Robust Queueing Network Analyzer (RQNA) traffic-variability equations.

Native Python implementation of the index-of-dispersion-for-counts (IDC) based
traffic equations of W. Whitt and W. You (2018), "A Robust Queueing Network
Analyzer Based on Indices of Dispersion", INFORMS J. on Computing.

Provides:
    npfqn_rqna_weight: canonical RBM correlation weight function w*(t)
    npfqn_traffic_idc: limiting variability equations + time-dependent IDC
                       equations with alpha/beta correction terms

References:
    MATLAB: matlab/src/api/npfqn/npfqn_rqna_weight.m
            matlab/src/api/npfqn/npfqn_traffic_idc.m
"""

import numpy as np
from math import sqrt, pi, exp, erfc


def npfqn_rqna_weight(t):
    """
    Canonical RBM correlation weight function w*(t) used by the RQNA.

    w*(t) = 1 - (1 - c*(t))/(2 t), where c*(t) is the correlation function of
    the stationary version of canonical reflected Brownian motion (drift -1,
    diffusion coefficient 1),

        c*(t) = 2(1 - 2t - t^2) Phi^c(sqrt(t)) + 2 sqrt(t) phi(sqrt(t)) (1 + t),

    with Phi^c the standard-normal complementary cdf and phi its density. The
    weight is monotonically increasing with w*(0)=0 and w*(Inf)=1.

    Reference: Whitt and You (2018), eqs. (24)-(25).

    Args:
        t: scalar or array of nonnegative time arguments

    Returns:
        Weight(s) w*(t) in [0,1], same shape as t (scalar in -> float out).
    """
    scalar_in = np.isscalar(t)
    ta = np.atleast_1d(np.asarray(t, dtype=np.float64))
    w = np.zeros(ta.shape)
    flat = ta.ravel()
    wf = w.ravel()
    for idx in range(flat.size):
        ti = flat[idx]
        if ti <= 0:
            wf[idx] = 0.0
            continue
        if not np.isfinite(ti):
            wf[idx] = 1.0
            continue
        st = sqrt(ti)
        phic = 0.5 * erfc(st / sqrt(2.0))       # complementary standard-normal cdf
        phi = exp(-ti / 2.0) / sqrt(2.0 * pi)   # standard-normal density at sqrt(t)
        cstar = 2.0 * (1.0 - 2.0 * ti - ti ** 2) * phic + 2.0 * st * phi * (1.0 + ti)
        if ti < 1e-6:
            # limit w*(t) -> 0 as t -> 0; series avoids catastrophic cancellation
            val = 0.0
        else:
            val = 1.0 - (1.0 - cstar) / (2.0 * ti)
        if val < 0.0:
            val = 0.0
        if val > 1.0:
            val = 1.0
        wf[idx] = val
    if scalar_in:
        return float(wf[0])
    return w


class _TrafficIdcContext:
    """Context returned by npfqn_traffic_idc holding the limiting variability
    solution and a callable IaFun(t) that returns the K-vector of total arrival
    IDCs I_{a,i}(t) for all internal flows at time t."""

    def __init__(self):
        pass

    def IaFun(self, t):
        return _local_idc_at(self, t)


def npfqn_traffic_idc(lambda0, P, c2a0, a0IdcFun, mu, cs2, sIdcFun, corrections=None):
    """
    Traffic variability equations for the RQNA. Assembles and solves:
      - the limiting variability equations (eq. 42/44) for the asymptotic
        total-arrival variability parameters c2_{a,i} = I_{a,i}(Inf);
      - a solver (ctx.IaFun) for the time-dependent IDC equations (eq. 40/43),
        returning I_{a,i}(t) for all internal arrival flows, using the default
        correction terms alpha_{i,j} (eq. 34) and beta_i (eqs. 38-39) and tuning
        function h(rho)=rho^2.

    Models a single-class open network of K single-server FCFS queues with
    Markovian routing P (P[i,j]=p_{i,j}).

    Args:
        lambda0: (K,) external arrival rate into each queue
        P: (K,K) routing matrix among queues
        c2a0: (K,) asymptotic IDC (SCV) of each external arrival process
        a0IdcFun: callable a0IdcFun(t) -> (K,) external arrival IDC I_{a,0,i}(t)
        mu: (K,) service rate at each queue
        cs2: (K,) service SCV c2_{s,i}
        sIdcFun: callable sIdcFun(t) -> (K,) service IDC I_{s,i}(t)
        corrections: optional dict with 'alpha'/'beta' bool toggles

    Returns:
        _TrafficIdcContext with lambda, rho, Xi, c2a, c2d, c2aij, c2x fields
        and method IaFun(t).
    """
    mu = np.asarray(mu, dtype=np.float64).ravel()
    K = mu.size
    lambda0 = np.asarray(lambda0, dtype=np.float64).ravel()
    cs2 = np.asarray(cs2, dtype=np.float64).ravel()
    c2a0 = np.asarray(c2a0, dtype=np.float64).ravel()
    P = np.asarray(P, dtype=np.float64).reshape(K, K)

    if corrections is None:
        corrections = {}
    useAlpha = corrections.get('alpha', True)
    useBeta = corrections.get('beta', True)

    # ----- traffic rate equations (eq. 20-21) -----
    Xi = np.linalg.inv(np.eye(K) - P.T)      # fundamental matrix (I-P')^{-1}
    lam = Xi @ lambda0                       # total arrival rate at each queue
    rho = lam / mu
    lam_ji = (lam.reshape(K, 1) * np.ones((1, K))) * P   # lam_ji[j,i]=lambda_j p_{j,i}

    # ----- correction terms (asymptotic, w*(Inf)=1) -----
    # alpha: c2alpha_{i,j} = 2 Xi_{i,j} p_{i,j} (1-p_{i,j})
    c2alpha = 2.0 * Xi * P * (1.0 - P)
    if not useAlpha:
        c2alpha = np.zeros((K, K))

    # beta: zeta_{j,i;k,i} from eq. (39), then c2beta_i = (2/lambda_i) sum_{j<k} zeta
    Sigma = []
    for l in range(K):
        pl = P[l, :]
        Sl = -(np.outer(pl, pl)) * lam[l]
        np.fill_diagonal(Sl, pl * (1.0 - pl) * lam[l])
        Sigma.append(Sl)
    Amat = np.diag(c2a0 * lambda0)           # external arrival Brownian variance-rate
    for l in range(K):
        Amat = Amat + Sigma[l]
    zetaAll = [np.zeros((K, K)) for _ in range(K)]
    c2beta = np.zeros(K)
    if useBeta:
        for i in range(K):
            nu = (P[:, i].reshape(K, 1) * np.ones((1, K))) * Xi   # nu[l,:]=p_{l,i} Xi[l,:]
            Z = nu @ Amat @ nu.T
            for j in range(K):
                for k in range(K):
                    Z[j, k] += nu[k, :] @ Sigma[j][:, i] + nu[j, :] @ Sigma[k][:, i]
            zetaAll[i] = Z
            s = 0.0
            for j in range(K):
                for k in range(j + 1, K):
                    s += Z[j, k]
            if lam[i] > 0:
                c2beta[i] = (2.0 / lam[i]) * s

    # ----- limiting variability equations (eq. 44): (E - Minf) c = binf -----
    Na = K
    Naij = K * K
    Nd = K
    N = Na + Naij + Nd

    def ia(i):
        return i

    def iaij(i, j):
        return Na + i * K + j

    def idd(i):
        return Na + Naij + i

    Minf = np.zeros((N, N))
    binf = np.zeros(N)
    for i in range(K):
        if lam[i] > 0:
            for j in range(K):
                Minf[ia(i), iaij(j, i)] = lam_ji[j, i] / lam[i]
            binf[ia(i)] = (lambda0[i] / lam[i]) * c2a0[i] + c2beta[i]
        for j in range(K):
            Minf[iaij(i, j), idd(i)] = P[i, j]
            binf[iaij(i, j)] = (1.0 - P[i, j]) + c2alpha[i, j]
        Minf[idd(i), ia(i)] = 1.0
    csol = np.linalg.solve(np.eye(N) - Minf, binf)
    c2a = csol[0:K]
    c2aij = csol[Na:Na + Naij].reshape(K, K)     # row-major -> c2aij[i,j]
    c2d = csol[Na + Naij:]
    c2x = c2a + cs2                              # c2_{x,i} = I_{a,i}(Inf)+c2_{s,i}

    ctx = _TrafficIdcContext()
    ctx.K = K
    ctx.lam = lam
    ctx._lambda = lam
    ctx.lambda0 = lambda0
    ctx.mu = mu
    ctx.rho = rho
    ctx.Xi = Xi
    ctx.P = P
    ctx.cs2 = cs2
    ctx.c2a = c2a
    ctx.c2aij = c2aij
    ctx.c2d = c2d
    ctx.c2x = c2x
    ctx.c2a0 = c2a0
    ctx.lam_ji = lam_ji
    ctx.zetaAll = zetaAll
    ctx.c2alpha = c2alpha
    ctx.a0IdcFun = a0IdcFun
    ctx.sIdcFun = sIdcFun
    return ctx


def _local_idc_at(ctx, t):
    """Solve the time-dependent IDC equations (43) at a single time t; return
    the (K,) vector of total arrival IDCs I_{a,i}(t)."""
    K = ctx.K
    P = ctx.P
    lam = ctx.lam
    rho = ctx.rho
    c2x = ctx.c2x
    lam_ji = ctx.lam_ji
    h = rho ** 2                              # tuning function h(rho)=rho^2

    # departure weights w_i(t)=w*((1-rho_i)^2 lambda_i t /(h_i c2x_i))
    warg = np.zeros(K)
    for i in range(K):
        if h[i] > 0 and c2x[i] > 0:
            warg[i] = (1.0 - rho[i]) ** 2 * lam[i] * t / (h[i] * c2x[i])
        else:
            warg[i] = np.inf
    w = npfqn_rqna_weight(warg)
    w = np.atleast_1d(w).astype(np.float64)

    Ia0 = np.asarray(ctx.a0IdcFun(t), dtype=np.float64).ravel()
    Is = np.asarray(ctx.sIdcFun(rho * t), dtype=np.float64).ravel()   # service IDC at rho*t

    # time-dependent correction terms
    alpha_t = ctx.c2alpha * (w.reshape(K, 1) * np.ones((1, K)))  # row i scaled by w[i]
    beta_t = np.zeros(K)
    for i in range(K):
        Z = ctx.zetaAll[i]
        wj = np.zeros(K)
        for j in range(K):
            if h[j] > 0 and c2x[j] > 0 and P[j, i] > 0:
                aj = (1.0 - rho[j]) ** 2 * P[j, i] * lam[j] * t / (h[j] * c2x[j])
                wj[j] = npfqn_rqna_weight(aj)
        s = 0.0
        for j in range(K):
            for k in range(K):
                if j != k:
                    s += Z[j, k] * wj[j]
        if lam[i] > 0:
            beta_t[i] = s / lam[i]

    # assemble (E - M(t)) I = b(t)
    Na = K
    Naij = K * K
    Nd = K
    N = Na + Naij + Nd

    def ia(i):
        return i

    def iaij(i, j):
        return Na + i * K + j

    def idd(i):
        return Na + Naij + i

    M = np.zeros((N, N))
    b = np.zeros(N)
    for i in range(K):
        if lam[i] > 0:
            for j in range(K):
                M[ia(i), iaij(j, i)] = lam_ji[j, i] / lam[i]
            b[ia(i)] = (ctx.lambda0[i] / lam[i]) * Ia0[i] + beta_t[i]
        for j in range(K):
            M[iaij(i, j), idd(i)] = P[i, j]
            b[iaij(i, j)] = (1.0 - P[i, j]) + alpha_t[i, j]
        M[idd(i), ia(i)] = w[i]
        b[idd(i)] = (1.0 - w[i]) * Is[i]
    sol = np.linalg.solve(np.eye(N) - M, b)
    return sol[0:K]
