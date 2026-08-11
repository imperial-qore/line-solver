"""
Multichain Convolution Algorithm for class-dependent (cdscaling) / FES Models.

Implements the convolution algorithm of Sauer (1983), Section 5.2,
"Computational Algorithms for State-Dependent Queueing Networks",
ACM TOCS, Vol. 1, No. 1, pp. 67-92.

Computes G(N) = (X_1 * X_2 * ... * X_M)(N) where X_m(n) is the
station factor at population vector n, and * denotes the multivariate
discrete convolution.

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_conv.m
    Original MATLAB: matlab/src/solvers/NC/solver_nc_conv.m
"""

import numpy as np
from typing import Optional, List, Tuple
from scipy.special import gammaln

from ..io.logging import line_warning


def _hashpop(n: np.ndarray, N: np.ndarray) -> int:
    """Map population vector to linear index (0-based)."""
    idx = 0
    R = len(N)
    for r in range(R):
        idx += int(np.prod(N[:r] + 1)) * int(n[r])
    return idx


def _cd_peak_scaling(beta, NK: np.ndarray, K: int) -> float:
    """
    Peak of the class-dependence callable over the reachable population lattice
    0 <= n[r] <= NK[r]. beta returns either a scalar (shared by every class) or
    a length-K vector, so the peak is taken over both the states and the
    classes: utilization is a per-station quantity, so the whole station shares
    one normalizer, as it does for max(lldscaling[ist, :]).
    """
    bmax = 0.0
    n = _pprod_init(NK)
    while n[0] >= 0:
        if int(np.sum(n)) > 0:
            v = np.atleast_1d(np.asarray(beta(n), dtype=float)).flatten()
            v = v[np.isfinite(v)]
            if v.size > 0:
                bmax = max(bmax, float(np.max(v)))
        n = _pprod_next(n, NK)
    return bmax


def _pprod_init(N: np.ndarray) -> np.ndarray:
    """Initialize population vector enumeration."""
    return np.zeros(len(N), dtype=int)


def _pprod_next(n: np.ndarray, N: np.ndarray) -> np.ndarray:
    """Advance to next population vector. Returns n[0]=-1 when done."""
    R = len(N)
    if np.all(n == N):
        return -np.ones(R, dtype=int)
    s = R - 1
    while s >= 0 and n[s] == N[s]:
        n[s] = 0
        s -= 1
    if s >= 0:
        n[s] += 1
    return n


def _fz(Z: np.ndarray, n: np.ndarray) -> float:
    """Delay server unnormalized probability factor.

    F = (Z[0]^n[0] / n[0]!) * ... * (Z[R-1]^n[R-1] / n[R-1]!)
    """
    if np.sum(n) == 0:
        return 1.0
    log_f = 0.0
    for r in range(len(n)):
        if Z[r] > 0:
            log_f += np.log(Z[r]) * n[r]
            log_f -= gammaln(1 + n[r])
        elif n[r] > 0:
            return 0.0
    return np.exp(log_f)


def pfqn_conv(
    L: np.ndarray,
    N: np.ndarray,
    Z: Optional[np.ndarray] = None,
    cdscaling: Optional[List] = None,
) -> Tuple[float, float]:
    """Multichain convolution for closed networks with class-dependent rates.

    Implements the convolution algorithm of Sauer (1983), Section 5.2. For
    class-dependent stations the station factor is built by eq. (40),

        X_m(n) = (u_km / mu_km(n)) * X_m(n - e_k),   X_m(0) = 1

    where mu_km(n) = (n_k/|n|) * beta_{m,k}(n) and beta is the DIMENSIONLESS
    class-dependent demand scaling supplied by the class-dependence handle:
    cdscaling[m] is a callable of the per-class population vector n at station m
    returning either a scalar (chain-independent) or a length-R array of
    per-class rates. Any saturation/cutoff is applied inside the handle.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (1 x R), must be finite
        Z: Think time vector (1 x R)
        cdscaling: List of M class-dependence callables beta_m(n); a None entry
            denotes a load-independent station

    Returns:
        (G, lG) where G is the normalizing constant and lG = log(G)
    """
    L = np.atleast_2d(L)
    M, R = L.shape
    N = np.asarray(N, dtype=int).flatten()

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z).flatten()

    if cdscaling is None:
        cdscaling = [None] * M

    # Total state space size
    state_space_size = int(np.prod(N + 1))

    # Identify which stations carry a class-dependence function
    is_cd = [ist < len(cdscaling) and cdscaling[ist] is not None for ist in range(M)]

    # --- Precompute X_m(n) tables for class-dependent stations ---
    Xm = [None] * M
    for ist in range(M):
        if is_cd[ist]:
            Xm[ist] = np.zeros(state_space_size)
            Xm[ist][0] = 1.0  # X_m(0) = 1

            n = _pprod_init(N)
            while n[0] >= 0:
                idx = _hashpop(n, N)
                if np.sum(n) == 0:
                    Xm[ist][idx] = 1.0
                else:
                    for r in range(R):
                        if n[r] > 0:
                            # see _kb/03-api-layer.md for rationale
                            bval = np.atleast_1d(
                                np.asarray(cdscaling[ist](n), dtype=float))
                            beta = float(bval[r]) if bval.size > 1 else float(bval[0])

                            # see _kb/03-api-layer.md for rationale
                            tot = float(np.sum(n))
                            nr = float(n[r])
                            n[r] -= 1
                            idx_prev = _hashpop(n, N)
                            n[r] += 1
                            if beta > 0:
                                Xm[ist][idx] = (tot / nr) * (L[ist, r] / beta) * Xm[ist][idx_prev]
                            else:
                                Xm[ist][idx] = 0.0
                            break
                n = _pprod_next(n, N)

    # --- Convolution ---
    G_curr = np.zeros(state_space_size)

    # Initialize G_0(n) = F_Z(n)
    n = _pprod_init(N)
    while n[0] >= 0:
        idx = _hashpop(n, N)
        G_curr[idx] = _fz(Z, n)
        n = _pprod_next(n, N)

    # Convolve one station at a time
    for ist in range(M):
        if is_cd[ist]:
            # class-dependent station: direct convolution sum
            G_old = G_curr.copy()
            G_curr = np.zeros(state_space_size)

            n = _pprod_init(N)
            while n[0] >= 0:
                idx_n = _hashpop(n, N)
                conv_sum = 0.0

                # Inner loop: enumerate all i from 0 to n
                i = _pprod_init(n.copy())
                while i[0] >= 0:
                    idx_i = _hashpop(i, N)
                    nmi = n - i
                    idx_nmi = _hashpop(nmi, N)
                    conv_sum += Xm[ist][idx_i] * G_old[idx_nmi]
                    i = _pprod_next(i, n.copy())

                G_curr[idx_n] = conv_sum
                n = _pprod_next(n, N)
        else:
            # LI station: efficient recurrence
            n = _pprod_init(N)
            while n[0] >= 0:
                idx_n = _hashpop(n, N)
                for r in range(R):
                    if n[r] >= 1:
                        n[r] -= 1
                        idx_n1r = _hashpop(n, N)
                        n[r] += 1
                        G_curr[idx_n] += L[ist, r] * G_curr[idx_n1r]
                n = _pprod_next(n, N)

    G = G_curr[-1]  # G(N) is at the last index
    lG = np.log(G) if G > 0 else -np.inf
    return G, lG


def solver_nc_conv(sn, options=None):
    """NC solver using convolution for class-dependent/FES models.

    Port of MATLAB solver_nc_conv.m.

    Args:
        sn: NetworkStruct
        options: Solver options

    Returns:
        SolverNCReturn with Q, U, R, T, X, lG, etc.
    """
    from ..solvers.nc.handler import SolverNCReturn
    import time
    start_time = time.time()

    M = sn.nstations
    K = sn.nclasses
    NK = sn.njobs.flatten().astype(int)
    nservers = sn.nservers.flatten()

    # Compute visits and service times
    V = sum(v for v in sn.visits.values() if v is not None)
    ST = np.zeros((M, K))
    rates = sn.rates
    for i in range(M):
        for k in range(K):
            if rates[i, k] > 0:
                ST[i, k] = 1.0 / rates[i, k]

    # Demands
    Ldemand = V * ST

    # Separate delay and queue stations
    is_delay = np.isinf(nservers)
    delay_idx = np.where(is_delay)[0]
    queue_idx = np.where(~is_delay)[0]
    n_queues = len(queue_idx)

    # Z: total delay demand per class
    Z_conv = np.zeros(K)
    for ist in delay_idx:
        Z_conv += Ldemand[ist, :]

    # L_conv: demands for queue stations only
    L_conv = Ldemand[queue_idx, :]

    # see _kb/03-api-layer.md for rationale
    cdscaling_conv = [None] * n_queues
    cd = getattr(sn, 'cdscaling', None)
    if cd is not None and len(cd) > 0:
        for qi in range(n_queues):
            ist = queue_idx[qi]
            if ist < len(cd) and cd[ist] is not None:
                cdscaling_conv[qi] = cd[ist]
    # Fold joint-dependence handles (sn.jdscaling, non-product-form eta_i) into
    # the same per-station handle used by the convolution recursion; cd and jd
    # are evaluated identically. When only one is present the product reproduces
    # the single-mechanism case (the other is treated as absent).
    jd = getattr(sn, 'jdscaling', None)
    if jd is not None and len(jd) > 0:
        for qi in range(n_queues):
            ist = queue_idx[qi]
            if ist < len(jd) and jd[ist] is not None:
                jd_h = jd[ist]
                cd_h = cdscaling_conv[qi]
                if cd_h is None:
                    cdscaling_conv[qi] = jd_h
                else:
                    cdscaling_conv[qi] = (lambda cf, jf: (lambda n: np.asarray(cf(n), dtype=float) * np.asarray(jf(n), dtype=float)))(cd_h, jd_h)

    # Compute G(N)
    G_N, lG = pfqn_conv(L_conv, NK, Z_conv, cdscaling_conv)

    # Compute G(N - e_k) for each class -> throughput
    XN = np.zeros(K)
    for k in range(K):
        if NK[k] > 0:
            NK_minus = NK.copy()
            NK_minus[k] -= 1
            G_Nk, _ = pfqn_conv(L_conv, NK_minus, Z_conv, cdscaling_conv)
            if G_N > 0:
                XN[k] = G_Nk / G_N

    # Per-station throughput
    TN = V * XN[np.newaxis, :]

    # Queue lengths
    QN = np.zeros((M, K))

    # Delay stations: Q = L * X
    for ist in delay_idx:
        for k in range(K):
            QN[ist, k] = Ldemand[ist, k] * XN[k]

    # Queue stations: use marginal distribution
    state_space_size = int(np.prod(NK + 1))

    for qi in range(n_queues):
        ist = queue_idx[qi]

        # Build X_m(n) for this station
        Xm_qi = np.zeros(state_space_size)
        Xm_qi[0] = 1.0
        is_cd_station = cdscaling_conv[qi] is not None

        n = _pprod_init(NK)
        while n[0] >= 0:
            idx = _hashpop(n, NK)
            if np.sum(n) > 0:
                if is_cd_station:
                    for r in range(K):
                        if n[r] > 0:
                            # see _kb/03-api-layer.md for rationale
                            bval = np.atleast_1d(
                                np.asarray(cdscaling_conv[qi](n), dtype=float))
                            beta = float(bval[r]) if bval.size > 1 else float(bval[0])
                            # X_m(n) = (|n|/n_r) * (L/beta) * X_m(n-e_r); at
                            # beta=1 this is the LI multinomial recurrence.
                            tot = float(np.sum(n))
                            nr = float(n[r])
                            n[r] -= 1
                            idx_prev = _hashpop(n, NK)
                            n[r] += 1
                            if beta > 0:
                                Xm_qi[idx] = (tot / nr) * (L_conv[qi, r] / beta) * Xm_qi[idx_prev]
                            break
                else:
                    for r in range(K):
                        if n[r] > 0:
                            n[r] -= 1
                            idx_prev = _hashpop(n, NK)
                            n[r] += 1
                            Xm_qi[idx] += L_conv[qi, r] * Xm_qi[idx_prev]
            n = _pprod_next(n, NK)

        # Build complement: all stations except qi
        L_comp = np.delete(L_conv, qi, axis=0)
        cd_comp = [cdscaling_conv[j] for j in range(n_queues) if j != qi]

        # Compute Q_m_k using marginal
        n = _pprod_init(NK)
        while n[0] >= 0:
            if np.any(n > 0):
                idx = _hashpop(n, NK)
                nmi = NK - n
                if np.all(nmi >= 0):
                    G_comp, _ = pfqn_conv(L_comp, nmi, Z_conv, cd_comp)
                    if G_N > 0:
                        prob = Xm_qi[idx] * G_comp / G_N
                        for k in range(K):
                            QN[ist, k] += n[k] * prob
            n = _pprod_next(n, NK)

    # Remaining metrics
    RN = np.zeros((M, K))
    with np.errstate(divide='ignore', invalid='ignore'):
        RN = np.where(TN > 0, QN / TN, 0.0)

    UN = TN * ST

    # see _kb/03-api-layer.md for rationale. Effective peak = product of the
    # class- and joint-dependence peaks declared at the station (missing = 1);
    # a jd-only station has NaN in cdscalingpeak, so guard each contribution.
    cd = getattr(sn, 'cdscaling', None)
    jd = getattr(sn, 'jdscaling', None)
    for qi in range(n_queues):
        if cdscaling_conv[qi] is None:
            continue
        ist = queue_idx[qi]
        has_cd_st = cd is not None and ist < len(cd) and cd[ist] is not None
        has_jd_st = jd is not None and ist < len(jd) and jd[ist] is not None
        for r in range(K):
            bmax = 1.0
            has_peak = False
            if has_cd_st:
                bmax *= sn.cdscalingpeak[ist, r]; has_peak = True
            if has_jd_st:
                bmax *= sn.jdscalingpeak[ist, r]; has_peak = True
            if has_peak and bmax > 0:
                UN[ist, r] = UN[ist, r] / bmax

    CN = np.where(XN > 0, NK / XN, 0.0)

    runtime = time.time() - start_time

    result = SolverNCReturn(
        Q=QN,
        U=UN,
        R=RN,
        T=TN,
        nchains=sn.nchains if hasattr(sn, 'nchains') else K,
        X=XN.reshape(1, -1),
        lG=lG,
        STeff=np.zeros((M, K)),
        it=1,
        runtime=runtime,
        method='conv',
    )
    return result
