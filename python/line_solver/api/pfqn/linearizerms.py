"""
Multi-server and Mixed Linearizer Algorithms.

Native Python implementations of the multi-server Linearizer algorithm
(Krzesinski/Conway/De Souza-Muntz). The mixed open/closed linearizer
(pfqn_linearizermx) lives in linearizermx.py.

Key functions:
    pfqn_linearizerms: Multi-server Linearizer
    pfqn_conwayms: Conway multi-server Linearizer

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_linearizerms.m
    Conway, 1989, "Fast Approximate Solution of Queueing Networks
    with Multi-Server Chain-Dependent FCFS Queues"
"""

import numpy as np
from typing import Tuple, Optional

from .linearizer import SchedStrategy, pfqn_linearizer, pfqn_gflinearizer, pfqn_egflinearizer
from .mva import pfqn_bs
from .utils import oner

def _get_max_finite_servers(nservers: np.ndarray) -> int:
    """Get maximum server count, excluding infinite values (Delay nodes)."""
    finite_servers = nservers[np.isfinite(nservers)]
    return int(np.max(finite_servers)) if len(finite_servers) > 0 else 1

def _is_fcfs(sched_val) -> bool:
    """Check if scheduling strategy is FCFS (handles multiple enum/string formats)."""
    if sched_val is None:
        return False
    if isinstance(sched_val, str):
        # String: check for 'FCFS' or numeric '0'
        upper_val = sched_val.upper()
        if 'FCFS' in upper_val:
            return True
        # Handle numeric string '0' (FCFS = 0 in IntEnum)
        try:
            return int(sched_val) == 0
        except ValueError:
            return False
    elif hasattr(sched_val, 'name'):
        # IntEnum or Enum: use name attribute
        return sched_val.name.upper() == 'FCFS'
    elif isinstance(sched_val, (int, float)):
        # Integer: FCFS is typically 0 in IntEnum definitions
        return int(sched_val) == 0
    else:
        # Fallback: convert to string and check
        return 'FCFS' in str(sched_val).upper()

def pfqn_linearizerms(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                      nservers: np.ndarray,
                      type_sched: Optional[np.ndarray] = None,
                      tol: float = 1e-8,
                      maxiter: int = 1000,
                      QN0: np.ndarray = None
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Multiserver Linearizer (Krzesinski/Conway/De Souza-Muntz).

    Extends the Linearizer algorithm to handle multi-server stations
    in product-form queueing networks.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        nservers: Number of servers per station (M,)
        type_sched: Scheduling strategy per station (M,), optional (default: PS)
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)

    Returns:
        Tuple of (Q, U, R, C, X, totiter):
            Q: Mean queue lengths (M, R)
            U: Utilization (M, R)
            R: Residence times (M, R)
            C: Cycle times (R,)
            X: System throughput (R,)
            totiter: Total iterations performed

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_linearizerms.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    nservers = np.asarray(nservers, dtype=float).flatten()

    M, R = L.shape

    if type_sched is None:
        type_sched = np.full(M, SchedStrategy.PS)

    if len(Z) == 0:
        Z = np.zeros(R)

    max_servers = _get_max_finite_servers(nservers)

    # Initialize Q, PB, P, Delta
    Q = np.zeros((M, R, 1 + R))
    PB = np.zeros((M, 1 + R))
    P = np.zeros((M, max_servers, 1 + R))
    Delta = np.zeros((M, R, R))
    # Linearizer corrections for the queue-length marginals; without them the
    # marginals stay at population N while the queue lengths are reduced to
    # N-e_s, which breaks Q + sum_j (m-1-j) p_j >= m-1 and lets W fall below the
    # mean service time.
    DeltaP = np.zeros((M, max_servers, R))
    DeltaPB = np.zeros((M, R))

    # see _kb/03-api-layer.md for rationale
    for r in range(R):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()
            if QN0 is None:
                result = pfqn_bs(L, N_1, Z)
            else:
                result = pfqn_bs(L, N_1, Z, tol, maxiter, QN0)   # warm start
            q = result[1]  # QN is at index 1 (XN, QN, UN, RN, it), shape (M, R)
            Q[:, r, s] = q[:, r]  # Extract column r

    for ist in range(M):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()
            pop = np.sum(N_1)
            # Only handle multiserver for finite server counts > 1 (skip Delay nodes with Inf)
            if nservers[ist] > 1 and np.isfinite(nservers[ist]):
                ns = int(nservers[ist])
                if pop == 0:
                    # empty network: the station is idle with probability one
                    P[ist, 1:ns, s] = 0.0
                    PB[ist, s] = 0.0
                    P[ist, 0, s] = 1.0
                    continue
                for j in range(1, ns):
                    P[ist, j, s] = 2 * np.sum(Q[ist, :, s]) / (pop * (pop + 1))

                if pop > nservers[ist] - 1:
                    PB[ist, s] = 2 * np.sum(Q[ist, :, s]) / (pop + 1 - nservers[ist]) / (pop * (pop + 1))
                else:
                    PB[ist, s] = 0.0

                P[ist, 0, s] = 1 - PB[ist, s] - np.sum(P[ist, 1:ns, s])

    totiter = 0

    # Main loop (2 iterations)
    for I in range(2):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()

            # Core iteration
            Q[:, :, s], _, _, P[:, :, s], PB[:, s], iter_count = _core_ms(
                L, M, R, N_1, Z, nservers, Q[:, :, s], P[:, :, s], PB[:, s],
                Delta, DeltaP, DeltaPB, type_sched, tol, maxiter - totiter
            )
            totiter += iter_count

        # Update Delta
        for ist in range(M):
            for r in range(R):
                for s in range(R):
                    Ns = oner(N, s)
                    if N[r] <= 0:
                        Delta[ist, r, s] = 0.0
                    elif Ns[r] > 0:
                        # Python stores full population at index R (unlike MATLAB which uses index 0)
                        Delta[ist, r, s] = Q[ist, r, s] / Ns[r] - Q[ist, r, R] / N[r]
                    else:
                        # Chandy-Neuse 0/0 convention: F_ir(N-e_s) = 0
                        Delta[ist, r, s] = -Q[ist, r, R] / N[r]

        # Update the marginal corrections. Probabilities do not scale with the
        # population, so the analogue of Delta is a plain difference.
        for ist in range(M):
            if nservers[ist] > 1 and np.isfinite(nservers[ist]):
                ns = int(nservers[ist])
                for s in range(R):
                    for j in range(ns):
                        DeltaP[ist, j, s] = P[ist, j, s] - P[ist, j, R]
                    DeltaPB[ist, s] = PB[ist, s] - PB[ist, R]

    # Final Core(N) - Python stores full population at index R
    Q_final, W, X, _, _, iter_count = _core_ms(
        L, M, R, N, Z, nservers, Q[:, :, R], P[:, :, R], PB[:, R],
        Delta, DeltaP, DeltaPB, type_sched, tol, maxiter - totiter
    )
    totiter += iter_count

    # Compute performance metrics
    U = np.zeros((M, R))
    for ist in range(M):
        for r in range(R):
            if np.isinf(nservers[ist]):
                # Delay node: utilization is 0 (infinite servers)
                U[ist, r] = 0.0
            elif nservers[ist] == 1:
                U[ist, r] = X[r] * L[ist, r]
            else:
                U[ist, r] = X[r] * L[ist, r] / nservers[ist]

    Q_out = Q_final
    C = N / X - Z
    R_out = W

    return Q_out, U, R_out, C, X, totiter

def _core_ms(L, M, R, N_1, Z, nservers, Q, P, PB, Delta, DeltaP, DeltaPB, type_sched, tol, maxiter):
    """Core iteration for multiserver linearizer."""
    max_servers = _get_max_finite_servers(nservers)
    iter_count = 0
    hasConverged = False

    while not hasConverged:
        iter_count += 1
        Qlast = Q.copy()

        # Estimate
        Q_1, P_1, PB_1 = _estimate_ms(M, R, N_1, nservers, Q, P, PB, Delta, DeltaP, DeltaPB)

        # Forward MVA
        Q, W, T, P, PB = _forward_mva_ms(L, M, R, N_1, Z, nservers, type_sched, Q_1, P_1, PB_1)

        if np.linalg.norm(Q - Qlast) < tol or iter_count >= maxiter:
            hasConverged = True

    return Q, W, T, P, PB, iter_count

def _estimate_ms(M, R, N_1, nservers, Q, P, PB, Delta, DeltaP, DeltaPB):
    """Estimate populations for linearizer."""
    max_servers = _get_max_finite_servers(nservers)

    # JIT dispatch for large problems
    P_1 = np.zeros((M, max_servers, 1 + R))
    PB_1 = np.zeros((M, 1 + R))
    Q_1 = np.zeros((M, R, 1 + R))

    for ist in range(M):
        # Only handle multiserver for finite server counts > 1 (skip Delay nodes with Inf)
        if nservers[ist] > 1 and np.isfinite(nservers[ist]):
            ns = int(nservers[ist])
            for j in range(ns):
                P_1[ist, j, 0] = P[ist, j]
                for s in range(R):
                    P_1[ist, j, s + 1] = P[ist, j] + DeltaP[ist, j, s]

            PB_1[ist, 0] = PB[ist]
            for s in range(R):
                PB_1[ist, s + 1] = PB[ist] + DeltaPB[ist, s]

        for r in range(R):
            for s in range(R):
                Ns = oner(N_1, s)
                if N_1[r] > 0:
                    # Store at s+1 to match MATLAB's 1+s indexing (1-based s=1:R -> 0-based indices 1:R)
                    Q_1[ist, r, s + 1] = Ns[r] * (Q[ist, r] / N_1[r] + Delta[ist, r, s])
                else:
                    Q_1[ist, r, s + 1] = 0.0

    return Q_1, P_1, PB_1

def _forward_mva_ms(L, M, R, N_1, Z, nservers, type_sched, Q_1, P_1, PB_1):
    """Forward MVA step for multiserver linearizer."""
    max_servers = _get_max_finite_servers(nservers)

    # JIT dispatch for large problems
    W = np.zeros((M, R))
    T = np.zeros(R)
    Q = np.zeros((M, R))
    P = np.zeros((M, max_servers))
    PB = np.zeros(M)

    for ist in range(M):
        for r in range(R):
            # For Delay nodes (Inf servers), W = L (pure delay, no queueing)
            if np.isinf(nservers[ist]):
                W[ist, r] = L[ist, r]
                continue

            W[ist, r] = L[ist, r] / nservers[ist]
            if L[ist, r] == 0:
                continue

            # Use r+1 for 3rd dim to match MATLAB's 1+r indexing (1-based r=1:R -> 0-based indices 1:R)
            if _is_fcfs(type_sched[ist]):
                for s in range(R):
                    W[ist, r] += (L[ist, s] / nservers[ist]) * Q_1[ist, s, r + 1]
            else:
                for s in range(R):
                    W[ist, r] += (L[ist, r] / nservers[ist]) * Q_1[ist, s, r + 1]

            # Partially-idle-server correction. It compensates the 1/m scaling of the
            # arriving job's OWN service, so it carries L[ist,r]/m and no sum over the
            # other classes: at N=e_r the terms must collapse to W = L[ist,r].
            if nservers[ist] > 1 and np.isfinite(nservers[ist]):
                ns = int(nservers[ist])
                for j in range(ns - 1):
                    W[ist, r] += (L[ist, r] / nservers[ist]) * (nservers[ist] - 1 - j) * P_1[ist, j, r + 1]

    for r in range(R):
        denom = Z[r] + np.sum(W[:, r])
        if denom > 0:
            T[r] = N_1[r] / denom
        else:
            T[r] = 0.0

        for ist in range(M):
            Q[ist, r] = T[r] * W[ist, r]

    # Queue-length marginals. The relations
    #   p_j = (A p_{j-1} + d_{j-1}) / j,  pB = (A (pB + p_{ns-1}) + dB) / ns,
    #   p_0 = 1 - pB - sum_j p_j
    # with A = sum_s X_s L_is the mean number of busy servers and d the population
    # corrections, are solved in closed form rather than iterated: as a Jacobi
    # iteration they amplify by A per sweep and diverge once A approaches ns.
    for ist in range(M):
        if nservers[ist] > 1 and np.isfinite(nservers[ist]):
            ns = int(nservers[ist])
            A = 0.0
            d = np.zeros(ns)
            dB = 0.0
            for s in range(R):
                a_s = L[ist, s] * T[s]
                A += a_s
                for j in range(ns):
                    d[j] += a_s * (P_1[ist, j, s + 1] - P_1[ist, j, 0])
                dB += a_s * (PB_1[ist, s + 1] - PB_1[ist, 0])
            if A >= ns:
                raise ValueError(
                    "pfqn_linearizerms: station %d offers %g busy servers out of %d; "
                    "the model is saturated and its queue-length marginals do not exist."
                    % (ist, A, ns))
            # p_j = alpha[j] * p_0 + beta[j]
            alpha = np.zeros(ns)
            beta = np.zeros(ns)
            alpha[0] = 1.0
            for j in range(1, ns):
                alpha[j] = A * alpha[j - 1] / j
                beta[j] = (A * beta[j - 1] + d[j - 1]) / j
            alphaB = A * alpha[ns - 1] / (ns - A)
            betaB = (A * beta[ns - 1] + dB + d[ns - 1]) / (ns - A)
            P[ist, 0] = (1 - np.sum(beta[1:ns]) - betaB) / (1 + np.sum(alpha[1:ns]) + alphaB)
            for j in range(1, ns):
                P[ist, j] = alpha[j] * P[ist, 0] + beta[j]
            PB[ist] = alphaB * P[ist, 0] + betaB

    return Q, W, T, P, PB

def sprod(R: int, n: int) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """
    Initialize state product space iterator.

    Used for enumerating states in multi-class queueing networks.

    Args:
        R: Number of classes
        n: Total population constraint

    Returns:
        Tuple of (s, nvec, SD, D):
            s: Current state index
            nvec: Current state vector
            SD: Upper bounds
            D: Direction vector
    """
    # The first state must be the one sprod_next unranks at 0, or the sweep both
    # repeats a state and misses another: with nvec[0] = n hardcoded here, R=2
    # started at (n,0) while rank 0 is (0,n), so (n,0) was counted twice and (0,n)
    # never. Both ends of the iterator now go through the same unranking.
    SD = np.full(R, n, dtype=int)
    D = np.arange(R, dtype=int)
    return 0, _sprod_unrank(0, n, R), SD, D


def _sprod_unrank(rank: int, n: int, R: int) -> np.ndarray:
    """The rank-th length-R nonnegative vector summing to n, in stars-and-bars order."""
    from math import comb

    nvec = np.zeros(R, dtype=int)
    remaining = n
    for i in range(R - 1):
        for k in range(remaining + 1):
            c = comb(remaining - k + R - i - 2, R - i - 2)
            if rank < c:
                nvec[i] = k
                remaining -= k
                break
            rank -= c
    nvec[R - 1] = remaining
    return nvec

def sprod_next(s: int, SD: np.ndarray, D: np.ndarray) -> Tuple[int, np.ndarray]:
    """
    Get next state in product space.

    The caller drives this as `while s >= 0: ...; s, nvec = sprod_next(s, SD, D)`,
    so the returned index IS the iterator: returning anything but the advanced
    rank makes the loop non-terminating. It used to return the literal 0, so `s`
    oscillated 0 -> 1 -> 0 and the caller re-decoded state 1 forever (reproducible
    at R=2, one station). The rank is therefore captured BEFORE the decoding loop,
    which consumes its own copy.

    Args:
        s: Current state index
        SD: Upper bounds
        D: Direction vector

    Returns:
        Tuple of (s, nvec) for next state, or (-1, nvec) if exhausted
    """
    from math import comb

    R = len(SD)
    n = int(SD[0])

    if s < 0:
        return s, np.zeros(R, dtype=int)

    # Number of length-R nonnegative vectors summing to n (stars and bars).
    # math.comb, not scipy.special.comb: this is an exact index bound, and the
    # float scipy returns loses integrality well before the loop does.
    total = comb(n + R - 1, R - 1)

    s = int(s) + 1
    if s >= total:
        return -1, np.zeros(R, dtype=int)

    return s, _sprod_unrank(s, n, R)

def multinomialln(n: np.ndarray) -> float:
    """Compute log of multinomial coefficient."""
    from scipy.special import gammaln
    return float(gammaln(np.sum(n) + 1) - np.sum(gammaln(np.asarray(n) + 1)))

def pfqn_conwayms(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                   nservers: np.ndarray,
                   type_sched: Optional[np.ndarray] = None,
                   tol: float = 1e-8,
                   maxiter: int = 1000,
                   QN0: np.ndarray = None
                   ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Conway (1989) multiserver Linearizer approximation for FCFS queues.

    Implements the algorithm from Conway (1989), "Fast Approximate Solution
    of Queueing Networks with Multi-Server Chain-Dependent FCFS Queues".

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        nservers: Number of servers per station (M,)
        type_sched: Scheduling strategy per station (M,), optional (default: FCFS)
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)

    Returns:
        Tuple of (Q, U, R, C, X, totiter):
            Q: Mean queue lengths (M, R)
            U: Utilization (M, R)
            R: Residence times (M, R)
            C: Cycle times (R,)
            X: System throughput (R,)
            totiter: Total iterations performed

    References:
        Conway, A. E., "Fast Approximate Solution of Queueing Networks with
        Multi-Server Chain-Dependent FCFS Queues", Performance Evaluation,
        Vol. 8, 1989, pp. 141-159.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    nservers = np.asarray(nservers, dtype=float).flatten()

    M, R = L.shape

    if type_sched is None:
        type_sched = np.full(M, SchedStrategy.FCFS)

    if len(Z) == 0:
        Z = np.zeros(R)
    Z = np.sum(Z.reshape(-1, R) if Z.ndim > 1 else Z.reshape(1, -1), axis=0)

    max_servers = _get_max_finite_servers(nservers)

    # Initialize Q, PB, P, Delta indexed by population reduction
    Q = np.zeros((M, R, 1 + R))
    PB = np.zeros((M, 1 + R))
    P = np.zeros((M, max_servers, 1 + R))
    Delta = np.zeros((M, R, R))

    # Initialize queue lengths
    for ist in range(M):
        for r in range(R):
            for s in range(R + 1):
                N_1 = oner(N, s) if s < R else N.copy()
                if QN0 is None:
                    Q[ist, r, s] = N_1[r] / M
                else:
                    Q[ist, r, s] = QN0[ist, r]   # warm start

    # Initialize probabilities
    for ist in range(M):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()
            pop = np.sum(N_1)
            # Only handle multiserver for finite server counts > 1 (skip Delay nodes with Inf)
            if nservers[ist] > 1 and np.isfinite(nservers[ist]):
                ns = int(nservers[ist])
                for j in range(1, ns):
                    if pop * (pop + 1) > 0:
                        P[ist, j, s] = 2 * np.sum(Q[ist, :, s]) / (pop * (pop + 1))

                denom = (pop + 1 - nservers[ist]) * pop * (pop + 1)
                if denom > 0:
                    PB[ist, s] = 2 * np.sum(Q[ist, :, s]) / denom
                else:
                    PB[ist, s] = 0.0

                P[ist, 0, s] = max(0, 1 - PB[ist, s] - np.sum(P[ist, 1:ns, s]))

    totiter = 0

    # Main loop (2 iterations for linearizer)
    for I in range(2):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()

            # Core iteration
            Q[:, :, s], W_temp, T_temp, P[:, :, s], PB[:, s], iter_count = _conway_core(
                L, M, R, N_1, Z, nservers, Q[:, :, s], P[:, :, s], PB[:, s],
                Delta, type_sched, tol, maxiter - totiter
            )
            totiter += iter_count

        # Update Delta
        for ist in range(M):
            for r in range(R):
                for s in range(R):
                    Ns = oner(N, s)
                    if N[s] > 2 and N[r] > 0 and Ns[r] > 0:
                        # Python stores full population at index R
                        Delta[ist, r, s] = Q[ist, r, s] / Ns[r] - Q[ist, r, R] / N[r]
                    else:
                        Delta[ist, r, s] = 0.0

    # Final Core(N) - Python stores full population at index R.
    # The budget is the FULL maxiter, as in MATLAB pfqn_conwayms.m:110, not the
    # remainder the two priming sweeps left. Charging this call for their
    # iterations truncates the only solve whose result is returned: with the
    # remainder it stopped after 106 iterations against MATLAB's 1035 and landed
    # 0.5% away.
    Q_final, W, X, P_final, PB_final, iter_count = _conway_core(
        L, M, R, N, Z, nservers, Q[:, :, R], P[:, :, R], PB[:, R],
        Delta, type_sched, tol, maxiter
    )
    totiter += iter_count

    # Compute performance metrics
    U = np.zeros((M, R))
    for ist in range(M):
        for r in range(R):
            if np.isinf(nservers[ist]):
                # Delay node: utilization is 0 (infinite servers)
                U[ist, r] = 0.0
            elif nservers[ist] == 1:
                U[ist, r] = X[r] * L[ist, r]
            else:
                U[ist, r] = X[r] * L[ist, r] / nservers[ist]

    Q_out = Q_final
    C = N / np.maximum(X, 1e-12) - Z
    R_out = W

    return Q_out, U, R_out, C, X, totiter

def _conway_core(L, M, R, N_1, Z, nservers, Q, P, PB, Delta, type_sched, tol, maxiter):
    """Core iteration for Conway multiserver linearizer."""
    max_servers = _get_max_finite_servers(nservers)
    hasConverged = False
    iter_count = 0
    # W seeds the first Estimate and is reset per Core call, as in MATLAB Core:
    # carrying it over from the previous population sweep makes the estimated
    # throughputs of sweep s depend on the sweep before it.
    W = L.copy()
    Wlast = None

    while not hasConverged:
        Qlast = Q.copy()

        # Estimate
        Q_1, P_1, PB_1, T_1 = _conway_estimate(M, R, N_1, nservers, Q, P, PB, Delta, W)

        # Forward MVA
        Q, W, T, P, PB = _conway_forward_mva(L, M, R, N_1, Z, nservers, type_sched, Q_1, P_1, PB_1, T_1)

        # W must enter the test: Q alone is satisfied on the FIRST sweep whenever Q
        # cannot move (M=1 seeds Q at its own fixed point), and the residence times
        # returned then are still the seed W=L, so T_1=Q/W is unbounded and the
        # throughput exceeds the station's own service capacity.
        moved = np.inf if Wlast is None else max(np.linalg.norm(Q - Qlast),
                                                 np.linalg.norm(W - Wlast))
        Wlast = W.copy()

        if moved < tol or iter_count >= maxiter:
            hasConverged = True

        iter_count += 1

    return Q, W, T, P, PB, iter_count

def _conway_estimate(M, R, N_1, nservers, Q, P, PB, Delta, W):
    """Estimate populations for Conway linearizer."""
    max_servers = _get_max_finite_servers(nservers)

    # JIT dispatch for large problems
    P_1 = np.zeros((M, max_servers, 1 + R))
    PB_1 = np.zeros((M, 1 + R))
    Q_1 = np.zeros((M, R, 1 + R))
    T_1 = np.zeros((R, 1 + R))

    for ist in range(M):
        # Only handle multiserver for finite server counts > 1 (skip Delay nodes with Inf)
        if nservers[ist] > 1 and np.isfinite(nservers[ist]):
            ns = int(nservers[ist])
            for j in range(ns):
                for s in range(R + 1):
                    P_1[ist, j, s] = P[ist, j]
            for s in range(R + 1):
                PB_1[ist, s] = PB[ist]

        for r in range(R):
            for s in range(R):
                Ns = oner(N_1, s)
                if N_1[r] > 0:
                    Q_1[ist, r, s] = Ns[r] * (Q[ist, r] / N_1[r] + Delta[ist, r, s])
                else:
                    Q_1[ist, r, s] = 0.0

    # T_1 is Little's law over the queueing part of the cycle, sum_i Q_1 / sum_i W,
    # and not the ratio at the FIRST station with a positive residence time: the
    # per-station estimates disagree, so picking one made the answer depend on the
    # station order. The demand matrix carries no order, so a model symmetric under
    # permuting classes and stations together must return equal class throughputs,
    # and with the single-station pick it did not.
    for r in range(R):
        for s in range(R):
            Nr = oner(N_1, r)
            num = 0.0
            den = 0.0
            for ist in range(M):
                if W[ist, s] > 0 and N_1[s] > 0:
                    # Delta is indexed (station, queued class, removed class), as the
                    # Q_1 loop above uses it: here class s queues and class r is removed.
                    num += Nr[s] * (Q[ist, s] / N_1[s] + Delta[ist, s, r])
                    den += W[ist, s]
            if den > 0:
                # a reduced-population throughput cannot be negative; a negative one
                # makes log(F) complex in the XR/XE sums of the forward step
                T_1[s, r] = max(0.0, num / den)

    return Q_1, P_1, PB_1, T_1

def _conway_forward_mva(L, M, R, N_1, Z, nservers, type_sched, Q_1, P_1, PB_1, T_1):
    """Forward MVA step for Conway multiserver linearizer."""
    max_servers = _get_max_finite_servers(nservers)
    W = np.zeros((M, R))
    T = np.zeros(R)
    Q = np.zeros((M, R))
    P = np.zeros((M, max_servers))
    PB = np.zeros(M)
    XR = np.zeros((M, R))
    XE = np.zeros((M, R, R))

    mu = 1.0 / np.maximum(L, 1e-12)

    # Compute F matrix for each class
    F = []
    for r in range(R):
        F_r = np.zeros((M, R))
        for ist in range(M):
            den = np.dot(L[ist, :], T_1[:, r])
            if den > 0:
                for c in range(R):
                    F_r[ist, c] = T_1[c, r] * L[ist, c] / den
        F.append(F_r)

    # Compute XR (multi-server specific)
    for ist in range(M):
        for r in range(R):
            # Only handle multiserver for finite server counts > 1 (skip Delay nodes with Inf)
            if nservers[ist] > 1 and np.isfinite(nservers[ist]):
                ns = int(nservers[ist])
                XR[ist, r] = 0.0
                C_val = 0.0

                s, nvec, SD, D = sprod(R, ns)
                while s >= 0:
                    Nr = oner(N_1, r)
                    if np.all(nvec <= Nr):
                        Ai = np.exp(multinomialln(nvec) + np.dot(nvec, np.log(np.maximum(F[r][ist, :], 1e-300))))
                        C_val += Ai
                        mu_sum = np.dot(mu[ist, :], nvec)
                        if mu_sum > 0:
                            XR[ist, r] += Ai / mu_sum

                    s, nvec = sprod_next(s, SD, D)

                if C_val > 0:
                    XR[ist, r] /= C_val

    # Compute XE
    for ist in range(M):
        for r in range(R):
            # Only handle multiserver for finite server counts > 1 (skip Delay nodes with Inf)
            if nservers[ist] > 1 and np.isfinite(nservers[ist]):
                ns = int(nservers[ist])
                for c in range(R):
                    XE[ist, r, c] = 0.0
                    Cx = 0.0

                    s, nvec, SD, D = sprod(R, ns)
                    while s >= 0:
                        Nr = oner(N_1, r)
                        if np.all(nvec <= Nr) and nvec[c] >= 1:
                            Aix = np.exp(multinomialln(nvec) + np.dot(nvec, np.log(np.maximum(F[r][ist, :], 1e-300))))
                            Cx += Aix
                            mu_sum = np.dot(mu[ist, :], nvec)
                            if mu_sum > 0:
                                XE[ist, r, c] += Aix / mu_sum

                        s, nvec = sprod_next(s, SD, D)

                    if Cx > 0:
                        XE[ist, r, c] /= Cx

    # Compute residence time
    for ist in range(M):
        for r in range(R):
            # Delay node (Inf servers): W = L (pure delay, no queueing)
            if np.isinf(nservers[ist]):
                W[ist, r] = L[ist, r]
            elif nservers[ist] == 1:
                if _is_fcfs(type_sched[ist]):
                    W[ist, r] = L[ist, r]
                    for c in range(R):
                        W[ist, r] += L[ist, c] * Q_1[ist, c, r]
                else:
                    W[ist, r] = L[ist, r]
                    for c in range(R):
                        W[ist, r] += L[ist, r] * Q_1[ist, c, r]
            else:
                # Multiserver (finite servers > 1)
                W[ist, r] = L[ist, r] + PB_1[ist, r] * XR[ist, r]
                for c in range(R):
                    W[ist, r] += XE[ist, r, c] * (Q_1[ist, c, r] - L[ist, c] * T_1[c, r])

    # Compute throughputs and queue lengths
    for r in range(R):
        denom = Z[r] + np.sum(W[:, r])
        if denom > 0:
            T[r] = N_1[r] / denom
        else:
            T[r] = 0.0

        for ist in range(M):
            Q[ist, r] = T[r] * W[ist, r]

    # Queue-length marginals. The relations
    #   p_j = A*p_{j-1}/j,  pB = A*(pB + p_{ms-1})/ms,  p_0 = 1 - pB - sum_j p_j
    # with A = sum_s X_s*L_is the mean number of busy servers are solved in closed
    # form rather than iterated. As a Jacobi iteration they amplify by A per sweep,
    # and since the convergence test watches Q and W but not P the routine returned
    # marginals whose mass had run to 334 behind the p_0 = max(0,1-...) floor.
    # _conway_estimate hands the same marginals to every reduced population, so the
    # population corrections that pfqn_linearizerms carries here are all zero.
    for ist in range(M):
        if nservers[ist] > 1 and np.isfinite(nservers[ist]):
            ms = int(nservers[ist])
            A = 0.0
            for s in range(R):
                A += L[ist, s] * T[s]
            P[ist, :] = 0
            if A >= ms:
                # Saturated: the closed form is singular and its limit is the
                # degenerate marginal, every server busy with probability one.
                # N = m with Z = 0 reaches it exactly, so this is a legal input.
                PB[ist] = 1.0
            else:
                alpha = np.zeros(ms)
                alpha[0] = 1.0
                for j in range(1, ms):
                    alpha[j] = A * alpha[j - 1] / j
                alphaB = A * alpha[ms - 1] / (ms - A)
                P[ist, 0] = 1.0 / (1.0 + np.sum(alpha[1:ms]) + alphaB)
                for j in range(1, ms):
                    P[ist, j] = alpha[j] * P[ist, 0]
                PB[ist] = alphaB * P[ist, 0]

    return Q, W, T, P, PB

__all__ = [
    'pfqn_linearizerms',
    'pfqn_conwayms',
]
