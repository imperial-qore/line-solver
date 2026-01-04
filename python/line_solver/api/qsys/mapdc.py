"""
MAP/D/c queue analysis using Non-Skip-Free (NSF) Markov chain methods.

Implements the Q-MAM algorithm for MAP/D/c queues with:
- Markovian Arrival Process (MAP) arrivals
- Deterministic service times
- c parallel FCFS servers

References:
    Van Houdt, B. "Q-MAM: A Toolbox for Numerical Analysis of MMAP[K]/G[K]/1
    Queues." ACM Transactions on Mathematical Software, 2011.
"""

import numpy as np
from scipy import linalg
from typing import Dict, Tuple, Optional
from dataclasses import dataclass

# Import CTMC solver from mc module
from ..mc import ctmc_solve


@dataclass
class MapDcResult:
    """Result container for MAP/D/c queue analysis."""
    mean_queue_length: float
    mean_waiting_time: float
    mean_sojourn_time: float
    utilization: float
    queue_length_dist: np.ndarray
    waiting_time_dist: np.ndarray
    analyzer: str


def qsys_mapdc(
    D0: np.ndarray,
    D1: np.ndarray,
    s: float,
    c: int,
    max_num_comp: int = 1000,
    num_steps: int = 1,
    verbose: int = 0
) -> Dict:
    """
    Analyze MAP/D/c queue (MAP arrivals, deterministic service, c servers).

    Uses Non-Skip-Free (NSF) Markov chain analysis embedding at deterministic
    service intervals. Multiple arrivals can occur per interval.

    Args:
        D0: MAP hidden transition matrix (n x n).
        D1: MAP arrival transition matrix (n x n).
        s: Deterministic service time (positive scalar).
        c: Number of servers.
        max_num_comp: Maximum number of queue length components (default 1000).
        num_steps: Number of waiting time distribution points per interval (default 1).
        verbose: Verbosity level (default 0).

    Returns:
        dict: Performance metrics including:
            - mean_queue_length: Mean number of customers in system
            - mean_waiting_time: Mean waiting time in queue
            - mean_sojourn_time: Mean sojourn time (waiting + service)
            - utilization: Server utilization (per server)
            - queue_length_dist: Queue length distribution P(Q=n)
            - waiting_time_dist: Waiting time CDF at discrete points
            - analyzer: Analyzer identifier
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    m = D0.shape[0]

    # Validate inputs
    assert D0.shape == (m, m) and D1.shape == (m, m), "D0 and D1 must be square matrices"
    assert s > 1e-14, "Service time s must be strictly positive"
    assert c >= 1, "Number of servers c must be at least 1"

    # Uniformization constant
    lam = np.max(np.abs(np.diag(D0)))
    P0 = D0 / lam + np.eye(m)
    P1 = D1 / lam

    # Test load
    Q = D0 + D1
    theta = ctmc_solve(Q)
    lambda_arr = theta @ D1 @ np.ones(m)
    epsilon = 1e-12

    load = lambda_arr * s / c
    if load >= 1 - epsilon:
        raise ValueError(f"The load {load} of the system exceeds one")

    # Compute NSF blocks P(k,s) - probability of k arrivals in time s
    # P(0,s) = exp(D0*s)
    P0s = linalg.expm(D0 * s)
    Ptot = P0s.copy()

    # Poisson terms for uniformization
    h = _compute_poisson_terms(lam * s, epsilon)
    Pterms = len(h)

    # Initialize K matrices
    Kold = [np.eye(m) if j == 0 else np.linalg.matrix_power(P0, j) for j in range(Pterms)]

    # Compute P(k,s) for k = 1, 2, ...
    Ps = []
    k = 1
    prob_cum = np.min(np.sum(Ptot, axis=1))

    while prob_cum < 1 - epsilon:
        K = [P1 @ Kold[0]]
        htemp = h[k-1] if k <= len(h) else _poisson_pdf(k, lam * s)
        Psk = htemp * K[0]

        for j in range(1, Pterms):
            K.append(P0 @ K[j-1] + P1 @ Kold[j])
            htemp = htemp * lam * s / (k + j)
            Psk = Psk + htemp * K[j]

        Ps.append(Psk)
        Ptot = Ptot + Psk
        Kold = K
        prob_cum = np.min(np.sum(Ptot, axis=1))
        k += 1

    # Build A matrix: A = [P0s Ps{1} Ps{2} ... ]
    num_blocks = 1 + len(Ps)
    A = np.zeros((m, m * num_blocks))
    A[:, :m] = P0s
    for i, Psk in enumerate(Ps):
        A[:, (i+1)*m:(i+2)*m] = Psk

    # Compute G using NSF-GHT algorithm
    G = _nsf_ght(A, c, max_num_it=10000, verbose=verbose)

    # Compute stationary distribution pi
    pi = _nsf_pi(A, G, c, max_num_comp=max_num_comp, verbose=verbose)

    # Queue length distribution
    num_levels = pi.shape[0] // m
    ql = np.zeros(num_levels)
    for i in range(num_levels):
        ql[i] = np.sum(pi[i*m:(i+1)*m])

    # Mean queue length
    mean_ql = np.sum(np.arange(len(ql)) * ql)

    # Waiting time distribution W(t) for t = {0, s, 2s, ...}
    w = []
    w0 = np.sum(pi[:min(c*m, len(pi))])
    w.append(w0)

    wt_accum = w0
    i = 2
    while wt_accum < 1 - 1e-10 and i * m * c < len(pi):
        start_idx = (i-1) * m * c
        end_idx = min(i * m * c, len(pi))
        w_level = w[-1] + np.sum(pi[start_idx:end_idx])
        w.append(w_level)
        wt_accum = w_level
        i += 1

    w = np.array(w)

    # Mean waiting time from CDF using survival function integration
    step_size = s / num_steps
    mean_wt = np.sum(1.0 - w[:-1]) * step_size if len(w) > 1 else 0.0

    # Mean sojourn time
    mean_st = mean_wt + s

    # Utilization
    rho = lambda_arr * s / c

    return {
        'mean_queue_length': float(mean_ql),
        'mean_waiting_time': float(mean_wt),
        'mean_sojourn_time': float(mean_st),
        'utilization': float(rho),
        'queue_length_dist': ql,
        'waiting_time_dist': w,
        'analyzer': f'Q-MAM:MAP/D/{c}'
    }


def qsys_mapd1(
    D0: np.ndarray,
    D1: np.ndarray,
    s: float,
    max_num_comp: int = 1000,
    num_steps: int = 1
) -> Dict:
    """
    Analyze MAP/D/1 queue (single-server convenience function).

    Args:
        D0: MAP hidden transition matrix.
        D1: MAP arrival transition matrix.
        s: Deterministic service time.
        max_num_comp: Maximum number of queue length components.
        num_steps: Number of waiting time points per interval.

    Returns:
        dict: Performance metrics (see qsys_mapdc).
    """
    return qsys_mapdc(D0, D1, s, 1, max_num_comp, num_steps)


def _compute_poisson_terms(lambda_s: float, epsilon: float) -> np.ndarray:
    """Compute Poisson probability terms for uniformization."""
    terms = []
    h = np.exp(-lambda_s) * lambda_s
    sum_h = h + np.exp(-lambda_s)
    terms.append(h)

    while sum_h < 1 - epsilon:
        next_h = terms[-1] * lambda_s / (len(terms) + 1)
        terms.append(next_h)
        sum_h += next_h

    return np.array(terms)


def _poisson_pdf(k: int, lam: float) -> float:
    """Compute Poisson PDF P(X=k) for parameter lambda."""
    result = np.exp(-lam)
    for i in range(1, k + 1):
        result *= lam / i
    return result


def _nsf_ght(
    A: np.ndarray,
    N: int,
    max_num_it: int = 10000,
    verbose: int = 0
) -> np.ndarray:
    """
    Non-Skip-Free Gail-Hantler-Taylor algorithm for computing G matrix.

    Computes the minimal nonnegative solution to:
    G = C0 + C1*G + C2*G^2 + ... + C_max*G^max

    Args:
        A: Block matrix [A0 A1 A2 ... A_max] with m rows
        N: Skip factor (e.g., c for MAP/D/c)
        max_num_it: Maximum iterations
        verbose: Verbosity level

    Returns:
        G matrix of dimension (m*N) x (m*N)
    """
    m = A.shape[0]
    K = A.shape[1] // m - 1

    # Initialize G = A[:, :m*N]
    G = A[:, :min(m*N, A.shape[1])].copy()
    if G.shape[1] < m * N:
        G_padded = np.zeros((m, m * N))
        G_padded[:, :G.shape[1]] = G
        G = G_padded

    check = 1.0
    numit = 0

    while check > 1e-14 and numit < max_num_it:
        numit += 1
        Gold = G.copy()
        temp = G.copy()

        # G = A[:, :m*N] + A[:, N*m:(N+1)*m] @ G
        AN = A[:, N*m:(N+1)*m] if (N+1)*m <= A.shape[1] else np.zeros((m, m))
        G = A[:, :m*N] + AN @ G

        # For j = N+1 to K
        for j in range(N + 1, K + 1):
            # temp = [zeros(m) temp[:, :(N-1)*m]] + temp[:, (N-1)*m:] @ Gold
            new_temp = np.zeros((m, m * N))
            new_temp[:, m:] = temp[:, :(N-1)*m]
            new_temp += temp[:, (N-1)*m:] @ Gold

            temp = new_temp

            # G = G + A[:, j*m:(j+1)*m] @ temp
            if (j+1)*m <= A.shape[1]:
                Aj = A[:, j*m:(j+1)*m]
                G = G + Aj @ temp

        check = np.max(np.abs(G - Gold))
        if verbose > 0 and numit % verbose == 0:
            print(f"Check after {numit} iterations: {check}")

    if numit == max_num_it and check > 1e-14:
        print(f"Warning: NSF_GHT - Maximum iterations {numit} reached")

    # Expand to full (m*N) x (m*N) if needed
    if N > 1:
        full_G = np.zeros((m * N, m * N))
        full_G[:m, :] = G

        # Compute remaining block rows
        for j in range(2, N + 1):
            prev_row = full_G[(j-2)*m:(j-1)*m, :]
            shifted = np.zeros((m, m * N))
            shifted[:, m:] = prev_row[:, :(N-1)*m]
            mult_part = prev_row[:, (N-1)*m:] @ full_G[:m, :]
            full_G[(j-1)*m:j*m, :] = shifted + mult_part

        G = full_G

    return G


def _nsf_pi(
    A: np.ndarray,
    G: np.ndarray,
    N: int,
    max_num_comp: int = 1000,
    verbose: int = 0
) -> np.ndarray:
    """
    Compute stationary distribution of NSF Markov chain.

    Args:
        A: Block matrix [A0 A1 ... A_max]
        G: G matrix from NSF-GHT
        N: Skip factor
        max_num_comp: Maximum components
        verbose: Verbosity level

    Returns:
        Stationary distribution vector
    """
    m = A.shape[0]
    K = A.shape[1] // m - 1

    # Construct reblocked C matrices
    extra_blocks_C = (N - 1) - ((K - 1) % N) if K > 0 else 0
    C_width = m * (N - 1) + A.shape[1] + m * extra_blocks_C
    C = np.zeros((N * m, C_width))

    # Fill C with shifted copies of A
    for row_block in range(N):
        row_start = row_block * m
        col_start = row_block * m
        for col in range(A.shape[1]):
            if col_start + col < C_width:
                C[row_start:row_start+m, col_start+col] = A[:, col].reshape(-1, 1).flatten()

    # Homogeneous case
    # Compute FirstRowSumsG
    first_row_sums_G = np.zeros((m, N * m))
    G1 = G[:m, :] if G.shape[0] >= m else G

    for i in range(m):
        for block in range(N):
            if block == 0:
                first_row_sums_G[i, block*m:(block+1)*m] = G1[i, block*m:(block+1)*m]
            else:
                first_row_sums_G[i, block*m:(block+1)*m] = (
                    first_row_sums_G[i, (block-1)*m:block*m] + G1[i, block*m:(block+1)*m]
                )

    # ghat = stat(FirstRowSumsG[:, (N-1)*m:])
    last_block = first_row_sums_G[:, (N-1)*m:N*m]
    # Compute stationary distribution of last_block
    try:
        eigenvalues, eigenvectors = linalg.eig(last_block.T)
        idx = np.argmin(np.abs(eigenvalues))
        ghat = np.real(eigenvectors[:, idx])
        ghat = ghat / np.sum(ghat)
        ghat = ghat.reshape(1, -1)
    except:
        ghat = np.ones((1, m)) / m

    # pi0 = ghat @ G1
    pi0 = ghat @ G1

    # Compute normalization factor
    g = ghat @ first_row_sums_G

    # Simplified normalization using series summation
    # For M/G/1-type, pi = pi0 * (I-R)^(-1) where R is derived from G
    # Use iterative computation for stability

    # Compute R matrix: R = sum(A_i * G^i)
    R = A[:, :m].copy()
    G_pow = G1.copy()
    for i in range(1, K + 1):
        if (i+1)*m <= A.shape[1]:
            R = R + A[:, i*m:(i+1)*m] @ G_pow
        G_pow = G_pow @ G1

    # Compute stationary vector using R
    pi_list = [pi0.flatten()]
    total_mass = np.sum(pi0)

    R_pow = R.copy()
    for i in range(1, max_num_comp):
        pi_i = pi0 @ R_pow
        mass = np.sum(pi_i)
        if mass < 1e-15:
            break
        pi_list.append(pi_i.flatten())
        total_mass += mass
        R_pow = R_pow @ R

    # Normalize
    pi = np.concatenate(pi_list)
    pi = pi / np.sum(pi)

    return pi
