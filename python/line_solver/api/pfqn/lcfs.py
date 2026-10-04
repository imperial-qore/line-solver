"""
LCFS Queueing Network algorithms for Product-Form Queueing Networks.

Native Python implementations for computing normalizing constants and
performance metrics in LCFS queueing networks.

Key functions:
    pfqn_lcfsqn_nc: Normalizing constant for LCFS queueing networks
    pfqn_lcfsqn_mva: Mean Value Analysis for LCFS queueing networks

References:
    G. Casale, "A family of multiclass LCFS queueing networks with
    order-dependent product-form solutions", QUESTA 2026.
"""

import numpy as np
from math import log, exp, comb, factorial
from typing import Tuple, List, Optional
from dataclasses import dataclass

from ..perm import compute_permanent


@dataclass
class LcfsqnNcResult:
    """Result of LCFS NC algorithm."""
    G: float
    Ax: List[np.ndarray]


@dataclass
class LcfsqnMvaResult:
    """Result of LCFS MVA algorithm."""
    T: np.ndarray  # Throughput vector (1 x R)
    Q: np.ndarray  # Queue length matrix (2 x R)
    U: np.ndarray  # Utilization matrix (2 x R)
    B: np.ndarray  # Back probability matrix (2 x R)


def _ryser_permanent(matrix: np.ndarray) -> float:
    """
    Compute the permanent of a matrix using Ryser's algorithm with Gray code.

    Uses the inclusion-exclusion principle with Gray code ordering for
    efficient computation.

    Args:
        matrix: Square matrix (n x n)

    Returns:
        Permanent value
    """
    n = matrix.shape[0]
    if n == 0:
        return 1.0
    if n == 1:
        return matrix[0, 0]

    # Generate Gray code bit modification order
    def bit_to_modify(m: int) -> List[int]:
        if m == 1:
            return [0]
        else:
            sub = bit_to_modify(m - 1)
            return sub + [m - 1] + sub

    permanent = 0.0
    row_sum = np.zeros(n)
    current_bit = np.zeros(n, dtype=bool)
    bit_order = bit_to_modify(n)

    for bit_index in bit_order:
        # Toggle bit
        current_bit[bit_index] = not current_bit[bit_index]

        # Update row sums
        if current_bit[bit_index]:
            row_sum += matrix[:, bit_index]
        else:
            row_sum -= matrix[:, bit_index]

        # Number of columns selected
        nb_col = np.sum(current_bit)

        # Calculate product of row sums
        product = np.prod(row_sum)

        # Add to permanent with appropriate sign
        permanent += ((-1) ** nb_col) * product

    permanent *= ((-1) ** n)
    return permanent


def _permanent_with_multiplicities(matrix: np.ndarray, multiplicities: np.ndarray) -> float:
    """
    Compute permanent of a matrix with repeated rows based on multiplicities.

    Uses the inclusion-exclusion principle with multinomial coefficients
    to efficiently handle matrices with duplicate rows.

    Args:
        matrix: Matrix (R x K) with unique rows
        multiplicities: (R,) array of row repetition counts

    Returns:
        Permanent value of the expanded matrix
    """
    R = len(multiplicities)
    n = int(np.sum(multiplicities))

    if n == 0:
        return 1.0

    # pprod iterator - generates all states from 0 to multiplicities
    def pprod_next(current: np.ndarray, bounds: np.ndarray) -> Optional[np.ndarray]:
        """Generate next state in product space."""
        result = current.copy()
        R_local = len(bounds)

        # Check if at maximum
        if np.all(result == bounds):
            return None

        s = R_local - 1
        while s >= 0 and result[s] == bounds[s]:
            result[s] = 0
            s -= 1

        if s < 0:
            return None

        result[s] += 1
        return result

    value = 0.0
    f = np.zeros(R, dtype=int)

    while True:
        term = ((-1) ** np.sum(f))

        # Multinomial coefficients
        for j in range(R):
            term *= comb(int(multiplicities[j]), int(f[j]))

        # Product term
        for i in range(n):
            sum_term = 0.0
            for k in range(R):
                if k < matrix.shape[0] and i < matrix.shape[1]:
                    sum_term += f[k] * matrix[k, i]
            term *= sum_term

        value += term

        f = pprod_next(f, multiplicities.astype(int))
        if f is None:
            break

    return ((-1) ** n) * value


def _perm_with_population(A: np.ndarray, N: np.ndarray) -> float:
    """
    Compute permanent of matrix A with rows repeated according to population N.

    Args:
        A: Matrix (R x K)
        N: Population vector (R,) specifying row repetitions

    Returns:
        Permanent of the expanded matrix
    """
    K = int(np.sum(N))
    R = len(N)

    if K == 0:
        return 1.0

    # Build expanded matrix with rows repeated according to N
    expanded_matrix = np.zeros((K, K))
    row_idx = 0
    for r in range(R):
        count = int(N[r])
        for _ in range(count):
            if row_idx < K and A.shape[1] >= K:
                expanded_matrix[row_idx, :] = A[r, :K]
            elif row_idx < K:
                expanded_matrix[row_idx, :A.shape[1]] = A[r, :]
            row_idx += 1

    # Routed through the shared permanent so the expansion picks the
    # better-conditioned orientation; see _kb/03-api-layer.md.
    return compute_permanent(expanded_matrix)


def pfqn_lcfsqn_nc(alpha: np.ndarray, beta: np.ndarray,
                   N: np.ndarray) -> LcfsqnNcResult:
    """
    Normalizing constant for multiclass LCFS queueing networks.

    Computes the normalizing constant for a 2-station closed queueing network with:
    - Station 1: LCFS (Last-Come-First-Served, non-preemptive)
    - Station 2: LCFS-PR (LCFS with Preemption-Resume)

    Args:
        alpha: Vector of inverse service rates at station 1 (LCFS).
               alpha(r) = 1/mu(1,r) for class r
        beta: Vector of inverse service rates at station 2 (LCFS-PR).
              beta(r) = 1/mu(2,r) for class r
        N: Population vector, N(r) = number of jobs of class r.

    Returns:
        LcfsqnNcResult containing normalizing constant G and A matrices
    """
    alpha = np.asarray(alpha, dtype=float).flatten()
    beta = np.asarray(beta, dtype=float).flatten()
    N = np.asarray(N, dtype=float).flatten()

    K = int(np.sum(N))
    R = len(N)

    if K == 0:
        return LcfsqnNcResult(G=1.0, Ax=[])

    def make_A(x: int) -> np.ndarray:
        """Construct A matrix for state x."""
        A = np.zeros((R, K))

        for i in range(R):
            # Columns 0 to x-1: alpha(i)^(j+1) for j=0..x-1
            for j in range(x):
                A[i, j] = alpha[i] ** (j + 1)

            # Columns x to K-1: alpha(i)^(x+j) * beta(i) for j=0..K-x-1
            for j in range(K - x):
                A[i, x + j] = (alpha[i] ** (x + j)) * beta[i]

        return A

    # Create array of A matrices for each state x=0:K
    Ax = [make_A(x) for x in range(K + 1)]

    # Sum over all states
    G = 0.0
    for x in range(K + 1):
        G += _perm_with_population(Ax[x], N)

    # see _kb/03-api-layer.md for rationale
    G = G / float(np.prod([factorial(int(n)) for n in N]))

    return LcfsqnNcResult(G=G, Ax=Ax)


def _hashpop_lcfs(n: np.ndarray, N: np.ndarray, prods: np.ndarray) -> int:
    """Compute hash index for LCFS population lattice."""
    idx = 0
    for r in range(len(N)):
        idx += int(prods[r] * n[r])
    return idx


def _oner(n: np.ndarray, r: int) -> np.ndarray:
    """Decrement the r-th component of vector n by 1."""
    result = n.copy()
    result[r] -= 1
    return result


def _logsumexp(terms: List[float], log_scales: List[float]) -> float:
    """
    Compute sum of scaled terms using log-sum-exp.

    Computes: result = sum_i (terms(i) * exp(log_scales(i)))
    """
    if not terms:
        return 0.0

    # Filter out invalid terms
    valid_pairs = [(t, s) for t, s in zip(terms, log_scales)
                   if t > 0 and np.isfinite(s)]

    if not valid_pairs:
        return 0.0

    # Compute log of each term
    log_values = [log(t) + s for t, s in valid_pairs]

    # Find maximum for numerical stability
    max_log = max(log_values)

    if not np.isfinite(max_log):
        return 0.0

    # log-sum-exp
    return exp(max_log) * sum(exp(lv - max_log) for lv in log_values)


def _pprod_init(N: np.ndarray) -> np.ndarray:
    """Initialize population product iterator."""
    return np.zeros_like(N, dtype=float)


def _pprod_next(n: np.ndarray, N: np.ndarray) -> Optional[np.ndarray]:
    """Get next population vector in lexicographic order."""
    R = len(N)
    result = n.copy()

    # Check if at maximum
    if np.all(result == N):
        return None

    s = R - 1
    while s >= 0 and result[s] == N[s]:
        result[s] = 0
        s -= 1

    if s < 0:
        return None

    result[s] += 1
    return result


def pfqn_lcfsqn_mva(alpha: np.ndarray, beta: np.ndarray,
                    N: Optional[np.ndarray] = None) -> LcfsqnMvaResult:
    """
    Mean Value Analysis for multiclass LCFS queueing networks.

    Computes performance metrics for a 2-station closed queueing network with:
    - Station 1: LCFS (Last-Come-First-Served, non-preemptive)
    - Station 2: LCFS-PR (LCFS with Preemption-Resume)

    This implementation uses log-space arithmetic to prevent numerical underflow.

    Args:
        alpha: Vector of inverse service rates at station 1 (LCFS).
               alpha(r) = 1/mu(1,r) for class r
        beta: Vector of inverse service rates at station 2 (LCFS-PR).
              beta(r) = 1/mu(2,r) for class r
        N: Population vector, N(r) = number of jobs of class r.
           Default: ones(R) - one job per class

    Returns:
        LcfsqnMvaResult containing throughputs, queue lengths, utilizations,
        and back probabilities
    """
    alpha = np.asarray(alpha, dtype=float).flatten()
    beta = np.asarray(beta, dtype=float).flatten()
    R = len(alpha)

    if N is None:
        N = np.ones(R)
    else:
        N = np.asarray(N, dtype=float).flatten()

    K = int(np.sum(N))

    if K == 0:
        return LcfsqnMvaResult(
            T=np.zeros((1, R)),
            Q=np.zeros((2, R)),
            U=np.zeros((2, R)),
            B=np.zeros((2, R))
        )

    # Precompute log values for numerical stability
    log_alpha = np.log(alpha)
    log_beta = np.log(beta)

    # Precompute products for hash computation
    prods = np.ones(R)
    for r in range(R):
        prod = 1.0
        for j in range(r):
            prod *= (N[j] + 1)
        prods[r] = prod

    # Calculate total size for arrays
    total_size = 1
    for r in range(R):
        total_size *= int(N[r] + 1)

    # Storage for intermediate results
    QN = [np.zeros((2, R)) for _ in range(total_size)]
    TN = [np.zeros((1, R)) for _ in range(total_size)]
    BN_scaled = [np.zeros((2, R)) for _ in range(total_size)]
    log_scale_BN = [[np.full(R, float('-inf')) for _ in range(2)]
                    for _ in range(total_size)]

    # Initialize population iterator
    n = _pprod_init(N)

    while n is not None:
        idx = _hashpop_lcfs(n, N, prods)
        sum_N = int(np.sum(n))

        if sum_N == 1:
            # Base case: single job
            for r in range(R):
                if n[r] > 0:
                    denom = alpha[r] + beta[r]
                    QN[idx][0, r] = alpha[r] / denom
                    QN[idx][1, r] = beta[r] / denom
                    BN_scaled[idx][0, r] = alpha[r] / denom
                    BN_scaled[idx][1, r] = beta[r] / denom
                    log_scale_BN[idx][0][r] = 0.0
                    log_scale_BN[idx][1][r] = 0.0
                    TN[idx][0, r] = 1.0 / denom

        elif sum_N > 1:
            # Recursive case: multiple jobs
            # log_scale_np = log(prod(alpha.^n))
            log_scale_np = np.dot(n, log_alpha)

            for k in range(R):
                if n[k] > 0:
                    n_minus_k = _oner(n, k)
                    idx_k = _hashpop_lcfs(n_minus_k, N, prods)

                    # Compute unscaled waiting time contributions
                    Wnp_unscaled = 1.0 + QN[idx_k][0, k]
                    Wpr_unscaled = 1.0 + QN[idx_k][1, k]

                    for r in range(R):
                        if r != k and n[r] > 0:
                            n_minus_r = _oner(n, r)
                            idx_r = _hashpop_lcfs(n_minus_r, N, prods)

                            # Compute B ratio in log-space
                            if (BN_scaled[idx_k][0, r] > 0 and
                                BN_scaled[idx_r][0, k] > 0):
                                log_ratio_np = (log(BN_scaled[idx_k][0, r]) +
                                               log_scale_BN[idx_k][0][r] -
                                               log(BN_scaled[idx_r][0, k]) -
                                               log_scale_BN[idx_r][0][k])
                                ratio_np = exp(log_ratio_np)
                                Wnp_unscaled += ((alpha[k] / alpha[r]) *
                                                 ratio_np * QN[idx_r][0, k])

                            if (BN_scaled[idx_k][1, r] > 0 and
                                BN_scaled[idx_r][1, k] > 0):
                                log_ratio_pr = (log(BN_scaled[idx_k][1, r]) +
                                               log_scale_BN[idx_k][1][r] -
                                               log(BN_scaled[idx_r][1, k]) -
                                               log_scale_BN[idx_r][1][k])
                                ratio_pr = exp(log_ratio_pr)
                                Wpr_unscaled += ((alpha[r] / alpha[k]) *
                                                 ratio_pr * QN[idx_r][1, k])

                    # log_scale_pr_k = log(alpha(k)^(sum(n)-1) * beta(k))
                    log_scale_pr_k = (sum_N - 1) * log_alpha[k] + log_beta[k]

                    # Compute Y(k) = n(k) / (Wnp + Wpr) using log-sum-exp
                    log_Wnp = log_scale_np + log(max(Wnp_unscaled, 1e-300))
                    log_Wpr = log_scale_pr_k + log(max(Wpr_unscaled, 1e-300))

                    # log-sum-exp for numerical stability
                    max_log = max(log_Wnp, log_Wpr)
                    log_sum_W = max_log + log(exp(log_Wnp - max_log) +
                                               exp(log_Wpr - max_log))

                    # B(1,k) = prod(alpha.^n) * n(k) / (Wnp + Wpr)
                    BN_scaled[idx][0, k] = n[k]
                    log_scale_BN[idx][0][k] = log_scale_np - log_sum_W

                    # B(2,k) = alpha(k)^(sum(n)-1) * beta(k) * n(k) / (Wnp + Wpr)
                    BN_scaled[idx][1, k] = n[k]
                    log_scale_BN[idx][1][k] = log_scale_pr_k - log_sum_W

            # Obtain queue lengths: Q(i,k) = B(i,k) + sum_r B(i,r) * Q_{n-e_r}(i,k)
            for k in range(R):
                if n[k] > 0:
                    for station in range(2):
                        terms = []
                        log_scales = []

                        # First term: B(station,k)
                        if BN_scaled[idx][station, k] > 0:
                            terms.append(BN_scaled[idx][station, k])
                            log_scales.append(log_scale_BN[idx][station][k])

                        # Sum over r: B(station,r) * Q_{n-e_r}(station,k)
                        for r in range(R):
                            if n[r] > 0:
                                n_minus_r = _oner(n, r)
                                idx_r = _hashpop_lcfs(n_minus_r, N, prods)
                                if (BN_scaled[idx][station, r] > 0 and
                                    QN[idx_r][station, k] > 0):
                                    terms.append(BN_scaled[idx][station, r] *
                                                QN[idx_r][station, k])
                                    log_scales.append(log_scale_BN[idx][station][r])

                        QN[idx][station, k] = _logsumexp(terms, log_scales)

            # Obtain throughput
            for k in range(R):
                if n[k] > 0:
                    n_minus_k = _oner(n, k)
                    idx_k = _hashpop_lcfs(n_minus_k, N, prods)

                    # Compute Unp = sum_r alpha(r) * T_{n-e_k}(r)
                    Unp = np.dot(alpha, TN[idx_k][0, :])

                    terms = []
                    log_scales = []

                    for r in range(R):
                        if n[r] > 0:
                            n_minus_r = _oner(n, r)
                            idx_r = _hashpop_lcfs(n_minus_r, N, prods)
                            if (BN_scaled[idx][0, r] > 0 and
                                TN[idx_r][0, k] > 0):
                                terms.append(BN_scaled[idx][0, r] *
                                            TN[idx_r][0, k])
                                log_scales.append(log_scale_BN[idx][0][r])

                    # Add (1/alpha(k)) * B(1,k) * (1-Unp) term
                    if BN_scaled[idx][0, k] > 0 and (1 - Unp) > 0:
                        terms.append((1.0 / alpha[k]) * BN_scaled[idx][0, k] *
                                    (1 - Unp))
                        log_scales.append(log_scale_BN[idx][0][k])

                    TN[idx][0, k] = _logsumexp(terms, log_scales)

        n = _pprod_next(n, N)

    # Get final results
    final_idx = total_size - 1
    Q = QN[final_idx]
    T = TN[final_idx]

    # Recover actual B values from scaled representation
    B = np.zeros((2, R))
    for r in range(R):
        if BN_scaled[final_idx][0, r] > 0:
            B[0, r] = BN_scaled[final_idx][0, r] * exp(log_scale_BN[final_idx][0][r])
        if BN_scaled[final_idx][1, r] > 0:
            B[1, r] = BN_scaled[final_idx][1, r] * exp(log_scale_BN[final_idx][1][r])

    # Compute utilizations
    U = np.zeros((2, R))
    for r in range(R):
        U[0, r] = T[0, r] * alpha[r]
        U[1, r] = T[0, r] * beta[r]

    return LcfsqnMvaResult(T=T, Q=Q, U=U, B=B)


def _pprod(n: np.ndarray, N: np.ndarray = None) -> np.ndarray:
    """
    Population product iterator - get next vector in population lattice.

    Args:
        n: Current population vector
        N: Maximum population (bounds)

    Returns:
        Next population vector, or vector with negative first element if exhausted
    """
    R = len(n)
    result = n.copy()

    if N is None:
        # Just increment last position
        result[-1] += 1
        return result

    # Check if at maximum
    if np.all(result >= N):
        return -np.ones(R)

    # Find rightmost position that can be incremented
    s = R - 1
    while s >= 0 and result[s] >= N[s]:
        result[s] = 0
        s -= 1

    if s < 0:
        return -np.ones(R)

    result[s] += 1
    return result


def _hashpop(n: np.ndarray, N: np.ndarray, R: int, prods: np.ndarray) -> int:
    """
    Compute hash index for population vector.

    Args:
        n: Current population vector
        N: Maximum population
        R: Number of classes
        prods: Precomputed products for indexing

    Returns:
        Linear index in flattened population lattice
    """
    idx = 0
    for r in range(R):
        idx += int(prods[r] * n[r])
    return int(idx)


def pfqn_lcfsqn_ca(
    alpha: np.ndarray,
    beta: np.ndarray,
    N: np.ndarray = None
) -> Tuple[float, float]:
    """
    Convolution algorithm for multiclass LCFS queueing networks.

    Computes the normalizing constant for a 2-station closed queueing network with:
    - Station 1: LCFS (Last-Come-First-Served, non-preemptive)
    - Station 2: LCFS-PR (LCFS with Preemption-Resume)

    This algorithm uses a recursive convolution approach that is more
    efficient than direct computation for networks with multiple job classes.

    Args:
        alpha: Vector of inverse service rates at station 1 (LCFS).
               alpha(r) = 1/mu(1,r) for class r
        beta: Vector of inverse service rates at station 2 (LCFS-PR).
              beta(r) = 1/mu(2,r) for class r
        N: Population vector, N(r) = number of jobs of class r.
           Default: ones(R) - one job per class

    Returns:
        Tuple of (G, V) where:
            G: Normalizing constant
            V: Auxiliary normalization term

    Reference:
        G. Casale, "A family of multiclass LCFS queueing networks with
        order-dependent product-form solutions", QUESTA 2026.

    Example:
        >>> G, V = pfqn_lcfsqn_ca([1.0, 2.0], [0.5, 0.5], [2, 3])
    """
    alpha = np.asarray(alpha, dtype=float).flatten()
    beta = np.asarray(beta, dtype=float).flatten()
    R = len(alpha)

    if N is None:
        N = np.ones(R)
    else:
        N = np.asarray(N, dtype=float).flatten()

    K = int(np.sum(N))

    # Empty system case
    if K == 0:
        return 1.0, 1.0

    # Precompute products for hashing
    prods = np.zeros(R)
    for r in range(R):
        prods[r] = np.prod(N[:r] + 1)

    # Allocate arrays for G and V indexed by population state
    total_size = int(np.prod(N + 1))
    G = np.zeros(total_size)
    V = np.zeros(total_size)

    # Initialize population iterator
    n = np.zeros(R)

    while n[0] >= 0:  # While not exhausted
        idx = _hashpop(n, N, R, prods)
        sum_n = int(np.sum(n))

        if sum_n == 0:
            # Base case: empty population
            G[idx] = 1.0
            V[idx] = 1.0
        else:
            # Recursive case
            V[idx] = 0.0
            G[idx] = 0.0

            for r in range(R):
                if n[r] > 0:
                    # Get index for n - e_r (one less job of class r)
                    n_minus_r = _oner(n, r)
                    idx_r = _hashpop(n_minus_r, N, R, prods)

                    # Recursive definition of V
                    V[idx] += V[idx_r]

                    # Recursive definition of G uses V and beta
                    G[idx] += (alpha[r] ** (sum_n - 1)) * beta[r] * G[idx_r]

            # Complete V computation
            V[idx] *= np.prod(alpha ** n)

            # Complete G computation
            G[idx] += V[idx]

        # Advance to next population
        n = _pprod(n, N)

    # Return final values (at population N)
    final_idx = total_size - 1
    return G[final_idx], V[final_idx]


@dataclass
class LcfsqnCaResult:
    """Result of LCFS CA algorithm."""
    G: float
    V: float


def pfqn_lcfsqn_mva_logsumexp(terms: np.ndarray, log_scales: np.ndarray) -> float:
    """
    Compute sum of scaled terms using log-sum-exp for numerical stability.

    Computes: result = sum_i (terms[i] * exp(log_scales[i]))

    This is mathematically exact - it computes the sum in a numerically
    stable way by factoring out the maximum exponent.

    Args:
        terms: Array of term values (coefficients)
        log_scales: Array of log scale factors

    Returns:
        Sum of scaled terms

    References:
        G. Casale, "A family of multiclass LCFS queueing networks with
        order-dependent product-form solutions", QUESTA 2026.
    """
    terms = np.atleast_1d(np.asarray(terms, dtype=float))
    log_scales = np.atleast_1d(np.asarray(log_scales, dtype=float))

    if len(terms) == 0:
        return 0.0

    # Filter out zero or negative terms
    valid = (terms > 0) & np.isfinite(log_scales)
    if not np.any(valid):
        return 0.0

    terms = terms[valid]
    log_scales = log_scales[valid]

    # Compute log of each term: log(term * exp(log_scale)) = log(term) + log_scale
    log_values = np.log(terms) + log_scales

    # Find maximum for numerical stability
    max_log = np.max(log_values)

    if not np.isfinite(max_log):
        return 0.0

    # log-sum-exp: log(sum(exp(x))) = max(x) + log(sum(exp(x - max(x))))
    # So: sum(exp(x)) = exp(max(x)) * sum(exp(x - max(x)))
    result = np.exp(max_log) * np.sum(np.exp(log_values - max_log))

    return float(result)


def pfqn_lcfsqn(
    L: np.ndarray,
    N: np.ndarray,
    lcfs_stat: int,
    lcfspr_stat: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    MVA for LCFS + LCFS-PR two-station closed queueing networks.

    Wrapper function that extracts service parameters from the demand matrix
    and calls the specialized LCFS MVA algorithm.

    Args:
        L: Demand matrix (M x K) where L[i, k] = service time * visits at station i for class k
        N: Population vector (K,) with number of jobs per class
        lcfs_stat: Index of LCFS station in L
        lcfspr_stat: Index of LCFS-PR station in L

    Returns:
        Tuple of (Q, U, R, T, C, X, lG) where:
            Q: Queue lengths (M x K)
            U: Utilizations (M x K)
            R: Response times (M x K)
            T: Throughputs (M x K)
            C: Cycle times (1 x K)
            X: System throughputs (1 x K)
            lG: Log normalizing constant

    References:
        G. Casale, "A family of multiclass LCFS queueing networks with
        order-dependent product-form solutions", QUESTA 2026.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()

    M = L.shape[0]
    K = L.shape[1]

    # see _kb/03-api-layer.md for rationale
    alpha = L[lcfs_stat, :].flatten()
    beta = L[lcfspr_stat, :].flatten()

    # Call the specialized MVA algorithm
    result = pfqn_lcfsqn_mva(alpha, beta, N)

    # Build full result matrices (M stations x K classes)
    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T_mat = np.zeros((M, K))

    # see _kb/03-api-layer.md for rationale
    Q[lcfs_stat, :] = result.Q[0, :]
    Q[lcfspr_stat, :] = result.Q[1, :]

    U[lcfs_stat, :] = result.U[0, :]
    U[lcfspr_stat, :] = result.U[1, :]

    # Throughputs at stations equal system throughput (closed network)
    T_mat[lcfs_stat, :] = result.T[0, :]
    T_mat[lcfspr_stat, :] = result.T[0, :]

    # Response times: R = Q / T (Little's law per station)
    for k in range(K):
        if result.T[0, k] > 0:
            R[lcfs_stat, k] = Q[lcfs_stat, k] / result.T[0, k]
            R[lcfspr_stat, k] = Q[lcfspr_stat, k] / result.T[0, k]

    # Cycle time is sum of response times
    C = np.zeros((1, K))
    for k in range(K):
        C[0, k] = R[lcfs_stat, k] + R[lcfspr_stat, k]

    # System throughput
    X = result.T.reshape(1, K)

    # Log normalizing constant (compute from NC result)
    nc_result = pfqn_lcfsqn_nc(alpha, beta, N)
    lG = np.log(nc_result.G) if nc_result.G > 0 else 0.0

    return Q, U, R, T_mat, C, X, lG


__all__ = [
    'pfqn_lcfsqn',
    'pfqn_lcfsqn_nc',
    'pfqn_lcfsqn_mva',
    'pfqn_lcfsqn_ca',
    'pfqn_lcfsqn_mva_logsumexp',
    'LcfsqnNcResult',
    'LcfsqnMvaResult',
    'LcfsqnCaResult',
]
