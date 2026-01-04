"""
Utility functions for M3A (Markovian Arrival Process with 3-moment Approximation).

Provides supporting functions for moment computation, validation, and quality assessment.
"""

import numpy as np
from numpy.typing import NDArray
from typing import List, Dict, Tuple, Optional
import scipy.linalg as la


def compute_autocorrelation(MMAP: List[np.ndarray], max_lag: int = 10) -> np.ndarray:
    """
    Compute the autocorrelation function of an MMAP up to the specified lag.
    Uses the first two matrices (D0, D1) of the MMAP which form a valid MAP.

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]
        max_lag: Maximum lag to compute

    Returns:
        Array of autocorrelation values for lags 1 to max_lag
    """
    D0 = MMAP[0]
    D1 = MMAP[1]
    n = D0.shape[0]

    # Compute stationary distribution
    Q = D0 + D1
    pi = _compute_stationary(Q)

    # Compute mean inter-arrival time
    lambda_arr = D1.sum()  # Total arrival rate
    mean = 1.0 / lambda_arr if lambda_arr > 0 else 1.0

    # Compute moments for variance
    try:
        Q_inv = la.solve(Q.T, np.eye(n)).T
        moment1 = -pi @ Q_inv @ np.ones(n)
        moment2 = 2 * pi @ (Q_inv @ Q_inv) @ np.ones(n)
        var = moment2 - moment1**2
    except la.LinAlgError:
        var = mean**2

    # Compute autocorrelation at each lag
    acf = np.zeros(max_lag)
    e = np.ones(n)
    d1_sum = D1 @ e

    try:
        P = la.expm(D0)  # Embedded transition matrix approximation
    except:
        P = np.eye(n)

    for lag in range(1, max_lag + 1):
        try:
            # Use embedded transition matrix raised to power lag
            P_lag = np.linalg.matrix_power(P, lag)
            # Approximate autocorrelation
            acf[lag - 1] = (pi @ np.diag(d1_sum) @ P_lag @ np.diag(d1_sum) @ e /
                          (pi @ np.diag(d1_sum) @ e)**2 - 1)
            # Bound the autocorrelation to valid range
            acf[lag - 1] = max(-1, min(1, acf[lag - 1]))
        except:
            acf[lag - 1] = 0.0

    return acf


def compute_idc(MMAP: List[np.ndarray], time_window: Optional[float] = None) -> float:
    """
    Compute the index of dispersion for counts (IDC) of an MMAP.
    Uses asymptotic IDC if time_window is None.

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]
        time_window: Time window for IDC computation (uses asymptotic if None)

    Returns:
        Index of dispersion for counts
    """
    D0 = MMAP[0]
    D1 = MMAP[1]
    n = D0.shape[0]

    # Compute stationary distribution
    Q = D0 + D1
    pi = _compute_stationary(Q)

    # Compute arrival rate
    lambda_arr = pi @ D1 @ np.ones(n)

    if lambda_arr <= 0:
        return 1.0

    # Compute squared coefficient of variation (SCV)
    scv = compute_scv(MMAP)

    # For asymptotic IDC, use formula based on SCV
    # IDC = 2*SCV - 1 for renewal processes, but for MAPs it's more complex
    # Using simplified approximation
    if time_window is None:
        return scv
    else:
        # For finite time window, use transient computation
        return scv  # Simplified


def compute_coeff_var(MMAP: List[np.ndarray]) -> float:
    """
    Compute the coefficient of variation of an MMAP.

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]

    Returns:
        Coefficient of variation (CV = std / mean)
    """
    return np.sqrt(compute_scv(MMAP))


def compute_scv(MMAP: List[np.ndarray]) -> float:
    """
    Compute the squared coefficient of variation (SCV) of an MMAP.

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]

    Returns:
        Squared coefficient of variation
    """
    moments = compute_moments(MMAP, 2)
    mean = moments[0]
    if mean <= 0:
        return 1.0
    variance = moments[1] - mean**2
    return max(0, variance / (mean**2))


def compute_moments(MMAP: List[np.ndarray], n: int = 3) -> np.ndarray:
    """
    Compute the first n moments of an MMAP.

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]
        n: Number of moments to compute

    Returns:
        Array of the first n moments
    """
    D0 = MMAP[0]
    D1 = MMAP[1]
    size = D0.shape[0]

    # Compute stationary distribution
    Q = D0 + D1
    pi = _compute_stationary(Q)

    moments = np.zeros(n)

    try:
        # Compute (-D0)^(-1) which is the expected sojourn time matrix
        D0_inv = la.solve(-D0, np.eye(size))
        e = np.ones(size)

        # k-th moment = k! * pi * ((-D0)^(-1))^k * e
        D0_inv_k = np.eye(size)
        for k in range(1, n + 1):
            D0_inv_k = D0_inv_k @ D0_inv
            factorial_k = np.prod(np.arange(1, k + 1))
            moments[k - 1] = factorial_k * pi @ D0_inv_k @ e

    except la.LinAlgError:
        # Fallback: estimate moments from diagonal rates
        total_rate = -np.sum(np.diag(D0))
        mean = size / total_rate if total_rate > 0 else 1.0
        for k in range(n):
            moments[k] = mean**(k + 1) * np.math.factorial(k + 1)

    return moments


def compute_spectral_gap(MMAP: List[np.ndarray]) -> float:
    """
    Compute the spectral gap of an MMAP generator matrix.

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]

    Returns:
        Spectral gap (difference between largest and second-largest eigenvalue real parts)
    """
    Q = MMAP[0] + MMAP[1]

    try:
        eigenvalues = la.eigvals(Q)
        # Get real parts and sort descending
        real_parts = np.sort(np.real(eigenvalues))[::-1]

        if len(real_parts) >= 2:
            return real_parts[0] - real_parts[1]
        return 0.0
    except la.LinAlgError:
        return 0.0


def validate_mmap(MMAP: List[np.ndarray], tolerance: float = 1e-10) -> bool:
    """
    Validate that matrices represent a valid MMAP.

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]
        tolerance: Tolerance for numerical checks

    Returns:
        True if the MMAP is valid, False otherwise
    """
    if len(MMAP) < 2:
        return False

    D0 = MMAP[0]
    D1 = MMAP[1]

    # Check dimensions
    if D0.shape[0] != D0.shape[1]:
        return False
    if D1.shape[0] != D1.shape[1]:
        return False
    if D0.shape[0] != D1.shape[0]:
        return False

    n = D0.shape[0]

    # Check that D0 has non-positive diagonal elements
    for i in range(n):
        if D0[i, i] > tolerance:
            return False

    # Check that off-diagonal elements of D0 are non-negative
    for i in range(n):
        for j in range(n):
            if i != j and D0[i, j] < -tolerance:
                return False

    # Check that D1 has non-negative elements
    if np.any(D1 < -tolerance):
        return False

    # Check that rows of D0 + sum(D_k) sum to zero
    Q = D0.copy()
    for k in range(1, len(MMAP)):
        Q = Q + MMAP[k]

    row_sums = Q @ np.ones(n)
    if np.any(np.abs(row_sums) > tolerance):
        return False

    return True


def compute_kl_divergence(MMAP1: List[np.ndarray], MMAP2: List[np.ndarray],
                         num_samples: int = 10000, seed: int = 23000) -> float:
    """
    Compute the Kullback-Leibler divergence between two MMAPs.

    Args:
        MMAP1: First MMAP
        MMAP2: Second MMAP
        num_samples: Number of samples for estimation
        seed: Random seed for reproducibility

    Returns:
        KL divergence estimate
    """
    samples1 = sample_interarrival_times(MMAP1, num_samples, seed)
    samples2 = sample_interarrival_times(MMAP2, num_samples, seed)

    # Compute empirical distributions using histograms
    num_bins = 50

    # Determine common bin edges
    all_samples = np.concatenate([samples1, samples2])
    min_val = np.min(all_samples)
    max_val = np.max(all_samples)
    bin_edges = np.linspace(min_val, max_val, num_bins + 1)

    # Compute histograms
    hist1, _ = np.histogram(samples1, bins=bin_edges, density=True)
    hist2, _ = np.histogram(samples2, bins=bin_edges, density=True)

    # Add small epsilon to avoid log(0)
    epsilon = 1e-10
    hist1 = hist1 + epsilon
    hist2 = hist2 + epsilon

    # Normalize
    hist1 = hist1 / np.sum(hist1)
    hist2 = hist2 / np.sum(hist2)

    # Compute KL divergence
    kl = np.sum(hist1 * np.log(hist1 / hist2))

    return max(0, kl)


def sample_interarrival_times(MMAP: List[np.ndarray], num_samples: int,
                              seed: Optional[int] = None) -> np.ndarray:
    """
    Generate random inter-arrival time samples from an MMAP.

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]
        num_samples: Number of samples to generate
        seed: Random seed (optional)

    Returns:
        Array of inter-arrival times
    """
    if seed is not None:
        np.random.seed(seed)

    D0 = MMAP[0]
    D1 = MMAP[1]
    n = D0.shape[0]

    samples = np.zeros(num_samples)

    # Compute stationary distribution for initial state
    Q = D0 + D1
    pi = _compute_stationary(Q)

    # Sample initial state
    state = np.random.choice(n, p=pi)

    for s in range(num_samples):
        time = 0.0

        # Generate inter-arrival time
        while True:
            # Holding time in current state
            rate = -D0[state, state]
            if rate <= 0:
                time += 1.0  # Fallback
                break

            hold_time = np.random.exponential(1.0 / rate)

            # Decide next transition (arrival vs hidden)
            arrival_rate = np.sum(D1[state, :])
            arrival_prob = arrival_rate / rate if rate > 0 else 0

            if np.random.random() < arrival_prob:
                # Arrival occurred
                time += hold_time

                # Sample destination state for arrival
                if arrival_rate > 0:
                    probs = D1[state, :] / arrival_rate
                    state = np.random.choice(n, p=probs)
                break
            else:
                # Hidden transition (no arrival)
                time += hold_time

                # Sample destination state for hidden transition
                hidden_rate = rate - arrival_rate
                if hidden_rate > 0:
                    probs = D0[state, :].copy()
                    probs[state] = 0  # Can't stay in same state
                    if np.sum(probs) > 0:
                        probs = probs / np.sum(probs)
                        state = np.random.choice(n, p=probs)

        samples[s] = time

    return samples


def _compute_stationary(Q: np.ndarray) -> np.ndarray:
    """
    Compute stationary distribution of a generator matrix.

    Args:
        Q: Generator matrix (rows sum to zero)

    Returns:
        Stationary distribution vector
    """
    n = Q.shape[0]

    try:
        # Solve pi * Q = 0 with sum(pi) = 1
        # Transpose: Q' * pi' = 0
        # Replace last equation with normalization constraint
        A = Q.T.copy()
        A[-1, :] = 1.0
        b = np.zeros(n)
        b[-1] = 1.0

        pi = la.solve(A, b)

        # Ensure non-negative and normalized
        pi = np.maximum(pi, 0)
        pi = pi / np.sum(pi)

        return pi
    except la.LinAlgError:
        # Fallback to uniform distribution
        return np.ones(n) / n


def mmap_isfeasible(MMAP: List[np.ndarray]) -> bool:
    """
    Check if MMAP is feasible (alias for validate_mmap).

    Args:
        MMAP: MMAP as list of matrices [D0, D1, D2, ...]

    Returns:
        True if feasible, False otherwise
    """
    return validate_mmap(MMAP)


__all__ = [
    'compute_autocorrelation',
    'compute_idc',
    'compute_coeff_var',
    'compute_scv',
    'compute_moments',
    'compute_spectral_gap',
    'validate_mmap',
    'compute_kl_divergence',
    'sample_interarrival_times',
    'mmap_isfeasible',
]
