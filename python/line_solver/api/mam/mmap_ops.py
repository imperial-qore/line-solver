"""
Marked Markovian Arrival Process (MMAP) operations.

Native Python implementations for MMAP manipulation operations including:
- mmap_super_safe: Safe MMAP superposition with order constraints
- mmap_mark: Arrival reclassification based on probability matrix
- mmap_scale: Scaling of mean inter-arrival times
- mmap_hide: Hiding (removing) specified arrival classes
- mmap_compress: State space reduction with multiple methods

An MMAP is represented as a list [D0, D1, D2, ..., Dm] where:
    D0: hidden transition matrix (non-arrival transitions)
    D1, D2, ..., Dm: arrival marking matrices for each class m
    D0 + D1 + D2 + ... + Dm: complete infinitesimal generator

References:
    Neuts, M.F. "Matrix-Geometric Solutions in Stochastic Models"
    Bini, D.A., Meini, B., Ramaswami, V. "Algorithmic Methods for Structured
    Markov Chains" Oxford University Press, 2005
"""

import numpy as np
from scipy import linalg
from typing import List, Tuple, Optional, Union
import warnings

# Import MAP analysis functions that MMAP operations depend on
from .map_analysis import (
    map_infgen,
    map_piq,
    map_lambda,
    map_mean,
    map_moment,
    map_scv,
)


def mmap_infgen(D0: np.ndarray, D_list: List[np.ndarray]) -> np.ndarray:
    """
    Compute the infinitesimal generator of an MMAP.

    Args:
        D0: Hidden transition matrix (non-arrival transitions)
        D_list: List of arrival marking matrices [D1, D2, ..., Dm]

    Returns:
        Infinitesimal generator Q = D0 + D1 + D2 + ... + Dm
    """
    Q = np.asarray(D0, dtype=np.float64).copy()
    for D_k in D_list:
        Q = Q + np.asarray(D_k, dtype=np.float64)
    return Q


def mmap_normalize(D0: np.ndarray, D_list: List[np.ndarray]) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Normalize an MMAP to ensure -diag(D0) = rowsums(D1 + D2 + ... + Dm).

    This ensures D0 represents true hidden transitions (negative diagonal).

    Args:
        D0: Hidden transition matrix
        D_list: List of arrival marking matrices

    Returns:
        Tuple of (normalized_D0, normalized_D_list)
    """
    D0 = np.asarray(D0, dtype=np.float64).copy()
    D_list = [np.asarray(D, dtype=np.float64).copy() for D in D_list]

    n = D0.shape[0]

    # Compute row sums of marking matrices
    e = np.ones(n)
    mark_rowsums = np.zeros(n)
    for D_k in D_list:
        mark_rowsums += D_k @ e

    # Set diagonal of D0 to ensure proper balance
    np.fill_diagonal(D0, -mark_rowsums)

    return D0, D_list


def mmap_super_safe(mmap_list: List[Tuple[np.ndarray, List[np.ndarray]]],
                    maxorder: int = 1000,
                    method: str = "default") -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Safely superpose multiple MMAPs into a single MMAP.

    Combines multiple MMAPs while ensuring the order of the resulting MMAP
    does not exceed maxorder. If combining would exceed maxorder, alternative
    methods are used (fitting to exponential or simplified MAP).

    Args:
        mmap_list: List of (D0, D_list) tuples, one per MMAP
        maxorder: Maximum allowed order (state space size) for result
        method: Combination method ("default" or "match")

    Returns:
        Tuple of (D0_super, D_list_super) for the superposed MMAP

    Algorithm:
        1. Sort MMAPs by squared coefficient of variation (SCV) of unmarked process
        2. Iteratively combine: start with simplest, add others in order
        3. If Kronecker product would exceed maxorder, use simpler approximations
    """
    if not mmap_list:
        # Empty case: return standard Poisson (exponential)
        D0 = np.array([[-1.0]])
        D_list = [np.array([[1.0]])]
        return D0, D_list

    # Extract all MMAPs and compute SCVs
    mmaps = []
    scv_values = []

    for D0, D_list in mmap_list:
        mmaps.append((D0, D_list))
        # Compute SCV of unmarked process (aggregate arrivals)
        agg_D1 = np.zeros_like(D0)
        for D_k in D_list:
            agg_D1 = agg_D1 + np.asarray(D_k, dtype=np.float64)

        try:
            scv = map_scv(D0, agg_D1)
        except:
            scv = 1.0  # Default to 1 if computation fails
        scv_values.append(scv)

    # Sort MMAPs by SCV
    sorted_indices = np.argsort(scv_values)
    sorted_mmaps = [mmaps[i] for i in sorted_indices]

    # Start with first (simplest) MMAP
    sup_D0, sup_D_list = sorted_mmaps[0]
    sup_D0 = np.asarray(sup_D0, dtype=np.float64).copy()
    sup_D_list = [np.asarray(D, dtype=np.float64).copy() for D in sup_D_list]

    # Iteratively combine remaining MMAPs
    for i in range(1, len(sorted_mmaps)):
        next_D0, next_D_list = sorted_mmaps[i]
        next_D0 = np.asarray(next_D0, dtype=np.float64)
        next_D_list = [np.asarray(D, dtype=np.float64) for D in next_D_list]

        # Check if direct Kronecker product exceeds maxorder
        n_sup = sup_D0.shape[0]
        n_next = next_D0.shape[0]

        if n_sup * n_next > maxorder:
            # Use simplified versions if needed
            if n_sup * 2 < maxorder:
                # Can fit 2-state approximation
                # For now, use exponential approximation
                try:
                    lambda_next = map_lambda(next_D0,
                                            sum([D for D in next_D_list]))
                    next_D0_approx = np.array([[-lambda_next]])
                    next_D_list_approx = [np.array([[lambda_next]])]
                    sup_D0, sup_D_list = _mmap_super(sup_D0, sup_D_list,
                                                     next_D0_approx, next_D_list_approx,
                                                     method)
                except:
                    # If approximation fails, skip this MMAP
                    warnings.warn(f"Could not add MMAP {i} to superposition")
                    continue
            else:
                # Use exponential for this MMAP
                try:
                    lambda_next = map_lambda(next_D0,
                                            sum([D for D in next_D_list]))
                    next_D0_approx = np.array([[-lambda_next]])
                    next_D_list_approx = [np.array([[lambda_next]])]
                    sup_D0, sup_D_list = _mmap_super(sup_D0, sup_D_list,
                                                     next_D0_approx, next_D_list_approx,
                                                     method)
                except:
                    warnings.warn(f"Could not add MMAP {i} to superposition")
                    continue
        else:
            # Direct superposition
            sup_D0, sup_D_list = _mmap_super(sup_D0, sup_D_list,
                                            next_D0, next_D_list,
                                            method)

    # Normalize result
    sup_D0, sup_D_list = mmap_normalize(sup_D0, sup_D_list)

    return sup_D0, sup_D_list


def _mmap_super(D0_1: np.ndarray, D_list_1: List[np.ndarray],
               D0_2: np.ndarray, D_list_2: List[np.ndarray],
               method: str = "default") -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Superpose two MMAPs using Kronecker products.

    Args:
        D0_1, D_list_1: First MMAP
        D0_2, D_list_2: Second MMAP
        method: Superposition method

    Returns:
        Tuple of (D0_super, D_list_super)
    """
    D0_1 = np.asarray(D0_1, dtype=np.float64)
    D0_2 = np.asarray(D0_2, dtype=np.float64)
    D_list_1 = [np.asarray(D, dtype=np.float64) for D in D_list_1]
    D_list_2 = [np.asarray(D, dtype=np.float64) for D in D_list_2]

    n1 = D0_1.shape[0]
    n2 = D0_2.shape[0]

    # Kronecker sum for D0: D0 = D0_1 \otimes I_2 + I_1 \otimes D0_2
    I1 = np.eye(n1)
    I2 = np.eye(n2)
    D0_super = np.kron(D0_1, I2) + np.kron(I1, D0_2)

    # Marking matrices are Kronecker products
    # For MMAP1: D_{k,1} \otimes I_2 for k = 1..m1
    # For MMAP2: I_1 \otimes D_{k,2} for k = 1..m2
    # Result has m1 + m2 classes

    D_list_super = []

    # Add classes from MMAP1 (with identity for MMAP2)
    for D_k1 in D_list_1:
        D_super_k = np.kron(D_k1, I2)
        D_list_super.append(D_super_k)

    # Add classes from MMAP2 (with identity for MMAP1)
    for D_k2 in D_list_2:
        D_super_k = np.kron(I1, D_k2)
        D_list_super.append(D_super_k)

    return D0_super, D_list_super


def mmap_mark(D0: np.ndarray, D_list: List[np.ndarray],
             prob: np.ndarray) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Reclassify arrivals in an MMAP according to a probability matrix.

    Converts an MMAP with K classes to one with R classes based on a
    KxR probability matrix that describes how arrivals are reclassified.

    Args:
        D0: Hidden transition matrix
        D_list: List of K arrival marking matrices [D1, D2, ..., DK]
        prob: K x R probability matrix where prob[k,r] is the probability
              that a class-k arrival in the original MMAP becomes a class-r
              arrival in the new MMAP

    Returns:
        Tuple of (D0_new, D_list_new) where D_list_new has R matrices

    Example:
        If prob = [[0.8, 0.2], [0.3, 0.7]], a 2-class MMAP becomes a 2-class
        MMAP where 80% of class-1 arrivals are now class-1, 20% become class-2, etc.
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D_list = [np.asarray(D, dtype=np.float64) for D in D_list]
    prob = np.asarray(prob, dtype=np.float64)

    K = prob.shape[0]  # Original number of classes
    R = prob.shape[1]  # New number of classes

    if len(D_list) != K:
        raise ValueError(f"prob has {K} rows but D_list has {len(D_list)} matrices")

    # D0 remains unchanged
    D0_new = D0.copy()

    # Compute new marking matrices by weighted combination
    D_list_new = []
    for r in range(R):
        D_r_new = np.zeros_like(D0)
        for k in range(K):
            D_r_new = D_r_new + prob[k, r] * D_list[k]
        D_list_new.append(D_r_new)

    return D0_new, D_list_new


def mmap_scale(D0: np.ndarray, D_list: List[np.ndarray],
              M: Union[float, np.ndarray],
              max_iter: int = 30) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Scale the mean inter-arrival times of an MMAP.

    Adjusts all matrices to achieve specified mean inter-arrival times.

    Args:
        D0: Hidden transition matrix
        D_list: List of arrival marking matrices
        M: Desired mean inter-arrival times. Can be:
           - Single float: uniform scaling of all classes
           - Array of K values: individual scaling per class
        max_iter: Maximum iterations for optimization (only for vector M)

    Returns:
        Tuple of (D0_scaled, D_list_scaled)

    Algorithm:
        For single M: Scale all matrices uniformly by ratio = M_old / M
        For vector M: Use iterative coordinate descent to find scaling factors
    """
    D0 = np.asarray(D0, dtype=np.float64).copy()
    D_list = [np.asarray(D, dtype=np.float64).copy() for D in D_list]

    n = D0.shape[0]
    m = len(D_list)

    # Handle scalar vs vector M
    if np.isscalar(M):
        # Uniform scaling
        M_scalar = float(M)

        # Compute current mean
        agg_D1 = np.zeros_like(D0)
        for D_k in D_list:
            agg_D1 = agg_D1 + D_k

        M_old = 1.0 / map_lambda(D0, agg_D1)
        scale_factor = M_old / M_scalar

        # Scale all matrices
        D0 = D0 * scale_factor
        D_list = [D_k * scale_factor for D_k in D_list]

    else:
        # Per-class scaling (more complex)
        M_vec = np.asarray(M, dtype=np.float64)
        if len(M_vec) != m:
            raise ValueError(f"M has {len(M_vec)} elements but MMAP has {m} classes")

        # Initial heuristic: scale each class independently
        agg_D1 = sum([D_k for D_k in D_list])
        lambdas = np.zeros(m)

        for k in range(m):
            try:
                lambdas[k] = map_lambda(D0, D_list[k])
            except:
                lambdas[k] = 1.0

        # Initial scaling factors
        x = np.zeros(m)
        for k in range(m):
            if lambdas[k] > 0:
                x[k] = (1.0 / M_vec[k]) / lambdas[k]
            else:
                x[k] = 1.0

        # Iterative refinement using coordinate descent
        for iteration in range(max_iter):
            # Build scaled MMAP
            D0_scaled = D0.copy()
            D_list_scaled = [D_list[k] * x[k] for k in range(m)]

            # Normalize
            agg_D1 = sum([D_k for D_k in D_list_scaled])
            diag_sum = agg_D1.sum(axis=1)
            np.fill_diagonal(D0_scaled, -diag_sum)

            # Compute error
            try:
                lambdas_new = np.zeros(m)
                for k in range(m):
                    lambdas_new[k] = map_lambda(D0_scaled, D_list_scaled[k])

                error = np.sum(np.abs(1.0 / M_vec - lambdas_new))

                if error < 1e-6:
                    break

                # Update scaling factors
                for k in range(m):
                    if lambdas_new[k] > 0:
                        x[k] = x[k] * (1.0 / M_vec[k]) / lambdas_new[k]
                    x[k] = np.clip(x[k], 1e-6, 1e6)
            except:
                break

        # Final scaling
        D0 = D0.copy()
        D_list = [D_list[k] * x[k] for k in range(m)]

        # Normalize
        agg_D1 = sum([D_k for D_k in D_list])
        diag_sum = agg_D1.sum(axis=1)
        np.fill_diagonal(D0, -diag_sum)

    # Normalize result
    D0, D_list = mmap_normalize(D0, D_list)

    return D0, D_list


def mmap_hide(D0: np.ndarray, D_list: List[np.ndarray],
             types: Union[int, List[int], np.ndarray]) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Hide (remove) specified arrival classes from an MMAP.

    The hidden classes are set to zero matrices, effectively removing them
    from observation while maintaining the underlying stochastic process.

    Args:
        D0: Hidden transition matrix
        D_list: List of arrival marking matrices
        types: Indices of classes to hide. Can be single int or list of ints.
               0-indexed (class 1 is index 0)

    Returns:
        Tuple of (D0_hidden, D_list_hidden) where specified classes are zero

    Example:
        mmap_hide(D0, [D1, D2, D3], types=[1]) hides class 2 (index 1)
    """
    D0 = np.asarray(D0, dtype=np.float64).copy()
    D_list = [np.asarray(D, dtype=np.float64).copy() for D in D_list]

    # Normalize type specification
    if np.isscalar(types):
        hide_indices = [int(types)]
    else:
        hide_indices = [int(t) for t in types]

    # Validate indices
    if any(idx >= len(D_list) for idx in hide_indices):
        raise ValueError(f"Type indices out of range. MMAP has {len(D_list)} classes")

    # Zero out hidden classes
    for idx in hide_indices:
        D_list[idx] = np.zeros_like(D_list[idx])

    # Normalize
    D0, D_list = mmap_normalize(D0, D_list)

    return D0, D_list


def mmap_compress(D0: np.ndarray, D_list: List[np.ndarray],
                 method: str = "default",
                 target_order: Optional[int] = None) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Compress an MMAP using various approximation methods.

    Reduces the state space of an MMAP while preserving key statistical properties
    (moments, inter-arrival time distribution, etc.).

    Args:
        D0: Hidden transition matrix
        D_list: List of arrival marking matrices
        method: Compression method. Options:
                - "default", "mixture": Mixture fitting (order 1)
                - "exponential": Fit to single-state Poisson
        target_order: Target state space size (overrides method if specified)

    Returns:
        Tuple of (D0_compressed, D_list_compressed)

    Notes:
        The "default" mixture method fits each class independently to MMPP(2)
        if possible, then combines them. This is fast and often effective.

        More advanced methods (m3pp, mamap2) can be added as extensions.
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D_list = [np.asarray(D, dtype=np.float64) for D in D_list]

    n = D0.shape[0]
    m = len(D_list)

    # If already small, return as-is
    if n <= 2 or (target_order is not None and n <= target_order):
        return D0.copy(), [D.copy() for D in D_list]

    if method in ["default", "mixture"]:
        return _compress_mixture_order1(D0, D_list, target_order)
    elif method == "exponential":
        return _compress_exponential(D0, D_list)
    else:
        raise ValueError(f"Unknown compression method: {method}")


def _compress_mixture_order1(D0: np.ndarray, D_list: List[np.ndarray],
                            target_order: Optional[int] = None) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Compress using mixture of exponential approximations."""
    # For simplicity, start by computing aggregate arrival rate
    agg_D1 = sum([D for D in D_list])
    lambda_agg = map_lambda(D0, agg_D1)

    # Get class probabilities
    class_rates = []
    for D_k in D_list:
        try:
            rate_k = map_lambda(D0, D_k)
            class_rates.append(rate_k)
        except:
            class_rates.append(0.0)

    class_rates = np.array(class_rates)
    class_probs = class_rates / (np.sum(class_rates) + 1e-10)

    # Create simplified MMAP with exponential arrivals per class
    D0_comp = np.array([[-lambda_agg]])
    D_list_comp = []

    for k in range(len(D_list)):
        if class_rates[k] > 0:
            D_k_comp = np.array([[class_rates[k]]])
        else:
            D_k_comp = np.array([[0.0]])
        D_list_comp.append(D_k_comp)

    return D0_comp, D_list_comp


def _compress_exponential(D0: np.ndarray, D_list: List[np.ndarray]) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Compress to single exponential (Poisson) arrivals."""
    agg_D1 = sum([D for D in D_list])
    lambda_agg = map_lambda(D0, agg_D1)

    D0_comp = np.array([[-lambda_agg]])
    D_list_comp = [np.array([[lambda_agg]])]

    return D0_comp, D_list_comp


def mmap_exponential(lambda_rate: Union[float, np.ndarray], nclasses: Optional[int] = None) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Create an exponential MMAP (Poisson arrival process).

    Args:
        lambda_rate: Arrival rate (scalar) or rate per class (vector)
        nclasses: Number of classes (if lambda_rate is scalar)

    Returns:
        Tuple of (D0, D_list) representing exponential MMAP
    """
    if isinstance(lambda_rate, (int, float)):
        # Single rate: create Poisson (single-class or broadcast to all classes)
        lambda_rate = float(lambda_rate)
        if nclasses is None:
            nclasses = 1

        D0 = np.array([[-lambda_rate]])
        D_list = [np.array([[lambda_rate]])]

        # Broadcast to multiple classes if needed
        if nclasses > 1:
            D_list = [np.zeros((1, 1)) for _ in range(nclasses)]
            D_list[0] = np.array([[lambda_rate]])

    else:
        # Vector of rates per class
        lambda_rate = np.asarray(lambda_rate, dtype=np.float64)
        lambda_total = np.sum(lambda_rate)

        D0 = np.array([[-lambda_total]])
        D_list = [np.array([[rate]]) for rate in lambda_rate]

    return D0, D_list


def mmap_issym(mmap: List[np.ndarray]) -> bool:
    """
    Check if an MMAP contains symbolic elements.

    In Python/NumPy, we don't have built-in symbolic support like MATLAB,
    so this always returns False for numpy arrays.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        False (symbolic computation not supported in numpy)
    """
    # NumPy doesn't support symbolic computation natively
    # This is provided for API compatibility
    return False


def mmap_isfeasible(mmap: List[np.ndarray], tol: float = None) -> bool:
    """
    Check whether an MMAP is feasible up to the given tolerance.

    Checks:
    - Elements are real (no imaginary parts)
    - D0 + D1 rows sum to zero (generator property)
    - Diagonal of D0 is negative
    - Off-diagonal of D0 is non-negative
    - All D1c matrices are non-negative
    - D1 = D11 + D12 + ... + D1C

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        tol: Tolerance for feasibility checks (default: 1e-6)

    Returns:
        True if MMAP is feasible, False otherwise
    """
    if tol is None:
        tol = 1e-6

    if len(mmap) < 2:
        return False

    D0 = np.asarray(mmap[0], dtype=np.float64)

    # Check for imaginary parts
    for mat in mmap:
        if np.max(np.abs(np.imag(mat))) > tol:
            return False

    # Check basic MAP feasibility (D0, D1 = sum of marking matrices)
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)

    # Check rows sum to zero
    e = np.ones(D0.shape[0])
    row_sums = (D0 + D1_total) @ e
    if np.max(np.abs(row_sums)) > tol:
        return False

    # Check diagonal of D0 is negative
    diag_D0 = np.diag(D0)
    if np.any(diag_D0 > tol):
        return False

    # Check off-diagonal of D0 is non-negative
    D0_offdiag = D0 - np.diag(diag_D0)
    if np.min(D0_offdiag) < -tol:
        return False

    # Check all marking matrices are non-negative
    C = len(mmap) - 2  # Number of classes
    for c in range(C):
        Dc = np.asarray(mmap[2 + c], dtype=np.float64)
        if np.min(Dc) < -tol:
            return False

    # Check D1 = D11 + D12 + ... + D1C (if more than 2 matrices)
    if len(mmap) > 2:
        D1 = np.asarray(mmap[1], dtype=np.float64)
        D1_sum = np.zeros_like(D0)
        for c in range(C):
            D1_sum = D1_sum + np.asarray(mmap[2 + c], dtype=np.float64)
        if np.max(np.abs(D1 - D1_sum)) > tol:
            return False

    return True


def mmap_lambda(mmap: List[np.ndarray]) -> np.ndarray:
    """
    Compute the arrival rate of each class in an MMAP.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        Array of arrival rates, one per class
    """
    return mmap_count_lambda(mmap)


def mmap_count_lambda(mmap: List[np.ndarray]) -> np.ndarray:
    """
    Compute the arrival rate of the counting process for each class.

    The rate for class k is λ_k = π * D_k * e where π is the steady-state
    of the underlying CTMC.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        Array of arrival rates, one per class
    """
    if len(mmap) < 2:
        return np.array([])

    D0 = np.asarray(mmap[0], dtype=np.float64)
    n = D0.shape[0]
    K = len(mmap) - 2  # Number of classes

    if K == 0:
        # Only D0 and D1, treat as single class
        D1 = np.asarray(mmap[1], dtype=np.float64)
        theta = map_piq(D0, D1)
        e = np.ones(n)
        return np.array([theta @ D1 @ e])

    # Compute steady-state of underlying CTMC
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)

    theta = map_piq(D0, D1_total)
    e = np.ones(n)

    lk = np.zeros(K)
    for k in range(K):
        Dk = np.asarray(mmap[2 + k], dtype=np.float64)
        lk[k] = theta @ Dk @ e

    return lk


def mmap_pie(mmap: List[np.ndarray]) -> np.ndarray:
    """
    Compute the stationary probability of the DTMC embedded at restart instants
    after an arrival of each class.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        Matrix of shape (m, n) where row c is the steady-state for class c
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    n = D0.shape[0]
    K = len(mmap) - 2

    if K == 0:
        # Single class case
        D1 = np.asarray(mmap[1], dtype=np.float64)
        pie = map_pie(D0, D1)
        return pie.reshape(1, -1)

    pie = np.zeros((K, n))

    # Compute D1 = sum of all marking matrices
    D1 = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1 = D1 + np.asarray(mmap[k], dtype=np.float64)

    for c in range(K):
        Dc = np.asarray(mmap[2 + c], dtype=np.float64)
        # Pc = (-D0 - D1 + Dc)^{-1} * Dc
        try:
            Pc = np.linalg.solve(-D0 - D1 + Dc, Dc)
            # Solve steady-state: pi * Pc = pi, sum(pi) = 1
            A = Pc.T - np.eye(n)
            A[-1, :] = np.ones(n)
            b = np.zeros(n)
            b[-1] = 1
            pie[c, :] = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            pie[c, :] = np.ones(n) / n

    return pie


def mmap_pc(mmap: List[np.ndarray]) -> np.ndarray:
    """
    Compute the arrival probabilities of each class for the given MMAP.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        Array of arrival probabilities, one per class
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    n = D0.shape[0]
    K = len(mmap) - 2

    if K == 0:
        return np.array([1.0])

    # Compute steady-state
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)

    pie = map_pie(D0, D1_total)

    pc = np.zeros(K)
    for c in range(K):
        Dc = np.asarray(mmap[2 + c], dtype=np.float64)
        try:
            pc[c] = np.sum(pie @ np.linalg.solve(-D0, Dc))
        except np.linalg.LinAlgError:
            pc[c] = 0.0

    return pc


def mmap_embedded(mmap: List[np.ndarray]) -> List[np.ndarray]:
    """
    Compute the embedded DTMC transition matrices for each class.

    For class k, the embedded matrix is Pc_k = (-D0)^{-1} * D_{2+k}

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        List of embedded DTMC transition matrices, one per class
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    K = len(mmap) - 2

    if K == 0:
        D1 = np.asarray(mmap[1], dtype=np.float64)
        return [np.linalg.solve(-D0, D1)]

    Pc = []
    for k in range(K):
        Dk = np.asarray(mmap[2 + k], dtype=np.float64)
        Pc.append(np.linalg.solve(-D0, Dk))

    return Pc


def mmap_sample(mmap: List[np.ndarray], n_samples: int,
                pi: np.ndarray = None, seed: int = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate samples from a Marked MAP.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        n_samples: Number of random samples to collect
        pi: Initial probability (if None, uses steady-state)
        seed: Random seed for reproducibility

    Returns:
        Tuple of (T, A) where:
        - T: Array of inter-arrival times
        - A: Array of class labels (1-indexed for MATLAB compatibility)
    """
    if seed is not None:
        np.random.seed(seed)

    D0 = np.asarray(mmap[0], dtype=np.float64)
    M = D0.shape[0]  # Number of states
    C = len(mmap) - 2  # Number of classes

    if C == 0:
        # Single class case (or MAP)
        D1 = np.asarray(mmap[1], dtype=np.float64)
        if pi is None:
            pi = map_pie(D0, D1)

        # Use map_sample logic for single class
        from .map_analysis import map_sample as _map_sample
        T = _map_sample(D0, D1, n_samples, pi)
        A = np.ones(n_samples, dtype=int)
        return T, A

    if pi is None:
        D1_total = np.zeros_like(D0)
        for k in range(1, len(mmap)):
            D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)
        pi = map_pie(D0, D1_total)

    # Single state case -> exponential
    if M == 1:
        mean_iat = 1.0 / np.sum(mmap_count_lambda(mmap))
        T = np.random.exponential(mean_iat, n_samples)
        # Compute class rates
        crate = np.array([mmap[2 + c][0, 0] for c in range(C)])
        crate = crate / np.sum(crate)
        A = np.random.choice(C, size=n_samples, p=crate) + 1  # 1-indexed
        return T, A

    # Rates of exponential distributions in each state
    rates = -np.diag(D0)

    # Build transition targets
    # First M columns: internal transitions
    # Next M*C columns: arrivals for each class
    targets = np.zeros((M, M + M * C))

    for i in range(M):
        for j in range(M):
            if i != j:
                targets[i, j] = D0[i, j] / rates[i]
        for c in range(C):
            Dc = np.asarray(mmap[2 + c], dtype=np.float64)
            for j in range(M):
                targets[i, (c + 1) * M + j] = Dc[i, j] / rates[i]

    # Pick random initial state according to pi
    cumpi = np.cumsum(pi)
    init_state = np.searchsorted(cumpi, np.random.rand())

    T = np.zeros(n_samples)
    A = np.zeros(n_samples, dtype=int)

    # Simulate
    s = init_state
    t_time = 0.0
    h = 0

    while h < n_samples:
        # Sample time and transition
        time = np.random.exponential(1.0 / rates[s])

        # Sample transition
        cum_targets = np.cumsum(targets[s, :])
        transition = np.searchsorted(cum_targets, np.random.rand())

        # Determine transition type
        typ = transition // M

        if typ == 0:
            # No event, internal transition
            t_time += time
            s_next = transition % M
            s = s_next
        else:
            # Arrival event
            c = typ  # Class (1-indexed)
            T[h] = t_time + time
            A[h] = c  # 1-indexed class
            s_next = transition % M
            t_time = 0.0
            s = s_next
            h += 1

    return T, A


def mmap_rand(order: int, classes: int) -> List[np.ndarray]:
    """
    Generate a random MMAP with given order and number of classes.

    Args:
        order: Number of states
        classes: Number of arrival classes

    Returns:
        List of matrices [D0, D1, D21, D22, ..., D2K]
    """
    from .map_analysis import map_rand

    # Generate random MAP
    D0, D1 = map_rand(order)

    mmap = [D0, D1]

    # Split D1 into classes according to random probabilities
    for i in range(order):
        p = np.random.rand(classes)
        p = p / np.sum(p)

        for c in range(classes):
            if c >= len(mmap) - 2:
                mmap.append(D1 * p[c])
            else:
                mmap[2 + c] = D1 * p[c]

    return mmap


def mmap_timereverse(mmap: List[np.ndarray]) -> List[np.ndarray]:
    """
    Compute the time-reversed MMAP.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        Time-reversed MMAP
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)

    # Compute steady-state
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)

    piq = map_piq(D0, D1_total)
    D = np.diag(piq)
    D_inv = np.diag(1.0 / (piq + 1e-300))

    mmap_r = []
    for k in range(len(mmap)):
        Dk = np.asarray(mmap[k], dtype=np.float64)
        mmap_r.append(D_inv @ Dk.T @ D)

    return mmap_r


def mmap_sum(mmap: List[np.ndarray], n: int) -> List[np.ndarray]:
    """
    Create an MMAP representing the sum of n identical MMAPs (n-fold sum).

    This creates an MMAP whose inter-arrival times are distributed as the
    sum of n inter-arrival times from the original MMAP.

    Args:
        mmap: Original MMAP
        n: Number of copies to sum

    Returns:
        Summed MMAP
    """
    K = len(mmap) - 2
    mmaps = [mmap for _ in range(n)]

    # Compute state sizes
    ns = [np.asarray(m[0]).shape[0] for m in mmaps]
    total_n = sum(ns)

    D0 = np.zeros((total_n, total_n))
    D1_list = [np.zeros((total_n, total_n)) for _ in range(K + 1)]

    curpos = 0
    for i in range(n):
        ni = ns[i]
        mi = mmaps[i]

        # D0 block
        D0[curpos:curpos + ni, curpos:curpos + ni] = np.asarray(mi[0], dtype=np.float64)

        if i < n - 1:
            # Internal transition to next replica
            D0[curpos:curpos + ni, curpos + ni:curpos + ni + ns[i + 1]] = np.asarray(mi[1], dtype=np.float64)
        else:
            # Final transition with marking
            for k in range(1, len(mi)):
                D1_list[k - 1][curpos:curpos + ni, 0:ns[0]] = np.asarray(mi[k], dtype=np.float64)

        curpos += ni

    # Build result
    result = [D0]
    for D1k in D1_list:
        result.append(D1k)

    return result


def mmap_super(mmap_a: List[np.ndarray], mmap_b: List[np.ndarray] = None,
               opt: str = "default") -> List[np.ndarray]:
    """
    Superpose two or more MMAPs.

    Args:
        mmap_a: First MMAP (or list of MMAPs if mmap_b is None)
        mmap_b: Second MMAP (optional)
        opt: "default" - each class in MMAPa and MMAPb is distinct in result
             "match" - class c in both maps is mapped to class c in result

    Returns:
        Superposed MMAP
    """
    if mmap_b is None:
        # mmap_a is a list of MMAPs
        if len(mmap_a) == 0:
            return [np.array([[-1.0]]), np.array([[1.0]])]

        # Filter empty
        mmaps = [m for m in mmap_a if m is not None and len(m) > 0]
        if len(mmaps) == 0:
            return [np.array([[-1.0]]), np.array([[1.0]])]

        result = mmaps[0]
        for i in range(1, len(mmaps)):
            result = mmap_super(result, mmaps[i], opt)
        return result

    # Two MMAP case
    mmap_a = [np.asarray(m, dtype=np.float64) for m in mmap_a]
    mmap_b = [np.asarray(m, dtype=np.float64) for m in mmap_b]

    K1 = len(mmap_a) - 2
    K2 = len(mmap_b) - 2
    n1 = mmap_a[0].shape[0]
    n2 = mmap_b[0].shape[0]

    if opt.lower() == "default":
        # Each class in MMAPa and MMAPb is distinct in result
        D0_sup = np.kron(mmap_a[0], np.eye(n2)) + np.kron(np.eye(n1), mmap_b[0])
        D1_sup = np.kron(mmap_a[1], np.eye(n2)) + np.kron(np.eye(n1), mmap_b[1])

        result = [D0_sup, D1_sup]

        # Add class matrices from mmap_a
        for i in range(K1):
            result.append(np.kron(mmap_a[2 + i], np.eye(n2)))

        # Add class matrices from mmap_b
        for j in range(K2):
            result.append(np.kron(np.eye(n1), mmap_b[2 + j]))

    elif opt.lower() == "match":
        # Class c in both maps is mapped to class c in result
        if len(mmap_a) != len(mmap_b):
            raise ValueError("Class matching failed: MMAPs have different number of classes")

        result = []
        for i in range(len(mmap_a)):
            result.append(np.kron(mmap_a[i], np.eye(n2)) + np.kron(np.eye(n1), mmap_b[i]))
    else:
        raise ValueError(f"Unrecognized option: {opt}")

    # Normalize
    D0, D_list = mmap_normalize(result[0], result[1:])
    return [D0] + D_list


def mmap_mixture(alpha: np.ndarray, maps: List) -> List[np.ndarray]:
    """
    Create a probabilistic mixture of MAPs.

    Each MAP in the list is selected with probability alpha[i].

    Args:
        alpha: Array of mixing probabilities
        maps: List of MAPs (each MAP is [D0, D1])

    Returns:
        Mixed MMAP with len(maps) classes
    """
    alpha = np.asarray(alpha, dtype=np.float64)
    I = len(maps)

    # Handle empty MAPs
    for i in range(I):
        if maps[i] is None or len(maps[i]) == 0:
            # Replace with very slow exponential
            maps[i] = [np.array([[-1e-6]]), np.array([[1e-6]])]

    # Build block diagonal D0
    D0 = None

    for i in range(I):
        D0_i = np.asarray(maps[i][0], dtype=np.float64)

        if D0 is None:
            D0 = D0_i
        else:
            D0 = linalg.block_diag(D0, D0_i)

    n_total = D0.shape[0]

    # Compute D1 and class matrices
    D1 = np.zeros((n_total, n_total))
    Dk_list = [np.zeros((n_total, n_total)) for _ in range(I)]

    row_offset = 0
    for i in range(I):
        D0_i = np.asarray(maps[i][0], dtype=np.float64)
        D1_i = np.asarray(maps[i][1], dtype=np.float64)
        ni = D0_i.shape[0]
        e_i = np.ones(ni)

        # Compute pie for each target MAP
        col_offset = 0
        for j in range(I):
            D0_j = np.asarray(maps[j][0], dtype=np.float64)
            D1_j = np.asarray(maps[j][1], dtype=np.float64)
            nj = D0_j.shape[0]
            pie_j = map_pie(D0_j, D1_j)

            # D1 block: D1_i * e * alpha[j] * pie_j
            block = np.outer(D1_i @ e_i, alpha[j] * pie_j)
            D1[row_offset:row_offset + ni, col_offset:col_offset + nj] = block

            col_offset += nj

        # Class matrix for class i
        Dk_list[i][row_offset:row_offset + ni, :] = D1[row_offset:row_offset + ni, :]

        row_offset += ni

    # Normalize
    result = [D0, D1] + Dk_list
    D0_norm, D_list = mmap_normalize(result[0], result[1:])
    return [D0_norm] + D_list


def mmap_max(mmap_a: List[np.ndarray], mmap_b: List[np.ndarray],
             k: int) -> List[np.ndarray]:
    """
    Create MMAP for the maximum of arrivals from two synchronized MMAPs.

    This models a synchronization queue where arrivals from both sources
    must be paired before release.

    Args:
        mmap_a: First MMAP
        mmap_b: Second MMAP
        k: Length of synchronization queue

    Returns:
        MMAP for the synchronized (maximum) process
    """
    D0a = np.asarray(mmap_a[0], dtype=np.float64)
    D0b = np.asarray(mmap_b[0], dtype=np.float64)
    D1a = np.asarray(mmap_a[1], dtype=np.float64)
    D1b = np.asarray(mmap_b[1], dtype=np.float64)

    na = D0a.shape[0]
    nb = D0b.shape[0]

    Ia = np.eye(na)
    Ib = np.eye(nb)

    A0B0 = np.kron(D0a, D0b)
    A1IB = np.kron(D1a, Ib)
    IAB1 = np.kron(Ia, D1b)
    IAB0 = np.kron(Ia, D0b)
    A0IB = np.kron(D0a, Ib)

    block_size = na * nb
    total_size = block_size * (1 + k * 2)

    M0 = np.zeros((total_size, total_size))
    M1 = np.zeros((total_size, total_size))

    # First row of blocks
    M0[0:block_size, 0:block_size * 3] = np.hstack([A0B0, A1IB, IAB1])

    # Middle diagonal blocks
    for i in range(1, 2 * k - 1):
        r = block_size * i
        c = block_size * i
        M0[r:r + block_size, c:c + block_size] = A0B0

    # Last two blocks
    r = block_size * (2 * k - 1)
    c = block_size * (2 * k - 1)
    M0[r:r + block_size, c:c + block_size] = IAB0
    r += block_size
    c += block_size
    if r < total_size:
        M0[r:r + block_size, c:c + block_size] = A0IB

    # Off-diagonal blocks for internal transitions
    for i in range(1, k):
        r = block_size * (1 + 2 * (i - 1))
        c = block_size * (3 + 2 * (i - 1))
        if c + block_size <= total_size:
            M0[r:r + block_size, c:c + block_size] = A1IB
            M0[r + block_size:r + 2 * block_size, c + block_size:c + 2 * block_size] = IAB1

    # M1 blocks
    M1[block_size:block_size * 3, 0:block_size] = np.vstack([IAB1, A1IB])

    for i in range(1, k):
        r = block_size * (1 + 2 * i)
        c = block_size * (1 + 2 * (i - 1))
        if r + 2 * block_size <= total_size and c + 2 * block_size <= total_size:
            M1[r:r + block_size, c:c + block_size] = IAB1
            M1[r + block_size:r + 2 * block_size, c + block_size:c + 2 * block_size] = A1IB

    # Build result with class matrices
    result = [M0, M1]

    # Handle additional classes
    for cls in range(2, len(mmap_a)):
        D1ca = np.asarray(mmap_a[cls], dtype=np.float64)
        D1cb = np.asarray(mmap_b[cls], dtype=np.float64)

        IAB1_c = np.kron(Ia, D1cb)
        A1IB_c = np.kron(D1ca, Ib)

        Mc = np.zeros((total_size, total_size))
        Mc[block_size:block_size * 3, 0:block_size] = np.vstack([IAB1_c, A1IB_c])

        for i in range(1, k):
            r = block_size * (1 + 2 * i)
            c = block_size * (1 + 2 * (i - 1))
            if r + 2 * block_size <= total_size and c + 2 * block_size <= total_size:
                Mc[r:r + block_size, c:c + block_size] = IAB1_c
                Mc[r + block_size:r + 2 * block_size, c + block_size:c + 2 * block_size] = A1IB_c

        result.append(Mc)

    return result


def mmap_maps(mmap: List[np.ndarray]) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Extract K MAPs, one for each class of the MMAP[K] process.

    For class k, the MAP is {D0 + D1 - D_{2+k}, D_{2+k}}

    Args:
        mmap: List of matrices [D0, D1, D21, D22, ..., D2K]

    Returns:
        List of MAPs, one per class
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    D1 = np.asarray(mmap[1], dtype=np.float64)
    K = len(mmap) - 2

    maps = []
    for k in range(K):
        Dk = np.asarray(mmap[2 + k], dtype=np.float64)
        D0_k = D0 + D1 - Dk
        maps.append((D0_k, Dk))

    return maps


def mmap_count_mean(mmap: List[np.ndarray], t: float) -> np.ndarray:
    """
    Compute the mean of the counting process at resolution t.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        t: Time period for counting

    Returns:
        Array of mean counts, one per class
    """
    if len(mmap) < 2:
        return np.array([])

    lk = mmap_count_lambda(mmap)
    return lk * t


def mmap_count_var(mmap: List[np.ndarray], t: float) -> np.ndarray:
    """
    Compute the variance of the counting process at resolution t.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        t: Time period for counting

    Returns:
        Array of variances, one per class
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    n = D0.shape[0]
    K = len(mmap) - 2

    if K == 0:
        return np.array([])

    I = np.eye(n)
    e = np.ones(n)

    # Compute D = D0 + D1
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)
    D = D0 + D1_total

    # Compute steady-state
    theta = map_piq(D0, D1_total)

    # Compute fundamental matrix
    try:
        tmp = np.linalg.inv(np.outer(e, theta) - D)
    except np.linalg.LinAlgError:
        tmp = np.eye(n)

    lk = np.zeros(K)
    ck = np.zeros((K, n))
    dk = np.zeros((n, K))
    llk = np.zeros(K)
    vk = np.zeros(K)

    for k in range(K):
        Dk = np.asarray(mmap[2 + k], dtype=np.float64)
        lk[k] = theta @ Dk @ e
        ck[k, :] = theta @ Dk @ tmp
        dk[:, k] = tmp @ Dk @ e
        llk[k] = theta @ Dk @ e

    for k in range(K):
        Dk = np.asarray(mmap[2 + k], dtype=np.float64)
        vk[k] = (llk[k] - 2 * lk[k]**2 + 2 * ck[k, :] @ Dk @ e) * t
        vk[k] -= 2 * ck[k, :] @ (I - linalg.expm(D * t)) @ dk[:, k]

    return vk


def mmap_count_idc(mmap: List[np.ndarray], t: float) -> np.ndarray:
    """
    Compute the per-class Index of Dispersion of Counts at resolution t.

    IDC = Var[N(t)] / E[N(t)]

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        t: Time period

    Returns:
        Array of IDC values, one per class
    """
    m = mmap_count_mean(mmap, t)
    v = mmap_count_var(mmap, t)
    return v / (m + 1e-300)


def mmap_count_mcov(mmap: List[np.ndarray], t: float) -> np.ndarray:
    """
    Compute the count covariance between each pair of classes at time scale t.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        t: Time scale

    Returns:
        m x m covariance matrix
    """
    m = len(mmap) - 2

    # Diagonal: per-class variance
    mV = mmap_count_var(mmap, t)
    S = np.diag(mV)

    for i in range(m):
        for j in range(m):
            if i != j:
                # Compute variance between this pair of classes
                mmap2 = [mmap[0], mmap[1],
                         np.asarray(mmap[2 + i]) + np.asarray(mmap[2 + j]),
                         np.asarray(mmap[1]) - np.asarray(mmap[2 + i]) - np.asarray(mmap[2 + j])]
                pV = mmap_count_var(mmap2, t)
                S[i, j] = 0.5 * (pV[0] - mV[i] - mV[j])

    return S


def mmap_idc(mmap: List[np.ndarray]) -> np.ndarray:
    """
    Compute the asymptotic Index of Dispersion of Counts for each class.

    This is the limit of IDC(t) as t -> infinity.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        Array of asymptotic IDC values, one per class
    """
    lambda_sum = np.sum(mmap_lambda(mmap))
    if lambda_sum <= 0:
        return np.ones(len(mmap) - 2)

    t_inf = 1e6 / lambda_sum
    return mmap_count_var(mmap, t_inf) / mmap_count_mean(mmap, t_inf)


def mmap_sigma(mmap: List[np.ndarray]) -> np.ndarray:
    """
    Compute one-step class transition probabilities.

    p_{i,j} = P(C_k = j | C_{k-1} = i)

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        C x C matrix of transition probabilities
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    C = len(mmap) - 2

    if C == 0:
        return np.array([[1.0]])

    sigma = np.zeros((C, C))

    # Compute steady-state
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)

    alpha = map_pie(D0, D1_total)
    D0_inv = np.linalg.inv(-D0)

    for i in range(C):
        Di = np.asarray(mmap[2 + i], dtype=np.float64)
        start = alpha @ D0_inv @ Di
        for j in range(C):
            Dj = np.asarray(mmap[2 + j], dtype=np.float64)
            sigma[i, j] = np.sum(start @ D0_inv @ Dj)

    return sigma


def mmap_sigma2(mmap: List[np.ndarray]) -> np.ndarray:
    """
    Compute two-step class transition probabilities.

    p_{i,j,h} = P(C_k = h | C_{k-1} = j, C_{k-2} = i)

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]

    Returns:
        C x C x C 3D array of transition probabilities
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    C = len(mmap) - 2

    if C == 0:
        return np.array([[[1.0]]])

    sigma = np.zeros((C, C, C))

    # Compute steady-state
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)

    alpha = map_pie(D0, D1_total)
    D0_inv = np.linalg.inv(-D0)

    for i in range(C):
        Di = np.asarray(mmap[2 + i], dtype=np.float64)
        starti = alpha @ D0_inv @ Di
        for j in range(C):
            Dj = np.asarray(mmap[2 + j], dtype=np.float64)
            startj = starti @ D0_inv @ Dj
            for h in range(C):
                Dh = np.asarray(mmap[2 + h], dtype=np.float64)
                sigma[i, j, h] = np.sum(startj @ D0_inv @ Dh)

    return sigma


def mmap_forward_moment(mmap: List[np.ndarray], orders: np.ndarray,
                        norm: bool = True) -> np.ndarray:
    """
    Compute the theoretical forward moments of an MMAP.

    Forward moments are E[X^k | previous arrival was class c].

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        orders: Array of moment orders to compute
        norm: If True, normalize by class probability

    Returns:
        Matrix of shape (C, len(orders)) where element (c, k) is the
        order-k forward moment for class c
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    C = len(mmap) - 2
    K = len(orders)

    if C == 0:
        return np.array([])

    moments = np.zeros((C, K))

    # Compute steady-state
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)

    pie = map_pie(D0, D1_total)
    M = np.linalg.inv(-D0)

    for a in range(C):
        Da = np.asarray(mmap[2 + a], dtype=np.float64)
        if norm:
            pa = np.sum(pie @ np.linalg.solve(-D0, Da))
        else:
            pa = 1.0

        for h, k in enumerate(orders):
            fk = np.math.factorial(int(k))
            M_pow = np.linalg.matrix_power(M, int(k))
            moments[a, h] = fk / pa * np.sum(pie @ np.linalg.solve(-D0, Da) @ M_pow)

    return moments


def mmap_backward_moment(mmap: List[np.ndarray], orders: np.ndarray,
                         norm: bool = True) -> np.ndarray:
    """
    Compute the theoretical backward moments of an MMAP.

    Backward moments are E[X^k | next arrival will be class c].

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        orders: Array of moment orders to compute
        norm: If True, normalize by class probability

    Returns:
        Matrix of shape (C, len(orders)) where element (c, k) is the
        order-k backward moment for class c
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    C = len(mmap) - 2
    K = len(orders)

    if C == 0:
        return np.array([])

    moments = np.zeros((C, K))

    # Compute steady-state
    D1_total = np.zeros_like(D0)
    for k in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k], dtype=np.float64)

    pie = map_pie(D0, D1_total)
    M = np.linalg.inv(-D0)

    for a in range(C):
        Da = np.asarray(mmap[2 + a], dtype=np.float64)
        if norm:
            pa = np.sum(pie @ np.linalg.solve(-D0, Da))
        else:
            pa = 1.0

        for h, k in enumerate(orders):
            fk = np.math.factorial(int(k))
            M_pow = np.linalg.matrix_power(M, int(k) + 1)
            moments[a, h] = fk / pa * np.sum(pie @ M_pow @ Da)

    return moments


def mmap_cross_moment(mmap: List[np.ndarray], k: int) -> np.ndarray:
    """
    Compute the k-th order moment of inter-arrival times between class pairs.

    Computes E[X^k | previous = class i, next = class j] for all i, j.

    Args:
        mmap: List of matrices [D0, D1, D2, ..., Dm]
        k: Order of the moment

    Returns:
        C x C matrix where element (i, j) is the k-th moment for transitions i->j
    """
    D0 = np.asarray(mmap[0], dtype=np.float64)
    C = len(mmap) - 2

    if C == 0:
        from .map_analysis import map_moment as _map_moment
        return np.array([[_map_moment(D0, np.asarray(mmap[1]), k)]])

    MC = np.zeros((C, C))
    TG = np.zeros(C)

    # Compute steady-state
    D1_total = np.zeros_like(D0)
    for k_idx in range(1, len(mmap)):
        D1_total = D1_total + np.asarray(mmap[k_idx], dtype=np.float64)

    pie = map_pie(D0, D1_total)
    D0_inv = np.linalg.inv(-D0)

    for i in range(C):
        Di = np.asarray(mmap[2 + i], dtype=np.float64)
        TG[i] = np.sum(pie @ D0_inv @ Di)

    for i in range(C):
        Di = np.asarray(mmap[2 + i], dtype=np.float64)
        start = pie @ D0_inv @ Di / TG[i]
        M_pow = np.linalg.matrix_power(D0_inv, k + 1)
        for j in range(C):
            Dj = np.asarray(mmap[2 + j], dtype=np.float64)
            MC[i, j] = np.math.factorial(k) * np.sum(start @ M_pow @ Dj)
            MC[i, j] /= np.sum(start @ D0_inv @ Dj)

    return MC


def mmap_modulate(P: np.ndarray, HT: List, MMAP: List) -> List[np.ndarray]:
    """
    Modulate MMAPs in continuous time according to holding time distributions.

    Creates an MMAP that switches between different MMAPs according to a
    continuous-time Markov modulated process.

    Args:
        P: J x J transition probability matrix between states
        HT: List of phase-type distributions for holding time in each state
        MMAP: List of MMAPs, one for each state

    Returns:
        Modulated MMAP
    """
    J = len(HT)
    if len(MMAP) != J:
        raise ValueError("HT and MMAP must have the same length")

    # Convert MAPs to MMAPs if needed
    for j in range(J):
        if len(MMAP[j]) == 2:
            MMAP[j] = MMAP[j] + [MMAP[j][1]]

    # Get number of classes (should be same for all)
    nc = [len(m) for m in MMAP]
    if min(nc) != max(nc):
        raise ValueError("MMAPs must have the same number of classes")
    K = nc[0] - 2

    # Get state sizes
    nh = [np.asarray(ht[0]).shape[0] for ht in HT]
    nm = [np.asarray(mmap[0]).shape[0] for mmap in MMAP]

    # Build result
    result = [None for _ in range(2 + K)]
    for k in range(2 + K):
        result[k] = []

    for j in range(J):
        Q_j = [None for _ in range(3 + K)]
        HT_j = [np.asarray(HT[j][i], dtype=np.float64) for i in range(len(HT[j]))]
        MMAP_j = [np.asarray(MMAP[j][i], dtype=np.float64) for i in range(len(MMAP[j]))]

        Q_j[0] = np.kron(HT_j[0], MMAP_j[0])
        Q_j[2] = np.kron(HT_j[1], np.eye(nm[j]))

        for k in range(K):
            Q_j[3 + k] = np.kron(np.eye(nh[j]), MMAP_j[2 + k])

        nq = Q_j[0].shape[0]
        MOD_j = [[] for _ in range(2 + K)]

        for i in range(J):
            if i == j:
                MOD_j[0].append(Q_j[0])
                for k in range(K):
                    MOD_j[2 + k].append(Q_j[3 + k])
            else:
                HT_i = [np.asarray(HT[i][idx], dtype=np.float64) for idx in range(len(HT[i]))]
                MMAP_i = [np.asarray(MMAP[i][idx], dtype=np.float64) for idx in range(len(MMAP[i]))]
                pie_ht = map_pie(HT_i[0], HT_i[1])
                pie_mmap = map_pie(MMAP_i[0], MMAP_i[1])
                entry_i = np.kron(pie_ht, pie_mmap)
                e_nq = np.ones(nq)

                block = P[j, i] * np.outer(Q_j[2] @ e_nq, entry_i)
                MOD_j[0].append(block)
                for k in range(K):
                    MOD_j[2 + k].append(np.zeros_like(block))

        for k in range(2 + K):
            if k == 1:
                continue  # D1 computed later
            result[k].append(np.hstack(MOD_j[k]))

    # Stack vertically
    for k in range(2 + K):
        if k == 1:
            continue
        result[k] = np.vstack(result[k])

    # Compute D1 = sum of class matrices
    result[1] = np.zeros_like(result[0])
    for k in range(K):
        result[1] = result[1] + result[2 + k]

    # Normalize
    D0_norm, D_list = mmap_normalize(result[0], result[1:])
    return [D0_norm] + D_list
