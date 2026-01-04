"""
Distribution conversion utilities for AoI solver.

Converts LINE distribution objects and other formats into Phase-Type (PH)
representations suitable for aoi-fluid solvers.
"""

import numpy as np
from typing import Tuple, Union, Any


def aoi_dist2ph(dist_proc: Any) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert LINE distribution to Phase-Type (PH) representation.

    Converts various distribution formats to (sigma, S) representation where:
    - sigma: Initial probability vector (row vector, shape (k,))
    - S: Sub-generator matrix (infinitesimal generator, shape (k, k))

    **Supported Distributions**:
    - Exponential: single phase with rate λ
    - HyperExponential: mixture of exponential phases
    - Erlang: sequential exponential phases
    - Coxian: general phase-type with feedback
    - PH/MAP: direct extraction from process definition

    Parameters
    ----------
    dist_proc : object
        Distribution object from LINE or BUTools

    Returns
    -------
    sigma : np.ndarray
        Initial probability vector (shape (k,))
        Represents starting phase distribution
    S : np.ndarray
        Sub-generator matrix (shape (k, k))
        Infinitesimal generator with S_ii < 0 and S_ij ≥ 0 for i ≠ j

    Raises
    ------
    ValueError
        If distribution type is not supported or cannot be converted

    Examples
    --------
    >>> from butools.dph import DPH
    >>> dph = DPH(alpha=[1, 0], A=[[-2, 2], [0, -3]])
    >>> sigma, S = aoi_dist2ph(dph)
    >>> sigma
    array([1., 0.])
    >>> S
    array([[-2.,  2.],
           [ 0., -3.]])

    >>> # For exponential distribution
    >>> from butools.dph import Exp
    >>> exp_dist = Exp(rate=2.0)
    >>> sigma, S = aoi_dist2ph(exp_dist)
    >>> sigma
    array([1.])
    >>> S
    array([[-2.]])
    """
    # Try to extract PH representation directly
    if hasattr(dist_proc, 'alpha') and hasattr(dist_proc, 'A'):
        # BUTools PH/Coxian format: has alpha (initial vector) and A (sub-generator)
        sigma = np.asarray(dist_proc.alpha, dtype=float)
        S = np.asarray(dist_proc.A, dtype=float)
        return sigma, S

    if hasattr(dist_proc, 'D0') and hasattr(dist_proc, 'D1'):
        # MAP/BMAP format: convert to PH by using D0 as sub-generator
        sigma = np.asarray(dist_proc.initial_state, dtype=float) if hasattr(dist_proc, 'initial_state') else None
        if sigma is None:
            # Default: uniform initial distribution
            n = dist_proc.D0.shape[0]
            sigma = np.ones(n) / n
        S = np.asarray(dist_proc.D0, dtype=float)
        return sigma, S

    # Handle LINE distribution formats
    proc_type = type(dist_proc).__name__

    if proc_type in ['Exponential', 'Exp']:
        # Exponential distribution: single phase with rate = 1/mean
        mean = dist_proc.mean if hasattr(dist_proc, 'mean') else 1.0
        rate = 1.0 / mean if mean > 0 else 1.0
        sigma = np.array([1.0])
        S = np.array([[-rate]])
        return sigma, S

    if proc_type in ['Erlang', 'Erla']:
        # Erlang(k, λ): k phases in series, each with rate kλ
        k = dist_proc.k if hasattr(dist_proc, 'k') else 2
        mean = dist_proc.mean if hasattr(dist_proc, 'mean') else 1.0
        rate = k / mean if mean > 0 else k

        # k sequential phases
        sigma = np.zeros(k)
        sigma[0] = 1.0  # Start in first phase

        S = np.zeros((k, k))
        for i in range(k):
            S[i, i] = -rate
            if i < k - 1:
                S[i, i + 1] = rate  # Forward transition

        return sigma, S

    if proc_type in ['HyperExponential', 'HyperExp', 'HExp']:
        # HyperExponential: mixture of exponentials
        # Try LINE native format first
        if hasattr(dist_proc, 'probs') and hasattr(dist_proc, 'rates'):
            probs = np.asarray(dist_proc.probs, dtype=float)
            rates = np.asarray(dist_proc.rates, dtype=float)
            n = len(probs)

            sigma = probs / probs.sum()  # Normalize
            S = -np.diag(rates)  # Diagonal sub-generator

            return sigma, S

        # Fallback: try with weights
        if hasattr(dist_proc, 'weights') and hasattr(dist_proc, 'rates'):
            weights = np.asarray(dist_proc.weights, dtype=float)
            rates = np.asarray(dist_proc.rates, dtype=float)
            n = len(weights)

            sigma = weights / weights.sum()  # Normalize
            S = -np.diag(rates)  # Diagonal sub-generator

            return sigma, S

    if proc_type in ['PH', 'Coxian', 'Cox']:
        # Already in PH form, try to extract alpha and A
        if hasattr(dist_proc, 'alpha') and hasattr(dist_proc, 'A'):
            sigma = np.asarray(dist_proc.alpha, dtype=float)
            S = np.asarray(dist_proc.A, dtype=float)
            return sigma, S

    # Fallback: Try to convert using mean and squared coefficient of variation
    if hasattr(dist_proc, 'mean'):
        mean = dist_proc.mean
        cv = dist_proc.cv if hasattr(dist_proc, 'cv') else 1.0
        scv = cv ** 2  # Squared coefficient of variation

        if mean <= 0:
            raise ValueError(f"Mean must be positive, got {mean}")

        # Use moment-based 2-phase Erlang-like distribution
        # For SCV ≈ 0.5: Erlang(2), SCV = 1: Exponential, SCV > 1: Hyperexponential
        if scv < 0.51:
            # Erlang(2) approximation
            k = 2
            rate = k / mean
            sigma = np.array([1.0, 0.0])
            S = np.array([[-rate, rate], [0.0, -rate]])
        else:
            # Exponential (single phase)
            rate = 1.0 / mean
            sigma = np.array([1.0])
            S = np.array([[-rate]])

        return sigma, S

    raise ValueError(f"Cannot convert distribution type {proc_type} to PH representation. "
                     f"Ensure distribution has 'mean', 'alpha'/'A', or 'D0'/'D1' attributes.")


def ph_validate(sigma: np.ndarray, S: np.ndarray, tolerance: float = 1e-10) -> bool:
    """
    Validate Phase-Type (PH) representation.

    Checks that (sigma, S) form a valid PH distribution:
    - sigma: probability vector (non-negative, sums to ≤ 1)
    - S: infinitesimal generator (S_ii < 0, S_ij ≥ 0 for i≠j)
    - S0 exit rates: 0 = -S*ones(k) (mass conservation)

    Parameters
    ----------
    sigma : np.ndarray
        Initial probability vector
    S : np.ndarray
        Sub-generator matrix
    tolerance : float
        Numerical tolerance for checks

    Returns
    -------
    valid : bool
        True if (sigma, S) is a valid PH representation
    """
    # Check dimensions
    if sigma.shape[0] != S.shape[0] or S.shape[0] != S.shape[1]:
        return False

    # Check sigma: non-negative
    if np.any(sigma < -tolerance):
        return False

    # Check sigma: sums to at most 1
    if np.sum(sigma) > 1.0 + tolerance:
        return False

    # Check S diagonal: negative
    if np.any(np.diag(S) > tolerance):
        return False

    # Check S off-diagonal: non-negative
    off_diag = S - np.diag(np.diag(S))
    if np.any(off_diag < -tolerance):
        return False

    # Check exit rates: -S*ones should be ≥ 0
    exit_rates = -S @ np.ones(S.shape[0])
    if np.any(exit_rates < -tolerance):
        return False

    return True
