"""
Station replica consolidation utilities for Product-Form Queueing Networks.

Provides functions to identify replicated stations (stations with identical demand
rows and service rates) and consolidate them into unique stations with multiplicity.
This optimization reduces computational cost for networks with symmetric station replicas.
"""

import numpy as np
from typing import Tuple, Optional, NamedTuple


class PfqnUniqueResult(NamedTuple):
    """Result class for pfqn_unique containing all output matrices and mapping information.

    Attributes:
        L_unique: Reduced demand matrix (M' x R) with M' <= M unique stations
        mu_unique: Reduced load-dependent rates (M' x Ntot), None if mu was empty
        gamma_unique: Reduced class-dependent rates (M' x R), None if gamma was empty
        mi: Multiplicity vector (1 x M'), mi[j] = count of stations mapping to unique station j
        mapping: Mapping vector (1 x M), mapping[i] = unique station index for original station i
    """
    L_unique: np.ndarray
    mu_unique: Optional[np.ndarray]
    gamma_unique: Optional[np.ndarray]
    mi: np.ndarray
    mapping: np.ndarray


def pfqn_unique(L: np.ndarray,
                mu: Optional[np.ndarray] = None,
                gamma: Optional[np.ndarray] = None,
                tol: float = 1e-14) -> PfqnUniqueResult:
    """
    Consolidate replicated stations into unique stations with multiplicity.

    Identifies stations with identical demand rows L[i,:] and (if present)
    identical load-dependent rates mu[i,:] or class-dependent rates gamma[i,:].
    Returns reduced matrices with only unique stations plus a multiplicity vector.

    Args:
        L: Service demand matrix (M x R)
        mu: Load-dependent rate matrix (M x Ntot), optional - pass None if not used
        gamma: Class-dependent service rate matrix (M x R), optional - pass None if not used
        tol: Tolerance for floating point comparison (default 1e-14)

    Returns:
        PfqnUniqueResult containing reduced matrices and mapping information
    """
    L = np.asarray(L, dtype=np.float64)
    if L.ndim == 1:
        L = L.reshape(-1, 1)

    M, R = L.shape

    # Build fingerprint matrix for comparison
    fingerprint = L.copy()
    if mu is not None:
        mu = np.asarray(mu, dtype=np.float64)
        fingerprint = np.hstack([fingerprint, mu])
    if gamma is not None:
        gamma = np.asarray(gamma, dtype=np.float64)
        fingerprint = np.hstack([fingerprint, gamma])

    # Find unique rows using tolerance-based comparison
    mapping = np.full(M, -1, dtype=int)
    unique_idx = []
    mi_list = []

    for i in range(M):
        if mapping[i] == -1:  # not yet assigned
            # This is a new unique station
            unique_idx.append(i)
            group_idx = len(unique_idx) - 1
            mapping[i] = group_idx
            count = 1

            # Find all stations identical to this one
            for j in range(i + 1, M):
                if mapping[j] == -1:  # not yet assigned
                    # Compute infinity norm of difference
                    max_diff = np.max(np.abs(fingerprint[i, :] - fingerprint[j, :]))
                    if max_diff < tol:
                        mapping[j] = group_idx
                        count += 1

            mi_list.append(count)

    # Extract unique rows
    M_unique = len(unique_idx)
    L_unique = L[unique_idx, :]

    mu_unique = None
    if mu is not None:
        mu_unique = mu[unique_idx, :]

    gamma_unique = None
    if gamma is not None:
        gamma_unique = gamma[unique_idx, :]

    # Create multiplicity vector
    mi = np.array(mi_list, dtype=np.float64).reshape(1, -1)

    return PfqnUniqueResult(L_unique, mu_unique, gamma_unique, mi, mapping)


def pfqn_expand(QN: np.ndarray,
                UN: np.ndarray,
                CN: np.ndarray,
                mapping: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Expand per-station metrics from reduced model to original dimensions.

    Expands performance metrics computed on a reduced model (with unique stations)
    back to the original model dimensions by replicating values according to mapping.

    Args:
        QN: Queue lengths from reduced model (M' x R)
        UN: Utilizations from reduced model (M' x R)
        CN: Cycle times from reduced model (M' x R)
        mapping: Mapping vector from pfqn_unique (length M), mapping[i] = unique station index

    Returns:
        Tuple of (QN_full, UN_full, CN_full) in original dimensions (M x R)
    """
    M_original = len(mapping)
    R = QN.shape[1] if QN.ndim > 1 else 1

    QN_full = np.zeros((M_original, R))
    UN_full = np.zeros((M_original, R))
    CN_full = np.zeros((M_original, R))

    for i in range(M_original):
        unique_idx = mapping[i]
        QN_full[i, :] = QN[unique_idx, :] if QN.ndim > 1 else QN[unique_idx]
        UN_full[i, :] = UN[unique_idx, :] if UN.ndim > 1 else UN[unique_idx]
        CN_full[i, :] = CN[unique_idx, :] if CN.ndim > 1 else CN[unique_idx]

    return QN_full, UN_full, CN_full


def pfqn_combine_mi(mi: np.ndarray,
                    mapping: np.ndarray,
                    M_unique: int) -> np.ndarray:
    """
    Combine user-provided multiplicity vector with detected replica multiplicity.

    For each unique station j, sums the mi values of all original stations mapping to it.

    Args:
        mi: User-provided multiplicity vector (1 x M_original or M_original,)
        mapping: Mapping vector from pfqn_unique (length M_original)
        M_unique: Number of unique stations

    Returns:
        Combined multiplicity vector (1 x M_unique)
    """
    mi = np.asarray(mi, dtype=np.float64).flatten()
    mi_combined = np.zeros(M_unique)

    for i in range(len(mapping)):
        mi_combined[mapping[i]] += mi[i]

    return mi_combined.reshape(1, -1)


__all__ = [
    'PfqnUniqueResult',
    'pfqn_unique',
    'pfqn_expand',
    'pfqn_combine_mi',
]
