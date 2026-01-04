"""
NPFQN Traffic Merging and Splitting Algorithms.

Implements traffic flow operations for non-product-form queueing networks
including merging, splitting, and class switching capabilities.

Key functions:
    npfqn_traffic_merge: Merge multiple MMAP traffic flows
    npfqn_traffic_merge_cs: Merge with class switching
    npfqn_traffic_split_cs: Split traffic with class switching

References:
    Casale, G., et al. "LINE: A unified library for queueing network modeling."
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Union
import warnings

from ..mam import (
    mmap_normalize,
    mmap_super_safe,
    mmap_mark,
    mmap_compress,
    map_isfeasible,
)


# Type alias for MMAP: List of numpy arrays [D0, D1, D2, ...]
MMAP = List[np.ndarray]


def _is_mmap_empty(mmap: Optional[MMAP]) -> bool:
    """Check if an MMAP is empty or None."""
    if mmap is None:
        return True
    if len(mmap) == 0:
        return True
    if mmap[0] is None or mmap[0].size == 0:
        return True
    return False


def _mmap_to_tuple(mmap: MMAP) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Convert MMAP list representation to (D0, D_list) tuple."""
    D0 = mmap[0]
    D_list = mmap[1:] if len(mmap) > 1 else []
    return D0, D_list


def _tuple_to_mmap(D0: np.ndarray, D_list: List[np.ndarray]) -> MMAP:
    """Convert (D0, D_list) tuple to MMAP list representation."""
    return [D0] + D_list


def npfqn_traffic_merge(
    MMAPa: Dict[int, Optional[MMAP]],
    config_merge: str = "default",
    config_compress: Optional[str] = None
) -> Optional[MMAP]:
    """
    Merge multiple MMAP traffic flows.

    Combines multiple MMAPs using specified aggregation strategies.
    Supports various merge configurations for different network topologies.

    Args:
        MMAPa: Dictionary of MMAP traffic flows to be merged.
               Keys are integer indices, values are MMAP lists [D0, D1, ...].
        config_merge: Merge configuration. Options:
            - "default", "super": Use MMAP superposition
            - "mixture": Apply mixture fitting after superposition
            - "interpos": Interposition method (falls back to super)
        config_compress: Compression configuration. Options:
            - None, "none": No compression
            - "default": Apply default compression

    Returns:
        Merged and normalized MMAP, or None if input is empty.

    Raises:
        RuntimeError: If unsupported merge configuration is provided.

    Example:
        >>> mmap1 = [np.array([[-1.0]]), np.array([[1.0]])]
        >>> mmap2 = [np.array([[-2.0]]), np.array([[2.0]])]
        >>> result = npfqn_traffic_merge({0: mmap1, 1: mmap2})
    """
    # Filter out empty/None MMAPs
    valid_mmaps = []
    for key in sorted(MMAPa.keys()):
        mmap = MMAPa[key]
        if not _is_mmap_empty(mmap):
            valid_mmaps.append(mmap)

    n = len(valid_mmaps)

    if n == 0:
        # Return empty/default MMAP
        D0 = np.array([[0.0]])
        D_list = [np.array([[0.0]])]
        D0_norm, D_list_norm = mmap_normalize(D0, D_list)
        return _tuple_to_mmap(D0_norm, D_list_norm)

    if n == 1:
        D0, D_list = _mmap_to_tuple(valid_mmaps[0])
        D0_norm, D_list_norm = mmap_normalize(D0, D_list)
        return _tuple_to_mmap(D0_norm, D_list_norm)

    # Validate merge configuration
    if config_merge not in ["default", "super", "mixture", "interpos"]:
        raise RuntimeError(f"Unsupported configuration for merge: {config_merge}")

    # Collect all (D0, D_list) tuples for superposition
    mmap_tuples = []
    for mmap in valid_mmaps:
        D0, D_list = _mmap_to_tuple(mmap)
        mmap_tuples.append((D0, D_list))

    # Perform merging based on configuration
    if config_merge in ["default", "super"]:
        # Use safe superposition
        D0_result, D_list_result = mmap_super_safe(mmap_tuples)

    elif config_merge == "mixture":
        # Superpose then apply mixture fitting
        D0_result, D_list_result = mmap_super_safe(mmap_tuples)

        # Check if result is feasible before mixture fitting
        try:
            if map_isfeasible(D0_result, sum(D_list_result)):
                # Use compression with mixture method as approximation
                # (mmap_mixture_fit_mmap is not available in native)
                D0_result, D_list_result = mmap_compress(D0_result, D_list_result, method="mixture")
        except Exception:
            pass  # Keep superposition result if fitting fails

    elif config_merge == "interpos":
        # Interposition method not yet ported - fall back to super
        warnings.warn(
            "npfqn_traffic_merge: 'interpos' merge mode requires m3pp2m functions "
            "that are not yet ported. Using 'super' merge instead.",
            RuntimeWarning
        )
        D0_result, D_list_result = mmap_super_safe(mmap_tuples)

    # Apply compression if requested
    if config_compress == "default":
        D0_result, D_list_result = mmap_compress(D0_result, D_list_result)

    # Normalize result
    D0_norm, D_list_norm = mmap_normalize(D0_result, D_list_result)

    return _tuple_to_mmap(D0_norm, D_list_norm)


def npfqn_traffic_split_cs(
    MMAP_input: MMAP,
    P: np.ndarray
) -> Dict[int, Optional[MMAP]]:
    """
    Split MMAP traffic flows with class switching.

    Decomposes a single MMAP into multiple MMAPs based on routing
    probabilities and class switching matrix.

    Args:
        MMAP_input: Input MMAP as list [D0, D1, D2, ...].
        P: Class switching probability matrix (R x J) where:
           - R is the number of arrival classes
           - J = M * R where M is the number of destinations
           - P[r, (jst-1)*R + s] = probability that class r arrival
             becomes class s at destination jst

    Returns:
        Dictionary mapping destination index to split MMAP.
        Keys are 0-indexed destination indices.

    Algorithm:
        1. Parse dimensions from P matrix
        2. For each destination, create MMAP with weighted marking matrices
        3. Normalize each result

    Example:
        >>> mmap = [np.array([[-3.0]]), np.array([[1.0]]), np.array([[2.0]])]
        >>> P = np.array([[0.8, 0.2, 0.5, 0.5], [0.3, 0.7, 0.6, 0.4]])
        >>> result = npfqn_traffic_split_cs(mmap, P)
    """
    P = np.asarray(P, dtype=np.float64)
    R = P.shape[0]  # Number of classes
    J = P.shape[1]  # Total columns = M * R
    M = int(np.round(J / R))  # Number of destinations

    D0 = np.asarray(MMAP_input[0], dtype=np.float64)
    n = D0.shape[0]  # State space size

    # Arrival matrices for each class (indices 1 to R in MMAP)
    D_arrivals = []
    for r in range(R):
        if r + 1 < len(MMAP_input):
            D_arrivals.append(np.asarray(MMAP_input[r + 1], dtype=np.float64))
        else:
            D_arrivals.append(np.zeros((n, n)))

    # Also need D1 (total arrivals) from original MMAP
    if len(MMAP_input) > 1:
        D1_total = np.asarray(MMAP_input[1], dtype=np.float64).copy()
    else:
        D1_total = np.zeros((n, n))

    SMMAP: Dict[int, Optional[MMAP]] = {}

    for jst in range(M):
        # Initialize split MMAP for this destination
        # D0_split = D0 + D1 (hidden + arrivals become hidden)
        D0_split = D0.copy() + D1_total.copy()
        D1_split = np.zeros((n, n))

        # Marking matrices for each class at this destination
        D_list_split = []
        for s in range(R):
            D_s = np.zeros((n, n))
            for r in range(R):
                # Weight by class switching probability
                # P[r, jst * R + s] = prob class r becomes class s at dest jst
                if jst * R + s < J:
                    weight = P[r, jst * R + s]
                else:
                    weight = 0.0
                D_s = D_s + weight * D_arrivals[r]

            D_list_split.append(D_s)
            D1_split = D1_split + D_s
            D0_split = D0_split + D_s

        # Normalize
        try:
            D0_norm, D_list_norm = mmap_normalize(D0_split, D_list_split)
            SMMAP[jst] = _tuple_to_mmap(D0_norm, D_list_norm)
        except Exception:
            SMMAP[jst] = None

    return SMMAP


def npfqn_traffic_merge_cs(
    MMAPs: Dict[int, MMAP],
    prob: np.ndarray,
    config: str = "default"
) -> Optional[MMAP]:
    """
    Merge MMAP traffic flows with class switching.

    Combines multiple MMAPs while applying class switching transformations
    based on the probability matrix.

    Args:
        MMAPs: Dictionary of MMAP traffic flows indexed by source.
        prob: Class switching probability matrix ((n*R) x R) where:
              - n is the number of sources
              - R is the number of classes
              - prob[(i-1)*R + r, s] = probability that class r from source i
                becomes class s in the merged stream
        config: Merge configuration ("default" or "super").

    Returns:
        Merged MMAP with class switching applied, or None if empty.

    Algorithm:
        1. Apply mmap_mark to each source MMAP to encode class switching
        2. Superpose all marked MMAPs
        3. Return result

    Example:
        >>> mmap1 = [np.array([[-1.0]]), np.array([[0.5]]), np.array([[0.5]])]
        >>> mmap2 = [np.array([[-2.0]]), np.array([[1.0]]), np.array([[1.0]])]
        >>> prob = np.array([[0.8, 0.2], [0.3, 0.7], [0.6, 0.4], [0.5, 0.5]])
        >>> result = npfqn_traffic_merge_cs({0: mmap1, 1: mmap2}, prob)
    """
    prob = np.asarray(prob, dtype=np.float64)
    n = len(MMAPs)
    R = prob.shape[1]  # Number of output classes

    if n == 0:
        return None

    # Create copies and apply class switching via mmap_mark
    MMAP_marked: Dict[int, Tuple[np.ndarray, List[np.ndarray]]] = {}

    for i, key in enumerate(sorted(MMAPs.keys())):
        mmap = MMAPs[key]
        D0, D_list = _mmap_to_tuple(mmap)

        # Build class switching matrix P for this source
        P_local = np.zeros((R, R))
        for r in range(R):
            for s in range(R):
                idx = i * R + r
                if idx < prob.shape[0]:
                    P_local[r, s] = prob[idx, s]

        # Apply marking with class switching
        D0_marked, D_list_marked = mmap_mark(D0, D_list, P_local)
        MMAP_marked[i] = (D0_marked, D_list_marked)

    # Handle single source case
    if n == 1:
        D0, D_list = list(MMAP_marked.values())[0]
        return _tuple_to_mmap(D0, D_list)

    # Merge using superposition
    if config in ["default", "super"]:
        mmap_tuples = list(MMAP_marked.values())
        D0_result, D_list_result = mmap_super_safe(mmap_tuples)
        return _tuple_to_mmap(D0_result, D_list_result)

    return None


__all__ = [
    'npfqn_traffic_merge',
    'npfqn_traffic_split_cs',
    'npfqn_traffic_merge_cs',
]
