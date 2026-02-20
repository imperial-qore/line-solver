"""
NC Solver probability handlers.

Native Python implementation of probability computation handlers for the NC solver.
Computes marginal and joint state probabilities for product-form queueing networks.

Port from:


"""

import numpy as np
from math import exp, log
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    sn_get_demands_chain,
)
from ...pfqn import pfqn_ncld


# Constants
FINE_TOL = 1e-12
ZERO = 1e-10


@dataclass
class StateMarginalStatistics:
    """Marginal statistics extracted from a state vector.

    Attributes:
        ni: Total number of jobs at the node (1,)
        nir: Number of jobs per class (1, K)
        sir: Number of jobs in service per class (1, K) - optional
        kir: Per-phase job distribution for each class - optional
    """
    ni: np.ndarray
    nir: np.ndarray
    sir: Optional[np.ndarray] = None
    kir: Optional[List[np.ndarray]] = None


@dataclass
class SolverNCMargReturn:
    """Result of marginal probability computation."""
    lPr: np.ndarray   # Log marginal probabilities (M, 1)
    G: float          # Normalizing constant
    lG: float         # Log normalizing constant
    runtime: float    # Runtime in seconds


@dataclass
class SolverNCJointReturn:
    """Result of joint probability computation."""
    Pr: float         # Joint probability
    G: float          # Normalizing constant
    lG: float         # Log normalizing constant
    runtime: float    # Runtime in seconds


def to_marginal_aggr(
    sn: NetworkStruct,
    ist: int,
    state_i: Optional[np.ndarray] = None
) -> StateMarginalStatistics:
    """
    Compute aggregated marginal statistics for a station.

    Extracts the number of jobs per class from the station state.
    This is a simplified version for aggregated (queue length) states.

    Args:
        sn: Network structure
        ist: Station index (0-based)
        state_i: State vector for the station (optional, uses sn.state if None)

    Returns:
        StateMarginalStatistics with job counts per class
    """
    K = sn.nclasses

    # Get state for this station
    if state_i is None:
        # Try to get state from sn.state (can be list or dict)
        if sn.state is not None:
            if isinstance(sn.state, dict):
                # Dict keyed by node index
                if sn.stationToStateful is not None and len(sn.stationToStateful) > ist:
                    isf = int(sn.stationToStateful[ist])
                    if sn.statefulToNode is not None and len(sn.statefulToNode) > isf:
                        stateful_node = int(sn.statefulToNode[isf])
                        state_i = sn.state.get(stateful_node, None)
            elif isinstance(sn.state, (list, np.ndarray)) and ist < len(sn.state):
                # List indexed by station
                state_i = sn.state[ist]

        if state_i is None:
            # Default: empty state (all zeros)
            state_i = np.zeros(K)

    # Ensure state_i is 1D or 2D
    if state_i.ndim == 1:
        state_i = state_i.reshape(1, -1)

    nrows = state_i.shape[0]
    ncols = state_i.shape[1] if state_i.ndim > 1 else len(state_i)

    # For aggregated state, assume state contains queue lengths per class
    # State format: [n_class1, n_class2, ..., n_classK, <optional phase info>]
    nir = np.zeros((nrows, K))

    for r in range(K):
        if r < ncols:
            nir[:, r] = state_i[:, r] if state_i.ndim > 1 else state_i[r]

    ni = np.sum(nir, axis=1, keepdims=True)

    return StateMarginalStatistics(
        ni=ni,
        nir=nir,
        sir=None,
        kir=None
    )


def solver_nc_margaggr(
    sn: NetworkStruct,
    options: Optional[Any] = None,
    lG: Optional[float] = None
) -> SolverNCMargReturn:
    """
    Compute aggregated marginal probabilities.

    Computes the marginal probability of observing the current queue length
    distribution at each station.

    Args:
        sn: Network structure with state information
        options: Solver options (optional)
        lG: Pre-computed log normalizing constant (optional)

    Returns:
        SolverNCMargReturn with marginal probabilities

    References:
        Port of MATLAB solver_nc_margaggr.m and Kotlin Solver_nc_margaggr.kt
    """
    start_time = time.time()

    M = sn.nstations
    K = sn.nclasses
    C = sn.nchains

    # Get server counts
    nservers = sn.nservers
    if nservers is None:
        nservers = np.ones(M)
    nservers = nservers.flatten()

    # Get chain demands and populations
    chain_result = sn_get_demands_chain(sn)
    Lchain = chain_result.Lchain
    STchain = chain_result.STchain
    Nchain = chain_result.Nchain.flatten()

    # Build visit matrix
    V = np.zeros((M, K))
    if sn.visits is not None:
        for c in range(C):
            if c in sn.visits and sn.visits[c] is not None:
                visit_mat = sn.visits[c]
                if visit_mat.shape == V.shape:
                    V += visit_mat
                else:
                    min_rows = min(visit_mat.shape[0], V.shape[0])
                    min_cols = min(visit_mat.shape[1], V.shape[1])
                    V[:min_rows, :min_cols] += visit_mat[:min_rows, :min_cols]

    # Compute service times
    rates = sn.rates
    if rates is None:
        rates = np.ones((M, K))
    with np.errstate(divide='ignore', invalid='ignore'):
        ST = np.where(rates > 0, 1.0 / rates, 0.0)
        ST = np.nan_to_num(ST, nan=0.0, posinf=0.0, neginf=0.0)

    # Build mu matrix for load-dependent scaling
    Nt = int(np.sum(Nchain))
    mu = np.zeros((M, max(1, Nt)))
    for ist in range(M):
        if np.isinf(nservers[ist]):  # Infinite server
            for j in range(Nt):
                mu[ist, j] = j + 1
        else:
            for j in range(Nt):
                mu[ist, j] = min(j + 1, nservers[ist])

    # Compute normalizing constant if not provided
    if lG is None or np.isnan(lG):
        Zchain = np.zeros(C)
        result = pfqn_ncld(Lchain, Nchain.reshape(1, -1), Zchain.reshape(1, -1), mu)
        lG = result.lG

    G = exp(lG)

    # Build chains matrix for aggregation
    chains = np.zeros((K, C))
    if sn.inchain is not None:
        for c in range(C):
            if c in sn.inchain:
                class_indices = sn.inchain[c].flatten().astype(int)
                for k in class_indices:
                    if k < K:
                        chains[k, c] = 1.0

    # Compute marginal probabilities for each station
    lPr = np.zeros((M, 1))

    for ist in range(M):
        # Get marginal statistics for this station
        marginal = to_marginal_aggr(sn, ist)
        nivec = marginal.nir

        if nivec is None or len(nivec) == 0:
            lPr[ist, 0] = 0.0
            continue

        # Check for negative population (invalid state)
        if np.any(nivec < 0):
            lPr[ist, 0] = np.nan
            continue

        # Convert to chain populations
        if nivec.ndim == 1:
            nivec = nivec.reshape(1, -1)
        nivec_chain = nivec @ chains  # (nrows, C)

        # Build reduced system (all stations except ist)
        Lchain_minus_i = np.delete(Lchain, ist, axis=0)
        mu_minus_i = np.delete(mu, ist, axis=0)

        # Reduced chain population
        Nchain_minus_i = Nchain - nivec_chain[0]

        # Compute normalizing constant for reduced system
        Zchain_minus_i = np.zeros(C)

        # Only compute if there's population left
        if np.all(Nchain_minus_i >= 0):
            if Lchain_minus_i.shape[0] > 0:
                result_minus_i = pfqn_ncld(
                    Lchain_minus_i,
                    Nchain_minus_i.reshape(1, -1),
                    Zchain_minus_i.reshape(1, -1),
                    mu_minus_i
                )
                lG_minus_i = result_minus_i.lG
            else:
                lG_minus_i = 0.0

            # Compute local normalizing constant for station ist
            ST_V_ist = (ST[ist:ist+1, :] * V[ist:ist+1, :])
            mu_ist = mu[ist:ist+1, :]
            Znivec = np.zeros_like(nivec)

            result_ist = pfqn_ncld(ST_V_ist, nivec, Znivec, mu_ist)
            lF_i = result_ist.lG

            # Marginal probability in log space
            lPr[ist, 0] = lF_i + lG_minus_i - lG
        else:
            lPr[ist, 0] = np.nan

    runtime = time.time() - start_time

    return SolverNCMargReturn(
        lPr=lPr,
        G=G,
        lG=lG,
        runtime=runtime
    )


def solver_nc_jointaggr(
    sn: NetworkStruct,
    options: Optional[Any] = None
) -> SolverNCJointReturn:
    """
    Compute aggregated joint state probability.

    Computes the joint probability of observing the current queue length
    distribution across all stations.

    Args:
        sn: Network structure with state information
        options: Solver options (optional)

    Returns:
        SolverNCJointReturn with joint probability

    References:
        Port of MATLAB solver_nc_jointaggr.m and Kotlin Solver_nc_jointaggr.kt
    """
    start_time = time.time()

    M = sn.nstations
    K = sn.nclasses
    C = sn.nchains

    # Get server counts
    nservers = sn.nservers
    if nservers is None:
        nservers = np.ones(M)
    nservers = nservers.flatten()

    # Build visit matrix
    V = np.zeros((M, K))
    if sn.visits is not None:
        for c in range(C):
            if c in sn.visits and sn.visits[c] is not None:
                visit_mat = sn.visits[c]
                if visit_mat.shape == V.shape:
                    V += visit_mat
                else:
                    min_rows = min(visit_mat.shape[0], V.shape[0])
                    min_cols = min(visit_mat.shape[1], V.shape[1])
                    V[:min_rows, :min_cols] += visit_mat[:min_rows, :min_cols]

    # Get chain demands and populations
    chain_result = sn_get_demands_chain(sn)
    Lchain = chain_result.Lchain
    Nchain = chain_result.Nchain.flatten()

    # Compute service times
    rates = sn.rates
    if rates is None:
        rates = np.ones((M, K))
    with np.errstate(divide='ignore', invalid='ignore'):
        ST = np.where(rates > 0, 1.0 / rates, 0.0)
        ST = np.nan_to_num(ST, nan=0.0, posinf=0.0, neginf=0.0)

    # Build mu matrix for load-dependent scaling
    Nt = int(np.sum(Nchain))
    mu = np.zeros((M, max(1, Nt)))
    for ist in range(M):
        if np.isinf(nservers[ist]):  # Infinite server
            for j in range(Nt):
                mu[ist, j] = j + 1
        else:
            for j in range(Nt):
                mu[ist, j] = min(j + 1, nservers[ist])

    # Compute normalizing constant
    Zchain = np.zeros(C)
    result = pfqn_ncld(Lchain, Nchain.reshape(1, -1), Zchain.reshape(1, -1), mu)
    lG = result.lG
    G = exp(lG)

    # Build chains matrix for aggregation
    chains = np.zeros((K, C))
    if sn.inchain is not None:
        for c in range(C):
            if c in sn.inchain:
                class_indices = sn.inchain[c].flatten().astype(int)
                for k in class_indices:
                    if k < K:
                        chains[k, c] = 1.0

    # Sum log probabilities over all stations
    lPr = 0.0

    for ist in range(M):
        # Get marginal statistics for this station
        marginal = to_marginal_aggr(sn, ist)
        nivec = marginal.nir

        if nivec is None or len(nivec) == 0:
            continue

        if nivec.ndim == 1:
            nivec = nivec.reshape(1, -1)

        # Get unique rows (for multi-row states)
        unique_nivec = _get_unique_rows(nivec)

        for row_idx in range(unique_nivec.shape[0]):
            current_nivec = unique_nivec[row_idx:row_idx+1, :]
            nivec_chain = current_nivec @ chains

            # Check if any population is positive
            if np.any(nivec_chain > 0):
                # Build service time matrix for this station
                ST_V_ist = ST[ist:ist+1, :] * V[ist:ist+1, :]
                mu_ist = mu[ist:ist+1, :]
                Znivec = np.zeros_like(current_nivec)

                result_ist = pfqn_ncld(ST_V_ist, current_nivec, Znivec, mu_ist)
                lF_i = result_ist.lG
                lPr += lF_i

    lPr -= lG
    Pr = exp(lPr) if np.isfinite(lPr) else 0.0

    runtime = time.time() - start_time

    return SolverNCJointReturn(
        Pr=Pr,
        G=G,
        lG=lG,
        runtime=runtime
    )


def _get_unique_rows(matrix: np.ndarray) -> np.ndarray:
    """
    Get unique rows from a matrix.

    Equivalent to MATLAB's unique(matrix, 'rows').

    Args:
        matrix: Input matrix (N, M)

    Returns:
        Matrix with unique rows
    """
    if matrix.ndim == 1:
        return matrix.reshape(1, -1)

    if matrix.shape[0] == 0:
        return matrix

    if matrix.shape[0] == 1:
        return matrix

    # Use numpy's unique with axis parameter for 2D arrays
    unique_rows, unique_indices = np.unique(matrix, axis=0, return_index=True)
    # Sort by original index to maintain order
    sorted_idx = np.argsort(unique_indices)
    return unique_rows[sorted_idx]


__all__ = [
    'StateMarginalStatistics',
    'SolverNCMargReturn',
    'SolverNCJointReturn',
    'to_marginal_aggr',
    'solver_nc_margaggr',
    'solver_nc_jointaggr',
]
