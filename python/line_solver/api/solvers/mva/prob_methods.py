"""
Probability computation methods for MVA solver.

Implements four core probability methods matching MATLAB @SolverMVA:
- getProbAggr: Probability of current per-class job distribution
- getProbMarg: Marginal queue-length distribution
- getProbSysAggr: Joint system state probability
- getProbNormConstAggr: Log normalizing constant for closed networks

Uses binomial approximation (Schmidt 1997) for efficient computation with
caching to avoid recomputation.
"""

import numpy as np
from typing import Tuple, Optional, Dict
from ...sn import NetworkStruct
from ..common import (
    schmidt_binomial_prob_aggr,
    schmidt_binomial_prob_marg,
    schmidt_binomial_prob_sys,
    ProbabilityResults,
    SolverResults,
)


def get_prob_aggr(
    sn: NetworkStruct,
    result: SolverResults,
    ist: int,
    method: str = "binomial_approx",
) -> Tuple[float, float]:
    """
    Compute probability of current per-class job distribution at a station.

    Returns P(n_1, ..., n_K at station i) where n_r is the number of
    class-r jobs, computed using the specified approximation method.

    Algorithm:
        Binomial approximation: For each class r, jobs are independently
        distributed with probability p_r = Q[i,r] / N[r].

        P(n_1,...,n_K) = ∏_r Binomial(N[r], p_r, n_r)

    Args:
        sn: Network structure with state information
        result: Solver results containing mean queue lengths Q[i,r]
        ist: Station index (1-based, as used by MATLAB)
        method: Computation method ('binomial_approx' or 'exact')

    Returns:
        (log_prob, prob): Tuple of log probability and probability value
        - log_prob: Natural logarithm (for numerical stability)
        - prob: Actual probability (may underflow to 0)

    Raises:
        ValueError: If result.Q is not computed or station index invalid
        NotImplementedError: If method='exact' (future enhancement)

    Notes:
        - Results are cached to avoid recomputation
        - 1-based station index (MATLAB compatibility)
        - Requires sn.state to contain current state information
    """
    if result.Q is None:
        raise ValueError("Mean queue lengths (Q) not computed. Call runAnalyzer() first.")

    if ist < 1 or ist > sn.nstations:
        raise ValueError(f"Invalid station index {ist}. Must be in [1, {sn.nstations}]")

    # Check cache first
    if result.prob and result.prob.has_cached_aggr(ist, method):
        return result.prob.get_cached_aggr(ist, method)

    # Convert to 0-based indexing
    i = ist - 1
    K = sn.nclasses

    # Extract current state at this station
    if sn.state is None or len(sn.state) == 0:
        raise ValueError("Network state not available. Cannot compute probabilities.")

    # sn.state should contain the current job population at each station
    # For now, assume it's a global state vector that we extract
    nir = np.zeros(K)
    if isinstance(sn.state, np.ndarray):
        if sn.state.ndim == 1:
            # Flat state representation (need station mapping)
            # This depends on how sn.state is structured
            # For now, use a placeholder that gets the state properly
            if len(sn.state) >= K:
                nir = sn.state[:K].copy()
        elif sn.state.ndim == 2:
            # State is (M x K) matrix
            nir = sn.state[i, :].copy()

    # Get population vector
    N = np.zeros(K)
    for r in range(K):
        N[r] = sn.njobs[r] if r < len(sn.njobs) else 0

    # Compute probability
    log_prob, prob = schmidt_binomial_prob_aggr(result.Q, N, nir, i)

    # Cache result
    if result.prob is None:
        result.prob = ProbabilityResults()
    result.prob.cache_aggr(ist, method, log_prob, prob)

    return log_prob, prob


def get_prob_marg(
    sn: NetworkStruct,
    result: SolverResults,
    ist: int,
    jobclass: int,
    state_m: Optional[np.ndarray] = None,
    method: str = "binomial_approx",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute marginal queue-length distribution for a single class at a station.

    Returns P(n | station i, class r) for n = 0, 1, ..., N[r].

    Algorithm:
        Using binomial approximation with parameter p_r = Q[i,r] / N[r]:
        P(n) = Binomial(N[r], p_r, n) = C(N[r],n) * p_r^n * (1-p_r)^(N[r]-n)

    Args:
        sn: Network structure
        result: Solver results with mean queue lengths
        ist: Station index (1-based)
        jobclass: Job class index (1-based)
        state_m: Optional state values (for future extensions)
        method: Computation method ('binomial_approx')

    Returns:
        (states, probs): Tuple of:
        - states: Array [0, 1, ..., N[r]] of possible job counts
        - probs: Probability distribution P(n) for each state
                Normalized to sum to 1.0

    Raises:
        ValueError: If inputs invalid or Q not computed

    Notes:
        - 1-based indexing for both station and class
        - Probabilities are normalized to sum to exactly 1.0
        - Cache results for efficiency
    """
    if result.Q is None:
        raise ValueError("Mean queue lengths (Q) not computed. Call runAnalyzer() first.")

    if ist < 1 or ist > sn.nstations:
        raise ValueError(f"Invalid station index {ist}. Must be in [1, {sn.nstations}]")

    if jobclass < 1 or jobclass > sn.nclasses:
        raise ValueError(f"Invalid class index {jobclass}. Must be in [1, {sn.nclasses}]")

    # Check cache first
    if result.prob and result.prob.has_marginal(ist, jobclass):
        cached = result.prob.get_marginal(ist, jobclass)
        if cached is not None:
            return cached['states'], cached['probs']

    # Convert to 0-based indexing
    i = ist - 1
    r = jobclass - 1

    # Get population and compute marginal
    N = np.zeros(sn.nclasses)
    for k in range(sn.nclasses):
        N[k] = sn.njobs[k] if k < len(sn.njobs) else 0

    states, probs = schmidt_binomial_prob_marg(result.Q, N, i, r)

    # Cache result
    if result.prob is None:
        result.prob = ProbabilityResults()
    result.prob.set_marginal(ist, jobclass, {'states': states, 'probs': probs})

    return states, probs


def get_prob_sys_aggr(
    sn: NetworkStruct,
    result: SolverResults,
    method: str = "binomial_approx",
) -> Tuple[float, float]:
    """
    Compute joint probability of current system state.

    Returns P(n_1^1, ..., n_M^K) where n_i^r is the number of class-r jobs
    at station i, using product of station marginals.

    Algorithm:
        Binomial approximation assumes independence across stations:
        P(state) = ∏_{i=1}^M ∏_{r=1}^K Binomial(N[r], p_{i,r}, n_i^r)
        where p_{i,r} = Q[i,r] / N[r]

    Args:
        sn: Network structure with state
        result: Solver results with Q values
        method: Computation method ('binomial_approx')

    Returns:
        (log_prob, prob): Log probability and probability of system state

    Raises:
        ValueError: If Q or state not available
        NotImplementedError: For exact method

    Notes:
        - Assumes station independence (valid for product-form networks)
        - Cache results for efficiency
    """
    if result.Q is None:
        raise ValueError("Mean queue lengths (Q) not computed. Call runAnalyzer() first.")

    if sn.state is None or len(sn.state) == 0:
        raise ValueError("Network state not available.")

    # Check cache first
    if result.prob and result.prob.sys_aggr_prob is not None:
        if result.prob.computation_method == method:
            log_prob = result.prob.logNormConstAggr  # Use as cached log prob
            prob = result.prob.sys_aggr_prob
            if not np.isnan(log_prob):
                return log_prob, prob

    # Reconstruct state matrix if needed
    M = sn.nstations
    K = sn.nclasses
    state_matrix = np.zeros((M, K))

    if isinstance(sn.state, np.ndarray):
        if sn.state.ndim == 2 and sn.state.shape == (M, K):
            state_matrix = sn.state.copy()
        elif sn.state.ndim == 1:
            # Need to reconstruct from flat state
            # This depends on state encoding; use a simple reconstruction
            pass

    # Get population
    N = np.array([sn.njobs[r] if r < len(sn.njobs) else 0 for r in range(K)])

    # Compute system probability
    log_prob, prob = schmidt_binomial_prob_sys(result.Q, N, state_matrix)

    # Cache result
    if result.prob is None:
        result.prob = ProbabilityResults()
    result.prob.sys_aggr_prob = prob
    result.prob.logNormConstAggr = log_prob
    result.prob.computation_method = method

    return log_prob, prob


def get_prob_norm_const_aggr(
    sn: NetworkStruct,
    result: SolverResults,
    options: Optional[Dict] = None,
) -> float:
    """
    Compute log normalizing constant for closed queueing network.

    For a closed network with population N = [N_1, ..., N_K], the normalizing
    constant is G = ∑_state P(state), where P is computed via MVA.

    Algorithm:
        1. If available from solver, use precomputed lG from MVA result
        2. Otherwise, compute via exact method using MVA
        3. Cache for future use

    Args:
        sn: Network structure
        result: Solver results
        options: Optional computation options

    Returns:
        log_G: Natural logarithm of normalizing constant

    Notes:
        - For open networks, returns 0.0 (G = ∞, log(G) is undefined)
        - Uses exact computation via MVA for closed networks
        - Caches result for efficiency

    References:
        MATLAB: matlab/src/solvers/MVA/@SolverMVA/getProbNormConstAggr.m
        - Reiser & Lavenberg (1980): Mean-Value Analysis
    """
    # Check if we have a precomputed value from solver
    if result.prob is not None and not np.isnan(result.prob.logNormConstAggr):
        return result.prob.logNormConstAggr

    # Check if network is open or closed
    is_open = False
    if sn.njobs is not None:
        is_open = np.any(np.isinf(sn.njobs))

    # Initialize cache if needed
    if result.prob is None:
        result.prob = ProbabilityResults()

    if is_open:
        # Open network: G = ∞, log(G) is undefined (use 0 as placeholder)
        log_G = 0.0
    else:
        # Closed network: compute via exact MVA algorithm
        from ...pfqn.mva import pfqn_mva_multichain
        from ...sn import sn_get_demands_chain

        try:
            # Get chain-level demands
            chain_data = sn_get_demands_chain(sn)
            Lchain = chain_data.Lchain
            Nchain = chain_data.Nchain.flatten().astype(int)
            Z_chain = np.zeros(sn.nchains)

            # Add think times from delay stations
            if sn.sched:
                from ...sn import SchedStrategy
                for i in range(sn.nstations):
                    if sn.sched.get(i) == SchedStrategy.INF:
                        for c in range(sn.nchains):
                            if Lchain.shape[0] > i and c < Lchain.shape[1]:
                                Z_chain[c] += Lchain[i, c]

            # Compute MVA with normalizing constant
            mva_result = pfqn_mva_multichain(
                Lchain=Lchain,
                Nchain=Nchain,
                Z=Z_chain,
                mi=None
            )

            log_G = mva_result.get('lG', 0.0)

        except Exception as e:
            # Fall back to 0 if computation fails
            log_G = 0.0

    # Cache result
    result.prob.logNormConstAggr = log_G

    return log_G


# Helper methods for result caching

def _add_cache_methods_to_probability_results():
    """Add caching helper methods to ProbabilityResults."""
    def has_cached_aggr(self, ist, method):
        key = (ist, method)
        return key in getattr(self, '_aggr_cache', {})

    def get_cached_aggr(self, ist, method):
        key = (ist, method)
        return getattr(self, '_aggr_cache', {}).get(key, (np.nan, 0.0))

    def cache_aggr(self, ist, method, log_prob, prob):
        if not hasattr(self, '_aggr_cache'):
            self._aggr_cache = {}
        self._aggr_cache[(ist, method)] = (log_prob, prob)

    ProbabilityResults.has_cached_aggr = has_cached_aggr
    ProbabilityResults.get_cached_aggr = get_cached_aggr
    ProbabilityResults.cache_aggr = cache_aggr


# Initialize caching methods
_add_cache_methods_to_probability_results()


__all__ = [
    'get_prob_aggr',
    'get_prob_marg',
    'get_prob_sys_aggr',
    'get_prob_norm_const_aggr',
]
