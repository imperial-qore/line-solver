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


@dataclass
class _MargAggrContext:
    """The model-wide quantities every aggregate marginal is written against.

    They depend on the network and on nothing else, so a caller that evaluates
    the same station at many occupancies -- `getProbMarg`'s enumeration over the
    per-class partitions of n -- builds this once instead of once per partition,
    which is also how `solver_nc_getprob_marg` in the C++ caches `lG`.
    """
    Lchain: np.ndarray
    mu: np.ndarray
    ST: np.ndarray
    V: np.ndarray
    chains: np.ndarray
    Nchain: np.ndarray
    lG: float


def _margaggr_context(sn: NetworkStruct, lG: Optional[float] = None) -> _MargAggrContext:
    """Assemble the chain demands, the load-dependent rates and lG."""
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

    # Build chains matrix for aggregation
    chains = np.zeros((K, C))
    if sn.inchain is not None:
        for c in range(C):
            if c in sn.inchain:
                class_indices = sn.inchain[c].flatten().astype(int)
                for k in class_indices:
                    if k < K:
                        chains[k, c] = 1.0

    return _MargAggrContext(Lchain=Lchain, mu=mu, ST=ST, V=V, chains=chains,
                            Nchain=Nchain, lG=float(lG))


def _margaggr_logp(ctx: _MargAggrContext, ist: int, nivec: np.ndarray) -> float:
    """log P(station `ist` holds `nivec` jobs per class), as lF_i + lG_{-i} - lG."""
    if nivec is None or len(nivec) == 0:
        return 0.0
    if np.any(nivec < 0):
        return float(np.nan)

    if nivec.ndim == 1:
        nivec = nivec.reshape(1, -1)
    nivec_chain = nivec @ ctx.chains  # (nrows, C)

    # Build reduced system (all stations except ist)
    Lchain_minus_i = np.delete(ctx.Lchain, ist, axis=0)
    mu_minus_i = np.delete(ctx.mu, ist, axis=0)

    # Reduced chain population
    Nchain_minus_i = ctx.Nchain - nivec_chain[0]
    if not np.all(Nchain_minus_i >= 0):
        return float(np.nan)

    # Compute normalizing constant for reduced system
    Zchain_minus_i = np.zeros(ctx.Nchain.size)
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
    ST_V_ist = (ctx.ST[ist:ist+1, :] * ctx.V[ist:ist+1, :])
    mu_ist = ctx.mu[ist:ist+1, :]
    Znivec = np.zeros_like(nivec)

    lF_i = pfqn_ncld(ST_V_ist, nivec, Znivec, mu_ist).lG

    return float(lF_i + lG_minus_i - ctx.lG)


def solver_nc_margaggr_state(
    sn: NetworkStruct,
    ist: int,
    nivec: np.ndarray,
    options: Optional[Any] = None,
    lG: Optional[float] = None
) -> Tuple[float, float]:
    """
    The aggregate marginal at ONE station for a GIVEN per-class occupancy.

    `solver_nc_margaggr` reads the occupancy off `sn.state`, which is what
    `getProbAggr` wants and what `getProbMarg` cannot use: the latter sums the
    same law over every per-class partition of a total n, none of which is the
    state the model is in. The `lG` returned beside the probability is the one
    the caller should hand back on the next partition.

    Args:
        sn: Network structure
        ist: Station index (0-based)
        nivec: Per-class job counts at that station, length nclasses
        options: Solver options (unused; kept for signature symmetry)
        lG: Pre-computed log normalizing constant (optional)

    Returns:
        (log probability, log normalizing constant)
    """
    ctx = _margaggr_context(sn, lG)
    return _margaggr_logp(ctx, int(ist), np.asarray(nivec, dtype=float)), ctx.lG


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
        Port of MATLAB solver_nc_margaggr.m
    """
    start_time = time.time()

    M = sn.nstations

    ctx = _margaggr_context(sn, lG)
    G = exp(ctx.lG)

    # Compute marginal probabilities for each station
    lPr = np.zeros((M, 1))
    for ist in range(M):
        # Get marginal statistics for this station
        lPr[ist, 0] = _margaggr_logp(ctx, ist, to_marginal_aggr(sn, ist).nir)

    runtime = time.time() - start_time

    return SolverNCMargReturn(
        lPr=lPr,
        G=G,
        lG=ctx.lG,
        runtime=runtime
    )


def solver_nc_marg(
    sn: NetworkStruct,
    options: Optional[Any] = None,
    lG: Optional[float] = None
) -> SolverNCMargReturn:
    """
    The DETAILED per-station marginal, `@SolverNC/getProb`'s own quantity.

    NOT `solver_nc_margaggr` to more places. That one evaluates the station's
    balance function with `pfqn_ncld` over the whole per-class vector, which is
    the AGGREGATE law; this one writes the balance function per class, so it
    carries the class-within-chain split that the aggregate sums out, and on a
    multichain model the two are different numbers.

    RETURNS LOGARITHMS, and the caller returns them unexponentiated: MATLAB's
    `solver_nc_marg.m` names its first output `lPr` and `@SolverNC/getProb.m`
    hands it straight back, as does `SolverNC.java`. Its sibling `getProbSys`
    returns a plain probability, so the inconsistency lives inside SolverNC and
    is reproduced here rather than quietly repaired in one codebase only.

    The discipline decides the balance function, as in the reference:

      FCFS      exponential service with one common mean; the state is weighted
                by prod_r V_ir^{n_ir} over the load-dependent denominator
      PS, INF   prod_r (V_ir S_ir)^{n_ir} / n_ir!, times |n_i|!, over the same
                denominator
      SIRO      refused: its weight needs the CLASS OF THE JOB IN SERVICE, which
                a per-class marginal does not carry
      other     left at zero, which is what the reference's switch does when no
                case matches

    References:
        Port of MATLAB solver_nc_marg.m; twin of the C++
        `line::nc::solver_nc_marg` in solver_nc_prob.h.
    """
    start_time = time.time()

    M = sn.nstations
    K = sn.nclasses
    C = sn.nchains

    nservers = sn.nservers
    if nservers is None:
        nservers = np.ones(M)
    nservers = np.asarray(nservers).flatten()

    chain_result = sn_get_demands_chain(sn)
    Lchain = chain_result.Lchain
    Nchain = chain_result.Nchain.flatten()

    V = np.zeros((M, K))
    if sn.visits is not None:
        for c in range(C):
            if c in sn.visits and sn.visits[c] is not None:
                visit_mat = sn.visits[c]
                rows = min(visit_mat.shape[0], V.shape[0])
                cols = min(visit_mat.shape[1], V.shape[1])
                V[:rows, :cols] += visit_mat[:rows, :cols]

    rates = sn.rates
    if rates is None:
        rates = np.ones((M, K))
    with np.errstate(divide='ignore', invalid='ignore'):
        ST = np.where(rates > 0, 1.0 / rates, 0.0)
        ST = np.nan_to_num(ST, nan=0.0, posinf=0.0, neginf=0.0)

    Nt = int(np.sum(Nchain))
    mu = np.zeros((M, max(1, Nt)))
    for ist in range(M):
        for j in range(max(1, Nt)):
            mu[ist, j] = (j + 1) if np.isinf(nservers[ist]) else min(j + 1, nservers[ist])

    if lG is None or np.isnan(lG):
        Zchain = np.zeros(C)
        lG = pfqn_ncld(Lchain, Nchain.reshape(1, -1), Zchain.reshape(1, -1), mu).lG

    chains = np.zeros((K, C))
    if sn.inchain is not None:
        for c in range(C):
            if c in sn.inchain:
                for k in np.asarray(sn.inchain[c]).flatten().astype(int):
                    if k < K:
                        chains[k, c] = 1.0

    def _is_exponential(ist, r):
        """One phase in the service law of (ist, r); a disabled class is vacuous."""
        try:
            if sn.phases is not None and np.asarray(sn.phases).size:
                return int(np.asarray(sn.phases)[ist, r]) <= 1
        except (IndexError, TypeError):
            pass
        return True

    lPr = np.zeros((M, 1))
    for ist in range(M):
        nivec = np.asarray(to_marginal_aggr(sn, ist).nir).flatten()
        if nivec.size == 0:
            continue
        if np.any(nivec < 0):
            lPr[ist, 0] = np.nan
            continue

        nivec_chain = nivec.reshape(1, -1) @ chains
        Nchain_minus_i = Nchain - nivec_chain[0]
        if not np.all(Nchain_minus_i >= 0):
            lPr[ist, 0] = np.nan
            continue

        if M > 1:
            lG_minus_i = pfqn_ncld(
                np.delete(Lchain, ist, axis=0),
                Nchain_minus_i.reshape(1, -1),
                np.zeros(C).reshape(1, -1),
                np.delete(mu, ist, axis=0)).lG
        else:
            lG_minus_i = 0.0

        ntot_i = int(np.sum(nivec))
        # sum_{n=1}^{|n_i|} log mu_i(n), the load-dependent denominator
        lmu = float(np.sum(np.log(mu[ist, :min(ntot_i, mu.shape[1])]))) if ntot_i > 0 else 0.0

        sched = sn.sched[ist] if sn.sched is not None else None
        sched_name = str(getattr(sched, 'name', sched)).upper()
        lF_i = 0.0
        if sched_name == 'FCFS':
            st_active = [ST[ist, r] for r in range(K) if not _disabled(sn, ist, r)]
            for r in range(K):
                if _disabled(sn, ist, r):
                    continue
                if not _is_exponential(ist, r):
                    raise ValueError(
                        "solver_nc_marg: the product-form state probability requires "
                        "exponential service times at FCFS nodes, and this station's class "
                        "%d is not exponential" % (r + 1))
            stmax = max(st_active) if st_active else 0.0
            for r in range(K):
                if _disabled(sn, ist, r) or nivec[r] == 0:
                    continue
                if abs(ST[ist, r] - stmax) > FINE_TOL:
                    raise ValueError(
                        "solver_nc_marg: the product-form state probability requires identical "
                        "service times across classes at FCFS nodes, and this station's class "
                        "%d differs" % (r + 1))
            if ntot_i > 0:
                for r in range(K):
                    if nivec[r] == 0:
                        continue
                    if not V[ist, r] > 0:
                        raise ValueError(
                            "solver_nc_marg: class %d holds jobs at a station it never visits"
                            % (r + 1))
                    lF_i += float(nivec[r]) * log(V[ist, r])
                lF_i -= lmu
        elif sched_name == 'SIRO':
            raise ValueError(
                "solver_nc_marg: the SIRO branch weighs the state by log(n_ci / sum n), which "
                "needs the CLASS OF THE JOB IN SERVICE; a per-class marginal does not carry it. "
                "Use getProbAggr, whose aggregate marginal has no such dependency")
        elif sched_name in ('PS', 'INF'):
            for r in range(K):
                if _disabled(sn, ist, r):
                    continue
                if not _is_exponential(ist, r):
                    raise ValueError(
                        "solver_nc_marg: a non-exponential service law at a %s station makes the "
                        "balance function depend on the PHASE-LEVEL occupancy, which a per-class "
                        "marginal does not carry" % ('PS' if sched_name == 'PS' else 'delay'))
                if nivec[r] == 0:
                    continue
                w = V[ist, r] * ST[ist, r]
                if not w > 0:
                    raise ValueError(
                        "solver_nc_marg: class %d holds jobs at a station whose demand for it is "
                        "zero" % (r + 1))
                lF_i += float(nivec[r]) * log(w)
                for q in range(2, int(nivec[r]) + 1):
                    lF_i -= log(float(q))
            for q in range(2, ntot_i + 1):
                lF_i += log(float(q))
            lF_i -= lmu

        lPr[ist, 0] = lF_i + lG_minus_i - lG

    return SolverNCMargReturn(lPr=lPr, G=exp(lG), lG=lG,
                              runtime=time.time() - start_time)


def _disabled(sn: NetworkStruct, ist: int, r: int) -> bool:
    """Whether class r has no service law at station ist."""
    try:
        if sn.rates is not None:
            v = np.asarray(sn.rates)[ist, r]
            return bool(np.isnan(v)) or not v > 0
    except (IndexError, TypeError):
        pass
    return False


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
        Port of MATLAB solver_nc_jointaggr.m
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

    # THE CONSTANT DEPENDS ON THE METHOD, and taking the load-dependent one
    # unconditionally made this disagree with the reference on every model with a
    # multiserver station. MATLAB's solver_nc_jointaggr switches: 'exact' takes
    # pfqn_ncld over the mu(n) = min(n, c) lattice, and every other method takes
    # solver_nc's own lG and effective service times, with the comment that it is
    # "unclear if this is correct as it doesn't consider the transformation to ld
    # model". Reproduced rather than improved: on the 2-job Delay -> PS -> PS(c=2)
    # chain MATLAB reports 0.1722498297 under the default method and this
    # function reported 0.18, so a caller comparing the codebases saw a 4.5% gap
    # that is a difference of normalizations, not of models.
    method = getattr(options, 'method', 'default') if options is not None else 'default'
    if method == 'exact':
        Zchain = np.zeros(C)
        lG = pfqn_ncld(Lchain, Nchain.reshape(1, -1), Zchain.reshape(1, -1), mu).lG
    else:
        from .handler import solver_nc, SolverOptions as _HandlerOptions
        # The handler reads knobs a SolverNCOptions does not carry (highvar), so
        # the caller's options are narrowed to the handler's own, exactly as
        # SolverNC.runAnalyzer narrows them before its own solve.
        nc_sol = solver_nc(sn, _HandlerOptions(
            method=getattr(options, 'method', 'default') or 'default',
            tol=getattr(options, 'tol', 1e-6),
            iter_max=int(getattr(options, 'iter_max', 1000) or 1000),
            iter_tol=getattr(options, 'iter_tol', 1e-4),
            verbose=0,
            samples=int(getattr(options, 'samples', 100000) or 100000),
            seed=getattr(options, 'seed', None),
        ) if options is not None else _HandlerOptions())
        lG = nc_sol.lG
        if nc_sol.STeff is not None and np.asarray(nc_sol.STeff).shape == ST.shape:
            ST = np.asarray(nc_sol.STeff, dtype=float)
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


@dataclass
class SolverNCJointMargReturn:
    """Result of a joint total-queue-length evaluation."""
    Pr: float         # Joint probability
    lPr: float        # Its logarithm
    lG: float         # Log normalizing constant, supplied or computed here
    runtime: float    # Runtime in seconds


def solver_nc_jointmarg(
    sn: NetworkStruct,
    options: Optional[Any] = None,
    nvec: np.ndarray = None,
    engine: str = 'exact',
    lG: Optional[float] = None
) -> SolverNCJointMargReturn:
    """
    Joint probability that station i holds nvec[i] jobs IN TOTAL, all classes
    summed out.

    This is NOT solver_nc_jointaggr, which fixes the per-class population of
    every station: each state here is the sum of jointaggr over the whole fibre
    of per-class tables with these row sums, and that fibre grows
    combinatorially. The permanent evaluates the sum in closed form
    (pfqn_jointmarg).

    Args:
        sn: Network structure
        options: Solver options
        nvec: (M,) per-station total job counts
        engine: Permanent engine, 'exact' by default
        lG: Log normalizing constant already known, the same one that
            SolverNC.getProbAggr caches as result.Prob.logNormConstAggr

    Returns:
        SolverNCJointMargReturn with the probability, its logarithm, the
        constant and the runtime
    """
    from ...pfqn import pfqn_ca, pfqn_jointmarg

    start_time = time.time()

    if engine is None or engine == '':
        engine = 'exact'

    M = sn.nstations
    nvec = np.asarray(nvec, dtype=float).ravel()
    if nvec.size != M:
        raise ValueError("solver_nc_jointmarg: the occupancy vector has %d entries but the model "
                         "has %d stations." % (nvec.size, M))

    chain_result = sn_get_demands_chain(sn)
    Lchain = np.atleast_2d(np.asarray(chain_result.Lchain, dtype=float)).copy()
    Lchain[~np.isfinite(Lchain)] = 0.0
    Nchain = np.asarray(chain_result.Nchain, dtype=float).ravel()

    _jointmarg_supports(sn, Nchain)

    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    infset = np.where(np.isinf(nservers))[0]

    if lG is not None and not np.isfinite(lG):
        lG = None

    Pr, lPr = pfqn_jointmarg(nvec, Lchain, Nchain, infset, lG, engine)

    if lG is None:
        # Recover the constant the API computed, so the caller can cache it and
        # the sweep over a whole lattice pays for it once.
        isinfrow = np.zeros(M, dtype=bool)
        isinfrow[infset] = True
        Z = np.sum(Lchain[isinfrow, :], axis=0) if np.any(isinfrow) else np.zeros(Lchain.shape[1])
        _, lG = pfqn_ca(Lchain[~isinfrow, :], Nchain.reshape(1, -1), Z.reshape(1, -1))

    runtime = time.time() - start_time
    return SolverNCJointMargReturn(Pr=float(Pr), lPr=float(lPr), lG=float(lG), runtime=runtime)


def _jointmarg_supports(sn: NetworkStruct, Nchain: np.ndarray) -> None:
    """
    The permanent identity supplies one n_i! per queueing station and none per
    infinite server. A multiserver or load-dependent station has neither, so it
    is refused by name rather than approximated.
    """
    njobs = np.asarray(sn.njobs, dtype=float).ravel() if sn.njobs is not None else np.zeros(0)
    if np.any(~np.isfinite(Nchain)) or np.any(np.isinf(njobs)):
        raise ValueError("getProbSysMarg requires a closed model: the joint law of the total "
                         "queue lengths is not defined when a class has an infinite population.")
    lld = getattr(sn, 'lldscaling', None)
    if lld is not None and np.size(lld) > 0:
        raise ValueError("getProbSysMarg does not support load-dependent stations (sn.lldscaling "
                         "is set): the permanent identity supplies exactly one n_i! per queueing "
                         "station.")
    cds = getattr(sn, 'cdscaling', None)
    if cds is not None and len(cds) > 0:
        raise ValueError("getProbSysMarg does not support class-dependent scaling "
                         "(sn.cdscaling is set).")
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    ms = np.where(np.isfinite(nservers) & (nservers > 1))[0]
    if ms.size > 0:
        raise ValueError("getProbSysMarg does not support the multiserver station %d (%d servers): "
                         "the permanent identity supplies exactly one n_i! per queueing station."
                         % (ms[0] + 1, int(nservers[ms[0]])))


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
