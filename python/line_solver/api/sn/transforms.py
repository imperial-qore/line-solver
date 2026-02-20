"""
SN transform functions.

Native Python implementations of network structure transformation and
parameter extraction functions.

Port from:
    - /matlab/src/api/sn/sn_get_*.m
    - /matlab/src/api/sn/sn_set_*.m
    - /matlab/src/api/sn/sn_refresh_*.m
"""

import numpy as np
from typing import Optional, Dict, Any, Tuple, List, NamedTuple
from dataclasses import dataclass

from .network_struct import NetworkStruct, NodeType, SchedStrategy


def get_chain_for_class(chains: np.ndarray, class_idx: int) -> int:
    """
    Get the chain ID that a class belongs to.

    Handles both 1D and 2D chain formats:
    - 1D: chains[class_idx] = chain_id
    - 2D: chains[chain_id, class_idx] > 0 if class in chain

    Args:
        chains: Chain membership array (1D or 2D)
        class_idx: Index of the class

    Returns:
        Chain ID, or -1 if not found
    """
    if chains is None or chains.size == 0:
        return -1

    chains_arr = np.asarray(chains)
    if chains_arr.ndim == 1:
        # 1D format: chains[k] = chain_id for class k
        if class_idx < len(chains_arr):
            return int(chains_arr[class_idx])
        return -1
    elif chains_arr.ndim == 2:
        # 2D format: chains[c, k] > 0 means class k is in chain c
        if class_idx < chains_arr.shape[1]:
            chain_ids = np.where(chains_arr[:, class_idx] > 0)[0]
            if len(chain_ids) > 0:
                return int(chain_ids[0])
        return -1
    return -1


class ProductFormParams(NamedTuple):
    """Result of sn_get_product_form_params calculation."""
    lam: np.ndarray       # Arrival rates for open classes (1, R)
    D: np.ndarray         # Service demands at queueing stations (Mq, R)
    N: np.ndarray         # Population vector (1, R)
    Z: np.ndarray         # Think times at delay stations (1, R)
    mu: np.ndarray        # Load-dependent scaling factors (Mq, Nmax)
    S: np.ndarray         # Number of servers at queueing stations (Mq,)
    V: np.ndarray         # Visit ratios (M, R)


def sn_get_product_form_params(sn: NetworkStruct) -> ProductFormParams:
    """
    Extract standard product-form parameters from the network structure.

    This function extracts class-level parameters from a network structure for
    use in product-form queueing network analysis.

    Args:
        sn: NetworkStruct object

    Returns:
        ProductFormParams containing:
            - lam: Arrival rates for open classes
            - D: Service demands at queueing stations
            - N: Population vector
            - Z: Think times (service demands at delay stations)
            - mu: Load-dependent service capacity scaling factors
            - S: Number of servers at queueing stations
            - V: Visit ratios

    References:
        MATLAB: matlab/src/api/sn/sn_get_product_form_params.m
    """
    R = sn.nclasses
    N = sn.njobs.flatten() if sn.njobs is not None else np.zeros(R)

    # Find node types
    queue_indices = []
    delay_indices = []
    source_index = None

    for i, nt in enumerate(sn.nodetype):
        if nt == NodeType.QUEUE:
            queue_indices.append(i)
        elif nt == NodeType.DELAY:
            delay_indices.append(i)
        elif nt == NodeType.SOURCE:
            source_index = i

    Mq = len(queue_indices)
    Mz = len(delay_indices)

    # Initialize outputs
    lam = np.zeros(R)
    S = np.ones(Mq)

    # Get arrival rates for open classes
    if source_index is not None and sn.rates is not None:
        source_station = int(sn.nodeToStation[source_index]) if source_index < len(sn.nodeToStation) else -1
        if source_station >= 0 and source_station < sn.rates.shape[0]:
            for r in range(R):
                if np.isinf(N[r]) and r < sn.rates.shape[1]:
                    lam[r] = sn.rates[source_station, r]

    # Get number of servers
    for i, qi in enumerate(queue_indices):
        station_idx = int(sn.nodeToStation[qi]) if qi < len(sn.nodeToStation) else -1
        if station_idx >= 0 and sn.nservers is not None and station_idx < len(sn.nservers):
            S[i] = sn.nservers.flatten()[station_idx]

    # Compute service demands D
    Nct = np.sum(N[np.isfinite(N)])
    max_servers = max(int(np.max(S[np.isfinite(S)])) if len(S) > 0 and np.any(np.isfinite(S)) else 1, 1)
    D = np.zeros((max(Mq, 1), R))
    mu = np.ones((max(Mq, 1), int(Nct) + max_servers))

    if sn.rates is not None and sn.visits:
        for ist in range(Mq):
            qi = queue_indices[ist]
            station_idx = int(sn.nodeToStation[qi]) if qi < len(sn.nodeToStation) else -1
            stateful_idx = int(sn.nodeToStateful[qi]) if qi < len(sn.nodeToStateful) else -1

            for r in range(R):
                # Find chain containing class r
                chain_id = get_chain_for_class(sn.chains, r)

                if chain_id >= 0 and chain_id in sn.visits:
                    visits = sn.visits[chain_id]
                    if (station_idx >= 0 and station_idx < sn.rates.shape[0] and
                            r < sn.rates.shape[1] and sn.rates[station_idx, r] > 0):
                        if stateful_idx >= 0 and stateful_idx < visits.shape[0] and r < visits.shape[1]:
                            visit_ratio = visits[stateful_idx, r]
                            rate = sn.rates[station_idx, r]
                            # Normalize by reference station visit ratio (like MATLAB)
                            ref_visit_ratio = 1.0
                            if sn.refclass is not None and sn.refstat is not None and sn.stationToStateful is not None:
                                refclass_c = int(sn.refclass.flatten()[chain_id]) if chain_id < len(sn.refclass.flatten()) else -1
                                if refclass_c > 0:
                                    refstat_r = int(sn.refstat.flatten()[r]) if r < len(sn.refstat.flatten()) else -1
                                    if refstat_r >= 0 and refstat_r < len(sn.stationToStateful):
                                        refstat_stateful = int(sn.stationToStateful[refstat_r])
                                        if refstat_stateful >= 0 and refstat_stateful < visits.shape[0] and refclass_c < visits.shape[1]:
                                            ref_visit_ratio = visits[refstat_stateful, refclass_c]
                            D[ist, r] = (visit_ratio / rate / ref_visit_ratio) if rate > 0 and ref_visit_ratio > 0 else 0

            # Set mu scaling for multi-server
            if station_idx >= 0 and sn.nservers is not None:
                nserv = sn.nservers.flatten()[station_idx] if station_idx < len(sn.nservers.flatten()) else 1
                for n in range(mu.shape[1]):
                    mu[ist, n] = min(n + 1, nserv)

    # Compute think times Z
    Z = np.zeros(R)
    if sn.rates is not None and sn.visits:
        for ist in range(Mz):
            di = delay_indices[ist]
            station_idx = int(sn.nodeToStation[di]) if di < len(sn.nodeToStation) else -1
            stateful_idx = int(sn.nodeToStateful[di]) if di < len(sn.nodeToStateful) else -1

            for r in range(R):
                chain_id = get_chain_for_class(sn.chains, r)

                if chain_id >= 0 and chain_id in sn.visits:
                    visits = sn.visits[chain_id]
                    if (station_idx >= 0 and station_idx < sn.rates.shape[0] and
                            r < sn.rates.shape[1] and sn.rates[station_idx, r] > 0):
                        if stateful_idx >= 0 and stateful_idx < visits.shape[0] and r < visits.shape[1]:
                            visit_ratio = visits[stateful_idx, r]
                            rate = sn.rates[station_idx, r]
                            # Normalize by reference station visit ratio (like MATLAB)
                            ref_visit_ratio = 1.0
                            if sn.refclass is not None and sn.refstat is not None and sn.stationToStateful is not None:
                                refclass_c = int(sn.refclass.flatten()[chain_id]) if chain_id < len(sn.refclass.flatten()) else -1
                                if refclass_c > 0:
                                    refstat_r = int(sn.refstat.flatten()[r]) if r < len(sn.refstat.flatten()) else -1
                                    if refstat_r >= 0 and refstat_r < len(sn.stationToStateful):
                                        refstat_stateful = int(sn.stationToStateful[refstat_r])
                                        if refstat_stateful >= 0 and refstat_stateful < visits.shape[0] and refclass_c < visits.shape[1]:
                                            ref_visit_ratio = visits[refstat_stateful, refclass_c]
                            Z[r] += (visit_ratio / rate / ref_visit_ratio) if rate > 0 and ref_visit_ratio > 0 else 0

    # Compute total visits
    V = np.zeros((sn.nstations, R))
    if sn.visits:
        for chain_id, visits in sn.visits.items():
            if isinstance(visits, np.ndarray):
                for sf in range(min(visits.shape[0], sn.nstateful)):
                    station_idx = int(sn.statefulToStation[sf]) if sf < len(sn.statefulToStation) else -1
                    if station_idx >= 0 and station_idx < V.shape[0]:
                        for r in range(min(visits.shape[1], R)):
                            V[station_idx, r] += visits[sf, r]

    # Clean up NaN values
    D = np.nan_to_num(D, nan=0.0)
    Z = np.nan_to_num(Z, nan=0.0)

    return ProductFormParams(lam, D, N, Z, mu, S, V)


def sn_get_residt_from_respt(
    sn: NetworkStruct,
    RN: np.ndarray,
    WH: Optional[Dict] = None
) -> np.ndarray:
    """
    Compute residence times from response times.

    This function converts response times to residence times by accounting
    for visit ratios at each station.

    Args:
        sn: NetworkStruct object
        RN: Average response times (M, K)
        WH: Residence time handles (optional)

    Returns:
        WN: Average residence times (M, K)

    References:
        MATLAB: matlab/src/api/sn/sn_get_residt_from_respt.m
    """
    M = sn.nstations
    K = sn.nclasses
    WN = np.zeros((M, K))

    # Compute total visits by summing across all chains
    # sn.visits[chain_id] is indexed by stateful nodes (nstateful x nclasses)
    # We need to convert to station indices using statefulToStation
    V = np.zeros((M, K))
    if sn.visits:
        for chain_id, visits in sn.visits.items():
            if isinstance(visits, np.ndarray):
                nstateful = visits.shape[0]
                for sf in range(nstateful):
                    # Get station index for this stateful node
                    station_idx = -1
                    if sn.statefulToStation is not None and sf < len(sn.statefulToStation):
                        station_idx = int(sn.statefulToStation[sf])
                    if station_idx >= 0 and station_idx < M:
                        for r in range(min(visits.shape[1], K)):
                            V[station_idx, r] += visits[sf, r]

    for ist in range(M):
        for k in range(K):
            if WH is not None and (ist, k) in WH and WH[(ist, k)].get('disabled', False):
                WN[ist, k] = np.nan
            elif RN is not None and ist < RN.shape[0] and k < RN.shape[1] and RN[ist, k] > 0:
                if RN[ist, k] < 1e-14:
                    WN[ist, k] = RN[ist, k]
                else:
                    # Find chain containing class k
                    chain_id = get_chain_for_class(sn.chains, k)

                    # Get reference station for class k
                    refstat_k = int(sn.refstat.flatten()[k]) if sn.refstat is not None and k < len(sn.refstat.flatten()) else 0

                    # MATLAB formula: WN(ist,k) = RN(ist,k) * V(ist,k) / sum(V(refstat(k), refclass))
                    # MATLAB uses sn.refclass(c) if defined, otherwise sn.inchain{c}
                    # refclass is the reference class for the chain (usually the class at refstat)

                    # Get refclass for this chain
                    refclass = -1
                    if chain_id >= 0 and hasattr(sn, 'refclass') and sn.refclass is not None:
                        refclass_arr = np.asarray(sn.refclass).flatten()
                        if chain_id < len(refclass_arr):
                            refclass = int(refclass_arr[chain_id])

                    # Determine which classes to sum visits for
                    # Note: In Python (0-indexed), refclass = -1 means "not set"
                    # (equivalent to MATLAB's refclass = 0 in 1-indexed arrays)
                    # so we use >= 0 to check if a valid refclass is set
                    # For transient class models, refclass may have 0 visits (transient)
                    # so we need to check if refclass has non-zero visits first
                    use_refclass = False
                    if refclass >= 0 and refclass < V.shape[1] and refstat_k < V.shape[0]:
                        # Check if refclass has non-zero visits at reference station
                        if V[refstat_k, refclass] > 1e-10:
                            use_refclass = True

                    if use_refclass:
                        # Use just the reference class (matches MATLAB when refclass > 0 and has visits)
                        refclass_list = [refclass]
                    elif chain_id is not None and sn.inchain is not None and chain_id in sn.inchain:
                        # Fallback to all classes in chain
                        refclass_list = list(sn.inchain[chain_id])
                    else:
                        refclass_list = [k]  # fallback to single class

                    # Sum visits at reference station for refclass(es)
                    ref_visits_sum = 0.0
                    if refstat_k < V.shape[0]:
                        for rc in refclass_list:
                            if rc < V.shape[1]:
                                ref_visits_sum += V[refstat_k, rc]

                    if ref_visits_sum > 0:
                        WN[ist, k] = RN[ist, k] * V[ist, k] / ref_visits_sum

    # Clean up
    WN = np.nan_to_num(WN, nan=0.0)
    WN[WN < 1e-12] = 0.0

    return WN


def sn_get_state_aggr(sn: NetworkStruct) -> Dict[int, np.ndarray]:
    """
    Get aggregated state representation.

    Args:
        sn: NetworkStruct object

    Returns:
        Dictionary mapping stateful node index to aggregated state

    References:
        MATLAB: matlab/src/api/sn/sn_get_state_aggr.m
    """
    state_aggr = {}
    if sn.state is None:
        return state_aggr

    for node_id, state in sn.state.items():
        if state is not None:
            # Aggregate state by summing across phases
            if isinstance(state, np.ndarray) and state.ndim > 1:
                state_aggr[node_id] = np.sum(state, axis=0)
            else:
                state_aggr[node_id] = state

    return state_aggr


# ============================================================================
# Set functions - modify NetworkStruct in place
# ============================================================================

def sn_set_arrival(
    sn: NetworkStruct,
    station_idx: int,
    class_idx: int,
    rate: float
) -> None:
    """
    Set arrival rate for a class at a station.

    Args:
        sn: NetworkStruct object (modified in place)
        station_idx: Station index (0-based)
        class_idx: Class index (0-based)
        rate: Arrival rate

    References:
        MATLAB: matlab/src/api/sn/sn_set_arrival.m
    """
    if sn.rates is None:
        sn.rates = np.zeros((sn.nstations, sn.nclasses))

    if station_idx < sn.rates.shape[0] and class_idx < sn.rates.shape[1]:
        sn.rates[station_idx, class_idx] = rate


def sn_set_service(
    sn: NetworkStruct,
    station_idx: int,
    class_idx: int,
    rate: float,
    scv: float = 1.0
) -> None:
    """
    Set service rate for a class at a station.

    Args:
        sn: NetworkStruct object (modified in place)
        station_idx: Station index (0-based)
        class_idx: Class index (0-based)
        rate: Service rate
        scv: Squared coefficient of variation (default 1.0 for exponential)

    References:
        MATLAB: matlab/src/api/sn/sn_set_service.m
    """
    if sn.rates is None:
        sn.rates = np.zeros((sn.nstations, sn.nclasses))
    if sn.scv is None:
        sn.scv = np.ones((sn.nstations, sn.nclasses))

    if station_idx < sn.rates.shape[0] and class_idx < sn.rates.shape[1]:
        sn.rates[station_idx, class_idx] = rate
        sn.scv[station_idx, class_idx] = scv


def sn_set_servers(
    sn: NetworkStruct,
    station_idx: int,
    nservers: int
) -> None:
    """
    Set number of servers at a station.

    Args:
        sn: NetworkStruct object (modified in place)
        station_idx: Station index (0-based)
        nservers: Number of servers

    References:
        MATLAB: matlab/src/api/sn/sn_set_servers.m
    """
    if sn.nservers is None:
        sn.nservers = np.ones(sn.nstations)

    nservers_flat = sn.nservers.flatten()
    if station_idx < len(nservers_flat):
        nservers_flat[station_idx] = nservers
        sn.nservers = nservers_flat.reshape(sn.nservers.shape)


def sn_set_population(
    sn: NetworkStruct,
    class_idx: int,
    njobs: float
) -> None:
    """
    Set population for a class.

    Args:
        sn: NetworkStruct object (modified in place)
        class_idx: Class index (0-based)
        njobs: Number of jobs (inf for open class)

    References:
        MATLAB: matlab/src/api/sn/sn_set_population.m
    """
    if sn.njobs is None:
        sn.njobs = np.zeros(sn.nclasses)

    njobs_flat = sn.njobs.flatten()
    if class_idx < len(njobs_flat):
        njobs_flat[class_idx] = njobs
        sn.njobs = njobs_flat.reshape(sn.njobs.shape)

        # Update nclosedjobs as sum of finite populations
        finite_jobs = njobs_flat[np.isfinite(njobs_flat)]
        sn.nclosedjobs = int(np.sum(finite_jobs))


def sn_set_priority(
    sn: NetworkStruct,
    class_idx: int,
    priority: int
) -> None:
    """
    Set priority for a class.

    Args:
        sn: NetworkStruct object (modified in place)
        class_idx: Class index (0-based)
        priority: Priority level (higher = more priority)

    References:
        MATLAB: matlab/src/api/sn/sn_set_priority.m
    """
    if sn.classprio is None:
        sn.classprio = np.zeros(sn.nclasses)

    prio_flat = sn.classprio.flatten()
    if class_idx < len(prio_flat):
        prio_flat[class_idx] = priority
        sn.classprio = prio_flat.reshape(sn.classprio.shape)


def sn_set_routing(
    sn: NetworkStruct,
    source_node: int,
    dest_node: int,
    source_class: int,
    dest_class: int,
    prob: float
) -> None:
    """
    Set routing probability between nodes and classes.

    Args:
        sn: NetworkStruct object (modified in place)
        source_node: Source node index (0-based)
        dest_node: Destination node index (0-based)
        source_class: Source class index (0-based)
        dest_class: Destination class index (0-based)
        prob: Routing probability

    References:
        MATLAB: matlab/src/api/sn/sn_set_routing.m
    """
    N = sn.nnodes
    K = sn.nclasses

    if sn.rt is None:
        sn.rt = np.zeros((N * K, N * K))

    source_idx = source_node * K + source_class
    dest_idx = dest_node * K + dest_class

    if source_idx < sn.rt.shape[0] and dest_idx < sn.rt.shape[1]:
        sn.rt[source_idx, dest_idx] = prob


# ============================================================================
# Refresh functions - recompute derived quantities
# ============================================================================

def sn_refresh_visits(sn: NetworkStruct) -> None:
    """
    Refresh visit ratios from routing matrix.

    This function solves traffic equations to compute visit ratios at each
    station from the routing probability matrix.

    Args:
        sn: NetworkStruct object (modified in place)

    References:
        MATLAB: matlab/src/api/sn/sn_refresh_visits.m
    """
    if sn.rt is None:
        return

    from ..mc.dtmc import dtmc_solve_reducible, dtmc_solve

    FINE_TOL = 1e-10
    M = sn.nstateful  # rt is stateful-indexed (matches MATLAB: M = sn.nstateful)
    K = sn.nclasses
    N = sn.nnodes

    # Use rt_visits (which includes Sink->Source routing for open classes) if available.
    # This matches MATLAB where rt = dtmc_stochcomp(rtnodes) is computed AFTER adding
    # Sink->Source routing (getRoutingMatrix.m lines 325-335).
    rt_for_visits = getattr(sn, 'rt_visits', None)
    if rt_for_visits is None:
        rt_for_visits = sn.rt

    # Initialize visits and nodevisits dictionaries
    sn.visits = {}
    sn.nodevisits = {}

    # Force all classes in a chain to have the same reference station
    # (matches MATLAB sn_refresh_visits.m lines 48-53)
    refstat = sn.refstat.flatten() if sn.refstat is not None else np.zeros(K, dtype=int)
    for c in range(sn.nchains):
        if c not in sn.inchain:
            continue
        classes_in_chain = np.array([int(k) for k in sn.inchain[c]])
        # Check if all classes have the same refstat
        if len(classes_in_chain) > 0:
            first_refstat = int(refstat[classes_in_chain[0]]) if classes_in_chain[0] < len(refstat) else 0
            for k in classes_in_chain:
                if k < len(refstat) and int(refstat[k]) != first_refstat:
                    refstat[k] = first_refstat
    # Update sn.refstat
    sn.refstat = refstat.reshape(sn.refstat.shape) if sn.refstat is not None else refstat

    # Process each chain
    for c in range(sn.nchains):
        if c not in sn.inchain:
            continue

        classes_in_chain = np.array([int(k) for k in sn.inchain[c]])
        nIC = len(classes_in_chain)

        # For open chains with zero total arrival rate (e.g., auxiliary classes
        # with Disabled arrival in MMT fork-join models), set visits to 0.
        # These classes have no jobs entering from Source and contribute nothing.
        # In MATLAB, this case produces NaN routing (0/0) which propagates through
        # the DTMC solver, effectively giving zero visits.
        chain_is_open = any(np.isinf(sn.njobs[k]) for k in classes_in_chain if k < len(sn.njobs))
        if chain_is_open and sn.rates is not None:
            # Find Source station index
            source_station_idx = None
            if hasattr(sn, 'nodetype') and sn.nodetype is not None:
                for _node_idx in range(len(sn.nodetype)):
                    if sn.nodetype[_node_idx] == NodeType.SOURCE:
                        if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
                            source_station_idx = int(sn.nodeToStation[_node_idx])
                        break
            if source_station_idx is not None and source_station_idx >= 0:
                chain_arv_rates = sn.rates[source_station_idx, classes_in_chain]
                chain_arv_rates = np.where(np.isnan(chain_arv_rates), 0, chain_arv_rates)
                if np.sum(chain_arv_rates) < FINE_TOL:
                    # All arrival rates are zero - set visits to 0
                    sn.visits[c] = np.zeros((M, K))
                    sn.nodevisits[c] = np.zeros((N, K))
                    continue

        # ========================================================================
        # STATION VISITS
        # ========================================================================

        # Extract chain-specific routing matrix
        # Pchain[i,j] = P[(ist-1)*nIC+ik, (ist-1)*nIC+ik'] for ist, ik, ist', ik'
        cols = np.zeros(M * nIC, dtype=int)
        for ist in range(M):
            for ik_idx, ik in enumerate(classes_in_chain):
                cols[(ist) * nIC + ik_idx] = (ist) * K + int(ik)

        if np.any(cols >= rt_for_visits.shape[1]):
            # Handle bounds checking
            cols = cols[cols < rt_for_visits.shape[1]]

        Pchain = rt_for_visits[np.ix_(cols, cols)] if len(cols) > 0 else np.eye(len(cols))
        visited = np.sum(Pchain, axis=1) > FINE_TOL

        # Normalize routing matrix for Fork-containing models
        # Fork nodes have row sums > 1 (sending to all branches with prob 1 each)
        # Record original row sums to correct visit ratios after DTMC solve.
        row_sums = np.ones(Pchain.shape[0])
        if any(nt == NodeType.FORK for nt in sn.nodetype):
            for row in range(Pchain.shape[0]):
                rs = np.sum(Pchain[row, :])
                row_sums[row] = rs
                if rs > FINE_TOL:
                    Pchain[row, :] = Pchain[row, :] / rs

        # Solve traffic equations using DTMC
        if np.sum(visited) > 0:
            Pchain_visited = Pchain[np.ix_(np.where(visited)[0], np.where(visited)[0])]

            # Use dtmc_solve as primary, fallback to dtmc_solve_reducible for chains
            # with transient states (matches MATLAB sn_refresh_visits.m lines 100-106)
            # CRITICAL: Do not change this order - it affects visit ratio computations
            try:
                alpha_visited = dtmc_solve(Pchain_visited)
            except Exception:
                alpha_visited = np.full(Pchain_visited.shape[0], np.nan)
            # Fallback to dtmc_solve_reducible if dtmc_solve fails (e.g., reducible chain)
            if np.all(alpha_visited == 0) or np.any(np.isnan(alpha_visited)):
                try:
                    alpha_visited = dtmc_solve_reducible(Pchain_visited)
                except Exception:
                    alpha_visited = np.zeros(Pchain_visited.shape[0])
        else:
            alpha_visited = np.ones(np.sum(visited)) / np.sum(visited)

        # Expand back to full visited set
        alpha = np.zeros(M * nIC)
        alpha[visited] = alpha_visited

        # Apply Fork fanout correction: multiply fork-scope station visits by fanout
        # Block propagation through Join nodes to prevent the correction from
        # leaking back through cycles (which would cancel out during normalization).
        if any(nt == NodeType.FORK for nt in sn.nodetype):
            n = Pchain.shape[0]
            # Build adjacency but block outgoing edges from Join nodes
            adj = (Pchain > FINE_TOL).astype(bool).copy()
            for isf in range(M):
                nd = int(sn.statefulToNode[isf]) if hasattr(sn, 'statefulToNode') and sn.statefulToNode is not None else isf
                if nd < len(sn.nodetype) and sn.nodetype[nd] == NodeType.JOIN:
                    for k in range(nIC):
                        adj[isf * nIC + k, :] = False
            # Compute transitive closure with join-blocked adjacency
            reachable = adj.copy()
            for _ in range(int(np.ceil(np.log2(max(n, 2))))):
                reachable = reachable | (reachable @ reachable > 0)

            # For each fork row, multiply only fork-scope reachable states by fanout
            for row in range(n):
                fanout = row_sums[row]
                if fanout > 1 + FINE_TOL:
                    for col in range(n):
                        if reachable[row, col]:
                            alpha[col] = alpha[col] * fanout

        # Create visit matrix
        visits = np.zeros((M, K))
        for ist in range(M):
            for ik_idx, ik in enumerate(classes_in_chain):
                visits[ist, int(ik)] = alpha[ist * nIC + ik_idx]

        # Normalize by reference station visit
        # refstat is a station index; convert to stateful index via stationToNode (stateful = node in most cases)
        refstat_station = int(sn.refstat.flatten()[classes_in_chain[0]])
        if refstat_station < len(sn.stationToNode):
            refstat_idx = int(sn.stationToNode[refstat_station])  # Convert to node/stateful index
        else:
            refstat_idx = refstat_station
        if refstat_idx < M:
            normSum = np.sum(visits[refstat_idx, classes_in_chain])
            if normSum > FINE_TOL:
                visits = visits / normSum

        # Remove numerical noise
        visits = np.abs(visits)
        sn.visits[c] = visits

        # ========================================================================
        # NODE VISITS
        # ========================================================================

        if sn.rtnodes is not None:
            # Extract chain-specific node routing matrix
            nodes_cols = np.zeros(N * nIC, dtype=int)
            for ind in range(N):
                for ik_idx, ik in enumerate(classes_in_chain):
                    nodes_cols[ind * nIC + ik_idx] = ind * K + int(ik)

            if np.any(nodes_cols >= sn.rtnodes.shape[1]):
                nodes_cols = nodes_cols[nodes_cols < sn.rtnodes.shape[1]]

            nodes_Pchain = sn.rtnodes[np.ix_(nodes_cols, nodes_cols)] if len(nodes_cols) > 0 else np.eye(len(nodes_cols))

            # Handle NaN values in routing matrix (e.g., from Cache class switching)
            # For visits calculation, replace NaN with equal probabilities
            # (matches MATLAB sn_refresh_visits.m lines 137-153)
            for row in range(nodes_Pchain.shape[0]):
                nan_cols = np.isnan(nodes_Pchain[row, :])
                if np.any(nan_cols):
                    non_nan_sum = np.nansum(nodes_Pchain[row, ~nan_cols])
                    remaining_prob = max(0.0, 1.0 - non_nan_sum)
                    n_nan = np.sum(nan_cols)
                    if n_nan > 0 and remaining_prob > 0:
                        nodes_Pchain[row, nan_cols] = remaining_prob / n_nan
                    else:
                        nodes_Pchain[row, nan_cols] = 0.0

            nodes_visited = np.sum(nodes_Pchain, axis=1) > FINE_TOL

            # Normalize for Fork nodes
            # Record original row sums to correct visit ratios after DTMC solve.
            nodes_row_sums = np.ones(nodes_Pchain.shape[0])
            if any(nt == NodeType.FORK for nt in sn.nodetype):
                for row in range(nodes_Pchain.shape[0]):
                    rs = np.sum(nodes_Pchain[row, :])
                    nodes_row_sums[row] = rs
                    if rs > FINE_TOL:
                        nodes_Pchain[row, :] = nodes_Pchain[row, :] / rs

            # Solve traffic equations
            if np.sum(nodes_visited) > 0:
                nodes_Pchain_visited = nodes_Pchain[np.ix_(np.where(nodes_visited)[0], np.where(nodes_visited)[0])]

                # Use dtmc_solve as primary, fallback to dtmc_solve_reducible for chains
                # with transient states (matches MATLAB sn_refresh_visits.m lines 167-173)
                # CRITICAL: Do not change this order - it affects visit ratio computations
                try:
                    nodes_alpha_visited = dtmc_solve(nodes_Pchain_visited)
                except Exception:
                    nodes_alpha_visited = np.full(nodes_Pchain_visited.shape[0], np.nan)
                # Fallback to dtmc_solve_reducible if dtmc_solve fails (e.g., reducible chain)
                if np.all(nodes_alpha_visited == 0) or np.any(np.isnan(nodes_alpha_visited)):
                    try:
                        nodes_alpha_visited = dtmc_solve_reducible(nodes_Pchain_visited)
                    except Exception:
                        nodes_alpha_visited = np.zeros(nodes_Pchain_visited.shape[0])
            else:
                nodes_alpha_visited = np.ones(np.sum(nodes_visited)) / np.sum(nodes_visited)

            # Expand back to full visited set
            nodes_alpha = np.zeros(N * nIC)
            nodes_alpha[nodes_visited] = nodes_alpha_visited

            # Apply Fork fanout correction for node visits
            # Block propagation through Join nodes to prevent the correction from
            # leaking back through cycles.
            if any(nt == NodeType.FORK for nt in sn.nodetype):
                n_nodes = nodes_Pchain.shape[0]
                # Build adjacency but block outgoing edges from Join nodes
                nodes_adj = (nodes_Pchain > FINE_TOL).astype(bool).copy()
                for nd in range(N):
                    if nd < len(sn.nodetype) and sn.nodetype[nd] == NodeType.JOIN:
                        for k in range(nIC):
                            idx = nd * nIC + k
                            if idx < n_nodes:
                                nodes_adj[idx, :] = False
                # Compute transitive closure with join-blocked adjacency
                nodes_reachable = nodes_adj.copy()
                for _ in range(int(np.ceil(np.log2(max(n_nodes, 2))))):
                    nodes_reachable = nodes_reachable | (nodes_reachable @ nodes_reachable > 0)

                for row in range(n_nodes):
                    fanout = nodes_row_sums[row]
                    if fanout > 1 + FINE_TOL:
                        for col in range(n_nodes):
                            if nodes_reachable[row, col]:
                                nodes_alpha[col] = nodes_alpha[col] * fanout

            # Create nodevisits matrix
            nodevisits = np.zeros((N, K))
            for ind in range(N):
                for ik_idx, ik in enumerate(classes_in_chain):
                    nodevisits[ind, int(ik)] = nodes_alpha[ind * nIC + ik_idx]

            # Normalize by reference node visit
            if hasattr(sn, 'statefulToNode') and sn.statefulToNode is not None:
                refstat_idx = int(sn.refstat.flatten()[classes_in_chain[0]])
                refnode_idx = int(sn.statefulToNode[refstat_idx])
                nodeNormSum = np.sum(nodevisits[refnode_idx, classes_in_chain])
                if nodeNormSum > FINE_TOL:
                    nodevisits = nodevisits / nodeNormSum

            # Clean up numerical noise
            nodevisits[nodevisits < 0] = 0
            nodevisits = np.nan_to_num(nodevisits, nan=0.0)
            sn.nodevisits[c] = nodevisits


# ============================================================================
# Fork/Join functions
# ============================================================================

def sn_set_fork_fanout(
    sn: NetworkStruct,
    fork_node_idx: int,
    fan_out: int
) -> NetworkStruct:
    """
    Set fork fanout (tasksPerLink) for a Fork node.

    Updates the fanOut field in nodeparam for a Fork node.

    Args:
        sn: NetworkStruct object
        fork_node_idx: Node index of the Fork node (0-based)
        fan_out: Number of tasks per output link (>= 1)

    Returns:
        Modified NetworkStruct

    Raises:
        ValueError: If the specified node is not a Fork node

    References:
        MATLAB: matlab/src/api/sn/sn_set_fork_fanout.m
    """
    # Verify it's a Fork node
    if sn.nodetype[fork_node_idx] != NodeType.FORK:
        raise ValueError(f'sn_set_fork_fanout: Node {fork_node_idx} is not a Fork node')

    # Initialize nodeparam if needed
    if sn.nodeparam is None:
        sn.nodeparam = [None] * sn.nnodes

    if sn.nodeparam[fork_node_idx] is None:
        sn.nodeparam[fork_node_idx] = {}

    # Update nodeparam
    sn.nodeparam[fork_node_idx]['fanOut'] = fan_out

    return sn


# ============================================================================
# Batch update functions
# ============================================================================

def sn_set_service_batch(
    sn: NetworkStruct,
    rates: np.ndarray,
    scvs: Optional[np.ndarray] = None,
    auto_refresh: bool = False
) -> NetworkStruct:
    """
    Set service rates for multiple station-class pairs.

    Batch update of service rates. NaN values are skipped (not updated).
    More efficient than calling sn_set_service multiple times.

    Args:
        sn: NetworkStruct object
        rates: Matrix of new rates (nstations x nclasses), NaN = skip
        scvs: Matrix of new SCVs (optional)
        auto_refresh: If True, refresh process fields (default False)

    Returns:
        Modified NetworkStruct

    References:
        MATLAB: matlab/src/api/sn/sn_set_service_batch.m
    """
    from .utils import sn_refresh_process_fields

    M = sn.nstations
    K = sn.nclasses

    rates = np.asarray(rates)

    # Track updated pairs for auto-refresh
    updated_pairs = []

    # Update rates
    for i in range(min(M, rates.shape[0])):
        for j in range(min(K, rates.shape[1])):
            if not np.isnan(rates[i, j]):
                if sn.rates is None:
                    sn.rates = np.zeros((M, K))
                sn.rates[i, j] = rates[i, j]
                updated_pairs.append((i, j))

    # Update SCVs if provided
    if scvs is not None:
        scvs = np.asarray(scvs)
        for i in range(min(M, scvs.shape[0])):
            for j in range(min(K, scvs.shape[1])):
                if not np.isnan(scvs[i, j]):
                    if sn.scv is None:
                        sn.scv = np.ones((M, K))
                    sn.scv[i, j] = scvs[i, j]

    # Auto-refresh if requested
    if auto_refresh:
        for ist, r in updated_pairs:
            sn = sn_refresh_process_fields(sn, ist, r)

    return sn


# ============================================================================
# Non-Markovian to PH conversion
# ============================================================================

def sn_nonmarkov_toph(
    sn: NetworkStruct,
    options: Optional[Dict[str, Any]] = None
) -> NetworkStruct:
    """
    Convert non-Markovian distributions to Phase-Type using approximation.

    This function scans all service and arrival processes in the network
    structure and converts non-Markovian distributions to Markovian Arrival
    Processes (MAPs) using the specified approximation method.

    Supported non-Markovian distributions:
    - GAMMA: Gamma distribution
    - WEIBULL: Weibull distribution
    - LOGNORMAL: Lognormal distribution
    - PARETO: Pareto distribution
    - UNIFORM: Uniform distribution
    - DET: Deterministic (converted to Erlang)

    Args:
        sn: NetworkStruct object (from getStruct())
        options: Solver options dict with fields:
            - config.nonmkv: Method for conversion ('none', 'bernstein')
            - config.nonmkvorder: Number of phases for approximation (default 20)
            - config.preserveDet: Keep deterministic distributions (for MAP/D/c)

    Returns:
        Modified NetworkStruct with converted processes

    References:
        MATLAB: matlab/src/api/sn/sn_nonmarkov_toph.m
    """
    from ...constants import ProcessType
    from ..mam import map_bernstein, map_scale, map_erlang, map_pie
    import warnings
    from scipy import stats

    if options is None:
        options = {}

    # Get non-Markovian conversion method from options
    config = options.get('config', {})
    nonmkv_method = config.get('nonmkv', 'bernstein')

    # If method is 'none', return without any conversion
    if nonmkv_method.lower() == 'none':
        return sn

    # Get number of phases from options (default 20)
    n_phases = config.get('nonmkvorder', 20)

    # Check if we should preserve deterministic distributions
    preserve_det = config.get('preserveDet', False)

    # Markovian ProcessType IDs (no conversion needed)
    markovian_types = {
        ProcessType.EXP, ProcessType.ERLANG, ProcessType.HYPEREXP,
        ProcessType.PH, ProcessType.APH, ProcessType.MAP,
        ProcessType.COXIAN, ProcessType.COX2, ProcessType.MMPP2,
        ProcessType.IMMEDIATE, ProcessType.DISABLED
    }

    M = sn.nstations
    K = sn.nclasses

    for ist in range(M):
        for r in range(K):
            # Get process type
            if sn.procid is None or ist >= sn.procid.shape[0] or r >= sn.procid.shape[1]:
                continue

            proc_type_val = sn.procid[ist, r]

            # Skip if procType is NaN
            if proc_type_val is None or (isinstance(proc_type_val, float) and np.isnan(proc_type_val)):
                continue

            # Convert to ProcessType enum if needed
            if isinstance(proc_type_val, ProcessType):
                proc_type = proc_type_val
            elif isinstance(proc_type_val, (int, float)):
                proc_type_list = list(ProcessType)
                idx = int(proc_type_val)
                if 0 <= idx < len(proc_type_list):
                    proc_type = proc_type_list[idx]
                else:
                    continue
            else:
                continue

            # Skip if already Markovian, disabled, or immediate
            if proc_type in markovian_types:
                continue

            # Get target mean from rates
            if sn.rates is None or ist >= sn.rates.shape[0] or r >= sn.rates.shape[1]:
                continue
            rate = sn.rates[ist, r]
            if rate <= 0 or np.isnan(rate) or np.isinf(rate):
                continue
            target_mean = 1.0 / rate

            # Check if we should skip Det conversion for exact MAP/D/c analysis
            if proc_type == ProcessType.DET and preserve_det:
                continue

            # Issue warning
            warnings.warn(
                f'Distribution {proc_type.name} at station {ist} class {r} is '
                f'non-Markovian and will be converted to PH ({n_phases} phases).',
                UserWarning
            )

            # Get original process parameters
            orig_proc = None
            if sn.proc is not None and ist < len(sn.proc) and sn.proc[ist] is not None:
                if r < len(sn.proc[ist]):
                    orig_proc = sn.proc[ist][r]

            # Define PDF function based on distribution type
            pdf_func = None

            if proc_type == ProcessType.GAMMA:
                if orig_proc is not None and len(orig_proc) >= 2:
                    shape = orig_proc[0]
                    scale = orig_proc[1]
                    pdf_func = lambda x, s=shape, sc=scale: stats.gamma.pdf(x, a=s, scale=sc)

            elif proc_type == ProcessType.WEIBULL:
                if orig_proc is not None and len(orig_proc) >= 2:
                    shape_param = orig_proc[0]  # r
                    scale_param = orig_proc[1]  # alpha
                    pdf_func = lambda x, c=shape_param, sc=scale_param: stats.weibull_min.pdf(x, c=c, scale=sc)

            elif proc_type == ProcessType.LOGNORMAL:
                if orig_proc is not None and len(orig_proc) >= 2:
                    mu = orig_proc[0]
                    sigma = orig_proc[1]
                    pdf_func = lambda x, m=mu, s=sigma: stats.lognorm.pdf(x, s=s, scale=np.exp(m))

            elif proc_type == ProcessType.PARETO:
                if orig_proc is not None and len(orig_proc) >= 2:
                    shape_param = orig_proc[0]  # alpha
                    scale_param = orig_proc[1]  # k (minimum value)
                    pdf_func = lambda x, a=shape_param, sc=scale_param: stats.pareto.pdf(x, b=a, scale=sc)

            elif proc_type == ProcessType.UNIFORM:
                if orig_proc is not None and len(orig_proc) >= 2:
                    min_val = orig_proc[0]
                    max_val = orig_proc[1]
                    pdf_func = lambda x, lo=min_val, hi=max_val: stats.uniform.pdf(x, loc=lo, scale=hi-lo)

            elif proc_type == ProcessType.DET:
                # Deterministic: use Erlang approximation
                MAP = map_erlang(target_mean, n_phases)
                sn = _update_sn_for_map(sn, ist, r, MAP, n_phases)
                continue

            # Apply Bernstein approximation if PDF function is defined
            if pdf_func is not None:
                try:
                    MAP = map_bernstein(pdf_func, n_phases)
                    MAP = map_scale(MAP, target_mean)
                except Exception:
                    # Fallback to Erlang approximation
                    MAP = map_erlang(target_mean, n_phases)
            else:
                # Generic fallback: Erlang approximation
                MAP = map_erlang(target_mean, n_phases)

            # Update the network structure for the converted MAP
            actual_phases = MAP[0].shape[0] if isinstance(MAP, (list, tuple)) else n_phases
            sn = _update_sn_for_map(sn, ist, r, MAP, actual_phases)

    return sn


def _update_sn_for_map(
    sn: NetworkStruct,
    ist: int,
    r: int,
    MAP: Tuple[np.ndarray, np.ndarray],
    n_phases: int
) -> NetworkStruct:
    """
    Update all network structure fields for converted MAP.

    Updates proc, procid, phases, phasessz, phaseshift, mu, phi, pie, nvars, state.

    Args:
        sn: NetworkStruct object
        ist: Station index
        r: Class index
        MAP: MAP representation (D0, D1)
        n_phases: Number of phases

    Returns:
        Modified NetworkStruct

    References:
        MATLAB: matlab/src/api/sn/sn_nonmarkov_toph.m (updateSnForMAP helper)
    """
    from ...constants import ProcessType
    from ..mam import map_pie

    # Save old phasessz before updating (needed for state expansion)
    old_phases = 1
    if sn.phasessz is not None and ist < sn.phasessz.shape[0] and r < sn.phasessz.shape[1]:
        old_phases = int(sn.phasessz[ist, r])

    # Update process representation
    if sn.proc is None:
        sn.proc = [[None] * sn.nclasses for _ in range(sn.nstations)]
    while len(sn.proc) <= ist:
        sn.proc.append([None] * sn.nclasses)
    while len(sn.proc[ist]) <= r:
        sn.proc[ist].append(None)
    sn.proc[ist][r] = MAP

    # Update procid
    if sn.procid is None:
        sn.procid = np.zeros((sn.nstations, sn.nclasses), dtype=object)
    sn.procid[ist, r] = ProcessType.MAP

    # Update phases
    if sn.phases is None:
        sn.phases = np.ones((sn.nstations, sn.nclasses))
    sn.phases[ist, r] = n_phases

    # Update phasessz
    if sn.phasessz is None:
        sn.phasessz = np.ones((sn.nstations, sn.nclasses))
    sn.phasessz[ist, r] = max(n_phases, 1)

    # Recompute phaseshift for this station (cumulative sum across classes)
    if sn.phaseshift is None:
        sn.phaseshift = np.zeros((sn.nstations, sn.nclasses + 1))
    sn.phaseshift[ist, :] = np.concatenate([[0], np.cumsum(sn.phasessz[ist, :])])

    # Update mu (rates from -diag(D0))
    D0 = MAP[0]
    D1 = MAP[1]
    if sn.mu is None:
        sn.mu = [[None] * sn.nclasses for _ in range(sn.nstations)]
    while len(sn.mu) <= ist:
        sn.mu.append([None] * sn.nclasses)
    while len(sn.mu[ist]) <= r:
        sn.mu[ist].append(None)
    sn.mu[ist][r] = -np.diag(D0)

    # Update phi (completion probabilities: sum(D1,2) / -diag(D0))
    D0_diag = -np.diag(D0)
    D1_rowsum = np.sum(D1, axis=1)
    if sn.phi is None:
        sn.phi = [[None] * sn.nclasses for _ in range(sn.nstations)]
    while len(sn.phi) <= ist:
        sn.phi.append([None] * sn.nclasses)
    while len(sn.phi[ist]) <= r:
        sn.phi[ist].append(None)
    with np.errstate(divide='ignore', invalid='ignore'):
        phi_val = D1_rowsum / D0_diag
        phi_val = np.nan_to_num(phi_val, nan=1.0, posinf=1.0, neginf=0.0)
    sn.phi[ist][r] = phi_val

    # Update pie (initial phase distribution)
    if sn.pie is None:
        sn.pie = [[None] * sn.nclasses for _ in range(sn.nstations)]
    while len(sn.pie) <= ist:
        sn.pie.append([None] * sn.nclasses)
    while len(sn.pie[ist]) <= r:
        sn.pie[ist].append(None)
    sn.pie[ist][r] = map_pie(MAP)

    # Update nvars for MAP local variable
    if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
        ind = int(sn.stationToNode[ist]) if ist < len(sn.stationToNode) else ist
        if sn.nvars is None:
            sn.nvars = np.zeros((sn.nnodes, sn.nclasses))
        if ind < sn.nvars.shape[0] and r < sn.nvars.shape[1]:
            sn.nvars[ind, r] = sn.nvars[ind, r] + 1

    return sn
