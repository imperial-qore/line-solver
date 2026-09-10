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
import scipy.sparse as sp
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
                                if refclass_c >= 0:
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
                                if refclass_c >= 0:
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

    # see _kb/03-api-layer.md for rationale
    V = np.zeros((M, K))
    visits_obj = getattr(sn, 'visits', None)
    if visits_obj is not None and len(visits_obj) > 0:
        if isinstance(visits_obj, dict):
            visit_matrices = list(visits_obj.values())
        elif isinstance(visits_obj, np.ndarray):
            if visits_obj.dtype == object:
                visit_matrices = [entry for entry in visits_obj.flat if entry is not None]
            else:
                visit_matrices = [visits_obj]
        elif isinstance(visits_obj, (list, tuple)):
            visit_matrices = list(visits_obj)
        else:
            visit_matrices = [visits_obj]

        stateful_to_station = getattr(sn, 'statefulToStation', None)
        if stateful_to_station is None or len(stateful_to_station) == 0:
            stateful_to_station = np.arange(M)
        else:
            stateful_to_station = np.asarray(stateful_to_station).flatten()

        for visits in visit_matrices:
            if visits is None:
                continue

            visits = np.asarray(visits)
            if visits.ndim == 1:
                visits = visits.reshape((-1, 1))

            for sf in range(visits.shape[0]):
                station_idx = -1
                if sf < len(stateful_to_station):
                    mapped_idx = int(stateful_to_station[sf])
                    if 0 <= mapped_idx < M:
                        station_idx = mapped_idx
                if 0 <= station_idx < M:
                    cols = min(visits.shape[1], K)
                    V[station_idx, :cols] += visits[sf, :cols]

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

                    # see _kb/03-api-layer.md for rationale

                    # Get refclass for this chain
                    refclass = -1
                    if chain_id >= 0 and hasattr(sn, 'refclass') and sn.refclass is not None:
                        refclass_arr = np.asarray(sn.refclass).flatten()
                        if chain_id < len(refclass_arr):
                            refclass = int(refclass_arr[chain_id])

                    # see _kb/03-api-layer.md for rationale
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

    # Clean up (preserve Inf: saturated open stations report unbounded times)
    WN[np.isnan(WN)] = 0.0
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
        priority: Priority level (lower = more priority; 0 is highest)

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

def _sn_has_server_types(sn: NetworkStruct, ist: int) -> bool:
    """True when the station declares heterogeneous server types, whose rates live in nodeparam."""
    param = None
    if getattr(sn, 'nodeparam', None) is not None:
        ind = int(sn.stationToNode[ist])
        try:
            param = sn.nodeparam[ind]
        except (KeyError, IndexError, TypeError):
            param = None
    if param is None:
        return False
    n = param.get('nservertypes', 0) if isinstance(param, dict) else getattr(param, 'nservertypes', 0)
    return bool(n) and int(n) > 0


def _csr_fill_nan_rows(P: sp.csr_matrix) -> None:
    """
    Replace NaN routing entries with equal probabilities, in place, on CSR data.

    Sparse twin of the MATLAB per-row NaN fill in sn_refresh_visits.m: a Cache
    leaves its hit/miss routing NaN until the cache itself is solved, so the
    probability mass the row is missing is spread equally over its NaN entries.
    Structural zeros are never NaN, so only the stored entries need scanning.
    """
    if not np.isnan(P.data).any():
        return
    for row in range(P.shape[0]):
        s, e = P.indptr[row], P.indptr[row + 1]
        row_data = P.data[s:e]
        nan_mask = np.isnan(row_data)
        n_nan = int(nan_mask.sum())
        if n_nan == 0:
            continue
        non_nan_sum = row_data[~nan_mask].sum()
        remaining_prob = max(0.0, 1.0 - non_nan_sum)
        row_data[nan_mask] = remaining_prob / n_nan if remaining_prob > 0 else 0.0
    P.eliminate_zeros()


def _csr_normalize_rows(P: sp.csr_matrix, tol: float) -> np.ndarray:
    """
    Scale each row of a CSR matrix by its own sum, in place, and return the sums.

    Fork nodes route to every branch with probability 1, so their rows sum above
    1 and the chain is not stochastic; the original sums are returned so the
    caller can detect the fork rows and correct the visit ratios afterwards.
    """
    row_sums = np.asarray(P.sum(axis=1)).ravel()
    for row in range(P.shape[0]):
        rs = row_sums[row]
        if rs > tol:
            s, e = P.indptr[row], P.indptr[row + 1]
            P.data[s:e] /= rs
    return row_sums


def _csr_threshold_mask(P: sp.csr_matrix, tol: float) -> sp.csr_matrix:
    """
    Return the sparsity pattern of the entries of P that exceed tol.

    Sparse twin of the dense `P > tol`, which scipy rejects as inefficient for a
    nonzero scalar; structural zeros never exceed a positive tol, so only the
    stored entries need testing.
    """
    mask = P.copy()
    mask.data = (mask.data > tol).astype(np.int8)
    mask.eliminate_zeros()
    return mask


def sn_refresh_cacheqn_visits(sn: NetworkStruct) -> None:
    """
    Relabel every Cache node's self-switch with the split standing on the node,
    then refresh the visits derived from it.

    RESTORING THE HIT/MISS SPLIT IS NOT ENOUGH, BECAUSE THE VISITS ARE DERIVED
    FROM IT. `link()` lays down a uniform hit/miss split before any cache has
    been analyzed; a solver that writes `actualhitprob` back onto the node
    without this step leaves `rtnodes` -- and so `nodevisits` -- carrying that
    guess, and the node-level ResidT is then RespT times the wrong visit. Every
    MATLAB solver that writes a hit probability follows it with `refreshChains`
    for this reason, and the native analyzers relabel and call
    `sn_refresh_visits` inline (`solver_nc_cacheqn_analyzer`).

    Intended for the delegating bridges (`lang='java'`, `lang='cpp'`), which
    take the split from a foreign engine and must reproduce that refresh here.
    A cache whose node carries no split is left alone rather than zeroed: absent
    means the engine reported none, and the offered routing is then all there is.

    Args:
        sn: NetworkStruct object (modified in place)
    """
    if sn.rtnodes is None or sn.nodeparam is None:
        return

    from ..mc.dtmc import dtmc_stochcomp

    K = int(sn.nclasses)
    I = int(sn.nnodes)
    touched = False

    for ind in range(I):
        if sn.nodetype is None or ind >= len(sn.nodetype) or sn.nodetype[ind] != NodeType.CACHE:
            continue
        cp = (sn.nodeparam.get(ind) if isinstance(sn.nodeparam, dict)
              else (sn.nodeparam[ind] if ind < len(sn.nodeparam) else None))
        if cp is None:
            continue
        hitprob = getattr(cp, 'actualhitprob', None)
        missprob = getattr(cp, 'actualmissprob', None)
        if hitprob is None or missprob is None:
            continue
        # A retrieval system routes the read class into its retrieval classes as
        # well, so a two-way hit/miss rewrite of that row would DELETE them; its
        # own analyzer (da_cacheqn_retrieval) owns that routing. Test for a class
        # actually declared, NOT for the array's size: `retrieval_classes` is an
        # (items x classes) block of -1 sentinels on every ordinary cache, so a
        # size test skips every model and the refresh never runs at all.
        rc = getattr(cp, 'retrieval_classes', None)
        if rc is not None and np.any(np.asarray(rc) >= 0):
            continue
        hitprob = np.atleast_1d(hitprob).flatten()
        missprob = np.atleast_1d(missprob).flatten()
        # A delayed hit is a hit that waited: it leaves by the hit class, so the
        # branch probability is the sum. SolverMVA's retrieval analyzer splits
        # the same way when it writes the node.
        delayed = getattr(cp, 'actualdelayedhitprob', None)
        if delayed is not None:
            delayed = np.atleast_1d(delayed).flatten()
        hitclass = np.atleast_1d(getattr(cp, 'hitclass', [])).flatten().astype(int)
        missclass = np.atleast_1d(getattr(cp, 'missclass', [])).flatten().astype(int)

        for r in range(K):
            if r >= len(hitclass) or hitclass[r] < 0:
                continue
            if r >= len(hitprob) or not np.isfinite(hitprob[r]):
                continue
            ph = float(hitprob[r])
            if delayed is not None and r < len(delayed) and np.isfinite(delayed[r]):
                ph += float(delayed[r])
            pm = float(missprob[r]) if r < len(missprob) and np.isfinite(missprob[r]) else 1.0 - ph
            sn.rtnodes[ind * K + r, :] = 0
            for jnd in range(I):
                if sn.connmatrix is None or ind >= sn.connmatrix.shape[0] \
                        or jnd >= sn.connmatrix.shape[1] or not sn.connmatrix[ind, jnd]:
                    continue
                hc = hitclass[r]
                mc = missclass[r] if r < len(missclass) else -1
                if 0 <= hc < K:
                    sn.rtnodes[ind * K + r, jnd * K + hc] = ph
                if 0 <= mc < K:
                    sn.rtnodes[ind * K + r, jnd * K + mc] = pm
            touched = True

    if not touched:
        return

    stateful_nodes = [i for i in range(I)
                      if sn.isstateful is None or (i < len(sn.isstateful) and sn.isstateful[i])]
    idx = np.array([ind * K + k for ind in stateful_nodes for k in range(K)], dtype=int)
    new_rt = dtmc_stochcomp(sn.rtnodes, idx)
    sn.rt = new_rt
    if getattr(sn, 'rt_visits', None) is not None:
        sn.rt_visits = new_rt.copy()
    sn_refresh_visits(sn)


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

    # see _kb/03-api-layer.md for rationale
    rt_for_visits = getattr(sn, 'rt_visits', None)
    if rt_for_visits is None:
        rt_for_visits = sn.rt

    # the per-chain routing matrices are slices of these, so sparsify once here
    # and keep sparse storage all the way to the DTMC solve
    rt_sparse = sp.csr_matrix(rt_for_visits)
    rtnodes_sparse = sp.csr_matrix(sn.rtnodes) if sn.rtnodes is not None else None

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

        # see _kb/03-api-layer.md for rationale
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

        if len(cols) > 0:
            Pchain = sp.csr_matrix(rt_sparse[cols, :][:, cols])
        else:
            Pchain = sp.csr_matrix((0, 0))

        # Match MATLAB sn_refresh_visits: replace NaN routing entries before DTMC solve.
        _csr_fill_nan_rows(Pchain)

        # the routing matrix carries a JMT-oriented uniform fill on DISABLED (node,class) pairs;
        # a class with no service at a station cannot be there -- see _kb/06-solver-catalog.md
        served = np.ones(Pchain.shape[0])
        # ABSENT IS THE EMPTY ARRAY HERE, NOT None: every index map on
        # NetworkStruct defaults to np.array([]), so an `is not None` test is
        # always true and the -1 fallback below was unreachable. A struct built
        # without rates has no disabled (station,class) pair to drop.
        if sn.rates is not None and np.size(sn.rates) > 0:
            stf2st = sn.statefulToStation
            have_stf2st = stf2st is not None and np.size(stf2st) > 0
            for ist in range(M):
                sti = int(stf2st[ist]) if (have_stf2st and ist < np.size(stf2st)) else -1
                if sti < 0 or sti >= sn.nstations:
                    continue
                # a Place, and a station declaring server types, carry NaN station rates by construction
                ntype = sn.nodetype[int(sn.stationToNode[sti])]
                if ntype in (NodeType.PLACE, NodeType.TRANSITION):
                    continue
                if _sn_has_server_types(sn, sti):
                    continue
                for ik_idx, ik in enumerate(classes_in_chain):
                    if not np.isnan(sn.rates[sti, int(ik)]):
                        continue
                    idx = ist * nIC + ik_idx
                    if idx < served.shape[0]:
                        served[idx] = 0.0
        if served.size and not served.all():
            D = sp.diags(served)
            Pchain = sp.csr_matrix(D @ Pchain @ D)
            Pchain.eliminate_zeros()

        visited = np.asarray(Pchain.sum(axis=1)).ravel() > 0

        # see _kb/03-api-layer.md for rationale
        row_sums = np.ones(Pchain.shape[0])
        if any(nt == NodeType.FORK for nt in sn.nodetype):
            row_sums = _csr_normalize_rows(Pchain, FINE_TOL)

        # Solve traffic equations using DTMC
        if np.sum(visited) > 0:
            vidx = np.where(visited)[0]
            Pchain_visited = Pchain[vidx, :][:, vidx]

            # see _kb/03-api-layer.md for rationale
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

        # see _kb/03-api-layer.md for rationale
        if any(nt == NodeType.FORK for nt in sn.nodetype) and np.any(row_sums > 1 + FINE_TOL):
            for idx in range(len(alpha)):
                if alpha[idx] > FINE_TOL:
                    alpha[idx] = 1

        # Create visit matrix
        visits = np.zeros((M, K))
        for ist in range(M):
            for ik_idx, ik in enumerate(classes_in_chain):
                visits[ist, int(ik)] = alpha[ist * nIC + ik_idx]

        # Normalize by reference station visit using the stateful mapping, like MATLAB.
        refstat_station = int(sn.refstat.flatten()[classes_in_chain[0]])
        if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None and refstat_station < len(sn.stationToStateful):
            refstat_idx = int(sn.stationToStateful[refstat_station])
        else:
            refstat_idx = refstat_station
        if refstat_idx < M:
            normSum = np.sum(visits[refstat_idx, classes_in_chain])
            if normSum > FINE_TOL:
                visits = visits / normSum

        # Remove numerical noise
        visits = np.abs(visits)
        sn.visits[c] = visits

        # see _kb/03-api-layer.md for rationale

        if sn.rtnodes is not None:
            # Extract chain-specific node routing matrix
            nodes_cols = np.zeros(N * nIC, dtype=int)
            for ind in range(N):
                for ik_idx, ik in enumerate(classes_in_chain):
                    nodes_cols[ind * nIC + ik_idx] = ind * K + int(ik)

            if np.any(nodes_cols >= sn.rtnodes.shape[1]):
                nodes_cols = nodes_cols[nodes_cols < sn.rtnodes.shape[1]]

            if len(nodes_cols) > 0:
                nodes_Pchain = sp.csr_matrix(rtnodes_sparse[nodes_cols, :][:, nodes_cols])
            else:
                nodes_Pchain = sp.csr_matrix((0, 0))

            # see _kb/03-api-layer.md for rationale
            _csr_fill_nan_rows(nodes_Pchain)

            # THE SAME DISABLED-PAIR MASK THE STATION BLOCK APPLIES ABOVE, and it
            # matters more here: at station level a (station,class) the class
            # cannot be served at is a dead end, while the node kernel keeps the
            # class-switch nodes between the stations, so the disabled states
            # close into a whole spurious CYCLE. A materialised LQN replica is
            # exactly that -- replica 2's stations still carry replica 1's classes
            # in rtnodes -- and dtmc_solve_reducible then splits the mass between
            # the real chain and the phantom one, giving every node of replica 2
            # a visit in replica 1's classes. See _kb/03-api-layer.md.
            nodes_served = np.ones(nodes_Pchain.shape[0])
            if sn.rates is not None and np.size(sn.rates) > 0:
                nd2st = sn.nodeToStation
                have_nd2st = nd2st is not None and np.size(nd2st) > 0
                for ind in range(N):
                    sti = int(nd2st[ind]) if (have_nd2st and ind < np.size(nd2st)) else -1
                    if sti < 0 or sti >= sn.nstations:
                        continue
                    # a Place, and a station declaring server types, carry NaN station rates by construction
                    if sn.nodetype[ind] in (NodeType.PLACE, NodeType.TRANSITION):
                        continue
                    if _sn_has_server_types(sn, sti):
                        continue
                    for ik_idx, ik in enumerate(classes_in_chain):
                        if not np.isnan(sn.rates[sti, int(ik)]):
                            continue
                        idx = ind * nIC + ik_idx
                        if idx < nodes_served.shape[0]:
                            nodes_served[idx] = 0.0
            if nodes_served.size and not nodes_served.all():
                nodes_D = sp.diags(nodes_served)
                nodes_Pchain = sp.csr_matrix(nodes_D @ nodes_Pchain @ nodes_D)
                nodes_Pchain.eliminate_zeros()

            nodes_visited = np.asarray(nodes_Pchain.sum(axis=1)).ravel() > 0

            # Normalize for Fork nodes
            # Record original row sums to correct visit ratios after DTMC solve.
            nodes_row_sums = np.ones(nodes_Pchain.shape[0])
            if any(nt == NodeType.FORK for nt in sn.nodetype):
                nodes_row_sums = _csr_normalize_rows(nodes_Pchain, FINE_TOL)

            # Solve traffic equations
            if np.sum(nodes_visited) > 0:
                nvidx = np.where(nodes_visited)[0]
                nodes_Pchain_visited = nodes_Pchain[nvidx, :][:, nvidx]

                # see _kb/03-api-layer.md for rationale
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

            # SPN-based fork correction for node visits: stations/Fork get visit=1,
            # Join nodes get visit = number of direct predecessors.
            if any(nt == NodeType.FORK for nt in sn.nodetype) and np.any(nodes_row_sums > 1 + FINE_TOL):
                for idx in range(len(nodes_alpha)):
                    if nodes_alpha[idx] > FINE_TOL:
                        nd = idx // nIC
                        if nd < len(sn.nodetype) and sn.nodetype[nd] == NodeType.JOIN:
                            r = int(classes_in_chain[idx % nIC])
                            col = nd * K + r
                            n_sources = int(np.sum(sn.rtnodes[:, col] > FINE_TOL))
                            nodes_alpha[idx] = n_sources
                        else:
                            nodes_alpha[idx] = 1

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
    from ..mam import map_bernstein, map_scale, map_erlang, map_pie, map_mean
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
    # Family used when a concrete distribution has to be replaced by a Markovian
    # surrogate: 'cme' fits a concentrated matrix exponential plus an exponential
    # tail, 'ph' keeps the Erlang. At a budget of n_phases the ME reaches an SCV
    # of O(1/n_phases^2) where the Erlang stops at 1/n_phases, so 'cme' is the
    # default; SSA, Fluid and JMT pass 'ph' because they cannot consume an ME.
    phfit = str(config.get('phfit', 'cme')).lower()

    # Check if we should preserve deterministic distributions
    preserve_det = config.get('preserveDet', False)

    # Markovian ProcessType IDs (no conversion needed)
    markovian_types = {
        ProcessType.EXP, ProcessType.ERLANG, ProcessType.HYPEREXP,
        ProcessType.PH, ProcessType.APH, ProcessType.MAP, ProcessType.MMAP,
        # see _kb/03-api-layer.md for rationale
        ProcessType.ME, ProcessType.RAP,
        ProcessType.COXIAN, ProcessType.COX2, ProcessType.MMPP2,
        ProcessType.IMMEDIATE, ProcessType.DISABLED,
        # see _kb/03-api-layer.md for rationale
        ProcessType.NHPP,
        ProcessType.MAPT, ProcessType.PHT
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

            # Set only by the families with a closed-form tail, and read only by
            # the long-tail fit below.
            ccdf_func = None

            if proc_type == ProcessType.GAMMA:
                if orig_proc is not None and len(orig_proc) >= 2:
                    shape = orig_proc[0]
                    scale = orig_proc[1]
                    pdf_func = lambda x, s=shape, sc=scale: stats.gamma.pdf(x, a=s, scale=sc)
                    ccdf_func = lambda x, s=shape, sc=scale: float(stats.gamma.sf(x, a=s, scale=sc))

            elif proc_type == ProcessType.WEIBULL:
                if orig_proc is not None and len(orig_proc) >= 2:
                    shape_param = orig_proc[0]  # r
                    scale_param = orig_proc[1]  # alpha
                    pdf_func = lambda x, c=shape_param, sc=scale_param: stats.weibull_min.pdf(x, c=c, scale=sc)
                    ccdf_func = lambda x, c=shape_param, sc=scale_param: float(stats.weibull_min.sf(x, c=c, scale=sc))

            elif proc_type == ProcessType.LOGNORMAL:
                if orig_proc is not None and len(orig_proc) >= 2:
                    mu = orig_proc[0]
                    sigma = orig_proc[1]
                    pdf_func = lambda x, m=mu, s=sigma: stats.lognorm.pdf(x, s=s, scale=np.exp(m))
                    ccdf_func = lambda x, m=mu, s=sigma: float(stats.lognorm.sf(x, s=s, scale=np.exp(m)))

            elif proc_type == ProcessType.PARETO:
                if orig_proc is not None and len(orig_proc) >= 2:
                    shape_param = orig_proc[0]  # alpha
                    scale_param = orig_proc[1]  # k (minimum value)
                    pdf_func = lambda x, a=shape_param, sc=scale_param: stats.pareto.pdf(x, b=a, scale=sc)
                    ccdf_func = lambda x, a=shape_param, sc=scale_param: float(stats.pareto.sf(x, b=a, scale=sc))

            elif proc_type == ProcessType.UNIFORM:
                if orig_proc is not None and len(orig_proc) >= 2:
                    min_val = orig_proc[0]
                    max_val = orig_proc[1]
                    pdf_func = lambda x, lo=min_val, hi=max_val: stats.uniform.pdf(x, loc=lo, scale=hi-lo)

            elif proc_type == ProcessType.DET:
                # Deterministic: the most concentrated surrogate the phase budget
                # allows. Erlang-20 only reaches SCV 0.05; the CME plus exponential
                # reaches 5.7e-3 at the same 20 phases.
                MAP, actual = _fit_concentrated_surrogate(target_mean, 0.0, n_phases, phfit)
                sn = _update_sn_for_map(sn, ist, r, MAP, actual)
                continue

            # A concentrated ME matching the first two moments EXACTLY reproduces
            # the Pollaczek-Khinchine mean, which the Bernstein fit does not: on
            # M/Gamma/1 at rho 0.5 the shape fit lands 2.7e-2 away from the exact
            # mean queue length while the two-moment ME lands on it. The Bernstein
            # path is kept for phfit='ph', where it carries shape information that
            # a two-moment fit cannot, and for the solvers that need a phase-type.
            # The long-tail fit, when it was asked for and the law has a tail to
            # fit. It matches the ccdf at points spread over decades rather than
            # matching two moments, so it is the route for a Pareto, a Weibull
            # with shape below one or a Lognormal with a large sigma; a
            # light-tailed law has nothing for it to do and falls through to the
            # fits below.
            if phfit == 'hyperexp' and ccdf_func is not None:
                MAP_lt, n_lt = _fit_long_tail_surrogate(ccdf_func, target_mean)
                if MAP_lt is not None:
                    sn = _update_sn_for_map(sn, ist, r, MAP_lt, n_lt)
                    continue

            target_scv = float(sn.scv[ist, r]) if sn.scv is not None else 1.0
            if phfit == 'cme' and 0.0 <= target_scv < 1.0:
                MAP, actual_me = _fit_concentrated_surrogate(
                    target_mean, target_scv, n_phases, phfit)
                sn = _update_sn_for_map(sn, ist, r, MAP, actual_me)
                continue

            # Apply Bernstein approximation if PDF function is defined
            if pdf_func is not None:
                MAP = map_bernstein(pdf_func, n_phases)
                # Rescale to the target mean: map_scale multiplies the rates
                # by factor, dividing the mean by factor.
                cur_mean = map_mean(MAP[0], MAP[1])
                MAP = map_scale(MAP[0], MAP[1], target_mean)
            else:
                # Generic fallback: same concentrated surrogate as the Det branch,
                # targeting the SCV recorded in sn.
                target_scv = float(sn.scv[ist, r]) if sn.scv is not None else 0.0
                MAP, _ = _fit_concentrated_surrogate(target_mean, target_scv, n_phases, phfit)

            # Update the network structure for the converted MAP
            actual_phases = MAP[0].shape[0] if isinstance(MAP, (list, tuple)) else n_phases
            sn = _update_sn_for_map(sn, ist, r, MAP, actual_phases)

    # see _kb/03-api-layer.md for rationale
    if hasattr(sn, 'nodeparam') and sn.nodeparam:
        from ...lang.base import NodeType as _NodeType
        transition_v = int(_NodeType.TRANSITION.value) if hasattr(_NodeType.TRANSITION, 'value') else int(_NodeType.TRANSITION)
        for ind, nparam in sn.nodeparam.items():
            if not hasattr(sn, 'nodetype') or sn.nodetype is None or ind >= len(sn.nodetype):
                continue
            nt = sn.nodetype[ind]
            nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
            if nt_val != transition_v:
                continue
            nmodes = int(getattr(nparam, 'nmodes', 0) if not isinstance(nparam, dict) else nparam.get('nmodes', 0))
            if nmodes <= 0:
                continue

            distributions = getattr(nparam, 'distributions', None) if not isinstance(nparam, dict) else nparam.get('distributions', None)
            firingproc = getattr(nparam, 'firingproc', None) if not isinstance(nparam, dict) else nparam.get('firingproc', None)
            firingpie = getattr(nparam, 'firingpie', None) if not isinstance(nparam, dict) else nparam.get('firingpie', None)
            firingphases = getattr(nparam, 'firingphases', None) if not isinstance(nparam, dict) else nparam.get('firingphases', None)
            firingprocid = getattr(nparam, 'firingprocid', None) if not isinstance(nparam, dict) else nparam.get('firingprocid', None)
            if distributions is None or firingproc is None:
                continue

            for m in range(nmodes):
                dist = distributions[m] if m < len(distributions) else None
                if dist is None:
                    continue
                # Skip if already PH (Markovian path populated firingproc).
                already_ph = (firingproc[m] is not None) and (firingphases is not None and m < len(firingphases) and not np.isnan(float(firingphases[m])))
                if already_ph:
                    continue

                # Resolve target mean.
                try:
                    target_mean = float(dist.getMean())
                except Exception:
                    target_mean = None
                if target_mean is None or target_mean <= 0 or not np.isfinite(target_mean):
                    continue

                # Identify class name for warning + DET branch.
                class_name = type(dist).__name__

                if class_name == 'Det' and not preserve_det:
                    MAP = map_erlang(target_mean, n_phases)
                else:
                    # Generic continuous distribution: use evalPDF directly.
                    if hasattr(dist, 'evalPDF'):
                        pdf_func = lambda x, _d=dist: float(_d.evalPDF(x))
                        MAP = map_bernstein(pdf_func, n_phases)
                        cur_mean = map_mean(MAP[0], MAP[1])
                        MAP = map_scale(MAP[0], MAP[1], target_mean)
                    else:
                        MAP = map_erlang(target_mean, n_phases)

                D0 = np.atleast_2d(np.asarray(MAP[0], dtype=float))
                D1 = np.atleast_2d(np.asarray(MAP[1], dtype=float))
                actual_phases = int(D0.shape[0])
                try:
                    pie = np.atleast_1d(np.asarray(map_pie(MAP), dtype=float)).ravel()
                except Exception:
                    pie = np.zeros(actual_phases, dtype=float)
                    if actual_phases > 0:
                        pie[0] = 1.0
                if pie.size < actual_phases:
                    pie = np.pad(pie, (0, actual_phases - pie.size))

                warnings.warn(
                    f'Firing distribution {class_name} at Transition node {ind} mode {m} is '
                    f'non-Markovian and will be converted to PH ({actual_phases} phases).',
                    UserWarning
                )

                firingproc[m] = (D0, D1)
                firingpie[m] = pie
                firingphases[m] = float(actual_phases)
                if firingprocid is not None and m < len(firingprocid):
                    firingprocid[m] = int(ProcessType.MAP.value) if hasattr(ProcessType.MAP, 'value') else int(ProcessType.MAP)

    return sn



def _fit_long_tail_surrogate(ccdf_func, target_mean):
    """
    A hyperexponential fitted to the ccdf across decades (Feldmann and Whitt
    1998), returned as its MAP pair and rescaled to the mean the struct carries.

    ``(None, 0)`` when the recursion declines the law: the components have to
    dominate one another at their own time scales, which a light-tailed law does
    not provide, and answering with a fit that does not hold is worse than
    falling through to the two-moment surrogate.
    """
    from ..mam.hyperexp_longtail import hyperexp_fit_longtail
    from ..mam import map_scale
    try:
        fit = hyperexp_fit_longtail(ccdf_func)
    except Exception:
        return None, 0
    p = np.asarray(fit['p'], dtype=float).ravel()
    lam = np.asarray(fit['lambda'], dtype=float).ravel()
    if p.size == 0 or not np.all(np.isfinite(lam)) or np.any(lam <= 0):
        return None, 0
    n = p.size
    D0 = -np.diag(lam)
    D1 = np.diag(lam) @ np.tile(p, (n, 1))
    MAP = map_scale(D0, D1, target_mean)
    return MAP, n


def _fit_concentrated_surrogate(target_mean, target_scv, n_phases, phfit):
    """Build the Markovian surrogate of a concrete distribution.

    With ``phfit='cme'`` the surrogate is a concentrated matrix exponential
    convolved with an exponential (see fit_me_mean_scv): under a budget of
    ``n_phases`` it reaches an SCV of O(1/n_phases^2), where an Erlang of the same
    order stops at 1/n_phases. With ``phfit='ph'`` the Erlang is kept, which is
    what SSA, Fluid and JMT need since they cannot consume a matrix exponential.

    Returns the (D0, D1) pair and its actual number of phases.
    """
    from ..mam import map_erlang

    if phfit == 'cme':
        from ...distributions.markovian import fit_me_mean_scv
        # A Det has SCV 0, which no ME attains; the budget-limited branch of the
        # fitter then returns the most concentrated member that fits.
        scv = max(float(target_scv), 1e-12)
        if scv < 1.0:
            fitted = fit_me_mean_scv(target_mean, min(scv, 1.0 - 1e-12), max_phases=n_phases)
            proc = fitted.getProcess()
            return (proc[0], proc[1]), fitted.getNumberOfPhases()

    MAP = map_erlang(target_mean, n_phases)
    return MAP, n_phases


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
    # The conversion methods (map_bernstein, map_erlang, the CME fit) always
    # produce a RENEWAL process, so it is tagged PH, not MAP, exactly as MATLAB
    # updateSnForMAP does. Tagging MAP made the MAM handler treat the service as
    # a correlated arrival process and fall back to the exponential rate: an
    # M/Gamma/1 at rho 0.5 returned the M/M/1 answer 1.0 instead of 0.8125.
    # A surrogate that is NOT a phase-type is tagged ME and recorded in sn.isph,
    # so the CTMC assembles a rational generator and the PH-only consumers refuse.
    from .utils import sn_is_phasetype
    is_ph = sn_is_phasetype(MAP)
    sn.procid[ist, r] = ProcessType.PH if is_ph else ProcessType.ME
    if sn.isph is None:
        sn.isph = np.ones((sn.nstations, sn.nclasses), dtype=bool)
    sn.isph[ist, r] = is_ph

    # Update phases
    if sn.phases is None:
        sn.phases = np.ones((sn.nstations, sn.nclasses))
    sn.phases[ist, r] = n_phases

    # Update phasessz (integer dtype: state-space slicing indexes with these)
    if sn.phasessz is None:
        sn.phasessz = np.ones((sn.nstations, sn.nclasses), dtype=int)
    elif not np.issubdtype(sn.phasessz.dtype, np.integer):
        sn.phasessz = sn.phasessz.astype(int)
    sn.phasessz[ist, r] = max(int(n_phases), 1)

    # Recompute phaseshift for this station (cumulative sum across classes)
    if sn.phaseshift is None:
        sn.phaseshift = np.zeros((sn.nstations, sn.nclasses + 1), dtype=int)
    elif not np.issubdtype(sn.phaseshift.dtype, np.integer):
        sn.phaseshift = sn.phaseshift.astype(int)
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

    # see _kb/03-api-layer.md for rationale

    return sn


def zero_source_metrics(QN, UN, sn):
    """
    Zero the queue length and utilization of every Source station.

    A Source holds no jobs and occupies no server, so both are zero BY
    DISCIPLINE rather than by sign. Applied to the stored result, not to the
    table, so a caller reading the solver's own arrays sees the same thing the
    table does. Solvers otherwise leave arbitrary values in that row (measured
    on an open M/M/1: NC U=1, MAM Q=1, CTMC Q=Inf), which the table layer only
    ever suppressed incidentally, through a threshold on RESPONSE TIME.

    Args:
        QN: (nstations, nclasses) queue-length matrix, or None.
        UN: (nstations, nclasses) utilization matrix, or None.
        sn: NetworkStruct providing nodetype and nodeToStation.

    Returns:
        Tuple ``(QN, UN)``. Both are copied, not mutated.
    """
    if sn is None:
        return QN, UN
    nodetype = getattr(sn, 'nodetype', None)
    node_to_station = getattr(sn, 'nodeToStation', None)
    if nodetype is None or node_to_station is None:
        return QN, UN
    nodetype = np.asarray(nodetype).flatten()
    node_to_station = np.asarray(node_to_station).flatten()
    src = [int(node_to_station[i]) for i in range(min(len(nodetype), len(node_to_station)))
           if nodetype[i] == NodeType.SOURCE and node_to_station[i] >= 0]
    if not src:
        return QN, UN
    out = []
    for M in (QN, UN):
        if M is None:
            out.append(M)
            continue
        M = np.array(M, dtype=float, copy=True)
        if M.ndim == 2:
            M[[i for i in src if i < M.shape[0]], :] = 0.0
        out.append(M)
    return out[0], out[1]


def cap_unstable_open_util(UN, TN, sn):
    """
    Cap the reported utilization of unstable open queueing stations at 1.0.

    A finite-server queueing station serving an open class is unstable when its
    offered load ``rho = sum_r T[i,r] / (nservers[i] * rate[i,r]) >= 1``. Such a
    station is fully saturated, so LINE reports its utilization as 1.0, split
    across classes in proportion to their offered load. ``rho`` is recomputed
    from throughput and service rate so the cap is independent of whatever value
    the solver algorithm left in ``U`` (some leave 0, others the raw rho). Source
    (EXT) and infinite-server / delay (INF) stations are excluded: their
    "utilization" is the mean number of busy servers and may legitimately
    exceed 1.

    Args:
        UN: (nstations, nclasses) per-class utilization matrix.
        TN: (nstations, nclasses) per-class throughput matrix.
        sn: NetworkStruct providing njobs, nservers, rates, sched.

    Returns:
        Tuple ``(UN_capped, any_unstable)``. ``UN`` is copied, not mutated.
    """
    UN = np.array(UN, dtype=float, copy=True)
    TN = np.asarray(TN, dtype=float)
    if UN.ndim != 2 or UN.shape != TN.shape:
        return UN, False

    any_unstable = False
    for i, rho, rho_tot, _ in _unstable_open_stations(TN, sn):
        any_unstable = True
        UN[i, :] = rho / rho_tot  # station total capped to 1.0
    return UN, any_unstable


def _station_reneges(sn, i) -> bool:
    """Whether any class at station ``i`` declares reneging, i.e. a waiting job
    may abandon. Such a station has a bounded queue at every offered load."""
    cls = getattr(sn, 'impatienceClass', None)
    if cls is None:
        return False
    arr = np.asarray(cls)
    if arr.ndim != 2 or i >= arr.shape[0]:
        return False
    from ...lang.base import ImpatienceType
    return bool(np.any(arr[i, :] == int(ImpatienceType.RENEGING)))


def _unstable_open_stations(TN, sn):
    """Yield ``(i, rho, rho_tot, open_cls)`` for every finite-server queueing
    station i whose offered load from open classes is >= 1 (saturated).
    ``rho`` is the per-class offered load vector recomputed from throughput and
    service rate; ``open_cls`` is the boolean open-class mask. Source (EXT) and
    infinite-server (INF) stations are skipped."""
    TN = np.asarray(TN, dtype=float)
    if TN.ndim != 2:
        return

    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    open_cls = ~np.isfinite(njobs)
    if not np.any(open_cls):
        return

    nservers = np.asarray(sn.nservers, dtype=float).flatten()
    rates = np.asarray(sn.rates, dtype=float)
    M, K = TN.shape
    for i in range(M):
        c = nservers[i] if i < len(nservers) else 1.0
        if not np.isfinite(c) or c <= 0:
            continue  # infinite-server / delay station: never capped
        try:
            sched_name = getattr(sn.sched[i], 'name', None)
        except Exception:
            sched_name = None
        if sched_name in ('INF', 'EXT'):
            continue  # delay (INF) or source (EXT) station
        # A STATION CUSTOMERS ABANDON CANNOT BE UNSTABLE, however heavily it is
        # offered: reneging bounds the queue, the excess leaving instead of
        # accumulating. rho = T/(c*rate) reaches exactly 1 there -- that is what
        # a saturated server WITH abandonment looks like -- so without this the
        # Erlang A queue length that qsys_erlanga and qsys_ggisgi_fluid compute
        # exactly would be overwritten with Inf.
        if _station_reneges(sn, i):
            continue
        # A load-dependent station serves faster than its nominal rate; without
        # the peak scaling the test reads the rate at population one and calls a
        # stable station saturated (e.g. the discrete-time p(n)=p*min(n,s)
        # server of Daduna's example 2.10).
        lldpeak = 1.0
        lld = getattr(sn, 'lldscaling', None)
        if lld is not None and np.size(lld) > 0:
            arr = np.asarray(lld, dtype=float)
            if arr.ndim == 2 and arr.shape[0] > i:
                row = arr[i, :]
                row = row[np.isfinite(row) & (row > 0)]
                if row.size > 0:
                    lldpeak = float(row.max())
        rho = np.zeros(K)
        has_open = False
        for r in range(K):
            if (i < rates.shape[0] and r < rates.shape[1]
                    and rates[i, r] > 0 and TN[i, r] > 0):
                rho[r] = TN[i, r] / (c * lldpeak * rates[i, r])
                if open_cls[r]:
                    has_open = True
        rho_open = float(rho[open_cls].sum())
        rho_tot = float(rho.sum())
        if has_open and rho_open >= 1.0 and rho_tot > 0:
            yield i, rho, rho_tot, open_cls


def saturate_unstable_open_metrics(QN, RN, TN, sn):
    """
    Sanitize per-station metrics at unstable (saturated) open stations.

    Analytical open-queue formulas diverge when the offered load rho >= 1:
    fixed-point algorithms leave overflow-scale garbage in QN/RN. At every
    saturated station the open classes have unbounded queue length and response
    time, reported Inf.

    THROUGHPUT IS LEFT AS THE ANALYZER COMPUTED IT, i.e. the offered rate, which
    is what getAvg.m and the JAR both report (measured on
    test_gallery_hyperl1_feedback: Tput 10, not the service capacity 2). Rescaling
    it to the capacity was a python-only rule and made the station's Tput and
    ArvR describe two different flows.

    Args:
        QN, RN, TN: (nstations, nclasses) metric matrices (copied, not mutated).
        sn: NetworkStruct providing njobs, nservers, rates, sched.

    Returns:
        Tuple ``(QN, RN, TN, any_unstable)``.
    """
    QN = np.array(QN, dtype=float, copy=True)
    RN = np.array(RN, dtype=float, copy=True)
    TN = np.array(TN, dtype=float, copy=True)
    if QN.ndim != 2 or QN.shape != TN.shape:
        return QN, RN, TN, False

    any_unstable = False
    for i, rho, _, open_cls in _unstable_open_stations(TN, sn):
        any_unstable = True
        # Open classes saturate the station: unbounded backlog.
        for r in range(QN.shape[1]):
            if open_cls[r] and rho[r] > 0:
                QN[i, r] = np.inf
                if RN.shape == QN.shape:
                    RN[i, r] = np.inf
    return QN, RN, TN, any_unstable


def sn_rt_stations(sn: NetworkStruct):
    """Station-to-station routing probabilities and per-station visits.

    ``sn.rt`` and ``sn.visits`` are indexed by STATEFUL node, so a solver that
    writes traffic equations over stations and indexes them by station index
    silently reads the wrong rows as soon as the model owns a stateful node
    that is not a station (Router, Cache, stateful class switch). The returned
    routing matrix absorbs those nodes,

        Pst = P_AA + P_AB * (I - P_BB)^-1 * P_BA,

    with A the station rows in station order and B the remaining stateful rows,
    which is exact because a non-station stateful node holds no jobs: it passes
    every arrival on instantaneously. When every stateful node is a station the
    result is ``sn.rt`` unchanged.

    Args:
        sn: NetworkStruct

    Returns:
        Tuple ``(rt_stations, visits_stations)`` of shapes (M*K, M*K) and (M, K).
    """
    K = int(sn.nclasses)
    M = int(sn.nstations)
    S = int(sn.nstateful)
    stationToStateful = np.asarray(sn.stationToStateful).ravel().astype(int)

    is_st = np.zeros(S, dtype=bool)
    is_st[stationToStateful] = True

    A = np.concatenate([stationToStateful[ist] * K + np.arange(K) for ist in range(M)]) \
        if M > 0 else np.zeros(0, dtype=int)
    B = np.concatenate([isf * K + np.arange(K) for isf in range(S) if not is_st[isf]]) \
        if S > int(is_st.sum()) else np.zeros(0, dtype=int)

    P = np.asarray(sn.rt, dtype=float)
    if B.size == 0:
        rtst = P[np.ix_(A, A)]
    else:
        rtst = P[np.ix_(A, A)] + P[np.ix_(A, B)] @ np.linalg.solve(
            np.eye(B.size) - P[np.ix_(B, B)], P[np.ix_(B, A)])

    Vall = np.zeros((S, K))
    if sn.visits:
        for c in sn.visits:
            Vall = Vall + np.asarray(sn.visits[c], dtype=float).reshape(S, K)
    Vst = Vall[stationToStateful, :]
    return rtst, Vst
