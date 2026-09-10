"""
Flow-Equivalent Server (FES) Analysis.

Native Python implementations for Flow-Equivalent Server analysis,
which allows decomposition of complex queueing networks by replacing
subnetworks with equivalent load-dependent servers.

Key functions:
    fes_aggregate: Aggregate a station subset into a Flow-Equivalent Server
    fes_aggregate_with_options: Aggregate with custom options
    fes_validate: Validate FES topology
    fes_compute_throughputs: Compute throughputs using FES decomposition
    fes_build_isolated: Build isolated subsystem for FES analysis

References:
    Original MATLAB: matlab/src/api/fes/fes_*.m
    JAR: jar/src/main/java/jline/api/fes/FESAggregator.java
    Chandy, Herzog, Woo, "Parametric analysis of queueing networks", 1975
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any, List
from dataclasses import dataclass
from collections import deque


@dataclass
class FesResult:
    """Result of FES throughput computation."""
    throughput: np.ndarray  # Throughput per chain
    queue_length: np.ndarray  # Queue lengths
    response_time: np.ndarray  # Response times


def fes_aggregate(model, station_subset):
    """
    Aggregate a station subset into a Flow-Equivalent Server.

    Replaces the specified stations in a closed product-form queueing network
    with a single FES station. The FES has per-class class-dependent service
    rates beta_{i,r}(n) computed from MVA throughputs of the isolated subnetwork
    for all reachable population vectors.

    Args:
        model: LINE Network model (Python native).
            Must be a closed product-form queueing network.
        station_subset: List of Station objects to aggregate.

    Returns:
        dict with keys:
            - fes_model: New Network with the FES station replacing the subset.
            - fes_station: Reference to the FES Queue station in the new model.
            - deagg_info (dict): Deaggregation information.
    """
    return fes_aggregate_with_options(model, station_subset)


def fes_aggregate_with_options(model, station_subset, solver='mva',
                                cutoffs=None, verbose=False):
    """
    Aggregate a station subset into a Flow-Equivalent Server with custom options.

    Args:
        model: LINE Network model (Python native).
        station_subset: List of Station objects to aggregate.
        solver (str): Solver for throughput computation. Default: 'mva'.
        cutoffs (numpy.ndarray or None): Per-class population cutoffs (1 x K).
            If None, uses the total population per class.
        verbose (bool): If True, print progress. Default: False.

    Returns:
        dict with keys:
            - fes_model: New Network with the FES station.
            - fes_station: Reference to the FES Queue station.
            - deagg_info (dict): Deaggregation information.
    """
    from ...lang.network import Network
    from ...lang.nodes import Queue, Delay
    from ...lang.classes import ClosedClass
    from ...lang.base import SchedStrategy
    from ...distributions.continuous import Exp, Disabled
    from ..pfqn.ljd import ljd_linearize
    from line_solver.lib.kpctoolbox.mc import dtmc_stochcomp

    # Validate inputs
    _validate_inputs(model, station_subset)

    # Get network structure
    sn = model.get_struct()
    M = sn.nstations
    K = sn.nclasses

    # Get station indices for subset
    model_stations = model.get_stations()
    subset_indices = np.zeros(len(station_subset), dtype=int)
    for i, sub_station in enumerate(station_subset):
        for j in range(M):
            if sub_station is model_stations[j]:
                subset_indices[i] = j
                break

    # Get complement indices
    all_indices = set(range(M))
    complement_indices = np.array(sorted(all_indices - set(subset_indices)), dtype=int)

    # Set cutoffs
    njobs = np.asarray(sn.njobs).flatten()
    if cutoffs is not None:
        cutoffs = np.asarray(cutoffs).flatten()
    else:
        cutoffs = njobs.copy()
    cutoffs = cutoffs.astype(int)

    # Build rt index sets (0-based in Python native)
    rt = np.asarray(sn.rt)
    station_to_stateful = np.asarray(sn.stationToStateful).astype(int)

    subset_rt_indices = []
    for i in subset_indices:
        isf = station_to_stateful[i]
        for k in range(K):
            subset_rt_indices.append(isf * K + k)

    complement_rt_indices = []
    for i in complement_indices:
        isf = station_to_stateful[i]
        for k in range(K):
            complement_rt_indices.append(isf * K + k)

    # Compute stochastic complements
    stoch_comp_subset = dtmc_stochcomp(rt, np.array(subset_rt_indices, dtype=int))
    stoch_comp_complement = dtmc_stochcomp(rt, np.array(complement_rt_indices, dtype=int))

    # Build isolated subnetwork
    isolated_model = _build_isolated_model(
        model, station_subset, subset_indices, stoch_comp_subset, sn
    )

    # Compute throughputs for all population states
    throughput_table = _compute_throughputs(
        isolated_model, cutoffs, sn, solver, verbose
    )

    # Create the FES model
    model_name = model.name
    fes_model = Network(f"{model_name}_FES")

    # Copy complement stations
    station_map = {}
    model_classes = model.get_classes()
    rates = np.asarray(sn.rates)

    for i in complement_indices:
        orig_station = model_stations[i]
        if isinstance(orig_station, Delay):
            new_station = Delay(fes_model, orig_station.name)
        elif isinstance(orig_station, Queue):
            sched = orig_station.get_sched_strategy()
            new_station = Queue(fes_model, orig_station.name, sched)
            nservers = orig_station.getNumberOfServers()
            if nservers < float('inf'):
                new_station.setNumberOfServers(nservers)
        else:
            raise ValueError(f"Unsupported station type: {type(orig_station)}")
        station_map[i] = new_station

    # Create FES station (PS, 1 server)
    fes_station = Queue(fes_model, "FES", SchedStrategy.PS)
    fes_station.setNumberOfServers(1)

    # Create job classes
    ref_station = station_map[complement_indices[0]] if len(complement_indices) > 0 else fes_station
    new_classes = []
    for k in range(K):
        orig_class = model_classes[k]
        population = int(njobs[k])
        class_name = orig_class.name
        new_class = ClosedClass(fes_model, class_name, population, ref_station, 0)
        new_classes.append(new_class)

    # Set service distributions for complement stations
    for i in complement_indices:
        new_station = station_map[i]
        for k in range(K):
            rate = rates[i, k]
            if rate > 0 and not np.isnan(rate):
                new_station.setService(new_classes[k], Exp(rate))
            else:
                new_station.setService(new_classes[k], Disabled())

    # Set FES base service rate Exp(1.0); the class-dependence handle provides
    # the actual rate
    for k in range(K):
        fes_station.setService(new_classes[k], Exp(1.0))

    # Store throughputs directly as the class-dependent rates (same as MATLAB).
    # The convolution algorithm uses X_m(n) = (L / mu_km(n)) * X_m(n-e_k)
    # where mu_km(n) = throughput of the subnetwork at population n.
    table_size = _compute_table_size(cutoffs)
    scaling_table = []
    for k in range(K):
        st_table = np.zeros(table_size)
        for idx in range(table_size):
            tput = throughput_table[k][idx]
            if tput > 1e-10:
                st_table[idx] = tput
            else:
                st_table[idx] = 1e-10  # Small positive value for zero throughput
        scaling_table.append(st_table)

    # Install the FES rates as a per-class class-dependence function
    # beta_{i,r}(n) = X_r(n), i.e. Sauer's chain-dependent service rate
    # mu_{r,i}(n) (Sauer 1983, eq. (40)). The handle indexes the linearized
    # table and clamps the population to the cutoffs it was built on.
    from ..pfqn.ljd import fes_beta_handle
    from ..pfqn.conv import _cd_peak_scaling
    fes_beta = fes_beta_handle(scaling_table, cutoffs)
    # Peak per-class FES rate = normalizer for Util = T*S/peak (the FES rates are
    # Sauer's chain-dependent service rates; their lattice peak plays the role of
    # the effective server count).
    import numpy as _np
    _cut = _np.asarray(cutoffs, dtype=int) if cutoffs is not None else _np.ones(len(scaling_table), dtype=int)
    fes_peak = _cd_peak_scaling(fes_beta, _cut, len(scaling_table))
    fes_station.set_class_dependence(fes_beta, fes_peak)

    # Build routing matrix
    nodes = fes_model.get_nodes()
    n_nodes = len(nodes)

    # Find FES node index
    fes_node_idx = 0
    for n in range(n_nodes):
        if nodes[n] is fes_station:
            fes_node_idx = n
            break

    # Build node index mapping for complement stations
    complement_node_map = {}
    for i in complement_indices:
        station = station_map[i]
        for n in range(n_nodes):
            if nodes[n] is station:
                complement_node_map[i] = n
                break

    # Set routing probabilities
    # Use ORIGINAL routing matrix (not stochastic complement) for the FES model.
    # The stochastic complement absorbs the subset into self-loops on the complement,
    # which is wrong for FES routing (we replace subset with FES, not absorb it).
    P = fes_model.init_routing_matrix()

    M_sub = len(subset_indices)
    M_comp = len(complement_indices)

    for k in range(K):
        Pk = np.zeros((n_nodes, n_nodes))

        # 1. Routes within complement: use ORIGINAL routing matrix
        for i in complement_indices:
            if i not in complement_node_map:
                continue
            i_node = complement_node_map[i]
            isf_i = station_to_stateful[i]

            for j in complement_indices:
                if j not in complement_node_map:
                    continue
                j_node = complement_node_map[j]
                isf_j = station_to_stateful[j]

                rt_idx_i = isf_i * K + k
                rt_idx_j = isf_j * K + k

                if rt_idx_i < rt.shape[0] and rt_idx_j < rt.shape[1]:
                    prob = rt[rt_idx_i, rt_idx_j]
                    if prob > 1e-10:
                        Pk[i_node, j_node] = prob

        # 2. Routes from complement to FES: sum of original routes to subset stations
        for i in complement_indices:
            if i not in complement_node_map:
                continue
            i_node = complement_node_map[i]
            isf_i = station_to_stateful[i]

            for j in subset_indices:
                isf_j = station_to_stateful[j]
                rt_idx_i = isf_i * K + k
                rt_idx_j = isf_j * K + k

                if rt_idx_i < rt.shape[0] and rt_idx_j < rt.shape[1]:
                    prob = rt[rt_idx_i, rt_idx_j]
                    if prob > 1e-10:
                        Pk[i_node, fes_node_idx] += prob

        # 3. Routes from FES to complement: absorption matrix approach
        # Build P_SS (subset→subset) and P_SC (subset→complement) for this class
        P_SS = np.zeros((M_sub, M_sub))
        P_SC = np.zeros((M_sub, M_comp))

        for si, i in enumerate(subset_indices):
            isf_i = station_to_stateful[i]
            for sj, j in enumerate(subset_indices):
                isf_j = station_to_stateful[j]
                rt_idx_i = isf_i * K + k
                rt_idx_j = isf_j * K + k
                if rt_idx_i < rt.shape[0] and rt_idx_j < rt.shape[1]:
                    P_SS[si, sj] = rt[rt_idx_i, rt_idx_j]

            for cj, j in enumerate(complement_indices):
                isf_j = station_to_stateful[j]
                rt_idx_i = isf_i * K + k
                rt_idx_j = isf_j * K + k
                if rt_idx_i < rt.shape[0] and rt_idx_j < rt.shape[1]:
                    P_SC[si, cj] = rt[rt_idx_i, rt_idx_j]

        # Absorption matrix: B = (I - P_SS)^{-1} * P_SC
        # B[i,j] = prob that starting at subset station i, job exits to complement station j
        I_sub = np.eye(M_sub)
        try:
            fundamental = np.linalg.solve(I_sub - P_SS, P_SC)
        except np.linalg.LinAlgError:
            fundamental = np.linalg.lstsq(I_sub - P_SS, P_SC, rcond=None)[0]

        # Entry weights: how much traffic enters each subset station from complement
        entry_weights = np.zeros(M_sub)
        for si, j in enumerate(subset_indices):
            isf_j = station_to_stateful[j]
            for i in complement_indices:
                isf_i = station_to_stateful[i]
                rt_idx_i = isf_i * K + k
                rt_idx_j = isf_j * K + k
                if rt_idx_i < rt.shape[0] and rt_idx_j < rt.shape[1]:
                    entry_weights[si] += rt[rt_idx_i, rt_idx_j]

        total_entry = np.sum(entry_weights)
        if total_entry > 1e-10:
            entry_weights /= total_entry

        # FES exit probability = weighted average of absorption probabilities
        fes_exit = entry_weights @ fundamental  # (M_comp,)

        for cj, j in enumerate(complement_indices):
            if j not in complement_node_map:
                continue
            j_node = complement_node_map[j]
            if fes_exit[cj] > 1e-10:
                Pk[fes_node_idx, j_node] = fes_exit[cj]

        # Normalize rows (should already sum to 1 but ensure numerical stability)
        for n in range(n_nodes):
            row_sum = np.sum(Pk[n, :])
            if row_sum > 1e-10:
                Pk[n, :] /= row_sum

        # Set routing for this class (same-class routing, no class switching)
        P.set(new_classes[k], new_classes[k], Pk)

    # Link the model
    fes_model.link(P)

    # Build deaggregation info
    deagg_info = {
        'original_model': model,
        'station_subset': station_subset,
        'subset_indices': subset_indices,
        'complement_indices': complement_indices,
        'throughput_table': throughput_table,
        'cutoffs': cutoffs,
        'stoch_comp_subset': stoch_comp_subset,
        'stoch_comp_complement': stoch_comp_complement,
        'isolated_model': isolated_model,
        'fes_node_idx': fes_node_idx,
    }

    return {
        'fes_model': fes_model,
        'fes_station': fes_station,
        'deagg_info': deagg_info,
    }


def _validate_inputs(model, station_subset):
    """Validate inputs for FES aggregation."""
    from ...lang.nodes import Queue, Delay

    if not station_subset:
        raise RuntimeError("Station subset cannot be empty")

    sn = model.get_struct()
    njobs = np.asarray(sn.njobs).flatten()

    # Check closed network
    for k in range(sn.nclasses):
        if np.isinf(njobs[k]):
            raise RuntimeError("FES aggregation only applies to closed queueing networks")

    model_stations = model.get_stations()
    if len(station_subset) >= len(model_stations):
        raise RuntimeError("Cannot aggregate all stations")

    # Validate each station
    for station in station_subset:
        if not isinstance(station, (Queue, Delay)):
            raise RuntimeError("FES aggregation only supports Queue and Delay stations")
        found = False
        for ms in model_stations:
            if station is ms:
                found = True
                break
        if not found:
            station_name = station.name
            raise RuntimeError(f"Station {station_name} does not belong to the model")


def _build_isolated_model(model, station_subset, subset_indices, stoch_comp_s, sn):
    """Build isolated subnetwork from station subset."""
    from ...lang.network import Network
    from ...lang.nodes import Queue, Delay
    from ...lang.classes import ClosedClass
    from ...lang.base import SchedStrategy
    from ...distributions.continuous import Exp, Disabled

    K = sn.nclasses
    M_sub = len(station_subset)
    rates = np.asarray(sn.rates)
    model_classes = model.get_classes()
    model_name = model.name

    isolated_model = Network(f"{model_name}_isolated")

    # Create new stations
    new_stations = []
    for i in range(M_sub):
        orig_station = station_subset[i]
        if isinstance(orig_station, Delay):
            new_station = Delay(isolated_model, orig_station.name)
        elif isinstance(orig_station, Queue):
            sched = orig_station.get_sched_strategy()
            new_station = Queue(isolated_model, orig_station.name, sched)
            nservers = orig_station.getNumberOfServers()
            if nservers < float('inf'):
                new_station.setNumberOfServers(nservers)
        else:
            raise RuntimeError("Unsupported station type")
        new_stations.append(new_station)

    # Create closed classes with population 1 (updated during enumeration)
    ref_station = new_stations[0]
    new_classes = []
    for k in range(K):
        orig_class = model_classes[k]
        class_name = orig_class.name
        new_class = ClosedClass(isolated_model, class_name, 1, ref_station, 0)
        new_classes.append(new_class)

    # Set service distributions
    for i in range(M_sub):
        orig_idx = subset_indices[i]
        new_station = new_stations[i]
        for k in range(K):
            rate = rates[orig_idx, k]
            if rate > 0 and not np.isnan(rate):
                new_station.setService(new_classes[k], Exp(rate))
            else:
                new_station.setService(new_classes[k], Disabled())

    # Build routing from stochastic complement
    nodes = isolated_model.get_nodes()
    n_nodes = len(nodes)

    station_to_node = np.zeros(M_sub, dtype=int)
    for i in range(M_sub):
        for n in range(n_nodes):
            if nodes[n] is new_stations[i]:
                station_to_node[i] = n
                break

    P = isolated_model.init_routing_matrix()

    for k in range(K):
        Pk = np.zeros((n_nodes, n_nodes))

        for i in range(M_sub):
            row_idx = i * K + k
            i_node = station_to_node[i]

            for j in range(M_sub):
                col_idx = j * K + k
                j_node = station_to_node[j]

                if row_idx < stoch_comp_s.shape[0] and col_idx < stoch_comp_s.shape[1]:
                    prob = stoch_comp_s[row_idx, col_idx]
                    if prob > 1e-10:
                        Pk[i_node, j_node] = prob

            # Normalize row
            row_sum = np.sum(Pk[i_node, :])
            if row_sum > 1e-10:
                Pk[i_node, :] /= row_sum
            else:
                # Self-loop if no routing
                Pk[i_node, i_node] = 1.0

        P.set(new_classes[k], new_classes[k], Pk)

    isolated_model.link(P)
    return isolated_model


def _compute_throughputs(isolated_model, cutoffs, sn, solver_name, verbose):
    """Compute throughputs for all population states using BFS enumeration."""
    from ..pfqn.ljd import ljd_linearize
    from ..pfqn.mvald import pfqn_mvams

    K = sn.nclasses
    table_size = _compute_table_size(cutoffs)

    # Initialize per-class throughput tables
    throughput_table = []
    for k in range(K):
        throughput_table.append(np.zeros(table_size))

    # Get isolated model structure for pfqn_mva
    iso_sn = isolated_model.get_struct()
    iso_rates = np.asarray(iso_sn.rates)
    iso_nservers = np.asarray(iso_sn.nservers).flatten()
    iso_M = iso_sn.nstations

    # Separate Delay (INF servers) from Queue stations for pfqn_mva
    # pfqn_mva takes L (demands at queues), N (population), Z (think times)
    from ...lang.base import SchedStrategy
    delay_indices = []
    queue_indices = []
    for i in range(iso_M):
        if iso_sn.sched.get(i) == SchedStrategy.INF or np.isinf(iso_nservers[i]):
            delay_indices.append(i)
        else:
            queue_indices.append(i)

    # BFS enumeration of all population states
    visited = set()
    queue = deque()
    zero_state = tuple([0] * K)
    queue.append(zero_state)

    while queue:
        nvec = queue.popleft()
        if nvec in visited:
            continue
        visited.add(nvec)

        nvec_arr = np.array(nvec, dtype=int)
        total_pop = sum(nvec)

        idx = ljd_linearize(nvec_arr, cutoffs)

        if total_pop == 0:
            # Zero population: throughput is 0
            for k in range(K):
                throughput_table[k][idx - 1] = 0.0  # ljd_linearize returns 1-based
        else:
            # Build demands and think times for the MVA solve. A failure is
            # NOT caught: an FES table silently filled with zeros is a wrong
            # aggregate, not a degraded one, and every downstream beta_r(n)
            # reads it as "the subnetwork serves nothing".
            # L = service demands at queue stations (M_queues x K)
            # Z = think times at delay stations (1 x K)
            N = nvec_arr.astype(float)

            # Compute visits for the isolated model
            # For a closed network, visits are derived from the routing matrix
            iso_sn_fresh = isolated_model.get_struct()
            iso_visits = None
            if hasattr(iso_sn_fresh, 'visits') and iso_sn_fresh.visits is not None:
                # visits is a list of per-chain visit matrices
                # For single-chain (no class switching), aggregate visits
                nchains = len(iso_sn_fresh.visits)
                iso_visits = np.zeros((iso_sn_fresh.nstateful, K))
                for c in range(nchains):
                    if iso_sn_fresh.visits[c] is not None:
                        v = np.asarray(iso_sn_fresh.visits[c])
                        iso_visits[:v.shape[0], :v.shape[1]] += v

            # Build L (demands) and Z (think times) from rates and visits
            if iso_visits is not None:
                # Demands = visits / rate (i.e., visits * service_time)
                Z = np.zeros(K)
                L_list = []
                mi_list = []
                for i in queue_indices:
                    isf = int(iso_sn_fresh.stationToStateful[i])
                    demands = np.zeros(K)
                    for k_cls in range(K):
                        rate = iso_rates[i, k_cls]
                        visit = iso_visits[isf, k_cls] if isf < iso_visits.shape[0] else 0.0
                        if rate > 0 and not np.isnan(rate) and visit > 0:
                            demands[k_cls] = visit / rate
                    L_list.append(demands)
                    mi_list.append(iso_nservers[i])

                for i in delay_indices:
                    isf = int(iso_sn_fresh.stationToStateful[i])
                    for k_cls in range(K):
                        rate = iso_rates[i, k_cls]
                        visit = iso_visits[isf, k_cls] if isf < iso_visits.shape[0] else 0.0
                        if rate > 0 and not np.isnan(rate) and visit > 0:
                            Z[k_cls] += visit / rate

                if len(L_list) > 0:
                    L = np.array(L_list)
                    mi = np.array(mi_list)
                else:
                    # All delays, no queues
                    L = np.zeros((1, K))
                    mi = np.array([1.0])
                    # With only delays, throughput = N / Z
                    for k_cls in range(K):
                        if Z[k_cls] > 0 and N[k_cls] > 0:
                            throughput_table[k_cls][idx - 1] = N[k_cls] / Z[k_cls]
                        else:
                            throughput_table[k_cls][idx - 1] = 0.0
                    # Skip pfqn_mva call
                    # Add neighboring states
                    for k_cls in range(K):
                        if nvec[k_cls] < cutoffs[k_cls]:
                            new_state = list(nvec)
                            new_state[k_cls] += 1
                            new_state_t = tuple(new_state)
                            if new_state_t not in visited:
                                queue.append(new_state_t)
                    continue

                # pfqn_mva's `mi` is NOT a server count -- it enters only as
                # the additive term of C(i,s)=L(i,s)*(mi(i)+Qarv), so passing
                # the real multiplicity INFLATES the residence time instead of
                # adding servers. pfqn_mvams forwards to pfqn_mva when every
                # station is a single server and to the load-dependent
                # recursion with mu(i,n)=min(n,S(i)) when one is not.
                # See _kb/03-api-layer.md (pfqn_mva: mi is not S).
                XN, QN, UN, CN, _lG = pfqn_mvams(
                    np.zeros(K), L, N, Z, np.ones(L.shape[0]), mi)

                # Extract system throughputs
                XN_flat = np.asarray(XN).flatten()
                for k_cls in range(K):
                    if nvec[k_cls] > 0 and k_cls < len(XN_flat):
                        throughput_table[k_cls][idx - 1] = XN_flat[k_cls]
                    else:
                        throughput_table[k_cls][idx - 1] = 0.0
            else:
                # Fallback: use rates directly as demands (visits = 1)
                Z = np.zeros(K)
                L_list = []
                mi_list = []
                for i in queue_indices:
                    demands = np.zeros(K)
                    for k_cls in range(K):
                        rate = iso_rates[i, k_cls]
                        if rate > 0 and not np.isnan(rate):
                            demands[k_cls] = 1.0 / rate
                    L_list.append(demands)
                    mi_list.append(iso_nservers[i])

                for i in delay_indices:
                    for k_cls in range(K):
                        rate = iso_rates[i, k_cls]
                        if rate > 0 and not np.isnan(rate):
                            Z[k_cls] += 1.0 / rate

                L = np.array(L_list) if L_list else np.zeros((1, K))
                mi = np.array(mi_list) if mi_list else np.array([1.0])

                # `mi` is the additive C=L*(mi+Qarv) term, not a server
                # count: multiservers go through pfqn_mvams (see above)
                XN, QN, UN, CN, _lG = pfqn_mvams(
                    np.zeros(K), L, N, Z, np.ones(L.shape[0]), mi)
                XN_flat = np.asarray(XN).flatten()
                for k_cls in range(K):
                    if nvec[k_cls] > 0 and k_cls < len(XN_flat):
                        throughput_table[k_cls][idx - 1] = XN_flat[k_cls]
                    else:
                        throughput_table[k_cls][idx - 1] = 0.0

            if verbose:
                print(f"  FES: nvec={nvec_arr}, X={[throughput_table[k][idx-1] for k in range(K)]}")


        # Add neighboring states
        for k_cls in range(K):
            if nvec[k_cls] < cutoffs[k_cls]:
                new_state = list(nvec)
                new_state[k_cls] += 1
                new_state_t = tuple(new_state)
                if new_state_t not in visited:
                    queue.append(new_state_t)

    return throughput_table


def _compute_table_size(cutoffs):
    """Compute total table size from cutoffs."""
    size = 1
    for k in range(len(cutoffs)):
        size *= (int(cutoffs[k]) + 1)
    return size


def fes_validate(sn) -> bool:
    """
    Validate that a network is suitable for FES analysis.

    Args:
        sn: Network structure

    Returns:
        True if network is valid for FES analysis
    """
    if hasattr(sn, 'classswitch') and sn.classswitch is not None:
        if np.any(sn.classswitch):
            return False

    if hasattr(sn, 'fj') and sn.fj is not None:
        if np.any(sn.fj):
            return False

    if hasattr(sn, 'nclosedjobs'):
        if sn.nclosedjobs == 0:
            return False

    return True


def fes_compute_throughputs(sn, mu_fes: np.ndarray,
                             population: Optional[np.ndarray] = None
                             ) -> FesResult:
    """
    Compute throughputs using FES decomposition.

    Args:
        sn: Network structure
        mu_fes: Load-dependent service rates of the FES (chain x population)
        population: Population vector (default: from sn)

    Returns:
        FesResult with throughput, queue_length, response_time
    """
    if population is None:
        if hasattr(sn, 'njobs'):
            population = sn.njobs
        else:
            raise ValueError("Population must be provided")

    population = np.asarray(population, dtype=int).flatten()
    mu_fes = np.atleast_2d(np.asarray(mu_fes, dtype=float))

    nchains = len(population)

    if nchains == 1:
        N = population[0]
        L = np.zeros(N + 1)
        W = np.zeros(N + 1)
        X = np.zeros(N + 1)

        for n in range(1, N + 1):
            mu_n = mu_fes[0, n - 1] if mu_fes.shape[1] >= n else mu_fes[0, -1]
            if mu_n > 0:
                W[n] = (1 + L[n - 1]) / mu_n
            else:
                W[n] = np.inf
            X[n] = n / W[n] if W[n] > 0 else 0
            L[n] = X[n] * W[n]

        return FesResult(
            throughput=np.array([X[N]]),
            queue_length=np.array([L[N]]),
            response_time=np.array([W[N]])
        )
    else:
        X = np.zeros(nchains)
        Q = np.zeros(nchains)
        W = np.zeros(nchains)

        for c in range(nchains):
            N_c = population[c]
            if N_c > 0 and mu_fes.shape[0] > c:
                mu_n = mu_fes[c, N_c - 1] if mu_fes.shape[1] >= N_c else mu_fes[c, -1]
                if mu_n > 0:
                    W[c] = 1.0 / mu_n
                    X[c] = N_c / W[c] if W[c] > 0 else 0
                    Q[c] = X[c] * W[c]

        return FesResult(throughput=X, queue_length=Q, response_time=W)


def fes_build_isolated(sn, subsystem_nodes: np.ndarray) -> Dict[str, Any]:
    """
    Build an isolated subsystem for FES analysis.

    Args:
        sn: Network structure
        subsystem_nodes: Node indices belonging to the subsystem

    Returns:
        Dictionary with subsystem parameters
    """
    subsystem_nodes = np.asarray(subsystem_nodes, dtype=int).flatten()

    stations = []
    for node in subsystem_nodes:
        if hasattr(sn, 'nodeToStation') and sn.nodeToStation[node] >= 0:
            stations.append(sn.nodeToStation[node])

    stations = np.array(stations, dtype=int)

    if len(stations) == 0:
        return {'D': np.array([]), 'Z': np.array([]), 'N': np.array([]), 'nservers': np.array([])}

    if hasattr(sn, 'proc') and sn.proc is not None:
        D = np.zeros((len(stations), sn.nclasses))
        for i, ist in enumerate(stations):
            for r in range(sn.nclasses):
                if sn.proc[ist, r] is not None:
                    D[i, r] = sn.proc[ist, r].get('mean', 0)
    else:
        D = np.zeros((len(stations), sn.nclasses))

    Z = np.zeros(sn.nclasses)
    if hasattr(sn, 'nodetype'):
        from ...lang.base import NodeType
        for node in subsystem_nodes:
            if sn.nodetype[node] == NodeType.Delay:
                ist = sn.nodeToStation[node]
                if ist >= 0 and hasattr(sn, 'proc'):
                    for r in range(sn.nclasses):
                        if sn.proc[ist, r] is not None:
                            Z[r] += sn.proc[ist, r].get('mean', 0)

    N = sn.njobs if hasattr(sn, 'njobs') else np.zeros(sn.nclasses)
    nservers = sn.nservers[stations] if hasattr(sn, 'nservers') else np.ones(len(stations))

    return {
        'D': D,
        'Z': Z,
        'N': N,
        'nservers': nservers,
        'stations': stations
    }


from .map_fes import (
    fes_map_levels,
    fes_map_interdeparture,
    fes_map_euler,
    fes_map_moments,
    fes_map_grid,
    fes_map_interp,
    fes_map_aggregate,
    fes_map_solve,
    fes_map_deaggregate,
)

__all__ = [
    'FesResult',
    'fes_map_levels',
    'fes_map_interdeparture',
    'fes_map_euler',
    'fes_map_moments',
    'fes_map_grid',
    'fes_map_interp',
    'fes_map_aggregate',
    'fes_map_solve',
    'fes_map_deaggregate',
    'fes_aggregate',
    'fes_aggregate_with_options',
    'fes_validate',
    'fes_compute_throughputs',
    'fes_build_isolated',
]


def fes_compute_metrics(isolated_model, cutoffs, K):
    """Per-station metrics of the ISOLATED subnetwork at every population state.

    The companion of :func:`_compute_throughputs`. That function returns the
    aggregate throughput X(n) that becomes the flow-equivalent server's rate,
    which is all the REDUCED model needs; it is not enough to report the
    COLLAPSED stations' own metrics. Those are recovered by conditioning on the
    FES population,

        E[Q_i] = sum_n P(N_fes = n) * Q_i(n),

    the Chandy-Herzog-Woo hierarchical decomposition, which is EXACT when the
    subnetwork is product-form. This supplies the Q_i(n) and U_i(n) the sum is
    taken over, indexed by ljd_linearize exactly as the throughput table is.
    Throughput needs no table: flow is fixed by the ROUTING, so the caller
    derives it from the original model's visit ratios.

    Args:
        isolated_model: the isolated subnetwork Network built by the transform
        cutoffs: per-class population cutoffs
        K: number of classes

    Returns:
        (Qtable, Utable): lists indexed by the 0-based linearized population
        state, each entry an (M_sub x K) array over the subnetwork's stations
        in their own order.

    References:
        MATLAB: matlab/src/api/fes/fes_compute_metrics.m
    """
    import numpy as _np

    from ..pfqn.ljd import ljd_linearize
    from ..pfqn.mvald import pfqn_mvams
    from ...lang.base import SchedStrategy

    cutoffs = _np.asarray(cutoffs, dtype=int).flatten()
    table_size = int(_np.prod(cutoffs + 1))

    iso_sn = isolated_model.get_struct()
    iso_rates = _np.asarray(iso_sn.rates)
    iso_nservers = _np.asarray(iso_sn.nservers).flatten()
    M_sub = iso_sn.nstations

    delay_indices, queue_indices = [], []
    for i in range(M_sub):
        if iso_sn.sched.get(i) == SchedStrategy.INF or _np.isinf(iso_nservers[i]):
            delay_indices.append(i)
        else:
            queue_indices.append(i)

    iso_visits = _np.zeros((iso_sn.nstateful, K))
    if getattr(iso_sn, 'visits', None) is not None:
        for c in range(len(iso_sn.visits)):
            if iso_sn.visits[c] is not None:
                v = _np.asarray(iso_sn.visits[c])
                iso_visits[:v.shape[0], :v.shape[1]] += v

    def _demand(i, k):
        isf = int(iso_sn.stationToStateful[i])
        rate = iso_rates[i, k]
        visit = iso_visits[isf, k] if isf < iso_visits.shape[0] else 0.0
        if rate > 0 and not _np.isnan(rate) and visit > 0:
            return visit / rate
        return 0.0

    L = _np.array([[_demand(i, k) for k in range(K)] for i in queue_indices]) \
        if queue_indices else _np.zeros((0, K))
    mi = _np.array([iso_nservers[i] for i in queue_indices]) if queue_indices else _np.zeros(0)
    Z = _np.zeros(K)
    for i in delay_indices:
        for k in range(K):
            Z[k] += _demand(i, k)

    Qtable = [_np.zeros((M_sub, K)) for _ in range(table_size)]
    Utable = [_np.zeros((M_sub, K)) for _ in range(table_size)]

    for idx0 in range(table_size):
        nvec = _ljd_delinearize(idx0 + 1, cutoffs)
        if nvec.sum() == 0:
            continue  # an empty subnetwork holds nothing and serves nothing
        # `mi` is the additive C=L*(mi+Qarv) term of pfqn_mva, not a server
        # count: multiservers go through pfqn_mvams. Its UN is per STATION on
        # the closed multiserver branch and per station-class elsewhere, so it
        # is not read here -- utilization is recomputed analytically as U=X*L/S,
        # the [0,1] convention LINE uses at every queueing station whatever its
        # multiplicity. See _kb/03-api-layer.md (pfqn_mva: mi is not S).
        XN, QN, _UN, _CN, _lG = pfqn_mvams(
            _np.zeros(K), L, nvec.astype(float), Z,
            _np.ones(L.shape[0]), mi)
        Q = _np.zeros((M_sub, K))
        U = _np.zeros((M_sub, K))
        QN = _np.atleast_2d(_np.asarray(QN))
        XN_row = _np.asarray(XN).flatten()
        for a, i in enumerate(queue_indices):
            if a < QN.shape[0]:
                Q[i, :QN.shape[1]] = QN[a, :]
            srv = max(1.0, float(mi[a])) if a < len(mi) else 1.0
            U[i, :] = XN_row[:K] * L[a, :] / srv
        # A Delay holds X_k * Z_i(k) jobs by Little's law on a station with no
        # queueing, and its INF "utilization" is that same population.
        XN_flat = _np.asarray(XN).flatten()
        for i in delay_indices:
            for k in range(K):
                Q[i, k] = XN_flat[k] * _demand(i, k) if k < len(XN_flat) else 0.0
                U[i, k] = Q[i, k]
        Qtable[idx0] = Q
        Utable[idx0] = U

    return Qtable, Utable


def _ljd_delinearize(idx, cutoffs):
    """Inverse of ljd_linearize: the population vector behind a linear index.

    The forward map is a mixed-radix numeral with class k in radix (Nk+1), so
    the inverse is the digit-by-digit division that reads it back.
    """
    import numpy as _np

    cutoffs = _np.asarray(cutoffs, dtype=int).flatten()
    nvec = _np.zeros(len(cutoffs), dtype=int)
    rem = int(idx) - 1
    for k in range(len(cutoffs)):
        radix = int(cutoffs[k]) + 1
        nvec[k] = rem % radix
        rem //= radix
    return nvec
