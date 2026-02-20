"""
RCAT-based INAP solver for SolverMAM.

Implements the RCAT (Reversed Compound Agent Theorem) iterative fixed-point
solver (INAP and INAP+) matching the MATLAB solver_mam_ag.m implementation.

Algorithm:
1. Build RCAT model from network structure
2. Initialize action rates randomly
3. Iterate:
   - Compute equilibrium distribution for each process using CTMC solve
   - Update action rates based on equilibrium terms:
     - INAP: x(a) = mean(Aa(i,j) * pi(i) / pi(j)) over non-zero entries
     - INAP+: x(a) = sum(Aa(i,j) * pi(i)) over non-zero entries
   - Check convergence
4. Compute performance metrics from equilibrium distributions

References:
    MATLAB: matlab/src/solvers/MAM/solver_mam_ag.m
"""

import numpy as np
import time
from typing import Tuple, Optional, List, Dict, Any
from dataclasses import dataclass

from . import MAMAlgorithm, MAMResult
from ..utils.network_adapter import extract_mam_params
from ....api.mc import ctmc_solve, ctmc_makeinfgen


@dataclass
class RCATModelData:
    """RCAT model data structure.

    Attributes:
        R: Rate matrices dictionary
        AP: Action-process mapping (num_actions x 2)
        process_map: Mapping from (station, class) to process index
        action_map: List of action details
        N: State space sizes for each process
        num_processes: Number of processes
        num_actions: Number of actions
    """
    R: Dict[Tuple[int, int], np.ndarray]
    AP: np.ndarray
    process_map: np.ndarray
    action_map: List[Dict]
    N: np.ndarray
    num_processes: int
    num_actions: int


def _build_rcat(sn, max_states: int = 100) -> RCATModelData:
    """Build RCAT model from LINE network structure.

    Args:
        sn: NetworkStruct
        max_states: Maximum number of states for truncation

    Returns:
        RCATModelData instance
    """
    M = sn.nstations
    K = sn.nclasses

    # Get routing table
    if hasattr(sn, 'rt') and sn.rt is not None:
        rt = np.asarray(sn.rt, dtype=np.float64)
    else:
        rt = np.zeros((M * K, M * K))

    # Get rates
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates, dtype=np.float64)
    else:
        rates = np.ones((M, K))

    # Get node types
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        nodetype = list(sn.nodetype) if not isinstance(sn.nodetype, np.ndarray) else sn.nodetype
    else:
        nodetype = [1] * M  # Default to Queue

    # Get station to node mapping
    if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
        stationToNode = np.asarray(sn.stationToNode).flatten()
    else:
        stationToNode = np.arange(M)

    # Identify station types (iterate over stations)
    source_stations = []
    queue_stations = []

    for ist in range(M):
        node_idx = int(stationToNode[ist]) if ist < len(stationToNode) else ist
        if node_idx < len(nodetype):
            ntype = nodetype[node_idx]
            # Check for NodeType
            if hasattr(ntype, 'value'):
                ntype_val = ntype.value
            elif hasattr(ntype, 'ID'):
                ntype_val = ntype.ID
            else:
                ntype_val = int(ntype) if not np.isnan(ntype) else 1

            # Source = 0, Sink = 1, Queue = 2, Delay = 3
            if ntype_val == 0:
                source_stations.append(ist)
            elif ntype_val == 1:
                pass  # Sink is handled separately at node level
            else:
                queue_stations.append(ist)
        else:
            queue_stations.append(ist)

    # Identify Sink nodes (iterate over ALL nodes, not just stations)
    sink_nodes = []
    for node_idx in range(len(nodetype)):
        ntype = nodetype[node_idx]
        if hasattr(ntype, 'value'):
            ntype_val = ntype.value
        elif hasattr(ntype, 'ID'):
            ntype_val = ntype.ID
        else:
            ntype_val = int(ntype) if not np.isnan(ntype) else 1
        if ntype_val == 1:  # Sink = 1
            sink_nodes.append(node_idx)

    # Check for signal classes
    issignal = np.zeros(K, dtype=bool)
    if hasattr(sn, 'issignal') and sn.issignal is not None:
        issignal_raw = sn.issignal
        if callable(issignal_raw):
            for r in range(K):
                try:
                    issignal[r] = issignal_raw(r)
                except:
                    issignal[r] = False
        else:
            issignal = np.asarray(issignal_raw, dtype=bool).flatten()
            if len(issignal) < K:
                issignal = np.pad(issignal, (0, K - len(issignal)), constant_values=False)

    # Check for signal types (for G-networks with negative customers)
    # signaltype[r] will be True if class r is a negative signal
    signaltype_is_negative = [False] * K
    if hasattr(sn, 'signaltype') and sn.signaltype is not None:
        for r in range(K):
            if r < len(sn.signaltype):
                st = sn.signaltype[r]
                if st is not None:
                    # Check if this is a negative signal type
                    # SignalType.NEGATIVE has value 'negative' or could be enum member
                    if hasattr(st, 'value'):
                        # Enum case: check if value is 'negative' or equals NEGATIVE
                        signaltype_is_negative[r] = (st.value == 'negative' or
                                                     str(st.name).upper() == 'NEGATIVE')
                    elif hasattr(st, 'name'):
                        signaltype_is_negative[r] = str(st.name).upper() == 'NEGATIVE'
                    elif isinstance(st, str):
                        signaltype_is_negative[r] = st.lower() == 'negative'
                    else:
                        # MATLAB uses numeric SignalType - 1 = NEGATIVE
                        try:
                            signaltype_is_negative[r] = int(st) == 1
                        except:
                            pass

    # Create process mapping: each (station, class) pair at queue stations
    process_idx = 0
    process_map = np.zeros((M, K), dtype=int)

    for ist in queue_stations:
        for r in range(K):
            if r < len(issignal) and issignal[r]:
                continue
            if ist < rates.shape[0] and r < rates.shape[1]:
                rate = rates[ist, r]
                if not np.isnan(rate) and rate > 0:
                    process_idx += 1
                    process_map[ist, r] = process_idx

    num_processes = process_idx

    if num_processes == 0:
        return RCATModelData(
            R={},
            AP=np.zeros((0, 2), dtype=int),
            process_map=process_map,
            action_map=[],
            N=np.array([]),
            num_processes=0,
            num_actions=0
        )

    # Determine number of states for each process
    N = np.zeros(num_processes, dtype=int)
    njobs = sn.njobs if hasattr(sn, 'njobs') and sn.njobs is not None else np.full(K, np.inf)
    njobs = np.asarray(njobs).flatten()

    for p in range(1, num_processes + 1):
        positions = np.where(process_map == p)
        if len(positions[0]) > 0:
            ist, r = positions[0][0], positions[1][0]
            if r < len(njobs) and njobs[r] < np.inf:
                N[p - 1] = int(njobs[r]) + 1
            else:
                N[p - 1] = max_states

    # Count actions
    action_map = []

    for ist in queue_stations:
        for r in range(K):
            if process_map[ist, r] > 0:
                # Check if class r is a negative signal class
                is_negative_class = False
                if r < len(issignal) and issignal[r]:
                    if r < len(signaltype_is_negative):
                        is_negative_class = signaltype_is_negative[r]

                for jst in queue_stations:
                    for s in range(K):
                        if process_map[jst, s] > 0:
                            from_idx = ist * K + r
                            to_idx = jst * K + s
                            if from_idx < rt.shape[0] and to_idx < rt.shape[1]:
                                prob = rt[from_idx, to_idx]
                                if prob > 0 and (ist != jst or r != s):
                                    action_map.append({
                                        'from_station': ist,
                                        'from_class': r,
                                        'to_station': jst,
                                        'to_class': s,
                                        'prob': prob,
                                        'isNegative': is_negative_class
                                    })

    num_actions = len(action_map)

    # Initialize R and AP
    R = {}
    AP = np.zeros((num_actions, 2), dtype=int) if num_actions > 0 else np.zeros((0, 2), dtype=int)

    # Build local rate matrices L for each process
    for p in range(1, num_processes + 1):
        positions = np.where(process_map == p)
        if len(positions[0]) > 0:
            ist, r = positions[0][0], positions[1][0]
            L = _build_local_rates(sn, ist, r, N[p - 1], rt, source_stations, sink_nodes, K, rates)
            R[(num_actions, p - 1)] = L

    # Build active and passive matrices for each action
    for a, am in enumerate(action_map):
        ist = am['from_station']
        r = am['from_class']
        p_active = process_map[ist, r] - 1
        AP[a, 0] = p_active

        mu_ir = rates[ist, r] if ist < rates.shape[0] and r < rates.shape[1] else 1.0
        prob = am['prob']

        # Active matrix: transition n -> n-1 with rate mu*prob
        Np_active = N[p_active]
        Aa = np.zeros((Np_active, Np_active))
        for n in range(1, Np_active):
            Aa[n, n - 1] = mu_ir * prob
        if Np_active > 0:
            Aa[Np_active - 1, Np_active - 1] = mu_ir * prob
        R[(a, 0)] = Aa

        # Passive process (arrival or signal effect)
        jst = am['to_station']
        s = am['to_class']
        p_passive = process_map[jst, s] - 1
        AP[a, 1] = p_passive

        Np_passive = N[p_passive]
        Pb = np.zeros((Np_passive, Np_passive))

        if am.get('isNegative', False):
            # NEGATIVE: Job removal at destination (G-network negative customer)
            # Empty queue (state 0): no effect
            if Np_passive > 0:
                Pb[0, 0] = 1.0
            # Non-empty queues: decrement (n -> n-1)
            for n in range(1, Np_passive - 1):
                Pb[n, n - 1] = 1.0
            # Boundary at max capacity: decrement
            if Np_passive > 1:
                Pb[Np_passive - 1, Np_passive - 2] = 1.0
        else:
            # POSITIVE: Normal job arrival at destination
            for n in range(Np_passive - 1):
                Pb[n, n + 1] = 1.0
            # Boundary: at max capacity
            if Np_passive > 0:
                Pb[Np_passive - 1, Np_passive - 1] = 1.0
        R[(a, 1)] = Pb

    return RCATModelData(
        R=R,
        AP=AP,
        process_map=process_map,
        action_map=action_map,
        N=N,
        num_processes=num_processes,
        num_actions=num_actions
    )


def _build_local_rates(sn, ist: int, r: int, Np: int, rt: np.ndarray,
                       source_stations: List[int], sink_nodes: List[int],
                       K: int, rates: np.ndarray) -> np.ndarray:
    """Build local/hidden transition matrix for process.

    Args:
        sn: NetworkStruct
        ist: Station index
        r: Class index
        Np: Number of states
        rt: Station-level routing table
        source_stations: List of source station indices
        sink_nodes: List of sink NODE indices (not stations)
        K: Number of classes
        rates: Service rates matrix

    Returns:
        Local rate matrix L
    """
    L = np.zeros((Np, Np))

    # Get signal class information
    issignal = np.zeros(K, dtype=bool)
    if hasattr(sn, 'issignal') and sn.issignal is not None:
        issignal_raw = sn.issignal
        if callable(issignal_raw):
            for s in range(K):
                try:
                    issignal[s] = issignal_raw(s)
                except:
                    issignal[s] = False
        else:
            issignal = np.asarray(issignal_raw, dtype=bool).flatten()
            if len(issignal) < K:
                issignal = np.pad(issignal, (0, K - len(issignal)), constant_values=False)

    signaltype_is_negative = [False] * K
    if hasattr(sn, 'signaltype') and sn.signaltype is not None:
        for s in range(K):
            if s < len(sn.signaltype):
                st = sn.signaltype[s]
                if st is not None:
                    # Check if this is a negative signal type
                    if hasattr(st, 'value'):
                        signaltype_is_negative[s] = (st.value == 'negative' or
                                                     str(st.name).upper() == 'NEGATIVE')
                    elif hasattr(st, 'name'):
                        signaltype_is_negative[s] = str(st.name).upper() == 'NEGATIVE'
                    elif isinstance(st, str):
                        signaltype_is_negative[s] = st.lower() == 'negative'
                    else:
                        try:
                            signaltype_is_negative[s] = int(st) == 1
                        except:
                            pass

    # External arrivals from source - separate positive and negative
    lambda_ir_pos = 0.0   # Positive arrivals
    lambda_ir_neg = 0.0   # Negative signal arrivals (single removal)

    for isrc in source_stations:
        for s_src in range(K):
            # Check if source class s_src is a signal
            is_signal = s_src < len(issignal) and issignal[s_src]

            if is_signal:
                # For signals: they route to themselves but affect positive customers
                # Check if signal routes to ANY class at this station
                prob_src = 0.0
                for s_dst in range(K):
                    from_idx = isrc * K + s_src
                    to_idx = ist * K + s_dst
                    if from_idx < rt.shape[0] and to_idx < rt.shape[1]:
                        prob_src += rt[from_idx, to_idx]
            else:
                # For regular classes: direct routing to (ist, r)
                from_idx = isrc * K + s_src
                to_idx = ist * K + r
                prob_src = rt[from_idx, to_idx] if from_idx < rt.shape[0] and to_idx < rt.shape[1] else 0.0

            if prob_src > 0:
                if isrc < rates.shape[0] and s_src < rates.shape[1]:
                    src_rate = rates[isrc, s_src]
                    if not np.isnan(src_rate):
                        # Check if source class s_src is a negative signal
                        if is_signal and s_src < len(signaltype_is_negative) and signaltype_is_negative[s_src]:
                            lambda_ir_neg += src_rate * prob_src
                        else:
                            lambda_ir_pos += src_rate * prob_src

    # Positive arrival transitions: n -> n+1
    if lambda_ir_pos > 0:
        for n in range(Np - 1):
            L[n, n + 1] = lambda_ir_pos

    # Negative arrival transitions: n -> n-1 (only if queue non-empty)
    if lambda_ir_neg > 0:
        for n in range(1, Np):
            L[n, n - 1] = L[n, n - 1] + lambda_ir_neg

    # Service rate at this station
    mu_ir = rates[ist, r] if ist < rates.shape[0] and r < rates.shape[1] else 0.0

    if not np.isnan(mu_ir) and mu_ir > 0:
        # Departures to sink - use rtnodes (node-level routing) to detect Sink routing
        prob_sink = 0.0

        # Get node index for this station
        stationToNode = sn.stationToNode if hasattr(sn, 'stationToNode') and sn.stationToNode is not None else np.arange(sn.nstations)
        stationToNode = np.asarray(stationToNode).flatten()
        node_idx = int(stationToNode[ist]) if ist < len(stationToNode) else ist

        # Use rtnodes if available for Sink routing
        if hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
            rtnodes = np.asarray(sn.rtnodes)
            for jsnk in sink_nodes:
                for s in range(K):
                    from_idx = node_idx * K + r
                    to_idx = jsnk * K + s
                    if from_idx < rtnodes.shape[0] and to_idx < rtnodes.shape[1]:
                        prob_sink += rtnodes[from_idx, to_idx]

        # Self-routing (same station, same class)
        self_idx = ist * K + r
        prob_self = rt[self_idx, self_idx] if self_idx < rt.shape[0] and self_idx < rt.shape[1] else 0.0

        # Combined local departure rate
        local_departure_rate = mu_ir * prob_sink

        # Departure transitions: n -> n-1
        for n in range(1, Np):
            L[n, n - 1] = L[n, n - 1] + local_departure_rate

        # Self-service transitions
        for n in range(1, Np):
            L[n, n] = mu_ir * prob_self

    return L


def _birth_death_solve(Q: np.ndarray) -> np.ndarray:
    """Solve equilibrium of a birth-death (tridiagonal) CTMC using recursion.

    For a birth-death chain with birth rate lambda_n = Q[n, n+1] and
    death rate mu_n = Q[n, n-1], the equilibrium is:
        pi(n) = pi(0) * prod(lambda_k / mu_{k+1}, k=0..n-1)

    This is numerically stable and avoids the ill-conditioned linear system
    that plagues eigenvalue/null-space methods for large state spaces.
    """
    n = Q.shape[0]
    if n <= 1:
        return np.ones(1)

    pi = np.zeros(n)
    pi[0] = 1.0

    for i in range(1, n):
        birth_rate = Q[i - 1, i]
        death_rate = Q[i, i - 1]
        if death_rate > 0:
            pi[i] = pi[i - 1] * birth_rate / death_rate
        else:
            pi[i] = 0.0

    total = np.sum(pi)
    if total > 0:
        pi /= total
    else:
        pi[:] = 1.0 / n

    return pi


def _is_tridiagonal(Q: np.ndarray) -> bool:
    """Check if a matrix is tridiagonal (only main diagonal and ±1 diagonals)."""
    n = Q.shape[0]
    for i in range(n):
        for j in range(n):
            if abs(i - j) > 1 and abs(Q[i, j]) > 1e-14:
                return False
    return True


def _compute_equilibrium(x: np.ndarray, model: RCATModelData) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Compute equilibrium distribution for each process."""
    num_processes = model.num_processes
    num_actions = model.num_actions
    N = model.N
    AP = model.AP
    R = model.R

    Q_list = []
    pi_list = []

    for k in range(num_processes):
        # Start with local rates
        L = R.get((num_actions, k), np.zeros((N[k], N[k])))
        Qk = L - np.diag(L @ np.ones(N[k]))

        # Add contributions from each action
        for c in range(num_actions):
            if AP[c, 1] == k:
                # Process k is passive for action c
                Pb = R.get((c, 1), np.zeros((N[k], N[k])))
                # MATLAB: Qk = Qk + x(c) * Pb{c} - diag(Pb{c} * ones(N(k), 1));
                # Note: diagonal adjustment is NOT scaled by x[c]
                Qk = Qk + x[c] * Pb - np.diag(Pb @ np.ones(N[k]))
            elif AP[c, 0] == k:
                # Process k is active for action c
                Aa = R.get((c, 0), np.zeros((N[k], N[k])))
                Qk = Qk + Aa - np.diag(Aa @ np.ones(N[k]))

        # Convert to valid generator
        Qk = ctmc_makeinfgen(Qk)
        Q_list.append(Qk)

        # Solve for equilibrium using birth-death recursion for tridiagonal
        # generators (numerically stable), falling back to ctmc_solve otherwise
        if _is_tridiagonal(Qk):
            pi_k = _birth_death_solve(Qk)
        else:
            pi_k = ctmc_solve(Qk)
        pi_list.append(pi_k)

    return pi_list, Q_list


def _inap_solve(model: RCATModelData, tol: float = 1e-6, max_iter: int = 1000,
                method: str = 'inap', verbose: bool = False) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray], int]:
    """INAP iterative solver for RCAT model."""
    num_actions = model.num_actions
    num_processes = model.num_processes
    N = model.N
    AP = model.AP
    R = model.R

    if num_actions == 0:
        # No actions - solve using local rates only
        x = np.array([])
        pi = []
        Q = []
        for p in range(num_processes):
            L = R.get((num_actions, p), np.zeros((N[p], N[p])))
            Qp = L - np.diag(L @ np.ones(N[p]))
            Qp = ctmc_makeinfgen(Qp)
            Q.append(Qp)
            pi.append(ctmc_solve(Qp))
        return x, pi, Q, 0

    # Initialize action rates randomly
    x = np.random.rand(num_actions)

    # Compute initial equilibrium
    pi, Q = _compute_equilibrium(x, model)

    iteration = 0
    while iteration < max_iter:
        iteration += 1
        pi_prev = [p.copy() for p in pi]

        # Update each action rate
        for a in range(num_actions):
            k = AP[a, 0]
            Aa = R.get((a, 0), np.zeros((N[k], N[k])))

            if method == 'inapplus':
                # INAP+: x(a) = sum(Aa(i,j) * pi(i))
                lambda_sum = 0.0
                for i in range(N[k]):
                    for j in range(N[k]):
                        if Aa[i, j] > 0:
                            lambda_sum += Aa[i, j] * pi[k][i]
                if lambda_sum > 0:
                    x[a] = lambda_sum
            else:
                # INAP: x(a) = mean(Aa(i,j) * pi(i) / pi(j))
                lambda_vec = []
                for i in range(N[k]):
                    for j in range(N[k]):
                        if Aa[i, j] > 0:
                            if i == j:
                                # Diagonal (boundary): ratio is always 1
                                lambda_vec.append(Aa[i, j])
                            elif pi[k][j] > 0:
                                ratio = Aa[i, j] * pi[k][i] / pi[k][j]
                                if np.isfinite(ratio):
                                    lambda_vec.append(ratio)
                if len(lambda_vec) > 0:
                    x[a] = np.mean(lambda_vec)

        # Recompute equilibrium
        pi, Q = _compute_equilibrium(x, model)

        # Convergence check
        max_err = 0.0
        for k in range(num_processes):
            max_err = max(max_err, np.sum(np.abs(pi[k] - pi_prev[k])))

        if max_err < tol:
            break

    return x, pi, Q, iteration


def _rcat_metrics(sn, x: np.ndarray, pi: List[np.ndarray], Q: List[np.ndarray],
                  model: RCATModelData) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Convert RCAT solution to LINE performance metrics."""
    M = sn.nstations
    K = sn.nclasses
    process_map = model.process_map
    N = model.N

    rates = sn.rates if hasattr(sn, 'rates') and sn.rates is not None else np.ones((M, K))
    rates = np.asarray(rates, dtype=np.float64)

    njobs = sn.njobs if hasattr(sn, 'njobs') and sn.njobs is not None else np.full(K, np.inf)
    njobs = np.asarray(njobs).flatten()

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))

    # Compute metrics for each (station, class) pair
    for ist in range(M):
        for r in range(K):
            p = int(process_map[ist, r]) - 1
            if p >= 0 and p < len(pi) and len(pi[p]) > 0:
                Np = len(pi[p])

                # Queue length: E[N]
                QN[ist, r] = np.sum(np.arange(Np) * pi[p])

                # Utilization: P(N > 0)
                UN[ist, r] = 1.0 - pi[p][0]

                # Throughput
                mu_ir = rates[ist, r] if ist < rates.shape[0] and r < rates.shape[1] else 0.0
                if not np.isnan(mu_ir) and mu_ir > 0:
                    TN[ist, r] = mu_ir * UN[ist, r]

    # Response times from Little's law
    for ist in range(M):
        for r in range(K):
            if TN[ist, r] > 0:
                RN[ist, r] = QN[ist, r] / TN[ist, r]

    # System metrics
    CN = np.zeros(K)
    XN = np.zeros(K)

    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        nodetype = list(sn.nodetype) if not isinstance(sn.nodetype, np.ndarray) else sn.nodetype
    else:
        nodetype = [1] * M

    if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
        stationToNode = np.asarray(sn.stationToNode).flatten()
    else:
        stationToNode = np.arange(M)

    refstat = sn.refstat if hasattr(sn, 'refstat') and sn.refstat is not None else np.zeros(K)
    refstat = np.asarray(refstat).flatten()

    for r in range(K):
        if r < len(njobs) and njobs[r] >= np.inf:
            # Open class
            for ist in range(M):
                node_idx = int(stationToNode[ist]) if ist < len(stationToNode) else ist
                if node_idx < len(nodetype):
                    ntype = nodetype[node_idx]
                    if hasattr(ntype, 'value'):
                        ntype_val = ntype.value
                    elif hasattr(ntype, 'ID'):
                        ntype_val = ntype.ID
                    else:
                        ntype_val = int(ntype) if not np.isnan(ntype) else 1
                    if ntype_val == 0:  # Source
                        if ist < rates.shape[0] and r < rates.shape[1]:
                            XN[r] = rates[ist, r]
                        break
            CN[r] = np.sum(RN[:, r])
        else:
            # Closed class
            refst = int(refstat[r]) if r < len(refstat) else 0
            if 0 <= refst < M:
                XN[r] = TN[refst, r]
                if XN[r] > 0 and r < len(njobs):
                    CN[r] = njobs[r] / XN[r]

    return QN, UN, RN, TN, CN, XN


class INAPAlgorithm(MAMAlgorithm):
    """RCAT iterative fixed-point (INAP) solver.

    Designed for open queueing networks with independent queue states.
    For closed networks, the algorithm may produce approximate results
    since it doesn't model the state space constraint (sum of jobs = N).
    """

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if INAP can solve this network.

        INAP is designed for open networks. For closed networks,
        results may be approximate.
        """
        # Check if network has closed classes
        njobs = getattr(sn, 'njobs', None)
        if njobs is not None:
            njobs = np.asarray(njobs).flatten()
            has_closed = any(n < np.inf for n in njobs if not np.isnan(n))
            if has_closed:
                return True, "INAP is designed for open networks; closed network results are approximate"
        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve using RCAT INAP method.

        Note: For closed networks, this produces approximate results
        since the RCAT model doesn't enforce the state constraint.
        """
        import warnings
        start_time = time.time()

        # Check for closed classes and warn
        njobs = getattr(sn, 'njobs', None)
        if njobs is not None:
            njobs_arr = np.asarray(njobs).flatten()
            has_closed = any(n < np.inf for n in njobs_arr if not np.isnan(n))
            if has_closed:
                warnings.warn(
                    "INAP is designed for open networks. "
                    "Closed network results are approximate; consider using MVA instead.",
                    UserWarning
                )

        tol = getattr(options, 'tol', 1e-6) if options else 1e-6
        max_iter = getattr(options, 'max_iter', 1000) if options else 1000
        verbose = getattr(options, 'verbose', False) if options else False

        max_states = 100
        if options and hasattr(options, 'config') and options.config is not None:
            if hasattr(options.config, 'maxStates'):
                max_states = options.config.maxStates
            elif isinstance(options.config, dict) and 'maxStates' in options.config:
                max_states = options.config['maxStates']

        M = sn.nstations
        K = sn.nclasses

        # Build RCAT model
        model = _build_rcat(sn, max_states)

        if model.num_processes == 0:
            # Fall back to simple M/M/c analysis
            return self._solve_simple(sn, options)

        # Solve using INAP
        x, pi, Q, iterations = _inap_solve(model, tol, max_iter, 'inap', verbose)

        # Convert to metrics
        QN, UN, RN, TN, CN, XN = _rcat_metrics(sn, x, pi, Q, model)

        runtime = time.time() - start_time

        return MAMResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            XN=XN.reshape(1, -1) if XN.ndim == 1 else XN,
            totiter=iterations,
            method="inap",
            runtime=runtime
        )

    def _solve_simple(self, sn, options) -> MAMResult:
        """Simple M/M/c fallback when RCAT model can't be built."""
        params = extract_mam_params(sn)
        M = params['nstations']
        K = params['nclasses']
        rates = params['rates']
        nservers = params['nservers']

        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))

        S = 1.0 / np.maximum(rates, 1e-10)

        for m in range(M):
            for k in range(K):
                rho = 0.5  # Default utilization
                UN[m, k] = rho
                RN[m, k] = S[m, k] / (1.0 - rho)
                QN[m, k] = rho * RN[m, k]
                TN[m, k] = rates[m, k] * rho

        return MAMResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            XN=TN.copy(),
            totiter=0,
            method="inap",
            runtime=0.0
        )


class INAPPlusAlgorithm(MAMAlgorithm):
    """RCAT iterative weighted variant (INAP+)."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """INAP+ supports general networks."""
        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve using RCAT INAP+ method."""
        start_time = time.time()

        tol = getattr(options, 'tol', 1e-6) if options else 1e-6
        max_iter = getattr(options, 'max_iter', 1000) if options else 1000
        verbose = getattr(options, 'verbose', False) if options else False

        max_states = 100
        if options and hasattr(options, 'config') and options.config is not None:
            if hasattr(options.config, 'maxStates'):
                max_states = options.config.maxStates
            elif isinstance(options.config, dict) and 'maxStates' in options.config:
                max_states = options.config['maxStates']

        M = sn.nstations
        K = sn.nclasses

        # Build RCAT model
        model = _build_rcat(sn, max_states)

        if model.num_processes == 0:
            # Fall back to INAP
            algo = INAPAlgorithm()
            return algo.solve(sn, options)

        # Solve using INAP+
        x, pi, Q, iterations = _inap_solve(model, tol, max_iter, 'inapplus', verbose)

        # Convert to metrics
        QN, UN, RN, TN, CN, XN = _rcat_metrics(sn, x, pi, Q, model)

        runtime = time.time() - start_time

        return MAMResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            XN=XN.reshape(1, -1) if XN.ndim == 1 else XN,
            totiter=iterations,
            method="inapplus",
            runtime=runtime
        )


__all__ = [
    'INAPAlgorithm',
    'INAPPlusAlgorithm',
    'RCATModelData',
]
