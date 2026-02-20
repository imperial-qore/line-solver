"""
SN Getter Functions for Parameter Extraction.

Native Python implementations for extracting parameters from
network structures including arrival rates, throughputs, and
product-form chain parameters.

Key functions:
    sn_get_arvr_from_tput: Compute arrival rates from throughputs
    sn_get_node_arvr_from_tput: Compute node arrival rates from throughputs
    sn_get_node_tput_from_tput: Compute node throughputs from station throughputs
    sn_get_product_form_chain_params: Extract chain-aggregated parameters

References:
    Original MATLAB: matlab/src/api/sn/sn_get_*.m
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass

from .network_struct import NetworkStruct, NodeType


@dataclass
class ChainParams:
    """Chain-aggregated product-form parameters."""
    lambda_vec: np.ndarray  # Chain arrival rates
    D: np.ndarray  # Chain service demands at queuing stations
    N: np.ndarray  # Chain populations
    Z: np.ndarray  # Chain think times
    mu: np.ndarray  # Load-dependent service capacity scaling
    S: np.ndarray  # Number of servers at queuing stations
    V: np.ndarray  # Chain visit ratios


def sn_get_arvr_from_tput(sn: NetworkStruct, TN: np.ndarray,
                          TH: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute average arrival rates at stations from throughputs.

    Calculates the average arrival rate at each station in steady-state
    from the station throughputs and routing matrix.

    Args:
        sn: Network structure
        TN: Average throughputs at stations (M x R)
        TH: Throughput handles (optional)

    Returns:
        AN: Average arrival rates at stations (M x R)

    References:
        Original MATLAB: matlab/src/api/sn/sn_get_arvr_from_tput.m
    """
    M = sn.nstations
    R = sn.nclasses

    if TN is None or len(TN) == 0:
        return np.array([])

    TN = np.atleast_2d(np.asarray(TN, dtype=float))
    AN = np.zeros((M, R))

    # Build mapping from stateful nodes to their position in rt matrix
    stateful_nodes = np.where(sn.isstateful)[0]
    n_stateful = len(stateful_nodes)

    # Build throughput vector for all stateful nodes (stations have TN, others need computation)
    TN_stateful = np.zeros((n_stateful, R))
    for sf, ind in enumerate(stateful_nodes):
        ist = sn.nodeToStation[ind]

        # Check if this is a Cache node - needs special handling regardless of station status
        # In some implementations (like Python), Cache may be treated as a station
        if sn.nodetype[ind] == NodeType.CACHE:
            # For Cache nodes, compute hit/miss class throughputs
            # from the reference station throughput and hit/miss probabilities
            if not hasattr(sn, 'nodeparam') or sn.nodeparam is None:
                # Fall back to station throughput if available
                if ist >= 0 and ist < TN.shape[0]:
                    TN_stateful[sf, :] = TN[ist, :]
                continue

            # nodeparam may be a dict or list
            if isinstance(sn.nodeparam, dict):
                nodeparam = sn.nodeparam.get(ind, None)
            else:
                nodeparam = sn.nodeparam[ind] if ind < len(sn.nodeparam) else None

            if nodeparam is None:
                # Fall back to station throughput if available
                if ist >= 0 and ist < TN.shape[0]:
                    TN_stateful[sf, :] = TN[ist, :]
                continue

            hitclass = getattr(nodeparam, 'hitclass', None)
            missclass = getattr(nodeparam, 'missclass', None)
            if hitclass is None or missclass is None:
                # Fall back to station throughput if available
                if ist >= 0 and ist < TN.shape[0]:
                    TN_stateful[sf, :] = TN[ist, :]
                continue

            # Get actual hit/miss probabilities if available
            actualhitprob = getattr(nodeparam, 'actualhitprob', None)
            actualmissprob = getattr(nodeparam, 'actualmissprob', None)
            if actualhitprob is None or actualmissprob is None:
                # Actual probabilities not yet computed - fall back to station throughput
                if ist >= 0 and ist < TN.shape[0]:
                    TN_stateful[sf, :] = TN[ist, :]
                continue

            hitclass = np.atleast_1d(hitclass)
            missclass = np.atleast_1d(missclass)
            actualhitprob = np.atleast_1d(actualhitprob)
            actualmissprob = np.atleast_1d(actualmissprob)

            # Find the chain and reference station for this cache
            if hasattr(sn, 'nchains') and hasattr(sn, 'inchain') and hasattr(sn, 'refstat'):
                for c in range(sn.nchains):
                    inchain = sn.inchain[c] if sn.inchain is not None else []
                    refstat = int(sn.refstat[c]) if sn.refstat is not None else 0
                    if refstat < TN.shape[0]:
                        totalTput = np.sum(TN[refstat, list(inchain)]) if len(inchain) > 0 else 0

                        # Set throughput for hit/miss classes
                        # Note: hitclass/missclass contain 0-indexed class indices
                        # Values > 0 are valid class indices; 0 or negative means no class
                        for origClass in range(len(hitclass)):
                            hc = int(hitclass[origClass])
                            mc = int(missclass[origClass])
                            if hc > 0 and hc < R:
                                TN_stateful[sf, hc] = totalTput * actualhitprob[origClass]
                            if mc > 0 and mc < R:
                                TN_stateful[sf, mc] = totalTput * actualmissprob[origClass]
        elif ist >= 0 and ist < TN.shape[0]:
            # This stateful node is a station - use station throughput
            TN_stateful[sf, :] = TN[ist, :]

    # Propagate throughputs to non-station, non-Cache stateful nodes
    # (e.g., Router) using already-known station throughputs and the rt matrix
    for sf in range(n_stateful):
        ind = stateful_nodes[sf]
        ist = sn.nodeToStation[ind]
        if ist < 0 and sn.nodetype[ind] != NodeType.CACHE:
            # Compute throughput from upstream stateful nodes
            for sf_src in range(n_stateful):
                for k in range(R):
                    for r in range(R):
                        from_idx = sf_src * R + r
                        to_idx = sf * R + k
                        if from_idx < sn.rt.shape[0] and to_idx < sn.rt.shape[1]:
                            TN_stateful[sf, k] += TN_stateful[sf_src, r] * sn.rt[from_idx, to_idx]

    # Compute arrival rates using stateful node throughputs and rt matrix
    for ist in range(M):
        ind_ist = sn.stationToNode[ist]
        if sn.nodetype[ind_ist] == NodeType.SOURCE:
            AN[ist, :] = 0
        else:
            # Find position of this station in stateful node list
            sf_ist_arr = np.where(stateful_nodes == ind_ist)[0]
            if len(sf_ist_arr) == 0:
                continue
            sf_ist = sf_ist_arr[0]

            for sf_jst in range(n_stateful):
                for k in range(R):
                    for r in range(R):
                        from_idx = sf_jst * R + r
                        to_idx = sf_ist * R + k
                        if from_idx < sn.rt.shape[0] and to_idx < sn.rt.shape[1]:
                            AN[ist, k] += TN_stateful[sf_jst, r] * sn.rt[from_idx, to_idx]

    # Fork-join special handling: delegate to node-level arrival rate computation
    # MATLAB's sn_get_arvr_from_tput lines 117-122
    if hasattr(sn, 'fj') and sn.fj is not None and np.any(sn.fj):
        ANn = sn_get_node_arvr_from_tput(sn, TN, TH, AN)
        for ist in range(M):
            ind = sn.stationToNode[ist]
            if ind >= 0 and ind < ANn.shape[0]:
                AN[ist, :] = ANn[ind, :]

    return AN


def sn_get_node_arvr_from_tput(sn: NetworkStruct, TN: np.ndarray,
                                TH: Optional[np.ndarray] = None,
                                AN: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute node arrival rates from station throughputs.

    This function handles:
    - Station nodes: Uses station arrival rates directly
    - Cache nodes: Only requesting classes arrive (not hit/miss classes)
    - Non-station nodes (ClassSwitch, Sink): Uses nodevisits-based computation

    Args:
        sn: Network structure
        TN: Station throughputs (M x R)
        TH: Throughput handles (optional)
        AN: Station arrival rates (optional, computed if not provided)

    Returns:
        ANn: Node arrival rates (I x R)

    References:
        Original MATLAB: matlab/src/api/sn/sn_get_node_arvr_from_tput.m
    """
    I = sn.nnodes
    M = sn.nstations
    C = sn.nchains
    R = sn.nclasses

    if AN is None:
        AN = sn_get_arvr_from_tput(sn, TN, TH)

    ANn = np.zeros((I, R))

    # Check if we have valid throughput data
    # Note: TH is optional and not used directly in the computation
    # We only need TN and AN (which may be computed from TN)
    if TN is None or TN.size == 0:
        return ANn

    # First, copy station arrival rates to station nodes
    for ist in range(M):
        ind = sn.stationToNode[ist]
        if ind >= 0 and ind < I:
            ANn[ind, :] = AN[ist, :]

    # Process non-station nodes
    for ind in range(I):
        if sn.nodetype is None or ind >= len(sn.nodetype):
            continue

        node_type = sn.nodetype[ind]

        # Skip Source nodes
        if node_type == NodeType.SOURCE:
            continue

        for c in range(C):
            if c not in sn.inchain:
                continue

            inchain = sn.inchain[c].flatten().astype(int)
            refstat_idx = int(sn.refstat[c]) if c < len(sn.refstat) else 0

            for r in inchain:
                if r >= R:
                    continue

                if node_type == NodeType.CACHE:
                    # For cache nodes, only the requesting class arrives
                    # Hit/miss classes don't arrive - they leave
                    hitclass = []
                    missclass = []

                    # Get hit/miss classes from nodeparam
                    if sn.nodeparam is not None and ind in sn.nodeparam:
                        node_param = sn.nodeparam[ind]
                        if hasattr(node_param, 'hitclass'):
                            hitclass = np.atleast_1d(node_param.hitclass).flatten()
                        if hasattr(node_param, 'missclass'):
                            missclass = np.atleast_1d(node_param.missclass).flatten()

                    # Check if this class is a hit or miss class
                    is_hit_or_miss = (r in hitclass) or (r in missclass)

                    if not is_hit_or_miss:
                        # This is a requesting class
                        # Use nodevisits-based calculation like MATLAB:
                        # ANn(ind, r) = (nodevisits{c}(ind,r) / sum(nodevisits{c}(refnode,:))) * sum(TN(refstat,:))
                        # Note: For Python where Cache is a Station, we must use nodevisits (not visits)
                        # for normalization because station-level visits are incorrect when routing
                        # goes through non-station nodes like ClassSwitch.
                        if c in sn.nodevisits:
                            nodevisits_c = sn.nodevisits[c]

                            # Get reference node index (node corresponding to refstat)
                            refnode_idx = int(sn.stationToNode[refstat_idx]) if refstat_idx < len(sn.stationToNode) else refstat_idx

                            if ind < nodevisits_c.shape[0] and r < nodevisits_c.shape[1]:
                                nodevisit_val = nodevisits_c[ind, r]

                                # Sum of nodevisits at refnode for all classes in chain
                                sum_nodevisits_refnode = 0.0
                                if refnode_idx < nodevisits_c.shape[0]:
                                    for rc in inchain:
                                        if rc < nodevisits_c.shape[1]:
                                            sum_nodevisits_refnode += nodevisits_c[refnode_idx, rc]

                                # Total throughput at refstat for all classes in chain
                                total_tput_refstat = 0.0
                                if refstat_idx < TN.shape[0]:
                                    for rc in inchain:
                                        if rc < TN.shape[1]:
                                            val = TN[refstat_idx, rc]
                                            if not np.isnan(val):
                                                total_tput_refstat += val

                                if sum_nodevisits_refnode > 0 and total_tput_refstat > 0:
                                    ANn[ind, r] = (nodevisit_val / sum_nodevisits_refnode) * total_tput_refstat
                                elif nodevisit_val > 0:
                                    # Fallback: if refstat has no throughput, try to get from any station
                                    # But skip stations with very high rates (1e6+) which are instant-service
                                    for ist in range(M):
                                        if r < TN.shape[1]:
                                            val = TN[ist, r]
                                            if val > 0 and not np.isnan(val) and val < 1e6:
                                                ANn[ind, r] = val
                                                break
                    # Hit/miss classes have 0 arrival rate at cache (they only depart)

                elif node_type == NodeType.CLASSSWITCH:
                    # For ClassSwitch nodes that follow Cache nodes
                    # Hit/miss classes arrive from the Cache
                    # Check if there's a Cache node that connects to this ClassSwitch
                    for cache_ind in range(I):
                        if cache_ind >= len(sn.nodetype):
                            continue
                        if sn.nodetype[cache_ind] != NodeType.CACHE:
                            continue
                        if sn.nodeparam is None or cache_ind not in sn.nodeparam:
                            continue

                        cache_param = sn.nodeparam[cache_ind]
                        hitclass = np.atleast_1d(cache_param.hitclass).flatten() if hasattr(cache_param, 'hitclass') else []
                        missclass = np.atleast_1d(cache_param.missclass).flatten() if hasattr(cache_param, 'missclass') else []
                        actual_hit_prob = np.atleast_1d(cache_param.actualhitprob).flatten() if hasattr(cache_param, 'actualhitprob') and cache_param.actualhitprob is not None else None
                        actual_miss_prob = np.atleast_1d(cache_param.actualmissprob).flatten() if hasattr(cache_param, 'actualmissprob') and cache_param.actualmissprob is not None else None

                        # Check if r is a hit or miss class for this cache
                        for orig_class in range(len(hitclass)):
                            if actual_hit_prob is None or orig_class >= len(actual_hit_prob):
                                continue

                            # Get throughput of requesting class
                            req_tput = 0.0
                            for ist in range(M):
                                if orig_class < TN.shape[1]:
                                    val = TN[ist, orig_class]
                                    if val > 0 and not np.isnan(val):
                                        req_tput = val
                                        break

                            if req_tput > 0:
                                if hitclass[orig_class] == r:
                                    # r is a hit class - arrival rate at ClassSwitch = hit throughput
                                    ANn[ind, r] = req_tput * actual_hit_prob[orig_class]
                                elif orig_class < len(missclass) and missclass[orig_class] == r:
                                    # r is a miss class - arrival rate at ClassSwitch = miss throughput
                                    miss_prob = actual_miss_prob[orig_class] if actual_miss_prob is not None and orig_class < len(actual_miss_prob) else (1 - actual_hit_prob[orig_class])
                                    ANn[ind, r] = req_tput * miss_prob

                elif node_type == NodeType.SINK:
                    # For Sink nodes, arrival rate equals total throughput of class r from upstream
                    # This handles Cache hit/miss classes correctly
                    # Check if class r is a hit/miss class from a Cache node
                    is_hit_miss_class = False
                    for cache_ind in range(I):
                        if cache_ind >= len(sn.nodetype):
                            continue
                        if sn.nodetype[cache_ind] != NodeType.CACHE:
                            continue
                        if sn.nodeparam is None or cache_ind not in sn.nodeparam:
                            continue

                        cache_param = sn.nodeparam[cache_ind]
                        hitclass = np.atleast_1d(cache_param.hitclass).flatten() if hasattr(cache_param, 'hitclass') else []
                        missclass = np.atleast_1d(cache_param.missclass).flatten() if hasattr(cache_param, 'missclass') else []
                        actual_hit_prob = np.atleast_1d(cache_param.actualhitprob).flatten() if hasattr(cache_param, 'actualhitprob') and cache_param.actualhitprob is not None else None
                        actual_miss_prob = np.atleast_1d(cache_param.actualmissprob).flatten() if hasattr(cache_param, 'actualmissprob') and cache_param.actualmissprob is not None else None

                        for orig_class in range(len(hitclass)):
                            if actual_hit_prob is None or orig_class >= len(actual_hit_prob):
                                continue

                            # Get throughput of requesting class
                            req_tput = 0.0
                            for ist in range(M):
                                if orig_class < TN.shape[1]:
                                    val = TN[ist, orig_class]
                                    if val > 0 and not np.isnan(val):
                                        req_tput = val
                                        break

                            if req_tput > 0:
                                if hitclass[orig_class] == r:
                                    # r is a hit class - arrival rate at Sink = hit throughput
                                    ANn[ind, r] = req_tput * actual_hit_prob[orig_class]
                                    is_hit_miss_class = True
                                    break
                                elif orig_class < len(missclass) and missclass[orig_class] == r:
                                    # r is a miss class - arrival rate at Sink = miss throughput
                                    miss_prob = actual_miss_prob[orig_class] if actual_miss_prob is not None and orig_class < len(actual_miss_prob) else (1 - actual_hit_prob[orig_class])
                                    ANn[ind, r] = req_tput * miss_prob
                                    is_hit_miss_class = True
                                    break
                        if is_hit_miss_class:
                            break

                    # If not a hit/miss class, use nodevisits-based computation
                    if not is_hit_miss_class:
                        if c in sn.nodevisits:
                            nodevisits_c = sn.nodevisits[c]
                            if ind < nodevisits_c.shape[0] and r < nodevisits_c.shape[1]:
                                nodevisit_val = nodevisits_c[ind, r]
                                refnode_idx = int(sn.stationToNode[refstat_idx]) if refstat_idx < len(sn.stationToNode) else refstat_idx
                                nodevisit_sum = 0.0
                                for s in inchain:
                                    if refnode_idx < nodevisits_c.shape[0] and s < nodevisits_c.shape[1]:
                                        nodevisit_sum += nodevisits_c[refnode_idx, s]
                                tput_sum = 0.0
                                for s in inchain:
                                    if refstat_idx < TN.shape[0] and s < TN.shape[1]:
                                        tput_sum += TN[refstat_idx, s]
                                if nodevisit_sum > 0:
                                    ANn[ind, r] = (nodevisit_val / nodevisit_sum) * tput_sum

                else:
                    # For other non-station nodes
                    # Check if this node has a station mapping
                    node_to_station = sn.nodeToStation[ind] if ind < len(sn.nodeToStation) else -1

                    if node_to_station < 0 or np.isnan(node_to_station):
                        # This is a non-station node (e.g., Router, Sink, ClassSwitch)
                        # Check if class r is a hit/miss class from a Cache node.
                        # If so, use actual cache hit/miss probabilities for accurate
                        # per-class arrival rates (nodevisits don't include cache splits).
                        is_cache_hit_miss = False
                        for cache_ind in range(I):
                            if cache_ind >= len(sn.nodetype):
                                continue
                            if sn.nodetype[cache_ind] != NodeType.CACHE:
                                continue
                            if sn.nodeparam is None or cache_ind not in sn.nodeparam:
                                continue

                            cache_param = sn.nodeparam[cache_ind]
                            hitclass = np.atleast_1d(cache_param.hitclass).flatten() if hasattr(cache_param, 'hitclass') else []
                            missclass = np.atleast_1d(cache_param.missclass).flatten() if hasattr(cache_param, 'missclass') else []
                            actual_hit_prob = np.atleast_1d(cache_param.actualhitprob).flatten() if hasattr(cache_param, 'actualhitprob') and cache_param.actualhitprob is not None else None
                            actual_miss_prob = np.atleast_1d(cache_param.actualmissprob).flatten() if hasattr(cache_param, 'actualmissprob') and cache_param.actualmissprob is not None else None

                            if actual_hit_prob is None:
                                continue

                            for orig_class in range(len(hitclass)):
                                if orig_class >= len(actual_hit_prob):
                                    continue

                                # Get throughput of requesting class
                                req_tput = 0.0
                                for ist_s in range(M):
                                    if orig_class < TN.shape[1]:
                                        val = TN[ist_s, orig_class]
                                        if val > 0 and not np.isnan(val):
                                            req_tput = val
                                            break

                                if req_tput > 0:
                                    if hitclass[orig_class] == r:
                                        ANn[ind, r] = req_tput * actual_hit_prob[orig_class]
                                        is_cache_hit_miss = True
                                        break
                                    elif orig_class < len(missclass) and missclass[orig_class] == r:
                                        miss_prob = actual_miss_prob[orig_class] if actual_miss_prob is not None and orig_class < len(actual_miss_prob) else (1 - actual_hit_prob[orig_class])
                                        ANn[ind, r] = req_tput * miss_prob
                                        is_cache_hit_miss = True
                                        break
                            if is_cache_hit_miss:
                                break

                        # If not a cache hit/miss class, use nodevisits-based computation
                        if not is_cache_hit_miss:
                            if c in sn.nodevisits:
                                nodevisits_c = sn.nodevisits[c]

                                if ind < nodevisits_c.shape[0] and r < nodevisits_c.shape[1]:
                                    nodevisit_val = nodevisits_c[ind, r]

                                    # Get reference node index
                                    refnode_idx = int(sn.stationToNode[refstat_idx]) if refstat_idx < len(sn.stationToNode) else refstat_idx

                                    # Sum nodevisits at refnode for all classes in chain
                                    nodevisit_sum = 0.0
                                    for s in inchain:
                                        if refnode_idx < nodevisits_c.shape[0] and s < nodevisits_c.shape[1]:
                                            nodevisit_sum += nodevisits_c[refnode_idx, s]

                                    # Sum throughput at refstat for all classes in chain
                                    tput_sum = 0.0
                                    for s in inchain:
                                        if refstat_idx < TN.shape[0] and s < TN.shape[1]:
                                            tput_sum += TN[refstat_idx, s]

                                    if nodevisit_sum > 0:
                                        ANn[ind, r] = (nodevisit_val / nodevisit_sum) * tput_sum

    # Replace NaN with 0
    ANn = np.nan_to_num(ANn, nan=0.0)

    return ANn


def sn_get_node_tput_from_tput(sn: NetworkStruct, TN: np.ndarray,
                                TH: Optional[np.ndarray] = None,
                                ANn: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute node throughputs from station throughputs.

    This function handles:
    - Station nodes: Uses station throughputs directly
    - Cache nodes: Uses actual hit/miss probabilities if available
    - Non-station nodes: Uses routing matrix (rtnodes) for computation

    Args:
        sn: Network structure
        TN: Station throughputs (M x R)
        TH: Throughput handles (optional)
        ANn: Node arrival rates (optional, computed if not provided)

    Returns:
        TNn: Node throughputs (I x R)

    References:
        Original MATLAB: matlab/src/api/sn/sn_get_node_tput_from_tput.m
    """
    I = sn.nnodes
    M = sn.nstations
    C = sn.nchains
    R = sn.nclasses

    TN = np.atleast_2d(np.asarray(TN, dtype=float))

    if ANn is None:
        ANn = sn_get_node_arvr_from_tput(sn, TN, TH)

    TNn = np.zeros((I, R))

    # Check if we have valid throughput data
    if TH is None or TN is None or TN.size == 0:
        return TNn

    # First pass: Process Cache nodes with actual hit/miss probabilities
    for ind in range(I):
        if sn.nodetype is None or ind >= len(sn.nodetype):
            continue

        node_type = sn.nodetype[ind]

        if node_type == NodeType.CACHE:
            # Get hit/miss class indices from nodeparam
            hitclass = np.array([])
            missclass = np.array([])
            actual_hit_prob = None
            actual_miss_prob = None

            if sn.nodeparam is not None and ind in sn.nodeparam:
                node_param = sn.nodeparam[ind]
                if hasattr(node_param, 'hitclass'):
                    hitclass = np.atleast_1d(node_param.hitclass).flatten()
                if hasattr(node_param, 'missclass'):
                    missclass = np.atleast_1d(node_param.missclass).flatten()
                if hasattr(node_param, 'actualhitprob') and node_param.actualhitprob is not None:
                    actual_hit_prob = np.atleast_1d(node_param.actualhitprob).flatten()
                if hasattr(node_param, 'actualmissprob') and node_param.actualmissprob is not None:
                    actual_miss_prob = np.atleast_1d(node_param.actualmissprob).flatten()

            # For Cache nodes with class switching, each class may be in a different chain
            # Find the total throughput from the chain with the requesting class (not hit/miss classes)
            # The requesting class is the one that has a non-negative hit/miss class mapping

            # First, get the station index for this cache node to check if TN already has
            # hit/miss class throughputs (e.g., from SSA which computes them directly)
            cache_ist = int(sn.nodeToStation[ind]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and ind < len(sn.nodeToStation) else -1

            for orig_class in range(len(hitclass)):
                h = int(hitclass[orig_class])
                m = int(missclass[orig_class])

                if h >= 0 and m >= 0:
                    # Check if TN already has hit/miss class throughputs at this cache station
                    # This happens when SSA computes them directly during simulation
                    if cache_ist >= 0 and cache_ist < TN.shape[0]:
                        t_hit_existing = TN[cache_ist, h] if h < TN.shape[1] else 0.0
                        t_miss_existing = TN[cache_ist, m] if m < TN.shape[1] else 0.0

                        if t_hit_existing > 0 or t_miss_existing > 0:
                            # Use existing hit/miss throughputs from TN directly
                            if h < R:
                                TNn[ind, h] = t_hit_existing
                            if m < R:
                                TNn[ind, m] = t_miss_existing
                            continue  # Skip recomputing for this orig_class

                    # orig_class is a requesting class with hit/miss classes
                    # Find the throughput for this requesting class
                    # Look at refstat (typically station 0)
                    refstat_idx = 0
                    if orig_class < C and orig_class in sn.inchain:
                        inchain_orig = sn.inchain[orig_class].flatten().astype(int)
                        if orig_class < len(sn.refstat):
                            refstat_idx = int(sn.refstat[orig_class])

                    # Get throughput of the requesting class at refstat
                    req_tput = 0.0
                    if refstat_idx < TN.shape[0] and orig_class < TN.shape[1]:
                        val = TN[refstat_idx, orig_class]
                        if not np.isnan(val):
                            req_tput = val

                    # If no throughput at refstat, try to find it from any station
                    # Skip stations with very high rates (1e6+) which are instant-service
                    if req_tput == 0:
                        for ist in range(TN.shape[0]):
                            if orig_class < TN.shape[1]:
                                val = TN[ist, orig_class]
                                if val > 0 and not np.isnan(val) and val < 1e6:
                                    req_tput = val
                                    break

                    if req_tput > 0 and actual_hit_prob is not None and orig_class < len(actual_hit_prob):
                        # Set hit class throughput at cache
                        if h < R:
                            TNn[ind, h] = req_tput * actual_hit_prob[orig_class]
                        # Set miss class throughput at cache
                        if m < R:
                            miss_prob = actual_miss_prob[orig_class] if actual_miss_prob is not None and orig_class < len(actual_miss_prob) else (1 - actual_hit_prob[orig_class])
                            TNn[ind, m] = req_tput * miss_prob

    # Second pass: Copy station throughputs directly to station nodes
    # Skip Cache nodes - their throughputs were computed in the first pass
    for ist in range(M):
        ind = sn.stationToNode[ist]
        if ind >= 0 and ind < I:
            # Don't overwrite Cache node throughputs - they have special handling
            if sn.nodetype is not None and ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE:
                continue
            TNn[ind, :] = TN[ist, :]

    # Third pass: Compute throughputs for non-station nodes using routing matrix
    for ind in range(I):
        if sn.nodetype is None or ind >= len(sn.nodetype):
            continue

        node_type = sn.nodetype[ind]

        # Skip Source, Sink, Join nodes and station nodes
        if node_type in [NodeType.SOURCE, NodeType.SINK, NodeType.JOIN]:
            continue

        # Check if this is a station node
        node_to_station = sn.nodeToStation[ind] if ind < len(sn.nodeToStation) else -1
        if node_to_station >= 0 and not np.isnan(node_to_station):
            continue  # Already handled above

        for c in range(C):
            if c not in sn.inchain:
                continue

            inchain = sn.inchain[c].flatten().astype(int)

            for r in inchain:
                if r >= R:
                    continue

                # Check if there are any visits for this class at any stateful node
                # (skip this check for ClassSwitch nodes which need to process hit/miss classes)
                if node_type != NodeType.CLASSSWITCH:
                    if c not in sn.visits:
                        continue

                    visits_c = sn.visits[c]
                    any_stateful = np.any(visits_c[:, r] > 0) if r < visits_c.shape[1] else False

                    if not any_stateful:
                        continue

                # For Cache nodes, compute throughput from arrival rate and routing
                # BUT skip if actual hit/miss probs were already used in the first pass
                if node_type == NodeType.CACHE:
                    # Check if first pass already handled this cache
                    has_actual_probs = False
                    if sn.nodeparam is not None and ind in sn.nodeparam:
                        node_param = sn.nodeparam[ind]
                        if hasattr(node_param, 'actualhitprob') and node_param.actualhitprob is not None:
                            actual_probs = np.atleast_1d(node_param.actualhitprob).flatten()
                            has_actual_probs = np.any(actual_probs > 0)

                    if has_actual_probs:
                        continue  # Skip - already computed with actual probs in first pass

                    for s in inchain:
                        if s >= R:
                            continue
                        for jnd in range(I):
                            if ind != jnd:
                                # Use rtnodes for routing probability
                                if sn.rtnodes is not None and sn.rtnodes.size > 0:
                                    from_idx = ind * R + r
                                    to_idx = jnd * R + s
                                    if from_idx < sn.rtnodes.shape[0] and to_idx < sn.rtnodes.shape[1]:
                                        TNn[ind, s] += ANn[ind, r] * sn.rtnodes[from_idx, to_idx]
                elif node_type == NodeType.CLASSSWITCH:
                    # For ClassSwitch nodes that follow Cache nodes
                    # The output class (JobClass) gets throughput = sum of hit + miss arrival rates
                    # Check if there's a Cache node that connects to this ClassSwitch
                    for cache_ind in range(I):
                        if cache_ind >= len(sn.nodetype):
                            continue
                        if sn.nodetype[cache_ind] != NodeType.CACHE:
                            continue
                        if sn.nodeparam is None or cache_ind not in sn.nodeparam:
                            continue

                        cache_param = sn.nodeparam[cache_ind]
                        hitclass = np.atleast_1d(cache_param.hitclass).flatten() if hasattr(cache_param, 'hitclass') else []
                        missclass = np.atleast_1d(cache_param.missclass).flatten() if hasattr(cache_param, 'missclass') else []

                        # Find the requesting class that has r as its hit or miss class
                        for orig_class in range(len(hitclass)):
                            if orig_class < len(hitclass) and hitclass[orig_class] == r:
                                # r is a hit class - output class is orig_class
                                # Throughput of output class = arrival rate of hit + miss classes
                                # which equals the throughput of the requesting class
                                req_tput = 0.0
                                for ist in range(M):
                                    if orig_class < TN.shape[1]:
                                        val = TN[ist, orig_class]
                                        if val > 0 and not np.isnan(val):
                                            req_tput = val
                                            break
                                if req_tput > 0 and orig_class < R:
                                    TNn[ind, orig_class] = req_tput
                            elif orig_class < len(missclass) and missclass[orig_class] == r:
                                # r is a miss class - output class is orig_class
                                req_tput = 0.0
                                for ist in range(M):
                                    if orig_class < TN.shape[1]:
                                        val = TN[ist, orig_class]
                                        if val > 0 and not np.isnan(val):
                                            req_tput = val
                                            break
                                if req_tput > 0 and orig_class < R:
                                    TNn[ind, orig_class] = req_tput
                else:
                    # For other non-station nodes (Router, etc.)
                    for s in inchain:
                        if s >= R:
                            continue
                        for jnd in range(I):
                            # Use rtnodes for routing probability
                            if sn.rtnodes is not None and sn.rtnodes.size > 0:
                                from_idx = ind * R + r
                                to_idx = jnd * R + s
                                if from_idx < sn.rtnodes.shape[0] and to_idx < sn.rtnodes.shape[1]:
                                    TNn[ind, s] += ANn[ind, r] * sn.rtnodes[from_idx, to_idx]

    # Handle Join nodes
    for ind in range(I):
        if sn.nodetype is None or ind >= len(sn.nodetype):
            continue

        node_type = sn.nodetype[ind]

        if node_type == NodeType.JOIN:
            for c in range(C):
                if c not in sn.inchain:
                    continue

                inchain = sn.inchain[c].flatten().astype(int)

                for r in inchain:
                    if r >= R:
                        continue

                    for s in inchain:
                        if s >= R:
                            continue
                        for jnd in range(I):
                            if sn.rtnodes is not None and sn.rtnodes.size > 0:
                                from_idx = ind * R + r
                                to_idx = jnd * R + s
                                if from_idx < sn.rtnodes.shape[0] and to_idx < sn.rtnodes.shape[1]:
                                    TNn[ind, s] += ANn[ind, r] * sn.rtnodes[from_idx, to_idx]

    # Replace NaN with 0
    TNn = np.nan_to_num(TNn, nan=0.0)

    return TNn


def sn_get_product_form_chain_params(sn: NetworkStruct) -> ChainParams:
    """
    Extract product-form parameters aggregated by chain.

    Extracts parameters from a network structure and aggregates them
    by chain for use in product-form analysis methods.

    Args:
        sn: Network structure

    Returns:
        ChainParams with lambda_vec, D, N, Z, mu, S, V

    References:
        Original MATLAB: matlab/src/api/sn/sn_get_product_form_chain_params.m
    """
    from .transforms import sn_get_product_form_params
    from .demands import sn_get_demands_chain

    # Get base product-form parameters
    params = sn_get_product_form_params(sn)

    # Get chain-aggregated demands
    demands = sn_get_demands_chain(sn)

    # Find queue and delay indices
    # Handle both list and numpy array nodetype
    nodetype = sn.nodetype if isinstance(sn.nodetype, np.ndarray) else np.array([nt.value if hasattr(nt, 'value') else nt for nt in sn.nodetype])
    queue_indices = np.where(nodetype == NodeType.QUEUE.value)[0]
    delay_indices = np.where(nodetype == NodeType.DELAY.value)[0]

    # Initialize chain parameters
    nchains = sn.nchains
    lambda_chains = np.zeros(nchains)

    for c in range(nchains):
        chain_classes = sn.inchain[c]
        lambda_chains[c] = np.nansum(params.lam[chain_classes])

    # Extract demands at queue and delay stations
    # Note: Lchain in Python corresponds to Dchain in MATLAB
    D_chains = demands.Lchain[sn.nodeToStation[queue_indices], :]
    Z_chains = demands.Lchain[sn.nodeToStation[delay_indices], :] if len(delay_indices) > 0 else np.zeros((0, nchains))

    # Number of servers at queuing stations
    S = sn.nservers[sn.nodeToStation[queue_indices]]

    # Visit ratios
    V = demands.Vchain.copy()
    ignore_indices = np.where((nodetype == NodeType.SOURCE.value) | (nodetype == NodeType.JOIN.value))[0]
    if len(ignore_indices) > 0:
        keep_stations = [sn.nodeToStation[i] for i in range(sn.nnodes) if i not in ignore_indices]
        V = V[keep_stations, :]

    if len(Z_chains) == 0:
        Z_chains = np.zeros((0, nchains))

    return ChainParams(
        lambda_vec=lambda_chains,
        D=D_chains,
        N=demands.Nchain,
        Z=np.sum(Z_chains, axis=0) if Z_chains.size > 0 else np.zeros(nchains),
        mu=params.mu,
        S=S,
        V=V
    )


def sn_set_routing_prob(sn: NetworkStruct, from_stateful: int, from_class: int,
                        to_stateful: int, to_class: int, prob: float,
                        auto_refresh: bool = False) -> NetworkStruct:
    """
    Set a routing probability between two stateful node-class pairs.

    Updates a single entry in the rt matrix.

    Args:
        sn: Network structure
        from_stateful: Source stateful node index (0-based)
        from_class: Source class index (0-based)
        to_stateful: Destination stateful node index (0-based)
        to_class: Destination class index (0-based)
        prob: Routing probability [0, 1]
        auto_refresh: If True, refresh visit ratios (default False)

    Returns:
        Modified network structure

    References:
        Original MATLAB: matlab/src/api/sn/sn_set_routing_prob.m
    """
    K = sn.nclasses

    # Calculate indices in rt matrix
    from_idx = from_stateful * K + from_class
    to_idx = to_stateful * K + to_class

    # Update rt matrix
    sn.rt[from_idx, to_idx] = prob

    # Auto-refresh visit ratios if requested
    if auto_refresh:
        from .transforms import sn_refresh_visits
        sn_refresh_visits(sn)

    return sn


__all__ = [
    'ChainParams',
    'sn_get_arvr_from_tput',
    'sn_get_node_arvr_from_tput',
    'sn_get_node_tput_from_tput',
    'sn_get_product_form_chain_params',
    'sn_set_routing_prob',
]
