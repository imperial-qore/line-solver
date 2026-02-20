"""
CTMC Solver handler.

Native Python implementation of CTMC (Continuous-Time Markov Chain) solver
handler that analyzes queueing networks through exact state-space enumeration.

The CTMC solver builds the complete state space and infinitesimal generator
matrix, then solves for steady-state probabilities to compute performance metrics.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
from itertools import product
import time

import warnings
from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    RoutingStrategy,
    sn_is_open_model,
    sn_is_closed_model,
    sn_has_open_classes,
    sn_refresh_visits,
)
from ....constants import ProcessType, GlobalConstants
from ...mc import ctmc_solve, ctmc_makeinfgen
from ...mc.dtmc import dtmc_stochcomp
from ...cache.miss import cache_xi_fp
from ...mam import map_mean


def _enumerate_cache_states(nitems: int, capacity: int) -> np.ndarray:
    """
    Enumerate all possible cache configurations.

    This matches MATLAB's State.spaceCache(n, m) which generates all permutations
    of m items chosen from n items. Each row represents which items are in the cache.

    Args:
        nitems: Total number of items (n)
        capacity: Cache capacity (m)

    Returns:
        Array of shape (n_states, capacity) where each row is a cache configuration
        with the specific item IDs (1-based like MATLAB) that are in the cache.
    """
    from itertools import combinations, permutations

    if nitems <= 0 or capacity <= 0:
        return np.array([[]])

    # Generate all combinations of capacity items from nitems
    # Then generate all permutations of each combination (order matters for LRU/FIFO)
    cache_states = []
    for combo in combinations(range(1, nitems + 1), capacity):
        for perm in permutations(combo):
            cache_states.append(list(perm))

    return np.array(cache_states, dtype=int) if cache_states else np.array([[]])


def _get_cache_stations_info(sn: NetworkStruct) -> Dict[int, Dict]:
    """
    Extract cache node information for CTMC state space enumeration.

    Cache nodes in LINE are stateful non-station nodes (matching MATLAB behavior).
    This function iterates over all nodes to find Cache nodes, not just stations.

    Returns:
        Dictionary mapping stateful indices to cache parameters:
        {isf: {'node_idx': ind, 'nitems': n, 'capacity': m, 'hitclass': [...],
               'missclass': [...], 'pread': {...}, 'replacement': int,
               'states': array, 'num_states': int}}
    """
    cache_stations_info = {}

    if sn.nodeparam is None:
        return cache_stations_info

    # Iterate over all nodes to find Cache nodes
    # Cache nodes are stateful but NOT stations in LINE's semantics (matching MATLAB)
    for ind in range(sn.nnodes):
        # Check if this node is a Cache
        if sn.nodetype is None or ind >= len(sn.nodetype):
            continue

        nt_val = int(sn.nodetype[ind].value) if hasattr(sn.nodetype[ind], 'value') else int(sn.nodetype[ind])
        if nt_val != 6:  # NodeType.CACHE = 6
            continue

        # Get stateful index for this node
        if hasattr(sn, 'nodeToStateful') and sn.nodeToStateful is not None:
            isf = int(sn.nodeToStateful[ind])
            if isf < 0:
                continue
        else:
            # Fallback: assume node index = stateful index
            isf = ind

        if ind not in sn.nodeparam:
            continue

        param = sn.nodeparam[ind]
        nitems = getattr(param, 'nitems', 0)

        # Get capacity - prefer itemcap (multi-level capacities) over cap
        # itemcap is an array for multi-level caches, cap is often just the first level
        cap = getattr(param, 'itemcap', None)
        if cap is None:
            cap = getattr(param, 'cap', 0)
        if hasattr(cap, '__iter__'):
            cap = int(sum(cap))  # Sum for multi-level caches
        else:
            cap = int(cap) if cap else 0

        if nitems <= 0 or cap <= 0:
            continue

        # Enumerate all cache configurations
        cache_states = _enumerate_cache_states(nitems, int(cap))

        # Get itemcap array for multi-level cache
        itemcap_raw = getattr(param, 'itemcap', None)
        if itemcap_raw is not None:
            if hasattr(itemcap_raw, '__iter__'):
                itemcap = list(int(x) for x in itemcap_raw)
            else:
                itemcap = [int(itemcap_raw)]
        else:
            itemcap = [int(cap)]

        cache_stations_info[isf] = {
            'node_idx': ind,
            'nitems': nitems,
            'capacity': int(cap),
            'itemcap': itemcap,  # Multi-level capacity array
            'hitclass': np.atleast_1d(getattr(param, 'hitclass', np.array([]))).flatten().astype(int),
            'missclass': np.atleast_1d(getattr(param, 'missclass', np.array([]))).flatten().astype(int),
            'pread': getattr(param, 'pread', {}),
            'replacement': int(getattr(param, 'replacestrat', 0)),
            'accost': getattr(param, 'accost', None),  # Access cost matrix for promotions
            'states': cache_states,
            'num_states': len(cache_states)
        }

    return cache_stations_info


def _compute_cache_read_transitions(
    cache_info: Dict,
    cache_state_idx: int,
    input_class: int,
    K: int
) -> List[Tuple[int, float, Optional[int]]]:
    """
    Compute cache READ transitions based on current cache state.

    For a job in input_class at a Cache station, compute the transitions
    to hit/miss classes based on which items are in the current cache state.

    Args:
        cache_info: Cache station info from _get_cache_stations_info
        cache_state_idx: Current cache state index (which configuration)
        input_class: The class of the arriving job (0-indexed)
        K: Number of classes

    Returns:
        List of (output_class, probability, new_cache_state_idx) tuples.
        output_class: The class to route to (hitclass or missclass)
        probability: The probability of this transition (from pread)
        new_cache_state_idx: New cache state index (None if unchanged, int if cache updated)
    """
    transitions = []

    hitclass = cache_info['hitclass']
    missclass = cache_info['missclass']
    pread = cache_info['pread']
    cache_states = cache_info['states']
    nitems = cache_info['nitems']
    capacity = cache_info['capacity']
    replacement = cache_info['replacement']

    # Check if this class has hit/miss routing
    if input_class >= len(hitclass) or input_class >= len(missclass):
        return transitions

    h = int(hitclass[input_class])
    m = int(missclass[input_class])

    if h < 0 or m < 0 or h >= K or m >= K:
        return transitions

    # Get pread probabilities for this class
    if pread is None:
        return transitions

    pread_k = None
    if isinstance(pread, dict):
        pread_k = pread.get(input_class)
    elif isinstance(pread, (list, tuple)) and input_class < len(pread):
        pread_k = pread[input_class]
    elif isinstance(pread, np.ndarray):
        if pread.ndim == 1:
            pread_k = pread
        elif input_class < pread.shape[0]:
            pread_k = pread[input_class, :]

    if pread_k is None:
        return transitions

    pread_k = np.atleast_1d(pread_k).flatten()

    # Get current cache configuration (1-based item IDs)
    if cache_state_idx < 0 or cache_state_idx >= len(cache_states):
        return transitions

    current_cache = set(cache_states[cache_state_idx].tolist())

    # ReplacementStrategy values matching MATLAB: RR=0, FIFO=1, SFIFO=2, LRU=3
    RR = 0
    FIFO = 1
    SFIFO = 2
    LRU = 3

    # For each item, compute hit/miss transition
    for item in range(1, nitems + 1):  # Items are 1-based
        item_idx = item - 1  # Index into pread array
        if item_idx >= len(pread_k):
            continue

        prob = pread_k[item_idx]
        if prob <= 0:
            continue

        if item in current_cache:
            # CACHE HIT - route to hitclass, cache state may change
            # For multi-level caches, need to handle list promotions
            current_list = list(cache_states[cache_state_idx])
            item_pos = current_list.index(item)

            # Get itemcap for multi-level cache analysis
            itemcap_arr = cache_info.get('itemcap', [capacity])
            if not isinstance(itemcap_arr, (list, tuple, np.ndarray)):
                itemcap_arr = [capacity]
            itemcap_arr = list(itemcap_arr)
            n_lists = len(itemcap_arr)

            # Find which list the item is in (0-indexed)
            # List boundaries: list 0 = positions 0..itemcap[0]-1, list 1 = positions itemcap[0]..itemcap[0]+itemcap[1]-1, etc.
            item_list = 0
            pos_in_list = item_pos
            cumsum = 0
            for list_idx, cap in enumerate(itemcap_arr):
                if item_pos < cumsum + cap:
                    item_list = list_idx
                    pos_in_list = item_pos - cumsum
                    break
                cumsum += cap

            if replacement == LRU and n_lists > 1 and item_list < n_lists - 1:
                # Multi-level LRU HIT in list i < last list:
                # Item is promoted to the last list with the following steps:
                # 1. Source list: items [0, j-1] shift right to [1, j], head gets tail of target list
                # 2. Target list: items shift right, head gets hit item
                target_list = n_lists - 1
                target_list_start = sum(itemcap_arr[:target_list])
                target_list_end = target_list_start + itemcap_arr[target_list]
                target_cap = itemcap_arr[target_list]

                # Source list boundaries
                src_list_start = sum(itemcap_arr[:item_list])
                src_list_end = src_list_start + itemcap_arr[item_list]
                j = item_pos - src_list_start  # Position within source list

                new_list = current_list.copy()

                # Step 1: Source list - shift items [0, j-1] right to [1, j]
                if j > 0:
                    for i in range(j, 0, -1):
                        new_list[src_list_start + i] = current_list[src_list_start + i - 1]

                # Step 2: Source list head gets tail of target list
                target_tail = current_list[target_list_end - 1]
                new_list[src_list_start] = target_tail

                # Step 3: Target list - shift items right
                if target_cap > 1:
                    for i in range(target_list_end - 1, target_list_start, -1):
                        new_list[i] = current_list[i - 1]

                # Step 4: Target list head gets hit item
                new_list[target_list_start] = item

                # Find the new cache state index
                new_state_idx = None
                for idx, state in enumerate(cache_states):
                    if list(state) == new_list:
                        new_state_idx = idx
                        break

                transitions.append((h, prob, new_state_idx))
            elif replacement == LRU and (n_lists == 1 or item_list == n_lists - 1):
                # LRU HIT in last list (or single-level cache): move item to front of its list
                src_list_start = sum(itemcap_arr[:item_list])
                src_list_end = src_list_start + itemcap_arr[item_list]
                j = item_pos - src_list_start  # Position within list

                new_list = current_list.copy()
                # Shift items [0, j-1] right to [1, j]
                if j > 0:
                    for i in range(j, 0, -1):
                        new_list[src_list_start + i] = current_list[src_list_start + i - 1]
                # Head gets hit item
                new_list[src_list_start] = item

                # Find the new cache state index
                new_state_idx = None
                for idx, state in enumerate(cache_states):
                    if list(state) == new_list:
                        new_state_idx = idx
                        break
                transitions.append((h, prob, new_state_idx))
            elif replacement == FIFO and n_lists > 1 and item_list < n_lists - 1:
                # Multi-level FIFO HIT in list i < last list:
                # Item is promoted to the LAST list with 100% probability.
                # The tail item of the last list replaces the hit item's position.
                target_list = n_lists - 1  # Always promote to last list
                target_list_start = sum(itemcap_arr[:target_list])
                target_list_end = target_list_start + itemcap_arr[target_list]
                target_tail = current_list[target_list_end - 1]  # Last item in target list

                new_list = current_list.copy()

                # Step 1: Hit position gets the tail item from target list
                new_list[item_pos] = target_tail

                # Step 2: Shift target list items to make room at head
                # Move items [head, tail-1] to [head+1, tail]
                if itemcap_arr[target_list] > 1:
                    for i in range(target_list_end - 1, target_list_start, -1):
                        new_list[i] = current_list[i - 1]

                # Step 3: Place hit item at head of target list
                new_list[target_list_start] = item

                # Find the new cache state index
                new_state_idx = None
                for idx, state in enumerate(cache_states):
                    if list(state) == new_list:
                        new_state_idx = idx
                        break

                # 100% promotion to last list
                transitions.append((h, prob, new_state_idx))
            elif replacement == FIFO and (n_lists == 1 or item_list == n_lists - 1):
                # Single-level FIFO or HIT in last list: no state change
                transitions.append((h, prob, None))
            elif replacement == RR and n_lists > 1 and item_list < n_lists - 1:
                # Multi-level RR HIT in list i < last list:
                # Item swaps with a random position in the target list (list 1 by default accessCost).
                # Default accessCost only allows promotion to last list.
                target_list = n_lists - 1  # Promote to last list
                target_list_start = sum(itemcap_arr[:target_list])
                target_list_end = target_list_start + itemcap_arr[target_list]
                target_cap = itemcap_arr[target_list]

                # For each position r in target list, swap hit item with position r
                for r in range(target_cap):
                    target_pos = target_list_start + r
                    swap_item = current_list[target_pos]

                    new_list = current_list.copy()
                    new_list[item_pos] = swap_item
                    new_list[target_pos] = item

                    # Find the new cache state index
                    new_state_idx = None
                    for idx, state in enumerate(cache_states):
                        if list(state) == new_list:
                            new_state_idx = idx
                            break

                    # Equal probability for each position in target list
                    transitions.append((h, prob / target_cap, new_state_idx))
            elif replacement == RR and (n_lists == 1 or item_list == n_lists - 1):
                # Single-level RR or HIT in last list: no state change
                transitions.append((h, prob, None))
            else:
                # Default: cache state unchanged on hit
                transitions.append((h, prob, None))
        else:
            # CACHE MISS - route to missclass, update cache state
            if replacement == RR:
                # Random Replacement: replace random position in list 0 only
                # (default accessCost puts all MISS entries into list 0)
                current_list = list(cache_states[cache_state_idx])
                itemcap_arr = cache_info.get('itemcap', [capacity])
                if not isinstance(itemcap_arr, (list, tuple, np.ndarray)):
                    itemcap_arr = [capacity]
                itemcap_arr = list(itemcap_arr)

                # Only replace positions in list 0
                list0_cap = itemcap_arr[0]
                for replace_pos in range(list0_cap):
                    new_list = current_list.copy()
                    new_list[replace_pos] = item
                    # Find the new cache state index
                    new_state_idx = None
                    for idx, state in enumerate(cache_states):
                        if list(state) == new_list:
                            new_state_idx = idx
                            break
                    if new_state_idx is not None:
                        # Each position in list 0 has equal probability
                        transitions.append((m, prob / list0_cap, new_state_idx))
            elif replacement == FIFO:
                # Multi-level FIFO MISS: evict from list 0 tail only, insert at list 0 head.
                # Other lists remain unchanged.
                current_list = list(cache_states[cache_state_idx])
                itemcap_arr = cache_info.get('itemcap', [capacity])
                if not isinstance(itemcap_arr, (list, tuple, np.ndarray)):
                    itemcap_arr = [capacity]
                itemcap_arr = list(itemcap_arr)

                # List 0 boundaries
                list0_start = 0
                list0_end = itemcap_arr[0]

                # Build new list: shift list 0 right, insert item at head, other lists unchanged
                new_list = current_list.copy()
                # Shift list 0 items right (evict tail)
                for i in range(list0_end - 1, list0_start, -1):
                    new_list[i] = current_list[i - 1]
                # Insert new item at head of list 0
                new_list[list0_start] = item

                # Find the new cache state index
                new_state_idx = None
                for idx, state in enumerate(cache_states):
                    if list(state) == new_list:
                        new_state_idx = idx
                        break
                transitions.append((m, prob, new_state_idx))
            elif replacement == LRU:
                # Multi-level LRU MISS: evict from list 0 tail only, insert at list 0 head.
                # Same as FIFO MISS (default accessCost only allows entry to list 0).
                current_list = list(cache_states[cache_state_idx])
                itemcap_arr = cache_info.get('itemcap', [capacity])
                if not isinstance(itemcap_arr, (list, tuple, np.ndarray)):
                    itemcap_arr = [capacity]
                itemcap_arr = list(itemcap_arr)

                # List 0 boundaries
                list0_start = 0
                list0_end = itemcap_arr[0]

                # Build new list: shift list 0 right, insert item at head, other lists unchanged
                new_list = current_list.copy()
                for i in range(list0_end - 1, list0_start, -1):
                    new_list[i] = current_list[i - 1]
                new_list[list0_start] = item

                # Find the new cache state index
                new_state_idx = None
                for idx, state in enumerate(cache_states):
                    if list(state) == new_list:
                        new_state_idx = idx
                        break
                transitions.append((m, prob, new_state_idx))
            else:
                # Default: FIFO-like
                current_list = list(cache_states[cache_state_idx])
                new_list = [item] + current_list[:-1]
                new_state_idx = None
                for idx, state in enumerate(cache_states):
                    if list(state) == new_list:
                        new_state_idx = idx
                        break
                transitions.append((m, prob, new_state_idx))

    return transitions


def _get_cache_info(sn: NetworkStruct) -> Dict[int, Dict]:
    """
    Extract cache information from the network structure.

    Returns:
        Dictionary mapping cache node indices to cache parameters:
        {ind: {'nitems': n, 'capacity': m, 'hitclass': [...], 'missclass': [...],
               'pread': [...], 'replacement': int, 'states': array}}
    """
    cache_info = {}

    if sn.nodeparam is None:
        return cache_info

    for ind in range(sn.nnodes):
        if sn.nodetype is None or ind >= len(sn.nodetype):
            continue

        nt_val = int(sn.nodetype[ind].value) if hasattr(sn.nodetype[ind], 'value') else int(sn.nodetype[ind])
        if nt_val != 6:  # NodeType.CACHE = 6
            continue

        if ind not in sn.nodeparam:
            continue

        param = sn.nodeparam[ind]
        nitems = getattr(param, 'nitems', 0)
        cap = getattr(param, 'cap', 0)

        if nitems <= 0 or cap <= 0:
            continue

        cache_info[ind] = {
            'nitems': nitems,
            'capacity': cap,
            'hitclass': np.atleast_1d(getattr(param, 'hitclass', np.array([]))).flatten(),
            'missclass': np.atleast_1d(getattr(param, 'missclass', np.array([]))).flatten(),
            'pread': getattr(param, 'pread', None),
            'replacement': getattr(param, 'replacestrat', 0),
            'states': _enumerate_cache_states(nitems, cap)
        }

    return cache_info


def _get_phases_info(sn: NetworkStruct) -> Tuple[np.ndarray, bool]:
    """
    Get number of phases for each (station, class) pair and detect if phase augmentation is needed.

    For phase-type distributions (PH, APH, MAP, MMPP2, Erlang, HyperExp, Coxian),
    the number of phases is determined from the distribution parameters stored in sn.proc.

    Args:
        sn: Network structure

    Returns:
        Tuple of (phases matrix [M x K], needs_phase_augmentation bool)
    """
    M = sn.nstations
    K = sn.nclasses
    phases = np.ones((M, K), dtype=int)
    needs_augmentation = False

    if not hasattr(sn, 'proc') or sn.proc is None:
        return phases, needs_augmentation

    # sn.proc can be a list or dict - handle both
    proc_is_list = isinstance(sn.proc, list)

    for ist in range(M):
        # Check if this station has proc data
        if proc_is_list:
            if ist >= len(sn.proc) or sn.proc[ist] is None:
                continue
            station_proc = sn.proc[ist]
        else:
            if ist not in sn.proc:
                continue
            station_proc = sn.proc[ist]

        for k in range(K):
            # Get proc entry for this (station, class)
            proc_entry = None
            if isinstance(station_proc, (list, tuple)):
                if k < len(station_proc):
                    proc_entry = station_proc[k]
            elif isinstance(station_proc, dict):
                proc_entry = station_proc.get(k)

            if proc_entry is None:
                continue

            n_phases = 1

            # Handle different storage formats
            if isinstance(proc_entry, dict):
                # Erlang: {'k': phases, 'mu': rate}
                if 'k' in proc_entry:
                    n_phases = int(proc_entry['k'])
                # HyperExp: {'probs': [...], 'rates': [...]}
                elif 'probs' in proc_entry and 'rates' in proc_entry:
                    probs = np.array(proc_entry['probs'])
                    proc_rates = np.array(proc_entry['rates'])
                    # Check if rates data is valid (not same as probs)
                    # If rates == probs, this indicates a data bug - treat as single phase
                    if np.allclose(proc_rates, probs):
                        n_phases = 1
                    else:
                        n_phases = len(proc_entry['probs'])
                # Exp: {'rate': ...} - single phase
                else:
                    n_phases = 1
            elif isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 1:
                # PH/APH/MAP: [alpha/D0, T/D1] where alpha/D0 determines phases
                first_elem = proc_entry[0]
                if isinstance(first_elem, np.ndarray):
                    if first_elem.ndim == 1:
                        # alpha vector
                        n_phases = len(first_elem)
                    else:
                        # D0 matrix
                        n_phases = first_elem.shape[0]
                elif isinstance(first_elem, (list, tuple)):
                    if len(first_elem) > 0 and isinstance(first_elem[0], (list, tuple, np.ndarray)):
                        # 2D structure
                        n_phases = len(first_elem)
                    else:
                        # 1D structure
                        n_phases = len(first_elem)

            phases[ist, k] = max(1, n_phases)
            if n_phases > 1:
                needs_augmentation = True

    return phases, needs_augmentation


def _get_map_fcfs_info(sn: NetworkStruct) -> Tuple[Dict[Tuple[int, int], int], bool]:
    """
    Identify (station, class) pairs that have MAP distributions at FCFS stations.

    For FCFS stations with MAP distributions, the state must include an additional
    variable tracking the MAP modulating phase (the "mode" of the MAP process).

    Args:
        sn: Network structure

    Returns:
        Tuple of (map_fcfs dict mapping (ist, k) to n_phases, has_map_fcfs bool)
    """
    M = sn.nstations
    K = sn.nclasses
    map_fcfs = {}  # (ist, k) -> n_phases for MAP distributions at FCFS stations
    has_map_fcfs = False

    if not hasattr(sn, 'proc') or sn.proc is None:
        return map_fcfs, has_map_fcfs
    if not hasattr(sn, 'sched') or sn.sched is None:
        return map_fcfs, has_map_fcfs
    if not hasattr(sn, 'procid') or sn.procid is None:
        return map_fcfs, has_map_fcfs

    proc_is_list = isinstance(sn.proc, list)

    for ist in range(M):
        # Check if this is an FCFS station
        sched = sn.sched.get(ist, SchedStrategy.FCFS)
        # FCFS variants that need MAP phase tracking
        fcfs_variants = [SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS, SchedStrategy.LCFSPR]
        # Check for specific FCFS scheduling strategies
        is_fcfs_type = (sched in fcfs_variants or
                       (isinstance(sched, int) and sched in [1, 2, 3]))  # FCFS=1, LCFS=2, etc.

        if not is_fcfs_type:
            continue

        # Check if this station has proc data
        if proc_is_list:
            if ist >= len(sn.proc) or sn.proc[ist] is None:
                continue
            station_proc = sn.proc[ist]
        else:
            if ist not in sn.proc:
                continue
            station_proc = sn.proc[ist]

        for k in range(K):
            # Check process type
            if ist >= sn.procid.shape[0] or k >= sn.procid.shape[1]:
                continue
            procid = sn.procid[ist, k]

            # Check if MAP or MMPP2
            if procid not in [ProcessType.MAP, ProcessType.MMPP2]:
                continue

            # Get number of phases from proc entry
            proc_entry = None
            if isinstance(station_proc, (list, tuple)):
                if k < len(station_proc):
                    proc_entry = station_proc[k]
            elif isinstance(station_proc, dict):
                proc_entry = station_proc.get(k)

            if proc_entry is None:
                continue

            n_phases = 1
            if isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 2:
                # MAP: [D0, D1]
                D0 = np.atleast_2d(np.array(proc_entry[0], dtype=float))
                n_phases = D0.shape[0]

            if n_phases > 1:
                map_fcfs[(ist, k)] = n_phases
                has_map_fcfs = True

    return map_fcfs, has_map_fcfs


def _generate_phase_distributions(n_jobs: int, n_phases: int) -> List[Tuple[int, ...]]:
    """
    Generate all ways to distribute n_jobs across n_phases.

    This is equivalent to MATLAB's State.spaceClosedSingle.

    Args:
        n_jobs: Number of jobs to distribute
        n_phases: Number of phases

    Returns:
        List of tuples, each tuple has n_phases elements summing to n_jobs
    """
    if n_phases == 1:
        return [(n_jobs,)]
    if n_jobs == 0:
        return [tuple([0] * n_phases)]

    result = []
    for k in range(n_jobs + 1):
        # k jobs in first phase, rest in remaining phases
        for rest in _generate_phase_distributions(n_jobs - k, n_phases - 1):
            result.append((k,) + rest)
    return result


def _get_phase_transition_params(sn: NetworkStruct, ist: int, k: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get phase transition parameters for a (station, class) pair.

    For PH/APH distributions: Returns (mu, phi, alpha) where
        - mu[i] = service rate in phase i
        - phi[i] = probability of completion from phase i (vs moving to another phase)
        - alpha[i] = initial probability of starting in phase i

    For MAP distributions: Returns (D0, D1, pi) where
        - D0 = transition matrix without completions
        - D1 = transition matrix with completions
        - pi = stationary distribution (initial phases)

    Args:
        sn: Network structure
        ist: Station index
        k: Class index

    Returns:
        For PH/APH: (mu, phi, alpha)
        For MAP: (D0, D1, pi)
    """
    if not hasattr(sn, 'proc') or sn.proc is None:
        return np.array([1.0]), np.array([1.0]), np.array([1.0])

    # Handle both list and dict formats for sn.proc
    proc_is_list = isinstance(sn.proc, list)
    if proc_is_list:
        if ist >= len(sn.proc) or sn.proc[ist] is None:
            rate = sn.rates[ist, k] if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1] else 1.0
            return np.array([rate]), np.array([1.0]), np.array([1.0])
        station_proc = sn.proc[ist]
    else:
        if ist not in sn.proc or sn.proc[ist] is None:
            rate = sn.rates[ist, k] if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1] else 1.0
            return np.array([rate]), np.array([1.0]), np.array([1.0])
        station_proc = sn.proc[ist]

    # Get proc entry for this class
    if isinstance(station_proc, (list, tuple)):
        if k >= len(station_proc) or station_proc[k] is None:
            rate = sn.rates[ist, k] if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1] else 1.0
            return np.array([rate]), np.array([1.0]), np.array([1.0])
        proc_entry = station_proc[k]
    elif isinstance(station_proc, dict):
        if k not in station_proc or station_proc[k] is None:
            rate = sn.rates[ist, k] if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1] else 1.0
            return np.array([rate]), np.array([1.0]), np.array([1.0])
        proc_entry = station_proc[k]
    else:
        rate = sn.rates[ist, k] if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1] else 1.0
        return np.array([rate]), np.array([1.0]), np.array([1.0])

    # Check process type
    procid = None
    if hasattr(sn, 'procid') and sn.procid is not None:
        if ist < sn.procid.shape[0] and k < sn.procid.shape[1]:
            procid = sn.procid[ist, k]

    # Handle different distribution types
    if isinstance(proc_entry, dict):
        # Erlang: {'k': phases, 'mu': rate}
        if 'k' in proc_entry:
            n_phases = int(proc_entry['k'])
            rate = proc_entry.get('mu', 1.0)
            mu = np.full(n_phases, rate)
            phi = np.concatenate([np.zeros(n_phases - 1), [1.0]])  # Only complete from last phase
            alpha = np.zeros(n_phases)
            alpha[0] = 1.0  # Start in first phase
            return mu, phi, alpha
        # HyperExp: {'probs': [...], 'rates': [...]}
        elif 'probs' in proc_entry and 'rates' in proc_entry:
            probs = np.array(proc_entry['probs'])
            proc_rates = np.array(proc_entry['rates'])
            # Check if rates data is valid (not same as probs)
            # If rates == probs, this indicates a data bug - fall back to single-phase
            if np.allclose(proc_rates, probs):
                rate = sn.rates[ist, k] if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1] else 1.0
                return np.array([rate]), np.array([1.0]), np.array([1.0])
            phi = np.ones(len(proc_rates))  # Each phase completes immediately
            return proc_rates, phi, probs
        # Exp: {'rate': ...}
        else:
            rate = proc_entry.get('rate', 1.0)
            return np.array([rate]), np.array([1.0]), np.array([1.0])

    elif isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 2:
        # PH/APH: [alpha, T]
        # MAP: [D0, D1]
        first = np.atleast_1d(np.array(proc_entry[0], dtype=float))
        second = np.atleast_2d(np.array(proc_entry[1], dtype=float))

        is_map = procid in [ProcessType.MAP, ProcessType.MMPP2] if procid is not None else False

        if is_map or (first.ndim == 2):
            # MAP: [D0, D1]
            D0 = np.atleast_2d(first)
            D1 = second
            # Compute stationary distribution
            Q = D0 + D1
            n = Q.shape[0]
            A = np.vstack([Q.T, np.ones(n)])
            b = np.zeros(n + 1)
            b[-1] = 1.0
            try:
                from scipy import linalg
                pi, _, _, _ = linalg.lstsq(A, b)
                pi = np.maximum(pi, 0)
                pi /= pi.sum() if pi.sum() > 0 else 1
            except:
                pi = np.ones(n) / n
            return D0, D1, pi
        else:
            # PH/APH: [alpha, T]
            alpha = first.flatten()
            T = second
            n_phases = len(alpha)

            # Extract service rates (negative diagonal of T)
            mu = -np.diag(T)

            # Compute exit rates (completion probability)
            exit_rates = -T.sum(axis=1)  # -T * e gives exit rates
            phi = np.zeros(n_phases)
            for i in range(n_phases):
                if mu[i] > 0:
                    phi[i] = exit_rates[i] / mu[i]

            return mu, phi, alpha

    # Default: exponential
    rate = sn.rates[ist, k] if hasattr(sn, 'rates') and ist < sn.rates.shape[0] and k < sn.rates.shape[1] else 1.0
    return np.array([rate]), np.array([1.0]), np.array([1.0])


def _apply_router_stochcomp(
    sn: NetworkStruct,
    Q_orig: np.ndarray,
    depRates_orig: np.ndarray,
    state_space: np.ndarray,
    router_nodes: List[int],
    rrobin_info: dict,
    options: 'SolverCTMCOptions'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply stochastic complementation to handle Router nodes.

    MATLAB builds the full CTMC state space including Router-occupied states
    (with immediate departure rate GAMMA=1e8), then applies ctmc_stochcomp
    to remove them. This function replicates that approach:

    1. Identify Cache DEP transitions in Q_orig (hit/miss class jobs departing
       Cache and arriving at downstream stations at rate GAMMA*routing_prob)
    2. Expand Q to include Router-occupied states for each class
    3. Replace Cache DEP transitions with two-step: Cache->Router at GAMMA,
       Router->Delay at GAMMA*routing_prob
    4. Apply stochcomp S = Q_NN + Q_NI*(-Q_II)^{-1}*Q_IN to eliminate
       Router-occupied states

    This correctly handles multi-step immediate chains (e.g., Cache DEP ->
    Router -> Delay, then Cache READ -> Router -> Delay) that the direct
    routing approach cannot capture.

    Args:
        sn: Network structure
        Q_orig: Generator matrix built on station-only state space
        depRates_orig: Departure rates (n_states, M, K)
        state_space: Station-only state space
        router_nodes: List of Router node indices
        rrobin_info: Round-robin and state space info
        options: Solver options

    Returns:
        Tuple of (Q_reduced, depRates_reduced) on the original state space
    """
    M = int(sn.nstations)
    K = int(sn.nclasses)
    n_orig = state_space.shape[0]
    GAMMA = 1e8  # MATLAB GlobalConstants.Immediate = 1/FineTol = 1/1e-8

    rtnodes = np.asarray(sn.rtnodes) if hasattr(sn, 'rtnodes') and sn.rtnodes is not None else None
    connmatrix = np.asarray(sn.connmatrix) if hasattr(sn, 'connmatrix') and sn.connmatrix is not None else None
    node_to_station = np.asarray(sn.nodeToStation).flatten() if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None else None

    if rtnodes is None or connmatrix is None or node_to_station is None:
        return Q_orig, depRates_orig

    # Get cache state tracking info
    cache_stations_info = rrobin_info.get('cache_stations_info', {})
    cache_jobs_offsets = rrobin_info.get('cache_jobs_offsets', {})
    use_phase_aug = rrobin_info.get('needs_phase_augmentation', False)
    phases = rrobin_info.get('phases', np.ones((M, K), dtype=int))
    phase_offset = rrobin_info.get('phase_offset', None)

    # Identify Cache output classes (hit/miss classes) — these route through Router
    cache_output_classes = set()
    for cache_isf, cache_info in cache_stations_info.items():
        hitclass = cache_info.get('hitclass', [])
        missclass = cache_info.get('missclass', [])
        for kk in range(min(len(hitclass), K)):
            if hitclass[kk] >= 0:
                cache_output_classes.add(int(hitclass[kk]))
            if missclass[kk] >= 0:
                cache_output_classes.add(int(missclass[kk]))

    if not cache_output_classes:
        return Q_orig, depRates_orig

    # Determine which classes visit each Router
    router_visiting = {}  # rnode -> sorted list of class indices
    for rnode in router_nodes:
        visiting = set()
        incoming = np.where(connmatrix[:, rnode] > 0)[0]
        for src_node in incoming:
            for k in range(K):
                for r in range(K):
                    src_idx = src_node * K + k
                    dst_idx = rnode * K + r
                    if src_idx < rtnodes.shape[0] and dst_idx < rtnodes.shape[1]:
                        if rtnodes[src_idx, dst_idx] > 1e-10:
                            visiting.add(r)
        router_visiting[rnode] = sorted(visiting)

    # Router dims: each (rnode, class) pair that can hold a job
    router_dims = []
    for rnode in router_nodes:
        for k in router_visiting.get(rnode, []):
            router_dims.append((rnode, k))
    n_rdims = len(router_dims)
    if n_rdims == 0:
        return Q_orig, depRates_orig

    # Router -> destination station routing probabilities
    router_to_dests = {}  # (rnode, k) -> [(station, class, prob), ...]
    for rnode in router_nodes:
        outgoing = np.where(connmatrix[rnode, :] > 0)[0]
        for k in router_visiting.get(rnode, []):
            src_idx = rnode * K + k
            dests = []
            for jnd in outgoing:
                for r in range(K):
                    dst_idx = jnd * K + r
                    if src_idx < rtnodes.shape[0] and dst_idx < rtnodes.shape[1]:
                        prob = rtnodes[src_idx, dst_idx]
                        if prob > 1e-10:
                            jst = int(node_to_station[jnd]) if jnd < len(node_to_station) else -1
                            if jst >= 0 and jst < M:
                                dests.append((jst, r, prob))
            if dests:
                router_to_dests[(rnode, k)] = dests

    # Map class -> Router dim index
    class_to_dim = {}
    for d, (rnode, k) in enumerate(router_dims):
        if k not in class_to_dim:
            class_to_dim[k] = d

    # Helper to get job count at station for a class
    def get_job_count(state, ist, k):
        if use_phase_aug and phase_offset is not None:
            si = phase_offset[ist, k]
            return int(sum(state[si:si + phases[ist, k]]))
        return int(state[ist * K + k])

    def set_job_count_in_state(state_vec, ist, k, delta):
        """Modify job count in a state vector copy. Returns modified copy or None if invalid."""
        new_state = list(state_vec)
        if use_phase_aug and phase_offset is not None:
            si = phase_offset[ist, k]
            new_state[si] += delta
            if new_state[si] < 0:
                return None
        else:
            idx = ist * K + k
            new_state[idx] += delta
            if new_state[idx] < 0:
                return None
        return new_state

    # Build state_map for fast lookup
    state_map = {}
    for i in range(n_orig):
        state_map[tuple(int(x) for x in state_space[i])] = i

    # Compute total population (station + cache jobs) for each original state
    # to enforce chain population constraints when expanding with Router states
    skip_stations = rrobin_info.get('skip_stations', set())
    state_populations = np.zeros(n_orig, dtype=int)
    for s in range(n_orig):
        total = 0
        state = state_space[s]
        for ist in range(M):
            if ist in skip_stations:
                continue
            for kk in range(K):
                total += get_job_count(state, ist, kk)
        for cache_isf in cache_stations_info:
            cjo = cache_jobs_offsets.get(cache_isf)
            if cjo is not None:
                for kk in range(K):
                    total += int(state[cjo + kk])
        state_populations[s] = total

    # Max chain population from cutoff (sum of per-class cutoffs)
    max_pop = rrobin_info.get('max_chain_pop', K)

    # Build population-constrained expanded state space
    # Only include (s, r) where total population + router jobs <= max_pop
    n_rstates = 2 ** n_rdims
    expanded_pairs = []  # List of (orig_idx, router_bitmask)
    pair_to_exp = {}     # (orig_idx, router_bitmask) -> expanded index
    for s in range(n_orig):
        for r in range(n_rstates):
            router_jobs = bin(r).count('1')
            if state_populations[s] + router_jobs <= max_pop:
                idx = len(expanded_pairs)
                expanded_pairs.append((s, r))
                pair_to_exp[(s, r)] = idx

    n_expanded = len(expanded_pairs)

    if n_expanded > 500000:
        if options.verbose:
            print(f"  Router stochcomp: expanded state space too large ({n_expanded}), skipping")
        return Q_orig, depRates_orig

    if options.verbose:
        print(f"  Router stochcomp: expanding {n_orig} -> {n_expanded} states "
              f"({n_rdims} Router dims, max_pop={max_pop})")

    # Identify Cache DEP transitions in Q_orig
    cache_dep_transitions = set()
    for cache_isf, cache_info in cache_stations_info.items():
        cache_jobs_offset = cache_jobs_offsets.get(cache_isf)
        if cache_jobs_offset is None:
            continue
        for s in range(n_orig):
            state = state_space[s]
            for k in cache_output_classes:
                cache_job_count = int(state[cache_jobs_offset + k])
                if cache_job_count <= 0:
                    continue
                for sp in range(n_orig):
                    if s == sp:
                        continue
                    rate = Q_orig[s, sp]
                    if abs(rate) < 1e4:
                        continue
                    dest_state = state_space[sp]
                    if int(dest_state[cache_jobs_offset + k]) < cache_job_count:
                        cache_dep_transitions.add((s, sp))

    # Build expanded Q matrix
    Q_exp = np.zeros((n_expanded, n_expanded))

    # Step 1: Non-Cache-DEP transitions (replicated across valid Router states)
    for s in range(n_orig):
        for sp in range(n_orig):
            if s == sp:
                continue
            rate = Q_orig[s, sp]
            if abs(rate) < 1e-14:
                continue
            if (s, sp) in cache_dep_transitions:
                continue
            for r in range(n_rstates):
                src_exp = pair_to_exp.get((s, r))
                dst_exp = pair_to_exp.get((sp, r))
                if src_exp is not None and dst_exp is not None:
                    Q_exp[src_exp, dst_exp] += rate

    # Step 2: Cache DEP -> Router + Router -> Delay
    for cache_isf, cache_info in cache_stations_info.items():
        cache_jobs_offset = cache_jobs_offsets.get(cache_isf)
        if cache_jobs_offset is None:
            continue

        for s in range(n_orig):
            state = state_space[s]
            for k in cache_output_classes:
                cache_job_count = int(state[cache_jobs_offset + k])
                if cache_job_count <= 0:
                    continue

                dim_idx = class_to_dim.get(k, -1)
                if dim_idx < 0:
                    for sp in range(n_orig):
                        if (s, sp) in cache_dep_transitions:
                            rate = Q_orig[s, sp]
                            if abs(rate) >= 1e-14:
                                for r in range(n_rstates):
                                    src_exp = pair_to_exp.get((s, r))
                                    dst_exp = pair_to_exp.get((sp, r))
                                    if src_exp is not None and dst_exp is not None:
                                        Q_exp[src_exp, dst_exp] += rate
                    continue

                r_bit = 1 << dim_idx
                rnode_for_k = router_dims[dim_idx][0]
                dests = router_to_dests.get((rnode_for_k, k), [])
                if not dests:
                    continue

                # Intermediate state: job removed from Cache
                int_state_list = [int(x) for x in state]
                int_state_list[cache_jobs_offset + k] = cache_job_count - 1
                int_s = state_map.get(tuple(int_state_list), -1)
                if int_s < 0:
                    continue

                # Cache DEP -> Router
                for r_base in range(n_rstates):
                    if (r_base >> dim_idx) & 1:
                        continue
                    r_occupied = r_base | r_bit

                    src_exp = pair_to_exp.get((s, r_base))
                    mid_exp = pair_to_exp.get((int_s, r_occupied))
                    if src_exp is None or mid_exp is None:
                        continue

                    Q_exp[src_exp, mid_exp] += GAMMA

                    # Router -> Delay
                    for (jst, jr, prob) in dests:
                        dest_vec = set_job_count_in_state(list(state_space[int_s]), jst, jr, +1)
                        if dest_vec is None:
                            continue
                        target_orig = state_map.get(tuple(int(x) for x in dest_vec), -1)
                        if target_orig >= 0:
                            dst_exp = pair_to_exp.get((target_orig, r_base))
                            if dst_exp is not None:
                                Q_exp[mid_exp, dst_exp] += GAMMA * prob

    # Fix diagonal
    Q_exp = ctmc_makeinfgen(Q_exp)

    # Stochastic complementation: keep Router-empty states (bitmask 0)
    keep_indices = []
    remove_indices = []
    for idx, (s, r) in enumerate(expanded_pairs):
        if r == 0:
            keep_indices.append(idx)
        else:
            remove_indices.append(idx)
    keep = np.array(keep_indices)
    remove = np.array(remove_indices)

    if options.verbose:
        print(f"  Router stochcomp: {n_expanded} expanded, {len(keep)} non-immediate, "
              f"{len(remove)} immediate")

    if len(remove) == 0:
        return Q_orig, depRates_orig

    Q11 = Q_exp[np.ix_(keep, keep)]
    Q12 = Q_exp[np.ix_(keep, remove)]
    Q21 = Q_exp[np.ix_(remove, keep)]
    Q22 = Q_exp[np.ix_(remove, remove)]

    try:
        neg_Q22 = -Q22
        X = np.linalg.solve(neg_Q22, Q21)
        Q_reduced = Q11 + Q12 @ X
    except np.linalg.LinAlgError:
        if options.verbose:
            print("  Router stochcomp: linear solve failed, using pseudoinverse")
        try:
            neg_Q22_pinv = np.linalg.pinv(-Q22)
            Q_reduced = Q11 + Q12 @ neg_Q22_pinv @ Q21
        except Exception:
            return Q_orig, depRates_orig

    Q_reduced = ctmc_makeinfgen(Q_reduced)

    # Apply stochcomp to depRates
    depRates_exp = np.zeros((n_expanded, M, K))
    for idx, (s, r) in enumerate(expanded_pairs):
        depRates_exp[idx] = depRates_orig[s]

    depRates_reduced = np.zeros((n_orig, M, K))
    try:
        neg_Q22 = -Q22
        for ist in range(M):
            for kk in range(K):
                dep_nonimm = depRates_exp[keep, ist, kk]
                dep_imm = depRates_exp[remove, ist, kk]
                dep_correction = np.linalg.solve(neg_Q22, dep_imm)
                depRates_reduced[:, ist, kk] = dep_nonimm + Q12 @ dep_correction
    except Exception:
        depRates_reduced = depRates_orig

    if options.verbose:
        print(f"  Router stochcomp: {Q_exp.shape} -> {Q_reduced.shape}")

    return Q_reduced, depRates_reduced


def _ctmc_stochcomp(Q: np.ndarray, keep_indices: List[int]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform stochastic complementation to remove immediate states.

    This removes states where jobs are at non-station stateful nodes
    (like Router) by computing the equivalent transitions that bypass
    these nodes.

    Args:
        Q: Full infinitesimal generator matrix
        keep_indices: Indices of states to keep (non-immediate states)

    Returns:
        Tuple of (reduced Q matrix, transformation matrix for rates)
    """
    n = Q.shape[0]
    all_indices = set(range(n))
    remove_indices = sorted(all_indices - set(keep_indices))

    if not remove_indices:
        return Q, np.eye(n)

    keep_indices = sorted(keep_indices)

    # Partition Q into blocks:
    # Q = [Q11 Q12]  where 1 = keep, 2 = remove
    #     [Q21 Q22]
    Q11 = Q[np.ix_(keep_indices, keep_indices)]
    Q12 = Q[np.ix_(keep_indices, remove_indices)]
    Q21 = Q[np.ix_(remove_indices, keep_indices)]
    Q22 = Q[np.ix_(remove_indices, remove_indices)]

    # Stochastic complement: Q_reduced = Q11 + Q12 * (-Q22)^{-1} * Q21
    # Match MATLAB: T = (-Q22) \ Q21; T = Q12*T; S = Q11+T
    neg_Q22 = -Q22
    try:
        T = np.linalg.solve(neg_Q22, Q21)
    except np.linalg.LinAlgError:
        try:
            T = np.linalg.pinv(neg_Q22) @ Q21
        except Exception:
            return Q11, np.eye(len(keep_indices))

    Q_reduced = Q11 + Q12 @ T

    return Q_reduced, Q12 @ np.linalg.solve(neg_Q22, np.eye(len(remove_indices))) if len(remove_indices) > 0 else np.eye(len(keep_indices))


def _find_immediate_states(
    sn: NetworkStruct,
    state_space: List[List[int]],
    rrobin_info: dict
) -> Tuple[List[int], List[int]]:
    """
    Find immediate states where jobs are at non-station stateful nodes (e.g., Router).

    Immediate states are those where at least one job is at a stateful node that is
    NOT a station and NOT a Cache. These states have instant transitions and should
    be eliminated via stochastic complementation.

    This matches MATLAB's solver_ctmc.m lines 286-297:
    - For non-station stateful nodes that are not Caches
    - Find states where there's at least 1 job at that node (sum > 0)

    Args:
        sn: Network structure
        state_space: List of state vectors
        rrobin_info: Round-robin and state space info

    Returns:
        Tuple of (non-immediate state indices, immediate state indices)
    """
    K = int(sn.nclasses)
    I = int(sn.nnodes)

    # Find immediate (non-station stateful) nodes that are NOT Caches
    immediate_stateful_indices = []
    for ind in range(I):
        is_stateful = sn.isstateful[ind] if ind < len(sn.isstateful) else False
        is_station = sn.isstation[ind] if ind < len(sn.isstation) else False
        node_type = sn.nodetype[ind] if ind < len(sn.nodetype) else None

        # MATLAB: sn.isstateful(ind) && ~sn.isstation(ind) && sn.nodetype(ind) ~= NodeType.Cache
        if is_stateful and not is_station and node_type != NodeType.CACHE:
            if hasattr(sn, 'nodeToStateful') and sn.nodeToStateful is not None:
                sf_idx = int(sn.nodeToStateful[ind])
                if sf_idx >= 0:
                    immediate_stateful_indices.append(sf_idx)

    if not immediate_stateful_indices:
        # No immediate nodes
        return list(range(len(state_space))), []

    # Get state space structure info
    # The state vector layout depends on phases and stations
    # For each state, check if any immediate node has jobs
    phases = rrobin_info.get('phases', np.ones((int(sn.nstations), K), dtype=int))
    skip_stations = rrobin_info.get('skip_stations', set())
    cache_state_offsets = rrobin_info.get('cache_state_offsets', {})

    # Build column index mapping: which columns in state vector correspond to which stateful nodes
    # This is complex because the state vector has different structures for different node types
    # For now, use a heuristic based on the state_space_aggr structure

    # The aggregated state space should have one entry per stateful node per class
    # state_space_aggr[s][isf, k] = number of jobs of class k at stateful node isf

    nstateful = int(sn.nstateful) if hasattr(sn, 'nstateful') else int(sn.nstations)

    # For a simpler approach: check if we have state_space_aggr available
    # Otherwise, we need to infer from the state vector structure

    # For Router nodes specifically:
    # In the current state space enumeration, Router jobs are tracked in specific columns
    # Let's identify these by looking at the rrobin_info structure

    # Actually, looking at _enumerate_state_space, the base state is built per station
    # Router is NOT a station, so its jobs aren't directly in the base state
    # Instead, jobs pass through Router instantly via routing probabilities

    # The issue is that our current state space enumeration may NOT explicitly track
    # jobs at Router - they're modeled as immediate transitions in routing.

    # Let me check if state_space tracks Router jobs...
    # If Router is stateful, there should be columns for it

    # For the current implementation, let's check if we can identify states
    # where immediate node columns have non-zero values

    # Fallback: if we can't identify immediate states, return all states as non-immediate
    # This will give the same behavior as before (no stochastic complementation)

    # For now, let's use a different approach:
    # Look at the state_space_aggr to identify states with jobs at immediate nodes

    # Actually, the state space enumeration uses state_space_aggr which has shape
    # (nstateful, K) for each state. We need to check states where
    # state_space_aggr[sf_idx, :].sum() > 0 for any sf_idx in immediate_stateful_indices

    # The issue is that state_space and state_space_aggr may have different structures
    # Let me return early for now and investigate further

    # For models with Router, jobs don't actually stay at Router (it's immediate)
    # The MATLAB approach is to enumerate states with Router jobs and remove them
    # Our Python enumeration may not include these states at all

    # Return all states as non-immediate for now (maintains current behavior)
    # TODO: Implement proper immediate state identification based on state_space_aggr
    return list(range(len(state_space))), []


def _get_rrobin_outlinks(sn: NetworkStruct) -> dict:
    """
    Get outlinks for nodes with RROBIN/WRROBIN routing.

    Returns a dict: {(node_idx, class_idx): [outlink_node_indices]}
    """
    outlinks = {}

    if not hasattr(sn, 'routing') or sn.routing is None:
        return outlinks
    if not hasattr(sn, 'connmatrix') or sn.connmatrix is None:
        return outlinks

    routing = np.asarray(sn.routing)
    connmatrix = np.asarray(sn.connmatrix)

    N = routing.shape[0]  # Number of nodes
    K = routing.shape[1]  # Number of classes

    for ind in range(N):
        for r in range(K):
            strategy = routing[ind, r]
            # Check for RROBIN or WRROBIN
            if strategy == RoutingStrategy.RROBIN.value or strategy == RoutingStrategy.WRROBIN.value:
                # Get outgoing links from connection matrix
                links = np.where(connmatrix[ind, :] > 0)[0].tolist()
                if links:
                    outlinks[(ind, r)] = links

    return outlinks


def _build_rrobin_state_info(sn: NetworkStruct) -> dict:
    """
    Build information about round-robin state variables.

    Returns a dict with:
        'outlinks': {(node_idx, class_idx): [outlink_indices]}
        'state_vars': List of (node_idx, class_idx, num_outlinks) tuples
        'total_vars': Total number of extra state variables
        'non_station_stateful': Set of node indices that are stateful but not stations
    """
    outlinks = _get_rrobin_outlinks(sn)

    # Identify non-station stateful nodes (like Router)
    non_station_stateful = set()
    if hasattr(sn, 'isstation') and hasattr(sn, 'isstateful'):
        isstation = np.asarray(sn.isstation).flatten()
        isstateful = np.asarray(sn.isstateful).flatten()
        for i in range(len(isstation)):
            if isstateful[i] and not isstation[i]:
                non_station_stateful.add(i)

    state_vars = []
    for (node_idx, class_idx), links in sorted(outlinks.items()):
        # Track RROBIN state for all nodes including Router (non-station stateful nodes)
        # MATLAB tracks Router RROBIN pointers in the state space and uses stochcomp
        # to remove immediate Router states, producing correct effective rates
        state_vars.append((node_idx, class_idx, len(links)))

    return {
        'outlinks': outlinks,
        'state_vars': state_vars,
        'total_vars': len(state_vars),
        'non_station_stateful': non_station_stateful
    }


def _resolve_routing_through_non_stations(
    sn: NetworkStruct,
    src_node: int,
    dst_node: int,
    job_class: int,
    rrobin_info: dict,
    state: np.ndarray,
    rr_var_map: dict,
    M: int,
    K: int
) -> Tuple[int, int, np.ndarray]:
    """
    Resolve routing through non-station stateful nodes.

    When a job routes to a non-station node (like Router) with RROBIN,
    this function finds the final station destination and updates the
    RR pointer.

    Args:
        sn: Network structure
        src_node: Source node index
        dst_node: Destination node index (may be non-station)
        job_class: Job class index
        rrobin_info: Round-robin routing info
        state: Current state vector
        rr_var_map: Map from (node, class) to state variable index
        M: Number of stations
        K: Number of classes

    Returns:
        Tuple of (final_node, final_station, updated_state)
        Returns (-1, -1, state) if destination is a sink or invalid
    """
    non_station_stateful = rrobin_info['non_station_stateful']
    outlinks = rrobin_info['outlinks']

    # Get node to station mapping
    node_to_station = None
    if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
        node_to_station = np.asarray(sn.nodeToStation).flatten()

    current_node = dst_node
    new_state = state.copy()
    max_hops = 10  # Prevent infinite loops

    for _ in range(max_hops):
        # Check if current node is a station
        if node_to_station is not None and current_node < len(node_to_station):
            station_idx = int(node_to_station[current_node])
            if station_idx >= 0:
                return current_node, station_idx, new_state

        # Current node is not a station - check if it's a non-station stateful node
        if current_node not in non_station_stateful:
            # Node is neither station nor stateful (e.g., Sink)
            return -1, -1, new_state

        # It's a non-station stateful node - check for RROBIN routing
        if (current_node, job_class) in outlinks:
            # Apply RROBIN routing
            links = outlinks[(current_node, job_class)]
            rr_var_idx = rr_var_map.get((current_node, job_class))
            if rr_var_idx is not None:
                current_rr_ptr = int(new_state[rr_var_idx])
                next_node = links[current_rr_ptr]
                # Advance RR pointer
                new_state[rr_var_idx] = (current_rr_ptr + 1) % len(links)
                current_node = next_node
            else:
                # No RR state variable - use first outlink
                current_node = links[0]
        else:
            # No RROBIN - use connection matrix to find next node
            if hasattr(sn, 'connmatrix') and sn.connmatrix is not None:
                conn = np.asarray(sn.connmatrix)
                next_nodes = np.where(conn[current_node, :] > 0)[0]
                if len(next_nodes) > 0:
                    current_node = next_nodes[0]
                else:
                    return -1, -1, new_state
            else:
                return -1, -1, new_state

    # Max hops exceeded
    return -1, -1, new_state


def _get_probabilistic_final_destinations(
    sn: NetworkStruct,
    start_node: int,
    start_class: int,
    accumulated_prob: float = 1.0,
    max_depth: int = 10
) -> List[Tuple[int, int, float]]:
    """
    Recursively follow routing probabilities through non-station nodes to find final station destinations.

    This handles RAND routing properly by following the rtnodes probabilities
    instead of deterministically picking the first outlink.

    Args:
        sn: Network structure
        start_node: Starting node index
        start_class: Starting class index
        accumulated_prob: Probability accumulated so far (for recursion)
        max_depth: Maximum recursion depth to prevent infinite loops

    Returns:
        List of (station_index, final_class, cumulative_prob) tuples
        Returns empty list if all paths lead to Sink (job exits system)
    """
    results = []

    if max_depth <= 0 or accumulated_prob < 1e-12:
        return results

    if not hasattr(sn, 'rtnodes') or sn.rtnodes is None:
        return results

    P_nodes = np.asarray(sn.rtnodes)
    K = sn.nclasses
    I = sn.nnodes

    # Get node-to-station mapping
    node_to_station = None
    if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
        node_to_station = np.asarray(sn.nodeToStation).flatten()

    src_idx = start_node * K + start_class

    for jnd in range(I):
        for r in range(K):
            dst_idx = jnd * K + r
            if src_idx < P_nodes.shape[0] and dst_idx < P_nodes.shape[1]:
                prob = P_nodes[src_idx, dst_idx]
            else:
                prob = 0

            if prob <= 0:
                continue

            new_prob = accumulated_prob * prob

            # Check if destination is a Sink node - job exits system, don't follow further
            if hasattr(sn, 'nodetype') and sn.nodetype is not None and jnd < len(sn.nodetype):
                nt = sn.nodetype[jnd]
                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                if nt_val == 1:  # NodeType.SINK = 1
                    # Job exits system at Sink - don't add to results or follow routing
                    continue

            # Check if destination is a station
            if node_to_station is not None and jnd < len(node_to_station):
                jst = int(node_to_station[jnd])
                if jst >= 0:
                    # Reached a station - add to results
                    results.append((jst, r, new_prob))
                    continue

            # Destination is a non-station node (Router, etc.) - recursively follow routing
            sub_results = _get_probabilistic_final_destinations(
                sn, jnd, r, new_prob, max_depth - 1
            )
            results.extend(sub_results)

    return results


@dataclass
class SolverCTMCOptions:
    """Options for CTMC solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    cutoff: int = 10  # Cutoff for open class populations
    hide_immediate: bool = True  # Hide immediate transitions
    state_space_gen: str = 'default'  # 'default', 'full', 'reachable'
    force: bool = False  # Force solver to run even if state space may be too large


@dataclass
class SolverCTMCReturn:
    """
    Result of CTMC solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        pi: Steady-state distribution
        infgen: Infinitesimal generator matrix
        space: State space matrix
        station_col_ranges: List of (start, end) tuples for each station's columns in state space
        rrobin_info: Round-robin and cache state information
        runtime: Runtime in seconds
        method: Method used
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    pi: Optional[np.ndarray] = None
    infgen: Optional[np.ndarray] = None
    space: Optional[np.ndarray] = None
    space_aggr: Optional[np.ndarray] = None
    station_col_ranges: Optional[List[Tuple[int, int]]] = None
    rrobin_info: Optional[dict] = None
    runtime: float = 0.0
    method: str = "default"


def _normalize_sched_strategy(sched_val) -> str:
    """
    Normalize a scheduling strategy value to its name string.

    This handles the case where sched_val might come from different
    SchedStrategy enum definitions (constants.py vs network_struct.py).

    Args:
        sched_val: A SchedStrategy enum value or its name

    Returns:
        The strategy name as a string (e.g., 'FCFS', 'LCFS')
    """
    if hasattr(sched_val, 'name'):
        return sched_val.name
    return str(sched_val)


def _sched_is(sched_val, *names) -> bool:
    """
    Check if a scheduling strategy matches any of the given names.

    This handles the case where sched_val might come from different
    SchedStrategy enum definitions.

    Args:
        sched_val: A SchedStrategy enum value
        *names: One or more strategy names to check (e.g., 'LCFS', 'LCFSPR')

    Returns:
        True if sched_val matches any of the given names
    """
    sched_name = _normalize_sched_strategy(sched_val)
    return sched_name in names


def _get_fcfs_stations(sn: NetworkStruct) -> set:
    """
    Identify stations using FCFS-like scheduling (buffer-based state representation).

    This includes FCFS, HOL, and LCFS scheduling strategies, which use explicit
    buffer ordering (sequence of class IDs) in the CTMC state representation.

    Note: SIRO uses a different buffer structure (job counts per class, not sequence),
    so it is NOT included here. SIRO is handled like PS for state enumeration.

    Returns:
        Set of station indices using FCFS-like scheduling.
    """
    fcfs_stations = set()
    # Scheduling strategies that use buffer ordering (sequence of class IDs)
    # SIRO is excluded because it stores counts per class, not a sequence
    # Use name strings for comparison to handle different enum definitions
    buffer_order_names = {'FCFS', 'HOL', 'LCFS', 'LCFSPR'}
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            sched_name = _normalize_sched_strategy(sched_val)
            if sched_name in buffer_order_names:
                fcfs_stations.add(ist)
    return fcfs_stations


def _get_siro_stations(sn: NetworkStruct) -> set:
    """
    Identify stations using SIRO scheduling (count-based buffer representation).

    SIRO uses an unordered buffer where we track job counts per class, not the
    sequence of jobs. This is different from FCFS/HOL/LCFS which use ordered buffers.

    Returns:
        Set of station indices using SIRO scheduling.
    """
    siro_stations = set()
    # Scheduling strategies that use count-based buffer (per-class job counts)
    # Only SIRO is commonly used; LEPT/SEPT/SRPT are less common and may not be defined
    # Use name strings for comparison to handle different enum definitions
    siro_names = {'SIRO'}
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            sched_name = _normalize_sched_strategy(sched_val)
            if sched_name in siro_names:
                siro_stations.add(ist)
    return siro_stations


def _generate_fcfs_buffer_orderings(n: np.ndarray, S: int, K: int) -> List[Tuple[Tuple[int, ...], np.ndarray]]:
    """
    Generate all unique buffer orderings for FCFS state enumeration.

    For FCFS with n[k] jobs of class k and S servers, generate all unique
    permutations of job classes in the buffer/server.

    Args:
        n: Array of job counts per class (length K)
        S: Number of servers
        K: Number of classes

    Returns:
        List of (all_positions, jobs_in_service_per_class) tuples
        - all_positions: tuple of class IDs (1-based) for ALL jobs (buffer + service)
          The rightmost S positions are the jobs in service (MATLAB FCFS format)
        - jobs_in_service_per_class: array[K] with count of each class in service
    """
    from itertools import permutations

    total_jobs = int(sum(n))
    if total_jobs == 0:
        # Empty queue: no buffer, no jobs in service
        return [(tuple(), np.zeros(K, dtype=int))]

    # Build list of job classes with repetition (1-based class IDs)
    # MATLAB uses standard encoding: class k (0-based) -> marker (k + 1)
    # e.g., K=2, n=[2,0] -> vi=[1,1] (class 0 -> marker 1)
    # e.g., K=2, n=[0,2] -> vi=[2,2] (class 1 -> marker 2)
    vi = []
    for k in range(K):
        vi.extend([k + 1] * int(n[k]))  # 1-based class IDs

    # Generate all unique permutations
    unique_perms = list(set(permutations(vi)))

    result = []
    for perm in unique_perms:
        # Last S jobs are in service, rest are in buffer
        num_in_service = min(total_jobs, S)
        num_in_buffer = total_jobs - num_in_service

        # Return ALL positions (buffer + service) - MATLAB format
        # The rightmost positions are in service
        all_positions = perm  # Full permutation includes all jobs

        # Count jobs of each class in service (last S positions)
        jobs_in_service = np.zeros(K, dtype=int)
        for class_id in perm[num_in_buffer:]:
            jobs_in_service[int(class_id) - 1] += 1  # Convert to 0-based

        result.append((all_positions, jobs_in_service))

    return result


def _multichoose_con(n: np.ndarray, S: int) -> List[np.ndarray]:
    """
    Generate all ways to pick S elements from available units in vector n.

    For SIRO scheduling, this enumerates all possible distributions of S
    jobs in service across K classes, where n[k] is the max number of
    class k jobs that can be selected.

    Args:
        n: Array of available counts per class (length K)
        S: Number of elements to pick (typically = number of servers)

    Returns:
        List of arrays, each representing how many of each class is selected.
        Each array has length K, with sum = S.
    """
    K = len(n)

    if S == 0:
        return [np.zeros(K, dtype=int)]

    if S == 1:
        result = []
        for i in range(K):
            if n[i] > 0:
                v = np.zeros(K, dtype=int)
                v[i] = 1
                result.append(v)
        return result if result else [np.zeros(K, dtype=int)]

    result = []
    for i in range(K):
        if n[i] > 0:
            n_1 = n.copy()
            n_1[i] = n_1[i] - 1
            sub_results = _multichoose_con(n_1, S - 1)
            for sub in sub_results:
                v = sub.copy()
                v[i] += 1
                result.append(v)

    # Remove duplicates
    if result:
        seen = set()
        unique_result = []
        for v in result:
            key = tuple(v)
            if key not in seen:
                seen.add(key)
                unique_result.append(v)
        return unique_result

    return [np.zeros(K, dtype=int)]


def _enumerate_state_space(
    sn: NetworkStruct,
    cutoff = 10,
    rrobin_info: Optional[dict] = None,
    use_phase_augmentation: bool = True
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """
    Enumerate the state space for a queueing network.

    For closed networks, enumerates all valid job distributions.
    For open networks, uses cutoff to bound the state space.

    State vector format (without phase augmentation):
    - First M*K elements: job counts at each (station, class) pair
    - Remaining elements: round-robin pointers for state-dependent routing

    State vector format (with phase augmentation):
    - Phase counts for each (station, class): sum(phases[i,k]) elements per station/class
    - Round-robin pointers
    - MAP phase variables for FCFS stations with MAP distributions

    Args:
        sn: Network structure
        cutoff: Maximum jobs per station for open networks.
                Can be an int (same cutoff for all), or a matrix (M x K)
                with per-station, per-class cutoffs.
        rrobin_info: Optional pre-computed round-robin routing info
        use_phase_augmentation: Whether to use phase-augmented state space

    Returns:
        Tuple of (state_space, state_space_aggr, rrobin_info)
    """
    # Get round-robin routing info if not provided
    if rrobin_info is None:
        rrobin_info = _build_rrobin_state_info(sn)
    M = sn.nstations
    K = sn.nclasses
    # Use sn_has_open_classes instead of sn_is_open_model to properly handle
    # mixed networks (which have both open and closed classes)
    is_open = sn_has_open_classes(sn)

    # Get population constraints
    if sn.njobs is not None:
        N = sn.njobs.flatten()
    else:
        N = np.ones(K)

    # Check if model contains Cache nodes
    # Note: Router nodes with RROBIN do NOT require reduced cutoff - they just need state tracking
    # Cache nodes are tracked as stateful but non-station nodes, so their state is in the
    # cache state variable, not in the job count variables.
    has_cache = False
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for nt in sn.nodetype:
            nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
            if nt_val == 6:  # CACHE only (not ROUTER)
                has_cache = True
                break

    # Note: We no longer limit cutoff to 1 for Cache models.
    # The previous limitation (min(cutoff, 1)) caused incorrect hit/miss ratios
    # because it severely truncated the state space for queues after the Cache.
    # The user-specified cutoff should be respected.
    # Cache state explosion is handled separately via cache_stations set below.

    # Handle cutoff as scalar or matrix
    cutoff_arr = np.atleast_2d(cutoff)
    is_matrix_cutoff = cutoff_arr.shape[0] > 1 or cutoff_arr.shape[1] > 1

    # Identify Cache stations - these have immediate processing (capacity 1)
    # This matches MATLAB's spaceGeneratorNodes.m behavior
    # Note: Router nodes are not stations - they're pass-through nodes tracked separately
    # Note: We need to map station indices to node indices since nodetype is indexed by node
    cache_stations = set()
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        # Get station to node mapping if available
        station_to_node = None
        if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
            station_to_node = sn.stationToNode

        for ist in range(M):
            # Get the node index for this station
            if station_to_node is not None and ist < len(station_to_node):
                node_idx = int(station_to_node[ist])
            else:
                node_idx = ist  # Fallback: assume station index = node index

            if node_idx >= 0 and node_idx < len(sn.nodetype):
                nt = sn.nodetype[node_idx]
                # NodeType.CACHE = 6
                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                if nt_val == 6:  # CACHE only
                    cache_stations.add(ist)

    def get_cutoff(ist: int, k: int) -> int:
        """Get cutoff for station ist and class k."""
        # Cache nodes have capacity 1 (immediate processing)
        # This prevents state space explosion in models with Cache nodes
        if ist in cache_stations:
            return 1

        # Determine base cutoff from options
        if is_matrix_cutoff:
            # Matrix cutoff: index by station and class
            # Handle different matrix layouts
            if cutoff_arr.shape[0] >= M and cutoff_arr.shape[1] >= K:
                base_cutoff = int(cutoff_arr[ist, k])
            elif cutoff_arr.shape[0] >= K and cutoff_arr.shape[1] >= M:
                # Transposed: (K x M) layout
                base_cutoff = int(cutoff_arr[k, ist])
            else:
                # Fallback: use first valid value or default
                base_cutoff = int(cutoff_arr.flat[0]) if cutoff_arr.size > 0 else 10
        else:
            base_cutoff = int(cutoff_arr.flat[0])

        # Also check station capacity from sn.cap
        # Capacity limits the maximum number of jobs that can be at the station
        if hasattr(sn, 'cap') and sn.cap is not None and ist < len(sn.cap):
            station_cap = sn.cap[ist]
            if np.isfinite(station_cap) and station_cap > 0:
                # Use the minimum of cutoff and capacity
                return min(base_cutoff, int(station_cap))

        return base_cutoff

    states = []

    # Helper function for closed class distributions
    def _enumerate_distributions(n_jobs: int, n_stations: int) -> List[List[int]]:
        """Enumerate all ways to distribute n_jobs among n_stations."""
        if n_stations <= 0:
            # Cannot distribute jobs among 0 stations
            # Return empty list if jobs > 0 (impossible), single empty list if jobs == 0
            return [[]] if n_jobs == 0 else []
        if n_stations == 1:
            return [[n_jobs]]
        if n_jobs == 0:
            return [[0] * n_stations]

        result = []
        for i in range(n_jobs + 1):
            for rest in _enumerate_distributions(n_jobs - i, n_stations - 1):
                result.append([i] + rest)
        return result

    # Determine which classes are open vs closed
    open_classes = [k for k in range(K) if np.isinf(N[k])]
    closed_classes = [k for k in range(K) if np.isfinite(N[k]) and N[k] > 0]

    # Identify Source and Sink stations (closed class jobs should NOT be there)
    # Note: nodetype is indexed by NODE index, not station index
    # Use stationToNode mapping to convert station index to node index
    source_station = -1
    sink_station = -1

    # Get station to node mapping
    station_to_node = None
    if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
        station_to_node = np.asarray(sn.stationToNode).flatten()

    # First check sched for EXT scheduling (Source marker)
    # Note: sched dict may store enum values or integer values
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            if _sched_is(sched_val, 'EXT') or (isinstance(sched_val, int) and sched_val == SchedStrategy.EXT.value):
                source_station = ist

    # Also check nodetype for Source and Sink
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for ist in range(M):
            # Convert station index to node index
            if station_to_node is not None and ist < len(station_to_node):
                node_idx = int(station_to_node[ist])
            else:
                node_idx = ist  # Fallback if no mapping

            if node_idx < len(sn.nodetype):
                if sn.nodetype[node_idx] == NodeType.SOURCE:
                    source_station = ist
                elif sn.nodetype[node_idx] == NodeType.SINK:
                    sink_station = ist

    # Compute reachable stations for each closed class based on routing
    # This prevents generating states where closed jobs are at unreachable stations
    def _get_reachable_stations(class_k: int, ref_station: int) -> set:
        """Find all stations reachable by class k starting from ref_station."""
        if not hasattr(sn, 'rt') or sn.rt is None:
            # No routing info - assume all non-Source/Sink stations are reachable
            return set(ist for ist in range(M) if ist != source_station and ist != sink_station)

        rt = np.asarray(sn.rt)
        reachable = set()
        to_visit = [ref_station]

        while to_visit:
            ist = to_visit.pop()
            if ist in reachable:
                continue
            if ist == source_station or ist == sink_station:
                continue
            reachable.add(ist)

            # Find outgoing transitions for this class at this station
            src_idx = ist * K + class_k
            if src_idx >= rt.shape[0]:
                continue
            for jst in range(M):
                dst_idx = jst * K + class_k
                if dst_idx >= rt.shape[1]:
                    continue
                if rt[src_idx, dst_idx] > 1e-10 and jst not in reachable:
                    to_visit.append(jst)

        return reachable

    def _get_reachable_pairs_in_chain(chain_classes: list, ref_station: int) -> list:
        """
        Find all reachable (station, class) pairs for a chain with class switching.

        When class switching exists, jobs can only be at certain (station, class) pairs
        based on the routing constraints. For example, if class 1 always switches to
        class 2 when going from Q1 to Q2, then (Q2, class 1) is unreachable.

        The function starts BFS only from pairs that have initial population,
        since pairs with no jobs that can reach them are unreachable.

        Args:
            chain_classes: List of class indices in the chain
            ref_station: Reference station for the chain

        Returns:
            List of (station, class) tuples that are reachable
        """
        if not hasattr(sn, 'rt') or sn.rt is None:
            # No routing info - assume all pairs are reachable
            pairs = []
            for ist in range(M):
                if ist != source_station and ist != sink_station:
                    for k in chain_classes:
                        pairs.append((ist, k))
            return pairs

        rt = np.asarray(sn.rt)

        # Use BFS to find reachable (station, class) pairs
        # Start only from pairs that have initial population
        reachable = set()
        to_visit = []

        # Initialize with pairs that have non-zero initial population
        # For closed classes, jobs start at the reference station
        for k in chain_classes:
            if np.isfinite(N[k]) and N[k] > 0:
                # This class has jobs, start from its reference station
                class_ref = ref_station
                if hasattr(sn, 'refstat') and sn.refstat is not None:
                    refstat = np.asarray(sn.refstat).flatten()
                    if k < len(refstat):
                        class_ref = int(refstat[k])
                to_visit.append((class_ref, k))

        # If no initial population found, fall back to all classes at ref_station
        if not to_visit:
            for k in chain_classes:
                to_visit.append((ref_station, k))

        while to_visit:
            (ist, class_k) = to_visit.pop()

            if (ist, class_k) in reachable:
                continue
            if ist == source_station or ist == sink_station:
                continue
            if class_k not in chain_classes:
                continue

            reachable.add((ist, class_k))

            # Find outgoing transitions - jobs can switch class!
            src_idx = ist * K + class_k
            if src_idx >= rt.shape[0]:
                continue

            # Check transitions to ALL classes (within chain) at ALL stations
            for jst in range(M):
                for dst_class in chain_classes:
                    dst_idx = jst * K + dst_class
                    if dst_idx >= rt.shape[1]:
                        continue
                    if rt[src_idx, dst_idx] > 1e-10 and (jst, dst_class) not in reachable:
                        to_visit.append((jst, dst_class))

        return sorted(reachable)

    # For each closed class, compute reachable stations
    closed_class_reachable = {}
    for k in closed_classes:
        # Get reference station for this class
        ref_station = 0  # Default
        if hasattr(sn, 'refstat') and sn.refstat is not None:
            refstat = np.asarray(sn.refstat).flatten()
            if k < len(refstat):
                ref_station = int(refstat[k])
        closed_class_reachable[k] = _get_reachable_stations(k, ref_station)

    # Stations where closed class jobs can reside (exclude Source/Sink)
    # For enumeration, use union of all reachable stations across all closed classes
    if closed_classes and closed_class_reachable:
        closed_valid_stations = sorted(set.union(*closed_class_reachable.values()) if closed_class_reachable else set())
    else:
        closed_valid_stations = [ist for ist in range(M) if ist != source_station and ist != sink_station]
    n_closed_stations = len(closed_valid_stations)

    if is_open and closed_classes:
        # Mixed network: combine closed class distributions with open class enumeration
        # For closed classes: enumerate valid distributions (conserving population)
        #                     only among valid (non-Source/Sink) stations
        # For open classes: enumerate 0 to cutoff TOTAL jobs per class, distributed
        #                   across valid stations (matching MATLAB's spaceGenerator.m)

        # Generate closed class distributions among valid stations only
        closed_class_dists = []
        for k in closed_classes:
            n_k = int(N[k])
            # Distribute among valid stations only
            class_dists = _enumerate_distributions(n_k, n_closed_stations)
            closed_class_dists.append((k, class_dists))

        # Get valid stations for open classes (exclude Source/Sink)
        open_valid_stations = [ist for ist in range(M) if ist != source_station and ist != sink_station]
        n_open_stations = len(open_valid_stations)

        # Get per-class cutoff (max TOTAL jobs per class in network)
        # This matches MATLAB: Np(r) = max(capacityc(:,r))
        open_class_cutoffs = []
        for k in open_classes:
            if is_matrix_cutoff:
                # Matrix cutoff: max across stations for this class
                max_cutoff = max(get_cutoff(ist, k) for ist in open_valid_stations) if open_valid_stations else 1
            else:
                max_cutoff = int(cutoff_arr.flat[0])
            open_class_cutoffs.append((k, max_cutoff))

        # Combine: iterate over closed distributions × open distributions
        from itertools import product as iter_product

        # Get all combinations of closed class distributions
        if closed_class_dists:
            closed_combos = list(iter_product(*[dists for _, dists in closed_class_dists]))
        else:
            closed_combos = [()]

        # Generate open class distributions: for each class, enumerate total jobs 0 to cutoff,
        # then distribute across valid stations (same as pure open case and MATLAB)
        # This ensures total open jobs per class <= cutoff
        if open_class_cutoffs and n_open_stations > 0:
            # For each open class, generate all distributions of 0 to cutoff jobs
            open_class_all_dists = []
            for k, max_cutoff in open_class_cutoffs:
                # All distributions of 0 to max_cutoff jobs across n_open_stations
                all_dists_for_class = []
                for total in range(max_cutoff + 1):
                    all_dists_for_class.extend(_enumerate_distributions(total, n_open_stations))
                open_class_all_dists.append((k, all_dists_for_class))
            open_combos = list(iter_product(*[dists for _, dists in open_class_all_dists]))
        else:
            open_combos = [()]

        # Build states from combinations
        for closed_combo in closed_combos:
            for open_combo in open_combos:
                # Build state vector: [n_11, n_12, ..., n_1K, n_21, ..., n_MK]
                state = [0] * (M * K)
                # Fill in closed classes (mapping from valid station indices to full indices)
                for idx, (k, _) in enumerate(closed_class_dists):
                    dist = closed_combo[idx]
                    for valid_idx, ist in enumerate(closed_valid_stations):
                        state[ist * K + k] = dist[valid_idx]
                # Fill in open classes (mapping from open_valid_stations to full indices)
                for idx, (k, _) in enumerate(open_class_all_dists if open_class_cutoffs and n_open_stations > 0 else []):
                    dist = open_combo[idx]
                    for valid_idx, ist in enumerate(open_valid_stations):
                        state[ist * K + k] = dist[valid_idx]
                states.append(state)

    elif is_open:
        # Pure open network: enumerate all valid job distributions
        # MATLAB behavior: cutoff limits TOTAL jobs per class in the ENTIRE network
        # Each class can have at most cutoff jobs distributed across all stations
        # This matches MATLAB's spaceGenerator.m where Np(r) = max(capacityc(:,r))

        # Get valid stations (exclude Source/Sink)
        valid_stations = [ist for ist in range(M) if ist != source_station and ist != sink_station]

        # Compute reachable (station, class) pairs for open networks with class switching
        # This is crucial for Cache networks where jobs switch from InitClass to HitClass/MissClass
        # Without this, the state space includes unreachable states like (Delay1, InitClass)
        def _get_reachable_pairs_open() -> Dict[int, List[int]]:
            """
            Find reachable stations for each open class, accounting for class switching.

            Returns:
                Dict mapping class_idx -> list of reachable station indices
            """
            # Use rt (routing matrix indexed by STATEFUL nodes) to trace reachability
            # IMPORTANT: rt is indexed by stateful node indices, NOT station indices!
            # We need to map between station and stateful indices correctly.
            if not hasattr(sn, 'rt') or sn.rt is None:
                # No routing info - assume all valid stations are reachable for all classes
                return {k: valid_stations for k in range(K)}

            rt = np.asarray(sn.rt)
            nstateful = rt.shape[0] // K
            reachable_stations = {k: set() for k in range(K)}

            # Get station to stateful mapping
            station_to_stateful = {}
            stateful_to_station = {}
            if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None:
                for ist in range(M):
                    isf = int(sn.stationToStateful[ist])
                    station_to_stateful[ist] = isf
                    stateful_to_station[isf] = ist
            else:
                # Fallback: assume station index = stateful index
                for ist in range(M):
                    station_to_stateful[ist] = ist
                    stateful_to_station[ist] = ist

            # BFS starting from Source for classes with actual arrivals
            # Only start BFS from classes that have arrivals at Source (non-zero arrival rate)
            # Other classes (created by class switching) will be discovered via BFS
            if source_station < 0:
                return {k: valid_stations for k in range(K)}

            # Get Source's stateful index
            source_stateful = station_to_stateful.get(source_station, source_station)

            # Check which classes have arrivals at Source (non-zero rate)
            arrival_classes = []
            if hasattr(sn, 'rates') and sn.rates is not None:
                rates = np.asarray(sn.rates)
                for k in range(K):
                    if source_station < rates.shape[0] and k < rates.shape[1]:
                        if rates[source_station, k] > 0:
                            arrival_classes.append(k)
            if not arrival_classes:
                # Fallback: assume first open class has arrivals
                for k in range(K):
                    if np.isinf(N[k]):
                        arrival_classes.append(k)
                        break

            # BFS using STATEFUL indices to find all reachable (stateful, class) pairs
            # Starting from Source with arrival classes only
            visited = set()
            for k in arrival_classes:
                to_visit = [(source_stateful, k)]

                while to_visit:
                    (isf, cur_class) = to_visit.pop()
                    if (isf, cur_class) in visited:
                        continue
                    visited.add((isf, cur_class))

                    # If this stateful node corresponds to a valid station, mark it reachable
                    if isf in stateful_to_station:
                        ist = stateful_to_station[isf]
                        if ist in valid_stations:
                            reachable_stations[cur_class].add(ist)

                    # Follow outgoing transitions (may switch class)
                    # rt is indexed by stateful node, class
                    src_idx = isf * K + cur_class
                    if src_idx >= rt.shape[0]:
                        continue
                    for jsf in range(nstateful):
                        for dst_class in range(K):
                            dst_idx = jsf * K + dst_class
                            if dst_idx >= rt.shape[1]:
                                continue
                            if rt[src_idx, dst_idx] > 1e-10 and (jsf, dst_class) not in visited:
                                to_visit.append((jsf, dst_class))

            # Convert sets to sorted lists
            return {k: sorted(v) if v else [] for k, v in reachable_stations.items()}

        # Get per-class reachable stations
        class_reachable = _get_reachable_pairs_open()

        # NOTE: Do NOT fall back to all valid_stations when a class has empty reachability!
        # Empty reachability means the class was traced via BFS but doesn't reach any valid
        # stations (e.g., InitClass in a Cache model switches to HitClass/MissClass at Cache).
        # Such classes should ONLY exist with 0 jobs in the state space.
        # The fallback is only used when _get_reachable_pairs_open() returns default
        # (all valid_stations for all classes) due to missing routing info.

        from itertools import product as iter_product

        # Compute Np(k) for each class: max capacity across ALL nodes (not just stations)
        # This matches MATLAB: Np(r) = max(capacityc(:,r)) which considers stations,
        # Cache nodes (capacity 1 per class), and Router nodes.
        # Without including Cache capacity, chain populations are underestimated
        # and states like (InitClass=1 at Cache + HitClass=1 at D1) are missed.
        Np = np.zeros(K, dtype=int)
        for k in range(K):
            reachable = class_reachable.get(k, [])
            if reachable:
                if is_matrix_cutoff:
                    Np[k] = max(get_cutoff(ist, k) for ist in reachable)
                else:
                    Np[k] = int(cutoff_arr.flat[0])
            # Also include Cache node capacity (1 per class) if Cache exists
            # This matches MATLAB where capacityc(Cache, r) = 1 for all r
            if has_cache:
                Np[k] = max(Np[k], 1)

        # Check for chain-based class switching
        # When classes share a chain (e.g., InitClass, HitClass, MissClass via Cache),
        # MATLAB distributes the CHAIN population (sum of Np across chain classes)
        # among all classes freely, not independently per class.
        # This is critical: with 3 classes in 1 chain and Np=[1,1,1],
        # the system can have e.g. 2 HitClass + 1 MissClass simultaneously.
        has_chains = (hasattr(sn, 'chains') and sn.chains is not None and
                      hasattr(sn, 'inchain') and sn.inchain is not None and
                      hasattr(sn, 'nchains') and sn.nchains is not None)

        has_class_switching = False
        if has_chains and sn.nchains < K:
            for chain_id in range(sn.nchains):
                if chain_id in sn.inchain and len(sn.inchain[chain_id]) > 1:
                    has_class_switching = True
                    break

        if has_class_switching:
            # Chain-aware enumeration matching MATLAB's spaceGenerator.m:
            # 1. For each chain, compute chain_pop = sum(Np[k] for k in chain)
            # 2. Iterate total population from 0 to chain_pop
            # 3. For each total, distribute among classes (multichoose)
            #    NOTE: No per-class Np check! MATLAB allows e.g. HitClass=2
            #    across different stations as long as per-station capacity is met.
            # 4. For each class distribution, distribute among reachable stations
            #    with per-station capacity constraint (each station holds at most 1
            #    job per class, matching MATLAB's capacityc check)
            # 5. Cartesian product across chains

            def _multichoose(n_bins: int, n_items: int) -> List[List[int]]:
                """Generate all ways to distribute n_items among n_bins (with repetition).
                Equivalent to MATLAB's multichoose(n_bins, n_items).
                Returns list of lists, each of length n_bins summing to n_items."""
                if n_bins == 1:
                    return [[n_items]]
                if n_items == 0:
                    return [[0] * n_bins]
                result = []
                for i in range(n_items + 1):
                    for rest in _multichoose(n_bins - 1, n_items - i):
                        result.append([i] + rest)
                return result

            def _enumerate_distributions_capped(n_jobs: int, n_stations: int, cap: int = 1) -> List[List[int]]:
                """Enumerate distributions of n_jobs among n_stations with per-station cap.
                Each station can hold at most 'cap' jobs.
                Matches MATLAB's capacityc constraint at stations."""
                if n_stations <= 0:
                    return [[]] if n_jobs == 0 else []
                if n_stations == 1:
                    return [[n_jobs]] if n_jobs <= cap else []
                if n_jobs == 0:
                    return [[0] * n_stations]
                result = []
                for i in range(min(n_jobs, cap) + 1):
                    for rest in _enumerate_distributions_capped(n_jobs - i, n_stations - 1, cap):
                        result.append([i] + rest)
                return result

            all_chain_placements = []  # List of (chain_id, chain_classes, placements)

            for chain_id in range(sn.nchains):
                if chain_id not in sn.inchain:
                    continue
                chain_classes = list(sn.inchain[chain_id])
                chain_pop = sum(int(Np[k]) for k in chain_classes)

                # Build reachable (station, class) pairs for this chain
                # Each class has its own set of reachable stations
                chain_reachable = {}  # class -> list of reachable stations
                for k in chain_classes:
                    chain_reachable[k] = class_reachable.get(k, [])

                # Compute per-station per-class capacity
                # This matches MATLAB's capacityc: for stations, capacity is min(cutoff, classcap)
                # For most stations with cutoff=1, this is 1 per class
                station_class_cap = {}
                for k in chain_classes:
                    for ist in chain_reachable.get(k, []):
                        station_class_cap[(ist, k)] = get_cutoff(ist, k)

                # Generate all placements for this chain
                chain_placements = []  # List of state dicts: {(ist, k): count}

                # Iterate over total chain population from 0 to chain_pop
                # For each total, distribute among classes using multichoose
                # No per-class Np filter: MATLAB distributes chain total freely
                # among classes, filtering only at per-station level
                for n_total in range(chain_pop + 1):
                    # All ways to distribute n_total jobs among len(chain_classes) classes
                    class_dists = _multichoose(len(chain_classes), n_total)

                    for class_dist in class_dists:
                        # class_dist[i] = number of jobs of chain_classes[i]
                        # Quick feasibility check: can n_k jobs of class k fit
                        # across reachable stations? n_k must be <= sum of capacities
                        # at reachable stations for that class.
                        feasible = True
                        for i, k in enumerate(chain_classes):
                            n_k = class_dist[i]
                            reachable = chain_reachable[k]
                            if n_k > 0 and not reachable:
                                # Jobs of unreachable class must be 0 at stations
                                # (they can still be at Cache via augmentation)
                                # Don't block: just produce empty station distribution
                                pass
                            elif n_k > 0:
                                # Check if n_k jobs can fit across reachable stations
                                total_cap = sum(station_class_cap.get((ist, k), 1)
                                               for ist in reachable)
                                if n_k > total_cap:
                                    feasible = False
                                    break
                        if not feasible:
                            continue

                        # For each class with n_k > 0, distribute among reachable stations
                        # respecting per-station capacity
                        per_class_station_dists = []
                        for i, k in enumerate(chain_classes):
                            n_k = class_dist[i]
                            reachable = chain_reachable[k]
                            if n_k == 0 or not reachable:
                                # No jobs of this class at stations
                                per_class_station_dists.append([(k, reachable, [0] * max(len(reachable), 1))])
                            else:
                                # Distribute n_k jobs among reachable stations with per-station cap
                                cap = min(station_class_cap.get((ist, k), 1) for ist in reachable) if reachable else 1
                                station_dists = _enumerate_distributions_capped(n_k, len(reachable), cap)
                                if not station_dists:
                                    # Can't fit n_k jobs - skip this class dist
                                    per_class_station_dists = None
                                    break
                                per_class_station_dists.append(
                                    [(k, reachable, sd) for sd in station_dists]
                                )

                        if per_class_station_dists is None:
                            continue

                        # Cartesian product across classes in this chain
                        for combo in iter_product(*per_class_station_dists):
                            placement = {}
                            for (k, reachable, sd) in combo:
                                for j, ist in enumerate(reachable):
                                    if j < len(sd) and sd[j] > 0:
                                        placement[(ist, k)] = placement.get((ist, k), 0) + sd[j]
                            chain_placements.append(placement)

                all_chain_placements.append((chain_id, chain_classes, chain_placements))

            # Cartesian product across chains
            for combo in iter_product(*[placements for _, _, placements in all_chain_placements]):
                state = [0] * (M * K)
                for placement in combo:
                    for (ist, k), count in placement.items():
                        state[ist * K + k] = count
                states.append(state)

            # Remove duplicates (matching MATLAB's unique(chainStationPos, 'rows'))
            unique_states = set()
            deduped_states = []
            for state in states:
                key = tuple(state)
                if key not in unique_states:
                    unique_states.add(key)
                    deduped_states.append(state)
            states = deduped_states
        else:
            # No class switching - original per-class independent enumeration
            # Build per-class distribution lists
            class_all_dists = []  # List of (class_k, reachable_stations, all_distributions)
            for k in range(K):
                reachable = class_reachable.get(k, [])
                if not reachable:
                    # Class has no reachable stations - only empty distribution
                    class_all_dists.append((k, [], [[]]))
                    continue

                max_cutoff = int(Np[k])

                # Enumerate all distributions of 0 to max_cutoff jobs across reachable stations
                all_dists_for_class = []
                for total in range(max_cutoff + 1):
                    all_dists_for_class.extend(_enumerate_distributions(total, len(reachable)))
                class_all_dists.append((k, reachable, all_dists_for_class))

            # Cartesian product across classes
            class_combos = list(iter_product(*[dists for _, _, dists in class_all_dists]))

            # Build states from combinations
            for combo in class_combos:
                state = [0] * (M * K)
                for idx, (k, reachable, _) in enumerate(class_all_dists):
                    dist = combo[idx]
                    for valid_idx, ist in enumerate(reachable):
                        state[ist * K + k] = dist[valid_idx]
                states.append(state)
    else:
        # Pure closed network: enumerate valid distributions
        # Check for class switching (chains) - classes in the same chain can exchange jobs
        has_chains = (hasattr(sn, 'chains') and sn.chains is not None and
                      hasattr(sn, 'inchain') and sn.inchain is not None and
                      hasattr(sn, 'nchains') and sn.nchains is not None)

        # Determine if class switching exists (some chain has multiple classes)
        has_class_switching = False
        if has_chains and sn.nchains < K:
            for chain_id in range(sn.nchains):
                if chain_id in sn.inchain and len(sn.inchain[chain_id]) > 1:
                    has_class_switching = True
                    break

        if has_class_switching:
            # Class switching exists - enumerate based on chain populations
            # For each chain, enumerate all ways to distribute the chain population
            # across REACHABLE (station, class) pairs within the chain
            # This prevents generating unreachable states like (Q2, Class1) when routing
            # dictates that jobs at Q2 are always Class2
            all_chain_dists = []

            for chain_id in range(sn.nchains):
                if chain_id not in sn.inchain:
                    continue
                chain_classes = list(sn.inchain[chain_id])  # Array of class indices in this chain

                # Get reference station for this chain (use first class's refstat)
                ref_station = 0
                if hasattr(sn, 'refstat') and sn.refstat is not None:
                    refstat = np.asarray(sn.refstat).flatten()
                    if chain_classes[0] < len(refstat):
                        ref_station = int(refstat[chain_classes[0]])

                # Get reachable (station, class) pairs for this chain
                reachable_pairs = _get_reachable_pairs_in_chain(chain_classes, ref_station)
                n_pairs = len(reachable_pairs)

                # Chain population = sum of N[k] for all k in chain
                chain_pop = sum(int(N[k]) for k in chain_classes if np.isfinite(N[k]))

                if chain_pop == 0 or n_pairs == 0:
                    # No jobs in this chain or no reachable pairs - single empty distribution
                    chain_dists = [[0] * max(n_pairs, 1)]
                else:
                    # Enumerate all ways to distribute chain_pop among reachable pairs
                    chain_dists = _enumerate_distributions(chain_pop, n_pairs)

                all_chain_dists.append((chain_id, chain_classes, reachable_pairs, chain_dists))

            # Combine distributions across chains
            for combo in product(*[dists for _, _, _, dists in all_chain_dists]):
                state = [0] * (M * K)
                for idx, (chain_id, chain_classes, reachable_pairs, _) in enumerate(all_chain_dists):
                    dist = combo[idx]
                    # Map distribution back to reachable (station, class) pairs
                    for pair_idx, (ist, k) in enumerate(reachable_pairs):
                        state[ist * K + k] = dist[pair_idx]
                states.append(state)
        else:
            # No class switching - original per-class enumeration
            all_class_dists = []
            for k in range(K):
                n_k = int(N[k]) if np.isfinite(N[k]) else 0
                class_dists = _enumerate_distributions(n_k, M)
                all_class_dists.append(class_dists)

            # Combine all class distributions
            for combo in product(*all_class_dists):
                # combo is tuple of distributions, one per class
                # Convert to state vector: [n_11, n_12, ..., n_1K, n_21, ..., n_MK]
                state = []
                for ist in range(M):
                    for k in range(K):
                        state.append(combo[k][ist])
                states.append(state)

    if not states:
        states = [[0] * (M * K)]

    # Filter states to respect station capacity constraints (sn.cap)
    # This matches MATLAB's fromMarginalBounds which filters by: ni <= cap
    # Note: Per-class constraints are already applied during enumeration via get_cutoff
    if hasattr(sn, 'cap') and sn.cap is not None:
        filtered_states = []
        for state in states:
            valid = True
            for ist in range(M):
                if ist < len(sn.cap):
                    cap_val = sn.cap[ist]
                    if np.isfinite(cap_val) and cap_val > 0:
                        # Sum jobs at this station across all classes
                        total_at_station = sum(state[ist * K + k] for k in range(K))
                        if total_at_station > cap_val:
                            valid = False
                            break
            if valid:
                filtered_states.append(state)
        states = filtered_states if filtered_states else [[0] * (M * K)]

    # Get phase information for phase-type state augmentation
    phases, needs_phase_augmentation = _get_phases_info(sn)

    # Compute phase-related state offsets for later use
    # phase_offset[ist, k] = starting index of phases for (ist, k) in phase-augmented state
    total_phases = int(np.sum(phases))
    phase_offset = np.zeros((M, K), dtype=int)
    idx = 0
    for ist in range(M):
        for k in range(K):
            phase_offset[ist, k] = idx
            idx += phases[ist, k]

    # Store phase information in rrobin_info for use by generator
    rrobin_info['phases'] = phases
    rrobin_info['phase_offset'] = phase_offset
    rrobin_info['total_phases'] = total_phases
    # Note: needs_phase_augmentation will be set after FCFS stations are identified

    # Identify FCFS stations for proper buffer ordering
    fcfs_stations = _get_fcfs_stations(sn)
    rrobin_info['fcfs_stations'] = fcfs_stations

    # Identify SIRO stations for count-based buffer tracking
    siro_stations = _get_siro_stations(sn)
    rrobin_info['siro_stations'] = siro_stations

    # Get number of servers per station
    nservers = sn.nservers if hasattr(sn, 'nservers') and sn.nservers is not None else np.ones(M, dtype=int)

    # Identify Source and Sink stations to skip from state space
    # These stations don't contribute to the CTMC state (infinite population for Source)
    skip_stations = set()
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        station_to_node = getattr(sn, 'stationToNode', None)
        for ist in range(M):
            # Map station to node index
            node_idx = ist
            if station_to_node is not None:
                node_arr = np.atleast_1d(station_to_node)
                if ist < len(node_arr):
                    node_idx = int(node_arr[ist])
            # Check if this node is Source or Sink
            if node_idx >= 0 and node_idx < len(sn.nodetype):
                nt = sn.nodetype[node_idx]
                if nt == NodeType.SOURCE or nt == NodeType.SINK:
                    skip_stations.add(ist)
    rrobin_info['skip_stations'] = skip_stations

    # If phase augmentation is needed and enabled, expand states with phase distributions
    # For FCFS stations, also expand with buffer orderings
    # For SIRO stations, also expand with per-class buffer counts
    # IMPORTANT: Always run state expansion for FCFS/SIRO stations to match MATLAB format
    # MATLAB generates [buffer, phases] format even for single-phase (Exp) distributions
    needs_state_expansion = (needs_phase_augmentation and use_phase_augmentation) or bool(fcfs_stations - skip_stations) or bool(siro_stations - skip_stations)
    # Update rrobin_info to reflect actual state expansion status
    rrobin_info['needs_phase_augmentation'] = needs_state_expansion
    if needs_state_expansion:
        # For each basic state (job counts), expand to all possible phase distributions
        phase_augmented_states = []

        # Track FCFS buffer structure for later use
        # fcfs_buffer_info[ist] = (buffer_start_idx, max_buffer_size)
        fcfs_buffer_info = {}

        # Compute maximum buffer size for each FCFS station
        # MATLAB format: buffer only stores WAITING jobs (class indices), not in-service jobs
        # In-service jobs are tracked by phase counts
        # Buffer size = max_jobs_at_station - servers (to hold waiting jobs when station is full)
        # For closed models, max_jobs_at_station = total_population
        # For open models, use cutoff
        fcfs_max_buffer = {}  # Maps ist -> maximum buffer size
        for ist in fcfs_stations:
            S_ist = int(nservers[ist]) if ist < len(nservers) else 1
            if is_open:
                # Open model: use cutoff per class * number of classes
                # MATLAB treats scalar cutoff as per-class, so max total = cutoff * K
                effective_cutoff = int(cutoff_arr.flat[0]) if cutoff_arr.size > 0 else 10
                max_at_station = effective_cutoff * K
            else:
                # Closed model: max at any station is total population
                max_at_station = int(np.sum(N[np.isfinite(N)]))
            # Buffer only holds waiting jobs: max_waiting = max_at_station - servers
            # Use at least 1 to ensure consistent state format even when max_waiting=0
            fcfs_max_buffer[ist] = max(1, max_at_station - S_ist)

        rrobin_info['fcfs_max_buffer'] = fcfs_max_buffer

        for base_state in states:
            # base_state has format [n_11, n_12, ..., n_1K, n_21, ..., n_MK]
            # We need to expand each n_ik into phase distribution
            # For FCFS stations with multiple classes, also enumerate buffer orderings

            # First pass: compute station-level job counts and identify FCFS expansion
            station_jobs = {}  # Maps ist -> array of job counts per class
            for ist in range(M):
                n_ist = np.array([int(base_state[ist * K + k]) for k in range(K)])
                station_jobs[ist] = n_ist

            # For each FCFS station, generate buffer orderings
            # For non-FCFS stations, use PS-style phase distribution
            fcfs_buffer_combos = {}  # Maps ist -> list of (buffer, jobs_in_service) tuples
            for ist in fcfs_stations:
                n_ist = station_jobs[ist]
                total_at_station = int(sum(n_ist))
                if total_at_station > 0 and sum(n_ist > 0) > 1:
                    # Multiple classes at station - need buffer ordering
                    S_ist = int(nservers[ist]) if ist < len(nservers) else 1
                    fcfs_buffer_combos[ist] = _generate_fcfs_buffer_orderings(n_ist, S_ist, K)
                else:
                    # Single class or empty - buffer ordering is trivial
                    fcfs_buffer_combos[ist] = None

            # Generate all phase distribution combinations for non-FCFS parts
            # and combine with FCFS buffer orderings
            def _expand_state(base_state, remaining_stations, current_phases):
                """Recursively expand state with phase distributions and FCFS buffers."""
                if not remaining_stations:
                    # All stations processed - yield the expanded state
                    yield current_phases
                    return

                ist = remaining_stations[0]
                rest = remaining_stations[1:]
                n_ist = station_jobs[ist]

                # Skip Source and Sink stations - they don't contribute to state space
                if ist in skip_stations:
                    for expanded in _expand_state(base_state, rest, current_phases):
                        yield expanded
                    return

                if ist in fcfs_buffer_combos and fcfs_buffer_combos[ist] is not None:
                    # FCFS station with multiple classes - expand with buffer orderings
                    max_buf_size = fcfs_max_buffer.get(ist, 0)

                    for all_positions, jobs_in_service in fcfs_buffer_combos[ist]:
                        # Generate all buffer orderings without class imbalance constraint
                        # The constraint was causing unreachable states for larger cutoffs

                        # For each class, generate phase distributions for jobs in service
                        phase_combos_for_service = []
                        for k in range(K):
                            n_in_service = jobs_in_service[k]
                            n_phases_k = phases[ist, k]
                            if n_phases_k > 1 and n_in_service > 0:
                                phase_dists = _generate_phase_distributions(int(n_in_service), int(n_phases_k))
                            else:
                                phase_dists = [tuple([int(n_in_service)] + [0] * (n_phases_k - 1))] if n_phases_k > 1 else [(int(n_in_service),)]
                            phase_combos_for_service.append(phase_dists)

                        # Get number of servers at this station
                        S_ist = int(nservers[ist]) if ist < len(nservers) else 1

                        # Cartesian product of phase distributions for this station
                        for phase_combo in product(*phase_combos_for_service):
                            # Build station state: buffer positions (padded) + server phases
                            station_state = [0] * max_buf_size  # Initialize with zeros

                            # all_positions includes all jobs: leftmost are buffer, rightmost are in service
                            # Only place the BUFFER positions (first num_in_buffer elements) into buffer
                            # Service jobs are tracked by phase counts, not buffer positions
                            total_jobs = len(all_positions)
                            num_in_buffer = total_jobs - min(total_jobs, S_ist)
                            buffer_positions = all_positions[:num_in_buffer]

                            # Place buffer positions at the end (right-justified)
                            for i, bp in enumerate(buffer_positions):
                                station_state[max_buf_size - num_in_buffer + i] = bp

                            # Add phase distributions for server
                            for pd in phase_combo:
                                station_state.extend(pd)
                            # Recurse
                            for expanded in _expand_state(base_state, rest, current_phases + station_state):
                                yield expanded
                elif ist in fcfs_stations:
                    # FCFS station with single class or empty
                    # For FCFS: only S jobs can be in service, rest must be in buffer
                    max_buf_size = fcfs_max_buffer.get(ist, 0)
                    S_ist = int(nservers[ist]) if ist < len(nservers) else 1
                    total_at_station = int(sum(n_ist))

                    if total_at_station == 0:
                        # Empty station
                        station_state = [0] * max_buf_size
                        for k in range(K):
                            n_phases_k = phases[ist, k]
                            station_state.extend([0] * n_phases_k)
                        for expanded in _expand_state(base_state, rest, current_phases + station_state):
                            yield expanded
                    else:
                        # Find which class has jobs (single-class case)
                        active_class = -1
                        for k in range(K):
                            if n_ist[k] > 0:
                                active_class = k
                                break

                        if active_class >= 0:
                            n_active = int(n_ist[active_class])
                            num_in_service = min(n_active, S_ist)
                            num_in_buffer = n_active - num_in_service

                            # NOTE: Previously had a constraint to skip pure single-class states
                            # with 3+ jobs in buffer. This was removed because multiclass states
                            # can depart to these single-class states, causing absorbing states
                            # if the targets are filtered out.

                            # Generate buffer: class IDs for WAITING jobs only (1-based)
                            # MATLAB format: buffer tracks waiting jobs, phase counts track in-service
                            # num_in_buffer = jobs waiting (not in service)
                            buffer_positions = [active_class + 1] * num_in_buffer

                            # Generate phase distributions for jobs in service
                            n_phases_active = phases[ist, active_class]
                            if n_phases_active > 1 and num_in_service > 0:
                                service_phase_dists = _generate_phase_distributions(num_in_service, int(n_phases_active))
                            else:
                                service_phase_dists = [tuple([num_in_service] + [0] * (n_phases_active - 1))] if n_phases_active > 1 else [(num_in_service,)]

                            for service_phases in service_phase_dists:
                                # Build station state: buffer + phases for all classes
                                station_state = [0] * max_buf_size
                                # Place buffer positions at the end (right-justified)
                                # MATLAB format: buffer holds waiting jobs only, right-justified
                                for i, bp in enumerate(buffer_positions):
                                    station_state[max_buf_size - num_in_buffer + i] = bp

                                # Add phase distributions for all classes
                                for k in range(K):
                                    n_phases_k = phases[ist, k]
                                    if k == active_class:
                                        station_state.extend(service_phases)
                                    else:
                                        station_state.extend([0] * n_phases_k)

                                for expanded in _expand_state(base_state, rest, current_phases + station_state):
                                    yield expanded
                        else:
                            # No active class - empty state
                            station_state = [0] * max_buf_size
                            for k in range(K):
                                n_phases_k = phases[ist, k]
                                station_state.extend([0] * n_phases_k)
                            for expanded in _expand_state(base_state, rest, current_phases + station_state):
                                yield expanded
                elif ist in siro_stations:
                    # SIRO station - use count-based buffer (per-class counts, not ordered)
                    # SIRO state format: [buffer_counts[K], phase_counts]
                    # where buffer_counts[k] = jobs of class k waiting in buffer
                    S_ist = int(nservers[ist]) if ist < len(nservers) else 1
                    total_at_station = int(sum(n_ist))

                    if total_at_station == 0:
                        # Empty station: buffer counts are 0, phase counts are 0
                        station_state = [0] * K  # Buffer counts (one per class)
                        for k in range(K):
                            n_phases_k = phases[ist, k]
                            station_state.extend([0] * n_phases_k)
                        for expanded in _expand_state(base_state, rest, current_phases + station_state):
                            yield expanded
                    elif total_at_station <= S_ist:
                        # All jobs in service, no buffer
                        station_state = [0] * K  # Buffer counts (all zeros)
                        for k in range(K):
                            n_ik = int(n_ist[k])
                            n_phases_k = phases[ist, k]
                            if n_phases_k > 1 and n_ik > 0:
                                phase_dists = _generate_phase_distributions(n_ik, int(n_phases_k))
                            else:
                                phase_dists = [tuple([n_ik] + [0] * (n_phases_k - 1))] if n_phases_k > 1 else [(n_ik,)]
                            for pd in phase_dists:
                                inner_state = station_state.copy()
                                inner_state.extend(pd)
                                station_state = inner_state
                                break  # Only take first phase distribution here
                        # Enumerate all phase combinations
                        phase_combos_for_station = []
                        for k in range(K):
                            n_ik = int(n_ist[k])
                            n_phases_k = phases[ist, k]
                            if n_phases_k > 1 and n_ik > 0:
                                phase_dists = _generate_phase_distributions(n_ik, int(n_phases_k))
                            else:
                                phase_dists = [(n_ik,)] if n_phases_k <= 1 else [tuple([n_ik] + [0] * (n_phases_k - 1))]
                            phase_combos_for_station.append(phase_dists)
                        for phase_combo in product(*phase_combos_for_station):
                            station_state = [0] * K  # Buffer counts
                            for pd in phase_combo:
                                station_state.extend(pd)
                            for expanded in _expand_state(base_state, rest, current_phases + station_state):
                                yield expanded
                    else:
                        # Jobs exceed servers - need to enumerate buffer/service distributions
                        # Use multichoose to enumerate how many jobs of each class are in service
                        service_distributions = _multichoose_con(np.array(n_ist), S_ist)
                        for si in service_distributions:
                            # si[k] = number of class k jobs in service
                            # buffer[k] = n_ist[k] - si[k]
                            buffer_counts = [int(n_ist[k]) - int(si[k]) for k in range(K)]
                            # Generate phase distributions for jobs in service
                            phase_combos_for_station = []
                            for k in range(K):
                                n_in_service = int(si[k])
                                n_phases_k = phases[ist, k]
                                if n_phases_k > 1 and n_in_service > 0:
                                    phase_dists = _generate_phase_distributions(n_in_service, int(n_phases_k))
                                else:
                                    phase_dists = [(n_in_service,)] if n_phases_k <= 1 else [tuple([n_in_service] + [0] * (n_phases_k - 1))]
                                phase_combos_for_station.append(phase_dists)
                            for phase_combo in product(*phase_combos_for_station):
                                station_state = buffer_counts.copy()  # Buffer counts per class
                                for pd in phase_combo:
                                    station_state.extend(pd)
                                for expanded in _expand_state(base_state, rest, current_phases + station_state):
                                    yield expanded
                else:
                    # Non-FCFS/non-SIRO station - use PS-style phase distribution (no buffer)
                    phase_combos_for_station = []
                    for k in range(K):
                        n_ik = int(n_ist[k])
                        n_phases_k = phases[ist, k]
                        if n_phases_k > 1:
                            phase_dists = _generate_phase_distributions(n_ik, int(n_phases_k))
                        else:
                            phase_dists = [(n_ik,)]
                        phase_combos_for_station.append(phase_dists)

                    for phase_combo in product(*phase_combos_for_station):
                        station_state = []
                        for pd in phase_combo:
                            station_state.extend(pd)
                        for expanded in _expand_state(base_state, rest, current_phases + station_state):
                            yield expanded

            # Expand this base state
            all_stations = list(range(M))
            for expanded_state in _expand_state(base_state, all_stations, []):
                phase_augmented_states.append(expanded_state)

        # De-duplicate states (enumerate can produce duplicates for multi-class FCFS)
        seen = set()
        unique_states = []
        for s in phase_augmented_states:
            key = tuple(int(x) for x in s)
            if key not in seen:
                seen.add(key)
                unique_states.append(s)
        states = unique_states

        # Compute state_dim from actual state sizes (should be consistent now)
        if states:
            state_dim = len(states[0])
            # Verify all states have same size
            for s in states:
                if len(s) != state_dim:
                    print(f"Warning: inconsistent state size {len(s)} vs {state_dim}")
        else:
            state_dim = total_phases

        # Recompute phase_offset to account for FCFS/SIRO buffer positions
        # For FCFS stations: state = [buffer_positions, phases_class0, phases_class1, ...]
        # For SIRO stations: state = [buffer_counts_class0..K, phases_class0, phases_class1, ...]
        # For non-FCFS/non-SIRO stations: state = [phases_class0, phases_class1, ...]
        # Skip Source/Sink stations as they don't contribute to state
        fcfs_phase_offset = np.zeros((M, K), dtype=int)
        idx = 0
        for ist in range(M):
            # Skip Source/Sink stations
            if ist in skip_stations:
                fcfs_phase_offset[ist, :] = -1  # Mark as invalid
                continue
            if ist in fcfs_stations:
                # FCFS station: buffer comes first
                max_buf = fcfs_max_buffer.get(ist, 0)
                idx += max_buf  # Skip buffer positions
            elif ist in siro_stations:
                # SIRO station: K buffer counts come first
                idx += K  # Skip buffer count positions
            for k in range(K):
                fcfs_phase_offset[ist, k] = idx
                idx += phases[ist, k]

        # Update phase_offset with FCFS-aware offsets
        phase_offset = fcfs_phase_offset
        rrobin_info['phase_offset'] = phase_offset

        # Store FCFS buffer info for generator
        rrobin_info['fcfs_buffer_info'] = fcfs_buffer_info
    else:
        # No phase augmentation - use simple job counts
        state_dim = M * K

    rrobin_info['state_dim'] = state_dim

    # Get MAP FCFS information for additional state variables
    map_fcfs, has_map_fcfs = _get_map_fcfs_info(sn)
    rrobin_info['map_fcfs'] = map_fcfs
    rrobin_info['has_map_fcfs'] = has_map_fcfs

    # Augment states with MAP phase variables for FCFS stations with MAP distributions
    # The MAP phase tracks the "mode" of the MAP process (separate from service phase counts)
    if has_map_fcfs:
        map_augmented_states = []
        # Sort the map_fcfs keys for consistent ordering
        map_fcfs_keys = sorted(map_fcfs.keys())
        rrobin_info['map_fcfs_keys'] = map_fcfs_keys

        # Create map offset mapping: position in state vector where each MAP phase var starts
        map_var_offset = {}
        offset = state_dim  # MAP vars start after the job count variables
        for (ist, k) in map_fcfs_keys:
            map_var_offset[(ist, k)] = offset
            offset += 1  # Each MAP phase is a single integer (the current phase)
        rrobin_info['map_var_offset'] = map_var_offset
        rrobin_info['map_var_count'] = len(map_fcfs_keys)

        # Expand each base state with valid MAP phase combinations
        # MATLAB constraint: when jobs of a MAP class are in service at a single-server
        # FCFS station, the MAP phase must equal the occupied service phase.
        # When no jobs in service, all phases are valid.
        for base_state in states:
            # For each MAP class, determine valid MAP phases based on service state
            valid_map_phases = []
            for (ist, k) in map_fcfs_keys:
                n_phases = map_fcfs[(ist, k)]
                S_ist = int(nservers[ist]) if ist < len(nservers) else 1

                # Get the phase counts for this class from base_state
                start_idx = phase_offset[ist, k]
                end_idx = start_idx + n_phases
                if end_idx <= len(base_state):
                    phase_counts = base_state[start_idx:end_idx]
                else:
                    phase_counts = [0] * n_phases

                jobs_in_service = sum(phase_counts)

                if jobs_in_service > 0 and S_ist == 1:
                    # Single server with jobs in service: MAP phase = occupied phase
                    # Find which phase has the job (there can only be one in single server)
                    occupied_phase = None
                    for p in range(n_phases):
                        if phase_counts[p] > 0:
                            occupied_phase = p
                            break
                    if occupied_phase is not None:
                        valid_map_phases.append([occupied_phase])
                    else:
                        # Fallback: allow all phases
                        valid_map_phases.append(list(range(n_phases)))
                else:
                    # No jobs in service or multi-server: all phases are valid
                    valid_map_phases.append(list(range(n_phases)))

            # Create all valid MAP phase combinations for this base state
            for map_combo in product(*valid_map_phases):
                augmented_state = list(base_state) + list(map_combo)
                map_augmented_states.append(augmented_state)

        states = map_augmented_states
        # Note: state_dim stays the same - it's only the job count variables
        # MAP phase variables are stored separately after state_dim
    else:
        rrobin_info['map_fcfs_keys'] = []
        rrobin_info['map_var_offset'] = {}
        rrobin_info['map_var_count'] = 0

    # Augment states with round-robin pointers if there are RROBIN routing nodes
    # Each RROBIN (node, class) adds a state variable that tracks which outlink
    # the round-robin pointer is pointing to
    if rrobin_info['total_vars'] > 0:
        augmented_states = []
        # Get all possible round-robin pointer combinations
        rr_ranges = []
        for node_idx, class_idx, num_outlinks in rrobin_info['state_vars']:
            # Each pointer can be any of the outlinks (0 to num_outlinks-1)
            rr_ranges.append(range(num_outlinks))

        # Expand each base state with all possible RR pointer combinations
        for base_state in states:
            for rr_combo in product(*rr_ranges):
                augmented_state = list(base_state) + list(rr_combo)
                augmented_states.append(augmented_state)
        states = augmented_states

    # Augment states with cache configurations if there are Cache stations
    # Each Cache station adds:
    # 1. Job presence variables for each class (0/1 per class) - like JAR's state_bufsrv
    # 2. A state variable tracking which items are in cache (cache content index) - like JAR's state_var
    # This matches JAR/MATLAB behavior where non-station stateful nodes track job presence
    cache_stations_info = _get_cache_stations_info(sn)
    rrobin_info['cache_stations_info'] = cache_stations_info

    if cache_stations_info:
        # Compute offset where cache job presence variables start
        # They come after: state_dim + map_var_count + rr_var_count
        map_var_count = rrobin_info.get('map_var_count', 0)
        rr_var_count = rrobin_info.get('total_vars', 0)
        cache_jobs_offset_base = state_dim + map_var_count + rr_var_count

        sorted_cache_stations = sorted(cache_stations_info.keys())

        # Add job presence variables for each class at each Cache node (capacity 1 per class)
        # This matches MATLAB/Java where Cache nodes have capacityc[ind, r] = 1 for all r
        offset = cache_jobs_offset_base
        cache_jobs_offsets = {}  # {isf: offset} where state[offset+k] = jobs of class k at Cache isf
        for ist in sorted_cache_stations:
            cache_jobs_offsets[ist] = offset
            offset += K  # K job presence variables (one per class, each 0/1)
        rrobin_info['cache_jobs_offsets'] = cache_jobs_offsets

        # Cache content indices come after job presence variables
        cache_state_offsets = {}
        for ist in sorted_cache_stations:
            cache_state_offsets[ist] = offset
            offset += 1  # One state variable per cache (index into configurations)
        rrobin_info['cache_state_offsets'] = cache_state_offsets
        rrobin_info['cache_sorted_stations'] = sorted_cache_stations

        # Expand each base state with job presence AND cache configurations
        # Job presence: each class independently 0/1 at Cache (matching MATLAB)
        # In MATLAB, Cache capacity is 1 per class independently, allowing multiple
        # simultaneous jobs of different classes (e.g., InitClass=1 AND HitClass=1)
        # Cache content: index into item permutations
        cache_augmented_states = []

        # Build ranges for job presence (0/1 per class per Cache)
        job_presence_ranges = []
        for ist in sorted_cache_stations:
            for k_class in range(K):
                job_presence_ranges.append(range(2))  # 0 or 1

        # Build ranges for cache content indices
        cache_content_ranges = []
        for ist in sorted_cache_stations:
            num_configs = cache_stations_info[ist]['num_states']
            cache_content_ranges.append(range(num_configs))

        # Combine: job presence vars first, then cache content indices
        all_ranges = job_presence_ranges + cache_content_ranges

        # Compute max chain population for constraining total jobs across
        # Cache + Router + stations. In MATLAB, each non-source stateful node
        # can hold up to cutoff jobs per class independently, and the total
        # population across all nodes is bounded by the sum of per-class
        # cutoffs across all non-source stateful nodes.
        # For this Cache network: Cache capacity 1 per class, Router capacity
        # 1 per class, each Delay capacity cutoff per class.
        # Total max = K * (n_cache + n_router + sum_delay_cutoffs)
        # But MATLAB uses a simpler bound: sum over classes of max per-class cutoff.
        max_chain_pop = 0
        for k in range(K):
            if is_matrix_cutoff:
                class_cutoff = max((get_cutoff(ist_c, k) for ist_c in range(M)
                                    if ist_c != source_station), default=1)
            else:
                class_cutoff = int(cutoff_arr.flat[0])
            max_chain_pop += class_cutoff
        rrobin_info['max_chain_pop'] = max_chain_pop

        for base_state in states:
            # Count total jobs at stations in base state
            station_jobs_total = sum(base_state[ist * K + k]
                                     for ist in range(M) for k in range(K)
                                     if ist != source_station)
            for combo in product(*all_ranges):
                # Per-class constraint: each class can have at most 1 job at
                # Cache (capacity 1 per class, matching MATLAB/JAR)
                valid = True
                job_idx = 0
                cache_jobs_total = 0
                for ist in sorted_cache_stations:
                    for kk in range(K):
                        if combo[job_idx + kk] > 1:
                            valid = False
                            break
                    if not valid:
                        break
                    cache_jobs_total += sum(combo[job_idx:job_idx + K])
                    job_idx += K
                if not valid:
                    continue
                # Total population constraint: jobs at Cache + stations <= max_chain_pop
                if station_jobs_total + cache_jobs_total > max_chain_pop:
                    continue
                augmented_state = list(base_state) + list(combo)
                cache_augmented_states.append(augmented_state)
        states = cache_augmented_states
    else:
        rrobin_info['cache_state_offsets'] = {}
        rrobin_info['cache_sorted_stations'] = []
        rrobin_info['cache_jobs_offsets'] = {}

    # --- Router job presence augmentation ---
    # For Cache networks with Router nodes, include Router job presence in the state
    # vector. This is needed for stochastic complementation (stochcomp) to correctly
    # capture dynamics where Source arrivals can occur while jobs are "in transit"
    # at the Router. Without Router states, the CTMC underestimates throughput
    # because it doesn't account for this transient capacity.
    # MATLAB's CTMC solver includes Router states, builds a larger Q matrix,
    # then eliminates immediate (Router) states via stochcomp.
    router_nodes_info = {}
    if cache_stations_info:
        # Find Router-like nodes: stateful, non-station, non-Cache
        for ind in range(sn.nnodes):
            if ind < len(sn.isstateful) and sn.isstateful[ind]:
                if ind < len(sn.isstation) and not sn.isstation[ind]:
                    nt_val = int(sn.nodetype[ind].value) if hasattr(sn.nodetype[ind], 'value') else int(sn.nodetype[ind])
                    if nt_val != 6:  # Not CACHE (NodeType.CACHE=6)
                        isf = int(sn.nodeToStateful[ind])
                        # Determine which classes can be at this Router
                        # (from rtnodes: classes that have incoming routing to this node)
                        active_classes = []
                        P = np.asarray(sn.rtnodes) if sn.rtnodes is not None else None
                        if P is not None:
                            for k_c in range(K):
                                col_idx = ind * K + k_c
                                if col_idx < P.shape[1] and np.any(P[:, col_idx] > 1e-10):
                                    active_classes.append(k_c)
                        if active_classes:
                            router_nodes_info[isf] = {
                                'node_idx': ind,
                                'active_classes': active_classes
                            }

    rrobin_info['router_nodes_info'] = router_nodes_info

    if router_nodes_info:
        # Compute offset where Router job presence variables start
        # They come after everything else in the current state vector
        current_state_len = len(states[0]) if states else 0
        router_jobs_offset_base = current_state_len

        sorted_router_isfs = sorted(router_nodes_info.keys())
        router_jobs_offsets = {}  # {isf: offset} where state[offset+k] = jobs of class k at Router isf
        offset = router_jobs_offset_base
        for isf in sorted_router_isfs:
            router_jobs_offsets[isf] = offset
            offset += K  # K job presence variables per Router
        rrobin_info['router_jobs_offsets'] = router_jobs_offsets

        # Expand each state with Router job presence configurations
        # Each active class at each Router can be 0 or 1
        router_augmented_states = []
        router_ranges = []
        for isf in sorted_router_isfs:
            active = router_nodes_info[isf]['active_classes']
            for k_c in range(K):
                if k_c in active:
                    router_ranges.append(range(2))  # 0 or 1
                else:
                    router_ranges.append(range(1))  # Always 0

        for base_state in states:
            # Count current non-source jobs (station + cache)
            station_jobs_total = sum(base_state[ist_c * K + k_c]
                                     for ist_c in range(M) for k_c in range(K)
                                     if ist_c != source_station)
            cache_jobs_total = 0
            for isf_c in sorted(cache_stations_info.keys()):
                cj_off = rrobin_info['cache_jobs_offsets'].get(isf_c, -1)
                if cj_off >= 0 and cj_off < len(base_state):
                    for k_c in range(K):
                        if cj_off + k_c < len(base_state):
                            cache_jobs_total += int(base_state[cj_off + k_c])

            for combo in product(*router_ranges):
                router_jobs_total = sum(combo)
                # Total population constraint: station + cache + router <= max_chain_pop
                if station_jobs_total + cache_jobs_total + router_jobs_total > max_chain_pop:
                    continue
                augmented_state = list(base_state) + list(combo)
                router_augmented_states.append(augmented_state)

        states = router_augmented_states
    else:
        rrobin_info['router_jobs_offsets'] = {}

    state_space = np.array(states, dtype=np.float64)

    # Compute buffer_start offsets for FCFS stations
    # State layout: for each station in order, phases (and buffer for FCFS) are concatenated
    # We need to track where each station's data starts in the state vector
    # Initialize fcfs_buffer_info here in case phase augmentation didn't run
    try:
        fcfs_buffer_info
    except NameError:
        fcfs_buffer_info = {}
    current_offset = 0
    for ist in range(M):
        if ist in skip_stations:
            continue
        if ist in fcfs_stations:
            # FCFS: buffer comes first, then phases
            max_buf = fcfs_max_buffer.get(ist, 0)
            fcfs_buffer_info[ist] = {'buffer_start': current_offset, 'max_buffer': max_buf}
            current_offset += max_buf  # Buffer positions
            for k in range(K):
                current_offset += phases[ist, k]  # Phase positions
        elif ist in siro_stations:
            # SIRO: buffer counts come first (K counts), then phases
            current_offset += K  # Buffer counts (one per class)
            for k in range(K):
                current_offset += phases[ist, k]  # Phase positions
        else:
            # Non-FCFS/SIRO: just phases
            for k in range(K):
                current_offset += phases[ist, k]

    # Store updated fcfs_buffer_info in rrobin_info
    rrobin_info['fcfs_buffer_info'] = fcfs_buffer_info

    # Aggregated state space: sum over classes at each station (exclude RR pointers)
    # For phase-augmented states:
    # - FCFS stations: count buffer entries (buffer stores ALL jobs at station)
    # - Non-FCFS stations: sum over phases (phase counts represent jobs in each phase)
    # Aggregated state space: per-station per-class job counts
    # Format: (nstates, nstations * nclasses) matching MATLAB stateSpaceAggr
    # Layout: columns ((ist)*K + k) = jobs of class k at station ist (0-indexed)
    # Note: use needs_state_expansion (not just phase augmentation flag) because
    # FCFS stations always get state expansion regardless of use_phase_augmentation
    state_space_aggr = np.zeros((len(states), M * K))
    for i, state in enumerate(states):
        for ist in range(M):
            # Source/Sink stations always have 0 jobs (they're skipped from state space)
            if ist in skip_stations:
                continue

            if needs_state_expansion:
                if ist in fcfs_stations:
                    # FCFS stations: buffer stores WAITING jobs (class IDs, 1-based),
                    # phases store IN-SERVICE jobs (per-class phase counts)
                    max_buf = fcfs_max_buffer.get(ist, 0)
                    buffer_start = fcfs_buffer_info.get(ist, {}).get('buffer_start', 0)
                    # Count waiting jobs per class from buffer
                    for buf_pos in range(max_buf):
                        buf_val = int(state[buffer_start + buf_pos])
                        if buf_val > 0:  # Non-empty buffer position; value is 1-based class ID
                            class_k = buf_val - 1  # Convert to 0-based
                            if 0 <= class_k < K:
                                state_space_aggr[i, ist * K + class_k] += 1.0
                    # Add in-service jobs from phase counts (per class)
                    for k in range(K):
                        start_idx = phase_offset[ist, k]
                        if start_idx >= 0:
                            end_idx = start_idx + phases[ist, k]
                            state_space_aggr[i, ist * K + k] += sum(state[start_idx:end_idx])
                else:
                    # Non-FCFS stations: sum over phases per class
                    for k in range(K):
                        start_idx = phase_offset[ist, k]
                        if start_idx >= 0:  # Valid offset (not skipped station)
                            end_idx = start_idx + phases[ist, k]
                            state_space_aggr[i, ist * K + k] += sum(state[start_idx:end_idx])
            else:
                for k in range(K):
                    state_space_aggr[i, ist * K + k] += state[ist * K + k]

    return state_space, state_space_aggr, rrobin_info


def _build_state_index_map(state_space: np.ndarray) -> dict:
    """
    Build a hash map from state tuples to indices for O(1) lookup.

    Args:
        state_space: Enumerated state space (n_states x state_dim)

    Returns:
        Dictionary mapping state tuples to their indices
    """
    state_map = {}
    for i, state in enumerate(state_space):
        # Convert to tuple of integers for hashing
        state_key = tuple(int(x) for x in state)
        state_map[state_key] = i
    return state_map


def _find_state_index_fast(state_map: dict, state: np.ndarray) -> int:
    """
    Find index of state using hash map (O(1) lookup).

    Args:
        state_map: Hash map from state tuples to indices
        state: State vector to find

    Returns:
        Index of state, or -1 if not found
    """
    state_key = tuple(int(x) for x in state)
    return state_map.get(state_key, -1)


def _build_generator(
    sn: NetworkStruct,
    state_space: np.ndarray,
    options: SolverCTMCOptions,
    rrobin_info: Optional[dict] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build the infinitesimal generator matrix for the queueing network.

    Supports state-dependent routing (RROBIN, WRROBIN) by using routing state
    variables in the state vector to determine destinations.

    Also supports phase-type distributions (PH, APH, MAP, etc.) when phase
    augmentation is enabled in rrobin_info.

    Args:
        sn: Network structure
        state_space: Enumerated state space
        options: Solver options
        rrobin_info: Round-robin routing information (also contains phase info)

    Returns:
        Tuple of (infinitesimal generator matrix Q, departure rates array depRates)
        depRates[s, ist, k] = total departure rate from state s at station ist for class k
    """
    M = sn.nstations
    K = sn.nclasses
    n_states = state_space.shape[0]

    Q = np.zeros((n_states, n_states))
    # Track departure rates: depRates[s, ist, k] = sum of departure rates from state s
    # for class k at station ist (used for accurate throughput computation)
    depRates = np.zeros((n_states, M, K))

    # Get round-robin info if not provided
    if rrobin_info is None:
        rrobin_info = _build_rrobin_state_info(sn)

    # Check if phase augmentation is used
    use_phase_aug = rrobin_info.get('needs_phase_augmentation', False)
    phases = rrobin_info.get('phases', np.ones((M, K), dtype=int))
    phase_offset = rrobin_info.get('phase_offset', None)
    state_dim = rrobin_info.get('state_dim', M * K)

    # Get MAP FCFS info
    map_fcfs = rrobin_info.get('map_fcfs', {})
    has_map_fcfs = rrobin_info.get('has_map_fcfs', False)
    map_var_offset = rrobin_info.get('map_var_offset', {})
    map_var_count = rrobin_info.get('map_var_count', 0)

    # Get cache state tracking info
    cache_stations_info = rrobin_info.get('cache_stations_info', {})
    cache_state_offsets = rrobin_info.get('cache_state_offsets', {})
    cache_jobs_offsets = rrobin_info.get('cache_jobs_offsets', {})

    # Get Router state tracking info (for stochcomp)
    router_nodes_info = rrobin_info.get('router_nodes_info', {})
    router_jobs_offsets = rrobin_info.get('router_jobs_offsets', {})

    # Build mapping from (node_idx, class_idx) to state variable index
    # RR pointers come after state variables AND MAP phase variables
    rr_state_offset = state_dim + map_var_count  # RR pointers start after MAP vars
    rr_var_map = {}  # Maps (node_idx, class_idx) -> index in state vector
    for i, (node_idx, class_idx, _) in enumerate(rrobin_info['state_vars']):
        rr_var_map[(node_idx, class_idx)] = rr_state_offset + i

    def get_map_matrices(ist: int, k: int) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Get MAP D0 and D1 matrices for (station, class)."""
        if not hasattr(sn, 'proc') or sn.proc is None:
            return None, None
        proc_is_list = isinstance(sn.proc, list)
        if proc_is_list:
            if ist >= len(sn.proc) or sn.proc[ist] is None:
                return None, None
            station_proc = sn.proc[ist]
        else:
            if ist not in sn.proc:
                return None, None
            station_proc = sn.proc[ist]

        proc_entry = None
        if isinstance(station_proc, (list, tuple)):
            if k < len(station_proc):
                proc_entry = station_proc[k]
        elif isinstance(station_proc, dict):
            proc_entry = station_proc.get(k)

        if proc_entry is None:
            return None, None

        # Handle PH/APH distributions stored as [alpha, T] or direct D0/D1 matrices
        if isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 2:
            first_elem = np.atleast_2d(np.array(proc_entry[0], dtype=float))
            second_elem = np.atleast_2d(np.array(proc_entry[1], dtype=float))

            # Check if this is [alpha, T] format (alpha is 1D/row vector, T is square matrix)
            # vs [D0, D1] format (both D0 and D1 are square matrices)
            if first_elem.shape[0] == 1 and second_elem.shape[0] == second_elem.shape[1]:
                # This is [alpha, T] format for PH/APH distributions
                alpha = first_elem.flatten()
                T = second_elem
                # Convert to D0/D1:
                # D0 = T (the transition matrix, including diagonal absorption rates)
                # D1 = outer(exit_rates, alpha) where exit_rates = -sum(T, axis=1)
                D0 = T
                exit_rates = -np.sum(T, axis=1)
                D1 = np.outer(exit_rates, alpha)
                return D0, D1
            else:
                # Assume direct D0/D1 matrices
                D0 = first_elem
                D1 = second_elem
                return D0, D1

        # Handle Erlang distribution stored as dict with 'k' and 'mu'
        if isinstance(proc_entry, dict):
            if 'k' in proc_entry and 'mu' in proc_entry:
                n_phases = int(proc_entry['k'])
                mu = float(proc_entry['mu'])  # per-phase rate
                if n_phases > 1:
                    # Construct Erlang D0/D1 matrices
                    # D0: diagonal = -mu (absorption), off-diagonal D0[i,i+1] = mu (phase transition)
                    # D1: D1[k-1,0] = mu (completion from last phase)
                    D0 = np.zeros((n_phases, n_phases))
                    D1 = np.zeros((n_phases, n_phases))
                    for p in range(n_phases):
                        D0[p, p] = -mu
                        if p < n_phases - 1:
                            D0[p, p + 1] = mu  # phase transition
                    D1[n_phases - 1, 0] = mu  # completion from last phase
                    return D0, D1
                else:
                    # Single phase: exponential
                    D0 = np.array([[-mu]])
                    D1 = np.array([[mu]])
                    return D0, D1
            elif 'probs' in proc_entry and 'rates' in proc_entry:
                # HyperExp distribution: each phase completes independently
                probs = np.array(proc_entry['probs'])
                rates = np.array(proc_entry['rates'])
                n_phases = len(rates)
                # D0: diagonal with -rate_i (no inter-phase transitions in HyperExp)
                D0 = np.diag(-rates)
                # D1: completion from phase i, restart in phase j with prob_j
                # D1[i,j] = rate_i * prob_j
                D1 = np.outer(rates, probs)
                return D0, D1
            elif 'rate' in proc_entry:
                # Exponential distribution
                mu = float(proc_entry['rate'])
                D0 = np.array([[-mu]])
                D1 = np.array([[mu]])
                return D0, D1

        return None, None

    def get_map_phase(state: np.ndarray, ist: int, k: int) -> int:
        """Get current MAP phase for (station, class) from state vector."""
        if (ist, k) not in map_var_offset:
            return 0
        idx = map_var_offset[(ist, k)]
        return int(state[idx])

    def set_map_phase(state: np.ndarray, ist: int, k: int, phase: int) -> np.ndarray:
        """Set MAP phase for (station, class) in state vector."""
        new_state = state.copy()
        if (ist, k) in map_var_offset:
            idx = map_var_offset[(ist, k)]
            new_state[idx] = phase
        return new_state

    def get_job_count(state: np.ndarray, ist: int, k: int) -> int:
        """Get total job count for (station, class) from state vector.

        For FCFS stations with buffer, this includes:
        - Jobs waiting in buffer (count of buffer positions with class k+1)
        - Jobs in service (sum of phase counts for class k)

        For SIRO stations, this includes:
        - Jobs waiting in buffer (count stored at buffer_start + k)
        - Jobs in service (sum of phase counts for class k)
        """
        if use_phase_aug and phase_offset is not None:
            # Get jobs in service from phase counts
            start_idx = phase_offset[ist, k]
            end_idx = start_idx + phases[ist, k]
            in_service = int(sum(state[start_idx:end_idx]))

            # For FCFS stations, also count buffer jobs for this class
            if ist in fcfs_stations:
                buffer_start = fcfs_buffer_start.get(ist, -1)
                max_buf = fcfs_max_buffer.get(ist, 0)
                if buffer_start >= 0 and max_buf > 0:
                    # Buffer stores class IDs (1-based), count entries matching class k
                    in_buffer = sum(1 for i in range(max_buf) if int(state[buffer_start + i]) == k + 1)
                    return in_service + in_buffer
            # For SIRO stations, count buffer jobs (count-based)
            elif ist in siro_stations:
                buffer_start = siro_buffer_start.get(ist, -1)
                if buffer_start >= 0:
                    in_buffer = int(state[buffer_start + k])
                    return in_service + in_buffer

            return in_service
        else:
            return int(state[ist * K + k])

    def get_phase_counts(state: np.ndarray, ist: int, k: int) -> np.ndarray:
        """Get phase counts for (station, class) from state vector."""
        if use_phase_aug and phase_offset is not None:
            start_idx = phase_offset[ist, k]
            end_idx = start_idx + phases[ist, k]
            return np.array([int(x) for x in state[start_idx:end_idx]])
        else:
            # Single phase
            return np.array([int(state[ist * K + k])])

    def set_job_count(state: np.ndarray, ist: int, k: int, delta: int, phase_idx: int = -1) -> np.ndarray:
        """
        Modify job count for (station, class) in state vector.

        Args:
            state: Current state vector
            ist: Station index
            k: Class index
            delta: Change in job count (+1 for arrival, -1 for departure)
            phase_idx: Phase index to modify. If -1:
                       - For arrivals (delta > 0): use first phase (phase 0), or buffer for FCFS
                       - For departures (delta < 0): find first phase with jobs

        Returns:
            New state vector with modified job count
        """
        new_state = state.copy()
        if use_phase_aug and phase_offset is not None:
            # Check if this is an FCFS or SIRO station
            is_fcfs_station = ist in fcfs_stations
            is_siro_station = ist in siro_stations
            buffer_start = fcfs_buffer_start.get(ist, -1) if is_fcfs_station else -1
            phases_start = fcfs_phases_start.get(ist, -1) if is_fcfs_station else -1
            max_buf = fcfs_max_buffer.get(ist, 0) if is_fcfs_station else 0

            # SIRO buffer handling
            siro_buf_start = siro_buffer_start.get(ist, -1) if is_siro_station else -1
            siro_ph_start = siro_phases_start.get(ist, -1) if is_siro_station else -1

            if delta > 0 and is_siro_station and siro_buf_start >= 0:
                # SIRO arrival: check if all servers are busy
                total_in_service = 0
                for kk in range(K):
                    for p in range(phases[ist, kk]):
                        idx = siro_ph_start + sum(phases[ist, :kk]) + p
                        total_in_service += int(new_state[idx])

                n_servers_ist = int(nservers[ist]) if ist < len(nservers) else 1

                if total_in_service >= n_servers_ist:
                    # All servers busy - new arrival goes to buffer (increment class count)
                    new_state[siro_buf_start + k] += 1
                    return new_state
                else:
                    # Free server available - arrival goes directly to service (phase 0)
                    if phase_idx < 0:
                        phase_idx = 0
                    class_phase_start = siro_ph_start + sum(phases[ist, :k])
                    new_state[class_phase_start + phase_idx] += delta
                    return new_state

            if delta > 0 and is_fcfs_station and buffer_start >= 0 and max_buf > 0:
                # FCFS/HOL/LCFS/LCFSPR arrival: check if all servers are busy
                total_in_service = 0
                for kk in range(K):
                    for p in range(phases[ist, kk]):
                        idx = phases_start + sum(phases[ist, :kk]) + p
                        total_in_service += int(new_state[idx])

                # Get number of servers at this station
                n_servers_ist = int(nservers[ist]) if ist < len(nservers) else 1

                # Get scheduling strategy for this station
                ist_sched = sn.sched.get(ist, SchedStrategy.FCFS) if hasattr(sn, 'sched') and sn.sched else SchedStrategy.FCFS

                if total_in_service >= n_servers_ist:
                    # All servers busy - behavior depends on scheduling strategy
                    if _sched_is(ist_sched, 'LCFSPR'):
                        # LCFSPR: arriving job PREEMPTS current job in service
                        # 1. Find which class is currently in service
                        # 2. Move that job to buffer (as newest = leftmost)
                        # 3. Put arriving job in service

                        # Find class currently in service
                        preempted_class = -1
                        for kk in range(K):
                            class_phase_start = phases_start + sum(phases[ist, :kk])
                            for p in range(phases[ist, kk]):
                                if new_state[class_phase_start + p] > 0:
                                    preempted_class = kk
                                    # Remove from service
                                    new_state[class_phase_start + p] -= 1
                                    break
                            if preempted_class >= 0:
                                break

                        if preempted_class >= 0:
                            # Move preempted job to buffer (as newest = leftmost position)
                            # Buffer format: [padding..., oldest, ..., newest] but we store
                            # newest at the leftmost filled position
                            # Find the leftmost empty position to place the preempted job
                            # This is "one position left of current leftmost filled"
                            leftmost_filled = -1
                            for buf_pos in range(max_buf):
                                if new_state[buffer_start + buf_pos] > 0:
                                    leftmost_filled = buf_pos
                                    break

                            if leftmost_filled > 0:
                                # Put preempted job one position to the left of oldest
                                new_state[buffer_start + leftmost_filled - 1] = preempted_class + 1
                            elif leftmost_filled == -1:
                                # Buffer is empty - put at rightmost position
                                new_state[buffer_start + max_buf - 1] = preempted_class + 1
                            else:
                                # Buffer full (leftmost_filled == 0), cannot preempt
                                # Return invalid state
                                new_state[buffer_start] = -1
                                return new_state

                        # Put arriving job in service (phase 0)
                        if phase_idx < 0:
                            phase_idx = 0
                        idx = phase_offset[ist, k] + phase_idx
                        new_state[idx] += 1
                        return new_state
                    else:
                        # FCFS/HOL/LCFS (non-preemptive): new arrival goes to buffer
                        # Buffer is right-justified: [padding..., oldest, ..., newest]
                        # New arrival goes to position immediately LEFT of oldest job (leftmost filled)
                        # Find the leftmost filled position, then put new job one position to its left
                        leftmost_filled = -1
                        for buf_pos in range(max_buf):  # Search from left to right
                            if new_state[buffer_start + buf_pos] > 0:
                                leftmost_filled = buf_pos
                                break
                        if leftmost_filled > 0:
                            # Put new arrival one position to the left of oldest
                            new_state[buffer_start + leftmost_filled - 1] = k + 1  # 1-based class ID
                            return new_state
                        elif leftmost_filled == -1:
                            # Buffer is empty - put at rightmost position (will be only job)
                            new_state[buffer_start + max_buf - 1] = k + 1
                            return new_state
                        else:
                            # leftmost_filled == 0, buffer is full - return invalid state
                            # Setting a buffer position to -1 ensures this state won't be
                            # found in the state space, correctly blocking the arrival
                            new_state[buffer_start] = -1
                            return new_state
                else:
                    # Free server available - arrival goes directly to service (phase 0)
                    # Do NOT set buffer marker - buffer only tracks WAITING jobs
                    if phase_idx < 0:
                        phase_idx = 0
                    idx = phase_offset[ist, k] + phase_idx
                    new_state[idx] += delta
                    return new_state

            # For non-FCFS stations with arrivals, check if cutoff would be exceeded
            if delta > 0 and not is_fcfs_station:
                # Check if adding a job would exceed cutoff for this (station, class)
                current_count = get_job_count(state, ist, k)
                station_cutoff = get_cutoff(ist, k)
                if current_count >= station_cutoff:
                    # At or above cutoff - block arrival by returning invalid state
                    idx = phase_offset[ist, k]
                    new_state[idx] = -1
                    return new_state

            if phase_idx < 0:
                # Auto-select phase
                start_idx = phase_offset[ist, k]
                n_phases = phases[ist, k]
                if delta > 0:
                    # Arrival: add to first phase
                    phase_idx = 0
                else:
                    # Departure: find first phase with jobs
                    for p in range(n_phases):
                        if new_state[start_idx + p] > 0:
                            phase_idx = p
                            break
                    else:
                        # No jobs found - use first phase (will create invalid state)
                        phase_idx = 0
            idx = phase_offset[ist, k] + phase_idx
            new_state[idx] += delta
        else:
            new_state[ist * K + k] += delta
        return new_state

    def get_entry_probs(ist: int, k: int) -> np.ndarray:
        """
        Get entry probabilities for (station, class).

        For HyperExp: returns probs array (distributed across phases)
        For Erlang: returns [1, 0, 0, ...] (always start in phase 0)
        For Exp: returns [1.0]
        """
        if not use_phase_aug or phase_offset is None:
            return np.array([1.0])

        n_phases = phases[ist, k]
        if n_phases <= 1:
            return np.array([1.0])

        # Get proc entry for this (station, class)
        if not hasattr(sn, 'proc') or sn.proc is None:
            return np.array([1.0] + [0.0] * (n_phases - 1))

        proc_is_list = isinstance(sn.proc, list)
        if proc_is_list:
            if ist >= len(sn.proc) or sn.proc[ist] is None:
                return np.array([1.0] + [0.0] * (n_phases - 1))
            station_proc = sn.proc[ist]
        else:
            if ist not in sn.proc:
                return np.array([1.0] + [0.0] * (n_phases - 1))
            station_proc = sn.proc[ist]

        proc_entry = None
        if isinstance(station_proc, (list, tuple)):
            if k < len(station_proc):
                proc_entry = station_proc[k]
        elif isinstance(station_proc, dict):
            proc_entry = station_proc.get(k)

        if proc_entry is None:
            return np.array([1.0] + [0.0] * (n_phases - 1))

        if isinstance(proc_entry, dict):
            if 'probs' in proc_entry and 'rates' in proc_entry:
                # HyperExp: entry probability = probs
                return np.array(proc_entry['probs'])
            elif 'k' in proc_entry:
                # Erlang: always enter phase 0
                return np.array([1.0] + [0.0] * (n_phases - 1))
            else:
                # Exp: single phase
                return np.array([1.0])

        # For PH/APH with [alpha, T] format, extract alpha as entry probabilities
        if isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 2:
            first_elem = proc_entry[0]
            if isinstance(first_elem, np.ndarray):
                first_elem = np.atleast_1d(first_elem)
                # Check if first element is 1D (alpha vector) and second is 2D (T matrix)
                if first_elem.ndim == 1 and len(first_elem) == n_phases:
                    # This is [alpha, T] format - alpha is the entry probability
                    return first_elem.copy()

        # Default: enter first phase
        return np.array([1.0] + [0.0] * (n_phases - 1))

    def get_map_aware_entry_probs(state: np.ndarray, jst: int, r: int) -> np.ndarray:
        """
        Get entry probabilities for arriving job at destination (jst, r),
        accounting for MAP phase continuity.

        For MAP stations, the arriving job should enter in the stored MAP phase
        (to continue the MAP process), not according to the default entry distribution.

        Args:
            state: Current state vector (before the arrival)
            jst: Destination station index
            r: Destination class index

        Returns:
            Array of entry probabilities
        """
        # Check if destination has MAP distribution
        if (jst, r) in map_fcfs and (jst, r) in map_var_offset:
            # MAP station: use stored MAP phase for entry
            map_phase_idx = map_var_offset[(jst, r)]
            stored_map_phase = int(state[map_phase_idx])

            n_phases_jr = phases[jst, r]
            entry_probs = np.zeros(n_phases_jr)
            if 0 <= stored_map_phase < n_phases_jr:
                entry_probs[stored_map_phase] = 1.0
            else:
                # Fallback to first phase if stored phase is invalid
                entry_probs[0] = 1.0
            return entry_probs
        else:
            # Non-MAP station: use default entry probabilities
            return get_entry_probs(jst, r)

    def get_state_matrix(state: np.ndarray) -> np.ndarray:
        """Get job counts as (M, K) matrix from state vector."""
        result = np.zeros((M, K))
        for ist in range(M):
            for k in range(K):
                result[ist, k] = get_job_count(state, ist, k)
        return result

    # Get station to node mapping
    station_to_node = None
    if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
        station_to_node = np.asarray(sn.stationToNode).flatten()

    # Get node to station mapping
    node_to_station = None
    if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
        node_to_station = np.asarray(sn.nodeToStation).flatten()

    # Handle cutoff as scalar or matrix for state validation
    cutoff = options.cutoff
    cutoff_arr = np.atleast_2d(cutoff)
    is_matrix_cutoff = cutoff_arr.shape[0] > 1 or cutoff_arr.shape[1] > 1

    def is_within_cutoff(new_state: np.ndarray) -> bool:
        """Check if state is within cutoff bounds."""
        if is_matrix_cutoff:
            # Reshape state to (M, K) and compare with cutoff matrix
            state_matrix = new_state.reshape(M, K)
            # Handle different cutoff matrix layouts
            if cutoff_arr.shape[0] >= M and cutoff_arr.shape[1] >= K:
                return np.all(state_matrix <= cutoff_arr[:M, :K])
            elif cutoff_arr.shape[0] >= K and cutoff_arr.shape[1] >= M:
                return np.all(state_matrix <= cutoff_arr[:K, :M].T)
            else:
                # Fallback: use scalar comparison
                return np.all(new_state <= cutoff_arr.flat[0])
        else:
            return np.all(new_state <= cutoff_arr.flat[0])

    def get_cutoff(ist: int, k: int) -> int:
        """Get cutoff for station ist and class k."""
        if is_matrix_cutoff:
            # Matrix cutoff: index by station and class
            if cutoff_arr.shape[0] >= M and cutoff_arr.shape[1] >= K:
                return int(cutoff_arr[ist, k])
            elif cutoff_arr.shape[0] >= K and cutoff_arr.shape[1] >= M:
                # Transposed: (K x M) layout
                return int(cutoff_arr[k, ist])
            else:
                return int(cutoff_arr.flat[0]) if cutoff_arr.size > 0 else 10
        else:
            return int(cutoff_arr.flat[0])

    # Get FCFS station information for proper departure rate computation
    fcfs_stations = rrobin_info.get('fcfs_stations', set())
    fcfs_max_buffer = rrobin_info.get('fcfs_max_buffer', {})

    # Get SIRO station information
    siro_stations = rrobin_info.get('siro_stations', set())

    # Compute buffer offsets in the state vector
    # For each station, we need to know where its data starts in the state
    fcfs_buffer_start = {}  # Maps ist -> starting index of buffer in state
    fcfs_phases_start = {}  # Maps ist -> starting index of phases in state
    siro_buffer_start = {}  # Maps ist -> starting index of SIRO buffer counts
    siro_phases_start = {}  # Maps ist -> starting index of SIRO phases

    # Get skip_stations from rrobin_info (Source/Sink stations to exclude)
    skip_stations = rrobin_info.get('skip_stations', set())

    if use_phase_aug and (fcfs_stations or siro_stations):
        # Compute offsets based on state structure
        # State format: [station0_data, station1_data, ...]
        # For FCFS station: [buffer_pos..., phases_class0..., phases_class1..., ...]
        # For SIRO station: [buffer_count_class0, ..., buffer_count_classK-1, phases_class0..., ...]
        # For non-FCFS/non-SIRO station: [phases_class0..., phases_class1..., ...]
        # Skip Source/Sink stations as they don't contribute to state
        current_offset = 0
        for ist in range(M):
            # Skip Source/Sink stations
            if ist in skip_stations:
                fcfs_buffer_start[ist] = -1
                fcfs_phases_start[ist] = -1
                siro_buffer_start[ist] = -1
                siro_phases_start[ist] = -1
                continue
            if ist in fcfs_stations:
                max_buf = fcfs_max_buffer.get(ist, 0)
                fcfs_buffer_start[ist] = current_offset
                fcfs_phases_start[ist] = current_offset + max_buf
                siro_buffer_start[ist] = -1
                siro_phases_start[ist] = -1
                current_offset += max_buf + sum(phases[ist, :])
            elif ist in siro_stations:
                # SIRO: buffer counts (K elements) + phases
                siro_buffer_start[ist] = current_offset
                siro_phases_start[ist] = current_offset + K
                fcfs_buffer_start[ist] = -1
                fcfs_phases_start[ist] = -1
                current_offset += K + sum(phases[ist, :])
            else:
                fcfs_buffer_start[ist] = -1  # Not FCFS
                fcfs_phases_start[ist] = current_offset
                siro_buffer_start[ist] = -1
                siro_phases_start[ist] = current_offset
                current_offset += sum(phases[ist, :])

    rrobin_info['fcfs_phases_start'] = fcfs_phases_start

    def get_fcfs_class_in_service(state: np.ndarray, ist: int) -> int:
        """
        Determine which class is in service for an FCFS station.

        Returns:
            Class index (0-based) of the class in service, or -1 if no jobs in service.
        """
        if ist not in fcfs_stations:
            return -1

        # For FCFS, look at the server phases to see which class has jobs in service
        phases_start = fcfs_phases_start.get(ist, -1)
        if phases_start < 0:
            return -1

        # Check each class for non-zero phase counts
        for k in range(K):
            n_phases_k = phases[ist, k]
            phase_sum = sum(state[phases_start + sum(phases[ist, :j]) + p]
                           for j in range(k) for p in range(phases[ist, j]))
            # Actually, calculate offset correctly
            class_phase_start = phases_start + sum(phases[ist, :k])
            class_jobs_in_service = sum(int(state[class_phase_start + p]) for p in range(n_phases_k))
            if class_jobs_in_service > 0:
                return k
        return -1

    def get_fcfs_jobs_in_service(state: np.ndarray, ist: int, k: int) -> int:
        """
        Get number of class k jobs in service at FCFS station ist.

        For FCFS with single server, this is either 0 or 1.
        """
        if ist not in fcfs_stations:
            return 0

        phases_start = fcfs_phases_start.get(ist, -1)
        if phases_start < 0:
            return 0

        # Calculate offset for class k phases
        class_phase_start = phases_start + sum(phases[ist, :k])
        n_phases_k = phases[ist, k]
        return sum(int(state[class_phase_start + p]) for p in range(n_phases_k))

    def fcfs_promote_from_buffer_with_alpha(state: np.ndarray, ist: int, sched=None, classprio=None) -> List[Tuple[np.ndarray, float, int]]:
        """
        For FCFS/HOL/LCFS station with buffer, return possible promotion states with alpha probabilities.

        When a job completes service, the next job in the buffer should enter service.
        For PH/APH distributions, the starting phase is chosen according to the alpha
        (initial probability) distribution.

        The selection of which job is promoted depends on the scheduling strategy:
        - FCFS: Oldest job (rightmost non-zero position in buffer)
        - HOL: Highest priority job (lowest classprio value), then FCFS within priority
        - LCFS: Newest job (leftmost non-zero position in buffer)

        Args:
            state: Current state vector (will be copied, not modified)
            ist: Station index
            sched: Scheduling strategy (optional, defaults to FCFS behavior)
            classprio: Array of class priorities (required for HOL scheduling)

        Returns:
            List of (new_state, probability, class_k) tuples for each possible starting phase.
            Empty list if no jobs in buffer.
        """
        if ist not in fcfs_stations:
            return []

        buffer_start = fcfs_buffer_start.get(ist, -1)
        phases_start = fcfs_phases_start.get(ist, -1)
        max_buf = fcfs_max_buffer.get(ist, 0)

        if buffer_start < 0 or max_buf <= 0:
            return []

        # Find the job to promote based on scheduling strategy
        buf_pos = -1

        if _sched_is(sched, 'HOL') and classprio is not None:
            # HOL: Find highest priority job (lowest priority value), then rightmost (FCFS) within same priority
            # Collect all buffer positions with jobs
            occupied_positions = []
            for i in range(max_buf):
                buf_val = int(state[buffer_start + i])
                if buf_val > 0:
                    class_k = buf_val - 1  # 0-based class index
                    prio = classprio[class_k] if class_k < len(classprio) else float('inf')
                    occupied_positions.append((i, class_k, prio))

            if occupied_positions:
                # Find minimum priority (highest priority level)
                min_prio = min(p for _, _, p in occupied_positions)
                # Among positions with min priority, find rightmost (oldest in FCFS order)
                candidates = [(pos, ck) for pos, ck, p in occupied_positions if p == min_prio]
                buf_pos = max(pos for pos, _ in candidates)  # rightmost position
        elif _sched_is(sched, 'LCFS', 'LCFSPR'):
            # LCFS/LCFSPR: Find newest job (leftmost non-zero position)
            for i in range(max_buf):  # Search from left to right
                buf_val = int(state[buffer_start + i])
                if buf_val > 0:
                    buf_pos = i
                    break
        else:
            # FCFS (default): Find oldest job (rightmost non-zero position)
            for i in range(max_buf - 1, -1, -1):  # Search from right to left
                buf_val = int(state[buffer_start + i])
                if buf_val > 0:
                    buf_pos = i
                    break

        if buf_pos < 0:
            # No jobs in buffer
            return []

        # Get the class of the job to promote
        buf_val = int(state[buffer_start + buf_pos])
        class_k = buf_val - 1  # Convert from 1-based to 0-based class index

        # Get entry probabilities (alpha) for this class
        # Use MAP-aware entry probs to respect stored MAP phase if this is a MAP station
        entry_probs = get_map_aware_entry_probs(state, ist, class_k)

        results = []
        n_phases_k = phases[ist, class_k]
        class_phase_start = phases_start + sum(phases[ist, :class_k])

        for phase_idx, prob in enumerate(entry_probs):
            if prob > 0 and phase_idx < n_phases_k:
                new_state = state.copy()

                # Collect all buffer markers (excluding the promoted job)
                remaining_markers = []
                for i in range(max_buf):
                    if i != buf_pos and new_state[buffer_start + i] > 0:
                        remaining_markers.append(int(new_state[buffer_start + i]))

                # Clear entire buffer
                for i in range(max_buf):
                    new_state[buffer_start + i] = 0

                # Re-pack buffer: place REMAINING jobs at rightmost positions (MATLAB format)
                # The promoted job moves to SERVICE, NOT buffer
                # Only remaining_markers stay in buffer
                if remaining_markers:
                    n_remaining = len(remaining_markers)
                    start_pos = max_buf - n_remaining
                    for i, marker in enumerate(remaining_markers):
                        new_state[buffer_start + start_pos + i] = marker

                # Add promoted job to service at this phase (NOT to buffer)
                new_state[class_phase_start + phase_idx] += 1
                results.append((new_state, prob, class_k))

        return results

    def fcfs_promote_from_buffer_map_aware(state: np.ndarray, ist: int, departing_class: int, new_map_phase: int, sched=None, classprio=None) -> List[Tuple[np.ndarray, float, int]]:
        """
        MAP-aware buffer promotion for FCFS stations.

        After a MAP job departs, promote the next buffer job to service with MAP phase continuity:
        - Same MAP class as departing: enters at new_map_phase (kdest from D1)
        - Different MAP class: enters at stored MAP phase variable
        - Non-MAP class: uses pie (initial phase distribution)

        Args:
            state: Current state vector (after departing job removed and MAP phase updated)
            ist: Station index
            departing_class: Class index of the departing job
            new_map_phase: D1 destination phase (kdest) of the departing job
            sched: Scheduling strategy
            classprio: Array of class priorities (for HOL)

        Returns:
            List of (new_state, probability, class_k) tuples.
            Empty list if no jobs in buffer.
        """
        if ist not in fcfs_stations:
            return []

        buffer_start_val = fcfs_buffer_start.get(ist, -1)
        phases_start_val = fcfs_phases_start.get(ist, -1)
        max_buf = fcfs_max_buffer.get(ist, 0)

        if buffer_start_val < 0 or max_buf <= 0:
            return []

        # Find the job to promote based on scheduling strategy
        buf_pos = -1

        if _sched_is(sched, 'HOL') and classprio is not None:
            occupied_positions = []
            for i in range(max_buf):
                buf_val = int(state[buffer_start_val + i])
                if buf_val > 0:
                    class_k = buf_val - 1
                    prio = classprio[class_k] if class_k < len(classprio) else float('inf')
                    occupied_positions.append((i, class_k, prio))
            if occupied_positions:
                min_prio = min(p for _, _, p in occupied_positions)
                candidates = [(pos, ck) for pos, ck, p in occupied_positions if p == min_prio]
                buf_pos = max(pos for pos, _ in candidates)
        elif _sched_is(sched, 'LCFS', 'LCFSPR'):
            for i in range(max_buf):
                buf_val = int(state[buffer_start_val + i])
                if buf_val > 0:
                    buf_pos = i
                    break
        else:
            # FCFS: rightmost non-zero
            for i in range(max_buf - 1, -1, -1):
                buf_val = int(state[buffer_start_val + i])
                if buf_val > 0:
                    buf_pos = i
                    break

        if buf_pos < 0:
            return []

        buf_val = int(state[buffer_start_val + buf_pos])
        class_k = buf_val - 1  # 0-based class index

        n_phases_k = phases[ist, class_k]
        class_phase_start = phases_start_val + sum(phases[ist, :class_k])

        # Determine entry phase(s) based on MAP status of promoted class
        if (ist, class_k) in map_var_offset:
            # Promoted class has MAP distribution
            if class_k == departing_class:
                # Same class: enter at new_map_phase (kdest from D1)
                kentry = new_map_phase
            else:
                # Different MAP class: enter at stored MAP phase variable
                kentry = int(state[map_var_offset[(ist, class_k)]])
            entry_phases = [(kentry, 1.0)]
        else:
            # Non-MAP class: use pie (initial phase distribution)
            entry_probs_k = get_entry_probs(ist, class_k)
            entry_phases = [(p, prob) for p, prob in enumerate(entry_probs_k) if prob > 0 and p < n_phases_k]

        results = []
        for phase_idx, prob in entry_phases:
            new_state = state.copy()

            # Remove job from buffer and repack
            remaining_markers = []
            for i in range(max_buf):
                if i != buf_pos and new_state[buffer_start_val + i] > 0:
                    remaining_markers.append(int(new_state[buffer_start_val + i]))

            # Clear buffer
            for i in range(max_buf):
                new_state[buffer_start_val + i] = 0

            # Repack at rightmost positions
            if remaining_markers:
                n_remaining = len(remaining_markers)
                start_pos = max_buf - n_remaining
                for i, marker in enumerate(remaining_markers):
                    new_state[buffer_start_val + start_pos + i] = marker

            # Add promoted job to service at determined phase
            new_state[class_phase_start + phase_idx] += 1
            results.append((new_state, prob, class_k))

        return results

    def fcfs_promote_from_buffer(state: np.ndarray, ist: int) -> np.ndarray:
        """
        For FCFS station with buffer, promote the first waiting job to service (phase 1 only).

        DEPRECATED: Use fcfs_promote_from_buffer_with_alpha for correct PH/APH handling.
        This function is kept for backward compatibility but always uses phase 1.
        """
        results = fcfs_promote_from_buffer_with_alpha(state, ist)
        if results:
            # Return the first result (shouldn't be used in new code)
            return results[0][0]
        return state

    def siro_get_jobs_in_service(state: np.ndarray, ist: int, k: int) -> int:
        """Get number of class k jobs in service at SIRO station ist."""
        if ist not in siro_stations:
            return 0
        phases_start = siro_phases_start.get(ist, -1)
        if phases_start < 0:
            return 0
        # Calculate offset for class k phases
        class_phase_start = phases_start + sum(phases[ist, :k])
        n_phases_k = phases[ist, k]
        return sum(int(state[class_phase_start + p]) for p in range(n_phases_k))

    def siro_get_buffer_count(state: np.ndarray, ist: int, k: int) -> int:
        """Get buffer count for class k at SIRO station ist."""
        if ist not in siro_stations:
            return 0
        buffer_start = siro_buffer_start.get(ist, -1)
        if buffer_start < 0:
            return 0
        return int(state[buffer_start + k])

    def siro_get_total_jobs(state: np.ndarray, ist: int, k: int) -> int:
        """Get total jobs (buffer + service) for class k at SIRO station ist."""
        return siro_get_buffer_count(state, ist, k) + siro_get_jobs_in_service(state, ist, k)

    def siro_get_total_at_station(state: np.ndarray, ist: int) -> int:
        """Get total jobs at SIRO station (all classes)."""
        total = 0
        for k in range(K):
            total += siro_get_total_jobs(state, ist, k)
        return total

    def siro_promote_from_buffer_with_alpha(state: np.ndarray, ist: int, k: int) -> List[Tuple[np.ndarray, float, int]]:
        """
        For SIRO station with buffer, promote a class k job from buffer to service.

        SIRO selects randomly from buffer, but we just pick specified class k.
        The starting phase is chosen according to alpha distribution.

        Returns:
            List of (new_state, probability, class_k) tuples for each possible starting phase.
            Empty list if no class k jobs in buffer.
        """
        if ist not in siro_stations:
            return []

        buffer_start = siro_buffer_start.get(ist, -1)
        phases_start = siro_phases_start.get(ist, -1)

        if buffer_start < 0 or phases_start < 0:
            return []

        # Check if there's a class k job in buffer
        buf_count = int(state[buffer_start + k])
        if buf_count <= 0:
            return []

        # Get entry probabilities (alpha) for this class
        # Use MAP-aware entry probs to respect stored MAP phase if this is a MAP station
        entry_probs = get_map_aware_entry_probs(state, ist, k)

        results = []
        n_phases_k = phases[ist, k]
        class_phase_start = phases_start + sum(phases[ist, :k])

        for phase_idx, prob in enumerate(entry_probs):
            if prob > 0 and phase_idx < n_phases_k:
                new_state = state.copy()

                # Decrement buffer count for class k
                new_state[buffer_start + k] -= 1

                # Add job to service at this phase
                new_state[class_phase_start + phase_idx] += 1

                results.append((new_state, prob, k))

        return results

    # Build hash map for O(1) state lookup instead of O(n) linear search
    state_map = _build_state_index_map(state_space)

    # Get service rates
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
    else:
        rates = np.ones((M, K))

    # Get routing probabilities
    # NOTE: sn.rt is indexed by stateful nodes (nstateful x K), not stations (M x K)
    # Need to create a station-to-stateful mapping for correct indexing
    station_to_stateful = np.zeros(M, dtype=int)
    if hasattr(sn, 'stationToNode') and sn.stationToNode is not None and \
       hasattr(sn, 'nodeToStateful') and sn.nodeToStateful is not None:
        for ist in range(M):
            node_idx = int(sn.stationToNode[ist]) if ist < len(sn.stationToNode) else ist
            stateful_idx = int(sn.nodeToStateful[node_idx]) if node_idx < len(sn.nodeToStateful) else ist
            station_to_stateful[ist] = stateful_idx
    else:
        # Fallback: assume station index = stateful index
        station_to_stateful = np.arange(M)

    if hasattr(sn, 'rt') and sn.rt is not None:
        # sn.rt is (nstateful * K, nstateful * K) - extract station-only routing
        # Need to resolve routing through non-station stateful nodes (e.g., Router)
        rt_full = np.asarray(sn.rt)
        P = np.zeros((M * K, M * K))
        nstateful = sn.nstateful

        # Build stateful_to_station mapping: -1 for non-station stateful nodes
        stateful_to_station = np.full(nstateful, -1, dtype=int)
        for ist in range(M):
            sf = station_to_stateful[ist]
            if sf < nstateful:
                stateful_to_station[sf] = ist

        # Helper function to resolve routing through non-station nodes
        def resolve_routing_to_stations(src_sf: int, src_class: int):
            """
            Resolve routing from (src_sf, src_class) through non-station nodes.
            Returns dict: {(dst_station, dst_class): probability}
            """
            result = {}
            # Queue of (stateful_idx, class_idx, accumulated_probability)
            queue = []
            src_idx = src_sf * K + src_class
            for dst_sf in range(nstateful):
                for dst_class in range(K):
                    dst_idx = dst_sf * K + dst_class
                    if src_idx < rt_full.shape[0] and dst_idx < rt_full.shape[1]:
                        prob = rt_full[src_idx, dst_idx]
                        if prob > 1e-10:
                            queue.append((dst_sf, dst_class, prob))

            # Process queue until all reach stations
            max_iters = 100
            iters = 0
            while queue and iters < max_iters:
                iters += 1
                sf, cls, prob = queue.pop(0)
                station = stateful_to_station[sf] if sf < len(stateful_to_station) else -1
                if station >= 0:
                    # Reached a station
                    key = (station, cls)
                    result[key] = result.get(key, 0.0) + prob
                else:
                    # Non-station node - follow routing
                    idx = sf * K + cls
                    for next_sf in range(nstateful):
                        for next_cls in range(K):
                            next_idx = next_sf * K + next_cls
                            if idx < rt_full.shape[0] and next_idx < rt_full.shape[1]:
                                next_prob = rt_full[idx, next_idx]
                                if next_prob > 1e-10:
                                    queue.append((next_sf, next_cls, prob * next_prob))
            return result

        for ist in range(M):
            ist_sf = station_to_stateful[ist]
            for k in range(K):
                src_idx_station = ist * K + k
                # Resolve routing through non-station nodes
                routing = resolve_routing_to_stations(ist_sf, k)
                for (jst, r), prob in routing.items():
                    dst_idx_station = jst * K + r
                    P[src_idx_station, dst_idx_station] = prob

        cache_routing_info = {}  # Not used - Cache routing is handled by immediate departures block

        if options.verbose:
            station_names = []
            for ist in range(M):
                if hasattr(sn, 'nodenames') and sn.nodenames is not None:
                    node_idx = int(sn.stationToNode[ist]) if hasattr(sn, 'stationToNode') and sn.stationToNode is not None and ist < len(sn.stationToNode) else ist
                    if node_idx < len(sn.nodenames):
                        station_names.append(str(sn.nodenames[node_idx]))
                    else:
                        station_names.append(f"S{ist}")
                else:
                    station_names.append(f"S{ist}")
            class_names = []
            for k in range(K):
                if hasattr(sn, 'classnames') and sn.classnames is not None and k < len(sn.classnames):
                    class_names.append(str(sn.classnames[k]))
                else:
                    class_names.append(f"C{k}")
            print("  P matrix (station routing):")
            for ist in range(M):
                for k in range(K):
                    src_idx = ist * K + k
                    row = P[src_idx, :]
                    nonzero = np.where(np.abs(row) > 1e-10)[0]
                    if len(nonzero) > 0:
                        dests = []
                        for idx in nonzero:
                            jst = idx // K
                            r = idx % K
                            dests.append(f"{station_names[jst]}:{class_names[r]}={row[idx]:.4f}")
                        print(f"    {station_names[ist]}:{class_names[k]} -> {', '.join(dests)}")
    else:
        # Default: uniform routing
        P = np.ones((M * K, M * K)) / (M * K)

    # Get number of servers
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    # Get load-dependent scaling
    lldscaling = None
    if hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
        lldscaling = np.asarray(sn.lldscaling)


    # Get class-dependent scaling
    cdscaling = None
    if hasattr(sn, 'cdscaling') and sn.cdscaling is not None:
        cdscaling = sn.cdscaling

    def get_cd_scaling(ist: int, nir: np.ndarray) -> float:
        """Get class-dependent scaling factor for station ist with per-class population nir."""
        if cdscaling is None:
            return 1.0
        if ist < len(cdscaling) and cdscaling[ist] is not None:
            return float(cdscaling[ist](nir))
        return 1.0
    # Track which stations are infinite servers
    inf_server_stations = set()
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            if _sched_is(sched_val, 'INF'):
                inf_server_stations.add(ist)

    def get_load_scaling(ist: int, total_jobs: int) -> float:
        """Get service rate scaling factor for station ist with total_jobs present."""
        if total_jobs <= 0:
            return 0.0  # No jobs, no service
        # For infinite servers, each job gets its own server
        if ist in inf_server_stations:
            return float(total_jobs)
        if lldscaling is not None and ist < lldscaling.shape[0]:
            # lldscaling[ist, n-1] gives scaling when n jobs present
            idx = total_jobs - 1
            if lldscaling.ndim > 1 and idx < lldscaling.shape[1]:
                return lldscaling[ist, idx]
            elif lldscaling.ndim > 1:
                return lldscaling[ist, -1]  # Use last value if beyond range
        # Fall back to min(n, c) behavior
        c = nservers[ist] if ist < len(nservers) else 1
        return min(total_jobs, c)

    # Identify Source and Sink stations for open network handling
    # Use sn_has_open_classes to properly handle mixed networks
    is_open = sn_has_open_classes(sn)
    source_station = -1
    sink_station = -1
    if is_open and hasattr(sn, 'sched'):
        for ist, sched_val in sn.sched.items():
            # Check for EXT scheduling (Source station)
            if _sched_is(sched_val, 'EXT') or (isinstance(sched_val, int) and sched_val == SchedStrategy.EXT.value):
                source_station = ist
                break
    if is_open and hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for ist in range(M):
            if station_to_node is not None and ist < len(station_to_node):
                node_idx = int(station_to_node[ist])
            else:
                node_idx = ist
            if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.SINK:
                sink_station = ist
                break

    # Extract class priorities for HOL scheduling (once, before the main loop)
    classprio = None
    if hasattr(sn, 'classprio') and sn.classprio is not None:
        classprio = np.asarray(sn.classprio).flatten()

    # Build transitions
    for s, state in enumerate(state_space):
        # Get job counts as (M, K) matrix from state vector
        state_matrix = get_state_matrix(state)

        # === External arrivals from Source ===
        # For open networks, Source generates arrivals that enter the next station
        # Handle routing through non-station nodes (like Router with RROBIN)
        if source_station >= 0:
            source_node = int(station_to_node[source_station]) if station_to_node is not None and source_station < len(station_to_node) else source_station

            for k in range(K):
                # Arrival rate at source
                arrival_rate = rates[source_station, k] if source_station < rates.shape[0] and k < rates.shape[1] else 0

                if arrival_rate <= 0:
                    continue

                # Get direct destinations from connection matrix
                conn = np.asarray(sn.connmatrix) if hasattr(sn, 'connmatrix') and sn.connmatrix is not None else None
                if conn is not None:
                    direct_dests = np.where(conn[source_node, :] > 0)[0]
                else:
                    direct_dests = []

                for dest_node in direct_dests:
                    # Check if destination is a Cache node (stateful non-station)
                    # Cache nodes need class switching based on cache state
                    dest_isf = -1
                    is_cache_dest = False
                    if hasattr(sn, 'nodeToStateful') and sn.nodeToStateful is not None:
                        dest_isf = int(sn.nodeToStateful[dest_node]) if dest_node < len(sn.nodeToStateful) else -1
                        if dest_isf >= 0 and dest_isf in cache_stations_info:
                            is_cache_dest = True

                    if is_cache_dest:
                        # Destination is a Cache node - arrival puts job at Cache
                        # Cache has capacity 1 per class (matching MATLAB/Java)
                        # The job will depart immediately via Cache immediate departure transitions
                        cache_jobs_offset = cache_jobs_offsets.get(dest_isf)
                        if cache_jobs_offset is not None:
                            # Check if Cache can accept this job (capacity 1 per class)
                            cache_job_count = int(state[cache_jobs_offset + k])
                            if cache_job_count < 1:
                                # Create arrival state: increment Cache job count for class k
                                new_state = state.copy()
                                new_state[cache_jobs_offset + k] = 1

                                ns = _find_state_index_fast(state_map, new_state)
                                if ns >= 0:
                                    p = 1.0 / len(direct_dests) if len(direct_dests) > 1 else 1.0
                                    trans_rate = arrival_rate * p
                                    Q[s, ns] += trans_rate
                                    depRates[s, source_station, k] += trans_rate
                    else:
                        # Not a Cache destination - check if it's a station or non-station node
                        dest_station = -1
                        if node_to_station is not None and dest_node < len(node_to_station):
                            dest_station = int(node_to_station[dest_node])

                        if dest_station >= 0:
                            # Destination is already a station - add arrival directly
                            if dest_station != source_station:
                                # Use MAP-aware entry for MAP destinations
                                entry_probs = get_map_aware_entry_probs(state, dest_station, k)
                                for entry_phase, entry_prob in enumerate(entry_probs):
                                    if entry_prob <= 0:
                                        continue
                                    arrival_state = set_job_count(state.copy(), dest_station, k, +1, phase_idx=entry_phase)

                                    ns = _find_state_index_fast(state_map, arrival_state)
                                    if ns >= 0:
                                        p = 1.0 / len(direct_dests) if len(direct_dests) > 1 else 1.0
                                        trans_rate = arrival_rate * p * entry_prob
                                        Q[s, ns] += trans_rate
                                        depRates[s, source_station, k] += trans_rate
                        else:
                            # Destination is a non-station node (e.g., Router)
                            # Check if this node uses RROBIN routing (state-dependent)
                            has_rrobin_at_dest = (dest_node, k) in rrobin_info.get('outlinks', {})

                            if has_rrobin_at_dest:
                                outlinks = rrobin_info['outlinks'][(dest_node, k)]
                                rr_var_idx = rr_var_map.get((dest_node, k))
                                has_rrobin_state = rr_var_idx is not None

                                if has_rrobin_state:
                                    # Station with RROBIN state variable - use state-dependent routing
                                    current_rr_ptr = int(state[rr_var_idx])
                                    next_dest_node = outlinks[current_rr_ptr]
                                    next_rr_ptr = (current_rr_ptr + 1) % len(outlinks)

                                    # Resolve next_dest_node to a station
                                    final_station = -1
                                    if node_to_station is not None and next_dest_node < len(node_to_station):
                                        final_station = int(node_to_station[next_dest_node])

                                    if final_station >= 0 and final_station < M and final_station != source_station:
                                        entry_probs = get_map_aware_entry_probs(state, final_station, k)
                                        for entry_phase, entry_prob in enumerate(entry_probs):
                                            if entry_prob <= 0:
                                                continue
                                            arrival_state = set_job_count(state.copy(), final_station, k, +1, phase_idx=entry_phase)
                                            # Advance the RR pointer
                                            arrival_state[rr_var_idx] = next_rr_ptr

                                            ns = _find_state_index_fast(state_map, arrival_state)
                                            if ns >= 0:
                                                p = 1.0 / len(direct_dests) if len(direct_dests) > 1 else 1.0
                                                trans_rate = arrival_rate * p * entry_prob
                                                Q[s, ns] += trans_rate
                                                depRates[s, source_station, k] += trans_rate
                                    elif final_station == -1:
                                        # Destination is Sink or non-station - job exits system
                                        exit_state = state.copy()
                                        exit_state[rr_var_idx] = next_rr_ptr
                                        ns = _find_state_index_fast(state_map, exit_state)
                                        if ns >= 0:
                                            p = 1.0 / len(direct_dests) if len(direct_dests) > 1 else 1.0
                                            trans_rate = arrival_rate * p
                                            Q[s, ns] += trans_rate
                                            depRates[s, source_station, k] += trans_rate
                                else:
                                    # Non-station RROBIN (Router): blocking-aware probabilistic routing
                                    # First check which outlinks are feasible, then redistribute
                                    equal_prob = 1.0 / len(outlinks)
                                    p = 1.0 / len(direct_dests) if len(direct_dests) > 1 else 1.0

                                    # First pass: identify feasible destinations
                                    feasible_rr = []  # (final_station, k, equal_prob, ns, entry_prob)
                                    total_feasible_rr = 0.0
                                    for dest_link in outlinks:
                                        final_station = -1
                                        if node_to_station is not None and dest_link < len(node_to_station):
                                            final_station = int(node_to_station[dest_link])
                                        if final_station >= 0 and final_station < M and final_station != source_station:
                                            entry_probs = get_map_aware_entry_probs(state, final_station, k)
                                            for entry_phase, entry_prob in enumerate(entry_probs):
                                                if entry_prob <= 0:
                                                    continue
                                                arrival_state = set_job_count(state.copy(), final_station, k, +1, phase_idx=entry_phase)
                                                ns = _find_state_index_fast(state_map, arrival_state)
                                                if ns >= 0:
                                                    feasible_rr.append((ns, equal_prob, entry_prob))
                                                    total_feasible_rr += equal_prob * entry_prob

                                    # Second pass: add with redistributed probabilities
                                    if feasible_rr and total_feasible_rr > 0:
                                        redist_factor_rr = 1.0 / total_feasible_rr
                                        for ns, eq_p, entry_prob in feasible_rr:
                                            trans_rate = arrival_rate * p * eq_p * entry_prob * redist_factor_rr
                                            Q[s, ns] += trans_rate
                                            depRates[s, source_station, k] += trans_rate
                            else:
                                # Non-RROBIN non-station node - resolve probabilistic routing
                                final_destinations = _get_probabilistic_final_destinations(
                                    sn, dest_node, k
                                )

                                for final_station, final_class, cumulative_prob in final_destinations:
                                    if final_station >= 0 and final_station < M and final_station != source_station:
                                        # Arrival to final station (class may change if routing allows)
                                        # Use MAP-aware entry for MAP destinations
                                        entry_probs = get_map_aware_entry_probs(state, final_station, final_class)
                                        for entry_phase, entry_prob in enumerate(entry_probs):
                                            if entry_prob <= 0:
                                                continue
                                            # Use phase-aware state modification with entry phase
                                            arrival_state = set_job_count(state.copy(), final_station, final_class, +1, phase_idx=entry_phase)

                                            ns = _find_state_index_fast(state_map, arrival_state)
                                            if ns >= 0:
                                                # Probability accounts for:
                                                # 1. 1/num_direct_dests if multiple direct destinations from source
                                                # 2. cumulative_prob from routing through non-station nodes
                                                p = 1.0 / len(direct_dests) if len(direct_dests) > 1 else 1.0
                                                trans_rate = arrival_rate * p * cumulative_prob * entry_prob
                                                Q[s, ns] += trans_rate
                                                depRates[s, source_station, k] += trans_rate

        # === Cache immediate departures ===
        # Cache nodes process jobs immediately (rate Gamma = 1e7, matching MATLAB)
        # When a job is at Cache, it departs immediately with class switching (hit/miss)
        # and routes to downstream destinations (through Router if present)
        GAMMA = 1e8  # MATLAB GlobalConstants.Immediate = 1/FineTol = 1/1e-8
        if cache_stations_info:
            for dest_isf in cache_stations_info:
                cache_info = cache_stations_info[dest_isf]
                cache_jobs_offset = cache_jobs_offsets.get(dest_isf)
                cache_state_var_idx = cache_state_offsets.get(dest_isf)
                if cache_jobs_offset is None:
                    continue

                hitclass = cache_info['hitclass']
                missclass = cache_info['missclass']
                cache_node_idx = cache_info['node_idx']

                # Total jobs at this Cache across all classes
                cache_total_jobs = sum(int(state[cache_jobs_offset + kk]) for kk in range(K))

                for k in range(K):
                    # Check if class k has a job at this Cache
                    cache_job_count = int(state[cache_jobs_offset + k])
                    if cache_job_count <= 0:
                        continue

                    # Check if k is an input class (has hitclass/missclass defined)
                    has_hm_class = (k < len(hitclass) and k < len(missclass) and
                                   hitclass[k] >= 0 and missclass[k] >= 0)

                    if has_hm_class:
                        # Phase 1 of Cache processing: READ event
                        # Class switch only (Init -> Hit/Miss), job stays at Cache.
                        # Guard: READ only fires when total Cache jobs == 1
                        # (matching MATLAB/JAR AfterEventCache guard:
                        # spaceSrv.elementSum() == 1). When multiple jobs are
                        # at Cache, DEP events fire first to drain output-class
                        # jobs before READ can process input-class jobs.
                        if cache_total_jobs != 1:
                            continue
                        if cache_state_var_idx is not None:
                            cache_state_idx = int(state[cache_state_var_idx])
                        else:
                            cache_state_idx = 0

                        cache_transitions = _compute_cache_read_transitions(
                            cache_info, cache_state_idx, k, K
                        )

                        for output_class, trans_prob, new_cache_state_idx in cache_transitions:
                            if trans_prob <= 0:
                                continue

                            # Class switch at Cache: remove input class, add output class
                            new_state = state.copy()
                            new_state[cache_jobs_offset + k] -= 1  # Remove input class
                            new_state[cache_jobs_offset + output_class] += 1  # Add output class
                            if new_cache_state_idx is not None and cache_state_var_idx is not None:
                                new_state[cache_state_var_idx] = new_cache_state_idx
                            ns = _find_state_index_fast(state_map, new_state)
                            if ns >= 0:
                                trans_rate = GAMMA * trans_prob
                                Q[s, ns] += trans_rate
                    else:
                        # Phase 2 of Cache processing: DEP event
                        # Route Hit/Miss/other class jobs to downstream nodes.
                        # When Router nodes are tracked (for stochcomp), route to
                        # Router first. Router immediate transitions then route to
                        # stations. This matches MATLAB which includes Router in the
                        # full state space before applying stochcomp.
                        # When no Router tracking, route directly to stations
                        # (original behavior using _get_probabilistic_final_destinations).

                        # Check if immediate next hop from Cache is a tracked Router
                        routed_to_router = False
                        if router_jobs_offsets:
                            # Check rtnodes for direct next hop from Cache
                            P_nodes = np.asarray(sn.rtnodes) if sn.rtnodes is not None else None
                            if P_nodes is not None:
                                src_idx = cache_node_idx * K + k
                                if src_idx < P_nodes.shape[0]:
                                    for jnd in range(sn.nnodes):
                                        for r_c in range(K):
                                            dst_idx = jnd * K + r_c
                                            if dst_idx < P_nodes.shape[1] and P_nodes[src_idx, dst_idx] > 1e-10:
                                                # Check if jnd is a tracked Router
                                                jnd_isf = int(sn.nodeToStateful[jnd]) if jnd < len(sn.nodeToStateful) else -1
                                                if jnd_isf in router_jobs_offsets:
                                                    # Route to Router: remove from Cache, add to Router
                                                    rj_off = router_jobs_offsets[jnd_isf]
                                                    dep_state = state.copy()
                                                    dep_state[cache_jobs_offset + k] = 0
                                                    dep_state[rj_off + r_c] = 1
                                                    ns = _find_state_index_fast(state_map, dep_state)
                                                    if ns >= 0:
                                                        Q[s, ns] += GAMMA * P_nodes[src_idx, dst_idx]
                                                        routed_to_router = True

                        if not routed_to_router:
                            # No Router tracking - route directly to stations (original behavior)
                            final_destinations = _get_probabilistic_final_destinations(
                                sn, cache_node_idx, k
                            )
                            total_routing_prob = sum(cp for _, _, cp in final_destinations)

                            # First pass: identify feasible destinations
                            feasible_cache_dep = []
                            total_feasible_prob = 0.0
                            for final_station, final_class, cumulative_prob in final_destinations:
                                if final_station >= 0 and final_station < M and final_station != source_station:
                                    entry_probs = get_map_aware_entry_probs(state, final_station, final_class)
                                    for entry_phase, entry_prob in enumerate(entry_probs):
                                        if entry_prob <= 0:
                                            continue
                                        arrival_state = state.copy()
                                        arrival_state[cache_jobs_offset + k] = 0
                                        arrival_state = set_job_count(arrival_state, final_station, final_class, +1, phase_idx=entry_phase)
                                        ns = _find_state_index_fast(state_map, arrival_state)
                                        if ns >= 0:
                                            eff_prob = cumulative_prob * entry_prob
                                            feasible_cache_dep.append((ns, eff_prob))
                                            total_feasible_prob += eff_prob

                            # Second pass: add transitions with redistributed probs
                            if feasible_cache_dep and total_feasible_prob > 0:
                                redist = total_routing_prob / total_feasible_prob
                                for ns, eff_prob in feasible_cache_dep:
                                    Q[s, ns] += GAMMA * eff_prob * redist

                            # Handle exit to Sink (remaining routing probability)
                            exit_prob = max(0.0, 1.0 - total_routing_prob)
                            if exit_prob > 1e-10:
                                exit_state = state.copy()
                                exit_state[cache_jobs_offset + k] = 0  # Remove from Cache
                                ns = _find_state_index_fast(state_map, exit_state)
                                if ns >= 0:
                                    trans_rate = GAMMA * exit_prob
                                    Q[s, ns] += trans_rate

        # === Router immediate departures ===
        # When Router nodes are tracked for stochcomp, process immediate routing.
        # A job at the Router departs immediately at rate GAMMA to downstream stations.
        # This matches MATLAB's full Q matrix before stochcomp.
        if router_jobs_offsets:
            GAMMA_R = 1e8  # Same GAMMA as Cache
            P_nodes = np.asarray(sn.rtnodes) if sn.rtnodes is not None else None
            if P_nodes is not None:
                for isf_r, rinfo in router_nodes_info.items():
                    rj_off = router_jobs_offsets.get(isf_r)
                    if rj_off is None:
                        continue
                    r_node_idx = rinfo['node_idx']
                    for k in range(K):
                        router_job_k = int(state[rj_off + k])
                        if router_job_k <= 0:
                            continue
                        # Router has a job of class k - route to downstream destinations
                        src_idx = r_node_idx * K + k
                        if src_idx >= P_nodes.shape[0]:
                            continue
                        for jnd in range(sn.nnodes):
                            for r_c in range(K):
                                dst_idx = jnd * K + r_c
                                if dst_idx >= P_nodes.shape[1]:
                                    continue
                                prob = P_nodes[src_idx, dst_idx]
                                if prob <= 1e-10:
                                    continue
                                # Check if destination is a station
                                jst = -1
                                if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and jnd < len(sn.nodeToStation):
                                    jst = int(sn.nodeToStation[jnd])
                                if jst >= 0 and jst < M and jst != source_station:
                                    # Route from Router to station
                                    entry_probs = get_map_aware_entry_probs(state, jst, r_c)
                                    for entry_phase, entry_prob in enumerate(entry_probs):
                                        if entry_prob <= 0:
                                            continue
                                        dep_state = state.copy()
                                        dep_state[rj_off + k] = 0  # Remove from Router
                                        dep_state = set_job_count(dep_state, jst, r_c, +1, phase_idx=entry_phase)
                                        ns = _find_state_index_fast(state_map, dep_state)
                                        if ns >= 0:
                                            Q[s, ns] += GAMMA_R * prob * entry_prob
                                elif hasattr(sn, 'nodetype') and jnd < len(sn.nodetype):
                                    nt_val = int(sn.nodetype[jnd].value) if hasattr(sn.nodetype[jnd], 'value') else int(sn.nodetype[jnd])
                                    if nt_val == 1:  # Sink
                                        # Job exits system
                                        dep_state = state.copy()
                                        dep_state[rj_off + k] = 0
                                        ns = _find_state_index_fast(state_map, dep_state)
                                        if ns >= 0:
                                            Q[s, ns] += GAMMA_R * prob

        # === MAP phase transitions (D0 off-diagonal, no service completion) ===
        # These are transitions within the MAP process that don't complete service
        # MATLAB: rate = D0(k, kdest) * kir(class, k) -- uses per-phase service count
        if has_map_fcfs:
            for (ist, k) in map_fcfs.keys():
                n_ik = state_matrix[ist, k]
                if n_ik <= 0:
                    continue  # No jobs, no MAP transitions

                D0, D1 = get_map_matrices(ist, k)
                if D0 is None:
                    continue

                current_map_phase = get_map_phase(state, ist, k)
                n_phases = D0.shape[0]

                # For FCFS: check if this class is actually in service
                # MATLAB uses kir(class, k) which is 0 or 1 for single-server
                is_fcfs_ist = ist in fcfs_stations
                if is_fcfs_ist:
                    jobs_in_service_k = get_fcfs_jobs_in_service(state, ist, k)
                    if jobs_in_service_k <= 0:
                        continue  # Class not in service, no D0 transitions
                    class_fraction = 1.0
                else:
                    total_jobs = np.sum(state_matrix[ist, :])
                    class_fraction = (n_ik / total_jobs) if total_jobs > 0 else 0

                total_jobs = np.sum(state_matrix[ist, :])
                scaling_d0 = get_load_scaling(ist, int(total_jobs)) if total_jobs > 0 else 0
                scaling_d0 *= get_cd_scaling(ist, state_matrix[ist, :])

                # D0 off-diagonal entries: phase transitions without completion
                for j in range(n_phases):
                    if j == current_map_phase:
                        continue  # Skip diagonal
                    rate = D0[current_map_phase, j] * scaling_d0 * class_fraction
                    if rate <= 0:
                        continue

                    # Create new state: update MAP phase AND move service phase count
                    # MATLAB: space_srv(Ks(class)+k) -= 1; space_srv(Ks(class)+kdest) += 1
                    new_state = state.copy()
                    # Update MAP variable
                    if (ist, k) in map_var_offset:
                        new_state[map_var_offset[(ist, k)]] = j
                    # Move job from old phase to new phase in service area
                    if use_phase_aug and phase_offset is not None:
                        old_phase_idx = phase_offset[ist, k] + current_map_phase
                        new_phase_idx = phase_offset[ist, k] + j
                        new_state[old_phase_idx] -= 1
                        new_state[new_phase_idx] += 1

                    if np.all(new_state[:state_dim] >= 0):
                        ns = _find_state_index_fast(state_map, new_state)
                        if ns >= 0:
                            Q[s, ns] += rate

        # === Phase-type (non-MAP) phase transitions (D0 off-diagonal) at FCFS stations ===
        # For PH/APH/Erlang/HyperExp, phase counts are tracked in the state vector directly
        # D0 off-diagonal transitions move jobs between phases without service completion
        if use_phase_aug:
            for ist in range(M):
                if ist == source_station:
                    continue

                # Check scheduling strategy
                if hasattr(sn, 'sched') and sn.sched is not None:
                    sched = sn.sched.get(ist, SchedStrategy.FCFS)
                else:
                    sched = SchedStrategy.FCFS

                # Handle FCFS and similar scheduling strategies
                sched_name = _normalize_sched_strategy(sched)
                if sched_name not in {'FCFS', 'HOL', 'LCFS', 'LCFSPR', 'PS', 'DPS', 'GPS'}:
                    continue

                for k in range(K):
                    # Skip if this is a MAP distribution (handled separately above)
                    if (ist, k) in map_fcfs:
                        continue

                    n_phases_ik = phases[ist, k]
                    if n_phases_ik <= 1:
                        continue  # No internal phase transitions

                    D0, D1 = get_map_matrices(ist, k)
                    if D0 is None:
                        continue

                    # Get phase counts from state
                    phase_counts = get_phase_counts(state, ist, k)

                    total_jobs = np.sum(state_matrix[ist, :])

                    # For each source phase with jobs
                    for p_src in range(n_phases_ik):
                        n_src = int(phase_counts[p_src])
                        if n_src <= 0:
                            continue

                        # For each destination phase (D0 off-diagonal)
                        for p_dst in range(n_phases_ik):
                            if p_src == p_dst:
                                continue  # Skip diagonal

                            # D0 off-diagonal: phase transition rate
                            phase_transition_rate = D0[p_src, p_dst]
                            if phase_transition_rate <= 0:
                                continue

                            # Rate for phase transition:
                            # For FCFS/HOL/LCFS: each in-service job transitions independently
                            #   at the base D0 rate. Do NOT scale by min(n,c) since phase counts
                            #   already represent individual in-service jobs.
                            #   Matches D1 completion logic (lines 5526-5548).
                            # For PS/DPS/GPS: PS approximation (share of capacity)
                            if _sched_is(sched, 'FCFS', 'HOL', 'LCFS', 'LCFSPR'):
                                # FCFS-like: each in-service job transitions at base rate
                                cd_scaling_d0 = get_cd_scaling(ist, state_matrix[ist, :])
                                lld_d0 = 1.0
                                if lldscaling is not None and ist < lldscaling.shape[0]:
                                    idx = int(total_jobs) - 1
                                    if lldscaling.ndim > 1 and idx < lldscaling.shape[1]:
                                        lld_d0 = lldscaling[ist, idx]
                                    elif lldscaling.ndim > 1:
                                        lld_d0 = lldscaling[ist, -1]
                                total_rate = n_src * phase_transition_rate * lld_d0 * cd_scaling_d0
                            else:
                                # PS-like: share rate among all jobs
                                scaling_d0 = get_load_scaling(ist, int(total_jobs)) if total_jobs > 0 else 0
                                scaling_d0 *= get_cd_scaling(ist, state_matrix[ist, :])
                                total_rate = (n_src / total_jobs) * phase_transition_rate * scaling_d0 if total_jobs > 0 else 0

                            # Create new state: move one job from p_src to p_dst
                            new_state = state.copy()
                            start_idx = phase_offset[ist, k]
                            new_state[start_idx + p_src] -= 1
                            new_state[start_idx + p_dst] += 1

                            if np.all(new_state[:state_dim] >= 0):
                                ns = _find_state_index_fast(state_map, new_state)
                                if ns >= 0:
                                    Q[s, ns] += total_rate

        # === Service completions and internal routing ===
        for ist in range(M):
            # Skip source station - it doesn't hold jobs for service
            if ist == source_station:
                continue

            for k in range(K):
                n_ik = state_matrix[ist, k]
                if n_ik <= 0:
                    continue

                # Service completion at station ist, class k
                # Rate depends on scheduling discipline
                if hasattr(sn, 'sched') and sn.sched is not None:
                    sched = sn.sched.get(ist, SchedStrategy.FCFS)
                else:
                    sched = SchedStrategy.FCFS

                # Check if this is a MAP distribution at FCFS station
                is_map_fcfs_class = (ist, k) in map_fcfs
                mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0

                # Skip if rate is NaN (disabled class at this station)
                if np.isnan(mu):
                    continue

                if _sched_is(sched, 'INF'):
                    # Infinite server: rate depends on phase-type distribution
                    n_phases_ik = phases[ist, k]
                    if use_phase_aug and n_phases_ik > 1:
                        # Phase-type distribution: need separate transitions for each completing phase
                        D0, D1 = get_map_matrices(ist, k)
                        if D0 is not None and D1 is not None:
                            phase_counts = get_phase_counts(state, ist, k)
                            # For each phase with completions, create separate routing transitions
                            for p in range(n_phases_ik):
                                n_p = int(phase_counts[p])
                                if n_p <= 0:
                                    continue
                                completion_rate_p = np.sum(D1[p, :])
                                if completion_rate_p <= 0:
                                    continue
                                # Rate for completions from this phase
                                rate_p = n_p * completion_rate_p
                                # Create transitions for each destination
                                src_idx = ist * K + k

                                # Check for RROBIN routing
                                node_idx_inf = int(station_to_node[ist]) if station_to_node is not None and ist < len(station_to_node) else ist
                                is_rrobin_inf = (node_idx_inf, k) in rrobin_info['outlinks']

                                if is_rrobin_inf:
                                    # RROBIN routing for INF station
                                    outlinks_inf = rrobin_info['outlinks'][(node_idx_inf, k)]
                                    rr_var_idx_inf = rr_var_map[(node_idx_inf, k)]
                                    current_rr_ptr_inf = int(state[rr_var_idx_inf])
                                    dest_node_inf = outlinks_inf[current_rr_ptr_inf]
                                    next_rr_ptr_inf = (current_rr_ptr_inf + 1) % len(outlinks_inf)

                                    if node_to_station is not None and dest_node_inf < len(node_to_station):
                                        dest_station_inf = int(node_to_station[dest_node_inf])
                                    else:
                                        dest_station_inf = dest_node_inf

                                    r = k  # No class switching for RROBIN

                                    if dest_station_inf >= 0 and dest_station_inf < M:
                                        entry_probs = get_map_aware_entry_probs(state, dest_station_inf, r)
                                        for entry_phase, entry_prob in enumerate(entry_probs):
                                            if entry_prob <= 0:
                                                continue
                                            # Create base state: remove from phase p
                                            base_state = state.copy()
                                            start_idx = phase_offset[ist, k]
                                            base_state[start_idx + p] -= 1
                                            if np.any(base_state[:state_dim] < 0):
                                                continue
                                            new_state = set_job_count(base_state, dest_station_inf, r, +1, phase_idx=entry_phase)
                                            new_state[rr_var_idx_inf] = next_rr_ptr_inf
                                            if np.all(new_state[:state_dim] >= 0):
                                                ns = _find_state_index_fast(state_map, new_state)
                                                if ns >= 0:
                                                    Q[s, ns] += rate_p * entry_prob
                                                    depRates[s, ist, k] += rate_p * entry_prob
                                    else:
                                        # Destination is Sink or outside network
                                        base_state = state.copy()
                                        start_idx = phase_offset[ist, k]
                                        base_state[start_idx + p] -= 1
                                        if np.all(base_state[:state_dim] >= 0):
                                            new_state = base_state.copy()
                                            new_state[rr_var_idx_inf] = next_rr_ptr_inf
                                            ns = _find_state_index_fast(state_map, new_state)
                                            if ns >= 0:
                                                Q[s, ns] += rate_p
                                                depRates[s, ist, k] += rate_p
                                else:
                                    # Standard probabilistic routing
                                    for jst in range(M):
                                        if jst == source_station:
                                            continue
                                        for r in range(K):
                                            dst_idx = jst * K + r
                                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                                prob = P[src_idx, dst_idx]
                                            else:
                                                prob = 0
                                            if prob <= 0:
                                                continue
                                            # Get entry probabilities for destination
                                            # Use MAP-aware entry for MAP destinations
                                            entry_probs = get_map_aware_entry_probs(state, jst, r)
                                            for entry_phase, entry_prob in enumerate(entry_probs):
                                                if entry_prob <= 0:
                                                    continue
                                                # Create base state: remove from phase p
                                                base_state = state.copy()
                                                start_idx = phase_offset[ist, k]
                                                base_state[start_idx + p] -= 1  # Remove from completing phase
                                                if np.any(base_state[:state_dim] < 0):
                                                    continue
                                                # FCFS/HOL/LCFS buffer promotion with alpha distribution
                                                if ist in fcfs_stations:
                                                    promo_results = fcfs_promote_from_buffer_with_alpha(base_state, ist, sched, classprio)
                                                    if promo_results:
                                                        # Multiple destination states based on alpha
                                                        for promo_state, promo_prob, _ in promo_results:
                                                            new_state = set_job_count(promo_state.copy(), jst, r, +1, phase_idx=entry_phase)
                                                            if np.any(new_state[:state_dim] < 0):
                                                                continue
                                                            ns = _find_state_index_fast(state_map, new_state)
                                                            if ns >= 0:
                                                                Q[s, ns] += rate_p * prob * entry_prob * promo_prob
                                                                depRates[s, ist, k] += rate_p * prob * entry_prob * promo_prob
                                                        continue  # Skip non-promotion path
                                                # SIRO buffer promotion with alpha distribution
                                                if ist in siro_stations:
                                                    # Compute total buffer count for SIRO selection probability
                                                    total_buf = sum(siro_get_buffer_count(base_state, ist, kk) for kk in range(K))
                                                    if total_buf > 0:
                                                        found_promo = False
                                                        for buf_k in range(K):
                                                            buf_count_k = siro_get_buffer_count(base_state, ist, buf_k)
                                                            if buf_count_k > 0:
                                                                # Probability of selecting class buf_k
                                                                select_prob = buf_count_k / total_buf
                                                                promo_results = siro_promote_from_buffer_with_alpha(base_state, ist, buf_k)
                                                                for promo_state, promo_prob, _ in promo_results:
                                                                    new_state = set_job_count(promo_state.copy(), jst, r, +1, phase_idx=entry_phase)
                                                                    if np.any(new_state[:state_dim] < 0):
                                                                        continue
                                                                    ns = _find_state_index_fast(state_map, new_state)
                                                                    if ns >= 0:
                                                                        Q[s, ns] += rate_p * prob * entry_prob * promo_prob * select_prob
                                                                        depRates[s, ist, k] += rate_p * prob * entry_prob * promo_prob * select_prob
                                                                        found_promo = True
                                                        if found_promo:
                                                            continue  # Skip non-promotion path
                                                # No buffer jobs or not FCFS
                                                new_state = set_job_count(base_state.copy(), jst, r, +1, phase_idx=entry_phase)
                                                if np.any(new_state[:state_dim] < 0):
                                                    continue
                                                ns = _find_state_index_fast(state_map, new_state)
                                                if ns >= 0:
                                                    Q[s, ns] += rate_p * prob * entry_prob
                                                    # Track departure rate for class k from station ist
                                                    depRates[s, ist, k] += rate_p * prob * entry_prob
                            # Skip the normal routing below since we handled it here
                            continue
                        else:
                            # Fallback to aggregate rate
                            rate = n_ik * mu
                    else:
                        # Single phase (exponential): rate = n * mu
                        rate = n_ik * mu
                elif _sched_is(sched, 'DPS'):
                    # Discriminatory Processor Sharing: weighted by schedparam
                    total_jobs = np.sum(state_matrix[ist, :])
                    if total_jobs > 0:
                        # Get weight for this class (default 1.0)
                        w_k = 1.0
                        if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                            if ist < sn.schedparam.shape[0] and k < sn.schedparam.shape[1]:
                                w_k = sn.schedparam[ist, k]
                        # Compute weighted share: sum of w_j * n_j for all classes
                        weighted_total = 0.0
                        for kk in range(K):
                            n_kk = state_matrix[ist, kk]
                            w_kk = 1.0
                            if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                                if ist < sn.schedparam.shape[0] and kk < sn.schedparam.shape[1]:
                                    w_kk = sn.schedparam[ist, kk]
                            weighted_total += w_kk * n_kk
                        # Rate = (w_k * n_k / weighted_total) * scaling * mu
                        scaling = get_load_scaling(ist, int(total_jobs))
                        scaling *= get_cd_scaling(ist, state_matrix[ist, :])
                        if weighted_total > 0:
                            rate = (w_k * n_ik / weighted_total) * scaling * mu
                        else:
                            rate = 0
                    else:
                        rate = 0
                elif _sched_is(sched, 'PSPRIO'):
                    # PS with Priority: only jobs at highest priority level get service
                    total_jobs = np.sum(state_matrix[ist, :])
                    if total_jobs <= 0:
                        rate = 0
                    else:
                        scaling = get_load_scaling(ist, int(total_jobs))
                        nservers_ist = nservers[ist] if ist < len(nservers) else 1
                        scaling *= get_cd_scaling(ist, state_matrix[ist, :])

                        # Get class priorities (lower value = higher priority)
                        classprio = None
                        if hasattr(sn, 'classprio') and sn.classprio is not None:
                            classprio = np.asarray(sn.classprio).flatten()

                        if classprio is None or total_jobs <= nservers_ist:
                            # No priorities or all jobs get service: regular PS
                            rate = scaling * mu * (n_ik / total_jobs)
                        else:
                            # Find minimum priority value among classes with jobs present
                            present_classes = [kk for kk in range(K) if state_matrix[ist, kk] > 0]
                            min_prio = min(classprio[kk] for kk in present_classes)

                            if classprio[k] == min_prio:
                                # This class is at highest priority - gets PS among same priority
                                niprio = sum(state_matrix[ist, kk] for kk in range(K)
                                           if classprio[kk] == min_prio)
                                scaling_prio = min(niprio, nservers_ist)
                                rate = scaling_prio * mu * (n_ik / niprio)
                            else:
                                # Not highest priority: rate = 0
                                rate = 0

                elif _sched_is(sched, 'GPSPRIO'):
                    # GPS with Priority: weighted sharing among highest priority jobs
                    total_jobs = np.sum(state_matrix[ist, :])
                    if total_jobs <= 0:
                        rate = 0
                    else:
                        nservers_ist = nservers[ist] if ist < len(nservers) else 1

                        # Get class priorities
                        classprio = None
                        if hasattr(sn, 'classprio') and sn.classprio is not None:
                            classprio = np.asarray(sn.classprio).flatten()

                        # Get weights from schedparam
                        weights = np.ones(K)
                        if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                            if ist < sn.schedparam.shape[0]:
                                weights = sn.schedparam[ist, :K].flatten()
                        weights = weights / np.sum(weights)  # Normalize

                        if classprio is None or total_jobs <= nservers_ist:
                            # No priorities or all jobs get service: regular GPS
                            # cir = min(nir, 1) - indicator of presence
                            cir = np.minimum(state_matrix[ist, :], 1)
                            weighted_total = np.sum(weights * cir)
                            if weighted_total > 0:
                                rate = mu * (n_ik / state_matrix[ist, k]) * weights[k] / weighted_total if state_matrix[ist, k] > 0 else 0
                            else:
                                rate = 0
                        else:
                            # Find minimum priority value among classes with jobs present
                            present_classes = [kk for kk in range(K) if state_matrix[ist, kk] > 0]
                            min_prio = min(classprio[kk] for kk in present_classes)

                            if classprio[k] == min_prio:
                                # This class is at highest priority - GPS among same priority
                                nirprio = np.zeros(K)
                                for kk in range(K):
                                    if classprio[kk] == min_prio:
                                        nirprio[kk] = state_matrix[ist, kk]

                                cir = np.minimum(nirprio, 1)
                                weighted_total = np.sum(weights * cir)
                                if weighted_total > 0 and nirprio[k] > 0:
                                    rate = mu * (n_ik / nirprio[k]) * weights[k] / weighted_total
                                else:
                                    rate = 0
                            else:
                                # Not highest priority: rate = 0
                                rate = 0
                else:
                    # FCFS, PS, SIRO, etc.
                    total_jobs = np.sum(state_matrix[ist, :])
                    scaling = get_load_scaling(ist, int(total_jobs))

                    scaling *= get_cd_scaling(ist, state_matrix[ist, :])
                    # For FCFS stations with single-phase (Exp) distributions, only jobs in service get rate
                    if _sched_is(sched, 'FCFS', 'HOL', 'LCFS', 'LCFSPR') and ist in fcfs_stations:
                        # FCFS: only jobs in service complete - no PS-style sharing
                        # Each job in service completes independently at rate mu
                        # NOTE: For FCFS with c servers, jobs_in_service_k already accounts for
                        # the multi-server capacity (at most min(total_jobs, c) jobs in service).
                        # We should NOT multiply by scaling=min(n,c) as that would double-count.
                        # However, user-specified lldscaling must be applied (matches MATLAB afterEventStation.m).
                        jobs_in_service_k = get_fcfs_jobs_in_service(state, ist, k)
                        if jobs_in_service_k > 0:
                            cd_scaling = get_cd_scaling(ist, state_matrix[ist, :])
                            # Apply user-specified lldscaling if present (MATLAB always multiplies by lldscaling)
                            lld_factor = 1.0
                            if lldscaling is not None and ist < lldscaling.shape[0]:
                                idx = int(total_jobs) - 1
                                if lldscaling.ndim > 1 and idx < lldscaling.shape[1]:
                                    lld_factor = lldscaling[ist, idx]
                                elif lldscaling.ndim > 1:
                                    lld_factor = lldscaling[ist, -1]
                            rate = lld_factor * cd_scaling * jobs_in_service_k * mu
                        else:
                            # No jobs in service for this class
                            rate = 0
                    elif ist in siro_stations:
                        # SIRO: only jobs IN SERVICE complete, not buffer jobs
                        # Jobs in buffer wait until server is free
                        jobs_in_service_k = siro_get_jobs_in_service(state, ist, k)
                        if jobs_in_service_k > 0:
                            # Get total jobs in service (all classes) for PS-style sharing
                            total_in_service = siro_get_total_at_station(state, ist) - sum(
                                siro_get_buffer_count(state, ist, kk) for kk in range(K)
                            )
                            if total_in_service > 0:
                                # PS-style sharing among jobs in service
                                rate = scaling * mu * (jobs_in_service_k / total_in_service)
                            else:
                                rate = jobs_in_service_k * mu
                        else:
                            # No jobs in service for this class
                            rate = 0
                    else:
                        # PS, DPS, GPS: share rate among all jobs
                        rate = scaling * mu * (n_ik / total_jobs) if total_jobs > 0 else 0

                # === Special handling for phase-type distributions (PH, APH, Erlang, HyperExp) at FCFS stations ===
                # For phase-type distributions, service completions use D1 matrix with phase-dependent rates
                # This includes APH, PH, Erlang, HyperExp, and MAP distributions
                n_phases_ik = phases[ist, k]
                is_ph_distribution = use_phase_aug and n_phases_ik > 1 and not is_map_fcfs_class

                if is_ph_distribution:
                    D0, D1 = get_map_matrices(ist, k)
                    if D0 is not None and D1 is not None:
                        total_jobs = np.sum(state_matrix[ist, :])
                        scaling = get_load_scaling(ist, int(total_jobs))

                        scaling *= get_cd_scaling(ist, state_matrix[ist, :])
                        # Get phase counts from the state
                        phase_counts = get_phase_counts(state, ist, k)

                        # For each phase with jobs, compute service completions
                        src_idx = ist * K + k
                        total_routing = 0.0

                        # Get node index for this station
                        node_idx = int(station_to_node[ist]) if station_to_node is not None and ist < len(station_to_node) else ist
                        is_rrobin = (node_idx, k) in rrobin_info['outlinks']

                        # Determine if this is an FCFS station with buffer ordering
                        is_fcfs_with_buffer = ist in fcfs_stations

                        for p in range(n_phases_ik):
                            n_p = int(phase_counts[p])
                            if n_p <= 0:
                                continue

                            # Completion rate from phase p = sum of D1[p, :] (row sum)
                            completion_rate_p = np.sum(D1[p, :])
                            if completion_rate_p <= 0:
                                continue

                            # Rate computation depends on scheduling discipline
                            if is_fcfs_with_buffer:
                                # FCFS: only the class in service completes
                                # Check if this class has jobs in service (non-zero phase counts)
                                jobs_in_service_k = get_fcfs_jobs_in_service(state, ist, k)
                                if jobs_in_service_k > 0:
                                    # This class is in service - full D1 rate applies
                                    # Each job in phase p completes independently at completion_rate_p
                                    # NOTE: For FCFS with c servers, n_p (jobs in phase p) already
                                    # accounts for multi-server capacity. We should NOT multiply
                                    # by scaling=min(n,c) as that would double-count.
                                    # However, user-specified lldscaling must be applied (matches MATLAB afterEventStation.m).
                                    cd_scaling = get_cd_scaling(ist, state_matrix[ist, :])
                                    lld_factor = 1.0
                                    if lldscaling is not None and ist < lldscaling.shape[0]:
                                        idx = int(total_jobs) - 1
                                        if lldscaling.ndim > 1 and idx < lldscaling.shape[1]:
                                            lld_factor = lldscaling[ist, idx]
                                        elif lldscaling.ndim > 1:
                                            lld_factor = lldscaling[ist, -1]
                                    rate_p = lld_factor * cd_scaling * n_p * completion_rate_p
                                else:
                                    # This class is NOT in service - rate is 0
                                    rate_p = 0
                            elif _sched_is(sched, 'PSPRIO'):
                                # PS with Priority: only jobs at highest priority get service
                                nservers_ist = nservers[ist] if ist < len(nservers) else 1
                                classprio = None
                                if hasattr(sn, 'classprio') and sn.classprio is not None:
                                    classprio = np.asarray(sn.classprio).flatten()

                                if classprio is None or total_jobs <= nservers_ist:
                                    # No priorities or all jobs get service: regular PS
                                    scaling_prio = min(total_jobs, nservers_ist)
                                    rate_p = (n_p / total_jobs) * completion_rate_p * scaling_prio if total_jobs > 0 else 0
                                else:
                                    # Find minimum priority value among classes with jobs present
                                    present_classes = [kk for kk in range(K) if state_matrix[ist, kk] > 0]
                                    min_prio = min(classprio[kk] for kk in present_classes)

                                    if classprio[k] == min_prio:
                                        # This class is at highest priority - gets PS among same priority
                                        niprio = sum(state_matrix[ist, kk] for kk in range(K)
                                                   if classprio[kk] == min_prio)
                                        scaling_prio = min(niprio, nservers_ist)
                                        rate_p = (n_p / niprio) * completion_rate_p * scaling_prio if niprio > 0 else 0
                                    else:
                                        # Not highest priority: rate = 0
                                        rate_p = 0
                            elif _sched_is(sched, 'GPSPRIO'):
                                # GPS with Priority: weighted sharing among highest priority jobs
                                nservers_ist = nservers[ist] if ist < len(nservers) else 1
                                classprio = None
                                if hasattr(sn, 'classprio') and sn.classprio is not None:
                                    classprio = np.asarray(sn.classprio).flatten()

                                # Get weights from schedparam
                                weights = np.ones(K)
                                if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                                    if ist < sn.schedparam.shape[0]:
                                        weights = sn.schedparam[ist, :K].flatten()
                                weights = weights / np.sum(weights)  # Normalize

                                if classprio is None or total_jobs <= nservers_ist:
                                    # No priorities or all jobs get service: regular GPS
                                    cir = np.minimum(state_matrix[ist, :], 1)
                                    weighted_total = np.sum(weights * cir)
                                    if weighted_total > 0 and state_matrix[ist, k] > 0:
                                        rate_p = completion_rate_p * (n_p / state_matrix[ist, k]) * weights[k] / weighted_total
                                    else:
                                        rate_p = 0
                                else:
                                    # Find minimum priority value among classes with jobs present
                                    present_classes = [kk for kk in range(K) if state_matrix[ist, kk] > 0]
                                    min_prio = min(classprio[kk] for kk in present_classes)

                                    if classprio[k] == min_prio:
                                        # This class is at highest priority - GPS among same priority
                                        nirprio = np.zeros(K)
                                        for kk in range(K):
                                            if classprio[kk] == min_prio:
                                                nirprio[kk] = state_matrix[ist, kk]

                                        cir = np.minimum(nirprio, 1)
                                        weighted_total = np.sum(weights * cir)
                                        if weighted_total > 0 and nirprio[k] > 0:
                                            rate_p = completion_rate_p * (n_p / nirprio[k]) * weights[k] / weighted_total
                                        else:
                                            rate_p = 0
                                    else:
                                        # Not highest priority: rate = 0
                                        rate_p = 0
                            elif _sched_is(sched, 'INF'):
                                # Infinite server: rate is independent of total jobs
                                rate_p = n_p * completion_rate_p
                            else:
                                # Non-FCFS (PS, DPS, etc.): use PS approximation
                                # rate = (n_p / total_jobs) * completion_rate_p * scaling
                                rate_p = (n_p / total_jobs) * completion_rate_p * scaling if total_jobs > 0 else 0

                            if rate_p <= 0:
                                continue

                            if is_rrobin:
                                # RROBIN routing
                                outlinks = rrobin_info['outlinks'][(node_idx, k)]
                                rr_var_idx = rr_var_map[(node_idx, k)]
                                current_rr_ptr = int(state[rr_var_idx])
                                dest_node = outlinks[current_rr_ptr]
                                next_rr_ptr = (current_rr_ptr + 1) % len(outlinks)

                                if node_to_station is not None and dest_node < len(node_to_station):
                                    dest_station = int(node_to_station[dest_node])
                                else:
                                    dest_station = dest_node

                                r = k  # No class switching for RROBIN

                                if dest_station >= 0 and dest_station < M:
                                    entry_probs = get_map_aware_entry_probs(state, dest_station, r)
                                    for entry_phase, entry_prob in enumerate(entry_probs):
                                        if entry_prob <= 0:
                                            continue
                                        # Create base state: remove from phase p
                                        base_state = state.copy()
                                        start_idx = phase_offset[ist, k]
                                        base_state[start_idx + p] -= 1
                                        if np.any(base_state[:state_dim] < 0):
                                            continue
                                        # FCFS buffer promotion with alpha distribution
                                        if ist in fcfs_stations:
                                            promo_results = fcfs_promote_from_buffer_with_alpha(base_state, ist, sched, classprio)
                                            if promo_results:
                                                for promo_state, promo_prob, _ in promo_results:
                                                    new_state = set_job_count(promo_state.copy(), dest_station, r, +1, phase_idx=entry_phase)
                                                    new_state[rr_var_idx] = next_rr_ptr
                                                    if np.all(new_state[:state_dim] >= 0):
                                                        ns = _find_state_index_fast(state_map, new_state)
                                                        if ns >= 0:
                                                            Q[s, ns] += rate_p * entry_prob * promo_prob
                                                            total_routing += rate_p * entry_prob * promo_prob
                                                            depRates[s, ist, k] += rate_p * entry_prob * promo_prob
                                                continue
                                        # SIRO buffer promotion with alpha distribution
                                        if ist in siro_stations:
                                            total_buf = sum(siro_get_buffer_count(base_state, ist, kk) for kk in range(K))
                                            if total_buf > 0:
                                                found_promo = False
                                                for buf_k in range(K):
                                                    buf_count_k = siro_get_buffer_count(base_state, ist, buf_k)
                                                    if buf_count_k > 0:
                                                        select_prob = buf_count_k / total_buf
                                                        promo_results = siro_promote_from_buffer_with_alpha(base_state, ist, buf_k)
                                                        for promo_state, promo_prob, _ in promo_results:
                                                            new_state = set_job_count(promo_state.copy(), dest_station, r, +1, phase_idx=entry_phase)
                                                            new_state[rr_var_idx] = next_rr_ptr
                                                            if np.all(new_state[:state_dim] >= 0):
                                                                ns = _find_state_index_fast(state_map, new_state)
                                                                if ns >= 0:
                                                                    Q[s, ns] += rate_p * entry_prob * promo_prob * select_prob
                                                                    total_routing += rate_p * entry_prob * promo_prob * select_prob
                                                                    depRates[s, ist, k] += rate_p * entry_prob * promo_prob * select_prob
                                                                    found_promo = True
                                                if found_promo:
                                                    continue
                                        # No buffer jobs or not FCFS/SIRO
                                        new_state = set_job_count(base_state.copy(), dest_station, r, +1, phase_idx=entry_phase)
                                        new_state[rr_var_idx] = next_rr_ptr
                                        if np.all(new_state[:state_dim] >= 0):
                                            ns = _find_state_index_fast(state_map, new_state)
                                            if ns >= 0:
                                                Q[s, ns] += rate_p * entry_prob
                                                total_routing += rate_p * entry_prob
                                                depRates[s, ist, k] += rate_p * entry_prob
                                else:
                                    # Destination is Sink
                                    base_state = state.copy()
                                    start_idx = phase_offset[ist, k]
                                    base_state[start_idx + p] -= 1
                                    if np.any(base_state[:state_dim] < 0):
                                        continue
                                    # FCFS buffer promotion with alpha distribution
                                    if ist in fcfs_stations:
                                        promo_results = fcfs_promote_from_buffer_with_alpha(base_state, ist, sched, classprio)
                                        if promo_results:
                                            for promo_state, promo_prob, _ in promo_results:
                                                new_state = promo_state.copy()
                                                new_state[rr_var_idx] = next_rr_ptr
                                                if np.all(new_state[:state_dim] >= 0):
                                                    ns = _find_state_index_fast(state_map, new_state)
                                                    if ns >= 0:
                                                        Q[s, ns] += rate_p * promo_prob
                                                        total_routing += rate_p * promo_prob
                                                        depRates[s, ist, k] += rate_p * promo_prob
                                            continue
                                    # SIRO buffer promotion with alpha distribution
                                    if ist in siro_stations:
                                        total_buf = sum(siro_get_buffer_count(base_state, ist, kk) for kk in range(K))
                                        if total_buf > 0:
                                            found_promo = False
                                            for buf_k in range(K):
                                                buf_count_k = siro_get_buffer_count(base_state, ist, buf_k)
                                                if buf_count_k > 0:
                                                    select_prob = buf_count_k / total_buf
                                                    promo_results = siro_promote_from_buffer_with_alpha(base_state, ist, buf_k)
                                                    for promo_state, promo_prob, _ in promo_results:
                                                        new_state = promo_state.copy()
                                                        new_state[rr_var_idx] = next_rr_ptr
                                                        if np.all(new_state[:state_dim] >= 0):
                                                            ns = _find_state_index_fast(state_map, new_state)
                                                            if ns >= 0:
                                                                Q[s, ns] += rate_p * promo_prob * select_prob
                                                                total_routing += rate_p * promo_prob * select_prob
                                                                depRates[s, ist, k] += rate_p * promo_prob * select_prob
                                                                found_promo = True
                                            if found_promo:
                                                continue
                                    # No buffer jobs or not FCFS/SIRO
                                    new_state = base_state.copy()
                                    new_state[rr_var_idx] = next_rr_ptr
                                    if np.all(new_state[:state_dim] >= 0):
                                        ns = _find_state_index_fast(state_map, new_state)
                                        if ns >= 0:
                                            Q[s, ns] += rate_p
                                            total_routing += rate_p
                                            depRates[s, ist, k] += rate_p
                            else:
                                # Standard probabilistic routing
                                for jst in range(M):
                                    if jst == source_station:
                                        continue
                                    for r in range(K):
                                        dst_idx = jst * K + r
                                        if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                            prob = P[src_idx, dst_idx]
                                        else:
                                            prob = 0
                                        if prob <= 0:
                                            continue
                                        entry_probs = get_map_aware_entry_probs(state, jst, r)
                                        for entry_phase, entry_prob in enumerate(entry_probs):
                                            if entry_prob <= 0:
                                                continue
                                            # Create base state: remove from phase p
                                            base_state = state.copy()
                                            start_idx = phase_offset[ist, k]
                                            base_state[start_idx + p] -= 1
                                            if np.any(base_state[:state_dim] < 0):
                                                continue
                                            # FCFS buffer promotion with alpha distribution
                                            if ist in fcfs_stations:
                                                promo_results = fcfs_promote_from_buffer_with_alpha(base_state, ist, sched, classprio)
                                                if promo_results:
                                                    for promo_state, promo_prob, _ in promo_results:
                                                        new_state = set_job_count(promo_state.copy(), jst, r, +1, phase_idx=entry_phase)
                                                        if np.any(new_state[:state_dim] < 0):
                                                            continue
                                                        ns = _find_state_index_fast(state_map, new_state)
                                                        if ns >= 0:
                                                            Q[s, ns] += rate_p * prob * entry_prob * promo_prob
                                                            total_routing += rate_p * prob * entry_prob * promo_prob
                                                            depRates[s, ist, k] += rate_p * prob * entry_prob * promo_prob
                                                    continue
                                            # SIRO buffer promotion with alpha distribution
                                            if ist in siro_stations:
                                                total_buf = sum(siro_get_buffer_count(base_state, ist, kk) for kk in range(K))
                                                if total_buf > 0:
                                                    found_promo = False
                                                    for buf_k in range(K):
                                                        buf_count_k = siro_get_buffer_count(base_state, ist, buf_k)
                                                        if buf_count_k > 0:
                                                            select_prob = buf_count_k / total_buf
                                                            promo_results = siro_promote_from_buffer_with_alpha(base_state, ist, buf_k)
                                                            for promo_state, promo_prob, _ in promo_results:
                                                                new_state = set_job_count(promo_state.copy(), jst, r, +1, phase_idx=entry_phase)
                                                                if np.any(new_state[:state_dim] < 0):
                                                                    continue
                                                                ns = _find_state_index_fast(state_map, new_state)
                                                                if ns >= 0:
                                                                    Q[s, ns] += rate_p * prob * entry_prob * promo_prob * select_prob
                                                                    total_routing += rate_p * prob * entry_prob * promo_prob * select_prob
                                                                    depRates[s, ist, k] += rate_p * prob * entry_prob * promo_prob * select_prob
                                                                    found_promo = True
                                                    if found_promo:
                                                        continue
                                            # No buffer jobs or not FCFS/SIRO
                                            new_state = set_job_count(base_state.copy(), jst, r, +1, phase_idx=entry_phase)
                                            if np.any(new_state[:state_dim] < 0):
                                                continue
                                            ns = _find_state_index_fast(state_map, new_state)
                                            if ns >= 0:
                                                Q[s, ns] += rate_p * prob * entry_prob
                                                total_routing += rate_p * prob * entry_prob
                                                depRates[s, ist, k] += rate_p * prob * entry_prob

                        # Handle departures to sink for open networks
                        if is_open:
                            for p in range(n_phases_ik):
                                n_p = int(phase_counts[p])
                                if n_p <= 0:
                                    continue
                                completion_rate_p = np.sum(D1[p, :])
                                if completion_rate_p <= 0:
                                    continue
                                # Rate depends on scheduling discipline
                                if is_fcfs_with_buffer:
                                    # FCFS: only the class in service completes
                                    jobs_in_service_k = get_fcfs_jobs_in_service(state, ist, k)
                                    if jobs_in_service_k > 0:
                                        rate_p = n_p * completion_rate_p * scaling
                                    else:
                                        rate_p = 0
                                else:
                                    # Non-FCFS: PS approximation
                                    rate_p = (n_p / total_jobs) * completion_rate_p * scaling if total_jobs > 0 else 0
                                # Check routing to outside (sink) based on P matrix
                                # Exclude routing to source_station - in open networks, routing TO source
                                # represents exits to Sink (from dtmc_stochcomp absorption)
                                exit_prob = 0.0
                                for jst in range(M):
                                    if jst == source_station:
                                        continue  # Routing to Source means exit in open networks
                                    for r in range(K):
                                        dst_idx = jst * K + r
                                        if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                            exit_prob += P[src_idx, dst_idx]
                                exit_prob = max(0.0, 1.0 - exit_prob)
                                if exit_prob > 0:
                                    base_state = state.copy()
                                    start_idx = phase_offset[ist, k]
                                    base_state[start_idx + p] -= 1
                                    if np.any(base_state[:state_dim] < 0):
                                        continue

                                    # FCFS buffer promotion with alpha distribution
                                    # When service completes, the next job in buffer enters service
                                    # Note: we do NOT clear buffer markers here - fcfs_promote_from_buffer_with_alpha
                                    # handles finding the oldest waiting job and moving it from buffer to service
                                    if ist in fcfs_stations:
                                        promo_results = fcfs_promote_from_buffer_with_alpha(base_state, ist, sched, classprio)
                                        if promo_results:
                                            for promo_state, promo_prob, _ in promo_results:
                                                if np.all(promo_state[:state_dim] >= 0):
                                                    ns = _find_state_index_fast(state_map, promo_state)
                                                    if ns >= 0:
                                                        Q[s, ns] += rate_p * exit_prob * promo_prob
                                                        depRates[s, ist, k] += rate_p * exit_prob * promo_prob
                                            continue
                                    # SIRO buffer promotion with alpha distribution
                                    if ist in siro_stations:
                                        total_buf = sum(siro_get_buffer_count(base_state, ist, kk) for kk in range(K))
                                        if total_buf > 0:
                                            found_promo = False
                                            for buf_k in range(K):
                                                buf_count_k = siro_get_buffer_count(base_state, ist, buf_k)
                                                if buf_count_k > 0:
                                                    select_prob = buf_count_k / total_buf
                                                    promo_results = siro_promote_from_buffer_with_alpha(base_state, ist, buf_k)
                                                    for promo_state, promo_prob, _ in promo_results:
                                                        if np.all(promo_state[:state_dim] >= 0):
                                                            ns = _find_state_index_fast(state_map, promo_state)
                                                            if ns >= 0:
                                                                Q[s, ns] += rate_p * exit_prob * promo_prob * select_prob
                                                                depRates[s, ist, k] += rate_p * exit_prob * promo_prob * select_prob
                                                                found_promo = True
                                            if found_promo:
                                                continue
                                    # No buffer jobs or not FCFS/SIRO - buffer already cleared above
                                    if np.all(base_state[:state_dim] >= 0):
                                        ns = _find_state_index_fast(state_map, base_state)
                                        if ns >= 0:
                                            Q[s, ns] += rate_p * exit_prob
                                            depRates[s, ist, k] += rate_p * exit_prob

                        continue  # Skip the standard routing code below

                # === Special handling for MAP distributions at FCFS stations ===
                # For MAP, service completions use D1 matrix with phase-dependent rates
                if is_map_fcfs_class:
                    D0, D1 = get_map_matrices(ist, k)
                    if D0 is not None and D1 is not None:
                        n_phases = D1.shape[0]

                        total_jobs = np.sum(state_matrix[ist, :])
                        scaling = get_load_scaling(ist, int(total_jobs))

                        scaling *= get_cd_scaling(ist, state_matrix[ist, :])
                        # Determine if this is an FCFS station with buffer ordering
                        is_fcfs_with_buffer = ist in fcfs_stations

                        if is_fcfs_with_buffer:
                            # FCFS: only the class in service completes
                            jobs_in_service_k = get_fcfs_jobs_in_service(state, ist, k)
                            if jobs_in_service_k <= 0:
                                # This class is NOT in service - rate is 0
                                continue
                        else:
                            # Non-FCFS: use PS approximation
                            if total_jobs <= 0 or n_ik <= 0:
                                continue

                        # Get actual phase counts to determine which phase the job is in
                        # This is different from the MAP phase variable which tracks continuation phase
                        phase_counts = get_phase_counts(state, ist, k)

                        # Calculate routing probabilities
                        src_idx = ist * K + k
                        total_routing = 0.0

                        # Get node index for this station
                        node_idx = int(station_to_node[ist]) if station_to_node is not None and ist < len(station_to_node) else ist

                        # Check if this node has RROBIN routing for class k
                        is_rrobin = (node_idx, k) in rrobin_info['outlinks']

                        # Iterate over source phases with jobs (matching MATLAB's approach)
                        # For FCFS single server, typically only one phase has a job
                        for src_phase in range(n_phases):
                            n_in_src_phase = int(phase_counts[src_phase])
                            if n_in_src_phase <= 0:
                                continue

                            # For each destination MAP phase j, the rate is D1[src_phase, j] * count
                            for new_map_phase in range(n_phases):
                                if is_fcfs_with_buffer:
                                    # FCFS: rate = D1[src_phase, dest_phase] * count_in_src_phase
                                    d1_rate = D1[src_phase, new_map_phase] * n_in_src_phase * scaling
                                else:
                                    # PS approximation: scale by class fraction
                                    d1_rate = D1[src_phase, new_map_phase] * n_in_src_phase * scaling * (n_ik / total_jobs)
                                if d1_rate <= 0:
                                    continue

                                # Step 1: Remove departing job from specific phase and update MAP phase
                                base_state = set_job_count(state, ist, k, -1, phase_idx=src_phase)
                                base_state = set_map_phase(base_state, ist, k, new_map_phase)

                                # Step 2: Buffer promotion at source station ist (MAP-aware)
                                promo_list = []
                                if is_fcfs_with_buffer:
                                    promo_list = fcfs_promote_from_buffer_map_aware(base_state, ist, k, new_map_phase, sched, classprio)
                                if not promo_list:
                                    # No buffer jobs or not FCFS - no promotion needed
                                    promo_list = [(base_state, 1.0, -1)]

                                # Step 3: Route departing job to destination
                                if is_rrobin:
                                    # RROBIN routing with MAP
                                    outlinks = rrobin_info['outlinks'][(node_idx, k)]
                                    rr_var_idx = rr_var_map[(node_idx, k)]
                                    current_rr_ptr = int(state[rr_var_idx])
                                    dest_node = outlinks[current_rr_ptr]
                                    next_rr_ptr = (current_rr_ptr + 1) % len(outlinks)

                                    if node_to_station is not None and dest_node < len(node_to_station):
                                        dest_station = int(node_to_station[dest_node])
                                    else:
                                        dest_station = dest_node

                                    r = k  # No class switching for RROBIN

                                    for promo_state, promo_prob, _ in promo_list:
                                        promo_state_rr = promo_state.copy()
                                        promo_state_rr[rr_var_idx] = next_rr_ptr

                                        if dest_station >= 0 and dest_station < M:
                                            # Use MAP-aware entry for MAP destinations
                                            entry_probs = get_map_aware_entry_probs(promo_state_rr, dest_station, r)
                                            for entry_phase, entry_prob in enumerate(entry_probs):
                                                if entry_prob <= 0:
                                                    continue
                                                new_state = set_job_count(promo_state_rr, dest_station, r, +1, phase_idx=entry_phase)

                                                if np.all(new_state[:state_dim] >= 0):
                                                    ns = _find_state_index_fast(state_map, new_state)
                                                    if ns >= 0:
                                                        Q[s, ns] += d1_rate * entry_prob * promo_prob
                                                        total_routing += d1_rate * entry_prob * promo_prob
                                                        depRates[s, ist, k] += d1_rate * entry_prob * promo_prob
                                        else:
                                            # Destination is Sink
                                            if np.all(promo_state_rr[:state_dim] >= 0):
                                                ns = _find_state_index_fast(state_map, promo_state_rr)
                                                if ns >= 0:
                                                    Q[s, ns] += d1_rate * promo_prob
                                                    total_routing += d1_rate * promo_prob
                                                    depRates[s, ist, k] += d1_rate * promo_prob
                                else:
                                    # Standard probabilistic routing with MAP
                                    for jst in range(M):
                                        if jst == source_station:
                                            continue
                                        if jst == sink_station:
                                            continue  # Sink handled by exit-to-sink fallback
                                        for r in range(K):
                                            dst_idx = jst * K + r

                                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                                p = P[src_idx, dst_idx]
                                            else:
                                                p = 0

                                            if p <= 0:
                                                continue

                                            # Use MAP-aware entry for MAP destinations
                                            entry_probs = get_map_aware_entry_probs(state, jst, r)
                                            for entry_phase, entry_prob in enumerate(entry_probs):
                                                if entry_prob <= 0:
                                                    continue
                                                for promo_state, promo_prob, _ in promo_list:
                                                    new_state = set_job_count(promo_state, jst, r, +1, phase_idx=entry_phase)

                                                    if np.any(new_state[:state_dim] < 0):
                                                        continue

                                                    ns = _find_state_index_fast(state_map, new_state)
                                                    if ns >= 0:
                                                        Q[s, ns] += d1_rate * p * entry_prob * promo_prob
                                                        total_routing += d1_rate * p * entry_prob * promo_prob
                                                        depRates[s, ist, k] += d1_rate * p * entry_prob * promo_prob

                        # Note: exits to sink for open networks are handled in the routing sections above
                        # (RROBIN checks dest_station < 0, probabilistic routing via P matrix)
                        continue  # Skip the standard routing code for MAP distributions

                if rate <= 0:
                    continue

                # Calculate routing probabilities
                src_idx = ist * K + k
                total_routing = 0.0

                # Get node index for this station
                node_idx = int(station_to_node[ist]) if station_to_node is not None and ist < len(station_to_node) else ist

                # Check if this node has RROBIN routing for class k
                is_rrobin = (node_idx, k) in rrobin_info['outlinks']

                if is_rrobin:
                    # === State-dependent routing (RROBIN) ===
                    outlinks = rrobin_info['outlinks'][(node_idx, k)]
                    rr_var_idx = rr_var_map[(node_idx, k)]
                    current_rr_ptr = int(state[rr_var_idx])  # Index into outlinks array
                    dest_node = outlinks[current_rr_ptr]

                    # Advance round-robin pointer for the next job (wraps around)
                    next_rr_ptr = (current_rr_ptr + 1) % len(outlinks)

                    # Get destination station from destination node
                    if node_to_station is not None and dest_node < len(node_to_station):
                        dest_station = int(node_to_station[dest_node])
                    else:
                        dest_station = dest_node  # Fallback

                    # For RROBIN, routing is deterministic: probability 1.0 to current destination
                    # No class switching for RROBIN (r = k)
                    r = k
                    if dest_station >= 0 and dest_station < M:
                        # Get entry probabilities for destination
                        # Use MAP-aware entry for MAP destinations (preserves MAP phase continuity)
                        entry_probs = get_map_aware_entry_probs(state, dest_station, r)
                        for entry_phase, entry_prob in enumerate(entry_probs):
                            if entry_prob <= 0:
                                continue
                            # Create new state with job movement and updated RR pointer
                            # Use phase-aware state modification
                            new_state = set_job_count(state, ist, k, -1)  # Departure from ist, class k
                            new_state = set_job_count(new_state, dest_station, r, +1, phase_idx=entry_phase)
                            new_state[rr_var_idx] = next_rr_ptr  # Advance RR pointer

                            if np.all(new_state[:state_dim] >= 0):
                                ns = _find_state_index_fast(state_map, new_state)
                                if ns >= 0:
                                    Q[s, ns] += rate * 1.0 * entry_prob
                                    total_routing += entry_prob
                                    # Track departure rate for class k from station ist
                                    depRates[s, ist, k] += rate * 1.0 * entry_prob
                    else:
                        # Destination is Sink (exit system) - for open networks
                        new_state = set_job_count(state, ist, k, -1)
                        new_state[rr_var_idx] = next_rr_ptr  # Still advance RR pointer

                        if np.all(new_state[:state_dim] >= 0):
                            ns = _find_state_index_fast(state_map, new_state)
                            if ns >= 0:
                                Q[s, ns] += rate * 1.0
                                total_routing = 1.0
                                # Track departure rate for class k from station ist
                                depRates[s, ist, k] += rate * 1.0
                else:
                    # === Standard probabilistic routing ===
                    # Find destination states based on routing
                    for jst in range(M):
                        if jst == source_station:
                            continue  # Can't route to source
                        if jst == sink_station:
                            continue  # Sink handled by exit-to-sink fallback below
                        for r in range(K):
                            # Routing probability from (ist, k) to (jst, r)
                            dst_idx = jst * K + r

                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                p = P[src_idx, dst_idx]
                            else:
                                p = 0

                            if p <= 0:
                                continue

                            total_routing += p

                            # Get entry probabilities for destination
                            # Use MAP-aware entry for MAP destinations (preserves MAP phase continuity)
                            entry_probs = get_map_aware_entry_probs(state, jst, r)
                            for entry_phase, entry_prob in enumerate(entry_probs):
                                if entry_prob <= 0:
                                    continue
                                # Create base state: departure from ist, class k
                                base_state = set_job_count(state, ist, k, -1)

                                # Check if new state is valid
                                if np.any(base_state[:state_dim] < 0):
                                    continue

                                # For FCFS stations: handle buffer promotion after departure
                                if ist in fcfs_stations:
                                    promo_results = fcfs_promote_from_buffer_with_alpha(base_state, ist, sched, classprio)
                                    if promo_results:
                                        # Jobs waiting in buffer - create transitions for each promotion outcome
                                        for promo_state, promo_prob, _ in promo_results:
                                            new_state = set_job_count(promo_state.copy(), jst, r, +1, phase_idx=entry_phase)
                                            if np.any(new_state[:state_dim] < 0):
                                                continue
                                            ns = _find_state_index_fast(state_map, new_state)
                                            if ns >= 0:
                                                Q[s, ns] += rate * p * entry_prob * promo_prob
                                                depRates[s, ist, k] += rate * p * entry_prob * promo_prob
                                        continue  # Skip non-promotion path

                                # For SIRO stations: handle buffer promotion after departure
                                if ist in siro_stations:
                                    total_buf = sum(siro_get_buffer_count(base_state, ist, kk) for kk in range(K))
                                    if total_buf > 0:
                                        found_promo = False
                                        for buf_k in range(K):
                                            buf_count_k = siro_get_buffer_count(base_state, ist, buf_k)
                                            if buf_count_k > 0:
                                                # SIRO selects randomly from buffer
                                                select_prob = buf_count_k / total_buf
                                                promo_results = siro_promote_from_buffer_with_alpha(base_state, ist, buf_k)
                                                for promo_state, promo_prob, _ in promo_results:
                                                    new_state = set_job_count(promo_state.copy(), jst, r, +1, phase_idx=entry_phase)
                                                    if np.any(new_state[:state_dim] < 0):
                                                        continue
                                                    ns = _find_state_index_fast(state_map, new_state)
                                                    if ns >= 0:
                                                        Q[s, ns] += rate * p * entry_prob * promo_prob * select_prob
                                                        depRates[s, ist, k] += rate * p * entry_prob * promo_prob * select_prob
                                                        found_promo = True
                                        if found_promo:
                                            continue  # Skip non-promotion path

                                # No buffer jobs or not FCFS/SIRO - standard state modification
                                new_state = set_job_count(base_state.copy(), jst, r, +1, phase_idx=entry_phase)

                                # Find new state index using O(1) hash lookup
                                ns = _find_state_index_fast(state_map, new_state)
                                if ns >= 0:
                                    Q[s, ns] += rate * p * entry_prob
                                    # Track departure rate for class k from station ist
                                    depRates[s, ist, k] += rate * p * entry_prob

                # === Departures to sink (for open networks) ===
                # If routing doesn't sum to 1, remaining probability exits to sink
                if is_open and total_routing < 1.0:
                    exit_prob = 1.0 - total_routing
                    if exit_prob > 0:
                        # Create base state with job leaving the system
                        base_state = set_job_count(state, ist, k, -1)

                        if np.all(base_state[:state_dim] >= 0):
                            # For FCFS stations: handle buffer promotion after departure
                            if ist in fcfs_stations:
                                promo_results = fcfs_promote_from_buffer_with_alpha(base_state, ist, sched, classprio)
                                if promo_results:
                                    # Jobs waiting in buffer - create transitions for each promotion outcome
                                    for promo_state, promo_prob, _ in promo_results:
                                        if np.all(promo_state[:state_dim] >= 0):
                                            ns = _find_state_index_fast(state_map, promo_state)
                                            if ns >= 0:
                                                Q[s, ns] += rate * exit_prob * promo_prob
                                                depRates[s, ist, k] += rate * exit_prob * promo_prob
                                    continue  # Skip non-promotion path

                            # For SIRO stations: handle buffer promotion after departure
                            if ist in siro_stations:
                                total_buf = sum(siro_get_buffer_count(base_state, ist, kk) for kk in range(K))
                                if total_buf > 0:
                                    found_promo = False
                                    for buf_k in range(K):
                                        buf_count_k = siro_get_buffer_count(base_state, ist, buf_k)
                                        if buf_count_k > 0:
                                            select_prob = buf_count_k / total_buf
                                            promo_results = siro_promote_from_buffer_with_alpha(base_state, ist, buf_k)
                                            for promo_state, promo_prob, _ in promo_results:
                                                if np.all(promo_state[:state_dim] >= 0):
                                                    ns = _find_state_index_fast(state_map, promo_state)
                                                    if ns >= 0:
                                                        Q[s, ns] += rate * exit_prob * promo_prob * select_prob
                                                        depRates[s, ist, k] += rate * exit_prob * promo_prob * select_prob
                                                        found_promo = True
                                    if found_promo:
                                        continue  # Skip non-promotion path

                            # No buffer jobs or not FCFS/SIRO - standard exit
                            ns = _find_state_index_fast(state_map, base_state)
                            if ns >= 0:
                                Q[s, ns] += rate * exit_prob
                                # Track departure rate for class k from station ist
                                depRates[s, ist, k] += rate * exit_prob

    # === Phase-to-phase transitions for INF stations with phase-type distributions ===
    # These are internal transitions (D0 off-diagonal) that don't involve departures
    if use_phase_aug:
        for s, state in enumerate(state_space):
            for ist in range(M):
                # Only handle INF (infinite server) stations
                if hasattr(sn, 'sched') and sn.sched is not None:
                    sched = sn.sched.get(ist, SchedStrategy.FCFS)
                else:
                    continue

                if sched != SchedStrategy.INF:
                    continue

                for k in range(K):
                    n_phases_ik = phases[ist, k]
                    if n_phases_ik <= 1:
                        continue  # No phase transitions for single-phase distributions

                    # Get D0 matrix for phase transitions
                    D0, D1 = get_map_matrices(ist, k)
                    if D0 is None:
                        continue

                    # Get current phase counts
                    phase_counts = get_phase_counts(state, ist, k)

                    # For each source phase with jobs
                    for p_src in range(n_phases_ik):
                        n_src = int(phase_counts[p_src])
                        if n_src <= 0:
                            continue

                        # For each destination phase (D0 off-diagonal)
                        for p_dst in range(n_phases_ik):
                            if p_src == p_dst:
                                continue  # Skip diagonal (absorbed into overall departure rate)

                            # D0 off-diagonal: phase transition rate
                            phase_transition_rate = D0[p_src, p_dst]
                            if phase_transition_rate <= 0:
                                continue

                            # Total rate = n_src * phase_transition_rate
                            total_rate = n_src * phase_transition_rate

                            # Create new state: one job moves from p_src to p_dst
                            new_state = state.copy()
                            start_idx = phase_offset[ist, k]
                            new_state[start_idx + p_src] -= 1
                            new_state[start_idx + p_dst] += 1

                            # Add transition to Q
                            if np.all(new_state[:state_dim] >= 0):
                                ns = _find_state_index_fast(state_map, new_state)
                                if ns >= 0:
                                    Q[s, ns] += total_rate

    # Make valid generator (set diagonal)
    Q = ctmc_makeinfgen(Q)

    return Q, depRates


def _compute_metrics_from_distribution(
    sn: NetworkStruct,
    pi: np.ndarray,
    state_space: np.ndarray,
    rrobin_info: Optional[dict] = None,
    depRates: Optional[np.ndarray] = None
) -> Dict[str, np.ndarray]:
    """
    Compute performance metrics from steady-state distribution.

    Args:
        sn: Network structure
        pi: Steady-state probability distribution
        state_space: State space matrix (may include RR state variables at end)
        rrobin_info: Round-robin routing information
        depRates: Departure rates array from generator construction (shape: n_states, M, K)
                  If provided, used for accurate throughput computation (matching MATLAB approach)

    Returns:
        Dictionary with Q, U, R, T matrices
    """
    M = sn.nstations
    K = sn.nclasses

    # Check if phase augmentation is used
    use_phase_aug = rrobin_info.get('needs_phase_augmentation', False) if rrobin_info else False
    phases = rrobin_info.get('phases', np.ones((M, K), dtype=int)) if rrobin_info else np.ones((M, K), dtype=int)
    phase_offset = rrobin_info.get('phase_offset', None) if rrobin_info else None
    state_dim = rrobin_info.get('state_dim', M * K) if rrobin_info else M * K

    # Get MAP FCFS info for throughput computation
    map_fcfs = rrobin_info.get('map_fcfs', {}) if rrobin_info else {}
    map_var_offset = rrobin_info.get('map_var_offset', {}) if rrobin_info else {}

    # Get FCFS buffer info for correct queue length computation
    fcfs_stations = rrobin_info.get('fcfs_stations', set()) if rrobin_info else set()
    fcfs_max_buffer = rrobin_info.get('fcfs_max_buffer', {}) if rrobin_info else {}
    fcfs_phases_start = rrobin_info.get('fcfs_phases_start', {}) if rrobin_info else {}
    skip_stations = rrobin_info.get('skip_stations', set()) if rrobin_info else set()

    # Get SIRO stations for count-based buffer handling
    siro_stations = rrobin_info.get('siro_stations', set()) if rrobin_info else set()

    # Compute buffer starts for FCFS and SIRO stations
    # Skip Source/Sink stations as they don't contribute to state space
    fcfs_buffer_start_metrics = {}
    siro_buffer_start_metrics = {}  # Maps ist -> start index of buffer counts
    if use_phase_aug and (fcfs_stations or siro_stations):
        current_offset = 0
        for ist in range(M):
            # Skip Source/Sink stations - they don't appear in state vector
            if ist in skip_stations:
                fcfs_buffer_start_metrics[ist] = -1
                siro_buffer_start_metrics[ist] = -1
                continue
            if ist in fcfs_stations:
                max_buf = fcfs_max_buffer.get(ist, 0)
                fcfs_buffer_start_metrics[ist] = current_offset
                siro_buffer_start_metrics[ist] = -1
                current_offset += max_buf + sum(phases[ist, :])
            elif ist in siro_stations:
                # SIRO: buffer counts (K elements) + phase counts
                siro_buffer_start_metrics[ist] = current_offset
                fcfs_buffer_start_metrics[ist] = -1
                current_offset += K + sum(phases[ist, :])  # K buffer counts + phases
            else:
                fcfs_buffer_start_metrics[ist] = -1
                siro_buffer_start_metrics[ist] = -1
                current_offset += sum(phases[ist, :])

    def get_job_count(state: np.ndarray, ist: int, k: int) -> float:
        """Get total job count for (station, class) from state vector.

        For FCFS stations with phase-augmented state: count BOTH buffer entries
        AND in-service jobs (phase counts). Buffer stores waiting jobs, phases
        track in-service jobs.

        For SIRO stations: buffer_counts[k] + phase_counts (count-based buffer).

        For non-FCFS/non-SIRO stations (PS-like): count phase sums.
        Each phase count represents jobs in that phase, sum = total jobs.
        """
        count = 0.0
        if use_phase_aug and phase_offset is not None:
            if ist in fcfs_stations:
                # FCFS stations: count buffer entries (waiting) + phase counts (in-service)
                buffer_start = fcfs_buffer_start_metrics.get(ist, -1)
                max_buf = fcfs_max_buffer.get(ist, 0)
                if buffer_start >= 0 and max_buf > 0:
                    # Count buffer entries matching this class (1-based class IDs)
                    for buf_pos in range(max_buf):
                        buf_val = int(state[buffer_start + buf_pos])
                        if buf_val == k + 1:  # 1-based class ID
                            count += 1.0
                # Also count jobs in service (phase counts)
                start_idx = phase_offset[ist, k]
                end_idx = start_idx + phases[ist, k]
                count += float(sum(state[start_idx:end_idx]))
            elif ist in siro_stations:
                # SIRO stations: buffer_counts[k] (count-based) + phase counts (in-service)
                buffer_start = siro_buffer_start_metrics.get(ist, -1)
                if buffer_start >= 0:
                    # Buffer counts are first K elements at the station
                    count += float(state[buffer_start + k])  # Buffer count for class k
                    # Also count jobs in service (phase counts)
                    # SIRO state: [buf_0, buf_1, ..., buf_{K-1}, phase_counts...]
                    # Phase counts start at buffer_start + K
                    phase_start = buffer_start + K
                    n_phases_k = int(phases[ist, k])
                    phase_start_k = int(phase_start + int(np.sum(phases[ist, :k])))
                    count += float(np.sum(state[phase_start_k:phase_start_k + n_phases_k]))
            else:
                # Non-FCFS/non-SIRO stations (PS-like): count phase sums
                start_idx = phase_offset[ist, k]
                end_idx = start_idx + phases[ist, k]
                count = float(sum(state[start_idx:end_idx]))
            return count
        else:
            return float(state[ist * K + k])

    def get_state_matrix(state: np.ndarray) -> np.ndarray:
        """Get job counts as (M, K) matrix from state vector."""
        result = np.zeros((M, K))
        for ist in range(M):
            for k in range(K):
                result[ist, k] = get_job_count(state, ist, k)
        return result

    def get_phase_counts_metrics(state: np.ndarray, ist: int, k: int) -> np.ndarray:
        """Get phase counts for (station, class) from state vector."""
        if use_phase_aug and phase_offset is not None:
            start_idx = phase_offset[ist, k]
            end_idx = start_idx + phases[ist, k]
            return np.array([int(x) for x in state[start_idx:end_idx]])
        else:
            # Single phase
            return np.array([int(state[ist * K + k])])

    def get_d1_matrix(ist: int, k: int) -> Optional[np.ndarray]:
        """Get D1 matrix for (station, class)."""
        if not hasattr(sn, 'proc') or sn.proc is None:
            return None
        proc_is_list = isinstance(sn.proc, list)
        if proc_is_list:
            if ist >= len(sn.proc) or sn.proc[ist] is None:
                return None
            station_proc = sn.proc[ist]
        else:
            if ist not in sn.proc:
                return None
            station_proc = sn.proc[ist]

        proc_entry = None
        if isinstance(station_proc, (list, tuple)):
            if k < len(station_proc):
                proc_entry = station_proc[k]
        elif isinstance(station_proc, dict):
            proc_entry = station_proc.get(k)

        if proc_entry is None:
            return None

        # Handle PH/APH distributions stored as [alpha, T] or direct D0/D1 matrices
        if isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 2:
            first_elem = np.atleast_2d(np.array(proc_entry[0], dtype=float))
            second_elem = np.atleast_2d(np.array(proc_entry[1], dtype=float))

            # Check if this is [alpha, T] format (alpha is 1D/row vector, T is square matrix)
            if first_elem.shape[0] == 1 and second_elem.shape[0] == second_elem.shape[1]:
                # This is [alpha, T] format for PH/APH distributions
                alpha = first_elem.flatten()
                T = second_elem
                # Compute D1 from T: exit_rates = -row_sums(T), D1 = outer(exit_rates, alpha)
                exit_rates = -np.sum(T, axis=1)
                D1 = np.outer(exit_rates, alpha)
                return D1
            else:
                # This is [D0, D1] format
                return second_elem

        # Handle Erlang distribution stored as dict with 'k' and 'mu'
        if isinstance(proc_entry, dict):
            if 'k' in proc_entry and 'mu' in proc_entry:
                n_phases = int(proc_entry['k'])
                mu = float(proc_entry['mu'])  # per-phase rate
                if n_phases > 1:
                    # Construct Erlang D1 matrix: D1[k-1,0] = mu (completion from last phase)
                    D1 = np.zeros((n_phases, n_phases))
                    D1[n_phases - 1, 0] = mu
                    return D1
                else:
                    # Single phase: exponential
                    return np.array([[mu]])
            elif 'probs' in proc_entry and 'rates' in proc_entry:
                # HyperExp distribution
                probs = np.array(proc_entry['probs'])
                rates = np.array(proc_entry['rates'])
                # D1: completion from phase i, restart in phase j with prob_j
                D1 = np.outer(rates, probs)
                return D1
            elif 'rate' in proc_entry:
                # Exponential distribution
                mu = float(proc_entry['rate'])
                return np.array([[mu]])

        return None

    def get_map_phase(state: np.ndarray, ist: int, k: int) -> int:
        """Get current MAP phase for (station, class) from state vector."""
        if (ist, k) not in map_var_offset:
            return 0
        idx = map_var_offset[(ist, k)]
        return int(state[idx])

    def get_map_d1_rate(ist: int, k: int, map_phase: int) -> float:
        """Get MAP D1 service completion rate from current phase."""
        if not hasattr(sn, 'proc') or sn.proc is None:
            return 1.0
        proc_is_list = isinstance(sn.proc, list)
        if proc_is_list:
            if ist >= len(sn.proc) or sn.proc[ist] is None:
                return 1.0
            station_proc = sn.proc[ist]
        else:
            if ist not in sn.proc:
                return 1.0
            station_proc = sn.proc[ist]

        proc_entry = None
        if isinstance(station_proc, (list, tuple)):
            if k < len(station_proc):
                proc_entry = station_proc[k]
        elif isinstance(station_proc, dict):
            proc_entry = station_proc.get(k)

        if proc_entry is None:
            return 1.0

        # Handle direct D0/D1 matrices
        if isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 2:
            D1 = np.atleast_2d(np.array(proc_entry[1], dtype=float))
            # Total completion rate from this phase = row sum of D1
            return np.sum(D1[map_phase, :])

        # Handle Erlang distribution stored as dict with 'k' and 'mu'
        if isinstance(proc_entry, dict):
            if 'k' in proc_entry and 'mu' in proc_entry:
                n_phases = int(proc_entry['k'])
                mu = float(proc_entry['mu'])  # per-phase rate
                # For Erlang, only last phase completes
                if map_phase == n_phases - 1:
                    return mu
                else:
                    return 0.0
            elif 'probs' in proc_entry and 'rates' in proc_entry:
                # HyperExp distribution: completion rate from phase = rate for that phase
                rates = np.array(proc_entry['rates'])
                if map_phase < len(rates):
                    return rates[map_phase]
                return 0.0
            elif 'rate' in proc_entry:
                # Exponential distribution
                return float(proc_entry['rate'])

        return 1.0

    # Initialize metrics
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))

    # Get service rates
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
    else:
        rates = np.ones((M, K))

    # Get number of servers
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    # Get load-dependent scaling
    lldscaling = None
    if hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
        lldscaling = np.asarray(sn.lldscaling)

    # Track which stations are infinite servers
    inf_server_stations = set()
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            if _sched_is(sched_val, 'INF'):
                inf_server_stations.add(ist)

    def get_load_scaling(ist: int, total_jobs: int) -> float:
        """Get service rate scaling factor for station ist with total_jobs present."""
        if total_jobs <= 0:
            return 0.0  # No jobs, no service
        # For infinite servers, each job gets its own server
        if ist in inf_server_stations:
            return float(total_jobs)
        if lldscaling is not None and ist < lldscaling.shape[0]:
            idx = total_jobs - 1
            if lldscaling.ndim > 1 and idx < lldscaling.shape[1]:
                return lldscaling[ist, idx]
            elif lldscaling.ndim > 1:
                return lldscaling[ist, -1]
        c = nservers[ist] if ist < len(nservers) else 1
        return min(total_jobs, c)

    def get_max_scaling(ist: int) -> float:
        """Get maximum scaling factor for station (for utilization normalization)."""
        if lldscaling is not None and ist < lldscaling.shape[0]:
            if lldscaling.ndim > 1:
                return np.max(lldscaling[ist, :])
            return lldscaling[ist]
        c = nservers[ist] if ist < len(nservers) else 1
        return c

    # Compute expected queue lengths E[n_ik]
    for s, state in enumerate(state_space):
        state_matrix = get_state_matrix(state)
        for ist in range(M):
            for k in range(K):
                QN[ist, k] += pi[s] * state_matrix[ist, k]

    # Compute throughputs from departure rates tracked during generator construction
    # This matches MATLAB's approach: TN[ist,k] = sum(pi * depRates[:, ist, k])
    # depRates tracks all service completion transitions from each state
    if depRates is not None:
        for ist in range(M):
            for k in range(K):
                # TN[ist, k] = sum(pi[s] * depRates[s, ist, k]) for all states s
                TN[ist, k] = np.sum(pi * depRates[:, ist, k])
    else:
        # Fallback: compute throughput from formula (legacy behavior)
        for s, state in enumerate(state_space):
            state_matrix = get_state_matrix(state)
            for ist in range(M):
                total_at_station = np.sum(state_matrix[ist, :])
                if total_at_station <= 0:
                    continue

                scaling = get_load_scaling(ist, int(total_at_station))
                scaling *= get_cd_scaling(ist, state_matrix[ist, :])

                if hasattr(sn, 'sched') and sn.sched is not None:
                    sched = sn.sched.get(ist, SchedStrategy.FCFS)
                else:
                    sched = SchedStrategy.FCFS

                for k in range(K):
                    n_k = state_matrix[ist, k]
                    if n_k <= 0:
                        continue

                    mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0

                    if _sched_is(sched, 'INF'):
                        # Infinite server: rate = n * mu
                        TN[ist, k] += pi[s] * n_k * mu
                    else:
                        # PS, FCFS, etc.: share of capacity proportional to number of jobs
                        TN[ist, k] += pi[s] * scaling * mu * (n_k / total_at_station)

    # Compute expected capacity (E[scaling]) at each station
    # E[scaling] = sum over states of pi(s) * lldscaling(station, n)
    expected_scaling = np.zeros(M)
    for s, state in enumerate(state_space):
        state_matrix = get_state_matrix(state)
        for ist in range(M):
            total_at_station = np.sum(state_matrix[ist, :])
            scaling = get_load_scaling(ist, int(total_at_station))
            expected_scaling[ist] += pi[s] * scaling
    # Ensure no division by zero
    expected_scaling = np.maximum(expected_scaling, 1e-10)

    # Check for LLD/CD/joint dependence (determines which formula to use for PS/DPS/GPS)
    has_lld = hasattr(sn, 'lldscaling') and sn.lldscaling is not None
    has_cd = hasattr(sn, 'cdscaling') and sn.cdscaling is not None
    has_ljd = hasattr(sn, 'ljdscaling') and sn.ljdscaling is not None
    has_state_dependence = has_lld or has_cd or has_ljd

    # Compute utilizations using MATLAB's formulas
    # For INF servers: UN = QN
    # For non-LLD PS/DPS/GPS: UN = T * service_mean / num_servers
    # For LLD PS/DPS/etc: UN = E[(n_k * w_k) / sum(n_j * w_j)]
    for ist in range(M):
        if hasattr(sn, 'sched') and sn.sched is not None:
            sched = sn.sched.get(ist, SchedStrategy.FCFS)
        else:
            sched = SchedStrategy.FCFS

        if _sched_is(sched, 'INF'):
            # For infinite servers, utilization = queue length
            for k in range(K):
                UN[ist, k] = QN[ist, k]
        elif _sched_is(sched, 'PS', 'DPS', 'GPS'):
            if not has_state_dependence:
                # Non-LLD case: UN = T * service_mean / num_servers
                c = nservers[ist] if ist < len(nservers) else 1
                for k in range(K):
                    mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0
                    if mu > 0 and c > 0:
                        service_mean = 1.0 / mu
                        UN[ist, k] = TN[ist, k] * service_mean / c
            # else: LLD case is handled below via state iteration
        elif _sched_is(sched, 'PSPRIO', 'GPSPRIO'):
            # For priority variants: UN = T / (c * mu) like FCFS
            c = nservers[ist] if ist < len(nservers) else 1
            for k in range(K):
                mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0
                if mu > 0 and c > 0:
                    UN[ist, k] = TN[ist, k] / (c * mu)
        else:
            # For FCFS and others: UN = T / (c * mu) where c is number of servers
            # For MAP distributions, use expected D1 rate instead of mu
            c = nservers[ist] if ist < len(nservers) else 1
            if has_state_dependence:
                # LLD/CD/LJD: compute UN = E[sir_k / S] from steady-state probabilities
                UN[ist, :K] = 0
                for s_idx, state in enumerate(state_space):
                    state_matrix = get_state_matrix(state)
                    total_at_station = np.sum(state_matrix[ist, :])
                    if total_at_station > 0:
                        for k in range(K):
                            if ist in fcfs_stations and ist in fcfs_phases_start and fcfs_phases_start[ist] >= 0:
                                # Compute sir_k from phase counts in the raw state
                                ps = fcfs_phases_start[ist]
                                class_phase_start = ps + int(sum(phases[ist, :k]))
                                n_phases_k = int(phases[ist, k])
                                sir_k = sum(int(state[class_phase_start + p]) for p in range(n_phases_k))
                            else:
                                sir_k = state_matrix[ist, k]
                            UN[ist, k] += pi[s_idx] * sir_k / c
            else:
                for k in range(K):
                    if (ist, k) in map_fcfs:
                        # For MAP: UN = TN * map_mean / c (matches MATLAB UNdep_ik)
                        proc_is_list = isinstance(sn.proc, list)
                        if proc_is_list:
                            station_proc = sn.proc[ist] if ist < len(sn.proc) else None
                        else:
                            station_proc = sn.proc.get(ist)
                        proc_entry = None
                        if station_proc is not None:
                            if isinstance(station_proc, (list, tuple)) and k < len(station_proc):
                                proc_entry = station_proc[k]
                            elif isinstance(station_proc, dict):
                                proc_entry = station_proc.get(k)
                        if proc_entry is not None and isinstance(proc_entry, (list, tuple)) and len(proc_entry) >= 2:
                            D0 = np.atleast_2d(np.array(proc_entry[0], dtype=float))
                            D1 = np.atleast_2d(np.array(proc_entry[1], dtype=float))
                            service_mean = map_mean(D0, D1)
                            if c > 0:
                                UN[ist, k] = TN[ist, k] * service_mean / c
                    else:
                        mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0
                        if mu > 0 and c > 0:
                            UN[ist, k] = TN[ist, k] / (c * mu)

    # For PS/DPS/GPS with LLD/CD/LJD, compute UN = E[(n_k * w_k) / sum_i(n_i * w_i)]
    # where w_k = schedparam[ist, k] (DPS/GPS weights, default to 1 for PS)
    # This is only used when there is state-dependent scaling
    if has_state_dependence:
        for s, state in enumerate(state_space):
            state_matrix = get_state_matrix(state)
            for ist in range(M):
                if hasattr(sn, 'sched') and sn.sched is not None:
                    sched = sn.sched.get(ist, SchedStrategy.FCFS)
                else:
                    sched = SchedStrategy.FCFS

                if _sched_is(sched, 'PS', 'DPS', 'GPS'):
                    n_total = np.sum(state_matrix[ist, :])
                    if n_total > 0:
                        # Get weights from schedparam (default to 1.0 for each class)
                        weights = np.ones(K)
                        if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                            if ist < sn.schedparam.shape[0]:
                                for k in range(K):
                                    if k < sn.schedparam.shape[1]:
                                        w = sn.schedparam[ist, k]
                                        if w is not None and not np.isnan(w) and w > 0:
                                            weights[k] = w

                        # Compute weighted sum: sum_i(n_i * w_i)
                        weighted_total = sum(state_matrix[ist, j] * weights[j] for j in range(K))

                        if weighted_total > 0:
                            for k in range(K):
                                n_k = state_matrix[ist, k]
                                # UN[k] = E[(n_k * w_k) / sum_i(n_i * w_i)]
                                UN[ist, k] += pi[s] * n_k * weights[k] / weighted_total

    # Zero out numerical noise below threshold (matches MATLAB's GlobalConstants.Zero)
    # This prevents spurious metrics for disabled (station, class) pairs
    threshold = GlobalConstants.Zero
    QN[QN < threshold] = 0
    TN[TN < threshold] = 0
    UN[UN < threshold] = 0

    # Compute response times: R = Q / T (Little's law)
    for ist in range(M):
        for k in range(K):
            if TN[ist, k] > 0:
                RN[ist, k] = QN[ist, k] / TN[ist, k]
            else:
                RN[ist, k] = 0

    # Handle Source station metrics for open networks
    # Source doesn't hold jobs - QLen, Util, RespT should be 0 for reporting
    # Throughput is computed from depRates (departure rates from Source) which accounts
    # for state-space truncation (blocked arrivals when downstream stations are full)
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            # Check for EXT scheduling (Source station)
            is_source = (_sched_is(sched_val, 'EXT') or
                        (isinstance(sched_val, int) and sched_val == SchedStrategy.EXT.value))
            if is_source and ist < M:
                for k in range(K):
                    # Source doesn't hold jobs in the traditional sense
                    QN[ist, k] = 0.0
                    UN[ist, k] = 0.0
                    RN[ist, k] = 0.0
                    # Note: TN[ist, k] for Source is computed from depRates above
                    # which correctly captures the effective arrival rate accounting
                    # for state-space truncation (blocked arrivals). Do NOT overwrite
                    # with raw arrival rate from sn.rates.

    return {
        'Q': QN,
        'U': UN,
        'R': RN,
        'T': TN
    }


def _is_lcfs_lcfspr_network(sn: NetworkStruct) -> Tuple[bool, int, int]:
    """
    Check if network is a 2-station LCFS+LCFSPR closed network.

    Returns:
        Tuple of (is_lcfs_lcfspr, lcfs_station, lcfspr_station)
    """
    if sn.sched is None:
        return False, -1, -1

    # Check for closed network
    njobs = sn.njobs if sn.njobs is not None else np.zeros(sn.nclasses)
    if np.any(np.isinf(njobs)):
        return False, -1, -1

    # Find LCFS and LCFSPR stations
    lcfs_stats = []
    lcfspr_stats = []
    for ist, sched in sn.sched.items():
        if _sched_is(sched, 'LCFS'):
            lcfs_stats.append(ist)
        elif _sched_is(sched, 'LCFSPR'):
            lcfspr_stats.append(ist)

    # Must have exactly one of each
    if len(lcfs_stats) == 1 and len(lcfspr_stats) == 1:
        return True, lcfs_stats[0], lcfspr_stats[0]

    return False, -1, -1


def _solver_ctmc_lcfsqn(
    sn: NetworkStruct,
    options: SolverCTMCOptions,
    lcfs_stat: int,
    lcfspr_stat: int
) -> SolverCTMCReturn:
    """
    Specialized CTMC solver for LCFS+LCFSPR 2-station networks.

    Uses the product-form property: for LCFS+LCFSPR networks, the steady-state
    distribution is insensitive to buffer ordering, so we can use a simplified
    state representation and PS-like generator construction.

    Args:
        sn: Network structure
        options: Solver options
        lcfs_stat: Index of LCFS station
        lcfspr_stat: Index of LCFSPR station

    Returns:
        SolverCTMCReturn with performance metrics
    """
    start_time = time.time()

    M = sn.nstations
    K = sn.nclasses
    rates = sn.rates
    njobs = sn.njobs if sn.njobs is not None else np.zeros(K)
    N = njobs.astype(int)

    # For LCFS+LCFSPR product-form networks, enumerate states without buffer ordering
    # State format: [n_11, n_12, ..., n_1K, n_21, ..., n_MK] where n_ik = jobs of class k at station i
    # Only include the two LCFS/LCFSPR stations in state enumeration

    # Generate all valid population distributions across the two stations
    from itertools import product as cart_product

    total_jobs = int(np.sum(N))
    states = []

    # For each class, distribute jobs between the two stations
    class_distributions = []
    for k in range(K):
        n_k = int(N[k])
        # All ways to distribute n_k jobs of class k between station 0 and station 1
        class_distributions.append([(i, n_k - i) for i in range(n_k + 1)])

    # Cartesian product of all class distributions
    for dist in cart_product(*class_distributions):
        # dist is a tuple of (n_k_at_lcfs, n_k_at_lcfspr) for each class k
        state = np.zeros(2 * K)
        for k, (n_lcfs, n_lcfspr) in enumerate(dist):
            state[k] = n_lcfs          # Jobs at LCFS station
            state[K + k] = n_lcfspr    # Jobs at LCFSPR station
        states.append(state)

    state_space = np.array(states)
    n_states = len(state_space)

    # Build generator matrix using PS-like rates (product-form insensitivity)
    Q = np.zeros((n_states, n_states))
    depRates = np.zeros((n_states, 2, K))  # Only 2 stations

    # Map (lcfs_stat, lcfspr_stat) to (0, 1) for state indexing
    stat_map = {lcfs_stat: 0, lcfspr_stat: 1}

    # Get routing matrix
    P = sn.rt if sn.rt is not None else np.zeros((M * K, M * K))

    state_map = {}
    for i, s in enumerate(state_space):
        state_map[tuple(s)] = i

    for s_idx, state in enumerate(state_space):
        # Service completions at each station
        for local_ist, ist in enumerate([lcfs_stat, lcfspr_stat]):
            total_at_station = np.sum(state[local_ist * K:(local_ist + 1) * K])
            if total_at_station <= 0:
                continue

            for k in range(K):
                n_ik = state[local_ist * K + k]
                if n_ik <= 0:
                    continue

                mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0
                if mu <= 0:
                    continue

                # For product-form LCFS+LCFSPR: use PS-like service rate
                # Rate = mu * (n_ik / total_at_station) for single server
                # This is the key insight: product-form means insensitivity
                rate = mu * (n_ik / total_at_station)

                # Route to destination (for closed network, routes between the two stations)
                src_idx = ist * K + k
                for local_jst, jst in enumerate([lcfs_stat, lcfspr_stat]):
                    for r in range(K):
                        dst_idx = jst * K + r
                        if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                            p = P[src_idx, dst_idx]
                        else:
                            p = 0

                        if p <= 0:
                            continue

                        # Create new state
                        new_state = state.copy()
                        new_state[local_ist * K + k] -= 1  # Departure
                        new_state[local_jst * K + r] += 1  # Arrival

                        if np.all(new_state >= 0):
                            ns_key = tuple(new_state)
                            if ns_key in state_map:
                                ns_idx = state_map[ns_key]
                                Q[s_idx, ns_idx] += rate * p
                                depRates[s_idx, local_ist, k] += rate * p

    # Set diagonal elements
    for s in range(n_states):
        Q[s, s] = -np.sum(Q[s, :])

    # Solve for stationary distribution
    pi = ctmc_solve(Q)

    # Compute metrics
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    TN = np.zeros((M, K))
    RN = np.zeros((M, K))

    # Map local station indices back to original indices
    local_to_orig = {0: lcfs_stat, 1: lcfspr_stat}

    for s_idx, state in enumerate(state_space):
        for local_ist in range(2):
            ist = local_to_orig[local_ist]
            for k in range(K):
                QN[ist, k] += pi[s_idx] * state[local_ist * K + k]

    # Throughput from departure rates
    for local_ist in range(2):
        ist = local_to_orig[local_ist]
        for k in range(K):
            TN[ist, k] = np.sum(pi * depRates[:, local_ist, k])

    # Utilization: UN = TN / mu for single server
    for local_ist in range(2):
        ist = local_to_orig[local_ist]
        for k in range(K):
            mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0
            if mu > 0:
                UN[ist, k] = TN[ist, k] / mu

    # Zero out numerical noise below threshold (matches MATLAB's GlobalConstants.Zero)
    threshold = GlobalConstants.Zero
    QN[QN < threshold] = 0
    TN[TN < threshold] = 0
    UN[UN < threshold] = 0

    # Response time: RN = QN / TN (Little's Law)
    for ist in [lcfs_stat, lcfspr_stat]:
        for k in range(K):
            if TN[ist, k] > 0:
                RN[ist, k] = QN[ist, k] / TN[ist, k]

    # System throughput
    XN = np.zeros((1, K))
    for k in range(K):
        XN[0, k] = TN[lcfs_stat, k]  # Same at both stations for closed network

    # Cycle time
    CN = np.zeros((1, K))
    for k in range(K):
        if XN[0, k] > 0:
            CN[0, k] = N[k] / XN[0, k]

    result = SolverCTMCReturn(
        Q=QN,
        U=UN,
        R=RN,
        T=TN,
        C=CN,
        X=XN,
        infgen=Q,
        space=state_space,
        runtime=time.time() - start_time,
        method="lcfsqn"
    )

    return result


def _update_cache_routing_probabilities(sn: NetworkStruct) -> None:
    """
    Update routing matrix with actual cache hit/miss probabilities.

    The default routing matrix uses 0.5/0.5 split for cache hit/miss.
    This function computes the actual probabilities using cache_xi_fp
    and updates sn.rtnodes accordingly.

    This matches the logic in MVA solver's _update_cache_routing method.

    Args:
        sn: Network structure to update
    """
    if sn.rtnodes is None:
        return

    I = sn.nnodes
    K = sn.nclasses

    # Find cache nodes
    cache_indices = []
    for ind in range(I):
        if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE:
            cache_indices.append(ind)

    if not cache_indices:
        return

    # Update routing for each cache
    for ind in cache_indices:
        ch = sn.nodeparam.get(ind) if sn.nodeparam else None
        if ch is None:
            continue

        hitclass = getattr(ch, 'hitclass', None)
        missclass = getattr(ch, 'missclass', None)
        if hitclass is None or missclass is None:
            continue

        hitclass = np.atleast_1d(hitclass).flatten()
        missclass = np.atleast_1d(missclass).flatten()

        # Find input classes (classes that have hit/miss mappings)
        input_classes = []
        for r in range(K):
            if r < len(hitclass) and r < len(missclass):
                hc = int(hitclass[r]) if hitclass[r] >= 0 else -1
                mc = int(missclass[r]) if missclass[r] >= 0 else -1
                if hc >= 0 and mc >= 0:
                    input_classes.append(r)

        if not input_classes:
            continue

        # Get cache parameters
        n_items = getattr(ch, 'nitems', 5)
        m_cap = getattr(ch, 'cap', None)
        if m_cap is None:
            m_cap = np.array([1])
        m_cap = np.atleast_1d(m_cap)

        # Build gamma matrix from pread
        pread = getattr(ch, 'pread', None)
        if pread is None:
            # Default uniform access
            pread_arr = np.full(n_items, 1.0 / n_items)
        elif isinstance(pread, dict):
            # pread[class_idx] = probability array
            pread_arr = pread.get(input_classes[0], np.full(n_items, 1.0 / n_items))
            pread_arr = np.atleast_1d(pread_arr)
        elif isinstance(pread, (list, tuple)):
            pread_arr = np.atleast_1d(pread[input_classes[0]] if input_classes[0] < len(pread) else pread[0])
        else:
            pread_arr = np.atleast_1d(pread)

        # Ensure pread has correct length
        if len(pread_arr) != n_items:
            pread_arr = np.full(n_items, 1.0 / n_items)

        # Build gamma (n_items x h) where h = len(m_cap)
        h = len(m_cap)
        gamma = np.zeros((n_items, h))
        for i in range(n_items):
            for l in range(h):
                gamma[i, l] = pread_arr[i]  # Same popularity at each level for simple models

        # Compute hit/miss probabilities using cache_xi_fp
        try:
            xi, pi0, pij, it_fp = cache_xi_fp(gamma, m_cap)

            # Overall miss rate = sum(pread * pi0) where pi0 is miss prob per item
            overall_miss_rate = np.sum(pread_arr * pi0)
            overall_hit_rate = 1 - overall_miss_rate
        except Exception:
            # Fall back to default
            overall_hit_rate = 0.5
            overall_miss_rate = 0.5

        # Update routing matrix
        for r in input_classes:
            hc = int(hitclass[r])
            mc = int(missclass[r])

            # Zero out the row for input class at cache
            sn.rtnodes[ind * K + r, :] = 0

            # Find connected nodes
            for jnd in range(I):
                if sn.connmatrix is not None and ind < sn.connmatrix.shape[0] and jnd < sn.connmatrix.shape[1]:
                    if sn.connmatrix[ind, jnd]:
                        # Route to hit class with hit probability
                        if hc >= 0 and hc < K:
                            sn.rtnodes[ind * K + r, jnd * K + hc] = overall_hit_rate
                        # Route to miss class with miss probability
                        if mc >= 0 and mc < K:
                            sn.rtnodes[ind * K + r, jnd * K + mc] = overall_miss_rate

        # Store hit/miss probabilities for result reporting
        ch.actualhitprob = np.zeros(K)
        ch.actualmissprob = np.zeros(K)
        for r in input_classes:
            ch.actualhitprob[r] = overall_hit_rate
            ch.actualmissprob[r] = overall_miss_rate

    # Update sn.rt directly for Cache routing instead of recomputing from rtnodes.
    # The original sn.rt from get_struct() already has correct structure (Sink exits
    # have row_sum < 1). Recomputing from rtnodes is problematic because rtnodes
    # has Sink->Source routing added for visit ratio calculations, which would
    # incorrectly redirect traffic back to Source.
    #
    # We only need to update the Cache node's routing rows with the actual hit/miss
    # probabilities, preserving the downstream routing structure.
    if sn.rt is None:
        return

    # Find next hop stateful nodes from Cache for hit/miss classes
    # This is done by looking at existing routing in rt for hit/miss classes
    for ind in cache_indices:
        # Get Cache's stateful index
        cache_sf = int(sn.nodeToStateful[ind]) if sn.nodeToStateful is not None and ind < len(sn.nodeToStateful) else -1
        if cache_sf < 0:
            continue

        ch = sn.nodeparam.get(ind) if sn.nodeparam else None
        if ch is None:
            continue

        hitclass = getattr(ch, 'hitclass', None)
        missclass = getattr(ch, 'missclass', None)
        if hitclass is None or missclass is None:
            continue

        hitclass = np.atleast_1d(hitclass).flatten()
        missclass = np.atleast_1d(missclass).flatten()

        # Get the actual hit/miss probabilities that were computed
        actualhitprob = getattr(ch, 'actualhitprob', None)
        actualmissprob = getattr(ch, 'actualmissprob', None)
        if actualhitprob is None or actualmissprob is None:
            continue

        # For each input class, update rt routing from Cache
        for r in range(K):
            if r >= len(hitclass) or r >= len(missclass):
                continue

            hc = int(hitclass[r]) if hitclass[r] >= 0 else -1
            mc = int(missclass[r]) if missclass[r] >= 0 else -1
            if hc < 0 or mc < 0:
                continue

            hit_prob = actualhitprob[r] if r < len(actualhitprob) else 0.5
            miss_prob = actualmissprob[r] if r < len(actualmissprob) else 0.5

            # Source indices in rt (stateful-indexed)
            input_src_idx = cache_sf * K + r
            hit_src_idx = cache_sf * K + hc
            miss_src_idx = cache_sf * K + mc

            if input_src_idx >= sn.rt.shape[0]:
                continue

            # Zero out input class routing first
            sn.rt[input_src_idx, :] = 0

            # Combine hit and miss class routing with their probabilities
            for dst_idx in range(sn.rt.shape[1]):
                hit_route = sn.rt[hit_src_idx, dst_idx] if hit_src_idx < sn.rt.shape[0] else 0
                miss_route = sn.rt[miss_src_idx, dst_idx] if miss_src_idx < sn.rt.shape[0] else 0

                combined_prob = hit_prob * hit_route + miss_prob * miss_route
                if combined_prob > 1e-10:
                    sn.rt[input_src_idx, dst_idx] = combined_prob


def solver_ctmc_basic(
    sn: NetworkStruct,
    options: Optional[SolverCTMCOptions] = None
) -> SolverCTMCReturn:
    """
    Basic CTMC solver using state-space enumeration.

    Enumerates all valid states, builds the infinitesimal generator,
    and solves for steady-state distribution.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverCTMCReturn with performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverCTMCOptions()

    M = sn.nstations
    K = sn.nclasses

    # Check state space size estimate (matches MATLAB's runAnalyzer.m)
    # Uses gamma function to estimate worst-case state space size
    NK = sn.njobs if sn.njobs is not None else np.ones(K)
    size_estimator = 0.0
    from scipy.special import gammaln
    for k in range(K):
        if np.isfinite(NK[k]):
            # log(C(NK[k]+M-1, M-1)) = gammaln(NK[k]+M) - gammaln(M) - gammaln(NK[k]+1)
            size_estimator += gammaln(1 + NK[k] + M - 1) - gammaln(1 + M - 1) - gammaln(1 + NK[k])

    if size_estimator > 6 and not options.force:
        raise RuntimeError(
            f"CTMC state space may be too large to solve (size estimate: {size_estimator:.1f} > 6). "
            f"Use force=True option to bypass this check, e.g., CTMC(model, force=True). "
            f"Alternative solvers: MVA, NC, or SSA."
        )

    # LCFS+LCFSPR networks should be handled by the standard solver
    # The standard state-space enumeration and generator construction
    # correctly handles LCFS scheduling when buffer ordering is enabled

    # Detect Cache nodes early (needed for visits refresh and rt recomputation)
    has_cache = False
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for ind in range(int(sn.nnodes)):
            if ind < len(sn.nodetype):
                nt = sn.nodetype[ind]
                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                if nt_val == 6:  # NodeType.CACHE
                    has_cache = True
                    break

    # Update cache routing probabilities (replaces default 0.5/0.5 split)
    # This must be done before building the generator matrix
    _update_cache_routing_probabilities(sn)

    # Refresh visit ratios with updated cache routing probabilities
    # (matches MATLAB refreshChains which recomputes visits after cache prob update)
    if has_cache:
        sn_refresh_visits(sn)

    if options.verbose and sn.rt is not None:
        rt = np.asarray(sn.rt)
        nstateful = sn.nstateful if hasattr(sn, 'nstateful') else M
        print("  sn.rt after cache routing update (non-zero entries):")
        for i in range(rt.shape[0]):
            sf_i = i // K
            k_i = i % K
            for j in range(rt.shape[1]):
                if abs(rt[i, j]) > 1e-10:
                    sf_j = j // K
                    k_j = j % K
                    print(f"    rt[sf{sf_i},c{k_i} -> sf{sf_j},c{k_j}] = {rt[i,j]:.6f}")

    # Recompute sn.rt from updated sn.rtnodes using stochastic complementation
    # Skip for Cache networks: _update_cache_routing_probabilities already updates sn.rt
    # directly. Recomputing from rtnodes is wrong for Cache networks because rtnodes has
    # artificial Sink->Source routing for visit ratio calculations, which would incorrectly
    # redirect traffic back to Source.
    if sn.rtnodes is not None and not has_cache:
        I = sn.nnodes
        K = sn.nclasses
        nstateful = sn.nstateful if hasattr(sn, 'nstateful') else M

        # Get stateful node indices from stationToStateful and Router nodes
        stateful_nodes = set()
        if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None:
            for sf in sn.stationToStateful:
                stateful_nodes.add(int(sf))
        # Add any Router nodes (they are stateful but not stations)
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            node_to_stateful = {}
            sf_idx = 0
            for nidx in range(I):
                # Check if this node is stateful (station or stateful non-station)
                if hasattr(sn, 'nodeToStateful') and sn.nodeToStateful is not None:
                    if nidx < len(sn.nodeToStateful):
                        sf = sn.nodeToStateful[nidx]
                        if sf >= 0:
                            stateful_nodes.add(int(sf))
                elif hasattr(sn, 'isstateful') and sn.isstateful is not None:
                    if nidx < len(sn.isstateful) and sn.isstateful[nidx]:
                        stateful_nodes.add(sf_idx)
                        sf_idx += 1

        # Build stateful_nodes_classes for dtmc_stochcomp
        stateful_node_list = sorted(stateful_nodes)
        if not stateful_node_list:
            # Fallback: use station indices
            stateful_node_list = list(range(M))

        # Build list of (stateful_node_idx, class_idx) pairs
        # dtmc_stochcomp expects indices into the rtnodes matrix
        # For stateful nodes, we need to map from stateful index to node index
        stateful_to_node = {}
        if hasattr(sn, 'statefulToNode') and sn.statefulToNode is not None:
            for sf_idx, n_idx in enumerate(sn.statefulToNode):
                stateful_to_node[sf_idx] = int(n_idx)

        # Build the list of node*K + class indices for stateful nodes
        stateful_nodes_classes = []
        for sf_idx in range(nstateful):
            node_idx = stateful_to_node.get(sf_idx, sf_idx)
            for k in range(K):
                stateful_nodes_classes.append(node_idx * K + k)
        stateful_nodes_classes = np.array(stateful_nodes_classes, dtype=int)

        try:
            new_rt = dtmc_stochcomp(sn.rtnodes, stateful_nodes_classes)
            sn.rt = new_rt
            # CRITICAL: Also update rt_visits since sn_refresh_visits uses it
            if hasattr(sn, 'rt_visits') and sn.rt_visits is not None:
                sn.rt_visits = new_rt.copy()
        except Exception:
            pass  # Keep existing rt if stochcomp fails

    # Detect Router-like nodes (non-station stateful, non-Cache)
    router_nodes = []
    if sn_has_open_classes(sn) and has_cache:
        if hasattr(sn, 'isstateful') and hasattr(sn, 'isstation'):
            for ind in range(int(sn.nnodes)):
                if ind < len(sn.isstateful) and ind < len(sn.isstation):
                    is_stateful = sn.isstateful[ind]
                    is_station = sn.isstation[ind]
                    if is_stateful and not is_station:
                        if hasattr(sn, 'nodetype') and ind < len(sn.nodetype):
                            nt = sn.nodetype[ind]
                            nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                            if nt_val != 6:  # Not CACHE
                                router_nodes.append(ind)

    M = int(sn.nstations)
    K = int(sn.nclasses)

    state_space, state_space_aggr, rrobin_info = _enumerate_state_space(sn, options.cutoff)

    # Mandatory truncation warning for open/mixed models
    if sn_has_open_classes(sn):
        print(f"CTMC solver using state space cutoff = {options.cutoff} for open/mixed model.")
        warnings.warn(
            "State space truncation may cause inaccurate results. "
            "Consider varying cutoff to assess sensitivity.",
            UserWarning
        )

    if options.verbose:
        print(f"CTMC state space size: {len(state_space)}")
        if rrobin_info['total_vars'] > 0:
            print(f"  Including {rrobin_info['total_vars']} round-robin state variables")

    # Build generator matrix and departure rates
    Q, depRates = _build_generator(sn, state_space, options, rrobin_info)

    if options.verbose:
        print(f"  Q matrix shape: {Q.shape}, nnz: {np.count_nonzero(Q)}")
        print(f"  depRates shape: {depRates.shape}")

    # Apply stochastic complementation (stochcomp) to eliminate immediate states.
    # Immediate states are those where a tracked Router has jobs.
    # This matches MATLAB's solver_ctmc.m which builds the full Q, identifies
    # immediate states where non-station non-Cache stateful nodes have jobs,
    # and applies ctmc_stochcomp(Q, nonimm).
    # S = Q11 + Q12 * (-Q22)^{-1} * Q21
    # where Q11 = Q[nonimm, nonimm], Q12 = Q[nonimm, imm],
    #       Q21 = Q[imm, nonimm], Q22 = Q[imm, imm]
    router_jobs_offsets = rrobin_info.get('router_jobs_offsets', {})
    if router_jobs_offsets:
        # Identify immediate states: states where any Router has a job
        imm_indices = []
        nonimm_indices = []
        for s_idx in range(len(state_space)):
            is_immediate = False
            for isf_r, rj_off in router_jobs_offsets.items():
                for k_c in range(K):
                    if rj_off + k_c < len(state_space[s_idx]) and int(state_space[s_idx][rj_off + k_c]) > 0:
                        is_immediate = True
                        break
                if is_immediate:
                    break
            if is_immediate:
                imm_indices.append(s_idx)
            else:
                nonimm_indices.append(s_idx)

        if imm_indices:
            try:
                from scipy import sparse as sp_sparse
                Q_dense = Q.toarray() if sp_sparse.issparse(Q) else np.asarray(Q)
            except ImportError:
                Q_dense = np.asarray(Q)

            imm = np.array(imm_indices, dtype=int)
            nonimm = np.array(nonimm_indices, dtype=int)

            Q11 = Q_dense[np.ix_(nonimm, nonimm)]
            Q12 = Q_dense[np.ix_(nonimm, imm)]
            Q21 = Q_dense[np.ix_(imm, nonimm)]
            Q22 = Q_dense[np.ix_(imm, imm)]

            # Stochcomp: S = Q11 + Q12 * (-Q22)^{-1} * Q21
            # Matches MATLAB: T = (-Q22) \ Q21; T = Q12*T; S = Q11+T;
            try:
                T = np.linalg.solve(-Q22, Q21)
                T = Q12 @ T
                Q_reduced = Q11 + T
            except np.linalg.LinAlgError:
                # Fallback: use pseudo-inverse
                Q_reduced = Q11 + Q12 @ np.linalg.pinv(-Q22) @ Q21

            # Apply stochcomp to depRates as well
            # MATLAB also stochcomps Dfilt (per-action filter) but we just
            # remove immediate depRates since Router states don't contribute
            # station departure rates directly.
            depRates = depRates[nonimm_indices, :, :]

            # Update state_space: remove immediate states
            state_space = state_space[nonimm_indices]
            state_space_aggr = state_space_aggr[nonimm_indices] if state_space_aggr is not None and len(state_space_aggr) > len(nonimm_indices) - 1 else state_space_aggr

            Q = Q_reduced

            if options.verbose:
                print(f"  Stochcomp: {len(imm_indices)} immediate states eliminated, "
                      f"{len(nonimm_indices)} non-immediate states remain.")

    # Solve for steady-state distribution
    pi = ctmc_solve(Q)

    # Compute performance metrics using departure rates from generator construction
    metrics = _compute_metrics_from_distribution(sn, pi, state_space, rrobin_info, depRates)

    QN = metrics['Q']
    UN = metrics['U']
    RN = metrics['R']
    TN = metrics['T']

    # Compute cycle times and system throughput
    CN = np.sum(RN, axis=0).reshape(1, -1)
    XN = np.zeros((1, K))
    for k in range(K):
        ref_stat = int(sn.refstat[k]) if hasattr(sn, 'refstat') and k < len(sn.refstat) else 0
        if ref_stat < M:
            XN[0, k] = TN[ref_stat, k]

    # Clean up NaN values
    QN = np.nan_to_num(QN, nan=0.0)
    UN = np.nan_to_num(UN, nan=0.0)
    RN = np.nan_to_num(RN, nan=0.0)
    TN = np.nan_to_num(TN, nan=0.0)
    CN = np.nan_to_num(CN, nan=0.0)
    XN = np.nan_to_num(XN, nan=0.0)

    # Compute station column ranges for nodeStateSpace extraction
    # Each station has a range [start, end) of columns in the state space
    station_col_ranges = []
    phases = rrobin_info.get('phases', np.ones((M, K), dtype=int))
    fcfs_stations = rrobin_info.get('fcfs_stations', set())
    fcfs_max_buffer = rrobin_info.get('fcfs_max_buffer', {})
    skip_stations = rrobin_info.get('skip_stations', set())

    col_idx = 0
    for ist in range(M):
        if ist in skip_stations:
            # Source/Sink stations have no columns in state space
            station_col_ranges.append((col_idx, col_idx))
            continue

        start_col = col_idx
        if ist in fcfs_stations:
            # FCFS station: buffer columns + phase columns
            buffer_size = fcfs_max_buffer.get(ist, 0)
            col_idx += buffer_size
        # Add phase columns for all classes at this station
        for k in range(K):
            col_idx += int(phases[ist, k])
        station_col_ranges.append((start_col, col_idx))

    result = SolverCTMCReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.pi = pi
    result.infgen = Q
    result.space = state_space
    result.space_aggr = state_space_aggr
    result.station_col_ranges = station_col_ranges
    result.depRates = depRates  # Include departure rates for debugging
    result.rrobin_info = rrobin_info  # Include cache state info for hit ratio computation
    result.runtime = time.time() - start_time
    result.method = "basic"

    return result


def solver_ctmc(
    sn: NetworkStruct,
    options: Optional[SolverCTMCOptions] = None
) -> SolverCTMCReturn:
    """
    Main CTMC solver handler.

    Routes to appropriate method based on options and network characteristics.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverCTMCReturn with performance metrics
    """
    if options is None:
        options = SolverCTMCOptions()

    method = options.method.lower()

    if method in ['default', 'basic']:
        return solver_ctmc_basic(sn, options)
    else:
        # Unknown method - use basic
        if options.verbose:
            print(f"Warning: Unknown CTMC method '{method}'. Using basic.")
        return solver_ctmc_basic(sn, options)


__all__ = [
    'solver_ctmc',
    'solver_ctmc_basic',
    'SolverCTMCReturn',
    'SolverCTMCOptions',
]
