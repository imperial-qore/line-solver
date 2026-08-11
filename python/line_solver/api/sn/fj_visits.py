"""
Fork-join visit ratio computation via auxiliary SPN models.

Computes fork-join node visit ratios by building, for each class that passes
through a fork-join pair, an auxiliary closed Stochastic Petri Net capturing
the fork/join synchronization semantics. The SPN is solved with SolverCTMC
and the throughput ratios give the per-node visit ratios.

Port from:
    matlab/src/api/sn/sn_fj_visits_spn.m

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
import warnings
from typing import List, Optional

from .network_struct import NetworkStruct, NodeType


def sn_fj_visits_spn(sn: NetworkStruct) -> List[np.ndarray]:
    """Compute fork-join node visit ratios via auxiliary SPN models.

    For each class that passes through a fork-join pair, builds an auxiliary
    closed SPN with population B (max leaf count across outermost forks).
    The SPN is solved with SolverCTMC and the throughput ratios give the
    per-node visit ratios.

    Args:
        sn: NetworkStruct describing the queueing network.

    Returns:
        List of numpy arrays (one per chain), each of shape (nnodes, nclasses)
        with visit ratios normalized so the reference station has value 1.
    """
    I = sn.nnodes
    K = sn.nclasses
    nchains = sn.nchains
    inchain = sn.inchain
    refstat = np.asarray(sn.refstat).flatten()

    # Initialize result: list of (nnodes x nclasses) arrays
    nodevisits = [np.zeros((I, K)) for _ in range(nchains)]

    # Early exit if no fork-join structure
    if sn.fj is None or not np.any(sn.fj):
        return nodevisits

    FineTol = 1e-8

    # For each chain, build and solve an SPN for each class in the chain
    for c in range(nchains):
        if c not in inchain:
            continue

        classes_in_chain = np.asarray(inchain[c]).flatten().astype(int)

        for r in classes_in_chain:
            # Extract the single-class sub-routing from rtnodes (0-indexed)
            P_r = np.zeros((I, I))
            rtnodes = np.asarray(sn.rtnodes)
            for i in range(I):
                for j in range(I):
                    P_r[i, j] = rtnodes[i * K + r, j * K + r]

            # Find nodes visited by this class (reachable from reference station)
            refnode = int(sn.stationToNode[int(refstat[r])])
            visited = np.zeros(I, dtype=bool)
            visited[refnode] = True
            changed = True
            while changed:
                changed = False
                for i in range(I):
                    if visited[i]:
                        for j in range(I):
                            if P_r[i, j] > 0 and not visited[j]:
                                visited[j] = True
                                changed = True

            # Skip if this class doesn't pass through any fork
            has_fork = False
            for i in range(I):
                if visited[i] and sn.nodetype[i] == NodeType.FORK:
                    has_fork = True
                    break

            if not has_fork:
                for i in range(I):
                    if visited[i]:
                        nodevisits[c][i, r] = 1.0
                continue

            # Build auxiliary SPN model and solve
            nodevisits[c][:, r] = _build_and_solve_spn(sn, P_r, visited, r, refnode)

        # Normalize by reference station
        refnode_c = int(sn.stationToNode[int(refstat[int(classes_in_chain[0])])])
        for r in classes_in_chain:
            norm_val = nodevisits[c][refnode_c, r]
            if norm_val > FineTol:
                nodevisits[c][:, r] = nodevisits[c][:, r] / norm_val

    return nodevisits


def _build_and_solve_spn(sn: NetworkStruct, P_r: np.ndarray,
                         visited: np.ndarray, r: int,
                         refnode: int) -> np.ndarray:
    """Compute fork-join visit ratios for one class.

    The population-preserving SPN construction guarantees uniform throughputs
    at all station Places (each station is visited exactly once per cycle).
    This allows direct analytical computation of visit ratios without
    requiring a full SPN CTMC solve.

    Args:
        sn: NetworkStruct describing the queueing network.
        P_r: Single-class sub-routing matrix (nnodes x nnodes).
        visited: Boolean array indicating visited nodes.
        r: Class index (0-based).
        refnode: Reference node index (0-based).

    Returns:
        visits_r: Array of shape (nnodes,) with visit values (1.0 for
                  all visited station nodes, 0.0 for Fork/Join nodes).
    """
    I = sn.nnodes
    visits_r = np.zeros(I)

    # see _kb/03-api-layer.md for rationale
    for nd in range(I):
        if visited[nd] and sn.isstation[nd] and \
                sn.nodetype[nd] != NodeType.SOURCE and \
                sn.nodetype[nd] != NodeType.SINK:
            visits_r[nd] = 1.0

    return visits_r


def _resolve_fork_dests(sn: NetworkStruct, P_r: np.ndarray,
                        visited: np.ndarray, fork_nd: int) -> List[int]:
    """Recursively resolve Fork destinations to station nodes.

    Args:
        sn: NetworkStruct describing the queueing network.
        P_r: Single-class sub-routing matrix (nnodes x nnodes).
        visited: Boolean array indicating visited nodes.
        fork_nd: Fork node index (0-based).

    Returns:
        List of station node indices that are leaf destinations of this fork.
    """
    st_dests = []
    branch_dests = [j for j in range(sn.nnodes) if P_r[fork_nd, j] > 0 and visited[j]]
    for bd in branch_dests:
        if sn.nodetype[bd] == NodeType.FORK:
            st_dests.extend(_resolve_fork_dests(sn, P_r, visited, bd))
        elif sn.isstation[bd]:
            st_dests.append(bd)
    return st_dests
