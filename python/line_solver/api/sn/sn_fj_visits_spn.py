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

    Evaluates the rule the reference's SPN solve produces instead of enumerating
    a state space exponential in B. The pre-fork transition consumes all B tokens
    at once and the Join returns them, so a station INSIDE a fork-join region
    fires once per B firings of the cycle: normalized on the reference station, a
    station outside the region carries 1, a station inside it carries 1/B, and a
    Fork or a Join, which holds no Place, carries 0. Verified against MATLAB and
    the JAR (both of which do solve the net) on two-branch, three-branch, nested,
    chained-branch, pre/post-fork-station and two-class models.

    Args:
        sn: NetworkStruct describing the queueing network.
        P_r: Single-class sub-routing matrix (nnodes x nnodes).
        visited: Boolean array indicating visited nodes.
        r: Class index (0-based).
        refnode: Reference node index (0-based).

    Returns:
        visits_r: Array of shape (nnodes,) carrying 1 outside the fork-join
                  region, 1/B inside it and 0 at a Fork or a Join.
    """
    I = sn.nnodes
    visits_r = np.zeros(I)
    B, in_region = _spn_structure(sn, P_r, visited, r)

    # see _kb/03-api-layer.md for rationale
    for nd in range(I):
        if visited[nd] and sn.isstation[nd] and \
                sn.nodetype[nd] not in (NodeType.SOURCE, NodeType.SINK,
                                        NodeType.FORK, NodeType.JOIN):
            visits_r[nd] = 1.0 / B if in_region[nd] else 1.0

    return visits_r


def _spn_structure(sn: NetworkStruct, P_r: np.ndarray,
                   visited: np.ndarray, r: int = 0):
    """Size the auxiliary net: its population B and its fork-join region.

    B is the largest EXPECTED leaf count over the OUTERMOST forks, a fork fed by
    a Join being a serial stage rather than an outer one. The region is
    everything an outermost fork opens, followed to the end of each branch and
    not only to its first station, because a chained branch station runs at the
    branch rate too.

    On a plain fork every link carries one certain task, so B is the leaf count
    it has always been and is integral. Under a variable forking level it is
    sum over links of P(branch fires) * E[tasks on it], which is generally
    FRACTIONAL -- and a fractional token population is not a net anyone can
    enumerate, which is why the closed form is the answer everywhere here.

    Returns:
        (B, in_region) with B a float >= 0 and in_region a boolean array.
    """
    I = sn.nnodes
    in_region = np.zeros(I, dtype=bool)
    B = 0
    for fnd in range(I):
        if not visited[fnd] or sn.nodetype[fnd] != NodeType.FORK:
            continue
        outermost = True
        for src in range(I):
            if P_r[src, fnd] > 0 and visited[src] and sn.nodetype[src] == NodeType.JOIN:
                outermost = False
                break
        if not outermost:
            continue
        frontier = [fnd]
        while frontier:
            nd = frontier.pop()
            for j in range(I):
                if P_r[nd, j] <= 0 or not visited[j]:
                    continue
                if sn.nodetype[j] == NodeType.JOIN or in_region[j]:
                    continue
                in_region[j] = True
                frontier.append(j)
        leaves, weights = _resolve_fork_dests(sn, P_r, visited, fnd, r)
        if not leaves:
            raise ValueError(
                'sn_fj_visits_spn: a Fork reaches no station on any branch, so '
                'the auxiliary net has nothing to synchronize')
        expected = float(sum(weights))
        if expected <= 0.0:
            raise ValueError(
                'sn_fj_visits_spn: a Fork emits no task in expectation, so its '
                'Join can never fire; at least one branch must be certain to '
                'emit at least one task')
        B = max(B, expected)
    return (B if B > 0 else 1.0), in_region


def _resolve_fork_dests(sn: NetworkStruct, P_r: np.ndarray,
                        visited: np.ndarray, fork_nd: int, r: int = 0,
                        w: float = 1.0):
    """Recursively resolve Fork destinations to station nodes.

    The second return is the EXPECTED number of tasks each leaf receives per
    firing of the outermost fork: P(branch fires) times E[tasks on that link],
    multiplied down through any nesting. A plain fork gives every link exactly
    1, so the weighted leaf count collapses to the leaf count it always was.

    Args:
        sn: NetworkStruct describing the queueing network.
        P_r: Single-class sub-routing matrix (nnodes x nnodes).
        visited: Boolean array indicating visited nodes.
        fork_nd: Fork node index (0-based).
        r: Class index (0-based).
        w: Expected tasks accumulated on the way to this fork.

    Returns:
        (st_dests, weights) for the leaf destinations of this fork.
    """
    st_dests, weights = [], []
    param = None
    if sn.nodeparam is not None and isinstance(sn.nodeparam.get(fork_nd), dict):
        param = sn.nodeparam[fork_nd]
    has_fan = param is not None and param.get('fanOutLink', None) is not None
    branch_dests = [j for j in range(sn.nnodes) if P_r[fork_nd, j] > 0 and visited[j]]
    for bd in branch_dests:
        w_bd = w
        if has_fan:
            w_bd = w * float(param['fanOutProb'][bd, r]) * float(param['fanOutLink'][bd, r])
        if sn.nodetype[bd] == NodeType.FORK:
            d2, w2 = _resolve_fork_dests(sn, P_r, visited, bd, r, w_bd)
            st_dests.extend(d2)
            weights.extend(w2)
        elif sn.isstation[bd]:
            st_dests.append(bd)
            weights.append(w_bd)
    return st_dests, weights
