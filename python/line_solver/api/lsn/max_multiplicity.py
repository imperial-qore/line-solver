"""
Layered Stochastic Network Maximum Multiplicity Computation.

Computes maximum multiplicity constraints for load sharing network (LSN)
analysis. Essential for determining feasible population bounds and
capacity constraints in layered queueing network models.

References:
    Casale, G., et al. "LINE: A unified library for queueing network modeling."
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, List
from enum import IntEnum


class LayeredNetworkElement(IntEnum):
    """Types of elements in a layered network."""
    TASK = 1
    ENTRY = 2
    ACTIVITY = 3
    PROCESSOR = 4
    HOST = 5


@dataclass
class LayeredNetworkStruct:
    """
    Structure representing a layered network for LSN analysis.

    Attributes:
        dag: Directed acyclic graph adjacency matrix (n x n).
             dag[i,j] > 0 indicates an edge from node i to node j.
        mult: Multiplicity (max concurrent instances) for each node (n,).
        type: Node type for each node (n,), using LayeredNetworkElement values.
        isref: Reference task flags (n,). Non-zero indicates a reference task.
    """
    dag: np.ndarray
    mult: np.ndarray
    type: np.ndarray
    isref: np.ndarray


def kahn_topological_sort(adjacency: np.ndarray) -> List[int]:
    """
    Perform Kahn's algorithm for topological sorting.

    Args:
        adjacency: Adjacency matrix where adjacency[i,j] > 0 means edge i -> j.

    Returns:
        List of node indices in topological order.

    Raises:
        ValueError: If the graph contains a cycle.
    """
    n = adjacency.shape[0]

    # Compute in-degree for each node
    in_degree = np.zeros(n, dtype=int)
    for j in range(n):
        for i in range(n):
            if adjacency[i, j] > 0:
                in_degree[j] += 1

    # Initialize queue with nodes having in-degree 0
    queue = [i for i in range(n) if in_degree[i] == 0]
    result = []

    while queue:
        # Remove a node with in-degree 0
        node = queue.pop(0)
        result.append(node)

        # Decrease in-degree of neighbors
        for j in range(n):
            if adjacency[node, j] > 0:
                in_degree[j] -= 1
                if in_degree[j] == 0:
                    queue.append(j)

    if len(result) != n:
        raise ValueError("Graph contains a cycle - cannot perform topological sort")

    return result


def lsn_max_multiplicity(lsn: LayeredNetworkStruct) -> np.ndarray:
    """
    Compute the maximum multiplicity for each task in a layered network.

    This function uses flow analysis based on Kahn's topological sorting algorithm
    to determine the maximum sustainable throughput for each task, considering
    both the incoming flow and the multiplicity constraints.

    Args:
        lsn: The layered network structure containing task dependencies and constraints.

    Returns:
        Matrix of maximum multiplicities for each task in the network (n x 1).

    Algorithm:
        1. Build binary adjacency graph from DAG
        2. Apply Kahn's topological sort to determine processing order
        3. Initialize inflow from reference tasks
        4. For each node in topological order:
           - outflow = min(inflow, multiplicity constraint)
           - Propagate outflow to downstream nodes
        5. Handle unreachable tasks (infinite multiplicity)

    Example:
        >>> lsn = LayeredNetworkStruct(
        ...     dag=np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]]),
        ...     mult=np.array([2, 3, 5]),
        ...     type=np.array([1, 1, 1]),
        ...     isref=np.array([1, 0, 0])
        ... )
        >>> max_mult = lsn_max_multiplicity(lsn)
    """
    dag = np.asarray(lsn.dag, dtype=np.float64)
    n = dag.shape[0]

    # Build binary adjacency graph
    ag = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            if dag[i, j] > 0:
                ag[i, j] = 1.0

    mult = np.asarray(lsn.mult, dtype=np.float64).flatten()
    node_type = np.asarray(lsn.type, dtype=np.float64).flatten()
    isref = np.asarray(lsn.isref, dtype=np.float64).flatten()

    # Get topological order
    order = kahn_topological_sort(ag)

    # Initialize inflow from reference tasks
    inflow = np.zeros(n)
    for ist in range(n):
        if node_type[ist] == LayeredNetworkElement.TASK and isref[ist] != 0:
            inflow[ist] = mult[ist] if ist < len(mult) else np.inf

    # Initialize outflow
    outflow = np.zeros(n)

    # Extend mult array if needed
    if len(mult) < n:
        mult_extended = np.full(n, np.inf)
        mult_extended[:len(mult)] = mult
        mult = mult_extended

    # Propagate through topological order
    for k in range(n):
        ist = order[k]
        inflow_val = inflow[ist]
        mult_val = mult[ist]
        outflow[ist] = min(inflow_val, mult_val)

        # Propagate to downstream nodes
        for jst in range(n):
            if jst != ist and ag[ist, jst] != 0:
                inflow[jst] = inflow[jst] + outflow[ist]

    # Handle non-reference tasks with infinite multiplicity
    for ist in range(n):
        if (node_type[ist] == LayeredNetworkElement.TASK and
            mult[ist] == np.inf and
            isref[ist] == 0):
            outflow[ist] = np.inf

    # Return as column vector
    return outflow.reshape(-1, 1)


__all__ = [
    'LayeredNetworkElement',
    'LayeredNetworkStruct',
    'kahn_topological_sort',
    'lsn_max_multiplicity',
]
