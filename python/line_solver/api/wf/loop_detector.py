"""
Workflow Loop Pattern Detection.

Implements algorithms for detecting loop and iterative patterns in workflow
traces. Identifies repetitive activity sequences, cyclic behaviors, and
iteration structures in business process analysis and workflow mining.

Based on AUTO_Loop_Detector.m from the MDN toolbox.

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/Wf_loop_detector.kt
"""

import numpy as np
from typing import List, Tuple, Dict, Any, Set, Optional
from collections import deque


def detect_loops(
    link_matrix: np.ndarray,
    service_nodes: List[int],
    router_nodes: List[int],
    join_nodes: Optional[List[int]] = None
) -> List[int]:
    """
    Detect loop patterns in a workflow network.

    Args:
        link_matrix: Matrix containing workflow link information:
                     - Column 0: start node IDs
                     - Column 1: end node IDs
                     - Column 2: transition probabilities
        service_nodes: List of service node IDs
        router_nodes: List of router node IDs
        join_nodes: List of join node IDs (optional, for complex loops)

    Returns:
        List of loop patterns, where each pattern contains the loop service node
    """
    if join_nodes is None:
        join_nodes = []

    loop_nodes = []
    router_set = set(router_nodes)

    # Build adjacency map
    adjacency = _build_adjacency_map(link_matrix)

    # Find simple loops: service -> router -> service
    for service_node in service_nodes:
        if _is_in_simple_loop(service_node, adjacency, router_set):
            loop_nodes.append(service_node)

    # Find complex loops involving join nodes
    if join_nodes:
        complex_loops = _find_complex_loops(link_matrix, service_nodes, router_nodes, join_nodes)
        loop_nodes.extend(complex_loops)

    return list(set(loop_nodes))


def _build_adjacency_map(link_matrix: np.ndarray) -> Dict[int, List[Tuple[int, float]]]:
    """Build adjacency map from link matrix."""
    adjacency: Dict[int, List[Tuple[int, float]]] = {}

    if link_matrix.size == 0:
        return adjacency

    for i in range(link_matrix.shape[0]):
        start = int(link_matrix[i, 0])
        end = int(link_matrix[i, 1])
        prob = float(link_matrix[i, 2])

        if start not in adjacency:
            adjacency[start] = []
        adjacency[start].append((end, prob))

    return adjacency


def _is_in_simple_loop(
    service_node: int,
    adjacency: Dict[int, List[Tuple[int, float]]],
    router_set: Set[int]
) -> bool:
    """Check if a service node is in a simple loop pattern."""
    neighbors = adjacency.get(service_node)
    if neighbors is None:
        return False

    # Check if service connects to a router
    for router_node, _ in neighbors:
        if router_node in router_set:
            # Check if router connects back to service
            router_neighbors = adjacency.get(router_node, [])
            for back_node, _ in router_neighbors:
                if back_node == service_node:
                    return True

    return False


def _find_complex_loops(
    link_matrix: np.ndarray,
    service_nodes: List[int],
    router_nodes: List[int],
    join_nodes: List[int]
) -> List[int]:
    """Find complex loops involving multiple nodes and join points."""
    complex_loops = []
    service_set = set(service_nodes)
    router_set = set(router_nodes)
    join_set = set(join_nodes)

    # Build graph for cycle detection
    graph = _build_directed_graph(link_matrix)

    # Find strongly connected components
    sccs = _find_strongly_connected_components(graph)

    for scc in sccs:
        if len(scc) > 1:
            # Check if SCC contains service nodes and forms a valid loop
            service_nodes_in_scc = [n for n in scc if n in service_set]
            has_router_or_join = any(n in router_set or n in join_set for n in scc)

            if service_nodes_in_scc and has_router_or_join:
                complex_loops.extend(service_nodes_in_scc)

    return complex_loops


def _build_directed_graph(link_matrix: np.ndarray) -> Dict[int, Set[int]]:
    """Build directed graph representation."""
    graph: Dict[int, Set[int]] = {}

    if link_matrix.size == 0:
        return graph

    for i in range(link_matrix.shape[0]):
        start = int(link_matrix[i, 0])
        end = int(link_matrix[i, 1])

        if start not in graph:
            graph[start] = set()
        graph[start].add(end)

    return graph


def _find_strongly_connected_components(graph: Dict[int, Set[int]]) -> List[List[int]]:
    """Find strongly connected components using Tarjan's algorithm."""
    sccs: List[List[int]] = []
    visited: Set[int] = set()
    stack: List[int] = []
    indices: Dict[int, int] = {}
    low_links: Dict[int, int] = {}
    on_stack: Set[int] = set()
    index = [0]  # Use list to allow modification in nested function

    def strong_connect(node: int) -> None:
        indices[node] = index[0]
        low_links[node] = index[0]
        index[0] += 1
        stack.append(node)
        on_stack.add(node)

        neighbors = graph.get(node, set())
        for neighbor in neighbors:
            if neighbor not in indices:
                strong_connect(neighbor)
                low_links[node] = min(low_links[node], low_links[neighbor])
            elif neighbor in on_stack:
                low_links[node] = min(low_links[node], indices[neighbor])

        if low_links[node] == indices[node]:
            scc = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                scc.append(w)
                if w == node:
                    break
            sccs.append(scc)

    # Run Tarjan's algorithm on all unvisited nodes
    for node in graph.keys():
        if node not in indices:
            strong_connect(node)
            visited.add(node)

    return sccs


def get_loop_probability(
    service_node: int,
    link_matrix: np.ndarray,
    router_nodes: List[int]
) -> float:
    """
    Get loop probability for a service node.

    Args:
        service_node: Service node ID
        link_matrix: Workflow link matrix
        router_nodes: List of router node IDs

    Returns:
        Loop probability (0.0 if no loop found)
    """
    router_set = set(router_nodes)

    if link_matrix.size == 0:
        return 0.0

    # Find the router that creates the loop
    for i in range(link_matrix.shape[0]):
        start = int(link_matrix[i, 0])
        end = int(link_matrix[i, 1])

        if start == service_node and end in router_set:
            # Found service -> router connection
            # Now find router -> service connection (loop back)
            for j in range(link_matrix.shape[0]):
                loop_start = int(link_matrix[j, 0])
                loop_end = int(link_matrix[j, 1])
                loop_prob = float(link_matrix[j, 2])

                if loop_start == end and loop_end == service_node:
                    return loop_prob

    return 0.0


def validate_loop_pattern(
    loop_node: int,
    link_matrix: np.ndarray,
    router_nodes: List[int]
) -> bool:
    """
    Validate loop pattern structure.

    Args:
        loop_node: Loop node ID
        link_matrix: Workflow link matrix
        router_nodes: List of router node IDs

    Returns:
        True if valid loop pattern exists
    """
    router_set = set(router_nodes)

    # Check if there's a path: loopNode -> router -> loopNode
    adjacency = _build_adjacency_map(link_matrix)
    neighbors = adjacency.get(loop_node)
    if neighbors is None:
        return False

    for router_node, _ in neighbors:
        if router_node in router_set:
            router_neighbors = adjacency.get(router_node, [])
            for back_node, _ in router_neighbors:
                if back_node == loop_node:
                    return True

    return False


def get_expected_loop_iterations(loop_probability: float) -> float:
    """
    Calculate expected number of loop iterations.

    Args:
        loop_probability: Probability of looping

    Returns:
        Expected number of iterations (Inf if probability >= 1.0)
    """
    if loop_probability >= 1.0:
        return np.inf
    return 1.0 / (1.0 - loop_probability)


def get_loop_stats(
    loop_nodes: List[int],
    link_matrix: np.ndarray,
    router_nodes: List[int]
) -> Dict[str, Any]:
    """
    Get loop pattern statistics.

    Args:
        loop_nodes: List of loop node IDs
        link_matrix: Workflow link matrix
        router_nodes: List of router node IDs

    Returns:
        Dictionary with statistics
    """
    stats: Dict[str, Any] = {}

    stats['numLoops'] = len(loop_nodes)

    probabilities = [get_loop_probability(n, link_matrix, router_nodes) for n in loop_nodes]
    stats['avgLoopProbability'] = np.mean(probabilities) if probabilities else 0.0
    stats['maxLoopProbability'] = max(probabilities) if probabilities else 0.0
    stats['minLoopProbability'] = min(probabilities) if probabilities else 0.0

    iterations = [
        get_expected_loop_iterations(p)
        for p in probabilities
        if np.isfinite(get_expected_loop_iterations(p))
    ]
    stats['avgExpectedIterations'] = np.mean(iterations) if iterations else 0.0
    stats['maxExpectedIterations'] = max(iterations) if iterations else 0.0

    return stats
