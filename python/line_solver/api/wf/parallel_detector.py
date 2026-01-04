"""
Workflow Parallel Pattern Detection.

Implements algorithms for detecting parallel execution patterns in workflow
traces. Identifies concurrent activities, synchronization points, and
parallel flow structures in business process analysis.

Based on AUTO_Parallel_Detector.m from the MDN toolbox.

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/Wf_parallel_detector.kt
"""

import numpy as np
from typing import List, Tuple, Dict, Any, Set, Optional
from collections import deque


def detect_parallel(
    link_matrix: np.ndarray,
    service_nodes: List[int],
    fork_nodes: List[int],
    join_nodes: List[int]
) -> List[List[int]]:
    """
    Detect parallel patterns in a workflow network.

    Args:
        link_matrix: Matrix containing workflow link information:
                     - Column 0: start node IDs
                     - Column 1: end node IDs
                     - Column 2: transition probabilities
        service_nodes: List of service node IDs
        fork_nodes: List of fork node IDs
        join_nodes: List of join node IDs

    Returns:
        List of parallel patterns, where each pattern is a list of parallel service nodes
    """
    parallel_patterns = []

    # Find all fork-join pairs
    fork_join_pairs = _find_fork_join_pairs(link_matrix, fork_nodes, join_nodes)

    for fork_node, join_node in fork_join_pairs:
        parallel_services = _find_parallel_services(link_matrix, service_nodes, fork_node, join_node)
        if len(parallel_services) > 1:
            parallel_patterns.append(parallel_services)

    return parallel_patterns


def _find_fork_join_pairs(
    link_matrix: np.ndarray,
    fork_nodes: List[int],
    join_nodes: List[int]
) -> List[Tuple[int, int]]:
    """Find all valid fork-join pairs in the workflow."""
    pairs = []
    fork_set = set(fork_nodes)
    join_set = set(join_nodes)

    # Build adjacency map for reachability analysis
    adjacency = _build_adjacency_map(link_matrix)

    for fork in fork_nodes:
        for join in join_nodes:
            if _is_valid_fork_join_pair(fork, join, adjacency, fork_set, join_set):
                pairs.append((fork, join))

    return pairs


def _build_adjacency_map(link_matrix: np.ndarray) -> Dict[int, List[int]]:
    """Build adjacency map from link matrix."""
    adjacency: Dict[int, List[int]] = {}

    if link_matrix.size == 0:
        return adjacency

    for i in range(link_matrix.shape[0]):
        start = int(link_matrix[i, 0])
        end = int(link_matrix[i, 1])

        if start not in adjacency:
            adjacency[start] = []
        adjacency[start].append(end)

    return adjacency


def _is_valid_fork_join_pair(
    fork: int,
    join: int,
    adjacency: Dict[int, List[int]],
    fork_set: Set[int],
    join_set: Set[int]
) -> bool:
    """Check if fork and join nodes form a valid pair."""
    # Use BFS to check if there are multiple paths from fork to join
    queue = deque([fork])
    visited: Set[int] = set()
    path_count: Dict[int, int] = {fork: 1}

    while queue:
        current = queue.popleft()
        if current in visited:
            continue
        visited.add(current)

        neighbors = adjacency.get(current, [])
        for neighbor in neighbors:
            if neighbor == join:
                # Found path to join
                current_paths = path_count.get(current, 0)
                path_count[join] = path_count.get(join, 0) + current_paths
            elif neighbor not in visited and neighbor not in fork_set and neighbor not in join_set:
                # Continue exploring (avoid other fork/join nodes)
                queue.append(neighbor)
                current_paths = path_count.get(current, 0)
                path_count[neighbor] = path_count.get(neighbor, 0) + current_paths

    # Valid if there are multiple paths to the join
    return path_count.get(join, 0) > 1


def _find_parallel_services(
    link_matrix: np.ndarray,
    service_nodes: List[int],
    fork_node: int,
    join_node: int
) -> List[int]:
    """Find service nodes that execute in parallel between fork and join."""
    parallel_services = []
    service_set = set(service_nodes)

    # Find all nodes reachable from fork
    reachable_from_fork = _find_reachable_nodes(link_matrix, fork_node, join_node)

    # Find all nodes that can reach join
    can_reach_join = _find_nodes_that_can_reach(link_matrix, join_node, fork_node)

    # Intersection gives nodes in parallel paths
    parallel_nodes = reachable_from_fork.intersection(can_reach_join)

    # Filter to only include service nodes
    for node in parallel_nodes:
        if node in service_set:
            parallel_services.append(node)

    return parallel_services


def _find_reachable_nodes(
    link_matrix: np.ndarray,
    start_node: int,
    end_node: int
) -> Set[int]:
    """Find all nodes reachable from a start node (stopping at end node)."""
    reachable: Set[int] = set()
    queue = deque([start_node])
    visited: Set[int] = set()

    if link_matrix.size == 0:
        return reachable

    while queue:
        current = queue.popleft()
        if current in visited or current == end_node:
            continue
        visited.add(current)

        for i in range(link_matrix.shape[0]):
            start = int(link_matrix[i, 0])
            end = int(link_matrix[i, 1])

            if start == current and end != end_node:
                reachable.add(end)
                queue.append(end)

    return reachable


def _find_nodes_that_can_reach(
    link_matrix: np.ndarray,
    target_node: int,
    start_node: int
) -> Set[int]:
    """Find all nodes that can reach a target node (starting from start node)."""
    can_reach: Set[int] = set()

    if link_matrix.size == 0:
        return can_reach

    # Build reverse adjacency map
    reverse_adj: Dict[int, List[int]] = {}
    for i in range(link_matrix.shape[0]):
        start = int(link_matrix[i, 0])
        end = int(link_matrix[i, 1])
        if end not in reverse_adj:
            reverse_adj[end] = []
        reverse_adj[end].append(start)

    # BFS backward from target
    queue = deque([target_node])
    visited: Set[int] = set()

    while queue:
        current = queue.popleft()
        if current in visited or current == start_node:
            continue
        visited.add(current)

        predecessors = reverse_adj.get(current, [])
        for pred in predecessors:
            if pred != start_node:
                can_reach.add(pred)
                queue.append(pred)

    return can_reach


def validate_parallel_pattern(
    pattern: List[int],
    link_matrix: np.ndarray,
    fork_nodes: List[int],
    join_nodes: List[int]
) -> bool:
    """
    Validate parallel pattern structure.

    Args:
        pattern: List of parallel node IDs
        link_matrix: Workflow link matrix
        fork_nodes: List of fork node IDs
        join_nodes: List of join node IDs

    Returns:
        True if valid parallel pattern
    """
    if len(pattern) < 2:
        return False

    # Find the fork and join nodes for this pattern
    fork_node = _find_common_source(pattern, link_matrix, fork_nodes)
    join_node = _find_common_target(pattern, link_matrix, join_nodes)

    return fork_node is not None and join_node is not None


def _find_common_source(
    pattern: List[int],
    link_matrix: np.ndarray,
    fork_nodes: List[int]
) -> Optional[int]:
    """Find common source (fork) for pattern nodes."""
    fork_set = set(fork_nodes)
    sources: Set[int] = set()

    if link_matrix.size == 0:
        return None

    for node in pattern:
        for i in range(link_matrix.shape[0]):
            start = int(link_matrix[i, 0])
            end = int(link_matrix[i, 1])

            if end == node and start in fork_set:
                sources.add(start)

    return next(iter(sources)) if len(sources) == 1 else None


def _find_common_target(
    pattern: List[int],
    link_matrix: np.ndarray,
    join_nodes: List[int]
) -> Optional[int]:
    """Find common target (join) for pattern nodes."""
    join_set = set(join_nodes)
    targets: Set[int] = set()

    if link_matrix.size == 0:
        return None

    for node in pattern:
        for i in range(link_matrix.shape[0]):
            start = int(link_matrix[i, 0])
            end = int(link_matrix[i, 1])

            if start == node and end in join_set:
                targets.add(end)

    return next(iter(targets)) if len(targets) == 1 else None


def get_parallel_stats(patterns: List[List[int]]) -> Dict[str, Any]:
    """
    Get parallel pattern statistics.

    Args:
        patterns: List of detected parallel patterns

    Returns:
        Dictionary with statistics
    """
    stats: Dict[str, Any] = {}

    stats['numPatterns'] = len(patterns)
    stats['totalParallelNodes'] = sum(len(p) for p in patterns)
    stats['avgParallelism'] = np.mean([len(p) for p in patterns]) if patterns else 0.0
    stats['maxParallelism'] = max(len(p) for p in patterns) if patterns else 0

    return stats
