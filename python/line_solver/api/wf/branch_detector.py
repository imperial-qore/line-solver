"""
Workflow Branch Pattern Detection.

Implements algorithms for detecting branching patterns in workflow traces.
Identifies decision points, conditional paths, and alternative execution
flows in business process analysis and workflow pattern recognition.

Based on AUTO_Branch_Detector.m from the MDN toolbox.

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/Wf_branch_detector.kt
"""

import numpy as np
from typing import List, Tuple, Dict, Any, Set, Optional
from dataclasses import dataclass
from collections import deque


@dataclass
class BranchPattern:
    """Data class to represent a branch pattern with probabilities."""
    branch_nodes: List[int]
    probabilities: List[float]
    fork_node: Optional[int]
    join_node: Optional[int]


def detect_branches(
    link_matrix: np.ndarray,
    service_nodes: List[int],
    join_nodes: List[int]
) -> List[BranchPattern]:
    """
    Detect branch patterns in a workflow network.

    Args:
        link_matrix: Matrix containing workflow link information:
                     - Column 0: start node IDs
                     - Column 1: end node IDs
                     - Column 2: transition probabilities
        service_nodes: List of service node IDs
        join_nodes: List of join node IDs

    Returns:
        List of branch patterns with their associated probabilities
    """
    branch_patterns = []
    service_set = set(service_nodes)
    join_set = set(join_nodes)

    # Build adjacency maps
    adjacency = _build_adjacency_map(link_matrix)
    reverse_adjacency = _build_reverse_adjacency_map(link_matrix)

    # Find nodes that have multiple outgoing edges (potential branch points)
    branch_points = _find_branch_points(adjacency, service_set, join_set)

    for branch_point in branch_points:
        pattern = _analyze_branch_pattern(
            branch_point,
            adjacency,
            reverse_adjacency,
            service_set,
            join_set
        )
        if pattern is not None and len(pattern.branch_nodes) > 1:
            branch_patterns.append(pattern)

    return branch_patterns


def _build_adjacency_map(link_matrix: np.ndarray) -> Dict[int, List[Tuple[int, float]]]:
    """Build forward adjacency map with probabilities."""
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


def _build_reverse_adjacency_map(link_matrix: np.ndarray) -> Dict[int, List[int]]:
    """Build reverse adjacency map."""
    reverse_adj: Dict[int, List[int]] = {}

    if link_matrix.size == 0:
        return reverse_adj

    for i in range(link_matrix.shape[0]):
        start = int(link_matrix[i, 0])
        end = int(link_matrix[i, 1])

        if end not in reverse_adj:
            reverse_adj[end] = []
        reverse_adj[end].append(start)

    return reverse_adj


def _find_branch_points(
    adjacency: Dict[int, List[Tuple[int, float]]],
    service_set: Set[int],
    join_set: Set[int]
) -> List[int]:
    """Find potential branch points (nodes with multiple outgoing edges)."""
    branch_points = []

    for node, neighbors in adjacency.items():
        if len(neighbors) > 1:
            # Check if this creates a valid branch pattern
            service_targets = sum(1 for target, _ in neighbors if target in service_set)
            if service_targets > 1:
                branch_points.append(node)

    return branch_points


def _analyze_branch_pattern(
    branch_point: int,
    adjacency: Dict[int, List[Tuple[int, float]]],
    reverse_adjacency: Dict[int, List[int]],
    service_set: Set[int],
    join_set: Set[int]
) -> Optional[BranchPattern]:
    """Analyze a potential branch pattern starting from a branch point."""
    neighbors = adjacency.get(branch_point)
    if neighbors is None:
        return None

    # Find service nodes that are direct targets of the branch
    branch_targets: List[Tuple[int, float]] = []
    for target, prob in neighbors:
        if target in service_set:
            branch_targets.append((target, prob))

    if len(branch_targets) < 2:
        return None

    # Find common join point for these branches
    common_join = _find_common_join_point(
        [t[0] for t in branch_targets],
        adjacency,
        reverse_adjacency,
        join_set
    )

    # Validate probabilities sum to 1.0 (or close to it)
    total_prob = sum(t[1] for t in branch_targets)
    if abs(total_prob - 1.0) > 0.01:
        # Not a valid probabilistic branch
        return None

    return BranchPattern(
        branch_nodes=[t[0] for t in branch_targets],
        probabilities=[t[1] for t in branch_targets],
        fork_node=branch_point,
        join_node=common_join
    )


def _find_common_join_point(
    branch_nodes: List[int],
    adjacency: Dict[int, List[Tuple[int, float]]],
    reverse_adjacency: Dict[int, List[int]],
    join_set: Set[int]
) -> Optional[int]:
    """Find the common join point for a set of branch nodes."""
    # Find all nodes reachable from each branch
    reachable_sets = [
        _find_reachable_nodes(branch, adjacency, join_set)
        for branch in branch_nodes
    ]

    # Find intersection of all reachable sets
    common_reachable = reachable_sets[0]
    for rs in reachable_sets[1:]:
        common_reachable = common_reachable.intersection(rs)

    # Prefer join nodes as common points
    join_points = common_reachable.intersection(join_set)
    if join_points:
        # Return the "closest" join point (heuristic: first in set)
        return next(iter(join_points))

    # If no join nodes, return any common reachable node
    return next(iter(common_reachable), None)


def _find_reachable_nodes(
    start_node: int,
    adjacency: Dict[int, List[Tuple[int, float]]],
    stop_set: Set[int]
) -> Set[int]:
    """Find all nodes reachable from a starting node."""
    reachable: Set[int] = set()
    queue = deque([start_node])
    visited: Set[int] = set()

    while queue:
        current = queue.popleft()
        if current in visited:
            continue
        visited.add(current)

        neighbors = adjacency.get(current, [])
        for neighbor, _ in neighbors:
            reachable.add(neighbor)
            if neighbor not in stop_set:
                queue.append(neighbor)

    return reachable


def validate_branch_pattern(pattern: BranchPattern, link_matrix: np.ndarray) -> bool:
    """
    Validate branch pattern structure.

    Args:
        pattern: BranchPattern to validate
        link_matrix: Workflow link matrix

    Returns:
        True if pattern is valid
    """
    # Check probability sum
    total_prob = sum(pattern.probabilities)
    if abs(total_prob - 1.0) > 0.01:
        return False

    # Check that all branch nodes are reachable from fork
    if pattern.fork_node is None:
        return False

    adjacency = _build_adjacency_map(link_matrix)
    fork_neighbors = adjacency.get(pattern.fork_node)
    if fork_neighbors is None:
        return False

    fork_targets = {t[0] for t in fork_neighbors}
    return all(bn in fork_targets for bn in pattern.branch_nodes)


def calculate_branch_diversity(pattern: BranchPattern) -> Dict[str, float]:
    """
    Calculate branch diversity metrics.

    Args:
        pattern: BranchPattern to analyze

    Returns:
        Dictionary with metrics:
            - entropy: Shannon entropy
            - normalizedEntropy: Normalized entropy
            - gini: Gini coefficient
            - balance: Balance measure
    """
    metrics: Dict[str, float] = {}

    probs = pattern.probabilities
    n = len(probs)

    # Shannon entropy
    entropy = -sum(p * np.log(p) if p > 0 else 0.0 for p in probs)
    metrics['entropy'] = float(entropy)

    # Normalized entropy
    metrics['normalizedEntropy'] = float(entropy / np.log(n)) if n > 1 else 0.0

    # Gini coefficient (inequality measure)
    sorted_probs = sorted(probs)
    gini = 0.0
    for i, p in enumerate(sorted_probs):
        gini += (2 * (i + 1) - n - 1) * p
    gini /= (n - 1) * sum(probs) if n > 1 and sum(probs) > 0 else 1
    metrics['gini'] = abs(float(gini))

    # Balance (inverse of max probability)
    max_prob = max(probs) if probs else 1.0
    metrics['balance'] = 1.0 / max_prob if max_prob > 0 else 0.0

    return metrics


def get_branch_stats(patterns: List[BranchPattern]) -> Dict[str, Any]:
    """
    Get branch pattern statistics.

    Args:
        patterns: List of detected branch patterns

    Returns:
        Dictionary with statistics
    """
    stats: Dict[str, Any] = {}

    stats['numPatterns'] = len(patterns)
    stats['totalBranchNodes'] = sum(len(p.branch_nodes) for p in patterns)

    branch_counts = [len(p.branch_nodes) for p in patterns]
    stats['avgBranches'] = np.mean(branch_counts) if branch_counts else 0.0
    stats['maxBranches'] = max(branch_counts) if branch_counts else 0
    stats['minBranches'] = min(branch_counts) if branch_counts else 0

    # Diversity statistics
    diversity_metrics = [calculate_branch_diversity(p) for p in patterns]
    stats['avgEntropy'] = np.mean([d['entropy'] for d in diversity_metrics]) if diversity_metrics else 0.0
    stats['avgBalance'] = np.mean([d['balance'] for d in diversity_metrics]) if diversity_metrics else 0.0

    return stats


def find_most_probable_branch(pattern: BranchPattern) -> Optional[Tuple[int, float]]:
    """
    Find the most probable branch in a pattern.

    Args:
        pattern: BranchPattern to analyze

    Returns:
        Tuple of (node_id, probability) or None if empty
    """
    if not pattern.branch_nodes:
        return None

    max_idx = np.argmax(pattern.probabilities)
    return (pattern.branch_nodes[max_idx], pattern.probabilities[max_idx])


def find_least_probable_branch(pattern: BranchPattern) -> Optional[Tuple[int, float]]:
    """
    Find the least probable branch in a pattern.

    Args:
        pattern: BranchPattern to analyze

    Returns:
        Tuple of (node_id, probability) or None if empty
    """
    if not pattern.branch_nodes:
        return None

    min_idx = np.argmin(pattern.probabilities)
    return (pattern.branch_nodes[min_idx], pattern.probabilities[min_idx])
