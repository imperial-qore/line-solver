"""
Workflow Sequence Pattern Detection.

Implements algorithms for detecting sequential patterns in workflow traces.
Identifies ordered sequences of activities and their temporal relationships
for workflow analysis and pattern recognition applications.

Based on AUTO_Sequence_Detector.m from the MDN toolbox.

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/Wf_sequence_detector.kt
"""

import numpy as np
from typing import List, Tuple, Dict, Any, Set


def detect_sequences(link_matrix: np.ndarray, service_nodes: List[int]) -> List[List[int]]:
    """
    Detect sequence patterns in a workflow network.

    Args:
        link_matrix: Matrix containing workflow link information:
                     - Column 0: start node IDs
                     - Column 1: end node IDs
                     - Column 2: transition probabilities
        service_nodes: List of service node IDs

    Returns:
        List of sequence chains, where each chain is a list of connected service nodes
    """
    chains = []

    # Step 1: Find connections between service nodes
    service_connections = _find_service_connections(link_matrix, service_nodes)
    if not service_connections:
        return chains  # No sequence structures found

    # Step 2: Count occurrences of each service node
    node_counts = _count_node_occurrences(service_connections)

    # Step 3: Detect sequences using node occurrence patterns
    num_sequences = sum(1 for count in node_counts.values() if count == 1) // 2

    remaining_connections = list(service_connections)

    for _ in range(num_sequences):
        if not remaining_connections:
            break

        sequence = _build_sequence_chain(remaining_connections)
        if sequence:
            chains.append(sequence)

    return chains


def _find_service_connections(link_matrix: np.ndarray, service_nodes: List[int]) -> List[Tuple[int, int]]:
    """Find all connections between service nodes in the workflow."""
    connections = []
    service_set = set(service_nodes)

    if link_matrix.size == 0:
        return connections

    for i in range(link_matrix.shape[0]):
        start_node = int(link_matrix[i, 0])
        end_node = int(link_matrix[i, 1])

        # Check if both nodes are service nodes
        if start_node in service_set and end_node in service_set:
            connections.append((start_node, end_node))

    return connections


def _count_node_occurrences(connections: List[Tuple[int, int]]) -> Dict[int, int]:
    """Count how many times each service node appears in connections."""
    counts: Dict[int, int] = {}

    for start, end in connections:
        counts[start] = counts.get(start, 0) + 1
        counts[end] = counts.get(end, 0) + 1

    return counts


def _build_sequence_chain(connections: List[Tuple[int, int]]) -> List[int]:
    """Build a sequence chain starting from the first available connection."""
    if not connections:
        return []

    sequence = []
    used_indices = []

    # Start with first connection
    first, last = connections[0]
    sequence.append(first)
    sequence.append(last)
    used_indices.append(0)

    found_extension = True
    while found_extension:
        found_extension = False
        current_size = len(sequence)

        # Try to extend forward or backward
        for i in range(1, len(connections)):
            if i in used_indices:
                continue

            start, end = connections[i]

            if start == last:
                # Extend forward
                last = end
                sequence.append(end)
                used_indices.append(i)
                found_extension = True
            elif end == first:
                # Extend backward
                first = start
                sequence.insert(0, start)
                used_indices.append(i)
                found_extension = True

        # Check if we made progress
        found_extension = found_extension and len(sequence) > current_size

    # Remove used connections from the list
    for index in sorted(used_indices, reverse=True):
        connections.pop(index)

    return sequence


def validate_sequence(sequence: List[int], link_matrix: np.ndarray) -> bool:
    """
    Validate that a sequence chain is properly connected.

    Args:
        sequence: List of node IDs in the sequence
        link_matrix: Workflow link matrix

    Returns:
        True if all consecutive pairs are connected in the link matrix
    """
    if len(sequence) < 2:
        return False

    connections: Set[Tuple[int, int]] = set()
    for i in range(link_matrix.shape[0]):
        start = int(link_matrix[i, 0])
        end = int(link_matrix[i, 1])
        connections.add((start, end))

    # Check all consecutive pairs in sequence
    for i in range(len(sequence) - 1):
        connection = (sequence[i], sequence[i + 1])
        if connection not in connections:
            return False

    return True


def get_sequence_stats(sequences: List[List[int]]) -> Dict[str, Any]:
    """
    Get sequence statistics for analysis.

    Args:
        sequences: List of detected sequences

    Returns:
        Dictionary with statistics:
            - numSequences: Number of sequences
            - totalNodes: Total nodes across all sequences
            - avgLength: Average sequence length
            - maxLength: Maximum sequence length
            - minLength: Minimum sequence length
    """
    stats: Dict[str, Any] = {}

    stats['numSequences'] = len(sequences)
    stats['totalNodes'] = sum(len(seq) for seq in sequences)
    stats['avgLength'] = np.mean([len(seq) for seq in sequences]) if sequences else 0.0
    stats['maxLength'] = max(len(seq) for seq in sequences) if sequences else 0
    stats['minLength'] = min(len(seq) for seq in sequences) if sequences else 0

    return stats
