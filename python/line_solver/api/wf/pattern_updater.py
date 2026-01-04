"""
Workflow Pattern Updater.

Implements dynamic pattern updating algorithms for workflow analysis.
Provides methods for refining and updating discovered workflow patterns
based on new observations and evolving process characteristics.

Based on AUTO_Pattern_Update.m from the MDN toolbox.

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/Wf_pattern_updater.kt
"""

import numpy as np
from typing import List, Tuple, Dict, Any, Set, Optional
from dataclasses import dataclass

from .sequence_detector import detect_sequences
from .parallel_detector import detect_parallel
from .loop_detector import detect_loops, get_loop_probability
from .branch_detector import detect_branches, BranchPattern


@dataclass
class ServiceParameters:
    """Data class to represent service parameters (phase-type distributions)."""
    alpha: np.ndarray  # Initial probability vector
    T: np.ndarray      # Transition rate matrix


@dataclass
class UpdatedWorkflow:
    """Data class to represent the updated workflow."""
    link_matrix: np.ndarray
    service_parameters: Dict[int, ServiceParameters]


def update_patterns(
    link_matrix: np.ndarray,
    service_nodes: List[int],
    fork_nodes: List[int],
    join_nodes: List[int],
    router_nodes: List[int],
    service_params: Dict[int, ServiceParameters]
) -> UpdatedWorkflow:
    """
    Update workflow network by simplifying detected patterns.

    Args:
        link_matrix: Original workflow link matrix
        service_nodes: List of service node IDs
        fork_nodes: List of fork node IDs
        join_nodes: List of join node IDs
        router_nodes: List of router node IDs
        service_params: Original service parameters (PH distributions)

    Returns:
        Updated workflow with simplified patterns
    """
    # Make mutable copies
    current_matrix = link_matrix.copy()
    current_params = dict(service_params)

    # Step 1: Simplify sequence patterns
    sequences = detect_sequences(current_matrix, service_nodes)
    current_matrix = _update_sequence_patterns(current_matrix, sequences, service_nodes, current_params)

    # Step 2: Simplify parallel patterns
    parallels = detect_parallel(current_matrix, service_nodes, fork_nodes, join_nodes)
    current_matrix = _update_parallel_patterns(current_matrix, parallels, service_nodes, current_params)

    # Step 3: Simplify loop patterns
    loops = detect_loops(current_matrix, service_nodes, router_nodes, join_nodes)
    current_matrix = _update_loop_patterns(current_matrix, loops, service_nodes, router_nodes, current_params)

    # Step 4: Simplify branch patterns
    branches = detect_branches(current_matrix, service_nodes, join_nodes)
    current_matrix = _update_branch_patterns(current_matrix, branches, service_nodes, current_params)

    return UpdatedWorkflow(current_matrix, current_params)


def _update_sequence_patterns(
    link_matrix: np.ndarray,
    sequences: List[List[int]],
    service_nodes: List[int],
    service_params: Dict[int, ServiceParameters]
) -> np.ndarray:
    """Update workflow by simplifying sequence patterns."""
    if not sequences:
        return link_matrix

    service_set = set(service_nodes)
    updated_matrix = link_matrix.copy()
    rows_to_remove = []

    # Find rows connecting service nodes in sequences
    for i in range(updated_matrix.shape[0]):
        start = int(updated_matrix[i, 0])
        end = int(updated_matrix[i, 1])

        if start in service_set and end in service_set:
            # Check if this connection is part of a sequence
            for sequence in sequences:
                for j in range(len(sequence) - 1):
                    if sequence[j] == start and sequence[j + 1] == end:
                        rows_to_remove.append(i)
                        break

    # Remove sequence connection rows
    final_matrix = _remove_matrix_rows(updated_matrix, rows_to_remove)

    # Update node references and parameters
    for sequence in sequences:
        if len(sequence) >= 2:
            first_node = sequence[0]
            last_node = sequence[-1]

            # Replace all occurrences of lastNode with firstNode
            _replace_node_references(final_matrix, last_node, first_node)

            # Convolve service parameters for sequence
            sequence_params = [service_params[n] for n in sequence if n in service_params]
            convolved_params = _convolve_sequence(sequence_params)
            service_params[first_node] = convolved_params

            # Remove parameters for other nodes in sequence
            for i in range(1, len(sequence)):
                if sequence[i] in service_params:
                    del service_params[sequence[i]]

    return final_matrix


def _update_parallel_patterns(
    link_matrix: np.ndarray,
    parallels: List[List[int]],
    service_nodes: List[int],
    service_params: Dict[int, ServiceParameters]
) -> np.ndarray:
    """Update workflow by simplifying parallel patterns."""
    if not parallels:
        return link_matrix

    updated_matrix = link_matrix.copy()

    for parallel in parallels:
        if len(parallel) >= 2:
            first_node = parallel[0]

            # Find fork and join nodes for this parallel pattern
            fork_node, join_node = _find_fork_join_for_parallel(updated_matrix, parallel)

            if fork_node is not None and join_node is not None:
                # Remove all connections involving parallel nodes
                rows_to_remove = _find_rows_involving_nodes(updated_matrix, parallel)
                updated_matrix = _remove_matrix_rows(updated_matrix, rows_to_remove)

                # Replace fork and join references with first parallel node
                _replace_node_references(updated_matrix, fork_node, first_node)
                _replace_node_references(updated_matrix, join_node, first_node)

                # Convolve service parameters for parallel execution
                parallel_params = [service_params[n] for n in parallel if n in service_params]
                convolved_params = _convolve_parallel(parallel_params)
                service_params[first_node] = convolved_params

                # Remove parameters for other nodes in parallel
                for i in range(1, len(parallel)):
                    if parallel[i] in service_params:
                        del service_params[parallel[i]]

    return updated_matrix


def _update_loop_patterns(
    link_matrix: np.ndarray,
    loops: List[int],
    service_nodes: List[int],
    router_nodes: List[int],
    service_params: Dict[int, ServiceParameters]
) -> np.ndarray:
    """Update workflow by simplifying loop patterns."""
    if not loops:
        return link_matrix

    updated_matrix = link_matrix.copy()
    router_set = set(router_nodes)
    ZERO = 1e-10

    for loop_node in loops:
        loop_prob = get_loop_probability(loop_node, updated_matrix, router_nodes)

        if loop_prob > ZERO:
            # Find router nodes involved in the loop
            involved_routers = _find_routers_for_loop(updated_matrix, loop_node, router_set)

            # Remove loop connections
            rows_to_remove = _find_loop_connections(updated_matrix, loop_node, involved_routers)
            updated_matrix = _remove_matrix_rows(updated_matrix, rows_to_remove)

            # Replace router references with loop node
            for router in involved_routers:
                _replace_node_references(updated_matrix, router, loop_node)

            # Update transition probabilities
            _update_loop_transition_probabilities(updated_matrix, loop_node)

            # Convolve service parameters with loop probability
            if loop_node in service_params:
                original_params = service_params[loop_node]
                looped_params = _convolve_loop(original_params, loop_prob)
                service_params[loop_node] = looped_params

    return updated_matrix


def _update_branch_patterns(
    link_matrix: np.ndarray,
    branches: List[BranchPattern],
    service_nodes: List[int],
    service_params: Dict[int, ServiceParameters]
) -> np.ndarray:
    """Update workflow by simplifying branch patterns."""
    if not branches:
        return link_matrix

    updated_matrix = link_matrix.copy()

    for branch in branches:
        if len(branch.branch_nodes) >= 2 and branch.fork_node is not None:
            first_node = branch.branch_nodes[0]

            # Remove branch connections
            rows_to_remove = _find_rows_involving_nodes(updated_matrix, branch.branch_nodes)
            updated_matrix = _remove_matrix_rows(updated_matrix, rows_to_remove)

            # Replace fork and join references
            _replace_node_references(updated_matrix, branch.fork_node, first_node)
            if branch.join_node is not None:
                _replace_node_references(updated_matrix, branch.join_node, first_node)

            # Convolve service parameters for branches
            branch_params = [service_params[n] for n in branch.branch_nodes if n in service_params]
            convolved_params = _convolve_branches(branch_params, branch.probabilities)
            service_params[first_node] = convolved_params

            # Remove parameters for other branch nodes
            for i in range(1, len(branch.branch_nodes)):
                if branch.branch_nodes[i] in service_params:
                    del service_params[branch.branch_nodes[i]]

    return updated_matrix


# Helper methods for matrix operations

def _remove_matrix_rows(matrix: np.ndarray, rows_to_remove: List[int]) -> np.ndarray:
    """Remove specified rows from matrix."""
    if not rows_to_remove or matrix.size == 0:
        return matrix

    sorted_rows = sorted(set(rows_to_remove), reverse=True)
    mask = np.ones(matrix.shape[0], dtype=bool)
    for row in sorted_rows:
        if 0 <= row < matrix.shape[0]:
            mask[row] = False

    return matrix[mask]


def _replace_node_references(matrix: np.ndarray, old_node: int, new_node: int) -> None:
    """Replace node references in matrix."""
    if matrix.size == 0:
        return

    for i in range(matrix.shape[0]):
        if int(matrix[i, 0]) == old_node:
            matrix[i, 0] = float(new_node)
        if int(matrix[i, 1]) == old_node:
            matrix[i, 1] = float(new_node)


def _find_rows_involving_nodes(matrix: np.ndarray, nodes: List[int]) -> List[int]:
    """Find rows involving specified nodes."""
    node_set = set(nodes)
    rows_to_remove = []

    if matrix.size == 0:
        return rows_to_remove

    for i in range(matrix.shape[0]):
        start = int(matrix[i, 0])
        end = int(matrix[i, 1])

        if start in node_set or end in node_set:
            rows_to_remove.append(i)

    return rows_to_remove


def _find_fork_join_for_parallel(
    matrix: np.ndarray,
    parallel: List[int]
) -> Tuple[Optional[int], Optional[int]]:
    """Find fork and join nodes for parallel pattern."""
    # Simplified implementation
    return (None, None)


def _find_routers_for_loop(
    matrix: np.ndarray,
    loop_node: int,
    router_set: Set[int]
) -> List[int]:
    """Find routers involved in loop."""
    routers = []

    if matrix.size == 0:
        return routers

    for i in range(matrix.shape[0]):
        start = int(matrix[i, 0])
        end = int(matrix[i, 1])

        if (start == loop_node and end in router_set) or \
           (end == loop_node and start in router_set):
            if start in router_set:
                routers.append(start)
            if end in router_set:
                routers.append(end)

    return list(set(routers))


def _find_loop_connections(
    matrix: np.ndarray,
    loop_node: int,
    routers: List[int]
) -> List[int]:
    """Find connections involved in loop."""
    connections = []
    router_set = set(routers)

    if matrix.size == 0:
        return connections

    for i in range(matrix.shape[0]):
        start = int(matrix[i, 0])
        end = int(matrix[i, 1])

        if (start == loop_node and end in router_set) or \
           (end == loop_node and start in router_set) or \
           (start in router_set and end == loop_node):
            connections.append(i)

    return connections


def _update_loop_transition_probabilities(matrix: np.ndarray, loop_node: int) -> None:
    """Update transition probabilities after loop simplification."""
    if matrix.size == 0:
        return

    for i in range(matrix.shape[0]):
        start = int(matrix[i, 0])
        if start == loop_node:
            matrix[i, 2] = 1.0  # Set probability to 1.0 after loop simplification


# Phase-type distribution convolution methods

def _convolve_sequence(params: List[ServiceParameters]) -> ServiceParameters:
    """Convolve service parameters for sequence (placeholder)."""
    if params:
        return params[0]  # Simplified - return first for now
    return ServiceParameters(np.ones((1, 1)), np.zeros((1, 1)))


def _convolve_parallel(params: List[ServiceParameters]) -> ServiceParameters:
    """Convolve service parameters for parallel execution (placeholder)."""
    if params:
        return params[0]  # Simplified - return first for now
    return ServiceParameters(np.ones((1, 1)), np.zeros((1, 1)))


def _convolve_loop(params: ServiceParameters, loop_prob: float) -> ServiceParameters:
    """Convolve service parameters for loops (placeholder)."""
    return params  # Simplified - return original for now


def _convolve_branches(params: List[ServiceParameters], probs: List[float]) -> ServiceParameters:
    """Convolve service parameters for branches (placeholder)."""
    if params:
        return params[0]  # Simplified - return first for now
    return ServiceParameters(np.ones((1, 1)), np.zeros((1, 1)))


def validate_updated_workflow(workflow: UpdatedWorkflow) -> bool:
    """
    Validate the updated workflow structure.

    Args:
        workflow: UpdatedWorkflow to validate

    Returns:
        True if valid
    """
    # Check that all referenced nodes have service parameters
    referenced_nodes: Set[int] = set()

    if workflow.link_matrix.size > 0:
        for i in range(workflow.link_matrix.shape[0]):
            referenced_nodes.add(int(workflow.link_matrix[i, 0]))
            referenced_nodes.add(int(workflow.link_matrix[i, 1]))

    # Check that all service nodes have parameters
    for node in referenced_nodes:
        if node not in workflow.service_parameters:
            return False

    return True


def get_update_stats(
    original_matrix: np.ndarray,
    updated_workflow: UpdatedWorkflow
) -> Dict[str, Any]:
    """
    Get statistics about the pattern update process.

    Args:
        original_matrix: Original workflow link matrix
        updated_workflow: Updated workflow

    Returns:
        Dictionary with statistics
    """
    stats: Dict[str, Any] = {}

    original_rows = original_matrix.shape[0] if original_matrix.size > 0 else 0
    updated_rows = updated_workflow.link_matrix.shape[0] if updated_workflow.link_matrix.size > 0 else 0

    stats['originalLinks'] = original_rows
    stats['updatedLinks'] = updated_rows
    stats['linksReduced'] = original_rows - updated_rows
    stats['serviceNodes'] = len(updated_workflow.service_parameters)

    reduction_ratio = (original_rows - updated_rows) / original_rows if original_rows > 0 else 0.0
    stats['reductionRatio'] = reduction_ratio

    return stats
