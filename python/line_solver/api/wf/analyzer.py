"""
Workflow Analysis Engine.

Provides comprehensive workflow analysis capabilities including pattern
recognition, performance analysis, and workflow model construction.
Integrates various detection algorithms for complete workflow characterization.

Based on Wf_analyzer.kt from the MDN toolbox.

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/Wf_analyzer.kt
"""

import numpy as np
from typing import List, Dict, Any, Optional, Set
from dataclasses import dataclass, field

from .sequence_detector import detect_sequences, get_sequence_stats, validate_sequence
from .parallel_detector import detect_parallel, get_parallel_stats, validate_parallel_pattern
from .loop_detector import detect_loops, get_loop_stats, validate_loop_pattern
from .branch_detector import detect_branches, get_branch_stats, validate_branch_pattern, BranchPattern, calculate_branch_diversity
from .pattern_updater import (
    update_patterns,
    get_update_stats,
    validate_updated_workflow,
    ServiceParameters,
    UpdatedWorkflow
)


@dataclass
class WorkflowRepresentation:
    """Data class to represent workflow in matrix form."""
    link_matrix: np.ndarray
    service_nodes: List[int]
    fork_nodes: List[int]
    join_nodes: List[int]
    router_nodes: List[int]
    service_parameters: Dict[int, ServiceParameters]


@dataclass
class DetectedPatterns:
    """Data class to hold all detected patterns."""
    sequences: List[List[int]]
    parallels: List[List[int]]
    loops: List[int]
    branches: List[BranchPattern]


@dataclass
class WorkflowAnalysis:
    """Data class to represent comprehensive workflow analysis results."""
    original_workflow: WorkflowRepresentation
    detected_patterns: DetectedPatterns
    optimized_workflow: UpdatedWorkflow
    statistics: Dict[str, Any]


class WfAnalyzer:
    """
    Main workflow analyzer that coordinates pattern detection and optimization.

    This class serves as the primary interface for workflow analysis, combining
    all pattern detectors and the pattern updater to provide comprehensive
    workflow optimization capabilities.

    Based on the AUTO solver functionality from the MDN toolbox.
    """

    def __init__(self, network=None, link_matrix: Optional[np.ndarray] = None):
        """
        Initialize workflow analyzer.

        Args:
            network: Optional LINE Network object
            link_matrix: Optional link matrix if network not provided
        """
        self._network = network
        self._link_matrix = link_matrix

    def analyze_workflow(self) -> WorkflowAnalysis:
        """
        Perform comprehensive workflow analysis and optimization.

        Returns:
            Complete workflow analysis results
        """
        # Step 1: Convert network to workflow representation
        workflow_rep = self._convert_network_to_workflow()

        # Step 2: Detect all patterns
        patterns = self._detect_all_patterns(workflow_rep)

        # Step 3: Optimize workflow by updating patterns
        optimized_workflow = update_patterns(
            workflow_rep.link_matrix,
            workflow_rep.service_nodes,
            workflow_rep.fork_nodes,
            workflow_rep.join_nodes,
            workflow_rep.router_nodes,
            workflow_rep.service_parameters
        )

        # Step 4: Generate statistics
        statistics = self._generate_analysis_statistics(workflow_rep, patterns, optimized_workflow)

        return WorkflowAnalysis(workflow_rep, patterns, optimized_workflow, statistics)

    def _convert_network_to_workflow(self) -> WorkflowRepresentation:
        """Convert LINE Network to workflow matrix representation."""
        if self._network is not None:
            link_matrix = self._build_link_matrix_from_network()
            node_classification = self._classify_nodes()
            service_params = self._extract_service_parameters(node_classification['service_nodes'])
        else:
            link_matrix = self._link_matrix if self._link_matrix is not None else np.zeros((0, 3))
            node_classification = {
                'service_nodes': [],
                'fork_nodes': [],
                'join_nodes': [],
                'router_nodes': []
            }
            service_params = {}

        return WorkflowRepresentation(
            link_matrix=link_matrix,
            service_nodes=node_classification['service_nodes'],
            fork_nodes=node_classification['fork_nodes'],
            join_nodes=node_classification['join_nodes'],
            router_nodes=node_classification['router_nodes'],
            service_parameters=service_params
        )

    def _build_link_matrix_from_network(self) -> np.ndarray:
        """Build link matrix from network routing matrix."""
        if self._network is None:
            return np.zeros((0, 3))

        links = []
        ZERO = 1e-10

        # Try to get routing matrix from network
        if hasattr(self._network, 'getLinkedRoutingMatrix'):
            routing_matrix = self._network.getLinkedRoutingMatrix()
            if routing_matrix:
                for from_class, class_routing in routing_matrix.items():
                    for to_class, matrix in class_routing.items():
                        for i in range(matrix.shape[0]):
                            for j in range(matrix.shape[1]):
                                prob = matrix[i, j]
                                if prob > ZERO:
                                    links.append((i, j, prob))
        elif hasattr(self._network, 'nodes'):
            # Create basic sequential connections
            nodes = self._network.nodes
            for i in range(len(nodes) - 1):
                links.append((i, i + 1, 1.0))

        # Create matrix representation
        if links:
            link_matrix = np.zeros((len(links), 3))
            for idx, (start, end, prob) in enumerate(links):
                link_matrix[idx, 0] = start
                link_matrix[idx, 1] = end
                link_matrix[idx, 2] = prob
        else:
            link_matrix = np.zeros((0, 3))

        return link_matrix

    def _classify_nodes(self) -> Dict[str, List[int]]:
        """Classify nodes by their type and function."""
        result = {
            'service_nodes': [],
            'fork_nodes': [],
            'join_nodes': [],
            'router_nodes': []
        }

        if self._network is None or not hasattr(self._network, 'nodes'):
            return result

        nodes = self._network.nodes

        for index, node in enumerate(nodes):
            node_type = type(node).__name__
            if node_type in ['Queue', 'Delay']:
                result['service_nodes'].append(index)
            elif node_type == 'Fork':
                result['fork_nodes'].append(index)
            elif node_type == 'Join':
                result['join_nodes'].append(index)
            elif node_type in ['Router', 'ClassSwitch']:
                result['router_nodes'].append(index)

        return result

    def _extract_service_parameters(
        self,
        service_nodes: List[int]
    ) -> Dict[int, ServiceParameters]:
        """Extract service parameters from service nodes."""
        params = {}

        if self._network is None or not hasattr(self._network, 'nodes'):
            return params

        nodes = self._network.nodes

        for node_index in service_nodes:
            if node_index < len(nodes):
                # Create default PH parameters (exponential distribution)
                params[node_index] = self._create_default_ph()

        return params

    def _create_default_ph(self) -> ServiceParameters:
        """Create default phase-type parameters (exponential distribution)."""
        alpha = np.ones((1, 1))
        T = np.array([[-1.0]])  # -1 rate matrix for exponential
        return ServiceParameters(alpha, T)

    def _detect_all_patterns(self, workflow: WorkflowRepresentation) -> DetectedPatterns:
        """Detect all workflow patterns."""
        sequences = detect_sequences(
            workflow.link_matrix,
            workflow.service_nodes
        )

        parallels = detect_parallel(
            workflow.link_matrix,
            workflow.service_nodes,
            workflow.fork_nodes,
            workflow.join_nodes
        )

        loops = detect_loops(
            workflow.link_matrix,
            workflow.service_nodes,
            workflow.router_nodes,
            workflow.join_nodes
        )

        branches = detect_branches(
            workflow.link_matrix,
            workflow.service_nodes,
            workflow.join_nodes
        )

        return DetectedPatterns(sequences, parallels, loops, branches)

    def _generate_analysis_statistics(
        self,
        original_workflow: WorkflowRepresentation,
        patterns: DetectedPatterns,
        optimized_workflow: UpdatedWorkflow
    ) -> Dict[str, Any]:
        """Generate comprehensive analysis statistics."""
        stats: Dict[str, Any] = {}

        # Pattern detection statistics
        stats['sequenceStats'] = get_sequence_stats(patterns.sequences)
        stats['parallelStats'] = get_parallel_stats(patterns.parallels)
        stats['loopStats'] = get_loop_stats(
            patterns.loops,
            original_workflow.link_matrix,
            original_workflow.router_nodes
        )
        stats['branchStats'] = get_branch_stats(patterns.branches)

        # Optimization statistics
        stats['updateStats'] = get_update_stats(
            original_workflow.link_matrix,
            optimized_workflow
        )

        # Overall complexity metrics
        stats['originalComplexity'] = self._calculate_workflow_complexity(original_workflow)
        stats['optimizedComplexity'] = self._calculate_optimized_complexity(optimized_workflow)

        return stats

    def _calculate_workflow_complexity(self, workflow: WorkflowRepresentation) -> Dict[str, Any]:
        """Calculate workflow complexity metrics."""
        complexity: Dict[str, Any] = {}

        complexity['totalNodes'] = (
            len(workflow.service_nodes) + len(workflow.fork_nodes) +
            len(workflow.join_nodes) + len(workflow.router_nodes)
        )
        complexity['totalLinks'] = workflow.link_matrix.shape[0] if workflow.link_matrix.size > 0 else 0
        complexity['serviceNodes'] = len(workflow.service_nodes)
        complexity['controlNodes'] = (
            len(workflow.fork_nodes) + len(workflow.join_nodes) + len(workflow.router_nodes)
        )

        # Calculate connectivity metrics
        node_set: Set[int] = set()
        if workflow.link_matrix.size > 0:
            for i in range(workflow.link_matrix.shape[0]):
                node_set.add(int(workflow.link_matrix[i, 0]))
                node_set.add(int(workflow.link_matrix[i, 1]))
        complexity['connectedNodes'] = len(node_set)

        # Calculate average degree
        degrees: Dict[int, int] = {}
        if workflow.link_matrix.size > 0:
            for i in range(workflow.link_matrix.shape[0]):
                start = int(workflow.link_matrix[i, 0])
                end = int(workflow.link_matrix[i, 1])
                degrees[start] = degrees.get(start, 0) + 1
                degrees[end] = degrees.get(end, 0) + 1
        complexity['avgDegree'] = np.mean(list(degrees.values())) if degrees else 0.0

        return complexity

    def _calculate_optimized_complexity(self, optimized_workflow: UpdatedWorkflow) -> Dict[str, Any]:
        """Calculate optimized workflow complexity."""
        complexity: Dict[str, Any] = {}

        complexity['totalNodes'] = len(optimized_workflow.service_parameters)
        complexity['totalLinks'] = (
            optimized_workflow.link_matrix.shape[0]
            if optimized_workflow.link_matrix.size > 0 else 0
        )

        # Calculate connectivity for optimized workflow
        node_set: Set[int] = set()
        if optimized_workflow.link_matrix.size > 0:
            for i in range(optimized_workflow.link_matrix.shape[0]):
                node_set.add(int(optimized_workflow.link_matrix[i, 0]))
                node_set.add(int(optimized_workflow.link_matrix[i, 1]))
        complexity['connectedNodes'] = len(node_set)

        return complexity

    def get_optimization_recommendations(self, analysis: WorkflowAnalysis) -> List[str]:
        """
        Get workflow optimization recommendations.

        Args:
            analysis: WorkflowAnalysis results

        Returns:
            List of optimization recommendations
        """
        recommendations = []

        patterns = analysis.detected_patterns
        update_stats = analysis.statistics.get('updateStats', {})

        # Sequence recommendations
        if patterns.sequences:
            recommendations.append(
                f"Found {len(patterns.sequences)} sequence patterns that can be simplified"
            )

        # Parallel recommendations
        if patterns.parallels:
            recommendations.append(
                f"Found {len(patterns.parallels)} parallel patterns for potential optimization"
            )

        # Loop recommendations
        if patterns.loops:
            recommendations.append(
                f"Found {len(patterns.loops)} loop patterns - consider loop unrolling for performance"
            )

        # Branch recommendations
        if patterns.branches:
            recommendations.append(
                f"Found {len(patterns.branches)} branch patterns - analyze probability distributions"
            )

            high_entropy_branches = [
                b for b in patterns.branches
                if calculate_branch_diversity(b).get('entropy', 0) > 1.0
            ]

            if high_entropy_branches:
                recommendations.append(
                    f"{len(high_entropy_branches)} branches have high entropy - consider load balancing"
                )

        # Optimization impact
        reduction_ratio = update_stats.get('reductionRatio', 0.0)
        if reduction_ratio > 0.1:
            recommendations.append(
                f"Workflow complexity reduced by {int(reduction_ratio * 100)}% through pattern optimization"
            )

        return recommendations

    def validate_analysis(self, analysis: WorkflowAnalysis) -> bool:
        """
        Validate workflow analysis results.

        Args:
            analysis: WorkflowAnalysis to validate

        Returns:
            True if valid
        """
        # Validate optimized workflow structure
        if not validate_updated_workflow(analysis.optimized_workflow):
            return False

        # Validate pattern consistency
        for sequence in analysis.detected_patterns.sequences:
            if not validate_sequence(sequence, analysis.original_workflow.link_matrix):
                return False

        for parallel in analysis.detected_patterns.parallels:
            if not validate_parallel_pattern(
                parallel,
                analysis.original_workflow.link_matrix,
                analysis.original_workflow.fork_nodes,
                analysis.original_workflow.join_nodes
            ):
                return False

        for loop in analysis.detected_patterns.loops:
            if not validate_loop_pattern(
                loop,
                analysis.original_workflow.link_matrix,
                analysis.original_workflow.router_nodes
            ):
                return False

        for branch in analysis.detected_patterns.branches:
            if not validate_branch_pattern(branch, analysis.original_workflow.link_matrix):
                return False

        return True
