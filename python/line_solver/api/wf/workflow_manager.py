"""
Workflow Management System.

Provides workflow management capabilities for queueing network analysis,
integrating pattern detection and automatic workflow construction with
LINE solver framework for complex system modeling and optimization.

Based on WorkflowManager.kt from the MDN toolbox.

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/WorkflowManager.kt
"""

import numpy as np
import time
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

from .analyzer import WfAnalyzer, WorkflowAnalysis, DetectedPatterns
from .auto_integration import WfAutoIntegration, ExtendedSolverRecommendation
from .branch_detector import calculate_branch_diversity


@dataclass
class WorkflowAnalysisResult:
    """Comprehensive workflow analysis result."""
    pattern_analysis: WorkflowAnalysis
    solver_recommendation: ExtendedSolverRecommendation
    optimization_insights: Dict[str, Any]
    performance_metrics: Dict[str, float]


class WorkflowManager:
    """
    Main facade for workflow management and optimization in LINE.

    This class provides a high-level interface for workflow analysis,
    pattern detection, optimization, and intelligent solver selection
    based on the AUTO workflow analysis algorithms from the MDN toolbox.

    Usage example:
        manager = WorkflowManager(network)
        analysis = manager.analyze_workflow()
        recommendations = manager.get_optimization_recommendations()
    """

    def __init__(self, network=None, options: Optional[Dict[str, Any]] = None):
        """
        Initialize WorkflowManager.

        Args:
            network: LINE Network object
            options: Optional solver options
        """
        self._network = network
        self._options = options or {}
        self._workflow_analyzer = WfAnalyzer(network)
        self._auto_integration = WfAutoIntegration(network, options)

    def analyze_workflow(self) -> WorkflowAnalysisResult:
        """
        Perform comprehensive workflow analysis.

        Returns:
            Complete analysis including pattern detection, solver recommendation,
            and optimization insights
        """
        # Perform pattern analysis
        pattern_analysis = self._workflow_analyzer.analyze_workflow()

        # Get solver recommendation with workflow insights
        solver_recommendation = self._auto_integration.recommend_solver_with_workflow_analysis()

        # Get optimization insights
        optimization_insights = self._auto_integration.get_optimization_insights()

        # Calculate performance metrics
        performance_metrics = self._calculate_performance_metrics(
            pattern_analysis,
            solver_recommendation
        )

        return WorkflowAnalysisResult(
            pattern_analysis,
            solver_recommendation,
            optimization_insights,
            performance_metrics
        )

    def get_pattern_analysis(self) -> DetectedPatterns:
        """
        Get detailed pattern detection results.

        Returns:
            Detected workflow patterns with statistics
        """
        analysis = self._workflow_analyzer.analyze_workflow()
        return analysis.detected_patterns

    def get_optimization_recommendations(self) -> List[str]:
        """
        Get workflow optimization recommendations.

        Returns:
            List of actionable optimization recommendations
        """
        analysis = self._workflow_analyzer.analyze_workflow()
        return self._workflow_analyzer.get_optimization_recommendations(analysis)

    def generate_complexity_report(self) -> Dict[str, Any]:
        """
        Generate a workflow complexity report.

        Returns:
            Detailed complexity analysis
        """
        analysis = self._workflow_analyzer.analyze_workflow()
        report: Dict[str, Any] = {}

        # Basic complexity metrics
        original_complexity = analysis.statistics.get('originalComplexity', {})
        optimized_complexity = analysis.statistics.get('optimizedComplexity', {})

        report['originalMetrics'] = original_complexity
        report['optimizedMetrics'] = optimized_complexity

        # Pattern complexity breakdown
        patterns = analysis.detected_patterns
        pattern_complexity: Dict[str, Any] = {}

        pattern_complexity['sequences'] = {
            'count': len(patterns.sequences),
            'totalNodes': sum(len(s) for s in patterns.sequences),
            'complexity': sum(len(s) ** 2 for s in patterns.sequences)  # Quadratic complexity estimate
        }

        pattern_complexity['parallels'] = {
            'count': len(patterns.parallels),
            'totalNodes': sum(len(p) for p in patterns.parallels),
            'complexity': sum(self._factorial(len(p)) for p in patterns.parallels)  # Exponential for fork-join
        }

        pattern_complexity['loops'] = {
            'count': len(patterns.loops),
            'complexity': len(patterns.loops) * 10  # Loop complexity multiplier
        }

        pattern_complexity['branches'] = {
            'count': len(patterns.branches),
            'totalBranches': sum(len(b.branch_nodes) for b in patterns.branches),
            'complexity': sum(len(b.branch_nodes) for b in patterns.branches)
        }

        report['patternComplexity'] = pattern_complexity

        # Overall complexity score
        original_nodes = original_complexity.get('totalNodes', 0)
        original_links = original_complexity.get('totalLinks', 0)
        complexity_score = self._calculate_complexity_score(original_nodes, original_links, patterns)

        report['overallComplexityScore'] = complexity_score
        if complexity_score < 10:
            report['complexityLevel'] = 'LOW'
        elif complexity_score < 50:
            report['complexityLevel'] = 'MEDIUM'
        elif complexity_score < 100:
            report['complexityLevel'] = 'HIGH'
        else:
            report['complexityLevel'] = 'VERY_HIGH'

        return report

    def benchmark_solvers(
        self,
        solvers: Optional[List[str]] = None
    ) -> Dict[str, Dict[str, Any]]:
        """
        Benchmark different solvers on this workflow.

        Args:
            solvers: List of solver names to benchmark

        Returns:
            Performance comparison results
        """
        if solvers is None:
            solvers = ['MVA', 'NC', 'SSA', 'FLUID']

        results: Dict[str, Dict[str, Any]] = {}

        for solver_name in solvers:
            try:
                solver = self._create_solver(solver_name)
                if solver is None:
                    results[solver_name] = {
                        'success': False,
                        'error': 'Solver not available'
                    }
                    continue

                start_time = time.time()

                # Attempt to solve
                solver_result = solver.getAvg() if hasattr(solver, 'getAvg') else None
                end_time = time.time()

                metrics: Dict[str, Any] = {}
                metrics['success'] = True
                metrics['solveTime'] = (end_time - start_time) * 1000  # ms
                metrics['hasResults'] = solver_result is not None

                if solver_result is not None and hasattr(solver_result, 'QN') and solver_result.QN is not None:
                    metrics['totalQueueLength'] = float(np.sum(solver_result.QN))
                    metrics['avgQueueLength'] = float(np.mean(solver_result.QN))

                results[solver_name] = metrics

            except Exception as e:
                results[solver_name] = {
                    'success': False,
                    'error': str(e)
                }

        return results

    def export_analysis(self, format: str = 'SUMMARY') -> str:
        """
        Export workflow analysis to different formats.

        Args:
            format: Export format ('JSON', 'CSV', 'SUMMARY')

        Returns:
            Formatted analysis data
        """
        analysis = self.analyze_workflow()

        format_upper = format.upper()
        if format_upper == 'JSON':
            return self._export_to_json(analysis)
        elif format_upper == 'CSV':
            return self._export_to_csv(analysis)
        elif format_upper == 'SUMMARY':
            return self._export_to_summary(analysis)
        else:
            raise ValueError(f"Unsupported export format: {format}")

    def validate_workflow(self) -> Dict[str, Any]:
        """
        Validate workflow structure and analysis results.

        Returns:
            Validation results with any issues found
        """
        validation: Dict[str, Any] = {}
        issues: List[str] = []

        try:
            # Check network structure
            if self._network is not None:
                if hasattr(self._network, 'nodes') and len(self._network.nodes) == 0:
                    issues.append("Network has no nodes")
            else:
                issues.append("No network provided")

            # Validate workflow analysis
            analysis = self._workflow_analyzer.analyze_workflow()
            if not self._workflow_analyzer.validate_analysis(analysis):
                issues.append("Workflow analysis validation failed")

            # Validate solver integration
            if not self._auto_integration.validate_workflow_enhancement():
                issues.append("Workflow-enhanced solver selection validation failed")

            validation['isValid'] = len(issues) == 0
            validation['issues'] = issues

            if self._network is not None and hasattr(self._network, 'nodes'):
                validation['nodeCount'] = len(self._network.nodes)
            else:
                validation['nodeCount'] = 0

            validation['hasRouting'] = False
            if self._network is not None and hasattr(self._network, 'getLinkedRoutingMatrix'):
                routing = self._network.getLinkedRoutingMatrix()
                validation['hasRouting'] = routing is not None and len(routing) > 0

        except Exception as e:
            validation['isValid'] = False
            validation['issues'] = [f"Validation error: {str(e)}"]

        return validation

    # Private helper methods

    def _calculate_performance_metrics(
        self,
        analysis: WorkflowAnalysis,
        recommendation: ExtendedSolverRecommendation
    ) -> Dict[str, float]:
        """Calculate performance metrics."""
        metrics: Dict[str, float] = {}

        # Complexity reduction metric
        update_stats = analysis.statistics.get('updateStats', {})
        reduction_ratio = update_stats.get('reductionRatio', 0.0)
        metrics['complexityReduction'] = float(reduction_ratio)

        # Solver confidence
        metrics['solverConfidence'] = recommendation.confidence

        # Pattern efficiency scores
        patterns = analysis.detected_patterns
        metrics['sequenceEfficiency'] = self._calculate_sequence_efficiency(patterns.sequences)
        metrics['parallelEfficiency'] = self._calculate_parallel_efficiency(patterns.parallels)
        metrics['loopEfficiency'] = self._calculate_loop_efficiency(
            patterns.loops,
            analysis.original_workflow.link_matrix
        )
        metrics['branchEfficiency'] = self._calculate_branch_efficiency(patterns.branches)

        return metrics

    def _calculate_sequence_efficiency(self, sequences: List[List[int]]) -> float:
        """Calculate sequence efficiency score."""
        if not sequences:
            return 1.0

        total_length = sum(len(s) for s in sequences)
        avg_length = total_length / len(sequences)

        # Longer sequences are more efficient to optimize
        return min(1.0, avg_length / 5.0)

    def _calculate_parallel_efficiency(self, parallels: List[List[int]]) -> float:
        """Calculate parallel efficiency score."""
        if not parallels:
            return 1.0

        total_parallelism = sum(len(p) for p in parallels)
        avg_parallelism = total_parallelism / len(parallels)

        # Higher parallelism can be more efficient but also more complex
        return max(0.1, 1.0 - (avg_parallelism - 2.0) / 10.0)

    def _calculate_loop_efficiency(
        self,
        loops: List[int],
        link_matrix: np.ndarray
    ) -> float:
        """Calculate loop efficiency score."""
        if not loops:
            return 1.0

        # Lower loop probabilities are more efficient
        avg_prob = 0.5  # Simplified estimate
        return 1.0 - avg_prob

    def _calculate_branch_efficiency(self, branches: List) -> float:
        """Calculate branch efficiency score."""
        if not branches:
            return 1.0

        entropies = [
            calculate_branch_diversity(b).get('normalizedEntropy', 0.0)
            for b in branches
        ]

        # Higher entropy means more balanced branches, which is more efficient
        return np.mean(entropies) if entropies else 0.0

    def _calculate_complexity_score(
        self,
        nodes: int,
        links: int,
        patterns: DetectedPatterns
    ) -> float:
        """Calculate overall complexity score."""
        score = float(nodes) + links * 0.5

        # Add pattern complexity
        score += sum(len(s) ** 2 for s in patterns.sequences) * 0.1
        score += sum(self._factorial(len(p)) for p in patterns.parallels) * 0.2
        score += len(patterns.loops) * 10
        score += sum(len(b.branch_nodes) for b in patterns.branches) * 2

        return score

    def _factorial(self, n: int) -> int:
        """Calculate factorial."""
        if n <= 1:
            return 1
        return n * self._factorial(n - 1)

    def _create_solver(self, solver_name: str):
        """Create solver based on name."""
        if self._network is None:
            return None

        # This would integrate with LINE's solver factory
        # Placeholder implementation
        return None

    def _export_to_json(self, analysis: WorkflowAnalysisResult) -> str:
        """Export analysis to JSON format."""
        import json

        data = {
            'solver_recommendation': analysis.solver_recommendation.recommended_solver,
            'confidence': analysis.solver_recommendation.confidence,
            'patterns': {
                'sequences': len(analysis.pattern_analysis.detected_patterns.sequences),
                'parallels': len(analysis.pattern_analysis.detected_patterns.parallels),
                'loops': len(analysis.pattern_analysis.detected_patterns.loops),
                'branches': len(analysis.pattern_analysis.detected_patterns.branches)
            }
        }
        return json.dumps(data, indent=2)

    def _export_to_csv(self, analysis: WorkflowAnalysisResult) -> str:
        """Export analysis to CSV format."""
        lines = [
            'Metric,Value',
            f'Recommended Solver,{analysis.solver_recommendation.recommended_solver}',
            f'Confidence,{analysis.solver_recommendation.confidence}',
            f'Sequences,{len(analysis.pattern_analysis.detected_patterns.sequences)}',
            f'Parallels,{len(analysis.pattern_analysis.detected_patterns.parallels)}',
            f'Loops,{len(analysis.pattern_analysis.detected_patterns.loops)}',
            f'Branches,{len(analysis.pattern_analysis.detected_patterns.branches)}'
        ]
        return '\n'.join(lines)

    def _export_to_summary(self, analysis: WorkflowAnalysisResult) -> str:
        """Export analysis to summary format."""
        lines = [
            '=== Workflow Analysis Summary ===',
            '',
            f'Recommended Solver: {analysis.solver_recommendation.recommended_solver}',
            f'Confidence: {analysis.solver_recommendation.confidence:.2f}',
            '',
            'Detected Patterns:',
            f'- Sequences: {len(analysis.pattern_analysis.detected_patterns.sequences)}',
            f'- Parallels: {len(analysis.pattern_analysis.detected_patterns.parallels)}',
            f'- Loops: {len(analysis.pattern_analysis.detected_patterns.loops)}',
            f'- Branches: {len(analysis.pattern_analysis.detected_patterns.branches)}',
            '',
            'Reasoning:'
        ]

        for i, reason in enumerate(analysis.solver_recommendation.reasoning, 1):
            lines.append(f'{i}. {reason}')

        return '\n'.join(lines)


def quick_analysis(network) -> str:
    """
    Quick analysis method for simple workflow inspection.

    Args:
        network: LINE Network object

    Returns:
        Summary string
    """
    manager = WorkflowManager(network)
    return manager.export_analysis('SUMMARY')


def get_optimal_solver(network, options: Optional[Dict[str, Any]] = None):
    """
    Get optimal solver for a network without detailed analysis.

    Args:
        network: LINE Network object
        options: Optional solver options

    Returns:
        Recommended solver name
    """
    manager = WorkflowManager(network, options)
    analysis = manager.analyze_workflow()
    return analysis.solver_recommendation.recommended_solver
