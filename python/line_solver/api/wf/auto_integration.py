"""
Workflow Automatic Integration.

Provides automatic integration capabilities for workflow analysis with
LINE solver framework. Enables seamless workflow model construction and
analysis through automated pattern detection and model generation.

Based on Wf_auto_integration.kt from the MDN toolbox.

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/Wf_auto_integration.kt
"""

import numpy as np
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

from .analyzer import WfAnalyzer, WorkflowAnalysis, DetectedPatterns
from .branch_detector import calculate_branch_diversity


@dataclass
class ExtendedSolverRecommendation:
    """Data class for extended solver recommendations."""
    recommended_solver: str
    confidence: float
    reasoning: List[str]
    workflow_features: Dict[str, Any]
    alternative_solvers: List[str]


class WfAutoIntegration:
    """
    Integration class connecting workflow analysis to the AUTO solver.

    This class extends the AUTO solver's capabilities by incorporating
    workflow pattern analysis and optimization techniques from the MDN toolbox.
    It provides intelligent solver selection based on detected workflow patterns.
    """

    def __init__(self, network=None, options: Optional[Dict[str, Any]] = None):
        """
        Initialize workflow auto integration.

        Args:
            network: LINE Network object
            options: Solver options dictionary
        """
        self._network = network
        self._options = options or {}
        self._workflow_analyzer = WfAnalyzer(network)

    def recommend_solver_with_workflow_analysis(self) -> ExtendedSolverRecommendation:
        """
        Perform workflow-aware solver selection.

        Returns:
            Extended solver recommendation with workflow analysis
        """
        # Perform workflow analysis
        workflow_analysis = self._workflow_analyzer.analyze_workflow()

        # Extract workflow features for solver selection
        workflow_features = self._extract_workflow_features(workflow_analysis)

        # Get base recommendation from model analyzer
        base_recommendation = self._get_base_recommendation()

        # Enhance recommendation with workflow insights
        enhanced_recommendation = self._enhance_recommendation_with_workflow(
            base_recommendation,
            workflow_features,
            workflow_analysis
        )

        return enhanced_recommendation

    def _extract_workflow_features(
        self,
        analysis: WorkflowAnalysis
    ) -> Dict[str, Any]:
        """Extract workflow features for solver selection."""
        features: Dict[str, Any] = {}

        patterns = analysis.detected_patterns
        original_complexity = analysis.statistics.get('originalComplexity', {})
        optimized_complexity = analysis.statistics.get('optimizedComplexity', {})

        # Pattern-based features
        features['hasSequencePatterns'] = len(patterns.sequences) > 0
        features['hasParallelPatterns'] = len(patterns.parallels) > 0
        features['hasLoopPatterns'] = len(patterns.loops) > 0
        features['hasBranchPatterns'] = len(patterns.branches) > 0

        features['numSequences'] = len(patterns.sequences)
        features['numParallels'] = len(patterns.parallels)
        features['numLoops'] = len(patterns.loops)
        features['numBranches'] = len(patterns.branches)

        # Complexity features
        features['originalNodeCount'] = original_complexity.get('totalNodes', 0)
        features['originalLinkCount'] = original_complexity.get('totalLinks', 0)
        features['optimizedNodeCount'] = optimized_complexity.get('totalNodes', 0)
        features['optimizedLinkCount'] = optimized_complexity.get('totalLinks', 0)

        # Pattern complexity metrics
        if patterns.sequences:
            seq_lengths = [len(s) for s in patterns.sequences]
            features['avgSequenceLength'] = np.mean(seq_lengths)
            features['maxSequenceLength'] = max(seq_lengths)

        if patterns.parallels:
            par_sizes = [len(p) for p in patterns.parallels]
            features['avgParallelism'] = np.mean(par_sizes)
            features['maxParallelism'] = max(par_sizes)

        if patterns.loops:
            loop_stats = analysis.statistics.get('loopStats', {})
            features['avgLoopProbability'] = loop_stats.get('avgLoopProbability', 0.0)
            features['maxLoopProbability'] = loop_stats.get('maxLoopProbability', 0.0)

        if patterns.branches:
            branch_stats = analysis.statistics.get('branchStats', {})
            features['avgBranches'] = branch_stats.get('avgBranches', 0.0)
            features['maxBranches'] = branch_stats.get('maxBranches', 0)

            # Calculate average entropy across all branches
            entropies = [
                calculate_branch_diversity(b).get('entropy', 0.0)
                for b in patterns.branches
            ]
            features['avgBranchEntropy'] = np.mean(entropies) if entropies else 0.0

        return features

    def _get_base_recommendation(self) -> str:
        """Get base solver recommendation using existing AUTO logic."""
        # Simplified recommendation logic
        if self._network is None:
            return 'MVA'

        # Check network characteristics
        has_product_form = self._has_product_form()
        has_single_chain = self._has_single_chain()
        has_multi_chain = self._has_multi_chain()
        total_jobs = self._get_total_jobs()
        avg_jobs_per_chain = self._get_avg_jobs_per_chain()

        if has_single_chain:
            return 'NC'  # Normalizing Constant for single chain
        elif has_multi_chain and has_product_form and total_jobs < 10:
            return 'NC'
        elif has_multi_chain and has_product_form and avg_jobs_per_chain < 30:
            return 'MVA'
        elif has_multi_chain and avg_jobs_per_chain > 30:
            return 'FLUID'
        else:
            return 'MVA'  # Default fallback

    def _has_product_form(self) -> bool:
        """Check if network has product form."""
        return True  # Simplified - would check actual network properties

    def _has_single_chain(self) -> bool:
        """Check if network has single chain."""
        if self._network is None:
            return False
        if hasattr(self._network, 'getNumberOfChains'):
            return self._network.getNumberOfChains() == 1
        return True

    def _has_multi_chain(self) -> bool:
        """Check if network has multiple chains."""
        if self._network is None:
            return False
        if hasattr(self._network, 'getNumberOfChains'):
            return self._network.getNumberOfChains() > 1
        return False

    def _get_total_jobs(self) -> int:
        """Get total number of jobs in network."""
        if self._network is None:
            return 0
        if hasattr(self._network, 'getTotalJobs'):
            return self._network.getTotalJobs()
        return 0

    def _get_avg_jobs_per_chain(self) -> float:
        """Get average jobs per chain."""
        total_jobs = self._get_total_jobs()
        num_chains = 1
        if self._network is not None and hasattr(self._network, 'getNumberOfChains'):
            num_chains = max(1, self._network.getNumberOfChains())
        return total_jobs / num_chains

    def _enhance_recommendation_with_workflow(
        self,
        base_recommendation: str,
        workflow_features: Dict[str, Any],
        analysis: WorkflowAnalysis
    ) -> ExtendedSolverRecommendation:
        """Enhance solver recommendation with workflow analysis insights."""
        reasoning: List[str] = []
        alternatives: List[str] = []
        final_recommendation = base_recommendation
        confidence = 0.7  # Base confidence

        # Analyze workflow patterns for solver selection
        has_sequences = workflow_features.get('hasSequencePatterns', False)
        has_parallels = workflow_features.get('hasParallelPatterns', False)
        has_loops = workflow_features.get('hasLoopPatterns', False)
        has_branches = workflow_features.get('hasBranchPatterns', False)

        # Sequence pattern analysis
        if has_sequences:
            reasoning.append("Detected sequence patterns - suitable for analytical methods")
            confidence += 0.1

            max_seq_length = workflow_features.get('maxSequenceLength', 0)
            if max_seq_length > 10:
                reasoning.append("Long sequences detected - consider FLUID approximation")
                if final_recommendation == 'MVA':
                    alternatives.append('FLUID')

        # Parallel pattern analysis
        if has_parallels:
            reasoning.append("Detected parallel patterns - fork-join structures present")

            max_parallelism = workflow_features.get('maxParallelism', 0)
            if max_parallelism > 5:
                reasoning.append(
                    "High parallelism detected - exact methods may be computationally expensive"
                )
                if final_recommendation in ['NC', 'MVA']:
                    final_recommendation = 'SSA'  # Switch to simulation
                    reasoning.append("Switching to SSA for high-parallelism workflow")
                alternatives.append('JMT')
            else:
                confidence += 0.05

        # Loop pattern analysis
        if has_loops:
            avg_loop_prob = workflow_features.get('avgLoopProbability', 0.0)
            max_loop_prob = workflow_features.get('maxLoopProbability', 0.0)

            reasoning.append(f"Detected loop patterns with avg probability {avg_loop_prob:.2f}")

            if max_loop_prob > 0.8:
                reasoning.append(
                    "High loop probability detected - may cause numerical instability"
                )
                confidence -= 0.1

                if final_recommendation in ['NC', 'MVA']:
                    alternatives.append('SSA')
                    alternatives.append('FLUID')
            elif max_loop_prob > 0.5:
                reasoning.append("Moderate loop probability - analytical methods suitable")
                confidence += 0.05

        # Branch pattern analysis
        if has_branches:
            avg_entropy = workflow_features.get('avgBranchEntropy', 0.0)
            max_branches = workflow_features.get('maxBranches', 0)

            reasoning.append(f"Detected branch patterns with avg entropy {avg_entropy:.2f}")

            if avg_entropy > 1.5:
                reasoning.append("High branching entropy - complex decision structure")
                if max_branches > 5:
                    reasoning.append("Many branches detected - consider simulation methods")
                    alternatives.append('JMT')
                    alternatives.append('SSA')

            if avg_entropy < 0.5:
                reasoning.append("Low branching entropy - deterministic-like behavior")
                confidence += 0.1

        # Complexity-based adjustments
        original_nodes = workflow_features.get('originalNodeCount', 0)
        optimized_nodes = workflow_features.get('optimizedNodeCount', 0)

        if original_nodes > optimized_nodes and original_nodes > 0:
            reduction = (original_nodes - optimized_nodes) / original_nodes
            reasoning.append(
                f"Workflow complexity reduced by {int(reduction * 100)}% through pattern optimization"
            )
            confidence += 0.1

        if optimized_nodes > 50:
            reasoning.append("Large optimized workflow - consider approximation methods")
            if final_recommendation == 'NC':
                final_recommendation = 'MVA'
                reasoning.append("Switching from NC to MVA for large workflow")
            alternatives.append('FLUID')

        # Ensure confidence is within bounds
        confidence = min(1.0, max(0.1, confidence))

        # Add standard alternatives if not already present
        standard_solvers = ['MVA', 'NC', 'SSA', 'FLUID', 'JMT', 'CTMC']
        for solver in standard_solvers:
            if solver != final_recommendation and solver not in alternatives:
                alternatives.append(solver)

        return ExtendedSolverRecommendation(
            recommended_solver=final_recommendation,
            confidence=confidence,
            reasoning=reasoning,
            workflow_features=workflow_features,
            alternative_solvers=alternatives[:3]  # Limit to top 3 alternatives
        )

    def get_optimization_insights(self) -> Dict[str, Any]:
        """
        Get workflow optimization insights for solver performance tuning.

        Returns:
            Dictionary with optimization insights
        """
        analysis = self._workflow_analyzer.analyze_workflow()
        insights: Dict[str, Any] = {}

        # Add optimization recommendations
        insights['recommendations'] = self._workflow_analyzer.get_optimization_recommendations(analysis)

        # Add pattern-specific insights
        insights['patternInsights'] = self._generate_pattern_insights(analysis.detected_patterns)

        # Add performance predictions
        insights['performancePredictions'] = self._generate_performance_predictions(analysis)

        return insights

    def _generate_pattern_insights(self, patterns: DetectedPatterns) -> Dict[str, Any]:
        """Generate insights specific to detected patterns."""
        insights: Dict[str, Any] = {}

        if patterns.sequences:
            insights['sequenceOptimization'] = (
                "Consider merging sequential services to reduce overhead"
            )

        if patterns.parallels:
            insights['parallelOptimization'] = (
                "Parallel patterns can benefit from resource pooling strategies"
            )

        if patterns.loops:
            insights['loopOptimization'] = (
                "High-probability loops may benefit from caching or memoization"
            )

        if patterns.branches:
            avg_branches = np.mean([len(b.branch_nodes) for b in patterns.branches])
            if avg_branches > 3:
                insights['branchOptimization'] = (
                    "Complex branching patterns - consider load balancing strategies"
                )

        return insights

    def _generate_performance_predictions(self, analysis: WorkflowAnalysis) -> Dict[str, Any]:
        """Generate performance predictions based on workflow analysis."""
        predictions: Dict[str, Any] = {}

        update_stats = analysis.statistics.get('updateStats', {})
        reduction_ratio = update_stats.get('reductionRatio', 0.0)

        if reduction_ratio > 0.1:
            predictions['complexityReduction'] = (
                f"Expected {int(reduction_ratio * 100)}% reduction in solve time"
            )

        patterns = analysis.detected_patterns
        if patterns.parallels:
            predictions['parallelizationPotential'] = (
                "High potential for parallel execution optimization"
            )

        if patterns.loops:
            predictions['convergenceConsiderations'] = (
                "Loop patterns may affect solver convergence rates"
            )

        return predictions

    def validate_workflow_enhancement(self) -> bool:
        """
        Validate that workflow analysis enhances solver selection.

        Returns:
            True if validation passes
        """
        try:
            recommendation = self.recommend_solver_with_workflow_analysis()
            analysis = self._workflow_analyzer.analyze_workflow()

            # Validate that recommendation is sensible
            return (
                recommendation.confidence > 0.1 and
                len(recommendation.reasoning) > 0 and
                self._workflow_analyzer.validate_analysis(analysis)
            )
        except Exception:
            return False
