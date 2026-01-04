/**
 * @file Workflow Automatic Integration
 * 
 * Provides automatic integration capabilities for workflow analysis with
 * LINE solver framework. Enables seamless workflow model construction and
 * analysis through automated pattern detection and model generation.
 * 
 * @since LINE 3.0
 */
package jline.api.wf

import jline.lang.Network
import jline.solvers.auto.SolverAUTO
import jline.solvers.auto.ModelAnalyzer
import jline.solvers.NetworkSolver
import jline.solvers.SolverOptions
import jline.solvers.fluid.SolverFluid
import java.util.*

/**
 * Integration class connecting workflow analysis to the AUTO solver.
 * 
 * This class extends the AUTO solver's capabilities by incorporating
 * workflow pattern analysis and optimization techniques from the MDN toolbox.
 * It provides intelligent solver selection based on detected workflow patterns.
 */
class Wf_auto_integration(
    private val network: Network,
    private val options: SolverOptions = SolverOptions()
) {
    
    private val workflowAnalyzer = Wf_analyzer(network)
    private val modelAnalyzer = ModelAnalyzer(network)
    
    /**
     * Data class for extended solver recommendations.
     */
    data class ExtendedSolverRecommendation(
        val recommendedSolver: String,
        val confidence: Double,
        val reasoning: List<String>,
        val workflowFeatures: Map<String, Any>,
        val alternativeSolvers: List<String>
    )
    
    /**
     * Perform workflow-aware solver selection.
     * 
     * @return Extended solver recommendation with workflow analysis
     */
    fun recommendSolverWithWorkflowAnalysis(): ExtendedSolverRecommendation {
        // Perform workflow analysis
        val workflowAnalysis = workflowAnalyzer.analyzeWorkflow()
        
        // Extract workflow features for solver selection
        val workflowFeatures = extractWorkflowFeatures(workflowAnalysis)
        
        // Get base recommendation from model analyzer
        val baseRecommendation = getBaseRecommendation()
        
        // Enhance recommendation with workflow insights
        val enhancedRecommendation = enhanceRecommendationWithWorkflow(
            baseRecommendation,
            workflowFeatures,
            workflowAnalysis
        )
        
        return enhancedRecommendation
    }
    
    /**
     * Extract workflow features for solver selection.
     */
    private fun extractWorkflowFeatures(
        analysis: Wf_analyzer.WorkflowAnalysis
    ): Map<String, Any> {
        val features = HashMap<String, Any>()
        
        val patterns = analysis.detectedPatterns
        val originalComplexity = analysis.statistics["originalComplexity"] as? Map<*, *>
        val optimizedComplexity = analysis.statistics["optimizedComplexity"] as? Map<*, *>
        
        // Pattern-based features
        features["hasSequencePatterns"] = patterns.sequences.isNotEmpty()
        features["hasParallelPatterns"] = patterns.parallels.isNotEmpty()
        features["hasLoopPatterns"] = patterns.loops.isNotEmpty()
        features["hasBranchPatterns"] = patterns.branches.isNotEmpty()
        
        features["numSequences"] = patterns.sequences.size
        features["numParallels"] = patterns.parallels.size
        features["numLoops"] = patterns.loops.size
        features["numBranches"] = patterns.branches.size
        
        // Complexity features
        features["originalNodeCount"] = originalComplexity?.get("totalNodes") ?: 0
        features["originalLinkCount"] = originalComplexity?.get("totalLinks") ?: 0
        features["optimizedNodeCount"] = optimizedComplexity?.get("totalNodes") ?: 0
        features["optimizedLinkCount"] = optimizedComplexity?.get("totalLinks") ?: 0
        
        // Pattern complexity metrics
        if (patterns.sequences.isNotEmpty()) {
            val avgSeqLength = patterns.sequences.map { it.size }.average()
            val maxSeqLength = patterns.sequences.maxByOrNull { it.size }?.size ?: 0
            features["avgSequenceLength"] = avgSeqLength
            features["maxSequenceLength"] = maxSeqLength
        }
        
        if (patterns.parallels.isNotEmpty()) {
            val avgParallelism = patterns.parallels.map { it.size }.average()
            val maxParallelism = patterns.parallels.maxByOrNull { it.size }?.size ?: 0
            features["avgParallelism"] = avgParallelism
            features["maxParallelism"] = maxParallelism
        }
        
        if (patterns.loops.isNotEmpty()) {
            val loopStats = analysis.statistics["loopStats"] as? Map<*, *>
            features["avgLoopProbability"] = loopStats?.get("avgLoopProbability") ?: 0.0
            features["maxLoopProbability"] = loopStats?.get("maxLoopProbability") ?: 0.0
        }
        
        if (patterns.branches.isNotEmpty()) {
            val branchStats = analysis.statistics["branchStats"] as? Map<*, *>
            features["avgBranches"] = branchStats?.get("avgBranches") ?: 0.0
            features["maxBranches"] = branchStats?.get("maxBranches") ?: 0
            
            // Calculate average entropy across all branches
            val avgEntropy = patterns.branches.map { 
                val diversity = Wf_branch_detector.calculateBranchDiversity(it)
                diversity["entropy"] as? Double ?: 0.0
            }.average()
            features["avgBranchEntropy"] = avgEntropy
        }
        
        return features
    }
    
    /**
     * Get base solver recommendation using existing AUTO logic.
     */
    private fun getBaseRecommendation(): String {
        // Use ModelAnalyzer to get base recommendation
        val hasProductForm = modelAnalyzer.hasProductForm()
        val hasSingleChain = modelAnalyzer.hasSingleChain()
        val hasMultiChain = modelAnalyzer.hasMultiChain()
        val totalJobs = modelAnalyzer.getTotalJobs()
        val avgJobsPerChain = modelAnalyzer.getAvgJobsPerChain()
        
        return when {
            hasSingleChain -> "NC"  // Normalizing Constant for single chain
            hasMultiChain && hasProductForm && totalJobs < 10 -> "NC"
            hasMultiChain && hasProductForm && avgJobsPerChain < 30 -> "MVA"
            hasMultiChain && avgJobsPerChain > 30 -> "FLUID"
            else -> "MVA"  // Default fallback
        }
    }
    
    /**
     * Enhance solver recommendation with workflow analysis insights.
     */
    private fun enhanceRecommendationWithWorkflow(
        baseRecommendation: String,
        workflowFeatures: Map<String, Any>,
        analysis: Wf_analyzer.WorkflowAnalysis
    ): ExtendedSolverRecommendation {
        val reasoning = ArrayList<String>()
        val alternatives = ArrayList<String>()
        var finalRecommendation = baseRecommendation
        var confidence = 0.7  // Base confidence
        
        // Analyze workflow patterns for solver selection
        val hasSequences = workflowFeatures["hasSequencePatterns"] as Boolean
        val hasParallels = workflowFeatures["hasParallelPatterns"] as Boolean
        val hasLoops = workflowFeatures["hasLoopPatterns"] as Boolean
        val hasBranches = workflowFeatures["hasBranchPatterns"] as Boolean
        
        // Sequence pattern analysis
        if (hasSequences) {
            reasoning.add("Detected sequence patterns - suitable for analytical methods")
            confidence += 0.1
            
            val maxSeqLength = workflowFeatures["maxSequenceLength"] as? Int ?: 0
            if (maxSeqLength > 10) {
                reasoning.add("Long sequences detected - consider FLUID approximation")
                if (finalRecommendation == "MVA") {
                    alternatives.add("FLUID")
                }
            }
        }
        
        // Parallel pattern analysis
        if (hasParallels) {
            reasoning.add("Detected parallel patterns - fork-join structures present")
            
            val maxParallelism = workflowFeatures["maxParallelism"] as? Int ?: 0
            if (maxParallelism > 5) {
                reasoning.add("High parallelism detected - exact methods may be computationally expensive")
                if (finalRecommendation in listOf("NC", "MVA")) {
                    finalRecommendation = "SSA"  // Switch to simulation
                    reasoning.add("Switching to SSA for high-parallelism workflow")
                }
                alternatives.add("JMT")
            } else {
                confidence += 0.05
            }
        }
        
        // Loop pattern analysis
        if (hasLoops) {
            val avgLoopProb = workflowFeatures["avgLoopProbability"] as? Double ?: 0.0
            val maxLoopProb = workflowFeatures["maxLoopProbability"] as? Double ?: 0.0
            
            reasoning.add("Detected loop patterns with avg probability $avgLoopProb")
            
            if (maxLoopProb > 0.8) {
                reasoning.add("High loop probability detected - may cause numerical instability")
                confidence -= 0.1
                
                if (finalRecommendation in listOf("NC", "MVA")) {
                    alternatives.add("SSA")
                    alternatives.add("FLUID")
                }
            } else if (maxLoopProb > 0.5) {
                reasoning.add("Moderate loop probability - analytical methods suitable")
                confidence += 0.05
            }
        }
        
        // Branch pattern analysis
        if (hasBranches) {
            val avgEntropy = workflowFeatures["avgBranchEntropy"] as? Double ?: 0.0
            val maxBranches = workflowFeatures["maxBranches"] as? Int ?: 0
            
            reasoning.add("Detected branch patterns with avg entropy $avgEntropy")
            
            if (avgEntropy > 1.5) {
                reasoning.add("High branching entropy - complex decision structure")
                if (maxBranches > 5) {
                    reasoning.add("Many branches detected - consider simulation methods")
                    alternatives.add("JMT")
                    alternatives.add("SSA")
                }
            }
            
            if (avgEntropy < 0.5) {
                reasoning.add("Low branching entropy - deterministic-like behavior")
                confidence += 0.1
            }
        }
        
        // Complexity-based adjustments
        val originalNodes = workflowFeatures["originalNodeCount"] as? Int ?: 0
        val optimizedNodes = workflowFeatures["optimizedNodeCount"] as? Int ?: 0
        
        if (originalNodes > optimizedNodes) {
            val reduction = (originalNodes - optimizedNodes).toDouble() / originalNodes
            reasoning.add("Workflow complexity reduced by ${(reduction * 100).toInt()}% through pattern optimization")
            confidence += 0.1
        }
        
        if (optimizedNodes > 50) {
            reasoning.add("Large optimized workflow - consider approximation methods")
            if (finalRecommendation == "NC") {
                finalRecommendation = "MVA"
                reasoning.add("Switching from NC to MVA for large workflow")
            }
            alternatives.add("FLUID")
        }
        
        // Ensure confidence is within bounds
        confidence = minOf(1.0, maxOf(0.1, confidence))
        
        // Add standard alternatives if not already present
        val standardSolvers = listOf("MVA", "NC", "SSA", "FLUID", "JMT", "CTMC")
        for (solver in standardSolvers) {
            if (solver != finalRecommendation && solver !in alternatives) {
                alternatives.add(solver)
            }
        }
        
        return ExtendedSolverRecommendation(
            recommendedSolver = finalRecommendation,
            confidence = confidence,
            reasoning = reasoning,
            workflowFeatures = workflowFeatures,
            alternativeSolvers = alternatives.take(3)  // Limit to top 3 alternatives
        )
    }
    
    /**
     * Create and configure solver based on workflow analysis.
     */
    fun createOptimalSolver(): NetworkSolver {
        val recommendation = recommendSolverWithWorkflowAnalysis()
        
        // Create solver based on recommendation
        return when (recommendation.recommendedSolver.uppercase()) {
            "MVA" -> jline.solvers.mva.SolverMVA(network, options)
            "NC" -> jline.solvers.nc.SolverNC(network, options)
            "SSA" -> jline.solvers.ssa.SolverSSA(network, options)
            "FLUID" -> jline.solvers.fluid.SolverFluid(network, options)
            "JMT" -> jline.solvers.jmt.SolverJMT(network, options)
            "CTMC" -> jline.solvers.ctmc.SolverCTMC(network, options)
            else -> SolverAUTO(network, options)  // Fallback to AUTO
        }
    }
    
    /**
     * Get workflow optimization insights for solver performance tuning.
     */
    fun getOptimizationInsights(): Map<String, Any> {
        val analysis = workflowAnalyzer.analyzeWorkflow()
        val insights = HashMap<String, Any>()
        
        // Add optimization recommendations
        insights["recommendations"] = workflowAnalyzer.getOptimizationRecommendations(analysis)
        
        // Add pattern-specific insights
        insights["patternInsights"] = generatePatternInsights(analysis.detectedPatterns)
        
        // Add performance predictions
        insights["performancePredictions"] = generatePerformancePredictions(analysis)
        
        return insights
    }
    
    /**
     * Generate insights specific to detected patterns.
     */
    private fun generatePatternInsights(patterns: Wf_analyzer.DetectedPatterns): Map<String, Any> {
        val insights = HashMap<String, Any>()
        
        if (patterns.sequences.isNotEmpty()) {
            insights["sequenceOptimization"] = "Consider merging sequential services to reduce overhead"
        }
        
        if (patterns.parallels.isNotEmpty()) {
            insights["parallelOptimization"] = "Parallel patterns can benefit from resource pooling strategies"
        }
        
        if (patterns.loops.isNotEmpty()) {
            val highProbLoops = patterns.loops.size  // Simplified - would check actual probabilities
            if (highProbLoops > 0) {
                insights["loopOptimization"] = "High-probability loops may benefit from caching or memoization"
            }
        }
        
        if (patterns.branches.isNotEmpty()) {
            val avgBranches = patterns.branches.map { it.branchNodes.size }.average()
            if (avgBranches > 3) {
                insights["branchOptimization"] = "Complex branching patterns - consider load balancing strategies"
            }
        }
        
        return insights
    }
    
    /**
     * Generate performance predictions based on workflow analysis.
     */
    private fun generatePerformancePredictions(analysis: Wf_analyzer.WorkflowAnalysis): Map<String, Any> {
        val predictions = HashMap<String, Any>()
        
        val updateStats = analysis.statistics["updateStats"] as? Map<*, *>
        val reductionRatio = updateStats?.get("reductionRatio") as? Double ?: 0.0
        
        if (reductionRatio > 0.1) {
            predictions["complexityReduction"] = "Expected ${(reductionRatio * 100).toInt()}% reduction in solve time"
        }
        
        val patterns = analysis.detectedPatterns
        if (patterns.parallels.isNotEmpty()) {
            predictions["parallelizationPotential"] = "High potential for parallel execution optimization"
        }
        
        if (patterns.loops.isNotEmpty()) {
            predictions["convergenceConsiderations"] = "Loop patterns may affect solver convergence rates"
        }
        
        return predictions
    }
    
    /**
     * Validate that workflow analysis enhances solver selection.
     */
    fun validateWorkflowEnhancement(): Boolean {
        try {
            val recommendation = recommendSolverWithWorkflowAnalysis()
            val analysis = workflowAnalyzer.analyzeWorkflow()
            
            // Validate that recommendation is sensible
            return recommendation.confidence > 0.1 && 
                   recommendation.reasoning.isNotEmpty() &&
                   workflowAnalyzer.validateAnalysis(analysis)
        } catch (e: Exception) {
            return false
        }
    }
}