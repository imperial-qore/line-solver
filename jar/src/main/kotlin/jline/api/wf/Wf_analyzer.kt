/**
 * @file Workflow Analysis Engine
 * 
 * Provides comprehensive workflow analysis capabilities including pattern
 * recognition, performance analysis, and workflow model construction.
 * Integrates various detection algorithms for complete workflow characterization.
 * 
 * @since LINE 3.0
 */
package jline.api.wf

import jline.GlobalConstants
import jline.util.matrix.Matrix
import jline.lang.Network
import jline.lang.nodes.*
import java.util.*
import kotlin.collections.ArrayList

/**
 * Main workflow analyzer that coordinates pattern detection and optimization.
 * 
 * This class serves as the primary interface for workflow analysis, combining
 * all pattern detectors and the pattern updater to provide comprehensive
 * workflow optimization capabilities.
 * 
 * Based on the AUTO solver functionality from the MDN toolbox.
 */
class Wf_analyzer(private val network: Network) {

    /**
     * Data class to represent comprehensive workflow analysis results.
     */
    data class WorkflowAnalysis(
        val originalWorkflow: WorkflowRepresentation,
        val detectedPatterns: DetectedPatterns,
        val optimizedWorkflow: Wf_pattern_updater.UpdatedWorkflow,
        val statistics: Map<String, Any>
    )
    
    /**
     * Data class to represent workflow in matrix form.
     */
    data class WorkflowRepresentation(
        val linkMatrix: Matrix,
        val serviceNodes: List<Int>,
        val forkNodes: List<Int>,
        val joinNodes: List<Int>,
        val routerNodes: List<Int>,
        val serviceParameters: Map<Int, Wf_pattern_updater.ServiceParameters>
    )
    
    /**
     * Data class to hold all detected patterns.
     */
    data class DetectedPatterns(
        val sequences: List<List<Int>>,
        val parallels: List<List<Int>>,
        val loops: List<Int>,
        val branches: List<Wf_branch_detector.BranchPattern>
    )
    
    /**
     * Perform comprehensive workflow analysis and optimization.
     * 
     * @return Complete workflow analysis results
     */
    fun analyzeWorkflow(): WorkflowAnalysis {
        // Step 1: Convert network to workflow representation
        val workflowRep = convertNetworkToWorkflow()
        
        // Step 2: Detect all patterns
        val patterns = detectAllPatterns(workflowRep)
        
        // Step 3: Optimize workflow by updating patterns
        val optimizedWorkflow = Wf_pattern_updater.updatePatterns(
            workflowRep.linkMatrix,
            workflowRep.serviceNodes,
            workflowRep.forkNodes,
            workflowRep.joinNodes,
            workflowRep.routerNodes,
            workflowRep.serviceParameters
        )
        
        // Step 4: Generate statistics
        val statistics = generateAnalysisStatistics(workflowRep, patterns, optimizedWorkflow)
        
        return WorkflowAnalysis(workflowRep, patterns, optimizedWorkflow, statistics)
    }
    
    /**
     * Convert LINE Network to workflow matrix representation.
     */
    private fun convertNetworkToWorkflow(): WorkflowRepresentation {
        val nodes = network.nodes
        val linkMatrix = buildLinkMatrix()
        val nodeClassification = classifyNodes(nodes)
        val serviceParams = extractServiceParameters(nodeClassification.serviceNodes)
        
        return WorkflowRepresentation(
            linkMatrix = linkMatrix,
            serviceNodes = nodeClassification.serviceNodes,
            forkNodes = nodeClassification.forkNodes,
            joinNodes = nodeClassification.joinNodes,
            routerNodes = nodeClassification.routerNodes,
            serviceParameters = serviceParams
        )
    }
    
    /**
     * Build link matrix from network routing matrix.
     */
    private fun buildLinkMatrix(): Matrix {
        val routingMatrix = network.getLinkedRoutingMatrix()
        val nodes = network.nodes
        val links = ArrayList<Triple<Int, Int, Double>>()
        
        // Extract non-zero transitions from the nested routing matrix structure
        if (routingMatrix.isNotEmpty()) {
            for ((fromClass, classRouting) in routingMatrix) {
                for ((toClass, matrix) in classRouting) {
                    // Iterate through the matrix to find non-zero transitions
                    for (i in 0 until matrix.getNumRows()) {
                        for (j in 0 until matrix.getNumCols()) {
                            val prob = matrix.get(i, j)
                            if (prob > GlobalConstants.Zero) {
                                // Use node indices as IDs
                                links.add(Triple(i, j, prob))
                            }
                        }
                    }
                }
            }
        }
        
        // If no routing found, create basic sequential connections
        if (links.isEmpty() && nodes.size > 1) {
            for (i in 0 until nodes.size - 1) {
                links.add(Triple(i, i + 1, 1.0))
            }
        }
        
        // Create matrix representation
        val linkMatrix = Matrix.zeros(maxOf(1, links.size), 3)
        for ((index, link) in links.withIndex()) {
            linkMatrix.set(index, 0, link.first.toDouble())
            linkMatrix.set(index, 1, link.second.toDouble())
            linkMatrix.set(index, 2, link.third)
        }
        
        return linkMatrix
    }
    
    /**
     * Data class for node classification results.
     */
    data class NodeClassification(
        val serviceNodes: List<Int>,
        val forkNodes: List<Int>,
        val joinNodes: List<Int>,
        val routerNodes: List<Int>
    )
    
    /**
     * Classify nodes by their type and function.
     */
    private fun classifyNodes(nodes: List<Node>): NodeClassification {
        val serviceNodes = ArrayList<Int>()
        val forkNodes = ArrayList<Int>()
        val joinNodes = ArrayList<Int>()
        val routerNodes = ArrayList<Int>()
        
        for ((index, node) in nodes.withIndex()) {
            when (node) {
                is jline.lang.nodes.Queue -> serviceNodes.add(index)
                is jline.lang.nodes.Delay -> serviceNodes.add(index)
                is jline.lang.nodes.Fork -> forkNodes.add(index)
                is jline.lang.nodes.Join -> joinNodes.add(index)
                is jline.lang.nodes.Router -> routerNodes.add(index)
                is jline.lang.nodes.ClassSwitch -> routerNodes.add(index)
                // Source and Sink are typically not modified in pattern analysis
            }
        }
        
        return NodeClassification(serviceNodes, forkNodes, joinNodes, routerNodes)
    }
    
    /**
     * Extract service parameters from service nodes.
     */
    private fun extractServiceParameters(
        serviceNodes: List<Int>
    ): Map<Int, Wf_pattern_updater.ServiceParameters> {
        val params = HashMap<Int, Wf_pattern_updater.ServiceParameters>()
        val nodes = network.nodes
        
        for (nodeIndex in serviceNodes) {
            if (nodeIndex < nodes.size) {
                val node = nodes[nodeIndex]
                
                // Extract phase-type parameters from service distributions
                val serviceParam = when (node) {
                    is jline.lang.nodes.Queue -> extractPHFromQueue(node)
                    is jline.lang.nodes.Delay -> extractPHFromDelay(node)
                    else -> createDefaultPH()
                }
                
                params[nodeIndex] = serviceParam
            }
        }
        
        return params
    }
    
    /**
     * Extract phase-type parameters from a Queue node.
     */
    private fun extractPHFromQueue(queue: jline.lang.nodes.Queue): Wf_pattern_updater.ServiceParameters {
        // Placeholder implementation - would extract actual PH parameters
        // from the queue's service distribution
        return createDefaultPH()
    }
    
    /**
     * Extract phase-type parameters from a Delay node.
     */
    private fun extractPHFromDelay(delay: jline.lang.nodes.Delay): Wf_pattern_updater.ServiceParameters {
        // Placeholder implementation - would extract actual PH parameters
        // from the delay's service distribution
        return createDefaultPH()
    }
    
    /**
     * Create default phase-type parameters (exponential distribution).
     */
    private fun createDefaultPH(): Wf_pattern_updater.ServiceParameters {
        val alpha = Matrix.ones(1, 1)
        val T = Matrix.ones(1, 1)
        T.set(0, 0, -1.0) // -1 rate matrix for exponential
        return Wf_pattern_updater.ServiceParameters(alpha, T)
    }
    
    /**
     * Detect all workflow patterns.
     */
    private fun detectAllPatterns(workflow: WorkflowRepresentation): DetectedPatterns {
        val sequences = Wf_sequence_detector.detectSequences(
            workflow.linkMatrix,
            workflow.serviceNodes
        )

        val parallels = Wf_parallel_detector.detectParallel(
            workflow.linkMatrix,
            workflow.serviceNodes,
            workflow.forkNodes,
            workflow.joinNodes
        )

        val loops = Wf_loop_detector.detectLoops(
            workflow.linkMatrix,
            workflow.serviceNodes,
            workflow.routerNodes,
            workflow.joinNodes
        )

        val branches = Wf_branch_detector.detectBranches(
            workflow.linkMatrix,
            workflow.serviceNodes,
            workflow.joinNodes
        )
        
        return DetectedPatterns(sequences, parallels, loops, branches)
    }
    
    /**
     * Generate comprehensive analysis statistics.
     */
    private fun generateAnalysisStatistics(
        originalWorkflow: WorkflowRepresentation,
        patterns: DetectedPatterns,
        optimizedWorkflow: Wf_pattern_updater.UpdatedWorkflow
    ): Map<String, Any> {
        val stats = HashMap<String, Any>()
        
        // Pattern detection statistics
        stats["sequenceStats"] = Wf_sequence_detector.getSequenceStats(patterns.sequences)
        stats["parallelStats"] = Wf_parallel_detector.getParallelStats(patterns.parallels)
        stats["loopStats"] = Wf_loop_detector.getLoopStats(
            patterns.loops,
            originalWorkflow.linkMatrix,
            originalWorkflow.routerNodes
        )
        stats["branchStats"] = Wf_branch_detector.getBranchStats(patterns.branches)

        // Optimization statistics
        stats["updateStats"] = Wf_pattern_updater.getUpdateStats(
            originalWorkflow.linkMatrix,
            optimizedWorkflow
        )
        
        // Overall complexity metrics
        stats["originalComplexity"] = calculateWorkflowComplexity(originalWorkflow)
        stats["optimizedComplexity"] = calculateOptimizedComplexity(optimizedWorkflow)
        
        return stats
    }
    
    /**
     * Calculate workflow complexity metrics.
     */
    private fun calculateWorkflowComplexity(workflow: WorkflowRepresentation): Map<String, Any> {
        val complexity = HashMap<String, Any>()
        
        complexity["totalNodes"] = workflow.serviceNodes.size + workflow.forkNodes.size + 
                                  workflow.joinNodes.size + workflow.routerNodes.size
        complexity["totalLinks"] = workflow.linkMatrix.getNumRows()
        complexity["serviceNodes"] = workflow.serviceNodes.size
        complexity["controlNodes"] = workflow.forkNodes.size + workflow.joinNodes.size + workflow.routerNodes.size
        
        // Calculate connectivity metrics
        val nodeSet = HashSet<Int>()
        for (i in 0 until workflow.linkMatrix.getNumRows()) {
            nodeSet.add(workflow.linkMatrix.get(i, 0).toInt())
            nodeSet.add(workflow.linkMatrix.get(i, 1).toInt())
        }
        complexity["connectedNodes"] = nodeSet.size
        
        // Calculate average degree
        val degrees = HashMap<Int, Int>()
        for (i in 0 until workflow.linkMatrix.getNumRows()) {
            val start = workflow.linkMatrix.get(i, 0).toInt()
            val end = workflow.linkMatrix.get(i, 1).toInt()
            degrees[start] = degrees.getOrDefault(start, 0) + 1
            degrees[end] = degrees.getOrDefault(end, 0) + 1
        }
        complexity["avgDegree"] = if (degrees.isNotEmpty()) degrees.values.average() else 0.0
        
        return complexity
    }
    
    /**
     * Calculate optimized workflow complexity.
     */
    private fun calculateOptimizedComplexity(
        optimizedWorkflow: Wf_pattern_updater.UpdatedWorkflow
    ): Map<String, Any> {
        val complexity = HashMap<String, Any>()
        
        complexity["totalNodes"] = optimizedWorkflow.serviceParameters.size
        complexity["totalLinks"] = optimizedWorkflow.linkMatrix.getNumRows()
        
        // Calculate connectivity for optimized workflow
        val nodeSet = HashSet<Int>()
        for (i in 0 until optimizedWorkflow.linkMatrix.getNumRows()) {
            nodeSet.add(optimizedWorkflow.linkMatrix.get(i, 0).toInt())
            nodeSet.add(optimizedWorkflow.linkMatrix.get(i, 1).toInt())
        }
        complexity["connectedNodes"] = nodeSet.size
        
        return complexity
    }
    
    /**
     * Get workflow optimization recommendations.
     */
    fun getOptimizationRecommendations(analysis: WorkflowAnalysis): List<String> {
        val recommendations = ArrayList<String>()
        
        val patterns = analysis.detectedPatterns
        val updateStats = analysis.statistics["updateStats"] as? Map<*, *>
        
        // Sequence recommendations
        if (patterns.sequences.isNotEmpty()) {
            recommendations.add("Found ${patterns.sequences.size} sequence patterns that can be simplified")
        }
        
        // Parallel recommendations
        if (patterns.parallels.isNotEmpty()) {
            recommendations.add("Found ${patterns.parallels.size} parallel patterns for potential optimization")
        }
        
        // Loop recommendations
        if (patterns.loops.isNotEmpty()) {
            recommendations.add("Found ${patterns.loops.size} loop patterns - consider loop unrolling for performance")
        }
        
        // Branch recommendations
        if (patterns.branches.isNotEmpty()) {
            recommendations.add("Found ${patterns.branches.size} branch patterns - analyze probability distributions")
            
            val highEntropyBranches = patterns.branches.filter {
                val diversity = Wf_branch_detector.calculateBranchDiversity(it)
                (diversity["entropy"] ?: 0.0) as Double > 1.0
            }
            
            if (highEntropyBranches.isNotEmpty()) {
                recommendations.add("${highEntropyBranches.size} branches have high entropy - consider load balancing")
            }
        }
        
        // Optimization impact
        val reductionRatio = updateStats?.get("reductionRatio") as? Double ?: 0.0
        if (reductionRatio > 0.1) {
            recommendations.add("Workflow complexity reduced by ${(reductionRatio * 100).toInt()}% through pattern optimization")
        }
        
        return recommendations
    }
    
    /**
     * Validate workflow analysis results.
     */
    fun validateAnalysis(analysis: WorkflowAnalysis): Boolean {
        // Validate optimized workflow structure
        if (!Wf_pattern_updater.validateUpdatedWorkflow(analysis.optimizedWorkflow)) {
            return false
        }

        // Validate pattern consistency
        for (sequence in analysis.detectedPatterns.sequences) {
            if (!Wf_sequence_detector.validateSequence(sequence, analysis.originalWorkflow.linkMatrix)) {
                return false
            }
        }

        for (parallel in analysis.detectedPatterns.parallels) {
            if (!Wf_parallel_detector.validateParallelPattern(
                parallel,
                analysis.originalWorkflow.linkMatrix,
                analysis.originalWorkflow.forkNodes,
                analysis.originalWorkflow.joinNodes
            )) {
                return false
            }
        }

        for (loop in analysis.detectedPatterns.loops) {
            if (!Wf_loop_detector.validateLoopPattern(
                loop,
                analysis.originalWorkflow.linkMatrix,
                analysis.originalWorkflow.routerNodes
            )) {
                return false
            }
        }

        for (branch in analysis.detectedPatterns.branches) {
            if (!Wf_branch_detector.validateBranchPattern(branch, analysis.originalWorkflow.linkMatrix)) {
                return false
            }
        }

        return true
    }
}