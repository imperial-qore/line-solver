/**
 * @file Workflow Pattern Updater
 * 
 * Implements dynamic pattern updating algorithms for workflow analysis.
 * Provides methods for refining and updating discovered workflow patterns
 * based on new observations and evolving process characteristics.
 * 
 * @since LINE 3.0
 */
package jline.api.wf

import jline.GlobalConstants
import jline.util.matrix.Matrix
import jline.lang.processes.APH
import java.util.*
import kotlin.collections.ArrayList

/**
 * Workflow pattern updater for simplifying and optimizing workflow networks.
 * 
 * This class is responsible for updating workflow networks after recognizing
 * basic workflow patterns (sequence, parallel, loop, branch) by combining
 * pattern detection with phase-type distribution convolution.
 * 
 * Based on AUTO_Pattern_Update.m from the MDN toolbox.
 */
object Wf_pattern_updater {
    
    /**
     * Data class to represent service parameters (phase-type distributions).
     */
    data class ServiceParameters(
        val alpha: Matrix,  // Initial probability vector
        val T: Matrix       // Transition rate matrix
    )
    
    /**
     * Data class to represent the updated workflow.
     */
    data class UpdatedWorkflow(
        val linkMatrix: Matrix,
        val serviceParameters: Map<Int, ServiceParameters>
    )
    
    /**
     * Update workflow network by simplifying detected patterns.
     * 
     * @param linkMatrix Original workflow link matrix
     * @param serviceNodes List of service node IDs
     * @param forkNodes List of fork node IDs
     * @param joinNodes List of join node IDs
     * @param routerNodes List of router node IDs
     * @param serviceParams Original service parameters (PH distributions)
     * @return Updated workflow with simplified patterns
     */
    @JvmStatic
    fun updatePatterns(
        linkMatrix: Matrix,
        serviceNodes: List<Int>,
        forkNodes: List<Int>,
        joinNodes: List<Int>,
        routerNodes: List<Int>,
        serviceParams: Map<Int, ServiceParameters>
    ): UpdatedWorkflow {
        // Make mutable copies
        var currentMatrix = linkMatrix.copy()
        val currentParams = HashMap(serviceParams)
        
        // Step 1: Simplify sequence patterns
        val sequences = Wf_sequence_detector.detectSequences(currentMatrix, serviceNodes)
        currentMatrix = updateSequencePatterns(currentMatrix, sequences, serviceNodes, currentParams)
        
        // Step 2: Simplify parallel patterns  
        val parallels = Wf_parallel_detector.detectParallel(currentMatrix, serviceNodes, forkNodes, joinNodes)
        currentMatrix = updateParallelPatterns(currentMatrix, parallels, serviceNodes, currentParams)
        
        // Step 3: Simplify loop patterns
        val loops = Wf_loop_detector.detectLoops(currentMatrix, serviceNodes, routerNodes, joinNodes)
        currentMatrix = updateLoopPatterns(currentMatrix, loops, serviceNodes, routerNodes, currentParams)
        
        // Step 4: Simplify branch patterns
        val branches = Wf_branch_detector.detectBranches(currentMatrix, serviceNodes, joinNodes)
        currentMatrix = updateBranchPatterns(currentMatrix, branches, serviceNodes, currentParams)
        
        return UpdatedWorkflow(currentMatrix, currentParams)
    }
    
    /**
     * Update workflow by simplifying sequence patterns.
     */
    private fun updateSequencePatterns(
        linkMatrix: Matrix,
        sequences: List<List<Int>>,
        serviceNodes: List<Int>,
        serviceParams: MutableMap<Int, ServiceParameters>
    ): Matrix {
        if (sequences.isEmpty()) return linkMatrix
        
        val serviceSet = serviceNodes.toSet()
        val updatedMatrix = linkMatrix.copy()
        val rowsToRemove = ArrayList<Int>()
        
        // Find rows connecting service nodes in sequences
        for (i in 0 until updatedMatrix.getNumRows()) {
            val start = updatedMatrix.get(i, 0).toInt()
            val end = updatedMatrix.get(i, 1).toInt()
            
            if (start in serviceSet && end in serviceSet) {
                // Check if this connection is part of a sequence
                for (sequence in sequences) {
                    for (j in 0 until sequence.size - 1) {
                        if (sequence[j] == start && sequence[j + 1] == end) {
                            rowsToRemove.add(i)
                            break
                        }
                    }
                }
            }
        }
        
        // Remove sequence connection rows
        val finalMatrix = removeMatrixRows(updatedMatrix, rowsToRemove)
        
        // Update node references and parameters
        for (sequence in sequences) {
            if (sequence.size >= 2) {
                val firstNode = sequence[0]
                val lastNode = sequence[sequence.size - 1]
                
                // Replace all occurrences of lastNode with firstNode
                replaceNodeReferences(finalMatrix, lastNode, firstNode)
                
                // Convolve service parameters for sequence
                val sequenceParams = sequence.map { serviceParams[it]!! }
                val convolvedParams = convolveSequence(sequenceParams)
                serviceParams[firstNode] = convolvedParams
                
                // Remove parameters for other nodes in sequence
                for (i in 1 until sequence.size) {
                    serviceParams.remove(sequence[i])
                }
            }
        }
        
        return finalMatrix
    }
    
    /**
     * Update workflow by simplifying parallel patterns.
     */
    private fun updateParallelPatterns(
        linkMatrix: Matrix,
        parallels: List<List<Int>>,
        serviceNodes: List<Int>,
        serviceParams: MutableMap<Int, ServiceParameters>
    ): Matrix {
        if (parallels.isEmpty()) return linkMatrix
        
        val updatedMatrix = linkMatrix.copy()
        
        for (parallel in parallels) {
            if (parallel.size >= 2) {
                val firstNode = parallel[0]
                
                // Find fork and join nodes for this parallel pattern
                val (forkNode, joinNode) = findForkJoinForParallel(updatedMatrix, parallel)
                
                if (forkNode != null && joinNode != null) {
                    // Remove all connections involving parallel nodes
                    val rowsToRemove = findRowsInvolvingNodes(updatedMatrix, parallel)
                    removeMatrixRows(updatedMatrix, rowsToRemove)
                    
                    // Replace fork and join references with first parallel node
                    replaceNodeReferences(updatedMatrix, forkNode, firstNode)
                    replaceNodeReferences(updatedMatrix, joinNode, firstNode)
                    
                    // Convolve service parameters for parallel execution
                    val parallelParams = parallel.map { serviceParams[it]!! }
                    val convolvedParams = convolveParallel(parallelParams)
                    serviceParams[firstNode] = convolvedParams
                    
                    // Remove parameters for other nodes in parallel
                    for (i in 1 until parallel.size) {
                        serviceParams.remove(parallel[i])
                    }
                }
            }
        }
        
        return updatedMatrix
    }
    
    /**
     * Update workflow by simplifying loop patterns.
     */
    private fun updateLoopPatterns(
        linkMatrix: Matrix,
        loops: List<Int>,
        serviceNodes: List<Int>,
        routerNodes: List<Int>,
        serviceParams: MutableMap<Int, ServiceParameters>
    ): Matrix {
        if (loops.isEmpty()) return linkMatrix
        
        val updatedMatrix = linkMatrix.copy()
        val routerSet = routerNodes.toSet()
        
        for (loopNode in loops) {
            val loopProb = Wf_loop_detector.getLoopProbability(loopNode, updatedMatrix, routerNodes)
            
            if (loopProb > GlobalConstants.Zero) {
                // Find router nodes involved in the loop
                val involvedRouters = findRoutersForLoop(updatedMatrix, loopNode, routerSet)
                
                // Remove loop connections
                val rowsToRemove = findLoopConnections(updatedMatrix, loopNode, involvedRouters)
                removeMatrixRows(updatedMatrix, rowsToRemove)
                
                // Replace router references with loop node
                for (router in involvedRouters) {
                    replaceNodeReferences(updatedMatrix, router, loopNode)
                }
                
                // Update transition probabilities
                updateLoopTransitionProbabilities(updatedMatrix, loopNode)
                
                // Convolve service parameters with loop probability
                val originalParams = serviceParams[loopNode]!!
                val loopedParams = convolveLoop(originalParams, loopProb)
                serviceParams[loopNode] = loopedParams
            }
        }
        
        return updatedMatrix
    }
    
    /**
     * Update workflow by simplifying branch patterns.
     */
    private fun updateBranchPatterns(
        linkMatrix: Matrix,
        branches: List<Wf_branch_detector.BranchPattern>,
        serviceNodes: List<Int>,
        serviceParams: MutableMap<Int, ServiceParameters>
    ): Matrix {
        if (branches.isEmpty()) return linkMatrix
        
        val updatedMatrix = linkMatrix.copy()
        
        for (branch in branches) {
            if (branch.branchNodes.size >= 2 && branch.forkNode != null) {
                val firstNode = branch.branchNodes[0]
                
                // Remove branch connections
                val rowsToRemove = findRowsInvolvingNodes(updatedMatrix, branch.branchNodes)
                removeMatrixRows(updatedMatrix, rowsToRemove)
                
                // Replace fork and join references
                replaceNodeReferences(updatedMatrix, branch.forkNode, firstNode)
                if (branch.joinNode != null) {
                    replaceNodeReferences(updatedMatrix, branch.joinNode, firstNode)
                }
                
                // Convolve service parameters for branches
                val branchParams = branch.branchNodes.map { serviceParams[it]!! }
                val convolvedParams = convolveBranches(branchParams, branch.probabilities)
                serviceParams[firstNode] = convolvedParams
                
                // Remove parameters for other branch nodes
                for (i in 1 until branch.branchNodes.size) {
                    serviceParams.remove(branch.branchNodes[i])
                }
            }
        }
        
        return updatedMatrix
    }
    
    // Helper methods for matrix operations
    
    private fun removeMatrixRows(matrix: Matrix, rowsToRemove: List<Int>): Matrix {
        val sortedRows = rowsToRemove.sorted().reversed()
        val result = matrix.copy()
        
        for (row in sortedRows) {
            if (row >= 0 && row < result.getNumRows()) {
                // Create new matrix without the specified row
                val newRows = result.getNumRows() - 1
                val newMatrix = Matrix.zeros(newRows, result.getNumCols())
                
                var newRowIndex = 0
                for (i in 0 until result.getNumRows()) {
                    if (i != row) {
                        for (j in 0 until result.getNumCols()) {
                            newMatrix.set(newRowIndex, j, result.get(i, j))
                        }
                        newRowIndex++
                    }
                }
                return newMatrix
            }
        }
        
        return result
    }
    
    private fun replaceNodeReferences(matrix: Matrix, oldNode: Int, newNode: Int) {
        for (i in 0 until matrix.getNumRows()) {
            if (matrix.get(i, 0).toInt() == oldNode) {
                matrix.set(i, 0, newNode.toDouble())
            }
            if (matrix.get(i, 1).toInt() == oldNode) {
                matrix.set(i, 1, newNode.toDouble())
            }
        }
    }
    
    private fun findRowsInvolvingNodes(matrix: Matrix, nodes: List<Int>): List<Int> {
        val nodeSet = nodes.toSet()
        val rowsToRemove = ArrayList<Int>()
        
        for (i in 0 until matrix.getNumRows()) {
            val start = matrix.get(i, 0).toInt()
            val end = matrix.get(i, 1).toInt()
            
            if (start in nodeSet || end in nodeSet) {
                rowsToRemove.add(i)
            }
        }
        
        return rowsToRemove
    }
    
    private fun findForkJoinForParallel(matrix: Matrix, parallel: List<Int>): Pair<Int?, Int?> {
        // Simplified implementation - would need more sophisticated logic
        return Pair(null, null)
    }
    
    private fun findRoutersForLoop(matrix: Matrix, loopNode: Int, routerSet: Set<Int>): List<Int> {
        val routers = ArrayList<Int>()
        
        for (i in 0 until matrix.getNumRows()) {
            val start = matrix.get(i, 0).toInt()
            val end = matrix.get(i, 1).toInt()
            
            if ((start == loopNode && end in routerSet) || 
                (end == loopNode && start in routerSet)) {
                if (start in routerSet) routers.add(start)
                if (end in routerSet) routers.add(end)
            }
        }
        
        return routers.distinct()
    }
    
    private fun findLoopConnections(matrix: Matrix, loopNode: Int, routers: List<Int>): List<Int> {
        val connections = ArrayList<Int>()
        val routerSet = routers.toSet()
        
        for (i in 0 until matrix.getNumRows()) {
            val start = matrix.get(i, 0).toInt()
            val end = matrix.get(i, 1).toInt()
            
            if ((start == loopNode && end in routerSet) || 
                (end == loopNode && start in routerSet) ||
                (start in routerSet && end == loopNode)) {
                connections.add(i)
            }
        }
        
        return connections
    }
    
    private fun updateLoopTransitionProbabilities(matrix: Matrix, loopNode: Int) {
        for (i in 0 until matrix.getNumRows()) {
            val start = matrix.get(i, 0).toInt()
            if (start == loopNode) {
                matrix.set(i, 2, 1.0) // Set probability to 1.0 after loop simplification
            }
        }
    }
    
    // Phase-type distribution convolution methods
    // These would typically call into specialized PH libraries
    
    private fun convolveSequence(params: List<ServiceParameters>): ServiceParameters {
        // Placeholder - would implement APH convolution for sequence
        return if (params.isNotEmpty()) {
            params[0] // Simplified - return first for now
        } else {
            ServiceParameters(Matrix.ones(1, 1), Matrix.zeros(1, 1))
        }
    }
    
    private fun convolveParallel(params: List<ServiceParameters>): ServiceParameters {
        // Placeholder - would implement APH convolution for parallel execution
        return if (params.isNotEmpty()) {
            params[0] // Simplified - return first for now
        } else {
            ServiceParameters(Matrix.ones(1, 1), Matrix.zeros(1, 1))
        }
    }
    
    private fun convolveLoop(params: ServiceParameters, loopProb: Double): ServiceParameters {
        // Placeholder - would implement APH convolution for loops
        return params // Simplified - return original for now
    }
    
    private fun convolveBranches(params: List<ServiceParameters>, probs: List<Double>): ServiceParameters {
        // Placeholder - would implement APH convolution for branches
        return if (params.isNotEmpty()) {
            params[0] // Simplified - return first for now
        } else {
            ServiceParameters(Matrix.ones(1, 1), Matrix.zeros(1, 1))
        }
    }
    
    /**
     * Validate the updated workflow structure.
     */
    @JvmStatic
    fun validateUpdatedWorkflow(workflow: UpdatedWorkflow): Boolean {
        // Check that all referenced nodes have service parameters
        val referencedNodes = HashSet<Int>()
        
        for (i in 0 until workflow.linkMatrix.getNumRows()) {
            referencedNodes.add(workflow.linkMatrix.get(i, 0).toInt())
            referencedNodes.add(workflow.linkMatrix.get(i, 1).toInt())
        }
        
        // Check that all service nodes have parameters
        for (node in referencedNodes) {
            if (node !in workflow.serviceParameters) {
                return false
            }
        }
        
        return true
    }
    
    /**
     * Get statistics about the pattern update process.
     */
    @JvmStatic
    fun getUpdateStats(
        originalMatrix: Matrix,
        updatedWorkflow: UpdatedWorkflow
    ): Map<String, Any> {
        val stats = HashMap<String, Any>()
        
        stats["originalLinks"] = originalMatrix.getNumRows()
        stats["updatedLinks"] = updatedWorkflow.linkMatrix.getNumRows()
        stats["linksReduced"] = originalMatrix.getNumRows() - updatedWorkflow.linkMatrix.getNumRows()
        stats["serviceNodes"] = updatedWorkflow.serviceParameters.size
        
        val reductionRatio = if (originalMatrix.getNumRows() > 0) {
            (originalMatrix.getNumRows() - updatedWorkflow.linkMatrix.getNumRows()).toDouble() / 
            originalMatrix.getNumRows()
        } else 0.0
        stats["reductionRatio"] = reductionRatio
        
        return stats
    }
}