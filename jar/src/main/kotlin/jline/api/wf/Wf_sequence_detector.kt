/**
 * @file Workflow Sequence Pattern Detection
 * 
 * Implements algorithms for detecting sequential patterns in workflow traces.
 * Identifies ordered sequences of activities and their temporal relationships
 * for workflow analysis and pattern recognition applications.
 * 
 * @since LINE 3.0
 */
package jline.api.wf

import jline.util.matrix.Matrix
import java.util.*
import kotlin.collections.ArrayList

/**
 * Automatic sequence structure detector for workflow networks.
 * 
 * This detector can automatically identify basic sequence structures 
 * in a workflow network, where service nodes are connected sequentially.
 * 
 * Based on AUTO_Sequence_Detector.m from the MDN toolbox.
 */
object Wf_sequence_detector {
    
    /**
     * Detect sequence patterns in a workflow network.
     * 
     * @param linkMatrix Matrix containing workflow link information:
     *                   - Column 1: start node IDs  
     *                   - Column 2: end node IDs
     *                   - Column 3: transition probabilities
     * @param serviceNodes List of service node IDs
     * @return List of sequence chains, where each chain is a list of connected service nodes
     */
    @JvmStatic
    fun detectSequences(linkMatrix: Matrix, serviceNodes: List<Int>): List<List<Int>> {
        val chains = ArrayList<List<Int>>()
        
        // Step 1: Find connections between service nodes
        val serviceConnections = findServiceConnections(linkMatrix, serviceNodes)
        if (serviceConnections.isEmpty()) {
            return chains // No sequence structures found
        }
        
        // Step 2: Count occurrences of each service node
        val nodeCounts = countNodeOccurrences(serviceConnections)
        
        // Step 3: Detect sequences using node occurrence patterns
        val numSequences = nodeCounts.values.count { it == 1 } / 2
        
        val remainingConnections = ArrayList(serviceConnections)
        
        for (seqIndex in 0 until numSequences) {
            if (remainingConnections.isEmpty()) break
            
            val sequence = buildSequenceChain(remainingConnections)
            if (sequence.isNotEmpty()) {
                chains.add(sequence)
            }
        }
        
        return chains
    }
    
    /**
     * Find all connections between service nodes in the workflow.
     */
    private fun findServiceConnections(linkMatrix: Matrix, serviceNodes: List<Int>): List<Pair<Int, Int>> {
        val connections = ArrayList<Pair<Int, Int>>()
        val serviceSet = serviceNodes.toSet()
        
        for (i in 0 until linkMatrix.getNumRows()) {
            val startNode = linkMatrix.get(i, 0).toInt()
            val endNode = linkMatrix.get(i, 1).toInt()
            
            // Check if both nodes are service nodes
            if (startNode in serviceSet && endNode in serviceSet) {
                connections.add(Pair(startNode, endNode))
            }
        }
        
        return connections
    }
    
    /**
     * Count how many times each service node appears in connections.
     */
    private fun countNodeOccurrences(connections: List<Pair<Int, Int>>): Map<Int, Int> {
        val counts = HashMap<Int, Int>()
        
        for ((start, end) in connections) {
            counts[start] = counts.getOrDefault(start, 0) + 1
            counts[end] = counts.getOrDefault(end, 0) + 1
        }
        
        return counts
    }
    
    /**
     * Build a sequence chain starting from the first available connection.
     */
    private fun buildSequenceChain(connections: MutableList<Pair<Int, Int>>): List<Int> {
        if (connections.isEmpty()) return emptyList()
        
        val sequence = ArrayList<Int>()
        val usedConnections = ArrayList<Int>()
        
        // Start with first connection
        var (first, last) = connections[0]
        sequence.add(first)
        sequence.add(last)
        usedConnections.add(0)
        
        var foundExtension = true
        while (foundExtension) {
            foundExtension = false
            val currentSize = sequence.size
            
            // Try to extend forward or backward
            for (i in 1 until connections.size) {
                if (i in usedConnections) continue
                
                val (start, end) = connections[i]
                
                when {
                    start == last -> {
                        // Extend forward
                        last = end
                        sequence.add(end)
                        usedConnections.add(i)
                        foundExtension = true
                    }
                    end == first -> {
                        // Extend backward
                        first = start
                        sequence.add(0, start)
                        usedConnections.add(i)
                        foundExtension = true
                    }
                }
            }
            
            // Check if we made progress
            foundExtension = foundExtension && sequence.size > currentSize
        }
        
        // Remove used connections from the list
        usedConnections.sortDescending()
        for (index in usedConnections) {
            connections.removeAt(index)
        }
        
        return sequence
    }
    
    /**
     * Validate that a sequence chain is properly connected.
     */
    @JvmStatic
    fun validateSequence(sequence: List<Int>, linkMatrix: Matrix): Boolean {
        if (sequence.size < 2) return false
        
        val connections = HashSet<Pair<Int, Int>>()
        for (i in 0 until linkMatrix.getNumRows()) {
            val start = linkMatrix.get(i, 0).toInt()
            val end = linkMatrix.get(i, 1).toInt()
            connections.add(Pair(start, end))
        }
        
        // Check all consecutive pairs in sequence
        for (i in 0 until sequence.size - 1) {
            val connection = Pair(sequence[i], sequence[i + 1])
            if (connection !in connections) {
                return false
            }
        }
        
        return true
    }
    
    /**
     * Get sequence statistics for analysis.
     */
    @JvmStatic
    fun getSequenceStats(sequences: List<List<Int>>): Map<String, Any> {
        val stats = HashMap<String, Any>()
        
        stats["numSequences"] = sequences.size
        stats["totalNodes"] = sequences.sumOf { it.size }
        stats["avgLength"] = if (sequences.isNotEmpty()) {
            sequences.map { it.size }.average()
        } else 0.0
        stats["maxLength"] = sequences.maxByOrNull { it.size }?.size ?: 0
        stats["minLength"] = sequences.minByOrNull { it.size }?.size ?: 0
        
        return stats
    }
}