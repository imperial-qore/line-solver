/**
 * @file Workflow Parallel Pattern Detection
 * 
 * Implements algorithms for detecting parallel execution patterns in workflow
 * traces. Identifies concurrent activities, synchronization points, and
 * parallel flow structures in business process analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.wf

import jline.util.matrix.Matrix
import java.util.*
import kotlin.collections.ArrayList

/**
 * Automatic parallel structure detector for workflow networks.
 * 
 * This detector identifies parallel execution patterns in workflow networks
 * where multiple service nodes execute concurrently between fork and join points.
 * 
 * Based on AUTO_Parallel_Detector.m from the MDN toolbox.
 */
object Wf_parallel_detector {
    
    /**
     * Detect parallel patterns in a workflow network.
     * 
     * @param linkMatrix Matrix containing workflow link information:
     *                   - Column 1: start node IDs
     *                   - Column 2: end node IDs  
     *                   - Column 3: transition probabilities
     * @param serviceNodes List of service node IDs
     * @param forkNodes List of fork node IDs
     * @param joinNodes List of join node IDs
     * @return List of parallel patterns, where each pattern is a list of parallel service nodes
     */
    @JvmStatic
    fun detectParallel(
        linkMatrix: Matrix,
        serviceNodes: List<Int>,
        forkNodes: List<Int>,
        joinNodes: List<Int>
    ): List<List<Int>> {
        val parallelPatterns = ArrayList<List<Int>>()
        
        // Find all fork-join pairs
        val forkJoinPairs = findForkJoinPairs(linkMatrix, forkNodes, joinNodes)
        
        for ((forkNode, joinNode) in forkJoinPairs) {
            val parallelServices = findParallelServices(linkMatrix, serviceNodes, forkNode, joinNode)
            if (parallelServices.size > 1) {
                parallelPatterns.add(parallelServices)
            }
        }
        
        return parallelPatterns
    }
    
    /**
     * Find all valid fork-join pairs in the workflow.
     */
    private fun findForkJoinPairs(
        linkMatrix: Matrix,
        forkNodes: List<Int>,
        joinNodes: List<Int>
    ): List<Pair<Int, Int>> {
        val pairs = ArrayList<Pair<Int, Int>>()
        val forkSet = forkNodes.toSet()
        val joinSet = joinNodes.toSet()
        
        // Build adjacency map for reachability analysis
        val adjacency = buildAdjacencyMap(linkMatrix)
        
        for (fork in forkNodes) {
            for (join in joinNodes) {
                if (isValidForkJoinPair(fork, join, adjacency, forkSet, joinSet)) {
                    pairs.add(Pair(fork, join))
                }
            }
        }
        
        return pairs
    }
    
    /**
     * Build adjacency map from link matrix.
     */
    private fun buildAdjacencyMap(linkMatrix: Matrix): Map<Int, List<Int>> {
        val adjacency = HashMap<Int, MutableList<Int>>()
        
        for (i in 0 until linkMatrix.getNumRows()) {
            val start = linkMatrix.get(i, 0).toInt()
            val end = linkMatrix.get(i, 1).toInt()
            
            adjacency.computeIfAbsent(start) { ArrayList() }.add(end)
        }
        
        return adjacency
    }
    
    /**
     * Check if fork and join nodes form a valid pair.
     */
    private fun isValidForkJoinPair(
        fork: Int,
        join: Int,
        adjacency: Map<Int, List<Int>>,
        forkSet: Set<Int>,
        joinSet: Set<Int>
    ): Boolean {
        // Use BFS to check if there are multiple paths from fork to join
        val queue = LinkedList<Int>()
        val visited = HashSet<Int>()
        val pathCount = HashMap<Int, Int>()
        
        queue.offer(fork)
        pathCount[fork] = 1
        
        while (queue.isNotEmpty()) {
            val current = queue.poll()
            if (current in visited) continue
            visited.add(current)
            
            val neighbors = adjacency[current] ?: emptyList()
            for (neighbor in neighbors) {
                if (neighbor == join) {
                    // Found path to join
                    val currentPaths = pathCount[current] ?: 0
                    pathCount[join] = pathCount.getOrDefault(join, 0) + currentPaths
                } else if (neighbor !in visited && neighbor !in forkSet && neighbor !in joinSet) {
                    // Continue exploring (avoid other fork/join nodes)
                    queue.offer(neighbor)
                    val currentPaths = pathCount[current] ?: 0
                    pathCount[neighbor] = pathCount.getOrDefault(neighbor, 0) + currentPaths
                }
            }
        }
        
        // Valid if there are multiple paths to the join
        return (pathCount[join] ?: 0) > 1
    }
    
    /**
     * Find service nodes that execute in parallel between fork and join.
     */
    private fun findParallelServices(
        linkMatrix: Matrix,
        serviceNodes: List<Int>,
        forkNode: Int,
        joinNode: Int
    ): List<Int> {
        val parallelServices = ArrayList<Int>()
        val serviceSet = serviceNodes.toSet()
        
        // Find all nodes reachable from fork
        val reachableFromFork = findReachableNodes(linkMatrix, forkNode, joinNode)
        
        // Find all nodes that can reach join
        val canReachJoin = findNodesThatCanReach(linkMatrix, joinNode, forkNode)
        
        // Intersection gives nodes in parallel paths
        val parallelNodes = reachableFromFork.intersect(canReachJoin)
        
        // Filter to only include service nodes
        for (node in parallelNodes) {
            if (node in serviceSet) {
                parallelServices.add(node)
            }
        }
        
        return parallelServices
    }
    
    /**
     * Find all nodes reachable from a start node (stopping at end node).
     */
    private fun findReachableNodes(linkMatrix: Matrix, startNode: Int, endNode: Int): Set<Int> {
        val reachable = HashSet<Int>()
        val queue = LinkedList<Int>()
        val visited = HashSet<Int>()
        
        queue.offer(startNode)
        
        while (queue.isNotEmpty()) {
            val current = queue.poll()
            if (current in visited || current == endNode) continue
            visited.add(current)
            
            for (i in 0 until linkMatrix.getNumRows()) {
                val start = linkMatrix.get(i, 0).toInt()
                val end = linkMatrix.get(i, 1).toInt()
                
                if (start == current && end != endNode) {
                    reachable.add(end)
                    queue.offer(end)
                }
            }
        }
        
        return reachable
    }
    
    /**
     * Find all nodes that can reach a target node (starting from start node).
     */
    private fun findNodesThatCanReach(linkMatrix: Matrix, targetNode: Int, startNode: Int): Set<Int> {
        val canReach = HashSet<Int>()
        
        // Build reverse adjacency map
        val reverseAdj = HashMap<Int, MutableList<Int>>()
        for (i in 0 until linkMatrix.getNumRows()) {
            val start = linkMatrix.get(i, 0).toInt()
            val end = linkMatrix.get(i, 1).toInt()
            reverseAdj.computeIfAbsent(end) { ArrayList() }.add(start)
        }
        
        // BFS backward from target
        val queue = LinkedList<Int>()
        val visited = HashSet<Int>()
        
        queue.offer(targetNode)
        
        while (queue.isNotEmpty()) {
            val current = queue.poll()
            if (current in visited || current == startNode) continue
            visited.add(current)
            
            val predecessors = reverseAdj[current] ?: emptyList()
            for (pred in predecessors) {
                if (pred != startNode) {
                    canReach.add(pred)
                    queue.offer(pred)
                }
            }
        }
        
        return canReach
    }
    
    /**
     * Validate parallel pattern structure.
     */
    @JvmStatic
    fun validateParallelPattern(
        pattern: List<Int>,
        linkMatrix: Matrix,
        forkNodes: List<Int>,
        joinNodes: List<Int>
    ): Boolean {
        if (pattern.size < 2) return false
        
        // Find the fork and join nodes for this pattern
        val forkNode = findCommonSource(pattern, linkMatrix, forkNodes)
        val joinNode = findCommonTarget(pattern, linkMatrix, joinNodes)
        
        return forkNode != null && joinNode != null
    }
    
    private fun findCommonSource(pattern: List<Int>, linkMatrix: Matrix, forkNodes: List<Int>): Int? {
        val forkSet = forkNodes.toSet()
        val sources = HashSet<Int>()
        
        for (node in pattern) {
            for (i in 0 until linkMatrix.getNumRows()) {
                val start = linkMatrix.get(i, 0).toInt()
                val end = linkMatrix.get(i, 1).toInt()
                
                if (end == node && start in forkSet) {
                    sources.add(start)
                }
            }
        }
        
        return if (sources.size == 1) sources.first() else null
    }
    
    private fun findCommonTarget(pattern: List<Int>, linkMatrix: Matrix, joinNodes: List<Int>): Int? {
        val joinSet = joinNodes.toSet()
        val targets = HashSet<Int>()
        
        for (node in pattern) {
            for (i in 0 until linkMatrix.getNumRows()) {
                val start = linkMatrix.get(i, 0).toInt()
                val end = linkMatrix.get(i, 1).toInt()
                
                if (start == node && end in joinSet) {
                    targets.add(end)
                }
            }
        }
        
        return if (targets.size == 1) targets.first() else null
    }
    
    /**
     * Get parallel pattern statistics.
     */
    @JvmStatic
    fun getParallelStats(patterns: List<List<Int>>): Map<String, Any> {
        val stats = HashMap<String, Any>()
        
        stats["numPatterns"] = patterns.size
        stats["totalParallelNodes"] = patterns.sumOf { it.size }
        stats["avgParallelism"] = if (patterns.isNotEmpty()) {
            patterns.map { it.size }.average()
        } else 0.0
        stats["maxParallelism"] = patterns.maxByOrNull { it.size }?.size ?: 0
        
        return stats
    }
}