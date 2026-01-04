/**
 * @file Workflow Branch Pattern Detection
 * 
 * Implements algorithms for detecting branching patterns in workflow traces.
 * Identifies decision points, conditional paths, and alternative execution
 * flows in business process analysis and workflow pattern recognition.
 * 
 * @since LINE 3.0
 */
package jline.api.wf

import jline.util.matrix.Matrix
import java.util.*
import kotlin.collections.ArrayList

/**
 * Automatic branch structure detector for workflow networks.
 * 
 * This detector identifies branch patterns in workflow networks where
 * execution splits into multiple paths with different probabilities
 * and later converges at join points.
 * 
 * Based on AUTO_Branch_Detector.m from the MDN toolbox.
 */
object Wf_branch_detector {
    
    /**
     * Data class to represent a branch pattern with probabilities.
     */
    data class BranchPattern(
        val branchNodes: List<Int>,
        val probabilities: List<Double>,
        val forkNode: Int?,
        val joinNode: Int?
    )
    
    /**
     * Detect branch patterns in a workflow network.
     * 
     * @param linkMatrix Matrix containing workflow link information:
     *                   - Column 1: start node IDs
     *                   - Column 2: end node IDs
     *                   - Column 3: transition probabilities
     * @param serviceNodes List of service node IDs
     * @param joinNodes List of join node IDs
     * @return List of branch patterns with their associated probabilities
     */
    @JvmStatic
    fun detectBranches(
        linkMatrix: Matrix,
        serviceNodes: List<Int>,
        joinNodes: List<Int>
    ): List<BranchPattern> {
        val branchPatterns = ArrayList<BranchPattern>()
        val serviceSet = serviceNodes.toSet()
        val joinSet = joinNodes.toSet()
        
        // Build adjacency maps
        val adjacency = buildAdjacencyMap(linkMatrix)
        val reverseAdjacency = buildReverseAdjacencyMap(linkMatrix)
        
        // Find nodes that have multiple outgoing edges (potential branch points)
        val branchPoints = findBranchPoints(adjacency, serviceSet, joinSet)
        
        for (branchPoint in branchPoints) {
            val pattern = analyzeBranchPattern(
                branchPoint, 
                adjacency, 
                reverseAdjacency, 
                serviceSet, 
                joinSet
            )
            if (pattern != null && pattern.branchNodes.size > 1) {
                branchPatterns.add(pattern)
            }
        }
        
        return branchPatterns
    }
    
    /**
     * Build forward adjacency map with probabilities.
     */
    private fun buildAdjacencyMap(linkMatrix: Matrix): Map<Int, List<Pair<Int, Double>>> {
        val adjacency = HashMap<Int, MutableList<Pair<Int, Double>>>()
        
        for (i in 0 until linkMatrix.getNumRows()) {
            val start = linkMatrix.get(i, 0).toInt()
            val end = linkMatrix.get(i, 1).toInt()
            val prob = linkMatrix.get(i, 2)
            
            adjacency.computeIfAbsent(start) { ArrayList() }.add(Pair(end, prob))
        }
        
        return adjacency
    }
    
    /**
     * Build reverse adjacency map.
     */
    private fun buildReverseAdjacencyMap(linkMatrix: Matrix): Map<Int, List<Int>> {
        val reverseAdj = HashMap<Int, MutableList<Int>>()
        
        for (i in 0 until linkMatrix.getNumRows()) {
            val start = linkMatrix.get(i, 0).toInt()
            val end = linkMatrix.get(i, 1).toInt()
            
            reverseAdj.computeIfAbsent(end) { ArrayList() }.add(start)
        }
        
        return reverseAdj
    }
    
    /**
     * Find potential branch points (nodes with multiple outgoing edges).
     */
    private fun findBranchPoints(
        adjacency: Map<Int, List<Pair<Int, Double>>>,
        serviceSet: Set<Int>,
        joinSet: Set<Int>
    ): List<Int> {
        val branchPoints = ArrayList<Int>()
        
        for ((node, neighbors) in adjacency) {
            if (neighbors.size > 1) {
                // Check if this creates a valid branch pattern
                val serviceTargets = neighbors.count { it.first in serviceSet }
                if (serviceTargets > 1) {
                    branchPoints.add(node)
                }
            }
        }
        
        return branchPoints
    }
    
    /**
     * Analyze a potential branch pattern starting from a branch point.
     */
    private fun analyzeBranchPattern(
        branchPoint: Int,
        adjacency: Map<Int, List<Pair<Int, Double>>>,
        reverseAdjacency: Map<Int, List<Int>>,
        serviceSet: Set<Int>,
        joinSet: Set<Int>
    ): BranchPattern? {
        val neighbors = adjacency[branchPoint] ?: return null
        
        // Find service nodes that are direct targets of the branch
        val branchTargets = ArrayList<Pair<Int, Double>>()
        for ((target, prob) in neighbors) {
            if (target in serviceSet) {
                branchTargets.add(Pair(target, prob))
            }
        }
        
        if (branchTargets.size < 2) return null
        
        // Find common join point for these branches
        val commonJoin = findCommonJoinPoint(
            branchTargets.map { it.first },
            adjacency,
            reverseAdjacency,
            joinSet
        )
        
        // Validate probabilities sum to 1.0 (or close to it)
        val totalProb = branchTargets.sumOf { it.second }
        if (kotlin.math.abs(totalProb - 1.0) > 0.01) {
            // Not a valid probabilistic branch
            return null
        }
        
        return BranchPattern(
            branchNodes = branchTargets.map { it.first },
            probabilities = branchTargets.map { it.second },
            forkNode = branchPoint,
            joinNode = commonJoin
        )
    }
    
    /**
     * Find the common join point for a set of branch nodes.
     */
    private fun findCommonJoinPoint(
        branchNodes: List<Int>,
        adjacency: Map<Int, List<Pair<Int, Double>>>,
        reverseAdjacency: Map<Int, List<Int>>,
        joinSet: Set<Int>
    ): Int? {
        // Find all nodes reachable from each branch
        val reachableSets = branchNodes.map { branch ->
            findReachableNodes(branch, adjacency, joinSet)
        }
        
        // Find intersection of all reachable sets
        val commonReachable = reachableSets.reduce { acc, set -> acc.intersect(set) }
        
        // Prefer join nodes as common points
        val joinPoints = commonReachable.intersect(joinSet)
        if (joinPoints.isNotEmpty()) {
            // Return the "closest" join point (heuristic: first in topological order)
            return joinPoints.first()
        }
        
        // If no join nodes, return any common reachable node
        return commonReachable.firstOrNull()
    }
    
    /**
     * Find all nodes reachable from a starting node.
     */
    private fun findReachableNodes(
        startNode: Int,
        adjacency: Map<Int, List<Pair<Int, Double>>>,
        stopSet: Set<Int>
    ): Set<Int> {
        val reachable = HashSet<Int>()
        val queue = LinkedList<Int>()
        val visited = HashSet<Int>()
        
        queue.offer(startNode)
        
        while (queue.isNotEmpty()) {
            val current = queue.poll()
            if (current in visited) continue
            visited.add(current)
            
            val neighbors = adjacency[current] ?: emptyList()
            for ((neighbor, _) in neighbors) {
                reachable.add(neighbor)
                if (neighbor !in stopSet) {
                    queue.offer(neighbor)
                }
            }
        }
        
        return reachable
    }
    
    /**
     * Validate branch pattern structure.
     */
    @JvmStatic
    fun validateBranchPattern(
        pattern: BranchPattern,
        linkMatrix: Matrix
    ): Boolean {
        // Check probability sum
        val totalProb = pattern.probabilities.sum()
        if (kotlin.math.abs(totalProb - 1.0) > 0.01) {
            return false
        }
        
        // Check that all branch nodes are reachable from fork
        if (pattern.forkNode == null) return false
        
        val adjacency = buildAdjacencyMap(linkMatrix)
        val forkNeighbors = adjacency[pattern.forkNode] ?: return false
        
        val forkTargets = forkNeighbors.map { it.first }.toSet()
        return pattern.branchNodes.all { it in forkTargets }
    }
    
    /**
     * Calculate branch diversity metrics.
     */
    @JvmStatic
    fun calculateBranchDiversity(pattern: BranchPattern): Map<String, Double> {
        val metrics = HashMap<String, Double>()
        
        val probs = pattern.probabilities
        val n = probs.size
        
        // Shannon entropy
        val entropy = -probs.sumOf { p -> 
            if (p > 0) p * kotlin.math.ln(p) else 0.0 
        }
        metrics["entropy"] = entropy
        
        // Normalized entropy
        metrics["normalizedEntropy"] = if (n > 1) entropy / kotlin.math.ln(n.toDouble()) else 0.0
        
        // Gini coefficient (inequality measure)
        val sortedProbs = probs.sorted()
        var gini = 0.0
        for (i in sortedProbs.indices) {
            gini += (2 * (i + 1) - n - 1) * sortedProbs[i]
        }
        gini /= (n - 1) * probs.sum()
        metrics["gini"] = kotlin.math.abs(gini)
        
        // Balance (inverse of max probability)
        metrics["balance"] = 1.0 / (probs.maxOrNull() ?: 1.0)
        
        return metrics
    }
    
    /**
     * Get branch pattern statistics.
     */
    @JvmStatic
    fun getBranchStats(patterns: List<BranchPattern>): Map<String, Any> {
        val stats = HashMap<String, Any>()
        
        stats["numPatterns"] = patterns.size
        stats["totalBranchNodes"] = patterns.sumOf { it.branchNodes.size }
        
        val branchCounts = patterns.map { it.branchNodes.size }
        stats["avgBranches"] = branchCounts.average()
        stats["maxBranches"] = branchCounts.maxOrNull() ?: 0
        stats["minBranches"] = branchCounts.minOrNull() ?: 0
        
        // Diversity statistics
        val diversityMetrics = patterns.map { calculateBranchDiversity(it) }
        stats["avgEntropy"] = diversityMetrics.mapNotNull { it["entropy"] }.average()
        stats["avgBalance"] = diversityMetrics.mapNotNull { it["balance"] }.average()
        
        return stats
    }
    
    /**
     * Find the most probable branch in a pattern.
     */
    @JvmStatic
    fun findMostProbableBranch(pattern: BranchPattern): Pair<Int, Double>? {
        if (pattern.branchNodes.isEmpty()) return null
        
        val maxIndex = pattern.probabilities.indices.maxByOrNull { pattern.probabilities[it] }
            ?: return null
        
        return Pair(pattern.branchNodes[maxIndex], pattern.probabilities[maxIndex])
    }
    
    /**
     * Find the least probable branch in a pattern.
     */
    @JvmStatic
    fun findLeastProbableBranch(pattern: BranchPattern): Pair<Int, Double>? {
        if (pattern.branchNodes.isEmpty()) return null
        
        val minIndex = pattern.probabilities.indices.minByOrNull { pattern.probabilities[it] }
            ?: return null
        
        return Pair(pattern.branchNodes[minIndex], pattern.probabilities[minIndex])
    }
}