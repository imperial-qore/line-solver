/**
 * @file Workflow Loop Pattern Detection
 * 
 * Implements algorithms for detecting loop and iterative patterns in workflow
 * traces. Identifies repetitive activity sequences, cyclic behaviors, and
 * iteration structures in business process analysis and workflow mining.
 * 
 * @since LINE 3.0
 */
package jline.api.wf

import jline.GlobalConstants.Inf
import jline.util.matrix.Matrix
import java.util.*
import kotlin.collections.ArrayList

/**
 * Automatic loop structure detector for workflow networks.
 * 
 * This detector identifies loop patterns in workflow networks where
 * service nodes are connected in cycles through router nodes.
 * 
 * Based on AUTO_Loop_Detector.m from the MDN toolbox.
 */
object Wf_loop_detector {
    
    /**
     * Detect loop patterns in a workflow network.
     * 
     * @param linkMatrix Matrix containing workflow link information:
     *                   - Column 1: start node IDs
     *                   - Column 2: end node IDs
     *                   - Column 3: transition probabilities
     * @param serviceNodes List of service node IDs
     * @param routerNodes List of router node IDs
     * @param joinNodes List of join node IDs (optional, for complex loops)
     * @return List of loop patterns, where each pattern contains the loop service node
     */
    @JvmStatic
    fun detectLoops(
        linkMatrix: Matrix,
        serviceNodes: List<Int>,
        routerNodes: List<Int>,
        joinNodes: List<Int> = emptyList()
    ): List<Int> {
        val loopNodes = ArrayList<Int>()
        val serviceSet = serviceNodes.toSet()
        val routerSet = routerNodes.toSet()
        
        // Build adjacency map
        val adjacency = buildAdjacencyMap(linkMatrix)
        
        // Find simple loops: service -> router -> service
        for (serviceNode in serviceNodes) {
            if (isInSimpleLoop(serviceNode, adjacency, routerSet)) {
                loopNodes.add(serviceNode)
            }
        }
        
        // Find complex loops involving join nodes
        if (joinNodes.isNotEmpty()) {
            val complexLoops = findComplexLoops(linkMatrix, serviceNodes, routerNodes, joinNodes)
            loopNodes.addAll(complexLoops)
        }
        
        return loopNodes.distinct()
    }
    
    /**
     * Build adjacency map from link matrix.
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
     * Check if a service node is in a simple loop pattern.
     */
    private fun isInSimpleLoop(
        serviceNode: Int,
        adjacency: Map<Int, List<Pair<Int, Double>>>,
        routerSet: Set<Int>
    ): Boolean {
        val neighbors = adjacency[serviceNode] ?: return false
        
        // Check if service connects to a router
        for ((routerNode, _) in neighbors) {
            if (routerNode in routerSet) {
                // Check if router connects back to service
                val routerNeighbors = adjacency[routerNode] ?: emptyList()
                for ((backNode, _) in routerNeighbors) {
                    if (backNode == serviceNode) {
                        return true
                    }
                }
            }
        }
        
        return false
    }
    
    /**
     * Find complex loops involving multiple nodes and join points.
     */
    private fun findComplexLoops(
        linkMatrix: Matrix,
        serviceNodes: List<Int>,
        routerNodes: List<Int>,
        joinNodes: List<Int>
    ): List<Int> {
        val complexLoops = ArrayList<Int>()
        val serviceSet = serviceNodes.toSet()
        val routerSet = routerNodes.toSet()
        val joinSet = joinNodes.toSet()
        
        // Build graph for cycle detection
        val graph = buildDirectedGraph(linkMatrix)
        
        // Find strongly connected components
        val sccs = findStronglyConnectedComponents(graph)
        
        for (scc in sccs) {
            if (scc.size > 1) {
                // Check if SCC contains service nodes and forms a valid loop
                val serviceNodesInSCC = scc.filter { it in serviceSet }
                val hasRouterOrJoin = scc.any { it in routerSet || it in joinSet }
                
                if (serviceNodesInSCC.isNotEmpty() && hasRouterOrJoin) {
                    complexLoops.addAll(serviceNodesInSCC)
                }
            }
        }
        
        return complexLoops
    }
    
    /**
     * Build directed graph representation.
     */
    private fun buildDirectedGraph(linkMatrix: Matrix): Map<Int, Set<Int>> {
        val graph = HashMap<Int, MutableSet<Int>>()
        
        for (i in 0 until linkMatrix.getNumRows()) {
            val start = linkMatrix.get(i, 0).toInt()
            val end = linkMatrix.get(i, 1).toInt()
            
            graph.computeIfAbsent(start) { HashSet() }.add(end)
        }
        
        return graph
    }
    
    /**
     * Find strongly connected components using Tarjan's algorithm.
     */
    private fun findStronglyConnectedComponents(graph: Map<Int, Set<Int>>): List<List<Int>> {
        val sccs = ArrayList<List<Int>>()
        val visited = HashSet<Int>()
        val stack = ArrayDeque<Int>()
        val indices = HashMap<Int, Int>()
        val lowLinks = HashMap<Int, Int>()
        val onStack = HashSet<Int>()
        var index = 0
        
        fun strongConnect(node: Int) {
            indices[node] = index
            lowLinks[node] = index
            index++
            stack.push(node)
            onStack.add(node)
            
            val neighbors = graph[node] ?: emptySet()
            for (neighbor in neighbors) {
                when {
                    neighbor !in indices -> {
                        strongConnect(neighbor)
                        lowLinks[node] = minOf(lowLinks[node]!!, lowLinks[neighbor]!!)
                    }
                    neighbor in onStack -> {
                        lowLinks[node] = minOf(lowLinks[node]!!, indices[neighbor]!!)
                    }
                }
            }
            
            if (lowLinks[node] == indices[node]) {
                val scc = ArrayList<Int>()
                var w: Int
                do {
                    w = stack.pop()
                    onStack.remove(w)
                    scc.add(w)
                } while (w != node)
                sccs.add(scc)
            }
        }
        
        // Run Tarjan's algorithm on all unvisited nodes
        for (node in graph.keys) {
            if (node !in visited) {
                strongConnect(node)
                visited.add(node)
            }
        }
        
        return sccs
    }
    
    /**
     * Get loop probability for a service node.
     */
    @JvmStatic
    fun getLoopProbability(
        serviceNode: Int,
        linkMatrix: Matrix,
        routerNodes: List<Int>
    ): Double {
        val routerSet = routerNodes.toSet()
        
        // Find the router that creates the loop
        for (i in 0 until linkMatrix.getNumRows()) {
            val start = linkMatrix.get(i, 0).toInt()
            val end = linkMatrix.get(i, 1).toInt()
            val prob = linkMatrix.get(i, 2)
            
            if (start == serviceNode && end in routerSet) {
                // Found service -> router connection
                // Now find router -> service connection (loop back)
                for (j in 0 until linkMatrix.getNumRows()) {
                    val loopStart = linkMatrix.get(j, 0).toInt()
                    val loopEnd = linkMatrix.get(j, 1).toInt()
                    val loopProb = linkMatrix.get(j, 2)
                    
                    if (loopStart == end && loopEnd == serviceNode) {
                        return loopProb
                    }
                }
            }
        }
        
        return 0.0
    }
    
    /**
     * Validate loop pattern structure.
     */
    @JvmStatic
    fun validateLoopPattern(
        loopNode: Int,
        linkMatrix: Matrix,
        routerNodes: List<Int>
    ): Boolean {
        val routerSet = routerNodes.toSet()
        
        // Check if there's a path: loopNode -> router -> loopNode
        val adjacency = buildAdjacencyMap(linkMatrix)
        val neighbors = adjacency[loopNode] ?: return false
        
        for ((routerNode, _) in neighbors) {
            if (routerNode in routerSet) {
                val routerNeighbors = adjacency[routerNode] ?: emptyList()
                for ((backNode, _) in routerNeighbors) {
                    if (backNode == loopNode) {
                        return true
                    }
                }
            }
        }
        
        return false
    }
    
    /**
     * Calculate expected number of loop iterations.
     */
    @JvmStatic
    fun getExpectedLoopIterations(loopProbability: Double): Double {
        return if (loopProbability >= 1.0) {
            Inf
        } else {
            1.0 / (1.0 - loopProbability)
        }
    }
    
    /**
     * Get loop pattern statistics.
     */
    @JvmStatic
    fun getLoopStats(
        loopNodes: List<Int>,
        linkMatrix: Matrix,
        routerNodes: List<Int>
    ): Map<String, Any> {
        val stats = HashMap<String, Any>()
        
        stats["numLoops"] = loopNodes.size
        
        val probabilities = loopNodes.map { getLoopProbability(it, linkMatrix, routerNodes) }
        stats["avgLoopProbability"] = probabilities.average()
        stats["maxLoopProbability"] = probabilities.maxOrNull() ?: 0.0
        stats["minLoopProbability"] = probabilities.minOrNull() ?: 0.0
        
        val iterations = probabilities.map { getExpectedLoopIterations(it) }
            .filter { it.isFinite() }
        stats["avgExpectedIterations"] = iterations.average()
        stats["maxExpectedIterations"] = iterations.maxOrNull() ?: 0.0
        
        return stats
    }
}