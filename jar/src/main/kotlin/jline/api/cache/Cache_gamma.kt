/**
 * @file Cache Access Factor Computation
 * 
 * Computes cache access factors from request arrival rates and routing matrices.
 * Access factors quantify the effective demand placed on each cache level
 * by different user streams and item types.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.io.Ret
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes access factors for the cache.
 *
 * @param lambda - MatrixCell representing request arrival rates from users to items of individual lists.
 * @param R      - MatrixCell of MatrixCell representing the reachability graph of a list for different streams and items.
 * @return cacheGamma - An object containing access factors (gamma), the number of users (u),
 * the number of items (n), and the number of lists (h).
 */
fun cache_gamma(lambda: MatrixCell, R: MatrixCell): Ret.cacheGamma {
    val u = lambda.size() // Number of users
    val n = lambda.get(0).numRows // Number of items
    val h = lambda.get(0).numCols - 1 // Number of lists
    val gamma = Matrix(n, h)

    for (i in 0..<n) { // for all items
        for (j in 0..<h) { // for all levels
            // Compute gamma(i,j)
            
            // Create a directed graph from R{1,i} (topology from stream 1)
            val RMatrixCell = R.get(0) as MatrixCell
            val graph = RMatrixCell.get(i)
            
            // Find shortest path from node 1 (0-indexed) to node j (0-indexed)
            val Pj = findShortestPath(graph, 0, j)
            
            if (Pj.isEmpty()) {
                gamma[i, j] = 0.0
            } else {
                // Initialize gamma(i,j) with sum of lambda(:,i,1) for level 1
                var gammaValue = 0.0
                for (v in 0..<u) {
                    gammaValue += lambda.get(v)[i, 0] // lambda(v, i, 1+0)
                }
                
                // For all levels up to the current one
                for (li in 1..<Pj.size) {
                    var y = 0.0
                    val l_1 = Pj[li - 1]
                    val l = Pj[li]
                    
                    // For all streams
                    for (v in 0..<u) {
                        val RvMatrixCell = R.get(v) as MatrixCell
                        y += lambda.get(v)[i, l_1] * RvMatrixCell.get(i)[l_1, l]
                    }
                    gammaValue *= y
                }
                gamma[i, j] = gammaValue
            }
        }
    }
    
    return Ret.cacheGamma(gamma, u, n, h)
}

/**
 * Finds the shortest path from source to destination in the graph represented by the adjacency matrix.
 * Uses BFS since all edges have equal weight in the cache reachability graph.
 *
 * @param adjacencyMatrix - Matrix representing the graph
 * @param source - Source node (0-indexed)
 * @param destination - Destination node (0-indexed)
 * @return List of nodes in the shortest path from source to destination
 */
private fun findShortestPath(adjacencyMatrix: Matrix, source: Int, destination: Int): List<Int> {
    val n = adjacencyMatrix.numRows
    if (source >= n || destination >= n || source < 0 || destination < 0) {
        return emptyList()
    }
    
    // BFS to find shortest path
    val visited = BooleanArray(n)
    val parent = IntArray(n) { -1 }
    val queue = ArrayDeque<Int>()
    
    queue.add(source)
    visited[source] = true
    
    while (queue.isNotEmpty()) {
        val current = queue.removeFirst()
        
        if (current == destination) {
            // Reconstruct path
            val path = mutableListOf<Int>()
            var node = destination
            while (node != -1) {
                path.add(0, node)
                node = parent[node]
            }
            return path
        }
        
        // Check all neighbors
        for (next in 0..<n) {
            if (!visited[next] && adjacencyMatrix[current, next] > 0) {
                visited[next] = true
                parent[next] = current
                queue.add(next)
            }
        }
    }
    
    return emptyList() // No path found
}
/**
 * Cache gamma algorithms
 */
@Suppress("unused")
class CacheGammaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}