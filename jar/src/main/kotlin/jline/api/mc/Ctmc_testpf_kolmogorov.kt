package jline.api.mc

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Test if a CTMC has product form using Kolmogorov's criteria
 * Checks if the product of transition rates along any cycle equals
 * the product of rates in the reverse direction
 *
 * @param Q infinitesimal generator matrix
 * @return true if the CTMC has product form, false otherwise
 */
fun ctmc_testpf_kolmogorov(Q: Matrix): Boolean {
    // Make Q stochastic (normalize if needed)
    val Q_norm = ctmc_makeStochastic(Q)
    
    // Solve for stationary distribution
    val pi = ctmc_solve(Q_norm)
    
    // Compute time-reversed CTMC
    val Qr = ctmc_timereverse(Q_norm)
    
    val n = Q.length()
    
    // Create adjacency matrix from Q (excluding diagonal)
    val A = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            if (i != j && Q_norm.get(i, j) > 0) {
                A.set(i, j, 1.0)
            }
        }
    }
    
    // Find all cycles and check Kolmogorov's criteria
    println("Finding cycles")
    
    for (start in 0 until n) {
        for (target in 0 until n) {
            if (target != start && A.get(target, start) > 0) {
                // Find all paths from start to target
                val unusedNodes = BooleanArray(n) { true }
                val emptyPath = mutableListOf<Int>()
                val cycles = findPaths(A, unusedNodes, emptyPath, start, target)
                
                // Check each cycle
                for (cycle in cycles) {
                    val completeCycle = cycle + start
                    val reverseCycle = completeCycle.reversed()
                    
                    // Calculate product of rates in forward direction
                    var q = 1.0
                    for (i in 0 until completeCycle.size - 1) {
                        q *= Q_norm.get(completeCycle[i], completeCycle[i + 1])
                    }
                    
                    // Calculate product of rates in reverse direction
                    var qr = 1.0
                    for (i in 0 until reverseCycle.size - 1) {
                        qr *= Qr.get(reverseCycle[i], reverseCycle[i + 1])
                    }
                    
                    println("Forward: $q, Reverse: $qr")
                    
                    // Check if products are equal (within tolerance)
                    if (abs(q - qr) / q > 1e-6) {
                        return false
                    }
                }
            }
        }
    }
    
    return true
}

/**
 * Helper function to make a CTMC generator stochastic
 * Ensures off-diagonal elements are non-negative and rows sum to zero
 *
 * @param Q infinitesimal generator matrix
 * @return normalized infinitesimal generator
 */
private fun ctmc_makeStochastic(Q: Matrix): Matrix {
    val n = Q.length()
    val result = Q.copy()
    
    // Ensure non-negative off-diagonal elements
    for (i in 0 until n) {
        for (j in 0 until n) {
            if (i != j && result.get(i, j) < 0) {
                result.set(i, j, 0.0)
            }
        }
    }
    
    // Adjust diagonal to ensure row sums are zero
    for (i in 0 until n) {
        var rowSum = 0.0
        for (j in 0 until n) {
            if (i != j) {
                rowSum += result.get(i, j)
            }
        }
        result.set(i, i, -rowSum)
    }
    
    return result
}

/**
 * Recursive function to find all paths from start to target
 *
 * @param adj adjacency matrix
 * @param nodes array tracking which nodes are still unused
 * @param currentPath current path being explored
 * @param start current node
 * @param target target node
 * @return list of paths from start to target
 */
private fun findPaths(
    adj: Matrix,
    nodes: BooleanArray,
    currentPath: MutableList<Int>,
    start: Int,
    target: Int
): List<List<Int>> {
    val paths = mutableListOf<List<Int>>()
    
    // Mark current node as used
    val newNodes = nodes.clone()
    newNodes[start] = false
    
    // Add current node to path
    val newPath = currentPath.toMutableList()
    newPath.add(start)
    
    // If we reached the target, add this path
    if (start == target) {
        paths.add(newPath)
        return paths
    }
    
    // Find all unvisited children
    val childList = mutableListOf<Int>()
    for (j in 0 until adj.getNumCols()) {
        if (adj.get(start, j) > 0 && newNodes[j]) {
            childList.add(j)
        }
    }
    
    // If no children or reached target, return
    if (childList.isEmpty()) {
        return paths
    }
    
    // Explore each child
    for (child in childList) {
        val childPaths = findPaths(adj, newNodes, newPath, child, target)
        paths.addAll(childPaths)
    }
    
    return paths
}
/**
 * CTMC testpf kolmogorov algorithms
 */
@Suppress("unused")
class CtmcTestpfKolmogorovAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}