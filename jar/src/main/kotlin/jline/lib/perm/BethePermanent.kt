package jline.lib.perm

import jline.util.matrix.Matrix
import kotlin.math.*

/**
 * Implementation of Sum Product Algorithm (SPA) to approximate the Bethe permanent.
 * 
 * The Bethe permanent approximation uses message passing between matrix elements
 * to compute an approximation of the permanent. This method is particularly useful
 * for large matrices where exact computation is computationally expensive.
 * 
 * @param matrix The matrix for which to compute the Bethe permanent approximation
 * @param epsilon Convergence threshold for the SPA algorithm (default: 0.001)
 * @param maxIteration Maximum number of iterations if convergence is not reached (default: 200000)
 * @param solve Whether to automatically run solve() after construction (default: false)
 */
class BethePermanent(
    matrix: Matrix,
    private val epsilon: Double = 0.001,
    private val maxIteration: Int = 200000,
    solve: Boolean = false
) : PermSolver(matrix) {
    
    private val matrixSqrt: Matrix
    private val minValue = 2.220446049250314e-16 // Machine epsilon equivalent
    
    init {
        // Ensure all matrix elements are at least minValue to avoid numerical issues
        val data = Array(n) { i ->
            DoubleArray(n) { j ->
                maxOf(matrix.get(i, j), minValue)
            }
        }
        matrixSqrt = Matrix(data).apply {
            // Compute square root of each element
            for (i in 0 until numRows) {
                for (j in 0 until numCols) {
                    set(i, j, sqrt(get(i, j)))
                }
            }
        }
        
        if (solve) {
            solve()
        }
    }
    
    override fun compute() {
        value = spa()
    }
    
    /**
     * Sum Product Algorithm implementation for Bethe permanent approximation.
     * 
     * @return The Bethe permanent approximation
     */
    private fun spa(): Double {
        var rPast = Matrix.ones(n, n)
        var lPast = Matrix.ones(n, n)
        
        var (r, l) = update(lPast)
        
        var iteration = 0
        while (convergence(rPast, lPast, r, l) > epsilon && iteration < maxIteration) {
            iteration++
            rPast = r.copy()
            lPast = l.copy()
            
            // Message passing
            val (newR, newL) = update(lPast)
            r = newR
            l = newL
        }
        
        return bethe(l, r)
    }
    
    /**
     * Update the right going message and the left going message.
     * 
     * @param l Left going message from previous iteration
     * @return Pair of (right going message, left going message) for current iteration
     */
    private fun update(l: Matrix): Pair<Matrix, Matrix> {
        val r1 = Matrix.zeros(n, n)
        val l1 = Matrix.zeros(n, n)
        
        // Compute right going messages
        for (i in 0 until n) {
            var denomSum = 0.0
            for (j in 0 until n) {
                if (i != j) {
                    denomSum += matrixSqrt.get(i, j) * l.get(i, j)
                }
            }
            for (j in 0 until n) {
                r1.set(i, j, matrixSqrt.get(i, j) / denomSum)
            }
        }
        
        // Compute left going messages
        for (j in 0 until n) {
            var denomSum = 0.0
            for (i in 0 until n) {
                if (i != j) {
                    denomSum += matrixSqrt.get(i, j) * r1.get(i, j)
                }
            }
            for (i in 0 until n) {
                l1.set(i, j, matrixSqrt.get(i, j) / denomSum)
            }
        }
        
        return Pair(r1, l1)
    }
    
    /**
     * Calculate the squared difference between past and present messages to measure convergence.
     * 
     * @param r0 Right going message from previous iteration
     * @param l0 Left going message from previous iteration
     * @param r1 Right going message from current iteration
     * @param l1 Left going message from current iteration
     * @return Convergence measure (squared differences sum)
     */
    private fun convergence(r0: Matrix, l0: Matrix, r1: Matrix, l1: Matrix): Double {
        var sum = 0.0
        for (i in 0 until n) {
            for (j in 0 until n) {
                val rDiff = r0.get(i, j) - r1.get(i, j)
                val lDiff = l0.get(i, j) - l1.get(i, j)
                sum += rDiff * rDiff + lDiff * lDiff
            }
        }
        return sum
    }
    
    /**
     * Compute the Bethe permanent from the left and right going messages.
     * 
     * @param l Left going message
     * @param r Right going message
     * @return The Bethe permanent approximation
     */
    private fun bethe(l: Matrix, r: Matrix): Double {
        // Compute term1: sum over rows
        val term1 = DoubleArray(n) { i ->
            var sum = 0.0
            for (j in 0 until n) {
                sum += matrixSqrt.get(i, j) * l.get(i, j)
            }
            maxOf(sum, minValue)
        }
        
        // Compute term2: sum over columns
        val term2 = DoubleArray(n) { j ->
            var sum = 0.0
            for (i in 0 until n) {
                sum += matrixSqrt.get(i, j) * r.get(i, j)
            }
            maxOf(sum, minValue)
        }
        
        // Compute term3: element-wise product plus 1
        val term3 = Array(n) { i ->
            DoubleArray(n) { j ->
                maxOf(r.get(i, j) * l.get(i, j) + 1.0, minValue)
            }
        }
        
        // Compute log terms
        var lterm1 = 0.0
        for (i in 0 until n) {
            lterm1 -= ln(term1[i])
        }
        
        var lterm2 = 0.0
        for (j in 0 until n) {
            lterm2 -= ln(term2[j])
        }
        
        var lterm3 = 0.0
        for (i in 0 until n) {
            for (j in 0 until n) {
                lterm3 += ln(term3[i][j])
            }
        }
        
        val result = exp(-(lterm1 + lterm2 + lterm3))
        return if (result.isNaN() || result.isInfinite()) 0.0 else result
    }
}