package jline.lib.perm

import jline.util.matrix.Matrix
import kotlin.math.*

/**
 * Heuristic approximation to the permanent of a positive matrix.
 *
 * This implementation uses Sinkhorn scaling to make the matrix approximately
 * doubly stochastic, then applies a mean-field approximation with van der Waerden
 * bounds combined with a Gurvits-like capacity bound.
 *
 * The algorithm:
 * 1. Uses Sinkhorn scaling to make A approximately doubly stochastic
 * 2. Applies a mean-field approximation with van der Waerden bounds
 * 3. Combines with a Gurvits-like capacity bound
 * 4. Scales the result back to the original matrix scale
 *
 * This is a heuristic approximation suitable for large matrices where exact
 * computation is computationally prohibitive. The approximation quality depends
 * on the structure of the input matrix.
 *
 * @param matrix The matrix for which to compute the permanent (must be strictly positive)
 * @param tolerance Convergence threshold for Sinkhorn scaling
 * @param maxIterations Maximum number of Sinkhorn iterations
 * @param solve Whether to automatically run solve() after construction
 * @throws IllegalArgumentException if matrix contains non-positive elements
 */
class HeuristicPermanent(
    matrix: Matrix,
    private val tolerance: Double,
    private val maxIterations: Int,
    solve: Boolean
) : PermSolver(matrix) {

    /**
     * Constructor with just matrix and solve parameters.
     * Uses default tolerance (1e-10) and maxIterations (1000).
     */
    constructor(matrix: Matrix, solve: Boolean) : this(matrix, 1e-10, 1000, solve)

    /**
     * Constructor with just matrix parameter.
     * Uses default tolerance (1e-10), maxIterations (1000), and solve (false).
     */
    constructor(matrix: Matrix) : this(matrix, 1e-10, 1000, false)

    init {
        // Validate that all matrix elements are non-negative
        // We'll add a small epsilon to zeros during computation
        for (i in 0 until matrix.getNumRows()) {
            for (j in 0 until matrix.getNumCols()) {
                if (matrix.get(i, j) < 0.0) {
                    throw IllegalArgumentException("Matrix must be non-negative. Found negative element at ($i, $j): ${matrix.get(i, j)}")
                }
            }
        }

        if (solve) {
            solve()
        }
    }

    override fun compute() {
        value = computeHeuristicPermanent()
    }

    /**
     * Computes the heuristic permanent approximation using Sinkhorn scaling
     * and mean-field approximation.
     *
     * @return The approximate permanent value
     */
    private fun computeHeuristicPermanent(): Double {
        // Create a working copy with small epsilon added to zeros to ensure positivity
        val epsilon = 1e-15
        val workingMatrix = matrix.copy()
        for (i in 0 until n) {
            for (j in 0 until n) {
                if (workingMatrix.get(i, j) == 0.0) {
                    workingMatrix.set(i, j, epsilon)
                }
            }
        }

        // Sinkhorn scaling to make matrix approximately doubly stochastic
        val (B, r, c) = sinkhornScaling(workingMatrix)

        // Compute approximate permanent of doubly stochastic matrix B
        // Mean-field approximation: product of row sums / n^n * n!
        val rowSums = DoubleArray(n) { i ->
            var sum = 0.0
            for (j in 0 until n) {
                sum += B.get(i, j)
            }
            sum
        }

        var rowProd = 1.0
        for (i in 0 until n) {
            rowProd *= rowSums[i]
        }

        val pMeanfield = factorial(n) * (rowProd / n.toDouble().pow(n))

        // Gurvits-like capacity bound (optional refinement)
        var logSumRowSums = 0.0
        for (i in 0 until n) {
            logSumRowSums += ln(rowSums[i])
        }
        val cap = exp(logSumRowSums / n)
        val pGurvits = factorial(n) * (cap / n).pow(n)

        // Combine (simple average)
        var pEst = 0.5 * (pMeanfield + pGurvits)

        // Undo scaling
        var scaleFactor = 1.0
        for (i in 0 until n) {
            scaleFactor *= (1.0 / r[i])
        }
        for (j in 0 until n) {
            scaleFactor *= (1.0 / c[j])
        }

        pEst *= scaleFactor

        return pEst
    }

    /**
     * Performs Sinkhorn scaling to make the matrix approximately doubly stochastic.
     *
     * @param inputMatrix The matrix to scale
     * @return Triple of (scaled matrix B, row scaling factors r, column scaling factors c)
     */
    private fun sinkhornScaling(inputMatrix: Matrix): Triple<Matrix, DoubleArray, DoubleArray> {
        val B = inputMatrix.copy()
        val r = DoubleArray(n) { 1.0 }
        val c = DoubleArray(n) { 1.0 }

        for (iter in 0 until maxIterations) {
            // Update row scaling: r = 1 / (B * c)
            for (i in 0 until n) {
                var sum = 0.0
                for (j in 0 until n) {
                    sum += B.get(i, j) * c[j]
                }
                r[i] = 1.0 / sum
            }

            // Update column scaling: c = 1 / (B' * r)
            for (j in 0 until n) {
                var sum = 0.0
                for (i in 0 until n) {
                    sum += B.get(i, j) * r[i]
                }
                c[j] = 1.0 / sum
            }

            // Check convergence: max(abs(r .* (B * c) - 1)) < tol
            var maxDiff = 0.0
            for (i in 0 until n) {
                var rowSum = 0.0
                for (j in 0 until n) {
                    rowSum += B.get(i, j) * c[j]
                }
                val diff = abs(r[i] * rowSum - 1.0)
                if (diff > maxDiff) {
                    maxDiff = diff
                }
            }

            if (maxDiff < tolerance) {
                break
            }
        }

        // Apply scaling to matrix: B = diag(r) * B * diag(c)
        val scaledB = Matrix(n, n)
        for (i in 0 until n) {
            for (j in 0 until n) {
                scaledB.set(i, j, r[i] * B.get(i, j) * c[j])
            }
        }

        return Triple(scaledB, r, c)
    }

    /**
     * Computes n! (factorial of n).
     * For large n, this uses Stirling's approximation to avoid overflow.
     *
     * @param n The value to compute factorial for
     * @return n! as a double
     */
    private fun factorial(n: Int): Double {
        if (n <= 1) return 1.0
        if (n <= 20) {
            // For small n, compute exactly
            var result = 1.0
            for (i in 2..n) {
                result *= i
            }
            return result
        } else {
            // For large n, use Stirling's approximation: n! ≈ sqrt(2*pi*n) * (n/e)^n
            return sqrt(2.0 * PI * n) * (n / E).pow(n)
        }
    }
}
