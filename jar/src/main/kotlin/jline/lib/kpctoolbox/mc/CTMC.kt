package jline.lib.kpctoolbox.mc

import jline.lib.kpctoolbox.basic.e
import jline.lib.kpctoolbox.basic.eye
import jline.lib.kpctoolbox.basic.zeros
import jline.util.matrix.Matrix
import org.apache.commons.math3.linear.LUDecomposition
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.util.FastMath
import java.util.Random

/**
 * Continuous-Time Markov Chain (CTMC) analysis functions.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/mc/ctmc_*.m
 */

/**
 * Normalizes a matrix to be a valid infinitesimal generator.
 * Sets diagonal elements such that row sums are zero.
 *
 * @param Q Input matrix
 * @return Valid infinitesimal generator matrix where rows sum to zero
 */
fun ctmc_makeinfgen(Q: Matrix): Matrix {
    val n = Q.numRows
    val result = Matrix(n, n)

    // Copy off-diagonal elements and compute row sums
    for (i in 0 until n) {
        var rowSum = 0.0
        for (j in 0 until n) {
            if (i != j) {
                val value = Q.get(i, j)
                result.set(i, j, value)
                rowSum += value
            }
        }
        // Set diagonal to make row sum zero
        result.set(i, i, -rowSum)
    }

    return result
}

/**
 * Generates a random infinitesimal generator matrix.
 *
 * @param n Size of the matrix (n x n)
 * @return Random infinitesimal generator
 */
fun ctmc_rand(n: Int): Matrix {
    val random = Random()
    val Q = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            Q.set(i, j, random.nextDouble())
        }
    }
    return ctmc_makeinfgen(Q)
}

/**
 * Result of connected component analysis.
 */
data class ConnectedComponents(
    val numComponents: Int,
    val componentAssignment: IntArray
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false
        other as ConnectedComponents
        if (numComponents != other.numComponents) return false
        if (!componentAssignment.contentEquals(other.componentAssignment)) return false
        return true
    }

    override fun hashCode(): Int {
        var result = numComponents
        result = 31 * result + componentAssignment.contentHashCode()
        return result
    }
}

/**
 * Finds weakly connected components in a directed graph.
 *
 * @param G Adjacency matrix of the graph
 * @return ConnectedComponents with number of components and component assignment vector
 */
fun weaklyconncomp(G: Matrix): ConnectedComponents {
    val n = G.numRows

    // Make symmetric (undirected) graph for weakly connected components
    val adj = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            val gij = G.get(i, j)
            val gji = G.get(j, i)
            // Match MATLAB: G(isnan(G)) = 1; NaN is treated as connected
            if (gij.isNaN() || gij != 0.0 || gji.isNaN() || gji != 0.0) {
                adj.set(i, j, 1.0)
                adj.set(j, i, 1.0)
            }
        }
    }

    // Simple BFS/DFS-based connected component finding
    val visited = BooleanArray(n)
    val component = IntArray(n)
    var numComponents = 0

    for (start in 0 until n) {
        if (!visited[start]) {
            numComponents++
            // BFS from start
            val queue = ArrayDeque<Int>()
            queue.add(start)
            visited[start] = true
            component[start] = numComponents

            while (queue.isNotEmpty()) {
                val current = queue.removeFirst()
                for (neighbor in 0 until n) {
                    if (!visited[neighbor] && adj.get(current, neighbor) != 0.0) {
                        visited[neighbor] = true
                        component[neighbor] = numComponents
                        queue.add(neighbor)
                    }
                }
            }
        }
    }

    return ConnectedComponents(numComponents, component)
}

/**
 * Result of CTMC solving.
 */
data class CTMCSolveResult(
    val equilibriumDistribution: DoubleArray,
    val generator: Matrix,
    val numComponents: Int,
    val componentAssignment: IntArray
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false
        other as CTMCSolveResult
        if (!equilibriumDistribution.contentEquals(other.equilibriumDistribution)) return false
        if (numComponents != other.numComponents) return false
        if (!componentAssignment.contentEquals(other.componentAssignment)) return false
        return true
    }

    override fun hashCode(): Int {
        var result = equilibriumDistribution.contentHashCode()
        result = 31 * result + numComponents
        result = 31 * result + componentAssignment.contentHashCode()
        return result
    }
}

/**
 * Computes the equilibrium distribution of a continuous-time Markov chain.
 *
 * @param Q Infinitesimal generator matrix
 * @return Equilibrium distribution vector
 */
fun ctmc_solve(Q: Matrix): DoubleArray {
    return ctmc_solveFull(Q).equilibriumDistribution
}

/**
 * Computes the equilibrium distribution with full details.
 *
 * @param Q Infinitesimal generator matrix
 * @return CTMCSolveResult with equilibrium distribution and component info
 */
fun ctmc_solveFull(Q: Matrix): CTMCSolveResult {
    val n = Q.numRows

    // Handle trivial case
    if (n == 1) {
        return CTMCSolveResult(doubleArrayOf(1.0), Q, 1, intArrayOf(1))
    }

    // Normalize to valid infinitesimal generator
    val normalizedQ = ctmc_makeinfgen(Q)

    // Check for connected components
    val symMatrix = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            if (FastMath.abs(normalizedQ.get(i, j)) + FastMath.abs(normalizedQ.get(j, i)) > 0) {
                symMatrix.set(i, j, 1.0)
            }
        }
    }

    val cc = weaklyconncomp(symMatrix)

    if (cc.numComponents > 1) {
        // Reducible generator - solve each component recursively
        val p = DoubleArray(n)
        for (c in 1..cc.numComponents) {
            val indices = (0 until n).filter { cc.componentAssignment[it] == c }
            val m = indices.size
            val Qc = Matrix(m, m)
            for ((ii, i) in indices.withIndex()) {
                for ((jj, j) in indices.withIndex()) {
                    Qc.set(ii, jj, normalizedQ.get(i, j))
                }
            }
            val pc = ctmc_solve(ctmc_makeinfgen(Qc))
            for ((ii, i) in indices.withIndex()) {
                p[i] = pc[ii]
            }
        }

        // Normalize
        val sum = p.sum()
        if (sum > 0) {
            for (i in 0 until n) {
                p[i] /= sum
            }
        }

        return CTMCSolveResult(p, normalizedQ, cc.numComponents, cc.componentAssignment)
    }

    // Check if all zero
    var allZero = true
    for (i in 0 until n) {
        for (j in 0 until n) {
            if (normalizedQ.get(i, j) != 0.0) {
                allZero = false
                break
            }
        }
        if (!allZero) break
    }

    if (allZero) {
        val p = DoubleArray(n) { 1.0 / n }
        return CTMCSolveResult(p, normalizedQ, 1, IntArray(n) { 1 })
    }

    // Iterative stripping of states with zero rows AND columns (matching MATLAB lines 91-118)
    var nnzel = (0 until n).toList().toIntArray()
    var Qnnz = normalizedQ
    var bnnz = DoubleArray(n)
    var Qnnz_prev = Qnnz
    var bnnz_prev = bnnz.copyOf()

    var goon = true
    while (goon) {
        val m = Qnnz.numRows
        // Find indices where both column sum and row sum of abs values are nonzero
        val keep = mutableListOf<Int>()
        for (idx in 0 until m) {
            var colSum = 0.0
            var rowSum = 0.0
            for (k in 0 until m) {
                colSum += FastMath.abs(Qnnz.get(k, idx))
                rowSum += FastMath.abs(Qnnz.get(idx, k))
            }
            if (colSum != 0.0 && rowSum != 0.0) {
                keep.add(idx)
            }
        }

        // Extract sub-matrix for kept indices
        val mKeep = keep.size
        val QnnzNew = Matrix(mKeep, mKeep)
        val bnnzNew = DoubleArray(mKeep)
        for ((ii, i) in keep.withIndex()) {
            bnnzNew[ii] = bnnz[i]
            for ((jj, j) in keep.withIndex()) {
                QnnzNew.set(ii, jj, Qnnz.get(i, j))
            }
        }

        val QnnzNorm = ctmc_makeinfgen(QnnzNew)

        // Map kept local indices back to original indices
        val newNnzel = IntArray(mKeep) { nnzel[keep[it]] }

        // Check convergence: sizes must match
        if (Qnnz_prev.numRows == QnnzNorm.numRows && bnnz_prev.size == bnnzNew.size) {
            goon = false
        } else {
            Qnnz_prev = QnnzNorm
            bnnz_prev = bnnzNew
            nnzel = newNnzel
        }

        Qnnz = QnnzNorm
        bnnz = bnnzNew
        nnzel = newNnzel
    }

    // If all states were stripped, return uniform
    if (Qnnz.numRows == 0) {
        val p = DoubleArray(n) { 1.0 / n }
        return CTMCSolveResult(p, normalizedQ, 1, IntArray(n) { 1 })
    }

    val mSolve = Qnnz.numRows

    // Solve the system: p * Q = 0, sum(p) = 1
    // Modified to: Q' * p' = b where b = [0,0,...,1] and last column of Q' is replaced with ones

    val Qmod = Matrix(mSolve, mSolve)
    for (i in 0 until mSolve) {
        for (j in 0 until mSolve) {
            Qmod.set(j, i, Qnnz.get(i, j)) // Transpose
        }
    }

    // Replace last column with ones for normalization constraint
    for (i in 0 until mSolve) {
        Qmod.set(i, mSolve - 1, 1.0)
    }

    // Build RHS vector
    val b = DoubleArray(mSolve)
    b[mSolve - 1] = 1.0

    // Solve using LU decomposition
    val realMatrix = MatrixUtils.createRealMatrix(mSolve, mSolve)
    for (i in 0 until mSolve) {
        for (j in 0 until mSolve) {
            realMatrix.setEntry(i, j, Qmod.get(i, j))
        }
    }

    val pSolve = try {
        val solver = LUDecomposition(realMatrix).solver
        solver.solve(MatrixUtils.createRealVector(b)).toArray()
    } catch (e: Exception) {
        // If singular, return NaN so callers can detect and fall back to dtmc_solve_reducible
        DoubleArray(mSolve) { Double.NaN }
    }

    // Map solution back to full state vector (zeros for stripped states)
    val p = DoubleArray(n)
    for ((ii, origIdx) in nnzel.withIndex()) {
        p[origIdx] = pSolve[ii]
    }

    return CTMCSolveResult(p, normalizedQ, cc.numComponents, cc.componentAssignment)
}

/**
 * Computes the time-reversed generator of a CTMC.
 *
 * @param Q Infinitesimal generator matrix
 * @return Infinitesimal generator matrix of the time-reversed process
 */
fun ctmc_timereverse(Q: Matrix): Matrix {
    val n = Q.numRows
    val Qrev = Matrix(n, n)
    val pie = ctmc_solve(Q)

    for (i in 0 until n) {
        for (j in 0 until n) {
            if (pie[j] != 0.0) {
                Qrev.set(j, i, Q.get(i, j) * pie[i] / pie[j])
            }
        }
    }

    // Ensure valid generator with zero row sums (safety for numerical noise)
    return ctmc_makeinfgen(Qrev)
}

/**
 * Applies uniformization (randomization) to transform a CTMC into a DTMC.
 *
 * @param Q Infinitesimal generator matrix
 * @param q Uniformization rate (if null, uses max|diag(Q)| + random)
 * @return Pair of (uniformized stochastic matrix, uniformization rate)
 */
fun ctmc_randomization(Q: Matrix, q: Double? = null): Pair<Matrix, Double> {
    val n = Q.numRows

    // Compute uniformization rate
    val rate = q ?: (maxAbsDiagonal(Q) + Random().nextDouble())

    // P = Q/q + I
    val P = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            P.set(i, j, Q.get(i, j) / rate)
        }
        P.set(i, i, P.get(i, i) + 1.0)
    }

    return Pair(dtmc_makestochastic(P), rate)
}

/**
 * Computes transient probabilities using uniformization method.
 *
 * @param pi0 Initial probability distribution vector
 * @param Q Infinitesimal generator matrix
 * @param t Time point for transient analysis
 * @param tol Error tolerance (default: 1e-12)
 * @param maxiter Maximum iterations (default: 100)
 * @return Pair of (probability distribution at time t, number of iterations)
 */
fun ctmc_uniformization(
    pi0: DoubleArray,
    Q: Matrix,
    t: Double,
    tol: Double = 1e-12,
    maxiter: Int = 100
): Pair<DoubleArray, Int> {
    val n = Q.numRows

    // Uniformization rate
    var maxDiag = 0.0
    for (i in 0 until n) {
        maxDiag = maxOf(maxDiag, FastMath.abs(Q.get(i, i)))
    }
    val q = 1.1 * maxDiag

    // Uniformized matrix: Qs = I + Q/q
    val Qs = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            Qs.set(i, j, Q.get(i, j) / q)
        }
        Qs.set(i, i, Qs.get(i, i) + 1.0)
    }

    // Find number of terms needed
    var k = 0
    var s = 1.0
    var r = 1.0
    var kmax = 1

    for (iter in 0 until maxiter) {
        k++
        r = r * (q * t) / k
        s += r
        if (1 - FastMath.exp(-q * t) * s <= tol) {
            kmax = k
            break
        }
    }

    // Compute transient probability
    var pi = pi0.map { it * FastMath.exp(-q * t) }.toDoubleArray()
    var P = pi0.copyOf()
    var ri = FastMath.exp(-q * t)

    for (j in 1..kmax) {
        // P = P * Qs
        val newP = DoubleArray(n)
        for (i in 0 until n) {
            var sum = 0.0
            for (l in 0 until n) {
                sum += P[l] * Qs.get(l, i)
            }
            newP[i] = sum
        }
        P = newP

        ri *= (q * t / j)
        for (i in 0 until n) {
            pi[i] += ri * P[i]
        }
    }

    return Pair(pi, kmax)
}

/**
 * Helper function to find maximum absolute diagonal value.
 */
private fun maxAbsDiagonal(M: Matrix): Double {
    val n = M.numRows
    var maxVal = 0.0
    for (i in 0 until n) {
        maxVal = maxOf(maxVal, FastMath.abs(M.get(i, i)))
    }
    return maxVal
}

/**
 * Computes the equilibrium distribution relative to a reference state.
 * Solves global balance equations with p(refstate) = 1 normalization.
 *
 * @param Q Infinitesimal generator matrix
 * @param refstate Index of reference state (0-based, default: 0)
 * @return Relative equilibrium distribution vector
 */
fun ctmc_relsolve(Q: Matrix, refstate: Int = 0): DoubleArray {
    val n = Q.numRows

    if (n == 1) {
        return doubleArrayOf(1.0)
    }

    val normalizedQ = ctmc_makeinfgen(Q)

    // Build Q^T (transpose of normalized generator)
    val Qmod = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            Qmod.set(j, i, normalizedQ.get(i, j)) // Transpose
        }
    }

    // Replace last row of Q^T (= last column of Q) with normalization constraint
    // MATLAB: Qnnz(:,end) = 0; Qnnz(refstate,end) = 1; bnnz(end) = 1; then solves Qnnz' \ bnnz
    val lastRow = n - 1
    for (j in 0 until n) {
        Qmod.set(lastRow, j, 0.0)
    }
    Qmod.set(lastRow, refstate, 1.0)

    // Build RHS vector
    val b = DoubleArray(n)
    b[lastRow] = 1.0

    // Solve using LU decomposition
    val realMatrix = MatrixUtils.createRealMatrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            realMatrix.setEntry(i, j, Qmod.get(i, j))
        }
    }

    return try {
        val solver = LUDecomposition(realMatrix).solver
        solver.solve(MatrixUtils.createRealVector(b)).toArray()
    } catch (e: Exception) {
        DoubleArray(n) { 1.0 / n }
    }
}
