package jline.lib.kpctoolbox.mc

import jline.util.matrix.Matrix
import org.apache.commons.math3.linear.LUDecomposition
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.util.FastMath
import java.util.Random

/**
 * Discrete-Time Markov Chain (DTMC) analysis functions.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/mc/dtmc_*.m
 */

/**
 * Normalizes a matrix to be a valid stochastic matrix.
 * Rescales rows to sum to 1. Rows with zero sum get 1 at diagonal.
 *
 * @param P Input matrix
 * @return Valid stochastic transition matrix
 */
fun dtmc_makestochastic(P: Matrix): Matrix {
    val n = P.numRows
    val result = Matrix(n, n)

    for (i in 0 until n) {
        // Compute row sum
        var rowSum = 0.0
        for (j in 0 until n) {
            rowSum += P.get(i, j)
        }

        if (rowSum > 0) {
            // Normalize row
            for (j in 0 until n) {
                result.set(i, j, P.get(i, j) / rowSum)
            }
            // Ensure diagonal adjustment for numerical stability
            var newRowSum = 0.0
            for (j in 0 until n) {
                if (j != i) {
                    newRowSum += result.get(i, j)
                }
            }
            result.set(i, i, minOf(maxOf(0.0, 1.0 - newRowSum), 1.0))
        } else {
            // Set absorbing state
            for (j in 0 until n) {
                result.set(i, j, 0.0)
            }
            result.set(i, i, 1.0)
        }
    }

    return result
}

/**
 * Checks feasibility of a stochastic matrix.
 * Verifies if row sums are close to 1 and elements are non-negative.
 *
 * @param P Matrix to check
 * @return Precision level (1-15) if feasible, or 0 if not feasible
 */
fun dtmc_isfeasible(P: Matrix): Int {
    val n = P.numRows

    // Compute row sums
    val rowSums = DoubleArray(n)
    var minElement = Double.MAX_VALUE

    for (i in 0 until n) {
        var sum = 0.0
        for (j in 0 until n) {
            val value = P.get(i, j)
            sum += value
            minElement = minOf(minElement, value)
        }
        rowSums[i] = sum
    }

    val minRowSum = rowSums.minOrNull() ?: 0.0
    val maxRowSum = rowSums.maxOrNull() ?: 0.0

    var result = 0
    for (tol in 1..15) {
        val tolerance = FastMath.pow(10.0, -tol.toDouble())
        if (minRowSum > 1 - tolerance && maxRowSum < 1 + tolerance && minElement > -tolerance) {
            result = tol
        }
    }

    return result
}

/**
 * Computes the equilibrium distribution of a discrete-time Markov chain.
 *
 * @param P Stochastic transition matrix
 * @return Equilibrium distribution vector
 */
fun dtmc_solve(P: Matrix): DoubleArray {
    val n = P.numRows

    // Convert DTMC to CTMC: Q = P - I
    val Q = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            Q.set(i, j, P.get(i, j))
        }
        Q.set(i, i, Q.get(i, i) - 1.0)
    }

    return ctmc_solve(Q)
}

/**
 * Generates a random stochastic transition matrix.
 *
 * @param n Size of the matrix (n x n)
 * @return Random stochastic transition matrix
 */
fun dtmc_rand(n: Int): Matrix {
    val (P, _) = ctmc_randomization(ctmc_rand(n))
    return P
}

/**
 * Simulates a trajectory of a discrete-time Markov chain.
 *
 * @param P Stochastic transition matrix
 * @param pi0 Initial probability distribution vector
 * @param nSteps Number of steps to simulate
 * @return Vector of state indices visited in the simulation (0-based)
 */
fun dtmc_simulate(P: Matrix, pi0: DoubleArray, nSteps: Int): IntArray {
    val random = Random()
    val states = IntArray(nSteps)
    val n = P.numRows

    // Sample initial state from pi0
    var rnd = random.nextDouble()
    var cumSum = 0.0
    var currentState = 0
    for (i in 0 until n) {
        cumSum += pi0[i]
        if (rnd <= cumSum) {
            currentState = i
            break
        }
    }

    // Precompute cumulative transition probabilities
    val cumP = Array(n) { i ->
        val row = DoubleArray(n)
        var sum = 0.0
        for (j in 0 until n) {
            sum += P.get(i, j)
            row[j] = sum
        }
        row
    }

    // Simulate trajectory
    for (step in 0 until nSteps) {
        states[step] = currentState

        // Check for absorbing state
        if (cumP[currentState][n - 1] == 0.0 || P.get(currentState, currentState) == 1.0) {
            // Fill remaining states with current state
            for (s in step + 1 until nSteps) {
                states[s] = currentState
            }
            break
        }

        // Sample next state
        rnd = random.nextDouble()
        for (j in 0 until n) {
            if (rnd <= cumP[currentState][j] && P.get(currentState, j) > 0) {
                currentState = j
                break
            }
        }
    }

    return states
}

/**
 * Computes the stochastic complement of a DTMC partition.
 *
 * @param P Stochastic transition matrix
 * @param I Indices of states to retain (0-based)
 * @return Stochastic complement matrix for the subset I
 */
fun dtmc_stochcomp(P: Matrix, I: IntArray? = null): Matrix {
    val n = P.numRows

    // Default: first half of states
    val indices = I ?: (0 until (n + 1) / 2).toList().toIntArray()

    // Complement of I
    val Ic = (0 until n).filter { !indices.contains(it) }.toIntArray()

    val m1 = indices.size
    val m2 = Ic.size

    // Extract submatrices
    // P11: transitions within I
    val P11 = Matrix(m1, m1)
    for ((ii, i) in indices.withIndex()) {
        for ((jj, j) in indices.withIndex()) {
            P11.set(ii, jj, P.get(i, j))
        }
    }

    // P12: transitions from I to Ic
    val P12 = Matrix(m1, m2)
    for ((ii, i) in indices.withIndex()) {
        for ((jj, j) in Ic.withIndex()) {
            P12.set(ii, jj, P.get(i, j))
        }
    }

    // P21: transitions from Ic to I
    val P21 = Matrix(m2, m1)
    for ((ii, i) in Ic.withIndex()) {
        for ((jj, j) in indices.withIndex()) {
            P21.set(ii, jj, P.get(i, j))
        }
    }

    // P22: transitions within Ic
    val P22 = Matrix(m2, m2)
    for ((ii, i) in Ic.withIndex()) {
        for ((jj, j) in Ic.withIndex()) {
            P22.set(ii, jj, P.get(i, j))
        }
    }

    // S = P11 + P12 * (I - P22)^{-1} * P21
    // MATLAB uses S2 \ P21 (linear solve), so we do the same
    val IminusP22 = Matrix(m2, m2)
    for (i in 0 until m2) {
        for (j in 0 until m2) {
            IminusP22.set(i, j, if (i == j) 1.0 - P22.get(i, j) else -P22.get(i, j))
        }
    }

    // Solve (I - P22) * X = P21 using LU decomposition (matches MATLAB S2 \ P21)
    val S2_real = MatrixUtils.createRealMatrix(m2, m2)
    for (i in 0 until m2) {
        for (j in 0 until m2) {
            S2_real.setEntry(i, j, IminusP22.get(i, j))
        }
    }

    val P21_real = MatrixUtils.createRealMatrix(m2, m1)
    for (i in 0 until m2) {
        for (j in 0 until m1) {
            P21_real.setEntry(i, j, P21.get(i, j))
        }
    }

    val solvedMatrix = try {
        LUDecomposition(S2_real).solver.solve(P21_real)
    } catch (e: Exception) {
        // Return P11 if solve fails
        return P11
    }

    // Compute P12 * solvedMatrix and add to P11
    // solvedMatrix is (m2 x m1) = (I-P22)^{-1} * P21
    val S = Matrix(m1, m1)
    for (i in 0 until m1) {
        for (j in 0 until m1) {
            var sum = P11.get(i, j)
            for (k in 0 until m2) {
                sum += P12.get(i, k) * solvedMatrix.getEntry(k, j)
            }
            S.set(i, j, sum)
        }
    }

    return S
}

/**
 * Computes the time-reversed transition matrix of a DTMC.
 *
 * @param P Stochastic transition matrix of the original process
 * @return Stochastic transition matrix of the time-reversed process
 */
fun dtmc_timereverse(P: Matrix): Matrix {
    val n = P.numRows
    val Prev = Matrix(n, n)
    val pie = dtmc_solve(P)

    for (i in 0 until n) {
        for (j in 0 until n) {
            if (pie[j] != 0.0) {
                Prev.set(j, i, P.get(i, j) * pie[i] / pie[j])
            }
        }
    }

    return Prev
}

/**
 * Computes transient probabilities for a DTMC using uniformization.
 *
 * @param pi0 Initial probability distribution vector
 * @param P Transition probability matrix
 * @param t Time point for transient analysis (default: 1e4)
 * @param tol Error tolerance (default: 1e-12)
 * @param maxiter Maximum iterations (default: 100)
 * @return Pair of (probability distribution at time t, number of iterations)
 */
fun dtmc_uniformization(
    pi0: DoubleArray,
    P: Matrix,
    t: Double = 1e4,
    tol: Double = 1e-12,
    maxiter: Int = 100
): Pair<DoubleArray, Int> {
    val n = P.numRows

    // Convert to CTMC generator: Q = P - I
    val Q = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            Q.set(i, j, P.get(i, j))
        }
        Q.set(i, i, Q.get(i, i) - 1.0)
    }

    return ctmc_uniformization(pi0, ctmc_makeinfgen(Q), t, tol, maxiter)
}
