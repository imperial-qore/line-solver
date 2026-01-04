/**
 * @file M/G/1-type Shift Technique
 *
 * Implements the shift technique for M/G/1-type Markov chains to improve
 * convergence of iterative solvers.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.util.matrix.Matrix

/**
 * Result of MG1 shift transformation
 */
data class MG1ShiftResult(
    val A: Matrix,
    val drift: Double,
    val tau: Double,
    val v: Matrix
)

/**
 * Applies the shift technique to an M/G/1-type block matrix.
 *
 * @param A Block matrix [A0 A1 A2 ... A_max]
 * @param shiftType Type of shift: "one", "tau", or "dbl"
 * @return Tuple of (shifted A, drift, tau, v)
 */
fun mg1_shifts(A: Matrix, shiftType: String = "one"): Quadruple<Matrix, Double, Double, Matrix> {
    val m = A.numRows
    val dega = A.numCols / m - 1

    // Compute sumA and beta for drift calculation
    var sumA = A.extractCols(dega * m, (dega + 1) * m)
    var beta = sumA.sumRows()

    for (i in dega - 1 downTo 1) {
        sumA = sumA.add(A.extractCols(i * m, (i + 1) * m))
        beta = beta.add(sumA.sumRows())
    }
    sumA = sumA.add(A.extractCols(0, m))

    // Compute stationary distribution theta
    val theta = stat(sumA)
    val drift = theta.mult(beta)[0, 0]

    // Compute v (left eigenvector) and tau
    var v = Matrix.ones(m, 1).scale(1.0 / m)
    var tau = 0.0

    if (drift > 1) {
        // Transient case: compute Neuts' tau and v
        // v is the left eigenvector of sum(i * A_i)
        var weightedSum = Matrix.zeros(m, m)
        for (i in 1..dega) {
            weightedSum = weightedSum.add(A.extractCols(i * m, (i + 1) * m).scale(i.toDouble()))
        }

        // Simple power iteration for left eigenvector
        v = Matrix.ones(1, m).scale(1.0 / m)
        for (iter in 0 until 100) {
            val vNew = v.mult(weightedSum)
            val norm = vNew.elementSum()
            if (norm > 0) {
                v = vNew.scale(1.0 / norm)
            }
        }
        v = v.transpose()

        // tau is the reciprocal of the spectral radius
        val spectralRadius = v.transpose().mult(weightedSum).mult(Matrix.ones(m, 1))[0, 0]
        if (spectralRadius > 0) {
            tau = 1.0 / spectralRadius
        }
    }

    // Apply shift transformation
    val shiftedA = Matrix(m, A.numCols)

    when (shiftType) {
        "one" -> {
            if (drift < 1) {
                // Shift with ones matrix
                val onesM = Matrix.ones(m, m).scale(1.0 / m)
                for (i in 0..dega) {
                    val Ai = A.extractCols(i * m, (i + 1) * m)
                    val shiftedAi = if (i == 0) {
                        Ai.sub(onesM)
                    } else {
                        Ai
                    }
                    for (r in 0 until m) {
                        for (c in 0 until m) {
                            shiftedA[r, i * m + c] = shiftedAi[r, c]
                        }
                    }
                }
            } else {
                // No shift needed
                for (r in 0 until m) {
                    for (c in 0 until A.numCols) {
                        shiftedA[r, c] = A[r, c]
                    }
                }
            }
        }

        "tau" -> {
            if (drift > 1) {
                // Shift with tau * v * ones
                val tauVones = v.mult(Matrix.ones(1, m)).scale(tau)
                for (i in 0..dega) {
                    val Ai = A.extractCols(i * m, (i + 1) * m)
                    val shiftedAi = if (i == dega) {
                        Ai.sub(tauVones)
                    } else {
                        Ai
                    }
                    for (r in 0 until m) {
                        for (c in 0 until m) {
                            shiftedA[r, i * m + c] = shiftedAi[r, c]
                        }
                    }
                }
            } else {
                for (r in 0 until m) {
                    for (c in 0 until A.numCols) {
                        shiftedA[r, c] = A[r, c]
                    }
                }
            }
        }

        "dbl" -> {
            // Apply both shifts if applicable
            val onesM = Matrix.ones(m, m).scale(1.0 / m)
            val tauVones = v.mult(Matrix.ones(1, m)).scale(tau)

            for (i in 0..dega) {
                var Ai = A.extractCols(i * m, (i + 1) * m)
                if (i == 0 && drift < 1) {
                    Ai = Ai.sub(onesM)
                }
                if (i == dega && drift > 1) {
                    Ai = Ai.sub(tauVones)
                }
                for (r in 0 until m) {
                    for (c in 0 until m) {
                        shiftedA[r, i * m + c] = Ai[r, c]
                    }
                }
            }
        }

        else -> {
            // No shift
            for (r in 0 until m) {
                for (c in 0 until A.numCols) {
                    shiftedA[r, c] = A[r, c]
                }
            }
        }
    }

    return Quadruple(shiftedA, drift, tau, v)
}

/**
 * Helper class for returning four values
 */
data class Quadruple<A, B, C, D>(
    val first: A,
    val second: B,
    val third: C,
    val fourth: D
)
