/**
 * @file M/G/1-type Functional Iteration solver
 *
 * Computes the minimal nonnegative solution to M/G/1-type matrix equations
 * using Functional Iteration methods (Neuts).
 *
 * Solves: G = A0 + A1*G + A2*G^2 + ... + A_max*G^max
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.io.line_warning
import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Options for MG1_FI solver
 */
data class MG1FIOptions(
    val mode: String = "U-Based",
    val maxNumIt: Int = 10000,
    val verbose: Int = 0,
    val shiftType: String = "one",
    val startValue: Matrix? = null,
    val nonZeroBlocks: IntArray? = null
)

/**
 * Functional Iterations for M/G/1-Type Markov Chains.
 *
 * Computes the minimal nonnegative solution to the matrix equation
 * G = A0 + A1*G + A2*G^2 + ... + A_max*G^max
 *
 * @param A Block matrix [A0 A1 A2 ... A_max] with m rows and m*(max+1) columns
 * @param options Solver options
 * @return G matrix (minimal nonnegative solution)
 */
fun mg1_fi(A: Matrix, options: MG1FIOptions = MG1FIOptions()): Matrix {
    val m = A.numRows
    val maxd = A.numCols / m - 1

    // Check whether G is known explicitly
    val explicitG = mg1_eg(A, options.verbose > 0)
    if (explicitG != null) {
        return explicitG
    }

    var numit = 0
    var check = 1.0
    var G = options.startValue ?: Matrix.zeros(m, m)

    // Working copy of A (may be shifted)
    var workA = A.copy()
    var drift = 0.0
    var tau = 0.0
    var v: Matrix? = null

    val useShift = options.mode.startsWith("Shift")
    if (useShift && options.nonZeroBlocks == null) {
        val shiftResult = mg1_shifts(workA, options.shiftType)
        workA = shiftResult.first
        drift = shiftResult.second
        tau = shiftResult.third
        v = shiftResult.fourth
    }

    val actualMode = if (useShift) options.mode.removePrefix("Shift") else options.mode

    if (options.nonZeroBlocks == null) {
        // Standard iteration without NonZeroBlocks
        when {
            actualMode.contains("Natural") -> {
                while (check > 1e-14 && numit < options.maxNumIt) {
                    val Gold = G.copy()
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy()
                    for (j in maxd - 1 downTo 0) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold))
                    }
                    check = G.sub(Gold).infinityNorm()
                    numit++
                }
            }

            actualMode.contains("Traditional") -> {
                while (check > 1e-14 && numit < options.maxNumIt) {
                    val Gold = G.copy()
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy()
                    for (j in maxd - 1 downTo 2) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold))
                    }
                    G = workA.extractCols(0, m).add(G.mult(Gold).mult(Gold))
                    val ImA1 = Matrix.eye(m).sub(workA.extractCols(m, 2 * m))
                    G = ImA1.inv().mult(G)
                    check = G.sub(Gold).infinityNorm()
                    numit++
                }
            }

            actualMode.contains("U-Based") -> {
                while (check > 1e-14 && numit < options.maxNumIt) {
                    val Gold = G.copy()
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy()
                    for (j in maxd - 1 downTo 1) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold))
                    }
                    val ImG = Matrix.eye(m).sub(G)
                    G = ImG.inv().mult(workA.extractCols(0, m))
                    check = G.sub(Gold).infinityNorm()
                    numit++
                }
            }
        }

        // Apply shift correction
        if (useShift && v != null) {
            when (options.shiftType) {
                "one" -> if (drift < 1) G = G.add(Matrix.ones(m, m).scale(1.0 / m))
                "tau" -> if (drift > 1) G = G.add(v.mult(Matrix.ones(1, m)).scale(tau))
                "dbl" -> {
                    if (drift < 1) G = G.add(Matrix.ones(m, m).scale(1.0 / m))
                    if (drift > 1) G = G.add(v.mult(Matrix.ones(1, m)).scale(tau))
                }
            }
        }
    } else {
        // NonZeroBlocks mode
        val vec = intArrayOf(0) + options.nonZeroBlocks
        val vecDiff = IntArray(vec.size - 1) { vec[it + 1] - vec[it] }

        if (vecDiff[0] > 1 && actualMode == "Traditional") {
            // If A1 = 0, Traditional = Natural
            return mg1_fi(A, options.copy(mode = "Natural"))
        }

        when {
            actualMode.contains("Natural") -> {
                while (check > 1e-14 && numit < options.maxNumIt) {
                    val Gold = G.copy()
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy()
                    for (j in maxd - 1 downTo 0) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold.pow(vecDiff[j + 1])))
                    }
                    check = G.sub(Gold).infinityNorm()
                    numit++
                }
            }

            actualMode.contains("Traditional") -> {
                while (check > 1e-14 && numit < options.maxNumIt) {
                    val Gold = G.copy()
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy()
                    for (j in maxd - 1 downTo 2) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold.pow(vecDiff[j + 1])))
                    }
                    G = workA.extractCols(0, m).add(G.mult(Gold.pow(1 + vecDiff[2])))
                    val ImA1 = Matrix.eye(m).sub(workA.extractCols(m, 2 * m))
                    G = ImA1.inv().mult(G)
                    check = G.sub(Gold).infinityNorm()
                    numit++
                }
            }

            actualMode.contains("U-Based") -> {
                while (check > 1e-14 && numit < options.maxNumIt) {
                    val Gold = G.copy()
                    G = workA.extractCols(maxd * m, (maxd + 1) * m).copy()
                    for (j in maxd - 1 downTo 1) {
                        G = workA.extractCols(j * m, (j + 1) * m).add(G.mult(Gold.pow(vecDiff[j + 1])))
                    }
                    val ImGpow = Matrix.eye(m).sub(G.mult(Gold.pow(vecDiff[0] - 1)))
                    G = ImGpow.inv().mult(workA.extractCols(0, m))
                    check = G.sub(Gold).infinityNorm()
                    numit++
                }
            }
        }
    }

    if (numit == options.maxNumIt) {
        line_warning("MG1_FI", "Maximum number of iterations %d reached", numit)
    }

    return G
}

/**
 * Compute matrix power
 */
private fun Matrix.pow(n: Int): Matrix {
    if (n == 0) return Matrix.eye(this.numRows)
    if (n == 1) return this.copy()

    var result = Matrix.eye(this.numRows)
    var base = this.copy()
    var exp = n

    while (exp > 0) {
        if (exp % 2 == 1) {
            result = result.mult(base)
        }
        base = base.mult(base)
        exp /= 2
    }
    return result
}
