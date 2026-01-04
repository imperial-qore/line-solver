/**
 * @file GI/M/1-type Stationary Distribution
 *
 * Computes the stationary distribution for GI/M/1-type Markov chains.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.io.line_warning
import jline.util.matrix.Matrix

/**
 * Options for GIM1_pi solver
 */
data class GIM1PiOptions(
    val boundary: Matrix? = null,
    val maxNumComp: Int = 500,
    val verbose: Int = 0
)

/**
 * Computes the stationary distribution of a GI/M/1-type Markov chain.
 *
 * The chain has a transition matrix of the form:
 *
 *     P = | B1  A0  0   0   0  ... |
 *         | B2  A1  A0  0   0  ... |
 *         | B3  A2  A1  A0  0  ... |
 *         | B4  A3  A2  A1  A0 ... |
 *         | ...                    |
 *
 * @param B Block matrix [B1; B2; ...; B_maxb] with maxb*m rows and m columns
 * @param R The minimal nonnegative solution to R = A0 + R*A1 + R^2*A2 + ... + R^max*A_max
 * @param options Solver options
 * @return Stationary distribution vector (1 x n)
 */
fun gim1_pi(B: Matrix, R: Matrix, options: GIM1PiOptions = GIM1PiOptions()): Matrix {
    val m = R.numRows

    // Check spectral radius of R
    val ImR = Matrix.eye(m).sub(R)
    val ImRinv = try {
        ImR.inv()
    } catch (e: Exception) {
        throw IllegalArgumentException("The spectral radius of R is not below 1: QBD is not positive recurrent")
    }

    // Verify (I-R)^{-1} is nonnegative
    for (i in 0 until m) {
        for (j in 0 until m) {
            if (ImRinv[i, j] < -100 * 1e-15) {
                throw IllegalArgumentException("The spectral radius of R is not below 1: QBD is not positive recurrent")
            }
        }
    }

    return if (options.boundary == null) {
        // Standard boundary case
        computePiStandard(B, R, ImRinv, m, options)
    } else {
        // General boundary case
        computePiGeneralBoundary(B, R, ImRinv, m, options)
    }
}

/**
 * Compute pi for standard boundary case
 */
private fun computePiStandard(
    B: Matrix,
    R: Matrix,
    ImRinv: Matrix,
    m: Int,
    options: GIM1PiOptions
): Matrix {
    val maxb = B.numRows / m

    // Compute BR = B1 + R*B2 + R^2*B3 + ... + R^(maxb-1)*B_maxb
    // Using Horner's method: BR = R*(...R*(R*B_maxb + B_{maxb-1}) + ... + B_2) + B_1
    var BR = B.extractRows((maxb - 1) * m, maxb * m)
    for (i in maxb - 1 downTo 1) {
        BR = R.mult(BR).add(B.extractRows((i - 1) * m, i * m))
    }

    // Compute pi_0 as stationary distribution of BR
    var pi0 = stat(BR)

    // Normalize: pi_0 = pi_0 / (pi_0 * (I-R)^{-1} * e)
    val e = Matrix.ones(m, 1)
    val normFactor = pi0.mult(ImRinv).mult(e)[0, 0]
    pi0 = pi0.scale(1.0 / normFactor)

    // Compute higher levels using pi_i = pi_{i-1} * R
    val piLevels = mutableListOf<Matrix>()
    piLevels.add(pi0)

    var sumPi = pi0.elementSum()
    var numit = 1

    while (sumPi < 1 - 1e-10 && numit < 1 + options.maxNumComp) {
        val piNext = piLevels.last().mult(R)
        piLevels.add(piNext)
        numit++
        sumPi += piNext.elementSum()

        if (options.verbose > 0 && numit % options.verbose == 0) {
            println("Accumulated mass after $numit iterations: $sumPi")
        }
    }

    if (numit == 1 + options.maxNumComp) {
        line_warning("GIM1_pi", "Maximum Number of Components %d reached", numit - 1)
    }

    // Reshape to row vector
    val totalCols = piLevels.size * m
    val result = Matrix(1, totalCols)
    for (i in piLevels.indices) {
        for (j in 0 until m) {
            result[0, i * m + j] = piLevels[i][0, j]
        }
    }

    return result
}

/**
 * Compute pi for general boundary case
 */
private fun computePiGeneralBoundary(
    B: Matrix,
    R: Matrix,
    ImRinv: Matrix,
    m: Int,
    options: GIM1PiOptions
): Matrix {
    val boundary = options.boundary!!
    val mb = B.numCols // number of states of boundary level
    val maxbm1 = (B.numRows - mb) / m // maxb - 1

    // Compute BR1 from B blocks (excluding first block)
    var BR1 = B.extractRows(mb + (maxbm1 - 1) * m, mb + maxbm1 * m)
    for (i in maxbm1 - 1 downTo 1) {
        BR1 = R.mult(BR1).add(B.extractRows(mb + (i - 1) * m, mb + i * m))
    }

    // Compute BR0 from boundary blocks
    val maxa = (boundary.numRows - mb) / m
    var BR0 = boundary.extractRows(mb + (maxa - 1) * m, mb + maxa * m)
    for (i in maxa - 1 downTo 1) {
        BR0 = R.mult(BR0).add(boundary.extractRows(mb + (i - 1) * m, mb + i * m))
    }

    // Construct combined matrix and compute pi0, pi1
    val B0 = B.extractRows(0, mb)
    val BoundaryTop = boundary.extractRows(0, mb)
    val combinedTop = Matrix.concatColumns(B0, BoundaryTop, null)
    val combinedBottom = Matrix.concatColumns(BR0, BR1, null)
    val combined = Matrix.concatRows(combinedTop, combinedBottom, null)

    var pi01 = stat(combined)

    // Normalize: pi0 * e + pi1 * (I-R)^{-1} * e = 1
    val e = Matrix.ones(m, 1)
    val emb = Matrix.ones(mb, 1)
    val pi0part = pi01.extractCols(0, mb)
    val pi1part = pi01.extractCols(mb, mb + m)
    val normFactor = pi0part.mult(emb)[0, 0] + pi1part.mult(ImRinv).mult(e)[0, 0]
    pi01 = pi01.scale(1.0 / normFactor)

    // Extract normalized pi0 and pi1
    val pi0 = pi01.extractCols(0, mb)
    var piCurrent = pi01.extractCols(mb, mb + m)

    // Compute higher levels
    val piLevels = mutableListOf<Matrix>()
    piLevels.add(piCurrent)

    var sumPi = pi0.elementSum() + piCurrent.elementSum()
    var numit = 1

    while (sumPi < 1 - 1e-10 && numit < options.maxNumComp) {
        val piNext = piLevels.last().mult(R)
        piLevels.add(piNext)
        numit++
        sumPi += piNext.elementSum()

        if (options.verbose > 0 && numit % options.verbose == 0) {
            println("Accumulated mass after $numit iterations: $sumPi")
        }
    }

    if (numit == options.maxNumComp) {
        line_warning("GIM1_pi", "Maximum Number of Components %d reached", numit)
    }

    // Reshape to row vector: [pi0 pi1 pi2 ...]
    val totalCols = mb + piLevels.size * m
    val result = Matrix(1, totalCols)

    // Copy pi0
    for (j in 0 until mb) {
        result[0, j] = pi0[0, j]
    }

    // Copy pi levels
    for (i in piLevels.indices) {
        for (j in 0 until m) {
            result[0, mb + i * m + j] = piLevels[i][0, j]
        }
    }

    return result
}
