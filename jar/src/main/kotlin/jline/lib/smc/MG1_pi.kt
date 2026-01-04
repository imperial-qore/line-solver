/**
 * @file M/G/1-type Stationary Distribution
 *
 * Computes the stationary distribution for M/G/1-type Markov chains.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.util.matrix.Matrix

/**
 * Options for MG1_pi solver
 */
data class MG1PiOptions(
    val boundary: Matrix? = null,
    val maxNumComp: Int = 500,
    val precision: Int = 200,
    val solver: String = "FI",
    val verbose: Boolean = false,
    val mode: String = "ShiftPWCR"
)

/**
 * Computes the stationary distribution of an M/G/1-type Markov chain.
 *
 * The chain is characterized by:
 * - B: boundary blocks [B0 B1 B2 ... B_maxb] with m rows
 * - A: repeating blocks [A0 A1 A2 ... A_maxa] with m rows
 *
 * @param B Boundary block matrix (or null to use first row of A)
 * @param A Repeating block matrix
 * @param options Solver options
 * @return Stationary distribution vector
 */
fun mg1_pi(B: Matrix?, A: Matrix, options: MG1PiOptions = MG1PiOptions()): Matrix {
    val m = A.numRows
    val dega = A.numCols / m - 1

    // Use boundary or default to A's first row structure
    val boundary = B ?: A

    // First, compute G using appropriate solver
    val G = when (options.solver.uppercase()) {
        "CR" -> mg1_cr(A, MG1CROptions(
            mode = options.mode,
            verbose = if (options.verbose) 1 else 0
        ))
        else -> mg1_fi(A, MG1FIOptions(
            verbose = if (options.verbose) 1 else 0
        ))
    }

    // Compute R matrix
    // R = sum(i=0 to dega) A_i * G^i
    var R = A.extractCols(0, m).copy()
    var Gpow = G.copy()
    for (i in 1..dega) {
        R = R.add(A.extractCols(i * m, (i + 1) * m).mult(Gpow))
        Gpow = Gpow.mult(G)
    }

    // Compute stationary distribution
    // pi_0 satisfies: pi_0 * B_hat = 0 where B_hat = sum(B_i * G^i) - I
    val degb = boundary.numCols / m - 1

    var Bhat = boundary.extractCols(0, m).copy()
    Gpow = G.copy()
    for (i in 1..degb) {
        Bhat = Bhat.add(boundary.extractCols(i * m, (i + 1) * m).mult(Gpow))
        Gpow = Gpow.mult(G)
    }

    // Compute pi_0 as stationary distribution of Bhat
    val pi0 = stat(Bhat)

    // Compute higher levels using pi_i = pi_0 * R^i
    val result = mutableListOf<Matrix>()
    result.add(pi0)

    var Rpow = R.copy()
    for (i in 1 until options.maxNumComp) {
        val pi_i = pi0.mult(Rpow)

        // Check if probability mass is negligible
        val mass = pi_i.elementSum()
        if (mass < 1e-15) {
            break
        }

        result.add(pi_i)
        Rpow = Rpow.mult(R)
    }

    // Normalize the distribution
    var totalMass = 0.0
    for (pi_i in result) {
        totalMass += pi_i.elementSum()
    }

    // Return as single row vector
    val totalSize = result.size * m
    val piVec = Matrix(1, totalSize)
    for (i in result.indices) {
        for (j in 0 until m) {
            piVec[0, i * m + j] = result[i][0, j] / totalMass
        }
    }

    return piVec
}

/**
 * Computes the G matrix using Cyclic Reduction for M/G/1-type chains.
 * Placeholder - calls FI for now.
 */
fun mg1_cr(A: Matrix, options: MG1CROptions = MG1CROptions()): Matrix {
    // For now, delegate to FI
    // TODO: Implement full CR algorithm with FFT
    return mg1_fi(A, MG1FIOptions(
        mode = if (options.mode.contains("Shift")) "ShiftU-Based" else "U-Based",
        verbose = options.verbose
    ))
}

/**
 * Options for MG1_CR solver
 */
data class MG1CROptions(
    val mode: String = "ShiftPWCR",
    val maxNumIt: Int = 50,
    val maxNumRoot: Int = 2048,
    val epsilonValue: Double = 1e-16,
    val verbose: Int = 0,
    val shiftType: String = "one"
)
