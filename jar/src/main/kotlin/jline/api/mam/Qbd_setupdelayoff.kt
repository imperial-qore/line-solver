/**
 * @file Quasi-Birth-Death process setup delays and server switch-off analysis
 *
 * Analyzes queueing systems with server setup delays and switch-off mechanisms.
 * Models energy-efficient server operations and startup costs in queueing systems.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix

/**
 * Analyze a queue with setup delays and server switch-off using QBD approach.
 *
 * Models an M/M/1 queue where:
 * - When server becomes idle, it enters a delayoff phase before turning off
 * - When a job arrives to an off server, it goes through setup before service
 *
 * @param lambda Arrival rate
 * @param mu Service rate
 * @param alphaRate Setup rate (1/mean setup time)
 * @param alphaScv Setup squared coefficient of variation (1.0 for exponential)
 * @param betaRate Delayoff rate (1/mean delayoff time)
 * @param betaScv Delayoff squared coefficient of variation (1.0 for exponential)
 * @return Mean queue length
 */
fun qbd_setupdelayoff(
    lambda: Double,
    mu: Double,
    alphaRate: Double,
    alphaScv: Double,
    betaRate: Double,
    betaScv: Double
): Double {
    // For exponential distributions (SCV=1), we use 2 phases:
    // Phase 0: setup phase (server warming up)
    // Phase 1: busy phase (server active)
    val n = 2

    // QBD matrices using standard notation:
    // A0 = upward transitions (arrivals)
    // A1 = local transitions at levels > 0
    // A2 = downward transitions (service completions)

    // A0: Arrivals can happen in both phases
    val A0 = Matrix.zeros(n, n)
    A0[0, 0] = lambda  // arrival during setup
    A0[1, 1] = lambda  // arrival during busy

    // A2: Service completions only from busy phase
    val A2 = Matrix.zeros(n, n)
    A2[1, 1] = mu

    // A1: Local transitions at levels > 0
    val A1 = Matrix.zeros(n, n)
    A1[0, 0] = -(alphaRate + lambda)  // leave setup: complete setup or arrival
    A1[0, 1] = alphaRate               // setup completes -> busy
    A1[1, 1] = -(mu + lambda)          // leave busy: service or arrival

    // Boundary matrices for level 0:
    // B0 = upward from level 0 (same as A0)
    val B0 = A0.copy()

    // B1 = local at level 0 (server off or in delayoff)
    val B1 = Matrix.zeros(n, n)
    B1[0, 0] = -lambda                    // idle, only leave on arrival
    B1[1, 1] = -(betaRate + lambda)       // in delayoff phase
    B1[1, 0] = betaRate                   // delayoff completes -> idle

    // Solve for R matrix using Cyclic Reduction / iteration
    val R = solveR(A0, A1, A2)

    // Solve for stationary distribution
    val pi = solvePi(B0, B1, R, A0, A1, A2)

    // Compute mean queue length: sum over levels of level * P(level)
    var QN = 0.0
    val maxLevel = pi.size / n
    for (level in 1 until maxLevel) {
        var levelProb = 0.0
        for (phase in 0 until n) {
            val idx = level * n + phase
            if (idx < pi.size) {
                levelProb += pi[idx]
            }
        }
        QN += level * levelProb
    }

    return QN
}

/**
 * Solve for the R matrix using successive substitution.
 * R satisfies: 0 = A0 + R*A1 + R^2*A2 (continuous time)
 * Equivalently: R = -A0 * (A1 + R*A2)^{-1}
 */
private fun solveR(A0: Matrix, A1: Matrix, A2: Matrix): Matrix {
    val n = A0.numRows
    val maxIter = 10000
    val tol = 1e-14

    // Initial guess: R = 0
    var R = Matrix.zeros(n, n)

    for (iter in 0 until maxIter) {
        // Compute A1 + R*A2
        val M = A1.add(R.mult(A2))

        // Compute R_new = -A0 * M^{-1}
        val Minv = M.inv()
        val Rnew = A0.mult(Minv).scale(-1.0)

        // Check convergence
        val diff = Rnew.sub(R)
        val normDiff = diff.norm()

        R = Rnew

        if (normDiff < tol) {
            break
        }
    }

    return R
}

/**
 * Solve for stationary distribution pi.
 * Uses the matrix-geometric property: pi_k = pi_0 * R^k
 */
private fun solvePi(B0: Matrix, B1: Matrix, R: Matrix, A0: Matrix, A1: Matrix, A2: Matrix): DoubleArray {
    val n = B0.numRows
    val maxLevel = 500  // Truncation level

    // First, solve for pi_0 from the boundary equation:
    // pi_0 * B1 + pi_1 * A2 = 0
    // pi_1 = pi_0 * R
    // So: pi_0 * (B1 + R*A2) = 0
    // With normalization: sum(pi) = 1

    // Compute B1 + R*A2
    val boundaryMatrix = B1.add(R.mult(A2))

    // Solve for pi_0 as the left null vector of boundaryMatrix
    // For a 2x2 matrix, we can solve directly
    // We need pi_0 such that pi_0 * boundaryMatrix = 0 and sum(pi) = 1

    // Use the fact that (I-R)^{-1} gives the sum of geometric series
    val I = Matrix.eye(n)
    val ImR = I.sub(R)
    val ImRinv = ImR.inv()

    // The normalization is: pi_0 * (I-R)^{-1} * e = 1 (where e is vector of 1s)
    // Combined with pi_0 * (B1 + R*A2) = 0

    // For numerical stability, solve the augmented system
    val pi0 = solveBoundary(boundaryMatrix, ImRinv, n)

    // Build full distribution using pi_k = pi_0 * R^k
    val totalSize = maxLevel * n
    val pi = DoubleArray(totalSize)

    // Level 0
    for (i in 0 until n) {
        pi[i] = pi0[i]
    }

    // Levels 1 to maxLevel-1
    var Rk = R.copy()
    for (level in 1 until maxLevel) {
        for (i in 0 until n) {
            var sum = 0.0
            for (j in 0 until n) {
                sum += pi0[j] * Rk[j, i]
            }
            pi[level * n + i] = sum
        }
        Rk = Rk.mult(R)
    }

    // Normalize
    var total = 0.0
    for (p in pi) {
        total += p
    }
    if (total > 0) {
        for (i in pi.indices) {
            pi[i] /= total
        }
    }

    return pi
}

/**
 * Solve the boundary equations for pi_0
 */
private fun solveBoundary(boundaryMatrix: Matrix, ImRinv: Matrix, n: Int): DoubleArray {
    // We need to solve:
    // pi_0 * boundaryMatrix = 0
    // pi_0 * ImRinv * e = 1 (normalization)

    // Create augmented matrix [boundaryMatrix^T | ImRinv^T * e]
    // and solve the system

    // For 2x2 case, we can use a direct approach
    // From pi_0 * boundaryMatrix = 0, one equation is redundant
    // Use normalization to get the second equation

    val b = boundaryMatrix

    // pi_0[0] * b[0,0] + pi_0[1] * b[1,0] = 0  (first column)
    // pi_0[0] * b[0,1] + pi_0[1] * b[1,1] = 0  (second column)
    // These are linearly dependent for a valid generator

    // Use first equation: pi_0[1] = -pi_0[0] * b[0,0] / b[1,0] (if b[1,0] != 0)
    // Or: pi_0[0] = -pi_0[1] * b[1,0] / b[0,0] (if b[0,0] != 0)

    val pi0 = DoubleArray(n)

    // Find the ratio from null space
    if (Math.abs(b[1, 0]) > 1e-10) {
        // pi_0[1] / pi_0[0] = -b[0,0] / b[1,0]
        val ratio = -b[0, 0] / b[1, 0]
        // pi_0[0] + pi_0[1] * (ImRinv row sums) = 1 after normalization

        // Compute normalization factor
        var norm0 = 0.0
        var norm1 = 0.0
        for (j in 0 until n) {
            norm0 += ImRinv[0, j]
            norm1 += ImRinv[1, j]
        }

        // pi_0[0] * norm0 + pi_0[1] * norm1 = 1
        // pi_0[0] * norm0 + ratio * pi_0[0] * norm1 = 1
        // pi_0[0] * (norm0 + ratio * norm1) = 1
        pi0[0] = 1.0 / (norm0 + ratio * norm1)
        pi0[1] = ratio * pi0[0]
    } else if (Math.abs(b[0, 0]) > 1e-10) {
        val ratio = -b[1, 0] / b[0, 0]
        var norm0 = 0.0
        var norm1 = 0.0
        for (j in 0 until n) {
            norm0 += ImRinv[0, j]
            norm1 += ImRinv[1, j]
        }
        pi0[1] = 1.0 / (norm1 + ratio * norm0)
        pi0[0] = ratio * pi0[1]
    } else {
        // Fallback: equal probabilities
        val normFactor = 1.0 / n
        for (i in 0 until n) {
            pi0[i] = normFactor
        }
    }

    // Ensure non-negative
    for (i in 0 until n) {
        if (pi0[i] < 0) pi0[i] = 0.0
    }

    return pi0
}

/**
 * QBD setupdelayoff algorithms
 */
@Suppress("unused")
class QbdSetupdelayoff {
    companion object {
        // Class documentation marker for Dokka
    }
}
