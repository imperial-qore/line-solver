/**
 * @file GI/M/1-type R Matrix Computation
 *
 * Computes the R matrix for GI/M/1-type Markov chains using various methods.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.io.line_warning
import jline.util.matrix.Matrix

/**
 * Options for GIM1_R solver
 */
data class GIM1ROptions(
    val mode: String = "FI",
    val maxNumIt: Int = 10000,
    val verbose: Int = 0
)

/**
 * Computes the R matrix for a GI/M/1-type Markov chain.
 *
 * R is the minimal nonnegative solution to:
 * R = A0 + R*A1 + R^2*A2 + ... + R^max*A_max
 *
 * @param A Block matrix [A0 A1 A2 ... A_max]
 * @param options Solver options
 * @return R matrix
 */
fun gim1_R(A: Matrix, options: GIM1ROptions = GIM1ROptions()): Matrix {
    val m = A.numRows
    val dega = A.numCols / m - 1

    // Use functional iteration by default
    var numit = 0
    var check = 1.0
    var R = Matrix.zeros(m, m)

    when (options.mode.uppercase()) {
        "FI" -> {
            // Functional Iteration: R(n+1) = A0 + R(n)*A1 + R(n)^2*A2 + ...
            while (check > 1e-14 && numit < options.maxNumIt) {
                val Rold = R.copy()

                // Compute new R
                R = A.extractCols(0, m).copy()  // Start with A0
                var Rpow = Rold.copy()

                for (i in 1..dega) {
                    R = R.add(Rpow.mult(A.extractCols(i * m, (i + 1) * m)))
                    Rpow = Rpow.mult(Rold)
                }

                check = R.sub(Rold).infinityNorm()
                numit++

            }
        }

        "NI" -> {
            // Newton Iteration for faster convergence
            while (check > 1e-14 && numit < options.maxNumIt) {
                val Rold = R.copy()

                // Compute A(R) and A'(R)
                var AR = A.extractCols(0, m).copy()
                var ARprime = Matrix.zeros(m, m)
                var Rpow = Matrix.eye(m)

                for (i in 1..dega) {
                    val Ai = A.extractCols(i * m, (i + 1) * m)
                    ARprime = ARprime.add(Rpow.mult(Ai).scale(i.toDouble()))
                    Rpow = Rpow.mult(Rold)
                    AR = AR.add(Rpow.mult(Ai))
                }

                // Newton step: R_new = R - (R - A(R)) * (I - A'(R))^{-1}
                val residual = Rold.sub(AR)
                val jacobian = Matrix.eye(m).sub(ARprime)

                try {
                    val correction = residual.mult(jacobian.inv())
                    R = Rold.sub(correction)
                } catch (e: Exception) {
                    // Fall back to functional iteration step
                    R = AR
                }

                check = R.sub(Rold).infinityNorm()
                numit++

            }
        }
    }

    if (numit == options.maxNumIt) {
        line_warning("GIM1_R", "Maximum number of iterations %d reached", numit)
    }

    return R
}
