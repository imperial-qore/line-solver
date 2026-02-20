/**
 * @file Markovian Arrival Process feasibility validation
 * 
 * Comprehensive validation of MAP matrices including stochastic properties, numerical stability,
 * and irreducibility checks. Essential for ensuring MAP representations are mathematically valid.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.Utils
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.linear.EigenDecomposition
import org.apache.commons.math3.linear.MatrixUtils
import kotlin.math.abs

/**
 * Check the feasibility of a MAP with detailed validation.
 * 
 * This function performs comprehensive validation of MAP matrices following the
 * same approach as the MATLAB implementation. It checks:
 * 1. Matrix dimensions and finite values
 * 2. D0 diagonal elements (must be ≤ 0)
 * 3. D0 off-diagonal elements (must be ≥ 0)
 * 4. D1 elements (must be ≥ 0)
 * 5. Generator property (row sums = 0)
 * 6. Irreducibility and numerical stability
 *
 * @param MAP the MatrixCell representing the MAP transition matrices
 * @param TOL the tolerance level for numerical stability checks
 * @return true if the MAP is feasible, false otherwise
 */
fun map_checkfeasible(MAP: MatrixCell, TOL: Double): Boolean {
    val n = MAP[0].length()
    val D0 = MAP[0]
    val D1 = MAP[1]
    if (D0.hasNaN()) {
        return false
    }
    var D0_has_inf = false
    var D1_has_inf = false

    for (i in 0..<D0.numRows) {
        for (j in 0..<D0.numCols) {
            if (Utils.isInf(D0[i, j])) {
                D0_has_inf = true
                break
            }
        }
    }

    for (i in 0..<D1.numRows) {
        for (j in 0..<D1.numCols) {
            if (Utils.isInf(D1[i, j])) {
                D1_has_inf = true
                break
            }
        }
    }

    if (D0.hasNaN() || D1.hasNaN() || D0_has_inf || D1_has_inf) {
        return false
    }

    val neg_D0 = D0.copy()
    neg_D0.scaleEq(-1.0)

    val P = neg_D0.inv().mult(D1)
    val Q = D0.add(1.0, D1)

    for (i in 0..<D0.numRows) {
        for (j in 0..<D0.numCols) {
            if (abs(D0[i, j]) < TOL) {
                D0[i, j] = 0
            }
        }
    }

    for (i in 0..<D1.numRows) {
        for (j in 0..<D1.numCols) {
            if (abs(D1[i, j]) < TOL) {
                D1[i, j] = 0
            }
        }
    }

    for (i in 0..<P.numCols) {
        for (j in 0..<P.numCols) {
            if (abs(P[i, j]) < TOL) {
                P[i, j] = 0
            }
        }
    }

    for (i in 0..<Q.numCols) {
        for (j in 0..<Q.numCols) {
            if (abs(Q[i, j]) < TOL) {
                Q[i, j] = 0
            }
        }
    }


    for (i in 0..<n) {
        for (j in 0..<n) {
            if (i != j && D0[i, j] < 0) {
                return false
            }
            if (i == j && D0[i, j] > 0) {
                return false
            }
            if (D1[i, j] < 0) {
                return false
            }
            if (i != j && Q[i, j] < 0) {
                return false
            }
            if (i == j && Q[i, j] > 0) {
                return false
            }
            if (P[i, j] < 0) {
                return false
            }
        }

        if (abs(P.sumRows(i)) < 1 - n * TOL) {
            return false
        }

        if (abs(P.sumRows(i)) > 1 + n * TOL) {
            return false
        }
        if (abs(Q.sumRows(i)) < 0 - n * TOL) {
            return false
        }
        if (abs(Q.sumRows(i)) > 0 + n * TOL) {
            return false
        }
    }

    if (n < map_largemap()) {
        val P_data = P.toArray2D()
        val P_matrix = MatrixUtils.createRealMatrix(P_data)
        val P_eigenDecomposition = EigenDecomposition(P_matrix)
        val P_eigenvalues = P_eigenDecomposition.realEigenvalues
        var P_sum = 0
        for (e in P_eigenvalues) {
            if (e > 1 - TOL) {
                P_sum++
            }
            if (P_sum > 1) {
                return false
            }
        }

        val Q_data = Q.toArray2D()
        val Q_matrix = MatrixUtils.createRealMatrix(Q_data)
        val Q_eigenDecomposition = EigenDecomposition(Q_matrix)
        val Q_eigenvalues = Q_eigenDecomposition.realEigenvalues
        var Q_sum = 0
        for (e in Q_eigenvalues) {
            if (e > 0 - TOL) {
                Q_sum++
            }
            if (Q_sum > 1) {
                return false
            }
        }
    }

    return true
}

/**
 * Check the feasibility of a MAP given separate D0 and D1 matrices.
 * 
 * @param D0 the hidden transition matrix (transitions without arrivals)
 * @param D1 the visible transition matrix (transitions with arrivals)
 * @param TOL the tolerance level for numerical stability checks (default: 1e-14)
 * @return true if the MAP is feasible, false otherwise
 */
@JvmOverloads
fun map_checkfeasible(D0: Matrix, D1: Matrix, TOL: Double = 1e-14): Boolean {
    val MAP = MatrixCell(2)
    MAP[0] = D0
    MAP[1] = D1
    return map_checkfeasible(MAP, TOL)
}
/**
 * MAP checkfeasible algorithms
 */
@Suppress("unused")
class MapCheckfeasibleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}