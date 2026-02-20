/**
 * @file Fixed-point iteration for Laguerre expansion method with think times
 * 
 * Implements the fixed-point iteration algorithm for the Laguerre expansion method
 * in closed queueing networks with think times (delays). Solves for both spatial and
 * temporal saddle points to enable accurate asymptotic approximation.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.GlobalConstants.Inf
import jline.io.Ret
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Fixed-point iteration used in the logistic expansion method in models with delays
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @return fixed point
 */
fun pfqn_le_fpiZ(L: Matrix, N: Matrix, Z: Matrix): Ret.pfqnLeFpiZ {
    val M = L.numRows
    val R = L.numCols
    val eta = N.elementSum() + M
    val u = Matrix(M, 1)
    for (i in 0..<M) {
        u[i] = 1.0 / M
    }
    var v = eta + 1
    var u_1 = Matrix(M, 1)
    u_1.fill(Inf)
    var v_1 = Inf
    var d: Matrix? = Matrix(0, M)
    val u_abs_diff = Matrix(M, 1)

    for (i in 0..<M) {
        u_abs_diff[i] = FastMath.abs(u[i] - u_1[i])
    }

    while (u_abs_diff.elementSum() > 1e-10) {
        u_1 = u.copy()
        v_1 = v
        for (i in 0..<M) {
            u[i] = 1 / eta
            for (r in 0..<R) {
                val L_col_r = Matrix(L.numRows, 1)
                Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
                u[i] =
                    (u[i] + N[r] / eta * (Z[r] + v * L[i, r]) * u_1[i] / (Z[r] + v * u_1.transpose().mult(L_col_r)[0]))
            }
        }

        val xi = Matrix(1, R)
        for (r in 0..<R) {
            val L_col_r = Matrix(L.numRows, 1)
            Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
            xi[r] = N[r] / (Z[r] + v * u_1.transpose().mult(L_col_r)[0])
        }
        v = eta + 1
        for (r in 0..<R) {
            v -= xi[r] * Z[r]
        }

        for (i in 0..<M) {
            u_abs_diff[i] = FastMath.abs(u[i] - u_1[i])
        }

        d = Matrix.concatRows(d, u_abs_diff.transpose().elementIncrease(FastMath.abs(v - v_1)), null)
    }

    return Ret.pfqnLeFpiZ(u, v, d)
}
/**
 * PFQN le fpiZ algorithms
 */
@Suppress("unused")
class PfqnLeFpizAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}