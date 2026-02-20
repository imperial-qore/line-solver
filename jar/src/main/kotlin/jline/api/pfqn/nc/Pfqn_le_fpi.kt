/**
 * @file Fixed-point iteration for Laguerre expansion method
 * 
 * Implements the fixed-point iteration algorithm used in the Laguerre expansion method
 * for computing normalizing constants. Iteratively solves for the saddle point of the
 * integrand to enable accurate asymptotic approximation.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.GlobalConstants.Inf
import jline.io.Ret
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Fixed-point iteration used in the logistic expansion method
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @return fixed point
 */
fun pfqn_le_fpi(L: Matrix, N: Matrix): Ret.pfqnLeFpi {
    val M = L.numRows
    val R = L.numCols
    val u = Matrix(M, 1)
    for (i in 0..<M) {
        u[i] = 1.0 / M
    }
    var u_1 = Matrix(M, 1)
    u_1.fill(Inf)
    var d: Matrix? = Matrix(0, M)
    val u_abs_diff = Matrix(M, 1)

    for (i in 0..<M) {
        u_abs_diff[i] = FastMath.abs(u[i] - u_1[i])
    }

    while (u_abs_diff.elementSum() > 1e-10) {
        u_1 = u.copy()
        for (i in 0..<M) {
            u[i] = 1 / (N.elementSum() + M)
            for (r in 0..<R) {
                val L_col_r = Matrix(L.numRows, 1)
                Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
                u[i] = (u[i] + N[r] / (N.elementSum() + M) * L[i, r] * u_1[i] / (u_1.transpose().mult(L_col_r)[0]))
            }
        }

        for (i in 0..<M) {
            u_abs_diff[i] = FastMath.abs(u[i] - u_1[i])
        }


        d = Matrix.concatRows(d, u_abs_diff.transpose(), null)
    }

    return Ret.pfqnLeFpi(u, d)
}
/**
 * PFQN le fpi algorithms
 */
@Suppress("unused")
class PfqnLeFpiAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}