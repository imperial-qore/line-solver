/**
 * @file Laguerre expansion method for normalizing constant computation
 * 
 * Implements the Laguerre expansion approach for computing normalizing constants in
 * closed product-form queueing networks. Uses asymptotic approximation with Hessian-based
 * corrections to achieve high accuracy for large population systems.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.GlobalConstants
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.sqrt

/**
 * Logistic expansion method to compute the normalizing constant
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @return normalizing constant and its logarithm
 */
fun pfqn_le(L: Matrix, N: Matrix, Z: Matrix = Matrix(1, L.numCols)): Ret.pfqnNc {
    val M = L.numRows
    val R = L.numCols
    val lGn: Double
    val Gn: Double

    if (L.isEmpty || N.isEmpty || N.elementSum() == 0.0 || L.elementSum() < GlobalConstants.CoarseTol) {
        val tmp = Matrix(1, Z.numCols)
        for (i in 0..<tmp.length()) {
            tmp[i] = FastMath.log(Z.sumCols(i))
        }
        lGn = (-Matrix.factln(N).elementSum() + N.elementMult(tmp, null).elementSum())
        Gn = FastMath.exp(lGn)
    } else if (Z.isEmpty) {
        val ret = pfqn_le_fpi(L, N)
        val umax = ret.u
        val A = pfqn_le_hessian(L, N, umax.transpose())
        var S = 0.0
        for (r in 0..<R) {
            val L_col_r = Matrix(L.numRows, 1)
            Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
            S += N[r] * FastMath.log(umax.transpose().mult(L_col_r)[0])
        }

        val tmp = Matrix(1, N.length() + 1)
        Matrix.extract(N, 0, 1, 0, N.length(), tmp, 0, 0)
        tmp[N.length()] = (M - 1).toDouble()
        val log_umax = umax.copy()
        for (i in 0..<log_umax.length()) {
            log_umax[i] = FastMath.log(log_umax[i])
        }
        lGn =
            ((Maths.multinomialln(tmp) + Maths.factln(M - 1) + (M - 1) * FastMath.log(sqrt(2 * FastMath.PI)) - FastMath.log(
                sqrt(A.det()))) + log_umax.elementSum() + S)
        Gn = FastMath.exp(lGn)
    } else {
        val ret = pfqn_le_fpiZ(L, N, Z)
        val umax = ret.u
        val vmax = ret.v
        val A = pfqn_le_hessianZ(L, N, Z, umax.transpose(), vmax)
        var S = 0.0
        for (r in 0..<R) {
            val L_col_r = Matrix(L.numRows, 1)
            Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
            S += N[r] * FastMath.log(Z[r] + vmax * umax.transpose().mult(L_col_r)[0])
        }
        val log_umax = umax.copy()
        for (i in 0..<log_umax.length()) {
            log_umax[i] = FastMath.log(log_umax[i])
        }
        lGn = ((-Matrix.factln(N)
            .elementSum() - vmax + M * FastMath.log(vmax) + M * FastMath.log(sqrt(2 * FastMath.PI)) - FastMath.log(sqrt(
            A.det()))) + log_umax.elementSum() + S)
        Gn = FastMath.exp(lGn)
    }
    return Ret.pfqnNc(Gn, lGn)
}
/**
 * PFQN le algorithms
 */
@Suppress("unused")
class PfqnLeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}