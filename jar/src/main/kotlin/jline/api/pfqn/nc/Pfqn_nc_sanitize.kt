/**
 * @file Parameter sanitization for product-form queueing network models
 * 
 * Sanitizes and preprocesses parameters for product-form queueing network models to
 * avoid numerical degeneracies and improve computational stability. Handles zero demands,
 * scaling normalization, and class reordering for robust normalizing constant computation.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.io.Ret
import jline.GlobalConstants
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.*

/**
 * Sanitizes product-form model parameters to avoid degeneracies
 *
 * @param lambda - arrival rates for open classes
 * @param L      - demands at all stations
 * @param N      - number of jobs for each class
 * @param Z      - think time for each class
 * @param atol   - absolute numerical tolerance
 * @return sanitized parameters
 */

fun pfqn_nc_sanitize(lambda: Matrix, L: Matrix, N: Matrix, Z: Matrix, atol: Double): Ret.pfqnNcSanitize {
    var L_new = L.copy()
    var Z_new = Z.copy()
    L_new.removeNaN()
    Z_new.removeNaN()
    var L_tmp = Matrix(L_new.numRows, 0)
    var N_tmp = Matrix(N.numRows, 0)
    var Z_tmp = Matrix(Z_new.numRows, 0)
    var lambda_tmp = Matrix(lambda.numRows, 0)
    for (i in 0..<N.length()) {
        if (FastMath.abs(N[i]) >= GlobalConstants.FineTol && !(L_new.sumCols(i) + Z_new.sumCols(i) < atol)) {
            val L_col = Matrix.extractColumn(L_new, i, null)
            val N_col = Matrix.extractColumn(N, i, null)
            val Z_col = Matrix.extractColumn(Z_new, i, null)
            val lambda_col = Matrix.extractColumn(lambda, i, null)
            L_tmp = Matrix.concatColumns(L_tmp, L_col, null)
            N_tmp = Matrix.concatColumns(N_tmp, N_col, null)
            Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null)
            lambda_tmp = Matrix.concatColumns(lambda_tmp, lambda_col, null)
        }
    }
    L_new = L_tmp
    Z_new = Z_tmp
    var N_new = N_tmp
    val lambda_new = lambda_tmp

    var lGremaind = 0.0

    var L_zeroDemand = Matrix(L_new.numRows, 0)
    var Z_zeroDemand = Matrix(Z_new.numRows, 0)
    var N_zeroDemand = Matrix(N_new.numRows, 0)
    L_tmp = Matrix(L_new.numRows, 0)
    N_tmp = Matrix(N_new.numRows, 0)
    Z_tmp = Matrix(Z_new.numRows, 0)

    for (i in 0..<L_new.numCols) {
        val L_col = Matrix.extractColumn(L_new, i, null)
        val N_col = Matrix.extractColumn(N_new, i, null)
        val Z_col = Matrix.extractColumn(Z_new, i, null)
        if (L_new[i] < atol) {
            L_zeroDemand = Matrix.concatColumns(L_zeroDemand, L_col, null)
            Z_zeroDemand = Matrix.concatColumns(Z_zeroDemand, Z_col, null)
            N_zeroDemand = Matrix.concatColumns(N_zeroDemand, N_col, null)
        } else {
            L_tmp = Matrix.concatColumns(L_tmp, L_col, null)
            N_tmp = Matrix.concatColumns(N_tmp, N_col, null)
            Z_tmp = Matrix.concatColumns(Z_tmp, Z_col, null)
        }
    }

    val log_Z_zeroDemand = Z_zeroDemand.copy()
    val log_N_zeroDemand = N_zeroDemand.copy()

    for (i in 0..<Z_zeroDemand.numRows) {
        for (j in 0..<Z_zeroDemand.numCols) {
            log_Z_zeroDemand[i, j] = FastMath.log(Z_zeroDemand[i, j])
        }
    }

    for (i in 0..<N_zeroDemand.numRows) {
        for (j in 0..<N_zeroDemand.numCols) {
            log_N_zeroDemand[i, j] = FastMath.log(N_zeroDemand[i, j])
        }
    }

    if (!L_zeroDemand.isEmpty) {
        lGremaind += (N_zeroDemand.mult(log_Z_zeroDemand.transpose())[0] - log_N_zeroDemand.elementSum())
    }
    L_new = L_tmp
    Z_new = Z_tmp
    N_new = N_tmp

    if (!L_new.isEmpty) {
        var Lmax = Matrix(1, L_new.numCols)
        for (i in 0..<Lmax.length()) {
            val L_col_i = Matrix.extractColumn(L_new, i, null)
            Lmax[i] = L_col_i.elementMax()
        }
        if (Lmax.isEmpty) {
            Lmax = Matrix(1, Z_new.numCols)
            Lmax.ones()
        }
        val repmat_Lmax_L = Lmax.repmat(L_new.numRows, 1)
        val repmat_Lmax_Z = Lmax.repmat(Z_new.numRows, 1)
        for (i in 0..<L_new.numRows) {
            for (j in 0..<L_new.numCols) {
                L_new[i, j] = L_new[i, j] / repmat_Lmax_L[i, j]
            }
        }

        for (i in 0..<Z_new.numRows) {
            for (j in 0..<Z_new.numCols) {
                Z_new[i, j] = Z_new[i, j] / repmat_Lmax_Z[i, j]
            }
        }

        val Lmax_log = Lmax.copy().transpose()
        for (i in 0..<Lmax_log.length()) {
            Lmax_log[i] = FastMath.log(Lmax_log[i])
        }
        lGremaind += N_new.mult(Lmax_log)[0]

        if (!Z_new.isEmpty) {
            val index = arrayOfNulls<Int>(Z_new.numCols)
            for (i in index.indices) {
                index[i] = i
            }

            val Z_sum_col = Matrix(1, Z_new.numCols)
            for (i in 0..<Z_sum_col.length()) {
                Z_sum_col[i] = Z_new.sumCols(i)
            }

            Arrays.sort(index) { i1, i2 -> Z_sum_col[i1!!].compareTo(Z_sum_col[i2!!]) }

            L_tmp = Matrix(L_new.numRows, 0)
            Z_tmp = Matrix(Z_new.numRows, 0)
            N_tmp = Matrix(N_new.numRows, 0)

            for (i in index.indices) {
                if (!L_new.isEmpty) {
                    L_tmp = Matrix.concatColumns(L_tmp, Matrix.extractColumn(L_new, index[i]!!, null), null)
                }
                Z_tmp = Matrix.concatColumns(Z_tmp, Matrix.extractColumn(Z_new, index[i]!!, null), null)
                N_tmp = Matrix.concatColumns(N_tmp, Matrix.extractColumn(N_new, index[i]!!, null), null)
            }
            L_new = L_tmp
            Z_new = Z_tmp
            N_new = N_tmp
        }
    }

    return Ret.pfqnNcSanitize(lambda_new, L_new, N_new, Z_new, lGremaind)
}
/**
 * PFQN nc sanitize algorithms
 */
@Suppress("unused")
class PfqnNcSanitizeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}