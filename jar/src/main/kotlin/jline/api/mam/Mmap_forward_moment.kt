/**
 * @file Marked Markovian Arrival Process forward moment analysis
 * 
 * Computes forward moments of MMAP inter-arrival times for each marked class.
 * Essential for statistical characterization and fitting of multiclass arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.CombinatoricsUtils

/**
 * Computes the forward moments of an MMAP for specified orders with normalization.
 *
 *
 * This method calculates the forward moments of the inter-arrival times for each type in the MMAP
 * according to the specified orders, with moments normalized by default.
 *
 * @param MMAP   the MMAP from which to compute the forward moments
 * @param ORDERS a matrix specifying the orders of the moments to be computed
 * @return a matrix where each row corresponds to a type and each column corresponds to a moment order
 */
@JvmOverloads
fun mmap_forward_moment(MMAP: MatrixCell, ORDERS: Matrix, NORM: Int = 1): Matrix {
    val C = MMAP.size() - 2
    val K = ORDERS.length()

    val MOMENTS = Matrix(C, K, C * K)
    val pie = map_pie(MMAP[0], MMAP[1])

    val neg_D0 = MMAP[0].copy()
    neg_D0.scaleEq(-1.0)
    val M = neg_D0.inv()

    for (a in 0..<C) {
        val pa = if (NORM == 1) {
            pie.mult(M.mult(MMAP[2 + a])).elementSum()
        } else {
            1.0
        }
        for (h in 0..<ORDERS.length()) {
            val k = ORDERS[h].toInt()
            val fk = CombinatoricsUtils.factorial(k).toDouble()
            MOMENTS[a, h] = fk / pa * pie.mult(M.mult(MMAP[2 + a])).mult(Matrix.pow(M, k)).elementSum()
        }
    }
    return MOMENTS
}
/**
 * MMAP forward moment algorithms
 */
@Suppress("unused")
class MmapForwardMomentAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}