/**
 * @file Marked Markovian Arrival Process cross-moment analysis
 * 
 * Computes cross-moment matrices between different marked classes in MMAP processes.
 * Essential for analyzing inter-class dependencies and correlation structures.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.CombinatoricsUtils
import org.apache.commons.math3.util.FastMath

/**
 * Computes the k-th cross-moment matrix for a given MMAP.
 *
 *
 * This method calculates the k-th cross-moments for each pair of types in the MMAP.
 * The cross-moment measures the expected value of the product of k inter-arrival times
 * for different types in the MMAP. The result is a CxC matrix, where C is the number of types.
 *
 * @param mmap the MMAP for which to compute the cross-moments
 * @param k    the order of the moment to compute
 * @return a matrix representing the k-th cross-moment for each pair of types
 */
fun mmap_cross_moment(mmap: MatrixCell, k: Int): Matrix {
    val C = mmap.size() - 2

    val TG = MatrixCell()
    val MC = Matrix(C, C, FastMath.pow(C.toDouble(), 2).toInt())

    for (i in 0..<C) {
        val a = mmap[0].inv().mult(mmap[2 + i])
        a.scaleEq(-1.0)
        TG[i] = map_pie(mmap[0], mmap[1]).mult(a).sumCols()
    }

    for (i in 0..<C) {
        val a = mmap[0].inv().mult(mmap[2 + i])
        a.scaleEq(-1.0)
        val start = map_pie(mmap[0], mmap[1]).mult(a).mult(TG[i].inv())
        for (j in 0..<C) {
            val neg_mmap0 = mmap[0].copy()
            neg_mmap0.scaleEq(-1.0)
            MC[i, j] =
                CombinatoricsUtils.factorial(k) * start.mult((Matrix.pow(neg_mmap0.inv(), k + 1)).mult(mmap[2 + j]))
                    .elementSum()
            MC[i, j] = MC[i, j] / start.mult(neg_mmap0.mult(mmap[2 + j].inv())).elementSum()
        }
    }
    return MC
}
/**
 * MMAP cross moment algorithms
 */
@Suppress("unused")
class MmapCrossMomentAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}