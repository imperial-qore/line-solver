/**
 * @file Marked Markovian Arrival Process counting process variance analysis
 * 
 * Computes variance vectors for counting processes of each marked class in MMAP.
 * Essential for analyzing variability and dispersion in multiclass arrival systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the variance of the count vector of events of different types in a Markovian Arrival Process with marked arrivals (MMAP) over a time period.
 *
 *
 * This method calculates the variance of the number of events of each type over a given time period `t` for a MMAP,
 * represented by a set of matrices. The MMAP contains multiple types of arrivals, each represented by a separate matrix
 * in the `MMAP` MatrixCell. The variance for each type of event is computed using the stationary distribution
 * of the underlying Markov chain (`theta`), the infinitesimal generator matrix (`Q`), and the event matrices.
 *
 * @param MMAP the MatrixCell containing the transition matrices of the MMAP, with D0, D1, ..., Dc representing different types of events
 * @param t    the time period over which to compute the variances of counts
 * @return a Matrix containing the variance vector of counts over time `t`
 */
fun mmap_count_var(MMAP: MatrixCell, t: Double): Matrix {
    val D0 = MMAP[0]
    val D1 = MMAP[1]
    val theta = map_piq(D0, D1)
    val C = MMAP.size() - 2
    val vt = Matrix(1, C)
    val Q = map_infgen(D0, D1)
    val n = D0.numRows
    val I = Matrix.eye(n)
    val e = Matrix.ones(n, 1)
    val tmp = Matrix.ones(n, 1).mult(theta).add(-1.0, Q).inv()
    val lk = Matrix(1, C)
    val ck = MatrixCell(C)
    val dk = MatrixCell(C)
    val llk = Matrix(1, C)
    for (c in 0..<C) {
        val Dc = MMAP[2 + c]
        lk[c] = theta.mult(Dc.mult(e)).value()
        ck[c] = theta.mult(Dc.mult(tmp))
        dk[c] = tmp.mult(Dc.mult(e))
        llk[c] = theta.mult(Dc.mult(e)).value()
    }
    for (c in 0..<C) {
        val Dc = MMAP[2 + c]
        var vc = llk[c] - 2 * lk[c] * lk[c]
        vc += ck[c].scale(2.0).mult(Dc.mult(e)).value()
        vc *= t
        vc -= 2 * ck[c].mult(I.sub(Maths.matrixExp(Q.scale(t)))).mult(dk[c]).value()
        vt[c] = vc
    }

    return vt
}
/**
 * MMAP count var algorithms
 */
@Suppress("unused")
class MmapCountVarAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}