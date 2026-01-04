package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the time-reversed version of a Markovian Arrival Process with marked arrivals (MMAP).
 *
 *
 * This method takes an MMAP and returns its time-reversed version. The time-reversed MMAP is computed by transposing
 * the transition matrices and adjusting them using the stationary distribution of the underlying Markov chain.
 * The resulting matrices represent the same process viewed in reverse time.
 *
 * @param mmap the MatrixCell containing the transition matrices of the original MMAP
 * @return a MatrixCell representing the time-reversed MMAP
 */

fun mmap_timereverse(mmap: MatrixCell): MatrixCell {
    val K = mmap.size()
    val piq = map_piq(mmap[0], mmap[1])
    val D = Matrix.diag(*piq.toArray1D())
    val iD = D.inv()
    val DK = arrayOfNulls<Matrix>(K)
    for (k in 0..<K) {
        DK[k] = Matrix(iD.mult(mmap[k].transpose()).mult(D))
    }
    return MatrixCell(DK)
}
/**
 * MMAP timereverse algorithms
 */
@Suppress("unused")
class MmapTimereverse {
    companion object {
        // Class documentation marker for Dokka
    }
}