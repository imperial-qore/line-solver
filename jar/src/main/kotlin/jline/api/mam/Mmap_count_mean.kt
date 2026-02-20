/**
 * @file Marked Markovian Arrival Process counting process mean analysis
 * 
 * Computes mean count vectors for each marked class in MMAP counting processes.
 * Essential for determining expected arrival rates per class in multiclass systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the mean count vector of events of different types in a Markovian Arrival Process with marked arrivals (MMAP) over a time period.
 *
 *
 * This method calculates the expected number of events of each type over a given time period `t` for a MMAP, which is represented by
 * a set of matrices. The MMAP contains multiple types of arrivals, each represented by a separate matrix in the `MMAP` MatrixCell.
 * The count mean for each type of event is computed using the stationary distribution of the underlying Markov chain (`theta`) and
 * the event matrices.
 *
 * @param MMAP the MatrixCell containing the transition matrices of the MMAP, with D0, D1, ..., Dc representing different types of events
 * @param t    the time period over which to compute the mean counts
 * @return a Matrix containing the mean count vector of events over time `t`
 */
fun mmap_count_mean(MMAP: MatrixCell, t: Double): Matrix {
    val theta = map_piq(MMAP[0], MMAP[1])
    val C = MMAP.size() - 2
    val et = Matrix(1, C)
    for (c in 0..<C) {
        et[c] = theta.mult(MMAP[2 + c]).elementSum() * t
    }
    return et
}
/**
 * MMAP count mean algorithms
 */
@Suppress("unused")
class MmapCountMeanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}