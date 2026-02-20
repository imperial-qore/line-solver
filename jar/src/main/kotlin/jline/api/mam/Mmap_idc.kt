/**
 * @file Marked Markovian Arrival Process asymptotic index of dispersion
 * 
 * Computes asymptotic IDC for each marked class in MMAP as time approaches infinity.
 * Essential for analyzing long-term variability characteristics in multiclass systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the asymptotic index of dispersion for counts (IDC) for a Markovian Arrival Process with marked arrivals (MMAP).
 *
 *
 * This method calculates the IDC for each type of event in the MMAP as t approaches infinity, using an approximation for large time `t`.
 * The IDC is computed as the variance of the counts divided by the mean of the counts for each type of event.
 *
 * @param MMAP the MatrixCell containing the transition matrices of the MMAP, with D0, D1, ..., Dc representing different types of events
 * @return a Matrix containing the IDC for each type of event as time approaches infinity
 */
fun mmap_idc(MMAP: MatrixCell): Matrix {
    val tinf = 1e6 / mmap_lambda(MMAP).elementSum()
    val m = mmap_count_mean(MMAP, tinf)
    val v = mmap_count_var(MMAP, tinf)
    return v.elementDiv(m)
}
/**
 * MMAP idc algorithms
 */
@Suppress("unused")
class MmapIdcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}