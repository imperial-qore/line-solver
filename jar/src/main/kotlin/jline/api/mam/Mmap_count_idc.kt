/**
 * @file Marked Markovian Arrival Process counting process index of dispersion
 * 
 * Computes IDC values for each marked class in MMAP counting processes.
 * Essential for analyzing variability and burstiness patterns in multiclass systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the index of dispersion for counts (IDC) for a Markovian Arrival Process with marked arrivals (MMAP) over a time period.
 *
 *
 * The IDC is calculated as the variance of the counts divided by the mean of the counts. This method computes the IDC for each type of event
 * in the MMAP over a specified time period `t`.
 *
 * @param MMAP the MatrixCell containing the transition matrices of the MMAP, with D0, D1, ..., Dc representing different types of events
 * @param t    the time period over which to compute the IDC
 * @return a Matrix containing the IDC for each type of event over time `t`
 */
fun mmap_count_idc(MMAP: MatrixCell, t: Double): Matrix {
    val m = mmap_count_mean(MMAP, t)
    val v = mmap_count_var(MMAP, t)
    return v.elementDiv(m)
}
/**
 * MMAP count idc algorithms
 */
@Suppress("unused")
class MmapCountIdcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}