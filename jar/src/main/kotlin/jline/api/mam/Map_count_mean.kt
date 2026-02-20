/**
 * @file Markovian Arrival Process counting process mean analysis
 * 
 * Computes mean number of arrivals in MAP counting processes over specified time intervals.
 * Fundamental for analyzing average arrival rates and capacity planning in queueing systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Computes the mean of the counting process over multiple specified interval lengths for a given Markovian Arrival Process (MAP).
 *
 *
 * This function calculates the mean number of events that occur in a MAP over each interval length provided in the array `t`.
 * It uses the arrival rate of the MAP, denoted as lambda (λ), and multiplies it by each interval length.
 * The results are returned in an array corresponding to the input intervals.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell.
 * @param t   An array of interval lengths over which to compute the mean.
 * @return An array of doubles, where each element represents the mean of the counting process over the corresponding interval length in `t`.
 */
fun map_count_mean(MAP: MatrixCell?, t: DoubleArray): DoubleArray {
    val ret = DoubleArray(t.size)
    val lambda = map_lambda(MAP!!)
    for (i in t.indices) {
        ret[i] = lambda * t[i]
    }
    return ret
}


/**
 * Computes the mean of the counting process over a specified interval length for a given Markovian Arrival Process (MAP).
 *
 *
 * This function calculates the mean number of events that occur in a MAP over an interval of length `t`.
 * It uses the arrival rate of the MAP, denoted as lambda (λ), and multiplies it by the interval length `t`.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell.
 * @param t   The length of the interval over which to compute the mean.
 * @return The mean of the counting process over the interval `t`.
 */
fun map_count_mean(MAP: MatrixCell?, t: Double): Double {
    return map_lambda(MAP!!) * t
}

/**
 * MAP counting process mean algorithms
 */
@Suppress("unused")
class MapCountMeanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}