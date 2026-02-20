/**
 * @file Markovian Arrival Process exponential distribution construction
 * 
 * Creates MAP representations of exponential inter-arrival time distributions with specified
 * means. The simplest MAP form with single-state Poisson arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Creates a Markovian Arrival Process (MAP) with an exponential inter-arrival time distribution.
 *
 *
 * The method constructs a MAP where the inter-arrival times are exponentially distributed with a specified mean.
 * This is achieved by setting the rate parameter (mu) of the exponential distribution as 1/mean.
 * The resulting MAP has a single state with transition matrices D0 and D1 that define the exponential arrival process.
 *
 * @param mean the desired mean of the exponential inter-arrival times
 * @return a MatrixCell containing the MAP transition matrices for the exponential distribution
 */

fun map_exponential(mean: Double): MatrixCell {
    val mu = 1.0 / mean
    val MAP = MatrixCell()
    val D0 = Matrix(1, 1, 1)
    D0[0, 0] = -mu
    val D1 = Matrix(1, 1, 1)
    D1[0, 0] = mu
    MAP[0] = D0
    MAP[1] = D1
    return MAP
}
/**
 * MAP exponential algorithms
 */
@Suppress("unused")
class MapExponentialAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}