/**
 * @file Markovian Arrival Process squared coefficient of variation analysis
 * 
 * Computes SCV of MAP inter-arrival times as normalized dispersion measure.
 * Fundamental metric for characterizing variability and burstiness in arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the squared coefficient of variation (SCV) of the inter-arrival times of a Markovian Arrival Process (MAP).
 *
 *
 * The SCV is a normalized measure of the dispersion of the inter-arrival time distribution. It is calculated as the
 * variance of the inter-arrival times divided by the square of the mean inter-arrival time. The MAP is represented by
 * two matrices: D0 and D1, where D0 is the hidden transition matrix and D1 is the visible transition matrix.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @return the squared coefficient of variation (SCV) of the inter-arrival times
 */

fun map_scv(D0: Matrix, D1: Matrix): Double {
    val e1 = map_moment(D0, D1, 1)
    val e2 = map_moment(D0, D1, 2)

    val `var` = e2 - e1 * e1
    val scv = `var` / e1 / e1
    return scv
}

/**
 * Computes the squared coefficient of variation (SCV) of the inter-arrival times of a MAP
 * stored in a MatrixCell that contains the MAP's transition matrices.
 *
 * @param MAP a MatrixCell containing the transition matrices D0 and D1 of the MAP
 * @return the squared coefficient of variation (SCV) of the inter-arrival times
 */
fun map_scv(MAP: MatrixCell): Double {
    return map_scv(MAP[0], MAP[1])
}
/**
 * MAP scv algorithms
 */
@Suppress("unused")
class MapScvAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}