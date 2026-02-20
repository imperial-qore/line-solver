/**
 * @file Markovian Arrival Process skewness computation
 * 
 * Computes skewness of MAP inter-arrival times measuring asymmetry in distributions.
 * Important for statistical characterization and shape analysis of arrival patterns.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import kotlin.math.sqrt

/**
 * Computes the skewness of the inter-arrival times for a MAP.
 *
 *
 * The skewness measures the asymmetry of the distribution of inter-arrival times. It is calculated using the third central moment
 * normalized by the cube of the standard deviation.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @return the skewness of the inter-arrival times
 */

fun map_skew(D0: Matrix, D1: Matrix): Double {
    val m: MutableMap<Int, Double> = HashMap()
    for (i in 1..3) {
        m[i] = map_moment(D0, D1, i)
    }
    val M3 = m[3]!! - 3 * m[2]!! * m[1]!! + 2 * FastMath.pow(m[1]!!, 3)
    return M3 / FastMath.pow(sqrt(map_scv(D0, D1)) * m[1]!!, 3)
}

/**
 * Computes the skewness of the inter-arrival times for a MAP using a MatrixCell.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the skewness of the inter-arrival times.
 *
 * @param MAP the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @return the skewness of the inter-arrival times
 */
fun map_skew(MAP: MatrixCell): Double {
    return map_skew(MAP[0], MAP[1])
}
/**
 * MAP skew algorithms
 */
@Suppress("unused")
class MapSkewAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}