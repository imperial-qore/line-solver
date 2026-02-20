package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * MAP variance computation algorithms.
 * 
 * Provides methods for computing the variance of inter-arrival times in Markovian Arrival 
 * Processes (MAP). The variance measures the variability in the arrival process and is
 * important for characterizing the burstiness of traffic patterns.
 *
 * @since LINE 3.0
 */

/**
 * Computes the variance of the inter-arrival times for a MAP.
 *
 *
 * The variance is calculated as the difference between the second moment and the square of the mean inter-arrival time.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @return the variance of the inter-arrival times
 */

fun map_var(D0: Matrix, D1: Matrix): Double {
    return map_moment(D0, D1, 2) - FastMath.pow(map_mean(D0, D1), 2)
}

/**
 * Computes the variance of the inter-arrival times for a MAP using a MatrixCell.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the variance of the inter-arrival times.
 *
 * @param MAP the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @return the variance of the inter-arrival times
 */
fun map_var(MAP: MatrixCell): Double {
    return map_var(MAP[0], MAP[1])
}

/**
 * MAP variance computation algorithms.
 * 
 * Provides methods for computing the variance of inter-arrival times in Markovian Arrival 
 * Processes (MAP). The variance measures the variability in the arrival process and is
 * important for characterizing the burstiness of traffic patterns.
 */
@Suppress("unused")
class MapVar {
    companion object {
        // Class documentation marker for Dokka
    }
}
