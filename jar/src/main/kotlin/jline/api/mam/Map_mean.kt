package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * MAP mean inter-arrival time computation algorithms.
 * 
 * Provides methods for computing the mean inter-arrival time of Markovian Arrival Processes (MAP).
 * The mean inter-arrival time is the reciprocal of the arrival rate and represents the expected
 * time between consecutive arrivals in the process.
 *
 * @since LINE 3.0
 */

/**
 * Computes the mean inter-arrival time of a Markovian Arrival Process (MAP).
 *
 *
 * The mean inter-arrival time is calculated as the inverse of the overall arrival rate (lambda) of the MAP.
 * This rate is derived from the hidden (D0) and visible (D1) transition matrices of the MAP.
 *
 *
 * @param D0 The hidden transition matrix of the MAP, representing transitions without visible events.
 * @param D1 The visible transition matrix of the MAP, representing transitions with visible events.
 * @return The mean inter-arrival time of the MAP.
 */

fun map_mean(D0: Matrix, D1: Matrix): Double {
    return 1.0 / map_lambda(D0, D1)
}

/**
 * Computes the mean inter-arrival time of a Markovian Arrival Process (MAP) using matrices stored in a MatrixCell.
 *
 *
 * This is a convenience method that extracts the D0 and D1 matrices from a given MAP stored in a MatrixCell and computes the
 * mean inter-arrival time, which is the inverse of the overall arrival rate (lambda) of the MAP.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices.
 * @return The mean inter-arrival time of the MAP.
 */

fun map_mean(MAP: MatrixCell): Double {
    return map_mean(MAP[0], MAP[1])
}

/**
 * MAP mean inter-arrival time computation algorithms.
 * 
 * Provides methods for computing the mean inter-arrival time of Markovian Arrival Processes (MAP).
 * The mean inter-arrival time is the reciprocal of the arrival rate and represents the expected
 * time between consecutive arrivals in the process.
 */
@Suppress("unused")
class MapMean {
    companion object {
        // Class documentation marker for Dokka
    }
}
