/**
 * @file Markovian Arrival Process temporal scaling operations
 * 
 * Rescales MAP inter-arrival time distributions to achieve specified mean values.
 * Essential for parameter adjustment and model calibration in queueing analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Rescales the mean inter-arrival time of a Markovian Arrival Process (MAP) to a specified new mean.
 *
 *
 * This method adjusts the transition matrices D0 and D1 of the MAP to achieve a desired mean inter-arrival time
 * while preserving the structure and relative rates of transitions. The MAP is represented by two matrices: D0 and D1,
 * where D0 is the hidden transition matrix and D1 is the visible transition matrix.
 *
 * @param D0      the hidden transition matrix of the MAP
 * @param D1      the visible transition matrix of the MAP
 * @param newMean the desired new mean inter-arrival time
 * @return a MatrixCell containing the scaled MAP transition matrices
 */

fun map_scale(D0: Matrix, D1: Matrix, newMean: Double): MatrixCell {
    val D = MatrixCell()
    D[0] = D0.copy()
    D[1] = D1.copy()
    map_normalize(D[0], D[1])
    val ratio = map_mean(D0, D1) / newMean
    D[0].scaleEq(ratio)
    D[1].scaleEq(ratio)
    return D
}

/**
 * Rescales the mean inter-arrival time of a MAP stored in a MatrixCell that contains the MAP's transition matrices.
 *
 * @param MAP     a MatrixCell containing the transition matrices D0 and D1 of the MAP
 * @param newMean the desired new mean inter-arrival time
 * @return a MatrixCell containing the scaled MAP transition matrices
 */
fun map_scale(MAP: MatrixCell, newMean: Double): MatrixCell {
    return map_scale(MAP[0], MAP[1], newMean)
}
/**
 * MAP scale algorithms
 */
@Suppress("unused")
class MapScaleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}