/**
 * @file Markovian Arrival Process matrix normalization and sanitization
 * 
 * Sanitizes MAP matrices by ensuring non-negativity constraints and proper diagonal adjustments.
 * Essential for maintaining mathematical validity and numerical stability in MAP algorithms.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Sanitizes the (D0, D1) matrices of a Markovian Arrival Process (MAP) by ensuring all elements are non-negative and adjusting diagonal elements.
 *
 *
 * This function modifies the given hidden (D0) and visible (D1) transition matrices of a MAP to ensure they adhere to the necessary constraints
 * for a valid MAP. Specifically, it sets any negative elements to zero and adjusts the diagonal elements of D0 such that each row sums to zero,
 * preserving the stochastic nature of the matrix.
 *
 *
 * @param D0 The hidden transition matrix of the MAP, representing transitions that do not result in visible events.
 * @param D1 The visible transition matrix of the MAP, representing transitions that result in visible events.
 * @return A MatrixCell containing the sanitized D0 and D1 matrices.
 */

fun map_normalize(D0: Matrix, D1: Matrix?): MatrixCell {
    val D = MatrixCell()
    // Create copies to avoid modifying the inputs
    D[0] = D0.copy()
    D[1] = D1?.copy()
    val nPhases = D0.length().toLong()
    for (i in 0..<nPhases) {
        for (j in 0..<nPhases) {
            if (D[0][i.toInt(), j.toInt()] < 0) {
                D[0].set(i.toInt(), j.toInt(), 0.0)
            }
            if (D[1][i.toInt(), j.toInt()] < 0) {
                D[1].set(i.toInt(), j.toInt(), 0.0)
            }
        }
    }

    for (i in 0..<nPhases) {
        var rowSum = 0.0
        for (j in 0..<nPhases) {
            if (j != i) {
                rowSum = rowSum + D[0][i.toInt(), j.toInt()]
            }
            rowSum = rowSum + D[1][i.toInt(), j.toInt()]
            D[0].set(i.toInt(), i.toInt(), -rowSum)
        }
    }
    return D
}

/**
 * Sanitizes the (D0, D1) matrices of a Markovian Arrival Process (MAP) stored in a MatrixCell.
 *
 *
 * This function applies the sanitization process to the hidden (D0) and visible (D1) transition matrices of a MAP stored in a MatrixCell.
 * It ensures that all elements are non-negative and adjusts the diagonal elements to maintain the stochastic nature of the matrices.
 * The function is a convenience method that extracts the D0 and D1 matrices from the MAP and passes them to the more detailed
 * sanitization function.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices.
 * @return A MatrixCell containing the sanitized (D0, D1) matrices.
 */

fun map_normalize(MAP: MatrixCell): MatrixCell {
    return map_normalize(MAP[0], MAP[1])
}
/**
 * MAP normalize algorithms
 */
@Suppress("unused")
class MapNormalizeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}