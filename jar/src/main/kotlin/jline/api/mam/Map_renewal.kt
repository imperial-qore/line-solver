/**
 * @file Markovian Arrival Process renewal process construction
 * 
 * Creates renewal MAP by removing correlations to obtain memoryless arrival processes.
 * Used for comparison analysis and studying impact of correlation in arrival streams.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Creates a renewal MAP by removing all correlations from the input MAP.
 *
 * This function takes a MAP and returns a renewal MAP process with the same cumulative
 * distribution function (CDF) as the input MAP, but with no correlations between
 * inter-arrival times. The resulting process is a renewal process.
 *
 * @param MAPIN The input MAP stored in a MatrixCell, containing the (D0, D1) matrices
 * @return A MatrixCell representing the renewal MAP with no correlations
 */
fun map_renewal(MAPIN: MatrixCell): MatrixCell {
    return map_renewal(MAPIN[0], MAPIN[1])
}

/**
 * Creates a renewal MAP by removing all correlations from the input MAP.
 *
 * This function takes a MAP defined by D0 and D1 matrices and returns a renewal MAP process
 * with the same cumulative distribution function (CDF) as the input MAP, but with no
 * correlations between inter-arrival times.
 *
 * @param D0 The hidden transition matrix of the input MAP
 * @param D1 The visible transition matrix of the input MAP
 * @return A MatrixCell representing the renewal MAP with no correlations
 */
fun map_renewal(D0: Matrix, D1: Matrix): MatrixCell {
    // Compute the pie vector (steady-state probabilities of the embedded DTMC)
    val pie = map_pie(D0, D1)
    
    // Create the output MAP
    val MAPOUT = MatrixCell()
    
    // D0 remains the same
    MAPOUT[0] = D0.copy()
    
    // D1 = D1 * ones(size(D1,1), 1) * pie
    val ones = Matrix.ones(D1.numRows, 1)
    val D1_ones = D1.mult(ones)
    val D1_new = D1_ones.mult(pie)
    
    MAPOUT[1] = D1_new
    
    return MAPOUT
}
/**
 * MAP renewal algorithms
 */
@Suppress("unused")
class MapRenewal {
    companion object {
        // Class documentation marker for Dokka
    }
}