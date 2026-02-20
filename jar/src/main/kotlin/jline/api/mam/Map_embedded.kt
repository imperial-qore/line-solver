/**
 * @file Markovian Arrival Process embedded discrete-time Markov chain analysis
 * 
 * Computes embedded DTMC matrices from MAP representations by extracting transition
 * probabilities at arrival epochs. Used for state transition analysis in MAP models.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the embedded discrete-time Markov chain (DTMC) matrix of a MAP.
 *
 *
 * This method calculates the embedded DTMC matrix by first inverting the negative of the D0 matrix and then multiplying it by the D1 matrix.
 * The resulting matrix represents the probabilities of transitioning between states in the embedded chain.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @return the embedded DTMC matrix
 */

fun map_embedded(D0: Matrix, D1: Matrix): Matrix {
    val neg_D0 = D0.copy()
    neg_D0.scaleEq(-1.0)
    return neg_D0.inv().mult(D1)
}

/**
 * Computes the embedded discrete-time Markov chain (DTMC) matrix of a MAP given as a MatrixCell.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the embedded DTMC matrix.
 *
 * @param MAP the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @return the embedded DTMC matrix
 */
fun map_embedded(MAP: MatrixCell): Matrix {
    return map_embedded(MAP[0], MAP[1])
}
/**
 * MAP embedded algorithms
 */
@Suppress("unused")
class MapEmbeddedAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}