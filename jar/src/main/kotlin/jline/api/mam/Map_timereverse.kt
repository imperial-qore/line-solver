/**
 * @file Markovian Arrival Process time reversal operations
 * 
 * Computes time-reversed MAP by adjusting transition rates based on stationary distributions.
 * Used for analyzing reversibility properties and theoretical MAP characterization.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the time-reversed MAP of a given MAP.
 *
 *
 * The time-reversed MAP has transition rates adjusted based on the stationary distribution of the original MAP.
 *
 * @param map the original MAP stored in a MatrixCell containing the D0 and D1 matrices
 * @return a MatrixCell representing the time-reversed MAP, containing the reversed D0 and D1 matrices
 */

fun map_timereverse(map: MatrixCell): MatrixCell {
    val piq = map_piq(map[0], map[1])
    val D = Matrix.diag(*piq.toArray1D())
    val iD = D.inv()
    val MAPr0 = iD.mult(map[0].transpose()).mult(D)
    val MAPr1 = iD.mult(map[1].transpose()).mult(D)
    return MatrixCell(MAPr0, MAPr1)
}
/**
 * MAP timereverse algorithms
 */
@Suppress("unused")
class MapTimereverseAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}