/**
 * @file Markovian Arrival Process feasibility checking interface
 * 
 * Provides convenient interface for MAP feasibility validation with configurable tolerance.
 * Wrapper around detailed feasibility checking algorithms for ease of use.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * Checks if the provided MAP is feasible based on the given tolerance.
 *
 *
 * This method evaluates whether the MAP transition matrices meet the necessary conditions for a valid MAP,
 * such as non-negativity of elements, proper row/column sums, etc.
 *
 * @param MAP the MatrixCell representing the MAP transition matrices
 * @param TOL the tolerance level for numerical stability checks
 * @return true if the MAP is feasible, false otherwise
 */
fun map_isfeasible(MAP: MatrixCell, TOL: Double): Boolean {
    return map_checkfeasible(MAP, TOL)
}

/**
 * Checks if the provided MAP is feasible using a default tolerance.
 *
 *
 * This method uses a standard tolerance level to determine the feasibility of the MAP.
 *
 * @param MAP the MatrixCell representing the MAP transition matrices
 * @return true if the MAP is feasible, false otherwise
 */

fun map_isfeasible(MAP: MatrixCell): Boolean {
    val TOLMAGNITUDE = 15
    for (k in TOLMAGNITUDE downTo 1) {
        val check = map_checkfeasible(MAP, FastMath.pow(10.0, -k))
        if (check) {
            return k > map_feastol()
        }
    }
    return false
}
/**
 * MAP isfeasible algorithms
 */
@Suppress("unused")
class MapIsfeasibleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}