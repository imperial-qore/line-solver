/**
 * @file Markovian Arrival Process summation operations
 * 
 * Computes MAP representations of sums of identical MAP processes for load scaling.
 * Essential for modeling parallel systems and aggregated arrival streams.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the Markovian Arrival Process (MAP) representing the sum of `n` identical MAPs.
 *
 *
 * This method constructs a new MAP with transition matrices `D0` and `D1` by summing `n` identical MAPs
 * specified by the input `MAP`. The result is a larger MAP with `n` times the number of states in the original MAP,
 * where transitions in the resulting MAP are structured to reflect the cumulative behavior of `n` MAPs.
 *
 * @param MAP the original MAP represented by a MatrixCell containing matrices `D0` and `D1`
 * @param n   the number of identical MAPs to sum
 * @return a MatrixCell containing the transition matrices `D0` and `D1` of the resulting summed MAP
 */
fun map_sum(MAP: MatrixCell, n: Int): MatrixCell {
    val order = MAP[0].numRows
    val D0 = Matrix.zeros(n * order, n * order)
    val D1 = Matrix.zeros(n * order, n * order)
    var curpos = 0
    for (i in 0..<n) {
        D0.setSliceEq(curpos, curpos + order, curpos, curpos + order, MAP[0])
        if (i < n - 1) {
            D0.setSliceEq(curpos, curpos + order, curpos + order, curpos + 2 * order, MAP[1])
        } else {
            D1.setSliceEq(curpos, curpos + order, 0, order, MAP[1])
        }
        curpos += order
    }
    return MatrixCell(D0, D1)
}
/**
 * MAP sum algorithms
 */
@Suppress("unused")
class MapSumAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}