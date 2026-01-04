/**
 * @file Marked Markovian Arrival Process feasibility validation
 * 
 * Validates mathematical feasibility of MMAP representations including stochastic
 * properties and marking consistency. Essential for ensuring valid multiclass models.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.GlobalConstants
import jline.util.matrix.MatrixCell

/**
 * Checks the feasibility of a Markovian Arrival Process with marked arrivals (MMAP).
 *
 *
 * This method verifies whether the given MMAP is feasible, which includes checking:
 * 1. The underlying MAP (formed by matrices D0 and D1) is feasible.
 * 2. The elements of each additional event matrix (D1c) are non-negative.
 * 3. The sum of the matrices representing different event types (D1c) does not exceed the visible transition matrix (D1) in absolute value.
 *
 * @param MMAP the MatrixCell containing the transition matrices of the MMAP, with D0, D1, ..., Dc representing different types of events
 * @return true if the MMAP is feasible, false otherwise
 */
fun mmap_isfeasible(MMAP: MatrixCell): Boolean {
    val MAP = MatrixCell()
    MAP[0] = MMAP[0]
    MAP[1] = MMAP[1]
    if (!map_isfeasible(MAP)) {
        return false
    }
    val C = MMAP.size() - 2
    // elements of D1c are >= 0
    var smallest: Double
    for (c in 0..<C) {
        smallest = MMAP[2 + c].elementMin()
        if (smallest < -GlobalConstants.FineTol) return false
    }
    val S = MMAP[1]
    for (c in 0..<C) {
        S.subEq(MMAP[2 + c])
    }
    return !(S.elementMaxAbs() > GlobalConstants.FineTol)
}
/**
 * MMAP isfeasible algorithms
 */
@Suppress("unused")
class MmapIsfeasibleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}