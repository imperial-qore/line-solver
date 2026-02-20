/**
 * @file Marked Markovian Arrival Process class hiding operations
 * 
 * Hides specified arrival classes in MMAP processes by removing observable events.
 * Used for model reduction and analyzing subsystems in multiclass arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Hides specified types of arrivals in a Markovian Arrival Process with marked arrivals (MMAP).
 *
 *
 * This method creates a new MMAP where the specified arrival types are hidden, meaning they are removed
 * from the observation but the underlying stochastic process remains the same. The hidden types are set
 * to zero matrices, effectively removing their impact from the observed data.
 *
 * @param MMAP  the original MMAP
 * @param types a matrix containing the indices of the types to be hidden
 * @return a new MMAP with the specified types hidden
 */

fun mmap_hide(MMAP: MatrixCell, types: Matrix): MatrixCell? {
    val mmap = MatrixCell()
    mmap[0] = MMAP[0]
    mmap[1] = MMAP[1]
    for (i in 0..<types.length()) {
        mmap[(2 + types[i]).toInt()] = Matrix(MMAP[0].numRows, MMAP[0].numRows)
    }

    return mmap_normalize(mmap)
}
/**
 * MMAP hide algorithms
 */
@Suppress("unused")
class MmapHideAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}