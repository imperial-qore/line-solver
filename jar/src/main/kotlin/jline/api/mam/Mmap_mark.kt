package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Converts a Markovian Arrival Process with marked arrivals (MMAP) into a new MMAP with redefined classes based on a given probability matrix.
 *
 *
 * This method reclassifies the arrivals in an MMAP according to a provided probability matrix. The probability matrix defines
 * the likelihood that an arrival of a particular type in the original MMAP is marked as a different type in the new MMAP.
 * The resulting MMAP has a new set of output types (classes) based on these probabilities.
 *
 * @param MMAP the original MMAP with K types
 * @param prob a KxR matrix describing the probability of a type-k arrival in the original MMAP being marked as a type-r arrival in the new MMAP
 * @return a new MMAP with R types (classes)
 */

fun mmap_mark(MMAP: MatrixCell, prob: Matrix): MatrixCell {
    val K = prob.numRows
    val R = prob.numCols
    val mmap = MatrixCell(2 + R)
    mmap[0] = MMAP[0].copy()
    mmap[1] = MMAP[1].copy()
    for (r in 0..<R) {
        mmap[2 + r] = Matrix(MMAP[0].length(), MMAP[0].length(), MMAP[0].length() * MMAP[0].length())
        for (k in 0..<K) {
            val a = MMAP[2 + k].copy()
            a.scaleEq(prob[k, r])
            mmap[2 + r] = mmap[2 + r].add(1.0, a)
        }
    }

    return mmap
}
/**
 * MMAP mark algorithms
 */
@Suppress("unused")
class MmapMark {
    companion object {
        // Class documentation marker for Dokka
    }
}