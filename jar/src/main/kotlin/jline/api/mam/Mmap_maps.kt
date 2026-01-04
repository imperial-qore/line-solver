/**
 * @file Marked Markovian Arrival Process class decomposition
 * 
 * Extracts individual MAP processes for each marked class from MMAP representations.
 * Used for analyzing per-class behavior and comparing multiclass vs single-class models.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Extracts K Markovian Arrival Processes (MAPs) from a given MMAP, one for each class.
 *
 *
 * This method creates a MAP for each class in the MMAP by separating the transitions related to that class.
 * Each resulting MAP has its own set of transition matrices, derived from the original MMAP.
 *
 * @param MMAP the original MMAP
 * @return a map containing K MAPs, each stored in a MatrixCell
 */
fun mmap_maps(MMAP: MatrixCell): MutableMap<Int, MatrixCell> {
    val K = MMAP.size() - 2
    val Maps: MutableMap<Int, MatrixCell> = HashMap()
    for (k in 0..<K) {
        Maps[k] = MatrixCell()
        Maps[k]!![0] = MMAP[0].add(1.0, MMAP[1]).add(-1.0, MMAP[2 + k])
        Maps[k]!![1] = MMAP[2 + k]
    }
    return Maps
}
/**
 * MMAP maps algorithms
 */
@Suppress("unused")
class MmapMapsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}