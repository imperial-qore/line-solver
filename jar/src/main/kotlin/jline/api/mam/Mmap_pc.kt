package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the proportion of counts (PC) for each type in a Markovian Arrival Process with marked arrivals (MMAP).
 *
 *
 * This method calculates the proportion of arrivals attributed to each type in the MMAP. It does so by computing the
 * stationary distribution of the underlying Markov chain and using it to weight the arrivals from each type.
 *
 * @param MMAP the MMAP from which to compute the proportions
 * @return a matrix where each element represents the proportion of counts for a type
 */
fun mmap_pc(MMAP: MatrixCell): Matrix {
    val m = MMAP.size() - 2

    val neg_D0 = MMAP[0].copy()
    neg_D0.scaleEq(-1.0)
    val PC = Matrix(m, 1, m)
    for (i in 0..<m) {
        PC[i] = map_pie(MMAP[0], MMAP[1]).mult(neg_D0.inv().mult(MMAP[2 + i])).elementSum()
    }
    return PC
}
/**
 * MMAP pc algorithms
 */
@Suppress("unused")
class MmapPcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}