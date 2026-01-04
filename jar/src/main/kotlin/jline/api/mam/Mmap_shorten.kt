package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Converts an MMAP representation from M3A format to BUTools format.
 *
 *
 * In the M3A format, an MMAP is represented as a MAP followed by multiple D1 matrices for different markings (D0, D1, D1a, D1b, ...).
 * In the BUTools format, the MMAP representation skips the initial D1 matrix, directly listing the marking matrices (D0, D1a, D1b, ...).
 * This method reorders the matrices accordingly.
 *
 * @param mmap the MatrixCell containing the MMAP representation in M3A format
 * @return a MatrixCell representing the MMAP in BUTools format
 */

fun mmap_shorten(mmap: MatrixCell): MatrixCell {
    val result = MatrixCell()
    for (i in 0..<mmap.size()) {
        if (i == 0) {
            result[0] = mmap[0]
        } else if (i != 1) {
            result[i - 1] = mmap[i]
        }
    }
    return result
}
/**
 * MMAP shorten algorithms
 */
@Suppress("unused")
class MmapShortenAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}