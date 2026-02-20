/**
 * @file Marked Markovian Arrival Process normalization and sanitization
 * 
 * Normalizes MMAP matrices to ensure feasibility and mathematical validity.
 * Essential for maintaining proper stochastic properties in multiclass models.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * Normalizes a Markovian Arrival Process with marked arrivals (MMAP) to ensure feasibility.
 *
 *
 * This method adjusts the MMAP by setting negative off-diagonal values in the D0 matrix to zero and ensuring that
 * all elements in the marking matrices are non-negative. It also recalculates the D1 matrix to reflect the sum of
 * the marking matrices. The diagonal elements of the D0 matrix are adjusted to ensure that each row sums to zero,
 * maintaining the properties of a valid generator matrix.
 *
 * @param MMAP the MatrixCell representing the MMAP to be normalized
 * @return the normalized MatrixCell, or null if the input MMAP is empty
 */

fun mmap_normalize(MMAP: MatrixCell): MatrixCell? {
    if (MMAP.isEmpty) {
        return null
    }
    val K = MMAP[0].numRows
    val C = MMAP.size() - 2

    for (i in 0..<K) {
        for (j in 0..<K) {
            if (i != j) {
                MMAP[0][i, j] = FastMath.max(MMAP[0][i, j], 0.0)
            }
        }
    }

    MMAP[1] = Matrix(MMAP[0].numRows, MMAP[0].numCols, MMAP[0].numRows * MMAP[0].numCols)

    for (c in 0..<C) {
        MMAP[2 + c].removeNegative()
        if (java.lang.Double.isNaN(MMAP[2 + c][0])) {
            MMAP[2 + c] = Matrix(MMAP[2 + c].numRows, MMAP[2 + c].numCols, MMAP[2 + c].numRows * MMAP[2 + c].numCols)
        }
        MMAP[1] = MMAP[1].add(1.0, MMAP[2 + c])
    }

    for (k in 0..<K) {
        MMAP[0][k, k] = 0
        MMAP[0][k, k] = -MMAP[0].sumRows(k) - MMAP[1].sumRows(k)
    }

    return MMAP
}
/**
 * MMAP normalize algorithms
 */
@Suppress("unused")
class MmapNormalizeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}