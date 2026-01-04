package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Checks if an MMAP is symmetric.
 * An MMAP is symmetric if all its matrices are symmetric.
 *
 * @param mmap The MMAP to check
 * @param tolerance Tolerance for numerical comparison (default: 1e-12)
 * @return True if the MMAP is symmetric, false otherwise
 */
fun mmap_issym(mmap: MatrixCell, tolerance: Double = 1e-12): Boolean {
    // Check each matrix in the MMAP
    for (i in 0 until mmap.size()) {
        val matrix = mmap[i]
        
        // Check if matrix is square
        if (matrix.getNumRows() != matrix.getNumCols()) {
            return false
        }
        
        // Check if matrix is symmetric
        for (row in 0 until matrix.getNumRows()) {
            for (col in 0 until matrix.getNumCols()) {
                val diff = Math.abs(matrix[row, col] - matrix[col, row])
                if (diff > tolerance) {
                    return false
                }
            }
        }
    }
    
    return true
}
/**
 * MMAP issym algorithms
 */
@Suppress("unused")
class MmapIssymAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}