/**
 * @file Marked Markovian Arrival Process maximum operations
 * 
 * Computes element-wise maximum of MMAP processes for synchronization analysis.
 * Used for worst-case timing analysis and parallel system modeling.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix

/**
 * Computes the element-wise maximum of two MMAPs.
 * The result has the same structure as the input MMAPs.
 *
 * @param mmap1 First MMAP
 * @param mmap2 Second MMAP  
 * @return MMAP with element-wise maximum values
 */
fun mmap_max(mmap1: MatrixCell, mmap2: MatrixCell): MatrixCell {
    require(mmap1.size() == mmap2.size()) {
        "MMAPs must have the same number of matrices"
    }
    
    val result = MatrixCell(mmap1.size())
    
    for (i in 0 until mmap1.size()) {
        val matrix1 = mmap1[i]
        val matrix2 = mmap2[i]
        
        require(matrix1.numRows == matrix2.numRows && matrix1.numCols == matrix2.numCols) {
            "Corresponding matrices must have the same dimensions"
        }
        
        val maxMatrix = Matrix(matrix1.numRows, matrix1.numCols)
        
        for (row in 0 until matrix1.numRows) {
            for (col in 0 until matrix1.numCols) {
                maxMatrix[row, col] = Math.max(matrix1[row, col], matrix2[row, col])
            }
        }
        
        result[i] = maxMatrix
    }
    
    return result
}

/**
 * Computes the element-wise maximum of multiple MMAPs.
 *
 * @param mmaps List of MMAPs
 * @return MMAP with element-wise maximum values across all inputs
 */
fun mmap_max_multiple(mmaps: List<MatrixCell>): MatrixCell {
    require(mmaps.isNotEmpty()) {
        "Must provide at least one MMAP"
    }
    
    if (mmaps.size == 1) {
        return mmaps[0]
    }
    
    var result = mmaps[0]
    for (i in 1 until mmaps.size) {
        result = mmap_max(result, mmaps[i])
    }
    
    return result
}
/**
 * MMAP max algorithms
 */
@Suppress("unused")
class MmapMaxAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}