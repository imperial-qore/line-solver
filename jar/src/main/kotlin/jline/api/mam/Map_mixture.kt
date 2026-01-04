/**
 * @file Markovian Arrival Process probabilistic mixture models
 * 
 * Creates probabilistic mixtures of MAP processes with specified mixture probabilities.
 * Used for modeling heterogeneous arrival patterns and traffic characterization.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Creates a probabilistic mixture of Markovian Arrival Processes (MAPs).
 *
 * This function creates a mixture MAP from a collection of MAPs with given mixture probabilities.
 * The resulting MAP combines the behavior of the input MAPs according to the specified probabilities.
 * Implementation follows the MATLAB map_mixture function.
 *
 * @param alpha The mixture probabilities (should sum to 1.0)
 * @param MAPs Array of MAPs to be mixed
 * @return A MatrixCell representing the mixture MAP
 */
fun map_mixture(alpha: DoubleArray, MAPs: Array<MatrixCell>): MatrixCell {
    if (alpha.size != MAPs.size) {
        throw IllegalArgumentException("Alpha array and MAPs array must have the same size")
    }
    
    // Verify that alpha sums to 1.0 (within tolerance)
    val alphaSum = alpha.sum()
    if (kotlin.math.abs(alphaSum - 1.0) > 1e-10) {
        throw IllegalArgumentException("Alpha probabilities must sum to 1.0, got $alphaSum")
    }
    
    // Pre-compute pie vectors for all MAPs
    val pies = Array(MAPs.size) { i -> map_pie(MAPs[i]) }
    
    // Build D0 as block diagonal matrix
    var D0: Matrix? = null
    for (i in MAPs.indices) {
        if (i == 0) {
            D0 = MAPs[i][0].copy()
        } else {
            D0 = blockDiag(D0!!, MAPs[i][0])
        }
    }
    
    // Build D1 matrix by vertically concatenating D1i blocks
    var D1: Matrix? = null
    for (i in MAPs.indices) {
        val mapSize = MAPs[i][0].numRows
        val ones = Matrix.ones(mapSize, 1)
        
        // Start with alpha[0] * MAPs[i]{2} * ones * pie[0]
        var D1i = MAPs[i][1].mult(ones.scale(alpha[0])).mult(pies[0])
        
        // Horizontally concatenate with alpha[j] * MAPs[i]{2} * ones * pie[j] for j=1..n-1
        for (j in 1 until MAPs.size) {
            val D1i_j = MAPs[i][1].mult(ones.scale(alpha[j])).mult(pies[j])
            D1i = Matrix.concatColumns(D1i, D1i_j, null)
        }
        
        // Vertically concatenate to build D1
        if (i == 0) {
            D1 = D1i
        } else {
            D1 = Matrix.concatRows(D1!!, D1i, null)
        }
    }
    
    return map_normalize(D0!!, D1!!)
}

/**
 * Creates a block diagonal matrix from two matrices.
 */
private fun blockDiag(A: Matrix, B: Matrix): Matrix {
    val result = Matrix(A.numRows + B.numRows, A.numCols + B.numCols)
    
    // Copy A to the top-left block
    for (i in 0 until A.numRows) {
        for (j in 0 until A.numCols) {
            result.set(i, j, A.get(i, j))
        }
    }
    
    // Copy B to the bottom-right block
    for (i in 0 until B.numRows) {
        for (j in 0 until B.numCols) {
            result.set(A.numRows + i, A.numCols + j, B.get(i, j))
        }
    }
    
    return result
}
/**
 * MAP mixture algorithms
 */
@Suppress("unused")
class MapMixture {
    companion object {
        // Class documentation marker for Dokka
    }
}