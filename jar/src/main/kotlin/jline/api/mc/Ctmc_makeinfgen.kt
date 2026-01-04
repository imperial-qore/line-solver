/**
 * @file CTMC infinitesimal generator construction and validation
 * 
 * Constructs and validates infinitesimal generator matrices for continuous-time 
 * Markov chains. Ensures row sums equal zero and non-positive diagonal elements, 
 * fundamental requirements for valid CTMC generator matrices.
 * 
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Converts a matrix into a valid infinitesimal generator for a CTMC.
 * An infinitesimal generator has row sums equal to zero and non-positive diagonal elements.
 *
 * @param Q candidate infinitesimal generator matrix
 * @return valid infinitesimal generator matrix
 */
fun ctmc_makeinfgen(Q: Matrix): Matrix {
    // Extract diagonal elements and create off-diagonal matrix
    val diagonalValues = DoubleArray(Q.length()) { i -> Q[i, i] }
    val offDiagonal = Q.sub(Matrix.diagMatrix(null, diagonalValues, 0, diagonalValues.size))
    
    // Create new diagonal from negative row sums to ensure row sums = 0
    val rowSums = offDiagonal.sumRows()
    val newDiagonal = Matrix(0, 0)
    Matrix.diagMatrix(newDiagonal, rowSums.nonZeroValues, 0, rowSums.nonZeroLength)
    
    // Combine off-diagonal elements with corrected diagonal
    val result = offDiagonal.sub(newDiagonal)
    result.removeZeros(0.0)
    return result
}
/**
 * CTMC makeinfgen algorithms
 */
@Suppress("unused")
class CtmcMakeinfgenAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}