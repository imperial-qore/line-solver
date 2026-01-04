/**
 * @file Markovian Arrival Process stochastic complementation
 * 
 * Performs state elimination through stochastic complementation while preserving MAP properties.
 * Used for model reduction and aggregation in large-scale MAP analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Performs stochastic complementation on a MAP by eliminating specified states.
 *
 * Stochastic complementation is a technique used to reduce the dimension of a MAP
 * by eliminating certain states while preserving the statistical properties of the
 * remaining states. This is useful for model reduction and analysis.
 *
 * @param MAP The input MAP stored in a MatrixCell, containing the (D0, D1) matrices
 * @param retainIdx Array of indices of states to retain (0-based indexing)
 * @return A MatrixCell representing the reduced MAP after stochastic complementation
 */
fun map_stochcomp(MAP: MatrixCell, retainIdx: IntArray): MatrixCell {
    return map_stochcomp(MAP[0], MAP[1], retainIdx)
}

/**
 * Performs stochastic complementation on a MAP by eliminating specified states.
 *
 * @param D0 The hidden transition matrix of the input MAP
 * @param D1 The visible transition matrix of the input MAP
 * @param retainIdx Array of indices of states to retain (0-based indexing)
 * @return A MatrixCell representing the reduced MAP after stochastic complementation
 */
fun map_stochcomp(D0: Matrix, D1: Matrix, retainIdx: IntArray): MatrixCell {
    val n = D0.numRows
    
    // Validate retain indices
    for (idx in retainIdx) {
        if (idx < 0 || idx >= n) {
            throw IllegalArgumentException("Retain index $idx is out of bounds [0, ${n-1}]")
        }
    }
    
    // Create the full Q matrix
    val Q = D0.add(1.0, D1)
    
    // Find eliminated indices
    val eliminatedIdx = (0 until n).filter { !retainIdx.contains(it) }.toIntArray()
    
    if (eliminatedIdx.isEmpty()) {
        // Nothing to eliminate, return original MAP
        return map_normalize(D0, D1)
    }
    
    // Extract submatrices
    val Q_RE = extractSubmatrix(Q, retainIdx, eliminatedIdx)
    val Q_EE = extractSubmatrix(Q, eliminatedIdx, eliminatedIdx)
    val Q_RR = extractSubmatrix(Q, retainIdx, retainIdx)
    val Q_ER = extractSubmatrix(Q, eliminatedIdx, retainIdx)
    
    // Compute the new Q matrix using stochastic complementation
    // QNew = Q_RR + Q_RE * (-Q_EE)^(-1) * Q_ER
    val minusQ_EE = Q_EE.scale(-1.0)
    val invMinusQ_EE = minusQ_EE.inv()
    val QNew = Q_RR.add(1.0, Q_RE.mult(invMinusQ_EE).mult(Q_ER))
    
    // Extract D1 submatrices
    val D1_RR = extractSubmatrix(D1, retainIdx, retainIdx)
    val D1_ER = extractSubmatrix(D1, eliminatedIdx, retainIdx)
    
    // Compute new D0 and D1 matrices
    val D0new = QNew.sub(1.0, D1_RR)
    val D1new = D1_RR.add(1.0, Q_RE.mult(invMinusQ_EE).mult(D1_ER))
    
    return map_normalize(D0new, D1new)
}

/**
 * Extracts a submatrix from the given matrix based on row and column indices.
 */
private fun extractSubmatrix(matrix: Matrix, rowIndices: IntArray, colIndices: IntArray): Matrix {
    val result = Matrix(rowIndices.size, colIndices.size)
    
    for (i in rowIndices.indices) {
        for (j in colIndices.indices) {
            result.set(i, j, matrix.get(rowIndices[i], colIndices[j]))
        }
    }
    
    return result
}
/**
 * MAP stochcomp algorithms
 */
@Suppress("unused")
class MapStochcompAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}