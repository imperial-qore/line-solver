/**
 * @file Markovian Arrival Process to Markov Modulated Poisson Process conversion
 * 
 * Converts MAP representations to MMPP format for compatibility with MMPP-specific algorithms.
 * Extracts generator and rate matrices for alternative arrival process formulations.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.io.line_warning
import jline.util.matrix.Matrix

/**
 * Convert a MAP to MMPP format by extracting generator matrix Q and rate matrix LAMBDA.
 * 
 * @param MAP Input MAP as Array<Matrix> where MAP[0] = D0 and MAP[1] = D1
 * @return Pair<Matrix, Matrix> containing (Q, LAMBDA) where:
 *         Q = D0 + D1 (generator matrix)
 *         LAMBDA = D1 (rate matrix)
 */
fun map2mmpp(MAP: Array<Matrix>): Pair<Matrix, Matrix> {
    if (MAP.size != 2) {
        throw IllegalArgumentException("MAP must contain exactly 2 matrices [D0, D1]")
    }
    
    val D0 = MAP[0]
    val D1 = MAP[1]
    
    // Check if D1 is diagonal (true MMPP requirement)
    // Extract diagonal elements from D1
    val D1RowVec = DoubleArray(D1.numRows) { i -> D1[i, i] }
    val D1_diag_matrix = Matrix.diag(*D1RowVec)
    val diff = D1.sub(D1_diag_matrix)
    
    val tolerance = 1e-10
    var maxDiff = 0.0
    for (i in 0..<diff.numRows) {
        for (j in 0..<diff.numCols) {
            val absVal = kotlin.math.abs(diff[i, j])
            if (absVal > maxDiff) {
                maxDiff = absVal
            }
        }
    }
    
    if (maxDiff > tolerance) {
        line_warning("map2mmpp", "The MAP is not a MMPP, LAMBDA is not diagonal")
    }
    
    // Q = MAP{1} + MAP{2} = D0 + D1
    val Q = D0.add(D1)
    
    // LAMBDA = MAP{2} = D1
    val LAMBDA = D1.copy()
    
    return Pair(Q, LAMBDA)
}
