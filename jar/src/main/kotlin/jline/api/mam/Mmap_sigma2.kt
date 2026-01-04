package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes two-step class transition probabilities for a Markovian Arrival Process (MMAP).
 * 
 * This function calculates the 3D matrix of class-transition probabilities:
 * p_{i,j,h} = P(C_k = h | C_{k-1} = j , C_{k-2} = i)
 * 
 * The computation follows the algorithm from the MATLAB implementation:
 * 1. Compute the steady-state probability vector using map_pie
 * 2. For each class combination (i,j,h), compute the conditional probability
 * 3. Use matrix inversions to compute the transition probabilities
 * 
 * @param mmap the MMAP represented as a MatrixCell containing {D0, D1, D2, ..., Dc}
 * @return a 3D array (Matrix[][][]) of class-transition probabilities
 */
fun mmap_sigma2(mmap: MatrixCell): Array<Array<Array<Double>>> {
    val C = mmap.size() - 2  // Number of classes (excluding D0 and D1)
    
    // Initialize the 3D sigma array
    val sigma = Array(C) { Array(C) { Array(C) { 0.0 } } }
    
    // Compute steady-state probability vector
    val alpha = map_pie(mmap)
    
    // Compute -D0^{-1} once for efficiency
    val negD0inv = mmap[0].scale(-1.0).inv()
    
    // Compute class transition probabilities
    for (i in 0 until C) {
        // starti = alpha * (-D0^{-1} * D_{i+2})
        val starti = alpha.mult(negD0inv).mult(mmap[i + 2])
        
        for (j in 0 until C) {
            // startj = starti * (-D0^{-1} * D_{j+2})
            val startj = starti.mult(negD0inv).mult(mmap[j + 2])
            
            for (h in 0 until C) {
                // sigma(i,j,h) = sum(startj * (-D0^{-1} * D_{h+2}))
                val result = startj.mult(negD0inv).mult(mmap[h + 2])
                // Sum all elements in the matrix
                var sum = 0.0
                for (row in 0 until result.numRows) {
                    for (col in 0 until result.numCols) {
                        sum += result[row, col]
                    }
                }
                sigma[i][j][h] = sum
            }
        }
    }
    
    return sigma
}

/**
 * Computes two-step class transition probabilities for a Markovian Arrival Process (MMAP).
 * 
 * This overload returns the result as a MatrixCell where each matrix represents
 * the transition probabilities for a specific initial class.
 * 
 * @param mmap the MMAP represented as a MatrixCell containing {D0, D1, D2, ..., Dc}
 * @return a MatrixCell containing the class-transition probability matrices
 */
fun mmap_sigma2_cell(mmap: MatrixCell): MatrixCell {
    val C = mmap.size() - 2  // Number of classes
    val sigma = mmap_sigma2(mmap)
    
    val result = MatrixCell(C)
    
    for (i in 0 until C) {
        val matrix = Matrix(C, C)
        for (j in 0 until C) {
            for (h in 0 until C) {
                matrix.set(j, h, sigma[i][j][h])
            }
        }
        result.set(i, matrix)
    }
    
    return result
}
/**
 * MMAP sigma2 algorithms
 */
@Suppress("unused")
class MmapSigma2Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}