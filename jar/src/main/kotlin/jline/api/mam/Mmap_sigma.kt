package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes one-step class transition probabilities for a Marked Markovian Arrival Process (MMAP).
 * 
 * The function computes the class transition probabilities p_{i,j} = P(C_k = j | C_{k-1} = i),
 * which represent the probability of transitioning from class i to class j in one step.
 *
 * For an MMAP with m classes, the algorithm computes:
 * sigma(i,j) = alpha * P_i * P_j * 1
 * where:
 * - alpha is the stationary probability vector of the underlying MAP
 * - P_i = (-D0)^{-1} * D_i is the probability of selecting class i given a transition
 *
 * @param MMAP the MMAP represented as a MatrixCell where MMAP[0] = D0, MMAP[1] = aggregate D1, 
 *             and MMAP[2+i] = D_{i+1} for i = 0, 1, ..., m-1
 * @return matrix of class transition probabilities where sigma(i,j) is the probability 
 *         of transitioning from class i to class j
 */
fun mmap_sigma(MMAP: MatrixCell): Matrix {
    val C = MMAP.size() - 2  // Number of classes
    val sigma = Matrix.zeros(C, C)
    
    // Get the stationary probability vector of the underlying MAP
    val alpha = map_pie(MMAP[0], MMAP[1])
    
    // Compute (-D0)^{-1} once for efficiency
    val invNegD0 = MMAP[0].mult(Matrix.singleton(-1.0)).inv()
    
    for (i in 0 until C) {
        // Compute alpha * (-D0)^{-1} * D_{i+1}
        val start = alpha.mult(invNegD0).mult(MMAP[2 + i])
        
        for (j in 0 until C) {
            // Compute start * (-D0)^{-1} * D_{j+1} * 1
            val result = start.mult(invNegD0).mult(MMAP[2 + j])
            sigma.set(i, j, result.elementSum())
        }
    }
    
    return sigma
}

/**
 * Computes one-step class transition probabilities for a Marked Markovian Arrival Process (MMAP).
 * Array-based overload for compatibility with functions expecting Array<Matrix>.
 *
 * @param mmap the MMAP represented as Array<Matrix> where mmap[0] = D0, mmap[1] = D1, 
 *             mmap[2] = D2, ..., mmap[m+1] = D_m for m classes
 * @return matrix of class transition probabilities where sigma(i,j) is the probability 
 *         of transitioning from class i to class j
 */
fun mmap_sigma(mmap: Array<Matrix>): Matrix {
    // Convert to MatrixCell and delegate to main implementation
    val mmapCell = MatrixCell(mmap.size)
    for (i in mmap.indices) {
        mmapCell[i] = mmap[i]
    }
    return mmap_sigma(mmapCell)
}
/**
 * MMAP sigma algorithms
 */
@Suppress("unused")
class MmapSigmaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}