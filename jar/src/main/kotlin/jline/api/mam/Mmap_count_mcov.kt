/**
 * @file Marked Markovian Arrival Process counting covariance analysis
 * 
 * Computes count covariance matrices between marked classes in MMAP processes.
 * Essential for analyzing dependencies and correlations between different arrival classes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the count covariance between each pair of classes at a given time scale.
 * 
 * For an MMAP with m classes, this function computes the m x m covariance matrix S
 * where S(i,j) represents the covariance between the count of class i events and
 * class j events over time period t.
 * 
 * The algorithm works as follows:
 * 1. For diagonal elements: S(i,i) = Var(N_i(t)) (per-class variance)
 * 2. For off-diagonal elements: S(i,j) = 1/2 * (Var(N_i(t) + N_j(t)) - Var(N_i(t)) - Var(N_j(t)))
 * 
 * This uses the property that for random variables X and Y:
 * Cov(X,Y) = 1/2 * (Var(X+Y) - Var(X) - Var(Y))
 *
 * @param MMAP the MMAP represented as a MatrixCell where MMAP[0] = D0, MMAP[1] = aggregate D1,
 *             and MMAP[2+i] = D_{i+1} for i = 0, 1, ..., m-1
 * @param t the time scale over which to compute covariances
 * @return the m x m covariance matrix between class pairs
 */
fun mmap_count_mcov(MMAP: MatrixCell, t: Double): Matrix {
    val m = MMAP.size() - 2  // Number of classes
    
    // Compute per-class variances (diagonal elements)
    val mV = mmap_count_var(MMAP, t)
    val S = Matrix.zeros(m, m)
    
    // Set diagonal elements to per-class variances
    for (i in 0 until m) {
        S.set(i, i, mV.get(0, i))
    }
    
    // Compute off-diagonal covariances
    for (i in 0 until m) {
        for (j in 0 until m) {
            if (i != j) {
                // Create temporary MMAP with combined classes i and j
                val mmap2 = MatrixCell(4)
                mmap2[0] = MMAP[0]  // D0 remains the same
                mmap2[1] = MMAP[1]  // Aggregate D1 remains the same
                mmap2[2] = MMAP[2 + i].add(MMAP[2 + j])  // Combined class matrix
                mmap2[3] = mmap2[1].sub(mmap2[2])  // Complement class matrix
                
                // Compute variance of the combined process
                val pV = mmap_count_var(mmap2, t)
                
                // Covariance = 1/2 * (Var(X+Y) - Var(X) - Var(Y))
                val covariance = 0.5 * (pV.get(0, 0) - mV.get(0, i) - mV.get(0, j))
                S.set(i, j, covariance)
            }
        }
    }
    
    return S
}

/**
 * Computes the count covariance between each pair of classes at a given time scale.
 * Array-based overload for compatibility with functions expecting Array<Matrix>.
 *
 * @param mmap the MMAP represented as Array<Matrix> where mmap[0] = D0, mmap[1] = D1, 
 *             mmap[2] = D2, ..., mmap[m+1] = D_m for m classes
 * @param t the time scale over which to compute covariances
 * @return the m x m covariance matrix between class pairs
 */
fun mmap_count_mcov(mmap: Array<Matrix>, t: Double): Matrix {
    // Convert to MatrixCell and delegate to main implementation
    val mmapCell = MatrixCell(mmap.size)
    for (i in mmap.indices) {
        mmapCell[i] = mmap[i]
    }
    return mmap_count_mcov(mmapCell, t)
}
/**
 * MMAP count mcov algorithms
 */
@Suppress("unused")
class MmapCountMcovAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}