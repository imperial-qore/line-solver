/**
 * @file Marked Markovian Arrival Process mixture modeling
 * 
 * Creates probabilistic mixtures of MMAP processes with specified weights.
 * Essential for modeling heterogeneous multiclass traffic patterns and aggregation.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Creates a mixture of MMAPs using the given weights (alpha) and MAPs.
 *
 *
 * This method combines multiple MAPs into a single MMAP, weighted by the vector `alpha`.
 * Each MAP is assumed to have its own transition matrices, and the combination is done
 * in a way that preserves the stochastic properties of the original processes.
 *
 * @param alpha a matrix (vector) representing the weights for each MAP in the mixture
 * @param MAPs  a map containing the individual MAPs to be combined, indexed by integer keys
 * @return a MatrixCell representing the combined MMAP
 */
fun mmap_mixture(alpha: Matrix?, MAPs: MutableMap<Int, MatrixCell>): MatrixCell {
    val Dk = MatrixCell()
    val I = MAPs.size
    
    // Initialize all matrices
    for (j in 0..<I + 2) {
        Dk[j] = Matrix(0, 0)
    }
    
    // Replace empty MAPs with exponential
    for (i in 0..<I) {
        if (MAPs[i]!!.isEmpty) {
            MAPs[i] = map_exponential(1e6)
        }
    }
    
    // Build D0 as block diagonal matrix
    for (i in 0..<I) {
        if (i == 0) {
            Dk[0] = MAPs[i]!![0].copy()
        } else {
            Dk[0] = Dk[0].createBlockDiagonal(MAPs[i]!![0])
        }
    }
    
    // Build arrival matrices
    for (i in 0..<I) {
        val mapI = MAPs[i]!!
        val D0i = mapI[0]
        val D1i_base = mapI[1]
        val numStatesI = D0i.numRows
        val e = Matrix.ones(numStatesI, 1)
        
        // Build D1i for MAP i
        var D1i = Matrix(numStatesI, 0)
        
        for (j in 0..<I) {
            val mapJ = MAPs[j]!!
            val pieJ = map_pie(mapJ)
            val numStatesJ = mapJ[0].numRows
            
            // alpha(j) * D1i_base * e * pie(MAP_j)
            val term = D1i_base.mult(e).mult(pieJ).scale(alpha!![j])
            D1i = Matrix.concatColumns(D1i, term, null)
        }
        
        // Add D1i to the total arrival matrix D1 (index 1)
        if (i == 0) {
            Dk[1] = D1i.copy()
        } else {
            Dk[1] = Matrix.concatRows(Dk[1], D1i, null)
        }
        
        // Add to class-specific matrices D{2+j}
        for (j in 0..<I) {
            if (i == j) {
                // Same class: use D1i
                if (i == 0) {
                    Dk[2 + j] = D1i.copy()
                } else {
                    Dk[2 + j] = Matrix.concatRows(Dk[2 + j], D1i, null)
                }
            } else {
                // Different class: use zeros
                val zeroMatrix = Matrix(D1i.numRows, D1i.numCols)
                zeroMatrix.zero()
                if (i == 0) {
                    Dk[2 + j] = zeroMatrix
                } else {
                    Dk[2 + j] = Matrix.concatRows(Dk[2 + j], zeroMatrix, null)
                }
            }
        }
    }
    
    return mmap_normalize(Dk)!!
}
/**
 * MMAP mixture algorithms
 */
@Suppress("unused")
class MmapMixtureAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}