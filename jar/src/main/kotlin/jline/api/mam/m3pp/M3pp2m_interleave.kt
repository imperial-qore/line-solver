/**
 * @file M3PP(2,m) interleaved MMAP construction
 * 
 * Implements interleaved superposition of multiple M3PP(2,m) processes to construct
 * a single MMAP. Performs matrix tensor operations and state space aggregation
 * for combining independent multi-class arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix
import kotlin.math.*

/**
 * Computes the interleaved MMAP obtained by multiple M3PP(2,m).
 * 
 * @param m3pps List of M3PP(2,m) processes to interleave
 * @return Interleaved MMAP
 */
fun m3pp2m_interleave(m3pps: List<MatrixCell>): MatrixCell {
    require(m3pps.isNotEmpty()) { "Cannot interleave empty list of M3PPs" }
    
    if (m3pps.size == 1) {
        return m3pps[0]
    }
    
    val L = m3pps.size
    
    // Check if any M3PP uses symbolic computation (not supported in Kotlin)
    val symbolic = false // Kotlin doesn't have symbolic computation like MATLAB
    
    // Compute transition rates r(1,L) and r(2,1:L)
    val r = Array(2) { DoubleArray(L) }
    
    // Compute r(1,i) for i = L down to 1
    r[0][L - 1] = m3pps[L - 1][0][0, 1] // D0(1,2) of last M3PP
    for (i in (L - 2) downTo 0) {
        r[0][i] = m3pps[i][0][0, 1] // D0(1,2) of i-th M3PP
        for (j in (i + 1) until L) {
            r[0][i] -= r[0][j] // Subtract already assigned rates
        }
    }
    
    // Compute r(2,i) for i = 1 to L
    r[1][0] = m3pps[0][0][1, 0] // D0(2,1) of first M3PP
    for (i in 1 until L) {
        r[1][i] = m3pps[i][0][1, 0] // D0(2,1) of i-th M3PP
        for (j in 0 until i) {
            r[1][i] -= r[1][j] // Subtract already assigned rates
        }
    }
    
    // Compute total number of class matrices M
    var M = 0
    for (i in 0 until L) {
        M += m3pps[i].size() - 2 // Subtract D0 and D1
    }
    
    // Create interleaved MMAP with n = 2 + (L-1) states
    val n = 2 + (L - 1)
    val interleavedMmap = MatrixCell(2 + M)
    
    // Initialize D0 matrix
    val D0 = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            if (j > i) {
                // Upper triangular: transitions from state i to state j
                D0[i, j] = r[0][j - 1]
            } else if (j < i) {
                // Lower triangular: transitions from state i to state j
                D0[i, j] = r[1][j]
            }
            // Diagonal elements will be set later
        }
    }
    interleavedMmap[0] = D0
    
    // Create class matrices D1c
    var classIndex = 0
    for (i in 0 until L) { // For each M3PP
        val currentM3pp = m3pps[i]
        val numClasses = currentM3pp.size() - 2
        
        for (j in 0 until numClasses) { // For each class in the i-th M3PP
            val Dic = Matrix(n, n)
            
            for (h in 0 until n) {
                // Set diagonal rates based on state
                if (h <= i) {
                    // States 0 to i: use rate from first component of M3PP
                    Dic[h, h] = currentM3pp[2 + j][0, 0]
                } else {
                    // States i+1 to n-1: use rate from second component of M3PP
                    Dic[h, h] = currentM3pp[2 + j][1, 1]
                }
            }
            
            interleavedMmap[2 + classIndex] = Dic
            classIndex++
        }
    }
    
    // Compute total D1 matrix
    val D1 = Matrix(n, n)
    for (i in 0 until M) {
        val classMatrix = interleavedMmap[2 + i]
        for (row in 0 until n) {
            for (col in 0 until n) {
                D1[row, col] += classMatrix[row, col]
            }
        }
    }
    interleavedMmap[1] = D1
    
    // Set diagonal elements of D0 to ensure stochastic property
    for (h in 0 until n) {
        var rowSum = 0.0
        
        // Sum off-diagonal elements of D0
        for (j in 0 until n) {
            if (h != j) {
                rowSum += D0[h, j]
            }
        }
        
        // Sum diagonal elements of D1 (arrival rates)
        rowSum += D1[h, h]
        
        // Set diagonal of D0 to make row sum zero
        D0[h, h] = -rowSum
    }
    
    return interleavedMmap
}

/**
 * Extension function to copy a Matrix
 */
private fun Matrix.copy(): Matrix {
    val copy = Matrix(this.getNumRows(), this.getNumCols())
    for (i in 0 until this.getNumRows()) {
        for (j in 0 until this.getNumCols()) {
            copy[i, j] = this[i, j]
        }
    }
    return copy
}
/**
 * M3Pp2M Interleave algorithms
 */
@Suppress("unused")
class M3pp2mInterleaveAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}