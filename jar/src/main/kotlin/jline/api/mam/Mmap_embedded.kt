/**
 * @file Marked Markovian Arrival Process embedded chain analysis
 * 
 * Computes embedded discrete-time Markov chain for MMAP processes.
 * Essential for analyzing state transition probabilities at marked arrival epochs.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell
import jline.util.matrix.Matrix

/**
 * Computes the embedded chain of an MMAP.
 * The embedded chain represents the transition probabilities between states
 * when arrivals occur.
 *
 * @param mmap The MMAP process
 * @return Embedded chain transition matrix
 */
fun mmap_embedded(mmap: MatrixCell): Matrix {
    val D0 = mmap[0]
    val order = D0.getNumRows()
    
    // Compute total arrival rates for each state
    val totalArrivalRates = DoubleArray(order)
    for (i in 1 until mmap.size()) {
        val Di = mmap[i]
        for (row in 0 until order) {
            for (col in 0 until order) {
                totalArrivalRates[row] += Di[row, col]
            }
        }
    }
    
    // Create embedded transition matrix
    val P = Matrix(order, order)
    
    for (i in 0 until order) {
        if (totalArrivalRates[i] > 1e-12) {
            // For each state with positive arrival rate
            for (j in 0 until order) {
                var transitionRate = 0.0
                
                // Sum arrival rates from state i to state j across all classes
                for (c in 1 until mmap.size()) {
                    val Dc = mmap[c]
                    transitionRate += Dc[i, j]
                }
                
                // Normalize by total arrival rate to get probability
                P[i, j] = transitionRate / totalArrivalRates[i]
            }
        } else {
            // If no arrivals from this state, use uniform distribution
            for (j in 0 until order) {
                P[i, j] = 1.0 / order
            }
        }
    }
    
    return P
}
/**
 * MMAP embedded algorithms
 */
@Suppress("unused")
class MmapEmbeddedAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}