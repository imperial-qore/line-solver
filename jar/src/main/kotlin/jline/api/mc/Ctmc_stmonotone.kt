package jline.api.mc

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Computes the stochastically monotone upper bound for a CTMC
 * Based on Forneau, Pekergin - An Algorithmic Approach to Stochastic Bounds, Performance 2002
 *
 * @param Q infinitesimal generator matrix
 * @return Qub stochastically monotone upper bound infinitesimal generator
 */
fun ctmc_stmonotone(Q: Matrix): Matrix {
    // Find the maximum absolute value in Q for randomization
    var maxAbsVal = 0.0
    for (i in 0 until Q.getNumRows()) {
        for (j in 0 until Q.getNumCols()) {
            maxAbsVal = max(maxAbsVal, abs(Q.get(i, j)))
        }
    }
    
    // Apply randomization to convert CTMC to DTMC
    val (P, _) = ctmc_randomization(Q, maxAbsVal)
    
    // Make sure P is stochastic
    val P_stochastic = dtmc_makestochastic(P)
    
    // Apply the stochastically monotone algorithm to the DTMC
    val Pub = dtmc_stmonotone(P_stochastic)
    
    // Convert back to infinitesimal generator
    val Qub = ctmc_makeinfgen(Pub)
    
    return Qub
}

/**
 * Implementation of the dtmc_stmonotone algorithm
 * Forneau, Pekergin - An Algorithmic Approach to Stochastic Bounds, Performance 2002 
 * Algorithm 1
 *
 * @param P transition probability matrix
 * @return Q stochastically monotone upper bound matrix
 */
fun dtmc_stmonotone(P: Matrix): Matrix {
    val n = P.length()
    val Q = Matrix(n, n)
    
    // Initialize last column
    Q.set(0, n - 1, P.get(0, n - 1))
    
    for (i in 1 until n) {
        Q.set(i, n - 1, max(Q.get(i - 1, n - 1), P.get(i, n - 1)))
    }
    
    // Fill remaining columns from right to left
    for (l in n - 2 downTo 0) {
        Q.set(0, l, P.get(0, l))
        
        for (i in 1 until n) {
            // Calculate sum(Q(i-1, l:n))
            var sumQ_prev = 0.0
            for (k in l until n) {
                sumQ_prev += Q.get(i - 1, k)
            }
            
            // Calculate sum(P(i, l:n))
            var sumP = 0.0
            for (k in l until n) {
                sumP += P.get(i, k)
            }
            
            // Calculate sum(Q(i, (l+1):n))
            var sumQ_right = 0.0
            for (k in l + 1 until n) {
                sumQ_right += Q.get(i, k)
            }
            
            Q.set(i, l, max(sumQ_prev, sumP) - sumQ_right)
        }
    }
    
    return Q
}
/**
 * CTMC stmonotone algorithms
 */
@Suppress("unused")
class CtmcStmonotoneAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}