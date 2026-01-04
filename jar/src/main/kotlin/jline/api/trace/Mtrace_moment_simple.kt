package jline.api.trace

import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Computes the k-th order moment of the inter-arrival time between an event
 * of class i and an event of class j, for all possible pairs of classes.
 * 
 * @param T inter-arrival times
 * @param L labels
 * @param k order of the moment
 * @return matrix where element (i,j) is the k-th order moment of the inter-arrival
 *         time between an event of class i and an event of class j
 */
fun mtrace_moment_simple(T: DoubleArray, L: IntArray, k: Int): Matrix {
    val marks = L.distinct().sorted()
    val C = marks.size
    
    val MC = Array(C) { DoubleArray(C) }
    val count = Array(C) { IntArray(C) }
    
    for (t in 1 until T.size) {
        for (i in 0 until C) {
            for (j in 0 until C) {
                if (L[t - 1] == marks[i] && L[t] == marks[j]) {
                    MC[i][j] += T[t].pow(k)
                    count[i][j]++
                }
            }
        }
    }
    
    // Normalize by counts
    for (i in 0 until C) {
        for (j in 0 until C) {
            if (count[i][j] > 0) {
                MC[i][j] /= count[i][j]
            } else {
                MC[i][j] = Double.NaN
            }
        }
    }
    
    return Matrix(MC)
}
/**
 * Mtrace Moment Simple algorithms
 */
@Suppress("unused")
class MtraceMomentSimpleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}