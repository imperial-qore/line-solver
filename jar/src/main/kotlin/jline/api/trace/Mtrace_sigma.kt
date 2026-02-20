package jline.api.trace

import jline.util.matrix.Matrix

/**
 * Computes the empirical probability of observing a specific 2-element
 * sequence of events, i.e. the one-step class transition probabilities.
 * 
 * @param T the inter-arrival times
 * @param L the event labels
 * @return matrix where (i,j) element is the probability of observing
 *         an event of class i followed by an event of class j
 */
fun mtrace_sigma(T: DoubleArray, L: IntArray): Matrix {
    val marks = L.distinct().sorted()
    val C = marks.size
    
    val sigma = Array(C) { DoubleArray(C) }
    
    for (i in 0 until C) {
        for (j in 0 until C) {
            var count = 0
            for (t in 0 until L.size - 1) {
                if (L[t] == marks[i] && L[t + 1] == marks[j]) {
                    count++
                }
            }
            sigma[i][j] = count.toDouble() / (L.size - 1)
        }
    }
    
    return Matrix(sigma)
}
/**
 * Mtrace Sigma algorithms
 */
@Suppress("unused")
class MtraceSigmaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}