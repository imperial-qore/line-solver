package jline.api.trace

/**
 * Computes the empirical probability of observing a specific 3-element
 * sequence of events, i.e. the two-step class transition probabilities.
 * 
 * @param T the inter-arrival times
 * @param L the event labels
 * @return 3D array where (i,j,h) element is the probability of observing
 *         an event of class i followed by class j followed by class h
 */
fun mtrace_sigma2(T: DoubleArray, L: IntArray): Array<Array<DoubleArray>> {
    val marks = L.distinct().sorted()
    val C = marks.size
    
    val sigma = Array(C) { Array(C) { DoubleArray(C) } }
    
    for (i in 0 until C) {
        for (j in 0 until C) {
            for (h in 0 until C) {
                var count = 0
                for (t in 0 until L.size - 2) {
                    if (L[t] == marks[i] && L[t + 1] == marks[j] && L[t + 2] == marks[h]) {
                        count++
                    }
                }
                sigma[i][j][h] = count.toDouble() / (L.size - 2)
            }
        }
    }
    
    return sigma
}
/**
 * Mtrace Sigma2 algorithms
 */
@Suppress("unused")
class MtraceSigma2Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}