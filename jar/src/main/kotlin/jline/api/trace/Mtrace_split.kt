package jline.api.trace

/**
 * Given a multi-class trace with inter-arrivals T and labels L,
 * creates the separate per-class traces. In each per-class trace
 * the inter-arrivals are inter-arrivals between events of the same class.
 * 
 * @param T vector of inter-arrival times
 * @param L vector of class labels
 * @return array where the c-th element contains the vector of inter-arrival
 *         times for class c
 */
fun mtrace_split(T: DoubleArray, L: IntArray): Array<DoubleArray> {
    val labels = L.distinct().sorted()
    val C = labels.size
    
    val TL = Array(C) { doubleArrayOf() }
    
    // Compute cumulative sum of T
    val TCUM = DoubleArray(T.size + 1)
    TCUM[0] = 0.0
    for (i in T.indices) {
        TCUM[i + 1] = TCUM[i] + T[i]
    }
    
    for (c in 0 until C) {
        val label = labels[c]
        
        // Find indices where L == label
        val indices = mutableListOf<Int>()
        indices.add(0) // Add initial 0
        for (i in L.indices) {
            if (L[i] == label) {
                indices.add(i + 1) // +1 because TCUM is offset by 1
            }
        }
        
        // Compute differences of cumulative times at those indices
        if (indices.size > 1) {
            val classInterArrivals = DoubleArray(indices.size - 1)
            for (i in 1 until indices.size) {
                classInterArrivals[i - 1] = TCUM[indices[i]] - TCUM[indices[i - 1]]
            }
            TL[c] = classInterArrivals
        }
    }
    
    return TL
}
/**
 * Mtrace Split algorithms
 */
@Suppress("unused")
class MtraceSplitAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}