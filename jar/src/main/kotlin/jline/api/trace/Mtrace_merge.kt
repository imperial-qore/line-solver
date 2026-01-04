package jline.api.trace

/**
 * Merges two traces in a single marked (multiclass) trace
 * 
 * @param t1 inter-arrival times of the first trace
 * @param t2 inter-arrival times of the second trace
 * @return pair of (inter-arrival times of marked process, labels of marked process)
 */
fun mtrace_merge(t1: DoubleArray, t2: DoubleArray): Pair<DoubleArray, IntArray> {
    // Create cumulative sums starting from 0
    val cum1 = DoubleArray(t1.size + 1)
    val cum2 = DoubleArray(t2.size + 1)
    
    cum1[0] = 0.0
    for (i in t1.indices) {
        cum1[i + 1] = cum1[i] + t1[i]
    }
    
    cum2[0] = 0.0
    for (i in t2.indices) {
        cum2[i + 1] = cum2[i] + t2[i]
    }
    
    // Combine and sort all cumulative times with indices
    val allTimes = mutableListOf<Pair<Double, Int>>()
    
    // Add times from first trace (indices 1 to t1.size)
    for (i in cum1.indices) {
        allTimes.add(Pair(cum1[i], if (i == 0) 0 else i))
    }
    
    // Add times from second trace (indices t1.size+2 to t1.size+t2.size+1)
    for (i in cum2.indices) {
        allTimes.add(Pair(cum2[i], if (i == 0) 0 else t1.size + 1 + i))
    }
    
    // Sort by time
    val sorted = allTimes.sortedBy { it.first }
    
    // Compute inter-arrival times
    val T = DoubleArray(sorted.size - 1)
    for (i in 1 until sorted.size) {
        T[i - 1] = sorted[i].first - sorted[i - 1].first
    }
    
    // Determine labels
    val L = IntArray(T.size)
    for (i in 1 until sorted.size) {
        val idx = sorted[i].second
        when {
            idx in 1..t1.size -> L[i - 1] = 1
            idx >= t1.size + 2 && idx <= t1.size + 1 + t2.size -> L[i - 1] = 2
            else -> L[i - 1] = 0
        }
    }
    
    return Pair(T, L)
}
/**
 * Mtrace Merge algorithms
 */
@Suppress("unused")
class MtraceMergeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}