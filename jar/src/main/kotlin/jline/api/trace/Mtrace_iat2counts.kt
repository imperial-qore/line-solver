package jline.api.trace

import jline.util.matrix.Matrix

/**
 * Computes the per-class counting processes of T, i.e., the counts after
 * "scale" units of time from an arrival.
 * 
 * @param T inter-arrival times
 * @param A class labels
 * @param scale time after an arrival
 * @return matrix where column k is the counting process for class k
 */
fun mtrace_iat2counts(T: DoubleArray, A: IntArray, scale: Double): Matrix {
    val n = T.size
    val CT = DoubleArray(n + 1) // cumulative sum with 0 at start
    CT[0] = 0.0
    for (i in 1..n) {
        CT[i] = CT[i - 1] + T[i - 1]
    }
    
    val K = A.distinct().sorted() // classes
    val C = Array(n - 1) { IntArray(K.size) }
    
    var previousCur = 1
    
    for (i in 0 until n - 1) {
        // Speedup loop by looking at the previous value
        var cur = if (i >= 1) maxOf(i, previousCur) else 1
        
        while (cur + 1 < n && CT[cur + 1] - CT[i + 1] <= scale) {
            cur++
        }
        
        // When the window first hits the end of the trace we return
        if (cur == n - 1) {
            for (j in K.indices) {
                C[i][j] = A.sliceArray(i + 1..cur).count { it == K[j] }
            }
            // Truncate and return
            val truncatedC = Array(i + 1) { IntArray(K.size) }
            for (idx in 0..i) {
                truncatedC[idx] = C[idx]
            }
            return Matrix(truncatedC.map { it.map { value -> value.toDouble() }.toDoubleArray() }.toTypedArray())
        }
        
        for (j in K.indices) {
            C[i][j] = A.sliceArray(i + 1..cur).count { it == K[j] }
        }
        
        previousCur = cur
    }
    
    return Matrix(C.map { it.map { value -> value.toDouble() }.toDoubleArray() }.toTypedArray())
}
/**
 * Mtrace Iat2Counts algorithms
 */
@Suppress("unused")
class MtraceIat2countsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}