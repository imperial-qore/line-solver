/**
 * @file Multi-class trace moment computation
 * 
 * Computes empirical class-dependent statistical moments for multi-class trace data. 
 * Enables higher-order statistical characterization of per-class arrival processes 
 * for advanced queueing model fitting and analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.trace

import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Computes the empirical class-dependent moments of a multi-class trace.
 * 
 * @param T vector of inter-arrival times
 * @param A vector of class labels
 * @param orders vector with the orders of the moments to compute
 * @param after 0 to compute moments of Horvath variables
 *              1 to compute moments of Bucholz variables
 * @param norm 0 to not normalize, M_i = sum M_{i,c}
 *             1 to normalize, M_i = sum M_{i,c} * g_c
 *             where M_i is the class independent moment of order i
 *                   M_{i,c} is the class c moment of order i
 *                   g_c is the fraction of arrivals of class c
 * @return matrix of moments
 */
fun mtrace_moment(T: DoubleArray, A: IntArray, orders: IntArray, after: Int = 0, norm: Int = 0): Matrix {
    val marks = A.distinct().sorted()
    val C = marks.size
    
    val M = Array(C) { DoubleArray(orders.size) }
    
    for (j in orders.indices) {
        val k = orders[j]
        for (c in 0 until C) {
            val mark = marks[c]
            
            if (after == 1) {
                // Bucholz variables: moments of T(2:end) where A(1:end-1) == mark
                var sum = 0.0
                var count = 0
                for (i in 0 until T.size - 1) {
                    if (A[i] == mark) {
                        sum += T[i + 1].pow(k)
                        count++
                    }
                }
                M[c][j] = if (count > 0) sum / count else 0.0
                
                if (norm == 1) {
                    val totalCount = A.sliceArray(0 until A.size - 1).size
                    val classCount = A.sliceArray(0 until A.size - 1).count { it == mark }
                    if (classCount > 0) {
                        M[c][j] = M[c][j] * totalCount.toDouble() / classCount
                    }
                }
            } else {
                // Horvath variables: moments of T where A == mark
                var sum = 0.0
                var count = 0
                for (i in T.indices) {
                    if (A[i] == mark) {
                        sum += T[i].pow(k)
                        count++
                    }
                }
                M[c][j] = if (count > 0) sum / count else 0.0
                
                if (norm == 1) {
                    val totalCount = A.size
                    val classCount = A.count { it == mark }
                    if (classCount > 0) {
                        M[c][j] = M[c][j] * totalCount.toDouble() / classCount
                    }
                }
            }
        }
    }
    
    return Matrix(M)
}
/**
 * Mtrace Moment algorithms
 */
@Suppress("unused")
class MtraceMomentAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}