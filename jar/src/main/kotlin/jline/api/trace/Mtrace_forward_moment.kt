package jline.api.trace

import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Computes the forward moments of a marked trace.
 * 
 * @param T the inter-arrival times
 * @param A the class labels
 * @param orders vector with the orders of the moments to compute
 * @param norm 0 to return F_{i,c}: M_i = sum F_{i,c}
 *             1 (default) to return F_{i,c}: M_i = sum F_{i,c} * p_c
 *             where:
 *               M_i is the class independent moment of order i
 *               F_{i,c} is the class-c forward moment of order i
 *               p_c is the probability of arrivals of class c
 * @return the forward moments as a matrix F(c,i) = F_{i,c}
 */
fun mtrace_forward_moment(T: DoubleArray, A: IntArray, orders: IntArray, norm: Int = 1): Matrix {
    val marks = A.distinct().sorted()
    val C = marks.size
    
    val M = Array(C) { DoubleArray(orders.size) }
    
    for (j in orders.indices) {
        val k = orders[j]
        for (c in 0 until C) {
            var sum = 0.0
            var count = 0
            
            for (i in 0 until T.size - 1) {
                if (A[i] == marks[c]) {
                    sum += T[i + 1].pow(k)
                    count++
                }
            }
            
            M[c][j] = if (count > 0) sum / count else 0.0
            
            if (norm != 0) {
                val countMarksInPrefix = A.sliceArray(0 until A.size - 1).count { it == marks[c] }
                if (countMarksInPrefix > 0) {
                    M[c][j] = M[c][j] * (T.size - 1).toDouble() / countMarksInPrefix
                }
            }
        }
    }
    
    return Matrix(M)
}
/**
 * Mtrace Forward Moment algorithms
 */
@Suppress("unused")
class MtraceForwardMomentAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}