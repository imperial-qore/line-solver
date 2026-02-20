package jline.api.trace

import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Given a multi-class trace, computes the empirical class-dependent joint
 * moments that estimate E[ ( X^(a)_j )^i(1) (X^(a)_(j+l) )^i(2) ]
 * for all classes a.
 * 
 * @param T the inter-event times
 * @param A the class of each event
 * @param i a vector specifying the power of each variable in the joint moment
 * @return the joint moment of each class
 */
fun mtrace_joint(T: DoubleArray, A: IntArray, i: IntArray): Matrix {
    // number of classes
    val C = A.maxOrNull() ?: 0
    
    // number of events
    val N = A.size
    
    // result
    val JM = DoubleArray(C)
    
    for (a in 1..C) {
        // Count events of class a, excluding the first and last event
        val Na = A.sliceArray(1 until N - 1).count { it == a }
        
        var tmp = 0.0
        for (j in 0 until N - 2) {
            if (A[j + 1] == a) {
                tmp += T[j].pow(i[0]) * T[j + 1].pow(i[1])
            }
        }
        
        JM[a - 1] = if (Na > 0) tmp / Na else 0.0
    }
    
    return Matrix(JM)
}
/**
 * Mtrace Joint algorithms
 */
@Suppress("unused")
class MtraceJointAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}