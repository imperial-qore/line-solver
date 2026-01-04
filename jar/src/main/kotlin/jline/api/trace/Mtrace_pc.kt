package jline.api.trace

import jline.util.matrix.Matrix

/**
 * Computes the probabilities of arrival for each class.
 * 
 * @param T the inter-arrival times (ignored, for orthogonality with other APIs)
 * @param C the class labels
 * @return column vector of probabilities for each class
 */
fun mtrace_pc(T: DoubleArray, C: IntArray): Matrix {
    val labels = C.distinct().sorted()
    val m = labels.size
    
    val pc = DoubleArray(m)
    
    for (i in 0 until m) {
        pc[i] = C.count { it == labels[i] }.toDouble() / C.size
    }
    
    return Matrix(pc)
}
/**
 * Mtrace Pc algorithms
 */
@Suppress("unused")
class MtracePcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}