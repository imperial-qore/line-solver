/**
 * @file Multi-trace summary statistics
 * 
 * Computes summary statistics for multiple trace analysis, providing
 * aggregate measures and descriptive statistics for multi-dimensional
 * trace data sets used in performance analysis and system characterization.
 * 
 * @since LINE 3.0
 */
package jline.api.trace

import jline.util.matrix.Matrix

/**
 * Data class representing a summary of multi-trace statistics
 */
data class MtraceSummary(
    val M: DoubleArray,           // Moments [M1, M2, M3, M4, M5]
    val ACF: DoubleArray,         // Autocorrelation function (lags 1-100)
    val F1: Matrix,               // Forward moment of order 1
    val F2: Matrix,               // Forward moment of order 2
    val B1: Matrix,               // Backward moment of order 1
    val B2: Matrix,               // Backward moment of order 2
    val C1: Matrix,               // Cross moment of order 1
    val C2: Matrix,               // Cross moment of order 2
    val Pc: Matrix,               // Class probabilities
    val Pab: Matrix               // Transition probabilities
)

/**
 * Computes comprehensive summary statistics for a multi-class trace.
 * 
 * @param T the inter-arrival times
 * @param C the class labels
 * @return summary object containing all statistics
 */
fun mtrace_summary(T: DoubleArray, C: IntArray): MtraceSummary {
    // Compute basic moments
    val M = DoubleArray(5)
    for (i in 1..5) {
        val poweredValues = DoubleArray(T.size) { Math.pow(T[it], i.toDouble()) }
        M[i - 1] = poweredValues.average()
    }
    
    // Compute autocorrelation function
    val lags = (1..100).toList().toIntArray()
    val ACF = trace_acf(T, lags)
    
    // Compute forward moments
    val F1 = mtrace_forward_moment(T, C, intArrayOf(1))
    val F2 = mtrace_forward_moment(T, C, intArrayOf(2))
    
    // Compute backward moments  
    val B1 = Matrix(mtrace_backward_moment(T, C, 1) as DoubleArray)
    val B2 = Matrix(mtrace_backward_moment(T, C, 2) as DoubleArray)
    
    // Compute cross moments
    val C1 = mtrace_cross_moment(T, C, 1)
    val C2 = mtrace_cross_moment(T, C, 2)
    
    // Compute class probabilities
    val Pc = mtrace_pc(T, C)
    
    // Compute transition probabilities
    val Pab = mtrace_sigma(T, C)
    
    return MtraceSummary(M, ACF, F1, F2, B1, B2, C1, C2, Pc, Pab)
}
/**
 * Mtrace Summary algorithms
 */
@Suppress("unused")
class MtraceSummaryAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}