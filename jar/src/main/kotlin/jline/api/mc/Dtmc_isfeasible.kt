package jline.api.mc

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Check if a matrix represents a feasible DTMC transition matrix
 *
 * @param P Transition matrix to check
 * @return Feasibility tolerance level (0 if not feasible, higher values indicate better precision)
 */
fun dtmc_isfeasible(P: Matrix): Int {
    val n = P.numRows
    val sP = DoubleArray(n)
    
    // Compute row sums
    for (i in 0 until n) {
        var sum = 0.0
        for (j in 0 until n) {
            sum += P[i, j]
        }
        sP[i] = sum
    }
    
    // Find minimum and maximum values
    var minSum = sP[0]
    var maxSum = sP[0]
    for (i in 1 until n) {
        minSum = FastMath.min(minSum, sP[i])
        maxSum = FastMath.max(maxSum, sP[i])
    }
    
    var minElement = P[0, 0]
    for (i in 0 until n) {
        for (j in 0 until n) {
            minElement = FastMath.min(minElement, P[i, j])
        }
    }
    
    var res = 0
    for (tol in 1..15) {
        val tolerance = FastMath.pow(10.0, -tol.toDouble())
        if (minSum > 1.0 - tolerance && 
            maxSum < 1.0 + tolerance && 
            minElement > -tolerance) {
            res = tol
        }
    }
    
    return res
}
/**
 * DTMC isfeasible algorithms
 */
@Suppress("unused")
class DtmcIsfeasibleAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}