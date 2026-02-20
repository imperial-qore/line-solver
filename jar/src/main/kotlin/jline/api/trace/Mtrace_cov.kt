package jline.api.trace

import jline.util.matrix.Matrix

/**
 * Computes the covariance matrix for multi-type traces.
 * 
 * @param T the trace data (inter-arrival times or values)
 * @param A the type labels for each element in T (1-indexed)
 * @return a matrix containing covariance values between different types
 */
fun mtrace_cov(T: DoubleArray, A: IntArray): Array<Array<Matrix>> {
    val C = A.maxOrNull() ?: 0
    val N = A.size
    
    val COV = Array(C) { Array(C) { Matrix.zeros(2, 2) } }
    
    for (c1 in 1..C) {
        for (c2 in 1..C) {
            val X0c1v = DoubleArray(N - 1)
            val X1c2v = DoubleArray(N - 1)
            
            for (i in 0 until N - 1) {
                val X0c1 = if (A[i] == c1) T[i] else 0.0
                val X1c2 = if (A[i + 1] == c2) T[i + 1] else 0.0
                X0c1v[i] = X0c1
                X1c2v[i] = X1c2
            }
            
            // Compute covariance matrix
            val cov = computeCovariance(X0c1v, X1c2v)
            COV[c1 - 1][c2 - 1] = cov
        }
    }
    
    return COV
}

/**
 * Computes the 2x2 covariance matrix between two vectors
 */
private fun computeCovariance(x: DoubleArray, y: DoubleArray): Matrix {
    val n = x.size
    val meanX = x.average()
    val meanY = y.average()
    
    var cov_xx = 0.0
    var cov_xy = 0.0
    var cov_yy = 0.0
    
    for (i in 0 until n) {
        val dx = x[i] - meanX
        val dy = y[i] - meanY
        cov_xx += dx * dx
        cov_xy += dx * dy
        cov_yy += dy * dy
    }
    
    cov_xx /= (n - 1)
    cov_xy /= (n - 1)
    cov_yy /= (n - 1)
    
    return Matrix(arrayOf(
        doubleArrayOf(cov_xx, cov_xy),
        doubleArrayOf(cov_xy, cov_yy)
    ))
}
/**
 * Mtrace Cov algorithms
 */
@Suppress("unused")
class MtraceCovAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}