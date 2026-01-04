/**
 * @file GI/M/1-type Caudal Characteristic
 *
 * Computes the spectral radius of R for GI/M/1-type Markov chains.
 * Uses binary search to find the dominant eigenvalue.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc

import jline.util.matrix.Matrix
import org.apache.commons.math3.linear.EigenDecomposition
import org.apache.commons.math3.linear.Array2DRowRealMatrix

/**
 * Result of GIM1 Caudal computation
 */
data class GIM1CaudalResult(
    val eta: Double,
    val v: Matrix?
)

/**
 * Computes the dominant eigenvalue of the R matrix for GI/M/1-type chains.
 *
 * Finds the unique solution eta of PF(A(z)) = z on (0,1), where PF denotes
 * the Perron-Frobenius eigenvalue and A(z) = A0 + A1*z + A2*z^2 + ...
 *
 * @param A Block matrix [A0 A1 A2 ... A_max]
 * @param dual When true, compute for Ramaswami dual
 * @param computeEigenvector When true, also compute the right eigenvector
 * @return GIM1CaudalResult containing eta and optionally v
 */
fun gim1_caudal(
    A: Matrix,
    dual: Boolean = false,
    computeEigenvector: Boolean = false
): GIM1CaudalResult {
    val m = A.numRows
    val dega = A.numCols / m - 1

    // Binary search for eta
    var etaMin = 0.0
    var etaMax = 1.0
    var eta = 0.5

    while (etaMax - etaMin > 1e-15) {
        // Compute temp = A(eta) = sum(A_i * eta^i)
        var temp = A.extractCols(dega * m, (dega + 1) * m).copy()
        for (i in dega - 1 downTo 0) {
            temp = temp.scale(eta).add(A.extractCols(i * m, (i + 1) * m))
        }

        // Compute spectral radius (max eigenvalue)
        val newEta = computeSpectralRadius(temp)

        if (newEta > eta) {
            etaMin = eta
        } else {
            etaMax = eta
        }
        eta = (etaMin + etaMax) / 2
    }

    // Compute eigenvector if requested
    var v: Matrix? = null
    if (computeEigenvector) {
        var temp = A.extractCols(dega * m, (dega + 1) * m).copy()
        for (i in dega - 1 downTo 0) {
            temp = temp.scale(eta).add(A.extractCols(i * m, (i + 1) * m))
        }
        v = computeDominantEigenvector(temp)
    }

    return GIM1CaudalResult(eta, v)
}

/**
 * Compute the spectral radius (largest eigenvalue magnitude) of a matrix
 */
private fun computeSpectralRadius(A: Matrix): Double {
    try {
        val m = A.numRows
        val array = Array(m) { i -> DoubleArray(m) { j -> A[i, j] } }
        val realMatrix = Array2DRowRealMatrix(array)
        val decomposition = EigenDecomposition(realMatrix)

        var maxEigenvalue = 0.0
        for (i in 0 until m) {
            val real = decomposition.realEigenvalues[i]
            val imag = decomposition.imagEigenvalues[i]
            val magnitude = kotlin.math.sqrt(real * real + imag * imag)
            if (real > maxEigenvalue) {
                maxEigenvalue = real
            }
        }
        return maxEigenvalue
    } catch (e: Exception) {
        // Fallback: use power iteration
        return powerIterationSpectralRadius(A)
    }
}

/**
 * Power iteration for spectral radius computation
 */
private fun powerIterationSpectralRadius(A: Matrix, maxIter: Int = 100): Double {
    val m = A.numRows
    var v = Matrix.ones(m, 1).scale(1.0 / m)

    for (iter in 0 until maxIter) {
        val Av = A.mult(v)
        val norm = Av.infinityNorm()
        if (norm > 0) {
            v = Av.scale(1.0 / norm)
        }
    }

    val Av = A.mult(v)
    return v.transpose().mult(Av)[0, 0] / v.transpose().mult(v)[0, 0]
}

/**
 * Compute the dominant eigenvector of a matrix
 */
private fun computeDominantEigenvector(A: Matrix): Matrix {
    try {
        val m = A.numRows
        val array = Array(m) { i -> DoubleArray(m) { j -> A[i, j] } }
        val realMatrix = Array2DRowRealMatrix(array)
        val decomposition = EigenDecomposition(realMatrix)

        // Find index of dominant eigenvalue
        var maxIdx = 0
        var maxReal = Double.NEGATIVE_INFINITY
        for (i in 0 until m) {
            val real = decomposition.realEigenvalues[i]
            if (real > maxReal) {
                maxReal = real
                maxIdx = i
            }
        }

        val eigenvector = decomposition.getEigenvector(maxIdx)
        val result = Matrix(m, 1)
        for (i in 0 until m) {
            result[i, 0] = eigenvector.getEntry(i)
        }
        return result
    } catch (e: Exception) {
        // Fallback: power iteration
        val m = A.numRows
        var v = Matrix.ones(m, 1).scale(1.0 / m)
        for (iter in 0 until 100) {
            val Av = A.mult(v)
            val norm = Av.infinityNorm()
            if (norm > 0) {
                v = Av.scale(1.0 / norm)
            }
        }
        return v
    }
}
