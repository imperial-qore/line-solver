package jline.lib.kpctoolbox.basic

import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex
import org.apache.commons.math3.linear.EigenDecomposition
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.linear.RealMatrix
import org.apache.commons.math3.util.FastMath
import java.util.Arrays

/**
 * Basic utility functions for the KPC-Toolbox.
 *
 * Ported from MATLAB: matlab/lib/kpctoolbox/basic/
 */

/**
 * Finds the position of the minimum value in a vector.
 *
 * @param v The input vector
 * @return The index (0-based) of the minimum value
 */
fun minpos(v: DoubleArray): Int {
    if (v.isEmpty()) {
        throw IllegalArgumentException("Vector cannot be empty")
    }
    var minIdx = 0
    var minVal = v[0]
    for (i in 1 until v.size) {
        if (v[i] < minVal) {
            minVal = v[i]
            minIdx = i
        }
    }
    return minIdx
}

/**
 * Finds the positions of the n smallest values in a vector.
 *
 * @param v The input vector
 * @param n The number of smallest elements to find
 * @return Array of indices (0-based) of the n smallest values, sorted by value
 */
fun minpos(v: DoubleArray, n: Int): IntArray {
    if (v.isEmpty()) {
        throw IllegalArgumentException("Vector cannot be empty")
    }
    val count = minOf(n, v.size)

    // Create array of indices paired with values
    val indexed = v.mapIndexed { idx, value -> Pair(idx, value) }

    // Sort by value and take first n
    return indexed.sortedBy { it.second }
        .take(count)
        .map { it.first }
        .toIntArray()
}

/**
 * Finds the position of the maximum value in a vector.
 *
 * @param v The input vector
 * @return The index (0-based) of the maximum value
 */
fun maxpos(v: DoubleArray): Int {
    if (v.isEmpty()) {
        throw IllegalArgumentException("Vector cannot be empty")
    }
    var maxIdx = 0
    var maxVal = v[0]
    for (i in 1 until v.size) {
        if (v[i] > maxVal) {
            maxVal = v[i]
            maxIdx = i
        }
    }
    return maxIdx
}

/**
 * Finds the positions of the n largest values in a vector.
 *
 * @param v The input vector
 * @param n The number of largest elements to find
 * @return Array of indices (0-based) of the n largest values, sorted by value descending
 */
fun maxpos(v: DoubleArray, n: Int): IntArray {
    if (v.isEmpty()) {
        throw IllegalArgumentException("Vector cannot be empty")
    }
    val count = minOf(n, v.size)

    // Create array of indices paired with values
    val indexed = v.mapIndexed { idx, value -> Pair(idx, value) }

    // Sort by value descending and take first n
    return indexed.sortedByDescending { it.second }
        .take(count)
        .map { it.first }
        .toIntArray()
}

/**
 * Generates logarithmically spaced integers in [a, b].
 *
 * @param a Lower bound
 * @param b Upper bound
 * @param points Number of points
 * @return Array of logarithmically spaced integers
 */
fun logspacei(a: Double, b: Double, points: Int): IntArray {
    if (points <= 0) {
        return IntArray(0)
    }
    if (points == 1) {
        return intArrayOf(FastMath.round(a).toInt())
    }

    val logA = FastMath.log10(a)
    val logB = FastMath.log10(b)
    val step = (logB - logA) / (points - 1)

    val result = IntArray(points)
    for (i in 0 until points) {
        var value = FastMath.round(FastMath.pow(10.0, logA + i * step)).toInt()
        // Clamp to [a, b]
        if (value < a) value = FastMath.round(a).toInt()
        if (value > b) value = FastMath.round(b).toInt()
        result[i] = value
    }
    return result
}

/**
 * Result of spectral decomposition.
 */
data class SpectralDecomposition(
    val spectrum: DoubleArray,
    val projectors: List<Matrix>,
    val eigenvectors: Matrix,
    val eigenvalueMatrix: Matrix
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false
        other as SpectralDecomposition
        if (!spectrum.contentEquals(other.spectrum)) return false
        if (projectors != other.projectors) return false
        return true
    }

    override fun hashCode(): Int {
        var result = spectrum.contentHashCode()
        result = 31 * result + projectors.hashCode()
        return result
    }
}

/**
 * Computes the spectral decomposition of a matrix.
 * Returns eigenvalues, eigenvectors, and projectors.
 *
 * @param A The input matrix
 * @return SpectralDecomposition containing spectrum, projectors, eigenvectors, and eigenvalue matrix
 */
fun spectd(A: Matrix): SpectralDecomposition {
    val n = A.numRows

    // Convert to Apache Commons Math matrix for eigendecomposition
    val realMatrix = MatrixUtils.createRealMatrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            realMatrix.setEntry(i, j, A.get(i, j))
        }
    }

    val eigen = EigenDecomposition(realMatrix)
    val realEigenvalues = eigen.realEigenvalues

    // Get eigenvector matrix V
    val vMatrix = eigen.v
    val V = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            V.set(i, j, vMatrix.getEntry(i, j))
        }
    }

    // Compute inverse of V
    val vInverse = MatrixUtils.inverse(vMatrix)
    val iV = Matrix(n, n)
    for (i in 0 until n) {
        for (j in 0 until n) {
            iV.set(i, j, vInverse.getEntry(i, j))
        }
    }

    // Build eigenvalue diagonal matrix
    val D = Matrix(n, n)
    for (i in 0 until n) {
        D.set(i, i, realEigenvalues[i])
    }

    // Compute projectors: P_k = V(:,k) * iV(k,:)
    val projectors = ArrayList<Matrix>()
    for (k in 0 until n) {
        val projector = Matrix(n, n)
        for (i in 0 until n) {
            for (j in 0 until n) {
                projector.set(i, j, V.get(i, k) * iV.get(k, j))
            }
        }
        projectors.add(projector)
    }

    return SpectralDecomposition(realEigenvalues, projectors, V, D)
}

/**
 * Creates a column vector of ones.
 *
 * @param n Size of the vector
 * @return Column vector of ones
 */
fun ones(n: Int): Matrix {
    val result = Matrix(n, 1)
    for (i in 0 until n) {
        result.set(i, 0, 1.0)
    }
    return result
}

/**
 * Creates a vector e(n) = [1, 1, ..., 1]^T used in matrix operations.
 *
 * @param n Size of the vector
 * @return Column vector of ones
 */
fun e(n: Int): Matrix {
    return ones(n)
}

/**
 * Creates an identity matrix.
 *
 * @param n Size of the matrix
 * @return Identity matrix of size n x n
 */
fun eye(n: Int): Matrix {
    return Matrix.eye(n)
}

/**
 * Creates a zero matrix.
 *
 * @param m Number of rows
 * @param n Number of columns
 * @return Zero matrix of size m x n
 */
fun zeros(m: Int, n: Int): Matrix {
    return Matrix(m, n)
}
