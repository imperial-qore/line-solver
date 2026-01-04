/**
 * @file Hessian matrix computation for Laguerre expansion method with think times
 * 
 * Computes the Hessian matrix used in the Laguerre expansion method for closed queueing
 * networks with think times (delays). Provides second-order derivatives for enhanced
 * accuracy in normalizing constant computation with delay stations.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.util.matrix.Matrix

/**
 * Auxiliary function to compute the Hessian used in the logistic expansion method in models with delays
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think times
 * @param u - term appearing in the integrand
 * @param v - term appearing in the integrand
 * @return Hessian matrix
 */
internal fun pfqn_le_hessianZ(L: Matrix, N: Matrix, Z: Matrix, u: Matrix, v: Double): Matrix {
    val K = L.numRows
    val R = L.numCols
    val Ntot = N.elementSum()
    var A = Matrix(K, K)
    A.fill(0.0)
    val csi = Matrix(1, R)
    for (r in 0..<R) {
        csi[r] = N[r] / (Z[r] + v * u.mult(Matrix.extractColumn(L, r, null))[0])
    }
    val Lhat = Matrix(K, R)
    Lhat.fill(0.0)
    for (k in 0..<K) {
        for (r in 0..<R) {
            Lhat[k, r] = Z[r] + v * L[k, r]
        }
    }
    val eta = Ntot + K
    for (i in 0..<K) {
        for (j in 0..<K) {
            if (i != j) {
                A[i, j] = -eta * u[i] * u[j]
                for (r in 0..<R) {
                    A[i, j] = (A[i, j] + csi[r] * csi[r] * Lhat[i, r] * Lhat[j, r] * u[i] * u[j] / N[r])
                }
            }
        }
    }
    for (i in 0..<K) {
        A[i, i] = -Matrix.allbut(Matrix.extractRows(A, i, i + 1, null), i).elementSum()
    }
    val tmp_A = Matrix(K, K)
    tmp_A.fill(0.0)
    Matrix.extract(A, 0, K - 1, 0, K - 1, tmp_A, 0, 0)
    A = tmp_A
    A[K - 1, K - 1] = 1.0

    for (r in 0..<R) {
        val L_col_r = Matrix(L.numRows, 1)
        Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
        A[K - 1, K - 1] = A[K - 1, K - 1] - (csi[r] * csi[r] / N[r]) * Z[r] * u.mult(L_col_r)[0]
    }

    A[K - 1, K - 1] = v * A[K - 1, K - 1]

    for (i in 0..<K - 1) {
        for (r in 0..<R) {
            val L_col_r = Matrix(L.numRows, 1)
            Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
            A[i, K - 1] =
                A[i, K - 1] + v * u[i] * (((csi[r] * csi[r] / N[r]) * Lhat[i, r] * (u.mult(L_col_r)[0])) - csi[r] * L[i, r])
        }
        A[K - 1, i] = A[i, K - 1]
    }

    return A
}
/**
 * PFQN le hessianZ algorithms
 */
@Suppress("unused")
class PfqnLeHessianzAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}