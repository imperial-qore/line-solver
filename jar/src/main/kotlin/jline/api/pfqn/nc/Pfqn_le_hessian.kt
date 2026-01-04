/**
 * @file Hessian matrix computation for Laguerre expansion method
 * 
 * Computes the Hessian matrix used in the Laguerre expansion method for second-order
 * corrections to the asymptotic approximation. Provides curvature information at the
 * saddle point for enhanced accuracy in normalizing constant computation.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.util.matrix.Matrix

/**
 * Auxiliary function to compute the Hessian used in the logistic expansion method
 *
 * @param L  - demands at all stations
 * @param N  - number of jobs for each class
 * @param u0 - term appearing in the integrand
 * @return normalizing constant and its logarithm
 */
fun pfqn_le_hessian(L: Matrix, N: Matrix, u0: Matrix): Matrix {
    val M = L.numRows
    val R = L.numCols
    val Ntot = N.elementSum()
    val hu = Matrix(M - 1, M - 1)
    hu.fill(0.0)

    for (i in 0..<M - 1) {
        for (j in 0..<M - 1) {
            if (i != j) {
                hu[i, j] = -(Ntot + M) * u0[i] * u0[j]
                for (r in 0..<R) {
                    val L_col_r = Matrix(L.numRows, 1)
                    Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
                    hu[i, j] =
                        (hu[i, j] + N[r] * L[i, r] * L[j, r] * u0[i] * u0[j] / (u0.mult(L_col_r)[0] * u0.mult(L_col_r)[0]))
                }
            } else {
                hu[i, j] = (Ntot + M) * u0[i] * (Matrix.allbut(u0, i).elementSum())
                for (r in 0..<R) {
                    val L_col_r = Matrix(L.numRows, 1)
                    Matrix.extract(L, 0, L.numRows, r, r + 1, L_col_r, 0, 0)
                    val tmp_L = Matrix.allbut(L_col_r, i).transpose()

                    hu[i, j] = (hu[i, j] - (N[r] * L[i, r] * u0[i] * (Matrix.allbut(u0, i).mult(tmp_L)[0])) / (u0.mult(
                        L_col_r)[0] * u0.mult(L_col_r)[0]))
                }
            }
        }
    }
    return hu
}
/**
 * PFQN le hessian algorithms
 */
@Suppress("unused")
class PfqnLeHessianAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}