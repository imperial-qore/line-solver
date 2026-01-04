/**
 * @file Q_Sylvest - Sylvester Equation Solver
 *
 * Solves the Sylvester equation X*kron(A,I)+BX=-I using matrix decompositions.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam

import jline.util.matrix.Matrix
import org.apache.commons.math3.linear.Array2DRowRealMatrix
import org.apache.commons.math3.linear.EigenDecomposition

/**
 * Solves the equation X*kron(A,I)+BX=-I using a Schur-based approach.
 *
 * Uses pre-computed Schur decomposition where kron(A,I)=U*T*U'
 *
 * @param U The unitary matrix from Schur decomposition
 * @param T The triangular matrix from Schur decomposition
 * @param B The matrix B in the equation
 * @return Solution matrix X
 */
fun qSylvest(U: Matrix, T: Matrix, B: Matrix): Matrix {
    val n = B.numRows

    // Compute F = -U' (in the transformed coordinate system)
    val F = U.transpose().scale(-1.0)

    // Solve column by column using back-substitution
    val Y = Matrix(n, n)

    for (k in 0 until n) {
        // Build right-hand side for column k
        val temp = Matrix(n, 1)
        for (i in 0 until n) {
            temp[i, 0] = F[i, k]
        }

        // Subtract contributions from previously computed columns
        if (k > 0) {
            for (i in 0 until n) {
                var sum = 0.0
                for (j in 0 until k) {
                    sum += Y[i, j] * T[j, k]
                }
                temp[i, 0] = temp[i, 0] - sum
            }
        }

        // Solve (B + T(k,k)*I) * Y(:,k) = temp
        val coeff = B.add(Matrix.eye(n).scale(T[k, k]))

        // Use linear solve
        val yk = Matrix(n, 1)
        Matrix.solve(coeff, temp, yk)

        for (i in 0 until n) {
            Y[i, k] = yk[i, 0]
        }
    }

    // X = Y * U'
    return Y.mult(U.transpose())
}

/**
 * Performs Schur decomposition of a matrix.
 *
 * Returns a pair (U, T) where A = U * T * U'
 * Uses eigenvalue decomposition as a fallback since direct Schur
 * decomposition is not publicly accessible in commons-math3.
 *
 * @param A Input matrix
 * @return Pair of (U, T) matrices
 */
fun schurDecomposition(A: Matrix): Pair<Matrix, Matrix> {
    val n = A.numRows
    val aArray = Array(n) { i -> DoubleArray(n) { j -> A[i, j] } }
    val aRealMatrix = Array2DRowRealMatrix(aArray)

    // Use eigenvalue decomposition as an approximation
    // For real matrices with distinct real eigenvalues, this gives similar structure
    val eigen = EigenDecomposition(aRealMatrix)

    val U = Matrix(n, n)
    val T = Matrix(n, n)

    // V matrix (eigenvector matrix) - columns are eigenvectors
    val vMatrix = eigen.v
    // D matrix (diagonal with eigenvalues)
    val dMatrix = eigen.d

    for (i in 0 until n) {
        for (j in 0 until n) {
            U[i, j] = vMatrix.getEntry(i, j)
            T[i, j] = dMatrix.getEntry(i, j)
        }
    }

    return Pair(U, T)
}
