/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.reptrans

import jline.util.matrix.Matrix

/**
 * Returns the matrix that transforms A1 to A2.
 *
 * @param A1 The smaller matrix (shape N,N)
 * @param A2 The larger matrix (shape M,M where M>=N)
 * @return The matrix B (shape N,M) satisfying A1*B = B*A2
 *
 * Note: For the existence of a (unique) solution the larger
 * matrix has to inherit the eigenvalues of the smaller one.
 */
fun similarityMatrix(A1: Matrix, A2: Matrix): Matrix {
    if (A1.numRows != A1.numCols || A2.numRows != A2.numCols) {
        throw IllegalArgumentException("SimilarityMatrix: The input matrices must be square!")
    }

    val N1 = A1.numRows
    val N2 = A2.numRows

    if (N1 > N2) {
        throw IllegalArgumentException("SimilarityMatrix: The first input matrix must be smaller than the second one!")
    }

    // Compute Schur decomposition - returns Map with "U" (orthogonal) and "T" (upper triangular)
    val schur1 = A1.schur()
    val Q1 = schur1["U"]!!
    val R1 = schur1["T"]!!

    val schur2 = A2.schur()
    val Q2 = schur2["U"]!!
    val R2 = schur2["T"]!!

    // c1 = sum(Q2', 2) = sum of columns of Q2
    val c1 = Matrix(N2, 1)
    for (i in 0 until N2) {
        var sum = 0.0
        for (j in 0 until N2) {
            sum += Q2[j, i]
        }
        c1[i, 0] = sum
    }

    // c2 = sum(Q1', 2) = sum of columns of Q1
    val c2 = Matrix(N1, 1)
    for (i in 0 until N1) {
        var sum = 0.0
        for (j in 0 until N1) {
            sum += Q1[j, i]
        }
        c2[i, 0] = sum
    }

    val I = Matrix.eye(N2)
    val X = Matrix.zeros(N1, N2)

    for (k in N1 - 1 downTo 0) {
        // M = R1(k,k)*I - R2
        val M = I.scale(R1[k, k]).sub(R2)

        // Calculate m vector
        val m = Matrix(1, N2)
        if (k < N1 - 1) {
            // m = -R1(k,k+1:end)*X(k+1:end,:)
            for (j in 0 until N2) {
                var sum = 0.0
                for (l in k + 1 until N1) {
                    sum += R1[k, l] * X[l, j]
                }
                m[0, j] = -sum
            }
        }

        // Solve the linear system [M, c1]' * x = [m, c2(k)]'
        // This is an overdetermined system - use least squares via normal equations:
        // A' * A * x = A' * b, where A = [M, c1]' and b = [m, c2(k)]'

        // Build A = [M, c1]' which is (N2+1) x N2
        val A = Matrix(N2 + 1, N2)
        for (i in 0 until N2) {
            for (j in 0 until N2) {
                A[i, j] = M[j, i]  // Transpose of M
            }
        }
        for (j in 0 until N2) {
            A[N2, j] = c1[j, 0]
        }

        // Build b = [m, c2(k)]' which is (N2+1) x 1
        val b = Matrix(N2 + 1, 1)
        for (j in 0 until N2) {
            b[j, 0] = m[0, j]
        }
        b[N2, 0] = c2[k, 0]

        // Solve using normal equations: (A' * A) * x = A' * b
        val At = A.transpose()
        val AtA = At.mult(A)
        val Atb = At.mult(b)

        // Solve AtA * solution = Atb
        val solution = Matrix(N2, 1)
        Matrix.solve(AtA, Atb, solution)

        for (j in 0 until N2) {
            X[k, j] = solution[j, 0]
        }
    }

    // B = real(Q1 * X * Q2')
    val result = Q1.mult(X).mult(Q2.transpose())

    return result
}
