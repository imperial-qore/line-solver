/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.reptrans

import jline.util.matrix.ComplexMatrix
import jline.util.matrix.Matrix
import org.apache.commons.math3.complex.Complex

/**
 * Returns the matrix that transforms A1 to A2.
 *
 * @param A1 The smaller matrix (shape N,N)
 * @param A2 The larger matrix (shape M,M where M>=N)
 * @return The matrix B (shape N,M) satisfying A1*B = B*A2
 *
 * Note: For the existence of a (unique) solution the larger
 * matrix has to inherit the eigenvalues of the smaller one.
 *
 * Uses the complex Schur decomposition (equivalent to MATLAB's schur(A,'complex'))
 * to ensure a truly upper-triangular T matrix, which is required for the
 * row-by-row backsolve algorithm.
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

    // Compute complex Schur decomposition: A = Q * R * Q^H
    val schur1 = A1.schurComplex()
    val Q1 = schur1["U"]!!
    val R1 = schur1["T"]!!

    val schur2 = A2.schurComplex()
    val Q2 = schur2["U"]!!
    val R2 = schur2["T"]!!

    // c1 = sum(Q2', 2) = row sums of Q2^H
    // Row i of Q2^H = conj(col i of Q2), so c1(i) = sum_j conj(Q2(j,i))
    val c1 = ComplexMatrix(N2, 1)
    for (i in 0 until N2) {
        var sumRe = 0.0
        var sumIm = 0.0
        for (j in 0 until N2) {
            val v = Q2.get(j, i)
            sumRe += v.real
            sumIm += -v.imaginary  // conjugate
        }
        c1.set(i, 0, Complex(sumRe, sumIm))
    }

    // c2 = sum(Q1', 2)
    val c2 = ComplexMatrix(N1, 1)
    for (i in 0 until N1) {
        var sumRe = 0.0
        var sumIm = 0.0
        for (j in 0 until N1) {
            val v = Q1.get(j, i)
            sumRe += v.real
            sumIm += -v.imaginary  // conjugate
        }
        c2.set(i, 0, Complex(sumRe, sumIm))
    }

    val I = ComplexMatrix.eye(N2)
    val X = ComplexMatrix.zeros(N1, N2)

    for (k in N1 - 1 downTo 0) {
        // M = R1(k,k)*I - R2
        val lambda_k = R1.get(k, k)  // Complex eigenvalue
        val M = I.scaleComplex(lambda_k).sub(R2)

        // Calculate m vector (1 x N2)
        val m = ComplexMatrix.zeros(1, N2)
        if (k < N1 - 1) {
            // m = -R1(k,k+1:end)*X(k+1:end,:)
            for (j in 0 until N2) {
                var sumRe = 0.0
                var sumIm = 0.0
                for (l in k + 1 until N1) {
                    val r = R1.get(k, l)
                    val x = X.get(l, j)
                    sumRe += r.real * x.real - r.imaginary * x.imaginary
                    sumIm += r.real * x.imaginary + r.imaginary * x.real
                }
                m.set(0, j, Complex(-sumRe, -sumIm))
            }
        }

        // Solve: MATLAB does X(k,:) = linsolve([M,c1]',[m,c2(k)]')'
        // Build A_sys = [M, c1]^H  ((N2+1) x N2)
        val A_sys = ComplexMatrix(N2 + 1, N2)
        for (i in 0 until N2) {
            for (j in 0 until N2) {
                val v = M.get(i, j)
                A_sys.set(j, i, v.conjugate())
            }
        }
        for (i in 0 until N2) {
            val v = c1.get(i, 0)
            A_sys.set(N2, i, v.conjugate())
        }

        // Build b_sys = [m, c2(k)]^H  ((N2+1) x 1)
        val b_sys = ComplexMatrix(N2 + 1, 1)
        for (j in 0 until N2) {
            val v = m.get(0, j)
            b_sys.set(j, 0, v.conjugate())
        }
        val c2k = c2.get(k, 0)
        b_sys.set(N2, 0, c2k.conjugate())

        // Solve A_sys * x = b_sys (overdetermined, least squares)
        val solution = A_sys.leftMatrixDivide(b_sys)

        for (j in 0 until N2) {
            // MATLAB: X(k,:) = linsolve(...)' — the trailing ' conjugate-transposes the solution
            X.set(k, j, solution.get(j, 0).conjugate())
        }
    }

    // B = real(Q1 * X * Q2^H)
    val resultComplex = Q1.mult(X).mult(Q2.conjugateTranspose())
    return resultComplex.real
}
