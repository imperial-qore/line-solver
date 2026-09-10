/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.reptrans;

import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

public final class SimilarityMatrix {
    private SimilarityMatrix() {}

    /**
     * Returns the matrix that transforms A1 to A2.
     */
    public static Matrix similarityMatrix(Matrix A1, Matrix A2) {
        if (A1.getNumRows() != A1.getNumCols() || A2.getNumRows() != A2.getNumCols()) {
            throw new IllegalArgumentException("SimilarityMatrix: The input matrices must be square!");
        }

        int N1 = A1.getNumRows();
        int N2 = A2.getNumRows();

        if (N1 > N2) {
            throw new IllegalArgumentException(
                    "SimilarityMatrix: The first input matrix must be smaller than the second one!");
        }

        Map<String, ComplexMatrix> schur1 = A1.schurComplex();
        ComplexMatrix Q1 = schur1.get("U");
        ComplexMatrix R1 = schur1.get("T");

        Map<String, ComplexMatrix> schur2 = A2.schurComplex();
        ComplexMatrix Q2 = schur2.get("U");
        ComplexMatrix R2 = schur2.get("T");

        ComplexMatrix c1 = new ComplexMatrix(N2, 1);
        for (int i = 0; i < N2; i++) {
            double sumRe = 0.0;
            double sumIm = 0.0;
            for (int j = 0; j < N2; j++) {
                Complex v = Q2.get(j, i);
                sumRe += v.getReal();
                sumIm += -v.getImaginary();
            }
            c1.set(i, 0, new Complex(sumRe, sumIm));
        }

        ComplexMatrix c2 = new ComplexMatrix(N1, 1);
        for (int i = 0; i < N1; i++) {
            double sumRe = 0.0;
            double sumIm = 0.0;
            for (int j = 0; j < N1; j++) {
                Complex v = Q1.get(j, i);
                sumRe += v.getReal();
                sumIm += -v.getImaginary();
            }
            c2.set(i, 0, new Complex(sumRe, sumIm));
        }

        ComplexMatrix I = ComplexMatrix.eye(N2);
        ComplexMatrix X = ComplexMatrix.zeros(N1, N2);

        for (int k = N1 - 1; k >= 0; k--) {
            Complex lambda_k = R1.get(k, k);
            ComplexMatrix M = I.scaleComplex(lambda_k).sub(R2);

            ComplexMatrix m = ComplexMatrix.zeros(1, N2);
            if (k < N1 - 1) {
                for (int j = 0; j < N2; j++) {
                    double sumRe = 0.0;
                    double sumIm = 0.0;
                    for (int l = k + 1; l < N1; l++) {
                        Complex r = R1.get(k, l);
                        Complex x = X.get(l, j);
                        sumRe += r.getReal() * x.getReal() - r.getImaginary() * x.getImaginary();
                        sumIm += r.getReal() * x.getImaginary() + r.getImaginary() * x.getReal();
                    }
                    m.set(0, j, new Complex(-sumRe, -sumIm));
                }
            }

            ComplexMatrix A_sys = new ComplexMatrix(N2 + 1, N2);
            for (int i = 0; i < N2; i++) {
                for (int j = 0; j < N2; j++) {
                    Complex v = M.get(i, j);
                    A_sys.set(j, i, v.conjugate());
                }
            }
            for (int i = 0; i < N2; i++) {
                Complex v = c1.get(i, 0);
                A_sys.set(N2, i, v.conjugate());
            }

            ComplexMatrix b_sys = new ComplexMatrix(N2 + 1, 1);
            for (int j = 0; j < N2; j++) {
                Complex v = m.get(0, j);
                b_sys.set(j, 0, v.conjugate());
            }
            Complex c2k = c2.get(k, 0);
            b_sys.set(N2, 0, c2k.conjugate());

            ComplexMatrix solution = A_sys.leftMatrixDivide(b_sys);

            for (int j = 0; j < N2; j++) {
                X.set(k, j, solution.get(j, 0).conjugate());
            }
        }

        ComplexMatrix resultComplex = Q1.mult(X).mult(Q2.conjugateTranspose());
        return resultComplex.real;
    }
}
