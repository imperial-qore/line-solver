package jline.lib.smc;

import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

public final class QBD_NI_Sylvest {
    private QBD_NI_Sylvest() {}

    /**
     * Solves the Sylvester equation AXB + CX = D using complex Schur and Hessenberg-triangular decomposition.
     *
     * <p>Mirrors MATLAB QBD_NI_Sylvest.m: the pencil (A,C) is reduced to
     * Hessenberg-triangular form W*A*V=LBAR, W*C*V=NBAR, B is reduced to complex
     * Schur form B=U*T*U', and the transformed equation is solved column by
     * column exploiting the triangularity of T, giving X=real(V*Y*U').</p>
     *
     * @param A coefficient matrix A
     * @param B coefficient matrix B
     * @param C coefficient matrix C
     * @param D right-hand side matrix D
     * @return solution matrix X satisfying AXB + CX = D
     */
    public static Matrix QBD_NI_Sylvest(Matrix A, Matrix B, Matrix C, Matrix D) {
        if (A == null || B == null || C == null || D == null) {
            throw new IllegalArgumentException("All input matrices must be non-null");
        }

        int n = A.getNumRows();
        int m = B.getNumCols();

        if (A.getNumCols() != n || C.getNumRows() != n || C.getNumCols() != n) {
            throw new IllegalArgumentException("Matrices A and C must be square and of the same size");
        }
        if (B.getNumRows() != m) {
            throw new IllegalArgumentException("Matrix B must be square");
        }
        if (D.getNumRows() != n || D.getNumCols() != m) {
            throw new IllegalArgumentException("Matrix D dimensions must match A rows and B columns");
        }

        // [LBAR,NBAR,W,V] = hess(A,C): W*A*V = LBAR, W*C*V = NBAR
        Map<String, Matrix> genHess = A.hess(C);
        Matrix lbar = genHess.get("LBAR");
        Matrix nbar = genHess.get("NBAR");
        Matrix w = genHess.get("W");
        Matrix v = genHess.get("V");

        // [U,T] = schur(B,'complex'): B = U*T*U'
        Map<String, ComplexMatrix> schur = B.schurComplex();
        ComplexMatrix uMat = schur.get("U");
        ComplexMatrix t = schur.get("T");

        // F = W*D*U
        ComplexMatrix f = new ComplexMatrix(w.mult(D)).mult(uMat);

        ComplexMatrix y = ComplexMatrix.zeros(n, m);
        ComplexMatrix tempmat = ComplexMatrix.zeros(n, m > 1 ? m - 1 : 1);

        for (int k = 0; k < m; k++) {
            ComplexMatrix temp = new ComplexMatrix(n, 1);
            if (k == 0) {
                for (int i = 0; i < n; i++) {
                    temp.set(i, 0, f.get(i, 0));
                }
            } else {
                // tempmat(:,k-1) = -LBAR*Y(:,k-1)
                for (int i = 0; i < n; i++) {
                    Complex acc = Complex.ZERO;
                    for (int p = 0; p < n; p++) {
                        acc = acc.add(y.get(p, k - 1).multiply(-lbar.get(i, p)));
                    }
                    tempmat.set(i, k - 1, acc);
                }
                // temp = F(:,k) + sum_{j<k} tempmat(:,j)*T(j,k)
                for (int i = 0; i < n; i++) {
                    Complex acc = f.get(i, k);
                    for (int j = 0; j < k; j++) {
                        acc = acc.add(tempmat.get(i, j).multiply(t.get(j, k)));
                    }
                    temp.set(i, 0, acc);
                }
            }

            // Y(:,k) = (NBAR + T(k,k)*LBAR) \ temp
            ComplexMatrix coeff = new ComplexMatrix(nbar)
                    .add(new ComplexMatrix(lbar).scaleComplex(t.get(k, k)));
            ComplexMatrix ySol = coeff.leftMatrixDivide(temp);
            for (int i = 0; i < n; i++) {
                y.set(i, k, ySol.get(i, 0));
            }
        }

        // X = real(V*Y*U')
        ComplexMatrix x = new ComplexMatrix(v).mult(y).mult(uMat.conjugateTranspose());
        return x.real;
    }
}
