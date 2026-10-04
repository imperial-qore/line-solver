/**
 * @file Q_Sylvest - Sylvester Equation Solver
 *
 * Solves the Sylvester equation X*kron(A,I)+BX=-I using matrix decompositions.
 *
 * Based on the Q-MAM library by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.qmam;

import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

public final class Q_Sylvest {
    private Q_Sylvest() {}

    /**
     * Solves the equation X*kron(A,I)+B*X=-I, where kron(A,I)=U*T*U' is given by
     * its complex Schur decomposition (U unitary, T upper triangular).
     *
     * <p>Mirrors MATLAB Q_Sylvest.m: B is reduced to upper Hessenberg form
     * NBAR=V'*B*V, the transformed equation Y*T+NBAR*Y=F with F=-V'*U is solved
     * column by column exploiting the triangularity of T, and the solution is
     * mapped back as X=real(V*Y*U').</p>
     *
     * @param U unitary factor of the complex Schur decomposition of kron(A,I)
     * @param T upper triangular factor of the complex Schur decomposition of kron(A,I)
     * @param B the coefficient matrix multiplying X from the left
     * @return the real solution X
     */
    public static Matrix qSylvest(ComplexMatrix U, ComplexMatrix T, Matrix B) {
        int n = B.getNumRows();

        // [V,NBAR] = hess(B), i.e. V orthogonal with V'*B*V = NBAR upper Hessenberg
        Map<String, Matrix> hessB = B.hess();
        Matrix V = hessB.get("V");
        Matrix NBAR = hessB.get("H");

        // F = -V'*U
        ComplexMatrix F = new ComplexMatrix(V.transpose().scale(-1.0)).mult(U);

        ComplexMatrix Y = ComplexMatrix.zeros(n, n);

        for (int k = 0; k < n; k++) {
            // temp = F(:,k) - Y(:,1:k-1)*T(1:k-1,k)
            ComplexMatrix temp = new ComplexMatrix(n, 1);
            for (int i = 0; i < n; i++) {
                Complex sum = F.get(i, k);
                for (int j = 0; j < k; j++) {
                    sum = sum.subtract(Y.get(i, j).multiply(T.get(j, k)));
                }
                temp.set(i, 0, sum);
            }

            // Y(:,k) = (NBAR + T(k,k)*I) \ temp
            ComplexMatrix coeff = new ComplexMatrix(NBAR)
                    .add(ComplexMatrix.eye(n).scaleComplex(T.get(k, k)));
            ComplexMatrix yk = coeff.leftMatrixDivide(temp);
            for (int i = 0; i < n; i++) {
                Y.set(i, k, yk.get(i, 0));
            }
        }

        // X = real(V*Y*U')
        ComplexMatrix X = new ComplexMatrix(V).mult(Y).mult(U.conjugateTranspose());
        return X.real;
    }

    /**
     * Computes the complex Schur decomposition of A, equivalent to MATLAB's
     * {@code [U,T] = schur(A,'complex')} with A = U*T*U'.
     *
     * @param A the matrix to decompose
     * @return a map with keys "U" (unitary) and "T" (upper triangular)
     */
    public static Map<String, ComplexMatrix> schurDecomposition(Matrix A) {
        return A.schurComplex();
    }

    /**
     * Computes kron(A, eye(m)) for a complex matrix A.
     *
     * @param A the complex matrix
     * @param m the order of the identity factor
     * @return the Kronecker product kron(A, eye(m))
     */
    public static ComplexMatrix kronEye(ComplexMatrix A, int m) {
        int rows = A.getNumRows();
        int cols = A.getNumCols();
        ComplexMatrix out = ComplexMatrix.zeros(rows * m, cols * m);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                Complex aij = A.get(i, j);
                if (aij.getReal() == 0.0 && aij.getImaginary() == 0.0) {
                    continue;
                }
                for (int d = 0; d < m; d++) {
                    out.set(i * m + d, j * m + d, aij);
                }
            }
        }
        return out;
    }
}
