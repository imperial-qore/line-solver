/**
 * @file Joint density inner product of two MAPs via recursive Sylvester equations
 *
 * Computes the inner product of lag-L joint densities of two MAPs using a recursive
 * scheme based on continuous Sylvester equations.
 *
 * Reference:
 *   G. Horvath, "Measuring the distance between MAPs and some
 *   applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
 *   https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;

public final class Map_exp_mul_int {
    private Map_exp_mul_int() {}

    /**
     * Computes the inner product of lag-L joint densities of two MAPs.
     *
     * Uses a recursive scheme based on solving continuous Sylvester equations:
     *   Z_1 = sylv(B0', A0, alB' * alA)
     *   Z_i = sylv(B0', A0, B1' * Z_{i-1} * A1)   for i = 2, ..., L
     * The result is exitB' * Z_L * exitA, where exitA = -A0 * e, exitB = -B0 * e.
     *
     * @param A0  hidden transition matrix of the first MAP
     * @param A1  visible transition matrix of the first MAP
     * @param B0  hidden transition matrix of the second MAP
     * @param B1  visible transition matrix of the second MAP
     * @param L   number of inter-arrival intervals (lag parameter)
     * @param alA stationary vector at arrivals of the first MAP
     * @param alB stationary vector at arrivals of the second MAP
     * @return the inner product value (scalar)
     */
    public static double map_exp_mul_int(Matrix A0, Matrix A1, Matrix B0, Matrix B1, int L, Matrix alA, Matrix alB) {
        Matrix P = B0.transpose();
        Matrix Z = contSylv(P, A0, alB.transpose().mult(alA));
        for (int i = 1; i < L; i++) {
            Z = contSylv(P, A0, B1.transpose().mult(Z).mult(A1));
        }
        Matrix aA = A0.copy();
        aA.scaleEq(-1.0);
        Matrix exitA = aA.sumRows();
        Matrix aB = B0.copy();
        aB.scaleEq(-1.0);
        Matrix exitB = aB.sumRows();
        return exitB.transpose().mult(Z).mult(exitA).toDouble();
    }

    /**
     * Solves the continuous Sylvester equation P*X + X*Q + C = 0 via the
     * Kronecker formulation (I (x) P + Q' (x) I) vec(X) = -vec(C). This
     * replaces the Schur-based Matrix.sylv, which returns an incorrect
     * solution on the forward (A' != B) path used by the MAP cross-distance,
     * yielding negative "squared" distances. Matches MATLAB/native-Python.
     */
    private static Matrix contSylv(Matrix P, Matrix Q, Matrix C) {
        int n = P.getNumRows();
        int m = Q.getNumRows();
        Matrix lhs = Matrix.eye(m).kron(P).add(1.0, Q.transpose().kron(Matrix.eye(n)));
        Matrix vecC = new Matrix(n * m, 1, n * m);
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
                vecC.set(j * n + i, 0, -C.get(i, j));
            }
        }
        Matrix vecX = Matrix.robustLeftDivide(lhs, vecC);
        Matrix X = new Matrix(n, m, n * m);
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
                X.set(i, j, vecX.get(j * n + i, 0));
            }
        }
        return X;
    }
}
