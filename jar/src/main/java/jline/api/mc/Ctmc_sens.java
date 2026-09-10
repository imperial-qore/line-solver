/**
 * @file Sensitivity of a CTMC steady-state distribution to a scalar parameter
 *
 * Given a generator Q, its derivative dQ = dQ/dtheta, and the steady-state vector pi,
 * returns dpi/dtheta by solving the differentiated balance equations. Port of the
 * MATLAB/Python ctmc_sens primitive.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

public final class Ctmc_sens {
    private Ctmc_sens() {}

    /**
     * Sensitivity of the steady-state distribution of a CTMC to a scalar parameter
     * theta, computing pi from Q if it is not supplied.
     *
     * @param Q  Generator matrix (n x n)
     * @param dQ Derivative of the generator with respect to theta (n x n)
     * @return Derivative of the steady-state distribution (1 x n)
     */
    public static Matrix ctmc_sens(Matrix Q, Matrix dQ) {
        return ctmc_sens(Q, dQ, null);
    }

    /**
     * Sensitivity of the steady-state distribution of a CTMC to a scalar parameter
     * theta, given the generator Q, its derivative dQ = dQ/dtheta, and the
     * steady-state vector pi.
     *
     * <p>Differentiating the balance equations pi*Q = 0 and pi*e = 1 with respect to
     * theta gives the linear system
     *
     * <pre>  (dpi/dtheta) * Q = -pi * (dQ/dtheta),   sum_i dpi_i/dtheta = 0,</pre>
     *
     * i.e. Trivedi and Bobbio (2017), Eq. (9.81). The system has the same coefficient
     * matrix as the steady-state solve itself, so obtaining a sensitivity costs one
     * extra solve against a matrix that is already assembled. The normalization
     * replaces one row of the singular Q', exactly as in the steady-state solve.</p>
     *
     * @param Q  Generator matrix (n x n)
     * @param dQ Derivative of the generator with respect to theta (n x n)
     * @param pi Steady-state distribution (1 x n); computed from Q if null
     * @return Derivative of the steady-state distribution (1 x n)
     */
    public static Matrix ctmc_sens(Matrix Q, Matrix dQ, Matrix pi) {
        int n = Q.length();
        if (dQ.getNumRows() != n || dQ.getNumCols() != n) {
            throw new RuntimeException("dQ must have the same size as Q");
        }
        if (pi == null) {
            pi = Ctmc_solve.ctmc_solve(Q);
        }

        // Right-hand side of Eq. (9.81): b = -pi * dQ, as a column vector.
        Matrix b = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double acc = 0.0;
            for (int j = 0; j < n; j++) {
                acc += pi.get(0, j) * dQ.get(j, i);
            }
            b.set(i, 0, -acc);
        }

        // Solve dpi * Q = b subject to sum(dpi) = 0. Transpose to column form and
        // replace the last equation by the normalization, mirroring ctmc_solve.
        Matrix A = Q.transpose();
        for (int j = 0; j < n; j++) {
            A.set(n - 1, j, 1.0);
        }
        b.set(n - 1, 0, 0.0);

        Matrix x = new Matrix(n, 1);
        Matrix.solveSafe(A, b, x);
        return x.transpose();
    }
}
