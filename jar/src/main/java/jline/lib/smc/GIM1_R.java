/**
 * @file GI/M/1-type R Matrix Computation
 *
 * Computes the R matrix for GI/M/1-type Markov chains using various methods.
 *
 * Based on the SMC Solver implementation by Benny Van Houdt.
 *
 * @since LINE 3.1.0
 */
package jline.lib.smc;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class GIM1_R {
    private GIM1_R() {}

    /**
     * Options for GIM1_R solver
     */
    public static final class GIM1ROptions {
        public final String mode;
        public final int maxNumIt;
        public final int verbose;

        public GIM1ROptions() {
            this("FI", 10000, 0);
        }

        public GIM1ROptions(String mode, int maxNumIt, int verbose) {
            this.mode = mode;
            this.maxNumIt = maxNumIt;
            this.verbose = verbose;
        }
    }

    public static Matrix gim1_R(Matrix A) {
        return gim1_R(A, new GIM1ROptions());
    }

    /**
     * Computes the R matrix for a GI/M/1-type Markov chain.
     */
    public static Matrix gim1_R(Matrix A, GIM1ROptions options) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        int numit = 0;
        double check = 1.0;
        Matrix R = Matrix.zeros(m, m);

        String mode = options.mode.toUpperCase();
        if ("FI".equals(mode)) {
            // Functional Iteration
            while (check > 1e-14 && numit < options.maxNumIt) {
                Matrix Rold = R.copy();

                R = A.extractCols(0, m).copy();
                Matrix Rpow = Rold.copy();

                for (int i = 1; i <= dega; i++) {
                    R = R.add(Rpow.mult(A.extractCols(i * m, (i + 1) * m)));
                    Rpow = Rpow.mult(Rold);
                }

                check = R.sub(Rold).infinityNorm();
                numit++;
            }
        } else if ("NI".equals(mode)) {
            // Newton Iteration
            while (check > 1e-14 && numit < options.maxNumIt) {
                Matrix Rold = R.copy();

                Matrix AR = A.extractCols(0, m).copy();
                Matrix ARprime = Matrix.zeros(m, m);
                Matrix Rpow = Matrix.eye(m);

                for (int i = 1; i <= dega; i++) {
                    Matrix Ai = A.extractCols(i * m, (i + 1) * m);
                    ARprime = ARprime.add(Rpow.mult(Ai).scale((double) i));
                    Rpow = Rpow.mult(Rold);
                    AR = AR.add(Rpow.mult(Ai));
                }

                Matrix residual = Rold.sub(AR);
                Matrix jacobian = Matrix.eye(m).sub(ARprime);

                try {
                    Matrix correction = residual.mult(jacobian.inv());
                    R = Rold.sub(correction);
                } catch (Exception e) {
                    R = AR;
                }

                check = R.sub(Rold).infinityNorm();
                numit++;
            }
        }

        if (numit == options.maxNumIt) {
            InputOutput.line_warning("GIM1_R", "Maximum number of iterations %d reached", numit);
        }

        return R;
    }
}
