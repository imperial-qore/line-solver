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

    /**
     * R of a GI/M/1-type chain through the G of its DUAL. Port of GIM1_R.m for
     * Dual 'A', 'R', 'B' and Algor 'FI', 'CR'.
     *
     * There is no cyclic reduction for R directly, so the chain is transposed
     * into an M/G/1-type one whose G carries the same information: the
     * Ramaswami dual for a transient chain, the Bright dual (which also
     * rescales block i by eta^(i-1), eta the caudal characteristic) for a
     * positive recurrent one. 'A' picks between them by the drift, and R comes
     * back by the inverse similarity, times eta in the Bright case.
     *
     * This is what GIM1_R_ETAQA asks for. The `gim1_R` above iterates on R
     * itself, which is a DIFFERENT algorithm, and is kept for its own callers.
     */
    public static Matrix gim1_R_dual(Matrix A, String dual, String algor) {
        int m = A.getNumRows();
        int dega = A.getNumCols() / m - 1;

        // drift > 1: positive recurrent GI/M/1; drift < 1: transient
        Matrix sumA = A.extractCols(dega * m, (dega + 1) * m).copy();
        Matrix beta = sumA.sumRows();  // column of row sums
        for (int i = dega - 1; i >= 1; i--) {
            sumA = sumA.add(A.extractCols(i * m, (i + 1) * m));
            beta = beta.add(sumA.sumRows());
        }
        sumA = sumA.add(A.extractCols(0, m));
        Matrix theta = Stat.stat(sumA);
        double drift = theta.mult(beta).get(0, 0);

        boolean ram = "R".equals(dual) || ("A".equals(dual) && drift <= 1.0);
        double eta = 1.0;
        Matrix work = new Matrix(m, m * (dega + 1));

        if (ram) {
            for (int b = 0; b <= dega; b++) {
                Matrix Ab = A.extractCols(b * m, (b + 1) * m);
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < m; j++) {
                        work.set(i, b * m + j, Ab.get(j, i) * theta.get(0, j) / theta.get(0, i));
                    }
                }
            }
        } else if ("A".equals(dual) || "B".equals(dual)) {
            eta = (drift > 1.0) ? GIM1_Caudal.gim1_caudal(A).getEta() : MG1_Decay.mg1_decay(A, false).getEta();
            Matrix sumAeta = A.extractCols(dega * m, (dega + 1) * m).scale(Math.pow(eta, (double) dega));
            for (int i = dega - 1; i >= 0; i--) {
                sumAeta = sumAeta.add(A.extractCols(i * m, (i + 1) * m).scale(Math.pow(eta, (double) i)));
            }
            Matrix shifted = sumAeta.copy();
            for (int i = 0; i < m; i++) {
                shifted.set(i, i, shifted.get(i, i) + (1.0 - eta));
            }
            theta = Stat.stat(shifted);
            for (int b = 0; b <= dega; b++) {
                Matrix Ab = A.extractCols(b * m, (b + 1) * m);
                double sc = Math.pow(eta, (double) b - 1.0);
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < m; j++) {
                        work.set(i, b * m + j,
                                sc * Ab.get(j, i) * theta.get(0, j) / theta.get(0, i));
                    }
                }
            }
        } else {
            throw new IllegalArgumentException(
                    "GIM1_R: Dual '" + dual + "' is not one of 'A', 'B', 'R'");
        }

        Matrix G;
        if ("FI".equals(algor)) {
            G = MG1_FI.mg1_fi(work);
        } else if ("CR".equals(algor)) {
            G = MG1_pi.mg1_cr(work);
        } else {
            throw new UnsupportedOperationException(
                    "GIM1_R: Algor '" + algor + "' is not ported; 'FI' (what GIM1_R_ETAQA asks "
                            + "for) and 'CR' are available");
        }

        Matrix R = new Matrix(m, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                R.set(i, j, G.get(j, i) * theta.get(0, j) / theta.get(0, i));
            }
        }
        if (!ram) {
            R = R.scale(eta);
        }
        return R;
    }
}
