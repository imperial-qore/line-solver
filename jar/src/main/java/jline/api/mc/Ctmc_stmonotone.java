package jline.api.mc;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Ctmc_stmonotone {
    private Ctmc_stmonotone() {}

    /**
     * Computes the stochastically monotone upper bound for a CTMC.
     */
    public static Matrix ctmc_stmonotone(Matrix Q) {
        double maxAbsVal = 0.0;
        for (int i = 0; i < Q.getNumRows(); i++) {
            for (int j = 0; j < Q.getNumCols(); j++) {
                maxAbsVal = Math.max(maxAbsVal, Math.abs(Q.get(i, j)));
            }
        }

        Pair<Matrix, Double> rand = Ctmc_randomization.ctmc_randomization(Q, maxAbsVal);
        Matrix P = rand.getLeft();

        Matrix P_stochastic = Dtmc_makestochastic.dtmc_makestochastic(P);

        Matrix Pub = dtmc_stmonotone(P_stochastic);

        Matrix Qub = Ctmc_makeinfgen.ctmc_makeinfgen(Pub);

        return Qub;
    }

    /**
     * Implementation of the dtmc_stmonotone algorithm.
     */
    public static Matrix dtmc_stmonotone(Matrix P) {
        int n = P.length();
        Matrix Q = new Matrix(n, n);

        Q.set(0, n - 1, P.get(0, n - 1));

        for (int i = 1; i < n; i++) {
            Q.set(i, n - 1, Math.max(Q.get(i - 1, n - 1), P.get(i, n - 1)));
        }

        for (int l = n - 2; l >= 0; l--) {
            Q.set(0, l, P.get(0, l));

            for (int i = 1; i < n; i++) {
                double sumQ_prev = 0.0;
                for (int k = l; k < n; k++) {
                    sumQ_prev += Q.get(i - 1, k);
                }

                double sumP = 0.0;
                for (int k = l; k < n; k++) {
                    sumP += P.get(i, k);
                }

                double sumQ_right = 0.0;
                for (int k = l + 1; k < n; k++) {
                    sumQ_right += Q.get(i, k);
                }

                Q.set(i, l, Math.max(sumQ_prev, sumP) - sumQ_right);
            }
        }

        return Q;
    }
}
