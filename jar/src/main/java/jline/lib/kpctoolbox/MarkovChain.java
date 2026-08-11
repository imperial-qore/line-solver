package jline.lib.kpctoolbox;

import jline.lib.kpctoolbox.mc.CTMC;
import jline.lib.kpctoolbox.mc.DTMC;
import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Facade class for Markov Chain functions.
 */
public final class MarkovChain {

    private MarkovChain() {}

    public static Matrix ctmcSteadyState(Matrix Q) {
        double[] pi = CTMC.ctmc_solve(Q);
        Matrix result = new Matrix(1, pi.length);
        for (int i = 0; i < pi.length; i++) {
            result.set(0, i, pi[i]);
        }
        return result;
    }

    public static Matrix dtmcSteadyState(Matrix P) {
        double[] pi = DTMC.dtmc_solve(P);
        Matrix result = new Matrix(1, pi.length);
        for (int i = 0; i < pi.length; i++) {
            result.set(0, i, pi[i]);
        }
        return result;
    }

    public static Matrix ctmcTransient(Matrix Q, double t) {
        int n = Q.getNumRows();
        double[] pi0 = new double[n];
        pi0[0] = 1.0;
        Pair<double[], ?> resultPair = CTMC.ctmc_uniformization(pi0, Q, t);
        double[] pi = resultPair.getFirst();
        Matrix result = new Matrix(1, pi.length);
        for (int i = 0; i < pi.length; i++) {
            result.set(0, i, pi[i]);
        }
        return result;
    }
}
