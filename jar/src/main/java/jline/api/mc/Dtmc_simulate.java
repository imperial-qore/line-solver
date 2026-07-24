/**
 * @file Discrete-time Markov chain Monte Carlo simulation
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.concurrent.ThreadLocalRandom;

import jline.util.matrix.Matrix;

public final class Dtmc_simulate {
    private Dtmc_simulate() {}

    /**
     * Simulate a discrete-time Markov chain trajectory.
     */
    public static int[] dtmc_simulate(Matrix P, Matrix pi0, int n) {
        int numStates = P.getNumRows();
        int[] sts = new int[n];

        double rnd0 = ThreadLocalRandom.current().nextDouble();
        double[] cpi0 = new double[numStates];
        double cumSum = 0.0;
        for (int i = 0; i < numStates; i++) {
            cumSum += pi0.get(i);
            cpi0[i] = cumSum;
        }

        int st = 0;
        for (int i = 0; i < numStates; i++) {
            if (rnd0 <= cpi0[i] && pi0.get(i) > 0) {
                st = i;
                break;
            }
        }

        Matrix F = Matrix.zeros(numStates, numStates);
        for (int i = 0; i < numStates; i++) {
            double cumulative = 0.0;
            for (int j = 0; j < numStates; j++) {
                cumulative += P.get(i, j);
                F.set(i, j, cumulative);
            }
        }

        for (int step = 0; step < n; step++) {
            sts[step] = st;

            if (F.get(st, numStates - 1) == 0.0 || P.get(st, st) == 1.0) {
                for (int k = step + 1; k < n; k++) {
                    sts[k] = st;
                }
                break;
            }

            double rnd = ThreadLocalRandom.current().nextDouble();
            int nextState = 0;
            for (int j = 0; j < numStates; j++) {
                if (rnd <= F.get(st, j) && P.get(st, j) > 0) {
                    nextState = j;
                    break;
                }
            }
            st = nextState;
        }

        return sts;
    }
}
