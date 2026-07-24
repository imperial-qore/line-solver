/**
 * @file Continuous-time Markov chain Monte Carlo simulation
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.io.Ret;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

import java.util.Random;

public final class Ctmc_simulate {
    private Ctmc_simulate() {}

    public static Ret.ctmcSimulation ctmc_simulate(Matrix Q, double[] pi0, int n) {
        return ctmc_simulate(Q, pi0, n, new Random());
    }

    public static Ret.ctmcSimulation ctmc_simulate(Matrix Q, double[] pi0, int n, Random random) {
        int numStates = Q.length();
        if (pi0 == null || pi0.length == 0) {
            pi0 = new double[numStates];
            double sum = 0.0;
            for (int i = 0; i < numStates; i++) {
                pi0[i] = random.nextDouble();
                sum += pi0[i];
            }
            for (int i = 0; i < numStates; i++) {
                pi0[i] /= sum;
            }
        }
        double cumulative = 0.0;
        double r = random.nextDouble();
        int st = 0;
        for (int i = 0; i < pi0.length; i++) {
            cumulative += pi0[i];
            if (r < cumulative) {
                st = i;
                break;
            }
        }

        Matrix F = new Matrix(numStates, numStates);
        for (int i = 0; i < numStates; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < numStates; j++) {
                if (i != j) {
                    rowSum += Q.get(i, j);
                    F.set(i, j, rowSum);
                }
            }
            for (int j = 0; j < numStates; j++) {
                F.set(i, j, F.get(i, j) / rowSum);
            }
        }

        Ret.ctmcSimulation result = new Ret.ctmcSimulation();
        result.states = new int[n];
        result.sojournTimes = new double[n];

        for (int i = 0; i < n; i++) {
            result.states[i] = st;
            Exp expDist = new Exp(-Q.get(st, st));
            result.sojournTimes[i] = expDist.sample(1, random)[0];

            r = random.nextDouble();
            for (int j = 0; j < numStates; j++) {
                if (r < F.get(st, j)) {
                    st = j;
                    break;
                }
            }
        }
        return result;
    }
}
