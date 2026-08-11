package jline.lib.kpctoolbox.smp;

import java.util.Random;

import jline.api.mam.Map_pie;
import jline.lib.kpctoolbox.mc.DTMC;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Pair;

/**
 * Deterministic (Semi-Markov) Process functions.
 * Ported from MATLAB: matlab/lib/kpctoolbox/smp/det/
 */
public final class DET {
    private DET() {}

    public static class Triple<A,B,C> {
        public final A first; public final B second; public final C third;
        public Triple(A a, B b, C c) { this.first = a; this.second = b; this.third = c; }
        public A getFirst() { return first; }
        public B getSecond() { return second; }
        public C getThird() { return third; }
    }

    public static Matrix det_embedded(MatrixCell DET) {
        Matrix D0 = DET.get(0);
        Matrix D1 = DET.get(1);
        int n = D0.getNumRows();
        Matrix P = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            double rate = -D0.get(i, i);
            if (rate > 0) {
                for (int j = 0; j < n; j++) P.set(i, j, D1.get(i, j) / rate);
            } else {
                P.set(i, i, 1.0);
            }
        }
        return P;
    }

    public static double[] det_moment(MatrixCell DET, int[] kset) {
        Matrix D0 = DET.get(0);
        int n = D0.getNumRows();
        Matrix P = det_embedded(DET);
        Matrix piMat = new Matrix(DTMC.dtmc_solve(P));
        Matrix invD0 = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            double rate = -D0.get(i, i);
            if (rate > 0) invD0.set(i, i, 1.0 / rate);
        }
        double[] moments = new double[kset.length];
        for (int idx = 0; idx < kset.length; idx++) {
            int k = kset[idx];
            Matrix invD0k = Matrix.eye(n);
            for (int p = 0; p < k; p++) {
                Matrix temp = new Matrix(n, n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        double sum = 0.0;
                        for (int l = 0; l < n; l++) sum += invD0k.get(i, l) * invD0.get(l, j);
                        temp.set(i, j, sum);
                    }
                }
                invD0k = temp;
            }
            double moment = 0.0;
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) sum += invD0k.get(i, j);
                moment += piMat.get(i) * sum;
            }
            moments[idx] = moment * factorial(k);
        }
        return moments;
    }

    private static double factorial(int n) {
        if (n <= 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; i++) result *= (double) i;
        return result;
    }

    public static double det_scv(MatrixCell DET) {
        double[] moments = det_moment(DET, new int[]{1, 2});
        double E1 = moments[0];
        double E2 = moments[1];
        return (E2 - E1 * E1) / (E1 * E1);
    }

    public static double[] det_acf(MatrixCell DET, int[] kset) {
        Matrix D0 = DET.get(0);
        int n = D0.getNumRows();
        Matrix P = det_embedded(DET);
        Matrix piMat = new Matrix(DTMC.dtmc_solve(P));
        double[] holdTimes = new double[n];
        for (int i = 0; i < n; i++) holdTimes[i] = 1.0 / (-D0.get(i, i));
        Matrix K = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) K.set(i, j, holdTimes[i] * P.get(i, j));
        }
        double E1 = det_moment(DET, new int[]{1})[0];
        double E2 = det_moment(DET, new int[]{2})[0];
        double variance = E2 - E1 * E1;
        double[] acf = new double[kset.length];
        for (int idx = 0; idx < kset.length; idx++) {
            int k = kset[idx];
            if (k <= 0) { acf[idx] = 1.0; continue; }
            Matrix Pk = Matrix.eye(n);
            for (int p = 0; p < k - 1; p++) {
                Matrix temp = new Matrix(n, n);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        double sum = 0.0;
                        for (int l = 0; l < n; l++) sum += Pk.get(i, l) * P.get(l, j);
                        temp.set(i, j, sum);
                    }
                }
                Pk = temp;
            }
            Matrix KPk = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double sum = 0.0;
                    for (int l = 0; l < n; l++) sum += K.get(i, l) * Pk.get(l, j);
                    KPk.set(i, j, sum);
                }
            }
            Matrix KPkK = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double sum = 0.0;
                    for (int l = 0; l < n; l++) sum += KPk.get(i, l) * K.get(l, j);
                    KPkK.set(i, j, sum);
                }
            }
            double joint = 0.0;
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) sum += KPkK.get(i, j);
                joint += piMat.get(i) * sum;
            }
            acf[idx] = (variance > 0) ? (joint - E1 * E1) / variance : 0.0;
        }
        return acf;
    }

    public static Triple<double[], int[], int[]> det_sample(MatrixCell DET, int nSamples, Integer initState) {
        Matrix D0 = DET.get(0);
        Matrix D1 = DET.get(1);
        int n = D0.getNumRows();
        Random random = new Random();
        int startState;
        if (initState != null) {
            startState = initState;
        } else {
            Matrix piMat = Map_pie.map_pie(DET);
            double[] cumPi = new double[n];
            cumPi[0] = piMat.get(0);
            for (int i = 1; i < n; i++) cumPi[i] = cumPi[i - 1] + piMat.get(i);
            double r = random.nextDouble();
            int state = 0;
            for (int i = 0; i < n; i++) {
                if (r <= cumPi[i]) { state = i; break; }
            }
            startState = state;
        }
        double[] samples = new double[nSamples];
        int[] lastStates = new int[nSamples];
        int[] firstStates = new int[nSamples];
        double[][] transitionProbs = new double[n][];
        for (int i = 0; i < n; i++) {
            double totalRate = -D0.get(i, i);
            double[] probs = new double[2 * n];
            if (totalRate > 0) {
                for (int j = 0; j < n; j++) {
                    probs[j] = D0.get(i, j) / totalRate;
                    probs[n + j] = D1.get(i, j) / totalRate;
                }
                probs[i] = 0.0;
            }
            transitionProbs[i] = probs;
        }
        double[][] cumProbs = new double[n][];
        for (int i = 0; i < n; i++) {
            double[] cum = new double[2 * n];
            cum[0] = transitionProbs[i][0];
            for (int j = 1; j < 2 * n; j++) cum[j] = cum[j - 1] + transitionProbs[i][j];
            cumProbs[i] = cum;
        }
        int currentState = startState;
        for (int s = 0; s < nSamples; s++) {
            firstStates[s] = currentState;
            double interarrivalTime = 0.0;
            while (true) {
                double rate = -D0.get(currentState, currentState);
                double holdTime = (rate > 0) ? 1.0 / rate : 1.0;
                double rnd = random.nextDouble();
                int nextDest = currentState;
                for (int j = 0; j < 2 * n; j++) {
                    if (rnd <= cumProbs[currentState][j]) { nextDest = j; break; }
                }
                if (nextDest >= n) {
                    interarrivalTime += holdTime;
                    currentState = nextDest - n;
                    break;
                } else {
                    interarrivalTime += holdTime;
                    currentState = nextDest;
                }
            }
            samples[s] = interarrivalTime;
            lastStates[s] = currentState;
        }
        return new Triple<double[], int[], int[]>(samples, lastStates, firstStates);
    }

    public static Triple<double[], int[], int[]> det_sample(MatrixCell DET, int nSamples) {
        return det_sample(DET, nSamples, null);
    }

    public static MatrixCell det_sum(MatrixCell DET1, MatrixCell DET2) {
        Matrix D0_1 = DET1.get(0);
        Matrix D1_1 = DET1.get(1);
        Matrix D0_2 = DET2.get(0);
        Matrix D1_2 = DET2.get(1);
        int n1 = D0_1.getNumRows();
        int n2 = D0_2.getNumRows();
        int n = n1 * n2;
        Matrix D0 = new Matrix(n, n);
        Matrix D1 = new Matrix(n, n);
        for (int i1 = 0; i1 < n1; i1++) {
            for (int j1 = 0; j1 < n1; j1++) {
                for (int i2 = 0; i2 < n2; i2++) {
                    for (int j2 = 0; j2 < n2; j2++) {
                        int i = i1 * n2 + i2;
                        int j = j1 * n2 + j2;
                        if (i2 == j2) D0.set(i, j, D0.get(i, j) + D0_1.get(i1, j1));
                        if (i1 == j1) D0.set(i, j, D0.get(i, j) + D0_2.get(i2, j2));
                        D1.set(i, j, D1_1.get(i1, j1) * D1_2.get(i2, j2));
                    }
                }
            }
        }
        MatrixCell result = new MatrixCell(2);
        result.set(0, D0);
        result.set(1, D1);
        return result;
    }
}
