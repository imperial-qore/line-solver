/**
 * @file M/M/1 Discriminatory Processor Sharing (DPS) queueing system analysis
 *
 * Numerically exact multiclass M/M/1-DPS solver via a truncated CTMC. Port of
 * the python-native qsys_mm1_dps and MATLAB qsys_mm1_dps.m.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.matrix.Matrix;

public final class Qsys_mm1_dps {
    private Qsys_mm1_dps() {}

    /**
     * Solves the M/M/1-DPS queue exactly on a truncated state space.
     *
     * The chain lives on the per-class population vector (n_1..n_K) with
     * arrival rates lambda_k and class-k completion rate
     * mu_k * n_k * w_k / sum_j n_j * w_j. The truncation level starts from the
     * geometric tail bound and doubles until the per-class mean counts are
     * stable, so conservation of the M/M/1 total holds for equal service rates
     * by construction.
     *
     * @param lambda per-class Poisson arrival rates (1,K)
     * @param mu     per-class exponential service rates (1,K)
     * @param w      per-class DPS weights (1,K), positive
     * @return per-class mean response times (1,K) via Little's law
     */
    public static Matrix qsys_mm1_dps(Matrix lambda, Matrix mu, Matrix w) {
        int K = lambda.getNumElements();
        double rho = 0.0;
        for (int k = 0; k < K; k++) {
            if (lambda.get(k) <= 0 || mu.get(k) <= 0 || w.get(k) <= 0) {
                throw new IllegalArgumentException("qsys_mm1_dps: lambda, mu, w must be positive");
            }
            rho += lambda.get(k) / mu.get(k);
        }
        if (rho >= 1.0) {
            throw new IllegalStateException("qsys_mm1_dps: system is unstable, rho >= 1");
        }

        // see _kb/03-api-layer.md for rationale
        double tolTail = 1e-8;
        int maxCutoff = (K >= 3) ? 48 : 200;
        int N = 16;
        while (N < maxCutoff && N * Math.pow(rho, N) / ((1 - rho) * (1 - rho)) > tolTail) {
            N += 8;
        }
        double[] en = solveTruncated(lambda, mu, w, K, N);
        Matrix T = new Matrix(1, K);
        for (int k = 0; k < K; k++) {
            T.set(0, k, en[k] / lambda.get(k));
        }
        return T;
    }

    private static double[] solveTruncated(Matrix lambda, Matrix mu, Matrix w, int K, int N) {
        // enumerate states with total population <= N
        List<int[]> states = new ArrayList<int[]>();
        enumerate(new int[K], 0, N, states);
        int n = states.size();
        Map<Long, Integer> index = new HashMap<Long, Integer>(2 * n);
        for (int i = 0; i < n; i++) {
            index.put(key(states.get(i), N), Integer.valueOf(i));
        }
        // sparse transition triplets and uniformization constant
        int[] src = new int[2 * K * n];
        int[] dst = new int[2 * K * n];
        double[] rate = new double[2 * K * n];
        double[] outSum = new double[n];
        int nnz = 0;
        for (int i = 0; i < n; i++) {
            int[] s = states.get(i);
            int tot = 0;
            double den = 0.0;
            for (int k = 0; k < K; k++) {
                tot += s[k];
                den += s[k] * w.get(k);
            }
            if (tot < N) {
                for (int k = 0; k < K; k++) {
                    s[k]++;
                    int j = index.get(key(s, N)).intValue();
                    s[k]--;
                    src[nnz] = i; dst[nnz] = j; rate[nnz] = lambda.get(k);
                    outSum[i] += lambda.get(k);
                    nnz++;
                }
            }
            if (tot > 0) {
                for (int k = 0; k < K; k++) {
                    if (s[k] > 0) {
                        double r = mu.get(k) * s[k] * w.get(k) / den;
                        s[k]--;
                        int j = index.get(key(s, N)).intValue();
                        s[k]++;
                        src[nnz] = i; dst[nnz] = j; rate[nnz] = r;
                        outSum[i] += r;
                        nnz++;
                    }
                }
            }
        }
        double unif = 0.0;
        for (int i = 0; i < n; i++) unif = Math.max(unif, outSum[i]);
        unif *= 1.05;
        // uniformized power iteration: pi <- pi * (I + Q/unif)
        double[] pi = new double[n];
        java.util.Arrays.fill(pi, 1.0 / n);
        double[] next = new double[n];
        for (int it = 0; it < 200000; it++) {
            for (int i = 0; i < n; i++) {
                next[i] = pi[i] * (1.0 - outSum[i] / unif);
            }
            for (int e = 0; e < nnz; e++) {
                next[dst[e]] += pi[src[e]] * rate[e] / unif;
            }
            double delta = 0.0, sum = 0.0;
            for (int i = 0; i < n; i++) {
                delta = Math.max(delta, Math.abs(next[i] - pi[i]));
                sum += next[i];
            }
            for (int i = 0; i < n; i++) pi[i] = next[i] / sum;
            if (delta < 1e-13) break;
        }
        double[] en = new double[K];
        for (int i = 0; i < n; i++) {
            int[] s = states.get(i);
            for (int k = 0; k < K; k++) {
                en[k] += pi[i] * s[k];
            }
        }
        return en;
    }

    private static void enumerate(int[] cur, int pos, int remaining, List<int[]> out) {
        if (pos == cur.length - 1) {
            for (int v = 0; v <= remaining; v++) {
                int[] s = cur.clone();
                s[pos] = v;
                out.add(s);
            }
            return;
        }
        for (int v = 0; v <= remaining; v++) {
            cur[pos] = v;
            enumerate(cur, pos + 1, remaining - v, out);
        }
        cur[pos] = 0;
    }

    private static long key(int[] s, int N) {
        long k = 0;
        for (int i = 0; i < s.length; i++) {
            k = k * (N + 1) + s[i];
        }
        return k;
    }
}
