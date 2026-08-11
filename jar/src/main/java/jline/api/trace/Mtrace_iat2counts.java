package jline.api.trace;

import jline.util.matrix.Matrix;

import java.util.TreeSet;

public final class Mtrace_iat2counts {
    private Mtrace_iat2counts() {}

    /**
     * Computes the per-class counting processes of T, i.e., the counts after
     * "scale" units of time from an arrival.
     *
     * @param T inter-arrival times
     * @param A class labels
     * @param scale time after an arrival
     * @return matrix where column k is the counting process for class k
     */
    public static Matrix mtrace_iat2counts(double[] T, int[] A, double scale) {
        int n = T.length;
        double[] CT = new double[n + 1];
        CT[0] = 0.0;
        for (int i = 1; i <= n; i++) {
            CT[i] = CT[i - 1] + T[i - 1];
        }

        TreeSet<Integer> uniq = new TreeSet<Integer>();
        for (int v : A) uniq.add(v);
        Integer[] K = uniq.toArray(new Integer[0]);
        int[][] C = new int[n - 1][K.length];

        int previousCur = 1;

        for (int i = 0; i < n - 1; i++) {
            // Speedup loop by looking at the previous value
            int cur = i >= 1 ? Math.max(i, previousCur) : 1;

            while (cur + 1 < n && CT[cur + 1] - CT[i + 1] <= scale) {
                cur++;
            }

            // When the window first hits the end of the trace we return
            if (cur == n - 1) {
                for (int j = 0; j < K.length; j++) {
                    int count = 0;
                    for (int idx = i + 1; idx <= cur; idx++) {
                        if (A[idx] == K[j].intValue()) count++;
                    }
                    C[i][j] = count;
                }
                // Truncate and return
                double[][] truncatedC = new double[i + 1][K.length];
                for (int idx = 0; idx <= i; idx++) {
                    for (int j = 0; j < K.length; j++) {
                        truncatedC[idx][j] = C[idx][j];
                    }
                }
                return new Matrix(truncatedC);
            }

            for (int j = 0; j < K.length; j++) {
                int count = 0;
                for (int idx = i + 1; idx <= cur; idx++) {
                    if (A[idx] == K[j].intValue()) count++;
                }
                C[i][j] = count;
            }

            previousCur = cur;
        }

        double[][] Cdouble = new double[n - 1][K.length];
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < K.length; j++) {
                Cdouble[i][j] = C[i][j];
            }
        }
        return new Matrix(Cdouble);
    }
}
