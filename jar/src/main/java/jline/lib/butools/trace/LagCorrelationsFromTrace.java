/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace;

import jline.util.matrix.Matrix;

public final class LagCorrelationsFromTrace {
    private LagCorrelationsFromTrace() {}

    /**
     * Returns the lag-k autocorrelation of a trace.
     *
     * @param trace The trace data
     * @param K The number of lags to compute
     * @return The lag-k autocorrelation function of the trace up to lag K
     */
    public static double[] lagCorrelationsFromTrace(double[] trace, int K) {
        int n = trace.length;
        double sum = 0.0;
        for (double sample : trace) {
            sum += sample;
        }
        double mean = sum / n;

        // Compute variance
        double sumSq = 0.0;
        for (double sample : trace) {
            sumSq += (sample - mean) * (sample - mean);
        }
        double variance = sumSq / (n - 1);

        double[] acf = new double[K];
        for (int k = 1; k <= K; k++) {
            double s = 0.0;
            for (int i = 0; i < n - k; i++) {
                s += trace[i] * trace[i + k];
            }
            double covariance = s / (n - k) - mean * mean;
            acf[k - 1] = covariance / variance;
        }

        return acf;
    }

    public static double[] lagCorrelationsFromTrace(double[] trace) {
        return lagCorrelationsFromTrace(trace, 3);
    }

    /**
     * Returns the lag-L joint moments of a trace.
     *
     * It is computed as Nm_{i,j} = (1/(N-L)) * sum_{k=0}^{N-L-1} x_k^i * x_{k+L}^j.
     *
     * @param trace The trace data
     * @param K The joint moments are computed up to order K
     * @param L The lag at which the joint moments are computed
     * @return Matrix of shape (K+1, K+1) containing the lag-L joint moments
     */
    public static Matrix lagkJointMomentsFromTrace(double[] trace, int K, int L) {
        int n = trace.length;
        int size = K + 1;
        Matrix Nm = new Matrix(size, size);

        int count = n - L;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                double sum = 0.0;
                for (int k = 0; k < count; k++) {
                    sum += Math.pow(trace[k], i) * Math.pow(trace[k + L], j);
                }
                Nm.set(i, j, sum / count);
            }
        }

        return Nm;
    }

    public static Matrix lagkJointMomentsFromTrace(double[] trace, int K) {
        return lagkJointMomentsFromTrace(trace, K, 1);
    }
}
