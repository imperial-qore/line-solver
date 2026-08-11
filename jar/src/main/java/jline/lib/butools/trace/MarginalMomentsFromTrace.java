/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace;

public final class MarginalMomentsFromTrace {
    private MarginalMomentsFromTrace() {}

    /**
     * Returns the marginal moments of a trace.
     *
     * @param trace The trace data
     * @param K The number of moments to compute (default 5)
     * @return The (raw) moments of the trace
     */
    public static double[] marginalMomentsFromTrace(double[] trace, int K) {
        double[] moms = new double[K];
        double n = (double) trace.length;

        for (int i = 0; i < K; i++) {
            double sum = 0.0;
            for (double sample : trace) {
                sum += Math.pow(sample, i + 1);
            }
            moms[i] = sum / n;
        }

        return moms;
    }

    public static double[] marginalMomentsFromTrace(double[] trace) {
        return marginalMomentsFromTrace(trace, 5);
    }

    /**
     * Returns the marginal moments of a weighted trace.
     *
     * @param trace The trace data
     * @param weights The weights for each sample
     * @param K The number of moments to compute (default 5)
     * @return The (raw) moments of the weighted trace
     */
    public static double[] marginalMomentsFromWeightedTrace(double[] trace, double[] weights, int K) {
        double[] moms = new double[K];
        double totalWeight = 0.0;
        for (double w : weights) totalWeight += w;

        for (int i = 0; i < K; i++) {
            double sum = 0.0;
            for (int j = 0; j < trace.length; j++) {
                sum += weights[j] * Math.pow(trace[j], i + 1);
            }
            moms[i] = sum / totalWeight;
        }

        return moms;
    }

    public static double[] marginalMomentsFromWeightedTrace(double[] trace, double[] weights) {
        return marginalMomentsFromWeightedTrace(trace, weights, 5);
    }
}
