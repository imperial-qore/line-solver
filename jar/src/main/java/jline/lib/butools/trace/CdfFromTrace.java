/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace;

import java.util.Arrays;

public final class CdfFromTrace {
    private CdfFromTrace() {}

    /**
     * Result class for CdfFromTrace containing x and y arrays.
     */
    public static final class CdfResult {
        public final double[] x;
        public final double[] y;

        public CdfResult(double[] x, double[] y) {
            this.x = x;
            this.y = y;
        }

        public double[] getX() { return x; }
        public double[] getY() { return y; }
    }

    /**
     * Returns the empirical distribution function of the trace.
     *
     * @param trace The trace data
     * @return CdfResult containing x (points) and y (values) of the empirical CDF
     */
    public static CdfResult cdfFromTrace(double[] trace) {
        double[] sorted = trace.clone();
        Arrays.sort(sorted);
        int n = trace.length;

        double[] x = sorted;
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = (double) i / (n - 1);
        }

        return new CdfResult(x, y);
    }

    /**
     * Returns the empirical distribution function of a weighted trace.
     *
     * @param trace The trace data
     * @param weights The weights for each sample
     * @return CdfResult containing x (points) and y (values) of the empirical CDF
     */
    public static CdfResult cdfFromWeightedTrace(double[] trace, double[] weights) {
        // Sort trace and weights together by trace values
        Integer[] indices = new Integer[trace.length];
        for (int i = 0; i < trace.length; i++) indices[i] = i;
        Arrays.sort(indices, (a, b) -> Double.compare(trace[a], trace[b]));

        double[] x = new double[trace.length];
        double[] sortedWeights = new double[trace.length];
        for (int i = 0; i < trace.length; i++) {
            x[i] = trace[indices[i]];
            sortedWeights[i] = weights[indices[i]];
        }

        double totalWeight = 0.0;
        for (double w : sortedWeights) totalWeight += w;
        double cumWeight = 0.0;
        double[] y = new double[trace.length];
        for (int i = 0; i < trace.length; i++) {
            cumWeight += sortedWeights[i];
            y[i] = cumWeight / totalWeight;
        }

        return new CdfResult(x, y);
    }
}
