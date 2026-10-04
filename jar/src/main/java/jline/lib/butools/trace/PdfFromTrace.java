/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.trace;

public final class PdfFromTrace {
    private PdfFromTrace() {}

    /**
     * Result class for PdfFromTrace containing x and y arrays.
     */
    public static final class PdfResult {
        public final double[] x;
        public final double[] y;

        public PdfResult(double[] x, double[] y) {
            this.x = x;
            this.y = y;
        }
    }

    /**
     * Returns the empirical density function of a trace.
     */
    public static PdfResult pdfFromTrace(double[] trace, double[] intBounds) {
        int numIntervals = intBounds.length - 1;
        int[] hist = new int[numIntervals];

        double traceMax = Double.NEGATIVE_INFINITY;
        for (double t : trace) {
            if (t > traceMax) traceMax = t;
        }

        for (double sample : trace) {
            for (int i = 0; i < numIntervals; i++) {
                if (sample >= intBounds[i] && sample < intBounds[i + 1]) {
                    hist[i]++;
                    break;
                }
            }
            if (trace.length > 0 && traceMax == intBounds[intBounds.length - 1]) {
                if (sample == intBounds[intBounds.length - 1]) {
                    hist[numIntervals - 1]++;
                }
            }
        }

        double n = (double) trace.length;
        double[] x = new double[numIntervals];
        double[] y = new double[numIntervals];

        for (int i = 0; i < numIntervals; i++) {
            double intLen = intBounds[i + 1] - intBounds[i];
            x[i] = (intBounds[i] + intBounds[i + 1]) / 2.0;
            y[i] = hist[i] / intLen / n;
        }

        return new PdfResult(x, y);
    }

    /**
     * Returns the empirical density function of a weighted trace.
     */
    public static PdfResult pdfFromWeightedTrace(double[] trace, double[] weights, double[] intBounds) {
        int numIntervals = intBounds.length - 1;
        double[] hist = new double[numIntervals];

        for (int j = 0; j < trace.length; j++) {
            double sample = trace[j];
            for (int i = 0; i < numIntervals; i++) {
                if (sample >= intBounds[i] && sample < intBounds[i + 1]) {
                    hist[i] += weights[j];
                    break;
                }
            }
        }

        double totalWeight = 0.0;
        for (double w : weights) totalWeight += w;
        double[] x = new double[numIntervals];
        double[] y = new double[numIntervals];

        for (int i = 0; i < numIntervals; i++) {
            double intLen = intBounds[i + 1] - intBounds[i];
            x[i] = (intBounds[i] + intBounds[i + 1]) / 2.0;
            y[i] = hist[i] / intLen / totalWeight;
        }

        return new PdfResult(x, y);
    }
}
