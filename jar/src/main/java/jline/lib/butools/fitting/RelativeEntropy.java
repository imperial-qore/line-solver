/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.fitting;

public final class RelativeEntropy {
    private RelativeEntropy() {}

    /**
     * Returns the relative entropy (aka Kullback-Leibler divergence) of two vectors.
     *
     * @param p1 The first vector
     * @param p2 The second vector
     * @return The relative entropy calculated as sum(p1_i * |log(p1_i/p2_i)|)
     */
    public static double relativeEntropy(double[] p1, double[] p2) {
        if (p1.length != p2.length) {
            throw new IllegalArgumentException("Vectors must have the same length");
        }

        double re = 0.0;
        for (int i = 0; i < p1.length; i++) {
            if (p1[i] > 0.0 && p2[i] > 0.0) {
                re += p1[i] * Math.abs(Math.log(p1[i] / p2[i]));
            }
        }
        return re;
    }

    /**
     * Returns the empirical relative entropy using trace data.
     *
     * @param p1 The first vector (from empirical data)
     * @param p2 The second vector (from model)
     * @return The relative entropy
     */
    public static double empiricalRelativeEntropy(double[] p1, double[] p2) {
        return relativeEntropy(p1, p2);
    }
}
