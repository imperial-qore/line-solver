/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.fitting;

public final class SquaredDifference {
    private SquaredDifference() {}

    /**
     * Returns the squared difference between two vectors.
     *
     * @param p1 The first vector
     * @param p2 The second vector
     * @return The squared difference calculated as sum((p1_i - p2_i)^2)
     */
    public static double squaredDifference(double[] p1, double[] p2) {
        if (p1.length != p2.length) {
            throw new IllegalArgumentException("Vectors must have the same length");
        }

        double sd = 0.0;
        for (int i = 0; i < p1.length; i++) {
            double diff = p1[i] - p2[i];
            sd += diff * diff;
        }
        return sd;
    }

    /**
     * Returns the empirical squared difference using trace data.
     *
     * @param p1 The first vector (from empirical data)
     * @param p2 The second vector (from model)
     * @return The squared difference
     */
    public static double empiricalSquaredDifference(double[] p1, double[] p2) {
        return squaredDifference(p1, p2);
    }
}
