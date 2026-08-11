/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.util.Pair;

public final class MAP2CorrelationBounds {
    private MAP2CorrelationBounds() {}

    /**
     * Returns the upper and lower correlation bounds for a MAP(2)
     * given the three marginal moments.
     *
     * @param moms First three marginal moments of the inter-arrival times
     * @return Pair of (lower bound, upper bound) for correlation
     */
    public static Pair<Double, Double> map2CorrelationBounds(double[] moms) {
        double m1 = moms[0];
        double m2 = moms[1];
        double m3 = moms[2];

        double h2 = m2 / (2.0 * m1 * m1) - 1;
        double h3 = m3 / (6.0 * m1 * m1 * m1) - m2 * m2 / (4.0 * m1 * m1 * m1 * m1);
        double cv2 = m2 / m1 / m1 - 1.0;

        double gub;
        if (h2 >= 0) {
            gub = h2;
        } else {
            gub = -(h2 + Math.sqrt(-h3)) * (h2 + Math.sqrt(-h3));
        }

        double glb;
        if (h2 <= 0 || h3 / h2 + h2 < 1) {
            glb = -h3 - h2 * h2;
        } else {
            double temp = h3 + h2 * h2 - h2;
            glb = h2 * (temp - Math.sqrt(temp * temp + 4.0 * h2 * h2 * h2))
                    / (temp + Math.sqrt(temp * temp + 4.0 * h2 * h2 * h2));
        }

        if (h2 >= 0) {
            return new Pair<Double, Double>(glb / cv2, gub / cv2);
        } else {
            return new Pair<Double, Double>(gub / cv2, glb / cv2);
        }
    }
}
