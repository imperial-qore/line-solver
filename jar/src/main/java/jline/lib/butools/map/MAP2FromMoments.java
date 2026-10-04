/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.ph.PH2From3Moments;
import jline.lib.butools.ph.PH2From3Moments.PH2Representation;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class MAP2FromMoments {
    private MAP2FromMoments() {}

    /**
     * Returns a MAP(2) which has the same 3 marginal moments
     * and lag-1 autocorrelation as given.
     *
     * @param moms First three marginal moments of the inter-arrival times
     * @param corr1 The lag-1 autocorrelation of the inter-arrival times
     * @return Pair of (D0, D1) matrices of the MAP(2)
     */
    public static Pair<Matrix, Matrix> map2FromMoments(double[] moms, double corr1) {
        double m1 = moms[0];
        double m2 = moms[1];

        double prec = 1e-12;

        // If we have an exponential distribution, we do not allow correlation
        if (Math.abs(m2 - 2.0 * m1 * m1) < prec && Math.abs(corr1) > prec) {
            throw new IllegalArgumentException("We do not allow correlation in case of exponentially distributed marginal");
        }

        // Perform PH fitting
        PH2Representation phResult = PH2From3Moments.ph2From3Moments(moms);
        Matrix tau = phResult.alpha;
        Matrix A = phResult.A;

        double l1 = -A.get(0, 0);
        double l2 = -A.get(1, 1);
        double p = tau.get(0, 0);
        double alpha = l1 / l2;

        // Check the feasibility of the correlation parameter
        Pair<Double, Double> bounds = MAP2CorrelationBounds.map2CorrelationBounds(moms);
        double corrl = bounds.getLeft();
        double corru = bounds.getRight();
        if (corr1 < corrl) {
            throw new IllegalArgumentException("The correlation parameter is too small!");
        }
        if (corr1 > corru) {
            throw new IllegalArgumentException("The correlation parameter is too large!");
        }

        double gamma = corr1 * (m2 - m1 * m1) / (m2 / 2.0 - m1 * m1);

        // Perform MAP fitting
        Matrix D0;
        Matrix D1;

        if (gamma > 0) {
            double discriminant = (1.0 + alpha * gamma - p * (1.0 - gamma)) * (1.0 + alpha * gamma - p * (1.0 - gamma)) - 4.0 * alpha * gamma;
            double a = (1.0 + alpha * gamma - p * (1.0 - gamma) - Math.sqrt(discriminant)) / (2.0 * alpha);
            double b = (1.0 + alpha * gamma - p * (1.0 - gamma) + Math.sqrt(discriminant)) / 2.0;

            D0 = new Matrix(2, 2);
            D0.set(0, 0, -l1);
            D0.set(0, 1, (1.0 - a) * l1);
            D0.set(1, 0, 0.0);
            D0.set(1, 1, -l2);

            D1 = new Matrix(2, 2);
            D1.set(0, 0, a * l1);
            D1.set(0, 1, 0.0);
            D1.set(1, 0, (1.0 - b) * l2);
            D1.set(1, 1, b * l2);
        } else if (gamma < 0) {
            double a = gamma / (alpha * gamma - p * (1.0 - gamma));
            double b = p * (1.0 - gamma) - alpha * gamma;

            D0 = new Matrix(2, 2);
            D0.set(0, 0, -l1);
            D0.set(0, 1, (1.0 - a) * l1);
            D0.set(1, 0, 0.0);
            D0.set(1, 1, -l2);

            D1 = new Matrix(2, 2);
            D1.set(0, 0, 0.0);
            D1.set(0, 1, a * l1);
            D1.set(1, 0, b * l2);
            D1.set(1, 1, (1.0 - b) * l2);
        } else {
            D0 = new Matrix(2, 2);
            D0.set(0, 0, -l1);
            D0.set(0, 1, l1);
            D0.set(1, 0, 0.0);
            D0.set(1, 1, -l2);

            D1 = new Matrix(2, 2);
            D1.set(0, 0, 0.0);
            D1.set(0, 1, 0.0);
            D1.set(1, 0, p * l2);
            D1.set(1, 1, (1.0 - p) * l2);
        }

        return new Pair<Matrix, Matrix>(D0, D1);
    }
}
