/**
 * @file Order statistics and G(K) bound factors for Fork-Join analysis
 *
 * Computes order statistics (CDF and expected value of the k-th smallest of K
 * i.i.d. random variables) and the G(K) scaling factors used in the
 * mean-variance approximation X_K^max ~ mu + sigma * G(K).
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.function.DoubleUnaryOperator;

import jline.util.Maths;

public final class FJ_order_stat {
    private FJ_order_stat() {}

    /**
     * Compute G(K) factors for expected maximum approximation.
     */
    public static double fj_gk_bound(int K, String type) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }

        String t = type.toLowerCase();
        if ("exp".equals(t)) {
            return FJ_harmonic.fj_harmonic(K) - 1.0;
        } else if ("uniform".equals(t)) {
            return Math.sqrt(3.0) * (K - 1) / (K + 1);
        } else if ("evd".equals(t)) {
            return Math.sqrt(6.0) * Math.log(K) / Math.PI;
        } else if ("bound".equals(t)) {
            return ((double) (K - 1)) / Math.sqrt(2.0 * K - 1.0);
        } else {
            throw new IllegalArgumentException("Unknown type: " + type + ". Valid: exp, uniform, evd, bound.");
        }
    }

    public static double fj_gk_bound(int K) {
        return fj_gk_bound(K, "exp");
    }

    /**
     * Compute all G(K) bound factors.
     */
    public static GKBoundResult fj_gk_bound_all(int K) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }

        return new GKBoundResult(
                K,
                FJ_harmonic.fj_harmonic(K) - 1.0,
                Math.sqrt(3.0) * (K - 1) / (K + 1),
                Math.sqrt(6.0) * Math.log(K) / Math.PI,
                ((double) (K - 1)) / Math.sqrt(2.0 * K - 1.0)
        );
    }

    /**
     * CDF of k-th order statistic of K i.i.d. random variables.
     */
    public static double fj_order_stat_cdf(double FXy, int k, int K) {
        if (k < 1 || k > K) {
            throw new IllegalArgumentException("k must satisfy 1 <= k <= K. Got k=" + k + ", K=" + K + ".");
        }

        if (k == K) {
            // Maximum (K-th order statistic)
            return Math.pow(FXy, K);
        }

        // General k-th order statistic
        double FYk = 0.0;
        for (int j = k; j <= K; j++) {
            FYk += Maths.binomialCoeff(K, j)
                    * Math.pow(FXy, j)
                    * Math.pow(1.0 - FXy, K - j);
        }
        return FYk;
    }

    /**
     * Expected value of the maximum of K i.i.d. random variables via numerical integration.
     */
    public static double fj_order_stat_expected_max(int K, DoubleUnaryOperator cdfFunc, double upperLimit) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }

        // Composite Simpson's rule
        int numPoints = 10001;
        double h = upperLimit / (numPoints - 1);
        double integral = 0.0;
        for (int i = 0; i < numPoints; i++) {
            double t = i * h;
            double cdfVal = cdfFunc.applyAsDouble(t);
            double integrandVal = 1.0 - Math.pow(cdfVal, K);

            double weight;
            if (i == 0 || i == numPoints - 1) {
                weight = 1.0;
            } else if (i % 2 == 1) {
                weight = 4.0;
            } else {
                weight = 2.0;
            }
            integral += weight * integrandVal;
        }
        integral *= h / 3.0;

        return integral;
    }

    /**
     * Expected value of the minimum of K i.i.d. random variables via numerical integration.
     */
    public static double fj_order_stat_expected_min(int K, DoubleUnaryOperator cdfFunc, double upperLimit) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }

        // Composite Simpson's rule
        int numPoints = 10001;
        double h = upperLimit / (numPoints - 1);
        double integral = 0.0;
        for (int i = 0; i < numPoints; i++) {
            double t = i * h;
            double cdfVal = cdfFunc.applyAsDouble(t);
            double integrandVal = Math.pow(1.0 - cdfVal, K);

            double weight;
            if (i == 0 || i == numPoints - 1) {
                weight = 1.0;
            } else if (i % 2 == 1) {
                weight = 4.0;
            } else {
                weight = 2.0;
            }
            integral += weight * integrandVal;
        }
        integral *= h / 3.0;

        return integral;
    }
}
