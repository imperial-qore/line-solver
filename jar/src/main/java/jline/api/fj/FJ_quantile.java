/**
 * @file Quantile approximation for maximum of K random variables
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import java.util.function.DoubleUnaryOperator;

public final class FJ_quantile {
    private FJ_quantile() {}

    /**
     * Quantile approximation for maximum of K random variables (standard Gumbel).
     */
    public static double fj_quantile(int K, double q) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (q <= 0.0 || q >= 1.0) {
            throw new IllegalArgumentException("Quantile q must satisfy 0 < q < 1. Got q=" + String.format("%.4f", q) + ".");
        }
        return Math.log((double) K) - Math.log(Math.log(1.0 / q));
    }

    /**
     * Quantile of maximum of K random variables using inverse CDF.
     */
    public static double fj_quantile(int K, double q, DoubleUnaryOperator Finv) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (q <= 0.0 || q >= 1.0) {
            throw new IllegalArgumentException("Quantile q must satisfy 0 < q < 1. Got q=" + String.format("%.4f", q) + ".");
        }
        return Finv.applyAsDouble(Math.pow(q, 1.0 / K));
    }

    /**
     * Quantile approximation for maximum of K random variables (array version).
     */
    public static double[] fj_quantile(int K, double[] q) {
        double[] result = new double[q.length];
        for (int i = 0; i < q.length; i++) {
            result[i] = fj_quantile(K, q[i]);
        }
        return result;
    }
}
