/**
 * @file Characteristic maximum M_K for order statistics in Fork-Join analysis
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import jline.util.Maths;

public final class FJ_char_max {
    private FJ_char_max() {}

    /**
     * Characteristic maximum M_K for exponential distribution.
     */
    public static CharMaxResult fj_char_max_exp(int K, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (mu <= 0.0) {
            throw new IllegalArgumentException("Rate mu must be positive.");
        }

        double mK = Math.log((double) K) / mu;
        double MK = FJ_harmonic.fj_harmonic(K) / mu;

        return new CharMaxResult(MK, mK);
    }

    /**
     * Characteristic maximum M_K for Erlang-k distribution.
     */
    public static CharMaxResult fj_char_max_erlang(int K, int kStages, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (kStages < 1) {
            throw new IllegalArgumentException("Erlang stages must be positive.");
        }
        if (mu <= 0.0) {
            throw new IllegalArgumentException("Rate mu must be positive.");
        }

        double target = 1.0 / K;
        double lo = 0.0;
        double hi = (double) kStages / mu + Math.log((double) K) / mu;
        while (erlangSurvival(hi, kStages, mu) > target) {
            hi *= 2.0;
        }

        for (int iter = 0; iter < 200; iter++) {
            double mid = (lo + hi) / 2.0;
            double sVal = erlangSurvival(mid, kStages, mu);
            if (Math.abs(sVal - target) < 1e-12) {
                break;
            }
            if (sVal > target) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double mK = (lo + hi) / 2.0;

        double MK = ((double) kStages / mu) * (1.0 + K * Math.exp(-mu * mK)
                * Math.pow(mu * mK, (double) kStages) / Maths.fact((double) kStages));

        return new CharMaxResult(MK, mK);
    }

    private static double erlangSurvival(double x, int k, double mu) {
        if (x <= 0.0) return 1.0;
        double S = 0.0;
        for (int i = 0; i < k; i++) {
            S += Math.pow(mu * x, (double) i) / Maths.fact((double) i);
        }
        return Math.exp(-mu * x) * S;
    }
}
