/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

/**
 * Fallback interval used when a QUEST stage test fails.
 *
 * <p>Three intervals are formed and the smallest interval containing all of them
 * is returned, which is the articles' prescription:
 *
 * <ul>
 * <li>Two symmetric intervals of half-width
 * {@code max(t_{1-alpha/2,K} sqrt(Ap/n), t_{1-alpha/2,K-1} sqrt(Np/n))}, one about
 * the full-sample quantile and one about the average batch quantile. Taking the
 * wider of the two variance components is deliberately conservative, since
 * neither can be trusted once a stage test has failed.</li>
 * <li>Willink's asymmetric interval, which corrects the batch quantiles for
 * skewness through the cube-root transform
 * {@code G(zeta) = ([1+6 gamma(zeta-gamma)]^(1/3)-1)/(2 gamma)} with
 * {@code gamma = skewness/(6 sqrt(K))}, evaluated at both t-quantiles so the two
 * arms differ.</li>
 * </ul>
 *
 * <p>The residual-autocorrelation factor {@code max(sqrt((1+phi1)/(1-phi1)), 1)}
 * is applied for {@link Sim_fquest}, where the batch quantiles come from one
 * sample path and can stay correlated, and omitted for {@link Sim_firquest},
 * where they come from independent replications and the article drops it.
 *
 * <p>Port of MATLAB sim_quest_heuristic_ci.m.
 *
 * <p>Reference: R. Willink, "A Confidence Interval and Test for the Mean of an
 * Asymmetric Distribution", Commun. Statist. Theory Methods 34, 2005.
 *
 * @since LINE 3.1.0
 */
final class SimQuestHeuristic {
    private SimQuestHeuristic() {}

    /**
     * Builds the fallback interval.
     *
     * @param bqe          the K batch quantile estimators, pooled over replications
     * @param centre       the full-sample empirical quantile
     * @param ap           batched STS area estimator
     * @param np           nonoverlapping batched quantile estimator
     * @param nstar        observations the estimators were computed from
     * @param alpha        significance level
     * @param useAutocorr  whether to apply the residual-autocorrelation factor
     * @return the two interval endpoints, lower then upper
     */
    static double[] interval(double[] bqe, double centre, double ap, double np,
                             int nstar, double alpha, boolean useAutocorr) {
        int K = bqe.length;
        if (K < 3) {
            throw new IllegalArgumentException(
                    "The heuristic interval needs at least 3 batch quantiles, got " + K);
        }

        double half = Math.max(SimDist.tinv(1.0 - alpha / 2.0, K) * Math.sqrt(ap / nstar),
                SimDist.tinv(1.0 - alpha / 2.0, K - 1) * Math.sqrt(np / nstar));

        double mean = 0.0;
        for (int j = 0; j < K; j++) {
            mean += bqe[j];
        }
        mean /= K;

        double s2 = 0.0;
        double s2tilde = 0.0;
        for (int j = 0; j < K; j++) {
            double d = bqe[j] - mean;
            s2 += d * d;
            double e = bqe[j] - centre;
            s2tilde += e * e;
        }
        s2 /= (K - 1);
        s2tilde /= (K - 1);

        double lower = Math.min(centre - half, mean - half);
        double upper = Math.max(centre + half, mean + half);
        if (!(s2 > 0.0)) {
            return new double[] {lower, upper};
        }

        double sd = Math.sqrt(s2);
        double skew = 0.0;
        for (int j = 0; j < K; j++) {
            double z = (bqe[j] - mean) / sd;
            skew += z * z * z;
        }
        skew *= (double) K / ((double) (K - 1) * (K - 2));
        double gamma = skew / (6.0 * Math.sqrt(K));

        double varphi = 1.0;
        if (useAutocorr) {
            double cov = 0.0;
            for (int j = 0; j < K - 1; j++) {
                cov += (bqe[j] - mean) * (bqe[j + 1] - mean);
            }
            double phi1 = cov / ((K - 1) * s2);
            if (Math.abs(phi1) < 1.0) {
                varphi = Math.max(Math.sqrt((1.0 + phi1) / (1.0 - phi1)), 1.0);
            }
        }

        double tq = SimDist.tinv(1.0 - alpha / 2.0, K - 1);
        double scale = varphi * Math.sqrt(s2tilde / K);
        double g1 = willink(tq, gamma) * scale;
        double g2 = willink(-tq, gamma) * scale;

        lower = Math.min(lower, Math.min(centre - g1, centre - g2));
        upper = Math.max(upper, Math.max(centre - g1, centre - g2));
        return new double[] {lower, upper};
    }

    private static double willink(double zeta, double gamma) {
        if (Math.abs(gamma) <= 0.001) {
            return zeta;
        }
        double arg = 1.0 + 6.0 * gamma * (zeta - gamma);
        // the cube root is taken on the reals, the argument may turn negative for a
        // strongly skewed and small batch sample
        double root = arg >= 0.0 ? Math.cbrt(arg) : -Math.cbrt(-arg);
        return (root - 1.0) / (2.0 * gamma);
    }
}
