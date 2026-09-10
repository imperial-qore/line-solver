/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;

/**
 * Two-moment approximation for the maximum of n iid non-negative variables.
 *
 * <p>THE SHAPE OF THE ANSWER. For a law with an exponential-like tail the
 * maximum of n samples grows like c~^2 (log n + ...): doubling n ADDS a
 * constant, it does not scale the answer. The two moments buy the SLOPE and an
 * offset: x_n(q) = c~^2[log(n eta) - log log(1/q)], E[M_n] = c~^2[log(n eta) +
 * gamma], with c~^2 = cs2 and eta = (cs2+1)/(2cs2^2) for cs2 &gt;= 1 (the H2
 * representative), and c~^2 = sqrt(cs2), eta = exp((1-sqrt(cs2))/sqrt(cs2)) below
 * (the shifted exponential).
 *
 * <p>WHEN NOT TO USE IT: n must pass n* ~ cs2/q, because with a highly variable
 * law only about n p of the samples can contend for the maximum. Measured
 * against exact maxima the closed form is within a few percent for n &gt;= 100 at
 * cs2 = 4 and 16, and useless at n = 10 for cs2 = 16 -- exactly what n* predicts.
 *
 * <p>AND WHEN TWO MOMENTS ARE NOT ENOUGH: below cs2 = 1 the maximum is genuinely
 * family-dependent, an Erlang and a shifted exponential with the same two
 * moments differing by tens of percent and diverging as n grows.
 *
 * <p>Port of MATLAB qsys_maxima_twomoment.m.
 *
 * <p>Reference: C. Crow, D. Goldberg, W. Whitt (2007). Two-moment
 * approximations for maxima. Operations Research 55(3), 532-548.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_maxima_twomoment {

    private Qsys_maxima_twomoment() {
    }

    private static final double EULER = 0.5772156649015329;

    /**
     * The mean of the maximum.
     *
     * @param n    the number of samples
     * @param mean the mean of the underlying law
     * @param cs2  its squared coefficient of variation
     * @return map with value, slope, eta, threshold, reliable, exactFittedValue
     */
    public static Map<String, Double> qsys_maxima_twomoment(int n, double mean, double cs2) {
        return qsys_maxima_twomoment(n, mean, cs2, Double.NaN, true);
    }

    /**
     * @param n           the number of samples
     * @param mean        the mean of the underlying law
     * @param cs2         its squared coefficient of variation
     * @param q           a quantile level in (0,1); NaN returns the mean
     * @param exactFitted also compute the maximum exactly from the fitted law
     * @return map with value, slope, eta, threshold, reliable (1 or 0) and,
     *         when requested, exactFittedValue
     */
    public static Map<String, Double> qsys_maxima_twomoment(int n, double mean, double cs2,
                                                             double q, boolean exactFitted) {
        if (n < 1) {
            throw new RuntimeException("qsys_maxima_twomoment: at least one sample is required");
        }
        if (mean <= 0) {
            throw new RuntimeException("qsys_maxima_twomoment: the mean must be positive");
        }
        if (cs2 <= 0) {
            throw new RuntimeException("qsys_maxima_twomoment: the SCV must be positive");
        }
        boolean wantMean = Double.isNaN(q);
        if (!wantMean && !(q > 0 && q < 1)) {
            throw new RuntimeException(
                    "qsys_maxima_twomoment: the quantile level must lie in (0,1)");
        }
        double ct;
        double eta;
        if (cs2 >= 1) {
            ct = cs2;
            eta = (cs2 + 1.0) / (2.0 * cs2 * cs2);
        } else {
            ct = Math.sqrt(cs2);
            eta = Math.exp((1.0 - Math.sqrt(cs2)) / Math.sqrt(cs2));
        }
        double inner = wantMean ? Math.log(n * eta) + EULER
                                : Math.log(n * eta) - Math.log(Math.log(1.0 / q));
        Map<String, Double> res = new HashMap<String, Double>();
        res.put("value", mean * ct * inner);
        res.put("slope", mean * ct);
        res.put("eta", eta);
        double threshold = cs2 / (wantMean ? 0.5 : q);
        res.put("threshold", threshold);
        res.put("reliable", n >= threshold ? 1.0 : 0.0);

        if (exactFitted) {
            // Fit the representative law and compute the maximum from F^n.
            double p1 = 0.0;
            double l1 = 0.0;
            double l2 = 0.0;
            double d = 0.0;
            double m = mean;
            double hi;
            if (cs2 >= 1) {
                p1 = 0.5 * (1.0 + Math.sqrt((cs2 - 1.0) / (cs2 + 1.0)));
                l1 = 2.0 * p1 / mean;
                l2 = 2.0 * (1.0 - p1) / mean;
                hi = 40.0 * mean * Math.max(cs2, 1.0);
            } else {
                d = mean * (1.0 - Math.sqrt(cs2));
                m = mean * Math.sqrt(cs2);
                hi = d + 40.0 * m;
            }
            int gn = 200000;
            double h = hi / gn;
            if (wantMean) {
                double acc = 0.0;
                for (int i = 0; i <= gn; i++) {
                    double t = i * h;
                    double cc = cs2 >= 1 ? p1 * Math.exp(-l1 * t) + (1 - p1) * Math.exp(-l2 * t)
                                         : (t <= d ? 1.0 : Math.exp(-(t - d) / m));
                    double v = 1.0 - Math.pow(1.0 - cc, n);
                    acc += (i == 0 || i == gn) ? v / 2.0 : v;
                }
                res.put("exactFittedValue", acc * h);
            } else {
                double val = hi;
                for (int i = 0; i <= gn; i++) {
                    double t = i * h;
                    double cc = cs2 >= 1 ? p1 * Math.exp(-l1 * t) + (1 - p1) * Math.exp(-l2 * t)
                                         : (t <= d ? 1.0 : Math.exp(-(t - d) / m));
                    if (Math.pow(1.0 - cc, n) >= q) {
                        val = t;
                        break;
                    }
                }
                res.put("exactFittedValue", val);
            }
        }
        return res;
    }
}
