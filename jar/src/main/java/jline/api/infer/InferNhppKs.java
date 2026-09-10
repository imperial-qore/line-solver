/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.infer;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;

/**
 * Kolmogorov-Smirnov tests for a non-homogeneous Poisson arrival process.
 *
 * <p>THE CONDITIONAL-UNIFORM TRANSFORMATION. Conditional on the number of
 * arrivals in the interval, the arrival times of an NHPP are the order
 * statistics of iid variables with cdf Lambda(t)/Lambda(T). Mapping the data
 * through that cdf turns ANY NHPP, whatever its rate, into iid uniforms, so one
 * KS test covers every rate function.
 *
 * <p>WHY THE PLAIN TEST IS WEAK, AND WHAT FIXES IT. The CU KS test has
 * "remarkably little power" against non-exponential interarrival times: it looks
 * at the POSITIONS of the points, and those stay nearly uniform for many
 * non-Poisson processes. Lewis (1965) applies the Durbin (1961) transformation
 * first -- reorder the GAPS ascending, rescale each by how many gaps remain,
 * cumulate -- which turns a difference in the gap DISTRIBUTION into a difference
 * in position. Measured on 400 replications of an Erlang-4 renewal process, the
 * CU test rejects at its own size while the Lewis test rejects essentially
 * always.
 *
 * <p>Port of MATLAB infer_nhpp_ks.m.
 *
 * <p>Reference: S.-H. Kim, W. Whitt (2014). Are call center and hospital
 * arrivals well modeled by nonhomogeneous Poisson processes? Manufacturing and
 * Service Operations Management 16(3), 464-480; J. Durbin (1961), Biometrika 48,
 * 41-55; P. A. W. Lewis (1965), JRSS B 27, 417-432.
 *
 * @since LINE 3.1.0
 */
public final class InferNhppKs {

    private InferNhppKs() {
    }

    /**
     * The Lewis (Durbin-transformed) test with a constant rate.
     *
     * @param times the arrival times
     * @param T     the right end of the observation interval
     * @return map with statistic, pvalue and n
     */
    public static Map<String, Double> infer_nhpp_ks(double[] times, double T) {
        return infer_nhpp_ks(times, T, null, "lewis", 0.0);
    }

    /**
     * @param times   the arrival times, within [T0,T]
     * @param T       the right end of the observation interval
     * @param cumRate the cumulative rate Lambda(t), or null for a constant rate
     * @param method  "cu" or "lewis"
     * @param T0      the left end of the observation interval
     * @return map with statistic, pvalue and n
     */
    public static Map<String, Double> infer_nhpp_ks(double[] times, double T,
                                                     DoubleUnaryOperator cumRate, String method,
                                                     double T0) {
        double[] t = Arrays.stream(times).filter(x -> x >= T0 && x <= T).sorted().toArray();
        int n = t.length;
        if (n < 2) {
            throw new RuntimeException("infer_nhpp_ks: at least two arrivals are needed to test");
        }
        double[] u = new double[n];
        if (cumRate == null) {
            for (int i = 0; i < n; i++) {
                u[i] = (t[i] - T0) / (T - T0);
            }
        } else {
            double lo = cumRate.applyAsDouble(T0);
            double hi = cumRate.applyAsDouble(T);
            if (hi <= lo) {
                throw new RuntimeException(
                        "infer_nhpp_ks: the cumulative rate must increase over the interval");
            }
            for (int i = 0; i < n; i++) {
                u[i] = (cumRate.applyAsDouble(t[i]) - lo) / (hi - lo);
            }
        }
        for (int i = 0; i < n; i++) {
            u[i] = Math.min(Math.max(u[i], 0.0), 1.0);
        }

        double[] s;
        String m = method == null ? "lewis" : method.toLowerCase();
        if ("cu".equals(m)) {
            s = u;
        } else if ("lewis".equals(m)) {
            // The Durbin (1961) transformation: gaps, sorted ascending, each
            // rescaled by how many gaps remain, then cumulated.
            double[] v = u.clone();
            Arrays.sort(v);
            double[] gaps = new double[n + 1];
            gaps[0] = v[0];
            for (int i = 1; i < n; i++) {
                gaps[i] = v[i] - v[i - 1];
            }
            gaps[n] = 1.0 - v[n - 1];
            Arrays.sort(gaps);
            s = new double[n];
            double prev = 0.0;
            double acc = 0.0;
            for (int i = 0; i <= n; i++) {
                double c = (n + 1 - i) * (gaps[i] - prev);
                prev = gaps[i];
                if (i < n) {
                    acc += c;
                    s[i] = Math.min(Math.max(acc, 0.0), 1.0);
                }
            }
        } else {
            throw new RuntimeException("infer_nhpp_ks: method must be 'cu' or 'lewis'");
        }

        double d = ksStatistic(s);
        Map<String, Double> res = new HashMap<String, Double>();
        res.put("statistic", d);
        res.put("pvalue", ksPvalue(d, n));
        res.put("n", (double) n);
        return res;
    }

    /** Two-sided KS distance between the sample and the uniform cdf. */
    private static double ksStatistic(double[] u) {
        double[] v = u.clone();
        Arrays.sort(v);
        int n = v.length;
        double d = 0.0;
        for (int i = 0; i < n; i++) {
            d = Math.max(d, (i + 1.0) / n - v[i]);
            d = Math.max(d, v[i] - (double) i / n);
        }
        return d;
    }

    /**
     * Asymptotic Kolmogorov p-value with the small-sample correction of
     * Stephens: the effective argument is (sqrt(n)+0.12+0.11/sqrt(n))D.
     */
    private static double ksPvalue(double d, int n) {
        if (n <= 0) {
            return 1.0;
        }
        double x = (Math.sqrt(n) + 0.12 + 0.11 / Math.sqrt(n)) * d;
        if (x <= 0) {
            return 1.0;
        }
        double q = 0.0;
        for (int k = 1; k <= 100; k++) {
            q += ((k % 2 == 1) ? 1 : -1) * Math.exp(-2.0 * k * k * x * x);
        }
        return Math.min(Math.max(2.0 * q, 0.0), 1.0);
    }
}
