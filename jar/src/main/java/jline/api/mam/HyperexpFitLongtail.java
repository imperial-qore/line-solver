/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mam;

import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;

/**
 * Fitting a hyperexponential to a long-tail distribution.
 *
 * <p>WHY MOMENTS ARE THE WRONG HANDLE. A Pareto law with tail index below 2 has
 * infinite variance, so no two- or three-moment fit exists at all; and even when
 * the moments are finite, matching them says nothing about the several ORDERS OF
 * MAGNITUDE of time scale over which a long-tail law acts. This procedure
 * matches the CCDF ITSELF at points spread across those decades.
 *
 * <p>THE RECURSION, with lambda_1 &lt; ... &lt; lambda_k. In the far tail only
 * the slowest component survives, so it is fitted there alone from the ccdf at
 * c_1 and b*c_1 (eqs. 4.4-4.5); subtract it and repeat one decade lower
 * (eqs. 4.6-4.11). The last component takes the remaining probability and its
 * rate follows from the ccdf at c_k (eqs. 4.12-4.14). This is Prony's method
 * applied to a ccdf.
 *
 * <p>DEFAULTS. (b, decade) = (1.5, 4) rather than the paper's illustrative
 * (2, 10): the fit is exact AT the fitting arguments and free between them, and
 * measured on a Weibull(0.3) the tighter grid cuts the worst between-point error
 * from about 54% to 12%, at the cost of more components.
 *
 * <p>Port of MATLAB hyperexp_fit_longtail.m.
 *
 * <p>Reference: A. Feldmann, W. Whitt (1998). Fitting mixtures of exponentials
 * to long-tail distributions to analyze network performance models. Performance
 * Evaluation 31, 245-279, Section 4.
 *
 * @since LINE 3.1.0
 */
public final class HyperexpFitLongtail {

    private HyperexpFitLongtail() {
    }

    /** The within-scale spacing b of the fitting pairs. */
    public static final double DEFAULT_B = 1.5;
    /** The ratio between successive fitting arguments. */
    public static final double DEFAULT_DECADE = 4.0;

    /**
     * The fit with the component count chosen automatically: one per decade
     * between the 0.9 quantile and the 1e-6 quantile, retrying with fewer when
     * the recursion runs out of probability near the body.
     *
     * @param ccdf F^c(t) = P(X &gt; t)
     * @return map with p, lambda, points (double[]), mean, targetMean,
     *         coverageLow, coverageHigh, maxRelError and maxRelErrorGrid
     */
    public static Map<String, Object> hyperexp_fit_longtail(DoubleUnaryOperator ccdf) {
        return hyperexp_fit_longtail(ccdf, DEFAULT_B, DEFAULT_DECADE);
    }

    /**
     * @param ccdf   F^c(t) = P(X &gt; t)
     * @param b      the within-scale spacing
     * @param decade the ratio between successive fitting arguments
     * @return the fit, keyed as in the MATLAB struct
     */
    public static Map<String, Object> hyperexp_fit_longtail(DoubleUnaryOperator ccdf, double b,
                                                             double decade) {
        double top = quantile(ccdf, 1e-6);
        double body = quantile(ccdf, 0.9);
        if (body <= 0 || top <= body) {
            throw new RuntimeException(
                    "hyperexp_fit_longtail: the ccdf gives no usable range of time scales");
        }
        int k0 = Math.max(2, (int) Math.round(Math.log(top / body) / Math.log(decade)) + 1);
        // The recursion needs each component to dominate at its own scale. Near
        // the body of a law with a lot of mass there (a Pareto, say) that fails
        // and the remaining probability runs out; back off one at a time.
        RuntimeException last = null;
        for (int k = k0; k >= 2; k--) {
            try {
                return hyperexp_fit_longtail_k(ccdf, k, top, b, decade);
            } catch (RuntimeException e) {
                last = e;
            }
        }
        throw new RuntimeException("hyperexp_fit_longtail: no component count admits the recursion; "
                + "the ccdf may not be long-tailed enough for this scheme"
                + (last == null ? "" : " (" + last.getMessage() + ")"));
    }

    /**
     * The recursion at a fixed component count.
     *
     * @param ccdf   F^c(t) = P(X &gt; t)
     * @param k      number of exponential components
     * @param c1     the largest fitting argument
     * @param b      the within-scale spacing
     * @param decade the ratio between successive fitting arguments
     * @return the fit, keyed as in the MATLAB struct
     */
    public static Map<String, Object> hyperexp_fit_longtail_k(DoubleUnaryOperator ccdf, int k,
                                                               double c1, double b, double decade) {
        if (k < 1) {
            throw new RuntimeException("hyperexp_fit_longtail: at least one component is required");
        }
        if (!(b > 1)) {
            throw new RuntimeException("hyperexp_fit_longtail: the spacing b must exceed 1");
        }
        if (decade <= b) {
            throw new RuntimeException("hyperexp_fit_longtail: the decade ratio must exceed the "
                    + "spacing b, or the fitting arguments would interleave");
        }
        double[] cs = new double[k];
        cs[0] = c1;
        for (int i = 1; i < k; i++) {
            cs[i] = cs[i - 1] / decade;
        }
        double[] p = new double[k];
        double[] lam = new double[k];
        for (int i = 0; i < k; i++) {
            double ci = cs[i];
            // Eqs. (4.6)-(4.7): what the already-fitted, slower components leave.
            double residC = ccdf.applyAsDouble(ci);
            double residBC = ccdf.applyAsDouble(b * ci);
            for (int j = 0; j < i; j++) {
                residC -= p[j] * Math.exp(-lam[j] * ci);
                residBC -= p[j] * Math.exp(-lam[j] * b * ci);
            }
            if (i + 1 < k) {
                if (residC <= 0 || residBC <= 0 || residC <= residBC) {
                    throw new RuntimeException("hyperexp_fit_longtail: the residual ccdf is not "
                            + "positive and decreasing at fitting argument " + ci + "; the "
                            + "recursion needs c_i/c_(i+1) >> b");
                }
                lam[i] = Math.log(residC / residBC) / ((b - 1.0) * ci);      // eq. (4.10)
                p[i] = residC * Math.exp(lam[i] * ci);                       // eq. (4.11)
            } else {
                // Eqs. (4.12)-(4.14): the last component takes the rest.
                double rest = 1.0;
                for (int j = 0; j < i; j++) {
                    rest -= p[j];
                }
                if (rest <= 0) {
                    throw new RuntimeException("hyperexp_fit_longtail: the fitted components "
                            + "already carry all the probability");
                }
                if (residC <= 0) {
                    throw new RuntimeException("hyperexp_fit_longtail: the residual ccdf has gone "
                            + "non-positive at the last fitting argument");
                }
                p[i] = rest;
                lam[i] = Math.log(p[i] / residC) / ci;                       // eq. (4.14)
            }
            if (lam[i] <= 0) {
                throw new RuntimeException("hyperexp_fit_longtail: a non-positive rate came out of "
                        + "the fit at argument " + ci);
            }
        }

        double coverageHigh = cs[0] * b;
        double mean = 0.0;
        for (int j = 0; j < k; j++) {
            mean += p[j] / lam[j];
        }
        // The target mean over the covered range: a k too small to reach the
        // body shows up here and nowhere else.
        int gn = 20000;
        double h = coverageHigh / gn;
        double acc = (ccdf.applyAsDouble(0.0) + ccdf.applyAsDouble(coverageHigh)) / 2.0;
        for (int i = 1; i < gn; i++) {
            acc += ccdf.applyAsDouble(i * h);
        }
        double maxErr = 0.0;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < 2; j++) {
                double t = j == 0 ? cs[i] : b * cs[i];
                double target = ccdf.applyAsDouble(t);
                if (target > 0) {
                    maxErr = Math.max(maxErr, Math.abs(fitted(p, lam, t) - target) / target);
                }
            }
        }
        // The fit is exact at the fitting arguments; this says whether it also
        // holds BETWEEN them.
        double gridErr = 0.0;
        double loLog = Math.log(cs[k - 1]);
        double hiLog = Math.log(coverageHigh);
        for (int i = 0; i < 200; i++) {
            double t = Math.exp(loLog + (hiLog - loLog) * i / 199.0);
            double target = ccdf.applyAsDouble(t);
            if (target > 1e-300) {
                gridErr = Math.max(gridErr, Math.abs(fitted(p, lam, t) - target) / target);
            }
        }

        Map<String, Object> res = new HashMap<String, Object>();
        res.put("p", p);
        res.put("lambda", lam);
        res.put("points", cs);
        res.put("mean", mean);
        res.put("targetMean", acc * h);
        res.put("coverageLow", cs[k - 1]);
        res.put("coverageHigh", coverageHigh);
        res.put("maxRelError", maxErr);
        res.put("maxRelErrorGrid", gridErr);
        return res;
    }

    /** The fitted ccdf at t. */
    private static double fitted(double[] p, double[] lam, double t) {
        double v = 0.0;
        for (int j = 0; j < p.length; j++) {
            v += p[j] * Math.exp(-lam[j] * t);
        }
        return v;
    }

    /** Smallest t with F^c(t) &lt;= prob, by doubling then bisection. */
    private static double quantile(DoubleUnaryOperator ccdf, double prob) {
        double hi = 1.0;
        while (ccdf.applyAsDouble(hi) > prob) {
            hi *= 2.0;
            if (hi > 1e15) {
                throw new RuntimeException(
                        "hyperexp_fit_longtail: the ccdf does not decay, so there is no tail to fit");
            }
        }
        double lo = 0.0;
        for (int i = 0; i < 200; i++) {
            double mid = 0.5 * (lo + hi);
            if (ccdf.applyAsDouble(mid) > prob) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }
}
