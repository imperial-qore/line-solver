/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.HashMap;
import java.util.Map;
import java.util.function.DoubleUnaryOperator;

import org.apache.commons.math3.special.Erf;

/**
 * Truncated Gaussian approximation (TGA-G) for the G/GI/n+GI queue.
 *
 * <p>A FLUID CENTRE PLUS A GAUSSIAN FLUCTUATION, TRUNCATED. In the
 * efficiency-driven regime (rho &gt; 1 fixed as n grows) the fluid limit gives
 * the centre -- every server busy, w = F^-1(1-1/rho), Q = lambda int_0^w F^c --
 * and the many-server central limit theorem gives a normal fluctuation of order
 * sqrt(n) around it, with
 *
 * <pre>  sigma_W^2 = [(ca^2-1) + (cs+1)rho] / (2 mu rho^2 f(w))</pre>
 *
 * (eq. 24). Adding fluid and fluctuation directly can produce negative queues
 * and negative waits, so BOTH ARE TRUNCATED at zero; that truncation is what
 * makes the formulas usable down to moderate overload, reportedly rho &gt; 1.02.
 *
 * <p>The three sources of variability enter separately, which is what lets the
 * exponential-service formula be generalized: the service law appears only as
 * the factor (cs+1)rho, which is 2rho at cs = 1.
 *
 * <p>Port of MATLAB qsys_ggingi_tga.m.
 *
 * <p>Reference: Y. Liu, W. Whitt, Y. Yu (2016). Approximations for
 * heavily-loaded G/GI/n+GI queues. Naval Research Logistics 63(3), 187-217.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_ggingi_tga {

    private Qsys_ggingi_tga() {
    }

    /**
     * @param lambda        arrival rate
     * @param mu            service rate of one server
     * @param n             number of servers
     * @param ca            coefficient of variation of the interarrival time
     * @param cs            coefficient of variation of the service time
     * @param patienceCcdf  F^c(x) = P(patience &gt; x)
     * @return the steady-state measures, keyed as in the MATLAB struct
     */
    public static Map<String, Double> qsys_ggingi_tga(double lambda, double mu, int n, double ca,
                                                       double cs, DoubleUnaryOperator patienceCcdf) {
        return qsys_ggingi_tga(lambda, mu, n, ca, cs, patienceCcdf, null, null);
    }

    /**
     * @param lambda       arrival rate
     * @param mu           service rate of one server
     * @param n            number of servers
     * @param ca           coefficient of variation of the interarrival time
     * @param cs           coefficient of variation of the service time
     * @param patienceCcdf F^c(x) = P(patience &gt; x)
     * @param patiencePdf  the patience density, or null to difference the ccdf
     * @param serviceCcdf  G^c(x), used only in the underloaded branch
     * @return map with regime (1 overloaded, 0 underloaded), trafficIntensity,
     *         fluidWait, fluidQueueLength, meanWait, varWait, meanQueueLength,
     *         varQueueLength, meanNumberInService, meanNumber, probDelay,
     *         probAbandon, sigmaW and sigmaX
     */
    public static Map<String, Double> qsys_ggingi_tga(double lambda, double mu, int n, double ca,
                                                       double cs, DoubleUnaryOperator patienceCcdf,
                                                       DoubleUnaryOperator patiencePdf,
                                                       DoubleUnaryOperator serviceCcdf) {
        if (lambda <= 0 || mu <= 0) {
            throw new RuntimeException("qsys_ggingi_tga: the arrival and service rates must be positive");
        }
        if (n < 1) {
            throw new RuntimeException("qsys_ggingi_tga: the number of servers n must be at least 1");
        }
        double ca2 = ca * ca;
        double rho = lambda / (n * mu);
        double lamPn = lambda / n;
        final DoubleUnaryOperator ccdf = patienceCcdf;
        DoubleUnaryOperator pdf = patiencePdf;
        if (pdf == null) {
            pdf = new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double x) {
                    double h = 1e-6;
                    return Math.max(0.0, (ccdf.applyAsDouble(Math.max(0.0, x - h))
                            - ccdf.applyAsDouble(x + h)) / (2 * h));
                }
            };
        }

        Map<String, Double> res = new HashMap<String, Double>();
        res.put("trafficIntensity", rho);
        if (rho <= 1.0) {
            // Underloaded: no queue in the limit; the content is normal with the
            // infinite-server variance (eq. 10).
            double omega = 0.5;
            if (serviceCcdf != null) {
                double hi = invCcdf(serviceCcdf, 1e-12);
                final DoubleUnaryOperator g = serviceCcdf;
                omega = simpson(new DoubleUnaryOperator() {
                    @Override
                    public double applyAsDouble(double x) {
                        double v = g.applyAsDouble(x);
                        return v * v;
                    }
                }, 0.0, hi) * mu;
            }
            res.put("regime", 0.0);
            res.put("fluidWait", 0.0);
            res.put("fluidQueueLength", 0.0);
            res.put("meanWait", 0.0);
            res.put("varWait", 0.0);
            res.put("meanQueueLength", 0.0);
            res.put("varQueueLength", 0.0);
            res.put("meanNumberInService", lambda / mu);
            res.put("meanNumber", lambda / mu);
            res.put("probDelay", 0.0);
            res.put("probAbandon", 0.0);
            res.put("sigmaW", 0.0);
            res.put("sigmaX", Math.sqrt((lambda / mu) * (1.0 + (ca2 - 1.0) * omega)));
            return res;
        }

        // Overloaded: the fluid centre of Theorem 2.1(b).
        final double w = invCcdf(patienceCcdf, 1.0 / rho);
        double fw = pdf.applyAsDouble(w);
        if (fw <= 0) {
            throw new RuntimeException("qsys_ggingi_tga: the patience density vanishes at the fluid "
                    + "waiting time, so the Gaussian correction is undefined there");
        }
        double qPerServer = lamPn * simpson(patienceCcdf, 0.0, w);

        // Eq. (24): the service law enters only through the (cs+1)rho term.
        double sigmaW2 = ((ca2 - 1.0) + (cs + 1.0) * rho) / (2.0 * mu * rho * rho * fw);
        final double fca2 = ca2;
        double sigmaX2 = mu * mu * sigmaW2 + lamPn * simpson(new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double x) {
                double v = ccdf.applyAsDouble(x);
                return v * (1.0 + (fca2 - 1.0) * v);
            }
        }, 0.0, w);
        double sigmaW = Math.sqrt(Math.max(sigmaW2, 0.0));
        double sigmaX = Math.sqrt(Math.max(sigmaX2, 0.0));

        final double aW = Math.sqrt(n) * w / sigmaW;              // eq. (21)
        double aX = Math.sqrt(n) * qPerServer / sigmaX;           // eq. (19)
        double vW = truncatedVariance(aW);
        double vX = truncatedVariance(aX);

        res.put("regime", 1.0);
        res.put("fluidWait", w);
        res.put("fluidQueueLength", n * qPerServer);
        res.put("meanWait", w * (Phi(aW) + phi(aW) / aW));
        res.put("varWait", (sigmaW * sigmaW / n) * vW);
        res.put("meanQueueLength", n * qPerServer * (Phi(aX) + phi(aX) / aX));
        res.put("varQueueLength", n * sigmaX * sigmaX * vX);
        // E[B] = E[min(X_n,n)]: every server busy but for the lower tail.
        double meanB = n - Math.sqrt(n) * sigmaX * (phi(aX) - aX * (1.0 - Phi(aX)));
        res.put("meanNumberInService", meanB);
        res.put("meanNumber", meanB + res.get("meanQueueLength"));
        res.put("probDelay", Phi(aW));                            // eq. (22)
        // Eq. (23): a customer abandons when its patience falls short of its wait.
        final DoubleUnaryOperator fpdf = pdf;
        double pa = simpson(new DoubleUnaryOperator() {
            @Override
            public double applyAsDouble(double x) {
                return (1.0 - Phi(aW * (x / w - 1.0))) * fpdf.applyAsDouble(x);
            }
        }, 0.0, Math.max(20 * w, w + 20));
        res.put("probAbandon", Math.min(Math.max(pa, 0.0), 1.0));
        res.put("sigmaW", sigmaW);
        res.put("sigmaX", sigmaX);
        return res;
    }

    /** Standard normal density. */
    private static double phi(double x) {
        return Math.exp(-x * x / 2.0) / Math.sqrt(2.0 * Math.PI);
    }

    /** Standard normal cdf, through erfc. */
    private static double Phi(double x) {
        return Erf.erfc(-x / Math.sqrt(2.0)) / 2.0;
    }

    /** Variance of max(Z,-a) for a standard normal Z. */
    private static double truncatedVariance(double a) {
        double m1 = phi(a) - a * (1.0 - Phi(a));
        double m2 = Phi(a) - a * phi(a) + a * a * (1.0 - Phi(a));
        return Math.max(0.0, m2 - m1 * m1);
    }

    /** Composite Simpson rule on a fixed even panel count. */
    private static double simpson(DoubleUnaryOperator f, double a, double b) {
        if (b <= a) {
            return 0.0;
        }
        int m = 2000;
        double h = (b - a) / m;
        double sum = f.applyAsDouble(a) + f.applyAsDouble(b);
        for (int i = 1; i < m; i++) {
            sum += (i % 2 == 1 ? 4.0 : 2.0) * f.applyAsDouble(a + i * h);
        }
        return h / 3.0 * sum;
    }

    /** Smallest w with F^c(w) = target, by doubling then bisection. */
    private static double invCcdf(DoubleUnaryOperator ccdf, double target) {
        double lo = 0.0;
        double hi = 1.0;
        while (ccdf.applyAsDouble(hi) > target) {
            hi *= 2.0;
            if (hi > 1e12) {
                throw new RuntimeException("qsys_ggingi_tga: the patience ccdf never falls to "
                        + "1/rho, so the overloaded model has no fluid equilibrium");
            }
        }
        while (hi - lo > 1e-12 * Math.max(1.0, hi)) {
            double mid = 0.5 * (lo + hi);
            if (ccdf.applyAsDouble(mid) > target) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }
}
