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
 * Diffusion approximation for the G/GI/n/m queue.
 *
 * <p>THE APPROXIMATION IS ONE DIFFUSION WITH TWO REGIONS. Below the staffing
 * level the queue behaves like an infinite-server system, whose limit is NORMAL
 * with variance-to-mean ratio the ASYMPTOTIC PEAKEDNESS
 *
 * <pre>  z = 1 + (ca^2 - 1) omega_G,  omega_G = int G^c(x)^2 dx / int G^c(x) dx</pre>
 *
 * (eqs. 1.6-1.7); above it the queue behaves like a single-server queue, whose
 * limit is EXPONENTIAL with variability v = (ca^2 + cs^2)/2 (eq. 3.7). The
 * steady-state law is a normal piece spliced to an exponential piece, and every
 * measure is an integral of that density (eq. 3.14).
 *
 * <p>WHAT z SAYS. The service-time distribution enters the delay probability
 * ONLY through omega_G, which is 1 for deterministic service, 1/2 for
 * exponential, and falls toward 0 as service gets more variable. At ca^2 = 1 the
 * delay probability does not depend on the service law at all (z = 1), which is
 * the long-standing M/GI/n-by-M/M/n approximation; away from ca^2 = 1 it does.
 *
 * <p>With m infinite the delay probability is alpha(beta/sqrt(z)) for the
 * Halfin-Whitt function alpha (eq. 3.10), so this generalizes
 * {@link Qsys_mmk_qed}.
 *
 * <p>Port of MATLAB qsys_ggnm_diffusion.m.
 *
 * <p>Reference: W. Whitt (2004). A diffusion approximation for the G/GI/n/m
 * queue. Operations Research 52(6), 922-941.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_ggnm_diffusion {

    private Qsys_ggnm_diffusion() {
    }

    /**
     * The queue with exponential service.
     *
     * @param lambda arrival rate
     * @param mu     service rate of one server
     * @param n      number of servers
     * @param m      extra waiting spaces, {@code Double.POSITIVE_INFINITY} if unbounded
     * @param ca     coefficient of variation of the interarrival time
     * @param cs     coefficient of variation of the service time
     * @return the steady-state measures
     */
    public static Map<String, Double> qsys_ggnm_diffusion(double lambda, double mu, int n, double m,
                                                           double ca, double cs) {
        return qsys_ggnm_diffusion(lambda, mu, n, m, ca, cs, null, 1e-12, 4000);
    }

    /**
     * @param lambda      arrival rate
     * @param mu          service rate of one server
     * @param n           number of servers
     * @param m           extra waiting spaces, infinite if unbounded
     * @param ca          coefficient of variation of the interarrival time
     * @param cs          coefficient of variation of the service time
     * @param serviceCcdf G^c(x) = P(S &gt; x), or null for the exponential of rate mu
     * @param tol         service-tail cut for the peakedness integral
     * @param panels      Simpson panels for it
     * @return map with beta, gamma, peakedness, peakednessWeight, variability,
     *         probDelay, probBlock, meanQueueLength, meanNumber, meanWait,
     *         utilization, throughput and trafficIntensity
     */
    public static Map<String, Double> qsys_ggnm_diffusion(double lambda, double mu, int n, double m,
                                                           double ca, double cs,
                                                           DoubleUnaryOperator serviceCcdf,
                                                           double tol, int panels) {
        if (lambda <= 0 || mu <= 0) {
            throw new RuntimeException(
                    "qsys_ggnm_diffusion: the arrival and service rates must be positive");
        }
        if (n < 1) {
            throw new RuntimeException("qsys_ggnm_diffusion: the number of servers n must be at least 1");
        }
        if (m < 0) {
            throw new RuntimeException(
                    "qsys_ggnm_diffusion: the number of extra waiting spaces m must be non-negative");
        }
        double ca2 = ca * ca;
        double cs2 = cs * cs;
        double ES = 1.0 / mu;
        double rho = lambda / (n * mu);
        double beta = Math.sqrt(n) * (1.0 - rho);                       // eq. (0.1)
        double gamma = Double.isInfinite(m) ? Double.POSITIVE_INFINITY : m / Math.sqrt(n);

        double omega = serviceCcdf == null ? 0.5 : omega(serviceCcdf, ES, tol, panels);
        double z = 1.0 + (ca2 - 1.0) * omega;                           // eq. (1.6)
        if (z <= 0) {
            throw new RuntimeException("qsys_ggnm_diffusion: the asymptotic peakedness came out "
                    + "non-positive; check ca and the service ccdf");
        }
        double v = (ca2 + cs2) / 2.0;                                   // eq. (3.7), weight w = 1
        double b = beta / Math.sqrt(z);
        double r = beta / v;

        // Mass on the exponential piece. The tail factor is negative together
        // with r when the queue is overloaded, so the ratio stays positive.
        double tail = Double.isInfinite(gamma) ? 1.0 : -Math.expm1(-r * gamma);
        double alpha;
        double meanAbove;
        double densityAtTop;
        if (Math.abs(r) < 1e-14) {
            // beta = 0: the exponential piece degenerates to a uniform on [0,gamma].
            if (Double.isInfinite(gamma)) {
                throw new RuntimeException("qsys_ggnm_diffusion: with beta = 0 the queue needs a "
                        + "finite waiting room to be stable");
            }
            alpha = 1.0 / (1.0 + Phi(b) / (phi(b) * gamma / Math.sqrt(z)));
            meanAbove = gamma / 2.0;
            densityAtTop = alpha / gamma;
        } else {
            alpha = 1.0 / (1.0 + b * Phi(b) / (phi(b) * tail));
            if (Double.isInfinite(gamma)) {
                meanAbove = 1.0 / r;
                densityAtTop = 0.0;
            } else {
                double e = Math.exp(-r * gamma);
                meanAbove = (1.0 / r - (gamma + 1.0 / r) * e) / tail;
                densityAtTop = alpha * r * e / tail;
            }
        }

        // Mean of the normal piece, N(-beta, z) conditioned below 0.
        double meanBelow = -beta - Math.sqrt(z) * phi(b) / Phi(b);
        double meanScaled = (1.0 - alpha) * meanBelow + alpha * meanAbove;

        // Eq. (7.5): the loss rate of the diffusion at the upper boundary,
        // divided by the arrival rate, is the density there times v / sqrt(n).
        double probBlock = Double.isInfinite(gamma) ? 0.0
                : Math.min(Math.max(densityAtTop * v / Math.sqrt(n), 0.0), 1.0);
        double meanQueue = Math.sqrt(n) * alpha * meanAbove;
        double throughput = lambda * (1.0 - probBlock);

        Map<String, Double> res = new HashMap<String, Double>();
        res.put("beta", beta);
        res.put("gamma", gamma);
        res.put("peakedness", z);
        res.put("peakednessWeight", omega);
        res.put("variability", v);
        res.put("probDelay", alpha);
        res.put("probBlock", probBlock);
        res.put("meanQueueLength", meanQueue);
        res.put("meanNumber", n + Math.sqrt(n) * meanScaled);
        res.put("meanWait", throughput > 0 ? meanQueue / throughput : 0.0);
        res.put("utilization", Math.min(rho, 1.0));
        res.put("throughput", throughput);
        res.put("trafficIntensity", rho);
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

    /**
     * omega_G = int G^c(x)^2 dx / int G^c(x) dx of eq. (1.7), by Simpson on a
     * grid cut where the ccdf is negligible. The denominator is E[S], so only
     * the numerator is integrated.
     */
    private static double omega(DoubleUnaryOperator ccdf, double ES, double tol, int panels) {
        double hi = 1.0;
        while (ccdf.applyAsDouble(hi) > tol) {
            hi *= 2.0;
            if (hi > 1e12) {
                throw new RuntimeException("qsys_ggnm_diffusion: the service ccdf does not decay, "
                        + "so its peakedness is undefined");
            }
        }
        double h = hi / panels;
        double g0 = ccdf.applyAsDouble(0.0);
        double gn = ccdf.applyAsDouble(hi);
        double sum = g0 * g0 + gn * gn;
        for (int i = 1; i < panels; i++) {
            double g = ccdf.applyAsDouble(i * h);
            sum += (i % 2 == 1 ? 4.0 : 2.0) * g * g;
        }
        return (h / 3.0 * sum) / ES;
    }
}
