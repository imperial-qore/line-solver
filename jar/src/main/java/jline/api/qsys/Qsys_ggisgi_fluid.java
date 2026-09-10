/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.function.DoubleUnaryOperator;

/**
 * Steady state of the G/GI/s+GI fluid model.
 *
 * <p>Scale the content by s and let s grow. Customers become quanta of fluid but
 * their sojourns do not shrink, so the ages survive the limit: the state is the
 * density b(x) of fluid that has been IN SERVICE for time x and the density q(x)
 * of fluid that has been WAITING for time x. With rho = lambda/(s*mu),
 *
 * <ul>
 *   <li>rho &lt;= 1: b(x) = rho G^c(x), q = 0, no wait and no abandonment;</li>
 *   <li>rho &gt; 1: b(x) = G^c(x), q(x) = rho F^c(x) on [0,w] and 0 beyond,</li>
 * </ul>
 *
 * <p>the queue boundary w solving F^c(w) = 1/rho (eq. 3.6). That one equation
 * carries the whole overloaded regime: fluid that survives its patience for w
 * enters service, so the surviving fraction F^c(w) must equal the fraction 1/rho
 * the servers can absorb.
 *
 * <p>WHAT THE DISTRIBUTIONS CONTRIBUTE (Corollary 3.1): the rates and the number
 * in service depend on G and F only through their means; the wait, the queue
 * content and its age profile depend on F beyond its mean but on G only through
 * its mean. Neither s nor anything about the arrival process beyond its rate
 * appears, which is why this model says everything about the overloaded regime
 * and nothing about the QED one.
 *
 * <p>Port of MATLAB qsys_ggisgi_fluid.m.
 *
 * <p>Reference: W. Whitt (2006). Fluid models for multiserver queues with
 * abandonments. Operations Research 54(1), 37-54.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_ggisgi_fluid {

    private Qsys_ggisgi_fluid() {
    }

    /** Bisection tolerance for the queue boundary w. */
    public static final double DEFAULT_TOL = 1e-12;

    /**
     * Steady state with the exponential service ccdf implied by mu.
     *
     * @param lambda       arrival rate
     * @param mu           service rate of one server
     * @param s            number of servers
     * @param patienceCcdf F^c(t) = P(patience &gt; t)
     * @return the fluid steady state
     */
    public static QsysFluidAbandonResult qsys_ggisgi_fluid(double lambda, double mu, int s,
                                                            DoubleUnaryOperator patienceCcdf) {
        return qsys_ggisgi_fluid(lambda, mu, s, patienceCcdf, null, null, DEFAULT_TOL,
                Double.NaN);
    }

    /**
     * Steady state, optionally with the two age densities.
     *
     * @param lambda       arrival rate
     * @param mu           service rate of one server
     * @param s            number of servers
     * @param patienceCcdf F^c(t) = P(patience &gt; t)
     * @param servingCcdf  G^c(x) = P(service &gt; x); null takes the exponential of rate mu
     * @param agePoints    ages at which to return the densities, or null
     * @param tol          bisection tolerance for w
     * @param maxTime      largest age searched for w; NaN grows the search automatically
     * @return the fluid steady state
     */
    public static QsysFluidAbandonResult qsys_ggisgi_fluid(double lambda, double mu, int s,
                                                            DoubleUnaryOperator patienceCcdf,
                                                            DoubleUnaryOperator servingCcdf,
                                                            double[] agePoints, double tol,
                                                            double maxTime) {
        if (lambda <= 0) {
            throw new RuntimeException("qsys_ggisgi_fluid: the arrival rate lambda must be positive");
        }
        if (mu <= 0) {
            throw new RuntimeException("qsys_ggisgi_fluid: the service rate mu must be positive");
        }
        if (s < 1) {
            throw new RuntimeException("qsys_ggisgi_fluid: the number of servers s must be at least 1");
        }
        if (patienceCcdf == null) {
            throw new RuntimeException("qsys_ggisgi_fluid: the patience ccdf must not be null");
        }
        final double fmu = mu;
        DoubleUnaryOperator gc = servingCcdf != null ? servingCcdf
                : new DoubleUnaryOperator() {
                    @Override
                    public double applyAsDouble(double x) {
                        return Math.exp(-fmu * x);
                    }
                };

        double rho = lambda / (s * mu);
        String regime;
        double w;
        double meanWait;
        double probAbandon;
        double meanWaitAbandon;
        if (rho <= 1.0) {
            // Underloaded and balanced, eq. (3.2): the queue is empty and the
            // model is the infinite-server fluid model.
            regime = Math.abs(rho - 1.0) <= Math.ulp(1.0) ? "balanced" : "underloaded";
            w = 0.0;
            meanWait = 0.0;
            probAbandon = 0.0;
            meanWaitAbandon = 0.0;
        } else {
            regime = "overloaded";
            w = invCcdf(patienceCcdf, 1.0 / rho, tol, maxTime);         // eq. (3.6)
            // Eq. (3.14): W = int_0^w F^c(t) dt = m_a F_e(w), over ALL fluid.
            meanWait = integral(patienceCcdf, 0.0, w);
            probAbandon = 1.0 - 1.0 / rho;
            // E[T | T <= w] = (W - w F^c(w)) / F(w) by parts, F^c(w) = 1/rho.
            meanWaitAbandon = (meanWait - w / rho) / probAbandon;
        }

        double queueLength = lambda * meanWait;                          // eq. (3.11), Little's law
        double inService = Math.min(lambda / mu, (double) s);
        double throughput = Math.min(lambda, s * mu);
        double[] serviceAge = null;
        double[] queueAge = null;
        double[] ages = null;
        if (agePoints != null && agePoints.length > 0) {
            ages = agePoints.clone();
            serviceAge = new double[ages.length];
            queueAge = new double[ages.length];
            double sigma = Math.min(rho, 1.0);          // rate into service, per server
            for (int i = 0; i < ages.length; i++) {
                serviceAge[i] = sigma * gc.applyAsDouble(ages[i]);
                queueAge[i] = (rho > 1.0 && ages[i] <= w)
                        ? rho * patienceCcdf.applyAsDouble(ages[i]) : 0.0;
            }
        }

        return new QsysFluidAbandonResult(regime, rho, w, meanWait, w, meanWaitAbandon,
                probAbandon, queueLength, inService, inService + queueLength, Math.min(rho, 1.0),
                throughput, lambda - throughput, ages, serviceAge, queueAge);
    }

    /**
     * Smallest w with F^c(w) = target, found by doubling then bisection. F^c is
     * non-increasing, so the doubling either brackets the crossing or proves
     * that the patience law never decays that far.
     */
    private static double invCcdf(DoubleUnaryOperator ccdf, double target, double tol,
                                  double maxTime) {
        if (ccdf.applyAsDouble(0.0) < target) {
            throw new RuntimeException(
                    "qsys_ggisgi_fluid: the patience ccdf is below 1/rho at t = 0, so it is not a ccdf");
        }
        double lo = 0.0;
        double hi;
        if (Double.isNaN(maxTime)) {
            hi = 1.0;
            while (ccdf.applyAsDouble(hi) > target) {
                hi *= 2.0;
                if (hi > 1e12) {
                    throw new RuntimeException("qsys_ggisgi_fluid: the patience ccdf never falls to "
                            + "1/rho, so the overloaded fluid model has no equilibrium: too little "
                            + "of the fluid is willing to abandon");
                }
            }
        } else {
            hi = maxTime;
            if (ccdf.applyAsDouble(hi) > target) {
                throw new RuntimeException(
                        "qsys_ggisgi_fluid: the patience ccdf is still above 1/rho at maxTime");
            }
        }
        while (hi - lo > tol * Math.max(1.0, hi)) {
            double mid = 0.5 * (lo + hi);
            if (ccdf.applyAsDouble(mid) > target) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return 0.5 * (lo + hi);
    }

    /**
     * Composite Simpson rule on a fixed fine grid: the integrand is a ccdf,
     * hence monotone and bounded, so a fixed grid is enough and is reproducible.
     */
    private static double integral(DoubleUnaryOperator f, double a, double b) {
        if (b <= a) {
            return 0.0;
        }
        int n = 2000;
        double h = (b - a) / n;
        double sum = f.applyAsDouble(a) + f.applyAsDouble(b);
        for (int i = 1; i < n; i++) {
            sum += (i % 2 == 1 ? 4.0 : 2.0) * f.applyAsDouble(a + i * h);
        }
        return h / 3.0 * sum;
    }
}
