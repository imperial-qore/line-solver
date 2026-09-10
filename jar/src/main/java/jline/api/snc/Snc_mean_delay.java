/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Upper bound on the mean delay, from integrating the delay tail bound.
 *
 * <p>For a nonnegative delay {@code E[D] = int_0^inf P{D>d} dd}, so integrating
 * the tail bound of {@link Snc_bound_delay} bounds the MEAN. At fixed theta the
 * bound is {@code K*exp(-a*d)} with {@code a = theta*rhoS} and
 * {@code K = exp(theta*(sigmaA+sigmaS))/(1-exp(-theta*(rhoS-rhoA)))}, so,
 * clipping the bound at 1 where it exceeds it, the integral is available in
 * CLOSED FORM: {@code (log(K)+1)/a} when K &gt;= 1 and {@code K/a} otherwise. No
 * quadrature is involved, so the result is a bound and not a bound plus a
 * discretization error.</p>
 *
 * <p>IT IS A LOOSE MEAN BOUND AND THAT IS INHERENT: on the M/M/1 read in job
 * units it returns 2.4x the exact 1/(mu-lambda) at rho = 0.1 and 10.4x at
 * rho = 0.95, because the prefactor of the tail bound, not its decay rate,
 * dominates an integral over the whole axis. Use {@link Snc_perc_delay} when
 * the quantile is what matters.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_mean_delay.m. This is the entry point
 * SolverBA calls for the {@code snc.upper} response-time column.</p>
 */
public final class Snc_mean_delay {
    private Snc_mean_delay() {}

    /**
     * @param arv arrival envelope
     * @param srv service envelope
     * @return the bound on E[D] and the minimizing theta
     */
    public static SncResult snc_mean_delay(SncEnvelope arv, SncEnvelope srv) {
        return snc_mean_delay(arv, srv, 1e3);
    }

    /**
     * @param arv      arrival envelope
     * @param srv      service envelope
     * @param thetamax upper end of the theta search
     * @return the bound on E[D] and the minimizing theta
     */
    public static SncResult snc_mean_delay(final SncEnvelope arv, final SncEnvelope srv,
                                           double thetamax) {
        return Snc_thetaopt.snc_thetaopt(theta -> {
            double[] p = pair(arv, srv, theta);
            if (p == null) {
                return Double.POSITIVE_INFINITY;
            }
            double logK = theta * (p[0] + p[2]) - Math.log(1.0 - Math.exp(-theta * (p[3] - p[1])));
            double a = theta * p[3];
            return logK >= 0.0 ? (logK + 1.0) / a : Math.exp(logK) / a;
        }, thetamax);
    }

    /** Both envelopes at one theta, or null when the composition is infeasible. */
    private static double[] pair(SncEnvelope arv, SncEnvelope srv, double theta) {
        double[] a = arv.eval(theta);
        double[] s = srv.eval(theta);
        if (!Double.isFinite(a[0]) || !Double.isFinite(a[1])
                || !Double.isFinite(s[0]) || !Double.isFinite(s[1]) || s[1] <= a[1]) {
            return null;
        }
        return new double[] {a[0], a[1], s[0], s[1]};
    }
}
