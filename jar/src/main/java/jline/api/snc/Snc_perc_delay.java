/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Delay quantile at a prescribed violation probability.
 *
 * <p>Inverts {@link Snc_bound_delay} in d: at fixed theta,</p>
 *
 * <pre>
 *   d(theta) = (sigmaA + sigmaS
 *               - log(eps*(1-exp(-theta*(rhoS-rhoA))))/theta) / rhoS,
 * </pre>
 *
 * <p>minimized over the feasible thetas. This is the deliverable of the domain:
 * a statistical delay guarantee, the quantity a service-level objective is
 * written against, as opposed to the mean delay returned by the
 * queueing-theoretic solvers.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_perc_delay.m.</p>
 */
public final class Snc_perc_delay {
    private Snc_perc_delay() {}

    /**
     * @param arv arrival envelope
     * @param srv service envelope
     * @param eps violation probability, 0 &lt; eps &lt; 1
     * @return the quantile and the minimizing theta
     */
    public static SncResult snc_perc_delay(SncEnvelope arv, SncEnvelope srv, double eps) {
        return snc_perc_delay(arv, srv, eps, 1e3);
    }

    /**
     * @param arv      arrival envelope
     * @param srv      service envelope
     * @param eps      violation probability, 0 &lt; eps &lt; 1
     * @param thetamax upper end of the theta search
     * @return the quantile and the minimizing theta
     */
    public static SncResult snc_perc_delay(final SncEnvelope arv, final SncEnvelope srv,
                                           final double eps, double thetamax) {
        if (!(eps > 0.0 && eps < 1.0)) {
            throw new IllegalArgumentException("snc_perc_delay: eps must lie in (0,1), got " + eps);
        }
        SncResult r = Snc_thetaopt.snc_thetaopt(theta -> {
            double[] p = pair(arv, srv, theta);
            if (p == null) {
                return Double.POSITIVE_INFINITY;
            }
            return (p[0] + p[2]
                    - Math.log(eps * (1.0 - Math.exp(-theta * (p[3] - p[1])))) / theta) / p[3];
        }, thetamax);
        return new SncResult(Math.max(r.value, 0.0), r.theta);
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
