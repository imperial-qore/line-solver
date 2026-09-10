/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Backlog quantile at a prescribed violation probability.
 *
 * <p>Inverts {@link Snc_bound_backlog} in b: at fixed theta the smallest level
 * for which the bound certifies P{B &gt; b} &lt;= eps is</p>
 *
 * <pre>
 *   b(theta) = sigmaA + sigmaS - log(eps*(1-exp(-theta*(rhoS-rhoA))))/theta,
 * </pre>
 *
 * <p>and the reported quantile is its minimum over the feasible thetas. The
 * minimizing theta differs from the one of the forward bound at a given level,
 * which is why the inversion is done in closed form and re-optimized rather
 * than by a search on the forward bound.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_perc_backlog.m.</p>
 */
public final class Snc_perc_backlog {
    private Snc_perc_backlog() {}

    /**
     * @param arv arrival envelope
     * @param srv service envelope
     * @param eps violation probability, 0 &lt; eps &lt; 1
     * @return the quantile and the minimizing theta
     */
    public static SncResult snc_perc_backlog(SncEnvelope arv, SncEnvelope srv, double eps) {
        return snc_perc_backlog(arv, srv, eps, 1e3);
    }

    /**
     * @param arv      arrival envelope
     * @param srv      service envelope
     * @param eps      violation probability, 0 &lt; eps &lt; 1
     * @param thetamax upper end of the theta search
     * @return the quantile and the minimizing theta
     */
    public static SncResult snc_perc_backlog(final SncEnvelope arv, final SncEnvelope srv,
                                             final double eps, double thetamax) {
        if (!(eps > 0.0 && eps < 1.0)) {
            throw new IllegalArgumentException("snc_perc_backlog: eps must lie in (0,1), got " + eps);
        }
        SncResult r = Snc_thetaopt.snc_thetaopt(theta -> {
            double[] p = pair(arv, srv, theta);
            if (p == null) {
                return Double.POSITIVE_INFINITY;
            }
            return p[0] + p[2] - Math.log(eps * (1.0 - Math.exp(-theta * (p[3] - p[1])))) / theta;
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
