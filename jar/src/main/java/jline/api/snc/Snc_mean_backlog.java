/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Upper bound on the mean backlog, from integrating the backlog tail bound.
 *
 * <p>The counterpart of {@link Snc_mean_delay} with {@code a = theta}: the
 * clipped integral of {@code K*exp(-theta*b)} is {@code (log(K)+1)/theta} when
 * K &gt;= 1 and {@code K/theta} otherwise. The unit of the answer is the unit of
 * the envelopes: jobs when the pair is {@link Snc_env_poisson} with
 * {@link Snc_srv_exp}, units of work when it is {@link Snc_env_cpoisson} with
 * {@link Snc_srv_rate}.</p>
 *
 * <p>SolverBA does NOT use this for its queue-length column: it applies
 * Little's law to the response-time bound instead, so that Q and R stay
 * consistent with the exact open-network throughput. The two are close but not
 * identical, since each optimizes its own theta.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_mean_backlog.m.</p>
 */
public final class Snc_mean_backlog {
    private Snc_mean_backlog() {}

    /**
     * @param arv arrival envelope
     * @param srv service envelope
     * @return the bound on E[B] and the minimizing theta
     */
    public static SncResult snc_mean_backlog(SncEnvelope arv, SncEnvelope srv) {
        return snc_mean_backlog(arv, srv, 1e3);
    }

    /**
     * @param arv      arrival envelope
     * @param srv      service envelope
     * @param thetamax upper end of the theta search
     * @return the bound on E[B] and the minimizing theta
     */
    public static SncResult snc_mean_backlog(final SncEnvelope arv, final SncEnvelope srv,
                                             double thetamax) {
        return Snc_thetaopt.snc_thetaopt(theta -> {
            double[] p = pair(arv, srv, theta);
            if (p == null) {
                return Double.POSITIVE_INFINITY;
            }
            double logK = theta * (p[0] + p[2]) - Math.log(1.0 - Math.exp(-theta * (p[3] - p[1])));
            return logK >= 0.0 ? (logK + 1.0) / theta : Math.exp(logK) / theta;
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
