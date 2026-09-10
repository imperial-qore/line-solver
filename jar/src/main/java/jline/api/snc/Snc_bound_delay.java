/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Violation probability of a delay target.
 *
 * <p>For a flow with arrival envelope (sigmaA,rhoA) served by an element with
 * service envelope (sigmaS,rhoS), the virtual delay of the stable station
 * obeys, for every theta &gt; 0,</p>
 *
 * <pre>
 *   P{D(t) &gt; d} &lt;= exp(-theta*(rhoS*d-sigmaA-sigmaS)) / (1-exp(-theta*(rhoS-rhoA))),
 * </pre>
 *
 * <p>the horizontal rather than vertical deviation between the arrival and
 * service envelopes. On the M/M/1 read in job units ({@link Snc_env_poisson}
 * with {@link Snc_srv_exp}) the optimal theta tends to log(mu/lambda), so the
 * bound reproduces the exact asymptotic decay rate exp(-(mu-lambda)*d).</p>
 *
 * <p>Port of matlab/src/api/snc/snc_bound_delay.m.</p>
 */
public final class Snc_bound_delay {
    private Snc_bound_delay() {}

    /**
     * @param arv arrival envelope
     * @param srv service envelope
     * @param d   delay target, slots
     * @return the bound on P{D &gt; d} in [0,1], and the minimizing theta
     */
    public static SncResult snc_bound_delay(SncEnvelope arv, SncEnvelope srv, double d) {
        return snc_bound_delay(arv, srv, d, 1e3);
    }

    /**
     * @param arv      arrival envelope
     * @param srv      service envelope
     * @param d        delay target
     * @param thetamax upper end of the theta search
     * @return the bound on P{D &gt; d} in [0,1], and the minimizing theta
     */
    public static SncResult snc_bound_delay(final SncEnvelope arv, final SncEnvelope srv,
                                            final double d, double thetamax) {
        if (d < 0) {
            throw new IllegalArgumentException("snc_bound_delay: d must be nonnegative, got " + d);
        }
        SncResult r = Snc_thetaopt.snc_thetaopt(theta -> {
            double[] p = pair(arv, srv, theta);
            if (p == null) {
                return Double.POSITIVE_INFINITY;
            }
            return Math.exp(-theta * (p[3] * d - p[0] - p[2]))
                    / (1.0 - Math.exp(-theta * (p[3] - p[1])));
        }, thetamax);
        double eps = (!Double.isFinite(r.value) || r.value > 1.0) ? 1.0 : r.value;
        return new SncResult(eps, r.theta);
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
