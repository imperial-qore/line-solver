/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Violation probability of a backlog level.
 *
 * <p>For a flow with arrival envelope (sigmaA,rhoA) served by an element with
 * service envelope (sigmaS,rhoS), the backlog of the stable station obeys, for
 * every theta &gt; 0,</p>
 *
 * <pre>
 *   P{B(t) &gt; b} &lt;= exp(-theta*(b-sigmaA-sigmaS)) / (1-exp(-theta*(rhoS-rhoA))),
 * </pre>
 *
 * <p>the union bound over the start of the backlogged period summed as a
 * geometric series on the unit-slot time axis. The returned value is the
 * infimum over theta, clipped at 1, and is an UPPER BOUND on the tail, never an
 * estimate of it: the decay rate is asymptotically exact and the prefactor is
 * loose.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_bound_backlog.m.</p>
 */
public final class Snc_bound_backlog {
    private Snc_bound_backlog() {}

    /**
     * @param arv arrival envelope
     * @param srv service envelope
     * @param b   backlog level, units of the envelopes
     * @return the bound on P{B &gt; b} in [0,1], and the minimizing theta
     */
    public static SncResult snc_bound_backlog(SncEnvelope arv, SncEnvelope srv, double b) {
        return snc_bound_backlog(arv, srv, b, 1e3);
    }

    /**
     * @param arv      arrival envelope
     * @param srv      service envelope
     * @param b        backlog level
     * @param thetamax upper end of the theta search
     * @return the bound on P{B &gt; b} in [0,1], and the minimizing theta
     */
    public static SncResult snc_bound_backlog(final SncEnvelope arv, final SncEnvelope srv,
                                              final double b, double thetamax) {
        if (b < 0) {
            throw new IllegalArgumentException("snc_bound_backlog: b must be nonnegative, got " + b);
        }
        SncResult r = Snc_thetaopt.snc_thetaopt(theta -> {
            double[] p = pair(arv, srv, theta);
            if (p == null) {
                return Double.POSITIVE_INFINITY;
            }
            return Math.exp(-theta * (b - p[0] - p[2])) / (1.0 - Math.exp(-theta * (p[3] - p[1])));
        }, thetamax);
        // no feasible theta, or the bound is vacuous at this level
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
