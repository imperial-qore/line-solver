/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.dqsys;

/**
 * Exact analysis of a state dependent Bernoulli server on a discrete time
 * scale.
 *
 * <p>Port of matlab/src/api/dqsys/dqsys_bernoulli1.m. Time advances in slots.
 * In the slot starting at t with n jobs present the job in service departs with
 * probability p(n) and an arrival occurs with probability b(n), independently;
 * both are recorded at the end of the slot with the departure resolved first
 * (Daduna's LA rule and D/A rule). The queue length at slot boundaries is a
 * discrete birth-death chain with</p>
 *
 * <pre>
 *   pi(n) = [prod_{m=0}^{n-1} b(m) / prod_{m=0}^{n} c(m)]
 *         * [prod_{m=1}^{n-1} q(m) / prod_{m=1}^{n} p(m)] / H
 * </pre>
 *
 * <p>c = 1-b and q = 1-p, which is theorem 2.3, and corollary 2.8 once b(n) = 0
 * above the capacity. For constant b and p it collapses to the Geo/Geo/1 law of
 * {@link jline.api.dqsys.Dqsys_geogeo1} under the LAS_DA convention.</p>
 *
 * <p>The law seen by an arriving customer, with himself not counted, is theorem
 * 2.11 and is returned in {@code arrivalPmf}. It is not the time-stationary
 * law: discrete time has no PASTA analogue, and the two differ even when the
 * arrival stream is a state independent Bernoulli process. In that state
 * independent case pi_1 is exactly the EAS-convention queue length law of
 * {@link jline.api.dqsys.Dqsys_geogeo1}, geometric with ratio
 * r = b(1-p)/(p(1-b)).</p>
 */
public final class Dqsys_bernoulli1 {

    private Dqsys_bernoulli1() {}

    /** Unbounded buffer with constant arrival and service probabilities. */
    public static Bernoulli1Result dqsys_bernoulli1(double b, double p) {
        if (b <= 0 || b > 1) {
            throw new IllegalArgumentException("arrival probabilities must be real and in [0,1]");
        }
        if (p <= 0 || p > 1) {
            throw new IllegalArgumentException("service probabilities must be real and in (0,1]");
        }
        if (b >= p) {
            throw new IllegalArgumentException(
                    "load b/p must be strictly less than 1 on an unbounded buffer");
        }
        GeoGeo1Result g = Dqsys_geogeo1.dqsys_geogeo1(b, p, GeoGeo1Convention.LAS_DA);
        return new Bernoulli1Result(Integer.MAX_VALUE, new double[]{b}, new double[]{p},
                null, null, g.getEmptyProb(), g.getUtilization(), g.getThroughput(), 0.0,
                g.getMeanQueueLength(), g.getMeanWaitingQueue(), g.getMeanSojournTime(),
                g.getMeanWaitingTime(), 1.0 / g.getEmptyProb());
    }

    /**
     * Finite buffer, state dependent probabilities. An arrival in a slot that
     * finds L jobs present is lost, which is the loss system of corollary 2.8.
     *
     * @param b offered arrival probabilities, {@code b[n]} for n = 0..L, or a
     *          single-entry array for a state independent stream
     * @param p service probabilities, {@code p[n-1] = p(n)} for n = 1..L, or a
     *          single-entry array for a state independent server
     * @param L buffer capacity in jobs
     */
    public static Bernoulli1Result dqsys_bernoulli1(double[] b, double[] p, int L) {
        if (L < 1) {
            throw new IllegalArgumentException("L must be a positive integer");
        }
        double[] boff = expandState(b, L + 1, "arrival");
        double[] pv = expandService(p, L);
        for (int i = 0; i <= L; i++) {
            if (boff[i] < 0 || boff[i] > 1) {
                throw new IllegalArgumentException("arrival probabilities must be real and in [0,1]");
            }
        }
        for (int i = 0; i < L; i++) {
            if (pv[i] <= 0 || pv[i] > 1) {
                throw new IllegalArgumentException("service probabilities must be real and in (0,1]");
            }
        }
        double[] badm = boff.clone();
        badm[L] = 0.0;                        // an arrival finding L jobs is lost
        for (int i = 0; i < L; i++) {
            if (badm[i] >= 1) {
                // c(n)=0 makes the weight of theorem 2.3 diverge at n; p(n)=1 is
                // fine and truncates the chain instead, which is example 2.9.
                throw new IllegalArgumentException(
                        "arrival probabilities below the capacity must be strictly less than one");
            }
        }

        double[] lw = new double[L + 1];      // log unnormalized pi, theorem 2.3
        double acc = -Math.log(1.0 - badm[0]);
        lw[0] = acc;
        for (int n = 1; n <= L; n++) {
            acc += log(badm[n - 1]) - Math.log(1.0 - badm[n]) - Math.log(pv[n - 1]);
            if (n >= 2) {
                acc += log(1.0 - pv[n - 2]);
            }
            lw[n] = acc;
        }
        double lwMax = Double.NEGATIVE_INFINITY;
        for (int n = 0; n <= L; n++) {
            if (lw[n] > lwMax) {
                lwMax = lw[n];
            }
        }
        double[] w = new double[L + 1];
        double H = 0.0;
        for (int n = 0; n <= L; n++) {
            w[n] = Math.exp(lw[n] - lwMax);
            H += w[n];
        }
        double[] pmf = new double[L + 1];
        for (int n = 0; n <= L; n++) {
            pmf[n] = w[n] / H;
        }

        double[] la = new double[L];          // log unnormalized pi_1, theorem 2.11
        acc = log(badm[0]) - Math.log(1.0 - badm[0]) - Math.log(1.0 - badm[1]);
        la[0] = acc;
        for (int n = 1; n < L; n++) {
            acc += log(badm[n]) - Math.log(1.0 - badm[n + 1])
                    + log(1.0 - pv[n - 1]) - Math.log(pv[n - 1]);
            la[n] = acc;
        }
        double laMax = Double.NEGATIVE_INFINITY;
        for (int n = 0; n < L; n++) {
            if (!Double.isInfinite(la[n]) && la[n] > laMax) {
                laMax = la[n];
            }
        }
        double[] arrivalPmf = new double[L];
        if (!Double.isInfinite(laMax)) {
            double sa = 0.0;
            for (int n = 0; n < L; n++) {
                arrivalPmf[n] = Double.isInfinite(la[n]) ? 0.0 : Math.exp(la[n] - laMax);
                sa += arrivalPmf[n];
            }
            for (int n = 0; n < L; n++) {
                arrivalPmf[n] /= sa;
            }
        }

        double meanQueueLength = 0.0;
        double throughput = 0.0;
        double offered = 0.0;
        for (int n = 0; n <= L; n++) {
            meanQueueLength += pmf[n] * n;
            offered += pmf[n] * boff[n];
            if (n >= 1) {
                throughput += pmf[n] * pv[n - 1];
            }
        }
        double utilization = 1.0 - pmf[0];
        double lossProb = offered > 0 ? pmf[L] * boff[L] / offered : 0.0;
        double meanWaitingQueue = meanQueueLength - utilization;
        double sojourn = throughput > 0 ? meanQueueLength / throughput : 0.0;
        double waiting = throughput > 0 ? meanWaitingQueue / throughput : 0.0;

        return new Bernoulli1Result(L, boff, pv, pmf, arrivalPmf, pmf[0], utilization,
                throughput, lossProb, meanQueueLength, meanWaitingQueue, sojourn,
                waiting, H * Math.exp(lwMax));
    }

    private static double log(double x) {
        return x > 0 ? Math.log(x) : Double.NEGATIVE_INFINITY;
    }

    private static double[] expandState(double[] x, int len, String what) {
        if (x == null || x.length == 0) {
            throw new IllegalArgumentException("the " + what + " probability vector must not be empty");
        }
        if (x.length == 1) {
            double[] v = new double[len];
            for (int i = 0; i < len; i++) {
                v[i] = x[0];
            }
            return v;
        }
        if (x.length == len) {
            return x.clone();
        }
        throw new IllegalArgumentException(String.format(
                "the %s probability vector must have %d entries, one per state 0..%d",
                what, len, len - 1));
    }

    private static double[] expandService(double[] x, int L) {
        if (x == null || x.length == 0) {
            throw new IllegalArgumentException("the service probability vector must not be empty");
        }
        if (x.length == 1) {
            double[] v = new double[L];
            for (int i = 0; i < L; i++) {
                v[i] = x[0];
            }
            return v;
        }
        if (x.length == L) {
            return x.clone();
        }
        if (x.length == L + 1) {
            // A vector of length L+1 is accepted with its first entry, which
            // would be p(0), ignored.
            double[] v = new double[L];
            System.arraycopy(x, 1, v, 0, L);
            return v;
        }
        throw new IllegalArgumentException(String.format(
                "the service probability vector must have %d or %d entries", L, L + 1));
    }
}
