/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Fixed-sample-size quantile interval from independent replications, FIRQUEST.
 *
 * <p>The replicated counterpart of {@link Sim_fquest}. It differs in four places.
 *
 * <ul>
 * <li>The warmup randomness test runs independently on each replicate path, and
 * the batch size it settles on may differ between replications.</li>
 * <li>Truncation removes the largest of those batch sizes from the front of every
 * replication, not just from one path. This is more aggressive on purpose: an
 * untruncated transient common to all replications biases every replicate
 * estimate the same way, and averaging cannot remove it.</li>
 * <li>The four stage tests act on the {@code R*b} signed areas and {@code R*b}
 * replicate batched quantile estimators pooled across replications, with the same
 * b in every replication and at least one batch from each.</li>
 * <li>The delivered interval is
 * {@code ytilde_p(N*) +- t_{1-alpha/2,2Rb-1} sqrt(Vtilde_p(w;R,b,m)/N*)},
 * {@code N* = R*b*m}, with the pooled estimators
 * <pre>
 *   A_p(w;R,b,m)    = (Rb)^-1 sum_j A_p(w;j,m)^2
 *   Ntilde_p(R,b,m) = m (Rb-1)^-1 sum_j (yhat_p(j,m) - ytilde_p(N*))^2
 *   Vtilde_p        = [Rb A_p + (Rb-1) Ntilde_p] / (2Rb-1).
 * </pre>
 * The heuristic fallback drops the residual-autocorrelation correction, since the
 * pooled batch quantiles come from independent paths.</li>
 * </ul>
 *
 * <p>Independent replications shorten the correlation the estimator has to fight,
 * and they parallelize, but they reintroduce initialization bias in every path, so
 * a short run length per replication is worse here than in {@link Sim_fquest}. The
 * article reports slight undercoverage at p = 0.99 when the total sample is under
 * 500000, down to 90.8%.
 *
 * <p>Port of MATLAB sim_firquest.m.
 *
 * <p>Reference: A. Lolos, C. Alexopoulos, D. Goldsman, K. D. Dingec,
 * A. C. Mokashi, J. R. Wilson, "A Fixed-Sample-Size Procedure for Estimating
 * Steady-State Quantiles Based on Independent Replications", Proc. Winter
 * Simulation Conference, 2025.
 *
 * @since LINE 3.1.0
 */
public final class Sim_firquest {
    private Sim_firquest() {}

    /**
     * Article default batch counts per replication, as a function of the
     * replication count.
     *
     * <p>Chosen so that {@code R*b} pooled statistics remain enough to test while
     * every replication still contributes at least one batch.
     *
     * @param R number of replications, at least 2
     * @return the descending batch counts
     */
    public static int[] defaultBatchCounts(int R) {
        if (R < 2) {
            throw new IllegalArgumentException("At least 2 replications are required, got " + R);
        }
        if (R == 2) {
            return new int[] {14, 11, 8, 5};
        }
        if (R == 3) {
            return new int[] {10, 8, 6, 4};
        }
        if (R == 4) {
            return new int[] {6, 5, 4, 3};
        }
        if (R < 10) {
            return new int[] {5, 4, 3, 2};
        }
        if (R < 17) {
            return new int[] {4, 3, 2, 1};
        }
        if (R < 23) {
            return new int[] {3, 2, 1};
        }
        if (R < 33) {
            return new int[] {2, 1};
        }
        return new int[] {1};
    }

    /**
     * Runs FIRQUEST at nominal 95% coverage with the default constants.
     *
     * @param y the replicate sample paths, {@code y[r]} being replication r, all
     *          of equal length
     * @param p quantile probability in (0,1)
     * @return the point estimate and interval
     */
    public static QuantileCIResult sim_firquest(double[][] y, double p) {
        return sim_firquest(y, p, 0.05, null);
    }

    /**
     * Runs FIRQUEST with the default constants.
     *
     * @param y     the replicate sample paths, all of equal length
     * @param p     quantile probability in (0,1)
     * @param alpha significance level in (0,1)
     * @return the point estimate and interval
     */
    public static QuantileCIResult sim_firquest(double[][] y, double p, double alpha) {
        return sim_firquest(y, p, alpha, null);
    }

    /**
     * Runs FIRQUEST.
     *
     * @param y       the replicate sample paths, all of equal length
     * @param p       quantile probability in (0,1)
     * @param alpha   significance level in (0,1)
     * @param options procedure constants, null for the defaults; when b0 or s are
     *                left at the {@link Sim_fquest} defaults they are replaced by
     *                25 and {@link #defaultBatchCounts(int)}
     * @return the point estimate and interval
     */
    public static QuantileCIResult sim_firquest(double[][] y, double p, double alpha,
                                               QuestOptions options) {
        if (y == null || y.length < 2) {
            throw new IllegalArgumentException(
                    "At least 2 replications are required, use Sim_fquest for a single path");
        }
        int R = y.length;
        int n = y[0] == null ? 0 : y[0].length;
        if (n < 1) {
            throw new IllegalArgumentException("The replicate paths must be nonempty");
        }
        for (int r = 0; r < R; r++) {
            if (y[r] == null || y[r].length != n) {
                throw new IllegalArgumentException(
                        "The replicate paths must all have length " + n);
            }
            for (int i = 0; i < n; i++) {
                if (!Double.isFinite(y[r][i])) {
                    throw new IllegalArgumentException("The sample paths must be finite");
                }
            }
        }
        if (!(p > 0.0) || !(p < 1.0)) {
            throw new IllegalArgumentException("p=" + p + " must lie in (0,1)");
        }
        if (!(alpha > 0.0) || !(alpha < 1.0)) {
            throw new IllegalArgumentException("alpha=" + alpha + " must lie in (0,1)");
        }

        QuestOptions opt = options;
        if (opt == null) {
            opt = new QuestOptions();
            opt.b0 = 25;
            opt.s = defaultBatchCounts(R);
        } else {
            QuestOptions dflt = new QuestOptions();
            if (opt.b0 == dflt.b0) {
                opt.b0 = 25;
            }
            if (Arrays.equals(opt.s, dflt.s)) {
                opt.s = defaultBatchCounts(R);
            }
        }
        opt.validate();
        int smallest = opt.s[opt.s.length - 1];
        if (R * smallest < 3) {
            throw new IllegalArgumentException("R*min(s) = " + (R * smallest)
                    + " pooled batches is below the 3 the stage tests need");
        }

        List<String> warnings = new ArrayList<String>();

        // warmup: one randomness loop per replicate path
        int b = opt.b0;
        int mStart = opt.m0;
        if (n < (long) b * mStart) {
            mStart = n / b;
        }
        if (mStart < 1) {
            throw new IllegalArgumentException("Each replication holds " + n
                    + " observations, too few for b0 = " + b + " batches");
        }
        int mMax = 0;
        boolean failed = false;
        for (int r = 0; r < R; r++) {
            int m = mStart;
            int ell = 1;
            boolean atMax = false;
            boolean passed = false;
            while (true) {
                StsQuantileStats stats = Sim_sts_quantile_areas.sim_sts_quantile_areas(
                        Arrays.copyOfRange(y[r], 0, b * m), b, m, p, opt.weight);
                double sig = opt.beta * Math.exp(-opt.eta * Math.pow(ell - 1, opt.theta));
                if (!Sim_vonneumann.sim_vonneumann(stats.areas, sig).reject) {
                    passed = true;
                    break;
                }
                if (atMax) {
                    break;
                }
                ell++;
                int mNext = (int) Math.round(m * Math.sqrt(2.0));
                if (n < (long) b * mNext && mNext != n / b) {
                    m = n / b;
                } else {
                    m = mNext;
                    if (n < (long) b * m) {
                        m = n / b;
                        atMax = true;
                    }
                }
                if (m < 1) {
                    break;
                }
            }
            failed = failed || !passed;
            mMax = Math.max(mMax, m);
        }
        if (failed) {
            warnings.add("the warmup randomness test could not be passed in every "
                    + "replication, the replicate paths are too short");
        }

        // truncation: delete the longest warmup batch from every replication
        int truncated = mMax;
        if (truncated >= n) {
            throw new IllegalArgumentException("The warmup batch size " + truncated
                    + " exhausts the replication length " + n);
        }
        double[][] yt = new double[R][];
        for (int r = 0; r < R; r++) {
            yt[r] = Arrays.copyOfRange(y[r], truncated, n);
        }
        int nStarRep = yt[0].length;

        // batch-count selection on the pooled statistics
        int v = 0;
        b = opt.s[v];
        int m = nStarRep / b;
        boolean ok = true;
        PooledStats pooled = null;
        for (int stage = 1; stage <= 4 && ok; stage++) {
            while (true) {
                if (m < 1) {
                    ok = false;
                    break;
                }
                pooled = pool(yt, b, m, p, opt.weight);
                double[] sample = stage <= 2 ? pooled.areas : pooled.bqe;
                boolean reject = stage % 2 == 1
                        ? Sim_vonneumann.sim_vonneumann(sample, opt.beta).reject
                        : Sim_shapirowilk.sim_shapirowilk(sample, opt.beta).reject;
                if (!reject) {
                    break;
                }
                v++;
                if (v >= opt.s.length) {
                    ok = false;
                    break;
                }
                b = opt.s[v];
                m = nStarRep / b;
            }
        }
        if (pooled == null || m < 1) {
            throw new IllegalArgumentException("The replicate paths are too short to form "
                    + smallest + " batches each");
        }

        int nStar = pooled.n;
        double estimate = pooled.quantile;
        double lower;
        double upper;
        double half;
        boolean heuristic;
        if (ok) {
            half = SimDist.tinv(1.0 - alpha / 2.0, 2.0 * R * b - 1)
                    * Math.sqrt(pooled.Vp / nStar);
            lower = estimate - half;
            upper = estimate + half;
            heuristic = false;
        } else {
            warnings.add("a randomness or normality test failed at b = " + smallest
                    + " per replication, the delivered interval is heuristic");
            heuristic = true;
            if (opt.force) {
                double[] ci = SimQuestHeuristic.interval(pooled.bqe, estimate, pooled.Ap,
                        pooled.Np, nStar, alpha, false);
                lower = ci[0];
                upper = ci[1];
                half = (upper - lower) / 2.0;
            } else {
                lower = Double.NaN;
                upper = Double.NaN;
                half = Double.NaN;
            }
        }

        return new QuantileCIResult(estimate, lower, upper, half, b, m, nStar, R, truncated,
                pooled.Ap, pooled.Np, pooled.Vp, heuristic, warnings);
    }

    private static PooledStats pool(double[][] yt, int b, int m, double p, double weight) {
        int R = yt.length;
        int nStarRep = yt[0].length;
        int K = R * b;
        double[] areas = new double[K];
        double[] bqe = new double[K];
        double[] all = new double[K * m];
        for (int r = 0; r < R; r++) {
            double[] seg = Arrays.copyOfRange(yt[r], nStarRep - b * m, nStarRep);
            StsQuantileStats stats = Sim_sts_quantile_areas.sim_sts_quantile_areas(
                    seg, b, m, p, weight);
            System.arraycopy(stats.areas, 0, areas, r * b, b);
            System.arraycopy(stats.bqe, 0, bqe, r * b, b);
            System.arraycopy(seg, 0, all, r * b * m, b * m);
        }

        int n = K * m;
        Arrays.sort(all);
        double quantile = all[(int) Math.ceil(n * p) - 1];

        double ap = 0.0;
        for (int j = 0; j < K; j++) {
            ap += areas[j] * areas[j];
        }
        ap /= K;

        double acc = 0.0;
        for (int j = 0; j < K; j++) {
            double d = bqe[j] - quantile;
            acc += d * d;
        }
        double np = m * acc / (K - 1);
        double vp = (K * ap + (K - 1) * np) / (2.0 * K - 1);

        return new PooledStats(areas, bqe, quantile, ap, np, vp, n);
    }

    private static final class PooledStats {
        final double[] areas;
        final double[] bqe;
        final double quantile;
        final double Ap;
        final double Np;
        final double Vp;
        final int n;

        PooledStats(double[] areas, double[] bqe, double quantile,
                    double Ap, double Np, double Vp, int n) {
            this.areas = areas;
            this.bqe = bqe;
            this.quantile = quantile;
            this.Ap = Ap;
            this.Np = Np;
            this.Vp = Vp;
            this.n = n;
        }
    }
}
