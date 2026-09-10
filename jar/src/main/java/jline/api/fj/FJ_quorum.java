/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fj;

import java.util.ArrayList;
import java.util.List;

/**
 * Completion time of a k-of-n (quorum) join over independent, not necessarily identically
 * distributed branches.
 *
 * <p>Each branch is summarised by its first two moments and expanded into a discrete step CDF
 * by a two-moment fit. The CDF of the k-th order statistic is then assembled by the
 * inclusion-exclusion identity</p>
 *
 * <pre>
 *   F_(k)(t) = sum_{i=k..n} (-1)^(i-k) * C(i-1, k-1) * e_i(F_1(t), ..., F_n(t))
 * </pre>
 *
 * <p>where e_i is the i-th elementary symmetric polynomial of the branch CDFs, evaluated by the
 * memoised recurrence in {@link #symmetric}. For k = n this collapses to the product of the branch
 * CDFs, i.e. the ordinary AND-join, and for k = 1 to 1 - prod(1 - F_i), i.e. the minimum.</p>
 *
 * <p>For identically distributed branches the result agrees with
 * {@link FJ_order_stat#fj_order_stat_cdf(double, int, int)}, which is the closed binomial form of
 * the same quantity; this class generalises it to heterogeneous branches, which is what an LQN
 * AND-join with a quorum requires since its branches carry different demands.</p>
 *
 * <p>Follows the formulation of Omari, Franks, Woodside and Pan, as implemented in LQNS 6.x
 * (randomvar.cc, DiscreteCDFs::quorumKofN and DiscretePoints::estimateCDF).</p>
 */
public class FJ_quorum {

    /**
     * A discrete step function on a finite, increasing time grid.
     *
     * <p>Holds parallel arrays t[0..m-1] of abscissae and A[0..m-1] of accumulated values, with
     * A implicitly 0 below t[0] and constant at A[m-1] above t[m-1]. When the values are a
     * probability distribution, A is its CDF; intermediate results of the inclusion-exclusion
     * sum are not distributions and may leave [0, 1].</p>
     */
    public static final class StepCDF {
        public final double[] t;
        public final double[] A;

        public StepCDF(double[] t, double[] A) {
            if (t.length != A.length) {
                throw new IllegalArgumentException("t and A must have equal length, got "
                        + t.length + " and " + A.length + ".");
            }
            this.t = t;
            this.A = A;
        }

        public boolean isEmpty() {
            return t.length == 0;
        }

        /** Value of the step function at time x. */
        public double at(double x) {
            double v = 0.0;
            for (int i = 0; i < t.length && t[i] <= x; i++) {
                v = A[i];
            }
            return v;
        }

        /** Scale the values, leaving the time grid untouched. */
        public StepCDF times(double c) {
            double[] B = new double[A.length];
            for (int i = 0; i < A.length; i++) {
                B[i] = A[i] * c;
            }
            return new StepCDF(t.clone(), B);
        }

        /** Mean of the distribution described by this CDF. */
        public double mean() {
            double m = 0.0;
            double prev = 0.0;
            for (int i = 0; i < t.length; i++) {
                m += (A[i] - prev) * t[i];
                prev = A[i];
            }
            return m;
        }

        /** Variance of the distribution described by this CDF, clamped at zero. */
        public double variance() {
            double m = mean();
            double v = 0.0;
            double prev = 0.0;
            for (int i = 0; i < t.length; i++) {
                double d = t[i] - m;
                v += (A[i] - prev) * d * d;
                prev = A[i];
            }
            return Math.max(v, 0.0);
        }
    }

    /**
     * Two-moment fit of a branch completion time to a three-point discrete distribution.
     *
     * <p>Matches the mean and variance exactly. A zero standard deviation degenerates to a single
     * deterministic point, and a zero mean yields the empty step function.</p>
     *
     * @param mean     branch mean completion time, must be non-negative
     * @param variance branch variance, must be non-negative
     * @return the fitted step CDF
     */
    public static StepCDF threePointFit(double mean, double variance) {
        if (mean < 0.0 || variance < 0.0) {
            throw new IllegalArgumentException("mean and variance must be non-negative, got mean="
                    + mean + ", variance=" + variance + ".");
        }
        if (mean == 0.0) {
            return new StepCDF(new double[0], new double[0]);
        }
        double sd = Math.sqrt(variance);
        if (sd == 0.0) {
            return new StepCDF(new double[]{mean}, new double[]{1.0});
        }

        double t1 = (mean > sd) ? mean - sd : 0.0;
        double t2 = mean;
        double t3 = (sd >= mean) ? mean + 2.0 * variance / mean : mean + 2.0 * sd;

        double delta = t1 * t1 * (t3 - t2) + t2 * t2 * (t1 - t3) + t3 * t3 * (t2 - t1);
        if (delta == 0.0) {
            // The three abscissae are not distinct, so the fit is not determined; fall back to
            // the deterministic point rather than dividing by zero.
            return new StepCDF(new double[]{mean}, new double[]{1.0});
        }
        double temp = variance + mean * mean;

        double a1 = (temp * (t3 - t2) + t2 * t2 * (mean - t3) + t3 * t3 * (t2 - mean)) / delta;
        double a3 = (t1 * t1 * (mean - t2) + t2 * t2 * (t1 - mean) + temp * (t2 - t1)) / delta;

        return new StepCDF(new double[]{t1, t2, t3}, new double[]{a1, 1.0 - a3, 1.0});
    }

    /**
     * Pointwise sum of two step functions on the union of their time grids.
     */
    public static StepCDF add(StepCDF x, StepCDF y) {
        return combine(x, y, true);
    }

    /**
     * Pointwise product of two step functions on the union of their time grids.
     *
     * <p>An empty operand annihilates the product, matching the convention that an absent branch
     * distribution contributes no mass.</p>
     */
    public static StepCDF multiply(StepCDF x, StepCDF y) {
        if (x.isEmpty() || y.isEmpty()) {
            return new StepCDF(new double[0], new double[0]);
        }
        return combine(x, y, false);
    }

    private static StepCDF combine(StepCDF x, StepCDF y, boolean isAdd) {
        int nx = x.t.length;
        int ny = y.t.length;
        List<Double> ts = new ArrayList<>(nx + ny);
        List<Double> as = new ArrayList<>(nx + ny);

        int i = 0;
        int j = 0;
        double xPrev = 0.0;
        double yPrev = 0.0;
        while (i < nx || j < ny) {
            double time;
            // Advance whichever grid holds the next time, stepping both when they coincide so a
            // shared abscissa is emitted once.
            if (j >= ny || (i < nx && x.t[i] < y.t[j])) {
                time = x.t[i];
                xPrev = x.A[i];
                i++;
            } else if (i >= nx || y.t[j] < x.t[i]) {
                time = y.t[j];
                yPrev = y.A[j];
                j++;
            } else {
                time = x.t[i];
                xPrev = x.A[i];
                yPrev = y.A[j];
                i++;
                j++;
            }
            ts.add(time);
            as.add(isAdd ? xPrev + yPrev : xPrev * yPrev);
        }

        double[] tOut = new double[ts.size()];
        double[] aOut = new double[as.size()];
        for (int k = 0; k < tOut.length; k++) {
            tOut[k] = ts.get(k);
            aOut[k] = as.get(k);
        }
        return new StepCDF(tOut, aOut);
    }

    /**
     * Largest branch count accepted by {@link #quorumKofN}.
     *
     * <p>The evaluation is cubic in the branch count, so this bounds the work at roughly 1e8
     * elementary operations. An LQN AND-join with more branches than this is degenerate; the
     * limit exists to turn a silent hang into an immediate, explicable failure.</p>
     */
    public static final int MAX_BRANCHES = 512;

    /**
     * The union of the branch time grids, sorted and de-duplicated.
     */
    private static double[] mergedGrid(StepCDF[] branches) {
        int total = 0;
        for (StepCDF b : branches) {
            total += b.t.length;
        }
        double[] all = new double[total];
        int p = 0;
        for (StepCDF b : branches) {
            for (double ti : b.t) {
                all[p++] = ti;
            }
        }
        java.util.Arrays.sort(all);
        int m = 0;
        for (int i = 0; i < all.length; i++) {
            if (i == 0 || all[i] != all[i - 1]) {
                all[m++] = all[i];
            }
        }
        return java.util.Arrays.copyOf(all, m);
    }

    /**
     * CDF of the k-th smallest of n independent branch completion times.
     *
     * <p>Evaluated pointwise on the union of the branch grids. At each time the number of
     * completed branches is Poisson-binomial, so its distribution is built by the recurrence
     * q_j &lt;- q_{j-1} * F_i + q_j * (1 - F_i) over branches i, and the k-th order statistic is
     * the upper tail sum_{j&gt;=k} q_j.</p>
     *
     * <p>This is the same quantity as the inclusion-exclusion identity
     * sum_{i=k..n} (-1)^(i-k) C(i-1,k-1) e_i(F_1..F_n) used by LQNS, but every term here is a
     * probability in [0, 1] and no term is subtracted, so it does not suffer the catastrophic
     * cancellation that the alternating binomial sum incurs as n grows.</p>
     *
     * @param branches the branch CDFs, one per branch
     * @param k        the quorum count, in [1, branches.length]
     * @return the step CDF of the k-th order statistic
     * @throws IllegalArgumentException if k is out of range or there are more than
     *                                  {@link #MAX_BRANCHES} branches
     */
    public static StepCDF quorumKofN(StepCDF[] branches, int k) {
        int n = branches.length;
        if (n == 0) {
            return new StepCDF(new double[0], new double[0]);
        }
        if (k < 1 || k > n) {
            throw new IllegalArgumentException("k must satisfy 1 <= k <= n. Got k=" + k + ", n=" + n + ".");
        }
        if (n > MAX_BRANCHES) {
            throw new IllegalArgumentException("quorum join has " + n + " branches, above the supported "
                    + "maximum of " + MAX_BRANCHES + "; evaluation is cubic in the branch count.");
        }

        // A branch with no mass completes instantly, so that it neither delays the join nor
        // suppresses the completion counts below.
        StepCDF[] b = new StepCDF[n];
        for (int i = 0; i < n; i++) {
            b[i] = branches[i].isEmpty() ? new StepCDF(new double[]{0.0}, new double[]{1.0}) : branches[i];
        }

        double[] grid = mergedGrid(b);
        int m = grid.length;

        // Tabulate each branch along the merged grid by a single monotone walk, so the
        // evaluation below is linear rather than quadratic in the grid size.
        double[][] fv = new double[n][m];
        for (int i = 0; i < n; i++) {
            StepCDF br = b[i];
            int p = 0;
            double cur = 0.0;
            for (int g = 0; g < m; g++) {
                while (p < br.t.length && br.t[p] <= grid[g]) {
                    cur = br.A[p];
                    p++;
                }
                fv[i][g] = cur;
            }
        }

        double[] out = new double[m];
        double[] q = new double[n + 1];
        for (int g = 0; g < m; g++) {
            java.util.Arrays.fill(q, 0.0);
            q[0] = 1.0;
            for (int i = 0; i < n; i++) {
                double f = fv[i][g];
                for (int j = Math.min(i + 1, n); j >= 1; j--) {
                    q[j] = q[j - 1] * f + q[j] * (1.0 - f);
                }
                q[0] = q[0] * (1.0 - f);
            }
            double tail = 0.0;
            for (int j = k; j <= n; j++) {
                tail += q[j];
            }
            out[g] = tail;
        }
        return new StepCDF(grid, out);
    }

    /**
     * Mean and variance of the completion time of a k-of-n join whose branches are given by
     * their first two moments.
     *
     * <p>This is the entry point used by the layered solver: each branch mean and variance comes
     * from the residence time and variance accumulated along that branch of the activity graph.</p>
     *
     * @param branchMeans     per-branch mean completion times
     * @param branchVariances per-branch variances, same length as branchMeans
     * @param k               the quorum count, in [1, branchMeans.length]
     * @return a two-element array holding the mean and the variance of the join completion time
     */
    public static double[] quorumMoments(double[] branchMeans, double[] branchVariances, int k) {
        if (branchMeans.length != branchVariances.length) {
            throw new IllegalArgumentException("branchMeans and branchVariances must have equal length, got "
                    + branchMeans.length + " and " + branchVariances.length + ".");
        }
        int n = branchMeans.length;
        if (n == 0) {
            return new double[]{0.0, 0.0};
        }
        StepCDF[] branches = new StepCDF[n];
        for (int i = 0; i < n; i++) {
            branches[i] = threePointFit(branchMeans[i], branchVariances[i]);
        }
        StepCDF join = quorumKofN(branches, k);
        return new double[]{join.mean(), join.variance()};
    }

}
