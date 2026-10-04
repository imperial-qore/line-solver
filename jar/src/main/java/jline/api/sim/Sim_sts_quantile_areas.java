/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

import java.util.Arrays;

/**
 * Standardized time series areas of the batched quantile process.
 *
 * <p>Splits {@code b*m} observations into b nonoverlapping batches of size m and
 * returns the signed standardized time series areas of the quantile-estimation
 * process, the batched quantile estimators, and the three variance-parameter
 * estimators {@link Sim_fquest} and {@link Sim_firquest} build intervals from.
 *
 * <p>With {@code yhat_p(j,m)} the empirical p-quantile of batch j and
 * {@code yhat_p(j,k)} that of its first k observations, the STS process of batch
 * j is
 *
 * <pre>
 *   T_{j,m}(k/m) = (k/sqrt(m)) (yhat_p(j,m) - yhat_p(j,k)),
 * </pre>
 *
 * <p>its signed area is {@code A_p(w;j,m) = m^-1 sum_k w(k/m) T_{j,m}(k/m)}, and
 * the three estimators of {@code sigma_p^2 = lim n Var(ytilde_p(n))} are
 *
 * <pre>
 *   A_p(w;b,m) = b^-1 sum_j A_p(w;j,m)^2                        (STS area)
 *   N_p(b,m)   = (b-1)^-1 m sum_j (yhat_p(j,m)-ytilde_p(n))^2   (NBQ)
 *   V_p(w;b,m) = [b A_p(w;b,m) + (b-1) N_p(b,m)] / (2b-1)       (combined)
 * </pre>
 *
 * <p>with {@code ytilde_p(n)} the full-sample empirical p-quantile. The first two
 * have limiting chi-square laws on b and b-1 degrees of freedom and are
 * asymptotically independent, so the combined estimator carries 2b-1 degrees of
 * freedom and is about sqrt(2) less variable than either component.
 *
 * <p>The requirement on a weight function w is that {@code int_0^1 w(t)B(t)dt} be
 * standard normal for a standard Brownian bridge B; for a constant {@code w = c}
 * that variance is {@code c^2/12}, so {@code sqrt(12)} is the normalizing choice
 * and the default.
 *
 * <p>The prefix quantiles are exact order statistics, obtained from a Fenwick
 * tree over the within-batch ranks, so the cost is O(b m log m). On i.i.d. Exp(1)
 * data, where {@code sigma_p^2 = p(1-p)/f(y_p)^2} is exact, both A_p and N_p are
 * unbiased to within 7%.
 *
 * <p>Port of MATLAB sim_sts_quantile_areas.m.
 *
 * <p>Reference: C. Alexopoulos, D. Goldsman, A. Lolos, K. D. Dingec,
 * J. R. Wilson, "Steady-State Quantile Estimation Using Standardized Time
 * Series", 2020/2023; A. Lolos et al., Proc. Winter Simulation Conference, 2023,
 * theorems 1 to 3.
 *
 * @since LINE 3.1.0
 */
public final class Sim_sts_quantile_areas {
    /** Constant STS weight function that normalizes the Brownian bridge area. */
    public static final double DEFAULT_WEIGHT = Math.sqrt(12.0);

    private Sim_sts_quantile_areas() {}

    /**
     * Computes the STS statistics with the default weight function.
     *
     * @param y exactly b*m observations, in time order
     * @param b batch count, at least 1
     * @param m batch size, at least 1
     * @param p quantile probability in (0,1)
     * @return the statistics
     */
    public static StsQuantileStats sim_sts_quantile_areas(double[] y, int b, int m, double p) {
        return sim_sts_quantile_areas(y, b, m, p, DEFAULT_WEIGHT);
    }

    /**
     * Computes the STS statistics.
     *
     * @param y      exactly b*m observations, in time order
     * @param b      batch count, at least 1
     * @param m      batch size, at least 1
     * @param p      quantile probability in (0,1)
     * @param weight constant STS weight function, nonzero
     * @return the statistics
     */
    public static StsQuantileStats sim_sts_quantile_areas(double[] y, int b, int m,
                                                         double p, double weight) {
        if (y == null) {
            throw new IllegalArgumentException("The sample path must not be null");
        }
        if (b < 1) {
            throw new IllegalArgumentException("The batch count b=" + b + " must be positive");
        }
        if (m < 1) {
            throw new IllegalArgumentException("The batch size m=" + m + " must be positive");
        }
        if (!(p > 0.0) || !(p < 1.0)) {
            throw new IllegalArgumentException("p=" + p + " must lie in (0,1)");
        }
        if (weight == 0.0 || !Double.isFinite(weight)) {
            throw new IllegalArgumentException("weight=" + weight + " must be nonzero and finite");
        }
        long need = (long) b * m;
        if (y.length != need) {
            throw new IllegalArgumentException("Y must hold exactly b*m = " + need
                    + " observations, got " + y.length);
        }
        for (int i = 0; i < y.length; i++) {
            if (!Double.isFinite(y[i])) {
                throw new IllegalArgumentException("The sample path must be finite");
            }
        }

        int n = (int) need;
        double[] areas = new double[b];
        double[] bqe = new double[b];
        int batchIdx = (int) Math.ceil(m * p);

        double[] batch = new double[m];
        double[] sorted = new double[m];
        int[] rank = new int[m];
        int[] fen = new int[m + 1];
        Integer[] order = new Integer[m];

        for (int j = 0; j < b; j++) {
            System.arraycopy(y, j * m, batch, 0, m);
            System.arraycopy(batch, 0, sorted, 0, m);
            Arrays.sort(sorted);

            // rank[i] is the 1-based position of batch[i] among the sorted values,
            // ties broken by time index so the map is a permutation
            for (int i = 0; i < m; i++) {
                order[i] = i;
            }
            final double[] key = batch;
            Arrays.sort(order, (l, r) -> {
                int c = Double.compare(key[l], key[r]);
                return c != 0 ? c : Integer.compare(l, r);
            });
            for (int r = 0; r < m; r++) {
                rank[order[r]] = r + 1;
            }

            double full = sorted[batchIdx - 1];
            bqe[j] = full;

            Arrays.fill(fen, 0);
            int log = 31 - Integer.numberOfLeadingZeros(m);
            double acc = 0.0;
            for (int k = 1; k <= m; k++) {
                for (int pos = rank[k - 1]; pos <= m; pos += pos & (-pos)) {
                    fen[pos]++;
                }
                int target = (int) Math.ceil(p * k);
                int pos = 0;
                int rem = target;
                for (int step = 1 << log; step >= 1; step >>= 1) {
                    int cand = pos + step;
                    if (cand <= m && fen[cand] < rem) {
                        pos = cand;
                        rem -= fen[cand];
                    }
                }
                acc += k * (full - sorted[pos]);
            }
            areas[j] = weight * acc / (m * Math.sqrt(m));
        }

        double[] all = y.clone();
        Arrays.sort(all);
        double quantile = all[(int) Math.ceil(n * p) - 1];

        double ap = 0.0;
        for (int j = 0; j < b; j++) {
            ap += areas[j] * areas[j];
        }
        ap /= b;

        double np;
        double vp;
        if (b >= 2) {
            double acc = 0.0;
            for (int j = 0; j < b; j++) {
                double d = bqe[j] - quantile;
                acc += d * d;
            }
            np = m * acc / (b - 1);
            vp = (b * ap + (b - 1) * np) / (2.0 * b - 1);
        } else {
            // a single batch carries no between-batch degrees of freedom; areas and
            // bqe stay valid and Sim_firquest pools them across replications instead
            np = Double.NaN;
            vp = Double.NaN;
        }

        return new StsQuantileStats(areas, bqe, quantile, ap, np, vp, b, m, n);
    }
}
