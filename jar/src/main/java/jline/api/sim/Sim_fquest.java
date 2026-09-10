/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Fixed-sample-size confidence interval for a steady-state quantile, FQUEST.
 *
 * <p>Takes a single simulation sample path of arbitrary fixed length and returns
 * a point estimate and an interval for the p-quantile of the steady-state
 * marginal distribution. No sequential control of the run length is needed, which
 * is the point: a dataset already in hand can be analyzed.
 *
 * <p>The procedure has four blocks.
 *
 * <ul>
 * <li><b>Warmup.</b> Starting from {@code b = b0} and {@code m = m0} it computes
 * the b signed STS areas of the batched quantile process and tests them for
 * randomness with von Neumann's ratio at the decaying significance
 * {@code beta*exp(-eta*(l-1)^theta)} on iteration l, growing m by sqrt(2)
 * whenever the test rejects. Passing means the areas are approximately
 * independent, so any initialization bias is confined to the first batch.</li>
 * <li><b>Truncation.</b> The first batch is deleted. That is the entire warmup
 * treatment; there is no separate transient detector.</li>
 * <li><b>Batch-count selection.</b> With b stepping down through {@code s} and
 * {@code m = nStar/b}, four tests must pass in order: von Neumann and Shapiro-Wilk
 * on the signed areas, then von Neumann and Shapiro-Wilk on the batched quantile
 * estimators. These check the asymptotic properties the interval rests on, that
 * both the areas and the batch quantiles behave like independent normal variates.
 * b only ever decreases and a failure at the smallest entry of s ends the
 * stage.</li>
 * <li><b>Delivery.</b> When all four pass the interval is
 * {@code ytilde_p(n*) +- t_{1-alpha/2,2b-1} sqrt(V_p(w;b,m)/n*)}. Otherwise the
 * sample was too small, {@link QuantileCIResult#heuristic} is true, and with
 * {@link QuestOptions#force} the interval is the conservative fallback of
 * {@link SimQuestHeuristic}.</li>
 * </ul>
 *
 * <p>Coverage was measured on the article's own test bed, the waiting-time
 * process of an M/M/1 queue with lambda = 0.8, mu = 1 started with 113 jobs in
 * system, over 500 independent replications at N = 200000, giving a standard
 * error near 1%: 95.2% at p = 0.5, 96.2% at p = 0.9 and 95.6% at p = 0.99
 * against a nominal 95%. The delivered half-width exceeds the empirically needed
 * one by factors of 1.16, 1.24 and 1.99 respectively, so the interval is
 * conservative and increasingly so into the tail, consistent with the half-widths
 * the article reports. On i.i.d. Exp(1) data, where
 * {@code sigma_p^2 = p(1-p)/f(y_p)^2} is exact, both A_p and N_p are unbiased to
 * within 7%.
 *
 * <p>Two properties are worth knowing. Coverage is not monotone in the sample
 * size over this range: a larger sample passes the stage tests more often and so
 * reaches the conservative fallback less. And a substantial fraction of runs
 * takes that fallback at all, 22% to 65% here and rising with p, so a delivered
 * interval may well be the heuristic one; {@link QuantileCIResult#heuristic} says
 * which. Fewer than about 100 replications cannot resolve a two-point difference
 * in coverage, so do not read a small experiment as a defect.
 *
 * <p>Applicability is a condition on the output process, not on the model that
 * produced it. The theory needs geometric moment contraction (Wu 2005), which
 * holds for ARMA series, a broad class of short-range-dependent linear and
 * nonlinear processes, many Markov chains, and was proved for M/M/1 and
 * non-heavy-tailed G/G/1 waiting times by Dingec et al. (2022); a density that is
 * positive and differentiable at the quantile of interest; short-range dependence
 * and an FCLT for the indicator process. M/M/1 is only the validation bed, chosen
 * because its exact quantiles are known.
 *
 * <p>Two practical exclusions follow. <b>Do not use this on integer-valued output
 * such as a queue length</b>: the marginal has no density, the density-regularity
 * condition fails, and the batched quantile has no Bahadur representation. Use it
 * on continuous output, that is response, waiting and sojourn times. And
 * heavy-tailed service, which can break geometric moment contraction and induce
 * long-range dependence, is outside the theory.
 *
 * <p>Port of MATLAB sim_fquest.m.
 *
 * <p>Reference: A. Lolos, C. Alexopoulos, D. Goldsman, K. D. Dingec,
 * A. C. Mokashi, J. R. Wilson, "A Fixed-Sample-Size Method for Estimating
 * Steady-State Quantiles", Proc. Winter Simulation Conference, 2023.
 *
 * @since LINE 3.1.0
 */
public final class Sim_fquest {
    private Sim_fquest() {}

    /**
     * Runs FQUEST at nominal 95% coverage with the default constants.
     *
     * @param y the sample path
     * @param p quantile probability in (0,1)
     * @return the point estimate and interval
     */
    public static QuantileCIResult sim_fquest(double[] y, double p) {
        return sim_fquest(y, p, 0.05, new QuestOptions());
    }

    /**
     * Runs FQUEST with the default constants.
     *
     * @param y     the sample path
     * @param p     quantile probability in (0,1)
     * @param alpha significance level in (0,1)
     * @return the point estimate and interval
     */
    public static QuantileCIResult sim_fquest(double[] y, double p, double alpha) {
        return sim_fquest(y, p, alpha, new QuestOptions());
    }

    /**
     * Runs FQUEST.
     *
     * @param y       the sample path
     * @param p       quantile probability in (0,1)
     * @param alpha   significance level in (0,1)
     * @param options procedure constants
     * @return the point estimate and interval
     */
    public static QuantileCIResult sim_fquest(double[] y, double p, double alpha,
                                             QuestOptions options) {
        if (y == null) {
            throw new IllegalArgumentException("The sample path must not be null");
        }
        if (!(p > 0.0) || !(p < 1.0)) {
            throw new IllegalArgumentException("p=" + p + " must lie in (0,1)");
        }
        if (!(alpha > 0.0) || !(alpha < 1.0)) {
            throw new IllegalArgumentException("alpha=" + alpha + " must lie in (0,1)");
        }
        QuestOptions opt = options == null ? new QuestOptions() : options;
        opt.validate();
        int smallest = opt.s[opt.s.length - 1];
        if (smallest < 3) {
            throw new IllegalArgumentException("min(s)=" + smallest
                    + ", but the stage tests need at least 3 batches");
        }
        int N = y.length;
        if (N < smallest * 2) {
            throw new IllegalArgumentException("The sample path holds " + N
                    + " observations, at least " + (smallest * 2) + " are needed");
        }
        for (int i = 0; i < N; i++) {
            if (!Double.isFinite(y[i])) {
                throw new IllegalArgumentException("The sample path must be finite");
            }
        }

        List<String> warnings = new ArrayList<String>();

        // warmup: grow the batch size until the signed areas look random
        int b = opt.b0;
        int m = opt.m0;
        if (N < (long) b * m) {
            m = N / b;
        }
        if (m < 1) {
            throw new IllegalArgumentException(
                    "The sample path is too short for the initial batch count b0 = " + b);
        }
        int ell = 1;
        boolean atMax = false;
        boolean passed = false;
        while (true) {
            StsQuantileStats stats = Sim_sts_quantile_areas.sim_sts_quantile_areas(
                    Arrays.copyOfRange(y, 0, b * m), b, m, p, opt.weight);
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
            if (N < (long) b * mNext && mNext != N / b) {
                m = N / b;
            } else {
                m = mNext;
                if (N < (long) b * m) {
                    m = N / b;
                    atMax = true;
                }
            }
            if (m < 1) {
                break;
            }
        }
        if (!passed) {
            warnings.add("the warmup randomness test could not be passed at the largest "
                    + "admissible batch size, the sample path is too short");
        }

        // truncation: delete the first batch
        int truncated = m;
        double[] yt = Arrays.copyOfRange(y, truncated, N);
        int nStar = yt.length;

        // batch-count selection: four tests in order, b only decreases
        int v = 0;
        b = opt.s[v];
        m = nStar / b;
        boolean ok = true;
        StsQuantileStats stats = null;
        for (int stage = 1; stage <= 4 && ok; stage++) {
            while (true) {
                if (m < 1) {
                    ok = false;
                    break;
                }
                stats = Sim_sts_quantile_areas.sim_sts_quantile_areas(
                        Arrays.copyOfRange(yt, nStar - b * m, nStar), b, m, p, opt.weight);
                double[] sample = stage <= 2 ? stats.areas : stats.bqe;
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
                m = nStar / b;
            }
        }
        if (stats == null || m < 1) {
            throw new IllegalArgumentException(
                    "The sample path is too short to form " + smallest + " batches");
        }

        int n = stats.n;
        double estimate = stats.quantile;
        double lower;
        double upper;
        double half;
        boolean heuristic;
        if (ok) {
            half = SimDist.tinv(1.0 - alpha / 2.0, 2.0 * b - 1) * Math.sqrt(stats.Vp / n);
            lower = estimate - half;
            upper = estimate + half;
            heuristic = false;
        } else {
            warnings.add("a randomness or normality test failed at b = " + smallest
                    + ", the delivered interval is heuristic");
            heuristic = true;
            if (opt.force) {
                double[] ci = SimQuestHeuristic.interval(stats.bqe, estimate, stats.Ap,
                        stats.Np, n, alpha, true);
                lower = ci[0];
                upper = ci[1];
                half = (upper - lower) / 2.0;
            } else {
                lower = Double.NaN;
                upper = Double.NaN;
                half = Double.NaN;
            }
        }

        return new QuantileCIResult(estimate, lower, upper, half, b, m, n, 1, truncated,
                stats.Ap, stats.Np, stats.Vp, heuristic, warnings);
    }
}
