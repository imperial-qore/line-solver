/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

/**
 * Procedure constants of {@link Sim_fquest} and {@link Sim_firquest}.
 *
 * <p>The defaults are the ones the articles report after their own
 * experimentation: {@code b0 = 50} gives the warmup randomness test enough
 * power, 32 batches suffice to estimate the variance parameter while fewer than
 * 10 make the interval unreliable, and the decaying warmup significance keeps the
 * batch size from growing so far that truncation eats a short sample. With these
 * values the fourth warmup iteration runs at
 * {@code beta*exp(-0.2*3^2.3) = 0.025}.
 *
 * <p>{@link Sim_firquest} substitutes {@code b0 = 25} and an {@code s} chosen from
 * the replication count; see {@link Sim_firquest#defaultBatchCounts(int)}.
 *
 * <p>Port of MATLAB sim_quest_options.m.
 *
 * @since LINE 3.1.0
 */
public class QuestOptions {
    /** Initial batch count for the warmup stage. */
    public int b0 = 50;
    /** Initial batch size for the warmup stage. */
    public int m0 = 500;
    /** Descending batch counts for the test stages. */
    public int[] s = {32, 24, 16, 10};
    /** Significance level of the stage tests. */
    public double beta = 0.30;
    /** Decay coefficient of the warmup significance. */
    public double eta = 0.2;
    /** Decay exponent of the warmup significance. */
    public double theta = 2.3;
    /** Constant STS weight function. */
    public double weight = Sim_sts_quantile_areas.DEFAULT_WEIGHT;
    /** Whether to deliver a heuristic interval when a stage test fails. */
    public boolean force = true;

    /**
     * Rejects inadmissible option values.
     *
     * @throws IllegalArgumentException when a field is out of range
     */
    public void validate() {
        if (b0 < 3) {
            throw new IllegalArgumentException("b0=" + b0 + " must be at least 3");
        }
        if (m0 < 1) {
            throw new IllegalArgumentException("m0=" + m0 + " must be positive");
        }
        if (s == null || s.length == 0) {
            throw new IllegalArgumentException("s must be nonempty");
        }
        for (int i = 0; i < s.length; i++) {
            if (s[i] < 1) {
                throw new IllegalArgumentException("s[" + i + "]=" + s[i] + " must be positive");
            }
            if (i > 0 && s[i] >= s[i - 1]) {
                throw new IllegalArgumentException("s must be strictly decreasing");
            }
        }
        if (!(beta > 0.0) || !(beta < 1.0)) {
            throw new IllegalArgumentException("beta=" + beta + " must lie in (0,1)");
        }
        if (!(eta >= 0.0)) {
            throw new IllegalArgumentException("eta=" + eta + " must be nonnegative");
        }
        if (!(theta > 0.0)) {
            throw new IllegalArgumentException("theta=" + theta + " must be positive");
        }
        if (weight == 0.0 || !Double.isFinite(weight)) {
            throw new IllegalArgumentException("weight=" + weight + " must be nonzero and finite");
        }
    }
}
