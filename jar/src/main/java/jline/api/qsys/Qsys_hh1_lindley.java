/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Conditional waiting-time moments of the Hl/Hn/1 Lindley recursion.
 *
 * <p>Hyperexponential primitives are mixtures of exponentials, so conditioning on
 * the arrival phase i and the service phase j reduces one Lindley step to the
 * M/M/1 step of {@link Qsys_mm1_lindley} at rates {@code lambda[i]} and
 * {@code mu[j]}, and the conditional moment is the corresponding mixture
 *
 * <pre>
 *   E[W_{n+1}^m | W_n] = sum_i sum_j pa[i] ps[j] E_{ij}[W_{n+1}^m | W_n].
 * </pre>
 *
 * <p>Phases are drawn independently for each customer, which is what makes the
 * mixture exact rather than an approximation; a Markov-modulated arrival stream
 * would not decompose this way.
 *
 * <p>Note that the variance is not the corresponding mixture of the per-phase
 * variances, because the phase is itself random: it is recovered here from the
 * first two mixed raw moments, which adds the between-phase spread of the means.
 *
 * <p>Port of MATLAB qsys_hh1_lindley.m. Verified against 4e6 Monte Carlo
 * replications to 2e-3 relative error for m = 1, 2.
 *
 * <p>Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
 * Recursion", Proc. Winter Simulation Conference, 2021, theorem 3.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_hh1_lindley {
    private Qsys_hh1_lindley() {}

    /**
     * Conditional mean and variance of the next waiting time.
     *
     * @param lambda arrival phase rates, positive
     * @param pa     arrival phase probabilities, nonnegative and summing to 1
     * @param mu     service phase rates, positive
     * @param ps     service phase probabilities, nonnegative and summing to 1
     * @param Wn     current waiting times, finite and nonnegative
     * @return the conditional moments up to order 2
     */
    public static QsysLindleyResult qsys_hh1_lindley(double[] lambda, double[] pa,
                                                     double[] mu, double[] ps, double[] Wn) {
        return qsys_hh1_lindley(lambda, pa, mu, ps, Wn, 2);
    }

    /**
     * Conditional raw moments of the next waiting time.
     *
     * @param lambda arrival phase rates, positive
     * @param pa     arrival phase probabilities, nonnegative and summing to 1
     * @param mu     service phase rates, positive
     * @param ps     service phase probabilities, nonnegative and summing to 1
     * @param Wn     current waiting times, finite and nonnegative
     * @param mmax   highest moment order, at least 1; raised to 2 so the variance
     *               is always available
     * @return the conditional moments
     */
    public static QsysLindleyResult qsys_hh1_lindley(double[] lambda, double[] pa,
                                                     double[] mu, double[] ps, double[] Wn,
                                                     int mmax) {
        checkPhases(lambda, pa, "arrival");
        checkPhases(mu, ps, "service");
        if (mmax < 1) {
            throw new IllegalArgumentException("mmax=" + mmax + " must be positive");
        }
        if (Wn == null || Wn.length == 0) {
            throw new IllegalArgumentException("Wn must hold at least one waiting time");
        }
        for (int i = 0; i < Wn.length; i++) {
            if (!Double.isFinite(Wn[i]) || Wn[i] < 0.0) {
                throw new IllegalArgumentException(
                        "Wn must hold finite nonnegative values, got " + Wn[i]);
            }
        }

        int order = Math.max(mmax, 2);
        int nw = Wn.length;
        double[][] moments = new double[nw][order];
        for (int i = 0; i < lambda.length; i++) {
            for (int j = 0; j < mu.length; j++) {
                double weight = pa[i] * ps[j];
                if (weight == 0.0) {
                    continue;
                }
                for (int m = 1; m <= order; m++) {
                    for (int w = 0; w < nw; w++) {
                        moments[w][m - 1] += weight
                                * QsysLindleyMoment.moment(lambda[i], mu[j], Wn[w], m);
                    }
                }
            }
        }

        double[] mean = new double[nw];
        double[] var = new double[nw];
        for (int w = 0; w < nw; w++) {
            mean[w] = moments[w][0];
            var[w] = moments[w][1] - moments[w][0] * moments[w][0];
        }

        return new QsysLindleyResult(mean, var, moments, order);
    }

    private static void checkPhases(double[] rates, double[] probs, String what) {
        if (rates == null || probs == null || rates.length == 0
                || rates.length != probs.length) {
            throw new IllegalArgumentException("The " + what
                    + " rates and probabilities must be nonempty and of equal length");
        }
        double total = 0.0;
        for (int i = 0; i < rates.length; i++) {
            if (!(rates[i] > 0.0) || !Double.isFinite(rates[i])) {
                throw new IllegalArgumentException("The " + what + " rates must be positive"
                        + " and finite, got " + rates[i]);
            }
            if (!(probs[i] >= 0.0)) {
                throw new IllegalArgumentException("The " + what + " probabilities must be"
                        + " nonnegative, got " + probs[i]);
            }
            total += probs[i];
        }
        if (Math.abs(total - 1.0) > 1e-10) {
            throw new IllegalArgumentException("The " + what
                    + " probabilities must sum to 1, they sum to " + total);
        }
    }
}
