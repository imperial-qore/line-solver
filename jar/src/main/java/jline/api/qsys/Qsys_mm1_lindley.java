/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Conditional waiting-time moments of the M/M/1 Lindley recursion.
 *
 * <p>One step of Lindley's recursion
 * {@code W_{n+1} = max(W_n + S_n - A_n, 0)} with {@code A_n ~ Exp(lambda)} and
 * {@code S_n ~ Exp(mu)}: given the waiting time of customer n, the exact
 * conditional moments of the waiting time of customer n+1. Unlike every other
 * {@code qsys_*} algorithm these are conditional on the current state rather than
 * stationary, so they are defined and finite for any load, including
 * {@code lambda >= mu}.
 *
 * <p>The m-th conditional moment is
 *
 * <pre>
 *   E[W_{n+1}^m | W_n] = lambda mu/(lambda+mu) [ S + T ],
 *   S = sum_{k=0}^{m} C(m,k) W_n^k (m-k)! / mu^(m-k+1),
 *   T = (-1)^m e^{-lambda W_n} (Gamma(m+1,-lambda W_n) - m!) / lambda^(m+1),
 * </pre>
 *
 * <p>where the density of {@code S_n - A_n} is the asymmetric Laplace density
 * {@code lambda mu/(lambda+mu)} times {@code e^{-mu x}} for {@code x > 0} and
 * {@code e^{lambda x}} for {@code x < 0}. Because {@code m+1} is a positive
 * integer, the upper incomplete gamma function admits the finite form
 * {@code Gamma(m+1,x) = m! e^{-x} sum_{k=0}^{m} x^k/k!}, valid at the negative
 * argument {@code -lambda W_n} needed here. Substituting it cancels the growing
 * exponential and leaves the numerically stable
 *
 * <pre>
 *   T = (-1)^m m! ( sum_{k=0}^{m} (-lambda W_n)^k/k! - e^{-lambda W_n} ) / lambda^(m+1),
 * </pre>
 *
 * <p>which is what this class evaluates. No incomplete gamma routine is needed.
 * The mean is returned from the equivalent explicit form
 * {@code W_n + (lambda-mu)/(lambda mu) + mu e^{-lambda W_n}/(lambda(lambda+mu))},
 * and the variance as the second moment less the squared mean.
 *
 * <p>Port of MATLAB qsys_mm1_lindley.m. Verified against 4e6 Monte Carlo
 * replications to 5e-4 relative error for m = 1, 2, 3.
 *
 * <p>Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
 * Recursion", Proc. Winter Simulation Conference, 2021, theorem 1 and
 * corollary 2.
 *
 * @since LINE 3.1.0
 */
public final class Qsys_mm1_lindley {
    private Qsys_mm1_lindley() {}

    /**
     * Conditional mean and variance of the next waiting time.
     *
     * @param lambda arrival rate, positive
     * @param mu     service rate, positive
     * @param Wn     current waiting times, finite and nonnegative
     * @return the conditional moments up to order 2
     */
    public static QsysLindleyResult qsys_mm1_lindley(double lambda, double mu, double[] Wn) {
        return qsys_mm1_lindley(lambda, mu, Wn, 2);
    }

    /**
     * Conditional raw moments of the next waiting time.
     *
     * @param lambda arrival rate, positive
     * @param mu     service rate, positive
     * @param Wn     current waiting times, finite and nonnegative
     * @param mmax   highest moment order, at least 1; raised to 2 so the variance
     *               is always available
     * @return the conditional moments
     */
    public static QsysLindleyResult qsys_mm1_lindley(double lambda, double mu, double[] Wn,
                                                     int mmax) {
        if (!(lambda > 0.0) || !Double.isFinite(lambda)) {
            throw new IllegalArgumentException("lambda=" + lambda + " must be positive and finite");
        }
        if (!(mu > 0.0) || !Double.isFinite(mu)) {
            throw new IllegalArgumentException("mu=" + mu + " must be positive and finite");
        }
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
        for (int m = 1; m <= order; m++) {
            for (int i = 0; i < nw; i++) {
                moments[i][m - 1] = QsysLindleyMoment.moment(lambda, mu, Wn[i], m);
            }
        }

        double[] mean = new double[nw];
        double[] var = new double[nw];
        for (int i = 0; i < nw; i++) {
            mean[i] = Wn[i] + (lambda - mu) / (lambda * mu)
                    + mu * Math.exp(-lambda * Wn[i]) / (lambda * (lambda + mu));
            var[i] = moments[i][1] - moments[i][0] * moments[i][0];
        }

        return new QsysLindleyResult(mean, var, moments, order);
    }
}
