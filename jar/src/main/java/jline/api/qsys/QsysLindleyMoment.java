/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * One conditional Lindley moment for exponential primitives.
 *
 * <p>Evaluates {@code E[max(Wn + S - A, 0)^m]} with {@code A ~ Exp(lambda)} and
 * {@code S ~ Exp(mu)}. This is the algorithm shared by {@link Qsys_mm1_lindley},
 * which calls it once per moment order, and {@link Qsys_hh1_lindley}, which mixes
 * it over the arrival and service phases. See {@link Qsys_mm1_lindley} for the
 * derivation and for why the upper incomplete gamma function reduces to a finite
 * sum here.
 *
 * <p>Port of MATLAB qsys_lindley_moment.m.
 *
 * @since LINE 3.1.0
 */
public final class QsysLindleyMoment {
    private QsysLindleyMoment() {}

    /**
     * Evaluates the conditional moment.
     *
     * @param lambda arrival rate, positive
     * @param mu     service rate, positive
     * @param Wn     current waiting time, finite and nonnegative
     * @param m      moment order, at least 1
     * @return {@code E[max(Wn + S - A, 0)^m]}
     */
    public static double moment(double lambda, double mu, double Wn, int m) {
        if (m < 1) {
            throw new IllegalArgumentException("The moment order m=" + m + " must be positive");
        }

        double sTerm = 0.0;
        for (int k = 0; k <= m; k++) {
            sTerm += binomial(m, k) * Math.pow(Wn, k) * factorial(m - k)
                    / Math.pow(mu, m - k + 1);
        }

        double x = -lambda * Wn;
        double inner = 0.0;
        for (int k = 0; k <= m; k++) {
            inner += Math.pow(x, k) / factorial(k);
        }
        double sign = (m % 2 == 0) ? 1.0 : -1.0;
        double tTerm = sign * factorial(m) * (inner - Math.exp(x)) / Math.pow(lambda, m + 1);

        return lambda * mu / (lambda + mu) * (sTerm + tTerm);
    }

    private static double factorial(int k) {
        double v = 1.0;
        for (int i = 2; i <= k; i++) {
            v *= i;
        }
        return v;
    }

    private static double binomial(int n, int k) {
        double v = 1.0;
        for (int i = 1; i <= k; i++) {
            v = v * (n - k + i) / i;
        }
        return v;
    }
}
