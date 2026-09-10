/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Exact sojourn-time moments of the multiclass M/M/1-PS queue.
 *
 * <p>Class j arrives in a Poisson stream of rate {@code lambda[j]} and requires
 * an exponential amount of service with rate {@code mu[j]}. The processor is
 * shared equally by all jobs in service, so the class of a job affects its
 * sojourn time both through its own service rate and through the mix of rates
 * of the jobs it shares the processor with. With
 * {@code alpha = 1 - sum_j lambda_j/mu_j} the unutilized fraction of the
 * processor, the moments of the sojourn time W_r of a tagged class-r job are</p>
 *
 * <pre>
 *   E[W_r]   = 1/(alpha*mu_r)
 *   E[W_r^2] = 2/(alpha*mu_r)^2 * [1 - sum_j lambda_j (mu_j-mu_r)/(mu_j(mu_j+mu_r))]
 *                               / [1 - sum_j lambda_j/(mu_j+mu_r)]
 * </pre>
 *
 * <p>which is equation (7) of Mitra and Morrison (1983). Both are exact, not
 * asymptotic: the open system is the N -&gt; infinity limit of the closed
 * terminal-driven system whose moments that paper expands in 1/N, and the
 * leading term of the expansion is exact in the limit. For a single class the
 * second moment reduces to the classical 4/(mu^2 (1-rho)^2 (2-rho)) of Coffman,
 * Muntz and Trotter (1970).</p>
 *
 * <p>Reference: D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of
 * the Waiting Time in Closed and Open Processor-Sharing Systems with Multiple
 * Job Classes", Adv. Appl. Prob. 15(4):813-839, 1983.</p>
 *
 * <p>Port of MATLAB qsys_mm1_ps.m.</p>
 */
public class Qsys_mm1_ps {

    private Qsys_mm1_ps() {
    }

    /**
     * Sojourn-time moments of the multiclass M/M/1-PS queue.
     *
     * @param lambda per-class Poisson arrival rates, non-negative
     * @param mu     per-class exponential service rates, positive
     * @return per-class mean and second moment of the sojourn time
     */
    public static QsysMm1PsResult qsys_mm1_ps(double[] lambda, double[] mu) {
        int R = lambda.length;
        if (mu.length != R) {
            throw new RuntimeException("qsys_mm1_ps: lambda and mu must have the same number of classes");
        }
        for (int j = 0; j < R; j++) {
            if (!isFinite(lambda[j]) || lambda[j] < 0) {
                throw new RuntimeException("qsys_mm1_ps: lambda must be finite and non-negative");
            }
            if (!isFinite(mu[j]) || mu[j] <= 0) {
                throw new RuntimeException("qsys_mm1_ps: mu must be finite and positive");
            }
        }
        double rho = 0.0;
        for (int j = 0; j < R; j++) {
            rho += lambda[j] / mu[j];
        }
        double alpha = 1.0 - rho;
        if (alpha <= 0) {
            throw new RuntimeException(String.format(
                    "qsys_mm1_ps: system is unstable: utilization %.6f >= 1", rho));
        }
        double[] W = new double[R];
        double[] W2 = new double[R];
        for (int r = 0; r < R; r++) {
            double mur = mu[r];
            double num = 1.0;
            double den = 1.0;
            for (int j = 0; j < R; j++) {
                num -= lambda[j] * (mu[j] - mur) / (mu[j] * (mu[j] + mur));
                den -= lambda[j] / (mu[j] + mur);
            }
            W[r] = 1.0 / (alpha * mur);
            W2[r] = 2.0 / ((alpha * mur) * (alpha * mur)) * num / den;
        }
        return new QsysMm1PsResult(W, W2, alpha);
    }

    private static boolean isFinite(double x) {
        return !Double.isNaN(x) && !Double.isInfinite(x);
    }
}
