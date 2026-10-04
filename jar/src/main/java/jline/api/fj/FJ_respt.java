/**
 * @file Fork-Join response time approximations
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import jline.util.Maths;

public final class FJ_respt {
    private FJ_respt() {}

    /**
     * Exact two-way Fork-Join response time R_2^{F/J}(rho).
     */
    public static double fj_respt_2way(double lambda, double mu) {
        double rho = lambda / mu;
        if (rho >= 1.0) {
            throw new IllegalArgumentException(
                    "System is unstable: rho = lambda/mu = " + String.format("%.4f", rho)
                            + " >= 1. Require lambda < mu.");
        }

        double R_rho = 1.0 / (mu - lambda);
        return (12.0 - rho) / 8.0 * R_rho;
    }

    /**
     * Nelson-Tantawi approximation for K-way F/J response time.
     */
    public static double fj_respt_nt(int K, double lambda, double mu) {
        if (K < 2) {
            throw new IllegalArgumentException("Nelson-Tantawi approximation requires K >= 2. Got K=" + K + ".");
        }

        double rho = lambda / mu;
        if (rho >= 1.0) {
            throw new IllegalArgumentException(
                    "System is unstable: rho = lambda/mu = " + String.format("%.4f", rho)
                            + " >= 1. Require lambda < mu.");
        }

        double H_K = FJ_harmonic.fj_harmonic(K);
        double H_2 = FJ_harmonic.fj_harmonic(2);

        double S_K = H_K / H_2 + (1.0 - H_K / H_2) * (4.0 * rho / 11.0);
        double R_2_factor = 1.5 - rho / 8.0;
        double R_rho = 1.0 / (mu - lambda);

        return S_K * R_2_factor * R_rho;
    }

    /**
     * Varma-Makowski approximation for K-way F/J response time.
     */
    public static double fj_respt_vm(int K, double lambda, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }

        double rho = lambda / mu;
        if (rho >= 1.0) {
            throw new IllegalArgumentException(
                    "System is unstable: rho = lambda/mu = " + String.format("%.4f", rho)
                            + " >= 1. Require lambda < mu.");
        }

        double H_K = FJ_harmonic.fj_harmonic(K);

        double A_K = 0.0;
        for (int i = 1; i <= K; i++) {
            double innerSum = 0.0;
            for (int m = 1; m <= i; m++) {
                innerSum += Maths.binomialCoeff(i, m) * Maths.fact((double) (m - 1))
                        / Math.pow((double) i, (double) (m + 1));
            }
            double sign = ((i - 1) % 2 == 0) ? 1.0 : -1.0;
            A_K += Maths.binomialCoeff(K, i) * sign * innerSum;
        }

        double R_rho = 1.0 / (mu - lambda);
        return (H_K + (A_K - H_K) * rho) * R_rho;
    }

    /**
     * Varki et al. approximation for K-way F/J response time.
     */
    public static double fj_respt_varki(int K, double lambda, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }

        double rho = lambda / mu;
        if (rho >= 1.0) {
            throw new IllegalArgumentException(
                    "System is unstable: rho = lambda/mu = " + String.format("%.4f", rho)
                            + " >= 1. Require lambda < mu.");
        }

        double H_K = FJ_harmonic.fj_harmonic(K);

        double S1 = 0.0;
        for (int i = 1; i <= K; i++) {
            S1 += 1.0 / (i - rho);
        }

        double S2 = 0.0;
        for (int i = 1; i <= K; i++) {
            S2 += 1.0 / ((double) i * (i - rho));
        }

        return (1.0 / mu) * (H_K + (rho / (2.0 * (1.0 - rho))) * (S1 + (1.0 - 2.0 * rho) * S2));
    }
}
