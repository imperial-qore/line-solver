/**
 * @file Maximum response time computations for Fork-Join systems
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import jline.util.Maths;

public final class FJ_rmax {
    private FJ_rmax() {}

    /**
     * Compute maximum response time R_K^max(rho) for K M/M/1 queues.
     */
    public static double fj_rmax(int K, double lambda, double mu) {
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
        double R_rho = 1.0 / (mu - lambda);

        return H_K * R_rho;
    }

    /**
     * Maximum response time R_K^max for Erlang service times.
     */
    public static double fj_rmax_erlang(int K, int k, double lambda, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (k < 1) {
            throw new IllegalArgumentException("k (number of stages) must be a positive integer. Got k=" + k + ".");
        }

        double meanService = (double) k / mu;
        double rho = lambda * meanService;
        if (rho >= 1.0) {
            throw new IllegalArgumentException("System is unstable: rho = " + String.format("%.4f", rho) + " >= 1.");
        }

        double cv2 = 1.0 / k;
        double R_single = meanService * (1.0 + rho * (1.0 + cv2) / (2.0 * (1.0 - rho)));

        if (K == 2) {
            double muResp = (double) k / R_single;
            double correction = 0.0;
            for (int m = 0; m < k; m++) {
                for (int n = 0; n < k; n++) {
                    correction += Maths.binomialCoeff(m + n, m)
                            * Math.pow(muResp, (double) m)
                            * Math.pow(muResp, (double) n)
                            / Math.pow(2.0 * muResp, (double) (m + n + 1));
                }
            }
            return 2.0 * R_single - correction;
        } else {
            double cv2Response = Math.max((cv2 + rho) / (1.0 + rho), 1.0 / 20.0);
            int kResp = Math.max(1, (int) Math.ceil(1.0 / cv2Response));
            double muResp = (double) kResp / R_single;

            double upperLimit = R_single * 20.0;
            int numPoints = 10001;
            double h = upperLimit / (numPoints - 1);

            double integral = 0.0;
            for (int i = 0; i < numPoints; i++) {
                double t = i * h;
                double cdfVal = erlangCdf(t, kResp, muResp);
                double integrandVal = 1.0 - Math.pow(cdfVal, (double) K);

                double weight;
                if (i == 0 || i == numPoints - 1) {
                    weight = 1.0;
                } else if (i % 2 == 1) {
                    weight = 4.0;
                } else {
                    weight = 2.0;
                }
                integral += weight * integrandVal;
            }
            integral *= h / 3.0;

            return integral;
        }
    }

    /**
     * Maximum response time using Extreme Value Distribution (EVD) approximation.
     */
    public static double fj_rmax_evd(int K, double R, double sigmaR, boolean calibrated) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (R <= 0.0) {
            throw new IllegalArgumentException("Mean response time R must be positive.");
        }
        if (sigmaR < 0.0) {
            throw new IllegalArgumentException("Standard deviation sigma_R must be non-negative.");
        }

        double correctionFactor = Math.sqrt(6.0) * Math.log((double) K) / Math.PI;
        if (calibrated) {
            correctionFactor /= 1.27;
        }

        return R + correctionFactor * sigmaR;
    }

    public static double fj_rmax_evd(int K, double R, double sigmaR) {
        return fj_rmax_evd(K, R, sigmaR, false);
    }

    /**
     * Compute Erlang CDF: F(t) = 1 - exp(-mu*t) * sum_{j=0}^{k-1} (mu*t)^j / j!
     */
    private static double erlangCdf(double t, int k, double mu) {
        if (t <= 0.0) {
            return 0.0;
        }
        double S = 0.0;
        for (int j = 0; j < k; j++) {
            S += Math.pow(mu * t, (double) j) / Maths.fact((double) j);
        }
        return 1.0 - Math.exp(-mu * t) * S;
    }
}
