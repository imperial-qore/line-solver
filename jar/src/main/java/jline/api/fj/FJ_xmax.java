/**
 * @file Expected maximum of K i.i.d. random variables for Fork-Join analysis.
 *
 * @since LINE 3.0
 */
package jline.api.fj;

import jline.util.Maths;

public final class FJ_xmax {
    private FJ_xmax() {}

    public static double fj_harmonic(int K) {
        double sum = 0.0;
        for (int i = 1; i <= K; i++) sum += 1.0 / i;
        return sum;
    }

    /**
     * Expected maximum of K i.i.d. exponential random variables.
     */
    public static double fj_xmax_exp(int K, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (mu <= 0.0) {
            throw new IllegalArgumentException("Service rate mu must be positive. Got mu=" + String.format("%.4f", mu) + ".");
        }
        return fj_harmonic(K) / mu;
    }

    /**
     * Expected maximum of 2 exponential random variables.
     */
    public static double fj_xmax_2(double lambda1, double lambda2) {
        if (lambda1 <= 0.0 || lambda2 <= 0.0) {
            throw new IllegalArgumentException(
                "Rates must be positive. Got lambda1=" + String.format("%.4f", lambda1)
                + ", lambda2=" + String.format("%.4f", lambda2) + ".");
        }
        return 1.0 / lambda1 + 1.0 / lambda2 - 1.0 / (lambda1 + lambda2);
    }

    /**
     * Expected maximum of 2 i.i.d. exponential random variables (same rate).
     */
    public static double fj_xmax_2(double lambda) {
        return fj_xmax_2(lambda, lambda);
    }

    /**
     * Expected maximum of K i.i.d. Erlang-k random variables.
     */
    public static double fj_xmax_erlang(int K, int k, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (k < 1) {
            throw new IllegalArgumentException("k (number of stages) must be a positive integer. Got k=" + k + ".");
        }
        if (mu <= 0.0) {
            throw new IllegalArgumentException("Rate mu must be positive. Got mu=" + String.format("%.4f", mu) + ".");
        }

        if (k == 2) {
            double outerSum = 0.0;
            for (int n = 1; n <= K; n++) {
                double innerSum = 0.0;
                for (int m = 1; m <= n; m++) {
                    innerSum += Maths.binomialCoeff(n, m) * Maths.fact((double) m) / (2.0 * Math.pow((double) n, (double) (m + 1)));
                }
                double sign = ((n - 1) % 2 == 0) ? 1.0 : -1.0;
                outerSum += Maths.binomialCoeff(K, n) * sign * innerSum;
            }
            return outerSum / mu;
        } else {
            double meanErlang = (double) k / mu;
            double upperLimit = meanErlang * 10.0 + 10.0 * Math.sqrt((double) k) / mu;

            int numPoints = 10001;
            double h = upperLimit / (numPoints - 1);
            double integral = 0.0;
            for (int i = 0; i < numPoints; i++) {
                double t = i * h;
                double cdfVal = erlangCdfGeneral(t, k, mu);
                double integrandVal = 1.0 - Math.pow(cdfVal, (double) K);
                double weight;
                if (i == 0 || i == numPoints - 1) weight = 1.0;
                else if (i % 2 == 1) weight = 4.0;
                else weight = 2.0;
                integral += weight * integrandVal;
            }
            integral *= h / 3.0;
            return integral;
        }
    }

    /**
     * Expected maximum of K i.i.d. Hyperexponential-2 random variables.
     */
    public static double fj_xmax_hyperexp(int K, double p1, double mu1, double mu2) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (p1 <= 0.0 || p1 >= 1.0) {
            throw new IllegalArgumentException("p1 must be in (0,1). Got p1=" + String.format("%.4f", p1) + ".");
        }
        if (mu1 <= 0.0 || mu2 <= 0.0) {
            throw new IllegalArgumentException("Rates mu1 and mu2 must be positive.");
        }
        double p2 = 1.0 - p1;
        double Xmax = 0.0;
        for (int n = 1; n <= K; n++) {
            double innerSum = 0.0;
            for (int m = 0; m <= n; m++) {
                double denominator = m * mu1 + (n - m) * mu2;
                if (denominator > 0.0) {
                    innerSum += Maths.binomialCoeff(n, m)
                            * Math.pow(p1, (double) m)
                            * Math.pow(p2, (double) (n - m))
                            / denominator;
                }
            }
            double sign = ((n + 1) % 2 == 0) ? 1.0 : -1.0;
            Xmax += sign * innerSum;
        }
        return Xmax;
    }

    /**
     * Result container for expected maximum of normal distribution.
     */
    public static final class FJXmaxNormalResult {
        public final double Xmax;
        public final double Vmax;
        public FJXmaxNormalResult(double Xmax, double Vmax) {
            this.Xmax = Xmax;
            this.Vmax = Vmax;
        }
        public double getXmax() { return Xmax; }
        public double getVmax() { return Vmax; }
    }

    /**
     * Expected maximum for normal distribution.
     */
    public static FJXmaxNormalResult fj_xmax_normal(int K, double mu, double sigma, String method) {
        if (K < 2) {
            throw new IllegalArgumentException("K must be at least 2 for normal approximation. Got K=" + K + ".");
        }
        if (sigma < 0.0) {
            throw new IllegalArgumentException("Standard deviation sigma must be non-negative.");
        }

        double gammaEM = 0.5772156649015329;
        double sqrt2lnK = Math.sqrt(2.0 * Math.log(K));

        double GK;
        String m = method.toLowerCase();
        if (m.equals("arnold")) {
            GK = sqrt2lnK;
        } else if (m.equals("johnson")) {
            double correction = (Math.log(Math.log(K)) - Math.log(4.0 * Math.PI) + 2.0 * gammaEM) / (2.0 * sqrt2lnK);
            GK = sqrt2lnK - correction;
        } else if (m.equals("corrected")) {
            double correction = (Math.log(Math.log(K)) - Math.log(4.0 * Math.PI) + 2.0 * gammaEM) / (2.0 * sqrt2lnK);
            double deltaK = 0.1727 * Math.pow((double) K, -0.2750);
            GK = sqrt2lnK - correction - deltaK;
        } else {
            throw new IllegalArgumentException("Unknown method: " + method + ". Valid: johnson, arnold, corrected.");
        }

        double Xmax = mu + sigma * GK;
        double Vmax = 1.64492 * sigma * sigma / (2.0 * Math.log(K));
        return new FJXmaxNormalResult(Xmax, Vmax);
    }

    public static FJXmaxNormalResult fj_xmax_normal(int K, double mu, double sigma) {
        return fj_xmax_normal(K, mu, sigma, "johnson");
    }

    /**
     * Result container for expected maximum of Pareto distribution.
     */
    public static final class FJXmaxParetoResult {
        public final double Xmax;
        public final double MK;
        public FJXmaxParetoResult(double Xmax, double MK) {
            this.Xmax = Xmax;
            this.MK = MK;
        }
        public double getXmax() { return Xmax; }
        public double getMK() { return MK; }
    }

    /**
     * Expected maximum for Pareto distribution.
     */
    public static double fj_xmax_pareto(int K, double beta, double k) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (beta <= 2.0) {
            throw new IllegalArgumentException("Shape parameter beta must be > 2 for finite moments. Got beta=" + String.format("%.4f", beta) + ".");
        }
        if (k <= 0.0) {
            throw new IllegalArgumentException("Scale parameter k must be positive.");
        }

        double upperLimit = k * Math.pow((double) K, 2.0 / beta) * 10.0;
        int numPoints = 10001;
        double h = upperLimit / (numPoints - 1);
        double integral = 0.0;
        for (int i = 0; i < numPoints; i++) {
            double x = i * h;
            double F = 1.0 - Math.pow(k / (k + Math.max(x, 0.0)), beta);
            double integrandVal = 1.0 - Math.pow(F, (double) K);
            double weight;
            if (i == 0 || i == numPoints - 1) weight = 1.0;
            else if (i % 2 == 1) weight = 4.0;
            else weight = 2.0;
            integral += weight * integrandVal;
        }
        integral *= h / 3.0;
        return integral;
    }

    public static double fj_xmax_pareto(int K, double beta) {
        return fj_xmax_pareto(K, beta, beta - 1.0);
    }

    /**
     * Characteristic maximum M_K for Pareto distribution.
     */
    public static double fj_xmax_pareto_char_max(int K, double beta, double k) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (beta <= 1.0) {
            return Double.POSITIVE_INFINITY;
        }
        double mK = k * (Math.pow((double) K, 1.0 / beta) - 1.0);
        double tailIntegral = Math.pow(k, beta) * Math.pow(k + mK, 1.0 - beta) / (beta - 1.0);
        return mK + K * tailIntegral;
    }

    /**
     * General approximation for expected maximum of K random variables.
     */
    public static double fj_xmax_approx(int K, double muX, double sigmaX, String distType) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (sigmaX < 0.0) {
            throw new IllegalArgumentException("Standard deviation sigma_X must be non-negative.");
        }

        double GK;
        String dt = distType.toLowerCase();
        if (dt.equals("exp")) {
            GK = fj_harmonic(K) - 1.0;
        } else if (dt.equals("uniform")) {
            GK = Math.sqrt(3.0) * (K - 1) / (K + 1);
        } else if (dt.equals("evd")) {
            GK = Math.sqrt(6.0) * Math.log(K) / Math.PI;
        } else if (dt.equals("bound")) {
            GK = (double) (K - 1) / Math.sqrt(2.0 * K - 1.0);
        } else {
            throw new IllegalArgumentException("Unknown distribution type: " + distType + ". Valid: exp, uniform, evd, bound.");
        }
        return muX + sigmaX * GK;
    }

    public static double fj_xmax_approx(int K, double muX, double sigmaX) {
        return fj_xmax_approx(K, muX, sigmaX, "exp");
    }

    /**
     * Expected maximum using EMMA.
     */
    public static double fj_xmax_emma(int K, double mu) {
        double PHI = 0.570376;
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (mu <= 0.0) {
            throw new IllegalArgumentException("Rate mu must be positive.");
        }
        return -(1.0 / mu) * Math.log(1.0 - Math.pow(PHI, 1.0 / K));
    }

    private static double erlangCdfGeneral(double t, int k, double mu) {
        if (t <= 0.0) return 0.0;
        double S = 0.0;
        for (int j = 0; j < k; j++) {
            S += Math.pow(mu * t, (double) j) / Maths.fact((double) j);
        }
        return 1.0 - Math.exp(-mu * t) * S;
    }
}
