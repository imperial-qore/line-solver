/**
 * @file D/M/1 FCFS Age of Information analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_fcfs_dm1 {
    private Aoi_fcfs_dm1() {}

    /**
     * Mean, variance, and peak AoI for D/M/1 FCFS queue.
     *
     * <p>Exact mean: E[A] = tau/2 + 1/(mu*(1-sigma)) with sigma the root of
     * sigma = exp(-mu*tau*(1-sigma)); E[Apeak] = tau + 1/(mu*(1-sigma)).</p>
     */
    public static AoiResult aoi_fcfs_dm1(double tau, double mu) {
        if (!(tau > 0)) throw new IllegalArgumentException("Interarrival time tau must be positive");
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");

        double lambda = 1.0 / tau;
        double rho = lambda / mu;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = 1/(tau*mu) = %.4f >= 1", rho));
        }

        double sigma = findSigmaDM1(mu, tau);

        double E_D = 1.0 / (mu * (1.0 - sigma));
        // Exact (as MATLAB aoi_fcfs_dm1): deterministic interarrivals kill the
        // Y-W correlation term, so E[A] = E[Y^2]/(2*E[Y]) + E[D]
        double meanAoI = tau / 2.0 + E_D;
        double peakAoI = tau + E_D;

        double E_D2 = 2.0 / ((mu * (1.0 - sigma)) * (mu * (1.0 - sigma)));
        double Var_D = E_D2 - E_D * E_D;
        double sigmaTerm = sigma / (mu * (1.0 - sigma));
        double varAoI = Math.max(0.0, Var_D + sigmaTerm * sigmaTerm);

        return new AoiResult(meanAoI, varAoI, peakAoI);
    }

    private static double findSigmaDM1(double mu, double tau) {
        double lo = 0.001;
        double hi = 0.999;
        int maxIter = 100;
        double tol = 1e-12;

        for (int iter = 0; iter < maxIter; iter++) {
            double mid = (lo + hi) / 2.0;
            double fMid = Math.exp(-mu * tau * (1.0 - mid)) - mid;
            if (Math.abs(fMid) < tol) return mid;
            double fLo = Math.exp(-mu * tau * (1.0 - lo)) - lo;
            if (fLo * fMid < 0) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        return (lo + hi) / 2.0;
    }
}
