/**
 * @file M/M/1 FCFS Age of Information analysis
 *
 * Computes mean, variance, and peak AoI for an M/M/1 queue with
 * First-Come First-Served (FCFS) discipline using closed-form formulas
 * from Inoue et al., IEEE Trans. IT, 2019.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_fcfs_mm1 {
    private Aoi_fcfs_mm1() {}

    /**
     * Mean, variance, and peak AoI for M/M/1 FCFS queue.
     *
     * <p>Exact (Inoue et al. 2019): E[A] = (1/mu)(1 + 1/rho + rho^2/(1-rho));
     * E[Apeak] = (1/mu)(1 + 1/rho + rho/(1-rho));
     * E[A^2] = (2/mu^2)(1 - rho - rho^3 + 4*rho^4 - 2*rho^5)/(rho^2*(1-rho)^2).</p>
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param mu Service rate (exponential service), must be positive
     * @return AoiResult containing meanAoI, varAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiResult aoi_fcfs_mm1(double lambda, double mu) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");

        double rho = lambda / mu;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = lambda/mu = %.4f >= 1", rho));
        }

        // Mean AoI (Proposition 1 / Corollary from Section III-A)
        double meanAoI = (1.0 / mu) * (1.0 + 1.0 / rho + rho * rho / (1.0 - rho));

        // Mean Peak AoI
        double peakAoI = (1.0 / mu) * (1.0 + 1.0 / rho + rho / (1.0 - rho));

        // Variance of AoI: exact second moment via the sawtooth identity
        // (as MATLAB aoi_fcfs_mm1)
        double E_A = meanAoI;
        double E_A2 = (2.0 / (mu * mu))
                * (1.0 - rho - Math.pow(rho, 3) + 4.0 * Math.pow(rho, 4) - 2.0 * Math.pow(rho, 5))
                / (rho * rho * (1.0 - rho) * (1.0 - rho));
        double varAoI = Math.max(0.0, E_A2 - E_A * E_A);

        return new AoiResult(meanAoI, varAoI, peakAoI);
    }
}
