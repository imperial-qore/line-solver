/**
 * @file M/M/1 preemptive LCFS Age of Information analysis
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfspr_mm1 {
    private Aoi_lcfspr_mm1() {}

    /**
     * Mean, variance, and peak AoI for M/M/1 preemptive LCFS queue.
     *
     * <p>Exact: E[A] = 1/lambda + 1/mu;
     * E[Apeak] = 1/(lambda+mu) + 1/lambda + 1/mu.</p>
     *
     * Note: LCFS-PR always achieves lower mean AoI than FCFS for M/M/1.
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param mu Service rate (exponential service), must be positive
     * @return AoiResult containing meanAoI, varAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiResult aoi_lcfspr_mm1(double lambda, double mu) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");

        double rho = lambda / mu;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = lambda/mu = %.4f >= 1", rho));
        }

        // Mean AoI for preemptive LCFS (Proposition 3 / Section IV)
        double meanAoI = (1.0 / mu) * (1.0 + 1.0 / rho);

        // Mean Peak AoI (exact, as MATLAB): E[S|success] + E[inter-delivery]
        double peakAoI = 1.0 / (lambda + mu) + 1.0 / lambda + 1.0 / mu;

        // Variance of AoI for preemptive LCFS
        double E_A = meanAoI;
        double E_A2 = 2.0 * (1.0 / (lambda * lambda) + 1.0 / (lambda * mu) + 1.0 / (mu * mu));
        double varAoI = Math.max(0.0, E_A2 - E_A * E_A);

        return new AoiResult(meanAoI, varAoI, peakAoI);
    }
}
