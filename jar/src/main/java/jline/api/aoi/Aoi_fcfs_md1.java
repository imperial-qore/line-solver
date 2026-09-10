/**
 * @file M/D/1 FCFS Age of Information analysis
 *
 * Computes mean, variance, and peak AoI for an M/D/1 queue with
 * First-Come First-Served (FCFS) discipline.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_fcfs_md1 {
    private Aoi_fcfs_md1() {}

    /**
     * Mean, variance, and peak AoI for M/D/1 FCFS queue.
     *
     * <p>Exact mean: E[A] = d*(1/2 + 1/(2*(1-rho)) + ((1-rho)/rho)*exp(rho)),
     * rho = lambda*d; E[Apeak] = E[T] + E[Y].</p>
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param d Deterministic service time, must be positive
     * @return AoiResult containing meanAoI, varAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiResult aoi_fcfs_md1(double lambda, double d) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(d > 0)) throw new IllegalArgumentException("Service time d must be positive");

        double rho = lambda * d;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = lambda*d = %.4f >= 1", rho));
        }

        // Service time moments (deterministic)
        double E_H = d;
        double E_H2 = d * d;

        // Mean waiting time (Pollaczek-Khinchine for M/G/1)
        double E_W = lambda * E_H2 / (2.0 * (1.0 - rho));

        // Mean system time (sojourn time)
        double E_T = E_W + E_H;

        // Mean interarrival time
        double E_Y = 1.0 / lambda;

        // Mean AoI for M/D/1 FCFS (exact, as MATLAB aoi_fcfs_md1)
        double meanAoI = d * (0.5 + 1.0 / (2.0 * (1.0 - rho))
                + ((1.0 - rho) / rho) * Math.exp(rho));

        // Mean Peak AoI
        double peakAoI = E_T + E_Y;

        // Variance of AoI
        double E_Y2 = 2.0 / (lambda * lambda); // Second moment of exponential
        double varAoI = Math.max(0.0,
                E_Y2 - E_Y * E_Y + 2.0 * E_W * E_H / (1.0 - rho)
                        + E_H2 * rho / ((1.0 - rho) * (1.0 - rho)));

        return new AoiResult(meanAoI, varAoI, peakAoI);
    }
}
