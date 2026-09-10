/**
 * @file M/GI/1 non-preemptive LCFS-S Age of Information analysis
 *
 * Computes mean AoI for an M/GI/1 queue with non-preemptive Last-Come
 * First-Served with Set-aside (LCFS-S) discipline.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfss_mgi1 {
    private Aoi_lcfss_mgi1() {}

    /**
     * Mean AoI for M/GI/1 non-preemptive LCFS-S queue.
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param H_lst LST of service time distribution (for reference; may be null)
     * @param E_H Mean service time (first moment), must be positive
     * @param E_H2 Second moment of service time, must be &gt;= E_H^2
     * @return AoiLstResult containing meanAoI, lstAoI (null), peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_lcfss_mgi1(double lambda, LstFunction H_lst, double E_H, double E_H2) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(E_H > 0)) throw new IllegalArgumentException("Mean service time E_H must be positive");
        if (!(E_H2 >= E_H * E_H)) throw new IllegalArgumentException("Second moment E_H2 must be >= E_H^2");

        double rho = lambda * E_H;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = lambda*E_H = %.4f >= 1", rho));
        }

        // Mean interarrival time
        double E_Y = 1.0 / lambda;

        // Busy period moments
        double E_B = E_H / (1.0 - rho);
        double E_B2 = E_H2 / ((1.0 - rho) * (1.0 - rho) * (1.0 - rho));

        // Mean AoI for LCFS-S (Proposition 5, simplified form)
        double meanAoI = E_Y + E_H + lambda * E_H2 / (2.0 * (1.0 - rho) * (1.0 - rho));

        // Mean Peak AoI (residual busy period)
        double peakAoI = E_Y + E_H + lambda * E_B2 / (2.0 * E_B);

        // LST is complex for LCFS-S; return null
        return new AoiLstResult(meanAoI, null, peakAoI);
    }
}
