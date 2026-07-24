/**
 * @file M/GI/1 non-preemptive LCFS-D Age of Information analysis
 *
 * Computes mean AoI for an M/GI/1 queue with non-preemptive Last-Come
 * First-Served with Discarding (LCFS-D) discipline.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfsd_mgi1 {
    private Aoi_lcfsd_mgi1() {}

    /**
     * Mean AoI for M/GI/1 non-preemptive LCFS-D queue.
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param H_lst LST of service time distribution (for reference; may be null)
     * @param E_H Mean service time (first moment), must be positive
     * @param E_H2 Second moment of service time, must be &gt;= E_H^2
     * @return AoiLstResult containing meanAoI, lstAoI (null), peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_lcfsd_mgi1(double lambda, LstFunction H_lst, double E_H, double E_H2) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(E_H > 0)) throw new IllegalArgumentException("Mean service time E_H must be positive");
        if (!(E_H2 >= E_H * E_H)) throw new IllegalArgumentException("Second moment E_H2 must be >= E_H^2");

        double rho = lambda * E_H;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = lambda*E_H = %.4f >= 1", rho));
        }

        // Mean interarrival time
        double E_Y = 1.0 / lambda;

        // Residual service time
        double E_H_residual = E_H2 / (2.0 * E_H);

        // Mean effective system time
        double E_T_eff = E_H + rho * E_H_residual;

        // Mean AoI for LCFS-D (Proposition 6, simplified form)
        double meanAoI = E_Y + E_H + rho * E_H2 / (2.0 * E_H) + rho * E_H / (1.0 + rho);

        // Mean Peak AoI
        double peakAoI = E_Y + E_T_eff;

        // LST is complex for LCFS-D; return null
        return new AoiLstResult(meanAoI, null, peakAoI);
    }
}
