/**
 * @file M/GI/1 non-preemptive LCFS-D Age of Information analysis
 *
 * Computes the mean AoI, its LST and the mean peak AoI of an M/GI/1 queue with
 * non-preemptive Last-Come First-Served with Discarding (LCFS-D, M/GI/1/2*):
 * a new arrival replaces the update waiting for the server, if any.
 *
 * Exact, Inoue, Masuyama, Takine and Tanaka (IEEE Trans. IT 65(12), 2019),
 * Section 3.3, NP-LCFS (C), with H*(s) the service LST:
 *   E[A]     = (lambda E[H^2]/2 + H*(lambda)/lambda - H*'(lambda))/(rho + H*(lambda))
 *              + (1 - H*(lambda))/lambda + H*'(lambda) + E[H]            (eq. 68)
 *   E[Apeak] = 1/lambda + H*'(lambda) + 2 E[H]
 *   A*(s)    = eq. (64)
 * Discarding keeps the system stable for every rho, so rho &gt;= 1 is accepted.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfsd_mgi1 {
    private Aoi_lcfsd_mgi1() {}

    /**
     * Mean AoI, LST and peak AoI for the M/GI/1 non-preemptive LCFS-D queue.
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param H_lst LST of service time distribution
     * @param E_H Mean service time (first moment), must be positive
     * @param E_H2 Second moment of service time, must be &gt;= E_H^2
     * @return AoiLstResult containing meanAoI, lstAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid
     */
    public static AoiLstResult aoi_lcfsd_mgi1(final double lambda, final LstFunction H_lst, final double E_H, double E_H2) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(E_H > 0)) throw new IllegalArgumentException("Mean service time E_H must be positive");
        if (!(E_H2 >= E_H * E_H)) throw new IllegalArgumentException("Second moment E_H2 must be >= E_H^2");
        if (H_lst == null) throw new IllegalArgumentException("The service LST H_lst is required");

        final double rho = lambda * E_H;
        // H*(lambda) = P(no arrival during a service), -H*'(lambda) = E[H exp(-lambda H)]
        final double hL = H_lst.evaluate(lambda);
        double hstep = 1e-6 * Math.max(1.0, lambda);
        double mdH = -(H_lst.evaluate(lambda + hstep) - H_lst.evaluate(lambda - hstep)) / (2.0 * hstep);

        double meanAoI = (lambda * E_H2 / 2.0 + hL / lambda + mdH) / (rho + hL)
                + (1.0 - hL) / lambda - mdH + E_H;
        double peakAoI = 1.0 / lambda - mdH + 2.0 * E_H;

        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                if (Math.abs(s) < 1e-12) return 1.0;  // removable singularity of the residual-service LST
                double Hs = H_lst.evaluate(s);
                double HsL = H_lst.evaluate(s + lambda);
                double hresS = (1.0 - Hs) / (s * E_H);
                double hresSL = (1.0 - HsL) / ((s + lambda) * E_H);
                return (hL + rho * hresSL) * Hs * (rho * hresS + HsL * lambda / (s + lambda)) / (rho + hL);
            }
        };
        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }
}
