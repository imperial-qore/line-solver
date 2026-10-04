/**
 * @file M/GI/1 non-preemptive LCFS-S Age of Information analysis
 *
 * Computes the mean AoI, its LST and the mean peak AoI of an M/GI/1 queue with
 * non-preemptive Last-Come First-Served with Set-aside (LCFS-S), i.e. LCFS
 * without discarding: the freshest waiting update is served next, and older
 * ones stay in queue and are still served later.
 *
 * Exact, Inoue, Masuyama, Takine and Tanaka (IEEE Trans. IT 65(12), 2019),
 * Section 3.3, NP-LCFS (D), with H*(s) the service LST:
 *   E[A]     = lambda E[H^2]/2 + ((1-rho)^2/(rho H*(lambda)) + 2) E[H]      (eq. 70)
 *   E[Apeak] = E[H] + E[W] + 1/(lambda (2 - rho - H*(lambda)))
 *   E[W]     = ((1 - H*(lambda))/lambda + H*'(lambda))/(2 - rho - H*(lambda))
 *   A*(s)    = eq. (66)
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfss_mgi1 {
    private Aoi_lcfss_mgi1() {}

    /**
     * Mean AoI, LST and peak AoI for the M/GI/1 non-preemptive LCFS-S queue.
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param H_lst LST of service time distribution
     * @param E_H Mean service time (first moment), must be positive
     * @param E_H2 Second moment of service time, must be &gt;= E_H^2
     * @return AoiLstResult containing meanAoI, lstAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_lcfss_mgi1(final double lambda, final LstFunction H_lst, final double E_H, double E_H2) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(E_H > 0)) throw new IllegalArgumentException("Mean service time E_H must be positive");
        if (!(E_H2 >= E_H * E_H)) throw new IllegalArgumentException("Second moment E_H2 must be >= E_H^2");
        if (H_lst == null) throw new IllegalArgumentException("The service LST H_lst is required");

        final double rho = lambda * E_H;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = lambda*E_H = %.4f >= 1", rho));
        }

        double hL = H_lst.evaluate(lambda);
        double hstep = 1e-6 * Math.max(1.0, lambda);
        double mdH = -(H_lst.evaluate(lambda + hstep) - H_lst.evaluate(lambda - hstep)) / (2.0 * hstep);

        double meanAoI = lambda * E_H2 / 2.0 + ((1.0 - rho) * (1.0 - rho) / (rho * hL) + 2.0) * E_H;
        // Informative arrivals: lambda_dagger = lambda (2 - rho - H*(lambda)); E[W] from eq. (98)
        double den = 2.0 - rho - hL;
        double E_W = ((1.0 - hL) / lambda - mdH) / den;
        double peakAoI = E_H + E_W + 1.0 / (lambda * den);

        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                if (Math.abs(s) < 1e-12) return 1.0;  // removable singularity of the residual-service LST
                double Hs = H_lst.evaluate(s);
                double HsL = H_lst.evaluate(s + lambda);
                double hresS = (1.0 - Hs) / (s * E_H);
                return lambda / (s + lambda) * Hs * (rho * hresS
                        + (1.0 - rho) * (s + lambda) * (1.0 - Hs + HsL) / (s + lambda * HsL));
            }
        };
        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }
}
