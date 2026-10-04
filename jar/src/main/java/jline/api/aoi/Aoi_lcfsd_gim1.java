/**
 * @file GI/M/1 non-preemptive LCFS-D Age of Information analysis
 *
 * Computes the mean AoI, its LST and the mean peak AoI of a GI/M/1 queue with
 * non-preemptive Last-Come First-Served with Discarding (LCFS-D): a new
 * arrival replaces the update waiting for the server, if any.
 *
 * Exact, Inoue, Masuyama, Takine and Tanaka (IEEE Trans. IT 65(12), 2019),
 * Section 3.3, NP-LCFS (C), with G*(s) the interarrival LST, rho = 1/(E_Y mu):
 *   E[A]     = 1/mu + E[G^2]/(2E[G])
 *              + rho (-G*'(mu) + mu G*(mu) G*''(mu)/(1 + mu G*'(mu)))      (eq. 69)
 *   E[Apeak] = P0 (E[G] + (1+G*(mu))/mu) + Pw (E[G]/(1-G*(mu)) + 1/mu),
 *              P0 = q/(q+G*(mu)), Pw = 1-P0, q = 1 + mu G*'(mu)/(1-G*(mu))  (eq. 91)
 *   A*(s)    = eq. (65)
 * Discarding keeps the system stable for every rho, so rho &gt;= 1 is accepted.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfsd_gim1 {
    private Aoi_lcfsd_gim1() {}

    /**
     * Mean AoI, LST and peak AoI for the GI/M/1 non-preemptive LCFS-D queue.
     *
     * @param Y_lst LST of interarrival time distribution
     * @param mu Service rate (exponential service), must be positive
     * @param E_Y Mean interarrival time (first moment), must be positive
     * @param E_Y2 Second moment of interarrival time, must be &gt;= E_Y^2
     * @return AoiLstResult containing meanAoI, lstAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid
     */
    public static AoiLstResult aoi_lcfsd_gim1(final LstFunction Y_lst, final double mu, final double E_Y, double E_Y2) {
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");
        if (!(E_Y > 0)) throw new IllegalArgumentException("Mean interarrival time E_Y must be positive");
        if (!(E_Y2 >= E_Y * E_Y)) throw new IllegalArgumentException("Second moment E_Y2 must be >= E_Y^2");
        if (Y_lst == null) throw new IllegalArgumentException("The interarrival LST Y_lst is required");

        double lambda = 1.0 / E_Y;
        final double rho = lambda / mu;

        final double gM = Y_lst.evaluate(mu);
        double hstep = 1e-6 * Math.max(1.0, mu);
        final double mdY = -(Y_lst.evaluate(mu + hstep) - Y_lst.evaluate(mu - hstep)) / (2.0 * hstep);
        // five-point stencil, O(h^4): truncation and roundoff both near 1e-10 at this step
        double h2 = 5e-3 * Math.max(1.0, mu);
        double d2Y = (-Y_lst.evaluate(mu + 2 * h2) + 16.0 * Y_lst.evaluate(mu + h2) - 30.0 * gM
                + 16.0 * Y_lst.evaluate(mu - h2) - Y_lst.evaluate(mu - 2 * h2)) / (12.0 * h2 * h2);

        double meanAoI = 1.0 / mu + E_Y2 / (2.0 * E_Y) + rho * (mdY + mu * gM * d2Y / (1.0 - mu * mdY));

        // Two-state chain "no wait"/"wait" of informative updates (eq. 91)
        double q = 1.0 - mu * mdY / (1.0 - gM);
        double P0 = q / (q + gM);
        double Pw = gM / (q + gM);
        double peakAoI = P0 * (E_Y + (1.0 + gM) / mu) + Pw * (E_Y / (1.0 - gM) + 1.0 / mu);

        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                if (Math.abs(s) < 1e-12) return 1.0;  // removable singularity of the residual-interarrival LST
                double ys = s + mu;
                double hs = 1e-6 * Math.max(1.0, Math.abs(ys));
                double mdYs = -(Y_lst.evaluate(ys + hs) - Y_lst.evaluate(ys - hs)) / (2.0 * hs);
                double gres = (1.0 - Y_lst.evaluate(s)) / (s * E_Y);
                return (gres + rho * mu / ys * (Y_lst.evaluate(ys)
                        - gM * (1.0 - mu * mdYs) / (1.0 - mu * mdY))) * mu / ys;
            }
        };
        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }
}
