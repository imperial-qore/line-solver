/**
 * @file GI/M/1 preemptive LCFS Age of Information analysis
 *
 * Computes mean AoI, LST, and peak AoI for a GI/M/1 queue with
 * preemptive Last-Come First-Served (LCFS-PR) discipline.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfspr_gim1 {
    private Aoi_lcfspr_gim1() {}

    /**
     * Mean AoI and LST for GI/M/1 preemptive LCFS queue.
     *
     * <p>Exact peak: E[Apeak] = E[S given success] + 1/(lambda*q) with
     * q = 1 - Y*(mu).</p>
     *
     * @param Y_lst LST of interarrival time distribution
     * @param mu Service rate (exponential service), must be positive
     * @param E_Y Mean interarrival time (first moment), must be positive
     * @param E_Y2 Second moment of interarrival time (for reference)
     * @return AoiLstResult containing meanAoI, lstAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_lcfspr_gim1(final LstFunction Y_lst, final double mu, double E_Y, double E_Y2) {
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");
        if (!(E_Y > 0)) throw new IllegalArgumentException("Mean interarrival time E_Y must be positive");

        double lambda = 1.0 / E_Y;
        double rho = lambda / mu;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = 1/(E_Y*mu) = %.4f >= 1", rho));
        }

        // Mean service time
        double E_S = 1.0 / mu;

        // Mean AoI for preemptive LCFS (Proposition 4)
        double meanAoI = E_Y + E_S;

        // Mean Peak AoI (exact, as MATLAB): q = 1 - Y*(mu),
        // E[S|success] = (1/mu + Y*'(mu) - Y*(mu)/mu)/q
        double hstep = 1e-6 * Math.max(1.0, mu);
        double yMu = Y_lst.evaluate(mu);
        double dYstar = (Y_lst.evaluate(mu + hstep) - Y_lst.evaluate(mu - hstep)) / (2.0 * hstep);
        double q = 1.0 - yMu;
        double esSucc = (1.0 / mu + dYstar - yMu / mu) / q;
        double peakAoI = esSucc + 1.0 / (lambda * q);

        // LST of AoI for preemptive LCFS
        // A*(s) = Y*(s) * (mu / (s + mu))
        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                return Y_lst.evaluate(s) * (mu / (s + mu));
            }
        };

        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }
}
