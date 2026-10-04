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
     * <p>Exact: the age is the backward recurrence time of the arrivals plus an
     * independent Exp(mu), so A*(s) = (mu/(s+mu))*lambda*(1 - Y*(s))/s and
     * E[A] = lambda*E[Y^2]/2 + 1/mu; E[Apeak] = E[S given success] + 1/(lambda*q)
     * with q = 1 - Y*(mu).</p>
     *
     * @param Y_lst LST of interarrival time distribution
     * @param mu Service rate (exponential service), must be positive
     * @param E_Y Mean interarrival time (first moment), must be positive
     * @param E_Y2 Second moment of interarrival time, must be at least E_Y^2
     * @return AoiLstResult containing meanAoI, lstAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_lcfspr_gim1(final LstFunction Y_lst, final double mu, double E_Y, double E_Y2) {
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");
        if (!(E_Y > 0)) throw new IllegalArgumentException("Mean interarrival time E_Y must be positive");
        if (E_Y2 < E_Y * E_Y) throw new IllegalArgumentException("Second moment E_Y2 must be >= E_Y^2");

        double lambda = 1.0 / E_Y;
        double rho = lambda / mu;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = 1/(E_Y*mu) = %.4f >= 1", rho));
        }

        // Mean AoI for preemptive LCFS: mean backward recurrence time plus 1/mu
        double meanAoI = lambda * E_Y2 / 2.0 + 1.0 / mu;

        // Mean Peak AoI (exact, as MATLAB): q = 1 - Y*(mu),
        // E[S|success] = (1/mu + Y*'(mu) - Y*(mu)/mu)/q
        double hstep = 1e-6 * Math.max(1.0, mu);
        double yMu = Y_lst.evaluate(mu);
        double dYstar = (Y_lst.evaluate(mu + hstep) - Y_lst.evaluate(mu - hstep)) / (2.0 * hstep);
        double q = 1.0 - yMu;
        double esSucc = (1.0 / mu + dYstar - yMu / mu) / q;
        double peakAoI = esSucc + 1.0 / (lambda * q);

        // LST of AoI for preemptive LCFS
        // A*(s) = (mu/(s+mu)) * lambda*(1 - Y*(s))/s, whose removable singularity at s = 0 is 1
        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                if (Math.abs(s) < 1e-12) return 1.0;
                return (mu / (s + mu)) * lambda * (1.0 - Y_lst.evaluate(s)) / s;
            }
        };

        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }
}
