/**
 * @file M/GI/1 preemptive LCFS Age of Information analysis
 *
 * Computes mean AoI, LST, and peak AoI for an M/GI/1 queue with
 * preemptive Last-Come First-Served (LCFS-PR) discipline.
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfspr_mgi1 {
    private Aoi_lcfspr_mgi1() {}

    /**
     * Mean AoI and LST for M/GI/1 preemptive LCFS queue.
     *
     * <p>Exact peak: E[Apeak] = -H*'(lambda)/H*(lambda)
     * + 1/(lambda*H*(lambda)).</p>
     *
     * @param lambda Arrival rate (Poisson arrivals), must be positive
     * @param H_lst LST of service time distribution
     * @param E_H Mean service time (first moment), must be positive
     * @param E_H2 Second moment of service time (for reference)
     * @return AoiLstResult containing meanAoI, lstAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_lcfspr_mgi1(final double lambda, final LstFunction H_lst, double E_H, double E_H2) {
        if (!(lambda > 0)) throw new IllegalArgumentException("Arrival rate lambda must be positive");
        if (!(E_H > 0)) throw new IllegalArgumentException("Mean service time E_H must be positive");

        double rho = lambda * E_H;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = lambda*E_H = %.4f >= 1", rho));
        }

        // Mean interarrival time
        double E_Y = 1.0 / lambda;

        // Mean AoI for preemptive LCFS (Proposition 3)
        double meanAoI = E_Y + E_H;

        // Mean Peak AoI (exact, as MATLAB): q = H*(lambda),
        // E[S|success] = -H*'(lambda)/H*(lambda)
        double hstep = 1e-6 * Math.max(1.0, lambda);
        double hLam = H_lst.evaluate(lambda);
        double dHstar = (H_lst.evaluate(lambda + hstep) - H_lst.evaluate(lambda - hstep)) / (2.0 * hstep);
        double peakAoI = -dHstar / hLam + 1.0 / (lambda * hLam);

        // LST of AoI for preemptive LCFS
        // A*(s) = (lambda / (s + lambda)) * H*(s)
        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                return (lambda / (s + lambda)) * H_lst.evaluate(s);
            }
        };

        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }
}
