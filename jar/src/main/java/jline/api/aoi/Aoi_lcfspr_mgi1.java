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
     * <p>Exact: A*(s) = lambda*H*(s+lambda)/(s + lambda*H*(s+lambda)),
     * E[A] = 1/(lambda*H*(lambda)), E[Apeak] = -H*'(lambda)/H*(lambda)
     * + 1/(lambda*H*(lambda)). The age looks back over Exp(lambda) gaps to the
     * first update whose service beat the next arrival.</p>
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

        // Mean AoI for preemptive LCFS: the delivery rate is lambda*H*(lambda)
        double meanAoI = 1.0 / (lambda * H_lst.evaluate(lambda));

        // Mean Peak AoI (exact, as MATLAB): q = H*(lambda),
        // E[S|success] = -H*'(lambda)/H*(lambda)
        double hstep = 1e-6 * Math.max(1.0, lambda);
        double hLam = H_lst.evaluate(lambda);
        double dHstar = (H_lst.evaluate(lambda + hstep) - H_lst.evaluate(lambda - hstep)) / (2.0 * hstep);
        double peakAoI = -dHstar / hLam + 1.0 / (lambda * hLam);

        // LST of AoI for preemptive LCFS
        // A*(s) = lambda*H*(s+lambda) / (s + lambda*H*(s+lambda)), A*(0) = 1 identically
        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                double f = lambda * H_lst.evaluate(s + lambda);
                return f / (s + f);
            }
        };

        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }
}
