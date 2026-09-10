/**
 * @file GI/M/1 non-preemptive LCFS-S Age of Information analysis
 *
 * Computes mean AoI for a GI/M/1 queue with non-preemptive Last-Come
 * First-Served with Set-aside (LCFS-S) discipline.
 *
 * In LCFS-S, when a new update arrives while the server is busy:
 *   - The new update waits in the queue
 *   - When service completes, the most recent update in queue is served next
 *   - The older update remains in queue (is "set aside")
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfss_gim1 {
    private Aoi_lcfss_gim1() {}

    /**
     * Mean AoI for GI/M/1 non-preemptive LCFS-S queue.
     *
     * @param Y_lst LST of interarrival time distribution
     * @param mu Service rate (exponential service), must be positive
     * @param E_Y Mean interarrival time (first moment), must be positive
     * @param E_Y2 Second moment of interarrival time
     * @return AoiLstResult containing meanAoI, lstAoI (null), peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_lcfss_gim1(LstFunction Y_lst, double mu, double E_Y, double E_Y2) {
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");
        if (!(E_Y > 0)) throw new IllegalArgumentException("Mean interarrival time E_Y must be positive");

        double lambda = 1.0 / E_Y;
        double rho = lambda / mu;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = 1/(E_Y*mu) = %.4f >= 1", rho));
        }

        // Mean service time
        double E_S = 1.0 / mu;

        // Find sigma: probability arriving customer finds server busy
        double sigma = findSigmaGIM1ForLcfss(Y_lst, mu, rho);

        // Mean system delay
        double E_D = 1.0 / (mu * (1.0 - sigma));

        // Mean AoI for LCFS-S (from Section V analysis adapted for GI/M/1)
        double meanAoI = E_Y + E_S + sigma * E_D;

        // Mean Peak AoI
        double peakAoI = E_Y + E_D;

        // LST is complex for LCFS-S; return null
        return new AoiLstResult(meanAoI, null, peakAoI);
    }

    /**
     * Find sigma for GI/M/1 LCFS-S: root of Y*(mu - mu*sigma) = sigma.
     */
    private static double findSigmaGIM1ForLcfss(LstFunction Y_lst, double mu, double rho) {
        double lo = 0.001;
        double hi = 0.999;
        int maxIter = 100;
        double tol = 1e-12;

        double fLo = Y_lst.evaluate(mu - mu * lo) - lo;
        double fHi = Y_lst.evaluate(mu - mu * hi) - hi;

        if (fLo * fHi < 0) {
            for (int iter = 0; iter < maxIter; iter++) {
                double mid = (lo + hi) / 2.0;
                double fMid = Y_lst.evaluate(mu - mu * mid) - mid;
                if (Math.abs(fMid) < tol) {
                    return mid;
                }
                double fLoNew = Y_lst.evaluate(mu - mu * lo) - lo;
                if (fLoNew * fMid < 0) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            return (lo + hi) / 2.0;
        } else {
            // Fallback
            return rho;
        }
    }
}
