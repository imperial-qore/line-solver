/**
 * @file GI/M/1 non-preemptive LCFS-D Age of Information analysis
 *
 * Computes mean AoI for a GI/M/1 queue with non-preemptive Last-Come
 * First-Served with Discarding (LCFS-D) discipline.
 *
 * In LCFS-D, when a new update arrives while the server is busy:
 *   - If there's an update waiting in queue, it is discarded
 *   - The new update takes its place in the queue
 *   - When service completes, the waiting update (if any) is served
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfsd_gim1 {
    private Aoi_lcfsd_gim1() {}

    /**
     * Mean AoI for GI/M/1 non-preemptive LCFS-D queue.
     */
    public static AoiLstResult aoi_lcfsd_gim1(LstFunction Y_lst, double mu, double E_Y, double E_Y2) {
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
        double sigma = findSigmaGIM1ForLcfsd(Y_lst, mu, rho);

        // Mean effective system time
        double E_T_eff = E_S + sigma * E_S;

        // Mean AoI for LCFS-D (from Section VI analysis adapted for GI/M/1)
        double meanAoI = E_Y + E_S * (1.0 + sigma) + sigma * E_S / (1.0 + sigma);

        // Mean Peak AoI
        double peakAoI = E_Y + E_T_eff;

        // LST is complex for LCFS-D; return null
        return new AoiLstResult(meanAoI, null, peakAoI);
    }

    /**
     * Find sigma for GI/M/1 LCFS-D: root of Y*(mu - mu*sigma) = sigma.
     */
    private static double findSigmaGIM1ForLcfsd(LstFunction Y_lst, double mu, double rho) {
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
