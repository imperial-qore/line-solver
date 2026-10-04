/**
 * @file GI/M/1 non-preemptive LCFS-S Age of Information analysis
 *
 * Computes the mean AoI, its LST and the mean peak AoI of a GI/M/1 queue with
 * non-preemptive Last-Come First-Served with Set-aside (LCFS-S), i.e. LCFS
 * without discarding: the freshest waiting update is served next, and older
 * ones stay in queue and are still served later.
 *
 * Exact, Inoue, Masuyama, Takine and Tanaka (IEEE Trans. IT 65(12), 2019),
 * Section 3.3, NP-LCFS (D), with G*(s) the interarrival LST, rho = 1/(E_Y mu)
 * and gamma the root of G*(mu - mu x) = x:
 *   E[A]     = 1/mu + E[G^2]/(2E[G]) + rho (-G*'(mu - mu gamma))            (eq. 71)
 *   E[Apeak] = P0 (E[G] + (1+G*(mu))/mu)
 *              + Pw (1/mu + (E[G] + G*'(mu) + (gamma-G*(mu))/(gamma mu))/(1-G*(mu))),
 *              P0 = (1-gamma)/(1-gamma G*(mu)), Pw = 1-P0              (eqs. 100-103)
 *   A*(s)    = eq. (67)
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

public final class Aoi_lcfss_gim1 {
    private Aoi_lcfss_gim1() {}

    /**
     * Mean AoI, LST and peak AoI for the GI/M/1 non-preemptive LCFS-S queue.
     *
     * @param Y_lst LST of interarrival time distribution
     * @param mu Service rate (exponential service), must be positive
     * @param E_Y Mean interarrival time (first moment), must be positive
     * @param E_Y2 Second moment of interarrival time, must be &gt;= E_Y^2
     * @return AoiLstResult containing meanAoI, lstAoI, peakAoI
     * @throws IllegalArgumentException if parameters are invalid or system is unstable
     */
    public static AoiLstResult aoi_lcfss_gim1(final LstFunction Y_lst, final double mu, final double E_Y, double E_Y2) {
        if (!(mu > 0)) throw new IllegalArgumentException("Service rate mu must be positive");
        if (!(E_Y > 0)) throw new IllegalArgumentException("Mean interarrival time E_Y must be positive");
        if (!(E_Y2 >= E_Y * E_Y)) throw new IllegalArgumentException("Second moment E_Y2 must be >= E_Y^2");
        if (Y_lst == null) throw new IllegalArgumentException("The interarrival LST Y_lst is required");

        double lambda = 1.0 / E_Y;
        final double rho = lambda / mu;
        if (!(rho < 1)) {
            throw new IllegalArgumentException(String.format("System unstable: rho = 1/(E_Y*mu) = %.4f >= 1", rho));
        }

        final double gam = findGamma(Y_lst, mu);
        double gM = Y_lst.evaluate(mu);
        double hstep = 1e-6 * Math.max(1.0, mu);
        double mdY = -(Y_lst.evaluate(mu + hstep) - Y_lst.evaluate(mu - hstep)) / (2.0 * hstep);
        double zg = mu - mu * gam;
        double hz = 1e-6 * Math.max(1.0, zg);
        double mdYg = -(Y_lst.evaluate(zg + hz) - Y_lst.evaluate(zg - hz)) / (2.0 * hz);

        double meanAoI = 1.0 / mu + E_Y2 / (2.0 * E_Y) + rho * mdYg;

        double P0 = (1.0 - gam) / (1.0 - gam * gM);
        double Pw = gam * (1.0 - gM) / (1.0 - gam * gM);
        double E0 = E_Y + (1.0 + gM) / mu;
        double Ew = 1.0 / mu + (E_Y - mdY + (gam - gM) / (gam * mu)) / (1.0 - gM);
        double peakAoI = P0 * E0 + Pw * Ew;

        LstFunction lstAoI = new LstFunction() {
            @Override
            public double evaluate(double s) {
                if (Math.abs(s) < 1e-12) return 1.0;  // removable singularity of the residual-interarrival LST
                double gres = (1.0 - Y_lst.evaluate(s)) / (s * E_Y);
                return (gres + rho * (Y_lst.evaluate(s + mu - mu * gam) - gam) * mu / (s + mu)) * mu / (s + mu);
            }
        };
        return new AoiLstResult(meanAoI, lstAoI, peakAoI);
    }

    /**
     * gamma: root in (0,1) of Y*(mu - mu x) = x, the probability that an arrival
     * finds the server busy, by bisection on MATLAB's bracket [0.001, 0.999].
     */
    private static double findGamma(LstFunction Y_lst, double mu) {
        double lo = 0.001;
        double hi = 0.999;
        double fLo = Y_lst.evaluate(mu - mu * lo) - lo;
        double fHi = Y_lst.evaluate(mu - mu * hi) - hi;
        if (fLo * fHi > 0) {
            throw new IllegalArgumentException("aoi_lcfss_gim1: the root of Y*(mu - mu*x) = x is not bracketed by [0.001, 0.999]");
        }
        for (int iter = 0; iter < 200; iter++) {
            double mid = 0.5 * (lo + hi);
            if (mid == lo || mid == hi) break;
            double fMid = Y_lst.evaluate(mu - mu * mid) - mid;
            if (fMid == 0) return mid;
            if ((fMid > 0) == (fLo > 0)) { lo = mid; fLo = fMid; } else { hi = mid; }
        }
        return 0.5 * (lo + hi);
    }
}
