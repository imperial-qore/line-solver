package jline.api.qsys;

import java.util.HashMap;

import jline.io.Ret;

public final class Qsys_gigk_approx_whitt {
    private Qsys_gigk_approx_whitt() {}

    /**
     * GI/G/k approximation of Whitt (1993), eqs. (2.16)-(2.25):
     * Wq = phi(rho,ca2,cs2,k) * ((ca2+cs2)/2) * Wq(M/M/k),
     * where phi interpolates the Cosmetatos M/D/k (phi1) and D/M/k (phi3)
     * correction factors. Exact for M/M/k; reduces to the Cosmetatos M/D/k
     * approximation for cs2=0. Implements eq. (2.25) as printed, which was
     * validated against the paper's Tables 5-7 (New column).
     *
     * <p>Reference: Whitt, W. (1993). Approximations for the GI/G/m queue.
     * Production and Operations Management 2(2), 114-161.
     *
     * @param lambda arrival rate
     * @param mu service rate per server
     * @param ca2 squared coefficient of variation of inter-arrival time
     * @param cs2 squared coefficient of variation of service time
     * @param k number of servers
     * @return HashMap containing L, W, Q, U
     */
    public static HashMap<String, Object> qsys_gigk_approx_whitt(
            double lambda, double mu, double ca2, double cs2, int k) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / (k * mu);

        // Exact M/M/k baseline (Erlang-C based)
        Qsys_mmk.qsys_mmk(lambda, mu, k);
        double W_mmk = Ret.qsys.W;
        double Wq_mmk = W_mmk - 1.0 / mu;

        // Cosmetatos correction, as modified by Whitt (1993), eq. (2.17)
        double gamma = Math.min(0.24, (1.0 - rho) * (k - 1) * (Math.sqrt(4.0 + 5.0 * k) - 2.0) / (16.0 * k * rho));
        double phi1 = 1.0 + gamma;                             // M/D/k factor, eq. (2.16)
        double phi2 = 1.0 - 4.0 * gamma;                       // eq. (2.18)
        double phi3 = phi2 * Math.exp(-2.0 * (1.0 - rho) / (3.0 * rho)); // D/M/k factor, eq. (2.20)
        double phi4 = Math.min(1.0, (phi1 + phi3) / 2.0);      // eq. (2.21)

        double c2 = (ca2 + cs2) / 2.0;
        double psi;
        if (c2 >= 1.0) {
            psi = 1.0;                                          // eq. (2.22)
        } else {
            psi = Math.pow(phi4, 2.0 * (1.0 - c2));
        }

        double phi;
        if (Math.abs(ca2 - cs2) < 1e-12) {
            phi = psi;                                          // eq. (2.25) reduces to psi
        } else if (ca2 > cs2) {
            phi = (4.0 * (ca2 - cs2) / (4.0 * ca2 - 3.0 * cs2)) * phi1 + (cs2 / (4.0 * ca2 - 3.0 * cs2)) * psi;
        } else {
            phi = ((cs2 - ca2) / (2.0 * (ca2 + cs2))) * phi3 + ((cs2 + 3.0 * ca2) / (2.0 * (ca2 + cs2))) * psi;
        }

        double Wq = phi * c2 * Wq_mmk;                          // eq. (2.24)
        double W = Wq + 1.0 / mu;
        double L = lambda * W;
        double Q = lambda * Wq;
        double U = rho;

        result.put("L", L);
        result.put("W", W);
        result.put("Q", Q);
        result.put("U", U);

        return result;
    }
}
