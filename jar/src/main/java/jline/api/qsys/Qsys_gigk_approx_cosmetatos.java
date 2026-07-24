package jline.api.qsys;

import jline.io.Ret;

import java.util.HashMap;

public final class Qsys_gigk_approx_cosmetatos {
    private Qsys_gigk_approx_cosmetatos() {}

    /**
     * GI/G/k approximation by interpolation of the M/M/k, M/D/k and D/M/k
     * queues (Cosmetatos 1982; Page 1982):
     * Wq = [ca2*cs2 + ca2*(1-cs2)*phi1/2 + (1-ca2)*cs2*phi3/2]*Wq(M/M/k),
     * where phi1 and phi3 are the Cosmetatos (1975) correction factors for
     * M/D/k and D/M/k, with the safeguards of Whitt (1993). The D/D/k corner
     * has Wq=0. The interpolation requires ca2&lt;=1 and cs2&lt;=1; outside
     * this region the Lee-Longton scaling Wq = ((ca2+cs2)/2)*Wq(M/M/k) is used.
     *
     * <p>References: Cosmetatos, G.P. (1975). Approximate explicit formulae
     * for the average queueing time in the processes (M/D/r) and (D/M/r).
     * INFOR 13, 328-331. Page, E. (1982). Tables of waiting times for M/M/n,
     * M/D/n and D/M/n and their use to give approximate waiting times in more
     * general queues. J. Opl. Res. Soc. 33, 453-473.
     *
     * @param lambda arrival rate
     * @param mu service rate per server
     * @param ca2 squared coefficient of variation of inter-arrival time
     * @param cs2 squared coefficient of variation of service time
     * @param k number of servers
     * @return HashMap containing L, W, Q, U
     */
    public static HashMap<String, Object> qsys_gigk_approx_cosmetatos(
            double lambda, double mu, double ca2, double cs2, int k) {
        HashMap<String, Object> result = new HashMap<String, Object>();

        double rho = lambda / (k * mu);

        // Exact M/M/k baseline waiting time (Erlang-C based)
        Qsys_mmk.qsys_mmk(lambda, mu, k);
        double W_mmk = Ret.qsys.W;
        double Wq_mmk = W_mmk - 1.0 / mu;

        double Wq;
        if (ca2 <= 1.0 && cs2 <= 1.0) {
            // Cosmetatos correction, as modified by Whitt (1993), eq. (2.17)
            double gamma = Math.min(0.24, (1.0 - rho) * (k - 1) * (Math.sqrt(4.0 + 5.0 * k) - 2.0) / (16.0 * k * rho));
            double phi1 = 1.0 + gamma;                                        // M/D/k factor
            double phi3 = (1.0 - 4.0 * gamma) * Math.exp(-2.0 * (1.0 - rho) / (3.0 * rho)); // D/M/k factor
            Wq = (ca2 * cs2 + ca2 * (1.0 - cs2) * phi1 / 2.0 + (1.0 - ca2) * cs2 * phi3 / 2.0) * Wq_mmk;
        } else {
            // Interpolation weights are invalid outside the unit box
            Wq = ((ca2 + cs2) / 2.0) * Wq_mmk;
        }
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
