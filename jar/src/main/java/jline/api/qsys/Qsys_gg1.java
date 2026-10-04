/**
 * @file G/G/1 queueing system analysis
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import java.util.HashMap;

import jline.io.Ret;

public final class Qsys_gg1 {
    private Qsys_gg1() {}

    /**
     * G/G/1 queue analysis (general arrivals and service).
     *
     * <p>Uses exact methods for the special cases (M/M/1, M/G/1, G/M/1) and the
     * Allen-Cunneen approximation for the general case. In the G/M/1 case the
     * interarrival-time distribution is fitted from (lambda, ca2) by a two-moment
     * renewal process (H2 with balanced means for ca2&gt;1, mixed Erlang for
     * ca2&lt;1) and sigma is the root of sigma = A*(mu*(1-sigma)), with A* the
     * interarrival-time LST.
     *
     * @param lambda arrival rate
     * @param mu service rate
     * @param ca2 squared coefficient of variation of inter-arrival time
     * @param cs2 squared coefficient of variation of service time
     * @return HashMap containing L, Lq, W, Wq, p0
     */
    public static HashMap<String, Object> qsys_gg1(double lambda, double mu, double ca2, double cs2) {
        HashMap<String, Object> result = new HashMap<String, Object>();
        double rho = lambda / mu;
        double tolerance = 1e-8;

        double W;
        if (Math.abs(ca2 - 1.0) < tolerance && Math.abs(cs2 - 1.0) < tolerance) {
            // M/M/1
            Qsys_mm1.qsys_mm1(lambda, mu);
            W = Ret.qsys.W;
        } else if (Math.abs(ca2 - 1.0) < tolerance) {
            // M/G/1 (qsys_mg1 expects the coefficient of variation, not the SCV)
            Qsys_mg1.qsys_mg1(lambda, mu, Math.sqrt(cs2));
            W = Ret.qsys.W;
        } else if (Math.abs(cs2 - 1.0) < tolerance) {
            // G/M/1: sigma from the two-moment fit of the interarrival-time LST
            double sigma = gm1SigmaTwoMoment(lambda, mu, ca2);
            Qsys_gm1.qsys_gm1(sigma, mu);
            W = Ret.qsys.W;
        } else {
            // General G/G/1: Allen-Cunneen (expects coefficients of variation)
            Qsys_gig1_approx_allencunneen.qsys_gig1_approx_allencunneen(lambda, mu, Math.sqrt(ca2), Math.sqrt(cs2));
            W = Ret.qsys.W;
        }

        double Wq = W - 1.0 / mu;
        double L = lambda * W;
        double Lq = lambda * Wq;
        double p0 = 1.0 - rho;
        result.put("L", L);
        result.put("Lq", Lq);
        result.put("W", W);
        result.put("Wq", Wq);
        result.put("p0", p0);
        return result;
    }

    /**
     * G/G/1 with state probability for k.
     */
    public static HashMap<String, Object> qsys_gg1(double lambda, double mu, double ca2, double cs2, int k) {
        HashMap<String, Object> result = qsys_gg1(lambda, mu, ca2, cs2);
        double rho = lambda / mu;
        double L = ((Double) result.get("L")).doubleValue();

        double pk;
        if (k == 0) {
            pk = ((Double) result.get("p0")).doubleValue();
        } else {
            double p0 = ((Double) result.get("p0")).doubleValue();
            double r = (L - rho) / L;
            pk = p0 * Math.pow(rho * (1.0 + r * (ca2 + cs2 - 2.0) / 2.0), (double) k);
        }

        result.put("pk", pk);
        return result;
    }

    /**
     * Root in (0,1) of sigma = A*(mu*(1-sigma)) for a two-moment fit of the
     * interarrival-time LST A*. The map T(x)=A*(mu*(1-x)) is increasing with
     * the queue root as its smallest fixed point, so fixed-point iterates
     * converge monotonically.
     */
    private static double gm1SigmaTwoMoment(double lambda, double mu, double ca2) {
        double sigma = lambda / mu;
        for (int it = 0; it < 100000; it++) {
            double signew = interarrivalLst(mu * (1.0 - sigma), lambda, ca2);
            if (Math.abs(signew - sigma) < 1e-13) {
                return signew;
            }
            sigma = signew;
        }
        return sigma;
    }

    private static double interarrivalLst(double s, double lambda, double ca2) {
        if (ca2 < 1e-6) {
            // deterministic interarrival times
            return Math.exp(-s / lambda);
        } else if (ca2 < 1.0) {
            // mixed Erlang(j-1,j) with common rate (Tijms, 1994)
            int j = (int) Math.ceil(1.0 / ca2);
            double p = (j * ca2 - Math.sqrt(j * (1.0 + ca2) - j * (double) j * ca2)) / (1.0 + ca2);
            double nu = (j - p) * lambda;
            return p * Math.pow(nu / (s + nu), j - 1) + (1.0 - p) * Math.pow(nu / (s + nu), j);
        } else {
            // hyperexponential H2 with balanced means
            double p1 = (1.0 + Math.sqrt((ca2 - 1.0) / (ca2 + 1.0))) / 2.0;
            double l1 = 2.0 * p1 * lambda;
            double l2 = 2.0 * (1.0 - p1) * lambda;
            return p1 * l1 / (s + l1) + (1.0 - p1) * l2 / (s + l2);
        }
    }
}
