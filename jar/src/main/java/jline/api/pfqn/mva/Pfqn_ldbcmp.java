/**
 * @file Anselmi-Cremonesi (2008) lower throughput bound for closed single-class BCMP
 *       networks with load-dependent stations
 *
 * Lower bound on system throughput via the asymptotic closed-open equivalence (their eq. 15),
 * refined to a monotone fixed point by Algorithm 1. Applicable when N >= Qhat. Ported at
 * parity from MATLAB pfqn_ldbcmp.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_ldbcmp {
    private Pfqn_ldbcmp() {}

    /**
     * Lower throughput bound (and response-time upper bound) for LD-BCMP closed networks.
     *
     * @param L service demand vector (M x 1); limiting demand for Heffes LD stations
     * @param N population
     * @param Z think time (non-bottleneck delay station)
     * @param c per-station Heffes coefficient (M x 1); c(i)=0 fixed-rate, c(i)>0 Heffes LD
     * @return {Xlo, Rhi, Qhat}; Xlo=Rhi=NaN if N < Qhat
     */
    public static double[] pfqn_ldbcmp(Matrix L, double N, double Z, Matrix c) {
        int M = L.length();
        double Dm = L.elementMax();
        int bmax = 0;
        for (int i = 0; i < M; i++) if (Math.abs(L.get(i) - Dm) <= 1e-12 * Dm) bmax++;
        double lambda = 1.0 / Dm;
        double Qhat = 0.0;
        for (int i = 0; i < M; i++) {
            if (Math.abs(L.get(i) - Dm) <= 1e-12 * Dm) continue;    // bottleneck
            double rho = lambda * L.get(i);
            if (rho >= 1) return new double[]{Double.NaN, Double.NaN, Double.NaN};
            double ci = (c != null && c.length() > i) ? c.get(i) : 0.0;
            Qhat += (ci + 1) * rho / (1 - rho);
        }
        Qhat += lambda * Z;                                          // delay station
        if (N - Qhat < 0) return new double[]{Double.NaN, Double.NaN, Qhat};

        double a = N - Qhat;
        double Xprime = 0.0, Xlo = 0.0;
        double tol = 1e-10;
        for (int it = 0; it < 10000; it++) {
            double Xprev = Xlo;
            double denom = Dm * (bmax + N - Qhat) - bmax * Math.pow(Dm * Xprime, N) * Dm;
            Xlo = a / denom;
            Xprime = Xlo;
            if (Xprev > 0 && Math.abs(Xprev - Xlo) / Xprev <= tol) break;
        }
        return new double[]{Xlo, N / Xlo, Qhat};
    }
}
