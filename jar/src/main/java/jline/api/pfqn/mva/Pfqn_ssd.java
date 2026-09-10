/**
 * @file Server-Station Disaggregation throughput bounds (Suri-Dallery 1986)
 *
 * O(K) throughput bounds for single-class closed networks with multiserver stations,
 * bracketing each C_k-server station between single-server disaggregations jointly with
 * ABA (Theorem 5). With Z&gt;0 the queueing terms carry the terminal-workload correction
 * of Lazowska et al. 1984, Table 5.2. Ported at parity from MATLAB pfqn_ssd.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_ssd {
    private Pfqn_ssd() {}

    /**
     * SSD multiserver throughput bounds.
     *
     * @param L        service demand vector (M x 1)
     * @param N        population
     * @param Z        think time
     * @param nservers per-station server counts C_k (M x 1)
     * @return {Xlo, Xhi} lower/upper throughput bounds
     */
    public static double[] pfqn_ssd(Matrix L, double N, double Z, Matrix nservers) {
        int K = L.length();
        double Rl = 0.0, Yl = 0.0, Ru = 0.0;
        int b = 0;
        double bestLc = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < K; k++) {
            double Ck = nservers.get(k);
            double Lc = L.get(k) / Ck;
            Rl += L.get(k);
            Ru += Lc;
            if (Lc > Yl) Yl = Lc;
            if (Lc > bestLc) { bestLc = Lc; b = k; }
        }
        double Yu = Ru / K;
        double Xlo = N / (Rl + Z + (N - 1) * Yl / (1 + Z / (N * Rl)));
        double xhi = N / (Ru + Z + (N - 1) * Yu / (1 + Z / Ru));  // Theorem 5 upper
        double aba1 = nservers.get(b) / L.get(b);         // ABA capacity bound
        double aba2 = N / (Rl + Z);                       // ABA population bound
        double Xhi = Math.min(xhi, Math.min(aba1, aba2));
        return new double[]{Xlo, Xhi};
    }
}
