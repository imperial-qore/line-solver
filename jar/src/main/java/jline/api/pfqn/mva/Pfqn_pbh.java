/**
 * @file Performance Bound Hierarchy (Eager-Sevcik 1983) for single-class closed networks
 *
 * Level-parameterized nested optimistic/pessimistic throughput bounds that converge to
 * exact MVA as level -> N. Ported at parity from MATLAB pfqn_pbh.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Pfqn_pbh {
    private Pfqn_pbh() {}

    /**
     * Level-`level` Performance Bound Hierarchy throughput bounds.
     *
     * @param L     service demand vector (M x 1)
     * @param N     population
     * @param Z     think time
     * @param level hierarchy level (>= 0), clamped to N
     * @return {Xlo, Xhi} lower/upper throughput bounds
     */
    public static double[] pfqn_pbh(Matrix L, double N, double Z, int level) {
        int K = L.length();
        double Lmax = L.elementMax();
        double[] Ro = pbhResidence(L, N, Z, level, true);   // optimistic
        double[] Rp = pbhResidence(L, N, Z, level, false);  // pessimistic
        double sumRo = 0.0, sumRp = 0.0, sumL = L.elementSum();
        for (int k = 0; k < K; k++) { sumRo += Ro[k]; sumRp += Rp[k]; }
        double RoC = Maths.max(sumRo, Maths.max(N * Lmax - Z, sumL));
        double Xhi = Maths.min(1.0 / Lmax, N / (Z + RoC));
        double Xlo = N / (Z + sumRp);
        return new double[]{Xlo, Xhi};
    }

    /** Per-station residence-time vector of the level-`level` PBH bound. */
    private static double[] pbhResidence(Matrix L, double N, double Z, int level, boolean opt) {
        int K = L.length();
        int b = 0;                        // bottleneck index (argmax L)
        for (int k = 1; k < K; k++) if (L.get(k) > L.get(b)) b = k;
        double Lb = L.get(b), sumL = L.elementSum();
        if (level > N) level = (int) N;
        double n0 = N - level;
        double[] Rk = new double[K];
        if (!opt) {
            // Pessimistic start, carried in QUEUE LENGTHS: all n0 customers at
            // the bottleneck (eq 13). Seeding a residence instead makes the
            // assumed population n0*Rb/(Z+Rtot) < n0 once Z > 0, so the
            // pessimism is diluted and Xlo stops being a bound (violated exact
            // on 14% of random delay models, worst 43%).
            double[] Q = new double[K];
            Q[b] = n0;
            for (int k = 0; k < K; k++) Rk[k] = L.get(k) * (1 + Q[k]);
            for (int n = (int) n0 + 1; n <= N; n++) {
                double Rtot = 0.0;
                for (int k = 0; k < K; k++) {
                    Rk[k] = L.get(k) * (1 + Q[k]);         // eq (5), MVA step
                    Rtot += Rk[k];
                }
                double X = n / (Z + Rtot);
                for (int k = 0; k < K; k++) Q[k] = X * Rk[k];
            }
            return Rk;
        }
        double v = Maths.max(n0 * Lb - Z, sumL) / K;       // eq (7), ABA opt
        for (int k = 0; k < K; k++) Rk[k] = v;
        if (n0 == 0) {
            for (int k = 0; k < K; k++) Rk[k] = 0.0;       // exact empty base
        }
        for (int n = (int) n0 + 1; n <= N; n++) {
            double Rtot = 0.0;
            for (int k = 0; k < K; k++) Rtot += Rk[k];
            if (n == 1 || Z + Rtot == 0) {
                for (int k = 0; k < K; k++) Rk[k] = L.get(k);   // empty-network base step
            } else {
                for (int k = 0; k < K; k++) {
                    Rk[k] = L.get(k) * (1 + (n - 1) * Rk[k] / (Z + Rtot));   // eq (5)
                }
            }
        }
        return Rk;
    }
}
