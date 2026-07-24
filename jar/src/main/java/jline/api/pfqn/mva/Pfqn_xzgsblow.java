/**
 * @file Lower Geometric Square-Root Bound (GSB) for throughput in single-class networks
 *
 * Computes the lower GSB for throughput in closed single-class queueing networks using
 * geometric mean analysis and square-root approximations. Provides tighter bounds than
 * basic asymptotic approximations through advanced geometric techniques.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Pfqn_xzgsblow {
    private Pfqn_xzgsblow() {}

    /**
     * Computes the lower Geometric Square-Root Bound (GSB) for the throughput of the given
     * closed single-class queueing networks.
     *
     * @param L service demand matrix
     * @param N population
     * @param Z think time
     * @return the lower GSB for the throughput
     */
    public static double pfqn_xzgsblow(Matrix L, double N, double Z) {
        int M = L.length();
        double maxL = L.elementMax();
        double R = Z + L.elementSum() + maxL * (N - 1);
        for (int i = 0; i < M; i++) {
            if (L.get(i) < maxL) {
                R += (L.get(i) - maxL) * Pfqn_qzgblow.pfqn_qzgblow(L, N - 1, Z, i);
            }
        }
        return 2 * N / (R + FastMath.sqrt(R * R - 4 * Z * maxL * (N - 1)));
    }
}
