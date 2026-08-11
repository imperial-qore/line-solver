/**
 * @file Upper Geometric Bound computation for queue lengths in single-class networks
 *
 * Computes the upper Geometric Bound (GB) for queue lengths in closed single-class queueing
 * networks. Provides tight upper bounds for performance analysis and bounding techniques
 * in approximate MVA algorithms.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_qzgbup {
    private Pfqn_qzgbup() {}

    /**
     * Computes the upper Geometric Bound (GB) for the queue length of the given closed
     * single-class queueing networks.
     *
     * @param L service demand matrix
     * @param N population
     * @param Z think time
     * @param i station index
     * @return the upper GB for the queue length
     */
    public static double pfqn_qzgbup(Matrix L, double N, double Z, int i) {
        double sumL = L.elementSum();
        double sumLSq = L.elementMult(L, null).elementSum();
        double sigma = sumLSq / sumL;
        double Yi = L.get(i) * Maths.min(1.0 / L.elementMax(),
                N / (Z + sumL + sigma * (N - 1 - Z * Pfqn_xzabaup.pfqn_xzabaup(L, N - 1, Z))));
        double Qgb;
        if (Yi < 1) {
            Qgb = Yi / (1 - Yi) - (Math.pow(Yi, N + 1) / (1 - Yi));
        } else {
            Qgb = N;
        }
        return Qgb;
    }
}
