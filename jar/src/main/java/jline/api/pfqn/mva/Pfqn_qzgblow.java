/**
 * @file Lower Geometric Bound computation for queue lengths in single-class networks
 *
 * Computes the lower Geometric Bound (GB) for queue lengths in closed single-class queueing
 * networks. Provides tight lower bounds for performance analysis and serves as input for
 * other approximation algorithms requiring queue length bounds.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_qzgblow {
    private Pfqn_qzgblow() {}

    /**
     * Computes the lower Geometric Bound (GB) for the queue length of the given closed single-class
     * queueing networks
     *
     * @param L - service demand matrix
     * @param N - population
     * @param Z - think time
     * @param i - station index
     * @return - the lower GB for the queue length
     */
    public static double pfqn_qzgblow(Matrix L, double N, double Z, int i) {
        double Yi = N * L.get(i) / (Z + L.elementSum() + L.elementMax() * N);
        double Qgb = Yi / (1 - Yi) - (Math.pow(Yi, N + 1) / (1 - Yi));
        return Qgb;
    }
}
