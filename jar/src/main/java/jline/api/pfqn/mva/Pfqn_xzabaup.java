/**
 * @file Upper Asymptotic Bound Approximation (ABA) for throughput in single-class networks
 *
 * Computes the upper ABA bound for throughput in closed single-class queueing networks.
 * Provides fundamental upper bounds based on bottleneck analysis and system capacity
 * constraints for performance evaluation.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_xzabaup {
    private Pfqn_xzabaup() {}

    /**
     * Computes the upper ABA for the throughput of the given closed single-class queueing networks
     *
     * @param L - service demand matrix
     * @param N - population
     * @param Z - think time
     * @return - the upper ABA for the throughput
     */
    public static double pfqn_xzabaup(Matrix L, double N, double Z) {
        double e1 = 1 / L.elementMax();
        double e2 = N / (L.elementSum() + Z);
        return Maths.min(e1, e2);
    }
}
