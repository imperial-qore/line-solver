package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

/**
 * Lower Asymptotic Bound Approximation (ABA) for throughput in single-class networks.
 *
 * <p>Computes the lower ABA bound for throughput in closed single-class queueing networks.
 * Provides fundamental lower bounds based on asymptotic behavior analysis for performance
 * evaluation and approximation algorithm validation.
 *
 * @since LINE 3.0
 */
public final class Pfqn_xzabalow {
    private Pfqn_xzabalow() {}

    /**
     * Computes the lower ABA for the throughput of the given closed single-class queueing networks
     *
     * <p>This is the textbook ABA lower bound N / (Z + N * sum(L)), obtained by
     * assuming every job queues behind all the others. It is NOT the tighter
     * Zahorjan-Balanced bound; see {@code Pfqn_xzgsblow} for that.
     *
     * @param L - service demand matrix
     * @param N - population
     * @param Z - think time
     * @return - the lower ABA for the throughput
     */
    public static double pfqn_xzabalow(Matrix L, double N, double Z) {
        return N / (L.elementSum() * N + Z);
    }
}
