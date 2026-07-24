/**
 * @file Split-Merge maximum throughput for Fork-Join systems
 *
 * Computes the maximum throughput for a K-way Split-Merge (SM) queueing system.
 * In an SM system, all tasks of a request must complete before the next request
 * can be issued, making the maximum throughput inversely proportional to the
 * expected maximum service time.
 *
 * @since LINE 3.0
 */
package jline.api.fj;

public final class FJ_sm_tput {
    private FJ_sm_tput() {}

    /**
     * Maximum throughput for Split-Merge queueing system.
     *
     *   lambda_K^SM = 1 / X_K^max = mu / H_K
     *
     * For comparison, the maximum throughput of a K-way F/J system is
     * lambda_K^FJ = mu, which exceeds lambda_K^SM by a factor H_K.
     *
     * @param K  Number of parallel servers (positive integer)
     * @param mu Service rate at each server
     * @return Maximum throughput lambda_K^SM = mu / H_K
     *
     * @throws IllegalArgumentException if K &lt; 1 or mu &lt;= 0
     *
     * Reference: Thomasian, page 17:2.
     */
    public static double fj_sm_tput(int K, double mu) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        if (mu <= 0.0) {
            throw new IllegalArgumentException("Service rate mu must be positive. Got mu="
                    + String.format("%.4f", mu) + ".");
        }

        return mu / FJ_harmonic.fj_harmonic(K);
    }
}
