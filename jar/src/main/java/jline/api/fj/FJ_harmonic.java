/**
 * @file Harmonic number computation for Fork-Join analysis
 *
 * Computes the K-th Harmonic number H_K = sum_{k=1}^{K} 1/k, which is a fundamental
 * quantity in Fork-Join queueing analysis. For large K, H_K ~ ln(K) + gamma,
 * where gamma ~ 0.57721 is the Euler-Mascheroni constant.
 *
 * @since LINE 3.0
 */
package jline.api.fj;

public final class FJ_harmonic {
    private FJ_harmonic() {}

    /**
     * Compute Harmonic sum H_K = sum(1/k) for k=1 to K.
     *
     * @param K Number of parallel servers (positive integer)
     * @return Harmonic sum H_K = 1 + 1/2 + 1/3 + ... + 1/K
     *
     * @throws IllegalArgumentException if K &lt; 1
     *
     * Reference: A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
     * ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
     */
    public static double fj_harmonic(int K) {
        if (K < 1) {
            throw new IllegalArgumentException("K must be a positive integer. Got K=" + K + ".");
        }
        double H = 0.0;
        for (int k = 1; k <= K; k++) {
            H += 1.0 / k;
        }
        return H;
    }
}
