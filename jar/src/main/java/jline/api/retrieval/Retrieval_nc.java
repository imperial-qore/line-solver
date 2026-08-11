/**
 * @file Retrieval_nc.java
 * @brief Exact normalizing constant of a delayed-hit (list-based) cache.
 *
 * Port of matlab/src/api/retrieval/retrieval_nc.m. Computes the delayed-hit
 * cache normalizing constant E(v,m) by the exact recurrence (paper eq. mainrec):
 *
 *   E(v,m) = (1 + lambda_k*eta_{0,k}) E_k(v,m)
 *          + sum_{s=1}^r lambda_k*eta_{s,k}*(v_s+1) E_k(v+1_s,m)
 *          + sum_{j=1}^h m_j*gamma_{k,j} E_k(v,m-1_j)
 *
 * The plain constant is E(m) = E(0,m).
 *
 * @since LINE 3.0
 */
package jline.api.retrieval;

public final class Retrieval_nc {
    private Retrieval_nc() {}

    /**
     * Normalizing constant E(v,m).
     *
     * @param v     moment-order vector for the PS stations (zeros(1,r) for the plain constant)
     * @param m     cache list capacities (length h)
     * @param lambda per-item arrival rates (length n)
     * @param eta   fetching demands eta[i][s], column 0 = IS aggregate, columns 1..r = PS stations (n x (r+1))
     * @param gamma access factors gamma[i][j] (n x h)
     * @return E(v,m)
     */
    public static double retrieval_nc(double[] v, double[] m, double[] lambda, double[][] eta, double[][] gamma) {
        return sub(v.clone(), m.clone(), lambda, eta, gamma, lambda.length);
    }

    private static double sub(double[] v, double[] m, double[] lambda, double[][] eta, double[][] gamma, int k) {
        int r = v.length;
        int h = m.length;
        double msum = 0;
        double mmin = Double.POSITIVE_INFINITY;
        for (int j = 0; j < h; j++) { msum += m[j]; if (m[j] < mmin) mmin = m[j]; }
        if (msum > k || mmin < 0) {
            return 0.0;
        }
        if (k == 0) {
            return 1.0;
        }
        // item k (1-indexed) -> lambda[k-1], eta[k-1][.], gamma[k-1][.]
        int ki = k - 1;
        double E = (1.0 + lambda[ki] * eta[ki][0]) * sub(v, m, lambda, eta, gamma, k - 1);

        // item k fetched at PS station s = 1..r
        for (int s = 1; s <= r; s++) {
            double[] vp = v.clone();
            vp[s - 1] = vp[s - 1] + 1;
            E += lambda[ki] * eta[ki][s] * (v[s - 1] + 1) * sub(vp, m, lambda, eta, gamma, k - 1);
        }

        // item k stored in cache list j = 1..h
        for (int j = 0; j < h; j++) {
            if (m[j] > 0) {
                double[] mp = m.clone();
                mp[j] = mp[j] - 1;
                E += gamma[ki][j] * m[j] * sub(v, mp, lambda, eta, gamma, k - 1);
            }
        }
        return E;
    }
}
