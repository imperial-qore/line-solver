/**
 * @file Random Replacement Model Mean Field ODE System
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;

public final class Cache_rrm_meanfield_ode {
    private Cache_rrm_meanfield_ode() {}

    /**
     * ODE system for RRM (Random Replacement Model) mean field equations.
     *
     * @param x      Current state matrix (n x (1+h)) where x[k,s] is probability of item k being in list s
     * @param lambda Request rates vector
     * @param m      Cache sizes per level
     * @param n      Number of items
     * @param h      Number of cache levels
     * @return Time derivative matrix dxdt
     */
    public static Matrix cache_rrm_meanfield_ode(Matrix x, Matrix lambda, Matrix m, int n, int h) {
        Matrix dxdt = Matrix.zeros(n, 1 + h);

        for (int k = 0; k < n; k++) {
            for (int s = 1; s <= h; s++) {
                int sIndex = s;

                // First term: promotion from list s-1
                double sum1 = 0.0;
                for (int k1 = 0; k1 < n; k1++) {
                    sum1 += lambda.get(k1) / m.get(s - 1) * x.get(k1, sIndex - 1) * x.get(k, sIndex);
                }

                // Second term: demotion from list s+1
                double sum2 = 0.0;
                if (s < h) {
                    for (int k1 = 0; k1 < n; k1++) {
                        sum2 += lambda.get(k1) / m.get(s) * x.get(k1, sIndex) * x.get(k, sIndex + 1);
                    }
                    sum2 -= lambda.get(k) * x.get(k, sIndex);
                }

                // Drift component
                dxdt.set(k, sIndex, lambda.get(k) * x.get(k, sIndex - 1) - sum1 + sum2);
            }

            // Case s=0: conservation law
            double sumDerivatives = 0.0;
            for (int s = 1; s <= h; s++) {
                sumDerivatives += dxdt.get(k, s);
            }
            dxdt.set(k, 0, -sumDerivatives);
        }

        return dxdt;
    }
}
