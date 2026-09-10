/**
 * @file Multi-server load-dependent scaling factor computation
 *
 * Computes load-dependent scaling factors for multi-server queueing stations with finite
 * server capacity. Implements recursive algorithms for computing state-dependent service rates
 * in multi-server environments with load-dependent service mechanisms.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import org.apache.commons.math3.util.FastMath;

import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_mu_ms {
    private Pfqn_mu_ms() {}

    public static Matrix pfqn_mu_ms(int N, int m, int c) {
        Matrix mu = Matrix.zeros(1, N);
        Matrix g = Matrix.zeros(m, N + 1); // table

        for (int n = 0; n <= N; n++) {
            for (int i = 0; i < m; i++) {
                g.set(i, n, pfqn_mu_ms_gnaux(n, i, c, g));
            }
        }

        for (int n = 1; n <= N; n++) {
            mu.set(0, n - 1, g.get(m - 1, n - 1) / g.get(m - 1, n));
        }

        return mu;
    }

    public static double pfqn_mu_ms_gnaux(int n, int m, int c, Matrix g) {
        if (n == 0) {
            return 1.0;
        } else {
            if (m == 0) {
                double prodmin = 1.0;
                for (int i = 1; i <= n; i++) {
                    prodmin *= FastMath.min(i, c);
                }
                return 1.0 / prodmin;
            } else {
                double gn = 0.0;
                for (int k = 0; k <= n; k++) {
                    double a;
                    if (k > 0) {
                        Matrix avec = new Matrix(1, k);
                        for (int t = 0; t < k; t++) {
                            avec.set(0, t, Maths.min((double) (1 + t), (double) c));
                        }
                        a = avec.elementMult();
                    } else {
                        a = 1.0;
                    }
                    double b = 1 / g.get(m - 1, n - k);
                    gn += 1.0 / (a * b);
                }
                return gn;
            }
        }
    }
}
