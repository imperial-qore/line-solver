/**
 * @file Complex-valued single-class load-dependent auxiliary normalizing constant computation
 *
 * Provides specialized auxiliary function for computing normalizing constants in single-class
 * load-dependent closed queueing networks with complex-valued service demands. Implements
 * efficient recursive computation optimized for complex parameter systems.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.ComplexMatrix;
import jline.util.matrix.Matrix;

public final class Pfqn_gldsingle_complex {
    private Pfqn_gldsingle_complex() {}

    /**
     * Auxiliary function used by pfqn_gld to compute the normalizing constant in a single-class load-dependent model
     * with complex demands.
     *
     * @param L       demands at all stations
     * @param N       number of jobs for each class
     * @param mu      load-dependent scaling factors
     * @param options solver options
     * @return normalizing constant (G) and its logarithm (lG)
     */
    public static Ret.pfqnNcComplex pfqn_gldsingle_complex(ComplexMatrix L, Matrix N, Matrix mu, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        if (R > 1) {
            throw new RuntimeException(
                    "pfqn_gldsingle_complex: multiclass model detected. pfqn_gldsingle_complex is for single class models.");
        }
        Map<Ret.pfqnGldIndex, Complex> g = new HashMap<Ret.pfqnGldIndex, Complex>();
        g.put(new Ret.pfqnGldIndex(1, 1, 1), Complex.valueOf(0.0));
        int n = 1;
        while (n <= N.get(0)) {
            g.put(new Ret.pfqnGldIndex(1, n + 1, 2), Complex.valueOf(0.0));
            n++;
        }
        for (int m = 1; m <= M; m++) {
            int tm = 1;
            while (tm <= N.get(0) + 1) {
                g.put(new Ret.pfqnGldIndex(m + 1, 1, tm + 1), Complex.valueOf(1.0));
                tm++;
            }
            int nn = 1;
            while (nn <= N.get(0)) {
                int tmm = 1;
                while (tmm <= N.get(0) - nn + 1) {
                    g.put(new Ret.pfqnGldIndex(m + 1, nn + 1, tmm + 1),
                            g.get(new Ret.pfqnGldIndex(m, nn + 1, 2)).add(
                                    L.get(m - 1).multiply(g.get(new Ret.pfqnGldIndex(m + 1, nn, tmm + 2)))
                                            .divide(mu.get(m - 1, tmm - 1))));
                    tmm++;
                }
                nn++;
            }
        }
        Complex G = g.get(new Ret.pfqnGldIndex(M + 1, (int) N.get(0) + 1, 2));
        Complex lG = G.log();
        return new Ret.pfqnNcComplex(G, lG);
    }
}
