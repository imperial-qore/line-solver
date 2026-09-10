/**
 * Linearizer Approximate MVA for Product-Form Networks
 *
 * Implements the linearizer approximate MVA method for large closed queueing networks
 * where exact MVA becomes computationally prohibitive. Provides near-exact accuracy
 * with significantly reduced computational complexity for multi-class systems.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

public final class Pfqn_linearizer {
    private Pfqn_linearizer() {}

    /**
     * Linearizer approximate mean value analysis algorithm
     */
    public static Ret.pfqnAMVA pfqn_linearizer(Matrix L,
                                                Matrix N,
                                                Matrix Z,
                                                SchedStrategy[] type,
                                                double tol,
                                                int maxiter) {
        return Pfqn_gflinearizer.pfqn_gflinearizer(L, N, Z, type, tol, maxiter, 1.0);
    }

    public static Ret.pfqnAMVA pfqn_linearizer(Matrix L,
                                                Matrix N,
                                                Matrix Z,
                                                SchedStrategy[] type,
                                                double tol,
                                                int maxiter,
                                                Matrix QN0) {
        return Pfqn_gflinearizer.pfqn_gflinearizer(L, N, Z, type, tol, maxiter, 1.0, QN0);
    }

    public static Ret.pfqnAMVA pfqn_linearizer(Matrix L,
                                                Matrix N,
                                                Matrix Z,
                                                SchedStrategy[] type,
                                                double tol) {
        return pfqn_linearizer(L, N, Z, type, tol, 1000);
    }

    public static Ret.pfqnAMVA pfqn_linearizer(Matrix L,
                                                Matrix N,
                                                Matrix Z,
                                                SchedStrategy[] type) {
        return pfqn_linearizer(L, N, Z, type, 1.0e-8, 1000);
    }
}
