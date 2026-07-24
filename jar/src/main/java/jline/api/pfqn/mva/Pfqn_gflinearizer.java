/**
 * @file General-form linearizer approximate MVA
 *
 * Implements the general-form linearizer approximation for closed queueing networks with
 * configurable linearization parameters. Provides wrapper for the extended general-form
 * linearizer algorithm with scalar linearization factor.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

public final class Pfqn_gflinearizer {
    private Pfqn_gflinearizer() {}

    /** General-form linearizer algorithm */
    public static Ret.pfqnAMVA pfqn_gflinearizer(Matrix L,
                                                 Matrix N,
                                                 Matrix Z,
                                                 SchedStrategy[] type,
                                                 double tol,
                                                 int maxiter,
                                                 double alpha) {
        return pfqn_gflinearizer(L, N, Z, type, tol, maxiter, alpha, null);
    }

    public static Ret.pfqnAMVA pfqn_gflinearizer(Matrix L,
                                                 Matrix N,
                                                 Matrix Z,
                                                 SchedStrategy[] type,
                                                 double tol,
                                                 int maxiter,
                                                 double alpha,
                                                 Matrix QN0) {
        Matrix alphaM = new Matrix(1, N.getNumCols());
        alphaM.fill(alpha);
        return Pfqn_egflinearizer.pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alphaM, QN0);
    }
}
