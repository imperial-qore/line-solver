/**
 * @file JMT-compatible Chow approximate Mean Value Analysis
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * JMT-compatible Chow approximate MVA.
 *
 * <p>JMT's Chow analyzer estimates every class's arrival-instant queue length
 * at a station by the aggregate queue length at the full population. This is
 * the Bard large-customer-population fixed point implemented by
 * {@link Pfqn_lcp}.</p>
 */
public final class Pfqn_chow {
    private Pfqn_chow() {}

    /** Legacy compatibility selector; JMT-compatible Chow ignores it. */
    public static final String FORWARD = "forward";
    /** Legacy compatibility selector; JMT-compatible Chow ignores it. */
    public static final String BACKWARD = "backward";

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N) {
        return Pfqn_lcp.pfqn_lcp(L, N);
    }

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N, Matrix Z) {
        return Pfqn_lcp.pfqn_lcp(L, N, Z);
    }

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0) {
        return Pfqn_lcp.pfqn_lcp(L, N, Z, tol, maxiter, QN0);
    }

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N, Matrix Z, double tol, int maxiter,
                                         Matrix QN0, SchedStrategy[] type) {
        return Pfqn_lcp.pfqn_lcp(L, N, Z, tol, maxiter, QN0, type);
    }

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N, Matrix Z, double tol, int maxiter,
                                         Matrix QN0, SchedStrategy[] type, String variant) {
        return Pfqn_lcp.pfqn_lcp(L, N, Z, tol, maxiter, QN0, type);
    }
}
