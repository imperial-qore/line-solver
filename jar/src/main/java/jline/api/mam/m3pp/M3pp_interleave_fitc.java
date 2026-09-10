/**
 * @file M3PP interleaved multi-process fitting
 *
 * Fits k second-order M3PP[m_j] and interleaves them into an M3PP[m] of order
 * k+1, matching a per-partition IDC at two time scales and a per-class
 * variance-plus-covariance.
 *
 * NOT AVAILABLE, and refused by name rather than approximated. The MATLAB
 * reference matlab/lib/m3a/m3a/m3pp/m3pp_interleave_fitc.m calls
 * compute_feasible_interleave and compute_feasible_interleave_reorder, which
 * are defined NOWHERE in the LINE tree -- neither in m3a nor in
 * matlab/src/api -- so the reference raises "Unrecognized function" on its
 * first call and cannot execute. There is therefore no ground truth to port,
 * and no way to check a port if one were written: the two missing functions
 * are exactly the step that decides which per-partition off-diagonal rates
 * admit an interleaving, which is the whole content of the algorithm.
 *
 * What stood here until this was written was an ad-hoc substitute rather than
 * a port. It set lambda2 = rate * (1 - IDC(t)) / 2, which is NEGATIVE for
 * every input the family admits (the reference's own compute_d refuses unless
 * IDC(t) > 1), so its D1 carried a negative arrival rate and the result was
 * not a MAP at all; it also ignored the binf, t and tinf arguments entirely,
 * truncated the interleaved order at ten, and rescaled the class matrices by a
 * factor that broke sum_c Dc = D1. Returning that silently is worse than
 * returning nothing.
 *
 * WHAT TO USE INSTEAD. The two composites of this family that ARE executable
 * in MATLAB are ported and cover the same ground:
 *  - {@link M3pp_superpos_fitc}, which fits one second-order process per class
 *    and superposes them, and
 *  - {@link M3pp22_interleave_fitc}, which solves for interleaving-admissible
 *    off-diagonal rates with a linear program before fitting the components,
 *    which is what the missing helpers were meant to do.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import jline.util.Pair;
import jline.util.matrix.MatrixCell;

import java.util.List;

public final class M3pp_interleave_fitc {
    private M3pp_interleave_fitc() {}

    private static final String UNPORTED =
            "m3pp_interleave_fitc is not available: its MATLAB reference "
            + "(matlab/lib/m3a/m3a/m3pp/m3pp_interleave_fitc.m) calls "
            + "compute_feasible_interleave and compute_feasible_interleave_reorder, which are "
            + "defined nowhere in the LINE tree, so the reference itself cannot execute and "
            + "there is no algorithm to port. Use m3pp22_interleave_fitc, which solves for "
            + "interleaving-admissible off-diagonal rates with a linear program, or "
            + "m3pp_superpos_fitc.";

    /**
     * Fits k second-order M3PP[m_j] and interleaves them into an M3PP[m] of order k+1.
     *
     * @throws UnsupportedOperationException always; see the file comment.
     */
    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc(
            double[] av,
            double[] btv,
            double[] binfv,
            double[][] acc,
            double[][] gtcc,
            double t,
            double tinf,
            int[] mapping,
            boolean reorder) {
        throw new UnsupportedOperationException(UNPORTED);
    }

    /**
     * Fits k second-order M3PP[m_j] and interleaves them, without reordering.
     *
     * @throws UnsupportedOperationException always; see the file comment.
     */
    public static Pair<MatrixCell, List<MatrixCell>> m3pp_interleave_fitc(
            double[] av, double[] btv, double[] binfv,
            double[][] acc, double[][] gtcc, double t, double tinf) {
        throw new UnsupportedOperationException(UNPORTED);
    }
}
