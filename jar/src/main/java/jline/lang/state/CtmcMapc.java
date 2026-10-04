package jline.lang.state;

import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * Pair form of a multiserver MAP (or MMPP2) service in the CTMC state space, set by
 * {@code Solver_ctmc_mapc.apply}; port of the MATLAB sn.ctmcmapc{ist,r} struct.
 *
 * <p>With V = (-D0)^-1 D1 a draw started in phase h ends in j with probability V(h,j). A busy
 * server is a pair (i,j): current phase i and landing j, fixed when the draw starts. The class
 * memory local variable (1-based) is the landing of the most recently started draw: a start
 * from h enters pair (h,j) w.p. V(h,j) and sets h := j; phase moves and completions keep it.</p>
 */
public final class CtmcMapc implements Serializable {
    private static final long serialVersionUID = 1L;

    /** Number of phases of the original MAP. */
    public final int p;
    /** Pair t = (pairs[t][0], pairs[t][1]), 0-based phases, i-major. */
    public final int[][] pairs;
    /** Landing probabilities V = (-D0)^-1 D1 of the original MAP. */
    public final Matrix V;
    /** Completion rate D1(i,j)/V(i,j) of each pair. */
    public final double[] done;

    public CtmcMapc(int p, int[][] pairs, Matrix V, double[] done) {
        this.p = p;
        this.pairs = pairs;
        this.V = V;
        this.done = done;
    }

    /** The pair form of class r at station ist, or null when the service is not in pair form. */
    public static CtmcMapc get(jline.lang.NetworkStruct sn, int ist, int r) {
        if (sn.ctmcmapc == null) {
            return null;
        }
        CtmcMapc[] row = sn.ctmcmapc.get(ist);
        return row == null || r < 0 || r >= row.length ? null : row[r];
    }
}
