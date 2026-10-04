/**
 * @file Normalizing constant of a multiclass LIMITED load-dependent closed model
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Same recursion, same arithmetic and the same result as {@link Pfqn_gld}, but
 * with the rate shift saturated at the limited load-dependence threshold, which
 * makes the recursion's state space finite and lets it be memoised. This is
 * what {@link Pfqn_lldsingle} does to {@link Pfqn_gldsingle}, one level up:
 * there the rate offset is an index into a table, here it is the shift that
 * {@link Pfqn_mushift} applies.
 *
 * <p>Pfqn_gld peels the last station and advances its rate lattice one job at a
 * time,
 *
 * <pre>
 *   g(m,n,j) = g(m-1,n,0) + sum_r L[m,r]/alpha_m(j+1) * g(m,n-e_r,j+1)
 * </pre>
 *
 * with j the number of shifts row m has taken, so that Pfqn_mushift's leading
 * element is alpha_m(j+1). Once j &gt;= s_m-1, where s_m is the population past
 * which alpha_m stays constant, every remaining entry of the row is
 * alpha_m(s_m) and a further shift LEAVES THE ROW UNCHANGED over the columns
 * the recursion can still read. Saturating j at s_m-1 therefore returns the
 * same value and makes the state (m, n, j) repeat, at which point one memo
 * answers what Pfqn_gld recomputes down an exponential tree.
 *
 * <p>COST. The state space is M * prod_r(N_r+1) * max_k s_k, against Pfqn_gld's
 * unmemoised recursion, which revisits the same states exponentially often.
 * Without the saturation a memo would still be bounded, but by
 * M * prod_r(N_r+1) * (|N|+1): the threshold is what replaces the population by
 * the server count, exactly as in Pfqn_lldsingle.
 *
 * <p>Every terminal case of Pfqn_gld is delegated back to it on the
 * materialised block, so the two agree to the last bit rather than to a
 * tolerance.
 */
public final class Pfqn_lld {
    private Pfqn_lld() {}

    /**
     * Normalizing constant of a multiclass limited load-dependent closed model.
     *
     * @param L       demands at all stations
     * @param N       number of jobs for each class
     * @param mu      load-dependent scaling factors
     * @param options solver options
     * @return normalizing constant (G) and its logarithm (lG)
     */
    public static Ret.pfqnNc pfqn_lld(Matrix L, Matrix N, Matrix mu, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        int Ntot0 = (int) Math.round(N.elementSum());

        if (mu == null) {
            mu = new Matrix(M, Math.max(Ntot0, 1));
            mu.fill(1.0);
        }
        int ncols = mu.getNumCols();

        if (M == 0 || Ntot0 <= 0 || R == 1 || M == 1) {
            // nothing to memoise: Pfqn_gld returns from a shortcut without recursing
            return Pfqn_gld.pfqn_gld(L, N, mu, options);
        }

        int[] s = thresholds(mu, M, ncols);
        Map<String, Double> memo = new HashMap<String, Double>();
        double G = node(L, N, mu, M, N, 0, s, memo, options, R, Ntot0, ncols);
        return new Ret.pfqnNc(G, FastMath.log(G));
    }

    /**
     * The value Pfqn_gld would return at demands L(0..m-1,:), population n and
     * row m-1 shifted j times. j is saturated by the caller.
     */
    private static double node(Matrix L, Matrix N0, Matrix mu, int m, Matrix n, int j,
                               int[] s, Map<String, Double> memo, SolverOptions options,
                               int R, int Ntot0, int ncols) {
        StringBuilder key = new StringBuilder();
        key.append(m).append('|').append(j).append('|');
        for (int r = 0; r < R; r++) key.append((int) Math.round(n.get(r))).append(',');
        Double hit = memo.get(key.toString());
        if (hit != null) return hit;

        int nsum = (int) Math.round(n.elementSum());
        // Rows 0..m-2 are never shifted, row m-1 is shifted j times, and every
        // row is truncated by one column per job already placed, exactly as
        // Pfqn_mushift leaves them.
        int cols = ncols - (Ntot0 - nsum);
        Matrix Lsub = Matrix.extractRows(L, 0, m, null);
        Matrix muSub = new Matrix(m, Math.max(cols, 0));
        if (cols > 0) {
            // j + cols <= ncols holds for an INTEGER population, since j never
            // exceeds the jobs already placed; the clamp only bites if a caller
            // passes a fractional one, where it keeps the read in range instead
            // of walking off the rate matrix
            int jsrc = Math.min(j, Math.max(ncols - cols, 0));
            for (int i = 0; i < m - 1; i++)
                for (int c = 0; c < cols; c++) muSub.set(i, c, mu.get(i, c));
            for (int c = 0; c < cols; c++) muSub.set(m - 1, c, mu.get(m - 1, jsrc + c));
        }

        double val;
        // Pfqn_gld's own cascade. Every arm below returns from a shortcut of
        // Pfqn_gld without recursing, so delegating keeps the last bit.
        if (m <= 1 || R == 1 || Lsub.isEmpty() || nsum == 0 || !isLoadDep(muSub, m, nsum)) {
            val = Pfqn_gld.pfqn_gld(Lsub, n, muSub, options).G;
            memo.put(key.toString(), val);
            return val;
        }

        // the recursion, memoised. The shift saturates at s[m-1]-1, past which
        // the row the child would see is the one it sees now.
        val = node(L, N0, mu, m - 1, n, 0, s, memo, options, R, Ntot0, ncols);
        for (int r = 0; r < R; r++) {
            if (n.get(r) > GlobalConstants.FineTol) {
                Matrix n1 = n.copy();
                n1.set(r, n1.get(r) - 1);
                val += L.get(m - 1, r) / muSub.get(m - 1, 0)
                        * node(L, N0, mu, m, n1, Math.min(j + 1, s[m - 1] - 1),
                               s, memo, options, R, Ntot0, ncols);
            }
        }
        memo.put(key.toString(), val);
        return val;
    }

    /** Pfqn_gld's own load-dependence scan, on the materialised block. */
    private static boolean isLoadDep(Matrix muSub, int m, int nsum) {
        for (int i = 0; i < m; i++) {
            boolean single = true;
            boolean delay = true;
            for (int c = 0; c < nsum; c++) {
                double v = c < muSub.getNumCols() ? muSub.get(i, c) : 1.0;
                if (FastMath.abs(v - 1.0) > GlobalConstants.FineTol) single = false;
                if (FastMath.abs(v - (c + 1)) > GlobalConstants.FineTol) delay = false;
            }
            if (!single && !delay) return true;
        }
        return false;
    }

    /**
     * Per-station limited load-dependence threshold: the smallest offset past
     * which the rate row is constant. The tolerance is confined to FINITE pairs,
     * since at tail = Inf the bound is itself Inf and every finite rate would
     * tie to it, collapsing the row on a false tie.
     */
    private static int[] thresholds(Matrix mu, int M, int ncols) {
        int[] s = new int[M];
        for (int m = 0; m < M; m++) {
            double tail = mu.get(m, ncols - 1);
            s[m] = ncols;
            for (int n = ncols - 1; n >= 1; n--) {
                double prev = mu.get(m, n - 1);
                if (prev == tail || (Double.isFinite(tail) && Double.isFinite(prev)
                        && Math.abs(prev - tail) <= Math.ulp(1.0) * Math.max(Math.abs(tail), 1.0))) {
                    s[m] = n;
                } else {
                    break;
                }
            }
            if (s[m] < 1) s[m] = 1;
        }
        return s;
    }
}
