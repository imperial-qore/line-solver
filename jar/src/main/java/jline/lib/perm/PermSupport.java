package jline.lib.perm;

import java.util.Arrays;

import jline.util.matrix.Matrix;

/**
 * Preconditions and helpers shared by the approximate permanent engines.
 *
 * Twin of MATLAB perm_require_support.m / perm_maxweight_assignment.m and of
 * the python line_solver.api.perm.base.require_full_support together with
 * HuberLawSampler._hungarian_assignment.
 */
public final class PermSupport {

    private PermSupport() {}

    /**
     * Refuse a matrix the permanent approximations cannot take.
     *
     * The four approximate engines -- BethePermanent, HeuristicPermanent,
     * HuberLawSampler and AdaPartSampler -- all rest, directly or through the
     * Sinkhorn scaling they share, on a strictly positive matrix.
     *
     * Why this is refused rather than floored. Each of these used to replace a
     * zero by a small constant eps before doing anything else, and that
     * substitution is not invertible: every permutation picks exactly one entry
     * from each row, so a matrix with an identically zero row has permanent 0
     * while the floored matrix has permanent n! eps times the permanent of the
     * rest. n! outruns eps quickly -- with eps = 2.22e-16 the fabricated value
     * passes 1% at n = 17 and 1 at n = 18, and at n = 20 the floor alone
     * manufactures a permanent of about 540 where the truth is exactly zero.
     * The order of the replicated demand matrix in Pfqn_jointmarg is sum(N), so
     * a closed model with 18 jobs and one structurally zero demand is already
     * in that regime.
     *
     * Positivity is sufficient but not necessary: the sharp precondition of the
     * Sinkhorn scaling is TOTAL SUPPORT (every positive entry lies on a
     * positive permutation), which a matrix with a strictly positive permanent
     * can still fail -- [[J3, 0], [J3, J3]] has permanent 36, no total support,
     * and the scaling exits on its tolerance rather than on convergence.
     * Positivity is used here because it is O(n^2) and is the contract the C++
     * header already states.
     *
     * @param matrix matrix to check
     * @param caller engine name, used in the message
     * @throws IllegalArgumentException if any entry is not strictly positive
     */
    public static void requireFullSupport(Matrix matrix, String caller) {
        int rows = matrix.getNumRows();
        int cols = matrix.getNumCols();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double v = matrix.get(i, j);
                if (!(v > 0.0)) {
                    throw new IllegalArgumentException("The '" + caller
                            + "' permanent approximation requires a strictly positive matrix: entry ("
                            + i + ", " + j + ") is " + v
                            + ", so the matrix has no full support. Flooring it would change the"
                            + " permanent by n!*eps, which is O(1) by n=18. Use the exact engine"
                            + " (Permanent / Pfqn_perm).");
                }
            }
        }
    }

    /**
     * Maximum-weight perfect assignment of a square weight matrix.
     *
     * This replaces the row-by-row greedy that HuberLawSampler called
     * hungarianAssignment. That routine was greedy despite the name and returns
     * a zero-weight assignment on inputs that admit a positive one: on
     * [[1, 2], [0, 3]] row 0 takes the larger entry in column 1, leaving row 1
     * with the zero in column 0.
     *
     * That matters because the weight is alpha3 in the sampler's rescale step,
     * and alpha3 sets the flooring level alpha1 = alpha3*delta/(3 n!) of the
     * Huber-Law bound -- the one principled zero-handling in the whole family.
     * A suboptimal assignment understates alpha3, hence alpha1, and weakens the
     * method's own guarantee. A correct assignment makes alpha3 &gt; 0 whenever
     * the matrix has a positive permanent.
     *
     * O(n^3) Hungarian algorithm by shortest augmenting paths, run on the
     * negated weights so that it minimizes. The MATLAB and python twins solve
     * the same problem and agree on the optimal VALUE, which is all alpha3
     * depends on; they need not return the same permutation when the optimum is
     * degenerate.
     *
     * @param weight square weight matrix to maximize over
     * @return assignment[i] is the column matched to row i
     */
    public static int[] maxWeightAssignment(Matrix weight) {
        final int n = weight.getNumRows();
        int[] assignment = new int[n];
        if (n == 0) {
            return assignment;
        }
        Arrays.fill(assignment, -1);

        // cost[i][j] = -weight, 1-indexed as the augmenting-path formulation needs
        final double[][] cost = new double[n + 1][n + 1];
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                cost[i][j] = -weight.get(i - 1, j - 1);
            }
        }

        final double[] u = new double[n + 1];
        final double[] v = new double[n + 1];
        final int[] p = new int[n + 1];
        final int[] way = new int[n + 1];

        for (int i = 1; i <= n; i++) {
            p[0] = i;
            int j0 = 0;
            final double[] minv = new double[n + 1];
            final boolean[] used = new boolean[n + 1];
            Arrays.fill(minv, Double.POSITIVE_INFINITY);
            do {
                used[j0] = true;
                final int i0 = p[j0];
                int j1 = -1;
                double delta = Double.POSITIVE_INFINITY;
                for (int j = 1; j <= n; j++) {
                    if (!used[j]) {
                        final double cur = cost[i0][j] - u[i0] - v[j];
                        if (cur < minv[j]) {
                            minv[j] = cur;
                            way[j] = j0;
                        }
                        if (minv[j] < delta) {
                            delta = minv[j];
                            j1 = j;
                        }
                    }
                }
                for (int j = 0; j <= n; j++) {
                    if (used[j]) {
                        u[p[j]] += delta;
                        v[j] -= delta;
                    } else {
                        minv[j] -= delta;
                    }
                }
                j0 = j1;
            } while (p[j0] != 0);
            do {
                final int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0 != 0);
        }

        for (int j = 1; j <= n; j++) {
            if (p[j] != 0) {
                assignment[p[j] - 1] = j - 1;
            }
        }
        return assignment;
    }
}
