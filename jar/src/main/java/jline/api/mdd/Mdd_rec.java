package jline.api.mdd;

/**
 * MDD-rec: the normalising constant of a product-form model whose reachable set
 * is held in a decision diagram.
 *
 * <p>S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant
 * for product-form models of distributed systems with synchronisation", Future
 * Generation Computer Systems 111 (2020) 475-490, Sec. 4.</p>
 *
 * <p>A product-form model has P(s) = (1/G) prod_k g_k(s_k) over its levels, and
 * G = sum_{s in S} prod_k g_k(s_k). Summing state by state is exponential and
 * numerically unstable. MDD-rec instead walks the diagram that already encodes
 * S, accumulating the unnormalised mass of each node ONCE (Def. 4.4,
 * Algorithm 1):</p>
 *
 * <pre>
 *   M(&lt;l.p&gt;) = sum_{v in S_l} g_l(v) * M(&lt;l.p&gt;[v]),  M(TRUE)=1, M(FALSE)=0
 * </pre>
 *
 * <p>so the cost is O(sum_l |nodes_l| * |S_l|) rather than O(|S|), and
 * G = M(root).</p>
 *
 * <p>FORMALISM-AGNOSTIC. Nothing here knows what a level is: the paper's
 * Appendix B shows that on the lattice sum_k s_k = n of a closed queueing
 * network this collapses to Buzen's convolution, and Sec. 5 that on an
 * S-invariant reachable Petri net it collapses to the Coleman-Henderson-Taylor
 * convolution ({@code Spn_conv}). Unlike either, it needs only that the
 * reachable set be finite and encoded -- no lattice, no S-invariant
 * reachability.</p>
 *
 * <p>THE MASK is how Sec. 5.3 computes measures. Restricting the sum at level l
 * to a subset of its local values gives the unnormalised mass of the
 * corresponding subset of S, so P(m_l = k) and P(e_j &gt;= k) are the same
 * recursion under a different mask rather than three separate algorithms.</p>
 *
 * <p>WHAT IS NOT HERE. The g_l themselves, and the test that the model has a
 * product form at all, are the caller's: the paper declares that out of scope
 * (Sec. 3.2). Passing g_l that do not describe a product-form model returns a
 * number that is not the normalising constant of anything, and nothing here can
 * detect it.</p>
 *
 * <p>MATLAB twin: {@code mdd_rec.m}. Python twin: {@code api/mdd/rec.py}.</p>
 */
public class Mdd_rec {

    private Mdd_rec() {}

    /** Shape agreement between the diagram, the factors and the mask. */
    private static void check(MddStruct mdds, double[][] g, boolean[][] mask) {
        if (g == null || g.length != mdds.K) {
            throw new RuntimeException("mdd_rec: one g_l per level is required");
        }
        for (int l = 0; l < mdds.K; l++) {
            if (g[l] == null || g[l].length != mdds.domain[l]) {
                throw new RuntimeException("mdd_rec: g_l must have one entry per local state");
            }
        }
        if (mask != null) {
            if (mask.length != mdds.K) {
                throw new RuntimeException("mdd_rec: the mask must have one row per level");
            }
            for (int l = 0; l < mdds.K; l++) {
                if (mask[l] == null || mask[l].length != mdds.domain[l]) {
                    throw new RuntimeException(
                            "mdd_rec: the mask must have one entry per local state");
                }
            }
        }
    }

    /**
     * Unnormalised mass of the masked subset of the reachable set (Algorithm 1).
     *
     * @param mdds the reachable set, in MDD orientation (level 0 is the root)
     * @param g g[l][v] is g_l(v), the per-level factor of the product form
     * @param mask per-level admissible values; null admits everything and
     *        returns the normalising constant G
     * @return the unnormalised mass of the masked subset
     */
    public static double mdd_rec_masked(MddStruct mdds, double[][] g, boolean[][] mask) {
        check(mdds, g, mask);
        if (mdds.root == MDD.TERM_FALSE) {
            return 0.0;
        }
        double[][] memo = new double[mdds.K][];
        for (int l = 0; l < mdds.K; l++) {
            memo[l] = new double[mdds.nnodes[l]];
        }
        // Levels are swept bottom-up, so every child mass a node needs is
        // already in memo when the node is reached; an explicit sweep avoids the
        // stack depth a level count in the hundreds would otherwise reach.
        for (int l = mdds.K - 1; l >= 0; l--) {
            for (int id = 0; id < mdds.nnodes[l]; id++) {
                int[] arcs = mdds.node[l][id];
                double acc = 0.0;
                for (int v = 0; v < mdds.domain[l]; v++) {
                    if (mask != null && !mask[l][v]) {
                        continue;
                    }
                    double gv = g[l][v];
                    if (gv == 0.0) {
                        continue;
                    }
                    int ch = arcs[v];
                    if (l == mdds.K - 1) {
                        if (ch != MDD.TERM_TRUE) {
                            continue;
                        }
                        acc += gv;
                    } else {
                        if (ch == MDD.TERM_FALSE) {
                            continue;
                        }
                        acc += gv * memo[l + 1][ch - 1];
                    }
                }
                memo[l][id] = acc;
            }
        }
        return memo[0][mdds.root - 1];
    }

    /** The normalising constant G = sum_{s in S} prod_l g_l(s_l). */
    public static double mdd_rec(MddStruct mdds, double[][] g) {
        return mdd_rec_masked(mdds, g, null);
    }

    /**
     * Unnormalised masses of {s in S : s_l = k}, one per local value k of level l.
     *
     * <p>Divided by G these are P(m_l = k) of Sec. 5.3: the mean occupancy of a
     * level is sum_k k * P(m_l = k), and its utilization 1 - P(m_l = 0).</p>
     *
     * @param mdds the reachable set
     * @param g per-level product-form factors
     * @param l level index, 0-based
     * @return one unnormalised mass per local value of level l
     */
    public static double[] mdd_rec_marginal(MddStruct mdds, double[][] g, int l) {
        if (l < 0 || l >= mdds.K) {
            throw new RuntimeException("mdd_rec_marginal: level index is out of range");
        }
        int d = mdds.domain[l];
        double[] out = new double[d];
        for (int k = 0; k < d; k++) {
            boolean[][] mask = new boolean[mdds.K][];
            for (int j = 0; j < mdds.K; j++) {
                mask[j] = new boolean[mdds.domain[j]];
                for (int v = 0; v < mdds.domain[j]; v++) {
                    mask[j][v] = true;
                }
            }
            for (int v = 0; v < d; v++) {
                mask[l][v] = (v == k);
            }
            out[k] = mdd_rec_masked(mdds, g, mask);
        }
        return out;
    }
}
