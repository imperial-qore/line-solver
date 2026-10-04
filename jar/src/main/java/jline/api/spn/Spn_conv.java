package jline.api.spn;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import jline.api.spn.Spn_sinvariants.SpnInvariants;

/**
 * Convolution algorithm for the normalising constant of an S-invariant reachable
 * product-form stochastic Petri net.
 *
 * <p>J. Coleman, W. Henderson, P. Taylor, "Product form equilibrium
 * distributions and a convolution algorithm for stochastic Petri nets",
 * Performance Evaluation 26(3), 1996, 159-180; presented as the point of
 * comparison for MDD-rec in S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020),
 * Sec. 5.1.</p>
 *
 * <p>With S the minimal-support S-invariant matrix and V = S m0 the load vector,
 * the reachability set of an S-INVARIANT REACHABLE net is exactly
 * {m &gt;= 0 : S m = V}, and conditioning on the marking of one place partitions
 * it (Lemma 5.1). Writing G_j(W) for the mass of the markings supported on the
 * first j places with S m = W,</p>
 *
 * <pre>
 *   G_0(W) = [W == 0],     G_j(W) = sum_i g_j(i) G_{j-1}(W - i S_j),
 * </pre>
 *
 * <p>and G = G_n(V). On a net whose only invariant is "the tokens are
 * conserved" this is Buzen's convolution for a closed queueing network, one
 * place per station.</p>
 *
 * <p>NO ILP IS SOLVED. The paper obtains the marking set M_p(P',W) from the
 * feasibility of an integer program (Prop. 5.2) so that the sum skips the terms
 * that contribute nothing. Here the sum simply runs over i whose residual
 * W - i S_j stays non-negative and the recursion returns zero on an infeasible
 * residual, which gives the same value: the ILP is an optimisation of the
 * enumeration, not part of the definition. Memoising on (j, W) keeps the walk
 * over the reachable residuals rather than over all of them.</p>
 *
 * <p>S-INVARIANT REACHABILITY IS NOT CHECKED, and cannot be cheaply: no
 * algorithm is known that decides it without generating the reachability set
 * (FGCS, Sec. 5.1). On a net that fails it, {m : S m = V} is strictly larger
 * than the reachable set and this returns a normalising constant over
 * unreachable markings too, which is why {@code Mdd_rec.mdd_rec} -- which walks
 * the reachable set itself -- is the general algorithm and this one the special
 * case. Compare the two on a new net before trusting this one on it.</p>
 *
 * <p>MATLAB twin: {@code spn_conv.m}. Python twin: {@code api/spn/conv.py}.</p>
 */
public class Spn_conv {

    private Spn_conv() {}

    /** Convolve straight off the invariant basis {@code Spn_sinvariants} returned. */
    public static double spn_conv(SpnInvariants inv, double[][] g) {
        return spn_conv(inv.S, inv.V, g);
    }

    /**
     * The normalising constant by convolution over the invariant load vector.
     *
     * @param S S[i][p], the minimal-support S-invariants, one row per invariant
     * @param V the load vector S m0, one entry per invariant
     * @param g g[p][i] is g_p(i), the product-form factor of i tokens in place
     *        level p; its length bounds the marking of that level
     * @return the normalising constant
     */
    public static double spn_conv(long[][] S, long[] V, double[][] g) {
        if (S == null || S.length == 0) {
            throw new RuntimeException("spn_conv: the net has no S-invariant to convolve over");
        }
        if (S.length != V.length) {
            throw new RuntimeException("spn_conv: one load-vector entry per invariant is required");
        }
        int n = g.length;
        for (int r = 0; r < S.length; r++) {
            if (S[r].length != n) {
                throw new RuntimeException(
                        "spn_conv: the invariant matrix and g must agree on the place-level count");
            }
            for (int p = 0; p < n; p++) {
                if (S[r][p] < 0) {
                    throw new RuntimeException("spn_conv: an S-invariant has a negative weight, so "
                            + "the residual recursion has no monotone bound on the marking");
                }
            }
        }
        Map<String, Double> memo = new HashMap<String, Double>();
        return rec(n, V.clone(), S, g, memo);
    }

    private static double rec(int j, long[] W, long[][] S, double[][] g, Map<String, Double> memo) {
        if (j == 0) {
            for (int r = 0; r < W.length; r++) {
                if (W[r] != 0) {
                    return 0.0;
                }
            }
            return 1.0;
        }
        String key = j + "|" + Arrays.toString(W);
        Double hit = memo.get(key);
        if (hit != null) {
            return hit.doubleValue();
        }
        int p = j - 1;
        double acc = 0.0;
        for (int i = 0; i < g[p].length; i++) {
            long[] rem = new long[W.length];
            boolean feasible = true;
            for (int r = 0; r < W.length && feasible; r++) {
                rem[r] = W[r] - ((long) i) * S[r][p];
                if (rem[r] < 0) {
                    feasible = false;
                }
            }
            if (!feasible) {
                break;                 // S is non-negative, so larger i only gets worse
            }
            if (g[p][i] == 0.0) {
                continue;
            }
            acc += g[p][i] * rec(j - 1, rem, S, g, memo);
        }
        memo.put(key, Double.valueOf(acc));
        return acc;
    }
}
