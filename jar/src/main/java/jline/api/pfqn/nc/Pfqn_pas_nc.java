/**
 * @file Pass-and-swap (P&amp;S) normalizing constant of one communicating class.
 *
 * With a non-empty swap graph the ordered-state chain of a P&amp;S network is
 * reducible (Comte and Dorsman, 2021, arXiv:2009.12299): the recurrent
 * communicating classes are the placement-order-adhering sets and the product
 * form holds per class. This routine returns the per-class constant G_C by a
 * microstate head-peeling walk over the feasible orderings.
 *
 * Port of matlab/src/api/pfqn/pfqn_pas_nc.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.io.Ret;
import org.apache.commons.math3.special.Gamma;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;

public final class Pfqn_pas_nc {
    private Pfqn_pas_nc() {}

    /**
     * Normalizing constant of a closed P&amp;S network with no placement order,
     * i.e. the plain OI constant. Prefer {@link Pfqn_ncoi}, which computes the
     * same value on the count lattice at far lower cost.
     *
     * @param Z  (R) think-time demand of the aggregated delay node.
     * @param N  (R) closed population, finite.
     * @param mu list of P&amp;S station rate handles.
     * @return G(N) and log G(N).
     */
    public static Ret.pfqnOiNc pfqn_pas_nc(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu) {
        return pfqn_pas_nc(Z, N, mu, null);
    }

    /**
     * Normalizing constant G_C of the communicating class selected by the
     * placement order.
     *
     * <p>Method: build station M's chain head-first; appending class r at chain
     * position k = sum(occ)+1 is admissible iff no class already placed at that
     * station must come after r, and contributes the reciprocal OI prefix rate
     * 1/mu_M(occ+e_r); the chain may be finalized (recursing to station M-1)
     * only when occ is a placement-order ideal at full multiplicity. When every
     * P&amp;S station has been peeled the residual population sits at the delay
     * node with the multinomial weight prod_r Z_r^{N_r}/N_r!.
     *
     * <p>This is a MICROSTATE routine: it walks the ordered chains position by
     * position, because with a placement order the reachable set is a set of
     * ORDERINGS that does not collapse onto the count lattice. Cost: one node
     * per feasible ordered prefix, which with an empty order is
     * sum_{b&lt;=N} C(|b|+M-1,M-1) |b|!/prod_r b_r!, factorial in sum(N).
     *
     * @param Z    (R) think-time demand of the aggregated delay node; null or
     *             empty for a network with no delay station.
     * @param N    (R) closed population, finite.
     * @param mu   list of P&amp;S station rate handles; each maps a per-class
     *             count vector to the total service rate of that station.
     * @param prec placement order, one (R x R) matrix per station, with
     *             prec.get(m)[i][j] true iff class i must be placed before
     *             class j at station m (the closure of
     *             {@code Pas_placement.pas_placement}, fed by the global DAG of
     *             {@code Pas_swap2order.pas_swap2order}). A single-entry list is
     *             broadcast to every station; null or empty means no order, so
     *             G is the plain OI constant. NOTE the orientation: around a
     *             cycle each downstream station traverses its chain in the
     *             opposite direction, so downstream stations take the TRANSPOSE
     *             of the upstream order. Passing the same matrix to both
     *             stations of a cycle silently returns a smaller, wrong G.
     * @return G_C and log G_C.
     */
    public static Ret.pfqnOiNc pfqn_pas_nc(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu,
                                           List<boolean[][]> prec) {
        int R = N.length;
        List<ToDoubleFunction<int[]>> muList = (mu == null) ? new ArrayList<ToDoubleFunction<int[]>>() : mu;
        int M = muList.size();

        if (Z == null || Z.length == 0) {
            Z = new double[R];
        }
        if (Z.length != R) {
            throw new RuntimeException("pfqn_pas_nc: Z and N must have the same number of classes.");
        }
        for (int r = 0; r < R; r++) {
            if (N[r] < 0) {
                throw new RuntimeException("pfqn_pas_nc requires finite (closed) populations.");
            }
        }

        // Normalize the placement order to one (R x R) matrix per station.
        List<boolean[][]> precList = new ArrayList<boolean[][]>();
        if (prec == null || prec.isEmpty()) {
            for (int m = 0; m < M; m++) {
                precList.add(new boolean[R][R]);
            }
        } else if (prec.size() == 1 && M > 1) {
            for (int m = 0; m < M; m++) {
                precList.add(prec.get(0) == null ? new boolean[R][R] : prec.get(0));
            }
        } else {
            if (prec.size() != M) {
                throw new RuntimeException("pfqn_pas_nc: prec must supply one precedence matrix per station.");
            }
            for (int m = 0; m < M; m++) {
                precList.add(prec.get(m) == null ? new boolean[R][R] : prec.get(m));
            }
        }
        for (int m = 0; m < M; m++) {
            boolean[][] p = precList.get(m);
            if (p.length != R) {
                throw new RuntimeException("pfqn_pas_nc: each precedence matrix must be R x R.");
            }
            for (int i = 0; i < R; i++) {
                if (p[i].length != R) {
                    throw new RuntimeException("pfqn_pas_nc: each precedence matrix must be R x R.");
                }
            }
        }

        double G = rec(Z, N.clone(), N, muList, precList, M, new int[R]);
        double lG = (G > 0) ? Math.log(G) : Double.NEGATIVE_INFINITY;
        return new Ret.pfqnOiNc(G, lG);
    }

    /** Head-peeling recursion over stations m, m-1, ..., 1 and the delay node. */
    private static double rec(double[] Z, int[] N, int[] Norig, List<ToDoubleFunction<int[]>> mu,
                              List<boolean[][]> prec, int m, int[] occ) {
        int R = N.length;

        if (m == 0) {
            // Base case: the residual population sits at the delay node with
            // unnormalized weight prod_r Z_r^{N_r} / N_r!.
            double logf = 0.0;
            for (int r = 0; r < R; r++) {
                if (N[r] > 0) {
                    if (Z[r] <= 0) {
                        return 0.0; // population but no delay demand: infeasible
                    }
                    logf += N[r] * Math.log(Z[r]) - Gamma.logGamma(N[r] + 1.0);
                }
            }
            return Math.exp(logf);
        }

        boolean[][] precm = prec.get(m - 1);

        // Step A: finalize station m's chain here and peel to station m-1, but
        // only if the accumulated occupancy is a placement-order ideal.
        double G = isIdeal(occ, precm, Norig) ? rec(Z, N, Norig, mu, prec, m - 1, new int[R]) : 0.0;

        // Step B: append one more class r at the next chain position of station m.
        ToDoubleFunction<int[]> active = mu.get(m - 1);
        for (int r = 0; r < R; r++) {
            if (N[r] <= 0) {
                continue;
            }
            boolean blocked = false;
            for (int j = 0; j < R; j++) {
                if (occ[j] > 0 && precm[r][j]) {
                    blocked = true;
                    break;
                }
            }
            if (blocked) {
                continue;
            }
            occ[r]++;
            double mu_r = active.applyAsDouble(occ.clone());
            if (mu_r > 0) {
                N[r]--;
                G += (1.0 / mu_r) * rec(Z, N, Norig, mu, prec, m, occ);
                N[r]++;
            }
            occ[r]--;
        }
        return G;
    }

    /**
     * True iff occ is a placement-order ideal at full multiplicity: for every
     * i prec j, occ[j] &gt; 0 requires occ[i] == Norig[i].
     */
    private static boolean isIdeal(int[] occ, boolean[][] precm, int[] Norig) {
        int R = occ.length;
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < R; j++) {
                if (precm[i][j] && occ[j] > 0 && occ[i] < Norig[i]) {
                    return false;
                }
            }
        }
        return true;
    }
}
