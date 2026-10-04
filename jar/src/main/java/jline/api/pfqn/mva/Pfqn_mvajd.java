/**
 * Joint-dependent name of {@link Pfqn_mvaoi}: the mean-value analysis of a closed
 * network whose station rates read the whole per-class occupancy vector.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.List;
import java.util.function.ToDoubleFunction;

/**
 * The two names denote the SAME routine because the recursion evaluates the rate
 * handle at a full occupancy vector, mu_i(s_i + e_r) with s_i the shift (the
 * occupancy already committed at the bottom of station i), and never inspects the
 * structure of mu_i. That is exactly the "third form" of the Conditional MVA of
 * Casale (QUESTA 2009), a rate depending on the full per-class occupancy vector,
 * so any joint-dependent scaling eta_i(n) (sn.jdscaling) is admissible.
 *
 * <p>Unlike the AMVA joint-dependence route (Solver_amvald with sn.jdscaling,
 * which evaluates eta at the MEAN arrival-instant vector 1 + E[Q] and therefore
 * collapses a support indicator to 1), this routine evaluates the rate at exact
 * integer occupancies and is exact for the balanced-fair station, at the cost of
 * walking prod_r C(N_r+K+1,K+1) states with K joint-dependent stations.
 */
public final class Pfqn_mvajd {
    private Pfqn_mvajd() {}

    /** @see Pfqn_mvaoi#pfqn_mvaoi(double[], int[], ToDoubleFunction) */
    public static Pfqn_mvaoi.Result pfqn_mvajd(double[] Z, int[] N, ToDoubleFunction<int[]> mu) {
        return Pfqn_mvaoi.pfqn_mvaoi(Z, N, mu);
    }

    /** @see Pfqn_mvaoi#pfqn_mvaoi(double[], int[], List, double[][]) */
    public static Pfqn_mvaoi.Result pfqn_mvajd(double[] Z, int[] N,
                                               List<ToDoubleFunction<int[]>> mu, double[][] Dli) {
        return Pfqn_mvaoi.pfqn_mvaoi(Z, N, mu, Dli);
    }

    /** @see Pfqn_mvaoi#pfqn_mvaoi(double[], int[], List, double[][], double[][]) */
    public static Pfqn_mvaoi.Result pfqn_mvajd(double[] Z, int[] N,
                                               List<ToDoubleFunction<int[]>> mu, double[][] Dli,
                                               double[][] visits) {
        return Pfqn_mvaoi.pfqn_mvaoi(Z, N, mu, Dli, visits);
    }
}
