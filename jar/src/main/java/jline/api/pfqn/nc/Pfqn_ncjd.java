/**
 * Joint-dependent name of {@link Pfqn_ncoi}: the balance-function convolution of
 * a closed network whose station rates read the whole per-class occupancy vector.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.io.Ret;

import java.util.List;
import java.util.function.ToDoubleFunction;

/**
 * The two names denote the SAME routine because the balanced-fairness recursion
 * Phi_i(0) = 1, mu_i(n) Phi_i(n) = sum_{r: n_r&gt;0} v_{i,r} Phi_i(n - e_r) never
 * inspects the structure of mu_i: it evaluates the handle at the full count
 * vector n. Order independence (mu_i constant on each support) is a modelling
 * restriction that buys insensitivity and a physical reading of Phi, not
 * something the convolution uses, so any joint-dependent scaling eta_i(n)
 * (sn.jdscaling) is admissible, with the proviso that the product form it
 * induces is the balanced-fair one matched to that rate.
 *
 * <p>Use {@link Pfqn_ncoi} when the model is genuinely order independent and the
 * name should say so, this class when the rate is a general joint dependence.
 * {@link Pfqn_clwjd} is the transform route and is NOT merely a renaming: it
 * needs the rate to saturate at a finite cutoff.
 */
public final class Pfqn_ncjd {
    private Pfqn_ncjd() {}

    /** @see Pfqn_ncoi#pfqn_ncoi(double[], int[], List) */
    public static Ret.pfqnOiNc pfqn_ncjd(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu) {
        return Pfqn_ncoi.pfqn_ncoi(Z, N, mu, null);
    }

    /** @see Pfqn_ncoi#pfqn_ncoi(double[], int[], List, double[][]) */
    public static Ret.pfqnOiNc pfqn_ncjd(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu,
                                         double[][] visits) {
        return Pfqn_ncoi.pfqn_ncoi(Z, N, mu, visits);
    }
}
