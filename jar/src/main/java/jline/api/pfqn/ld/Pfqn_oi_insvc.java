/**
 * @file Conditional mean number of in-service jobs at an order-independent (OI) station.
 *
 * This is the quantity underlying the LINE utilization convention at OI
 * stations, U_r = E[sir_r] / c, with c the number of servers and sir_r the
 * number of class-r jobs receiving a strictly positive service rate.
 *
 * Port of matlab/src/api/pfqn/pfqn_oi_insvc.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.function.ToDoubleFunction;

/**
 * Conditional mean number of in-service jobs per class at an OI station.
 *
 * <p>In an OI station the state is the ordered list c = (c_1,...,c_n) of job
 * classes (position 1 = head) and the job in position p is served at the rank
 * rate increment Delta_p(c) = mu(c_1..c_p) - mu(c_1..c_{p-1}), so that the
 * total rate telescopes to mu(c). Position p is IN SERVICE when Delta_p(c) &gt; 0,
 * and sir_r(c) = #{p : c_p = r, Delta_p(c) &gt; 0}. Note that sir_r counts JOBS,
 * not servers: a single job served concurrently by several compatible servers
 * counts once. This matches the definition used by the exact CTMC solver
 * (State.toMarginal, PAS branch) and by LDES.
 *
 * <p>Because mu is permutation-invariant (a function of the count vector), the
 * unnormalized weight of an ordering c of the multiset n factorizes over its
 * prefixes as w(c) = prod_{p=1}^{|n|} 1/mu(n(c_1..c_p)), and the OI balance
 * function Phi(n) = sum_{orderings c of n} w(c) obeys the standard
 * balanced-fairness recursion (condition on the tail element c_{|n|}):
 *
 * <pre>
 *   Phi(0) = 1,   Phi(n) = (1/mu(n)) sum_{r: n_r&gt;0} Phi(n - e_r).
 * </pre>
 *
 * <p>Conditioning the same way on the tail element and using
 * sir_r(c) = sir_r(c_1..c_{|n|-1}) + [c_{|n|} = r] * 1{mu(n) &gt; mu(n - e_r)}
 * gives the companion recursion for the sir-weighted balance
 * Xi_r(n) = sum_{orderings c of n} w(c) sir_r(c):
 *
 * <pre>
 *   Xi_r(0) = 0,
 *   Xi_r(n) = (1/mu(n)) [ sum_{s: n_s&gt;0} Xi_r(n - e_s)
 *                         + 1{n_r &gt; 0} 1{mu(n) &gt; mu(n - e_r)} Phi(n - e_r) ].
 * </pre>
 *
 * <p>Given n, every ordering carries the same class-weight factor, so the
 * conditional law of the ordering is w(c)/Phi(n) and
 * E[sir_r | n] = Xi_r(n) / Phi(n) =: g_r(n), a function of the count vector
 * alone. The station's in-service mean then follows from the count marginal pM
 * by E[sir_r] = sum_n pM(n) g_r(n), or, in normalizing-constant form, from the
 * functional-server identity of {@link Pfqn_oi_fnc} applied to f(n) = g_r(n)
 * (note g_r(0) = 0, as required).
 *
 * @see Pfqn_oi_fnc
 */
public class Pfqn_oi_insvc {

    /** Result of {@link #pfqn_oi_insvc}. */
    public static final class Result {
        /** g[i][r] = E[sir_r | n], column-major over the lattice 0 &lt;= n &lt;= N. */
        public final double[][] g;
        /** Xi[i][r], the sir-weighted balance function. */
        public final double[][] Xi;
        /** Phi[i], the OI balance function. */
        public final double[] Phi;

        public Result(double[][] g, double[][] Xi, double[] Phi) {
            this.g = g;
            this.Xi = Xi;
            this.Phi = Phi;
        }
    }

    /**
     * Conditional mean number of in-service jobs per class over the population
     * lattice 0 &lt;= n &lt;= N.
     *
     * @param oirate rank rate mu(n) on a per-class count vector; mu(0) is taken as 0
     * @param N      (1 x R) closed population vector, finite and nonnegative
     * @return the tables g, Xi and Phi, flat column-major over the lattice
     */
    public static Result pfqn_oi_insvc(ToDoubleFunction<int[]> oirate, int[] N) {
        int R = N.length;
        int[] shp = new int[R];
        int[] stride = new int[R];
        int total = 1;
        for (int d = 0; d < R; d++) {
            if (N[d] < 0) {
                throw new IllegalArgumentException("pfqn_oi_insvc requires nonnegative populations.");
            }
            shp[d] = N[d] + 1;
            total *= shp[d];
        }
        stride[0] = 1;
        for (int d = 1; d < R; d++) {
            stride[d] = stride[d - 1] * shp[d - 1];
        }

        // Decode the lattice once and tabulate the rank rate mu(n).
        int[][] subs = new int[total][R];
        double[] muv = new double[total];
        for (int i = 0; i < total; i++) {
            int li = i;
            int tot = 0;
            for (int d = 0; d < R; d++) {
                subs[i][d] = li % shp[d];
                li /= shp[d];
                tot += subs[i][d];
            }
            if (tot > 0) {
                muv[i] = oirate.applyAsDouble(subs[i]);
            }
        }

        double[] Phi = new double[total];
        double[][] Xi = new double[total][R];
        for (int i = 0; i < total; i++) {
            int[] n = subs[i];
            int tot = 0;
            for (int d = 0; d < R; d++) {
                tot += n[d];
            }
            if (tot == 0) {
                Phi[i] = 1.0;
                continue;
            }
            double mun = muv[i];
            if (mun <= 0) {
                // Unreachable state (no server can serve this composition).
                continue;
            }
            double sPhi = 0;
            double[] sXi = new double[R];
            for (int s = 0; s < R; s++) {
                if (n[s] > 0) {
                    int j = i - stride[s];
                    sPhi += Phi[j];
                    for (int r = 0; r < R; r++) {
                        sXi[r] += Xi[j][r];
                    }
                }
            }
            Phi[i] = sPhi / mun;
            for (int r = 0; r < R; r++) {
                double acc = sXi[r];
                if (n[r] > 0) {
                    int j = i - stride[r];
                    if (mun > muv[j]) {
                        acc += Phi[j];      // the tail class-r job is in service
                    }
                }
                Xi[i][r] = acc / mun;
            }
        }

        double[][] g = new double[total][R];
        for (int i = 0; i < total; i++) {
            if (Phi[i] > 0) {
                for (int r = 0; r < R; r++) {
                    g[i][r] = Xi[i][r] / Phi[i];
                }
            }
        }
        return new Result(g, Xi, Phi);
    }
}
