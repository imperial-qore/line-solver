/**
 * @file Order-independent (OI) balanced-fairness normalizing constant.
 *
 * Computes the exact normalizing constant G(N) of a closed product-form
 * queueing network comprising an arbitrary number of order-independent (OI) /
 * pass-and-swap stations with empty swap graph and a single aggregated
 * infinite-server (delay) node. The OI stations are analyzed by the
 * balanced-fairness recursion of Bonald and Proutiere (2003), combined with
 * the multichain convolution over stations on the count lattice.
 *
 * Port of matlab/src/api/pfqn/pfqn_ncoi.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.io.Ret;
import org.apache.commons.math3.special.Gamma;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;

public final class Pfqn_ncoi {
    private Pfqn_ncoi() {}

    /**
     * Normalizing constant of a closed OI + single-delay product-form network.
     *
     * <p>The balance function of each OI station is the balanced-fairness
     * recursion Phi(0) = 1, Phi(n) = (1/mu(n)) sum_{r: n_r&gt;0} Phi(n - e_r),
     * and G(N) is the convolution of those balance functions with the
     * multinomial delay factor F_Z(n) = prod_r Z_r^{n_r}/n_r!.
     *
     * <p>This is a MACROSTATE routine: everything is tabulated over the count
     * lattice 0 &lt;= n &lt;= N, never over orderings, which is legitimate
     * because an OI rate is permutation-invariant so Phi closes on the count
     * vector. With a non-empty swap graph that closure fails and
     * {@link Pfqn_pas_nc} must be used instead. Cost: O(M R L) for the balance
     * functions and O(M prod_r (N_r+1)(N_r+2)/2) for the convolutions, with
     * L = prod_r (N_r+1).
     *
     * @param Z  (R) think-time demand of the aggregated delay node.
     * @param N  (R) closed population, finite.
     * @param mu list of OI station rate handles; each maps a per-class count
     *           vector n (length R) to the total service rate at that station.
     *           A state with non-positive rate is unreachable and gets a zero
     *           balance value. May be empty/null to model a pure delay network.
     * @return G(N) and log G(N).
     */
    public static Ret.pfqnOiNc pfqn_ncoi(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu) {
        return pfqn_ncoi(Z, N, mu, null);
    }

    // General per-OI-station class visit ratios: visits[m] is the 1xR visit
    // vector of OI station m, entering the v-weighted balanced-fairness balance
    // Phi^v(n)=(1/mu(n)) sum_r v_r Phi^v(n-e_r). null -> unit visits.
    public static Ret.pfqnOiNc pfqn_ncoi(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu, double[][] visits) {
        int R = N.length;
        List<ToDoubleFunction<int[]>> muList = (mu == null) ? new ArrayList<ToDoubleFunction<int[]>>() : mu;

        if (Z == null || Z.length == 0) {
            Z = new double[R];
        }
        if (Z.length != R) {
            throw new RuntimeException("pfqn_ncoi: Z and N must have the same number of classes.");
        }
        for (int r = 0; r < R; r++) {
            if (N[r] < 0) {
                throw new RuntimeException("pfqn_ncoi requires finite (closed) populations.");
            }
        }
        if (R == 0) {
            return new Ret.pfqnOiNc(1.0, 0.0, new double[]{1.0});
        }

        // Count lattice 0 <= n <= N, flattened column-major; lin(v) = sum(v*strides)
        // is linear, so lin(x+y) = lin(x) + lin(y), which the convolution exploits.
        int[] dims = new int[R];
        int[] strides = new int[R];
        int ngrid = 1;
        for (int r = 0; r < R; r++) {
            dims[r] = N[r] + 1;
            strides[r] = ngrid;
            ngrid *= dims[r];
        }

        int[][] counts = new int[ngrid][R];
        int total = 0;
        for (int r = 0; r < R; r++) {
            total += N[r];
        }
        // Lattice points bucketed by total population: the balanced-fairness
        // recursion must be swept in population-increasing order.
        List<List<Integer>> byPop = new ArrayList<List<Integer>>();
        for (int s = 0; s <= total; s++) {
            byPop.add(new ArrayList<Integer>());
        }
        for (int k = 0; k < ngrid; k++) {
            int rest = k;
            int sum = 0;
            for (int r = 0; r < R; r++) {
                counts[k][r] = rest % dims[r];
                rest /= dims[r];
                sum += counts[k][r];
            }
            byPop.get(sum).add(k);
        }

        // Delay balance function: the multinomial factor F_Z(n). A class with
        // population but no delay demand makes the state infeasible.
        double[] g = new double[ngrid];
        for (int k = 0; k < ngrid; k++) {
            double logf = 0.0;
            boolean feas = true;
            for (int r = 0; r < R; r++) {
                if (counts[k][r] > 0) {
                    if (Z[r] <= 0) {
                        feas = false;
                        break;
                    }
                    logf += counts[k][r] * Math.log(Z[r]) - Gamma.logGamma(counts[k][r] + 1.0);
                }
            }
            if (feas) {
                g[k] = Math.exp(logf);
            }
        }

        // Convolve in one OI station at a time.
        int[] rem = new int[R];
        int[] y = new int[R];
        for (int m = 0; m < muList.size(); m++) {
            double[] vism = (visits != null && m < visits.length) ? visits[m] : null;
            double[] phi = balance(counts, byPop, strides, muList.get(m), R, ngrid, vism);
            double[] gnext = new double[ngrid];
            for (int kx = 0; kx < ngrid; kx++) {
                double px = phi[kx];
                if (px == 0.0) {
                    continue;
                }
                int base = 0;
                for (int r = 0; r < R; r++) {
                    rem[r] = N[r] - counts[kx][r];
                    y[r] = 0;
                    base += counts[kx][r] * strides[r];
                }
                // Odometer over the sub-box 0 <= y <= rem.
                while (true) {
                    int ylin = 0;
                    for (int r = 0; r < R; r++) {
                        ylin += y[r] * strides[r];
                    }
                    gnext[base + ylin] += px * g[ylin];
                    int d = 0;
                    while (d < R && y[d] == rem[d]) {
                        y[d] = 0;
                        d++;
                    }
                    if (d == R) {
                        break;
                    }
                    y[d]++;
                }
            }
            g = gnext;
        }

        // g now holds G(n) for EVERY n on the lattice, not just n = N.
        double G = g[ngrid - 1];
        double lG = (G > 0) ? Math.log(G) : Double.NEGATIVE_INFINITY;
        return new Ret.pfqnOiNc(G, lG, g);
    }

    /**
     * Balanced-fairness recursion Phi(n) = (1/mu(n)) sum_{r: n_r&gt;0} Phi(n-e_r)
     * over the count lattice, swept in population-increasing order.
     */
    private static double[] balance(int[][] counts, List<List<Integer>> byPop, int[] strides,
                                    ToDoubleFunction<int[]> murate, int R, int ngrid, double[] vis) {
        double[] phi = new double[ngrid];
        for (int s = 0; s < byPop.size(); s++) {
            List<Integer> bucket = byPop.get(s);
            for (int b = 0; b < bucket.size(); b++) {
                int k = bucket.get(b);
                if (s == 0) {
                    phi[k] = 1.0;
                    continue;
                }
                double mun = murate.applyAsDouble(counts[k]);
                if (!(mun > 0)) {
                    continue; // unreachable station state: zero balance value
                }
                double acc = 0.0;
                for (int r = 0; r < R; r++) {
                    if (counts[k][r] > 0) {
                        double vr = (vis == null) ? 1.0 : vis[r];
                        acc += vr * phi[k - strides[r]];
                    }
                }
                phi[k] = acc / mun;
            }
        }
        return phi;
    }
}
