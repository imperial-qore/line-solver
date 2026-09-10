package jline.api.sn;

import java.util.function.Function;

/**
 * Whittle balance check for a globally state-dependent rate scaling.
 *
 * <p>For every state n of the lattice 0..cutoffs and every pair of stations
 * (s,t) populated in n, the balance property requires
 *
 * <pre>  phi_s(n) phi_t(n - e_s) = phi_t(n) phi_s(n - e_t).</pre>
 *
 * <p>When it holds, the chain is reversible with pi(n) ~ Phi(n) prod rho^n for
 * the balance function Phi implied by phi, and the stationary law is insensitive
 * to the service-time distribution beyond its mean. When it fails, the model is
 * still solvable by SolverCTMC but has no product form and is sensitive.
 *
 * <p>Twin of MATLAB matlab/src/api/sn/sn_gd_balance.m and of python
 * line_solver.api.sn.sn_gd_balance.
 *
 * <p>Reference: P. Whittle, "Partial balance and insensitivity", J. Appl. Prob.
 * 22(1), 1985; T. Bonald, A. Proutiere, "Insensitivity in processor-sharing
 * networks", Perf. Eval. 49, 2002.
 */
public class SnGdBalance {

    private SnGdBalance() {
    }

    /**
     * Worst relative violation of the balance property over the given lattice.
     *
     * @param phi     scaling evaluated on an (nstations) population vector, returning
     *                a scalar (broadcast) or one entry per station
     * @param cutoffs per-station lattice bound, one entry per station
     * @return the worst relative violation; 0 (to rounding) when phi is balanced
     */
    public static double sn_gd_balance(Function<double[], double[]> phi, int[] cutoffs) {
        if (phi == null) {
            throw new IllegalArgumentException("phi must not be null.");
        }
        final int S = cutoffs.length;
        if (S < 2) {
            throw new IllegalArgumentException(
                    "cutoffs must have one entry per station (at least two stations are needed for a balance pair).");
        }
        long total = 1;
        for (int s = 0; s < S; s++) total *= (cutoffs[s] + 1);

        double viol = 0;
        double[] n = new double[S];
        for (long idx = 0; idx < total; idx++) {
            long rem = idx;
            for (int s = 0; s < S; s++) {
                n[s] = rem % (cutoffs[s] + 1);
                rem /= (cutoffs[s] + 1);
            }
            for (int s = 0; s < S; s++) {
                if (n[s] == 0) continue;
                for (int t = s + 1; t < S; t++) {
                    if (n[t] == 0) continue;
                    final double[] xn = eval(phi, n, S);
                    n[s] -= 1;
                    final double[] xs = eval(phi, n, S);
                    n[s] += 1;
                    n[t] -= 1;
                    final double[] xt = eval(phi, n, S);
                    n[t] += 1;
                    final double lhs = xn[s] * xs[t];
                    final double rhs = xn[t] * xt[s];
                    final double scale = Math.max(Math.abs(lhs), Math.abs(rhs));
                    if (scale > 0) {
                        final double v = Math.abs(lhs - rhs) / scale;
                        if (v > viol) viol = v;
                    }
                }
            }
        }
        return viol;
    }

    private static double[] eval(Function<double[], double[]> phi, double[] n, int S) {
        final double[] v = phi.apply(n.clone());
        if (v.length == 1) {
            final double[] out = new double[S];
            for (int s = 0; s < S; s++) out[s] = v[0];
            return out;
        }
        if (v.length != S) {
            throw new IllegalArgumentException("phi must return a scalar or a vector of length " + S + ".");
        }
        return v;
    }
}
