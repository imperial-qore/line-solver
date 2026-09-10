/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.nc;

import java.util.List;
import java.util.Random;
import java.util.function.ToDoubleFunction;

/**
 * Importance-sampling (IS) estimate of the normalizing constant of a SINGLE
 * communicating class of a closed two-station pass-and-swap (P&S) /
 * order-independent (OI) tandem with swap graph H (Casale, Comte and Dorsman,
 * 2026). Monte-Carlo counterpart of the exact convolution: it estimates the
 * same per-communicating-class constant G_C but scales to populations where the
 * exact count-lattice convolution becomes expensive.
 *
 * <p>Two OI/P&S stations (1 = upstream, 2 = downstream) hold all N jobs. With a
 * non-empty swap graph the ordered-state chain is reducible; the recurrent
 * communicating class is the set of splits of the orderings that are the linear
 * extensions of the placement partial order induced by H. Writing D for that
 * set of orderings and Phi_m for the balanced-fairness balance of station m,
 *   G_C = sum_{c in D} sum_{k=0..ell} Phi_1(c_{1..k}) Phi_2(c_{ell..k+1}),
 * with Phi_m(q) = prod_p 1/mu_m(n(q_{1..p})), n(.) the per-class COUNT vector of
 * the prefix (OI property P1 makes mu permutation-invariant, i.e. a function of the
 * counts -- NOT of the support alone, which differs once a class holds two jobs).
 *
 * <p>Auto-normalized IS (notebook generator IS_3): orderings c are drawn from D
 * by placing, at each step, a uniformly random placement-order-minimal present
 * class; p(c) is the product of reciprocal branching factors. Then
 *   G_C[xi] = E_{C~p}[ (sum_k xi(C,k) Phi_1 Phi_2) / p(C) ]
 * and E[xi] = G_C[xi]/G_C[1] reuses the same samples for numerator and
 * denominator. xi = number of class-r jobs in the prefix gives the mean queue
 * length of class r at station 1.
 *
 * <p>Port of matlab/src/api/pfqn/pfqn_pas_is.m.
 */
public final class Pfqn_pas_is {
    private Pfqn_pas_is() {}

    /** Result of {@link #pfqn_pas_is}: G, log G, and (2 x R) mean queue lengths. */
    public static final class Result {
        public final double G;
        public final double lG;
        public final double[][] Q;   // Q[0] station 1, Q[1] station 2

        public Result(double G, double lG, double[][] Q) {
            this.G = G;
            this.lG = lG;
            this.Q = Q;
        }
    }

    /**
     * @param N        (1 x R) closed population vector (macrostate).
     * @param mu       list of exactly two OI rank-rate handles; mu.get(m)(n)
     *                 returns the total service rate of station m for the
     *                 per-class count vector n (permutation-invariant, but not
     *                 in general a function of supp(n) alone).
     * @param H        (R x R) placement-order DAG; H[i][j] != 0 iff i precedes j.
     *                 Null/all-zero => pure OI (every ordering feasible).
     * @param nsamples number of IS samples.
     * @param seed     RNG seed (for reproducibility / common random numbers).
     * @param verbose  print progress when true.
     */
    public static Result pfqn_pas_is(int[] N, List<ToDoubleFunction<int[]>> mu, int[][] H,
                                    long nsamples, long seed, boolean verbose) {
        return pfqn_pas_is(N, mu, H, nsamples, seed, verbose, true);
    }

    /**
     * As {@link #pfqn_pas_is(int[], List, int[][], long, long, boolean)}, with
     * the queue-length coefficients made optional.
     *
     * @param wantQlen estimate the per-class queue lengths as well as the
     *                 constant. False estimates ONLY G: the xi = n_{1,r}
     *                 coefficients are not accumulated and Q comes back zero.
     *                 The ordering is drawn from the same stream either way, so
     *                 G is unchanged to the last bit -- this is for the callers
     *                 that want G(N - e_r) and read nothing else from it.
     */
    public static Result pfqn_pas_is(int[] N, List<ToDoubleFunction<int[]>> mu, int[][] H,
                                    long nsamples, long seed, boolean verbose, boolean wantQlen) {
        if (mu == null || mu.size() != 2) {
            throw new RuntimeException("pfqn_pas_is models a two-station pass-and-swap tandem: mu must have exactly two rate functions.");
        }
        int R = N.length;
        for (int r = 0; r < R; r++) {
            if (!Double.isFinite(N[r])) {
                throw new RuntimeException("pfqn_pas_is requires finite (closed) populations.");
            }
        }
        int ell = 0;
        for (int r = 0; r < R; r++) {
            ell += N[r];
        }

        // Placement-order logic isolated in Pas_placement.
        int[][] P = Pas_placement.closure(H);

        int nCoef = wantQlen ? R + 1 : 1;           // xi = [1, n_{1,1}, ..., n_{1,R}]
        double[] Q0 = new double[R];
        double[][] Q = new double[2][R];
        if (ell == 0) {
            return new Result(1.0, 0.0, Q);
        }

        final ToDoubleFunction<int[]> mu1 = mu.get(0);
        final ToDoubleFunction<int[]> mu2 = mu.get(1);
        Random rng = new Random(seed);

        double[] accum = new double[nCoef];
        int[] x = new int[R];
        int[] c = new int[ell];
        int[] occ2 = new int[R];
        double[] Phi1 = new double[ell + 1];
        double[] Phi2cut = new double[ell + 1];
        int[][] cnt1 = new int[ell + 1][R];
        long report = Math.max(1L, nsamples / 10L);

        for (long s = 0; s < nsamples; s++) {
            // ---- draw an ordering c from D (auto-normalized IS, generator IS_3)
            System.arraycopy(N, 0, x, 0, R);
            double logp = 0.0;
            for (int ppos = 0; ppos < ell; ppos++) {
                int[] avail = Pas_placement.placeable(x, P);
                int na = avail.length;
                if (na == 0) {
                    throw new RuntimeException("swap graph induces no feasible ordering (cyclic placement order).");
                }
                int pick = avail[rng.nextInt(na)];
                c[ppos] = pick;
                logp -= Math.log(na);
                x[pick]--;
            }
            double p_c = Math.exp(logp);

            // ---- prefix balance Phi_1 and class counts up to each cut ----------
            for (int r = 0; r < R; r++) {
                cnt1[0][r] = 0;
            }
            Phi1[0] = 1.0;
            double phi = 1.0;
            for (int k = 1; k <= ell; k++) {
                int cls = c[k - 1];
                for (int r = 0; r < R; r++) {
                    cnt1[k][r] = cnt1[k - 1][r];
                }
                cnt1[k][cls]++;
                // the rank rate is evaluated at the prefix COUNT vector, not at its
                // support: OI property P1 makes mu permutation-invariant, i.e. a
                // function of the counts (an INF station has mu(n)=sum_r n_r sigma_r)
                phi /= mu1.applyAsDouble(cnt1[k]);
                Phi1[k] = phi;
            }

            // see _kb/03-api-layer.md for rationale
            for (int r = 0; r < R; r++) {
                occ2[r] = 0;
            }
            phi = 1.0;
            Phi2cut[ell] = 1.0;   // placeholder; overwritten below for indices 1..ell
            for (int k = ell; k >= 1; k--) {
                int cls = c[k - 1];
                occ2[cls]++;
                phi /= mu2.applyAsDouble(occ2);
                Phi2cut[k] = phi;
            }

            // sample_val(coef) = sum_{k=0..ell} xi_coef(k) Phi1[k] Phi2cut-at-k
            for (int k = 0; k <= ell; k++) {
                double phi2 = (k >= ell) ? 1.0 : Phi2cut[k + 1];
                double w = Phi1[k] * phi2;
                accum[0] += w / p_c;
                if (wantQlen && k > 0) {
                    for (int r = 0; r < R; r++) {
                        if (cnt1[k][r] != 0) {
                            accum[1 + r] += (w * cnt1[k][r]) / p_c;
                        }
                    }
                }
            }

            if (verbose && ((s + 1) % report == 0)) {
                System.out.println("pfqn_pas_is: " + (s + 1) + "/" + nsamples + " samples");
            }
        }

        double G = accum[0] / nsamples;
        double lG = (G > 0) ? Math.log(G) : Double.NEGATIVE_INFINITY;
        if (wantQlen) {
            if (G > 0) {
                for (int r = 0; r < R; r++) {
                    Q0[r] = (accum[1 + r] / nsamples) / G;
                }
            }
            for (int r = 0; r < R; r++) {
                Q[0][r] = Q0[r];
                Q[1][r] = N[r] - Q0[r];
            }
        }
        return new Result(G, lG, Q);
    }
}
