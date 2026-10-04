/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.nc;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

import java.util.Random;

/**
 * Importance-sampling (IS) estimate of the normalizing constant of a closed
 * LOAD-DEPENDENT product-form queueing network. Load-dependent counterpart of
 * {@link Pfqn_pas_is} / {@link Pfqn_oi_is}: the same sample-an-ordering
 * estimator, with the order-independent rank rate replaced by the
 * load-dependent capacity.
 *
 * <p>Identity. Every product-form station's balance function is the sum, over
 * the orderings q of a given per-class count vector n, of an ordered product of
 * a per-position factor:
 * <pre>
 *   F_i(n) = |n|!/prod_r(n_r!) * prod_r L(i,r)^{n_r} / prod_{k=1}^{|n|} mu_i(k)
 *          = sum_{q: |q|=n} prod_{p=1}^{|n|} L(i,q_p) / mu_i(p),
 * </pre>
 * since the multiset has |n|!/prod_r(n_r!) orderings and each contributes the
 * same ordered product. The delay (infinite-server) node is the special case
 * mu_Z(k) = k, giving F_Z(n) = prod_r Z_r^{n_r}/n_r!; a single-server queue is
 * mu_i(k) = 1; a c-server queue is mu_i(k) = min(k,c).
 *
 * <p>Consequently, writing ell = sum(N) and letting a "cut vector" split an
 * ordering c of all ell jobs into S contiguous segments (one per station),
 * <pre>
 *   G(N) = sum_{c} sum_{cuts} prod_{m=1}^{S} w_m(seg_m),
 *   w_m(q) = prod_{p=1}^{|q|} L(m,q_p) / mu_m(p),
 * </pre>
 * because summing over the orderings of each segment independently reproduces
 * prod_m F_m(n_m), and each count split (n_1,...,n_S) is realized exactly once.
 *
 * <p>Estimator. An ordering c is drawn by placing, at each step, a uniformly
 * random present class; p(c) is the product of the reciprocal branching
 * factors. For the sampled c the inner sum over ALL cut vectors is computed
 * exactly by the dynamic program
 * <pre>
 *   A_0(0) = 1,   A_m(k) = sum_{j=0}^{k} A_{m-1}(j) * w_m(c_{j+1..k}),
 * </pre>
 * so S(c) = A_S(ell) in O(S*ell^2) time (no cut enumeration). Then
 * G = E_{C~p}[ S(C) / p(C) ] is unbiased and is estimated by the sample mean.
 *
 * <p>Port of matlab/src/api/pfqn/pfqn_ld_is.m.
 *
 * @see Pfqn_is
 * @see Pfqn_oi_is
 * @see Pfqn_pas_is
 */
public final class Pfqn_ld_is {
    private Pfqn_ld_is() {}

    /**
     * Importance-sampling estimate of the load-dependent normalizing constant.
     *
     * @param L       (M x R) per-class service demands at the M queueing stations.
     * @param N       (1 x R) closed population vector, finite.
     * @param Z       (1 x R) aggregated think time (delay) demand; null or zeros if none.
     * @param mu      (M x ell) load-dependent capacities, mu(i,k) the capacity of
     *                station i holding k jobs; null for the load-independent case
     *                mu(i,k)=1 (see {@link Pfqn_is}). Shorter rows are extended
     *                with their last capacity.
     * @param options solver options; uses options.samples (number of IS samples)
     *                and options.seed (RNG seed for reproducibility).
     * @return G and lG = log(G).
     */
    public static Ret.pfqnNc pfqn_ld_is(Matrix L, Matrix N, Matrix Z, Matrix mu, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        int[] Nv = new int[R];
        for (int r = 0; r < R; r++) {
            double nr = N.get(r);
            if (Double.isInfinite(nr) || Double.isNaN(nr)) {
                throw new RuntimeException("pfqn_ld_is requires finite (closed) populations.");
            }
            Nv[r] = (int) Math.round(nr);
        }
        if (N.length() != R) {
            throw new RuntimeException("L must have as many columns as N has classes.");
        }

        double[] Zv = new double[R];
        if (Z != null && !Z.isEmpty()) {
            // accept a (1 x R) row or an (M x R) block; aggregate over rows
            Matrix Zs = (Z.getNumRows() > 1) ? Z.sumCols() : Z;
            for (int r = 0; r < R; r++) {
                Zv[r] = Zs.get(r);
            }
        }

        int ell = 0;
        for (int r = 0; r < R; r++) {
            ell += Nv[r];
        }

        int nsamples = 10000;
        if (options != null && options.samples > 0) {
            nsamples = options.samples;
        }
        Random rng = (options != null) ? new Random(options.seed) : new Random();

        if (ell == 0) {
            return new Ret.pfqnNc(1.0, 0.0);
        }

        // ---- assemble the station list: M queues, plus the delay as mu_Z(k)=k --
        // D(m,r) is the per-class demand of station m; B(m,k) its capacity at k jobs.
        boolean hasZ = false;
        for (int r = 0; r < R; r++) {
            if (Zv[r] > 0) {
                hasZ = true;
                break;
            }
        }
        int S = M + (hasZ ? 1 : 0);
        double[][] D = new double[S][R];
        double[][] B = new double[S][ell];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                D[i][r] = L.get(i, r);
            }
            if (mu == null || mu.isEmpty()) {
                for (int k = 0; k < ell; k++) {
                    B[i][k] = 1.0;              // load-independent single server
                }
            } else {
                int muCols = mu.getNumCols();
                int lim = Math.min(ell, muCols);
                for (int k = 0; k < lim; k++) {
                    B[i][k] = mu.get(i, k);
                }
                for (int k = lim; k < ell; k++) {
                    B[i][k] = mu.get(i, muCols - 1);   // extend with the last capacity
                }
            }
        }
        if (hasZ) {
            for (int r = 0; r < R; r++) {
                D[S - 1][r] = Zv[r];
            }
            for (int k = 0; k < ell; k++) {
                B[S - 1][k] = k + 1;           // delay: mu_Z(k) = k
            }
        }
        for (int m = 0; m < S; m++) {
            for (int k = 0; k < ell; k++) {
                if (B[m][k] <= 0) {
                    throw new RuntimeException("load-dependent capacities must be strictly positive.");
                }
            }
        }

        double acc = 0.0;
        int[] x = new int[R];
        int[] c = new int[ell];
        int[] avail = new int[R];
        double[] A = new double[ell + 1];
        double[] Anew = new double[ell + 1];

        for (int s = 0; s < nsamples; s++) {
            // ---- draw an ordering c (uniformly random present class at each step)
            System.arraycopy(Nv, 0, x, 0, R);
            double logp = 0.0;
            for (int ppos = 0; ppos < ell; ppos++) {
                int na = 0;
                for (int r = 0; r < R; r++) {
                    if (x[r] > 0) {
                        avail[na] = r;
                        na++;
                    }
                }
                int pick = avail[rng.nextInt(na)];
                c[ppos] = pick;
                logp -= Math.log((double) na);
                x[pick]--;
            }

            // see _kb/03-api-layer.md for rationale
            java.util.Arrays.fill(A, 0.0);
            A[0] = 1.0;                        // A[k] indexes k jobs placed
            for (int m = 0; m < S; m++) {
                java.util.Arrays.fill(Anew, 0.0);
                for (int j = 0; j <= ell; j++) {
                    if (A[j] == 0.0) {
                        continue;
                    }
                    double w = 1.0;
                    Anew[j] += A[j];                             // empty segment
                    for (int k = j + 1; k <= ell; k++) {
                        w = w * D[m][c[k - 1]] / B[m][k - j - 1];   // position within segment
                        if (w == 0.0) {
                            break;
                        }
                        Anew[k] += A[j] * w;
                    }
                }
                System.arraycopy(Anew, 0, A, 0, ell + 1);
            }
            acc += A[ell] * Math.exp(-logp);
        }

        double G = acc / nsamples;
        return new Ret.pfqnNc(G, Math.log(G));
    }
}
