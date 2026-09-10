/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.retrieval;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;

/**
 * Exact MVA-style recursion for delayed-hit (list-based) cache metrics.
 *
 * <p>Port of {@code matlab/src/api/retrieval/retrieval_mva.m}, twin of
 * {@code cpp/include/line/api/retrieval/retrieval_mva.h} and of the native
 * Python {@code api.retrieval.retrieval_mva}.
 *
 * <p>This is the exact recursion that {@link Retrieval_fpi} approximates.
 * Writing phi^(k) for the delayed-hit probability in the system WITHOUT item k,
 *
 * <pre>
 *   theta_ij(m) = gamma_ij / (1 + lambda_i eta_0i
 *                   + sum_s lambda_i eta_si (1 + sum_{k!=i} phi^(i)_sk(m - 1_j)))
 *   xi_j(m)     = m_j / sum_i theta_ij(m) (1 - pihit_i(m - 1_j))
 *   pi_ij(m)    = theta_ij(m) xi_j(m) (1 - pihit_i(m - 1_j))
 *   pi_i0(m)    = (1 - pihit_i(m)) / (1 + lambda_i eta_0i
 *                   + sum_s lambda_i eta_si (1 + sum_{k!=i} phi^(i)_sk(m)))
 *   phi_sk(m)   = lambda_k pi_k0(m) eta_sk (1 + sum_{i!=k} phi^(k)_si(m))
 *   phi_0k(m)   = lambda_k eta_0k pi_k0(m)
 * </pre>
 *
 * memoized over (item subset, capacity vector). The recursion bottoms out at the
 * empty item set and at any capacity able to hold every remaining item, where
 * the items are permanently cached: hit probability one, no miss and no fetch.
 *
 * <p>Cost is O(2^n n^2 h r prod_j (1+m_j)) in time and memory, so this is a
 * small-case ORACLE; use {@link Retrieval_fpi} beyond that.
 */
public final class Retrieval_mva {

    private Retrieval_mva() {
    }

    /** Mirrors the [pmiss, phit, pdh] return list of the MATLAB function. */
    public static final class Result {
        /** 1 x n, miss ratios pi_i0. */
        public Matrix pmiss;
        /** h x n, hit ratios pi_ij. */
        public Matrix phit;
        /** (r+1) x n, delayed-hit probabilities phi_si, s = 0..r. */
        public Matrix pdh;
    }

    /**
     * @param m      1 x h cache list capacities
     * @param lambda 1 x n per-item arrival rates
     * @param eta    n x (r+1) fetching demands; column 0 is the IS station,
     *               columns 1..r the PS stations
     * @param gamma  n x h access factors
     */
    public static Result retrieval_mva(int[] m, double[] lambda, Matrix eta, Matrix gamma) {
        int n = lambda.length;
        int h = m.length;
        if (eta.getNumRows() != n || gamma.getNumRows() != n) {
            line_error("retrieval_mva", "eta/gamma and lambda disagree on the item count.");
        }
        if (gamma.getNumCols() != h) {
            line_error("retrieval_mva", "gamma and m disagree on the number of lists.");
        }
        if (eta.getNumCols() == 0) {
            line_error("retrieval_mva", "eta has no columns.");
        }
        for (int j = 0; j < h; j++) {
            if (m[j] < 0) {
                line_error("retrieval_mva", "negative list capacity.");
            }
        }
        if (n > 30) {
            line_error("retrieval_mva", "too many items for a bitmask subset enumeration: " + n);
        }
        return new Solver(m, lambda, eta, gamma).run();
    }

    /**
     * Memo tables and the recursion. Held in an object rather than in statics so
     * the function is reentrant and leaves no state behind.
     */
    private static final class Solver {
        private final int[] m;
        private final double[] lambda;
        private final Matrix eta;
        private final Matrix gamma;
        private final int n;
        private final int h;
        private final int r;
        private final int[] radix;
        private final int ncap;
        private final int nmask;
        private final boolean[] done;
        private final double[] pi0;
        private final double[] pihit;
        private final double[] pij;
        private final double[] phi;

        Solver(int[] m, double[] lambda, Matrix eta, Matrix gamma) {
            this.m = m;
            this.lambda = lambda;
            this.eta = eta;
            this.gamma = gamma;
            this.n = lambda.length;
            this.h = m.length;
            this.r = eta.getNumCols() - 1;
            this.radix = new int[h];
            int cap = 1;
            for (int j = 0; j < h; j++) {
                radix[j] = m[j] + 1;
                cap *= radix[j];
            }
            this.ncap = cap;
            this.nmask = 1 << n;
            int cells = nmask * ncap;
            this.done = new boolean[cells];
            this.pi0 = new double[cells * n];
            this.pihit = new double[cells * n];
            this.pij = new double[cells * n * h];
            this.phi = new double[cells * n * (r + 1)];
        }

        Result run() {
            int full = nmask - 1;
            solve(full, m.clone());
            int ci = capidx(m);
            Result out = new Result();
            out.pmiss = new Matrix(1, n);
            out.phit = new Matrix(h, n);
            out.pdh = new Matrix(r + 1, n);
            for (int i = 0; i < n; i++) {
                out.pmiss.set(0, i, pi0[(full * ncap + ci) * n + i]);
                for (int j = 0; j < h; j++) {
                    out.phit.set(j, i, pij[((full * ncap + ci) * n + i) * h + j]);
                }
                for (int s = 0; s <= r; s++) {
                    out.pdh.set(s, i, phi[((full * ncap + ci) * n + i) * (r + 1) + s]);
                }
            }
            return out;
        }

        private int capidx(int[] c) {
            int idx = 0;
            int mul = 1;
            for (int j = 0; j < h; j++) {
                idx += c[j] * mul;
                mul *= radix[j];
            }
            return idx;
        }

        private double getPi0(int mask, int ci, int i) {
            return pi0[(mask * ncap + ci) * n + i];
        }

        private double getPihit(int mask, int ci, int i) {
            return pihit[(mask * ncap + ci) * n + i];
        }

        private double getPhi(int mask, int ci, int i, int s) {
            return phi[((mask * ncap + ci) * n + i) * (r + 1) + s];
        }

        private void solve(int mask, int[] c) {
            int ci = capidx(c);
            if (done[mask * ncap + ci]) {
                return;
            }
            if (mask == 0) {
                done[mask * ncap + ci] = true;
                return;
            }
            List<Integer> active = new ArrayList<Integer>();
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    active.add(Integer.valueOf(i));
                }
            }
            long csum = 0;
            for (int j = 0; j < h; j++) {
                csum += c[j];
            }
            if (csum >= active.size()) {
                // The cache holds every active item, so all of them are cached
                // permanently: no miss, no fetch, hit probability one.
                for (int idx = 0; idx < active.size(); idx++) {
                    int i = active.get(idx).intValue();
                    pihit[(mask * ncap + ci) * n + i] = 1.0;
                }
                done[mask * ncap + ci] = true;
                return;
            }

            // dependencies at m - 1_j, with and without each active item
            for (int j = 0; j < h; j++) {
                if (c[j] > 0) {
                    int[] cj = c.clone();
                    cj[j] -= 1;
                    solve(mask, cj);
                    for (int idx = 0; idx < active.size(); idx++) {
                        solve(mask & ~(1 << active.get(idx).intValue()), cj);
                    }
                }
            }
            for (int idx = 0; idx < active.size(); idx++) {
                solve(mask & ~(1 << active.get(idx).intValue()), c);
            }

            // theta, xi and pi_ij, all evaluated on the cache recursion m - 1_j
            for (int j = 0; j < h; j++) {
                if (c[j] == 0) {
                    continue;
                }
                int[] cj = c.clone();
                cj[j] -= 1;
                int cjx = capidx(cj);
                double[] theta = new double[n];
                for (int idx = 0; idx < active.size(); idx++) {
                    int i = active.get(idx).intValue();
                    int maski = mask & ~(1 << i);
                    double acc = 0.0;
                    for (int s = 0; s < r; s++) {
                        double sphi = 0.0;
                        for (int q = 0; q < active.size(); q++) {
                            int k = active.get(q).intValue();
                            if (k != i) {
                                sphi += getPhi(maski, cjx, k, s + 1);
                            }
                        }
                        acc += lambda[i] * eta.get(i, s + 1) * (1.0 + sphi);
                    }
                    theta[i] = gamma.get(i, j) / (1.0 + lambda[i] * eta.get(i, 0) + acc);
                }
                double sden = 0.0;
                for (int idx = 0; idx < active.size(); idx++) {
                    int i = active.get(idx).intValue();
                    sden += theta[i] * (1.0 - getPihit(mask, cjx, i));
                }
                if (sden == 0.0) {
                    line_error("retrieval_mva",
                            "degenerate list occupancy (zero denominator).");
                }
                double xij = ((double) c[j]) / sden;
                for (int idx = 0; idx < active.size(); idx++) {
                    int i = active.get(idx).intValue();
                    pij[((mask * ncap + ci) * n + i) * h + j] =
                            theta[i] * xij * (1.0 - getPihit(mask, cjx, i));
                }
            }

            for (int idx = 0; idx < active.size(); idx++) {
                int i = active.get(idx).intValue();
                double ph = 0.0;
                for (int j = 0; j < h; j++) {
                    ph += pij[((mask * ncap + ci) * n + i) * h + j];
                }
                pihit[(mask * ncap + ci) * n + i] = ph;
            }

            for (int idx = 0; idx < active.size(); idx++) {
                int i = active.get(idx).intValue();
                int maski = mask & ~(1 << i);
                double acc = 0.0;
                for (int s = 0; s < r; s++) {
                    double sphi = 0.0;
                    for (int q = 0; q < active.size(); q++) {
                        int k = active.get(q).intValue();
                        if (k != i) {
                            sphi += getPhi(maski, ci, k, s + 1);
                        }
                    }
                    acc += lambda[i] * eta.get(i, s + 1) * (1.0 + sphi);
                }
                pi0[(mask * ncap + ci) * n + i] =
                        (1.0 - getPihit(mask, ci, i)) / (1.0 + lambda[i] * eta.get(i, 0) + acc);
            }

            for (int idx = 0; idx < active.size(); idx++) {
                int k = active.get(idx).intValue();
                int maskk = mask & ~(1 << k);
                double pi0k = getPi0(mask, ci, k);
                phi[((mask * ncap + ci) * n + k) * (r + 1)] = lambda[k] * eta.get(k, 0) * pi0k;
                for (int s = 0; s < r; s++) {
                    double sphi = 0.0;
                    for (int q = 0; q < active.size(); q++) {
                        int i = active.get(q).intValue();
                        if (i != k) {
                            sphi += getPhi(maskk, ci, i, s + 1);
                        }
                    }
                    phi[((mask * ncap + ci) * n + k) * (r + 1) + s + 1] =
                            lambda[k] * pi0k * eta.get(k, s + 1) * (1.0 + sphi);
                }
            }

            done[mask * ncap + ci] = true;
        }
    }
}
