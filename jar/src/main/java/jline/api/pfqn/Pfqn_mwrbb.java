/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn;

import jline.util.matrix.Matrix;

/**
 * Majumdar-Woodside robust box bounds on throughput for closed multiclass
 * queueing networks with mixed scheduling disciplines.
 *
 * <p>Computes distribution-insensitive (NBUE) upper and lower bounds on the
 * per-class system throughput of a closed multiclass queueing network, per
 * S. Majumdar and C.M. Woodside, "Robust bounds and throughput guarantees for
 * closed multiclass queueing networks", Performance Evaluation 32 (1998)
 * 101-136. The upper bound intersects the no-contention bound (eq. 2) with the
 * utilization-based bound (eq. 3) and is independent of the scheduling
 * discipline. The lower bound is the multiclass throughput guarantee of
 * Theorem 2 (eq. 15): X_c &gt;= N_c / (Z_c + sum_k V_kc (S_kc + d_kc+)), where
 * the per-visit queueing delay bound d_kc+ depends on the discipline at
 * station k -- FIFO (Theorem 1 / eq. 6, via Lemma 1), processor sharing
 * (Lemma 2), preemptive priority (Lemma 3) and non-preemptive priority
 * (Lemmas 4-5). The coupled inequalities are resolved by the interval-
 * narrowing fixed point reproducing the BNR-Prolog robust box bounds; for a
 * single FIFO class it reduces to the Muntz-Wong asymptotic bounds.</p>
 *
 * <p>Bounds are insensitive to the service-time distributions (only NBUE is
 * assumed) and to routing dependencies. The think time Z aggregates the
 * pure-delay (infinite-server) stations; only queueing stations are passed in
 * V,S.</p>
 *
 * @since LINE 3.0
 */
public final class Pfqn_mwrbb {

    /** Discipline codes for the {@code sched} argument. */
    public static final int FIFO = 0;
    public static final int PS = 1;
    public static final int NPPRIO = 2;   // non-preemptive priority
    public static final int PPPRIO = 3;   // preemptive priority
    public static final int ABA = 4;      // ABA full-contention (discipline-independent)

    private Pfqn_mwrbb() {
    }

    /** Result holder for {@link #pfqn_mwrbb}. */
    public static final class Result {
        /** (1 x C) lower bound on class throughput (Theorem 2). */
        public final Matrix Xlo;
        /** (1 x C) upper bound on class throughput (eqs. 2-3). */
        public final Matrix Xup;
        /** (K x C) per-visit residence time consistent with the lower bound. */
        public final Matrix Wlo;

        public Result(Matrix Xlo, Matrix Xup, Matrix Wlo) {
            this.Xlo = Xlo;
            this.Xup = Xup;
            this.Wlo = Wlo;
        }
    }

    /** FIFO, equal-priority convenience overload. */
    public static Result pfqn_mwrbb(Matrix V, Matrix S, Matrix N, Matrix Z) {
        return pfqn_mwrbb(V, S, N, Z, null, null);
    }

    /**
     * @param V     (K x C) mean visits of class c at queueing station k
     * @param S     (K x C) mean service demand per visit of class c at station k
     * @param N     (1 x C) or (C x 1) population of class c
     * @param Z     (1 x C) or (C x 1) think time of class c (may be null)
     * @param sched (K x 1) discipline code per station (0=FIFO, 1=PS,
     *              2=non-preemptive priority, 3=preemptive priority); null =&gt; all FIFO
     * @param prio  (1 x C) class priority, lower value = higher priority; null =&gt; all equal
     * @return robust box bounds and per-visit residence times
     */
    public static Result pfqn_mwrbb(Matrix V, Matrix S, Matrix N, Matrix Z,
                                    Matrix sched, Matrix prio) {
        int K = V.getNumRows();
        int C = V.getNumCols();

        double[] n = new double[C];
        double[] z = new double[C];
        double[] pr = new double[C];
        int[] sc = new int[K];
        for (int c = 0; c < C; c++) {
            n[c] = (N.getNumRows() == 1) ? N.get(0, c) : N.get(c, 0);
            if (Z == null || Z.isEmpty()) {
                z[c] = 0.0;
            } else {
                z[c] = (Z.getNumRows() == 1) ? Z.get(0, c) : Z.get(c, 0);
            }
            if (prio == null || prio.isEmpty()) {
                pr[c] = 0.0;
            } else {
                pr[c] = (prio.getNumRows() == 1) ? prio.get(0, c) : prio.get(c, 0);
            }
        }
        for (int k = 0; k < K; k++) {
            sc[k] = (sched == null || sched.isEmpty()) ? FIFO
                    : (int) Math.round((sched.getNumRows() == 1) ? sched.get(0, k) : sched.get(k, 0));
        }

        double[][] Va = to2d(V, K, C);
        double[][] Sa = to2d(S, K, C);

        // no-contention upper bound on the cycle rate f_c = X_c/N_c (eqs. 1-2)
        double[] fup = new double[C];
        double[] flo = new double[C];
        for (int c = 0; c < C; c++) {
            double dem = z[c];
            for (int k = 0; k < K; k++) {
                dem += Va[k][c] * Sa[k][c];
            }
            fup[c] = 1.0 / dem;
        }

        int maxiter = 20000;
        double tol = 1e-13;
        for (int it = 0; it < maxiter; it++) {
            double maxdelta = 0.0;

            // utilization-based narrowing of the upper bounds (eq. 3)
            for (int c = 0; c < C; c++) {
                double cap = fup[c];
                for (int k = 0; k < K; k++) {
                    double other = 0.0;
                    for (int m = 0; m < C; m++) {
                        if (m != c) {
                            other += n[m] * Va[k][m] * Sa[k][m] * flo[m];
                        }
                    }
                    double denomk = n[c] * Va[k][c] * Sa[k][c];
                    if (denomk > 0) {
                        cap = Math.min(cap, (1.0 - other) / denomk);
                    }
                }
                double newfup = Math.min(fup[c], Math.max(cap, 0.0));
                maxdelta = Math.max(maxdelta, Math.abs(newfup - fup[c]));
                fup[c] = newfup;
            }

            // lower-bound narrowing (Theorem 2, eq. 15). Higher-priority delay
            // carries a 1/f_c factor and is isolated as Bh: f_c = (1-Bh)/DEN.
            for (int c = 0; c < C; c++) {
                double[] db = denom(c, Va, Sa, n, z, fup, flo, sc, pr, K, C);
                double val = (1.0 - db[1]) / db[0];
                if (val < 0) {
                    val = 0.0;
                }
                double newflo = Math.max(flo[c], val);
                maxdelta = Math.max(maxdelta, Math.abs(newflo - flo[c]));
                flo[c] = newflo;
            }

            if (maxdelta < tol) {
                break;
            }
        }

        Matrix Xlo = new Matrix(1, C);
        Matrix Xup = new Matrix(1, C);
        Matrix Wlo = new Matrix(K, C);
        for (int c = 0; c < C; c++) {
            Xlo.set(0, c, n[c] * flo[c]);
            Xup.set(0, c, n[c] * fup[c]);
            for (int k = 0; k < K; k++) {
                if (Va[k][c] == 0) {
                    continue;
                }
                Wlo.set(k, c, residence(k, c, Va, Sa, n, fup, flo, sc, pr, C));
            }
        }
        return new Result(Xlo, Xup, Wlo);
    }

    // lower-bound denominator DEN (index 0) and isolated higher-priority work
    // Bh (index 1) for class c.
    private static double[] denom(int c, double[][] V, double[][] S, double[] n,
                                  double[] z, double[] fup, double[] flo,
                                  int[] sc, double[] pr, int K, int C) {
        double DEN = z[c];
        double Bh = 0.0;
        double fc = flo[c];
        for (int k = 0; k < K; k++) {
            double Vkc = V[k][c];
            if (Vkc == 0) {
                continue;
            }
            if (sc[k] == PPPRIO || sc[k] == NPPRIO) {
                for (int m = 0; m < C; m++) {
                    if (pr[m] < pr[c]) {
                        Bh += n[m] * fup[m] * V[k][m] * S[k][m];
                    }
                }
            }
            DEN += Vkc * wrest(k, c, V, S, n, fup, fc, sc, pr, C);
        }
        return new double[]{DEN, Bh};
    }

    // per-visit residence at station k for class c EXCLUDING the isolated
    // higher-priority (1/f_c) term; includes own service and bounded delays.
    private static double wrest(int k, int c, double[][] V, double[][] S,
                                double[] n, double[] fup, double fc,
                                int[] sc, double[] pr, int C) {
        double Vkc = V[k][c];
        double Skc = S[k][c];
        int d = sc[k];
        if (d == FIFO) {
            double s = 0.0;
            for (int m = 0; m < C; m++) {
                double pcm = (fc * Vkc == 0) ? 1.0
                        : Math.min(1.0, (fup[m] * V[k][m]) / (fc * Vkc));
                s += n[m] * S[k][m] * pcm;
            }
            return s;   // own service is the m=c term (= N_c S_kc)
        } else if (d == PS) {
            double dp = 0.0;
            for (int m = 0; m < C; m++) {
                double Ncont = (m == c) ? n[c] - 1 : n[m];
                double term = (fc * Vkc == 0) ? Skc
                        : Math.min(Skc, (fup[m] * V[k][m] * S[k][m]) / (fc * Vkc));
                dp += Ncont * term;
            }
            return Skc + dp;
        } else if (d == ABA) {   // ABA full-contention (P_cm = 1)
            double s = 0.0;
            for (int m = 0; m < C; m++) {
                s += n[m] * S[k][m];   // wait behind full service of all customers
            }
            return s;
        } else {   // NPPRIO or PPPRIO
            double dp = 0.0;
            for (int m = 0; m < C; m++) {
                if (pr[m] == pr[c]) {   // equal priority (includes c)
                    double Ncont = (m == c) ? n[c] - 1 : n[m];
                    double pcm = (fc * Vkc == 0) ? 1.0
                            : Math.min(1.0, (fup[m] * V[k][m]) / (fc * Vkc));
                    dp += Ncont * S[k][m] * pcm;
                }
                // higher priority handled via Bh in denom()
            }
            if (d == NPPRIO) {   // lower-priority water-filling (Lemma 4)
                Integer[] L = lowerSortedByService(k, c, S, pr, C);
                double budget = 1.0;
                for (int idx = 0; idx < L.length; idx++) {
                    int l = L[idx];
                    double capr = (fc * Vkc == 0) ? Double.POSITIVE_INFINITY
                            : (fup[l] * V[k][l]) / (fc * Vkc);
                    double al = (n[l] <= 0) ? 0.0 : Math.min(budget / n[l], capr);
                    if (al < 0) {
                        al = 0.0;
                    }
                    dp += n[l] * al * S[k][l];
                    budget -= n[l] * al;
                    if (budget < 0) {
                        budget = 0.0;
                    }
                }
            }
            return Skc + dp;
        }
    }

    // full per-visit residence (including higher-priority delay) for Q.
    private static double residence(int k, int c, double[][] V, double[][] S,
                                    double[] n, double[] fup, double[] flo,
                                    int[] sc, double[] pr, int C) {
        double fc = flo[c];
        double Vkc = V[k][c];
        double W = wrest(k, c, V, S, n, fup, fc, sc, pr, C);
        int d = sc[k];
        if ((d == NPPRIO || d == PPPRIO) && fc * Vkc > 0) {
            for (int m = 0; m < C; m++) {
                if (pr[m] < pr[c]) {
                    W += n[m] * fup[m] * V[k][m] * S[k][m] / (fc * Vkc);
                }
            }
        }
        return W;
    }

    // lower-priority classes (prio > prio(c)) sorted by decreasing service.
    private static Integer[] lowerSortedByService(int k, int c, double[][] S,
                                                  double[] pr, int C) {
        java.util.List<Integer> L = new java.util.ArrayList<Integer>();
        for (int m = 0; m < C; m++) {
            if (pr[m] > pr[c]) {
                L.add(m);
            }
        }
        final int kk = k;
        final double[][] Sf = S;
        L.sort(new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Double.compare(Sf[kk][b], Sf[kk][a]);
            }
        });
        return L.toArray(new Integer[0]);
    }

    private static double[][] to2d(Matrix M, int rows, int cols) {
        double[][] a = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                a[i][j] = M.get(i, j);
            }
        }
        return a;
    }
}
