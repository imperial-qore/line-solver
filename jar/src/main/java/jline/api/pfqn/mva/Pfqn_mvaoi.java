/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

import jline.api.pfqn.ld.Pfqn_oi_insvc;

/**
 * Mean-value analysis of a closed product-form queueing network composed of an
 * aggregated infinite-server (delay) node, any number of load-independent (LI)
 * single-server product-form queues, and any number of order-independent (OI) /
 * pass-and-swap stations with empty swap graph.
 *
 * This is the mean-value counterpart of {@link jline.api.pfqn.nc.Pfqn_ncoi} and
 * the marginal-distribution form {@link Pfqn_mvaoi_marg}: it returns the same
 * exact per-class throughput and queue-lengths but WITHOUT computing any
 * normalizing constant or joint marginal, using only mean quantities. It is the
 * composition-dependent generalization of the Conditional MVA (CMVA) of Casale,
 * "A Note on Stable Flow-Equivalent Aggregation in Closed Networks" (QUESTA
 * 2009), extended to MULTIPLE OI stations by carrying one rate-shift vector s_i
 * per OI station i (row i of the shift matrix S).
 *
 * Throughout, r and s index job classes; i indexes OI stations; j indexes LI
 * queues. State (S, Nn) is processed by increasing sum(Nn); each OI station keeps
 * its own D^i, rho^i and Q^i recursions driven by the common throughput
 * X^{(S)}(Nn), and the population conservation aggregates every station's
 * contribution:
 *   Nn_r = X_r Z_r + sum_j Q^{(j)}_r + sum_i Q^{(i)}_r,
 * with the LI queue term Q^{(j)}_r = X_r D_{j,r}(1 + sum_s Q^{(j)}_s(Nn - e_r)).
 *
 * Port of matlab/src/api/pfqn/pfqn_mvaoi.m.
 */
public final class Pfqn_mvaoi {
    private Pfqn_mvaoi() {}

    /**
     * Result: per-class throughput X, OI queue-lengths Qoi (K x R), LI queue-lengths
     * Qli (J x R), delay queue-length Qdelay (R).
     */
    public static final class Result {
        public final double[] X;
        public final double[][] Qoi;
        public final double[][] Qli;
        public final double[] Qdelay;
        /**
         * (K x R) per-class mean number of IN-SERVICE jobs at each OI station,
         * i.e. E[sir_r] with sir_r the count of class-r jobs receiving a strictly
         * positive rank rate (see Pfqn_oi_insvc); the utilization of OI station i
         * is Soi[i][r]/c_i. Unlike X/Qoi/Qli, which are pure mean-value
         * quantities, Soi is a distributional statistic and is therefore obtained
         * from the OI count marginal
         *   pM_i(n|k) = (1/mu_i(n)) sum_r X_r(k) pM_i(n-e_r|k-e_r),
         *   pM_i(0|k) = 1 - sum_{n != 0} pM_i(n|k),
         * assembled from the zero-shift throughputs X^{(0)}(k) already cached by
         * the mean-value recursion (no normalizing constant is formed).
         *
         * Computed ON DEMAND, mirroring MATLAB's `if nargout >= 5` guard
         * (pfqn_mvaoi.m:277). This is not a micro-optimization: insvcMeans
         * allocates a (total x total) table with total = prod_r (N_r + 1), so
         * N = [20 20 20] costs a 9261 x 9261 matrix (about 686 MB) and ~2.6e8
         * inner iterations. Paying that on a call that only wanted X reads as a
         * hang rather than an error.
         */
        private double[][] soi;
        private final Ctx ctx;
        private final int[][] zeroS;

        Result(double[] X, double[][] Qoi, double[][] Qli, double[] Qdelay, Ctx ctx, int[][] zeroS) {
            this.X = X; this.Qoi = Qoi; this.Qli = Qli; this.Qdelay = Qdelay;
            this.ctx = ctx; this.zeroS = zeroS;
        }

        /** (K x R) per-class mean number of in-service jobs; see the field doc. */
        public double[][] getSoi() {
            if (soi == null) soi = ctx.insvcMeans(zeroS);
            return soi;
        }
    }

    /** Single-OI, delay-only convenience overload. */
    public static Result pfqn_mvaoi(double[] Z, int[] N, ToDoubleFunction<int[]> mu) {
        List<ToDoubleFunction<int[]>> muList = new ArrayList<ToDoubleFunction<int[]>>();
        muList.add(mu);
        return pfqn_mvaoi(Z, N, muList, new double[0][]);
    }

    /**
     * @param Z   (R) think-time demand vector of the aggregated delay node.
     * @param N   (R) closed population vector, finite.
     * @param mu  list {mu_1,...,mu_K} of OI total-service-rate functions of the
     *            per-class occupancy (count) vector n (same convention as Pfqn_ncoi).
     * @param Dli (J x R) per-class demand matrix of the LI single-server queues;
     *            null or empty when J = 0.
     */
    public static Result pfqn_mvaoi(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu, double[][] Dli) {
        return pfqn_mvaoi(Z, N, mu, Dli, null);
    }

    // General per-OI-station class visit ratios: visits[i] is the 1xR visit
    // vector of OI station i, entering the N_r=1 demand base case
    // theta_{i,r}=v_{i,r}/mu_i(...). ms-promoted stations pass unit visits
    // (their visits are already folded into the rate handle). null -> unit.
    public static Result pfqn_mvaoi(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu, double[][] Dli, double[][] visits) {
        // see _kb/03-api-layer.md for rationale
        if (mu == null || mu.isEmpty())
            throw new IllegalArgumentException(
                    "pfqn_mvaoi: mu must be a nonempty list of OI rate handles");
        return new Ctx(Z, N, mu, Dli, visits).run();
    }

    // ------------------------------------------------------------------
    private static final class Ctx {
        final int R, K, J;
        final double[] Z;
        final int[] N;
        final List<ToDoubleFunction<int[]>> mu;
        final double[][] Dli;
        final double[][] visits;   // per-OI-station class visit ratios (null -> unit)
        final Map<String, double[]> Xc = new HashMap<String, double[]>();   // X^{(S)}(Nn)
        final Map<String, double[][]> Qlc = new HashMap<String, double[][]>(); // Qli^{(S)}(Nn)
        final List<Map<String, double[]>> Dc = new ArrayList<Map<String, double[]>>();
        final List<Map<String, double[]>> Qc = new ArrayList<Map<String, double[]>>();
        final List<Map<String, double[]>> Rc = new ArrayList<Map<String, double[]>>();

        Ctx(double[] Z, int[] N, List<ToDoubleFunction<int[]>> mu, double[][] Dli, double[][] visits) {
            this.R = N.length;
            this.Z = Z;
            this.N = N;
            this.mu = mu;
            this.K = mu.size();
            this.visits = visits;
            this.Dli = (Dli == null) ? new double[0][] : Dli;
            this.J = this.Dli.length;
            for (int i = 0; i < K; i++) {
                Dc.add(new HashMap<String, double[]>());
                Qc.add(new HashMap<String, double[]>());
                Rc.add(new HashMap<String, double[]>());
            }
        }

        String key(int[][] S, int[] Nn) {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < K; i++) {
                for (int r = 0; r < R; r++) sb.append(S[i][r]).append('_');
            }
            for (int r = 0; r < R; r++) sb.append(Nn[r]).append('_');
            return sb.toString();
        }

        /** Return S with the shift row of OI station i incremented in class r. */
        static int[][] shiftPlus(int[][] S, int i, int r) {
            int[][] Sp = new int[S.length][];
            for (int k = 0; k < S.length; k++) Sp[k] = S[k].clone();
            Sp[i][r]++;
            return Sp;
        }

        static int[] minus(int[] a, int r) { int[] o = a.clone(); o[r]--; return o; }

        double rhoF(int i, int r, int[][] S, int[] M) {
            String rkey = key(S, M);
            double[] rrow = Rc.get(i).get(rkey);
            if (rrow != null && !Double.isNaN(rrow[r])) {
                return rrow[r];
            }
            if (rrow == null) {
                rrow = new double[R];
                Arrays.fill(rrow, Double.NaN);
            }
            int tot = 0;
            for (int v : M) tot += v;
            if (tot == 0) {
                rrow[r] = 1.0;
                Rc.get(i).put(rkey, rrow);
                return 1.0;
            }
            int s = -1;
            for (int ss = 0; ss < R; ss++) {
                if (ss != r && M[ss] > 0) { s = ss; break; }
            }
            double xu = Xc.get(key(S, M))[s];
            double xu2 = Xc.get(key(shiftPlus(S, i, r), M))[s];
            double ratio = (xu2 > 0) ? xu / xu2 : 0.0;
            double val = rhoF(i, r, S, minus(M, s)) * ratio;
            rrow[r] = val;
            Rc.get(i).put(rkey, rrow);
            return val;
        }

        Result run() {
            int[][] zeroS = new int[K][R];
            String k0 = key(zeroS, new int[R]);
            Xc.put(k0, new double[R]);
            Qlc.put(k0, new double[J][R]);
            for (int i = 0; i < K; i++) {
                Dc.get(i).put(k0, new double[R]);
                Qc.get(i).put(k0, new double[R]);
            }

            List<int[][]> states = enumStates();
            states.sort(new Comparator<int[][]>() {
                public int compare(int[][] a, int[][] b) {
                    int sa = 0, sb = 0;
                    int[] na = a[K], nb = b[K];
                    for (int v : na) sa += v;
                    for (int v : nb) sb += v;
                    return Integer.compare(sa, sb);
                }
            });

            for (int[][] st : states) {
                int[][] S = new int[K][];
                for (int i = 0; i < K; i++) S[i] = st[i];
                int[] Nn = st[K];
                String kk = key(S, Nn);
                int ntot = 0;
                for (int v : Nn) ntot += v;
                if (ntot == 0) {
                    Xc.put(kk, new double[R]);
                    Qlc.put(kk, new double[J][R]);
                    for (int i = 0; i < K; i++) {
                        Dc.get(i).put(kk, new double[R]);
                        Qc.get(i).put(kk, new double[R]);
                    }
                    continue;
                }

                double[][] Dt = new double[K][R];
                double[][][] Qsub = new double[K][R][R];
                for (int i = 0; i < K; i++) {
                    for (int r = 0; r < R; r++) {
                        if (Nn[r] == 0) continue;
                        if (Nn[r] == 1) {
                            double vfac = (visits == null) ? 1.0 : visits[i][r];
                            double mur = mu.get(i).applyAsDouble(shiftRow(S, i, r));
                            if (mur > 0) {
                                Dt[i][r] = (vfac / mur) * rhoF(i, r, S, minus(Nn, r));
                            }
                        } else {
                            int[] Nmr = minus(Nn, r);
                            double xs = Xc.get(key(S, Nmr))[r];
                            double xs2 = Xc.get(key(shiftPlus(S, i, r), Nmr))[r];
                            if (xs2 > 0) {
                                Dt[i][r] = (xs / xs2) * Dc.get(i).get(key(S, Nmr))[r];
                            }
                        }
                    }
                    for (int s = 0; s < R; s++) {
                        if (Nn[s] > 0) {
                            Qsub[i][s] = Qc.get(i).get(key(shiftPlus(S, i, s), minus(Nn, s)));
                        }
                    }
                }

                double[][] betaLI = new double[J][R];
                for (int r = 0; r < R; r++) {
                    if (Nn[r] == 0) continue;
                    double[][] QsubLI = Qlc.get(key(S, minus(Nn, r)));
                    for (int j = 0; j < J; j++) {
                        double qsum = 0;
                        for (int s = 0; s < R; s++) qsum += QsubLI[j][s];
                        betaLI[j][r] = Dli[j][r] * (1.0 + qsum);
                    }
                }

                List<Integer> idx = new ArrayList<Integer>();
                for (int r = 0; r < R; r++) if (Nn[r] > 0) idx.add(r);
                int m = idx.size();
                double[][] A = new double[m][m];
                for (int a = 0; a < m; a++) {
                    int r = idx.get(a);
                    for (int b = 0; b < m; b++) {
                        int s = idx.get(b);
                        double val = 0;
                        if (s == r) {
                            val = Z[r];
                            for (int j = 0; j < J; j++) val += betaLI[j][r];
                            for (int i = 0; i < K; i++) val += Dt[i][r] * (1.0 + Qsub[i][r][r]);
                        } else {
                            for (int i = 0; i < K; i++) val += Dt[i][s] * Qsub[i][s][r];
                        }
                        A[a][b] = val;
                    }
                }
                // see _kb/03-api-layer.md for rationale
                double[] Xk = new double[R];
                for (int a = 0; a < m; a++) {
                    int r = idx.get(a);
                    double denom = A[a][a];
                    for (int b = 0; b < m; b++) {
                        if (b != a) {
                            int s = idx.get(b);
                            double[] Xner = Xc.get(key(S, minus(Nn, r)));  // X(S, Nn-e_r)
                            double[] Xnes = Xc.get(key(S, minus(Nn, s)));  // X(S, Nn-e_s)
                            if (Xnes[r] > 0) denom += A[a][b] * (Xner[s] / Xnes[r]);
                        }
                    }
                    if (denom > 0) Xk[r] = Nn[r] / denom;
                }

                double[][] QkLI = new double[J][R];
                for (int r = 0; r < R; r++) {
                    if (Nn[r] == 0) continue;
                    for (int j = 0; j < J; j++) QkLI[j][r] = Xk[r] * betaLI[j][r];
                }
                for (int i = 0; i < K; i++) {
                    double[] U = new double[R];
                    for (int r = 0; r < R; r++) U[r] = Dt[i][r] * Xk[r];
                    double[] Qi = new double[R];
                    for (int r = 0; r < R; r++) {
                        double acc = U[r];
                        for (int s = 0; s < R; s++) acc += U[s] * Qsub[i][s][r];
                        Qi[r] = acc;
                    }
                    Qc.get(i).put(kk, Qi);
                    Dc.get(i).put(kk, Dt[i].clone());
                }
                Xc.put(kk, Xk);
                Qlc.put(kk, QkLI);
            }

            String keyN = key(zeroS, N);
            double[] X = Xc.get(keyN);
            double[][] Qoi = new double[K][R];
            for (int i = 0; i < K; i++) Qoi[i] = Qc.get(i).get(keyN);
            double[][] Qli = Qlc.get(keyN);
            double[] Qdelay = new double[R];
            for (int r = 0; r < R; r++) Qdelay[r] = X[r] * Z[r];
            return new Result(X, Qoi, Qli, Qdelay, this, zeroS);
        }

        /**
         * Mean number of in-service jobs per class at each OI station, E[sir_r],
         * from the OI count marginal pM_i(n|k) built on the cached zero-shift
         * throughputs X^{(0)}(k). Exact because in product form
         * pM_i(n|k) = Phi_i(n) G_{-i}(k-n)/G(k) and X_r(k) = G(k-e_r)/G(k), so
         * the recursion reproduces the balanced-fairness identity for Phi_i.
         */
        double[][] insvcMeans(int[][] zeroS) {
            int[] shp = new int[R];
            int[] stride = new int[R];
            int total = 1;
            for (int d = 0; d < R; d++) { shp[d] = N[d] + 1; total *= shp[d]; }
            stride[0] = 1;
            for (int d = 1; d < R; d++) stride[d] = stride[d - 1] * shp[d - 1];

            int[][] subs = new int[total][R];
            Integer[] order = new Integer[total];
            final int[] sums = new int[total];
            for (int i = 0; i < total; i++) {
                int li = i, tot = 0;
                for (int d = 0; d < R; d++) { subs[i][d] = li % shp[d]; li /= shp[d]; tot += subs[i][d]; }
                sums[i] = tot;
                order[i] = i;
            }
            Arrays.sort(order, new Comparator<Integer>() {
                @Override
                public int compare(Integer a, Integer b) { return sums[a] - sums[b]; }
            });

            // Cache X^{(0)}(k) over the lattice.
            double[][] Xk = new double[total][];
            for (int i = 0; i < total; i++) Xk[i] = Xc.get(key(zeroS, subs[i]));

            double[][] Soi = new double[K][R];
            for (int m = 0; m < K; m++) {
                double[][] gm = Pfqn_oi_insvc.pfqn_oi_insvc(mu.get(m), N).g;
                double[] muv = new double[total];
                for (int i = 0; i < total; i++) {
                    if (sums[i] > 0) muv[i] = mu.get(m).applyAsDouble(subs[i]);
                }
                // pMv[a][b] = pM_m(n_a | k_b), filled for n_a <= k_b.
                double[][] pMv = new double[total][total];
                pMv[0][0] = 1.0;                       // pM(0|0) = 1
                for (int bb = 0; bb < total; bb++) {
                    int b = order[bb];
                    if (sums[b] == 0) continue;
                    int[] k = subs[b];
                    double acc0 = 0;
                    for (int aa = 0; aa < total; aa++) {
                        int a = order[aa];
                        if (sums[a] == 0 || muv[a] <= 0) continue;
                        int[] n = subs[a];
                        boolean fits = true;
                        for (int d = 0; d < R; d++) { if (n[d] > k[d]) { fits = false; break; } }
                        if (!fits) continue;
                        double acc = 0;
                        for (int r = 0; r < R; r++) {
                            if (n[r] > 0) acc += Xk[b][r] * pMv[a - stride[r]][b - stride[r]];
                        }
                        pMv[a][b] = acc / muv[a];
                        acc0 += pMv[a][b];
                    }
                    pMv[0][b] = 1 - acc0;              // empty state by complement
                }
                int idxN = 0;
                for (int d = 0; d < R; d++) idxN += N[d] * stride[d];
                for (int r = 0; r < R; r++) {
                    double acc = 0;
                    for (int i = 0; i < total; i++) acc += pMv[i][idxN] * gm[i][r];
                    Soi[m][r] = acc;
                }
            }
            return Soi;
        }

        /** Per-class occupancy row s_i + e_r of OI station i. */
        int[] shiftRow(int[][] S, int i, int r) {
            int[] row = S[i].clone();
            row[r]++;
            return row;
        }

        /**
         * All reachable (S, Nn): per class r the K+1 buckets (K OI shifts, then the
         * free bucket Nn) sum to <= N_r. Obtained as the first K+1 entries of each
         * (K+2)-part composition of N_r, whose dropped (K+2)-th part is the slack
         * bucket; this mirrors MATLAB sprod(K+2, N) with the slack row removed.
         */
        List<int[][]> enumStates() {
            List<List<int[]>> perclass = new ArrayList<List<int[]>>();
            for (int r = 0; r < R; r++) {
                List<int[]> full = new ArrayList<int[]>();
                compositions(N[r], K + 2, full);   // (K+2)-part compositions of N_r
                List<int[]> rows = new ArrayList<int[]>();
                for (int[] comp : full) {
                    int[] row = new int[K + 1];     // drop the slack (last) bucket
                    System.arraycopy(comp, 0, row, 0, K + 1);
                    rows.add(row);
                }
                perclass.add(rows);
            }
            int[] counts = new int[R];
            long total = 1;
            for (int r = 0; r < R; r++) { counts[r] = perclass.get(r).size(); total *= counts[r]; }
            List<int[][]> out = new ArrayList<int[][]>();
            int[] odo = new int[R];
            for (long p = 0; p < total; p++) {
                int[][] st = new int[K + 1][];
                for (int i = 0; i < K; i++) st[i] = new int[R];
                st[K] = new int[R];
                for (int r = 0; r < R; r++) {
                    int[] row = perclass.get(r).get(odo[r]);
                    for (int i = 0; i < K; i++) st[i][r] = row[i];
                    st[K][r] = row[K];
                }
                out.add(st);
                for (int r = 0; r < R; r++) {
                    if (++odo[r] < counts[r]) break;
                    odo[r] = 0;
                }
            }
            return out;
        }

        void compositions(int mm, int p, List<int[]> out) {
            compRec(mm, p, new int[p], 0, out);
        }

        void compRec(int rem, int p, int[] cur, int pos, List<int[]> out) {
            if (pos == p - 1) {
                cur[pos] = rem;
                out.add(cur.clone());
                return;
            }
            for (int v = 0; v <= rem; v++) {
                cur[pos] = v;
                compRec(rem - v, p, cur, pos + 1, out);
            }
        }
    }

}
