/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

/**
 * Exact marginal load-dependent MVA for a closed product-form network of
 * infinite-server (delay) and load-independent (single-server, product-form)
 * stations plus ANY number of order-independent (OI) stations. Carries, for each
 * OI station, its joint count-vector marginal distribution pM_i(n | k):
 *
 *   pM_i(n | k) = (1/mu_i(n)) sum_r X_r(k) pM_i(n - e_r | k - e_r),  n != 0
 *   pM_i(0 | k) = 1 - sum_{n != 0} pM_i(n | k)
 *
 * The per-class throughput X_r(k) is closed at each population level by
 *   X_r(k) A_r(k) + sum_i QM_ir(k;X) = k_r,  A_r(k) = sum_{i not OI} R_ir(k).
 * This is the marginal-distribution counterpart of {@link Pfqn_mvaoi} (mean-value
 * CMVA form). Port of matlab/src/api/pfqn/pfqn_mvaoi_marg.m.
 */
public final class Pfqn_mvaoi_marg {
    private Pfqn_mvaoi_marg() {}

    /** Result: per-class throughput XN (R) and per-station queue-lengths QN (M x R). */
    public static final class Result {
        public final double[] XN;
        public final double[][] QN;
        public Result(double[] XN, double[][] QN) { this.XN = XN; this.QN = QN; }
    }

    /**
     * @param D       (M x R) per-class demand at every station (OI rows ignored).
     * @param N       (R) closed population vector, finite.
     * @param isDelay (M) true for infinite-server (delay) stations.
     * @param mu      (M) OI rate functions of the count vector n; null for non-OI.
     */
    public static Result pfqn_mvaoi_marg(double[][] D, int[] N, boolean[] isDelay,
                                         List<ToDoubleFunction<int[]>> mu) {
        int M = D.length;
        int R = N.length;
        List<Integer> oiList = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) if (mu.get(i) != null) oiList.add(i);
        int nOI = oiList.size();
        if (nOI == 0) {
            throw new RuntimeException("pfqn_mvaoi_marg requires at least one OI station.");
        }

        Map<String, double[]> Xc = new HashMap<String, double[]>();
        Map<String, double[][]> Qc = new HashMap<String, double[][]>();
        List<Map<String, Map<String, Double>>> pM = new ArrayList<Map<String, Map<String, Double>>>();
        for (int o = 0; o < nOI; o++) pM.add(new HashMap<String, Map<String, Double>>());

        int[] zero = new int[R];
        String zeroKey = key(zero);
        Xc.put(zeroKey, new double[R]);
        Qc.put(zeroKey, new double[M][R]);
        for (int o = 0; o < nOI; o++) {
            Map<String, Double> pm0 = new HashMap<String, Double>();
            pm0.put(zeroKey, 1.0);
            pM.get(o).put(zeroKey, pm0);
        }

        List<int[]> pops = enumVecs(N);
        Collections.sort(pops, new Comparator<int[]>() {
            public int compare(int[] a, int[] b) {
                int sa = 0, sb = 0;
                for (int v : a) sa += v;
                for (int v : b) sb += v;
                if (sa != sb) return Integer.compare(sa, sb);
                for (int i = 0; i < a.length; i++) if (a[i] != b[i]) return Integer.compare(a[i], b[i]);
                return 0;
            }
        });

        for (int[] k : pops) {
            int ktot = 0;
            for (int v : k) ktot += v;
            if (ktot == 0) continue;
            String kkey = key(k);

            double[][] Rfix = new double[M][R];
            double[] A = new double[R];
            for (int r = 0; r < R; r++) {
                if (k[r] == 0) continue;
                int[] kr = k.clone(); kr[r]--;
                double[][] Qkr = Qc.get(key(kr));
                for (int i = 0; i < M; i++) {
                    if (mu.get(i) != null) continue;
                    if (isDelay[i]) {
                        Rfix[i][r] = D[i][r];
                    } else {
                        double qsum = 0;
                        for (int s = 0; s < R; s++) qsum += Qkr[i][s];
                        Rfix[i][r] = D[i][r] * (1 + qsum);
                    }
                    A[r] += Rfix[i][r];
                }
            }

            double[] Xk = new double[R];
            for (int r = 0; r < R; r++) if (k[r] > 0) Xk[r] = k[r] / (A[r] + 1.0);
            List<Map<String, Double>> pmk = new ArrayList<Map<String, Double>>();
            for (int o = 0; o < nOI; o++) pmk.add(null);
            for (int it = 0; it < 2000; it++) {
                double[] QMtot = new double[R];
                for (int o = 0; o < nOI; o++) {
                    pmk.set(o, oiMarginal(k, Xk, R, mu.get(oiList.get(o)), pM.get(o)));
                    double[] qm = marginalMeans(pmk.get(o), R);
                    for (int r = 0; r < R; r++) QMtot[r] += qm[r];
                }
                double[] Xnew = new double[R];
                double diff = 0;
                for (int r = 0; r < R; r++) {
                    if (k[r] > 0 && A[r] > 0) Xnew[r] = Math.max((k[r] - QMtot[r]) / A[r], 0.0);
                    diff = Math.max(diff, Math.abs(Xnew[r] - Xk[r]));
                }
                if (diff < 1e-13) { Xk = Xnew; break; }
                for (int r = 0; r < R; r++) Xk[r] = 0.5 * Xk[r] + 0.5 * Xnew[r];
            }

            for (int o = 0; o < nOI; o++) {
                pmk.set(o, oiMarginal(k, Xk, R, mu.get(oiList.get(o)), pM.get(o)));
            }
            double[][] Qk = new double[M][R];
            for (int o = 0; o < nOI; o++) {
                double[] qm = marginalMeans(pmk.get(o), R);
                Qk[oiList.get(o)] = qm;
            }
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    if (mu.get(i) != null) continue;
                    Qk[i][r] = Xk[r] * Rfix[i][r];
                }
            }
            Xc.put(kkey, Xk);
            Qc.put(kkey, Qk);
            for (int o = 0; o < nOI; o++) pM.get(o).put(kkey, pmk.get(o));
        }

        String Nkey = key(N);
        return new Result(Xc.get(Nkey), Qc.get(Nkey));
    }

    private static Map<String, Double> oiMarginal(int[] k, double[] Xk, int R,
            ToDoubleFunction<int[]> mu, Map<String, Map<String, Double>> pM) {
        Map<String, Double> pmk = new HashMap<String, Double>();
        List<int[]> vecs = enumVecs(k);
        double psum = 0;
        for (int[] n : vecs) {
            int tot = 0;
            for (int v : n) tot += v;
            if (tot == 0) continue;
            double rate = mu.applyAsDouble(n);
            if (rate <= 0) continue;
            double acc = 0;
            for (int r = 0; r < R; r++) {
                if (n[r] >= 1 && k[r] >= 1) {
                    int[] nr = n.clone(); nr[r]--;
                    int[] kr = k.clone(); kr[r]--;
                    Map<String, Double> prev = pM.get(key(kr));
                    if (prev != null) {
                        Double pv = prev.get(key(nr));
                        if (pv != null) acc += Xk[r] * pv;
                    }
                }
            }
            double p = acc / rate;
            pmk.put(key(n), p);
            psum += p;
        }
        pmk.put(key(new int[R]), 1.0 - psum);
        return pmk;
    }

    private static double[] marginalMeans(Map<String, Double> pmk, int R) {
        double[] QM = new double[R];
        for (Map.Entry<String, Double> e : pmk.entrySet()) {
            int[] n = unkey(e.getKey(), R);
            double p = e.getValue();
            for (int r = 0; r < R; r++) QM[r] += n[r] * p;
        }
        return QM;
    }

    private static List<int[]> enumVecs(int[] bound) {
        List<int[]> out = new ArrayList<int[]>();
        int[] cur = new int[bound.length];
        enumRec(bound, 0, cur, out);
        return out;
    }

    private static void enumRec(int[] bound, int idx, int[] cur, List<int[]> out) {
        if (idx == bound.length) { out.add(cur.clone()); return; }
        for (int v = 0; v <= bound[idx]; v++) {
            cur[idx] = v;
            enumRec(bound, idx + 1, cur, out);
        }
        cur[idx] = 0;
    }

    private static String key(int[] v) {
        StringBuilder sb = new StringBuilder();
        for (int x : v) sb.append(x).append('_');
        return sb.toString();
    }

    private static int[] unkey(String s, int R) {
        String[] parts = s.split("_");
        int[] v = new int[R];
        for (int r = 0; r < R && r < parts.length; r++) v[r] = Integer.parseInt(parts[r]);
        return v;
    }
}
