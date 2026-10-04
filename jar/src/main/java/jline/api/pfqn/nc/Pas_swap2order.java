/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.pfqn.nc;

import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Derive the GLOBAL placement-order DAG H of a closed two-station pass-and-swap
 * (P&S) tandem 1-&gt;2-&gt;1 directly from its swap graph, for use with
 * {@link Pfqn_pas_is}.
 *
 * <p>With a non-empty swap graph the ordered-state chain is reducible (Comte and
 * Dorsman, 2021): the recurrent communicating class is the set of splits of the
 * orderings that are the linear extensions of a single placement partial order.
 * {@link Pfqn_pas_is} samples those orderings from H, so it needs exactly this
 * global order. The order is a class-level property (multiplicity independent),
 * extracted from the single-job-per-class instance: enumerate the reachable
 * communicating class from the all-in-queue-1 initial state; each reachable
 * state (l1;l2) exposes the full ordering c = [l1, reverse(l2)] in D; then
 * H[i][j] = 1 iff class i precedes class j in EVERY c in D (forced precedence).
 *
 * <p>Port of matlab/src/api/pfqn/pas_swap2order.m.
 */
public final class Pas_swap2order {
    private Pas_swap2order() {}

    /**
     * @param G1       swap graph of queue 1 (R x R, 0-based; G(a,b) != 0 means a
     *                 chases b). Null treated as all-zero.
     * @param G2       swap graph of queue 2.
     * @param svc1     ordered-list service rate function of queue 1.
     * @param svc2     ordered-list service rate function of queue 2.
     * @param N0       minimal probing population (typically ones(1,R)).
     * @return (R x R) global placement-order DAG H; H[i][j]=1 iff i precedes j.
     */
    public static int[][] pas_swap2order(Matrix G1, Matrix G2,
                                         SerializableFunction<Matrix, Double> svc1,
                                         SerializableFunction<Matrix, Double> svc2,
                                         int[] N0) {
        int R = N0.length;
        int[][] H = new int[R][R];

        boolean empty1 = (G1 == null) || !hasNonZero(G1);
        boolean empty2 = (G2 == null) || !hasNonZero(G2);
        if (empty1 && empty2) {
            return H;   // pure OI: no placement constraint
        }

        Matrix[] swap = {G1, G2};
        @SuppressWarnings("unchecked")
        SerializableFunction<Matrix, Double>[] listRate = new SerializableFunction[]{svc1, svc2};

        // minimal single-job-per-class initial placement, all jobs at queue 1
        List<Integer> initList = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            for (int t = 0; t < N0[r]; t++) {
                initList.add(r);
            }
        }
        int[] initState = toArray(initList);
        int[][] init = {initState, new int[0]};

        // breadth-first enumeration of the reachable (communicating) class
        Set<String> seen = new HashSet<String>();
        List<int[][]> classStates = new ArrayList<int[][]>();
        Deque<int[][]> frontier = new ArrayDeque<int[][]>();
        seen.add(encode(init));
        frontier.push(init);
        while (!frontier.isEmpty()) {
            int[][] st = frontier.pop();
            classStates.add(st);
            for (int m = 0; m < 2; m++) {
                int[] c = st[m];
                int next = (m + 1) % 2;
                for (int pp = 0; pp < c.length; pp++) {
                    // marginal OI completion rate of the pp-th customer; the
                    // empty-prefix rate is 0 by definition (guard pp==0).
                    double prevRate = (pp == 0) ? 0.0 : applyPrefix(listRate[m], c, pp);
                    double rate = applyPrefix(listRate[m], c, pp + 1) - prevRate;
                    if (rate <= 1e-12) {
                        continue;
                    }
                    SwapResult sr = swapLocal(c, pp, swap[m]);
                    int[][] stn = new int[2][];
                    stn[m] = sr.cnew;
                    stn[next] = append(st[next], sr.dep);
                    String k = encode(stn);
                    if (!seen.contains(k)) {
                        seen.add(k);
                        frontier.push(stn);
                    }
                }
            }
        }

        // full orderings D: every reachable state exposes c = [l1, reverse(l2)]
        Set<String> dseen = new HashSet<String>();
        List<int[]> D = new ArrayList<int[]>();
        for (int s = 0; s < classStates.size(); s++) {
            int[] l1 = classStates.get(s)[0];
            int[] l2 = classStates.get(s)[1];
            int[] c = new int[l1.length + l2.length];
            System.arraycopy(l1, 0, c, 0, l1.length);
            for (int i = 0; i < l2.length; i++) {
                c[l1.length + i] = l2[l2.length - 1 - i];   // reversed suffix
            }
            String kc = encodeArr(c);
            if (!dseen.contains(kc)) {
                dseen.add(kc);
                D.add(c);
            }
        }

        // forced precedence: i precedes j iff in every c in D where both present,
        // every copy of i is older than every copy of j (position-wise).
        for (int a = 0; a < R; a++) {
            for (int b = 0; b < R; b++) {
                if (a == b) {
                    continue;
                }
                boolean both = false;
                boolean forced = true;
                for (int s = 0; s < D.size(); s++) {
                    int[] c = D.get(s);
                    int maxA = -1;
                    int minB = Integer.MAX_VALUE;
                    boolean hasA = false, hasB = false;
                    for (int q = 0; q < c.length; q++) {
                        if (c[q] == a) {
                            hasA = true;
                            if (q > maxA) maxA = q;
                        } else if (c[q] == b) {
                            hasB = true;
                            if (q < minB) minB = q;
                        }
                    }
                    if (!hasA || !hasB) {
                        continue;
                    }
                    both = true;
                    if (!(maxA < minB)) {
                        forced = false;
                        break;
                    }
                }
                H[a][b] = (both && forced) ? 1 : 0;
            }
        }
        return H;
    }

    // ======================================================================
    private static final class SwapResult {
        final int[] cnew;
        final int dep;

        SwapResult(int[] cnew, int dep) {
            this.cnew = cnew;
            this.dep = dep;
        }
    }

    /**
     * Pass-and-swap scan: the customer at position p (0-based) completes; it
     * chases along the swap graph, classes shift one step along the chain, and
     * the last ejected customer departs; the served slot is removed.
     */
    private static SwapResult swapLocal(int[] c, int p, Matrix G) {
        int n = c.length;
        List<Integer> chain = new ArrayList<Integer>();
        chain.add(p);
        int moving = c[p];
        int cur = p;
        while (true) {
            int q = -1;
            for (int j = cur + 1; j < n; j++) {
                if (G != null && G.get(moving, c[j]) != 0) {
                    q = j;
                    break;
                }
            }
            if (q == -1) {
                break;
            }
            chain.add(q);
            moving = c[q];
            cur = q;
        }
        int dep = c[chain.get(chain.size() - 1)];
        int[] tmp = c.clone();
        for (int i = 0; i < chain.size() - 1; i++) {
            tmp[chain.get(i + 1)] = c[chain.get(i)];
        }
        // remove position chain.get(0)
        int removeAt = chain.get(0);
        int[] cnew = new int[n - 1];
        int idx = 0;
        for (int i = 0; i < n; i++) {
            if (i == removeAt) {
                continue;
            }
            cnew[idx++] = tmp[i];
        }
        return new SwapResult(cnew, dep);
    }

    /** listRate on the length-len prefix of c, as a (1 x len) 0-based Matrix. */
    private static double applyPrefix(SerializableFunction<Matrix, Double> fun, int[] c, int len) {
        Matrix pref = new Matrix(1, len);
        for (int i = 0; i < len; i++) {
            pref.set(0, i, c[i]);
        }
        return fun.apply(pref);
    }

    private static int[] append(int[] arr, int v) {
        int[] out = new int[arr.length + 1];
        System.arraycopy(arr, 0, out, 0, arr.length);
        out[arr.length] = v;
        return out;
    }

    private static int[] toArray(List<Integer> l) {
        int[] out = new int[l.size()];
        for (int i = 0; i < l.size(); i++) {
            out[i] = l.get(i);
        }
        return out;
    }

    private static String encode(int[][] st) {
        StringBuilder sb = new StringBuilder();
        for (int m = 0; m < st.length; m++) {
            for (int v : st[m]) {
                sb.append(v).append(',');
            }
            sb.append('|');
        }
        return sb.toString();
    }

    private static String encodeArr(int[] c) {
        StringBuilder sb = new StringBuilder();
        for (int v : c) {
            sb.append(v).append(',');
        }
        return sb.toString();
    }

    private static boolean hasNonZero(Matrix m) {
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                if (m.get(i, j) != 0) {
                    return true;
                }
            }
        }
        return false;
    }
}
