/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Groups closed chains by similarity of their per-station SERVICE DEMAND, for
 * the aggregation {@link Solver_mam_bgchain} applies to the chains it is not
 * carrying exactly.
 *
 * <p>WHY DEMAND IS THE RIGHT CRITERION. The aggregate that replaces a group
 * carries the flow-weighted mean of its members' service times and routing, so
 * the group aggregates EXACTLY when its members place the same demand at every
 * station, and distorts both quantities in proportion to how far apart they are.
 * The distance is therefore the symmetric relative L1 gap between the demand
 * vectors,</p>
 *
 * <pre>dist(a,b) = sum_i |D[i][a] - D[i][b]| / ((sum_i D[i][a] + sum_i D[i][b])/2),</pre>
 *
 * <p>which is scale-relative rather than absolute: it separates two chains whose
 * demand PROFILE across the stations differs and two chains whose profile agrees
 * but whose magnitude does not, and being dimensionless it groups a model the
 * same way whatever its time unit.</p>
 *
 * <p>WHY COMPLETE LINKAGE. The clustering is agglomerative from singletons,
 * merging at each step the pair of clusters whose WORST member-to-member
 * distance is smallest. The aggregation error inside a group is driven by its
 * worst mismatch and not by its average one, so complete linkage is the
 * criterion that bounds what the aggregation actually costs; average or single
 * linkage would let one distant chain ride along inside an otherwise tight
 * group.</p>
 *
 * <p>DETERMINISM. Ties are broken by the lexicographically smallest pair of
 * cluster indices and the groups are relabelled by their smallest member, so the
 * same input gives the same grouping in MATLAB, the JAR, Python and C++.</p>
 *
 * @see Solver_mam_bgchain
 */
public final class Mam_bgchain_groups {
    private Mam_bgchain_groups() {}

    private static final double ZERO = 1e-14;

    /**
     * @param D (Mc x n) per-station demand, one column per chain
     * @param G number of groups wanted, clamped to [1, n]
     * @return group index in 0..G-1 of each chain
     */
    public static int[] mam_bgchain_groups(double[][] D, int G) {
        final int Mc = D.length;
        final int n = (Mc == 0) ? 0 : D[0].length;
        int[] grp = new int[n];
        if (n == 0) return grp;
        if (G < 1) G = 1;
        if (G > n) G = n;

        double[] tot = new double[n];
        for (int c = 0; c < n; c++) {
            for (int i = 0; i < Mc; i++) tot[c] += D[i][c];
        }
        double[][] dist = new double[n][n];
        for (int a = 0; a < n; a++) {
            for (int b = a + 1; b < n; b++) {
                double den = (tot[a] + tot[b]) / 2.0;
                double d = 0.0;
                if (den > ZERO) {
                    double acc = 0.0;
                    for (int i = 0; i < Mc; i++) acc += Math.abs(D[i][a] - D[i][b]);
                    d = acc / den;
                }
                dist[a][b] = d;
                dist[b][a] = d;
            }
        }

        List<List<Integer>> clusters = new ArrayList<List<Integer>>();
        for (int a = 0; a < n; a++) {
            List<Integer> one = new ArrayList<Integer>();
            one.add(Integer.valueOf(a));
            clusters.add(one);
        }
        boolean[] active = new boolean[n];
        java.util.Arrays.fill(active, true);
        int nactive = n;
        while (nactive > G) {
            double best = Double.POSITIVE_INFINITY;
            int bp = -1, bq = -1;
            for (int p = 0; p < n; p++) {
                if (!active[p]) continue;
                for (int q = p + 1; q < n; q++) {
                    if (!active[q]) continue;
                    double d = 0.0;
                    for (Integer x : clusters.get(p)) {
                        for (Integer y : clusters.get(q)) {
                            double v = dist[x.intValue()][y.intValue()];
                            if (v > d) d = v;
                        }
                    }
                    if (d < best - ZERO) {
                        best = d;
                        bp = p;
                        bq = q;
                    }
                }
            }
            if (bp < 0) break;
            clusters.get(bp).addAll(clusters.get(bq));
            Collections.sort(clusters.get(bp));
            clusters.get(bq).clear();
            active[bq] = false;
            nactive--;
        }

        // Relabel by smallest member, so the group numbering is canonical
        List<Integer> firsts = new ArrayList<Integer>();
        for (int p = 0; p < n; p++) {
            if (active[p]) firsts.add(Integer.valueOf(p));
        }
        Collections.sort(firsts, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Integer.compare(clusters.get(a.intValue()).get(0).intValue(),
                        clusters.get(b.intValue()).get(0).intValue());
            }
        });
        for (int g = 0; g < firsts.size(); g++) {
            for (Integer x : clusters.get(firsts.get(g).intValue())) {
                grp[x.intValue()] = g;
            }
        }
        return grp;
    }
}
