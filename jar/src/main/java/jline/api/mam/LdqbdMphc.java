/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.matrix.Matrix;

/**
 * Exact level-dependent QBD blocks of an M/PH/c queue.
 *
 * <p>Port of matlab/src/api/mam/ldqbd_mphc.m and ph_multisets.m. The level of
 * the chain is the number of jobs at the station; the coordinate INSIDE a level
 * is the MULTISET of the phases the min(n,c) busy servers sit in, which is what
 * makes the construction exact for phase-type service at any number of servers.
 *
 * <p>The collapsed alternative -- one PH process run at min(n,c) times its
 * speed -- gets the aggregate service rate right but forgets which phase each
 * busy server is in, turning the c servers into one fast server whose remaining
 * work is a single phase-type variable.
 *
 * <p>References: S. Asmussen and J.R. Moller, "Calculation of the steady state
 * waiting time distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems
 * 37(1):9-29, 2001; M. F. Neuts, "Matrix-geometric solutions in stochastic
 * models", Johns Hopkins University Press, 1981.
 */
public final class LdqbdMphc {
    private LdqbdMphc() {}

    /**
     * The repeating level is the widest one, so it is the size worth guarding:
     * the LD-QBD recursion inverts one matrix of that order per level.
     */
    public static final int MAX_CONFIGS = 2000;

    /** The three block lists of a level-dependent QBD, in the order {@link Ldqbd} takes them. */
    public static final class Blocks {
        public final List<Matrix> Q0;   // size Nlev   : upward, level n -> n+1
        public final List<Matrix> Q1;   // size Nlev+1 : local
        public final List<Matrix> Q2;   // size Nlev   : downward, Q2.get(n-1) leaves level n

        Blocks(List<Matrix> Q0, List<Matrix> Q1, List<Matrix> Q2) {
            this.Q0 = Q0; this.Q1 = Q1; this.Q2 = Q2;
        }
    }

    /**
     * Configurations of k identical servers over p service phases.
     *
     * <p>Rows are the compositions of k into p nonnegative parts: entry (r,i) is
     * the number of the k busy servers sitting in phase i. There are
     * nchoosek(k+p-1, p-1) of them, the multiset count of Asmussen and Moller
     * (2001) -- identical servers are exchangeable, so only the phase COUNTS
     * carry information and the ordered space of size p^k collapses onto this.
     *
     * <p>The order is fixed and shared by every caller, so a configuration index
     * means the same thing in each of them: the first part descends. k == 1
     * therefore yields the identity rows e_1 ... e_p in phase order, which is
     * what makes the c == 1 case of {@link #ldqbd_mphc} coincide with plain
     * phase indexing.
     */
    public static int[][] ph_multisets(int p, int k) {
        if (k == 0) {
            return new int[][]{new int[p]};
        }
        if (p == 1) {
            return new int[][]{new int[]{k}};
        }
        List<int[]> rows = new ArrayList<int[]>();
        for (int first = k; first >= 0; first--) {
            int[][] rest = ph_multisets(p - 1, k - first);
            for (int r = 0; r < rest.length; r++) {
                int[] row = new int[p];
                row[0] = first;
                System.arraycopy(rest[r], 0, row, 1, p - 1);
                rows.add(row);
            }
        }
        return rows.toArray(new int[rows.size()][]);
    }

    /**
     * Block-tridiagonal generator of an M/PH/c queue with level-dependent arrivals.
     *
     * @param D0      service sub-generator (p x p), phase changes without completion
     * @param D1      service completion block (p x p); D1 = (-D0*1)*alpha for a PH
     * @param alpha   1 x p vector a server starts each new job in
     * @param c       number of identical servers (&gt;= 1; capped at the top level)
     * @param arrRate length Nlev+1; arrRate[n] is the arrival rate out of level n
     * @param sf      optional length-Nlev multiplier on the station's TOTAL service
     *                rate at level n (load dependence); each busy server then runs at
     *                sf[n-1]/min(n,c) of nominal, so sf[n-1] == min(n,c) reproduces
     *                the unscaled queue exactly. Pass null for no scaling.
     *
     * <p>Level sizes grow over the boundary levels 0..c and repeat above them, so the
     * blocks joining differently sized neighbours are rectangular; Ldqbd, its rate
     * matrices and its stationary vector all accept that heterogeneity.
     */
    public static Blocks ldqbd_mphc(Matrix D0, Matrix D1, Matrix alpha, double c,
                                    double[] arrRate, double[] sf) {
        final int p = D0.getNumRows();
        final int Nlev = arrRate.length - 1;

        if (Nlev < 1) {
            throw new RuntimeException("ldqbd_mphc needs at least one level above the empty one.");
        }
        if (D0.getNumCols() != p || D1.getNumRows() != p || D1.getNumCols() != p
                || alpha.length() != p) {
            throw new RuntimeException("D0, D1 and alpha must all have the same order.");
        }
        if (sf != null && sf.length < Nlev) {
            throw new RuntimeException(
                    "sf must give one total-service-rate factor per level 1..Nlev.");
        }

        // Servers that can never be busy do not need a coordinate: above the top
        // level there is nothing left to serve.
        final int cmax = !Double.isFinite(c) ? Nlev
                : Math.min(Math.max(1, (int) Math.round(c)), Nlev);

        List<int[][]> cfg = new ArrayList<int[][]>();
        List<Map<String, Integer>> pos = new ArrayList<Map<String, Integer>>();
        int[] nCfg = new int[cmax + 1];
        for (int k = 0; k <= cmax; k++) {
            int[][] Ck = ph_multisets(p, k);
            cfg.add(Ck);
            Map<String, Integer> ix = new HashMap<String, Integer>();
            for (int r = 0; r < Ck.length; r++) {
                ix.put(cfgKey(Ck[r]), Integer.valueOf(r));
            }
            pos.add(ix);
            nCfg[k] = Ck.length;
        }

        if (nCfg[cmax] > MAX_CONFIGS) {
            throw new RuntimeException(String.format(
                    "the exact M/PH/c chain needs nchoosek(%d+%d-1,%d-1) = %d configurations per "
                    + "level for %d servers and %d service phases, above the %d the level-by-level "
                    + "inverses can carry. Use fewer phases (a lower-order fit), fewer servers, or "
                    + "SolverCTMC/SolverLDES on this model.",
                    cmax, p, p, nCfg[cmax], cmax, p, MAX_CONFIGS));
        }

        // Completion rate out of each phase, summed over targets.
        double[] t = new double[p];
        for (int i = 0; i < p; i++) {
            double s = 0.0;
            for (int j = 0; j < p; j++) s += D1.get(i, j);
            t[i] = s;
        }

        // Structural blocks per busy-server count: LOC the within-level phase
        // changes (with the full outflow on its diagonal), UP the entry of a newly
        // busy server, DN a completion that leaves a server idle.
        List<Matrix> LOC = new ArrayList<Matrix>();
        List<Matrix> UP = new ArrayList<Matrix>();
        List<Matrix> DN = new ArrayList<Matrix>();
        for (int k = 0; k <= cmax; k++) {
            int[][] Ck = cfg.get(k);
            int nk = nCfg[k];
            Matrix Lk = new Matrix(nk, nk);
            Lk.zero();
            for (int row = 0; row < nk; row++) {
                int[] m = Ck[row];
                for (int i = 0; i < p; i++) {
                    if (m[i] == 0) continue;
                    for (int j = 0; j < p; j++) {
                        if (j == i) continue;
                        int col = pos.get(k).get(cfgKey(move(m, i, j))).intValue();
                        Lk.set(row, col, Lk.get(row, col) + m[i] * D0.get(i, j));
                    }
                    // D0(i,i) is the total outflow of phase i, completions included
                    Lk.set(row, row, Lk.get(row, row) + m[i] * D0.get(i, i));
                }
            }
            LOC.add(Lk);

            if (k < cmax) {
                Matrix Uk = new Matrix(nk, nCfg[k + 1]);
                Uk.zero();
                for (int row = 0; row < nk; row++) {
                    int[] m = Ck[row];
                    for (int j = 0; j < p; j++) {
                        int col = pos.get(k + 1).get(cfgKey(plus(m, j))).intValue();
                        Uk.set(row, col, Uk.get(row, col) + alpha.get(j));
                    }
                }
                UP.add(Uk);
            } else {
                UP.add(null);
            }

            if (k > 0) {
                Matrix Dk = new Matrix(nk, nCfg[k - 1]);
                Dk.zero();
                for (int row = 0; row < nk; row++) {
                    int[] m = Ck[row];
                    for (int i = 0; i < p; i++) {
                        if (m[i] == 0) continue;
                        int col = pos.get(k - 1).get(cfgKey(minus(m, i))).intValue();
                        Dk.set(row, col, Dk.get(row, col) + m[i] * t[i]);
                    }
                }
                DN.add(Dk);
            } else {
                DN.add(null);
            }
        }

        // A completion at a full server bank takes the next waiting job at once, so
        // the server stays busy and only its phase moves: the repeating down block.
        int[][] Cc = cfg.get(cmax);
        int nc = nCfg[cmax];
        Matrix CDEP = new Matrix(nc, nc);
        CDEP.zero();
        for (int row = 0; row < nc; row++) {
            int[] m = Cc[row];
            for (int i = 0; i < p; i++) {
                if (m[i] == 0) continue;
                for (int j = 0; j < p; j++) {
                    int col = pos.get(cmax).get(cfgKey(move(m, i, j))).intValue();
                    CDEP.set(row, col, CDEP.get(row, col) + m[i] * D1.get(i, j));
                }
            }
        }

        // Per-server speed. Without load dependence every busy server runs at its
        // nominal rate; with it, the aggregate sf(n) is shared over the busy
        // servers. sf(n) == min(n,c) is passed through as exactly 1 so the unscaled
        // chain is reproduced bit for bit.
        double[] speed = new double[Nlev];
        for (int n = 1; n <= Nlev; n++) {
            speed[n - 1] = 1.0;
            if (sf != null) {
                int b = Math.min(n, cmax);
                if (sf[n - 1] != b) speed[n - 1] = sf[n - 1] / b;
            }
        }

        List<Matrix> Q0 = new ArrayList<Matrix>();
        List<Matrix> Q1 = new ArrayList<Matrix>();
        List<Matrix> Q2 = new ArrayList<Matrix>();

        Matrix Q1_0 = new Matrix(1, 1);          // level 0: arrivals only
        Q1_0.set(0, 0, -arrRate[0]);
        Q1.add(Q1_0);
        for (int n = 1; n <= Nlev; n++) {
            int b = Math.min(n, cmax);
            Q1.add(LOC.get(b).scale(speed[n - 1]).sub(Matrix.eye(nCfg[b]).scale(arrRate[n])));
        }

        for (int n = 0; n < Nlev; n++) {
            if (n < cmax) {
                Q0.add(UP.get(n).scale(arrRate[n]));            // a free server takes the job
            } else {
                Q0.add(Matrix.eye(nCfg[cmax]).scale(arrRate[n]));  // it waits, phases unchanged
            }
        }

        for (int n = 1; n <= Nlev; n++) {
            if (n <= cmax) {
                Q2.add(DN.get(n).scale(speed[n - 1]));          // the server falls idle
            } else {
                Q2.add(CDEP.scale(speed[n - 1]));               // it takes the next job
            }
        }

        return new Blocks(Q0, Q1, Q2);
    }

    private static int[] move(int[] m, int i, int j) {
        int[] mm = m.clone();
        mm[i]--;
        mm[j]++;
        return mm;
    }

    private static int[] plus(int[] m, int j) {
        int[] mm = m.clone();
        mm[j]++;
        return mm;
    }

    private static int[] minus(int[] m, int i) {
        int[] mm = m.clone();
        mm[i]--;
        return mm;
    }

    /** Key of a configuration row, for the index maps above. */
    private static String cfgKey(int[] m) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < m.length; i++) {
            sb.append(m[i]).append(',');
        }
        return sb.toString();
    }
}
