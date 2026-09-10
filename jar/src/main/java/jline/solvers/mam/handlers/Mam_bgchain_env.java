/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;

/**
 * Lumps the background modulating chain onto the number of closed jobs held by
 * one station, giving the Markovian environment that the open classes at that
 * station see.
 *
 * <p>Station i does not observe the whole closed population vector, only how
 * many closed jobs compete with the open ones for its server. The lumped
 * generator is the stationary-weighted aggregation of the chain's generator over
 * the level sets {s : totocc(s,i) = e},</p>
 *
 * <pre>A(e,e') = sum_{s in e} pi(s) sum_{s' in e'} Q(s,s') / sum_{s in e} pi(s),</pre>
 *
 * <p>which is exact when the partition is lumpable in the Kemeny-Snell sense and
 * is the standard exact-aggregation approximation otherwise. The diagonal is set
 * from the off-diagonal row sums, so A is a proper generator whatever the
 * lumping error is. Environment states of zero stationary probability are
 * unreachable and are dropped, so the support need not be 0..N.</p>
 *
 * @see Mam_bgchain_ctmc
 * @see Solver_mam_bgchain
 */
public final class Mam_bgchain_env {
    private Mam_bgchain_env() {}

    /** Lumped environment of one station. */
    public static final class Env {
        /** Lumped generator, (me x me). */
        public Matrix A;
        /** Stationary probability of each environment state. */
        public double[] phi;
        /** Closed jobs each environment state stands for, ascending. */
        public int[] esup;
    }

    /**
     * @param bg the solved background chain
     * @param i  index of the station within the chain's support
     * @return the lumped environment
     */
    public static Env mam_bgchain_env(Mam_bgchain_ctmc.Result bg, int i) {
        int nstates = bg.nstates;
        TreeSet<Integer> distinct = new TreeSet<Integer>();
        for (int s = 0; s < nstates; s++) {
            distinct.add(Integer.valueOf(bg.totocc[s][i]));
        }
        List<Integer> levels = new ArrayList<Integer>(distinct);
        Collections.sort(levels);

        int meAll = levels.size();
        double[] wAll = new double[meAll];
        for (int s = 0; s < nstates; s++) {
            int pos = Collections.binarySearch(levels, Integer.valueOf(bg.totocc[s][i]));
            wAll[pos] += bg.pi.get(0, s);
        }

        List<Integer> keptLevels = new ArrayList<Integer>();
        List<Double> keptW = new ArrayList<Double>();
        for (int e = 0; e < meAll; e++) {
            if (wAll[e] > GlobalConstants.Zero) {
                keptLevels.add(levels.get(e));
                keptW.add(Double.valueOf(wAll[e]));
            }
        }

        Env env = new Env();
        if (keptLevels.isEmpty()) {
            // degenerate chain: the station never holds a closed job
            env.A = new Matrix(1, 1);
            env.phi = new double[] { 1.0 };
            env.esup = new int[] { 0 };
            return env;
        }

        int me = keptLevels.size();
        env.esup = new int[me];
        double[] w = new double[me];
        double wsum = 0;
        for (int e = 0; e < me; e++) {
            env.esup[e] = keptLevels.get(e).intValue();
            w[e] = keptW.get(e).doubleValue();
            wsum += w[e];
        }

        // level index of each state, -1 when its level was dropped
        int[] lvl = new int[nstates];
        for (int s = 0; s < nstates; s++) {
            lvl[s] = Collections.binarySearch(keptLevels, Integer.valueOf(bg.totocc[s][i]));
        }

        Matrix A = new Matrix(me, me);
        if (me > 1) {
            java.util.Iterator<jline.util.matrix.MatrixEntry> it = bg.Q.nonZeroIterator();
            while (it.hasNext()) {
                jline.util.matrix.MatrixEntry entry = it.next();
                int s = entry.row;
                int sp = entry.col;
                if (s == sp) continue;
                int e = lvl[s];
                int ep = lvl[sp];
                if (e < 0 || ep < 0 || e == ep) continue;
                double p = bg.pi.get(0, s);
                if (p == 0) continue;
                A.set(e, ep, A.get(e, ep) + p * entry.value);
            }
            for (int e = 0; e < me; e++) {
                if (w[e] <= 0) continue;
                double diag = 0;
                for (int ep = 0; ep < me; ep++) {
                    if (ep == e) continue;
                    double v = A.get(e, ep) / w[e];
                    A.set(e, ep, v);
                    diag += v;
                }
                A.set(e, e, -diag);
            }
        }
        env.A = A;
        env.phi = new double[me];
        for (int e = 0; e < me; e++) {
            env.phi[e] = wsum > 0 ? w[e] / wsum : 1.0 / me;
        }
        return env;
    }
}
