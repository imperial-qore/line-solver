/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.lqn;

import jline.api.pfqn.Pfqn_mwrbb;
import jline.lang.constant.CallType;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.LayeredNetworkStruct;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Majumdar-Woodside robust box bounds on throughput for a layered queueing
 * network (LQN), computed on the processor-contention model.
 *
 * <p>Builds a closed multiclass queueing network from a {@link LayeredNetworkStruct}
 * in which the stations are the processors (hosts) and the classes are the
 * reference-task call chains, then applies the Majumdar-Woodside robust box
 * bounds ({@link Pfqn_mwrbb}). The per-chain demand at each processor is the
 * total host demand executed on that processor during one cycle of the
 * reference task, obtained by traversing the activity/call graph and scaling by
 * the mean number of synchronous calls. The reference-task multiplicity is the
 * class population and its think time is the class think time.</p>
 *
 * <p>This generalizes the classical LQN "Type 1 throughput bound"
 * X_ref &lt;= mult/(Z + D_total) (the no-contention upper bound computed by
 * lqns -b) with the processor-utilization upper bound and the Majumdar-Woodside
 * lower bound (Theorem 2).</p>
 *
 * <p>Scope: processors are the queueing resources; reference-task chains are the
 * classes. Software (finite-thread task) bottlenecks and non-deterministic
 * activity precedence are not modeled (sequential execution assumed).</p>
 *
 * @since LINE 3.0
 */
public final class Lqn_boxbounds {

    private Lqn_boxbounds() {
    }

    /** Result holder. Vectors are 1-based (index 0 unused), length nidx+1. */
    public static final class Result {
        public int[] refidx;
        public double[] Xlo;
        public double[] Xup;
        public double[] TN_lo;
        public double[] TN_up;
        public double[] UN_lo;
        public double[] UN_up;
    }

    public static Result compute(LayeredNetworkStruct lqn) {
        int nidx = lqn.nidx;
        int nH = lqn.nhosts;

        List<Integer> refs = new ArrayList<Integer>();
        for (int t = 1; t <= lqn.ntasks; t++) {
            int tidx = lqn.tshift + t;
            if ((int) lqn.isref.get(tidx) != 0) {
                refs.add(tidx);
            }
        }
        int R = refs.size();

        double[][] D = new double[nH][R];               // 0-based host row, class col
        double[][] Vis = new double[nidx + 1][R];       // 1-based element row
        double[] Nref = new double[R];
        double[] Zref = new double[R];
        for (int r = 0; r < R; r++) {
            int tidx = refs.get(r);
            double m = lqn.mult.get(0, tidx);
            Nref[r] = (Double.isInfinite(m) || Double.isNaN(m)) ? 1.0 : m;
            Double z = lqn.think_mean.get(tidx);
            Zref[r] = (z == null || Double.isNaN(z)) ? 0.0 : z;
            if (lqn.entriesof.get(tidx) != null) {
                for (int eidx : lqn.entriesof.get(tidx)) {
                    visitEntry(lqn, eidx, 1.0, nH, r, D, Vis);
                }
            }
        }

        Matrix V = new Matrix(nH, R);
        Matrix S = new Matrix(nH, R);
        for (int h = 0; h < nH; h++) {
            for (int r = 0; r < R; r++) {
                S.set(h, r, D[h][r]);
                V.set(h, r, D[h][r] > 0 ? 1.0 : 0.0);
            }
        }
        Matrix Nm = new Matrix(1, R);
        Matrix Zm = new Matrix(1, R);
        Matrix schedm = new Matrix(nH, 1);
        Matrix priom = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            Nm.set(0, r, Nref[r]);
            Zm.set(0, r, Zref[r]);
        }
        for (int h = 0; h < nH; h++) {
            int hidx = lqn.hshift + h + 1;              // host absolute index
            schedm.set(h, 0, discCode(lqn.sched.get(hidx)));
        }

        Pfqn_mwrbb.Result rbb = Pfqn_mwrbb.pfqn_mwrbb(V, S, Nm, Zm, schedm, priom);

        Result out = new Result();
        out.refidx = new int[R];
        out.Xlo = new double[R];
        out.Xup = new double[R];
        for (int r = 0; r < R; r++) {
            out.refidx[r] = refs.get(r);
            out.Xlo[r] = rbb.Xlo.get(0, r);
            out.Xup[r] = rbb.Xup.get(0, r);
        }
        out.TN_lo = new double[nidx + 1];
        out.TN_up = new double[nidx + 1];
        out.UN_lo = new double[nidx + 1];
        out.UN_up = new double[nidx + 1];
        for (int i = 0; i <= nidx; i++) {
            out.TN_lo[i] = Double.NaN;
            out.TN_up[i] = Double.NaN;
            out.UN_lo[i] = Double.NaN;
            out.UN_up[i] = Double.NaN;
        }
        for (int i = 1; i <= nidx; i++) {
            double tlo = 0, tup = 0;
            boolean visited = false;
            for (int r = 0; r < R; r++) {
                if (Vis[i][r] > 0) {
                    visited = true;
                    tlo += out.Xlo[r] * Vis[i][r];
                    tup += out.Xup[r] * Vis[i][r];
                }
            }
            if (visited) {
                out.TN_lo[i] = tlo;
                out.TN_up[i] = tup;
            }
        }
        for (int h = 0; h < nH; h++) {
            int hidx = lqn.hshift + h + 1;
            double ulo = 0, uup = 0;
            for (int r = 0; r < R; r++) {
                ulo += out.Xlo[r] * D[h][r];
                uup += out.Xup[r] * D[h][r];
            }
            out.UN_lo[hidx] = ulo;
            out.UN_up[hidx] = uup;
        }
        return out;
    }

    private static void visitEntry(LayeredNetworkStruct lqn, int eidx, double mult,
                                   int nH, int r, double[][] D, double[][] Vis) {
        Vis[eidx][r] += mult;
        int tidx = (int) lqn.parent.get(0, eidx);       // task hosting this entry
        if (tidx >= 1 && tidx < Vis.length) {
            Vis[tidx][r] += mult;
        }
        if (lqn.actsof.get(eidx) == null) {
            return;
        }
        for (int aidx : lqn.actsof.get(eidx)) {
            if ((int) lqn.parent.get(0, aidx) == tidx) {   // activity of this entry
                visitActivity(lqn, aidx, mult, nH, r, D, Vis);
            }
        }
    }

    private static void visitActivity(LayeredNetworkStruct lqn, int aidx, double mult,
                                      int nH, int r, double[][] D, double[][] Vis) {
        Vis[aidx][r] += mult;
        int tidx = (int) lqn.parent.get(0, aidx);
        int hidx = (int) lqn.parent.get(0, tidx);         // host absolute index
        Double hd = lqn.hostdem_mean.get(aidx);
        double dem = (hd == null || Double.isNaN(hd)) ? 0.0 : hd;
        int hrow = hidx - lqn.hshift - 1;                 // 0-based host row
        if (hrow >= 0 && hrow < nH) {
            D[hrow][r] += mult * dem;
        }
        if (lqn.callsof.get(aidx) == null) {
            return;
        }
        for (int cidx : lqn.callsof.get(aidx)) {
            if (lqn.calltype.get(cidx) == CallType.SYNC) {
                Double cm = lqn.callproc_mean.get(cidx);
                double cmean = (cm == null || Double.isNaN(cm)) ? 0.0 : cm;
                int callee = (int) lqn.callpair.get(cidx, 2);
                visitEntry(lqn, callee, mult * cmean, nH, r, D, Vis);
            }
        }
    }

    private static int discCode(SchedStrategy s) {
        if (s == null) {
            return Pfqn_mwrbb.FIFO;
        }
        switch (s) {
            case FCFS:
                return Pfqn_mwrbb.FIFO;
            case PS:
            case DPS:
            case GPS:
            case PSPRIO:
            case DPSPRIO:
            case GPSPRIO:
                return Pfqn_mwrbb.PS;
            case HOL:
            case FCFSPRIO:
                return Pfqn_mwrbb.NPPRIO;
            case FCFSPRPRIO:
            case LCFSPRPRIO:
                return Pfqn_mwrbb.PPPRIO;
            default:
                return Pfqn_mwrbb.ABA;
        }
    }
}
