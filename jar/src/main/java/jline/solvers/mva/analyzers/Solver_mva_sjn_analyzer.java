package jline.solvers.mva.analyzers;

import jline.api.pfqn.mva.Pfqn_amvasjn;
import jline.api.pfqn.mva.Pfqn_mvasjn;
import jline.api.pfqn.mva.SjnOptions;
import jline.api.pfqn.mva.SjnStarvationException;
import jline.api.sn.SnDeaggregateChainResults;
import jline.api.sn.SnGetDemandsChain;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Closed networks with non-preemptive shortest-job-next (SJF) stations, wired as a SolverMVA
 * dispatch path. The station is modelled by the conditional waiting time equation of K. Kant,
 * "MVA approximations for SJN scheduling", Performance Evaluation 15(1):41-61, 1992, evaluated
 * either over the full population lattice ({@link Pfqn_mvasjn}) or through its Schweitzer fixed
 * point ({@link Pfqn_amvasjn}).
 *
 * <p>The lattice costs prod(N+1) steps, so 'default' switches to the fixed point once the lattice
 * exceeds the configured limit. Ask for 'exact' or 'mva' to force the lattice, 'amva' to force the
 * fixed point.</p>
 *
 * <p>Ported from matlab/src/solvers/MVA/solver_mva_sjn_analyzer.m.</p>
 */
public final class Solver_mva_sjn_analyzer {

    private Solver_mva_sjn_analyzer() {}

    /** Lattice size above which 'default' prefers the fixed point. */
    private static final double LATTICE_MAX = 1e5;

    public static MVAResult solver_mva_sjn_analyzer(NetworkStruct sn, SolverOptions options) {
        long startTime = System.nanoTime();
        int M = sn.nstations;
        int C = sn.nchains;

        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(sn);
        Matrix Lchain = dem.Dchain;
        Matrix STchain = dem.STchain;
        Matrix Vchain = dem.Vchain;
        Matrix alpha = dem.alpha;
        Matrix Nchain = dem.Nchain;
        Matrix SCVchain = dem.SCVchain;

        for (int c = 0; c < C; c++) {
            if (Double.isInfinite(Nchain.get(c))) {
                throw new RuntimeException("SJN scheduling is supported by SolverMVA only in closed"
                        + " models, the open case has no population recursion.");
            }
        }

        // the scheduling strategy alone selects delay against queue, nservers never does
        boolean[] isDelay = new boolean[M];
        boolean[] isSjn = new boolean[M];
        List<Integer> rows = new ArrayList<>();
        List<Integer> sjnList = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            double ns = sn.nservers.get(i);
            if (s == SchedStrategy.EXT) {
                continue;
            } else if (s == SchedStrategy.INF) {
                isDelay[i] = true;
            } else if (s == SchedStrategy.SJF) {
                if (ns != 1) {
                    throw new RuntimeException("SJN scheduling at station " + (i + 1) + " requires a"
                            + " single server, the response time equation is a single-server one.");
                }
                isSjn[i] = true;
                rows.add(i);
                sjnList.add(rows.size() - 1);
            } else if (s == SchedStrategy.PS || s == SchedStrategy.LCFSPR
                    || s == SchedStrategy.FCFS || s == SchedStrategy.SIRO) {
                if (ns != 1) {
                    throw new RuntimeException("station " + (i + 1) + " has "
                            + (Double.isFinite(ns) ? String.valueOf((int) ns) : "Inf")
                            + " servers, the SJN analyzer solves the remaining stations with the"
                            + " single-server MVA equation.");
                }
                rows.add(i);
            } else {
                throw new RuntimeException("The SJN analyzer does not support " + s
                        + " scheduling at the other stations.");
            }
        }

        int Mq = rows.size();
        int[] queueRows = new int[Mq];
        for (int j = 0; j < Mq; j++) {
            queueRows[j] = rows.get(j);
        }
        int[] sjnrows = new int[sjnList.size()];
        for (int j = 0; j < sjnrows.length; j++) {
            sjnrows[j] = sjnList.get(j);
        }

        Matrix L = new Matrix(Mq, C);
        Matrix V = new Matrix(Mq, C);
        Matrix scv = new Matrix(Mq, C);
        for (int j = 0; j < Mq; j++) {
            int i = queueRows[j];
            for (int c = 0; c < C; c++) {
                L.set(j, c, STchain.get(i, c) * Vchain.get(i, c));
                V.set(j, c, Vchain.get(i, c));
                double v = 1.0;
                if (isSjn[i] && Double.isFinite(SCVchain.get(i, c)) && SCVchain.get(i, c) > 0) {
                    v = SCVchain.get(i, c);
                }
                scv.set(j, c, v);
            }
        }
        Matrix Z = new Matrix(1, C);
        Matrix Nrow = new Matrix(1, C);
        for (int c = 0; c < C; c++) {
            double z = 0;
            for (int i = 0; i < M; i++) {
                if (isDelay[i]) {
                    z += STchain.get(i, c) * Vchain.get(i, c);
                }
            }
            Z.set(0, c, z);
            Nrow.set(0, c, Nchain.get(c));
        }

        SjnOptions sjnopt = new SjnOptions();
        if (options.iter_tol > 0) {
            sjnopt.tol = options.iter_tol;
        }
        if (options.iter_max > 0) {
            sjnopt.iterMax = options.iter_max;
        }
        if (options.config != null) {
            if (options.config.sjn_ns != null) {
                sjnopt.ns = options.config.sjn_ns;
            }
            if (options.config.sjn_lfactor != null) {
                sjnopt.Lfactor = options.config.sjn_lfactor;
            }
            if (options.config.sjn_umax != null) {
                sjnopt.umax = options.config.sjn_umax;
            }
        }
        // SJN applies within a class and the classes are then non-preemptively prioritised;
        // without distinct priorities the jobs of every class are compared by size directly
        if (sn.nchains == sn.nclasses && hasDistinctPriorities(sn)) {
            int[] prio = new int[C];
            for (int c = 0; c < C; c++) {
                prio[c] = (int) sn.classprio.get(c);
            }
            sjnopt.prio = prio;
        }

        double lattice = 1;
        for (int c = 0; c < C; c++) {
            lattice *= Nchain.get(c) + 1;
        }
        double latticemax = LATTICE_MAX;
        if (options.config != null && options.config.sjn_lattice_max != null) {
            latticemax = options.config.sjn_lattice_max.doubleValue();
        }
        boolean uselattice;
        String method = options.method == null ? "default" : options.method;
        if (method.equals("amva") || method.equals("bs") || method.equals("sjn.amva")) {
            uselattice = false;
        } else if (method.equals("exact") || method.equals("mva") || method.equals("sjn.mva")) {
            uselattice = true;
        } else {
            uselattice = lattice <= latticemax;
        }

        Pfqn_mvasjn.Result sjn;
        String actualmethod;
        if (uselattice) {
            try {
                sjn = Pfqn_mvasjn.pfqn_mvasjn(L, Nrow, Z, scv, sjnrows, V, sjnopt);
                actualmethod = "sjn.mva";
            } catch (SjnStarvationException e) {
                if (!method.equals("default")) {
                    throw e;
                }
                sjn = Pfqn_amvasjn.pfqn_amvasjn(L, Nrow, Z, scv, sjnrows, V, sjnopt);
                actualmethod = "sjn.amva";
            }
        } else {
            sjn = Pfqn_amvasjn.pfqn_amvasjn(L, Nrow, Z, scv, sjnrows, V, sjnopt);
            actualmethod = "sjn.amva";
        }

        Matrix Xchain = new Matrix(1, C);
        Matrix Qchain = new Matrix(M, C);
        Matrix Uchain = new Matrix(M, C);
        Matrix Rchain = new Matrix(M, C);
        Matrix Tchain = new Matrix(M, C);
        for (int c = 0; c < C; c++) {
            Xchain.set(0, c, sjn.X.get(0, c));
        }
        for (int j = 0; j < Mq; j++) {
            int i = queueRows[j];
            for (int c = 0; c < C; c++) {
                Qchain.set(i, c, sjn.Q.get(j, c));
                Uchain.set(i, c, sjn.U.get(j, c));
            }
        }
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < C; c++) {
                Tchain.set(i, c, Xchain.get(0, c) * Vchain.get(i, c));
            }
        }
        for (int i = 0; i < M; i++) {
            if (isDelay[i]) {
                for (int c = 0; c < C; c++) {
                    double ts = Tchain.get(i, c) * STchain.get(i, c);
                    Qchain.set(i, c, ts);
                    Uchain.set(i, c, ts);
                }
            }
        }
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < C; c++) {
                double t = Tchain.get(i, c);
                if (t > 0) {
                    Rchain.set(i, c, Qchain.get(i, c) / t);
                }
            }
        }

        for (int c = 0; c < C; c++) {
            if (!Double.isFinite(Xchain.get(0, c))) {
                Xchain.set(0, c, 0);
            }
            for (int i = 0; i < M; i++) {
                if (!Double.isFinite(Qchain.get(i, c))) {
                    Qchain.set(i, c, 0);
                }
                if (!Double.isFinite(Uchain.get(i, c))) {
                    Uchain.set(i, c, 0);
                }
                if (!Double.isFinite(Rchain.get(i, c))) {
                    Rchain.set(i, c, 0);
                }
            }
            // an empty chain carries no jobs, so every one of its metrics is zero
            if (Nchain.get(c) == 0) {
                Xchain.set(0, c, 0);
                for (int i = 0; i < M; i++) {
                    Qchain.set(i, c, 0);
                    Uchain.set(i, c, 0);
                    Rchain.set(i, c, 0);
                    Tchain.set(i, c, 0);
                }
            }
        }

        // MATLAB passes [] here and lets the deaggregation rebuild Q and U from Rchain and alpha
        Ret.snDeaggregateChainResults dre = SnDeaggregateChainResults.snDeaggregateChainResults(
                sn, Lchain, null, STchain, Vchain, alpha, null, null, Rchain, Tchain, null,
                Xchain);

        MVAResult res = new MVAResult();
        res.QN = dre.Q;
        res.UN = dre.U;
        res.RN = dre.R;
        res.TN = dre.T;
        res.CN = dre.C;
        res.XN = dre.X;
        res.AN = new Matrix(0, 0);
        res.WN = new Matrix(0, 0);
        res.logNormConstAggr = Double.NaN;
        res.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        res.iter = sjn.iter;
        res.method = actualmethod;
        return res;
    }

    /** True when every class carries a different priority level, as Kant's method A requires. */
    private static boolean hasDistinctPriorities(NetworkStruct sn) {
        if (sn.classprio == null || sn.classprio.isEmpty()) {
            return false;
        }
        int R = sn.nclasses;
        for (int i = 0; i < R; i++) {
            for (int j = i + 1; j < R; j++) {
                if (sn.classprio.get(i) == sn.classprio.get(j)) {
                    return false;
                }
            }
        }
        return true;
    }
}
