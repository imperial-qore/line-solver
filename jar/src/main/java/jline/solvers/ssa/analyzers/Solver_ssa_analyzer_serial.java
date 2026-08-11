package jline.solvers.ssa.analyzers;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.api.mam.Map_mean;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.EventCache;
import jline.solvers.SolverOptions;
import jline.solvers.ssa.SSAResult;
import jline.solvers.ssa.SSAValues;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.ssa.handlers.Solver_ssa;
import jline.util.matrix.Matrix;

public final class Solver_ssa_analyzer_serial {
    private Solver_ssa_analyzer_serial() {}

    public static SSAResult solver_ssa_analyzer_serial(NetworkStruct sn, boolean hash,
                                                       Map<StatefulNode, Matrix> init_state,
                                                       SolverOptions options, SolverSSA solverSSA) {
        int M = sn.nstations;
        int K = sn.nclasses;

        Matrix S = sn.nservers;
        Matrix NK = sn.njobs.transpose();
        Map<jline.lang.nodes.Station, SchedStrategy> schedid = sn.sched;

        Map<jline.lang.nodes.Station, Map<jline.lang.JobClass, jline.util.matrix.MatrixCell>> PH = sn.proc;
        Map<Integer, Matrix> tranSysState = new HashMap<Integer, Matrix>();
        Matrix tranSync = new Matrix(0, 0);

        Matrix XN = new Matrix(1, K);
        XN.fill(Double.NaN);
        Matrix UN = new Matrix(M, K);
        UN.fill(Double.NaN);
        Matrix QN = new Matrix(M, K);
        QN.fill(Double.NaN);
        Matrix RN = new Matrix(M, K);
        RN.fill(Double.NaN);
        Matrix TN = new Matrix(M, K);
        TN.fill(Double.NaN);
        Matrix CN = new Matrix(1, K);
        CN.fill(Double.NaN);

        options.samples++;
        solverSSA.eventCache = new EventCache(false, options.config.eventcache);
        // see _kb/06-solver-catalog.md for rationale
        Matrix userCap = sn.cap.copy();
        Matrix userClasscap = sn.classcap.copy();
        SSAValues result = Solver_ssa.solver_ssa(sn, solverSSA.eventCache, init_state, options, solverSSA);
        Matrix probSysState = result.pi;
        Matrix StateSpaceAggr = result.SSq;
        Map<Integer, Matrix> arvRates = result.arvRates;
        Map<Integer, Matrix> depRates = result.depRates;
        tranSysState = result.tranSysState;
        tranSync = result.tranSync;

        for (int k = 0; k < K; k++) {
            int refsf = (int) sn.stationToStateful.get((int) sn.refstat.get(k));
            Matrix departure = depRates.get(k);
            Matrix dep_wset_refsf = new Matrix(StateSpaceAggr.getNumRows(), 1);
            for (int i = 0; i < StateSpaceAggr.getNumRows(); i++) {
                dep_wset_refsf.set(i, 0, departure.get(i, refsf));
            }
            // toDouble call since 1xn mult nx1
            XN.set(k, probSysState.mult(dep_wset_refsf).toDouble());
        }

        for (int i = 0; i < M; i++) {
            int isf = (int) sn.stationToStateful.get(i);
            for (int k = 0; k < K; k++) {
                Matrix departure = depRates.get(k);
                Matrix dep_wset_isf = new Matrix(StateSpaceAggr.getNumRows(), 1);
                for (int j = 0; j < StateSpaceAggr.getNumRows(); j++) {
                    dep_wset_isf.set(j, 0, departure.get(j, isf));
                }
                TN.set(i, k, probSysState.mult(dep_wset_isf).toDouble());

                Matrix ssaggr_wset_isf =
                        Matrix.extract(StateSpaceAggr, 0, StateSpaceAggr.getNumRows(), i * K + k, i * K + k + 1);
                QN.set(i, k, probSysState.mult(ssaggr_wset_isf).toDouble());
            }
            // see _kb/06-solver-catalog.md for rationale
            boolean stationCapFinite = !Double.isInfinite(userCap.get(i)) && userCap.get(i) < Integer.MAX_VALUE;
            boolean[] canDropClass = new boolean[K];
            for (int r = 0; r < K; r++) {
                boolean classCapFinite = !Double.isInfinite(userClasscap.get(i, r)) && userClasscap.get(i, r) < Integer.MAX_VALUE;
                canDropClass[r] = Double.isInfinite(sn.njobs.get(r)) && (stationCapFinite || classCapFinite);
            }
            SchedStrategy sched = schedid.get(sn.stations.get(i));
            if (sched == SchedStrategy.INF) {
                int k = 0;
                while (k < K) {
                    UN.set(i, k, QN.get(i, k));
                    k++;
                }
            } else if (sched == SchedStrategy.PS || sched == SchedStrategy.DPS || sched == SchedStrategy.GPS
                    || sched == SchedStrategy.PSPRIO || sched == SchedStrategy.DPSPRIO || sched == SchedStrategy.GPSPRIO) {
                if ((sn.lldscaling == null || sn.lldscaling.isEmpty())
                        && (sn.cdscaling == null || sn.cdscaling.isEmpty())
                        && (sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                    int k = 0;
                    while (k < K) {
                        // FJ tag-augmented structs: skip zero-visit classes (see
                        // the FCFS branch below for the rationale)
                        boolean fjVisitZeroPs = false;
                        if (sn.isfjaugmented) {
                            for (int cc = 0; cc < sn.chains.getNumRows(); cc++) {
                                if (sn.chains.get(cc, k) > 0) {
                                    fjVisitZeroPs = (sn.visits.get(cc).get(isf, k) == 0.0);
                                    break;
                                }
                            }
                        }
                        if (!fjVisitZeroPs && !PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).isEmpty()) {
                            Matrix arrival = arvRates.get(k);
                            Matrix arv_wset_isf = new Matrix(StateSpaceAggr.getNumRows(), 1);
                            int c = 0;
                            while (c < StateSpaceAggr.getNumRows()) {
                                arv_wset_isf.set(c, 0, arrival.get(c, isf));
                                c++;
                            }
                            // see _kb/06-solver-catalog.md for rationale
                            if (canDropClass[k] || sn.isfjaugmented) {
                                UN.set(i, k, TN.get(i, k) / sn.rates.get(i, k) / S.get(i));
                            } else {
                                UN.set(i, k, probSysState.mult(arv_wset_isf).toDouble() / sn.rates.get(i, k) / S.get(i));
                            }
                        }
                        k++;
                    }
                } else {
                    // see _kb/06-solver-catalog.md for rationale
                    double ceffLd = S.get(i);
                    if (sn.lldscaling != null && !sn.lldscaling.isEmpty() && i < sn.lldscaling.getNumRows()) {
                        for (int cidx = 0; cidx < sn.lldscaling.getNumCols(); cidx++) {
                            ceffLd = Math.max(ceffLd, sn.lldscaling.get(i, cidx));
                        }
                    }
                    Matrix cdPeakVec = (sn.cdscaling != null && sn.cdscaling.get(sn.stations.get(i)) != null
                            && sn.cdscalingpeak != null) ? sn.cdscalingpeak.get(sn.stations.get(i)) : null;
                    Matrix jdPeakVec = (sn.jdscaling != null && sn.jdscaling.get(sn.stations.get(i)) != null
                            && sn.jdscalingpeak != null) ? sn.jdscalingpeak.get(sn.stations.get(i)) : null;
                    int col = 0;
                    while (col < K) {
                        if (!PH.get(sn.stations.get(i)).get(sn.jobclasses.get(col)).isEmpty()) {
                            double meanSvcLd = Map_mean.map_mean(
                                    PH.get(sn.stations.get(i)).get(sn.jobclasses.get(col)).get(0),
                                    PH.get(sn.stations.get(i)).get(sn.jobclasses.get(col)).get(1));
                            // effective peak = product of the declared class- and
                            // joint-dependence peaks; fall back to ceffLd only when neither is set.
                            double cdiv;
                            if (cdPeakVec != null || jdPeakVec != null) {
                                cdiv = 1.0;
                                if (cdPeakVec != null) cdiv *= cdPeakVec.get(0, col);
                                if (jdPeakVec != null) cdiv *= jdPeakVec.get(0, col);
                            } else {
                                cdiv = ceffLd;
                            }
                            if (cdiv > 0) {
                                UN.set(i, col, TN.get(i, col) * meanSvcLd / cdiv);
                            } else {
                                UN.set(i, col, 0.0);
                            }
                        }
                        col++;
                    }
                }
            } else {
                if ((sn.lldscaling == null || sn.lldscaling.isEmpty())
                        && (sn.cdscaling == null || sn.cdscaling.isEmpty())
                        && (sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                    int k = 0;
                    while (k < K) {
                        // see _kb/06-solver-catalog.md for rationale
                        boolean fjVisitZero = false;
                        if (sn.isfjaugmented) {
                            for (int cc = 0; cc < sn.chains.getNumRows(); cc++) {
                                if (sn.chains.get(cc, k) > 0) {
                                    fjVisitZero = (sn.visits.get(cc).get(isf, k) == 0.0);
                                    break;
                                }
                            }
                        }
                        if (!fjVisitZero && !PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).isEmpty()) {
                            Matrix arrival = arvRates.get(k);
                            Matrix arv_wset_isf = new Matrix(StateSpaceAggr.getNumRows(), 1);
                            int c = 0;
                            while (c < StateSpaceAggr.getNumRows()) {
                                arv_wset_isf.set(c, 0, arrival.get(c, isf));
                                c++;
                            }
                            double map_mean = Map_mean.map_mean(
                                    PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(0),
                                    PH.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).get(1)) / S.get(i);
                            if (canDropClass[k] || sn.isfjaugmented) {
                                UN.set(i, k, TN.get(i, k) * map_mean);
                            } else {
                                UN.set(i, k, probSysState.mult(arv_wset_isf).toDouble() * map_mean);
                            }
                        }
                        k++;
                    }
                } else {
                    // see _kb/06-solver-catalog.md for rationale
                    double ceffLd = S.get(i);
                    if (sn.lldscaling != null && !sn.lldscaling.isEmpty() && i < sn.lldscaling.getNumRows()) {
                        for (int cidx = 0; cidx < sn.lldscaling.getNumCols(); cidx++) {
                            ceffLd = Math.max(ceffLd, sn.lldscaling.get(i, cidx));
                        }
                    }
                    Matrix cdPeakVec = (sn.cdscaling != null && sn.cdscaling.get(sn.stations.get(i)) != null
                            && sn.cdscalingpeak != null) ? sn.cdscalingpeak.get(sn.stations.get(i)) : null;
                    Matrix jdPeakVec = (sn.jdscaling != null && sn.jdscaling.get(sn.stations.get(i)) != null
                            && sn.jdscalingpeak != null) ? sn.jdscalingpeak.get(sn.stations.get(i)) : null;
                    int col = 0;
                    while (col < K) {
                        if (!PH.get(sn.stations.get(i)).get(sn.jobclasses.get(col)).isEmpty()) {
                            double meanSvcLd = Map_mean.map_mean(
                                    PH.get(sn.stations.get(i)).get(sn.jobclasses.get(col)).get(0),
                                    PH.get(sn.stations.get(i)).get(sn.jobclasses.get(col)).get(1));
                            // effective peak = product of the declared class- and
                            // joint-dependence peaks; fall back to ceffLd only when neither is set.
                            double cdiv;
                            if (cdPeakVec != null || jdPeakVec != null) {
                                cdiv = 1.0;
                                if (cdPeakVec != null) cdiv *= cdPeakVec.get(0, col);
                                if (jdPeakVec != null) cdiv *= jdPeakVec.get(0, col);
                            } else {
                                cdiv = ceffLd;
                            }
                            if (cdiv > 0) {
                                UN.set(i, col, TN.get(i, col) * meanSvcLd / cdiv);
                            } else {
                                UN.set(i, col, 0.0);
                            }
                        }
                        col++;
                    }
                }
            }
        }

        for (int k = 0; k < K; k++) {
            for (int i = 0; i < M; i++) {
                if (TN.get(i, k) > 0) {
                    RN.set(i, k, QN.get(i, k) / TN.get(i, k));
                } else {
                    RN.set(i, k, 0);
                }
            }
            CN.set(k, NK.get(k) / XN.get(k));
        }

        // update routing probabilities in nodes with state-dependent routing
        Matrix TNcache = new Matrix(sn.nstateful, K);
        Matrix XNcache = new Matrix(sn.nstateful, K);
        for (int k = 0; k < K; k++) {
            for (int isf = 0; isf < sn.nstateful; isf++) {
                if (sn.nodetype.get(isf) == NodeType.Cache) {
                    double TNcacheValue = probSysState.mult(depRates.get(k).getColumn(isf)).get(0);
                    double XNcacheValue = probSysState.mult(arvRates.get(k).getColumn(isf)).get(0);
                    TNcache.set(isf, k, TNcacheValue);
                    XNcache.set(isf, k, XNcacheValue);
                }
            }
        }

        boolean retrievalLatencyWarned = false;
        for (int k = 0; k < K; k++) {
            for (int isf = 0; isf < sn.nstateful; isf++) {
                // index the struct's own stateful list (on FJ tag-augmented
                // structs it differs from the original model's)
                StatefulNode statefulNode = sn.stateful.get(isf);
                if (statefulNode instanceof Cache) {
                    CacheNodeParam cacheNp = (CacheNodeParam) sn.nodeparam.get(statefulNode);
                    if (cacheNp.hitclass.getNumCols() > k) {
                        int h = (int) cacheNp.hitclass.get(k);
                        int m = (int) cacheNp.missclass.get(k);

                        if (cacheNp.actualhitprob == null) {
                            cacheNp.actualhitprob = new Matrix(1, K);
                            cacheNp.actualhitprob.fill(Double.NaN);
                        }
                        if (cacheNp.actualmissprob == null) {
                            cacheNp.actualmissprob = new Matrix(1, K);
                            cacheNp.actualmissprob.fill(Double.NaN);
                        }

                        double actualmissprobValue = Double.NaN;
                        double actualhitprobValue = Double.NaN;
                        if (h != -1 && m != -1) {
                            actualhitprobValue = TNcache.get(isf, h) / (TNcache.get(isf, h) + TNcache.get(isf, m));
                            actualmissprobValue = TNcache.get(isf, m) / (TNcache.get(isf, h) + TNcache.get(isf, m));
                        }
                        cacheNp.actualhitprob.set(k, actualhitprobValue);
                        cacheNp.actualmissprob.set(k, actualmissprobValue);

                        // see _kb/06-solver-catalog.md for rationale
                        double actualresidt = Double.NaN;
                        List<Integer> retrievalSystemQueueIndices =
                                (cacheNp.retrievalSystemQueueIndices != null)
                                        ? cacheNp.retrievalSystemQueueIndices.get(k) : null;
                        if (retrievalSystemQueueIndices != null && !retrievalSystemQueueIndices.isEmpty()) {
                            if (!retrievalLatencyWarned) {
                                InputOutput.line_warning("solver_ssa_analyzer_serial",
                                        "Retrieval-system expected latency is not currently implemented; "
                                        + "reporting NaN.");
                                retrievalLatencyWarned = true;
                            }
                        }

                        if (cacheNp.actualresidt == null) {
                            cacheNp.actualresidt = new Matrix(1, K);
                            cacheNp.actualresidt.fill(Double.NaN);
                        }
                        cacheNp.actualresidt.set(k, actualresidt);
                    }
                }
            }
        }

        // matrices QN, CN, RN, UN, XN, TN, where they are Double.isNan set to 0
        QN.apply(Double.NaN, 0.0, "equal");
        CN.apply(Double.NaN, 0.0, "equal");
        RN.apply(Double.NaN, 0.0, "equal");
        UN.apply(Double.NaN, 0.0, "equal");
        XN.apply(Double.NaN, 0.0, "equal");
        TN.apply(Double.NaN, 0.0, "equal");

        return new SSAResult(QN, UN, RN, TN, CN, XN, tranSysState, tranSync, sn);
    }
}
