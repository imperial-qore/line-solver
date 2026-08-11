/**
 * @file Parallel SSA analyzer
 *
 * @since LINE 3.0
 */
package jline.solvers.ssa.analyzers;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.util.FastMath;

import jline.api.mam.*;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.EventCache;
import jline.solvers.SolverOptions;
import jline.solvers.ssa.SSAResult;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.ssa.handlers.Solver_ssa;
import jline.solvers.ssa.SSAValues;

import jline.util.RandomManager;
import jline.util.matrix.Matrix;

public final class Solver_ssa_analyzer_parallel {
    private Solver_ssa_analyzer_parallel() {}

    private static Matrix averageMatrices(Collection<Matrix> matrices, int numThreads) {
        if (matrices.isEmpty()) {
            throw new IllegalArgumentException("Cannot average empty collection of matrices");
        }
        Matrix first = matrices.iterator().next();
        Matrix result = new Matrix(first.getNumRows(), first.getNumCols());
        result.zero();
        for (Matrix m : matrices) {
            result.addEq(m);
        }
        return Matrix.scaleMult(result, 1.0 / (double) numThreads);
    }

    public static SSAResult solver_ssa_analyzer_parallel(final NetworkStruct sn,
                                                         final Map<StatefulNode, Matrix> init_state,
                                                         final SolverOptions options,
                                                         final SolverSSA solverSSA) {
        final int M = sn.nstations;
        final int K = sn.nclasses;
        final Map<?, ?> PH = sn.proc;
        final Matrix S = sn.nservers;
        final Matrix NK = sn.njobs.transpose();

        final Map<Integer, Matrix> QNs = new ConcurrentHashMap<Integer, Matrix>();
        final Map<Integer, Matrix> UNs = new ConcurrentHashMap<Integer, Matrix>();
        final Map<Integer, Matrix> RNs = new ConcurrentHashMap<Integer, Matrix>();
        final Map<Integer, Matrix> TNs = new ConcurrentHashMap<Integer, Matrix>();
        final Map<Integer, Matrix> CNs = new ConcurrentHashMap<Integer, Matrix>();
        final Map<Integer, Matrix> XNs = new ConcurrentHashMap<Integer, Matrix>();

        options.samples = (int) FastMath.ceil(options.samples / (double) solverSSA.numThreads);
        solverSSA.eventCache = new EventCache(true, options.config.eventcache);

        // Snapshot the USER capacities before Solver_ssa's preamble rewrites
        // sn.cap/classcap in place (it folds the state-space cutoff into them,
        // making every open class look capacity-constrained). Unbounded capacity
        // is stored as Integer.MAX_VALUE (Station default).
        final Matrix userCap = sn.cap.copy();
        final Matrix userClasscap = sn.classcap.copy();

        for (int t = 0; t < solverSSA.numThreads; t++) {
            final int threadIndex = t;
            solverSSA.threadPool.submit(new Runnable() {
                @Override
                public void run() {
                    java.util.Random threadRandom = RandomManager.getParallelRandomAsRandom("SSA_PARALLEL", threadIndex);
                    Matrix probSysState = new Matrix(0, 0);
                    Matrix SSq = new Matrix(0, 0);
                    Map<Integer, Matrix> arvRates = new HashMap<Integer, Matrix>();
                    Map<Integer, Matrix> depRates = new HashMap<Integer, Matrix>();

                    if ("ssa.parallel".equals(options.method) || "parallel".equals(options.method)) {
                        SSAValues result = Solver_ssa.solver_ssa(sn, solverSSA.eventCache, init_state, options, solverSSA);
                        probSysState = result.pi;
                        SSq = result.SSq;
                        arvRates = result.arvRates;
                        depRates = result.depRates;
                    }

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

                    for (int k = 0; k < K; k++) {
                        int refsf = (int) sn.stationToStateful.get((int) sn.refstat.get(k));
                        Matrix departure = depRates.get(k);
                        Matrix dep_wset_refsf = Matrix.extractColumn(departure, refsf, null);
                        XN.set(k, probSysState.mult(dep_wset_refsf).toDouble());
                        for (int i = 0; i < M; i++) {
                            int isf = (int) sn.stationToStateful.get(i);

                            Matrix dep_isf = Matrix.extractColumn(departure, isf, null);
                            TN.set(i, k, probSysState.mult(dep_isf).toDouble());

                            Matrix ssq_extracted = Matrix.extractColumn(SSq, i * K + k, null);
                            QN.set(i, k, probSysState.mult(ssq_extracted).toDouble());

                            SchedStrategy sched = (SchedStrategy) ((Map<?, ?>) sn.sched).get(sn.stations.get(i));
                            if (sched == SchedStrategy.INF) {
                                UN.set(i, k, QN.get(i, k));
                            } else {
                                Map<?, ?> phStation = (Map<?, ?>) PH.get(sn.stations.get(i));
                                Object phEntry = phStation == null ? null : phStation.get(sn.jobclasses.get(k));
                                if (phEntry != null && !((jline.util.matrix.MatrixCell) phEntry).isEmpty()) {
                                    jline.util.matrix.MatrixCell ph = (jline.util.matrix.MatrixCell) phEntry;
                                    double map_mean = Map_mean.map_mean(ph.get(0), ph.get(1)) / S.get(i);
                                    // see _kb/06-solver-catalog.md for rationale
                                    boolean stationCapFinite = !Double.isInfinite(userCap.get(i)) && userCap.get(i) < Integer.MAX_VALUE;
                                    boolean classCapFinite = !Double.isInfinite(userClasscap.get(i, k)) && userClasscap.get(i, k) < Integer.MAX_VALUE;
                                    if (Double.isInfinite(sn.njobs.get(k)) && (stationCapFinite || classCapFinite)) {
                                        UN.set(i, k, TN.get(i, k) * map_mean);
                                    } else {
                                        Matrix arrival = arvRates.get(k);
                                        Matrix arv_ik = Matrix.extractColumn(arrival, i, null);
                                        UN.set(i, k, probSysState.mult(arv_ik).toDouble() * map_mean);
                                    }
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
                    QN.apply(Double.NaN, 0.0, "equal");
                    CN.apply(Double.NaN, 0.0, "equal");
                    RN.apply(Double.NaN, 0.0, "equal");
                    UN.apply(Double.NaN, 0.0, "equal");
                    XN.apply(Double.NaN, 0.0, "equal");
                    TN.apply(Double.NaN, 0.0, "equal");

                    QNs.put(threadIndex, QN);
                    CNs.put(threadIndex, CN);
                    RNs.put(threadIndex, RN);
                    UNs.put(threadIndex, UN);
                    XNs.put(threadIndex, XN);
                    TNs.put(threadIndex, TN);
                }
            });
        }
        solverSSA.threadPool.shutdown();
        try {
            solverSSA.threadPool.awaitTermination(600, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }

        Matrix XN = averageMatrices(XNs.values(), solverSSA.numThreads);
        Matrix UN = averageMatrices(UNs.values(), solverSSA.numThreads);
        Matrix QN = averageMatrices(QNs.values(), solverSSA.numThreads);
        Matrix RN = averageMatrices(RNs.values(), solverSSA.numThreads);
        Matrix TN = averageMatrices(TNs.values(), solverSSA.numThreads);
        Matrix CN = averageMatrices(CNs.values(), solverSSA.numThreads);

        SSAResult res = new SSAResult(QN, UN, RN, TN, CN, XN, null, null, sn);
        res.method = "parallel";
        return res;
    }
}
