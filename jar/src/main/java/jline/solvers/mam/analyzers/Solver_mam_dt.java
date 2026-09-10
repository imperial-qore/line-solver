package jline.solvers.mam.analyzers;

import java.util.ArrayList;
import java.util.List;

import jline.api.mam.Dmap_batch;
import jline.api.mam.Dph_from_dist;
import jline.api.mam.Dph_to_dmap;
import jline.api.mam.Mg1_dt_queue;
import jline.api.sn.SnIsOpenModel;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Station;
import jline.lib.qmam.DTQueueResult;
import jline.lib.qmam.MAPMAP1Options;
import jline.lib.qmam.Q_DT_MAP_MAP_1;
import jline.lib.qmam.Q_DT_PH_PH_1;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAMResult;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Discrete-time (slotted) analysis of an open network whose interarrival and
 * service laws all live on the slot lattice.
 *
 * <p>A single queueing station is solved EXACTLY by the Q-MAM discrete-time
 * algorithms, {@link Q_DT_PH_PH_1} when both laws are renewal discrete
 * phase-type and {@link Q_DT_MAP_MAP_1} when either side is a DMAP. Several
 * stations are solved by a discrete-time parametric decomposition, which is an
 * approximation.
 *
 * <p>Time is measured in slots internally and converted back on exit, so QN and
 * UN are dimensionless, TN is per time unit and RN is in time units.
 *
 * <p>Convention: late arrival system with delayed access (LAS-DA), matching the
 * Q-MAM discrete-time queues and the LDES slotted engine.
 *
 * <p>MATLAB twin: solver_mam_dt.m
 */
public final class Solver_mam_dt {
    private Solver_mam_dt() {}

    public static MAMResult solver_mam_dt(NetworkStruct sn, SolverOptions options, double slotLength) {
        assertScope(sn);

        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix RN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        Matrix CN = new Matrix(1, K);
        Matrix XN = new Matrix(1, K);

        MatrixCell[] law = new MatrixCell[M];
        for (int ist = 0; ist < M; ist++) {
            law[ist] = stationLaw(sn, ist, 0, slotLength);
        }

        int sourceIdx = -1;
        List<Integer> queueIdx = new ArrayList<Integer>();
        for (int ist = 0; ist < M; ist++) {
            if (sn.sched.get(sn.stations.get(ist)) == SchedStrategy.EXT) {
                sourceIdx = ist;
            } else {
                queueIdx.add(Integer.valueOf(ist));
            }
        }
        if (sourceIdx < 0) {
            throw new RuntimeException("The discrete-time path requires an open model with a Source.");
        }

        int maxNumComp = intConfig(options, "dt_maxlevel", 1000);
        int spaceMax = intConfig(options, "space_max", 128);

        double[] QNq = new double[queueIdx.size()];
        double[] UNq = new double[queueIdx.size()];
        double[] TNq = new double[queueIdx.size()];
        String method;
        int totiter = 1;

        if (queueIdx.size() == 1) {
            DTQueueResult r = solveSingleStation(law[sourceIdx], law[queueIdx.get(0)], maxNumComp);
            QNq[0] = r.getMeanQueueLength();
            UNq[0] = r.getUtilization();
            TNq[0] = Dmap_batch.dmap_lambda(law[sourceIdx]);
            method = "dt.qmam";
        } else {
            totiter = solveNetwork(sn, law, sourceIdx, queueIdx, options, spaceMax, maxNumComp,
                    QNq, UNq, TNq);
            method = "dt.dec";
        }

        double lambdaSlot = Dmap_batch.dmap_lambda(law[sourceIdx]);
        for (int idx = 0; idx < queueIdx.size(); idx++) {
            int ist = queueIdx.get(idx).intValue();
            QN.set(ist, 0, QNq[idx]);
            UN.set(ist, 0, UNq[idx]);
            TN.set(ist, 0, TNq[idx] / slotLength);
            if (TNq[idx] > 0) {
                RN.set(ist, 0, QNq[idx] / TNq[idx] * slotLength);
            }
        }
        TN.set(sourceIdx, 0, lambdaSlot / slotLength);
        XN.set(0, 0, lambdaSlot / slotLength);
        CN.set(0, 0, QN.elementSum() / (lambdaSlot / slotLength));

        MAMResult result = new MAMResult();
        result.QN = QN;
        result.UN = UN;
        result.RN = RN;
        result.TN = TN;
        result.CN = CN;
        result.XN = XN;
        result.iter = totiter;
        result.method = method;
        return result;
    }

    /** Rejects the model features the discrete-time path cannot represent. */
    private static void assertScope(NetworkStruct sn) {
        if (!SnIsOpenModel.snIsOpenModel(sn)) {
            throw new RuntimeException("The discrete-time path supports open models only. A closed "
                    + "slotted model needs a level-dependent discrete chain, which the Q-MAM "
                    + "discrete-time catalogue does not cover.");
        }
        if (sn.nclasses > 1) {
            throw new RuntimeException("The discrete-time path supports one class only. Independent "
                    + "per-class lattice sources fire in the same slot with positive probability, and "
                    + "a batch of simultaneous arrivals of different classes is not an MMAP[K], which "
                    + "is what Q_DT_MMAPK_PHK_1 consumes.");
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            if (sched == SchedStrategy.EXT) {
                continue;
            }
            if (sched != SchedStrategy.FCFS) {
                throw new RuntimeException("Station " + ist + " uses scheduling " + sched
                        + ". The discrete-time path supports FCFS single-server stations and the Source only.");
            }
            if (sn.nservers.get(ist, 0) > 1) {
                throw new RuntimeException("Station " + ist + " has " + sn.nservers.get(ist, 0)
                        + " servers. The discrete-time path models one server per station: a slotted "
                        + "multiserver queue needs the level-dependent boundary of Geo/Geo/c.");
            }
        }
    }

    /** Discrete-time law of a station, expressed in slots. */
    private static MatrixCell stationLaw(NetworkStruct sn, int ist, int r, double slotLength) {
        Station station = sn.stations.get(ist);
        ProcessType procType = sn.procid.get(station).get(sn.jobclasses.get(r));
        double meanSlots = 1.0 / (sn.rates.get(ist, r) * slotLength);

        if (procType == ProcessType.DMAP) {
            if (slotLength != 1.0) {
                throw new RuntimeException("A DMAP is defined on its own slot, so it cannot be "
                        + "combined with config.slotlength=" + slotLength + ".");
            }
            MatrixCell proc = sn.proc.get(station).get(sn.jobclasses.get(r));
            return new MatrixCell(proc.get(0).copy(), proc.get(1).copy());
        }
        Pair<Matrix, Matrix> dph = Dph_from_dist.dph_from_dist(procType, meanSlots, sn.scv.get(ist, r));
        return Dph_to_dmap.dph_to_dmap(dph.getLeft(), dph.getRight());
    }

    /** Exact single-station analysis through the Q-MAM discrete-time queues. */
    private static DTQueueResult solveSingleStation(MatrixCell ARV, MatrixCell SVC, int maxNumComp) {
        MAPMAP1Options opts = new MAPMAP1Options("CR", maxNumComp, 0);
        boolean arvRenewal = Dph_to_dmap.dmap_is_renewal(ARV.get(0), ARV.get(1));
        boolean svcRenewal = Dph_to_dmap.dmap_is_renewal(SVC.get(0), SVC.get(1));
        if (arvRenewal && svcRenewal) {
            Matrix alpha = Dph_to_dmap.dmap_to_dph_alpha(ARV.get(0), ARV.get(1));
            Matrix beta = Dph_to_dmap.dmap_to_dph_alpha(SVC.get(0), SVC.get(1));
            return Q_DT_PH_PH_1.qDtPhPh1(alpha, ARV.get(0), beta, SVC.get(0), opts);
        }
        return Q_DT_MAP_MAP_1.qDtMapMap1(ARV.get(0), ARV.get(1), SVC.get(0), SVC.get(1), opts);
    }

    /**
     * Discrete-time parametric decomposition over several stations. Each station
     * is solved as a DBMAP/DMAP/1 queue given its arrival stream, and its
     * departure stream is extracted from the truncated stationary chain and
     * split by the routing probabilities. Superposing discrete streams produces
     * batches, which is why the station solve is M/G/1-type rather than a QBD.
     */
    private static int solveNetwork(NetworkStruct sn, MatrixCell[] law, int sourceIdx,
                                    List<Integer> queueIdx, SolverOptions options, int spaceMax,
                                    int maxNumComp, double[] QN, double[] UN, double[] TN) {
        int nq = queueIdx.size();
        Matrix P = routing(sn, sourceIdx, queueIdx);
        double lambda = Dmap_batch.dmap_lambda(law[sourceIdx]);

        int iterMax = options.iter_max > 0 ? options.iter_max : 100;
        double iterTol = options.iter_tol > 0 ? options.iter_tol : 1e-3;

        MatrixCell[] DEP = new MatrixCell[nq];
        for (int idx = 0; idx < nq; idx++) {
            double p = lambda * visitRatio(sn, queueIdx.get(idx).intValue());
            Matrix d0 = new Matrix(1, 1);
            d0.set(0, 0, 1 - p);
            Matrix d1 = new Matrix(1, 1);
            d1.set(0, 0, p);
            DEP[idx] = new MatrixCell(d0, d1);
        }

        double[] QNprev = new double[nq];
        double[] QNprev2 = new double[nq];
        double[] UNprev = new double[nq];
        double[] TNprev = new double[nq];
        int totiter = 0;

        for (int it = 1; it <= iterMax; it++) {
            totiter = it;
            for (int idx = 0; idx < nq; idx++) {
                MatrixCell ARV = arrivals(law[sourceIdx], DEP, P, idx, spaceMax);
                Mg1_dt_queue.Result r = Mg1_dt_queue.mg1_dt_queue(ARV,
                        law[queueIdx.get(idx).intValue()], maxNumComp, true);
                QN[idx] = r.QN;
                UN[idx] = r.UN;
                TN[idx] = r.TN;
                DEP[idx] = Dmap_batch.dmap_compress(r.dep, spaceMax);
            }
            if (it > 1 && maxRelChange(QN, QNprev) < iterTol) {
                break;
            }
            if (it > 2 && maxRelChange(QN, QNprev2) < iterTol) {
                // Feedback loops settle into a period-two cycle rather than a
                // point: re-solving a station with the departure process it just
                // produced moves it back. The cycle amplitude sits far below the
                // error of the decomposition itself, so the midpoint is reported
                // instead of burning iter_max sweeps on an orbit that will not close.
                for (int idx = 0; idx < nq; idx++) {
                    QN[idx] = (QN[idx] + QNprev[idx]) / 2;
                    UN[idx] = (UN[idx] + UNprev[idx]) / 2;
                    TN[idx] = (TN[idx] + TNprev[idx]) / 2;
                }
                break;
            }
            System.arraycopy(QNprev, 0, QNprev2, 0, nq);
            System.arraycopy(QN, 0, QNprev, 0, nq);
            System.arraycopy(UN, 0, UNprev, 0, nq);
            System.arraycopy(TN, 0, TNprev, 0, nq);
        }
        return totiter;
    }

    /** Arrival stream of a queue: source share plus thinned upstream departures. */
    private static MatrixCell arrivals(MatrixCell SRC, MatrixCell[] DEP, Matrix P, int idx, int spaceMax) {
        MatrixCell ARV = null;
        if (P.get(0, idx) > 0) {
            ARV = Dmap_batch.dmap_thin(SRC, P.get(0, idx));
        }
        for (int j = 0; j < DEP.length; j++) {
            double p = P.get(1 + j, idx);
            if (p <= 0) {
                continue;
            }
            MatrixCell contrib = Dmap_batch.dmap_thin(DEP[j], p);
            if (ARV == null) {
                ARV = contrib;
            } else {
                ARV = Dmap_batch.dmap_super(ARV, contrib);
                ARV = Dmap_batch.dmap_compress_batch(ARV, spaceMax);
            }
        }
        if (ARV == null) {
            throw new RuntimeException("Queue " + idx
                    + " receives no arrivals in the discrete-time routing matrix.");
        }
        return ARV;
    }

    /** Station-to-station routing probabilities, source first, sink stripped. */
    private static Matrix routing(NetworkStruct sn, int sourceIdx, List<Integer> queueIdx) {
        int I = sn.nnodes;
        int K = sn.nclasses;
        Matrix Pn = new Matrix(I, I);
        for (int a = 0; a < I; a++) {
            for (int b = 0; b < I; b++) {
                Pn.set(a, b, sn.rtnodes.get(a * K, b * K));
            }
        }
        // the Sink feeds back into the Source to keep rt stochastic; an open
        // traffic equation must not see that edge
        List<Integer> sinkNodes = new ArrayList<Integer>();
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) == NodeType.Sink) {
                sinkNodes.add(Integer.valueOf(ind));
                for (int b = 0; b < I; b++) {
                    Pn.set(ind, b, 0.0);
                }
            }
        }

        // linear indexing: NetworkStruct documents stationToNode as nstations x 1
        // but Network.refreshStruct builds it 1 x nstations, so (ist,0) is out of
        // bounds for every station past the first
        List<Integer> stationNodes = new ArrayList<Integer>();
        stationNodes.add(Integer.valueOf((int) sn.stationToNode.get(sourceIdx)));
        for (int idx = 0; idx < queueIdx.size(); idx++) {
            stationNodes.add(Integer.valueOf((int) sn.stationToNode.get(queueIdx.get(idx).intValue())));
        }
        List<Integer> interNodes = new ArrayList<Integer>();
        for (int ind = 0; ind < I; ind++) {
            if (!stationNodes.contains(Integer.valueOf(ind)) && !sinkNodes.contains(Integer.valueOf(ind))) {
                interNodes.add(Integer.valueOf(ind));
            }
        }

        int ns = stationNodes.size();
        Matrix Pss = sub(Pn, stationNodes, stationNodes);
        Matrix full;
        if (interNodes.isEmpty()) {
            full = Pss;
        } else {
            // censor the intermediate nodes: routers and class switches carry no
            // service, so their transit collapses into (I-Pnn)^-1
            Matrix Psn = sub(Pn, stationNodes, interNodes);
            Matrix Pnn = sub(Pn, interNodes, interNodes);
            Matrix Pns = sub(Pn, interNodes, stationNodes);
            Matrix inv = Matrix.eye(interNodes.size()).add(-1.0, Pnn).inv();
            full = Pss.add(1.0, Psn.mult(inv).mult(Pns));
        }
        // column 0 is the source, which receives nothing
        Matrix P = new Matrix(ns, ns - 1);
        for (int a = 0; a < ns; a++) {
            for (int b = 1; b < ns; b++) {
                P.set(a, b - 1, full.get(a, b));
            }
        }
        return P;
    }

    private static Matrix sub(Matrix src, List<Integer> rows, List<Integer> cols) {
        Matrix out = new Matrix(rows.size(), cols.size());
        for (int a = 0; a < rows.size(); a++) {
            for (int b = 0; b < cols.size(); b++) {
                out.set(a, b, src.get(rows.get(a).intValue(), cols.get(b).intValue()));
            }
        }
        return out;
    }

    private static double visitRatio(NetworkStruct sn, int ist) {
        double v = 0;
        for (Integer chain : sn.visits.keySet()) {
            Matrix vm = sn.visits.get(chain);
            if (vm != null && ist < vm.getNumRows()) {
                for (int r = 0; r < vm.getNumCols(); r++) {
                    v += vm.get(ist, r);
                }
            }
        }
        return v;
    }

    private static double maxRelChange(double[] a, double[] b) {
        double worst = 0;
        for (int i = 0; i < a.length; i++) {
            double denom = Math.max(Math.abs(b[i]), 1e-14);
            worst = Math.max(worst, Math.abs(a[i] - b[i]) / denom);
        }
        return worst;
    }

    private static int intConfig(SolverOptions options, String key, int fallback) {
        if (options != null && options.config != null) {
            Object v = options.config.get(key);
            if (v instanceof Number) {
                int val = ((Number) v).intValue();
                if (val > 0) {
                    return val;
                }
            }
        }
        return fallback;
    }
}
