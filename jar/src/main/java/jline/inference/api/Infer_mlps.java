/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import java.util.function.Function;

import jline.api.mam.Map_pdf;
import jline.inference.util.OptimUtils;
import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.RoutingMatrix;
import jline.lang.constant.EventType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.lang.Sync;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.ResultCTMC;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ctmc.handlers.Solver_ctmc;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Infer_mlps {
    private Infer_mlps() {}

    /**
     * Pre-built CTMC model information for a unique (tagClass, arrivalQueue) pair.
     */
    private static final class PrebuiltCTMC {
        final NetworkStruct sn;
        final int queueStIdx;
        final List<Integer> taggedDepIdx;
        final int[] subset;
        final Matrix SSqueue;
        final int[] N;
        final int tagClass;

        PrebuiltCTMC(NetworkStruct sn, int queueStIdx, List<Integer> taggedDepIdx,
                     int[] subset, Matrix SSqueue, int[] N, int tagClass) {
            this.sn = sn;
            this.queueStIdx = queueStIdx;
            this.taggedDepIdx = taggedDepIdx;
            this.subset = subset;
            this.SSqueue = SSqueue;
            this.N = N;
            this.tagClass = tagClass;
        }
    }

    private static final class CachedResult {
        final Matrix A;
        final Matrix SSqueue;
        final int[] N;

        CachedResult(Matrix A, Matrix SSqueue, int[] N) {
            this.A = A;
            this.SSqueue = SSqueue;
            this.N = N;
        }
    }

    public static double[] infer_mlps(Network model, Queue node, final double[] rt,
                                       final int[] classVec, final Matrix ql) {
        final NetworkStruct sn = model.getStruct(false);
        final int R = sn.nclasses;
        int nCores = node.getNumberOfServers();

        // Get delay rates from model
        int delayIdx = -1;
        for (int i = 0; i < sn.nstations; i++) {
            Station station = sn.stations.get(i);
            if (sn.sched.get(station) == SchedStrategy.INF) {
                delayIdx = i;
                break;
            }
        }
        final double[] muZ = new double[R];
        for (int k = 0; k < R; k++) {
            Station station = sn.stations.get(delayIdx);
            JobClass jobclass = sn.jobclasses.get(k);
            muZ[k] = sn.mu.get(station).get(jobclass).get(0, 0);
        }

        // Initial point estimate
        double meanQLSum = 0.0;
        int meanQLCnt = 0;
        for (int i = 0; i < ql.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < ql.getNumCols(); j++) {
                rowSum += ql.get(i, j);
            }
            meanQLSum += rowSum;
            meanQLCnt++;
        }
        double meanQL = meanQLSum / meanQLCnt;
        double Vtilde = Math.min(meanQL, (double) nCores);
        double[] x0 = new double[R];
        for (int j = 0; j < R; j++) {
            List<Double> classRT = new ArrayList<Double>();
            for (int idx2 = 0; idx2 < rt.length; idx2++) {
                if (classVec[idx2] == j) classRT.add(Double.valueOf(rt[idx2]));
            }
            if (!classRT.isEmpty()) {
                double avg = 0.0;
                for (Double d : classRT) avg += d.doubleValue();
                avg /= classRT.size();
                x0[j] = Vtilde * avg / meanQL;
            } else {
                x0[j] = 1e-3;
            }
        }

        double[] xLB = new double[R];
        double[] xUB = new double[R];
        double rtMax = rt[0];
        for (double v : rt) if (v > rtMax) rtMax = v;
        for (int i = 0; i < R; i++) {
            xLB[i] = 1e-10;
            xUB[i] = rtMax;
        }

        final int newR = R + 1;
        final double[] augMuZ = new double[newR];
        for (int k = 0; k < R; k++) augMuZ[k] = muZ[k];

        // Collect unique (tagClass, arrivalQueue) pairs
        java.util.Set<Integer> tcSet = new java.util.TreeSet<Integer>();
        for (int v : classVec) tcSet.add(Integer.valueOf(v));
        List<Integer> uniqueTC = new ArrayList<Integer>(tcSet);
        Collections.sort(uniqueTC);

        LinkedHashSet<String> uniqueQLRows = new LinkedHashSet<String>();
        for (int i = 0; i < ql.getNumRows(); i++) {
            int[] row = new int[R];
            for (int j = 0; j < R; j++) row[j] = (int) ql.get(i, j);
            uniqueQLRows.add(intArrayToCsv(row));
        }

        SolverOptions ctmcOpts = null;
        final Map<String, PrebuiltCTMC> prebuilt = new HashMap<String, PrebuiltCTMC>();

        for (Integer tcBox : uniqueTC) {
            int tc = tcBox.intValue();
            augMuZ[newR - 1] = muZ[tc];
            for (String qlRowStr : uniqueQLRows) {
                String[] parts = qlRowStr.split(",");
                int[] aq = new int[parts.length];
                for (int i = 0; i < parts.length; i++) aq[i] = Integer.parseInt(parts[i]);

                // Population: move 1 job from tagClass to tagged class
                int[] N = new int[newR];
                for (int k = 0; k < R; k++) N[k] = aq[k];
                N[tc] = N[tc] - 1;
                N[newR - 1] = 1;

                if (N[tc] < 0) continue;

                Network augModel = new Network("mlps_aug");
                Delay augDelay = new Delay(augModel, "Think");
                Queue augQueue = new Queue(augModel, "Queue1", SchedStrategy.PS);
                augQueue.setNumberOfServers(nCores);

                ClosedClass[] augClasses = new ClosedClass[newR];
                for (int r = 0; r < newR; r++) {
                    augClasses[r] = new ClosedClass(augModel, "Class" + (r + 1), N[r], augDelay, 0);
                    augDelay.setService(augClasses[r], new Exp(augMuZ[r]));
                    augQueue.setService(augClasses[r], new Exp(1.0));
                }

                RoutingMatrix P = augModel.initRoutingMatrix();
                List<Node> nodeList = new ArrayList<Node>();
                nodeList.add(augDelay);
                nodeList.add(augQueue);
                for (int r = 0; r < newR; r++) {
                    P.set(augClasses[r], Network.serialRouting(nodeList));
                }
                augModel.link(P);

                SolverCTMC solver = new SolverCTMC(augModel);
                if (ctmcOpts == null) ctmcOpts = solver.getOptions();

                NetworkStruct augSn = augModel.getStruct(false);

                ResultCTMC ctmcResult = Solver_ctmc.solver_ctmc(augSn, ctmcOpts);
                Matrix stateSpaceAggr = ctmcResult.getStateSpaceAggr();

                int queueStIdx = augQueue.getStationIdx();
                int queueNodeIdx = augQueue.getNodeIndex();
                List<Integer> taggedDepIdx = new ArrayList<Integer>();
                Map<Integer, Sync> sync = augSn.sync;
                if (sync != null) {
                    for (int e = 0; e < sync.size(); e++) {
                        Sync syncEntry = sync.get(Integer.valueOf(e));
                        if (syncEntry == null) continue;
                        jline.lang.Event activeEvent = syncEntry.active.get(Integer.valueOf(0));
                        if (activeEvent == null) continue;
                        if (activeEvent.getNode() == queueNodeIdx
                                && activeEvent.getJobClass() == newR - 1
                                && activeEvent.getEvent() == EventType.DEP) {
                            taggedDepIdx.add(Integer.valueOf(e));
                        }
                    }
                }

                int taggedColAtQueue = queueStIdx * newR + (newR - 1);
                List<Integer> subsetList = new ArrayList<Integer>();
                if (stateSpaceAggr != null) {
                    for (int row = 0; row < stateSpaceAggr.getNumRows(); row++) {
                        if (taggedColAtQueue < stateSpaceAggr.getNumCols()
                                && (int) stateSpaceAggr.get(row, taggedColAtQueue) == 1) {
                            subsetList.add(Integer.valueOf(row));
                        }
                    }
                }
                int[] subset = new int[subsetList.size()];
                for (int i = 0; i < subset.length; i++) subset[i] = subsetList.get(i).intValue();

                int[] queueCols = new int[newR];
                for (int r = 0; r < newR; r++) queueCols[r] = queueStIdx * newR + r;
                Matrix SSqueue = new Matrix(subset.length, newR);
                for (int si = 0; si < subset.length; si++) {
                    for (int r = 0; r < newR; r++) {
                        if (queueCols[r] < stateSpaceAggr.getNumCols()) {
                            SSqueue.set(si, r, stateSpaceAggr.get(subset[si], queueCols[r]));
                        }
                    }
                }

                String cacheKey = tc + "," + intArrayToCsv(aq);
                prebuilt.put(cacheKey, new PrebuiltCTMC(augSn, queueStIdx, taggedDepIdx, subset, SSqueue, N, tc));
            }
        }

        final SolverOptions finalCtmcOpts = (ctmcOpts != null) ? ctmcOpts : SolverCTMC.defaultOptions();

        Function<double[], Double> objFun = new Function<double[], Double>() {
            @Override
            public Double apply(double[] x) {
                final double TOL = 1e-6;
                double[] rates = new double[R];
                for (int i = 0; i < R; i++) rates[i] = 1.0 / x[i];

                Map<String, CachedResult> cache = new HashMap<String, CachedResult>();

                for (Map.Entry<String, PrebuiltCTMC> entry : prebuilt.entrySet()) {
                    String key = entry.getKey();
                    PrebuiltCTMC pb = entry.getValue();

                    double[] augRates = new double[newR];
                    for (int k = 0; k < R; k++) augRates[k] = rates[k];
                    augRates[newR - 1] = rates[pb.tagClass];

                    NetworkStruct snUpd = pb.sn;
                    for (int cc = 0; cc < newR; cc++) {
                        Sn_set_service_coc.sn_set_service_coc(snUpd, pb.queueStIdx, cc, augRates[cc]);
                    }

                    ResultCTMC ctmcResult = Solver_ctmc.solver_ctmc(snUpd, finalCtmcOpts);
                    Matrix infGen = ctmcResult.getQ();

                    MatrixCell Dfilt = ctmcResult.getDfilt();
                    int nStates = infGen.getNumRows();
                    Matrix D1 = new Matrix(nStates, nStates);
                    if (Dfilt != null) {
                        for (Integer diBox : pb.taggedDepIdx) {
                            int di = diBox.intValue();
                            if (di < Dfilt.size()) {
                                Matrix dfiltMat = Dfilt.get(di);
                                if (dfiltMat == null) continue;
                                int rr = Math.min(nStates, dfiltMat.getNumRows());
                                int cc = Math.min(nStates, dfiltMat.getNumCols());
                                for (int ii = 0; ii < rr; ii++) {
                                    for (int jj = 0; jj < cc; jj++) {
                                        D1.set(ii, jj, D1.get(ii, jj) + dfiltMat.get(ii, jj));
                                    }
                                }
                            }
                        }
                    }

                    int subSize = pb.subset.length;
                    if (subSize > 0) {
                        Matrix MAPQ1 = new Matrix(nStates, nStates);
                        for (int ii = 0; ii < nStates; ii++) {
                            for (int jj = 0; jj < nStates; jj++) {
                                MAPQ1.set(ii, jj, infGen.get(ii, jj) - D1.get(ii, jj));
                            }
                        }
                        Matrix A = new Matrix(subSize, subSize);
                        for (int ii = 0; ii < subSize; ii++) {
                            for (int jj = 0; jj < subSize; jj++) {
                                A.set(ii, jj, MAPQ1.get(pb.subset[ii], pb.subset[jj]));
                            }
                        }
                        cache.put(key, new CachedResult(A, pb.SSqueue, pb.N));
                    }
                }

                double totalLogLike = 0.0;
                for (int i = 0; i < rt.length; i++) {
                    int tc = classVec[i];
                    int[] aqRow = new int[R];
                    for (int j = 0; j < R; j++) aqRow[j] = (int) ql.get(i, j);
                    String cacheKey = tc + "," + intArrayToCsv(aqRow);
                    CachedResult cached = cache.get(cacheKey);
                    if (cached != null) {
                        double like = evalMlpsLikelihood(cached.A, cached.SSqueue, cached.N, rt[i]);
                        totalLogLike += Math.log(TOL + like);
                    } else {
                        totalLogLike += Math.log(TOL);
                    }
                }
                return Double.valueOf(-totalLogLike);
            }
        };

        Pair<double[], Double> opt = OptimUtils.fmincon(objFun, x0, xLB, xUB, 10000);
        return opt.getLeft();
    }

    private static String intArrayToCsv(int[] arr) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < arr.length; i++) {
            if (i > 0) sb.append(',');
            sb.append(arr[i]);
        }
        return sb.toString();
    }

    private static double evalMlpsLikelihood(Matrix A, Matrix SSqueue, int[] N, double Rsam) {
        int nStates = A.getNumRows();
        if (nStates == 0) return 0.0;

        Matrix pie = new Matrix(1, nStates);
        for (int row = 0; row < SSqueue.getNumRows(); row++) {
            boolean matches = true;
            int upTo = Math.min(SSqueue.getNumCols(), N.length);
            for (int col = 0; col < upTo; col++) {
                if ((int) SSqueue.get(row, col) != N[col]) {
                    matches = false;
                    break;
                }
            }
            if (matches) {
                pie.set(0, row, 1.0);
                break;
            }
        }

        Matrix D0 = A;
        Matrix rowSums = new Matrix(nStates, 1);
        for (int i = 0; i < nStates; i++) {
            double sum = 0.0;
            for (int j = 0; j < nStates; j++) {
                sum += A.get(i, j);
            }
            rowSums.set(i, 0, -sum);
        }

        Matrix D1 = new Matrix(nStates, nStates);
        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
                D1.set(i, j, rowSums.get(i, 0) * pie.get(0, j));
            }
        }

        MatrixCell MAP = new MatrixCell(2);
        MAP.set(0, D0);
        MAP.set(1, D1);
        return Map_pdf.map_pdf(MAP, Rsam);
    }
}
