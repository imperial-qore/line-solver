package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SignalType;
import jline.lang.processes.DiscreteDistribution;
import jline.util.matrix.Matrix;

public final class Solver_mam_build_ag {
    private Solver_mam_build_ag() {}

    public static RCATModel solver_mam_build_ag(NetworkStruct sn, int maxStates) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix rt = sn.rt;

        List<Integer> sourceStations = new ArrayList<Integer>();
        List<Integer> queueStations = new ArrayList<Integer>();
        List<Integer> sinkNodes = new ArrayList<Integer>();
        for (int nodeIdx = 0; nodeIdx < sn.nodetype.size(); nodeIdx++) {
            if (sn.nodetype.get(nodeIdx) == NodeType.Sink) sinkNodes.add(nodeIdx);
        }
        for (int ist = 0; ist < M; ist++) {
            int nodeIdx = (int) sn.stationToNode.get(ist);
            NodeType nt = sn.nodetype.get(nodeIdx);
            if (nt == NodeType.Source) sourceStations.add(ist);
            else queueStations.add(ist);
        }

        int processIdx = 0;
        Matrix processMap = new Matrix(M, K);
        processMap.fill(-1.0);
        for (int ist : queueStations) {
            for (int r = 0; r < K; r++) {
                if (sn.issignal != null && sn.issignal.get(r, 0) > 0) continue;
                double rate = sn.rates.get(ist, r);
                if (!Double.isNaN(rate) && rate > 0) {
                    processMap.set(ist, r, processIdx);
                    processIdx++;
                }
            }
        }
        int numProcesses = processIdx;
        if (numProcesses == 0) {
            return new RCATModel(new Matrix[1][1], new Matrix(0, 2), processMap, new ArrayList<ActionInfo>(), new int[0]);
        }

        int[] N = new int[numProcesses];
        for (int p = 0; p < numProcesses; p++) {
            outerLoop:
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < K; r++) {
                    if ((int) processMap.get(ist, r) == p) {
                        double njobs = sn.njobs.get(r);
                        N[p] = Double.isInfinite(njobs) ? maxStates : ((int) njobs + 1);
                        break outerLoop;
                    }
                }
            }
        }

        List<ActionInfo> actionMap = new ArrayList<ActionInfo>();
        for (int ist : queueStations) {
            for (int r = 0; r < K; r++) {
                if ((int) processMap.get(ist, r) >= 0) {
                    // Removal signal: NEGATIVE or CATASTROPHE (the two are
                    // distinct SignalType values, so both must be tested).
                    boolean isNegativeClass = sn.issignal != null
                            && sn.issignal.get(r, 0) > 0
                            && sn.signaltype != null
                            && (sn.signaltype.get(r) == SignalType.NEGATIVE
                                || sn.signaltype.get(r) == SignalType.CATASTROPHE);
                    boolean isCatastropheClass = (sn.iscatastrophe != null
                            && sn.iscatastrophe.get(r, 0) > 0)
                            || (sn.signaltype != null && sn.signaltype.get(r) == SignalType.CATASTROPHE);
                    DiscreteDistribution removalDist = null;
                    if (sn.signalremdist != null && r < sn.signalremdist.size()) {
                        removalDist = sn.signalremdist.get(r);
                    }
                    for (int jst : queueStations) {
                        for (int s = 0; s < K; s++) {
                            if ((int) processMap.get(jst, s) >= 0) {
                                double probIJRS = rt.get(ist * K + r, jst * K + s);
                                if (probIJRS > 0 && (ist != jst || r != s)) {
                                    actionMap.add(new ActionInfo(ist, r, jst, s, probIJRS, isNegativeClass, isCatastropheClass, removalDist));
                                }
                            }
                        }
                    }
                }
            }
        }
        int numActions = actionMap.size();

        Matrix[][] R = new Matrix[numActions + 1][Math.max(numProcesses, 2)];
        Matrix AP = new Matrix(Math.max(numActions, 1), 2);

        for (int p = 0; p < numProcesses; p++) {
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < K; r++) {
                    if ((int) processMap.get(ist, r) == p) {
                        R[numActions][p] = buildLocalRates(sn, ist, r, N[p], rt, sourceStations, sinkNodes, K);
                        break;
                    }
                }
                if (R[numActions][p] != null) break;
            }
        }

        for (int a = 0; a < numActions; a++) {
            ActionInfo am = actionMap.get(a);
            int ist = am.fromStation;
            int r = am.fromClass;
            int pActive = (int) processMap.get(ist, r);
            AP.set(a, 0, (double) pActive);

            double muIr = sn.rates.get(ist, r);
            double prob = am.prob;

            Matrix Aa = new Matrix(N[pActive], N[pActive]);
            for (int n = 1; n < N[pActive]; n++) Aa.set(n, n - 1, muIr * prob);
            // see _kb/06-solver-catalog.md for rationale
            if (N[pActive] > 0 && !Double.isInfinite(sn.njobs.get(r))) {
                Aa.set(N[pActive] - 1, N[pActive] - 1, muIr * prob);
            }
            R[a][0] = Aa;

            int jst = am.toStation;
            int s = am.toClass;
            int pPassive = (int) processMap.get(jst, s);
            AP.set(a, 1, (double) pPassive);

            Matrix Pb = new Matrix(N[pPassive], N[pPassive]);
            if (am.isNegative) {
                if (am.isCatastrophe) {
                    for (int n = 0; n < N[pPassive]; n++) Pb.set(n, 0, 1.0);
                } else if (am.removalDistribution != null) {
                    DiscreteDistribution dist = am.removalDistribution;
                    for (int n = 0; n < N[pPassive]; n++) {
                        if (n == 0) {
                            Pb.set(0, 0, 1.0);
                        } else {
                            for (int mIdx = 0; mIdx <= n; mIdx++) {
                                int kRem = n - mIdx;
                                if (mIdx > 0) {
                                    double prRem = dist.evalPMF((double) kRem);
                                    if (prRem > 0) Pb.set(n, mIdx, Pb.get(n, mIdx) + prRem);
                                } else {
                                    double cdfNMinus1 = 0.0;
                                    for (int j = 0; j < n; j++) cdfNMinus1 += dist.evalPMF((double) j);
                                    double probAtLeastN = 1.0 - cdfNMinus1;
                                    if (probAtLeastN > 0) Pb.set(n, 0, Pb.get(n, 0) + probAtLeastN);
                                }
                            }
                        }
                    }
                } else {
                    Pb.set(0, 0, 1.0);
                    for (int n = 1; n < N[pPassive] - 1; n++) Pb.set(n, n - 1, 1.0);
                    if (N[pPassive] > 1) Pb.set(N[pPassive] - 1, N[pPassive] - 2, 1.0);
                }
            } else {
                for (int n = 0; n < N[pPassive] - 1; n++) Pb.set(n, n + 1, 1.0);
                if (N[pPassive] > 0) Pb.set(N[pPassive] - 1, N[pPassive] - 1, 1.0);
            }
            R[a][1] = Pb;
        }
        return new RCATModel(R, AP, processMap, actionMap, N);
    }

    private static Matrix buildLocalRates(NetworkStruct sn, int ist, int r, int Np, Matrix rt,
                                          List<Integer> sourceStations, List<Integer> sinkNodes, int K) {
        Matrix L = new Matrix(Np, Np);
        double lambdaIrPos = 0.0;
        double lambdaIrNegSingle = 0.0;
        double lambdaIrCatastrophe = 0.0;
        List<double[]> batchArrivals = new ArrayList<double[]>();
        List<DiscreteDistribution> batchDists = new ArrayList<DiscreteDistribution>();

        for (int isrc : sourceStations) {
            for (int sSrc = 0; sSrc < K; sSrc++) {
                boolean isSignalSrc = sn.issignal != null && sn.issignal.get(sSrc, 0) > 0;
                double probSrc;
                if (isSignalSrc) {
                    double pSum = 0.0;
                    for (int sDst = 0; sDst < K; sDst++) pSum += rt.get(isrc * K + sSrc, ist * K + sDst);
                    probSrc = pSum;
                } else {
                    probSrc = rt.get(isrc * K + sSrc, ist * K + r);
                }
                double srcRate = sn.rates.get(isrc, sSrc);
                if (probSrc > 0 && !Double.isNaN(srcRate)) {
                    // see _kb/06-solver-catalog.md for rationale
                    boolean isRemovalSrc = sn.issignal != null
                            && sn.issignal.get(sSrc, 0) > 0
                            && sn.signaltype != null
                            && (sn.signaltype.get(sSrc) == SignalType.NEGATIVE
                                || sn.signaltype.get(sSrc) == SignalType.CATASTROPHE);
                    if (isRemovalSrc) {
                        boolean isCat = (sn.iscatastrophe != null && sn.iscatastrophe.get(sSrc, 0) > 0)
                                || sn.signaltype.get(sSrc) == SignalType.CATASTROPHE;
                        if (isCat) {
                            lambdaIrCatastrophe += srcRate * probSrc;
                        } else {
                            DiscreteDistribution removalDist = null;
                            if (sn.signalremdist != null && sSrc < sn.signalremdist.size()) {
                                removalDist = sn.signalremdist.get(sSrc);
                            }
                            if (removalDist != null) {
                                batchArrivals.add(new double[]{srcRate * probSrc});
                                batchDists.add(removalDist);
                            } else {
                                lambdaIrNegSingle += srcRate * probSrc;
                            }
                        }
                    } else {
                        lambdaIrPos += srcRate * probSrc;
                    }
                }
            }
        }
        if (lambdaIrPos > 0) {
            for (int n = 0; n < Np - 1; n++) L.set(n, n + 1, lambdaIrPos);
        }
        if (lambdaIrCatastrophe > 0) {
            for (int n = 1; n < Np; n++) L.set(n, 0, L.get(n, 0) + lambdaIrCatastrophe);
        }
        for (int b = 0; b < batchArrivals.size(); b++) {
            double rate = batchArrivals.get(b)[0];
            DiscreteDistribution dist = batchDists.get(b);
            for (int n = 1; n < Np; n++) {
                for (int mIdx = 0; mIdx < n; mIdx++) {
                    int kRem = n - mIdx;
                    double pr;
                    if (mIdx > 0) pr = dist.evalPMF((double) kRem);
                    else {
                        double cdfNMinus1 = 0.0;
                        for (int j = 0; j < n; j++) cdfNMinus1 += dist.evalPMF((double) j);
                        pr = 1.0 - cdfNMinus1;
                    }
                    if (pr > 0) L.set(n, mIdx, L.get(n, mIdx) + rate * pr);
                }
            }
        }
        if (lambdaIrNegSingle > 0) {
            for (int n = 1; n < Np; n++) L.set(n, n - 1, L.get(n, n - 1) + lambdaIrNegSingle);
        }
        double muIr = sn.rates.get(ist, r);
        if (!Double.isNaN(muIr) && muIr > 0) {
            int nodeIdx = (int) sn.stationToNode.get(ist);
            double probSink = 0.0;
            if (sn.rtnodes != null && sn.rtnodes.getNumRows() > 0) {
                for (int jsnk : sinkNodes) {
                    for (int s = 0; s < K; s++) {
                        int fromIdx = nodeIdx * K + r;
                        int toIdx = jsnk * K + s;
                        if (fromIdx < sn.rtnodes.getNumRows() && toIdx < sn.rtnodes.getNumCols()) {
                            probSink += sn.rtnodes.get(fromIdx, toIdx);
                        }
                    }
                }
            }
            double probSelf = rt.get(ist * K + r, ist * K + r);
            double localDepartureRate = muIr * probSink;
            for (int n = 1; n < Np; n++) L.set(n, n - 1, L.get(n, n - 1) + localDepartureRate);
            for (int n = 1; n < Np; n++) L.set(n, n, muIr * probSelf);
        }
        return L;
    }
}
