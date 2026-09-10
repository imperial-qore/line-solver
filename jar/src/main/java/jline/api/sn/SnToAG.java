package jline.api.sn;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import jline.io.Ret;
import jline.lang.Event;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.constant.SignalType;
import jline.util.Utils;
import jline.util.matrix.Matrix;

/**
 * Convert LINE network structure to Agent (RCAT) format for SolverAG.
 * Converts a LINE queueing network to the R, AP format required by INAP and AUTOCAT algorithms.
 */
public final class SnToAG {
    private SnToAG() {}

    public static Ret.snToAG snToAG(NetworkStruct sn) {
        return snToAG(sn, 100);
    }

    public static Ret.snToAG snToAG(NetworkStruct sn, int maxStates) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Map<Integer, Sync> sync = sn.sync;
        int local = sn.nnodes + 1;

        List<Integer> sourceStations = new ArrayList<Integer>();
        List<Integer> sinkStations = new ArrayList<Integer>();
        List<Integer> queueStations = new ArrayList<Integer>();

        for (int ist = 0; ist < M; ist++) {
            int nodeIdx = (int) sn.stationToNode.get(ist, 0);
            NodeType nt = sn.nodetype.get(nodeIdx);
            if (nt == NodeType.Source) sourceStations.add(ist);
            else if (nt == NodeType.Sink) sinkStations.add(ist);
            else queueStations.add(ist);
        }

        int processIdx = 0;
        Matrix processMap = new Matrix(M, K);
        for (int ist : queueStations) {
            for (int r = 0; r < K; r++) {
                if (sn.issignal != null && sn.issignal.get(r, 0) != 0.0) continue;
                double rate = sn.rates.get(ist, r);
                if (!Double.isNaN(rate) && rate > 0) {
                    processIdx++;
                    processMap.set(ist, r, (double) processIdx);
                }
            }
        }
        int numProcesses = processIdx;

        if (numProcesses == 0) {
            return new Ret.snToAG(new Matrix[0][0], new Matrix(1, 2), processMap, new ArrayList<Ret.ActionMapEntry>(), new int[0]);
        }

        int[] N = new int[numProcesses];
        for (int p = 0; p < numProcesses; p++) {
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < K; r++) {
                    if ((int) processMap.get(ist, r) == p + 1) {
                        N[p] = Utils.isInf(sn.njobs.get(r)) ? maxStates : ((int) sn.njobs.get(r) + 1);
                    }
                }
            }
        }

        List<Ret.ActionMapEntry> actionMapList = new ArrayList<Ret.ActionMapEntry>();
        if (sync != null) {
            for (Integer s : sync.keySet()) {
                Sync syncEntry = sync.get(s);
                if (syncEntry == null) continue;
                if (syncEntry.active.isEmpty() || syncEntry.passive.isEmpty()) continue;
                Event activeEvent = syncEntry.active.get(0);
                Event passiveEvent = syncEntry.passive.get(0);
                if (activeEvent == null || passiveEvent == null) continue;
                if (activeEvent.getEvent() != EventType.DEP || passiveEvent.getEvent() != EventType.ARV) continue;
                int nodeA = activeEvent.getNode();
                int nodeP = passiveEvent.getNode();
                if (nodeP == local) continue;
                if (sn.isstation.get(nodeA, 0) == 0.0 || sn.isstation.get(nodeP, 0) == 0.0) continue;
                int ist = (int) sn.nodeToStation.get(nodeA, 0);
                int jst = (int) sn.nodeToStation.get(nodeP, 0);
                if (sourceStations.contains(ist) || sinkStations.contains(ist)) continue;
                if (sourceStations.contains(jst) || sinkStations.contains(jst)) continue;
                int r = activeEvent.getJobClass();
                int sClass = passiveEvent.getJobClass();
                if (processMap.get(ist, r) == 0.0 || processMap.get(jst, sClass) == 0.0) continue;
                double prob;
                if (passiveEvent.getProbFun() != null) prob = 1.0;
                else {
                    double rawProb = passiveEvent.getProb();
                    prob = Double.isNaN(rawProb) ? 1.0 : rawProb;
                }
                if (prob <= 0) continue;
                boolean isNegativeClass = false;
                if (sn.issignal != null && sn.issignal.get(r, 0) != 0.0) {
                    if (sn.signaltype != null && r < sn.signaltype.size()) {
                        if (sn.signaltype.get(r) == SignalType.NEGATIVE) isNegativeClass = true;
                    }
                }
                actionMapList.add(new Ret.ActionMapEntry(ist, r, jst, sClass, prob, isNegativeClass));
            }
        }
        int numActions = actionMapList.size();

        int numCols = (numProcesses < 2) ? 2 : numProcesses;
        Matrix[][] R = new Matrix[numActions + 1][numCols];
        Matrix AP = new Matrix(numActions < 1 ? 1 : numActions, 2);

        for (int p = 0; p < numProcesses; p++) {
            int pStation = -1;
            int pClass = -1;
            outerLoop:
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < K; r++) {
                    if ((int) processMap.get(ist, r) == p + 1) {
                        pStation = ist;
                        pClass = r;
                        break outerLoop;
                    }
                }
            }
            if (pStation >= 0 && pClass >= 0) {
                R[numActions][p] = buildLocalRates(sn, pStation, pClass, N[p], sourceStations, K);
            }
        }

        for (int a = 0; a < numActions; a++) {
            Ret.ActionMapEntry am = actionMapList.get(a);
            int ist = am.fromStation;
            int r = am.fromClass;
            int pActive = (int) processMap.get(ist, r) - 1;
            AP.set(a, 0, pActive + 1);

            double muIr = sn.rates.get(ist, r);
            double prob = am.prob;
            Matrix Aa = Matrix.zeros(N[pActive], N[pActive]);
            for (int n = 1; n < N[pActive]; n++) Aa.set(n, n - 1, muIr * prob);
            R[a][0] = Aa;

            int jst = am.toStation;
            int sClass = am.toClass;
            int pPassive = (int) processMap.get(jst, sClass) - 1;
            AP.set(a, 1, pPassive + 1);

            Matrix Pb = Matrix.zeros(N[pPassive], N[pPassive]);
            if (am.isNegative) {
                Pb.set(0, 0, 1.0);
                for (int n = 1; n < N[pPassive]; n++) Pb.set(n, n - 1, 1.0);
            } else {
                for (int n = 0; n < N[pPassive] - 1; n++) Pb.set(n, n + 1, 1.0);
                Pb.set(N[pPassive] - 1, N[pPassive] - 1, 1.0);
            }
            R[a][1] = Pb;
        }

        Matrix[][] RResult = new Matrix[R.length][];
        for (int i = 0; i < R.length; i++) {
            RResult[i] = new Matrix[R[i].length];
            for (int j = 0; j < R[i].length; j++) {
                RResult[i][j] = (R[i][j] != null) ? R[i][j] : Matrix.zeros(1, 1);
            }
        }
        return new Ret.snToAG(RResult, AP, processMap, actionMapList, N);
    }

    private static Matrix buildLocalRates(NetworkStruct sn, int ist, int r, int Np, List<Integer> sourceStations, int K) {
        Matrix L = Matrix.zeros(Np, Np);
        Matrix rt = sn.rt;

        double lambdaIr = 0.0;
        for (int isrc : sourceStations) {
            for (int sSrc = 0; sSrc < K; sSrc++) {
                int isfSrc = (int) sn.nodeToStateful.get((int) sn.stationToNode.get(isrc, 0), 0);
                int isfIst = (int) sn.nodeToStateful.get((int) sn.stationToNode.get(ist, 0), 0);
                double probSrc = rt.get(isfSrc * K + sSrc, isfIst * K + r);
                double srcRate = sn.rates.get(isrc, sSrc);
                if (probSrc > 0 && !Double.isNaN(srcRate)) {
                    boolean isNegativeSignal = false;
                    if (sn.issignal != null && sn.issignal.get(sSrc, 0) != 0.0) {
                        if (sn.signaltype != null && sSrc < sn.signaltype.size()) {
                            if (sn.signaltype.get(sSrc) == SignalType.NEGATIVE) isNegativeSignal = true;
                        }
                    }
                    if (!isNegativeSignal) lambdaIr += srcRate * probSrc;
                }
            }
        }
        if (lambdaIr > 0) {
            for (int n = 0; n < Np - 1; n++) L.set(n, n + 1, lambdaIr);
        }

        double muIr = sn.rates.get(ist, r);
        if (!Double.isNaN(muIr) && muIr > 0) {
            double probSink = 0.0;
            int isfIst = (int) sn.nodeToStateful.get((int) sn.stationToNode.get(ist, 0), 0);
            for (int jsrc : sourceStations) {
                int jsf = (int) sn.nodeToStateful.get((int) sn.stationToNode.get(jsrc, 0), 0);
                for (int s = 0; s < K; s++) probSink += rt.get(isfIst * K + r, jsf * K + s);
            }
            if (probSink > 0) {
                for (int n = 1; n < Np; n++) L.set(n, n - 1, L.get(n, n - 1) + muIr * probSink);
            }
        }
        return L;
    }

    /** Stochastic network ToAG algorithms. */
    public static final class SntoagKt {}
}
