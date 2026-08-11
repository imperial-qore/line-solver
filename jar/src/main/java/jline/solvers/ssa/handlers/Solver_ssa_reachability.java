package jline.solvers.ssa.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.Triple;

import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.state.State;
import jline.lang.nodes.StatefulNode;
import jline.lang.Sync;
import jline.solvers.SolverOptions;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;

/**
 * Kotlin migration of solver_ssa_reachability.m
 * Computes the reachable state space for SSA analysis
 */
public final class Solver_ssa_reachability {
    private Solver_ssa_reachability() {}

    public static Triple<Matrix, Matrix, NetworkStruct> solver_ssa_reachability(NetworkStruct sn, SolverOptions options) {
        int nstateful = sn.nstateful;
        int R = sn.nclasses;
        Matrix N = sn.njobs.transpose();
        Map<Integer, Sync> sync = sn.sync;
        Matrix csmask = sn.csmask;
        List<List<Matrix>> stack = new ArrayList<List<Matrix>>();
        List<Integer> stackIndex = new ArrayList<Integer>();

        // Initialize with starting state
        List<Matrix> initialStateCell = new ArrayList<Matrix>();
        for (int i = 0; i < nstateful; i++) {
            initialStateCell.add(sn.state.get(sn.stateful.get(i)));
        }
        stack.add(initialStateCell);

        Matrix SSq = null;
        int A = sync.size();
        boolean isSimulation = false;
        int local = sn.nnodes + 1;

        // Pre-compute sync node and class information
        int[] nodeA = new int[A];
        int[] nodeP = new int[A];
        int[] classA = new int[A];
        int[] classP = new int[A];
        EventType[] eventA = new EventType[A];
        EventType[] eventP = new EventType[A];
        for (int act = 0; act < A; act++) {
            nodeA[act] = sync.get(act).active.get(0).getNode();
            nodeP[act] = sync.get(act).passive.get(0).getNode();
            classA[act] = sync.get(act).active.get(0).getJobClass();
            classP[act] = sync.get(act).passive.get(0).getJobClass();
            eventA[act] = sync.get(act).active.get(0).getEvent();
            eventP[act] = sync.get(act).passive.get(0).getEvent();
        }

        List<Matrix> space = new ArrayList<Matrix>();
        for (int i = 0; i < nstateful; i++) {
            space.add(sn.state.get(sn.stateful.get(i)).copy());
        }

        Matrix SSh = new Matrix(1, nstateful);
        SSh.fill(1.0);

        int ih = 1;
        stackIndex.add(Integer.valueOf(1));
        int[] maxstatesz = new int[nstateful];

        while (!stack.isEmpty()) {
            if (stack.isEmpty()) {
                int totalCols = 0;
                for (Matrix m : space) totalCols += m.getNumCols();
                SSq = new Matrix(SSh.getNumRows(), totalCols);
                for (int i = 0; i < SSh.getNumRows(); i++) {
                    int colctr = 0;
                    for (int j = 0; j < nstateful; j++) {
                        int stateIdx = (int) SSh.get(i, j) - 1;
                        Matrix stateRow = space.get(j).getRow(stateIdx);
                        for (int k = 0; k < stateRow.getNumCols(); k++) {
                            SSq.set(i, colctr + k, stateRow.get(0, k));
                        }
                        colctr += space.get(j).getNumCols();
                    }
                }
                Map<StatefulNode, Matrix> spaceMap = new HashMap<StatefulNode, Matrix>();
                for (int i = 0; i < space.size(); i++) {
                    spaceMap.put(sn.stateful.get(i), space.get(i));
                }
                sn.space = spaceMap;
                State.buildSpaceHash(sn);
                return new Triple<Matrix, Matrix, NetworkStruct>(SSq, SSh, sn);
            }

            // Pop state from stack
            List<Matrix> stateCell = stack.remove(stack.size() - 1);
            ih = stackIndex.remove(stackIndex.size() - 1).intValue();

            List<List<Matrix>> newStateCell = new ArrayList<List<Matrix>>();
            for (int act = 0; act < A; act++) {
                newStateCell.add(new ArrayList<Matrix>(stateCell));
            }

            List<Integer> enabledSync = new ArrayList<Integer>();
            List<Double> enabledRates = new ArrayList<Double>();

            // Process each synchronization action
            for (int act = 0; act < A; act++) {
                boolean updateCondA = true;
                if (updateCondA) {
                    int isf = (int) sn.nodeToStateful.get(nodeA[act]);
                    Ret.EventResult activeResult = State.afterEvent(sn, nodeA[act], stateCell.get(isf), eventA[act], classA[act], isSimulation);

                    if (activeResult.outspace.isEmpty() || activeResult.outrate.isEmpty()) {
                        continue;
                    }

                    newStateCell.get(act).set((int) sn.nodeToStateful.get(nodeA[act]), activeResult.outspace);
                    Matrix rateA = activeResult.outrate;

                    for (int ia = 0; ia < activeResult.outspace.getNumRows(); ia++) {
                        if (activeResult.outspace.getRow(ia).elementSum() == -1.0) {
                            continue;
                        }

                        if (rateA.get(ia, 0) > 0) {
                            double probSyncP = 1.0;

                            if (nodeP[act] != local) {
                                if (nodeP[act] == nodeA[act]) {
                                    Ret.EventResult passiveResult = State.afterEvent(sn, nodeP[act],
                                            newStateCell.get(act).get((int) sn.nodeToStateful.get(nodeA[act])),
                                            eventP[act], classP[act], isSimulation);
                                    newStateCell.get(act).set((int) sn.nodeToStateful.get(nodeP[act]), passiveResult.outspace);
                                } else {
                                    Ret.EventResult passiveResult = State.afterEvent(sn, nodeP[act],
                                            stateCell.get((int) sn.nodeToStateful.get(nodeP[act])),
                                            eventP[act], classP[act], isSimulation);
                                    newStateCell.get(act).set((int) sn.nodeToStateful.get(nodeP[act]), passiveResult.outspace);
                                }

                                if (!newStateCell.get(act).get((int) sn.nodeToStateful.get(nodeP[act])).isEmpty()) {
                                    if (sn.isstatedep.get(nodeA[act], 2) == 1.0) {
                                        probSyncP = sync.get(act).passive.get(0).getProb();
                                    } else {
                                        probSyncP = sync.get(act).passive.get(0).getProb();
                                    }
                                } else {
                                    probSyncP = 0.0;
                                }
                            }

                            if (!newStateCell.get(act).get((int) sn.nodeToStateful.get(nodeA[act])).isEmpty()) {
                                if (nodeP[act] == local) {
                                    probSyncP = 1.0;
                                }

                                if (!Double.isNaN(rateA.get(ia, 0))) {
                                    boolean allNonEmpty = true;
                                    for (Matrix stateMatrix : newStateCell.get(act)) {
                                        if (stateMatrix.isEmpty()) {
                                            allNonEmpty = false;
                                            break;
                                        }
                                    }
                                    if (allNonEmpty) {
                                        if (nodeP[act] < local && csmask.get(classA[act], classP[act]) == 0.0
                                                && sn.nodetype.get(nodeP[act]) != NodeType.Source
                                                && (rateA.get(ia, 0) * probSyncP > 0)) {
                                            throw new RuntimeException("Error: state-dependent routing at node "
                                                    + nodeA[act] + " violates the class switching mask");
                                        }

                                        enabledRates.add(Double.valueOf(rateA.get(ia, 0) * probSyncP));
                                        enabledSync.add(Integer.valueOf(act));
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Process enabled transitions
            for (int firingCtr = 0; firingCtr < enabledRates.size(); firingCtr++) {
                double firingRate = enabledRates.get(firingCtr).doubleValue();
                int act = enabledSync.get(firingCtr).intValue();
                List<Matrix> netstates = newStateCell.get(act);

                boolean allNetStatesNonEmpty = true;
                for (Matrix stateMatrix : netstates) {
                    if (stateMatrix.isEmpty()) {
                        allNetStatesNonEmpty = false;
                        break;
                    }
                }
                if (firingRate > 0 && allNetStatesNonEmpty) {
                    Matrix nvec = new Matrix(1, netstates.size());
                    for (int i = 0; i < netstates.size(); i++) {
                        nvec.set(0, i, (double) (netstates.get(i).getNumRows() - 1));
                    }
                    Matrix n = PopulationLattice.pprod(nvec);

                    while (n != null && !n.isEmpty()) {
                        List<Matrix> newstatec = new ArrayList<Matrix>();

                        for (int i = 0; i < netstates.size(); i++) {
                            int nIdx = (int) n.get(0, i);
                            maxstatesz[i] = Math.max(maxstatesz[i], netstates.get(i).getRow(nIdx).getNumCols());
                            Matrix paddedState = new Matrix(1, maxstatesz[i]);
                            paddedState.fill(0.0);
                            Matrix originalState = netstates.get(i).getRow(nIdx);
                            for (int j = 0; j < originalState.getNumCols(); j++) {
                                paddedState.set(0, j + (maxstatesz[i] - originalState.getNumCols()), originalState.get(0, j));
                            }
                            newstatec.add(paddedState);
                        }

                        Matrix hashednewstate = new Matrix(1, nstateful);
                        for (int i = 0; i < nstateful; i++) {
                            int matchIdx = findMatchingRow(space.get(i), newstatec.get(i));
                            hashednewstate.set(0, i, (double) (matchIdx + 1));
                        }

                        int jh = findMatchingRow(SSh, hashednewstate);
                        if (jh == -1) {
                            for (int i = 0; i < nstateful; i++) {
                                if (hashednewstate.get(0, i) == 0.0) {
                                    Matrix newSpace = new Matrix(space.get(i).getNumRows() + 1, space.get(i).getNumCols());
                                    for (int r = 0; r < space.get(i).getNumRows(); r++) {
                                        for (int c = 0; c < space.get(i).getNumCols(); c++) {
                                            newSpace.set(r, c, space.get(i).get(r, c));
                                        }
                                    }
                                    for (int c = 0; c < newstatec.get(i).getNumCols(); c++) {
                                        newSpace.set(space.get(i).getNumRows(), c, newstatec.get(i).get(0, c));
                                    }
                                    space.set(i, newSpace);
                                    hashednewstate.set(0, i, (double) space.get(i).getNumRows());
                                }
                            }

                            Matrix newSSh = new Matrix(SSh.getNumRows() + 1, SSh.getNumCols());
                            for (int r = 0; r < SSh.getNumRows(); r++) {
                                for (int c = 0; c < SSh.getNumCols(); c++) {
                                    newSSh.set(r, c, SSh.get(r, c));
                                }
                            }
                            for (int c = 0; c < hashednewstate.getNumCols(); c++) {
                                newSSh.set(SSh.getNumRows(), c, hashednewstate.get(0, c));
                            }
                            stack.add(newstatec);
                            stackIndex.add(Integer.valueOf(newSSh.getNumRows()));
                        }
                        n = PopulationLattice.pprod(n, nvec);
                    }
                }
            }
        }

        return new Triple<Matrix, Matrix, NetworkStruct>(SSq, SSh, sn);
    }

    private static int findMatchingRow(Matrix matrix, Matrix targetRow) {
        if (targetRow.getNumRows() != 1) {
            throw new IllegalArgumentException("Target must be a single row");
        }

        for (int i = 0; i < matrix.getNumRows(); i++) {
            if (matrix.getNumCols() == targetRow.getNumCols()) {
                boolean match = true;
                for (int j = 0; j < matrix.getNumCols(); j++) {
                    if (matrix.get(i, j) != targetRow.get(0, j)) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    return i;
                }
            }
        }
        return -1;
    }
}
