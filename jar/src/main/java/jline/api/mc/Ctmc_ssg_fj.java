/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.mc;

import jline.io.Ret;
import jline.lang.FJSync;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.EventType;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.AfterFJEvent;
import jline.lang.state.EventCache;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Reachability-based CTMC state space generation for FJ tag-augmented
 * structs (ModelAdapter.fjtag). Explores the global state space forward
 * from the initial state through the regular synchronizations (sn.sync)
 * and the fork firing synchronizations (sn.fjsync). Required for
 * fork-join models, whose fork firings break per-chain population
 * conservation so the population-lattice enumeration of ctmc_ssg cannot
 * be used.
 *
 * Mirrors matlab/src/lang/+State/reachableSpaceGenerator.m plus its
 * fjsync exploration block. Hash values in the returned hashed state
 * space are 0-based row indices into sn.space, consistent with
 * Solver_ctmc.
 */
public final class Ctmc_ssg_fj {
    private Ctmc_ssg_fj() {}

    private static Matrix padLeft(Matrix row, int width) {
        if (row.getNumCols() >= width) {
            return row;
        }
        Matrix out = new Matrix(1, width);
        out.zero();
        int shift = width - row.getNumCols();
        for (int c = 0; c < row.getNumCols(); c++) {
            out.set(0, shift + c, row.get(0, c));
        }
        return out;
    }

    public static CtmcSsgReachabilityResult ctmc_ssg_fj(NetworkStruct sn, jline.solvers.SolverOptions options) {
        int nstateful = sn.nstateful;
        Map<Integer, Sync> sync = sn.sync;
        int A = sync.size();
        int local = sn.nnodes; // 0-based: passive node == nnodes means local action
        EventCache eventCache = new EventCache(false, false);

        // initial raw per-node state rows
        List<Matrix> init = new ArrayList<Matrix>();
        for (int isf = 0; isf < nstateful; isf++) {
            StatefulNode statefulNode = sn.stateful.get(isf);
            Matrix st = sn.state.get(statefulNode);
            if (st.getNumRows() > 1 && st.getNumCols() == 1) {
                st = st.transpose();
            }
            if (st.getNumRows() > 1) {
                st = Matrix.extractRows(st, 0, 1, null);
            }
            init.add(st.copy());
        }

        // per-node spaces, growing during exploration
        List<Matrix> space = new ArrayList<Matrix>();
        for (int isf = 0; isf < nstateful; isf++) {
            space.add(init.get(isf).copy());
        }

        List<int[]> SSh = new ArrayList<int[]>();
        Map<String, Integer> seen = new HashMap<String, Integer>();
        List<List<Matrix>> stack = new ArrayList<List<Matrix>>();
        int[] h0 = new int[nstateful];
        SSh.add(h0);
        seen.put(hashKey(h0), 0);
        stack.add(init);

        while (!stack.isEmpty()) {
            List<Matrix> stateCell = stack.remove(stack.size() - 1);

            // regular sync actions
            for (int act = 0; act < A; act++) {
                Sync syncA = sync.get(act);
                int nodeA = syncA.active.get(0).getNode();
                int classA = syncA.active.get(0).getJobClass();
                EventType eventA = syncA.active.get(0).getEvent();
                if (sn.isstateful.get(nodeA, 0) == 0) {
                    continue;
                }
                int isfA = (int) sn.nodeToStateful.get(nodeA);
                Ret.EventResult resA = State.afterEvent(sn, nodeA, stateCell.get(isfA), eventA, classA, false, eventCache);
                if (resA.outspace == null || resA.outspace.isEmpty()) {
                    continue;
                }
                for (int ia = 0; ia < resA.outspace.getNumRows(); ia++) {
                    double rateIa = (ia < resA.outrate.getNumRows()) ? resA.outrate.get(ia, 0) : resA.outrate.get(0, 0);
                    if (Double.isNaN(rateIa) || rateIa <= 0) {
                        continue;
                    }
                    List<Matrix> newCell = copyCell(stateCell);
                    newCell.set(isfA, Matrix.extractRows(resA.outspace, ia, ia + 1, null));

                    int nodeP = syncA.passive.get(0).getNode();
                    if (nodeP < sn.nnodes && sn.isstateful.get(nodeP, 0) == 1) {
                        int isfP = (int) sn.nodeToStateful.get(nodeP);
                        int classP = syncA.passive.get(0).getJobClass();
                        EventType eventP = syncA.passive.get(0).getEvent();
                        Matrix baseP = (nodeP == nodeA) ? newCell.get(isfP) : stateCell.get(isfP);
                        Ret.EventResult resP = State.afterEvent(sn, nodeP, baseP, eventP, classP, false, eventCache);
                        if (resP.outspace == null || resP.outspace.isEmpty()) {
                            continue;
                        }
                        for (int ip = 0; ip < resP.outspace.getNumRows(); ip++) {
                            List<Matrix> nc2 = copyCell(newCell);
                            nc2.set(isfP, Matrix.extractRows(resP.outspace, ip, ip + 1, null));
                            pushState(sn, nc2, space, SSh, seen, stack);
                        }
                    } else {
                        pushState(sn, newCell, space, SSh, seen, stack);
                    }
                }
            }

            // fork firing synchronizations
            if (sn.fjsync != null && !sn.fjsync.isEmpty()) {
                for (FJSync entry : sn.fjsync.values()) {
                    AfterFJEvent.AfterFJEventResult fjRes = AfterFJEvent.afterFJEvent(sn, entry, stateCell, false, eventCache);
                    for (int io = 0; io < fjRes.outGlobalStates.size(); io++) {
                        if (fjRes.outprob.get(io, 0) > 0) {
                            pushState(sn, fjRes.outGlobalStates.get(io), space, SSh, seen, stack);
                        }
                    }
                }
            }
        }

        // assemble outputs
        int nStates = SSh.size();
        Matrix stateSpaceHashed = new Matrix(nStates, nstateful);
        for (int s = 0; s < nStates; s++) {
            for (int isf = 0; isf < nstateful; isf++) {
                stateSpaceHashed.set(s, isf, SSh.get(s)[isf]);
            }
        }
        int totalCols = 0;
        for (int isf = 0; isf < nstateful; isf++) {
            totalCols += space.get(isf).getNumCols();
        }
        Matrix stateSpace = new Matrix(nStates, totalCols);
        for (int s = 0; s < nStates; s++) {
            int col = 0;
            for (int isf = 0; isf < nstateful; isf++) {
                Matrix sp = space.get(isf);
                for (int c = 0; c < sp.getNumCols(); c++) {
                    stateSpace.set(s, col + c, sp.get(SSh.get(s)[isf], c));
                }
                col += sp.getNumCols();
            }
        }

        Map<StatefulNode, Matrix> nodeStateSpace = new HashMap<StatefulNode, Matrix>();
        for (int isf = 0; isf < nstateful; isf++) {
            nodeStateSpace.put(sn.stateful.get(isf), space.get(isf));
        }
        sn.space = nodeStateSpace;

        // aggregated state space (per-station per-class counts)
        int nclasses = sn.nclasses;
        Matrix stateSpaceAggr = new Matrix(nStates, sn.nstations * nclasses);
        stateSpaceAggr.zero();
        for (int s = 0; s < nStates; s++) {
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.isstateful.get(ind, 0) == 1 && sn.isstation.get(ind, 0) == 1) {
                    int isf = (int) sn.nodeToStateful.get(ind);
                    int ist = (int) sn.nodeToStation.get(ind);
                    Matrix row = space.get(isf).getRow(SSh.get(s)[isf]);
                    State.StateMarginalStatistics stats = ToMarginal.toMarginal(sn, ind, row, null, null, null, null, null);
                    for (int c = 0; c < Math.min(nclasses, (int) stats.nir.length()); c++) {
                        stateSpaceAggr.set(s, ist * nclasses + c, stats.nir.get(c));
                    }
                }
            }
        }

        return new CtmcSsgReachabilityResult(stateSpace, stateSpaceAggr, stateSpaceHashed, nodeStateSpace, sn);
    }

    private static List<Matrix> copyCell(List<Matrix> cell) {
        List<Matrix> out = new ArrayList<Matrix>();
        for (Matrix m : cell) {
            out.add(m.copy());
        }
        return out;
    }

    private static String hashKey(int[] h) {
        StringBuilder sb = new StringBuilder();
        for (int v : h) {
            sb.append(v).append(',');
        }
        return sb.toString();
    }

    private static void pushState(NetworkStruct sn, List<Matrix> stateCell, List<Matrix> space,
                                  List<int[]> SSh, Map<String, Integer> seen, List<List<Matrix>> stack) {
        int nstateful = space.size();
        int[] hashed = new int[nstateful];
        for (int isf = 0; isf < nstateful; isf++) {
            Matrix row = stateCell.get(isf);
            Matrix sp = space.get(isf);
            if (row.getNumCols() > sp.getNumCols()) {
                // widen the space (left-zero padding, e.g. FCFS buffer growth)
                Matrix grown = new Matrix(sp.getNumRows(), row.getNumCols());
                grown.zero();
                int shift = row.getNumCols() - sp.getNumCols();
                for (int rr = 0; rr < sp.getNumRows(); rr++) {
                    for (int c = 0; c < sp.getNumCols(); c++) {
                        grown.set(rr, shift + c, sp.get(rr, c));
                    }
                }
                sp = grown;
                space.set(isf, sp);
            }
            Matrix prow = padLeft(row, sp.getNumCols());
            int h = Matrix.matchrow(sp, prow);
            if (h < 0) {
                Matrix newSp = Matrix.concatRows(sp, prow, null);
                space.set(isf, newSp);
                h = newSp.getNumRows() - 1;
            }
            hashed[isf] = h;
        }
        String key = hashKey(hashed);
        if (!seen.containsKey(key)) {
            seen.put(key, SSh.size());
            SSh.add(hashed);
            List<Matrix> stored = new ArrayList<Matrix>();
            for (int isf = 0; isf < nstateful; isf++) {
                stored.add(padLeft(stateCell.get(isf), space.get(isf).getNumCols()).copy());
            }
            stack.add(stored);
        }
    }
}
