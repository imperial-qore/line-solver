/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.mam.Mmap_compress;
import jline.api.mam.Mmap_lambda;
import jline.api.mam.Mmap_max;
import jline.api.mam.Mmap_normalize;
import jline.api.mam.Mmap_super;
import jline.api.mc.Dtmc_stochcomp;
import jline.api.npfqn.Npfqn_traffic_merge;
import jline.api.npfqn.Npfqn_traffic_split_cs;
import jline.api.sn.SnBuildFjSyncMap;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * FJ-aware traffic solver extending solver_mam_traffic with mmap_max
 * synchronization at join points.
 *
 * <p>Port of matlab/src/solvers/MAM/solver_mam_traffic_mmap.m. DEP is indexed by
 * node (ind, r) rather than by station.</p>
 */
public final class Solver_mam_traffic_mmap {
    private Solver_mam_traffic_mmap() {}

    /**
     * @param sn        the network struct
     * @param DEP       departure processes, DEP.get(ind).get(r) in (D0,D1) format
     * @param config    solver configuration
     * @param fjSyncMap fork-join synchronization map
     * @return ARV.get(ind), the arrival MMAP at node ind (empty cell when none)
     */
    public static Map<Integer, MatrixCell> solver_mam_traffic_mmap(NetworkStruct sn,
                                                                   Map<Integer, Map<Integer, MatrixCell>> DEP,
                                                                   SolverOptions.Config config,
                                                                   SnBuildFjSyncMap.FjSyncMap fjSyncMap) {
        int I = sn.nnodes;
        int R = sn.nclasses;

        // MATLAB defaults fj_sync_q_len to 2 only when the field is absent; an
        // explicitly configured value is used as given.
        int fjSyncQLen = 2;
        Object fjq = config.get("fj_sync_q_len");
        if (fjq instanceof Integer) {
            fjSyncQLen = ((Integer) fjq).intValue();
        }

        // Index over all non-ClassSwitch nodes
        List<Integer> non_cs_classes = new ArrayList<Integer>();
        boolean[] isNCS = new boolean[I];
        int[] nodeToNCS = new int[I];
        int ncsCount = 0;
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) != NodeType.ClassSwitch) {
                for (int i = 0; i < R; i++) {
                    non_cs_classes.add(Integer.valueOf(ind * R + i));
                }
                isNCS[ind] = true;
                ncsCount++;
                nodeToNCS[ind] = ncsCount; // 1-based, as in solver_mam_traffic
            } else {
                isNCS[ind] = false;
            }
        }

        // Hide the nodes that are class switches
        Matrix rtncs = Dtmc_stochcomp.dtmc_stochcomp(sn.rtnodes, non_cs_classes);
        int Inc = ncsCount;

        // DEP is indexed by node (ind, r) - convert to MMAP format
        Map<Integer, Map<Integer, MatrixCell>> MMAP = new HashMap<Integer, Map<Integer, MatrixCell>>();
        for (int ind = 0; ind < I; ind++) {
            Map<Integer, MatrixCell> row = new HashMap<Integer, MatrixCell>();
            for (int r = 0; r < R; r++) {
                MatrixCell src = (DEP.get(ind) != null) ? DEP.get(ind).get(r) : null;
                MatrixCell cell;
                if (src == null || src.isEmpty() || src.get(0).hasNaN()) {
                    // no arrivals from this class
                    cell = new MatrixCell(3);
                    cell.set(0, new Matrix(1, 1, 0));
                    cell.set(1, new Matrix(1, 1, 0));
                    cell.set(2, new Matrix(1, 1, 0));
                } else {
                    cell = new MatrixCell(src);
                    cell.set(2, cell.get(1));
                }
                row.put(r, cell);
            }
            MMAP.put(ind, row);
        }

        Map<Integer, MatrixCell> ARV = new HashMap<Integer, MatrixCell>();
        Map<Integer, MatrixCell> DEP_NCS = new HashMap<Integer, MatrixCell>();
        Map<Integer, Map<Integer, MatrixCell>> LINKS = new HashMap<Integer, Map<Integer, MatrixCell>>();

        // Build the nodeSync matrix in NCS indexing
        Matrix nodeSyncNCS = new Matrix(Inc, Inc);
        for (int ind = 0; ind < I; ind++) {
            if (!isNCS[ind]) continue;
            int inc = nodeToNCS[ind];
            for (int jnd = 0; jnd < I; jnd++) {
                if (!isNCS[jnd]) continue;
                int jnc = nodeToNCS[jnd];
                double g = fjSyncMap.nodeSync.get(ind, jnd);
                if (g > 0) {
                    nodeSyncNCS.set(inc - 1, jnc - 1, g);
                }
            }
        }

        // First determine all outgoing flows from all nodes
        for (int ind = 0; ind < I; ind++) {
            if (!isNCS[ind]) continue;
            NodeType nt = sn.nodetype.get(ind);
            if (nt != NodeType.Source && nt != NodeType.Delay && nt != NodeType.Queue
                    && nt != NodeType.Fork && nt != NodeType.Join) {
                continue;
            }
            int inc = nodeToNCS[ind];

            MatrixCell dep;
            if (R > 1) {
                // see _kb/06-solver-catalog.md for rationale
                dep = MMAP.get(ind).get(0);
                for (int rr = 1; rr < R; rr++) {
                    dep = Mmap_super.mmap_super(dep, MMAP.get(ind).get(rr));
                    if (dep.get(0).getNumRows() > config.space_max) {
                        dep = compress(dep, config.compress);
                    }
                }
            } else {
                dep = MMAP.get(ind).get(0);
            }
            DEP_NCS.put(inc, dep);

            Matrix Psplit = new Matrix(R, Inc * R, R * Inc * R);
            for (int r = 0; r < R; r++) {
                for (int jnd = 0; jnd < I; jnd++) {
                    if (!isNCS[jnd]) continue;
                    int jnc = nodeToNCS[jnd];
                    for (int s = 0; s < R; s++) {
                        Psplit.set(r, (jnc - 1) * R + s,
                                rtncs.get((inc - 1) * R + r, (jnc - 1) * R + s));
                    }
                }
            }

            Map<Integer, MatrixCell> Fsplit =
                    Npfqn_traffic_split_cs.npfqn_traffic_split_cs(DEP_NCS.get(inc), Psplit);
            LINKS.put(inc, new HashMap<Integer, MatrixCell>());
            for (int jnc = 0; jnc < Inc; jnc++) {
                LINKS.get(inc).put(jnc, Mmap_normalize.mmap_normalize(Fsplit.get(jnc)));
            }
        }

        // Then determine all incoming flows, with FJ synchronization
        for (int ind = 0; ind < I; ind++) {
            if (!isNCS[ind] || sn.nodetype.get(ind) == NodeType.Source) {
                ARV.put(ind, new MatrixCell());
                continue;
            }
            int inc = nodeToNCS[ind];

            // see _kb/06-solver-catalog.md for rationale
            List<MatrixCell> independentFlows = new ArrayList<MatrixCell>();
            Map<Integer, List<MatrixCell>> syncFlows = new LinkedHashMap<Integer, List<MatrixCell>>();
            List<Integer> syncGroupsAtNode = new ArrayList<Integer>();
            for (int jnc = 0; jnc < Inc; jnc++) {
                int gid = (int) nodeSyncNCS.get(inc - 1, jnc);
                if (gid > 0 && !syncGroupsAtNode.contains(Integer.valueOf(gid))) {
                    syncGroupsAtNode.add(Integer.valueOf(gid));
                }
            }
            java.util.Collections.sort(syncGroupsAtNode);

            for (int jnc = 1; jnc <= Inc; jnc++) {
                Map<Integer, MatrixCell> outLinks = LINKS.get(jnc);
                MatrixCell flow = (outLinks != null) ? outLinks.get(inc - 1) : null;
                if (flow == null || flow.isEmpty()
                        || Mmap_lambda.mmap_lambda(flow).elementSum() <= GlobalConstants.FineTol) {
                    continue;
                }
                int gid = (int) nodeSyncNCS.get(inc - 1, jnc - 1);
                if (gid == 0) {
                    independentFlows.add(flow);
                } else {
                    if (!syncFlows.containsKey(Integer.valueOf(gid))) {
                        syncFlows.put(Integer.valueOf(gid), new ArrayList<MatrixCell>());
                    }
                    syncFlows.get(Integer.valueOf(gid)).add(flow);
                }
            }

            // Process synchronized flows: apply mmap_max iteratively within each group
            List<MatrixCell> syncResults = new ArrayList<MatrixCell>();
            for (int gi = 0; gi < syncGroupsAtNode.size(); gi++) {
                Integer gid = syncGroupsAtNode.get(gi);
                List<MatrixCell> groupFlows = syncFlows.get(gid);
                if (groupFlows == null || groupFlows.isEmpty()) {
                    continue;
                }
                MatrixCell syncedFlow = groupFlows.get(0);
                for (int f = 1; f < groupFlows.size(); f++) {
                    syncedFlow = Mmap_max.mmap_max(syncedFlow, groupFlows.get(f), fjSyncQLen);
                    // mmap_max builds the synchronization state space but does not
                    // enforce MMAP feasibility; normalize before compressing.
                    syncedFlow = Mmap_normalize.mmap_normalize(syncedFlow);
                    if (syncedFlow.get(0).getNumRows() > config.space_max) {
                        syncedFlow = compress(syncedFlow, config.compress);
                    }
                }
                syncResults.add(syncedFlow);
            }

            // Merge synced flows with independent flows
            List<MatrixCell> allFlows = new ArrayList<MatrixCell>();
            allFlows.addAll(syncResults);
            allFlows.addAll(independentFlows);

            if (allFlows.size() > 1) {
                Map<Integer, MatrixCell> FLOWS = new HashMap<Integer, MatrixCell>();
                for (int i = 0; i < allFlows.size(); i++) {
                    FLOWS.put(Integer.valueOf(i), allFlows.get(i));
                }
                ARV.put(ind, Npfqn_traffic_merge.npfqn_traffic_merge(FLOWS, config.merge, config.compress));
            } else if (allFlows.size() == 1) {
                ARV.put(ind, allFlows.get(0));
            } else {
                // No flows: take the FIRST non-empty link (MATLAB breaks on match)
                MatrixCell fallback = null;
                for (int jnc = 1; jnc <= Inc; jnc++) {
                    Map<Integer, MatrixCell> outLinks = LINKS.get(jnc);
                    MatrixCell cand = (outLinks != null) ? outLinks.get(inc - 1) : null;
                    if (cand != null && !cand.isEmpty()) {
                        fallback = cand;
                        break;
                    }
                }
                if (fallback == null) {
                    fallback = new MatrixCell(3);
                    fallback.set(0, new Matrix(1, 1, 0));
                    fallback.set(1, new Matrix(1, 1, 0));
                    fallback.set(2, new Matrix(1, 1, 0));
                }
                ARV.put(ind, fallback);
            }
        }

        return ARV;
    }

    /**
     * Adapts Mmap_compress (which operates on Matrix[]) to MatrixCell.
     */
    static MatrixCell compress(MatrixCell cell, String method) {
        Matrix[] arr = new Matrix[cell.size()];
        for (int i = 0; i < cell.size(); i++) {
            arr[i] = cell.get(i);
        }
        Matrix[] out = Mmap_compress.mmap_compress(arr, method);
        MatrixCell res = new MatrixCell(out.length);
        for (int i = 0; i < out.length; i++) {
            res.set(i, out[i]);
        }
        return res;
    }
}
