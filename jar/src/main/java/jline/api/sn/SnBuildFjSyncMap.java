/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

/**
 * Builds a fork-join synchronization map from LINE's sn structure.
 *
 * <p>Port of matlab/src/lang/sn/sn_build_fj_sync_map.m. For each (Fork, Join)
 * pair it identifies which source nodes feed into the join and must therefore
 * be synchronized with mmap_max rather than merged as independent flows.</p>
 */
public final class SnBuildFjSyncMap {
    private SnBuildFjSyncMap() {}

    /**
     * Result of the fork-join synchronization scan.
     */
    public static final class FjSyncMap {
        /**
         * nodeSync(joinIdx, srcIdx) = groupId; groupId &gt; 0 means srcIdx belongs
         * to sync group groupId at joinIdx, 0 means srcIdx is an independent flow.
         */
        public Matrix nodeSync;
        /** forkOfGroup.get(g) = fork node index of group g+1 (0-based node index). */
        public List<Integer> forkOfGroup;
        /** joinOfGroup.get(g) = join node index of group g+1 (0-based node index). */
        public List<Integer> joinOfGroup;
        /** total number of sync groups. */
        public int nGroups;
    }

    public static FjSyncMap sn_build_fj_sync_map(NetworkStruct sn) {
        int I = sn.nnodes;
        int K = sn.nclasses;

        Matrix nodeSync = new Matrix(I, I);
        List<Integer> forkOfGroup = new ArrayList<Integer>();
        List<Integer> joinOfGroup = new ArrayList<Integer>();
        int groupId = 0;

        if (sn.fj != null) {
            for (int forkIdx = 0; forkIdx < sn.fj.getNumRows(); forkIdx++) {
                // MATLAB: forkIndices = find(any(sn.fj,2)) then joinIdx = find(sn.fj(forkIdx,:))
                for (int jnd = 0; jnd < sn.fj.getNumCols(); jnd++) {
                    if (sn.fj.get(forkIdx, jnd) <= 0) {
                        continue;
                    }
                    groupId = groupId + 1;
                    forkOfGroup.add(Integer.valueOf(forkIdx));
                    joinOfGroup.add(Integer.valueOf(jnd));

                    // A node is on the parallel path if the fork routes to it and
                    // it routes to the join, for at least one class.
                    for (int ind = 0; ind < I; ind++) {
                        if (ind == forkIdx || ind == jnd) {
                            continue;
                        }
                        boolean forkRoutesToNode = false;
                        boolean nodeRoutesToJoin = false;
                        for (int k = 0; k < K; k++) {
                            if (sn.rtnodes.get(forkIdx * K + k, ind * K + k) > 0) {
                                forkRoutesToNode = true;
                            }
                            if (sn.rtnodes.get(ind * K + k, jnd * K + k) > 0) {
                                nodeRoutesToJoin = true;
                            }
                        }
                        if (forkRoutesToNode && nodeRoutesToJoin) {
                            nodeSync.set(jnd, ind, groupId);
                        }
                    }
                }
            }
        }

        FjSyncMap out = new FjSyncMap();
        out.nodeSync = nodeSync;
        out.forkOfGroup = forkOfGroup;
        out.joinOfGroup = joinOfGroup;
        out.nGroups = groupId;
        return out;
    }
}
