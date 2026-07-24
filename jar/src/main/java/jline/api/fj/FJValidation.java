/**
 * Fork-Join topology validation
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.api.fj;

import java.util.ArrayList;
import java.util.List;

import jline.util.Pair;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;

public final class FJValidation {
    private FJValidation() {}

    /**
     * Checks for a single fork-join pair with homogeneous parallel branches.
     *
     * <p>This is NOT a test for the presence of Fork/Join nodes: most fork-join
     * models fail this predicate. It tests membership in the homogeneous class on
     * which the FJ_codes tail approximation of Qiu, Perez and Harrison (IFIP
     * Performance 2015) is defined: a single fork-join pair, K homogeneous
     * single-server parallel queues between them, open classes only, and service
     * in {Exp, HyperExp(2), Erlang(2), MAP(2)}.
     *
     * @param sn the network structure to test
     * @return the verdict paired with the extracted fork-join indexes (null if the
     *         model is outside the class)
     */
    public static Pair<Boolean, FJInfo> isHomogeneous(NetworkStruct sn) {
        try {
            if (sn.nclosedjobs > 0) {
                return new Pair<Boolean, FJInfo>(Boolean.FALSE, null);
            }

            int sourceIdx = -1;
            int sinkIdx = -1;
            int forkIdx = -1;
            int joinIdx = -1;

            for (int i = 0; i < sn.nnodes; i++) {
                Node node = sn.nodes.get(i);
                if (node instanceof Source) {
                    sourceIdx = i;
                } else if (node instanceof Sink) {
                    sinkIdx = i;
                } else if (node instanceof Fork) {
                    forkIdx = i;
                } else if (node instanceof Join) {
                    joinIdx = i;
                }
            }

            if (sourceIdx < 0 || sinkIdx < 0 || forkIdx < 0 || joinIdx < 0) {
                return new Pair<Boolean, FJInfo>(Boolean.FALSE, null);
            }

            List<Integer> queueIndices = new ArrayList<Integer>();

            for (int i = 0; i < sn.nnodes; i++) {
                Node node = sn.nodes.get(i);
                if (node instanceof Queue || node instanceof Delay) {
                    queueIndices.add(Integer.valueOf(i));
                }
            }

            int K = queueIndices.size();
            if (K < 1) {
                return new Pair<Boolean, FJInfo>(Boolean.FALSE, null);
            }

            if (K > 1) {
                int firstQueueStationIdx = (int) sn.nodeToStation.get(queueIndices.get(0).intValue());
                for (int r = 0; r < sn.nclasses; r++) {
                    double refRate = sn.rates.get(firstQueueStationIdx, r);
                    java.util.Map<JobClass, ProcessType> refMap = sn.procid.get(sn.stations.get(firstQueueStationIdx));
                    ProcessType refProcId = (refMap == null) ? null : refMap.get(sn.jobclasses.get(r));
                    for (int q = 1; q < K; q++) {
                        int qStationIdx = (int) sn.nodeToStation.get(queueIndices.get(q).intValue());
                        double qRate = sn.rates.get(qStationIdx, r);
                        java.util.Map<JobClass, ProcessType> qMap = sn.procid.get(sn.stations.get(qStationIdx));
                        ProcessType qProcId = (qMap == null) ? null : qMap.get(sn.jobclasses.get(r));
                        if (Math.abs(refRate - qRate) > 1e-10 * Math.max(1.0, Math.abs(refRate))) {
                            return new Pair<Boolean, FJInfo>(Boolean.FALSE, null);
                        }
                        if (refProcId == null ? qProcId != null : !refProcId.equals(qProcId)) {
                            return new Pair<Boolean, FJInfo>(Boolean.FALSE, null);
                        }
                    }
                }
            }

            for (int i = 0; i < queueIndices.size(); i++) {
                int qIdx = queueIndices.get(i).intValue();
                int stationIdx = (int) sn.nodeToStation.get(qIdx);
                Station station = sn.stations.get(stationIdx);
                SchedStrategy sched = sn.sched.get(station);
                if (sched != SchedStrategy.FCFS && sched != SchedStrategy.PS) {
                    return new Pair<Boolean, FJInfo>(Boolean.FALSE, null);
                }
            }

            int[] queueIndicesArr = new int[queueIndices.size()];
            for (int i = 0; i < queueIndices.size(); i++) {
                queueIndicesArr[i] = queueIndices.get(i).intValue();
            }

            FJInfo fjInfo = new FJInfo(K, forkIdx, joinIdx, queueIndicesArr, sourceIdx, sinkIdx, true);

            return new Pair<Boolean, FJInfo>(Boolean.TRUE, fjInfo);

        } catch (Exception e) {
            return new Pair<Boolean, FJInfo>(Boolean.FALSE, null);
        }
    }
}
