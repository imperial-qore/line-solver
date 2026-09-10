/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import jline.api.mam.Map_exponential;
import jline.api.polling.Polling_qsys_1limited;
import jline.api.polling.Polling_qsys_decrementing;
import jline.api.polling.Polling_qsys_exhaustive;
import jline.api.polling.Polling_qsys_gated;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.lang.nodes.Station;
import jline.lang.constant.NodeType;
import jline.lang.constant.PollingType;
import jline.lang.NodeParam;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mva_polling_analyzer {
    private Solver_mva_polling_analyzer() {}

    /**
     * Converts a Distribution to a MAP representation.
     * Currently supports Exponential and Immediate distributions.
     */
    private static MatrixCell distributionToMAP(Distribution dist) {
        if (dist instanceof Exp) {
            return Map_exponential.map_exponential(dist.getMean());
        }
        // Use very small mean (1e-10) for immediate/null distributions to avoid division by zero
        // in polling analysis equations. This matches Python native behavior.
        if (dist instanceof Immediate) {
            return Map_exponential.map_exponential(1e-10);
        }
        if (dist == null) {
            return Map_exponential.map_exponential(1e-10);
        }
        return Map_exponential.map_exponential(dist.getMean());
    }

    /**
     * MVA Polling System analyzer
     */
    public static MVAResult solver_mva_polling_analyzer(NetworkStruct sn, SolverOptions options) {
        long startTime = System.nanoTime();
        String method = options.method;
        int totiter = 1;

        Matrix QN = new Matrix(sn.nstations, sn.nclasses);
        Matrix UN = new Matrix(sn.nstations, sn.nclasses);
        Matrix RN = new Matrix(sn.nstations, sn.nclasses);
        Matrix TN = new Matrix(sn.nstations, sn.nclasses);
        Matrix CN = new Matrix(sn.nstations, sn.nclasses);
        Matrix AN = new Matrix(sn.nstations, sn.nclasses);
        Matrix WN = new Matrix(sn.nstations, sn.nclasses);
        Matrix XN = new Matrix(sn.nstations, sn.nclasses);
        double lG = 0.0;

        int source_ist = -1;
        int queue_ist = -1;
        int queue_node = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                source_ist = (int) sn.nodeToStation.get(i);
            } else if (sn.nodetype.get(i) == NodeType.Queue) {
                queue_ist = (int) sn.nodeToStation.get(i);
                queue_node = i;
            }
        }

        if (source_ist == -1 || queue_ist == -1) {
            throw new RuntimeException("Polling analyzer requires a Source and a Queue node");
        }

        Matrix lambda = new Matrix(1, sn.nclasses);
        Matrix mu = new Matrix(1, sn.nclasses);
        int k = (int) sn.nservers.get(queue_ist);

        for (int r = 0; r < sn.nclasses; r++) {
            // Get arrival rate directly from rates matrix (matches Python native behavior)
            lambda.set(0, r, sn.rates.get(source_ist, r));
            mu.set(0, r, sn.rates.get(queue_ist, r));
        }

        Station sourceStation = sn.stations.get(source_ist);
        Station queueStation = sn.stations.get(queue_ist);

        if (sourceStation == null || queueStation == null) {
            throw new RuntimeException("Source or queue station not found in stations map");
        }

        MatrixCell[] arvMAPs = new MatrixCell[sn.nclasses];
        for (int r = 0; r < sn.nclasses; r++) {
            JobClass jobClass = sn.jobclasses.get(r);
            java.util.Map<JobClass, MatrixCell> procMap = sn.proc.get(sourceStation);
            MatrixCell mc = (procMap != null) ? procMap.get(jobClass) : null;
            if (mc == null) {
                throw new RuntimeException("No arrival MAP for class " + r);
            }
            arvMAPs[r] = mc;
        }

        MatrixCell[] svcMAPs = new MatrixCell[sn.nclasses];
        for (int r = 0; r < sn.nclasses; r++) {
            JobClass jobClass = sn.jobclasses.get(r);
            java.util.Map<JobClass, MatrixCell> procMap = sn.proc.get(queueStation);
            MatrixCell mc = (procMap != null) ? procMap.get(jobClass) : null;
            if (mc == null) {
                throw new RuntimeException("No service MAP for class " + r);
            }
            svcMAPs[r] = mc;
        }

        MatrixCell[] switchMAPs = new MatrixCell[sn.nclasses];
        for (int r = 0; r < sn.nclasses; r++) {
            JobClass jobClass = sn.jobclasses.get(r);
            Node queueNodeRef = sn.nodes.get(queue_node);
            NodeParam nodeParam = sn.nodeparam.get(queueNodeRef);

            if (nodeParam instanceof QueueNodeParam
                    && ((QueueNodeParam) nodeParam).switchoverTime != null
                    && ((QueueNodeParam) nodeParam).switchoverTime.containsKey(jobClass)) {
                Distribution switchoverDist = ((QueueNodeParam) nodeParam).switchoverTime.get(jobClass);
                switchMAPs[r] = distributionToMAP(switchoverDist);
            } else {
                switchMAPs[r] = distributionToMAP(null);
            }
        }

        Node queueNodeRef = sn.nodes.get(queue_node);
        NodeParam nodeParam = sn.nodeparam.get(queueNodeRef);
        PollingType pollingType = PollingType.EXHAUSTIVE;
        int pollingK = 1;

        if (nodeParam instanceof QueueNodeParam && ((QueueNodeParam) nodeParam).pollingType != null) {
            pollingType = ((QueueNodeParam) nodeParam).pollingType;
            if (pollingType == PollingType.KLIMITED && ((QueueNodeParam) nodeParam).pollingPar != null) {
                pollingK = ((QueueNodeParam) nodeParam).pollingPar.intValue();
            }
        }

        if ("exact".equalsIgnoreCase(method)) {
            Matrix ca = new Matrix(1, sn.nclasses);
            for (int r = 0; r < sn.nclasses; r++) {
                ca.set(0, r, Math.sqrt(sn.scv.get(source_ist, r)));
            }
            if (k == 1 && ca.elementMax() == 1.0 && ca.elementMin() == 1.0
                    && (pollingType == PollingType.EXHAUSTIVE || pollingType == PollingType.GATED)) {
                method = "stationtime";
            } else {
                throw new RuntimeException("MVA exact method unavailable for this model.");
            }
        }

        if ("default".equals(method)) {
            method = "stationtime";
        }

        Matrix W = new Matrix(1, sn.nclasses);
        Matrix R = new Matrix(1, sn.nclasses);

        if ("stationtime".equals(method)) {
            double[] waitTimes;
            if (pollingType == PollingType.EXHAUSTIVE) {
                waitTimes = Polling_qsys_exhaustive.polling_qsys_exhaustive(arvMAPs, svcMAPs, switchMAPs);
            } else if (pollingType == PollingType.GATED) {
                waitTimes = Polling_qsys_gated.polling_qsys_gated(arvMAPs, svcMAPs, switchMAPs);
            } else if (pollingType == PollingType.KLIMITED) {
                if (pollingK == 1) {
                    waitTimes = Polling_qsys_1limited.polling_qsys_1limited(arvMAPs, svcMAPs, switchMAPs);
                } else {
                    throw new RuntimeException("MVA method unavailable for K-limited polling with K>1.");
                }
            } else if (pollingType == PollingType.DECREMENTING) {
                waitTimes = Polling_qsys_decrementing.polling_qsys_decrementing(arvMAPs, svcMAPs, switchMAPs);
            } else {
                throw new RuntimeException("Unsupported polling type: " + pollingType);
            }
            for (int r = 0; r < sn.nclasses; r++) {
                W.set(0, r, waitTimes[r]);
                R.set(0, r, waitTimes[r] + 1.0 / mu.get(0, r));
            }
        } else {
            throw new RuntimeException("Unsupported polling solution method: " + method);
        }

        for (int r = 0; r < sn.nclasses; r++) {
            // In polling systems, each class visits the queue exactly once
            // No need to multiply by visit ratio (matches Python native behavior)
            RN.set(queue_ist, r, R.get(0, r));
            CN.set(queue_ist, r, R.get(0, r));
            XN.set(queue_ist, r, lambda.get(0, r));
            UN.set(queue_ist, r, lambda.get(0, r) / mu.get(0, r) / k);
            TN.set(source_ist, r, lambda.get(0, r));
            TN.set(queue_ist, r, lambda.get(0, r));
            QN.set(queue_ist, r, XN.get(queue_ist, r) * RN.get(queue_ist, r));
        }

        MVAResult res = new MVAResult();
        res.method = method;
        res.QN = QN;
        res.RN = RN;
        res.XN = XN;
        res.UN = UN;
        res.TN = TN;
        res.CN = CN;
        res.AN = AN;
        res.WN = WN;
        res.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        res.iter = totiter;
        res.logNormConstAggr = lG;

        return res;
    }
}
