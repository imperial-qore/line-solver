/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.NodeType;

/**
 * The single-station M/M/1/K loss system.
 *
 * <p>Port of {@code matlab/src/api/sn/sn_is_mm1k_loss.m}. True for a
 * single-class open Source-Queue-Sink system whose queue is a single-server
 * exponential M/M/1/K with tail drop ({@code DropStrategy.DROP}). This is the
 * exact regime of the closed-form loss scripts {@code qsys_mm1k_loss}
 * (probability-based, SolverNC) and {@code qsys_mg1k_loss_mgs} (moment-based,
 * SolverMVA), and the one truncated shape that keeps a product form over its
 * single station, hence the exemption in {@link SnHasBlocking}.
 */
public final class SnIsMm1kLoss {

    private SnIsMm1kLoss() {
    }

    /**
     * @param sn network structure
     * @return true if the model is a single-station M/M/1/K with tail drop
     */
    public static boolean snIsMm1kLoss(NetworkStruct sn) {
        if (sn.nclasses != 1 || sn.nclosedjobs != 0) {
            return false;
        }
        if (sn.nodetype == null || sn.nodetype.size() != 3) {
            return false;
        }
        int qnode = -1;
        int snode = -1;
        int nsink = 0;
        for (int a = 0; a < sn.nodetype.size(); a++) {
            NodeType t = sn.nodetype.get(a);
            if (t == NodeType.Queue) {
                if (qnode >= 0) {
                    return false;
                }
                qnode = a;
            } else if (t == NodeType.Source) {
                if (snode >= 0) {
                    return false;
                }
                snode = a;
            } else if (t == NodeType.Sink) {
                nsink++;
            } else {
                return false;
            }
        }
        if (qnode < 0 || snode < 0 || nsink != 1) {
            return false;
        }
        int qist = (int) sn.nodeToStation.get(qnode);
        int sist = (int) sn.nodeToStation.get(snode);
        if (qist < 0 || sist < 0) {
            return false;
        }
        if (sn.nservers == null || sn.nservers.get(qist) != 1) {
            return false;
        }
        if (sn.droprule == null || sn.stations == null || sn.jobclasses == null) {
            return false;
        }
        java.util.Map<jline.lang.JobClass, DropStrategy> rules = sn.droprule.get(sn.stations.get(qist));
        if (rules == null || rules.get(sn.jobclasses.get(0)) != DropStrategy.Drop) {
            return false;
        }
        // the JAR writes an unbounded capacity as Integer.MAX_VALUE where MATLAB writes Inf
        if (sn.cap == null || sn.cap.get(qist) <= 0
                || sn.cap.get(qist) >= (double) jline.GlobalConstants.MaxInt) {
            return false;
        }
        if (sn.scv == null || sn.scv.getNumRows() <= Math.max(qist, sist)) {
            return false;
        }
        if (Math.abs(sn.scv.get(sist, 0) - 1) > 1e-6 || Math.abs(sn.scv.get(qist, 0) - 1) > 1e-6) {
            return false;
        }
        // A DECLARED STATE DEPENDENCE IS NOT AN M/M/1/K. Both closed forms this
        // predicate gates read ONE service rate -- qsys_mm1k_loss a scalar mu,
        // qsys_mg1k_loss_mgs the first two moments of one service law -- so a
        // station carrying lldscaling, cdscaling or jdscaling would be answered
        // with the UNSCALED queue: a wrong number, not a coarse one, and silent.
        // The model is single-station by the tests above, so a model-wide test
        // IS a station test. See _kb/06-solver-catalog.md (NC, class-dependent
        // scaling).
        if (sn.lldscaling != null && !sn.lldscaling.isEmpty()) {
            return false;
        }
        if (sn.cdscaling != null && !sn.cdscaling.isEmpty()) {
            return false;
        }
        if (sn.jdscaling != null && !sn.jdscaling.isEmpty()) {
            return false;
        }
        return true;
    }
}
