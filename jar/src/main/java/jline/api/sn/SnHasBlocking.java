/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;

/**
 * Finite-buffer blocking or loss anywhere in the network.
 *
 * <p>Port of {@code matlab/src/api/sn/sn_has_blocking.m}. True when some
 * station can refuse a job, either because its own buffer BINDS (Kendall's K
 * below the population that can reach it, whatever the drop rule: WAITQ, DROP,
 * BAS, BBS, RSRD) or because a finite capacity region caps a set of stations
 * jointly. Such a network is not product form: the truncation couples the
 * station occupancies, so no BCMP factorization of the equilibrium
 * distribution exists.
 *
 * <p>Only a buffer that can actually BIND counts, which is what
 * {@link SnGetBufferSize} decides: refreshCapacity derives a finite classcap
 * (the chain population) at every station of every closed model, so a plain
 * finiteness test would call every closed model blocking.
 *
 * <p>Two shapes are exempt. A Cache builds its own capped retrieval queues
 * (classCap = 1), which the cache analyzers solve rather than treat as a buffer
 * constraint. And the single-station M/M/1/K loss system keeps the truncated
 * geometric distribution, a product form over its one station, evaluated in
 * closed form by the loss branches of SolverMVA and SolverNC.
 */
public final class SnHasBlocking {

    private SnHasBlocking() {
    }

    /**
     * @param sn network structure
     * @return true if the network has binding finite buffers or capacity regions
     */
    public static boolean snHasBlocking(NetworkStruct sn) {
        // a finite capacity region caps a SET of stations, which no per-station
        // capacity can express and no product form survives
        if (sn.nregions > 0) {
            return true;
        }
        if (sn.nodetype != null && sn.nodetype.contains(NodeType.Cache)) {
            return false;
        }
        if (SnIsMm1kLoss.snIsMm1kLoss(sn)) {
            return false;
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            if (Double.isFinite(SnGetBufferSize.snGetBufferSize(sn, ist))) {
                return true;
            }
        }
        return false;
    }
}
