/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import jline.api.da.Da_cacheqn_retrieval;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;

/**
 * NC analyzer for a CLOSED integrated cache-queueing model whose Cache node has a
 * delayed-hit retrieval system. Delegates to Da_cacheqn_retrieval with an NC
 * network solve (load-dependent ncld when the coalescing fetch station is
 * present, otherwise plain nc). Returns the true cache hit/miss probabilities
 * (hit = P(item cached), miss = 1 - hit); the delayed-hit fraction folds into
 * miss (delayedprob = 0). Java port of
 * matlab/src/solvers/NC/solver_nc_cacheqn_retrieval_analyzer.m.
 */
public final class Solver_nc_cacheqn_retrieval_analyzer {

    public static NCResult solver_nc_cacheqn_retrieval_analyzer(NetworkStruct sn, final SolverOptions options) {
        long t0 = System.nanoTime();
        int K = sn.nclasses;

        // Find the cache node before the driver relabels it to a ClassSwitch.
        Cache cache = null;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Cache) { cache = (Cache) sn.nodes.get(i); break; }
        }

        final Da_cacheqn_retrieval.NetSolve netfun = new Da_cacheqn_retrieval.NetSolve() {
            @Override public SolverResult solve(NetworkStruct snit) {
                if ((snit.lldscaling != null && !snit.lldscaling.isEmpty())
                        || (snit.cdscaling != null && !snit.cdscaling.isEmpty())
                        || (snit.jdscaling != null && !snit.jdscaling.isEmpty())) {
                    return Solver_ncld_analyzer.solver_ncld_analyzer(snit, options.copy());
                }
                return Solver_nc_analyzer.solver_nc_analyzer(snit, options.copy());
            }
        };

        Da_cacheqn_retrieval.Result dr = Da_cacheqn_retrieval.da_cacheqn_retrieval(sn, netfun, options);
        SolverResult r = dr.res;

        NCResult res = new NCResult();
        res.QN = r.QN; res.UN = r.UN; res.RN = r.RN; res.TN = r.TN;
        res.XN = r.XN; res.CN = r.CN; res.AN = r.AN; res.WN = r.WN;
        if (r instanceof NCResult) res.lG = ((NCResult) r).lG;
        res.hitProb = dr.hitprob;
        res.missProb = dr.missprob;
        res.it = dr.it;
        res.method = "fpi";

        if (cache != null) {
            int h = ((jline.lang.nodeparam.CacheNodeParam) sn.nodeparam.get(cache)).itemcap.length();
            Matrix hitProbList = new Matrix(K, h); hitProbList.fill(Double.NaN);
            Matrix latency = new Matrix(1, K); latency.fill(Double.NaN);
            cache.setResultHitProb(dr.hitprob);
            cache.setResultMissProb(dr.missprob);
            cache.setResultDelayedHitProb(dr.delayedprob);
            cache.setResultHitProbList(hitProbList);
            cache.setResultResidT(latency);
        }

        res.runtime = (System.nanoTime() - t0) / 1e9;
        return res;
    }
}
