/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.da;

import jline.api.cache.Cache_gamma_lp;
import jline.api.cache.Cache_miss_fpi;
import jline.api.mc.Dtmc_stochcomp;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.api.sn.SnRefreshVisits;

import java.util.ArrayList;
import java.util.List;

/**
 * Decomposition-aggregation driver for a CLOSED integrated cache-queueing model
 * whose Cache node has a delayed-hit retrieval system (Cache.setRetrievalSystem).
 * Java port of matlab/src/api/da/da_cacheqn_retrieval.m.
 *
 * Adapts the integrated cacheqn decomposition: the cache is relabelled as a
 * ClassSwitch (read -> hit / miss) and the finite-population delayed-hit
 * coalescing is made to emerge from the closed AMVA by giving the retrieval
 * (fetch) station a load-dependent COALESCING service rate lldscaling(k)=k/d(k),
 * d(k)=n_eff*(1-(1-1/n_eff)^k), n_eff=nitems-totalcap. The distinct-fetch
 * throughput saturates at 1/F and the coalescing benefit emerges from the finite
 * population.
 *
 * LIMITATIONS (EXPERIMENTAL - see the MATLAB reference): hitprob/missprob are the
 * true cache probabilities (hit=P(cached), miss=1-hit) and accurate; the delayed
 * fraction is NOT recovered (delayedprob=0, folds into miss); the coalescing
 * throughput benefit is captured only in direction and understated; single fetch
 * station (single-backend) only.
 *
 * NO CLOSED-RETRIEVAL EXAMPLE ships in the suite (every retrieval_* example is
 * OPEN). On a hand-built closed model this driver readily hits a reducible /
 * singular routing and a zero read-rate denominator (NaN read rate), so a plain
 * closed model can return an empty or unreliable table. Treat the closed path as
 * experimental and validate against LDES. The C++ SolverMVA (line-mp)
 * deliberately REFUSES this closed path by name; only the OPEN retrieval
 * analyzer is ported there. Verified 2026-07-24 that the open path matches
 * across codebases while the closed path lacks a clean, exampled reference.
 */
public class Da_cacheqn_retrieval {

    /** Network solver used for the aggregation step (mvald / ncld). */
    public interface NetSolve {
        SolverResult solve(NetworkStruct sn);
    }

    public static class Result {
        public SolverResult res;
        public Matrix hitprob;      // 1 x K on the read class
        public Matrix missprob;     // 1 x K
        public Matrix delayedprob;  // 1 x K
        public int it;
    }

    public static Result da_cacheqn_retrieval(NetworkStruct sn, NetSolve netfun, SolverOptions options) {
        final int I = sn.nnodes;
        final int K = sn.nclasses;

        // flat (node,class) indices of the stateful nodes (for stochastic complement)
        final List<Integer> statefulNodesClasses = new ArrayList<Integer>();
        for (int ind = 0; ind < I; ind++) {
            if (sn.isstateful.get(ind) != 0.0) {
                for (int r = 0; r < K; r++) statefulNodesClasses.add(ind * K + r);
            }
        }

        int ciTmp = -1;
        for (int i = 0; i < I; i++) if (sn.nodetype.get(i) == NodeType.Cache) { ciTmp = i; break; }
        if (ciTmp < 0) throw new RuntimeException("da_cacheqn_retrieval requires exactly one Cache node.");
        final int ci = ciTmp;
        final CacheNodeParam ch = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ci));

        // retrieval configuration (single read class, single fetch station)
        Integer readKey = ch.retrievalSystemQueueIndices.keySet().iterator().next();
        final int readClass = readKey;                    // 0-indexed read class
        List<Integer> queueNodes = ch.retrievalSystemQueueIndices.get(readKey);
        if (queueNodes.size() != 1) {
            throw new RuntimeException("da_cacheqn_retrieval currently supports a single-station (single-backend) retrieval system.");
        }
        final int fetchNode = queueNodes.get(0);
        final int fetchStation = (int) sn.nodeToStation.get(fetchNode);
        final int nitems = ch.nitems;
        double totcap = 0.0;
        for (int j = 0; j < ch.itemcap.length(); j++) totcap += ch.itemcap.get(j);
        final double n_eff = Math.max(1.0, nitems - totcap);

        // closed population
        int Npop = 0;
        for (int r = 0; r < K; r++) {
            double nj = sn.njobs.get(r);
            if (!Double.isInfinite(nj)) Npop += (int) Math.round(nj);
        }
        final int Nfin = Math.max(1, Npop);

        // mean fetch service F of the read class at the fetch station
        final int rcls0 = (int) ch.retrievalClasses.get(0, readClass);
        final double F = 1.0 / sn.rates.get(fetchStation, rcls0);

        // load-dependent COALESCING rate on the fetch station: lldscaling(k)=k/d(k)
        double[] alpha = new double[Nfin];
        for (int k = 1; k <= Nfin; k++) {
            double d = n_eff * (1.0 - Math.pow(1.0 - 1.0 / n_eff, k));
            alpha[k - 1] = k / d;
        }
        if (sn.lldscaling == null || sn.lldscaling.isEmpty()) {
            sn.lldscaling = new Matrix(sn.nstations, Nfin);
            sn.lldscaling.fill(1.0);
        } else if (sn.lldscaling.getNumCols() < Nfin) {
            Matrix grown = new Matrix(sn.nstations, Nfin);
            grown.fill(1.0);
            for (int i = 0; i < sn.lldscaling.getNumRows(); i++)
                for (int j = 0; j < sn.lldscaling.getNumCols(); j++)
                    grown.set(i, j, sn.lldscaling.get(i, j));
            sn.lldscaling = grown;
        }
        for (int k = 0; k < Nfin; k++) sn.lldscaling.set(fetchStation, k, alpha[k]);

        // relabel the cache as a class switch
        sn.nodetype.set(ci, NodeType.ClassSwitch);

        final Matrix hitClass = ch.hitclass;
        final Matrix missClass = ch.missclass;
        final Matrix retrievalClasses = ch.retrievalClasses;
        final double[] pread = new double[nitems];
        double preadSum = 0.0;
        for (int i = 0; i < nitems; i++) { pread[i] = ch.pread.get(readClass).get(i); preadSum += pread[i]; }
        for (int i = 0; i < nitems; i++) pread[i] /= preadSum;

        final int h = ch.itemcap.length();
        final Matrix mMat = ch.itemcap;

        final Matrix hitprob = new Matrix(1, K);
        final Matrix missprob = new Matrix(1, K);
        final Matrix delayedprob = new Matrix(1, K);
        final SolverResult[] resHolder = new SolverResult[1];

        double[] lambda0 = new double[K];
        lambda0[readClass] = 1.0;

        Da_fpi.Options<double[]> fpopts = new Da_fpi.Options<double[]>(
                options.iter_max, options.iter_tol,
                new Da_fpi.Norm<double[]>() {
                    @Override public double eval(double[] xnew, double[] xref) {
                        double d = 0.0;   // 1-norm
                        for (int i = 0; i < xnew.length; i++) d += Math.abs(xnew[i] - xref[i]);
                        return d;
                    }
                });

        Da_fpi.Result<double[]> fp = Da_fpi.run(new Da_fpi.Sweep<double[]>() {
            @Override public Da_fpi.SweepResult<double[]> sweep(double[] x, int it) {
                double[] lambda = x.clone();

                // isolated cache occupancy -> per-item uncached (would-be-miss) prob pi0
                Matrix[] lambda3d = new Matrix[1];
                lambda3d[0] = new Matrix(nitems, h + 1);
                for (int i = 0; i < nitems; i++)
                    for (int l = 0; l < h + 1; l++)
                        lambda3d[0].set(i, l, lambda[readClass] * pread[i]);
                Matrix[][] R = ch.accost;
                if (R == null || R.length == 0) {
                    R = new Matrix[1][nitems];
                    for (int k = 0; k < nitems; k++) {
                        Matrix Rmat = new Matrix(h + 1, h + 1);
                        for (int l = 0; l < h; l++) Rmat.set(l, l + 1, 1.0);
                        Rmat.set(h, h, 1.0);
                        R[0][k] = Rmat;
                    }
                }
                Matrix gamma = Cache_gamma_lp.cache_gamma_lp(lambda3d, R).gamma;
                Matrix pi0 = Cache_miss_fpi.cache_miss_fpi(gamma, mMat).pi0;   // n x 1

                double nonhit = 0.0;
                for (int i = 0; i < nitems; i++) nonhit += pread[i] * pi0.get(i, 0);
                double hp = 1.0 - nonhit;

                // see _kb/03-api-layer.md for rationale
                int r = readClass;
                for (int j = 0; j < sn.rtnodes.getNumCols(); j++) sn.rtnodes.set(ci * K + r, j, 0.0);
                sn.rtnodes.set(ci * K + r, ci * K + (int) hitClass.get(r), hp);
                for (int i = 0; i < nitems; i++) {
                    int rcls = (int) retrievalClasses.get(i, r);
                    if (rcls >= 0) {
                        sn.rtnodes.set(ci * K + r, fetchNode * K + rcls, pread[i] * pi0.get(i, 0));
                        for (int j = 0; j < sn.rtnodes.getNumCols(); j++) sn.rtnodes.set(fetchNode * K + rcls, j, 0.0);
                        sn.rtnodes.set(fetchNode * K + rcls, ci * K + (int) missClass.get(r), 1.0);
                        for (int j = 0; j < sn.rtnodes.getNumCols(); j++) sn.rtnodes.set(ci * K + rcls, j, 0.0);
                    }
                }
                sn.rt = Dtmc_stochcomp.dtmc_stochcomp(sn.rtnodes, statefulNodesClasses);

                SnRefreshVisits.snRefreshVisits(sn, sn.chains, sn.rt, sn.rtnodes);   // mutates sn in place

                SolverResult res = netfun.solve(sn);
                resHolder[0] = res;

                // aggregate nodevisits over chains
                Matrix nodevisits = null;
                for (Object key : sn.nodevisits.keySet()) {
                    if (nodevisits == null) nodevisits = sn.nodevisits.get(key);
                    else nodevisits = nodevisits.add(1.0, sn.nodevisits.get(key));
                }

                // throughput -> new read arrival rate at the cache
                int c = 0;
                while (c < sn.chains.getNumRows() && sn.chains.get(c, r) == 0.0) c++;
                double sumXN = 0.0;
                for (int j = 0; j < sn.chains.getNumCols(); j++)
                    if (sn.chains.get(c, j) != 0.0) sumXN += res.XN.get(j);
                int refnode = (int) sn.stationToNode.get((int) sn.refstat.get(r));
                double denom;
                if (sn.refclass.get(c) > -1) denom = nodevisits.get(refnode, (int) sn.refclass.get(c));
                else denom = nodevisits.get(refnode, r);
                lambda[r] = sumXN * nodevisits.get(ci, r) / denom;

                hitprob.set(0, r, hp);
                missprob.set(0, r, nonhit);
                delayedprob.set(0, r, 0.0);

                return new Da_fpi.SweepResult<double[]>(lambda, x);
            }
        }, lambda0, fpopts);

        Result out = new Result();
        out.res = resHolder[0];
        out.hitprob = hitprob;
        out.missprob = missprob;
        out.delayedprob = delayedprob;
        out.it = fp.it;
        return out;
    }
}
