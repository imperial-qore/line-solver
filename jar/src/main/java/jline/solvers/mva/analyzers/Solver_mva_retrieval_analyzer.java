/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import jline.api.retrieval.Cache_retrieval_inputs;
import jline.api.retrieval.Retrieval_fpi;
import jline.api.retrieval.Retrieval_fpi_latency;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.solvers.nc.analyzers.Solver_nc_retrieval_analyzer;
import jline.util.matrix.Matrix;

/**
 * Approximate analysis of a delayed-hit (retrieval-system) cache via the FPI
 * algorithms Retrieval_fpi (hit/miss/delayed-hit) and Retrieval_fpi_latency
 * (expected latency Z). Port of matlab/src/solvers/MVA/solver_mva_retrieval_analyzer.m.
 */
public final class Solver_mva_retrieval_analyzer {
    private Solver_mva_retrieval_analyzer() {}

    public static MVAResult solver_mva_retrieval_analyzer(NetworkStruct sn, SolverOptions options) {
        long t0 = System.nanoTime();
        MVAResult res = new MVAResult();
        int K = sn.nclasses;
        Matrix XN = new Matrix(1, K);
        Matrix QN = new Matrix(sn.nnodes, K);
        Matrix UN = new Matrix(sn.nnodes, K);
        Matrix RN = new Matrix(sn.nnodes, K);
        Matrix TN = new Matrix(sn.nnodes, K);
        Matrix CN = new Matrix(1, K);
        Matrix AN = new Matrix(sn.nnodes, K);
        Matrix WN = new Matrix(sn.nnodes, K);

        Cache_retrieval_inputs.Inputs in = Cache_retrieval_inputs.cache_retrieval_inputs(sn);
        int n = in.lambda.length;

        Retrieval_fpi.Result mm = Retrieval_fpi.retrieval_fpi(in.m, in.lambda, in.eta, in.gamma);
        Retrieval_fpi_latency.Result lat = Retrieval_fpi_latency.retrieval_fpi_latency(
                in.m, in.lambda, in.gamma, in.alpha, in.T, in.R, in.stationType);

        double tot = 0; for (double x : in.lambda) tot += x;
        double hitAgg = 0, missAgg = 0, delayedAgg = 0;
        for (int i = 0; i < n; i++) {
            double w = in.lambda[i] / tot;
            double ph = 0; for (double[] row : mm.phit) ph += row[i];
            double pd = 0; for (double[] row : mm.pdh) pd += row[i];
            hitAgg += w * ph; missAgg += w * mm.pmiss[i]; delayedAgg += w * pd;
        }

        Cache cache = null;
        for (int i = 0; i < sn.nodetype.size(); i++) if (sn.nodetype.get(i) == NodeType.Cache) { cache = (Cache) sn.nodes.get(i); break; }
        CacheNodeParam ch = (CacheNodeParam) sn.nodeparam.get(cache);
        Matrix hitProb = new Matrix(1, K), missProb = new Matrix(1, K), latency = new Matrix(1, K);
        Matrix delayedProb = new Matrix(1, K);
        hitProb.fill(Double.NaN); missProb.fill(Double.NaN); latency.fill(Double.NaN); delayedProb.fill(Double.NaN);
        hitProb.set(0, in.jobinClass, hitAgg);
        missProb.set(0, in.jobinClass, missAgg);
        delayedProb.set(0, in.jobinClass, delayedAgg);
        latency.set(0, in.jobinClass, lat.Z);
        // per-list (per-level) hit fractions for the read class: phit is (h x n);
        // access-weighted over items gives the per-list hit probability.
        int h = in.m.length;
        Matrix hitProbList = new Matrix(K, h);
        hitProbList.fill(Double.NaN);
        for (int l = 0; l < h; l++) {
            double acc = 0;
            for (int i = 0; i < n; i++) acc += (in.lambda[i] / tot) * mm.phit[l][i];
            hitProbList.set(in.jobinClass, l, acc);
        }
        // per-item occupancy [n x (h+1)]: column 0 = miss, columns 1.. = per-list.
        Matrix itemProb = new Matrix(n, h + 1);
        for (int i = 0; i < n; i++) {
            itemProb.set(i, 0, mm.pmiss[i]);
            for (int l = 0; l < h; l++) itemProb.set(i, l + 1, mm.phit[l][i]);
        }
        cache.setResultHitProb(hitProb);
        cache.setResultMissProb(missProb);
        cache.setResultDelayedHitProb(delayedProb);
        cache.setResultHitProbList(hitProbList);
        cache.setResultItemProb(itemProb);
        cache.setResultResidT(latency);
        res.hitProb = hitProb; res.missProb = missProb;

        int hc = (int) ch.hitclass.get(in.jobinClass);
        int mc = (int) ch.missclass.get(in.jobinClass);
        if (hc >= 0) XN.set(0, hc, in.sourceRate * (hitAgg + delayedAgg));
        if (mc >= 0) XN.set(0, mc, in.sourceRate * missAgg);

        int sourceNode = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) if (sn.nodetype.get(i) == NodeType.Source) { sourceNode = i; break; }
        // see _kb/06-solver-catalog.md for rationale
        int sourceStation = sourceNode >= 0 ? (int) sn.nodeToStation.get(sourceNode) : -1;
        if (sourceStation >= 0) for (int k = 0; k < K; k++) { double rt = sn.rates.get(sourceStation, k); if (!Double.isNaN(rt)) TN.set(sourceStation, k, rt); }

        Solver_nc_retrieval_analyzer.setStationMetrics(sn, in, mm.pdh, mm.pmiss, QN, UN, RN, TN);

        res.QN = QN; res.UN = UN; res.RN = RN; res.TN = TN; res.CN = CN; res.XN = XN; res.AN = AN; res.WN = WN;
        res.method = "fpi";
        res.runtime = (System.nanoTime() - t0) / 1e9;
        return res;
    }
}
