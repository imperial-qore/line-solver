/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import jline.api.retrieval.Cache_retrieval_inputs;
import jline.api.retrieval.Retrieval_metrics;
import jline.api.retrieval.Retrieval_nc;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;

/**
 * Exact analysis of a delayed-hit (retrieval-system) cache via the product-form
 * algorithms Retrieval_nc (normalizing constant) and Retrieval_metrics (hit/miss/
 * delayed-hit). Latency is left to SolverMVA (Retrieval_fpi_latency); it is NaN here.
 * Port of matlab/src/solvers/NC/solver_nc_retrieval_analyzer.m.
 */
public final class Solver_nc_retrieval_analyzer {
    private Solver_nc_retrieval_analyzer() {}

    public static NCResult solver_nc_retrieval_analyzer(NetworkStruct sn, SolverOptions options) {
        long t0 = System.nanoTime();
        NCResult res = new NCResult();
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
        int r = in.eta[0].length - 1;

        // exact normalizing constant E(m) = retrieval_nc(0,m,...)
        double E = Retrieval_nc.retrieval_nc(new double[r], in.m, in.lambda, in.eta, in.gamma);
        res.lG = Math.log(E);

        Retrieval_metrics.Result mm = Retrieval_metrics.retrieval_metrics(in.m, in.lambda, in.eta, in.gamma);

        double tot = 0; for (double x : in.lambda) tot += x;
        double hitAgg = 0, missAgg = 0, delayedAgg = 0;
        for (int i = 0; i < n; i++) {
            double w = in.lambda[i] / tot;
            double ph = 0; for (double[] row : mm.phit) ph += row[i];
            double pd = 0; for (double[] row : mm.pdh) pd += row[i];
            hitAgg += w * ph; missAgg += w * mm.pmiss[i]; delayedAgg += w * pd;
        }

        // cache node hit/miss/latency (per-class vectors; read-class entry set, rest NaN)
        Cache cache = null; int ci = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) if (sn.nodetype.get(i) == NodeType.Cache) { cache = (Cache) sn.nodes.get(i); ci = i; break; }
        CacheNodeParam ch = (CacheNodeParam) sn.nodeparam.get(cache);
        Matrix hitProb = new Matrix(1, K), missProb = new Matrix(1, K), latency = new Matrix(1, K);
        Matrix delayedProb = new Matrix(1, K);
        hitProb.fill(Double.NaN); missProb.fill(Double.NaN); latency.fill(Double.NaN); delayedProb.fill(Double.NaN);
        hitProb.set(0, in.jobinClass, hitAgg);
        missProb.set(0, in.jobinClass, missAgg);
        delayedProb.set(0, in.jobinClass, delayedAgg);
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
        cache.setResultResidT(latency);   // NaN: NC has no latency algorithm
        res.hitProb = hitProb; res.missProb = missProb;

        // throughputs: misses -> missclass, hits+delayed -> hitclass
        int hc = (int) ch.hitclass.get(in.jobinClass);
        int mc = (int) ch.missclass.get(in.jobinClass);
        if (hc >= 0) XN.set(0, hc, in.sourceRate * (hitAgg + delayedAgg));
        if (mc >= 0) XN.set(0, mc, in.sourceRate * missAgg);

        // source throughput
        int sourceNode = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) if (sn.nodetype.get(i) == NodeType.Source) { sourceNode = i; break; }
        // Station index on BOTH sides: the rate is read at the Source's station
        // row, so the throughput must be written there too. They coincide only
        // while the Source is station 0.
        int sourceStation = sourceNode >= 0 ? (int) sn.nodeToStation.get(sourceNode) : -1;
        if (sourceStation >= 0) for (int k = 0; k < K; k++) { double rt = sn.rates.get(sourceStation, k); if (!Double.isNaN(rt)) TN.set(sourceStation, k, rt); }

        // retrieval-station occupancy (QN=UN=phi_s) and throughput
        setStationMetrics(sn, in, mm.pdh, mm.pmiss, QN, UN, RN, TN);

        res.QN = QN; res.UN = UN; res.RN = RN; res.TN = TN; res.CN = CN; res.XN = XN; res.AN = AN; res.WN = WN;
        res.method = "exact";
        res.runtime = (System.nanoTime() - t0) / 1e9;
        return res;
    }

    /**
     * Shared station-metric population (QN=UN=phi_s, TN=fetch throughput, RN=Little).
     *
     * @param pmiss per-item miss probability pi_{i,0}: only a MISS is fetched
     *              through a retrieval station, so it weights the throughput.
     */
    public static void setStationMetrics(NetworkStruct sn, Cache_retrieval_inputs.Inputs in, double[][] pdh,
                                  double[] pmiss, Matrix QN, Matrix UN, Matrix RN, Matrix TN) {
        int n = in.lambda.length;
        int S = in.stationType.length;
        double tot = 0; for (double x : in.lambda) tot += x;
        // psIdx mapping (pdh row 1+p)
        int r = 0; int[] psIdx = new int[S];
        for (int s = 0; s < S; s++) if (!in.stationType[s].equals("IS")) psIdx[r++] = s;
        for (int s = 0; s < S; s++) {
            // QN/UN/RN/TN are indexed by STATION, not by node. MATLAB's
            // solver_nc_retrieval_analyzer.m converts with
            // sst = sn.nodeToStation(queueNodes(s)) before writing, and
            // Cache_retrieval_inputs already carries that conversion in
            // queueStation. Writing at the node index put the retrieval
            // station's metrics outside the station range whenever a non-station
            // node (a Cache, here) precedes it, so the row silently vanished
            // from the table and the model reported only its Source.
            int qstation = in.queueStation[s];
            double phi_s = 0;
            if (in.stationType[s].equals("IS")) {
                for (int i = 0; i < n; i++) phi_s += pdh[0][i];
            } else {
                int p = -1; for (int q = 0; q < r; q++) if (psIdx[q] == s) { p = q; break; }
                if (p >= 0 && 1 + p < pdh.length) for (int i = 0; i < n; i++) phi_s += pdh[1 + p][i];
            }
            double tput = fetchTput(in, pmiss, s, S);
            if (qstation < 0) continue;   // not a station: nothing to tabulate
            QN.set(qstation, in.jobinClass, phi_s);
            UN.set(qstation, in.jobinClass, phi_s);
            TN.set(qstation, in.jobinClass, tput);
            if (tput > 0) RN.set(qstation, in.jobinClass, phi_s / tput);
        }
    }

    private static double fetchTput(Cache_retrieval_inputs.Inputs in, double[] pmiss, int s, int S) {
        // MATLAB solver_nc_retrieval_analyzer.m:
        //   tput_s += sourceRate(jobinClass) * (lambda(i)/sum(lambda)) * pi0(i) * vis(s)
        // Only a MISS is fetched through a retrieval station, so the arrival
        // rate of fetches is the request rate times the ACCESS-WEIGHTED MISS
        // probability, not the request rate itself. Omitting pi0 (and the
        // sourceRate and the lambda normalisation) reported every request as a
        // fetch, which on a model whose access probabilities already sum to one
        // returns exactly the source rate: 1.0 where MATLAB gives 0.462797.
        int n = in.lambda.length;
        double tot = 0;
        for (int i = 0; i < n; i++) tot += in.lambda[i];
        double tput = 0;
        for (int i = 0; i < n; i++) {
            double[] visits = visitsOf(in.R[i], S);
            double w = tot > 0 ? in.lambda[i] / tot : 0.0;
            double miss = (pmiss != null && i < pmiss.length) ? pmiss[i] : 1.0;
            tput += in.sourceRate * w * miss * visits[s];
        }
        return tput;
    }

    private static double[] visitsOf(double[][] R, int S) {
        double[] a = new double[S];
        double[][] At = new double[S][S];
        for (int s = 0; s < S; s++) { a[s] = R[0][s + 1]; for (int sp = 0; sp < S; sp++) At[sp][s] = (s == sp ? 1.0 : 0.0) - R[s + 1][sp + 1]; }
        return jsolve(At, a);
    }

    private static double[] jsolve(double[][] Ain, double[] bin) {
        int N = bin.length; double[][] A = new double[N][N]; double[] b = bin.clone();
        for (int i = 0; i < N; i++) A[i] = Ain[i].clone();
        for (int col = 0; col < N; col++) {
            int piv = col; double best = Math.abs(A[col][col]);
            for (int rr = col + 1; rr < N; rr++) if (Math.abs(A[rr][col]) > best) { best = Math.abs(A[rr][col]); piv = rr; }
            if (piv != col) { double[] tr = A[piv]; A[piv] = A[col]; A[col] = tr; double tb = b[piv]; b[piv] = b[col]; b[col] = tb; }
            double d = A[col][col];
            for (int rr = col + 1; rr < N; rr++) { double f = A[rr][col] / d; if (f == 0) continue; for (int c = col; c < N; c++) A[rr][c] -= f * A[col][c]; b[rr] -= f * b[col]; }
        }
        double[] x = new double[N];
        for (int rr = N - 1; rr >= 0; rr--) { double sum = b[rr]; for (int c = rr + 1; c < N; c++) sum -= A[rr][c] * x[c]; x[rr] = sum / A[rr][rr]; }
        return x;
    }
}
