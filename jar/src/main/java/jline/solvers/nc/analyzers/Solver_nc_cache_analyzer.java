/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import java.util.List;

import jline.api.cache.Cache_cost;
import jline.api.cache.Cache_cost_pathcheck;
import jline.api.cache.Cache_gamma_lp;
import jline.api.cache.Cache_miss_is;
import jline.api.cache.Cache_miss_spm;
import jline.api.cache.Cache_prob_erec;
import jline.api.cache.Cache_spm_size;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_nc_cache_analyzer {
    private Solver_nc_cache_analyzer() {}

    public static NCResult solver_nc_cache_analyzer(NetworkStruct sn, SolverOptions options) {
        long startTime = System.nanoTime();
        Matrix QN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix UN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix RN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix TN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix CN = new Matrix(1, sn.nclasses);
        Matrix AN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix WN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix XN = new Matrix(1, sn.nclasses);
        double lG = Double.NaN;

        int source_ist = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                source_ist = (int) sn.nodeToStation.get(i);
                break;
            }
        }
        Matrix sourceRate = new Matrix(1, sn.rates.getNumCols());
        for (int i = 0; i < sn.rates.getNumCols(); i++) {
            if (!Double.isNaN(sn.rates.get(source_ist, i))) {
                sourceRate.set(i, sn.rates.get(source_ist, i));
            }
        }
        TN = new Matrix(sn.nodetype.size(), sn.nclasses);
        for (int i = 0; i < sourceRate.getNumCols(); i++) {
            TN.set(source_ist, i, sourceRate.get(i));
        }

        Cache cache = null;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Cache) {
                cache = (Cache) sn.nodes.get(i);
                break;
            }
        }
        CacheNodeParam ch = (CacheNodeParam) sn.nodeparam.get(cache);

        Matrix m = ch.itemcap;
        int n = ch.nitems; // Number of items

        // MATLAB carries a commented-out n<m+2 guard here (never enabled); not ported
        int h = m.length();          // Number of lists
        int u = sn.nclasses;         // Number of users
        Matrix[] lambda = new Matrix[u];
        for (int i = 0; i < u; i++) {
            lambda[i] = new Matrix(n, h + 1);
        }

        for (int v = 0; v < u; v++) {
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < h + 1; l++) {
                    if (ch.pread != null && ch.pread.get(v) != null) {
                        lambda[v].set(k, l, sourceRate.get(v) * ch.pread.get(v).get(k));
                    }
                }
            }
        }

        MatrixCell lambdaCell = new MatrixCell(lambda);

        Matrix[][] R = ch.accost;
        if (R == null) {
            // Default linear cache routing: items flow from list l to list l+1
            R = new Matrix[u][n];
            for (int v = 0; v < u; v++) {
                for (int k = 0; k < n; k++) {
                    R[v][k] = new Matrix(h + 1, h + 1);
                    for (int l = 0; l < h; l++) {
                        R[v][k].set(l, l + 1, 1.0);
                    }
                    R[v][k].set(h, h, 1.0);
                }
            }
        }
        Ret.cacheGamma gammaRet = Cache_gamma_lp.cache_gamma_lp(lambda, R);
        Matrix gamma = gammaRet.gamma;

        // per-item storage costs and per-list cost caps (ton21cache Sec. IX)
        Matrix sigma = ch.itemsize;
        Matrix costcap = ch.costcap;
        if (costcap != null && !costcap.isEmpty()) {
            if (sigma == null || sigma.isEmpty()) {
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "Storage cost caps require per-item sizes; call Cache.setItemSizes first.");
            }
            if (sigma.length() != n) {
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "The item size vector must have one entry per item.");
            }
            if (costcap.length() != h) {
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "The cost cap vector must have one entry per cache list.");
            }
            List<Cache_cost_pathcheck.BlockedPair> viol =
                    Cache_cost_pathcheck.cache_cost_pathcheck(gamma, sigma, costcap, gammaRet.parent);
            if (!viol.isEmpty()) {
                Cache_cost_pathcheck.BlockedPair first = viol.get(0);
                InputOutput.line_warning(InputOutput.mfilename(new Object()),
                        "Storage cost caps block the promotion path of item %d into list %d at list %d (and %d further pairs). The exact recursion normalizes over all size-feasible states, which is then a strict superset of the states the cache can reach; cross-check with SolverLDES.",
                        first.item + 1, first.list + 1, first.blockingList + 1, viol.size() - 1);
            }
        }

        Matrix pij = null;
        Matrix missRate = new Matrix(1, u);

        String cacheMethod = options.method;
        // 'rayint' is an alias of 'spm': on a cache both name the SPM saddle point,
        // and the method name stays live for Solver_nc_retrieval_analyzer's delayed-hit
        // expansion.
        if ("rayint".equals(cacheMethod) || "spm".equals(cacheMethod)) {
            cacheMethod = "default";
        }
        // The SPM family serves its size-tilted form (Cache_spm_size) once the items
        // carry storage costs. The saddle escapes to infinity at sum(m) = n, so the
        // size-free saddle point takes over there rather than the exact recursion,
        // which would refuse every replacement policy outside RR/FIFO.
        double msum = 0;
        for (int j = 0; j < h; j++) {
            msum += m.get(j);
        }
        boolean useSpmSize = "default".equals(cacheMethod)
                && sigma != null && !sigma.isEmpty() && msum < n;
        if (costcap != null && !costcap.isEmpty() && !useSpmSize
                && !"exact".equals(cacheMethod) && !"sampling".equals(cacheMethod)) {
            // The size-free SPM and the mean-field methods have no cost-capped
            // counterpart: the k - sigma_i 1_j argument couples item sizes into the
            // recursion graph, which only Cache_spm_size and the exact path carry.
            double lattice = n;
            for (int j = 0; j < h; j++) {
                lattice *= (m.get(j) + 1) * (costcap.get(j) + 1);
            }
            cacheMethod = (lattice <= 1e6) ? "exact" : "sampling";
            InputOutput.line_warning(InputOutput.mfilename(new Object()),
                    "Method '%s' does not support storage cost caps; using '%s' instead.", options.method, cacheMethod);
        }

        NCResult res = new NCResult();
        if ("exact".equals(cacheMethod)) {
            res.method = "exact";
            ReplacementStrategy rs = cache.getReplacementStrategy();
            if (rs == ReplacementStrategy.RR || rs == ReplacementStrategy.FIFO) {
                // cache_prob_erec already returns [n x (h+1)] with col 0 = miss
                pij = Cache_prob_erec.cache_prob_erec(gamma, m, sigma, costcap);
                int v = 0;
                while (v < u) {
                    missRate.set(v,
                            Matrix.extractColumn(lambda[v], 0, null).transpose()
                                    .mult(Matrix.extractColumn(pij, 0, null)).get(0));
                    v++;
                }
            } else {
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "NC does not support exact solution of the specified cache replacement policy.");
            }
        } else if (useSpmSize) {
            // Size-tilted SPM: a 2h Newton solve whose cost does not grow with the
            // (m,k) lattice the exact recursion walks. O(1/n), so it wants room
            // between the occupancies and n.
            res.method = "spm.size";
            double[] raysigma = new double[n];
            for (int i = 0; i < n; i++) {
                raysigma[i] = sigma.get(i);
            }
            double[] raycap = new double[h];
            if (costcap != null && !costcap.isEmpty()) {
                for (int j = 0; j < h; j++) {
                    raycap[j] = costcap.get(j);
                }
            } else {
                // Sizes but no caps: cap each list at the dearest load it can hold,
                // which is exactly slack, so the cost coordinate leaves the saddle
                // and the expansion degenerates to the size-free one.
                double[] srt = raysigma.clone();
                java.util.Arrays.sort(srt);
                for (int j = 0; j < h; j++) {
                    double top = 0;
                    for (int a = 0; a < (int) m.get(j); a++) {
                        top += srt[n - 1 - a];
                    }
                    raycap[j] = top;
                }
            }
            double[][] graw = new double[n][h];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < h; j++) {
                    graw[i][j] = gamma.get(i, j);
                }
            }
            double[] mraw = new double[h];
            for (int j = 0; j < h; j++) {
                mraw[j] = m.get(j);
            }
            Cache_spm_size.Result ray =
                    Cache_spm_size.cache_spm_size(graw, mraw, raysigma, raycap);
            pij = new Matrix(n, h + 1);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= h; j++) {
                    pij.set(i, j, ray.pij[i][j]);
                }
            }
            int v = 0;
            while (v < u) {
                missRate.set(v,
                        Matrix.extractColumn(lambda[v], 0, null).transpose()
                                .mult(Matrix.extractColumn(pij, 0, null)).get(0));
                v++;
            }
        } else if ("sampling".equals(cacheMethod)) {
            // importance sampling over the exchangeable-family measure, any replacement policy: see _kb/09-ldes-and-cache.md
            res.method = "sampling";
            int samples = options.samples;
            Ret.cacheMissSpm cacheMissIs = Cache_miss_is.cache_miss_is(gamma, m, lambdaCell, samples, sigma, costcap);
            double[] missRateList = cacheMissIs.getMU();
            int v = 0;
            while (v < u) {
                missRate.set(v, missRateList[v]);
                v++;
            }
        } else {
            // Size-free SPM expansion, the default/spm/rayint branch with no item
            // sizes: see _kb/09-ldes-and-cache.md ("Solver_nc_cache_analyzer's
            // non-exact methods share one SPM expansion")
            Ret.cacheMissSpm cacheMissSpm = Cache_miss_spm.cache_miss_spm(gamma, m, lambdaCell);
            double[] missRateList = cacheMissSpm.getMU();
            int v = 0;
            while (v < u) {
                missRate.set(v, missRateList[v]);
                v++;
            }
            cacheMissSpm.getlE(); // never used
            res.method = "spm";
        }

        // per-item occupancy [n x (h+1)] (col 0 = miss, cols 1.. = per-list)
        // via the exact recursion (cache_prob_erec). The spm/is algorithms
        // approximate each per-list column independently and do not form a
        // proper distribution for h>1, so the per-item table is always taken
        // from the exact product-form algorithm. The exact recursion is only
        // tractable for small item sets, so it is skipped (NaN, with a warning)
        // for caches with more than 10 items.
        if (cache != null) {
            Matrix itemProb;
            if (n > 10) {
                InputOutput.line_warning(InputOutput.mfilename(new Object()),
                        "Per-item cache occupancy (getAvgItemTable) requires the exact algorithm and is skipped for caches with more than 10 items (%d items); reporting NaN.", n);
                itemProb = new Matrix(n, m.length() + 1);
                itemProb.fill(Double.NaN);
            } else {
                itemProb = Cache_prob_erec.cache_prob_erec(gamma, m, sigma, costcap);
            }
            cache.setResultItemProb(itemProb);
            // mean storage cost held by each list, K_j = sum_i sigma_i pi_ij
            if (sigma != null && !sigma.isEmpty()
                    && itemProb.getNumRows() == n && itemProb.getNumCols() == h + 1) {
                Matrix listCost = Cache_cost.cache_cost(gamma, m, sigma, costcap, itemProb);
                cache.setResultListCost(listCost);
                ch.actuallistcost = listCost;
            }
        }

        for (int r = 0; r < sn.nclasses; r++) {
            if (ch.hitclass.length() > r && ch.missclass.get(r) > -1 && ch.hitclass.get(r) > -1) {
                XN.set((int) ch.missclass.get(r), XN.get((int) ch.missclass.get(r)) + missRate.get(r));
                XN.set((int) ch.hitclass.get(r), XN.get((int) ch.hitclass.get(r)) + sourceRate.get(r) - missRate.get(r));
            }
        }

        long endTime = System.nanoTime();
        res.QN = QN;
        res.UN = UN;
        res.RN = RN;
        res.TN = TN;
        res.CN = CN;
        res.XN = XN;
        res.AN = AN;
        res.WN = WN;
        res.lG = lG; // ???
        res.runtime = (endTime - startTime) / 1000000000.0;
        // res.iter = iter;
        return res;
    }
}
