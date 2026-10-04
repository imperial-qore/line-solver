/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import org.apache.commons.math3.util.FastMath;

import jline.api.cache.Cache_gamma_lp;
import jline.api.cache.Cache_mva;
import jline.api.cache.Cache_prob_fpi;
import jline.api.cache.Cache_ttl_lrua;
import jline.api.cache.Cache_ttl_lrum_map;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * MVA Analyzer class for non-rentrant caches
 */
public final class Solver_mva_cache_analyzer {
    private Solver_mva_cache_analyzer() {}

    public static MVAResult solver_mva_cache_analyzer(NetworkStruct sn, SolverOptions options) {
        MVAResult res = new MVAResult();
        long startTime = System.nanoTime();
        String method = options.method;
        Matrix QN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix UN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix RN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix TN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix CN = new Matrix(1, sn.nclasses);
        Matrix AN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix WN = new Matrix(sn.nnodes, sn.nclasses);
        Matrix XN = new Matrix(1, sn.nclasses);
        double lG = Double.NaN;
        int iter = 1;

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
        int h = m.length(); // Number of lists
        int u = sn.nclasses; // Number of users
        Matrix[] lambda = new Matrix[u];
        for (int v = 0; v < u; v++) {
            lambda[v] = new Matrix(n, h + 1);
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
        Matrix gamma = Cache_gamma_lp.cache_gamma_lp(lambda, R).gamma;
        Matrix pij = null;
        Matrix missRate = new Matrix(1, u);
        if ("exact".equals(method)) {
            ReplacementStrategy rs = cache.getReplacementStrategy();
            if (rs == ReplacementStrategy.RR || rs == ReplacementStrategy.FIFO) {
                Ret.cacheMVA mvaRes = Cache_mva.cache_mva(gamma, m);
                pij = mvaRes.pij;
                Matrix newPij = new Matrix(pij.getNumRows(), 1 + pij.getNumCols());
                int i = 0;
                while (i < pij.getNumRows()) {
                    int j = 0;
                    while (j < newPij.getNumCols()) {
                        if (j < 1) {
                            newPij.set(i, j, FastMath.abs(1 - pij.sumRows(i)));
                        } else {
                            newPij.set(i, j, pij.get(i, j - 1));
                        }
                        j++;
                    }
                    i++;
                }
                pij = newPij;
                missRate = new Matrix(1, u);
                int v = 0;
                while (v < u) {
                    missRate.set(v,
                            Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null))
                                    .get(0));
                    v++;
                }
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "MVA does not support exact solution of the specified cache replacement policy.");
            } else {
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "MVA does not support exact solution of the specified cache replacement policy.");
            }
        } else {
            ReplacementStrategy rs = cache.getReplacementStrategy();
            if (rs == ReplacementStrategy.RR || rs == ReplacementStrategy.FIFO) {
                pij = Cache_prob_fpi.cache_prob_fpi(gamma, m); // FPI method
            } else if (rs == ReplacementStrategy.LRU) {
                // see _kb/06-solver-catalog.md for rationale
                int carrier = -1;
                if (sn.markidx != null && source_ist < sn.markidx.getNumRows()) {
                    for (int r0 = 0; r0 < u; r0++) {
                        if (sn.markidx.get(source_ist, r0) > 0) {
                            carrier = r0;
                            break;
                        }
                    }
                }
                if (carrier >= 0) {
                    MatrixCell Dcell = sn.proc.get(sn.stations.get(source_ist))
                            .get(sn.jobclasses.get(carrier));
                    Matrix D0 = Dcell.get(0);
                    Matrix D1agg = Dcell.get(1);
                    int d = D0.getNumRows();
                    MatrixCell[] D0c = new MatrixCell[n];
                    MatrixCell[] D1c = new MatrixCell[n];
                    boolean allmarked = true;
                    for (int k = 0; k < n; k++) {
                        Matrix D1k = new Matrix(d, d);
                        for (int v = 0; v < u; v++) {
                            if (ch.pread != null && ch.pread.get(v) != null) {
                                double pk = ch.pread.get(v).get(k);
                                if (sn.markidx.get(source_ist, v) > 0) {
                                    int mk = (int) sn.markidx.get(source_ist, v);
                                    D1k = D1k.add(pk, Dcell.get(1 + mk));
                                } else if (sourceRate.get(v) > 0 && pk > 0) {
                                    allmarked = false; // unmarked reader mixed in
                                }
                            }
                        }
                        Matrix D0k = D0.add(1.0, D1agg).add(-1.0, D1k);
                        D0c[k] = new MatrixCell(h);
                        D1c[k] = new MatrixCell(h);
                        for (int l = 0; l < h; l++) {
                            D0c[k].set(l, D0k);
                            D1c[k].set(l, D1k);
                        }
                    }
                    if (allmarked) {
                        pij = Cache_ttl_lrum_map.cache_ttl_lrum_map_pij(D0c, D1c, m);
                    }
                }
                if (pij == null) {
                    // TTL tree approximation (allows trees and access costs)
                    pij = Cache_ttl_lrua.cache_ttl_lrua(lambda, R, m);
                }
            } else if (rs == ReplacementStrategy.HLRU) {
                // h-LRU / LRU(m) characteristic-time approximation (linear
                // list topology; mirrors MATLAB solver_mva_cache_analyzer)
                pij = jline.api.cache.Cache_ttl_hlru.cache_ttl_hlru(lambda, m);
            } else {
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "MVA does not support approximate solution of the specified cache replacement policy.");
            }
            missRate = new Matrix(1, u);
            int v = 0;
            while (v < u) {
                missRate.set(v,
                        Matrix.extractColumn(lambda[v], 0, null).transpose().mult(Matrix.extractColumn(pij, 0, null))
                                .get(0));
                v++;
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        Matrix itemProb;
        ReplacementStrategy rsItem = cache.getReplacementStrategy();
        if (!"exact".equals(method) && (rsItem == ReplacementStrategy.RR || rsItem == ReplacementStrategy.FIFO)) {
            if (n > 10) {
                InputOutput.line_warning(InputOutput.mfilename(new Object()),
                        "Per-item cache occupancy (getAvgItemTable) requires the exact algorithm for RR/FIFO and is skipped for caches with more than 10 items (%d items); reporting NaN.", n);
                itemProb = new Matrix(n, h + 1);
                itemProb.fill(Double.NaN);
            } else {
                Matrix pe = Cache_mva.cache_mva(gamma, m).pij;
                itemProb = new Matrix(pe.getNumRows(), pe.getNumCols() + 1);
                for (int i = 0; i < pe.getNumRows(); i++) {
                    itemProb.set(i, 0, FastMath.abs(1 - pe.sumRows(i)));
                    for (int j = 0; j < pe.getNumCols(); j++) {
                        itemProb.set(i, j + 1, pe.get(i, j));
                    }
                }
            }
        } else {
            itemProb = pij;
        }
        if (cache != null && itemProb != null && itemProb.getNumCols() == h + 1) {
            cache.setResultItemProb(itemProb);
        }
        if (!"ttl.map".equals(method)) {
            for (int r = 0; r < sn.nclasses; r++) {
                if (ch.hitclass.length() > r && ch.missclass.get(r) > -1 && ch.hitclass.get(r) > -1) {
                    XN.set((int) ch.missclass.get(r), XN.get((int) ch.missclass.get(r)) + missRate.get(r));
                    XN.set((int) ch.hitclass.get(r),
                            XN.get((int) ch.hitclass.get(r)) + sourceRate.get(r) - missRate.get(r));
                }
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
        res.logNormConstAggr = lG;
        res.runtime = (endTime - startTime) / 1000000000.0;
        res.iter = iter;

        // Set the actual method used
        if ("exact".equals(method)) {
            res.method = "exact";
        } else {
            ReplacementStrategy rs = cache.getReplacementStrategy();
            if (rs == ReplacementStrategy.RR || rs == ReplacementStrategy.FIFO) {
                res.method = "fpi";
            } else if (rs == ReplacementStrategy.LRU || rs == ReplacementStrategy.HLRU) {
                res.method = "ttl";
            } else {
                res.method = method;
            }
        }

        return res;
    }
}
