/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.analyzers;

import jline.api.cache.Cache_gamma_lp;
import jline.api.cache.Cache_miss_rmf;
import jline.api.cache.Cache_miss_fifo_rmf;
import jline.api.cache.Cache_miss_sfifo_rmf;
import jline.api.cache.CacheMissRmfTranResult;
import jline.lang.constant.ReplacementStrategy;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Transient refined-mean-field cache trajectory for an (open) integrated
 * cache-queueing network. Java analog of the MATLAB
 * {@code solver_fld_cacheqn_tran.m}.
 *
 * <p>For the open-cache topology used by the random-environment mean-field
 * cache aggregation (Source -&gt; Cache -&gt; Sink), each read class arrives at
 * the cache at its known source rate, so no arrival-rate fixed point is needed:
 * the per-cache isolated inputs (gamma, m, lambda_cache) are built directly from
 * the stage struct, and the RMF drift is integrated over the window from an
 * optional per-cache initial occupancy (used to carry the mean occupancy across
 * environment switches).</p>
 */
public final class FluidCacheTran {
    private FluidCacheTran() {}

    /**
     * Integrate the transient RMF cache trajectory for all cache nodes.
     *
     * @param sn             stage NetworkStruct (must contain at least one Cache).
     * @param time           end time of the transient window (start is 0).
     * @param nPoints        number of uniform grid points (>= 2).
     * @param xinitPerCache  optional per-cache initial occupancy; entry may be
     *                       null for the default initial state, and the array
     *                       itself may be null.
     * @return {@link CacheTranResult}.
     */
    public static CacheTranResult tran(NetworkStruct sn, double time, int nPoints, double[][] xinitPerCache) {
        int K = sn.nclasses;

        List<Integer> caches = new ArrayList<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Cache) {
                caches.add(i);
            }
        }
        int ncaches = caches.size();
        int[] cacheArr = new int[ncaches];

        double[] tOut = null;
        double[][][] hitprob = new double[ncaches][][];
        double[][][] missprob = new double[ncaches][][];
        double[][] arate = new double[ncaches][K];
        double[][][] xocc = new double[ncaches][][];

        for (int cIdx = 0; cIdx < ncaches; cIdx++) {
            int ind = caches.get(cIdx);
            cacheArr[cIdx] = ind;
            NodeParam param = sn.nodeparam.get(sn.nodes.get(ind));
            if (param == null || !(param instanceof CacheNodeParam)) {
                continue;
            }
            CacheNodeParam ch = (CacheNodeParam) param;
            Matrix hitClass = ch.hitclass;
            Matrix m = ch.itemcap;
            int n = ch.nitems;
            int h = m.length();

            // Input (read) classes: those with a defined hit class.
            List<Integer> inputClass = new ArrayList<Integer>();
            for (int i = 0; i < hitClass.getNumCols(); i++) {
                if (hitClass.get(0, i) != -1.0) {
                    inputClass.add(i);
                }
            }

            // Open-cache arrival rate per read class = its source (reference
            // station) arrival rate.
            Matrix lambda = new Matrix(1, K);
            for (int r : inputClass) {
                int rs = (int) sn.refstat.get(r);
                double arr = rs >= 0 ? sn.rates.get(rs, r) : 0.0;
                if (!Double.isFinite(arr)) {
                    arr = 0.0;
                }
                lambda.set(0, r, arr);
            }

            // Build per-user per-item per-list arrival rates.
            Matrix[] lambda_cache = new Matrix[K];
            for (int v = 0; v < K; v++) {
                lambda_cache[v] = new Matrix(n, h + 1);
                List<Double> pread_v = ch.pread.get(v);
                if (pread_v != null) {
                    for (int k = 0; k < n; k++) {
                        for (int l = 0; l < h + 1; l++) {
                            if (k < pread_v.size()) {
                                lambda_cache[v].set(k, l, lambda.get(0, v) * pread_v.get(k));
                            }
                        }
                    }
                }
            }

            // Access cost matrix R (default linear list-to-list flow).
            Matrix[][] R = ch.accost;
            if (R == null) {
                R = new Matrix[K][n];
                for (int v = 0; v < K; v++) {
                    for (int k = 0; k < n; k++) {
                        R[v][k] = new Matrix(h + 1, h + 1);
                        for (int l = 0; l < h; l++) {
                            R[v][k].set(l, l + 1, 1.0);
                        }
                        R[v][k].set(h, h, 1.0);
                    }
                }
            }

            Ret.cacheGamma gammaResult = Cache_gamma_lp.cache_gamma_lp(lambda_cache, R);
            Matrix gamma = gammaResult.gamma;

            double[] xinit = (xinitPerCache != null && cIdx < xinitPerCache.length) ? xinitPerCache[cIdx] : null;
            // RANDOM(m) and FIFO(m) share the steady state (Gast15 Thm 1) but NOT
            // the transient, so FIFO uses its own position-resolved drift; strict
            // FIFO(m) likewise. LRU/HLRU/CLIMB/QLRU have no drift-based transient.
            CacheMissRmfTranResult tr;
            if (ch.replacestrat == ReplacementStrategy.FIFO) {
                tr = Cache_miss_fifo_rmf.cache_miss_fifo_rmf_tran(gamma, m, lambda_cache, time, nPoints, xinit, ch.accost);
            } else if (ch.replacestrat == ReplacementStrategy.SFIFO) {
                tr = Cache_miss_sfifo_rmf.cache_miss_sfifo_rmf_tran(gamma, m, lambda_cache, time, nPoints, xinit, ch.accost);
            } else if (ch.replacestrat == ReplacementStrategy.RR) {
                tr = Cache_miss_rmf.cache_miss_rmf_tran(gamma, m, lambda_cache, time, nPoints, xinit);
            } else {
                throw new RuntimeException("Transient cache analysis is only available for "
                        + "RANDOM(m)/FIFO(m) and strict FIFO(m) replacement via a drift-based "
                        + "mean field; replacement strategy " + ch.replacestrat + " has none.");
            }
            if (tOut == null) {
                tOut = tr.t;
            }
            int nt = tr.t.length;

            for (int v = 0; v < K; v++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) {
                    double val = lambda_cache[v].get(i, 0);
                    if (Double.isFinite(val)) {
                        s += val;
                    }
                }
                arate[cIdx][v] = s;
            }

            hitprob[cIdx] = new double[K][nt];
            missprob[cIdx] = new double[K][nt];
            for (int v = 0; v < K; v++) {
                double a = arate[cIdx][v];
                if (a > 0) {
                    for (int t = 0; t < nt; t++) {
                        double mp = tr.MU_t[v][t] / a;
                        mp = Math.min(1.0, Math.max(0.0, mp));
                        missprob[cIdx][v][t] = mp;
                        hitprob[cIdx][v][t] = 1.0 - mp;
                    }
                }
            }
            xocc[cIdx] = tr.xtraj;
        }

        return new CacheTranResult(tOut, cacheArr, hitprob, missprob, arate, xocc);
    }
}
