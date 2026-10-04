/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import java.util.HashMap;
import java.util.Map;

import jline.api.cache.Cache_gamma_lp;
import jline.api.cache.Cache_miss_spm;
import jline.api.cache.Cache_prob_erec;
import jline.api.cache.Cache_ttl_lrua;
import jline.api.mc.Dtmc_stochcomp;
import jline.api.sn.SnRefreshVisits;
import jline.io.InputOutput;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.solvers.SolverOptions;
import jline.solvers.nc.NCResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * NC Analyzer class for solver_nc_cacheqn_analyzer.
 */
public final class Solver_nc_cache_qn_analyzer {
    private Solver_nc_cache_qn_analyzer() {}

    public static NCResult solver_nc_cache_qn_analyzer(NetworkStruct sn, SolverOptions options) {
        NetworkStruct snorig = sn;
        NCResult res = new NCResult();
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(sn);
            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream in = new ObjectInputStream(bis);
            snorig = (NetworkStruct) in.readObject();
        } catch (IOException e) {
            InputOutput.line_error(InputOutput.mfilename(new Object()),
                    "Could not create a copy of the NetworkStruct in SolverNCCacheQNAnalyzer");
        } catch (ClassNotFoundException e) {
            InputOutput.line_error(InputOutput.mfilename(new Object()),
                    "Could not create a copy of the NetworkStruct in SolverNCCacheQNAnalyzer");
        }
        sn = snorig;
        Random random = new Random((long) options.seed);
        int I = sn.nnodes;
        int K = sn.nclasses;
        List<Integer> statefulNodes = new ArrayList<Integer>();
        for (int i = 0; i < sn.isstateful.length(); i++) {
            if (sn.isstateful.get(i) == 1.0) {
                statefulNodes.add(i);
            }
        }
        Matrix statefulNodeClasses = new Matrix(1, statefulNodes.size() * K);
        int idx = 0;
        for (int ind : statefulNodes) {
            for (int i = 0; i < K; i++) {
                statefulNodeClasses.set(idx, (double) (ind * K + i));
                idx++;
            }
        }
        Matrix lambda = new Matrix(1, K);
        Matrix lambda_1 = new Matrix(1, K);
        List<Integer> caches = new ArrayList<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Cache) {
                caches.add(i);
            }
        }
        Matrix hitprob = new Matrix(sn.nodetype.size(), K);
        Matrix missprob = new Matrix(sn.nodetype.size(), K);
        // converged isolated-cache inputs, kept for the per-item occupancy law
        Map<Integer, Matrix> cacheGamma = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> cacheM = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix[]> cacheLambda = new HashMap<Integer, Matrix[]>();
        Map<Integer, Matrix[][]> cacheR = new HashMap<Integer, Matrix[][]>();
        Map<Integer, ReplacementStrategy> cacheStrat = new HashMap<Integer, ReplacementStrategy>();

        for (int it = 1; it <= options.iter_max; it++) {
            List<Integer> inputClass = new ArrayList<Integer>();
            for (int ind : caches) {
                CacheNodeParam ch = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                Matrix hitClass = ch.hitclass;
                Matrix missClass = ch.missclass;
                for (int i = 0; i < hitClass.length(); i++) {
                    if (hitClass.get(i) != -1.0) {
                        inputClass.add(i);
                    }
                }

                // Solution of the isolated cache
                Matrix m = ch.itemcap;
                int n = ch.nitems;
                int h = m.length();
                int u = lambda.length();

                if (it == 1) {
                    if (n < m.elementSum() + 2) {
                        InputOutput.line_error(InputOutput.mfilename(new Object()),
                                "NC requires the number of items to exceed the cache capacity at least by 2.");
                    }
                    for (int i = 0; i < inputClass.size(); i++) {
                        lambda_1.set(inputClass.get(i), random.nextDouble());
                    }
                    lambda = new Matrix(lambda_1);
                    sn.nodetype.set(ind, NodeType.ClassSwitch);
                }

                Matrix missrate = new Matrix(sn.nodetype.size(), u);

                Matrix[] lambda_cache = new Matrix[u];
                for (int i = 0; i < u; i++) {
                    lambda_cache[i] = new Matrix(n, h + 1);
                }

                for (int v = 0; v < u; v++) {
                    for (int kk = 0; kk < n; kk++) {
                        for (int l = 0; l < h + 1; l++) {
                            if (ch.pread != null && ch.pread.get(v) != null) {
                                lambda_cache[v].set(kk, l, lambda.get(v) * ch.pread.get(v).get(kk));
                            }
                        }
                    }
                }

                Matrix[][] R = ch.accost;
                if (R == null) {
                    // Default linear cache routing: items flow from list l to list l+1
                    R = new Matrix[u][n];
                    for (int v = 0; v < u; v++) {
                        for (int kk = 0; kk < n; kk++) {
                            R[v][kk] = new Matrix(h + 1, h + 1);
                            for (int l = 0; l < h; l++) {
                                R[v][kk].set(l, l + 1, 1.0);
                            }
                            R[v][kk].set(h, h, 1.0);
                        }
                    }
                }
                Matrix gamma = Cache_gamma_lp.cache_gamma_lp(lambda_cache, R).gamma;
                cacheGamma.put(ind, gamma);
                cacheM.put(ind, m);
                cacheLambda.put(ind, lambda_cache);
                cacheR.put(ind, R);
                cacheStrat.put(ind, ch.replacestrat);

                if ("exact".equals(options.method)) {
                    Matrix pij = Cache_prob_erec.cache_prob_erec(gamma, m);
                    res.method = "exact";
                    for (int i = 0; i < missprob.getNumCols(); i++) {
                        missprob.set(ind, i, 0);
                    }
                    for (int v = 0; v < u; v++) {
                        Matrix A = Matrix.extractColumn(lambda_cache[v], 0, null).transpose();
                        Matrix B = Matrix.extractColumn(pij, 0, null);
                        missrate.set(ind, v, A.mult(B).get(0));
                    }
                } else {
                    // Use cache_miss_spm for non-exact method (consistent with MATLAB)
                    MatrixCell lambdaCell = new MatrixCell(lambda_cache);
                    Ret.cacheMissSpm cacheMissResult = Cache_miss_spm.cache_miss_spm(gamma, m, lambdaCell);
                    res.method = "spm";
                    for (int i = 0; i < missprob.getNumCols(); i++) {
                        missprob.set(ind, i, 0);
                    }
                    // Extract miss rates directly from MU field
                    for (int v = 0; v < u; v++) {
                        missrate.set(ind, v, cacheMissResult.getMU()[v]);
                    }
                }
                for (int i = 0; i < missprob.getNumCols(); i++) {
                    double missValue = missrate.get(ind, i) / lambda.get(i);
                    if (Double.isNaN(missValue)) {
                        missprob.set(ind, i, 0);
                    } else {
                        missprob.set(ind, i, missValue);
                    }
                }
                for (int i = 0; i < hitprob.getNumCols(); i++) {
                    double hitValue = 1 - missrate.get(ind, i) / lambda.get(i);
                    if (Double.isNaN(hitValue)) {
                        hitprob.set(ind, i, 0);
                    } else {
                        hitprob.set(ind, i, hitValue);
                    }
                }
                // Bring back isolated model results into the queueing model
                for (int i = 0; i < inputClass.size(); i++) {
                    int r = inputClass.get(i);
                    for (int j = 0; j < sn.rtnodes.getNumCols(); j++) {
                        sn.rtnodes.set(ind * K + r, j, 0);
                    }
                    for (int jnd = 0; jnd < I; jnd++) {
                        if (sn.connmatrix.get(ind, jnd) == 1.0) {
                            sn.rtnodes.set(ind * K + r, (int) (jnd * K + hitClass.get(r)), hitprob.get(ind, r));
                            sn.rtnodes.set(ind * K + r, (int) (jnd * K + missClass.get(r)), missprob.get(ind, r));
                        }
                    }
                }
                List<Integer> statefulNodeClassesList = new ArrayList<Integer>(statefulNodeClasses.length());
                for (int i = 0; i < statefulNodeClasses.length(); i++) {
                    statefulNodeClassesList.add((int) statefulNodeClasses.get(i));
                }
                sn.rt = Dtmc_stochcomp.dtmc_stochcomp(sn.rtnodes, statefulNodeClassesList);
            }
            sn = SnRefreshVisits.snRefreshVisits(sn, sn.chains, sn.rt, sn.rtnodes);

            // The network solve below REPLACES res, so the isolated-cache method
            // decided above is lost with it unless it is carried across. Without
            // this the banner reported the queueing analyzer's method ("exact")
            // for a run whose cache was solved by the spm approximation.
            String cacheMethod = res.method;

            // Default branch: choose ld vs nc based on scaling
            if (!(sn.lldscaling == null || sn.lldscaling.isEmpty())
                    || !(sn.cdscaling == null || sn.cdscaling.isEmpty())
                    || !(sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                res = Solver_ncld_analyzer.solver_ncld_analyzer(sn, options);
            } else {
                res = Solver_nc_analyzer.solver_nc_analyzer(sn, options);
            }

            if (cacheMethod != null) {
                res.method = cacheMethod;
            }

            Matrix nodevisits = null;
            for (Integer key : sn.nodevisits.keySet()) {
                if (nodevisits == null) {
                    nodevisits = sn.nodevisits.get(key);
                } else {
                    nodevisits = nodevisits.add(1.0, sn.nodevisits.get(key));
                }
            }
            for (int ind : caches) {
                for (int i = 0; i < inputClass.size(); i++) {
                    int r = inputClass.get(i);
                    int c = 0;
                    while (c < sn.chains.getNumRows() && sn.chains.get(c, r) == 0.0) {
                        c++;
                    }
                    List<Integer> inchain = new ArrayList<Integer>();
                    for (int j = 0; j < sn.chains.getNumCols(); j++) {
                        if (sn.chains.get(c, j) != 0.0) {
                            inchain.add(j);
                        }
                    }
                    double sumXN = 0.0;
                    for (int ix : inchain) {
                        sumXN += res.XN.get(ix);
                    }
                    if (sn.refclass.get(c) > -1) {
                        lambda.set(r,
                                sumXN * nodevisits.get(ind, r)
                                        / nodevisits.get((int) sn.stationToNode.get((int) sn.refstat.get(r)),
                                                (int) sn.refclass.get(c)));
                    } else {
                        lambda.set(r,
                                sumXN * nodevisits.get(ind, r)
                                        / nodevisits.get((int) sn.stationToNode.get((int) sn.refstat.get(r)), r));
                    }
                }
            }
            if (lambda.sub(lambda_1).elementSum() < options.iter_tol) {
                res.it = it;
                break;
            }
            lambda_1 = lambda;
        }
        res.hitProb = hitprob;
        res.missProb = missprob;

        // Per-item occupancy from the converged access factors, as SolverMVA reports it.
        for (Integer ind : caches) {
            Matrix gamma = cacheGamma.get(ind);
            if (gamma == null) {
                continue;
            }
            Matrix m = cacheM.get(ind);
            int ni = gamma.getNumRows();
            int hi = m.length();
            Matrix itemProb;
            if (cacheStrat.get(ind) == ReplacementStrategy.LRU) {
                itemProb = Cache_ttl_lrua.cache_ttl_lrua(cacheLambda.get(ind), cacheR.get(ind), m);
            } else if (ni > 10) {
                InputOutput.line_warning(InputOutput.mfilename(new Object() {}),
                        "Per-item cache occupancy (getAvgItemTable) requires the exact algorithm for RR/FIFO and is skipped for caches with more than 10 items (%d items); reporting NaN.", ni);
                itemProb = new Matrix(ni, hi + 1);
                itemProb.fill(Double.NaN);
            } else {
                itemProb = Cache_prob_erec.cache_prob_erec(gamma, m);
            }
            res.cacheItemProb.put(ind, itemProb);
        }
        return res;
    }
}
