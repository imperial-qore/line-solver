/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.analyzers;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import jline.api.cache.Cache_gamma_lp;
import jline.api.cache.Cache_mva;
import jline.api.cache.Cache_prob_erec;
import jline.api.cache.Cache_prob_fpi;
import jline.api.cache.Cache_ttl_lrua;
import jline.api.mc.Dtmc_stochcomp;
import jline.api.sn.SnRefreshVisits;
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

public final class Solver_mva_cacheqn_analyzer {
    private Solver_mva_cacheqn_analyzer() {}

    public static MVAResult solver_mva_cacheqn_analyzer(NetworkStruct sn, SolverOptions options) {
        MVAResult res = new MVAResult();
        NetworkStruct snorig = sn;
        String method = options.method;
        try {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();
            ObjectOutputStream out = new ObjectOutputStream(bos);
            out.writeObject(sn);
            ByteArrayInputStream bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream in = new ObjectInputStream(bis);
            snorig = (NetworkStruct) in.readObject();
        } catch (IOException e) {
            InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                    "Could not create a copy of the NetworkStruct in SolverMVACacheQNAnalyzer");
        } catch (ClassNotFoundException e) {
            InputOutput.line_error(InputOutput.mfilename(new Object() {}),
                    "Could not create a copy of the NetworkStruct in SolverMVACacheQNAnalyzer");
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
        for (Integer ind : statefulNodes) {
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

        // per-cache converged inputs, used to report the per-item occupancy
        // through getAvgItemTable (computed once after convergence)
        Map<Integer, Matrix> cacheGamma = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix> cacheM = new HashMap<Integer, Matrix>();
        Map<Integer, Matrix[]> cacheLambda = new HashMap<Integer, Matrix[]>();
        Map<Integer, Matrix[][]> cacheR = new HashMap<Integer, Matrix[][]>();
        Map<Integer, ReplacementStrategy> cacheStrat = new HashMap<Integer, ReplacementStrategy>();

        List<Integer> inputClass = new ArrayList<Integer>();
        for (int it = 1; it <= options.iter_max; it++) {
            inputClass = new ArrayList<Integer>();
            for (Integer ind : caches) {
                CacheNodeParam ch = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                Matrix hitClass = ch.hitclass;
                Matrix missClass = ch.missclass;
                for (int i = 0; i < hitClass.getNumCols(); i++) {
                    if (hitClass.get(0, i) != -1.0) {
                        inputClass.add(i);
                    }
                }
                if (it == 1) {
                    for (int i = 0; i < inputClass.size(); i++) {
                        lambda_1.set(inputClass.get(i), random.nextDouble());
                    }
                    lambda = new Matrix(lambda_1);
                    sn.nodetype.set(ind, NodeType.ClassSwitch);
                }

                Matrix m = ch.itemcap;
                int n = ch.nitems;
                int h = m.length();
                int u = lambda.length();

                Matrix missrate = new Matrix(sn.nodetype.size(), u);

                Matrix[] lambda_cache = new Matrix[u];
                for (int i = 0; i < u; i++) {
                    lambda_cache[i] = new Matrix(n, h + 1);
                }

                for (int v = 0; v < u; v++) {
                    for (int k = 0; k < n; k++) {
                        for (int l = 0; l < h + 1; l++) {
                            if (ch.pread.getOrDefault(v, null) != null) {
                                lambda_cache[v].set(k, l, lambda.get(v) * ch.pread.get(v).get(k));
                            }
                        }
                    }
                }

                Matrix[][] R = ch.accost;
                if (R == null) {
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
                Matrix gamma = Cache_gamma_lp.cache_gamma_lp(lambda_cache, R).gamma;
                cacheGamma.put(ind, gamma);
                cacheM.put(ind, m);
                cacheLambda.put(ind, lambda_cache);
                cacheR.put(ind, R);
                cacheStrat.put(ind, ((Cache) sn.nodes.get(ind)).getReplacementStrategy());
                Matrix pij;
                if ("exact".equals(method)) {
                    pij = Cache_mva.cache_mva(gamma, m).pij;
                    Matrix newpij = new Matrix(pij.getNumRows(), pij.getNumCols() + 1);
                    for (int i = 0; i < newpij.getNumRows(); i++) {
                        for (int j = 0; j < newpij.getNumCols(); j++) {
                            if (j == 0) {
                                newpij.set(i, j, FastMath.abs(1 - pij.sumRows(i)));
                            } else {
                                newpij.set(i, j, pij.get(i, j - 1));
                            }
                        }
                    }
                    pij = newpij;
                    res.method = "exact";
                } else {
                    pij = Cache_prob_fpi.cache_prob_fpi(gamma, m);
                    res.method = "default".equals(method) ? "default/fpi" : "fpi";
                }
                for (int i = 0; i < missprob.getNumCols(); i++) {
                    missprob.set(ind, i, 0);
                }
                for (int v = 0; v < u; v++) {
                    Matrix A = Matrix.extractColumn(lambda_cache[v], 0, null).transpose();
                    Matrix B = Matrix.extractColumn(pij, 0, null);
                    missrate.set(ind, v, A.mult(B).get(0));
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
                for (int i = 0; i < inputClass.size(); i++) {
                    int r = inputClass.get(i);
                    for (int j = 0; j < sn.rtnodes.getNumCols(); j++) {
                        sn.rtnodes.set(ind * K + r, j, 0);
                    }
                    for (int jnd = 0; jnd < I; jnd++) {
                        if (sn.connmatrix.get(ind, jnd) == 1.0) {
                            sn.rtnodes.set(ind * K + r,
                                    (int) (jnd * K + hitClass.get(r)), hitprob.get(ind, r));
                            sn.rtnodes.set(ind * K + r,
                                    (int) (jnd * K + missClass.get(r)), missprob.get(ind, r));
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

            String cacheMethod = res.method;

            if ("aba.upper".equals(method) || "aba.lower".equals(method)
                    || "bjb.upper".equals(method) || "bjb.lower".equals(method)
                    || "pb.upper".equals(method) || "pb.lower".equals(method)
                    || "gb.upper".equals(method) || "gb.lower".equals(method)
                    || "sb.upper".equals(method) || "sb.lower".equals(method)) {
                res = Solver_mva_bound_analyzer.solver_mva_bound_analyzer(sn, options.copy());
            } else if (!(sn.lldscaling == null || sn.lldscaling.isEmpty())
                    || !(sn.cdscaling == null || sn.cdscaling.isEmpty())
                    || !(sn.jdscaling == null || sn.jdscaling.isEmpty())) {
                res = Solver_mvald_analyzer.solver_mvald_analyzer(sn, options.copy());
            } else {
                res = Solver_mva_analyzer.solver_mva_analyzer(sn, options.copy());
            }

            if (cacheMethod != null) {
                res.method = cacheMethod;
            }

            Matrix nodevisits = null;
            for (Object key : sn.nodevisits.keySet()) {
                if (nodevisits == null) {
                    nodevisits = sn.nodevisits.get(key);
                } else {
                    nodevisits = nodevisits.add(1.0, sn.nodevisits.get(key));
                }
            }
            for (Integer ind : caches) {
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
                    for (Integer ix : inchain) {
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
                res.iter = it;
                break;
            }
            lambda_1 = lambda;
        }
        res.hitProb = hitprob;
        res.missProb = missprob;

        // see _kb/06-solver-catalog.md for rationale
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
                Matrix pe = Cache_prob_erec.cache_prob_erec(gamma, m);
                itemProb = pe;
            }
            res.cacheItemProb.put(ind, itemProb);
        }

        return res;
    }
}
