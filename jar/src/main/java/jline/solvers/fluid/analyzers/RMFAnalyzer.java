/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.analyzers;

import jline.api.cache.Cache_gamma_lp;
import jline.api.cache.Cache_miss_rmf;
import jline.api.cache.Cache_miss_sfifo_rmf;
import jline.api.cache.Cache_miss_fifo_rmf;
import jline.api.mc.Dtmc_stochcomp;
import jline.api.sn.SnRefreshVisits;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.lang.constant.NodeType;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.FluidResult;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Integrated cache-queueing network analyzer using fluid approximation.
 *
 * <p>Iterates between isolated cache analysis and fluid ODE solution
 * of the surrounding queueing network until arrival rates converge.
 * This matches the MATLAB solver_fld_cacheqn_analyzer implementation.</p>
 *
 * @see CacheNodeParam
 * @see MatrixMethodAnalyzer
 */
public class RMFAnalyzer implements FluidAnalyzer {

    private Matrix xvecIt;

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        int I = sn.nnodes;
        int K = sn.nclasses;
        int M = sn.nstations;

        // Build statefulNodesClasses indices
        List<Integer> statefulNodes = new ArrayList<Integer>();
        for (int i = 0; i < sn.isstateful.length(); i++) {
            if (sn.isstateful.get(i) == 1.0) {
                statefulNodes.add(i);
            }
        }
        List<Integer> statefulNodeClassesList = new ArrayList<Integer>();
        for (int ind : statefulNodes) {
            for (int k = 0; k < K; k++) {
                statefulNodeClassesList.add(ind * K + k);
            }
        }

        // Find cache nodes
        List<Integer> caches = new ArrayList<Integer>();
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Cache) {
                caches.add(i);
            }
        }

        Matrix lambda = new Matrix(1, K);
        Matrix lambda_1 = new Matrix(1, K);
        Matrix hitprob = new Matrix(caches.size(), K);
        Matrix missprob = new Matrix(caches.size(), K);
        Matrix missrate = new Matrix(caches.size(), K);

        Random random = new Random(options.seed);
        int convergedIter = options.iter_max;

        for (int it = 1; it <= options.iter_max; it++) {
            for (int cIdx = 0; cIdx < caches.size(); cIdx++) {
                int ind = caches.get(cIdx);
                NodeParam param = sn.nodeparam.get(sn.nodes.get(ind));
                if (param == null || !(param instanceof CacheNodeParam)) {
                    continue;
                }
                CacheNodeParam ch = (CacheNodeParam) param;
                Matrix hitClass = ch.hitclass;
                Matrix missClass = ch.missclass;

                // Find input classes (non-negative hitClass entries)
                List<Integer> inputClass = new ArrayList<Integer>();
                for (int i = 0; i < hitClass.getNumCols(); i++) {
                    if (hitClass.get(0, i) != -1.0) {
                        inputClass.add(i);
                    }
                }

                Matrix m = ch.itemcap;
                int n = ch.nitems;
                int h = m.length();
                int u = lambda.getNumCols();

                if (it == 1) {
                    // Initial random arrival rates
                    for (int i : inputClass) {
                        lambda_1.set(0, i, random.nextDouble());
                    }
                    for (int i = 0; i < K; i++) {
                        lambda.set(0, i, lambda_1.get(0, i));
                    }
                    // Convert cache to ClassSwitch for QN solving
                    sn.nodetype.set(ind, NodeType.ClassSwitch);
                }

                // Build lambda_cache[u][n x (h+1)]
                Matrix[] lambda_cache = new Matrix[u];
                for (int v = 0; v < u; v++) {
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

                // Get access cost matrix R
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

                // Solve isolated cache
                Ret.cacheGamma gammaResult = Cache_gamma_lp.cache_gamma_lp(lambda_cache, R);
                Matrix gamma = gammaResult.gamma;

                // see _kb/06-solver-catalog.md (JAR-only implementation notes: RMFAnalyzer miss_isolated RANDOM(m) vs FPI)
                // RANDOM(m) and FIFO(m) share the refined mean field (Gast15
                // Thm 1: pi_FIFO(m) = pi_RAND(m)); strict FIFO(m) uses its own
                // position-resolved mean field. Every other strategy (LRU/HLRU/
                // CLIMB/QLRU) has no drift-based fluid model and is refused
                // rather than served a non-fluid FPI/characteristic-time fixed point.
                // RR/FIFO/SFIFO honour a custom access graph (accost) via their
                // general drift; the linear default keeps the refined RAND (FIFO)
                // / linear position-resolved (SFIFO) path. FIFO(m) == RANDOM(m)
                // only for the linear chain (Gast15 Thm 1), so a non-linear graph
                // uses FIFO's own position-resolved drift.
                boolean nonlinear = !accostIsLinear(ch.accost, h);
                if (ch.replacestrat == ReplacementStrategy.RR) {
                    Matrix MU = Cache_miss_rmf.cache_miss_rmf(gamma, m, lambda_cache, ch.accost).MU;
                    for (int v = 0; v < u; v++) {
                        missrate.set(cIdx, v, MU.get(v, 0));
                    }
                } else if (ch.replacestrat == ReplacementStrategy.FIFO) {
                    Matrix MU = nonlinear
                            ? Cache_miss_fifo_rmf.cache_miss_fifo_rmf(gamma, m, lambda_cache, ch.accost).MU
                            : Cache_miss_rmf.cache_miss_rmf(gamma, m, lambda_cache).MU;
                    for (int v = 0; v < u; v++) {
                        missrate.set(cIdx, v, MU.get(v, 0));
                    }
                } else if (ch.replacestrat == ReplacementStrategy.SFIFO) {
                    Matrix MU = Cache_miss_sfifo_rmf.cache_miss_sfifo_rmf(gamma, m, lambda_cache, ch.accost).MU;
                    for (int v = 0; v < u; v++) {
                        missrate.set(cIdx, v, MU.get(v, 0));
                    }
                } else {
                    throw new RuntimeException("SolverFluid supports only RANDOM(m)/FIFO(m) "
                            + "(refined mean field) and strict FIFO(m) (position-resolved "
                            + "mean field) cache replacement; replacement strategy "
                            + ch.replacestrat + " has no drift-based fluid model. Use "
                            + "SolverNC/SolverMVA or SolverLDES for this cache.");
                }

                // Compute hit/miss probabilities for input classes only
                for (int i : inputClass) {
                    double lam = lambda.get(0, i);
                    if (lam > 0) {
                        missprob.set(cIdx, i, missrate.get(cIdx, i) / lam);
                        hitprob.set(cIdx, i, 1.0 - missprob.get(cIdx, i));
                    } else {
                        missprob.set(cIdx, i, 0);
                        hitprob.set(cIdx, i, 0);
                    }
                }
                // Clean NaN values
                for (int i = 0; i < K; i++) {
                    if (Double.isNaN(hitprob.get(cIdx, i))) {
                        hitprob.set(cIdx, i, 0);
                    }
                    if (Double.isNaN(missprob.get(cIdx, i))) {
                        missprob.set(cIdx, i, 0);
                    }
                }

                // Update routing matrix with hit/miss probabilities
                for (int r : inputClass) {
                    // Zero the row for this class at this node
                    for (int j = 0; j < sn.rtnodes.getNumCols(); j++) {
                        sn.rtnodes.set(ind * K + r, j, 0);
                    }
                    // Set hit/miss routing to connected nodes
                    for (int jnd = 0; jnd < I; jnd++) {
                        if (sn.connmatrix.get(ind, jnd) == 1.0) {
                            sn.rtnodes.set(ind * K + r, (int) (jnd * K + hitClass.get(r)), hitprob.get(cIdx, r));
                            sn.rtnodes.set(ind * K + r, (int) (jnd * K + missClass.get(r)), missprob.get(cIdx, r));
                        }
                    }
                }

                // Stochastic complement to get station-level routing
                sn.rt = Dtmc_stochcomp.dtmc_stochcomp(sn.rtnodes, statefulNodeClassesList);
            }

            // Refresh visits
            sn = SnRefreshVisits.snRefreshVisits(sn, sn.chains, sn.rt, sn.rtnodes);

            // Solve the queueing network using the fluid matrix method
            SolverOptions fluidOptions = options.copy();
            fluidOptions.method = "matrix";
            fluidOptions.init_sol = computeFluidInitSol(sn);

            MatrixMethodAnalyzer matrixAnalyzer = new MatrixMethodAnalyzer();
            matrixAnalyzer.analyze(sn, fluidOptions, result);

            // Compute system throughputs XN
            Matrix XN = new Matrix(1, K);
            for (int k = 0; k < K; k++) {
                int refstat = (int) sn.refstat.get(k);
                if (refstat >= 0) {
                    XN.set(0, k, result.TN.get(refstat, k));
                }
            }

            // Update arrival rates to the cache
            // Sum nodevisits across chains
            Matrix nodevisits = null;
            for (Object key : sn.nodevisits.keySet()) {
                Matrix nv = sn.nodevisits.get(key);
                if (nodevisits == null) {
                    nodevisits = new Matrix(nv);
                } else {
                    nodevisits = nodevisits.add(1.0, nv);
                }
            }

            for (int cIdx = 0; cIdx < caches.size(); cIdx++) {
                int ind = caches.get(cIdx);
                NodeParam param = sn.nodeparam.get(sn.nodes.get(ind));
                if (param == null || !(param instanceof CacheNodeParam)) {
                    continue;
                }
                CacheNodeParam ch = (CacheNodeParam) param;
                Matrix hitClass = ch.hitclass;
                List<Integer> inputClass = new ArrayList<Integer>();
                for (int i = 0; i < hitClass.getNumCols(); i++) {
                    if (hitClass.get(0, i) != -1.0) {
                        inputClass.add(i);
                    }
                }

                for (int r : inputClass) {
                    // Find which chain class r belongs to
                    int c = 0;
                    while (c < sn.chains.getNumRows() && sn.chains.get(c, r) == 0.0) {
                        c++;
                    }
                    // Find classes in this chain
                    List<Integer> inchain = new ArrayList<Integer>();
                    for (int j = 0; j < sn.chains.getNumCols(); j++) {
                        if (sn.chains.get(c, j) != 0.0) {
                            inchain.add(j);
                        }
                    }
                    // Sum throughput across chain
                    double sumXN = 0.0;
                    for (int ix : inchain) {
                        sumXN += XN.get(0, ix);
                    }
                    // Compute lambda
                    int refNode = (int) sn.stationToNode.get((int) sn.refstat.get(r));
                    int refClass;
                    if (sn.refclass.get(c) > -1) {
                        refClass = (int) sn.refclass.get(c);
                    } else {
                        refClass = r;
                    }
                    double nvRef = nodevisits.get(refNode, refClass);
                    double nvCache = nodevisits.get(ind, r);
                    if (nvRef > 0) {
                        lambda.set(0, r, sumXN * nvCache / nvRef);
                    }
                }
            }

            // Check convergence
            double diff = 0.0;
            for (int i = 0; i < K; i++) {
                diff += Math.abs(lambda.get(0, i) - lambda_1.get(0, i));
            }
            if (diff < options.iter_tol) {
                convergedIter = it;
                break;
            }
            for (int i = 0; i < K; i++) {
                lambda_1.set(0, i, lambda.get(0, i));
            }
        }

        // Compute CN
        result.CN = new Matrix(1, K);
        Matrix XN = new Matrix(1, K);
        for (int k = 0; k < K; k++) {
            int refstat = (int) sn.refstat.get(k);
            if (refstat >= 0) {
                XN.set(0, k, result.TN.get(refstat, k));
                if (XN.get(0, k) > 0) {
                    result.CN.set(0, k, sn.njobs.get(k) / XN.get(0, k));
                }
            }
        }
        result.XN = XN;

        // Store hit/miss probabilities in the result
        if (result instanceof FluidResult) {
            // hitprob is (numCaches x K), expand to (nnodes x K) for consistency with MVA
            Matrix hitprobFull = new Matrix(sn.nnodes, K);
            Matrix missprobFull = new Matrix(sn.nnodes, K);
            for (int cIdx = 0; cIdx < caches.size(); cIdx++) {
                int ind = caches.get(cIdx);
                for (int k = 0; k < K; k++) {
                    hitprobFull.set(ind, k, hitprob.get(cIdx, k));
                    missprobFull.set(ind, k, missprob.get(cIdx, k));
                }
            }
            ((FluidResult) result).hitProb = hitprobFull;
            ((FluidResult) result).missProb = missprobFull;
        }

        xvecIt = (result instanceof FluidResult) ? ((FluidResult) result).odeStateVec : null;
    }

    @Override
    public Matrix getXVecIt() {
        return xvecIt;
    }

    /**
     * Compute fluid initial solution from network struct.
     * Matches MATLAB solver_fluid_initsol: iterates over stateful station nodes,
     * extracts per-class per-phase populations, skips non-station nodes.
     */
    private static Matrix computeFluidInitSol(NetworkStruct sn) {
        Matrix initSol = new Matrix(1, 0);
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.isstateful.get(ind, 0) != 1) {
                continue;
            }
            int ist = (int) sn.nodeToStation.get(ind);
            // Skip non-station stateful nodes (e.g., Router, ClassSwitch from cache conversion)
            if (ist < 0) {
                continue;
            }
            int isf = (int) sn.nodeToStateful.get(ind);
            Matrix state_i = sn.state.get(sn.nodes.get(ind));
            if (state_i == null) {
                continue;
            }

            State.StateMarginalStatistics stats = ToMarginal.toMarginal(sn, ind, state_i, null, null, null, null, null);
            Matrix nir = stats.nir;
            java.util.List<Matrix> kir_i = stats.kir;
            if (kir_i == null || kir_i.isEmpty()) {
                continue;
            }

            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            int rMax = kir_i.get(0).getNumCols();

            for (int r = 0; r < rMax; r++) {
                int kMax = sn.mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).length();
                for (int k = 0; k < kMax; k++) {
                    if (Double.isNaN(sn.rates.get(ist, r))) {
                        continue;
                    }
                    initSol.expandMatrix(1, initSol.getNumCols() + 1, initSol.getNumElements() + 1);
                    if (sched == SchedStrategy.EXT) {
                        // Source: use kir directly
                        initSol.set(0, initSol.getNumCols() - 1, kir_i.get(k).get(0, r));
                    } else {
                        // Station: phase 0 = waiting buffer jobs, phases 1+ = in-service
                        if (k == 0) {
                            double sumKir = 0;
                            for (int m = 1; m < kir_i.size(); m++) {
                                sumKir += kir_i.get(m).get(0, r);
                            }
                            initSol.set(0, initSol.getNumCols() - 1, nir.get(0, r) - sumKir);
                        } else {
                            initSol.set(0, initSol.getNumCols() - 1, kir_i.get(k).get(0, r));
                        }
                    }
                }
            }
        }
        return initSol;
    }

    /**
     * True when every per-(user,item) access graph is the standard linear chain
     * (miss -> list 1, hit in list a -> list a+1, self-loop on the top list).
     */
    private static boolean accostIsLinear(Matrix[][] accost, int h) {
        if (accost == null) {
            return true;
        }
        double[][] lin = new double[h + 1][h + 1];
        lin[0][1] = 1.0;
        for (int a = 1; a < h; a++) {
            lin[a][a + 1] = 1.0;
        }
        lin[h][h] = 1.0;
        for (int v = 0; v < accost.length; v++) {
            if (accost[v] == null) {
                continue;
            }
            for (int k = 0; k < accost[v].length; k++) {
                Matrix g = accost[v][k];
                if (g == null) {
                    continue;
                }
                for (int a = 0; a <= h; a++) {
                    for (int b = 0; b <= h; b++) {
                        if (Math.abs(g.get(a, b) - lin[a][b]) > 1e-9) {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }
}
