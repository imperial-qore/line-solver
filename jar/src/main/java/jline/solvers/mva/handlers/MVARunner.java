/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva.handlers;

import jline.io.Ret;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.NodeParam;
import jline.GlobalConstants;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodeparam.ForkNodeParam;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.AvgHandle;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.SolverAvgHandles;
import jline.solvers.mva.MVAResult;
import jline.solvers.mva.SolverMVA;
import jline.lang.FeatureSet;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.api.sn.SnGetResidTFromRespT.snGetResidTFromRespT;
import static jline.api.sn.SnHasPolling.snHasPolling;
import static jline.io.InputOutput.*;
import static jline.solvers.mva.analyzers.Solver_mva_analyzer.solver_mva_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_marie_analyzer.solver_mva_marie_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_cache_analyzer.solver_mva_cache_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_cacheqn_analyzer.solver_mva_cacheqn_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_retrieval_analyzer.solver_mva_retrieval_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_polling_analyzer.solver_mva_polling_analyzer;
import static jline.solvers.mva.analyzers.Solver_mva_qsys_analyzer.solver_mva_qsys_analyzer;
import static jline.solvers.mva.analyzers.Solver_mvald_analyzer.solver_mvald_analyzer;
import jline.solvers.mva.SolverMVAOIAnalyzer;
import jline.lang.ModelAdapter;
import static jline.lang.ModelAdapter.*;
import jline.solvers.fj.FJFixedPoint;
import static jline.util.Utils.isInf;

public class MVARunner {

    private final Network model;
    private final boolean enableChecks;
    protected NetworkStruct sn;
    protected SolverOptions options;
    protected MVAResult res;
    /**
     * MMT fork-join transformation carried in from (and handed back to) the owning
     * SolverMVA, so that it can be reused across outer iterations. Null means no
     * usable cache; see runAnalyzer for the validity test.
     */
    private Ret.FJApprox mmtCache;

    /**
     * Auxiliary-class arrival rates of the fork-join (MMT) fixed point, carried in
     * from (and handed back to) the owning SolverMVA so that an outer iteration
     * such as SolverLN resumes the fixed point instead of restarting it from
     * GlobalConstants.FineTol. Null means no retained iterate.
     */
    private Matrix fjForkLambda;

    public MVARunner(Network model, SolverOptions options, boolean enableChecks) {
        this(model, options, enableChecks, null, null);
    }

    public MVARunner(Network model, SolverOptions options, boolean enableChecks, Ret.FJApprox mmtCache) {
        this(model, options, enableChecks, mmtCache, null);
    }

    /** The MMT transformation this run built or reused, for the owning solver to keep. */
    public Ret.FJApprox getMmtCache() {
        return this.mmtCache;
    }

    /** The MMT fixed-point iterate this run ended on, for the owning solver to keep. */
    public Matrix getFjForkLambda() {
        return this.fjForkLambda;
    }

    public MVARunner(Network model, SolverOptions options, boolean enableChecks, Ret.FJApprox mmtCache,
                     Matrix fjForkLambda) {
        this.mmtCache = mmtCache;
        this.fjForkLambda = fjForkLambda;
        this.model = model;
        this.options = options;
        this.enableChecks = enableChecks;
        this.sn = this.model.getStruct(false);
        this.res = null;
    }

    /**
     * runAnalyzer() method from LINE.
     *
     * @param avgHandles - the average handles for the model
     * @return - the performance measures corresponding to the given network
     */
    public SolverResult runAnalyzer(SolverAvgHandles avgHandles) {
        long T0 = System.nanoTime();
        int iter = 0;

        // Case 'java' or 'jline.amva' can be ignored, so we can remove the switch
        this.sn = this.model.getStruct(false);

        if (this.enableChecks) {
            // see _kb/06-solver-catalog.md for rationale
            String gateMethod = this.options.method;
            if ("default".equals(gateMethod)) {
                boolean allOpen = true;
                for (int r = 0; r < this.sn.nclasses; r++) {
                    if (!Double.isInfinite(this.sn.njobs.get(r))) { allOpen = false; break; }
                }
                if ((this.sn.nclasses == 1) && allOpen
                        && jline.api.sn.SnHasBurstyArrival.snHasBurstyArrival(this.sn)) {
                    gateMethod = "rqna";
                }
            }
            String reason = FeatureSet.supportsReason(
                    SolverMVA.methodFeatureSet(gateMethod),
                    this.model.getUsedLangFeatures());
            if (!reason.isEmpty()) {
                line_error(mfilename(new Object() {
                }), "This model contains features not supported by the solver. " + reason);
                return null;
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        Matrix rtOrig = this.sn.rt != null ? this.sn.rt.copy() : null;
        // see _kb/06-solver-catalog.md for rationale
        FJFixedPoint.FJState fjState = new FJFixedPoint.FJState(this.mmtCache, this.fjForkLambda);
        FJFixedPoint.FJOutcome fjOut = FJFixedPoint.run(this.model, this.sn, this.options, fjState,
                new FJFixedPoint.InnerSolve() {
                    @Override
                    public MVAResult solve(NetworkStruct snIn, SolverOptions opts) {
                        return MVARunner.this.dispatch(snIn, opts);
                    }
                }, T0);
        MVAResult ret = fjOut.ret;
        Matrix QN = fjOut.QN;
        iter += fjOut.iter;
        this.sn = fjOut.sn;
        this.mmtCache = fjOut.state.mmtCache;
        this.fjForkLambda = fjOut.state.fjForkLambda;
        this.sn = this.model.getStruct(true);

        // see _kb/06-solver-catalog.md for rationale
        for (int ind = 0; ind < this.sn.nnodes; ind++) {
            if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                CacheNodeParam cacheParam = (CacheNodeParam) this.sn.nodeparam.get(cacheNode);
                if (cacheParam != null) {
                    Matrix hitProb = cacheNode.getHitRatio();
                    Matrix missProb = cacheNode.getMissRatio();
                    Matrix delayedProb = cacheNode.getDelayedHitRatio();
                    Matrix hitProbList = cacheNode.getHitRatioByList();
                    if (hitProb != null && !hitProb.isEmpty()) {
                        cacheParam.actualhitprob = hitProb;
                    }
                    if (missProb != null && !missProb.isEmpty()) {
                        cacheParam.actualmissprob = missProb;
                    }
                    if (delayedProb != null && !delayedProb.isEmpty()) {
                        cacheParam.actualdelayedhitprob = delayedProb;
                    }
                    if (hitProbList != null && !hitProbList.isEmpty()) {
                        cacheParam.actualhitproblist = hitProbList;
                    }
                }
            }
        }

        // Compute average arrival rate at steady-state
        AvgHandle T = avgHandles.getAvgTputHandles();
        Matrix rtRefreshed = sn.rt;
        if (rtOrig != null) {
            sn.rt = rtOrig;
        }
        Matrix AN = snGetArvRFromTput(sn, ret.TN, T);
        if (rtRefreshed != null) {
            sn.rt = rtRefreshed;
        }

        AvgHandle W = avgHandles.getAvgResidTHandles();
        Matrix WN;
        if (ret.RN != null && !ret.RN.isEmpty()) {
            WN = snGetResidTFromRespT(sn, ret.RN, W);
        } else {
            // Handle case where RN is null or empty by creating a zero matrix
            WN = new Matrix(sn.nstations, sn.nclasses);
        }

        this.res = new MVAResult();
        this.res.method = ret.method;
        this.res.QN = ret.QN; //
        this.res.RN = ret.RN; //
        this.res.XN = ret.XN;
        this.res.UN = ret.UN; //
        this.res.TN = ret.TN; //
        this.res.CN = ret.CN;
        this.res.AN = AN;
        this.res.WN = WN;
        this.res.runtime = (System.nanoTime() - T0) / 1000000000.0;
        this.res.iter = iter;
        this.res.logNormConstAggr = ret.logNormConstAggr;
        if (Solver.timeExceeded(T0, this.options.timeout)) {
            this.res.timedOut = true;
            line_warning(mfilename(new Object() {
            }), "Solver exceeded the wall-clock time budget (options.timeout=" + this.options.timeout + "s); returning the interim solution.");
        }

        return this.res;
    }

    /**
     * Check if the node types match exactly the expected types (order independent)
     */
    private boolean checkNodeTypes(List<NodeType> nodeTypes, NodeType... expectedTypes) {
        if (nodeTypes.size() != expectedTypes.length) {
            return false;
        }
        
        // Count occurrences of each node type
        Map<NodeType, Integer> actualCounts = new HashMap<>();
        Map<NodeType, Integer> expectedCounts = new HashMap<>();
        
        for (NodeType type : nodeTypes) {
            actualCounts.put(type, actualCounts.getOrDefault(type, 0) + 1);
        }
        
        for (NodeType type : expectedTypes) {
            expectedCounts.put(type, expectedCounts.getOrDefault(type, 0) + 1);
        }

        return actualCounts.equals(expectedCounts);
    }

    /**
     * True if the (Source-Queue-Sink) model is a single-server HOL priority
     * queue with Poisson (ca=1) arrivals for every class, so the exact
     * non-preemptive Cobham formula applies.
     */
    /**
     * True if the (Source-Queue-Sink) model is a single-server DPS queue with
     * Poisson (ca=1) arrivals and exponential (cs=1) service for every class,
     * so the numerically-exact M/M/1-DPS solver applies.
     */
    private boolean snIsDpsQsys(NetworkStruct sn) {
        int source_ist = -1;
        int queue_ist = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                source_ist = (int) sn.nodeToStation.get(i);
            } else if (sn.nodetype.get(i) == NodeType.Queue) {
                queue_ist = (int) sn.nodeToStation.get(i);
            }
        }
        if (source_ist < 0 || queue_ist < 0) {
            return false;
        }
        if (sn.sched.get(sn.stations.get(queue_ist)) != SchedStrategy.DPS) {
            return false;
        }
        double s = sn.nservers.get(queue_ist);
        if (Double.isInfinite(s) || (int) s != 1) {
            return false;
        }
        for (int r = 0; r < sn.nclasses; r++) {
            double scvA = sn.scv.get(source_ist, r);
            if (Double.isFinite(scvA) && Math.abs(scvA - 1.0) > 1e-6) {
                return false;
            }
            double scvS = sn.scv.get(queue_ist, r);
            if (Double.isFinite(scvS) && Math.abs(scvS - 1.0) > 1e-6) {
                return false;
            }
        }
        return true;
    }

    private boolean snIsHolQsys(NetworkStruct sn) {
        int source_ist = -1;
        int queue_ist = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                source_ist = (int) sn.nodeToStation.get(i);
            } else if (sn.nodetype.get(i) == NodeType.Queue) {
                queue_ist = (int) sn.nodeToStation.get(i);
            }
        }
        if (source_ist < 0 || queue_ist < 0) {
            return false;
        }
        if (sn.sched.get(sn.stations.get(queue_ist)) != SchedStrategy.HOL) {
            return false;
        }
        double s = sn.nservers.get(queue_ist);
        if (Double.isInfinite(s) || (int) s != 1) {
            return false;
        }
        for (int r = 0; r < sn.nclasses; r++) {
            double scv = sn.scv.get(source_ist, r);
            if (Double.isFinite(scv) && Math.abs(scv - 1.0) > 1e-6) {
                return false;
            }
        }
        return true;
    }


    /**
     * One inner solve of the MVA analyzer dispatch: pick the analyzer that fits
     * SN and return its metrics. This is the callback that
     * {@link FJFixedPoint} drives on each pass of the fork-join fixed point; on
     * a model without forks it runs exactly once. The body is the dispatch that
     * used to sit inline in runAnalyzer.
     *
     * @param snIn the struct to solve, the transformed one under a fork-join model
     * @param opts the solver options
     * @return the metrics of that solve
     */
    private MVAResult dispatch(NetworkStruct snIn, SolverOptions opts) {
        this.sn = snIn;
        MVAResult ret = new MVAResult();

            if (this.options.method.equals("exact") && !this.model.hasProductFormSolution()) {
                line_error(mfilename(new Object() {
                }), "The exact method requires the model to have a product-form solution. This model does not have one.");
            }
            if (this.options.method.equals("mva") && !this.model.hasProductFormSolution()) {
                line_warning(mfilename(new Object() {
                }), "The exact method requires the model to have a product-form solution. This model does not have one. SolverMVA will return an approximation generated by an exact MVA algorithm.");
            }
            String method = this.options.method;
            line_debug(this.options.verbose, String.format("MVA analyzer starting: method=%s, nclasses=%d, nclosed=%d, nnodes=%d",
                method, this.sn.nclasses, this.sn.nclosedjobs, this.sn.nnodes));

            if (this.sn.nclasses == 1 && this.sn.nclosedjobs == 0 && this.sn.nodetype.size() == 3 && checkNodeTypes(this.sn.nodetype, NodeType.Source, NodeType.Queue, NodeType.Sink)) {
                // Single-class open queueing system
                line_debug(this.options.verbose, "Detected open queueing system (Source-Queue-Sink), calling solver_mva_qsys_analyzer");
                ret = solver_mva_qsys_analyzer(this.sn, this.options.copy());
            } else if (this.sn.nclasses > 1 && this.sn.nclosedjobs == 0 && this.sn.nodetype.size() == 3 && checkNodeTypes(this.sn.nodetype, NodeType.Source, NodeType.Queue, NodeType.Sink) && snHasPolling(this.sn)) {
                // Multi-class open polling system
                line_debug(this.options.verbose, "Detected multi-class polling system (Source-Queue-Sink with POLLING scheduling), calling solver_mva_polling_analyzer");
                ret = solver_mva_polling_analyzer(this.sn, this.options.copy());
            } else if (this.sn.nclasses > 1 && this.sn.nclosedjobs == 0 && this.sn.nodetype.size() == 3 && checkNodeTypes(this.sn.nodetype, NodeType.Source, NodeType.Queue, NodeType.Sink) && snIsHolQsys(this.sn)) {
                // see _kb/06-solver-catalog.md for rationale
                line_debug(this.options.verbose, "Detected multi-class HOL priority queue (Source-Queue-Sink), calling solver_mva_qsys_prio_analyzer");
                ret = jline.solvers.mva.analyzers.Solver_mva_qsys_analyzer.solver_mva_qsys_prio_analyzer(this.sn, this.options.copy());
            } else if (this.sn.nclasses > 1 && this.sn.nclasses <= 3 && this.sn.nclosedjobs == 0 && this.sn.nodetype.size() == 3 && checkNodeTypes(this.sn.nodetype, NodeType.Source, NodeType.Queue, NodeType.Sink) && snIsDpsQsys(this.sn)) {
                // see _kb/06-solver-catalog.md for rationale
                line_debug(this.options.verbose, "Detected multi-class M/M/1-DPS queue (Source-Queue-Sink), calling solver_mva_qsys_dps_analyzer");
                ret = jline.solvers.mva.analyzers.Solver_mva_qsys_analyzer.solver_mva_qsys_dps_analyzer(this.sn, this.options.copy());
            } else if (jline.solvers.nc.SolverNC.hasRetrievalCache(this.sn)) {
                // Delayed-hit cache with a retrieval system: open (Source) -> FPI
                // analyzer; closed integrated -> da_cacheqn_retrieval driver.
                boolean hasSource = false;
                for (int i = 0; i < this.sn.nodetype.size(); i++) if (this.sn.nodetype.get(i) == NodeType.Source) { hasSource = true; break; }
                if (hasSource) {
                    line_debug(this.options.verbose, "Detected open delayed-hit retrieval cache, calling solver_mva_retrieval_analyzer");
                    ret = solver_mva_retrieval_analyzer(this.sn, this.options.copy());
                } else {
                    line_debug(this.options.verbose, "Detected closed integrated delayed-hit retrieval cache, calling solver_mva_cacheqn_retrieval_analyzer");
                    ret = jline.solvers.mva.analyzers.Solver_mva_cacheqn_retrieval_analyzer.solver_mva_cacheqn_retrieval_analyzer(this.sn, this.options.copy());
                }
            } else if (this.sn.nclosedjobs == 0 && this.sn.nodetype.size() == 3 && checkNodeTypes(this.sn.nodetype, NodeType.Source, NodeType.Cache, NodeType.Sink)) {
                // Non-rentrant cache
                // Random initialisation
                for (int ind = 0; ind < this.sn.nnodes; ind++) {
                    if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                        Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                        Matrix prob = new Matrix(cacheNode.getHitClass());
                        for (int i = 0; i < prob.getNumRows(); i++) {
                            for (int j = 0; j < prob.getNumCols(); j++) {
                                if (prob.get(i, j) > 0) {
                                    prob.set(i, j, 0.5);
                                }
                            }
                        }
                        cacheNode.setResultHitProb(prob);
                        Matrix missProb = new Matrix(prob.getNumRows(), prob.getNumCols());
                        for (int i = 0; i < prob.getNumRows(); i++) {
                            for (int j = 0; j < prob.getNumCols(); j++) {
                                missProb.set(i, j, 1 - prob.get(i, j));
                            }
                        }
                        cacheNode.setResultMissProb(missProb);
                    }
                }
                this.model.refreshChains(true);
                // Start iteration
                line_debug(this.options.verbose, "Detected cache system (Source-Cache-Sink), calling solver_mva_cache_analyzer");
                ret = solver_mva_cache_analyzer(this.sn, this.options.copy());

                for (int ind = 0; ind < this.sn.nnodes; ind++) {
                    if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                        Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                        Matrix hitClass = cacheNode.getHitClass();
                        Matrix missClass = cacheNode.getMissClass();
                        Matrix hitProb = new Matrix(1, hitClass.length());
                        for (int k = 0; k < hitClass.length(); k++) {
                            int chain_k = 0;
                            for (; chain_k < this.sn.chains.getNumRows(); chain_k++) {
                                if (this.sn.chains.get(chain_k, k) > 0)
                                    break;
                            }
                            Matrix inchain = new Matrix(1, this.sn.chains.getNumCols());
                            for (int i = 0; i < inchain.getNumCols(); i++) {
                                inchain.set(0, i, this.sn.chains.get(chain_k, i) > 0 ? 1 : 0);
                            }
                            int h = (int) hitClass.get(k);
                            int m = (int) missClass.get(k);
                            if (h > -1 && m > -1) {
                                double sumXN = 0;
                                for (int i = 0; i < inchain.getNumCols(); i++) {
                                    if (inchain.get(i) > 0 && !Double.isNaN(ret.XN.get(i))) {
                                        sumXN += ret.XN.get(i);
                                    }
                                }
                                hitProb.set(k, ret.XN.get(h) / sumXN);
                            }
                        }
                        Matrix missProb = new Matrix(1, hitClass.length());
                        for (int i = 0; i < hitClass.length(); i++) {
                            missProb.set(i, 1 - hitProb.get(i));
                        }
                        cacheNode.setResultHitProb(hitProb);
                        cacheNode.setResultMissProb(missProb);
                    }
                }
                this.model.refreshStruct(true);
            } else {
                // see _kb/06-solver-catalog.md for rationale
                int noi_idx = -1;
                boolean closedNet = true;
                for (int r = 0; r < this.sn.nclasses; r++) {
                    if (Double.isInfinite(this.sn.njobs.get(0, r))) { closedNet = false; break; }
                }
                if (closedNet && this.sn.nodeparam != null) {
                    for (int i = 0; i < this.sn.nstations; i++) {
                        jline.lang.nodes.Station st = this.sn.stations.get(i);
                        jline.lang.constant.SchedStrategy s = this.sn.sched.get(st);
                        if (s != jline.lang.constant.SchedStrategy.PAS
                                && s != jline.lang.constant.SchedStrategy.OI) {
                            continue;
                        }
                        Object param = this.sn.nodeparam.get(st);
                        if (!(param instanceof jline.lang.nodeparam.QueueNodeParam)) {
                            continue;
                        }
                        jline.lang.nodeparam.QueueNodeParam qp = (jline.lang.nodeparam.QueueNodeParam) param;
                        if (qp.svcRateFun == null || qp.swapGraph == null) {
                            continue;
                        }
                        if (qp.swapGraph.isEmpty() || qp.swapGraph.elementMaxAbs() == 0) {
                            noi_idx = i;
                            break;
                        }
                    }
                }

                boolean cachePresent = false;
                for (NodeType t : this.sn.nodetype) {
                    if (t == NodeType.Cache) {
                        cachePresent = true;
                        break;
                    }
                }

                // An OI/PAS station is admissible only when every other station is
                // product-form, which is exactly the NC-oi gate.
                boolean hasOIStation = false;
                for (int i = 0; i < this.sn.nstations; i++) {
                    jline.lang.constant.SchedStrategy s = this.sn.sched.get(this.sn.stations.get(i));
                    if (s == jline.lang.constant.SchedStrategy.PAS
                            || s == jline.lang.constant.SchedStrategy.OI) {
                        hasOIStation = true;
                        break;
                    }
                }
                boolean oiMethod = "exact".equals(method) || "default".equals(method);
                boolean oiExact = noi_idx >= 0 && oiMethod
                        && jline.solvers.nc.handlers.Solver_nc_oi.nc_is_oi_model(this.sn);

                if (oiExact) {
                    // Order-independent queueing network
                    line_debug(this.options.verbose, "Detected order-independent network, routing to SolverMVAOIAnalyzer");
                    SolverMVAOIAnalyzer analyzer = new SolverMVAOIAnalyzer(this.sn, this.options);
                    SolverMVAOIAnalyzer.AnalysisResults results = analyzer.analyze();
                    ret = new MVAResult();
                    ret.QN = new Matrix(results.QN);
                    ret.UN = new Matrix(results.UN);
                    ret.RN = new Matrix(results.RN);
                    ret.TN = new Matrix(results.TN);
                    ret.CN = new Matrix(results.CN);
                    ret.XN = new Matrix(results.XN);
                    ret.runtime = results.runtime;
                    ret.iter = results.iter;
                    ret.method = results.method;
                } else if (hasOIStation) {
                    // An OI/PAS station carries a rank-rate function mu(n) of the whole
                    // per-class occupancy. The AMVA iteration only ever sees sn.rates
                    // (the single-job rates), so it cannot represent such a station and
                    // would silently return a zero queue-length there. MVA therefore
                    // supports OI/PAS stations only via the exact path above.
                    throw new RuntimeException(String.format(
                            "SolverMVA supports order-independent (OI) and pass-and-swap (PAS) stations only%n"
                            + "through its exact order-independent analyzer, which requires method 'default'%n"
                            + "or 'exact' (got '%s'), an empty/zero swap graph at every OI/PAS station, a closed%n"
                            + "model, and every other station to be product-form (INF, PS, LCFS-PR, SIRO, or%n"
                            + "class-independent-rate FCFS). Use SolverCTMC or SolverLDES for this model.", method));
                } else if (cachePresent) {
                    // Integrated Cache Queueing
                    line_debug(this.options.verbose, "Detected cache+queueing network, calling solver_mva_cacheqn_analyzer");
                    ret = solver_mva_cacheqn_analyzer(this.sn, this.options.copy());
                    for (int ind = 0; ind < this.sn.nnodes; ind++) {
                        if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                            Cache cache = (Cache) this.model.getNodes().get(ind);
                            cache.setResultHitProb(Matrix.extractRows(ret.hitProb, ind, ind + 1, null));
                            cache.setResultMissProb(Matrix.extractRows(ret.missProb, ind, ind + 1, null));
                            if (ret.cacheItemProb != null && ret.cacheItemProb.containsKey(ind)) {
                                cache.setResultItemProb(ret.cacheItemProb.get(ind));
                            }
                        }
                    }
                    this.model.refreshStruct(true);
                } else {
                    // Ordinary queueing network
                    switch (method) {
                        case "marie":
                        case "amva.marie":
                            line_debug(this.options.verbose, "Using Marie aggregation-decomposition analyzer");
                            ret = solver_mva_marie_analyzer(this.sn, this.options.copy());
                            break;
                        case "aba.upper":
                        case "aba.lower":
                        case "bjb.upper":
                        case "bjb.lower":
                        case "pb.upper":
                        case "pb.lower":
                        case "gb.upper":
                        case "gb.lower":
                        case "sb.upper":
                        case "sb.lower":
                        case "harel.lower":
                        case "harel.upper":
                        case "mwba.upper":
                        case "mwba.lower":
                            // Bounds are served by the dedicated SolverBA solver
                            // (matching MATLAB, where bounds were moved out of
                            // SolverMVA). Redirect rather than silently answer.
                            throw new RuntimeException("Bound method '" + method + "' is served by SolverBA, "
                                    + "not SolverMVA. Use SolverBA(model, \"method\", \"" + method + "\").");
                        default: // this is the main solution block for standard queueing networks
                            if ((this.sn.lldscaling != null && !this.sn.lldscaling.isEmpty()) ||
                                    (this.sn.cdscaling != null && !this.sn.cdscaling.isEmpty()) ||
                                    (this.sn.jdscaling != null && !this.sn.jdscaling.isEmpty())) {
                                line_debug(this.options.verbose, "Detected load-/class-/joint-dependent scaling, calling solver_mvald_analyzer");
                                ret = solver_mvald_analyzer(this.sn, this.options.copy());
                            } else {
                                // see _kb/06-solver-catalog.md for rationale
                                SolverOptions subopts = this.options.copy();
                                if (this.model.hasFork() && "default".equals(subopts.method)) {
                                    subopts.method = "amva";
                                }
                                line_debug(this.options.verbose, "Using standard MVA method, calling solver_mva_analyzer");
                                ret = solver_mva_analyzer(this.sn, subopts);
                            }
                    }
                }
            }
        return ret;
    }
}
