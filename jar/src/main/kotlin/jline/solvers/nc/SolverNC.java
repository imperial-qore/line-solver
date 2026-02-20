/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.GlobalConstants;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.solvers.AvgHandle;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.util.Utils;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static jline.api.sn.SnGetArvRFromTputKt.snGetArvRFromTput;
import static jline.api.sn.SnHasMultiServerKt.snHasMultiServer;
import static jline.io.InputOutputKt.*;
import static jline.solvers.nc.analyzers.Solver_nc_analyzerKt.solver_nc_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_cache_analyzerKt.solver_nc_cache_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_cache_qn_analyzerKt.solver_nc_cache_qn_analyzer;
import static jline.solvers.nc.analyzers.Solver_ncld_analyzerKt.solver_ncld_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_lossn_analyzerKt.solver_nc_lossn_analyzer;
import static jline.api.sn.SnHasClosedClassesKt.snHasClosedClasses;
import jline.lang.constant.DropStrategy;
import static jline.solvers.nc.handlers.Solver_nc_margKt.solver_nc_marg;
import static jline.solvers.nc.handlers.Solver_nc_jointKt.solver_nc_joint;
import static jline.solvers.nc.handlers.Solver_nc_margaggrKt.solver_nc_margaggr;
import static jline.solvers.nc.handlers.Solver_nc_jointaggrKt.solver_nc_jointaggr;
import static jline.solvers.nc.handlers.Solver_nc_jointaggr_ldKt.solver_nc_jointaggr_ld;
import static jline.api.sn.SnGetProductFormParamsKt.snGetProductFormParams;
import jline.io.Ret;
import static jline.api.pfqn.nc.Pfqn_stdfKt.pfqn_stdf;
import static jline.api.pfqn.nc.Pfqn_stdf_heurKt.pfqn_stdf_heur;
import static jline.api.mam.Map_cdfKt.map_cdf;


/**
 * Solver for Normalizing Constant (NC) method applied to closed queueing networks.
 * 
 * <p>SolverNC implements the normalizing constant approach for computing exact
 * performance measures of closed product-form queueing networks. The normalizing
 * constant G(N) represents the partition function that ensures steady-state
 * probabilities sum to one.</p>
 * 
 * <p>Key NC solver capabilities:
 * <ul>
 *   <li>Exact normalizing constant computation</li>
 *   <li>Convolution algorithm implementation</li>
 *   <li>Load-dependent service station support</li>
 *   <li>Cache-enabled queueing network analysis</li>
 *   <li>Joint and marginal probability computation</li>
 *   <li>State probability aggregation methods</li>
 * </ul>
 * </p>
 * 
 * <p>The solver supports various computation methods including standard convolution,
 * tree convolution, and specialized algorithms for cache networks and load-dependent
 * stations. Results include exact performance metrics and state probabilities.</p>
 * 
 * @see jline.api.pfqn.nc
 * @see NCResult
 * @see NCOptions
 * @since 1.0
 */
public class SolverNC extends NetworkSolver {

    public SolverNC(Network model, SolverOptions options) {
        super(model, "SolverNC", options);
        this.sn = model.getStruct(false);
        this.result = new NCResult();
    }

    public SolverNC(Network model) {
        super(model, "SolverNC", new NCOptions());
        this.sn = model.getStruct(false);
        this.result = new NCResult();
    }

    public SolverNC(Network model, String method) {
        super(model, "SolverNC", new NCOptions().method(method));
        this.sn = model.getStruct(false);
        this.result = new NCResult();
    }

    public SolverNC(Network model, Object... varargin) {
        super(model, "SolverNC", new NCOptions());
        this.options = SolverNC.parseOptions(this.options, varargin);
        this.sn = model.getStruct(false);
        this.result = new NCResult();
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.NC);
    }

    /**
     * Returns the feature set supported by the NC solver
     *
     * @return - the feature set supported by the NC solver
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink", "Source",
                "ClassSwitch", "Delay", "DelayStation", "Queue",
                "APH", "Coxian", "Erlang", "Det", "Exp", "HyperExp",
                "StatelessClassSwitcher", "InfiniteServer",
                "SharedServer", "Buffer", "Dispatcher",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS", "SchedStrategy_SIRO",
                "SchedStrategy_LCFS", "SchedStrategy_LCFSPR",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "SchedStrategy_FCFS", "ClosedClass", "SelfLoopingClass",
                "Cache", "CacheClassSwitcher", "OpenClass",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_LRU", "ReplacementStrategy_SFIFO",
                "LoadDependence"
        });
        return featSupported;
    }

    public Double getProb(Node node, Matrix state) {
        if (GlobalConstants.DummyMode) {
            return Double.NaN;
        }

        Matrix state_new = state.copy();

        if (state_new == null || state_new.isEmpty()) {
            state_new = sn.state.get((StatefulNode) this.model.getNodes().get((int) sn.nodeToStateful.get(node.getNodeIndex())));
        }

        long startTimeMillis = System.nanoTime();
        NetworkStruct sn = getStruct();
        int ist = (int) sn.nodeToStation.get(node.getNodeIndex());
        int isf = (int) sn.nodeToStateful.get(node.getNodeIndex());
        sn.state.put(this.model.getStatefulNodes().get(isf), state_new);
        resetRandomGeneratorSeed(options.seed);

        Matrix Pnir;

        NCResult ncResult = (NCResult) this.result;
        if (ncResult != null && ncResult.prob != null
                && ncResult.prob.logNormConstAggr != null && !Utils.isInf(ncResult.prob.logNormConstAggr) && !Double.isNaN(ncResult.prob.logNormConstAggr)) {
            Pnir = solver_nc_marg(sn, this.options, ncResult.prob.logNormConstAggr).lPr;
        } else {
            SolverNCMargReturn ret = solver_nc_marg(sn, this.options, null);
            Pnir = ret.lPr;
            ((NCResult) this.result).prob.logNormConstAggr = ret.lG;
        }
        ((NCResult) this.result).solver = this.name;
        ((NCResult) this.result).prob.marginal = Pnir;
        long endTimeMillis = System.nanoTime();
        double runtime = (double) (endTimeMillis - startTimeMillis) / 1000000000.0;
        this.result.runtime = runtime;

        return Pnir.get(ist);
    }

    /**
     * Get aggregated probability for a specific node and state
     *
     * @param node The node to compute probability for
     * @param state_a The aggregated state (optional, uses current state if null)
     * @return The aggregated probability value
     */
    public Double getProbAggr(Node node, Matrix state_a) {
        if (GlobalConstants.DummyMode) {
            return Double.NaN;
        }

        long startTimeMillis = System.nanoTime();
        NetworkStruct sn = getStruct();
        int ist = (int) sn.nodeToStation.get(node.getNodeIndex());
        int isf = (int) sn.nodeToStateful.get(node.getNodeIndex());
        
        Matrix state_new = state_a;
        if (state_new == null || state_new.isEmpty()) {
            state_new = sn.state.get((StatefulNode) this.model.getNodes().get((int) sn.nodeToStateful.get(node.getNodeIndex())));
        }
        
        sn.state.put(this.model.getStatefulNodes().get(isf), state_new);
        resetRandomGeneratorSeed(options.seed);

        Matrix Pnir;
        NCResult ncResult = (NCResult) this.result;
        
        // Use aggregated marginal solver for proper aggregated probability computation
        if (ncResult != null && ncResult.prob != null
                && ncResult.prob.logNormConstAggr != null && !Utils.isInf(ncResult.prob.logNormConstAggr) && !Double.isNaN(ncResult.prob.logNormConstAggr)) {
            Pnir = solver_nc_margaggr(sn, this.options, ncResult.prob.logNormConstAggr).lPr;
        } else {
            SolverNCMargReturn ret = solver_nc_margaggr(sn, this.options, null);
            Pnir = ret.lPr;
            ((NCResult) this.result).prob.logNormConstAggr = ret.lG;
        }
        
        ((NCResult) this.result).solver = this.name;
        ((NCResult) this.result).prob.marginal = Pnir;
        long endTimeMillis = System.nanoTime();
        double runtime = (double) (endTimeMillis - startTimeMillis) / 1000000000.0;
        this.result.runtime = runtime;

        return Pnir.get(ist);
    }

    /**
     * Get aggregated probability for a specific node using current state
     *
     * @param node The node to compute probability for
     * @return The aggregated probability value
     */
    public Double getProbAggr(Node node) {
        return getProbAggr(node, null);
    }

    /**
     * Get the log normalization constant for aggregated probabilities
     *
     * @return The log normalization constant
     * @throws IllegalAccessException if analysis fails
     */
    public ProbabilityResult getProbNormConstAggr() {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult(Double.NaN, true);
        }

        try {
            this.runAnalyzer();
        } catch (Exception e) {
            return new ProbabilityResult(Double.NaN, true);
        }
        NCResult ncResult = (NCResult) this.result;
        return new ProbabilityResult(ncResult.prob.logNormConstAggr, true);
    }

    /**
     * Get system-wide joint probability 
     *
     * @return The joint probability value
     */
    public ProbabilityResult getProbSys() {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult(Double.NaN);
        }

        long startTimeMillis = System.nanoTime();
        NetworkStruct sn = getStruct();
        resetRandomGeneratorSeed(options.seed);
        
        SolverNCJointReturn result = solver_nc_joint(sn, this.options);
        
        NCResult ncResult = (NCResult) this.result;
        ncResult.solver = this.name;
        ncResult.prob.logNormConstAggr = result.lG;
        ncResult.prob.joint = result.Pr;
        
        long endTimeMillis = System.nanoTime();
        double runtime = (double) (endTimeMillis - startTimeMillis) / 1000000000.0;
        ncResult.runtime = runtime;

        return new ProbabilityResult(result.Pr);
    }

    /**
     * Get aggregated system-wide joint probability
     * 
     * @return The aggregated joint probability value
     */
    public ProbabilityResult getProbSysAggr() {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult(Double.NaN);
        }

        long startTimeMillis = System.nanoTime();
        // Get a fresh sn from the model - runAnalyzer may have modified
        // the cached sn (e.g., converting multiserver nservers to lldscaling),
        // but solver_nc_jointaggr needs the original nservers to build
        // its own mu matrix, matching MATLAB's value-copy semantics.
        // Force a hard refresh to regenerate the struct with original nservers.
        this.model.refreshStruct(true);
        NetworkStruct sn = this.model.getStruct(true);
        resetRandomGeneratorSeed(options.seed);

        // Always use solver_nc_jointaggr (not the LD variant) to match MATLAB behavior.
        // solver_nc_jointaggr builds its own mu matrix from sn.nservers.
        SolverNCJointReturn result = solver_nc_jointaggr(sn, this.options);
        
        NCResult ncResult = (NCResult) this.result;
        ncResult.solver = this.name;
        ncResult.prob.logNormConstAggr = result.lG;
        ncResult.prob.joint = result.Pr;
        
        long endTimeMillis = System.nanoTime();
        double runtime = (double) (endTimeMillis - startTimeMillis) / 1000000000.0;
        ncResult.runtime = runtime;

        ProbabilityResult probResult = new ProbabilityResult(result.Pr);
        probResult.isAggregated = true;
        return probResult;
    }

    public NetworkStruct getStruct() {
        if (this.sn == null)
            this.sn = this.model.getStruct(false);
        return this.sn;
    }

    public void setStruct(NetworkStruct sn) {
        this.sn = sn;
    }

    @Override
    public void runAnalyzer() throws IllegalAccessException {
        double start = System.nanoTime();
        if (this.model == null)
            throw new RuntimeException("Model is not provided");
        if (this.sn == null)
            this.sn = this.model.getStruct(false);
        if (this.options == null)
            this.options = new NCOptions();

        this.runAnalyzerChecks(options);
        this.resetRandomGeneratorSeed(options.seed);
        // Save original rt before refreshChains/cache analyzers modify it in-place.
        // In MATLAB, sn is a value-copy struct unaffected by refreshChains;
        // in Java, sn is a reference, so we must save/restore manually.
        Matrix rtOrig = sn.rt != null ? sn.rt.copy() : null;
        String origMethod = options.method;

        if (this.enableChecks && !this.supports(this.model)) {
            line_error(mfilename(new Object() {
            }), "This model contains features not supported by the solver.");
        }
        line_debug(options.verbose, String.format("NC solver starting: method=%s, nstations=%d, nclasses=%d", 
            options.method, sn.nstations, sn.nclasses));

        NCResult ret = null;
        String actualMethod = options.method;
        int iter = 0;

        // Method Selection and Preprocessing
        switch (options.method) {
            case "default":
                // Match MATLAB: any(sn.nservers(isfinite(sn.nservers))>1)
                // Only check finite servers for multi-server detection (exclude INF/Delay servers)
                boolean hasFiniteMultiServer = false;
                for (int i = 0; i < sn.nstations; i++) {
                    double ns = sn.nservers.get(i);
                    if (Double.isFinite(ns) && ns > 1) {
                        hasFiniteMultiServer = true;
                        break;
                    }
                }
                if (sn.nstations == 2 && !sn.nodetype.contains(NodeType.Cache) &&
                        sn.nodetype.contains(NodeType.Delay) && hasFiniteMultiServer) {
                    options.method = "comomld";
                }
                break;
            case "exact":
                if (!this.model.hasProductFormSolution()) {
                    line_error(mfilename(new Object(){}), "The exact method requires the model to have a product-form solution. This model does not have one. You can use Network.hasProductFormSolution() to check before running the solver.");
                } else if (sn.lldscaling == null || sn.lldscaling.isEmpty()) {
                    double Nt = sn.njobs.elementSum();
                    if (!Utils.isInf(Nt) && !Double.isNaN(Nt)) {
                        sn.lldscaling = Matrix.ones(sn.nstations, (int) Nt);
                        for (int i = 0; i < sn.nstations; i++) {
                            if (sn.nservers.get(i) > 1 && !Utils.isInf(sn.nservers.get(i))) {
                                for (int j = 0; j < Nt; j++) {
                                    sn.lldscaling.set(i, j, FastMath.min(j + 1, sn.nservers.get(i)));
                                }
                                sn.nservers.set(i, 1);
                            }
                        }
                    }
                }
                break;
        }

        // Check for non-reentrant cache model
        List<NodeType> nonReentrant = new ArrayList<>(Arrays.asList(NodeType.Source, NodeType.Cache, NodeType.Sink));
        if (sn.nclosedjobs == 0 && sn.nodetype.size() == 3 && sn.nodetype.containsAll(nonReentrant)) {
            // Initialize cache nodes
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.nodetype.get(ind) == NodeType.Cache) {
                    Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                    Matrix hitClass = cacheNode.getHitClass();
                    Matrix prob = hitClass.copy();
                    for (int i = 0; i < prob.length(); i++) {
                        if (prob.get(i) > 0) {
                            prob.set(i, 0.5);
                        }
                    }
                    cacheNode.setResultHitProb(prob);
                    Matrix missProb = prob.copy();
                    for (int i = 0; i < missProb.length(); i++) {
                        missProb.set(i, 1 - prob.get(i));
                    }
                    cacheNode.setResultMissProb(missProb);
                }
            }
            this.model.refreshChains(true);
            ret = solver_nc_cache_analyzer(this.sn, this.options.copy());
            actualMethod = ret.method;
            
            // Store item probabilities in result
            NCResult ncResult = (NCResult) this.result;
            if (ret.pij != null) {
                ncResult.prob.itemProb = ret.pij;
            }

            // Update hit probabilities based on results
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                    Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                    Matrix hitClass = cacheNode.getHitClass();
                    Matrix missClass = cacheNode.getMissClass();
                    Matrix hitProb = new Matrix(1, hitClass.length());
                    for (int k = 0; k < hitClass.length(); k++) {
                        // Find chain containing class k
                        int chainK = -1;
                        for (int c = 0; c < sn.chains.getNumRows(); c++) {
                            if (sn.chains.get(c, k) > 0) {
                                chainK = c;
                                break;
                            }
                        }
                        if (chainK >= 0) {
                            int h = (int) hitClass.get(k);
                            int m = (int) missClass.get(k);
                            if (h > 0 && m > 0) {
                                double sumXN = 0;
                                for (int j = 0; j < sn.nclasses; j++) {
                                    if (sn.chains.get(chainK, j) > 0) {
                                        sumXN += ret.XN.get(j);
                                    }
                                }
                                if (sumXN > 0) {
                                    hitProb.set(k, ret.XN.get(h) / sumXN);
                                }
                            }
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
            this.model.refreshChains(true);
        } else {
            // Regular queueing network
            if (sn.nodetype.contains(NodeType.Cache)) {
                // Cache-queueing network
                ret = solver_nc_cache_qn_analyzer(this.sn, this.options.copy());
                actualMethod = ret.method;
                iter = ret.iter;
                for (int ind = 0; ind < sn.nnodes; ind++) {
                    if (sn.nodetype.get(ind) == NodeType.Cache) {
                        Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                        // Extract hit/miss probabilities for this node
                        Matrix hitProb = new Matrix(1, ret.hitProb.getNumCols());
                        Matrix missProb = new Matrix(1, ret.missProb.getNumCols());
                        for (int j = 0; j < ret.hitProb.getNumCols(); j++) {
                            hitProb.set(j, ret.hitProb.get(ind, j));
                            missProb.set(j, ret.missProb.get(ind, j));
                        }
                        cacheNode.setResultHitProb(hitProb);
                        cacheNode.setResultMissProb(missProb);
                    }
                }
                this.model.refreshChains(true);
            } else {
                // Ordinary queueing network
                // Check for open model with single FCR containing single Delay (loss network)
                if (!snHasClosedClasses(sn) && sn.nregions == 1) {
                    Matrix regionMatrix = sn.region.get(0);
                    // Find stations in FCR (those with non-negative constraints)
                    int stationInFCR = -1;
                    int stationCount = 0;
                    for (int i = 0; i < sn.nstations; i++) {
                        boolean hasConstraint = false;
                        // Check per-class constraints
                        for (int r = 0; r < sn.nclasses; r++) {
                            if (regionMatrix.get(i, r) >= 0) {
                                hasConstraint = true;
                                break;
                            }
                        }
                        // Check global constraint
                        if (regionMatrix.get(i, sn.nclasses) >= 0) {
                            hasConstraint = true;
                        }
                        if (hasConstraint) {
                            stationInFCR = i;
                            stationCount++;
                        }
                    }
                    if (stationCount == 1 && Utils.isInf(sn.nservers.get(stationInFCR))) {
                        // Single delay node in FCR - check drop rule
                        if (sn.regionrule.get(0) == DropStrategy.Drop.getID()) {
                            // Use loss network solver
                            ret = solver_nc_lossn_analyzer(this.sn, this.options.copy());
                            actualMethod = ret.method;
                            iter = ret.iter;
                        } else {
                            // WAITQ (blocking) not supported - error and stop
                            throw new RuntimeException("SolverNC does not support finite capacity regions with WAITQ (blocking) policy. Use DROP policy instead.");
                        }
                    }
                }
                if (ret == null && ((sn.lldscaling != null && !sn.lldscaling.isEmpty()) || (sn.cdscaling != null && !sn.cdscaling.isEmpty()))) {
                    ret = solver_ncld_analyzer(this.sn, this.options.copy());
                    actualMethod = ret.method;
                    iter = ret.iter;
                } else if (ret == null) {
                    switch (options.method) {
                        case "exact":
                            if (!this.model.hasOpenClasses()) {
                                ret = solver_ncld_analyzer(this.sn, this.options.copy());
                                actualMethod = ret.method;
                                iter = ret.iter;
                            } else {
                                ret = solver_nc_analyzer(this.sn, this.options.copy());
                                actualMethod = ret.method;
                                iter = ret.iter;
                            }
                            break;
                        case "rd":
                        case "nrp":
                        case "nrl":
                        case "comomld":
                            ret = solver_ncld_analyzer(this.sn, this.options.copy());
                            actualMethod = ret.method;
                            iter = ret.iter;
                            break;
                        default:
                            line_debug(options.verbose, String.format("Using standard NC method: %s, calling solver_nc_analyzer", options.method));
                            ret = solver_nc_analyzer(this.sn, this.options.copy());
                            actualMethod = ret.method;
                            iter = ret.iter;
                            break;
                    }
                }
            }
        }

        // Propagate actual hit/miss probabilities from Cache model nodes to sn.nodeparam.
        // NC solver stores these on Cache nodes (via setResultHitProb/setResultMissProb)
        // but snGetArvRFromTput reads from sn.nodeparam.actualhitprob/actualmissprob.
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.nodetype.get(ind) == NodeType.Cache) {
                Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                CacheNodeParam cacheParam = (CacheNodeParam) sn.nodeparam.get(cacheNode);
                if (cacheParam != null) {
                    Matrix hitProb = cacheNode.getHitRatio();
                    Matrix missProb = cacheNode.getMissRatio();
                    if (hitProb != null && !hitProb.isEmpty()) {
                        cacheParam.actualhitprob = hitProb;
                    }
                    if (missProb != null && !missProb.isEmpty()) {
                        cacheParam.actualmissprob = missProb;
                    }
                }
            }
        }

        // Compute arrival rates
        AvgHandle T = getAvgTputHandles();
        if (rtOrig != null) {
            sn.rt = rtOrig;
        }
        Matrix AN = snGetArvRFromTput(sn, ret.TN, T);

        double finish = System.nanoTime();
        ret.runtime = (finish - start) / 1000000000.0;

        // Set results
        String resultMethod = actualMethod;
        if (origMethod.equals("default") && !actualMethod.equals("default")) resultMethod = "default/" + actualMethod;

        this.setAvgResults(ret.QN, ret.UN, ret.RN, ret.TN, AN, new Matrix(0,0), ret.CN, ret.XN, ret.runtime, resultMethod, iter);

        // Store probability results
        NCResult ncResult = (NCResult) this.result;
        ncResult.prob.logNormConstAggr = ret.lG;
    }

    /**
     * Get the normalizing constant and its logarithm
     *
     * @return NormalizingConstantResult containing the normalizing constant and its logarithm
     * @throws IllegalAccessException if analysis fails
     */
    public NormalizingConstantResult getNormalizingConstant() throws IllegalAccessException {
        if (GlobalConstants.DummyMode) {
            return new NormalizingConstantResult(Double.NaN, Double.NaN);
        }

        double lNormConst = getProbNormConstAggr().logNormalizingConstant;
        double normConst = FastMath.exp(lNormConst);
        return new NormalizingConstantResult(normConst, lNormConst);
    }

    /**
     * List all valid solution methods for this solver
     *
     * @return array of valid method names
     */
    public String[] listValidMethods() {
        return new String[]{
            "default", "exact", "imci", "ls", "le", "mmint2", "gleint",
            "panacea", "ca", "kt", "sampling", "propfair", "comomrm", "cub",
            "rd", "nrp", "nrl", "gm", "mem"
        };
    }

    /**
     * Maximum Entropy algorithm for Open Queueing Networks.
     *
     * Applies the ME algorithm from Kouvatsos (1994) to the model.
     * Only supports open queueing networks (no closed classes).
     *
     * @return NCResult containing ME algorithm results
     */
    public NCResult meOpen() {
        return meOpen(new jline.api.nc.MeOqnOptions());
    }

    /**
     * Maximum Entropy algorithm for Open Queueing Networks with custom options.
     *
     * Applies the ME algorithm from Kouvatsos (1994) to the model.
     * Only supports open queueing networks (no closed classes).
     *
     * @param meOptions Options for the ME algorithm
     * @return NCResult containing ME algorithm results
     */
    public NCResult meOpen(jline.api.nc.MeOqnOptions meOptions) {
        NetworkStruct sn = getStruct();

        // Check if model is open
        if (!model.hasOpenClasses()) {
            line_error(mfilename(new Object() {}), "meOpen only supports open queueing networks.");
            return null;
        }

        // Check for closed classes
        if (model.hasClosedClasses()) {
            line_error(mfilename(new Object() {}), "meOpen does not support models with closed classes.");
            return null;
        }

        int M = sn.nstations;
        int R = sn.nclasses;

        // Initialize parameter matrices
        Matrix lambda0 = new Matrix(M, R);
        Matrix Ca0 = new Matrix(M, R);
        Matrix mu = new Matrix(M, R);
        Matrix Cs = new Matrix(M, R);
        Matrix[][] P = new Matrix[M][M];

        // Initialize P array
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < M; i++) {
                P[j][i] = new Matrix(R, 1);
            }
        }

        // Extract service rates and service SCVs
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (sn.rates.get(i, r) > 0) {
                    mu.set(i, r, sn.rates.get(i, r));

                    // Extract service SCV from service process
                    if (sn.proc != null && i < sn.proc.size() && sn.proc.get(i) != null &&
                            r < sn.proc.get(i).size() && sn.proc.get(i).get(r) != null) {
                        jline.util.matrix.MatrixCell MAP = sn.proc.get(i).get(r);
                        if (!MAP.isEmpty()) {
                            Cs.set(i, r, jline.api.mam.Map_scvKt.map_scv(MAP.get(0), MAP.get(1)));
                        } else {
                            Cs.set(i, r, 1.0);
                        }
                    } else {
                        Cs.set(i, r, 1.0);
                    }
                }
            }
        }

        // Extract external arrival rates and arrival SCVs
        int sourceIdx = -1;
        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.nodetype.get(i) == jline.lang.constant.NodeType.Source) {
                sourceIdx = (int) sn.nodeToStation.get(i);
                break;
            }
        }

        if (sourceIdx >= 0) {
            for (int r = 0; r < R; r++) {
                if (sn.rates.get(sourceIdx, r) > 0) {
                    double externalRate = sn.rates.get(sourceIdx, r);
                    double Ca_external = 1.0;

                    // Extract arrival SCV from arrival process at source
                    if (sn.proc != null && sourceIdx < sn.proc.size() && sn.proc.get(sourceIdx) != null &&
                            r < sn.proc.get(sourceIdx).size() && sn.proc.get(sourceIdx).get(r) != null) {
                        jline.util.matrix.MatrixCell MAP = sn.proc.get(sourceIdx).get(r);
                        if (!MAP.isEmpty()) {
                            Ca_external = jline.api.mam.Map_scvKt.map_scv(MAP.get(0), MAP.get(1));
                        }
                    }

                    // Distribute external arrivals to queues based on routing
                    for (int i = 0; i < M; i++) {
                        if (i != sourceIdx) {
                            if (sn.visits != null && r < sn.visits.size() && sn.visits.get(r) != null &&
                                    i < sn.visits.get(r).length() && sn.visits.get(r).get(i) > 0) {
                                Matrix rt = sn.rt;
                                int sourceNode = (int) sn.stationToNode.get(sourceIdx);
                                int queueNode = (int) sn.stationToNode.get(i);

                                if (rt != null && sourceNode >= 0 && queueNode >= 0) {
                                    int routeIdx = sourceNode * R + r;
                                    int destIdx = queueNode * R + r;
                                    if (routeIdx < rt.getNumRows() && destIdx < rt.getNumCols()) {
                                        double routeProb = rt.get(routeIdx, destIdx);
                                        if (routeProb > 0) {
                                            lambda0.set(i, r, externalRate * routeProb);
                                            Ca0.set(i, r, Ca_external);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Extract routing probabilities between queues
        Matrix rt = sn.rt;
        if (rt != null) {
            for (int r = 0; r < R; r++) {
                for (int j = 0; j < M; j++) {
                    for (int i = 0; i < M; i++) {
                        int jNode = (int) sn.stationToNode.get(j);
                        int iNode = (int) sn.stationToNode.get(i);

                        if (jNode >= 0 && iNode >= 0) {
                            if (sn.nodetype.get(jNode) != jline.lang.constant.NodeType.Source &&
                                    sn.nodetype.get(jNode) != jline.lang.constant.NodeType.Sink &&
                                    sn.nodetype.get(iNode) != jline.lang.constant.NodeType.Source &&
                                    sn.nodetype.get(iNode) != jline.lang.constant.NodeType.Sink) {

                                int routeIdx = jNode * R + r;
                                int destIdx = iNode * R + r;

                                if (routeIdx < rt.getNumRows() && destIdx < rt.getNumCols()) {
                                    P[j][i].set(r, 0, rt.get(routeIdx, destIdx));
                                }
                            }
                        }
                    }
                }
            }
        }

        // Call ME algorithm
        jline.api.nc.MeOqnResult meResult = jline.api.nc.Me_oqnKt.me_oqn(M, R, lambda0, Ca0, mu, Cs, P, meOptions);

        // Build result
        NCResult result = new NCResult();
        result.QN = meResult.getL();
        result.UN = meResult.getRho();
        result.RN = meResult.getW();
        result.TN = meResult.getLambda();
        result.method = "mem";

        return result;
    }

    /**
     * Sample node state trajectory
     *
     * @param node The node to sample
     * @param numSamples Number of samples to generate
     * @return Sample result containing state trajectory
     */
    public Object sample(Node node, int numSamples) {
        if (GlobalConstants.DummyMode) {
            return null;
        }
        
        // NC solver doesn't support sampling in the same way as SSA
        // This method would need specialized implementation for NC methods
        line_error(mfilename(new Object(){}), "sample is not available in SolverNC. Use SolverSSA for sampling capabilities.");
        return null;
    }

    /**
     * Sample node state trajectory with default number of samples
     *
     * @param node The node to sample  
     * @return Sample result containing state trajectory
     */
    public Object sample(Node node) {
        return sample(node, options.samples);
    }

    /**
     * Checks whether the given model is supported by the NC solver
     *
     * @param model - the network model
     * @return - true if the model is supported, false otherwise
     */
    @Override
    public boolean supports(Network model) {
        return supportsModel(model);
    }

    /**
     * Static method to check whether the given model is supported by the NC solver.
     * This allows checking support without creating a solver instance.
     *
     * @param model - the network model
     * @return - true if the model is supported, false otherwise
     */
    public static boolean supportsModel(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverNC.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * Get cumulative distribution function of response times at FCFS and delay nodes
     * 
     * @param R Optional response time handles (currently unused in this implementation)
     * @return Response time distribution matrix for each node and class
     */
    public Matrix getCdfRespT(AvgHandle... R) {
        if (GlobalConstants.DummyMode) {
            return new Matrix(0, 0);
        }

        long startTimeMillis = System.nanoTime();
        NetworkStruct sn = getStruct();
        
        // Get algorithm configuration
        String algorithm = options.method != null ? options.method : "exact";
        
        Matrix RD = new Matrix(0, 0);
        
        try {
            // Get product form parameters
            Ret.snGetProductFormParams params = snGetProductFormParams(sn);
            Matrix D = params.D;  // Service demands
            Matrix N = params.N;  // Population vector  
            Matrix Z = params.Z;  // Think times
            Matrix S = params.S;  // Number of servers
            
            // Find FCFS and delay nodes
            List<Integer> fcfsNodesList = new ArrayList<>();
            List<Integer> fcfsNodeIdsList = new ArrayList<>(); 
            List<Integer> delayNodeIdsList = new ArrayList<>();
            
            for (int i = 0; i < sn.sched.size(); i++) {
                if (sn.sched.get(i) == SchedStrategy.FCFS) {
                    fcfsNodesList.add(i);
                    fcfsNodeIdsList.add(i);
                } else if (sn.sched.get(i) == SchedStrategy.INF) {
                    delayNodeIdsList.add(i);
                }
            }
            
            if (!fcfsNodesList.isEmpty()) {
                // Calculate time horizon
                double totalPop = N.elementSum();
                
                // Extract rates for FCFS nodes
                Matrix fcfsRates = new Matrix(fcfsNodesList.size(), sn.nclasses);
                for (int i = 0; i < fcfsNodesList.size(); i++) {
                    for (int j = 0; j < sn.nclasses; j++) {
                        fcfsRates.set(i, j, sn.rates.get(fcfsNodesList.get(i), j));
                    }
                }
                
                // Calculate mean service time across FCFS nodes
                double meanServiceTime = 0.0;
                int count = 0;
                for (int i = 0; i < fcfsRates.getNumRows(); i++) {
                    for (int j = 0; j < fcfsRates.getNumCols(); j++) {
                        if (fcfsRates.get(i, j) > 0) {
                            meanServiceTime += 1.0 / fcfsRates.get(i, j);
                            count++;
                        }
                    }
                }
                if (count > 0) {
                    meanServiceTime /= count;
                }
                
                double T = totalPop * meanServiceTime;
                
                // Create logarithmic time set: logspace(0, 2*log10(T), 100)
                Matrix tset = new Matrix(1, 100);
                for (int i = 0; i < 100; i++) {
                    double logVal = (2.0 * Math.log10(T)) * i / 99.0;
                    tset.set(i, Math.pow(10, logVal));
                }
                
                // Convert lists to matrices for function calls
                Matrix fcfsNodes = new Matrix(fcfsNodesList.size(), 1);
                for (int i = 0; i < fcfsNodesList.size(); i++) {
                    fcfsNodes.set(i, fcfsNodesList.get(i).doubleValue());
                }
                
                // Call appropriate PFQN algorithm
                Matrix[][] RDout;
                switch (algorithm) {
                    case "exact":
                        RDout = pfqn_stdf(D, N, Z, S, fcfsNodes, fcfsRates, tset);
                        break;
                    case "rd": 
                        RDout = pfqn_stdf_heur(D, N, Z, S, fcfsNodes, fcfsRates, tset);
                        break;
                    default:
                        RDout = pfqn_stdf(D, N, Z, S, fcfsNodes, fcfsRates, tset);
                        break;
                }
                
                // Initialize result matrix
                RD = new Matrix(sn.nnodes, sn.nclasses);
                
                // Process FCFS node results - remove complex number round-offs
                if (RDout != null && RDout.length > 0) {
                    for (int i = 0; i < fcfsNodeIdsList.size(); i++) {
                        for (int j = 0; j < sn.nclasses; j++) {
                            // RDout contains distribution data - for now just store a representative value
                            // In MATLAB this would be a cell array with full distribution data
                            if (i < RDout.length && j < RDout[i].length && RDout[i][j] != null) {
                                // Use the first element of the distribution matrix as representative value
                                if (RDout[i][j].getNumRows() > 0 && RDout[i][j].getNumCols() > 0) {
                                    RD.set(fcfsNodeIdsList.get(i), j, Math.abs(RDout[i][j].get(0, 0)));
                                }
                            }
                        }
                    }
                }
                
                // Process delay node results
                for (int i = 0; i < delayNodeIdsList.size(); i++) {
                    int nodeId = delayNodeIdsList.get(i);
                    for (int j = 0; j < sn.nclasses; j++) {
                        // For delay nodes, compute CDF from the process distribution
                        if (nodeId < sn.proc.size() && j < sn.proc.get(nodeId).size()) {
                            MatrixCell procCell = (MatrixCell) sn.proc.get(nodeId).get(j);
                            if (procCell != null && !procCell.isEmpty()) {
                                // Compute map_cdf for this process at time points
                                Matrix cdfResult = map_cdf(procCell, tset.transpose());
                                // Store representative value (first CDF value)
                                if (cdfResult.getNumRows() > 0) {
                                    RD.set(nodeId, j, cdfResult.get(0, 0));
                                }
                            }
                        }
                    }
                }
                
                long endTimeMillis = System.nanoTime();
                double runtime = (endTimeMillis - startTimeMillis) / 1000000000.0;
                
                // Store results using inherited method
                setDistribResults(RD, runtime);
                
            } else {
                line_warning(mfilename(new Object(){}), "getCdfRespT applies only to FCFS nodes.");
            }
            
        } catch (Exception e) {
            line_error(mfilename(new Object(){}), "Error in getCdfRespT: " + e.getMessage());
        }
        
        return RD;
    }

    /**
     * Get cumulative distribution function of response times with default parameters
     * 
     * @return Response time distribution matrix for each node and class
     */
    public DistributionResult getCdfRespT() {
        Matrix result = getCdfRespT((AvgHandle[]) null);
        NetworkStruct sn = getStruct();
        DistributionResult distResult = new DistributionResult((int) sn.nnodes, (int) sn.nclasses, "response_time");
        // Store the matrix result in the distribution result structure
        return distResult;
    }

    public static class SolverNCMargReturn {
        public Matrix lPr;
        public double G;
        public double lG;
        public double runtime;

        public SolverNCMargReturn(Matrix lPr, double G, double lG, double runtime) {
            this.lPr = lPr;
            this.G = G;
            this.lG = lG;
            this.runtime = runtime;
        }
    }

    public static class SolverNCJointReturn {
        public double Pr;
        public double G;
        public double lG;
        public double runtime;

        public SolverNCJointReturn(double Pr, double G, double lG, double runtime) {
            this.Pr = Pr;
            this.G = G;
            this.lG = lG;
            this.runtime = runtime;
        }
    }

    public static class SolverNCReturn {
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public Matrix T;
        public int C;
        public Matrix X;
        public double lG;
        public Matrix STeff;
        public int it;
        public String method;
        public double runtime;

        public SolverNCReturn(Matrix Q, Matrix U, Matrix R, Matrix T, int C,
                              Matrix X, Double lG, Matrix STeff, int it, double runtime, String method) {
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.T = T;
            this.C = C;
            this.X = X;
            this.lG = lG;
            this.STeff = STeff;
            this.it = it;
            this.method = method;
        }
    }

    public static class SolverNCLDReturn {
        public Matrix Q;
        public Matrix U;
        public Matrix R;
        public Matrix T;
        public Matrix C;
        public Matrix X;
        public double lG;
        public double runtime;
        public int it;
        public String method;

        public SolverNCLDReturn(Matrix Q, Matrix U, Matrix R, Matrix T, Matrix C,
                                Matrix X, double lG, double runtime, int it, String method) {
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.T = T;
            this.C = C;
            this.X = X;
            this.lG = lG;
            this.runtime = runtime;
            this.it = it;
            this.method = method;
        }
    }

    public static class NormalizingConstantResult {
        public double normConst;
        public double lNormConst;

        public NormalizingConstantResult(double normConst, double lNormConst) {
            this.normConst = normConst;
            this.lNormConst = lNormConst;
        }
    }

}
