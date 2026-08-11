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
import jline.lang.NodeParam;
import jline.lang.constant.SchedStrategy;
import jline.api.sn.SnHasProductForm;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.state.FromMarginal;
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

import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import jline.solvers.fj.FJFixedPoint;
import jline.solvers.mva.MVAResult;
import static jline.api.sn.SnHasMultiServer.snHasMultiServer;
import static jline.io.InputOutput.*;
import static jline.solvers.nc.analyzers.Solver_nc_analyzer.solver_nc_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_cache_analyzer.solver_nc_cache_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_cache_qn_analyzer.solver_nc_cache_qn_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_retrieval_analyzer.solver_nc_retrieval_analyzer;
import static jline.solvers.nc.analyzers.Solver_ncld_analyzer.solver_ncld_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_lossn_analyzer.solver_nc_lossn_analyzer;
import static jline.api.sn.SnHasClosedClasses.snHasClosedClasses;
import jline.lang.constant.DropStrategy;
import jline.solvers.nc.handlers.Solver_nc_pas_is;
import static jline.solvers.nc.handlers.Solver_nc_marg.solver_nc_marg;
import static jline.solvers.nc.handlers.Solver_nc_joint.solver_nc_joint;
import static jline.solvers.nc.handlers.Solver_nc_margaggr.solver_nc_margaggr;
import static jline.solvers.nc.handlers.Solver_nc_jointaggr.solver_nc_jointaggr;
import static jline.solvers.nc.handlers.Solver_nc_jointaggr_ld.solver_nc_jointaggr_ld;
import static jline.api.sn.SnGetProductFormParams.snGetProductFormParams;
import jline.io.Ret;
import static jline.api.pfqn.nc.Pfqn_stdf.pfqn_stdf;
import static jline.api.pfqn.nc.Pfqn_stdf_heur.pfqn_stdf_heur;
import static jline.api.pfqn.nc.Pfqn_procomom.pfqn_procomom;
import static jline.api.sn.SnGetDemandsChain.snGetDemandsChain;
import static jline.api.mam.Map_cdf.map_cdf;


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

    /**
     * The normalizing-constant solver is exact on the same product-form class that
     * pfqn_sens differentiates, so getSensitivityTable uses the analytic branch.
     *
     * @return true
     */
    @Override
    public boolean supportsExactSensitivity() {
        return true;
    }

    /** True if the model contains a Cache equipped with a delayed-hit retrieval system. */
    public static boolean hasRetrievalCache(NetworkStruct sn) {
        if (sn.nodeparam == null) return false;
        for (NodeParam np : sn.nodeparam.values()) {
            if (np instanceof CacheNodeParam && ((CacheNodeParam) np).retrievalSystemCapacity > 0) return true;
        }
        return false;
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
                // Finite capacity regions: NC solves the OPEN single-Delay
                // loss-network case exactly (Erlang fixed point,
                // Solver_nc_lossn_analyzer). It cannot do an FCR on queueing
                // stations -- a boolean feature cannot express that split, so
                // runAnalyzer keeps its residual-FCR check for that case.
                "Region",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS", "SchedStrategy_SIRO",
                "SchedStrategy_LCFS", "SchedStrategy_LCFSPR", "SchedStrategy_OI",
                "SchedStrategy_PAS",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "SchedStrategy_FCFS", "ClosedClass", "SelfLoopingClass",
                "Cache", "CacheClassSwitcher", "CacheRetrieval", "OpenClass",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO",
                "ReplacementStrategy_HLRU",
                "LoadDependence", "ClassDependence", "JointDependence",
                // Fork-join through the MMT/HT transformation, driven by
                // jline.solvers.fj.FJFixedPoint (as in SolverMVA)
                "Fork", "Forker", "Join", "Joiner"
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
        // getStruct(true) syncs node states into sn.state (a plain getStruct
        // would query the PREVIOUS state, giving a one-behind lag).
        NetworkStruct sn = this.model.getStruct(true);
        int ist = (int) sn.nodeToStation.get(node.getNodeIndex());
        int isf = (int) sn.nodeToStateful.get(node.getNodeIndex());

        Matrix state_new;
        if (state_a == null || state_a.isEmpty()) {
            state_new = sn.state.get((StatefulNode) this.model.getNodes().get((int) sn.nodeToStateful.get(node.getNodeIndex())));
        } else {
            // state_a is a per-class marginal count vector; encode it into the
            // internal state under the STATEFUL index expected by
            // solver_nc_margaggr (raw counts would be misinterpreted).
            state_new = FromMarginal.fromMarginal(sn, node.getNodeIndex(), state_a);
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
     * Get marginal queue-length probability distribution at a node.
     *
     * Returns P(n total jobs) for n=0,1,...,N, summing over all class combinations.
     * When method is "comom", uses pfqn_procomom directly for efficiency.
     * Otherwise falls back to enumeration via getProbAggr.
     *
     * @param node The node to compute marginal probability for
     * @return Matrix of marginal probabilities where element j = P(j total jobs at node)
     */
    /**
     * Marginal queue-length distribution at a station, routing the CLI's
     * (station, class) call to the exact procomom/enumeration path in
     * {@link #getProbMarg(Node)}. The class index is not used: NC returns the
     * per-station total-jobs marginal.
     *
     * @param ist      station index
     * @param jobclass class index (unused; marginal is over total jobs)
     * @param state_m  optional state (unused)
     * @return the exact marginal probability distribution
     */
    @Override
    public ProbabilityResult getProbMarg(int ist, int jobclass, Matrix state_m) {
        NetworkStruct sn = getStruct();
        int nodeIdx = (int) sn.stationToNode.get(ist);
        return getProbMarg(model.getNodes().get(nodeIdx));
    }

    public ProbabilityResult getProbMarg(Node node) {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult();
        }

        NetworkStruct sn = getStruct();
        int ist = (int) sn.nodeToStation.get(node.getNodeIndex());

        // Check that all classes are closed
        Matrix N = sn.njobs;
        for (int r = 0; r < sn.nclasses; r++) {
            if (Utils.isInf(N.get(r))) {
                line_error(mfilename(new Object(){}), "getProbMarg not yet implemented for models with open classes.");
                return new ProbabilityResult();
            }
        }

        int R = sn.nclasses;
        int Ntotal = (int) N.elementSum();

        // Exact marginal via pfqn_procomom (default path, mirrored by native).
        if (true) {
            int M = sn.nstations;
            int C = sn.nchains;
            Matrix nservers = sn.nservers;

            Ret.snGetDemands ret = snGetDemandsChain(sn);
            Matrix Lchain = ret.Dchain;
            Matrix Nchain = ret.Nchain;

            // Replace non-finite values with 0
            for (int i = 0; i < Lchain.getNumRows(); i++) {
                for (int j = 0; j < Lchain.getNumCols(); j++) {
                    if (!Double.isFinite(Lchain.get(i, j))) {
                        Lchain.set(i, j, 0.0);
                    }
                }
            }

            // Separate queue vs delay stations (matching solver_nc.m / MATLAB getProbMarg)
            Matrix Lms = new Matrix(M, C);
            Lms.fill(0.0);
            Matrix Ztotal = new Matrix(1, C);
            Ztotal.fill(0.0);
            Matrix Zms = new Matrix(1, C);
            Zms.fill(0.0);
            List<Integer> queueStations = new ArrayList<>();

            for (int i = 0; i < M; i++) {
                if (Utils.isInf(nservers.get(i))) {
                    // Delay station: accumulate into Ztotal
                    for (int c = 0; c < C; c++) {
                        Ztotal.set(0, c, Ztotal.get(0, c) + Lchain.get(i, c));
                    }
                } else {
                    queueStations.add(i);
                    for (int c = 0; c < C; c++) {
                        Lms.set(i, c, Lchain.get(i, c) / nservers.get(i));
                        Zms.set(0, c, Zms.get(0, c) + Lchain.get(i, c) * (nservers.get(i) - 1) / nservers.get(i));
                    }
                }
            }

            // Build L_queues (only queue station rows from Lms)
            int Mq = queueStations.size();
            Matrix L_queues = new Matrix(Mq, C);
            for (int qi = 0; qi < Mq; qi++) {
                int origIdx = queueStations.get(qi);
                for (int c = 0; c < C; c++) {
                    L_queues.set(qi, c, Lms.get(origIdx, c));
                }
            }

            // Z_total = Ztotal + Zms
            Matrix Z_total = new Matrix(1, C);
            for (int c = 0; c < C; c++) {
                Z_total.set(0, c, Ztotal.get(0, c) + Zms.get(0, c));
            }

            // Call pfqn_procomom
            Ret.pfqnProcomom procomom = pfqn_procomom(L_queues, Nchain, Z_total);
            Matrix Pr = procomom.Pr;

            // Find queue index for the requested station
            int queueIdx = queueStations.indexOf(ist);
            if (queueIdx >= 0) {
                int sumNchain = (int) Nchain.elementSum();
                Matrix Pmarg = new Matrix(1, Ntotal + 1);
                Pmarg.fill(0.0);
                int len = Math.min(sumNchain + 1, Ntotal + 1);
                for (int j = 0; j < len; j++) {
                    Pmarg.set(0, j, Pr.get(queueIdx, j));
                }
                return new ProbabilityResult(Pmarg);
            } else {
                // Delay station: fall through to enumeration
                line_warning(mfilename(new Object(){}), "comom method does not directly support delay stations, using enumeration.");
                // Fall through to enumeration below
            }
        }

        // Enumeration-based fallback: sum getProbAggr over all class partitions
        Matrix Pmarg = new Matrix(1, Ntotal + 1);
        Pmarg.fill(0.0);

        for (int n = 0; n <= Ntotal; n++) {
            List<int[]> partitions = generatePartitions(n, R, N);
            double probN = 0.0;
            for (int[] partition : partitions) {
                Matrix stateVec = new Matrix(1, R);
                for (int r = 0; r < R; r++) {
                    stateVec.set(0, r, partition[r]);
                }
                Double prob = getProbAggr(node, stateVec);
                if (prob != null && prob > 0) {
                    probN += prob;
                }
            }
            Pmarg.set(0, n, probN);
        }

        // Normalize
        double totalProb = 0.0;
        for (int j = 0; j <= Ntotal; j++) {
            totalProb += Pmarg.get(0, j);
        }
        if (totalProb > 0 && FastMath.abs(totalProb - 1.0) > 1e-10) {
            for (int j = 0; j <= Ntotal; j++) {
                Pmarg.set(0, j, Pmarg.get(0, j) / totalProb);
            }
        }

        return new ProbabilityResult(Pmarg);
    }

    /**
     * Generate all partitions of n into R non-negative integers,
     * each bounded by the corresponding element of Nmax.
     */
    private List<int[]> generatePartitions(int n, int R, Matrix Nmax) {
        List<int[]> result = new ArrayList<>();
        generatePartitionsHelper(n, R, 0, Nmax, new int[R], result);
        return result;
    }

    private void generatePartitionsHelper(int remaining, int R, int idx, Matrix Nmax, int[] current, List<int[]> result) {
        if (idx == R - 1) {
            if (remaining <= (int) Nmax.get(idx)) {
                current[idx] = remaining;
                result.add(current.clone());
            }
            return;
        }
        int maxVal = (int) Math.min(remaining, Nmax.get(idx));
        for (int v = 0; v <= maxVal; v++) {
            current[idx] = v;
            generatePartitionsHelper(remaining - v, R, idx + 1, Nmax, current, result);
        }
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
        if (this.model == null)
            throw new RuntimeException("Model is not provided");
        if (this.sn == null)
            this.sn = this.model.getStruct(false);
        if (this.options == null)
            this.options = new NCOptions();

        // Finite station/class capacity: a product-form method has no
        // representation of a finite buffer, so it would silently return the
        // unconstrained answer (QLen=4 instead of the M/M/1/2 value 0.8525).
        // Same defect as BUG-39 on the MVA side.
        //
        // method='mem' is the one exception, and only on a single-class open
        // model that memUnsupportedReason has cleared: MEM does represent the
        // buffer, as a censored GE/GE/c/0;N queue, and handles both a lost
        // arrival (DROP) and a job held in the upstream server (BAS).
        if (!memFiniteBufferPath(this.sn, this.options)) {
            String capacityReason = NetworkSolver.bindingCapacityReason(this.model, this.sn, "SolverNC");
            if (capacityReason != null) {
                throw new RuntimeException(capacityReason);
            }
        }

        // MATLAB's runAnalyzer takes sn and options by value, so the
        // multiserver-to-lldscaling conversion and the method downgrade it may
        // perform are scoped to the call. Here both are references into state
        // shared with the model and outlive it, and SolverLN re-solves each layer
        // once per iteration: a leaked lldscaling made the next solve read it as
        // a user-supplied load dependence, misclassify the layer as
        // non-product-form, downgrade the method to 'comom', and then route that
        // method to the ncld analyzer, which does not implement it.
        Matrix lldOrig = this.sn.lldscaling != null ? this.sn.lldscaling.copy() : null;
        String methodOrig = this.options.method;
        try {
            runAnalyzerBody();
        } finally {
            this.sn.lldscaling = lldOrig;
            this.options.method = methodOrig;
        }
    }

    private void runAnalyzerBody() throws IllegalAccessException {
        double start = System.nanoTime();

        this.runAnalyzerChecks(options);
        if (sn.immfeed != null && sn.immfeed.elementSum() > 0) {
            line_warning(mfilename(new Object(){}), "SolverNC does not handle immediate feedback (immfeed); the solver will treat self-loops as class-switching with re-queueing.");
        }
        this.resetRandomGeneratorSeed(options.seed);
        // Save original rt before refreshChains/cache analyzers modify it in-place.
        // In MATLAB, sn is a value-copy struct unaffected by refreshChains;
        // in Java, sn is a reference, so we must save/restore manually.
        Matrix rtOrig = sn.rt != null ? sn.rt.copy() : null;
        String origMethod = options.method;

        // Warn (debug level) when the requested method is unrecognized. With solver
        // checks disabled (e.g. when SolverNC is used as a SolverLN layer backend,
        // where enableChecks=false bypasses runAnalyzerChecks) an unknown method name
        // silently falls through to the default normalizing-constant analyzer instead
        // of raising an error, which masks typos such as 'adaptive'.
        boolean knownMethod = "comomld".equals(origMethod); // comomld is auto-selected internally from 'default'
        for (String vm : listValidMethods()) {
            if (vm.equals(origMethod)) {
                knownMethod = true;
                break;
            }
        }
        if (!knownMethod && this.enableChecks) {
            // Mirrors MATLAB NetworkSolver.runAnalyzerChecks, which rejects any method
            // outside listValidMethods before applying the feature gate below. This
            // used to be a debug-level message only, so an unadvertised method (e.g.
            // 'mom', which SolverNC does not implement) fell through to the default
            // analyzer and returned ALL-ZERO queue lengths while still reporting a
            // completed analysis -- and NetworkAvgTable then hid the zero rows, so the
            // caller saw an empty table rather than an error. Refuse instead.
            line_error(mfilename(new Object() {
            }), String.format("The '%s' method is unsupported by this solver.", origMethod));
        }
        if (!knownMethod) {
            // Only reachable with enableChecks=false, e.g. when SolverNC is used as a
            // SolverLN layer backend, which bypasses the checks above by design.
            line_debug(options.verbose, String.format("NC: unrecognized method '%s', falling back to the default normalizing-constant analyzer (nc_analyzer/comom).", origMethod));
        }

        if (this.enableChecks) {
            // Method-aware feature gate: method='mem' is gated by its structural
            // applicability (memUnsupportedReason); all other methods use the
            // coarse product-form feature set. Mirrors Solver.supportsModelMethod.
            String reason = this.supportsModelMethod(options.method);
            if (!reason.isEmpty()) {
                line_error(mfilename(new Object() {
                }), "This model contains features not supported by the solver. " + reason);
                return;
            }
        }
        line_debug(options.verbose, String.format("NC solver starting: method=%s, nstations=%d, nclasses=%d", 
            options.method, sn.nstations, sn.nclasses));

        // Fork-join: the MMT/HT transformation rewrites the model as a plain
        // network (forks -> routers, joins -> zero-service delays, parallelism
        // carried by auxiliary open classes), which the normalizing-constant
        // analyzer solves directly. The fixed point itself is solver-agnostic
        // and lives in jline.solvers.fj.FJFixedPoint; NC drives it with
        // ncDispatch as the inner solve, exactly as SolverMVA drives it with
        // its analyzer dispatch. Mirrors matlab/src/solvers/NC/@SolverNC.
        boolean hasFork = false;
        for (NodeType nt : sn.nodetype) {
            if (nt == NodeType.Fork) {
                hasFork = true;
                break;
            }
        }
        if (hasFork) {
            FJFixedPoint.FJState fjState = new FJFixedPoint.FJState(null, null);
            FJFixedPoint.FJOutcome fjOut = FJFixedPoint.run(this.model, this.sn, this.options, fjState,
                    new FJFixedPoint.InnerSolve() {
                        @Override
                        public MVAResult solve(NetworkStruct snIn, SolverOptions opts) {
                            return SolverNC.this.ncDispatch(snIn, opts);
                        }
                    }, (long) start);
            MVAResult fjRet = fjOut.ret;
            // The driver leaves the transformed struct behind; the result tail
            // below indexes the ORIGINAL stations and classes, so recompile it.
            this.sn = this.model.getStruct(true);
            AvgHandle Tfj = getAvgTputHandles();
            Matrix ANfj = snGetArvRFromTput(this.sn, fjRet.TN, Tfj);
            double runtimeFj = (System.nanoTime() - start) / 1000000000.0;
            this.setAvgResults(fjRet.QN, fjRet.UN, fjRet.RN, fjRet.TN, ANfj, new Matrix(0, 0),
                    fjRet.CN, fjRet.XN, runtimeFj, options.method, fjOut.iter);
            ((NCResult) this.result).prob.logNormConstAggr = fjRet.logNormConstAggr;
            return;
        }

        NCResult ret = null;
        String actualMethod = options.method;
        int iter = 0;

        // Maximum Entropy Method (Kouvatsos 1994): explicit request only. The
        // method='default' path routes to the native normalizing-constant
        // analyzer, as it did before MEM was introduced.
        boolean useMem = "mem".equals(options.method);
        if (useMem) {
            NCResult memRet;
            if (model.hasClosedClasses() && !model.hasOpenClasses()) {
                line_debug(options.verbose, "NC method=mem, routing to meClosed (Maximum Entropy, closed QN)");
                memRet = this.meClosed();
            } else if (model.hasClosedClasses() && model.hasOpenClasses()) {
                line_debug(options.verbose, "NC method=mem, routing to meMixed (Maximum Entropy, mixed QN)");
                memRet = this.meMixed();
            } else {
                line_debug(options.verbose, "NC method=mem, routing to meOpen (Maximum Entropy, open QN)");
                memRet = this.meOpen();
            }
            AvgHandle Tmem = getAvgTputHandles();
            Matrix ANmem = snGetArvRFromTput(sn, memRet.TN, Tmem);
            memRet.runtime = (System.nanoTime() - start) / 1000000000.0;
            this.setAvgResults(memRet.QN, memRet.UN, memRet.RN, memRet.TN, ANmem,
                    new Matrix(0, 0), memRet.CN, memRet.XN,
                    memRet.runtime, "mem", memRet.iter);
            return;
        }

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
                    // 2-station Delay+multiserver (the topology of every SolverLN
                    // layer submodel). Product-form models are solved exactly via
                    // load-dependent CoMoM (convert the multiserver to
                    // mu(n)=min(n,c) and route to the ncld analyzer), which matches
                    // CTMC/MVA exactly. Non-product-form models (e.g. LN layers
                    // with heterogeneous per-class FCFS rates) have no exact
                    // CoMoM-LD solution and fall back to Seidmann's approximation
                    // (comom).
                    if (this.model.hasProductFormSolution() && (sn.lldscaling == null || sn.lldscaling.isEmpty())) {
                        double Nt = sn.njobs.elementSum();
                        if (!Utils.isInf(Nt) && !Double.isNaN(Nt)) {
                            sn.lldscaling = Matrix.ones(sn.nstations, (int) Nt);
                            for (int i = 0; i < sn.nstations; i++) {
                                if (sn.nservers.get(i) > 1 && !Utils.isInf(sn.nservers.get(i))) {
                                    // The queueing solve reads mu(n) from lldscaling, so the
                                    // server count must be kept: utilization is the fraction of
                                    // the c servers busy, and c is not recoverable from
                                    // lldscaling once Nt<c (min(1:Nt,c) is then just 1:Nt).
                                    // Zeroing it to 1 made Solver_ncld normalize U by
                                    // max(lldscaling)=min(Nt,c) instead of c, overstating U by
                                    // c/Nt whenever the population is below the servers.
                                    for (int j = 0; j < Nt; j++) {
                                        sn.lldscaling.set(i, j, FastMath.min(j + 1, sn.nservers.get(i)));
                                    }
                                }
                            }
                        }
                        line_debug(options.verbose, "NC: default method for 2-station multiserver Delay product-form network, using exact load-dependent comomld");
                    } else {
                        options.method = "comom";
                        line_debug(options.verbose, "NC: default method for 2-station multiserver Delay non-product-form network, switching to comom");
                    }
                }
                break;
            case "is":
                // 'is' (importance sampling) on an ORDINARY product-form model is a
                // normalizing-constant estimator and needs the same model as
                // "exact": convert multiserver stations to load-dependent here too,
                // so it routes to solver_ncld_analyzer -> Pfqn_ncld -> Pfqn_ld_is
                // rather than falling through to Seidmann's approximation, which
                // would faithfully estimate a DIFFERENT (approximated) model.
                //
                // OI / pass-and-swap models are the exception: an OI station is
                // multiserver but is not a plain min(n,c) load-dependent station,
                // and a P&S tandem has only per-communicating-class product form
                // (so hasProductFormSolution() is false). Leave them untouched for
                // Solver_nc to route to Pfqn_pas_is / Pfqn_oi_is.
                if (Solver_nc_pas_is.nc_is_pas_model(sn)) {
                    break;   // handled by Solver_nc (Pfqn_pas_is)
                }
                // fall through to the "exact" preprocessing
            case "panaceald":
                // 'panaceald' is a load-dependent normalizing-constant expansion
                // and needs the same multiserver conversion as "exact"
            case "exact":
                if (!this.model.hasProductFormSolution()) {
                    line_error(mfilename(new Object(){}), "The " + options.method + " method requires the model to have a product-form solution. This model does not have one. You can use Network.hasProductFormSolution() to check before running the solver.");
                } else if ((sn.lldscaling == null || sn.lldscaling.isEmpty()) && ncHasFiniteMultiserver(sn)) {
                    // Only convert to load-dependent when a genuine multiserver
                    // station is present. Setting lldscaling=ones on an
                    // all-single-server model is a semantic no-op but forces the
                    // ncld path (solver_ncld_analyzer), whose Pfqn_nc_sanitize
                    // crashes ("Outside of matrix bounds") on a single-station-
                    // confined (self-looping) closed chain. Mirrors the MATLAB
                    // @SolverNC/runAnalyzer.m fix.
                    double Nt = sn.njobs.elementSum();
                    if (!Utils.isInf(Nt) && !Double.isNaN(Nt)) {
                        // transform multi-server nodes into lld nodes; the server count is
                        // kept so that utilization stays normalized by c (see the
                        // "default" branch above)
                        sn.lldscaling = Matrix.ones(sn.nstations, (int) Nt);
                        for (int i = 0; i < sn.nstations; i++) {
                            if (sn.nservers.get(i) > 1 && !Utils.isInf(sn.nservers.get(i))) {
                                for (int j = 0; j < Nt; j++) {
                                    sn.lldscaling.set(i, j, FastMath.min(j + 1, sn.nservers.get(i)));
                                }
                            }
                        }
                    }
                }
                break;
        }

        // Check for delayed-hit retrieval cache model (Cache with a retrieval system) first
        List<NodeType> nonReentrant = new ArrayList<>(Arrays.asList(NodeType.Source, NodeType.Cache, NodeType.Sink));
        if (hasRetrievalCache(this.sn)) {
            boolean hasSource = false;
            for (int i = 0; i < this.sn.nodetype.size(); i++) if (this.sn.nodetype.get(i) == NodeType.Source) { hasSource = true; break; }
            if (hasSource) {
                line_debug(options.verbose, "NC: detected open delayed-hit retrieval cache, calling solver_nc_retrieval_analyzer");
                ret = solver_nc_retrieval_analyzer(this.sn, this.options.copy());
            } else {
                line_debug(options.verbose, "NC: detected closed integrated delayed-hit retrieval cache, calling solver_nc_cacheqn_retrieval_analyzer");
                ret = jline.solvers.nc.analyzers.Solver_nc_cacheqn_retrieval_analyzer.solver_nc_cacheqn_retrieval_analyzer(this.sn, this.options.copy());
            }
            actualMethod = ret.method;
        } else if (sn.nclosedjobs == 0 && sn.nodetype.size() == 3 && sn.nodetype.containsAll(nonReentrant)) {
            line_debug(options.verbose, "NC: detected non-reentrant cache model (Source-Cache-Sink), calling solver_nc_cache_analyzer");
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
                line_debug(options.verbose, "NC: detected cache-queueing network, calling solver_nc_cache_qn_analyzer");
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
                            line_debug(options.verbose, "NC: detected loss network (single Delay in FCR with Drop), calling solver_nc_lossn_analyzer");
                            ret = solver_nc_lossn_analyzer(this.sn, this.options.copy());
                            actualMethod = ret.method;
                            iter = ret.iter;
                        } else {
                            // WAITQ (blocking) not supported - error and stop
                            throw new RuntimeException("SolverNC does not support finite capacity regions with WAITQ (blocking) policy. Use DROP policy instead.");
                        }
                    }
                }
                // Residual FCR (not the single-Delay loss-network case dispatched
                // above): NC does not enforce the aggregate region limit and would
                // silently return the unconstrained answer. Reject.
                if (ret == null && sn.nregions > 0) {
                    throw new RuntimeException("This model uses a Finite Capacity Region (addRegion) "
                            + "on queueing stations, which is not supported by SolverNC (only the "
                            + "single-Delay loss-network case is). Use SolverJMT, or setCapacity "
                            + "for a single-station limit.");
                }
                if (ret == null && ((sn.lldscaling != null && !sn.lldscaling.isEmpty()) || (sn.cdscaling != null && !sn.cdscaling.isEmpty()) || (sn.jdscaling != null && !sn.jdscaling.isEmpty()))) {
                    line_debug(options.verbose, "NC: detected load-/class-/joint-dependent scaling, calling solver_ncld_analyzer");
                    ret = solver_ncld_analyzer(this.sn, this.options.copy());
                    actualMethod = ret.method;
                    iter = ret.iter;
                } else if (ret == null) {
                    switch (options.method) {
                        case "exact":
                            // This branch is reached only with EMPTY lldscaling
                            // (the test above routes any load-dependent model to
                            // ncld). Genuine multiserver stations were converted
                            // to lldscaling above, so an empty lldscaling here
                            // means an all-single-server closed model, which the
                            // standard normalizing-constant path solves exactly.
                            // Routing it to solver_ncld_analyzer instead crashes
                            // (Pfqn_nc_sanitize, "Outside of matrix bounds") on a
                            // single-station-confined (self-looping) closed chain.
                            // Mirrors the MATLAB @SolverNC/runAnalyzer.m fix.
                            line_debug(options.verbose, "NC: exact method for single-server closed/open model, calling solver_nc_analyzer");
                            ret = solver_nc_analyzer(this.sn, this.options.copy());
                            actualMethod = ret.method;
                            iter = ret.iter;
                            break;
                        case "rd":
                        case "nrp":
                        case "nrl":
                        case "comomld":
                        case "panaceald":
                            line_debug(options.verbose, String.format("NC: load-dependent method=%s, calling solver_ncld_analyzer", options.method));
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

        // Compute arrival rates
        AvgHandle T = getAvgTputHandles();
        Matrix rtRefreshed = sn.rt;
        if (rtOrig != null) {
            sn.rt = rtOrig;
        }
        Matrix AN = snGetArvRFromTput(sn, ret.TN, T);
        if (rtRefreshed != null) {
            sn.rt = rtRefreshed;
        }

        double finish = System.nanoTime();
        ret.runtime = (finish - start) / 1000000000.0;

        // Set results
        String resultMethod = actualMethod;
        if (origMethod.equals("default") && !actualMethod.equals("default")) resultMethod = "default/" + actualMethod;

        line_debug(options.verbose, String.format("NC solver completed: method=%s, runtime=%.3fs", resultMethod, ret.runtime));
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
            "default", "exact", "erlangfp", "mci", "imci", "ls", "le", "mmint2", "gleint",
            "panacea", "panaceald", "ca", "clw", "kt", "sampling", "is", "propfair", "comom", "cub",
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
     * True when any finite arrival or service scv in the model deviates from
     * one, i.e. the model carries non-Markovian variability that product-form
     * normalizing-constant methods would silently exponentialize.
     */
    private static boolean snHasNonUnitScv(NetworkStruct sn) {
        return snHasNonUnitScv(sn, false);
    }

    /** True if any station has a finite server count greater than 1. */
    private static boolean ncHasFiniteMultiserver(NetworkStruct sn) {
        for (int i = 0; i < sn.nstations; i++) {
            double c = sn.nservers.get(i);
            if (!Utils.isInf(c) && c > 1) {
                return true;
            }
        }
        return false;
    }

    private static boolean snHasNonUnitScv(NetworkStruct sn, boolean productForm) {
        if (sn.scv == null) {
            return false;
        }
        for (int i = 0; i < sn.scv.getNumRows(); i++) {
            if (productForm && i < sn.nstations && isInsensitiveStation(sn, i)) {
                continue; // insensitive station: service distribution irrelevant
            }
            for (int r = 0; r < sn.scv.getNumCols(); r++) {
                double v = sn.scv.get(i, r);
                if (!Double.isNaN(v) && !Double.isInfinite(v) && Math.abs(v - 1.0) > 1e-8) {
                    return true;
                }
            }
        }
        return false;
    }

    private static boolean isInsensitiveStation(NetworkStruct sn, int ist) {
        SchedStrategy sk = sn.sched.get(sn.stations.get(ist));
        return sk == SchedStrategy.INF || sk == SchedStrategy.PS || sk == SchedStrategy.LCFSPR;
    }

    /**
     * Returns null when the Maximum Entropy Method (Kouvatsos 1994) supports
     * the features of the model described by sn (node types, class switching,
     * scheduling, source presence); otherwise a message describing the first
     * unsupported feature. Open/closed class membership is checked separately
     * on the model; this overload assumes the open-network variant.
     */
    public static String memUnsupportedReason(NetworkStruct sn) {
        return memUnsupportedReason(sn, false);
    }

    /**
     * Returns null when the Maximum Entropy Method (Kouvatsos 1994) supports
     * the features of the model described by sn; otherwise a message
     * describing the first unsupported feature. Open models (Section 3.2)
     * allow Source, Queue, Delay and Sink nodes with GE/GE/1, GE/GE/c and
     * GE/GE/inf building blocks; closed models (Section 3.3) allow Queue and
     * Delay nodes with G/G/1 and G/G/inf building blocks only, so finite
     * multiserver stations are rejected.
     */
    public static String memUnsupportedReason(NetworkStruct sn, boolean closed) {
        return memUnsupportedReason(sn, !closed, closed);
    }

    /**
     * Feature check for the Maximum Entropy Method given the class
     * composition of the model: hasOpen/hasClosed flag the presence of
     * open and closed classes (both true for mixed models). Returns null
     * when supported, otherwise the first unsupported feature found.
     */
    public static String memUnsupportedReason(NetworkStruct sn, boolean hasOpen, boolean hasClosed) {
        for (int ind = 0; ind < sn.nnodes; ind++) {
            NodeType nt = sn.nodetype.get(ind);
            if (nt != NodeType.Source && nt != NodeType.Sink
                    && nt != NodeType.Queue && nt != NodeType.Delay) {
                return "MEM supports only Source, Queue, Delay and Sink nodes.";
            }
            if (!hasOpen && (nt == NodeType.Source || nt == NodeType.Sink)) {
                return "MEM supports only Queue and Delay nodes in closed models.";
            }
        }

        // Class switching is not part of the Kouvatsos (1994) network model
        if (sn.csmask != null) {
            for (int r = 0; r < sn.nclasses; r++) {
                for (int s = 0; s < sn.nclasses; s++) {
                    if (r != s && sn.csmask.get(r, s) > 0) {
                        return "MEM does not support class switching.";
                    }
                }
            }
        }

        // Absorbing self-loops (p_ii=1) make the routing reducible and the
        // geometric feedback transform 1/(1-p_ii) degenerate
        for (int ist = 0; ist < sn.nstations; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            for (int r = 0; r < sn.nclasses; r++) {
                if (sn.rtnodes.get(ind * sn.nclasses + r, ind * sn.nclasses + r) >= 1 - 1e-9) {
                    return "MEM does not support absorbing self-loop routing (reducible network).";
                }
            }
        }

        // Only non-priority disciplines are supported; the PR/HOL constraint
        // formulae are not given in Kouvatsos (1994)
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            if (sched != SchedStrategy.EXT && sched != SchedStrategy.INF
                    && sched != SchedStrategy.FCFS && sched != SchedStrategy.PS
                    && sched != SchedStrategy.SIRO && sched != SchedStrategy.LCFS
                    && sched != SchedStrategy.LCFSPR) {
                return "MEM does not support the " + sched + " scheduling strategy.";
            }
        }

        if (hasOpen) {
            // A Source node must be present for the external arrival extraction
            boolean hasSource = false;
            for (int ind = 0; ind < sn.nnodes; ind++) {
                if (sn.nodetype.get(ind) == NodeType.Source) {
                    hasSource = true;
                    break;
                }
            }
            if (!hasSource) {
                return "MEM requires a Source node when open classes are present.";
            }
        }
        if (hasClosed) {
            // Closed classes build on G/G/1 and G/G/inf queues only (Section 3.3)
            for (int ist = 0; ist < sn.nstations; ist++) {
                double ns = sn.nservers.get(ist);
                if (!Double.isInfinite(ns) && ns > 1) {
                    return "MEM does not support multiserver stations in closed or mixed models.";
                }
            }
        }

        // Finite station buffers. The censored GE/GE/c/0;N building block of
        // Section 4.1 and the transfer-blocking expansion built on it are
        // single class, so a finite buffer is admissible only in a
        // single-class open model.
        boolean anyCapped = false;
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            if (sched == SchedStrategy.EXT) {
                continue;
            }
            if (!Double.isInfinite(memBufferSize(sn, ist))) {
                anyCapped = true;
                break;
            }
        }
        if (anyCapped) {
            if (!hasOpen || hasClosed) {
                return "MEM supports finite station buffers only in open models.";
            }
            if (sn.nclasses > 1) {
                return "MEM supports finite station buffers only in single-class models: the censored GE/GE/c/0;N building block is single class.";
            }
            for (int ist = 0; ist < sn.nstations; ist++) {
                SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                if (sched == SchedStrategy.EXT) {
                    if (sn.scv != null && ist < sn.scv.getNumRows()) {
                        double v = sn.scv.get(ist, 0);
                        if (!Double.isNaN(v) && !Double.isInfinite(v) && v < 1 - 1e-12) {
                            return "MEM with finite buffers needs an external interarrival scv of at least 1: the GE distribution is not defined below 1.";
                        }
                    }
                    continue;
                }
                if (Double.isInfinite(memBufferSize(sn, ist))) {
                    continue;
                }
                double ns = sn.nservers.get(ist);
                if (Double.isInfinite(ns) || ns < 1) {
                    return "MEM cannot apply a finite buffer to an infinite-server station.";
                }
                if (sched != SchedStrategy.FCFS) {
                    return "MEM supports finite station buffers only under FCFS scheduling; a station uses " + sched + ".";
                }
                DropStrategy dr = memDropRule(sn, ist);
                if (dr != DropStrategy.Drop && dr != DropStrategy.BlockingAfterService) {
                    return "MEM supports the DROP and BAS drop rules at a finite buffer; a station uses " + dr + ".";
                }
                if (sn.scv != null && ist < sn.scv.getNumRows()) {
                    double v = sn.scv.get(ist, 0);
                    if (!Double.isNaN(v) && !Double.isInfinite(v) && v < 1 - 1e-12) {
                        return "MEM with finite buffers needs a service scv of at least 1: the GE distribution is not defined below 1.";
                    }
                }
            }
        }
        return null;
    }

    /**
     * Physical buffer size of a station in jobs, in service included: the
     * tighter of the station capacity and the sum of the per-class
     * capacities, infinite when the station is unbounded.
     * <p>
     * Only a buffer that can actually BIND is reported. refreshCapacity
     * derives a FINITE classcap (the chain population) for EVERY closed
     * model, so a plain finiteness test would report a buffer at every
     * station of every closed model; a capacity at least as large as the
     * total population can never refuse a job and is returned as infinite.
     * The population sum is infinite as soon as one class is open, so any
     * finite capacity reachable by an open class binds.
     *
     * @param sn  network structure
     * @param ist station index
     * @return buffer size in jobs
     */
    public static double memBufferSize(NetworkStruct sn, int ist) {
        double N = Double.POSITIVE_INFINITY;
        if (sn.cap != null && ist < sn.cap.length() && sn.cap.get(ist) >= 0
                && sn.cap.get(ist) < Integer.MAX_VALUE) {
            // sn.cap and sn.classcap encode "unbounded" as Integer.MAX_VALUE
            // here, not as Inf as MATLAB does, so the sentinel has to be
            // decoded rather than compared as a number: see Station.hasFiniteCap.
            N = Math.min(N, sn.cap.get(ist));
        }
        if (sn.classcap != null && ist < sn.classcap.getNumRows()) {
            double tot = 0.0;
            boolean any = false;
            boolean unbounded = false;
            for (int r = 0; r < sn.classcap.getNumCols(); r++) {
                double v = sn.classcap.get(ist, r);
                if (v > 0) {
                    if (Double.isInfinite(v) || v >= Integer.MAX_VALUE) {
                        unbounded = true; // one unbounded class unbounds the station
                        break;
                    }
                    tot += v;
                    any = true;
                }
            }
            if (any && !unbounded) {
                N = Math.min(N, tot);
            }
        }
        if (sn.njobs != null) {
            double totalJobs = 0.0;
            for (int r = 0; r < sn.nclasses; r++) {
                totalJobs += sn.njobs.get(0, r);
            }
            if (N >= totalJobs) {
                N = Double.POSITIVE_INFINITY; // declared but unreachable
            }
        }
        return N;
    }

    /**
     * Drop rule declared at a station for the first class, defaulting to
     * Drop, which is what refreshCapacity assigns to a finite buffer
     * reachable by an open class.
     *
     * @param sn  network structure
     * @param ist station index
     * @return the drop strategy in force at the station
     */
    public static DropStrategy memDropRule(NetworkStruct sn, int ist) {
        if (sn.droprule != null && ist < sn.stations.size()) {
            java.util.Map<jline.lang.JobClass, DropStrategy> perClass = sn.droprule.get(sn.stations.get(ist));
            if (perClass != null && !sn.jobclasses.isEmpty()) {
                DropStrategy dr = perClass.get(sn.jobclasses.get(0));
                if (dr != null) {
                    return dr;
                }
            }
        }
        return DropStrategy.Drop;
    }

    /**
     * True when the model is a single-class open network carrying a finite
     * buffer that MEM represents explicitly, as a censored GE/GE/c/0;N queue.
     * Callers use it to bypass the product-form capacity gate, which exists
     * because the other NC methods have no representation of a buffer.
     *
     * @param sn network structure
     * @return true when the finite-buffer MEM path applies
     */
    public static boolean memHasFiniteBuffers(NetworkStruct sn) {
        if (sn.nclasses != 1) {
            return false;
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
            if (sched == SchedStrategy.EXT) {
                continue;
            }
            if (!Double.isInfinite(memBufferSize(sn, ist))) {
                return true;
            }
        }
        return false;
    }

    /**
     * Maximum Entropy algorithm for Closed Queueing Networks
     * (Kouvatsos 1994, Section 3.3): pseudo-open decomposition followed by a
     * convolution over the population lattice, iterated on the flow (work
     * rate) equations.
     *
     * @return NCResult containing the closed ME algorithm results
     */
    public NCResult meClosed() {
        return meClosed(new jline.api.nc.MeOqnOptions());
    }

    /**
     * Maximum Entropy algorithm for Closed Queueing Networks with custom
     * options (Kouvatsos 1994, Section 3.3).
     *
     * @param meOptions options for the ME algorithm
     * @return NCResult containing the closed ME algorithm results
     */
    public NCResult meClosed(jline.api.nc.MeOqnOptions meOptions) {
        NetworkStruct sn = getStruct();

        if (!model.hasClosedClasses() || model.hasOpenClasses()) {
            line_error(mfilename(new Object() {}), "meClosed only supports closed queueing networks.");
            return null;
        }

        int M = sn.nstations;
        int R = sn.nclasses;

        // Validate the model against the closed MEM feature set
        String memReason = memUnsupportedReason(sn, true);
        if (memReason != null) {
            line_error(mfilename(new Object() {}), memReason);
            return null;
        }

        // Populations and reference stations
        Matrix N = new Matrix(1, R);
        int[] refstat = new int[R];
        for (int r = 0; r < R; r++) {
            N.set(0, r, sn.njobs.get(r));
            refstat[r] = (int) sn.refstat.get(r);
        }

        // Service rates, scvs and server counts
        Matrix mu = new Matrix(M, R);
        Matrix Cs = Matrix.ones(M, R);
        Matrix nservers = Matrix.ones(M, 1);
        for (int ist = 0; ist < M; ist++) {
            nservers.set(ist, 0, sn.nservers.get(ist));
            for (int r = 0; r < R; r++) {
                double rate = sn.rates.get(ist, r);
                if (!Double.isNaN(rate) && !Double.isInfinite(rate) && rate > 0) {
                    mu.set(ist, r, rate);
                    if (sn.scv != null && ist < sn.scv.getNumRows() && r < sn.scv.getNumCols()) {
                        double v = sn.scv.get(ist, r);
                        if (!Double.isNaN(v) && !Double.isInfinite(v) && v > 0) {
                            Cs.set(ist, r, v);
                        }
                    }
                }
            }
        }

        // Routing probabilities between stations
        Matrix[][] P = new Matrix[M][M];
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < M; k++) {
                P[j][k] = new Matrix(R, 1);
            }
        }
        for (int r = 0; r < R; r++) {
            for (int j = 0; j < M; j++) {
                int jNode = (int) sn.stationToNode.get(j);
                for (int k = 0; k < M; k++) {
                    int iNode = (int) sn.stationToNode.get(k);
                    P[j][k].set(r, 0, sn.rtnodes.get(jNode * R + r, iNode * R + r));
                }
            }
        }

        // Run the two-stage Maximum Entropy algorithm
        boolean[] insens = new boolean[M];
        for (int ist = 0; ist < M; ist++) {
            SchedStrategy sk = sn.sched.get(sn.stations.get(ist));
            insens[ist] = (sk == SchedStrategy.PS || sk == SchedStrategy.LCFSPR);
        }
        jline.api.nc.MeCqnResult meResult = jline.api.nc.Me_cqn.me_cqn(M, R, N, mu, Cs, P, nservers, refstat, insens, meOptions);

        NCResult result = new NCResult();
        result.QN = meResult.getL().copy();
        result.UN = meResult.getRho().copy();
        result.RN = meResult.getW().copy();
        result.TN = meResult.getLambda().copy();
        result.CN = new Matrix(1, R);
        result.XN = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double xr = meResult.getX().get(0, r);
            result.XN.set(0, r, xr);
            if (xr > 0) {
                result.CN.set(0, r, N.get(0, r) / xr); // class cycle time
            }
        }
        result.iter = meResult.getIter();
        result.method = "mem";

        return result;
    }

    /**
     * Maximum Entropy algorithm for Mixed Queueing Networks: composition of
     * the open (Section 3.2) and closed (Section 3.3) algorithms with
     * product-form-style conditioning.
     *
     * @return NCResult containing the mixed ME algorithm results
     */
    public NCResult meMixed() {
        return meMixed(new jline.api.nc.MeOqnOptions());
    }

    /**
     * Maximum Entropy algorithm for Mixed Queueing Networks with custom
     * options.
     *
     * @param meOptions options for the ME algorithm
     * @return NCResult containing the mixed ME algorithm results
     */
    public NCResult meMixed(jline.api.nc.MeOqnOptions meOptions) {
        NetworkStruct sn = getStruct();

        if (!model.hasClosedClasses() || !model.hasOpenClasses()) {
            line_error(mfilename(new Object() {}), "meMixed only supports mixed open/closed queueing networks.");
            return null;
        }

        int M = sn.nstations;
        int R = sn.nclasses;

        // Validate the model against the mixed MEM feature set
        String memReason = memUnsupportedReason(sn, true, true);
        if (memReason != null) {
            line_error(mfilename(new Object() {}), memReason);
            return null;
        }

        boolean[] openCls = new boolean[R];
        Matrix N = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            N.set(0, r, sn.njobs.get(r));
            openCls[r] = Double.isInfinite(sn.njobs.get(r));
        }

        // Locate the source station (external arrivals of the open classes)
        int sourceIdx = -1;
        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                sourceIdx = (int) sn.nodeToStation.get(i);
                break;
            }
        }

        // Queueing/delay stations (all stations except the source)
        int Mq = M - 1;
        int[] qs = new int[Mq];
        int pos = 0;
        for (int ist = 0; ist < M; ist++) {
            if (ist != sourceIdx) {
                qs[pos] = ist;
                pos++;
            }
        }

        Matrix mu = new Matrix(Mq, R);
        Matrix Cs = Matrix.ones(Mq, R);
        Matrix nservers = Matrix.ones(Mq, 1);
        int[] refstat = new int[R];
        for (int k = 0; k < Mq; k++) {
            int ist = qs[k];
            nservers.set(k, 0, sn.nservers.get(ist));
            for (int r = 0; r < R; r++) {
                double rate = sn.rates.get(ist, r);
                if (!Double.isNaN(rate) && !Double.isInfinite(rate) && rate > 0) {
                    mu.set(k, r, rate);
                    if (sn.scv != null && ist < sn.scv.getNumRows() && r < sn.scv.getNumCols()) {
                        double v = sn.scv.get(ist, r);
                        if (!Double.isNaN(v) && !Double.isInfinite(v) && v > 0) {
                            Cs.set(k, r, v);
                        }
                    }
                }
            }
        }
        for (int r = 0; r < R; r++) {
            if (!openCls[r]) {
                int ref = (int) sn.refstat.get(r);
                for (int k = 0; k < Mq; k++) {
                    if (qs[k] == ref) {
                        refstat[r] = k;
                        break;
                    }
                }
            }
        }

        // External arrivals of the open classes along the source routing
        Matrix lambda0 = new Matrix(Mq, R);
        Matrix Ca0 = new Matrix(Mq, R);
        int sourceNode = (int) sn.stationToNode.get(sourceIdx);
        for (int r = 0; r < R; r++) {
            double extRate = sn.rates.get(sourceIdx, r);
            if (openCls[r] && !Double.isNaN(extRate) && !Double.isInfinite(extRate) && extRate > 0) {
                double Ca_ext = 1.0;
                if (sn.scv != null && sourceIdx < sn.scv.getNumRows() && r < sn.scv.getNumCols()) {
                    double v = sn.scv.get(sourceIdx, r);
                    if (!Double.isNaN(v) && !Double.isInfinite(v) && v > 0) {
                        Ca_ext = v;
                    }
                }
                for (int k = 0; k < Mq; k++) {
                    int destNode = (int) sn.stationToNode.get(qs[k]);
                    double routeProb = sn.rtnodes.get(sourceNode * R + r, destNode * R + r);
                    if (routeProb > 0) {
                        lambda0.set(k, r, extRate * routeProb);
                        Ca0.set(k, r, Ca_ext);
                    }
                }
            }
        }

        // Routing probabilities between queueing stations
        Matrix[][] P = new Matrix[Mq][Mq];
        for (int j = 0; j < Mq; j++) {
            for (int k = 0; k < Mq; k++) {
                P[j][k] = new Matrix(R, 1);
            }
        }
        for (int r = 0; r < R; r++) {
            for (int j = 0; j < Mq; j++) {
                int jNode = (int) sn.stationToNode.get(qs[j]);
                for (int k = 0; k < Mq; k++) {
                    int iNode = (int) sn.stationToNode.get(qs[k]);
                    P[j][k].set(r, 0, sn.rtnodes.get(jNode * R + r, iNode * R + r));
                }
            }
        }

        // Run the composed Maximum Entropy algorithm
        boolean[] insens = new boolean[Mq];
        for (int k = 0; k < Mq; k++) {
            SchedStrategy sk = sn.sched.get(sn.stations.get(qs[k]));
            insens[k] = (sk == SchedStrategy.PS || sk == SchedStrategy.LCFSPR);
        }
        jline.api.nc.MeCqnResult meResult = jline.api.nc.Me_mqn.me_mqn(Mq, R, openCls, lambda0, Ca0, N, mu, Cs, P, nservers, refstat, insens, meOptions);

        NCResult result = new NCResult();
        result.QN = new Matrix(M, R);
        result.UN = new Matrix(M, R);
        result.RN = new Matrix(M, R);
        result.TN = new Matrix(M, R);
        for (int k = 0; k < Mq; k++) {
            int ist = qs[k];
            for (int r = 0; r < R; r++) {
                result.QN.set(ist, r, meResult.getL().get(k, r));
                result.UN.set(ist, r, meResult.getRho().get(k, r));
                result.RN.set(ist, r, meResult.getW().get(k, r));
                result.TN.set(ist, r, meResult.getLambda().get(k, r));
            }
        }

        // Cap utilization of unstable stations at 1 (LINE convention)
        for (int k = 0; k < Mq; k++) {
            int ist = qs[k];
            if (!Double.isInfinite(sn.nservers.get(ist))) {
                double utot = 0.0;
                for (int r = 0; r < R; r++) {
                    utot += result.UN.get(ist, r);
                }
                if (utot > 1) {
                    for (int r = 0; r < R; r++) {
                        result.UN.set(ist, r, result.UN.get(ist, r) / utot);
                    }
                }
            }
        }

        result.CN = new Matrix(1, R);
        result.XN = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double xr = meResult.getX().get(0, r);
            result.XN.set(0, r, xr);
            if (openCls[r]) {
                result.TN.set(sourceIdx, r, xr);
                if (xr > 0) {
                    double qtot = 0.0;
                    for (int k = 0; k < Mq; k++) {
                        qtot += result.QN.get(qs[k], r);
                    }
                    result.CN.set(0, r, qtot / xr);
                }
            } else if (xr > 0) {
                result.CN.set(0, r, N.get(0, r) / xr); // class cycle time
            }
        }
        result.iter = meResult.getIter();
        result.method = "mem";

        return result;
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

        // Validate the model against the MEM feature set
        String memReason = memUnsupportedReason(sn);
        if (memReason != null) {
            line_error(mfilename(new Object() {}), memReason);
            return null;
        }

        // Locate the source station (external arrivals)
        int sourceIdx = -1;
        for (int i = 0; i < sn.nnodes; i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                sourceIdx = (int) sn.nodeToStation.get(i);
                break;
            }
        }

        // Queueing/delay stations (all stations except the source)
        int Mq = M - 1;
        int[] qs = new int[Mq];
        int pos = 0;
        for (int ist = 0; ist < M; ist++) {
            if (ist != sourceIdx) {
                qs[pos] = ist;
                pos++;
            }
        }

        // Service rates, scvs and server counts (rates/scv are NaN for
        // classes not served at a station)
        Matrix mu = new Matrix(Mq, R);
        Matrix Cs = Matrix.ones(Mq, R);
        Matrix nservers = Matrix.ones(Mq, 1);
        for (int k = 0; k < Mq; k++) {
            int ist = qs[k];
            nservers.set(k, 0, sn.nservers.get(ist));
            for (int r = 0; r < R; r++) {
                double rate = sn.rates.get(ist, r);
                if (!Double.isNaN(rate) && !Double.isInfinite(rate) && rate > 0) {
                    mu.set(k, r, rate);
                    if (sn.scv != null && ist < sn.scv.getNumRows() && r < sn.scv.getNumCols()) {
                        double v = sn.scv.get(ist, r);
                        if (!Double.isNaN(v) && !Double.isInfinite(v) && v > 0) {
                            Cs.set(k, r, v);
                        }
                    }
                }
            }
        }

        // External arrivals: distribute the source output along its routing
        Matrix lambda0 = new Matrix(Mq, R);
        Matrix Ca0 = new Matrix(Mq, R);
        int sourceNode = (int) sn.stationToNode.get(sourceIdx);
        for (int r = 0; r < R; r++) {
            double extRate = sn.rates.get(sourceIdx, r);
            if (!Double.isNaN(extRate) && !Double.isInfinite(extRate) && extRate > 0) {
                double Ca_ext = 1.0;
                if (sn.scv != null && sourceIdx < sn.scv.getNumRows() && r < sn.scv.getNumCols()) {
                    double v = sn.scv.get(sourceIdx, r);
                    if (!Double.isNaN(v) && !Double.isInfinite(v) && v > 0) {
                        Ca_ext = v;
                    }
                }
                for (int k = 0; k < Mq; k++) {
                    int destNode = (int) sn.stationToNode.get(qs[k]);
                    double routeProb = sn.rtnodes.get(sourceNode * R + r, destNode * R + r);
                    if (routeProb > 0) {
                        lambda0.set(k, r, extRate * routeProb);
                        Ca0.set(k, r, Ca_ext);
                    }
                }
            }
        }

        // Routing probabilities between queueing stations
        Matrix[][] P = new Matrix[Mq][Mq];
        for (int j = 0; j < Mq; j++) {
            for (int k = 0; k < Mq; k++) {
                P[j][k] = new Matrix(R, 1);
            }
        }
        for (int r = 0; r < R; r++) {
            for (int j = 0; j < Mq; j++) {
                int jNode = (int) sn.stationToNode.get(qs[j]);
                for (int k = 0; k < Mq; k++) {
                    int iNode = (int) sn.stationToNode.get(qs[k]);
                    P[j][k].set(r, 0, sn.rtnodes.get(jNode * R + r, iNode * R + r));
                }
            }
        }

        // Finite buffers: censored GE/GE/c/0;N building blocks, with the
        // holding-node expansion where the drop rule is BAS (Section 4.1 of
        // the source plus Tahilramani, Manjunath and Bose 1999)
        if (memHasFiniteBuffers(sn)) {
            Matrix Nbuf = new Matrix(Mq, 1);
            Matrix Pflat = new Matrix(Mq, Mq);
            int[] blockrule = new int[Mq];
            for (int k = 0; k < Mq; k++) {
                Nbuf.set(k, 0, memBufferSize(sn, qs[k]));
                blockrule[k] = memDropRule(sn, qs[k]) == DropStrategy.BlockingAfterService
                        ? jline.api.nc.Me_oqn_blk.RULE_BAS : jline.api.nc.Me_oqn_blk.RULE_LOSS;
                for (int l = 0; l < Mq; l++) {
                    Pflat.set(k, l, P[k][l].get(0, 0));
                }
            }
            Matrix lam0 = new Matrix(Mq, 1);
            Matrix ca0 = new Matrix(Mq, 1);
            Matrix mu1 = new Matrix(Mq, 1);
            Matrix cs1 = new Matrix(Mq, 1);
            for (int k = 0; k < Mq; k++) {
                lam0.set(k, 0, lambda0.get(k, 0));
                ca0.set(k, 0, Ca0.get(k, 0) > 0 ? Ca0.get(k, 0) : 1.0);
                mu1.set(k, 0, mu.get(k, 0));
                cs1.set(k, 0, Cs.get(k, 0));
            }
            jline.api.nc.MeOqnBlkResult blkResult = jline.api.nc.Me_oqn_blk.me_oqn_blk(
                    Mq, lam0, ca0, mu1, cs1, Pflat, nservers, Nbuf, blockrule, meOptions, 0.5);
            NCResult blk = new NCResult();
            blk.QN = new Matrix(M, R);
            blk.UN = new Matrix(M, R);
            blk.RN = new Matrix(M, R);
            blk.TN = new Matrix(M, R);
            blk.CN = new Matrix(1, R);
            blk.XN = new Matrix(1, R);
            for (int k = 0; k < Mq; k++) {
                blk.QN.set(qs[k], 0, blkResult.getQ().get(k, 0));
                blk.UN.set(qs[k], 0, blkResult.getU().get(k, 0));
                blk.RN.set(qs[k], 0, blkResult.getW().get(k, 0));
                blk.TN.set(qs[k], 0, blkResult.getT().get(k, 0));
            }
            // The source emits at its nominal rate; the jobs lost at a full
            // buffer never reach a station, so the carried flow reported per
            // station is below it by the loss probability at the entry
            // stations.
            double extRate0 = sn.rates.get(sourceIdx, 0);
            if (!Double.isNaN(extRate0) && !Double.isInfinite(extRate0) && extRate0 > 0) {
                blk.TN.set(sourceIdx, 0, extRate0);
                blk.XN.set(0, 0, extRate0);
                double accepted = extRate0;
                for (int k = 0; k < Mq; k++) {
                    if (lam0.get(k, 0) > 0) {
                        accepted -= lam0.get(k, 0) * blkResult.getPBa().get(k, 0);
                    }
                }
                double qtot = 0.0;
                for (int k = 0; k < Mq; k++) {
                    qtot += blk.QN.get(qs[k], 0);
                }
                if (accepted > 0) {
                    blk.CN.set(0, 0, qtot / accepted);
                }
            }
            blk.iter = blkResult.getIter();
            // Reported under its own name so that solver.citations() reaches
            // the transfer-blocking reference on top of the base MEM one.
            blk.method = "mem.blocking";
            return blk;
        }

        // Run the Maximum Entropy fixed-point algorithm
        boolean[] insens = new boolean[Mq];
        for (int k = 0; k < Mq; k++) {
            SchedStrategy sk = sn.sched.get(sn.stations.get(qs[k]));
            insens[k] = (sk == SchedStrategy.PS || sk == SchedStrategy.LCFSPR);
        }
        jline.api.nc.MeOqnResult meResult = jline.api.nc.Me_oqn.me_oqn(Mq, R, lambda0, Ca0, mu, Cs, P, nservers, insens, meOptions);

        // Map results back to station-indexed LINE outputs
        NCResult result = new NCResult();
        result.QN = new Matrix(M, R);
        result.UN = new Matrix(M, R);
        result.RN = new Matrix(M, R);
        result.TN = new Matrix(M, R);
        for (int k = 0; k < Mq; k++) {
            int ist = qs[k];
            for (int r = 0; r < R; r++) {
                result.QN.set(ist, r, meResult.getL().get(k, r));
                result.UN.set(ist, r, meResult.getRho().get(k, r));
                result.RN.set(ist, r, meResult.getW().get(k, r));
                result.TN.set(ist, r, meResult.getLambda().get(k, r));
            }
        }

        // Cap utilization of unstable stations at 1 (LINE convention)
        for (int k = 0; k < Mq; k++) {
            int ist = qs[k];
            if (!Double.isInfinite(sn.nservers.get(ist))) {
                double utot = 0.0;
                for (int r = 0; r < R; r++) {
                    utot += result.UN.get(ist, r);
                }
                if (utot > 1) {
                    for (int r = 0; r < R; r++) {
                        result.UN.set(ist, r, result.UN.get(ist, r) / utot);
                    }
                }
            }
        }

        // Source station: report the external arrival rates as throughputs
        // and derive system metrics by Little's law
        result.CN = new Matrix(1, R);
        result.XN = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double extRate = sn.rates.get(sourceIdx, r);
            if (!Double.isNaN(extRate) && !Double.isInfinite(extRate) && extRate > 0) {
                result.TN.set(sourceIdx, r, extRate);
                result.XN.set(0, r, extRate);
                double qtot = 0.0;
                for (int k = 0; k < Mq; k++) {
                    qtot += result.QN.get(qs[k], r);
                }
                result.CN.set(0, r, qtot / extRate);
            }
        }

        result.iter = meResult.getIter();
        result.method = "mem";

        return result;
    }

    /**
     * Sample node state trajectory
     *
     * @param node The node to sample
     * @param numEvents Number of samples to generate
     * @return Sample result containing state trajectory
     */
    public Object sample(Node node, int numEvents) {
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
        if (!supportsModel(model)) {
            return false;
        }
        // Registry inclusion cannot see finite capacity (there is no feature
        // name for it), so apply the structural gate as well, unless the
        // finite-buffer MEM path covers the model.
        NetworkStruct s = model.getStruct(false);
        if (!memFiniteBufferPath(s, this.options)) {
            String reason = NetworkSolver.bindingCapacityReason(model, s, "SolverNC");
            if (reason != null) {
                line_warning(mfilename(new Object() {}), reason);
                return false;
            }
        }
        return true;
    }

    /**
     * True when method='mem' has been requested on a model whose finite
     * buffers MEM carries explicitly. Callers use it to skip the
     * product-form capacity gate, which exists only because the other NC
     * methods cannot represent a buffer.
     *
     * @param sn      network structure
     * @param options solver options, whose method selects MEM
     * @return true when the finite-buffer MEM path applies
     */
    public static boolean memFiniteBufferPath(NetworkStruct sn, SolverOptions options) {
        if (sn == null || options == null || !"mem".equals(options.method)) {
            return false;
        }
        if (!memHasFiniteBuffers(sn)) {
            return false;
        }
        return memUnsupportedReason(sn, true, false) == null;
    }

    /**
     * Feature-driven resolution of method='default': an open network with
     * non-Markovian (non-unit SCV) variability within the MEM feature set is
     * solved by the Maximum Entropy Method by default, since the
     * normalizing-constant path would silently exponentialize it. Mirrors the
     * dispatch in runAnalyzer and the MATLAB/Python SolverNC.resolveMethod.
     */
    @Override
    public String resolveMethod(SolverOptions options) {
        String method = options == null ? "default" : options.method;
        if ("default".equals(method) && this.model != null
                && model.hasOpenClasses() && !model.hasClosedClasses()) {
            NetworkStruct s = this.sn != null ? this.sn : this.model.getStruct(false);
            if (memUnsupportedReason(s) == null && snHasNonUnitScv(s)) {
                return "mem";
            }
        }
        return method;
    }

    /**
     * Method-aware gate. MEM (Kouvatsos maximum entropy) has structural
     * applicability rules beyond a flat feature set (open-only, Source/Queue/
     * Delay/Sink, non-priority scheduling); delegate to memUnsupportedReason,
     * which returns a precise reason (null when supported). All other NC methods
     * use the coarse product-form feature gate.
     */
    @Override
    public String supportsModelMethod(String method) {
        if ("mem".equals(method)) {
            NetworkStruct s = this.sn != null ? this.sn : this.model.getStruct(false);
            // MEM now supports both open (Section 3.2) and closed (Section 3.3)
            // networks, so validate against the model's actual class composition.
            String r = memUnsupportedReason(s, model.hasOpenClasses(), model.hasClosedClasses());
            return r == null ? "" : r;
        }
        return supports((Network) this.model) ? "" : "Some features are not supported by the chosen solver.";
    }

    /**
     * NC is deterministic except for the Monte Carlo integration methods
     * (mci/imci), logistic sampling (ls), and the sampling method, whose
     * estimates depend on the random seed. Method names are tokenized so that
     * runtime-resolved names such as "default/imci" and prefixed names such as
     * "nc.ls" classify correctly.
     *
     * @param method the method name to classify
     * @return true if the method returns stochastic estimates
     */
    @Override
    public boolean isStochasticMethod(String method) {
        if (method == null || method.isEmpty()) {
            return false;
        }
        String[] tokens = method.toLowerCase().split("[./]");
        for (String tok : tokens) {
            if (tok.equals("mci") || tok.equals("imci") || tok.equals("ls") || tok.equals("sampling") || tok.equals("is")) {
                return true;
            }
        }
        return false;
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

    /**
     * One inner solve of the normalizing-constant analyzer, as the callback of
     * the fork-join fixed point. Mirrors matlab/src/solvers/NC/@SolverNC/ncDispatch.m.
     *
     * <p>No method override is applied: the transformed model carries auxiliary
     * open classes at a vanishing rate, which the default route already resolves.
     *
     * @param snIn the struct of the transformed model
     * @param opts the solver options
     * @return the metrics of that solve, in the neutral carrier the driver consumes
     */
    private MVAResult ncDispatch(NetworkStruct snIn, SolverOptions opts) {
        this.sn = snIn;
        NCResult r;
        boolean ld = (snIn.lldscaling != null && !snIn.lldscaling.isEmpty())
                || (snIn.cdscaling != null && !snIn.cdscaling.isEmpty())
                || (snIn.jdscaling != null && !snIn.jdscaling.isEmpty());
        if (ld) {
            r = solver_ncld_analyzer(snIn, opts.copy());
        } else {
            r = solver_nc_analyzer(snIn, opts.copy());
        }
        MVAResult out = new MVAResult();
        out.QN = r.QN;
        out.UN = r.UN;
        out.RN = r.RN;
        out.TN = r.TN;
        out.CN = r.CN;
        out.XN = r.XN;
        out.logNormConstAggr = r.lG;
        out.runtime = r.runtime;
        out.iter = r.it;
        out.method = opts.method;
        return out;
    }

}
