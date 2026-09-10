/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import java.util.HashMap;
import java.util.Map;
import jline.lang.FeatureSet;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.GlobalConstants;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.NodeParam;
import jline.lang.constant.SchedStrategy;
import jline.api.sn.SnHasProductForm;
import jline.api.sn.SnGetDemandsChain;
import jline.api.sn.SnRtStations;
import jline.api.pfqn.Pfqn_busyp;
import jline.io.Ret;
import jline.util.Pair;
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
import static jline.api.sn.SnHasDPS.snHasDPS;
import static jline.api.sn.SnHasMultiServer.snHasMultiServer;
import static jline.io.InputOutput.*;
import static jline.solvers.nc.analyzers.Solver_nc_analyzer.solver_nc_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_cache_analyzer.solver_nc_cache_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_cache_qn_analyzer.solver_nc_cache_qn_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_retrieval_analyzer.solver_nc_retrieval_analyzer;
import static jline.solvers.nc.analyzers.Solver_ncld_analyzer.solver_ncld_analyzer;
import static jline.solvers.nc.analyzers.Solver_nc_lossn_analyzer.solver_nc_lossn_analyzer;
import jline.solvers.nc.analyzers.Solver_nc_dt_analyzer;
import jline.solvers.nc.analyzers.Solver_nc_dps_analyzer;
import jline.solvers.nc.analyzers.Solver_nc_sdr_analyzer;
import static jline.api.sn.SnHasClosedClasses.snHasClosedClasses;
import static jline.api.sn.SnHasOpenClasses.snHasOpenClasses;
import jline.lang.constant.DropStrategy;
import jline.solvers.nc.handlers.Solver_nc_pas_is;
import static jline.solvers.nc.handlers.Solver_nc_marg.solver_nc_marg;
import static jline.solvers.nc.handlers.Solver_nc_joint.solver_nc_joint;
import static jline.solvers.nc.handlers.Solver_nc_margaggr.solver_nc_margaggr;
import static jline.solvers.nc.handlers.Solver_nc_jointaggr.solver_nc_jointaggr;
import jline.solvers.nc.handlers.Solver_nc_jointmarg;
import static jline.solvers.nc.handlers.Solver_nc_jointmarg.solver_nc_jointmarg;
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
                // Geometric is admitted for the discrete-time route only
                // (options.config.slotted, Solver_nc_dt_analyzer); on the
                // continuous-time routes it is treated by its mean and SCV
                "Geometric",
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
                // DPS is served ONLY in Morrison's closed think+DPS shape
                // (Solver_nc_dps_analyzer.nc_is_dps_model). A boolean feature cannot
                // express that restriction, so runAnalyzer keeps an imperative check
                // for every other DPS model, the same pattern as "Region".
                "SchedStrategy_DPS",
                "SchedStrategy_LCFS", "SchedStrategy_LCFSPR", "SchedStrategy_OI",
                "SchedStrategy_PAS",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND", "RoutingStrategy_SDR",
                "SchedStrategy_FCFS", "ClosedClass", "SelfLoopingClass",
                "Cache", "CacheClassSwitcher", "CacheRetrieval", "CacheItemSize", "OpenClass",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO",
                "ReplacementStrategy_HLRU",
                "LoadDependence", "ClassDependence", "JointDependence",
                // Fork-join through the MMT/HT transformation, driven by
                // jline.solvers.fj.FJFixedPoint (as in SolverMVA)
                "Fork", "Forker", "Join", "Joiner",
                // quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
                "JoinPartial",
                // Petri nets: the "rec" route (Solver_nc_spn_analyzer) walks the
                // reachable set in a decision diagram, so a Place is a token
                // container rather than a station with a service process. A
                // queueing Place is refused by Spn_pf, which is where the
                // product-form class is decided.
                "Place", "Transition", "Linkage", "Enabling", "Inhibiting", "Timing",
                "Firing", "Storage",
                // c-server stations: every route but "divdiff" carries the count
                // (see methodFeatureSet). FiniteCapacity is deliberately NOT here:
                // mem, default and exact are granted it per method, the rest solve
                // a buffer away.
                "MultiServer"
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

    /**
     * Joint probability that station i holds nvec(i) jobs IN TOTAL, all classes
     * summed out.
     *
     * <p>Compare with {@link #getProbSysAggr()}, which fixes the PER-CLASS
     * population of every station and is a product form; each value returned
     * here is the sum of that one over every per-class table with these row
     * sums. Compare also with getProbMarg, which is the one-station marginal of
     * this law. The quantity is a permanent of the demand matrix replicated once
     * per job (Ryser 1963), so it needs no enumeration of that fibre.</p>
     *
     * @param nvec per-station total job counts, summing to the closed population
     * @return the joint probability and its logarithm
     */
    public ProbabilityResult getProbSysMarg(Matrix nvec) {
        return getProbSysMarg(nvec, "exact");
    }

    /**
     * Joint probability of the per-station total queue lengths, evaluated with
     * a chosen permanent engine.
     *
     * @param nvec   per-station total job counts
     * @param engine "exact" (default), "spm", "bethe", "heur", "huberlaw" or "adapart".
     *               Only "exact" is exact; the others are refused on a demand
     *               matrix with a structural zero rather than having it floored,
     *               since they need full support.
     * @return the joint probability and its logarithm
     */
    public ProbabilityResult getProbSysMarg(Matrix nvec, String engine) {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult(Double.NaN);
        }

        long startTimeMillis = System.nanoTime();
        this.model.refreshStruct(true);
        NetworkStruct sn = this.model.getStruct(true);
        resetRandomGeneratorSeed(options.seed);

        NCResult ncResult = (NCResult) this.result;
        // Reuse the constant when a previous getter already paid for it:
        // sweeping the whole lattice of total states otherwise recomputes G once
        // per state.
        Double lGin = null;
        if (ncResult != null && ncResult.prob != null && ncResult.prob.logNormConstAggr != null
                && !Utils.isInf(ncResult.prob.logNormConstAggr)
                && !Double.isNaN(ncResult.prob.logNormConstAggr)) {
            lGin = ncResult.prob.logNormConstAggr;
        }

        Solver_nc_jointmarg.Ret_jointmarg ret =
                solver_nc_jointmarg(sn, this.options, nvec, engine, lGin);

        this.lastPermEngine = (engine == null || engine.isEmpty()) ? "exact" : engine.toLowerCase();
        ncResult.solver = this.name;
        ncResult.prob.logNormConstAggr = ret.lG;
        ncResult.prob.joint = ret.Pr;

        long endTimeMillis = System.nanoTime();
        ncResult.runtime = (double) (endTimeMillis - startTimeMillis) / 1000000000.0;

        ProbabilityResult probResult = new ProbabilityResult(ret.Pr);
        probResult.logNormalizingConstant = ret.lG;
        probResult.isAggregated = true;
        probResult.state = nvec.copy();
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

        // MODEL TRANSFORMATION, opt-in through options.config.transform. Mirrors
        // the branch MATLAB puts in the shared runAnalyzerPreamble, so a strategy
        // written once serves NC as well as MVA and CTMC. Note this is NOT the
        // 'lc' METHOD method name, which selects the pfqn_bklc kernel and stays the
        // fast path; the two are different computations.
        if (jline.solvers.tr.TransformSolve.isRequested(options)) {
            runTransformAnalyzer();
            return;
        }

        // Finite station/class capacity: a product-form method has no
        // representation of a finite buffer, so it would silently return the
        // unconstrained answer (QLen=4 instead of the M/M/1/2 value 0.8525).
        // Same defect as BUG-39 on the MVA side.
        //
        // method='mem' is the one exception, and only on a single-class open
        // model that memUnsupportedReason has cleared: MEM does represent the
        // buffer, as a censored GE/GE/c/0;N queue, and handles both a lost
        // arrival (DROP) and a job held in the upstream server (BAS).
        // The discrete-time route is the other exception: a finite buffer on a
        // Bernoulli server is the loss system of Daduna's corollary 2.8, which
        // Solver_nc_dt_analyzer solves exactly. The third is the single-station
        // M/M/1/K with tail drop, answered exactly by the Qsys_mm1k_loss branch
        // of runAnalyzerBody; the names that branch does NOT serve are refused
        // by ncMethodRefusal instead, so nothing reaches the recursion.
        if (!memFiniteBufferPath(this.sn, this.options)
                && !Solver_nc_dt_analyzer.isSlotted(this.options)
                && !jline.api.sn.SnIsMm1kLoss.snIsMm1kLoss(this.sn)) {
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
        boolean knownMethod = false;
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
            // forReport=false: this gate sits on the RUN path, and the run asks what
            // the reference DOES rather than what the report should offer. The two
            // differ for "mmint2"/"gleint", which the reference performs by name.
            String reason = this.supportsModelMethod(options.method, false);
            if (!reason.isEmpty()) {
                line_error(mfilename(new Object() {
                }), "This model contains features not supported by the solver. " + reason);
                return;
            }
        }
        line_debug(options.verbose, String.format("NC solver starting: method=%s, nstations=%d, nclasses=%d", 
            options.method, sn.nstations, sn.nclasses));

        // THE STRUCTURAL METHOD GATE, asked once and in one place.
        //
        // ncMethodRefusal holds every rule of the form "this method has no route
        // on this model": the DPS shape, state-dependent routing, the Petri net,
        // the order-independent rank rate, the cache and loss-network tokens,
        // PANACEA's normal usage. supportsModelMethod asks the SAME method, which
        // is what keeps Network.findSolver from offering a pair that would raise
        // here. It runs unconditionally, so it also holds for the
        // enableChecks=false SolverLN layer backend.
        // forReport=false: this is the RUN, and it asks what the reference DOES rather
        // than what the report should offer. The two answers differ for
        // "mmint2"/"gleint", which Pfqn_nc answers with an empty constant and a zero
        // table; see ncMethodRefusal.
        String ncRefusal = ncMethodRefusal(this.sn, options.method, this.options, false);
        if (!ncRefusal.isEmpty()) {
            throw new RuntimeException(ncRefusal);
        }

        // Closed think+DPS network -> Solver_nc_dps_analyzer (Morrison's heavy-usage
        // generating-function expansion), the DEFAULT for that shape. Intercepted
        // FIRST, ahead of every other route: "SchedStrategy_DPS" is declared in
        // featSupported, which opens all of them to a DPS model, and each would
        // silently drop the weights and answer with the egalitarian-PS network. Not a
        // product-form route: lG is NaN. See _kb/06-solver-catalog.md (NC section).
        if (Solver_nc_dps_analyzer.nc_is_dps_model(this.sn)) {
            if ("default".equalsIgnoreCase(options.method) || "morrison".equalsIgnoreCase(options.method)) {
                line_debug(options.verbose, "NC: closed think+DPS network, routing to Solver_nc_dps_analyzer (Morrison)");
                NCResult dpsRet = Solver_nc_dps_analyzer.solver_nc_dps_analyzer(this.sn, this.options);
                AvgHandle Tdps = getAvgTputHandles();
                Matrix ANdps = snGetArvRFromTput(this.sn, dpsRet.TN, Tdps);
                this.setAvgResults(dpsRet.QN, dpsRet.UN, dpsRet.RN, dpsRet.TN, ANdps, new Matrix(0, 0),
                        dpsRet.CN, dpsRet.XN, dpsRet.runtime, dpsRet.method, dpsRet.it);
                ((NCResult) this.result).prob.logNormConstAggr = dpsRet.lG;
                return;
            }
            // The three refusal arms that used to follow -- another method on a DPS
            // model, a DPS station outside Morrison's shape, and "morrison" on a model
            // with no DPS station at all -- moved into ncMethodRefusal above, with
            // their wording unchanged.
        }

        // Krzesinski state-dependent routing: the model has its own product form
        // (eq. 16), so it is intercepted before the standard convolution and MVA
        // analyzers, which assume state-independent routing; see
        // _kb/16-state-dependent-routing.md
        if (this.sn.sdr != null) {
            line_debug(options.verbose, "NC: state-dependent routing, routing to Solver_nc_sdr_analyzer");
            NCResult sdrRet = Solver_nc_sdr_analyzer.solver_nc_sdr_analyzer(this.sn, this.options);
            AvgHandle Tsdr = getAvgTputHandles();
            Matrix ANsdr = snGetArvRFromTput(this.sn, sdrRet.TN, Tsdr);
            this.setAvgResults(sdrRet.QN, sdrRet.UN, sdrRet.RN, sdrRet.TN, ANsdr, new Matrix(0, 0),
                    sdrRet.CN, sdrRet.XN, sdrRet.runtime, sdrRet.method, sdrRet.it);
            ((NCResult) this.result).prob.logNormConstAggr = sdrRet.lG;
            return;
        }
        // The arm that used to follow -- "sdr"/"sdr.mva" on a model declaring no
        // state-dependent routing -- moved into ncMethodRefusal, wording unchanged.

        // Discrete-time (slotted) route: explicit request only, and an error
        // rather than a fallback when the model is outside the discrete-time
        // product form; see _kb/06-solver-catalog.md (NC section)
        if (Solver_nc_dt_analyzer.isSlotted(options)) {
            line_debug(options.verbose, "NC: slotted model, routing to Solver_nc_dt_analyzer");
            NCResult dtRet = Solver_nc_dt_analyzer.solver_nc_dt_analyzer(this.sn, this.options);
            AvgHandle Tdt = getAvgTputHandles();
            Matrix ANdt = snGetArvRFromTput(this.sn, dtRet.TN, Tdt);
            this.setAvgResults(dtRet.QN, dtRet.UN, dtRet.RN, dtRet.TN, ANdt, new Matrix(0, 0),
                    dtRet.CN, dtRet.XN, dtRet.runtime, dtRet.method, dtRet.it);
            ((NCResult) this.result).prob.logNormConstAggr = dtRet.lG;
            return;
        }

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
                        public MVAResult solve(jline.lang.Network net, NetworkStruct snIn, SolverOptions opts) {
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

        // A stochastic Petri net takes the MDD-rec route: the reachable set lives
        // in a decision diagram and the product form supplies the rates, so none
        // of the queueing-network branches below apply to it.
        if (sn.nodetype.contains(NodeType.Place)) {
            line_debug(options.verbose, "NC: detected stochastic Petri net, calling "
                    + "solver_nc_spn_analyzer (MDD-rec)");
            ret = jline.solvers.nc.analyzers.Solver_nc_spn_analyzer.solver_nc_spn_analyzer(
                    this.model, this.sn, this.options.copy());
            this.result = ret;
            this.result.solver = this.getName();
            setAvgResults(ret.QN, ret.UN, ret.RN, ret.TN,
                    snGetArvRFromTput(sn, ret.TN, getAvgTputHandles()), null,
                    ret.CN, ret.XN, ret.runtime, ret.method, 1);
            return;
        }

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

        // How this model's finite multiserver stations are represented: Seidmann's
        // approximation or the exact mu(n)=min(n,c) lattice. The shipped "default"
        // reproduces the historical dispatch exactly, so no result moves unless
        // config.multiserver is set. See _kb/06-solver-catalog.md (NC section)
        final String ncMultiserverPolicy = ncMultiserverPolicy(options);

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
                    // Exact LD CoMoM enumerates the per-chain population lattice,
                    // so its cost is unbounded in N while the branch condition
                    // tests only the topology. The budget is the same one
                    // SolverCTMC uses for exact enumeration (6000 states),
                    // applied to prod(1+Nchain); above it Seidmann's comom is the
                    // only affordable option.
                    final double exactLatticeMax = 6000;
                    double latticeSize = 1;
                    if (sn.chains != null && sn.chains.getNumRows() > 0) {
                        for (int c = 0; c < sn.chains.getNumRows(); c++) {
                            double popc = 0;
                            for (int r = 0; r < sn.nclasses && r < sn.chains.getNumCols(); r++) {
                                if (sn.chains.get(c, r) > 0) {
                                    double v = sn.njobs.get(r);
                                    if (!Utils.isInf(v) && !Double.isNaN(v)) {
                                        popc += v;
                                    }
                                }
                            }
                            latticeSize *= (1 + popc);
                        }
                    } else {
                        for (int r = 0; r < sn.nclasses; r++) {
                            double v = sn.njobs.get(r);
                            if (!Utils.isInf(v) && !Double.isNaN(v)) {
                                latticeSize *= (1 + v);
                            }
                        }
                    }
                    if (latticeSize > exactLatticeMax) {
                        options.method = "comom";
                        line_debug(options.verbose, String.format(
                                "NC: default method, population lattice %g exceeds %g, switching to comom",
                                latticeSize, exactLatticeMax));
                    } else if (this.model.hasProductFormSolution() && (sn.lldscaling == null || sn.lldscaling.isEmpty())) {
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
                } else if (ncMultiserverPolicy.equals("lld") && this.model.hasProductFormSolution()) {
                    // config.multiserver="lld" generalizes the exact load-dependent
                    // lattice of the branch above to any closed product-form model,
                    // under the same 6000-state enumeration budget. Off unless asked
                    // for: with the shipped "default" policy this branch never runs
                    // and the model keeps Seidmann's approximation, as it always has
                    Matrix lldFromServers = ncLldFromNservers(sn, 6000.0);
                    if (lldFromServers != null) {
                        sn.lldscaling = lldFromServers;
                        line_debug(options.verbose, "NC: default method, config.multiserver=lld, converted multiserver stations to load-dependent");
                    } else {
                        line_debug(options.verbose, "NC: default method, config.multiserver=lld not applicable (no finite multiserver, non-closed model, or lattice over budget), keeping Seidmann");
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
                // nc_is_oi_model too, not only nc_is_pas_model: the latter demands that
                // BOTH stations be OI/PAS, so a Delay + OI cycle -- the canonical
                // topology Pfqn_oi_is exists to sample -- fell to the guard below.
                if (Solver_nc_pas_is.nc_is_pas_model(sn)
                        || jline.solvers.nc.handlers.Solver_nc_oi.nc_is_oi_model(sn)) {
                    break;   // handled by Solver_nc (Pfqn_pas_is / Pfqn_oi_is)
                }
                // fall through to the "exact" preprocessing
            case "panald":
                // 'panald' is a load-dependent normalizing-constant expansion
                // and needs the same multiserver conversion as "exact"
            case "exact":
                if (!this.model.hasProductFormSolution()) {
                    line_error(mfilename(new Object(){}), "The " + options.method + " method requires the model to have a product-form solution. This model does not have one. You can use Network.hasProductFormSolution() to check before running the solver.");
                } else if (ncMultiserverPolicy.equals("seidmann") && ncHasFiniteMultiserver(sn)) {
                    // config.multiserver="seidmann" asks for Seidmann's approximation
                    // on this arm too, so the conversion below is skipped and the
                    // model goes to the plain nc analyzer. Off by default
                    line_debug(options.verbose, "NC: exact method, config.multiserver=seidmann, keeping Seidmann approximation for multiserver stations");
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
                        if (ret.cacheItemProb != null && ret.cacheItemProb.containsKey(ind)) {
                            cacheNode.setResultItemProb(ret.cacheItemProb.get(ind));
                        }
                    }
                }
                this.model.refreshChains(true);
            } else if (jline.api.sn.SnIsMm1kLoss.snIsMm1kLoss(this.sn)
                    && ("default".equalsIgnoreCase(options.method)
                        || "exact".equalsIgnoreCase(options.method))) {
                // Single-station M/M/1/K with tail drop: exact probability-based
                // loss analysis off the M/M/1/K stationary distribution; see
                // Qsys_mm1k_loss. "mem" asked by name goes past it to the
                // censored GE/GE/1/N block of Solver_nc_mem (mem.blocking).
                // Every other name is refused upstream by ncMethodRefusal, since
                // the closed form reads none.
                line_debug(options.verbose, "NC: single-station M/M/1/K with tail drop, using the qsys_mm1k_loss closed form");
                int queueIst = -1;
                int sourceIst = -1;
                for (int ind = 0; ind < sn.nodetype.size(); ind++) {
                    if (sn.nodetype.get(ind) == NodeType.Queue) {
                        queueIst = (int) sn.nodeToStation.get(ind);
                    } else if (sn.nodetype.get(ind) == NodeType.Source) {
                        sourceIst = (int) sn.nodeToStation.get(ind);
                    }
                }
                int qStateful = (int) sn.stationToStateful.get(queueIst);
                double Vq = sn.visits.get(0).get(qStateful);
                double Kcap = sn.cap.get(queueIst);
                double lambda = sn.rates.get(sourceIst) * Vq;
                double mu = sn.rates.get(queueIst);
                double rho = lambda / mu;
                double Ploss = (Double) jline.api.qsys.Qsys_mm1k_loss
                        .qsys_mm1k_loss(lambda, mu, (int) Math.round(Kcap)).get("lossprob");
                double Tq = lambda * (1.0 - Ploss);   // carried throughput
                double Lsys;
                if (Math.abs(rho - 1.0) < 1e-10) {
                    Lsys = Kcap / 2.0;                // L'Hopital limit at rho=1
                } else {
                    double rKp1 = FastMath.pow(rho, Kcap + 1.0);
                    Lsys = rho / (1.0 - rho) - (Kcap + 1.0) * rKp1 / (1.0 - rKp1);
                }
                ret = new NCResult();
                ret.QN = new Matrix(sn.nstations, sn.nclasses);
                ret.UN = new Matrix(sn.nstations, sn.nclasses);
                ret.RN = new Matrix(sn.nstations, sn.nclasses);
                ret.TN = new Matrix(sn.nstations, sn.nclasses);
                ret.XN = new Matrix(1, sn.nclasses);
                ret.CN = new Matrix(1, sn.nclasses);
                double Rq = Lsys / Tq;                // per-visit response time, by Little
                ret.RN.set(queueIst, 0, Rq);
                ret.QN.set(queueIst, 0, Lsys);
                ret.UN.set(queueIst, 0, Tq / mu);     // single-server utilization
                ret.TN.set(queueIst, 0, Tq);          // carried (effective) rate
                ret.TN.set(sourceIst, 0, lambda);     // offered arrival rate
                ret.XN.set(0, 0, Tq);                 // system throughput = carried rate
                ret.CN.set(0, 0, Rq * Vq);
                ret.lG = 0;
                ret.it = 1;
                iter = 1;
                actualMethod = "mm1k.loss";
                ret.method = actualMethod;
            } else {
                // Ordinary queueing network.
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
                        // Single delay node in FCR - check drop rule. EVERY class
                        // must be dropped, as the reference tests: a region that
                        // discards one class and holds another back is a mixed
                        // system whose blocked class occupies the region while it
                        // waits, so the per-class loss probabilities the loss
                        // network implies are not the ones the model implies.
                        boolean allDrop = true;
                        for (int r = 0; r < sn.nclasses; r++) {
                            if (sn.regionrule.get(0, r) != DropStrategy.Drop.getID()) {
                                allDrop = false;
                            }
                        }
                        if (allDrop) {
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
                // The token gates that used to stand here -- "ms"/"erlangfp" and
                // "rec" off a loss network, "rayint"/"spm" with no Cache node, and
                // the residual Finite Capacity Region on queueing stations -- moved
                // into ncMethodRefusal, which decides them from the same struct
                // before the dispatch begins and which the support gate asks too;
                // their wording is unchanged. The six load-dependent evaluators on
                // an OPEN chain are now refused by the feature set instead
                // (methodFeatureSet drops OpenClass from them): "closed population
                // only" is a rule the registry CAN name, and naming it there is what
                // makes Network.findSolver drop the row rather than report it
                // runnable.
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
                        case "nre":
                        case "comomld":
                        case "panald":
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
        // "ms" names the Manjunath-Sikdar transform of the loss-network analyzer,
        // which is the only place it is admissible; it is listed because
        // solver_nc_lossn_analyzer branches on the token and would otherwise be
        // unreachable.
        return new String[]{
            // "rayint" and "spm" both name the SPM saddle point on a cache, which
            // serves Cache_spm_size once the items carry storage costs. On a retrieval
            // model "rayint" is instead the ray/WKB delayed-hit expansion, admissible
            // only with an infinite-server fetch system; Solver_nc_retrieval_analyzer
            // branches on the token and warns and falls back to "exact" anywhere else.
            // "divdiff" is the divided-difference closed form of Casale (SIGMETRICS
            // 2017); it needs no think time, since a delay would ask for the integral
            // form of Cor. 3.4, and Pfqn_nc refuses one by name. Load-dependent rates
            // ARE served: Pfqn_ncld substitutes the limited load-dependent kernel of
            // Casale-Harrison-Ong (Perform. Eval. 2021), Thm. 1, and reports itself as
            // "divdiff.ld/...".
            "default", "exact", "divdiff", "rayint", "spm", "ms", "erlangfp", "mci", "imci", "ls", "le", "ble", "aghq", "mmint2", "gleint",
            "pana", "panald", "ca", "clw", "kt", "bkt", "lekt", "bk", "bkue", "lc", "lc.ue", "sampling", "is",
            // Chen-O'Cinneide regularization; a Markov chain Monte Carlo estimator of
            // the throughput RATIOS G(N-e_r)/G(N), which supplies no constant of its own
            "mcmc",
            "propfair", "comom", "comomld", "cub",
            // "rgf" (recursion by generating functions, single-class) was the one
            // name in the reference's list with neither a dispatch arm nor an
            // entry here, although Pfqn_rgf sat at the API layer unreached.
            "rgf", "rd", "nrp", "nrl", "nre", "gm", "mem", "ger", "sdr", "sdr.mva",
            // "morrison" is the heavy-usage asymptotic expansion of the generating
            // function for a closed think+DPS network (Npfqn_dps_morrison,
            // Solver_nc_dps_analyzer). It is the DEFAULT on that shape and inadmissible
            // anywhere else, where runAnalyzer refuses it: nothing else in NC can see
            // the DPS weights. Non-product-form, so it returns no lG.
            "morrison",
            // "rec" is the MDD-rec route: on a loss network it is the exact
            // normalizing constant without the residue transform's integrality
            // demand, and on a Petri net it is the only admissible method.
            "rec"
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

    /**
     * Resolves options.config.multiserver into the multiserver handling SolverNC
     * implements: "default" (the historical dispatch), "seidmann" (Seidmann's
     * approximation everywhere, including on "exact") or "lld" (the exact
     * mu(n)=min(n,c) lattice everywhere it is admissible, including on "default").
     *
     * config.multiserver belongs to the general SolverOptions, shared with
     * SolverMVA, which implements approximations SolverNC has no counterpart for
     * ("softmin", "conway", "krzesinski", "suri", "erlang"). Those warn and fall
     * back to "default" rather than erroring, because one options object is
     * commonly reused across solvers.
     */
    private static String ncMultiserverPolicy(SolverOptions options) {
        if (options == null || options.config == null || options.config.multiserver == null) {
            return "default";
        }
        String requested = options.config.multiserver.toLowerCase();
        switch (requested) {
            case "":
            case "default":
                return "default";
            case "seidmann":
                return "seidmann";
            case "lld":
            case "exact":
            case "loaddep":
            case "load-dependent":
                return "lld";
            default:
                line_warning(mfilename(new Object(){}), String.format(
                        "SolverNC does not implement config.multiserver='%s' (it is a SolverMVA "
                        + "approximation); using 'default'. SolverNC accepts 'default', 'seidmann' "
                        + "and 'lld'.", requested));
                return "default";
        }
    }

    /**
     * Exact mu(n)=min(n,c) lattice for a closed model's multiserver stations, the
     * form Solver_ncld consumes. Returns null when the conversion does not apply --
     * no finite multiserver station, an open or mixed model, an lldscaling already
     * installed, or a per-chain population lattice above latticeMax -- in which
     * case the caller keeps Seidmann's approximation.
     */
    private static Matrix ncLldFromNservers(NetworkStruct sn, Double latticeMax) {
        if (sn.lldscaling != null && !sn.lldscaling.isEmpty()) {
            return null;
        }
        if (!ncHasFiniteMultiserver(sn)) {
            return null;
        }
        for (int r = 0; r < sn.nclasses; r++) {
            double v = sn.njobs.get(r);
            if (Utils.isInf(v) || Double.isNaN(v)) {
                return null;
            }
        }
        double Nt = sn.njobs.elementSum();
        if (Utils.isInf(Nt) || Double.isNaN(Nt) || Nt < 1) {
            return null;
        }
        if (latticeMax != null) {
            double latticeSize = 1;
            if (sn.chains != null && sn.chains.getNumRows() > 0) {
                for (int c = 0; c < sn.chains.getNumRows(); c++) {
                    double popc = 0;
                    for (int r = 0; r < sn.nclasses && r < sn.chains.getNumCols(); r++) {
                        if (sn.chains.get(c, r) > 0) {
                            double v = sn.njobs.get(r);
                            if (!Utils.isInf(v) && !Double.isNaN(v)) {
                                popc += v;
                            }
                        }
                    }
                    latticeSize *= (1 + popc);
                }
            } else {
                for (int r = 0; r < sn.nclasses; r++) {
                    double v = sn.njobs.get(r);
                    if (!Utils.isInf(v) && !Double.isNaN(v)) {
                        latticeSize *= (1 + v);
                    }
                }
            }
            if (latticeSize > latticeMax) {
                return null;
            }
        }
        Matrix lldscaling = Matrix.ones(sn.nstations, (int) Nt);
        for (int i = 0; i < sn.nstations; i++) {
            if (sn.nservers.get(i) > 1 && !Utils.isInf(sn.nservers.get(i))) {
                for (int j = 0; j < Nt; j++) {
                    lldscaling.set(i, j, FastMath.min(j + 1, sn.nservers.get(i)));
                }
            }
        }
        return lldscaling;
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
        // finite-buffer MEM path, the slotted route or the single-station
        // M/M/1/K closed form covers the model.
        NetworkStruct s = model.getStruct(false);
        if (!memFiniteBufferPath(s, this.options)
                && !Solver_nc_dt_analyzer.isSlotted(this.options)
                && !jline.api.sn.SnIsMm1kLoss.snIsMm1kLoss(s)) {
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
     * Per-method feature deltas applied to the base NC envelope.
     *
     * ONLY THE RESTRICTIONS A FEATURE NAME CAN CARRY LIVE HERE. A feature set
     * declares what the method ACCEPTS, so it can refuse a model for HAVING a
     * construct and never for lacking one: "closed population only" and "no
     * think time" are expressible by dropping OpenClass and SchedStrategy_INF,
     * while "requires a cache" or "requires a loss network" are not and belong
     * to {@link #ncMethodRefusal}, which supportsModelMethod consults next.
     * Mirrors the MATLAB/python SolverNC.getMethodFeatureSet and the C++
     * nc_feature_set.
     *
     * @param method the concrete method name
     * @return the per-method FeatureSet
     */
    public static FeatureSet methodFeatureSet(String method) {
        FeatureSet featSupported = SolverNC.getFeatureSet();
        if (method == null) {
            return featSupported;
        }
        String m = method.toLowerCase();
        if ("divdiff".equals(m)) {
            // The divided-difference closed form of Casale (SIGMETRICS 2017),
            // Eqs. (15)-(16), covers load-independent queues; a think time would
            // ask for the integral form of Cor. 3.4, which is not implemented, so
            // Pfqn_nc and Pfqn_ncld both refuse one by name. An infinite server is
            // where a think time comes from, so the envelope drops it.
            featSupported.setFalse(new String[]{"SchedStrategy_INF"});
        } else if ("rd".equals(m) || "nrp".equals(m) || "nrl".equals(m) || "nre".equals(m)
                || "comomld".equals(m) || "panald".equals(m)) {
            // The load-dependent normalizing-constant evaluators are reached by
            // Solver_ncld only on its CLOSED branch, where Pfqn_ncld reads the
            // method name. An open chain sends the model to the mixed route
            // (Pfqn_mvaldmx), which never reads it, so every one of these names
            // silently became "ncldmx".
            featSupported.setFalse(new String[]{"OpenClass"});
        } else if ("is".equals(m)) {
            // The sample-an-ordering estimator of Pfqn_is integrates over a closed
            // population simplex; there is no open-class form of it, and Solver_nc
            // refuses one by name. Use "sampling" (Pfqn_mci/Pfqn_ls) for an open
            // or mixed model.
            featSupported.setFalse(new String[]{"OpenClass"});
        }
        // MULTISERVER (registry name since 2026-09-05): the divided-difference
        // closed form of "divdiff" covers load-independent single-server
        // queues, and ncMethodRefusal keeps wording why (a c-server station
        // enters the constant as Seidmann's surrogate delay); every other route
        // folds the count into its own kernel.
        if ("divdiff".equals(m)) {
            featSupported.setFalse(new String[]{"MultiServer"});
        }
        // FINITECAPACITY (registry name since 2026-09-05) is NOT in the base
        // envelope: the product-form routes solve a buffer away, which is what
        // the binding-capacity gate in supportsModelMethod refuses. Two arms
        // honour one: "mem" represents it as a GE/GE/c/0;N queue
        // (memFiniteBufferPath), and the single-station M/M/1/K with tail drop
        // is solved in closed form (Qsys_mm1k_loss, the branch runAnalyzerBody
        // takes on that shape) under "default" and "exact". The shape half of
        // each rule stays structural.
        if ("mem".equals(m) || "default".equals(m) || "exact".equals(m)) {
            featSupported.setTrue(new String[]{"FiniteCapacity"});
        }
        return featSupported;
    }

    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        if (!(this.model instanceof Network)) {
            return null;
        }
        return SolverNC.methodFeatureSet(method);
    }

    /**
     * May METHOD run on this model? Empty string when it may, otherwise the
     * reason it may not, in the words the analyzer refuses with.
     *
     * ONE PREDICATE, TWO CALLERS. runAnalyzer asks it once, ahead of the
     * dispatch, and turns a non-empty answer into an error;
     * {@link #supportsModelMethod} asks it so that Network.findSolver never
     * offers a (solver, method) pair that would raise, and so that SolverAUTO
     * never delegates to one. Two copies of these rules is precisely how the
     * report and the run drift apart, which is the failure this method exists to
     * prevent, so a new rule goes here and not at a call site.
     *
     * Only what the feature registry cannot name lives here; see
     * {@link #methodFeatureSet} for the rules that do have a feature name.
     *
     * @param sn      the network struct
     * @param method  the concrete method name
     * @param options the solver options, read for the discrete-time route only
     * @return empty string if the method may run, else the offending reason
     */
    public static String ncMethodRefusal(NetworkStruct sn, String method, SolverOptions options) {
        return ncMethodRefusal(sn, method, options, true);
    }

    /**
     * The same predicate, told WHICH QUESTION IS BEING ASKED. For two method names
     * the two questions have different answers:
     *
     * <ul>
     *   <li>{@code forReport = true} -- "should model.help() offer this pair?" A pair
     *       that comes back as a table of zeros must not be offered, so the answer is
     *       no.</li>
     *   <li>{@code forReport = false} -- "what does the reference DO when asked for it
     *       by name?" For "mmint2" and "gleint" outside their shape the reference
     *       deliberately WARNS AND RETURNS A ZERO TABLE (pfqn_nc.m, case
     *       {'mmint2','gleint'}: lG = [] and return, unconditionally), and a caller who
     *       names the method keeps that answer.</li>
     * </ul>
     *
     * <p>THE ASYMMETRY IS A RULING, NOT AN OVERSIGHT (2026-07-25, reaffirmed when this
     * gate was added): the report answers "should this be offered" and the run answers
     * "what does the reference do". "comomld" is NOT in that bucket --
     * Pfqn_comomrm_ld raises "The solver accepts at most a single queueing station."
     * natively -- so it is refused on both paths.</p>
     *
     * @param sn        the network struct
     * @param method    the concrete method name
     * @param options   the solver options, read for the discrete-time route only
     * @param forReport true when the caller is the report, false when it is the run
     * @return empty string if the method may run, else the offending reason
     */
    public static String ncMethodRefusal(NetworkStruct sn, String method, SolverOptions options,
                                         boolean forReport) {
        if (sn == null) {
            return "";
        }
        String m = (method == null || method.isEmpty()) ? "default" : method.toLowerCase();

        // The discrete-time route answers for itself: ncIsDtModel decides
        // admissibility on the slot lattice, and every gate below is written
        // about a continuous-time queueing network.
        if (options != null && Solver_nc_dt_analyzer.isSlotted(options)) {
            return "";
        }

        // -- discriminatory processor sharing ------------------------------
        // Morrison's heavy-usage expansion is the ONLY NC route that can see the
        // DPS weights; every other method builds a product-form normalizing
        // constant that silently drops them and answers with the egalitarian-PS
        // network, which is a wrong number rather than a coarse one.
        if (Solver_nc_dps_analyzer.nc_is_dps_model(sn)) {
            if (!"default".equals(m) && !"morrison".equals(m)) {
                return String.format("Method '%s' cannot represent the DPS weights of a "
                        + "discriminatory processor-sharing station; it would return the "
                        + "egalitarian-PS network. Use method 'default' or 'morrison' "
                        + "(Npfqn_dps_morrison), SolverMVA, SolverFLD or SolverCTMC.", method);
            }
            return "";
        }
        if (snHasDPS(sn)) {
            // A DPS station outside Morrison's shape. SchedStrategy_DPS is
            // declared in the feature set because a boolean feature cannot
            // express "this shape only"; this is that imperative half.
            return "SolverNC analyzes a discriminatory processor-sharing station only in "
                    + "the shape Morrison's expansion is derived for: a CLOSED network of exactly "
                    + "two stations, one infinite-server (think) station and one single-server DPS "
                    + "station, exponential service, each class visiting the two equally often. Use "
                    + "SolverMVA, SolverFLD or SolverCTMC for any other DPS model.";
        }
        if ("morrison".equals(m)) {
            // The method named on a model that is not the shape at all -- not
            // even a DPS station in it. Left ungated it reaches no route of its
            // own and falls through to the ordinary normalizing-constant path,
            // which would answer the product-form model UNDER THE CALLER'S LABEL.
            return "Method 'morrison' is the heavy-usage expansion of a CLOSED network of "
                    + "exactly two stations, one infinite-server (think) station and one "
                    + "single-server DPS station with exponential service, which this model is not. "
                    + "Remove the method option to let SolverNC choose, or use SolverMVA, SolverFLD "
                    + "or SolverCTMC.";
        }

        // -- Krzesinski state-dependent routing -----------------------------
        // An SDR model is intercepted by Solver_nc_sdr_analyzer whatever the
        // method says, so reaching the second test means the model declares none.
        if (sn.sdr != null) {
            return "";
        }
        if ("sdr".equals(m) || "sdr.mva".equals(m)) {
            return "Method " + method + " requires state-dependent routing, which this model "
                    + "does not declare.";
        }

        // -- stochastic Petri net -------------------------------------------
        // A net is served only by the MDD-rec route, and none of the gates below
        // -- written about stations, capacities and the queueing-network product
        // form -- says anything about a net. Spn_pf decides its product-form
        // class, by name.
        if (sn.nodetype != null && sn.nodetype.contains(NodeType.Place)) {
            if (!"default".equals(m) && !"rec".equals(m)) {
                return String.format("a stochastic Petri net is solved by the MDD-rec route; "
                        + "method '%s' is a normalizing-constant algorithm for queueing networks. "
                        + "Use 'rec' or 'default'", method);
            }
            return "";
        }

        // -- order-independent stations --------------------------------------
        // Every method other than the four listed reads sn.rates, which holds
        // only the single-job rate mu([r]) of an OI station: the rank rate mu(n)
        // is silently dropped and the answer is that of an ordinary queue.
        if (jline.solvers.nc.handlers.Solver_nc_oi.nc_is_oi_model(sn)) {
            if (!"default".equals(m) && !"exact".equals(m) && !"is".equals(m)
                    && !"sampling".equals(m)) {
                return String.format("Method '%s' cannot represent the rank rate mu(n) of an "
                        + "order-independent station; use method 'default' or 'exact' (Pfqn_ncoi), "
                        + "'is', SolverMVA, or SolverCTMC.", method);
            }
            return "";
        }

        // -- caches -----------------------------------------------------------
        // "rayint" and "spm" both name the SPM saddle point of a cache (and, on a
        // retrieval model, the ray/WKB delayed-hit expansion), so they are
        // admissible here and nowhere else.
        if (sn.nodetype != null && sn.nodetype.contains(NodeType.Cache)) {
            if ("exact".equals(m) && ncIsNoreentrantCache(sn)) {
                CacheNodeParam cp = ncFirstCacheParam(sn);
                // Cache_prob_erec is exact for the exchangeable (RR/FIFO) family
                // only; anything else has to take the approximate route.
                if (cp != null && cp.replacestrat != null
                        && cp.replacestrat != jline.lang.constant.ReplacementStrategy.RR
                        && cp.replacestrat != jline.lang.constant.ReplacementStrategy.FIFO) {
                    return "NC does not support exact solution of the specified cache replacement "
                            + "policy; use the default (approximate) method or SolverCTMC.";
                }
            }
            return "";
        }
        if ("rayint".equals(m) || "spm".equals(m)) {
            return "SolverNC: method " + method + " names the SPM saddle point of a cache and, on "
                    + "a retrieval model, the ray/WKB delayed-hit expansion; this model declares no "
                    + "Cache node.";
        }

        // -- single-station M/M/1/K with tail drop ----------------------------
        // Answered exactly by the probability-based Qsys_mm1k_loss branch under
        // "default" and "exact", and by the censored GE/GE/1/N block of the
        // maximum entropy route under "mem" (runAnalyzerBody lets that name past
        // the closed form). No other name has a route: the closed form reads
        // none, and letting one through would report the UNCONSTRAINED
        // product-form answer under the caller's name.
        if (jline.api.sn.SnIsMm1kLoss.snIsMm1kLoss(sn)) {
            if (!"default".equals(m) && !"exact".equals(m) && !"mem".equals(m)) {
                return String.format("Method '%s' has no route on a single-station M/M/1/K with "
                        + "tail drop, which is answered by the closed form of Qsys_mm1k_loss "
                        + "under 'default' and 'exact' and by the censored GE/GE/1/N block "
                        + "under 'mem'.", method);
            }
            return "";
        }

        // -- loss networks and finite capacity regions -------------------------
        int lossn = ncLossnKind(sn);
        if (lossn == 2) {
            return "";  // "ms", "erlangfp", "rec" and "default" all have a route here
        }
        if (lossn == 1) {
            return "SolverNC does not support finite capacity regions with WAITQ (blocking) "
                    + "policy. Use DROP policy instead.";
        }
        if ("ms".equals(m) || "erlangfp".equals(m)) {
            return "Method " + method + " is admissible only on a loss network (open model, one "
                    + "DROP region holding a single Delay).";
        }
        if ("rec".equals(m)) {
            return "SolverNC: method rec is the MDD-rec route, admissible on a stochastic Petri "
                    + "net or on a loss network (open model, one DROP region holding a single "
                    + "Delay); this model is neither.";
        }
        if (sn.nregions > 0) {
            // NC does not enforce an aggregate region limit on queueing stations;
            // refuse rather than silently return the unconstrained answer.
            return "This model uses a Finite Capacity Region (addRegion) on queueing stations, "
                    + "which is not supported by SolverNC (only the single-Delay loss-network case "
                    + "is). Use SolverJMT, or setCapacity for a single-station limit.";
        }

        // -- PANACEA's domain ---------------------------------------------------
        // Normal usage is a property of the demands rather than of a declared
        // construct, so it has no feature name; an open chain is refused earlier
        // by the closed-population feature set of the load-dependent evaluators.
        //
        // BOTH TOKENS ARE GATED, because Pfqn_ncld evaluates "pana" and
        // "panald" with the SAME Pfqn_panaceald -- its case label is
        // {"pana", "panald"} -- so on a model carrying a rate lattice the
        // load-INDEPENDENT name reaches the load-dependent expansion and throws
        // with it. Off that lattice "pana" takes its own Pfqn_nc arm, which
        // warns and returns an empty constant rather than throwing, so it is left
        // alone there. Class- or joint-dependent scaling diverts the whole model
        // to Solver_nc_conv, which never reads the method at all.
        if (("pana".equals(m) || "panald".equals(m)) && !snHasOpenClasses(sn)) {
            boolean divertedToConv = (sn.cdscaling != null && !sn.cdscaling.isEmpty())
                    || (sn.jdscaling != null && !sn.jdscaling.isEmpty());
            boolean reachesLdKernel = "panald".equals(m)
                    || (sn.lldscaling != null && !sn.lldscaling.isEmpty());
            if (!divertedToConv && reachesLdKernel && !ncIsNormalUsage(sn)) {
                String why = "The model is not in normal usage, so the 'panald' asymptotic "
                        + "expansion does not apply. Use 'exact', 'clw' or an approximate "
                        + "load-dependent method instead.";
                if ("pana".equals(m)) {
                    why = "Method 'pana' reaches the load-dependent kernel on this model, where "
                            + "Pfqn_ncld evaluates it as 'panald'. " + why;
                }
                return why;
            }
        }

        // -- the single-queueing-station recursions ------------------------------
        // Two families are stated for a model with a delay and ONE queueing
        // station, and neither can say so with a feature name: it is a COUNT, and
        // a feature set has no arithmetic. Pfqn_nc states it for "mmint2"/"gleint"
        // in those words and Pfqn_comomrm_ld raises "The solver accepts at most a
        // single queueing station."
        //
        // The count is taken over the CLOSED chains only, and the rule is inactive
        // without a closed population, because Pfqn_nc answers an open network
        // with the exact open formulas BEFORE its method switch -- the method name is
        // never read there, so a purely open model with three queues runs these
        // names correctly today and must go on doing so.
        // "mmint2" and "gleint" are gated for the REPORT ONLY: Pfqn_nc answers them
        // with an empty constant and the caller renders a table of zeros, which is a
        // pair the report must not offer and a run the reference nonetheless
        // performs. See the forReport parameter above.
        if ("comomld".equals(m) || (forReport && ("mmint2".equals(m) || "gleint".equals(m)))) {
            int nq = ncClosedQueueingStations(sn);
            if (nq > 1) {
                if ("comomld".equals(m)) {
                    return "Method 'comomld' is the load-dependent CoMoM recursion, and "
                            + "Pfqn_comomrm_ld accepts at most a single queueing station; this "
                            + "model has " + nq + ".";
                }
                return "The '" + method + "' method requires a model with a delay and a single "
                        + "queueing station; this model has " + nq + ".";
            }
        }

        // -- "exact" outside its domain -----------------------------------------
        if ("exact".equals(m)) {
            boolean multiserver = ncHasFiniteMultiserver(sn);
            boolean hasOpen = snHasOpenClasses(sn);
            if (multiserver && hasOpen) {
                return "NC solver cannot provide exact solutions for open or mixed queueing "
                        + "networks. Remove the 'exact' option.";
            }
            boolean scaling = (sn.lldscaling != null && !sn.lldscaling.isEmpty())
                    || (sn.cdscaling != null && !sn.cdscaling.isEmpty())
                    || (sn.jdscaling != null && !sn.jdscaling.isEmpty());
            boolean fractional = false;
            for (int r = 0; r < sn.njobs.getNumElements(); r++) {
                double v = sn.njobs.get(r);
                if (!Utils.isInf(v) && FastMath.abs(v - FastMath.floor(v)) > GlobalConstants.FineTol) {
                    fractional = true;
                }
            }
            if ((scaling || multiserver) && fractional) {
                // The load-dependent analyzer interpolates a fractional population
                // between the two integer neighbours, which is an approximation, so
                // it refuses the exactness the caller asked for by name.
                return "NC load-dependent solver cannot provide exact solutions for fractional "
                        + "populations.";
            }
        }
        return "";
    }

    /**
     * How many queueing (non-infinite-server) stations carry demand from a CLOSED
     * chain, which is the row count L reaches Pfqn_nc and Pfqn_comomrm_ld with
     * once the delay rows have been folded into Z and the zero-demand rows
     * dropped. Zero when the model has no closed population at all.
     *
     * @param sn the network struct
     * @return the number of closed-demand queueing stations
     */
    private static int ncClosedQueueingStations(NetworkStruct sn) {
        Ret.snGetDemands dem = snGetDemandsChain(sn);
        Matrix Lchain = dem.Dchain;
        Matrix Nchain = dem.Nchain;
        if (Lchain == null || Nchain == null) {
            return 0;
        }
        int C = Nchain.getNumElements();
        boolean[] closed = new boolean[C];
        boolean anyClosed = false;
        for (int c = 0; c < C; c++) {
            double v = Nchain.get(c);
            closed[c] = !Utils.isInf(v) && !Double.isNaN(v) && v > 0;
            if (closed[c]) {
                anyClosed = true;
            }
        }
        if (!anyClosed) {
            return 0;
        }
        int nq = 0;
        for (int i = 0; i < sn.nstations; i++) {
            if (Utils.isInf(sn.nservers.get(i))) {
                continue;
            }
            for (int c = 0; c < C; c++) {
                if (closed[c] && FastMath.abs(Lchain.get(i, c)) > GlobalConstants.FineTol) {
                    nq++;
                    break;
                }
            }
        }
        return nq;
    }

    /** The Source-Cache-Sink model Solver_nc_cache_analyzer serves. */
    private static boolean ncIsNoreentrantCache(NetworkStruct sn) {
        if (sn.nclosedjobs != 0 || sn.nodetype == null || sn.nodetype.size() != 3) {
            return false;
        }
        return sn.nodetype.contains(NodeType.Source) && sn.nodetype.contains(NodeType.Cache)
                && sn.nodetype.contains(NodeType.Sink);
    }

    /** The parameters of the first Cache node, or null when the model has none. */
    private static CacheNodeParam ncFirstCacheParam(NetworkStruct sn) {
        if (sn.nodeparam == null || sn.nodes == null) {
            return null;
        }
        for (int ind = 0; ind < sn.nodetype.size(); ind++) {
            if (sn.nodetype.get(ind) == NodeType.Cache && ind < sn.nodes.size()) {
                NodeParam np = sn.nodeparam.get(sn.nodes.get(ind));
                if (np instanceof CacheNodeParam) {
                    return (CacheNodeParam) np;
                }
            }
        }
        return null;
    }

    /**
     * The loss-network verdict: 0 for neither, 1 for the shape without the DROP
     * rule, 2 for the loss network Solver_nc_lossn_analyzer solves.
     *
     * The shape is an OPEN model with a single finite capacity region whose only
     * member is an infinite server. It is a LOSS network when the admission rule
     * DROPS every class: a region of that shape under WAITQ (or any blocking
     * rule) holds the arrival back instead of discarding it, which keeps the job
     * in the region while it waits and is a queueing phenomenon the Erlang loss
     * model has no state for. The two verdicts are separate because they have
     * different remedies -- switching the rule to DROP makes the first solvable
     * here, while a region on queueing stations needs a solver that carries the
     * region as state.
     *
     * @param sn the network struct
     * @return 0, 1 or 2 as above
     */
    public static int ncLossnKind(NetworkStruct sn) {
        if (sn.nregions != 1 || snHasClosedClasses(sn) || sn.region == null) {
            return 0;
        }
        Matrix regionMatrix = sn.region.get(0);
        if (regionMatrix == null) {
            return 0;
        }
        int stationInFCR = -1;
        int stationCount = 0;
        for (int i = 0; i < sn.nstations; i++) {
            boolean hasConstraint = false;
            for (int r = 0; r < sn.nclasses; r++) {
                if (regionMatrix.get(i, r) >= 0) {
                    hasConstraint = true;
                    break;
                }
            }
            if (regionMatrix.getNumCols() > sn.nclasses && regionMatrix.get(i, sn.nclasses) >= 0) {
                hasConstraint = true;
            }
            if (hasConstraint) {
                stationInFCR = i;
                stationCount++;
            }
        }
        if (stationCount != 1 || !Utils.isInf(sn.nservers.get(stationInFCR))) {
            return 0;
        }
        if (sn.regionrule == null) {
            return 1;
        }
        // ALL classes, not merely one: a region that discards one class and holds
        // another back is a mixed system whose blocked class occupies the region
        // while it waits, so the per-class loss probabilities the Erlang fixed
        // point returns would not be the ones the model implies.
        for (int r = 0; r < sn.nclasses; r++) {
            if (sn.regionrule.get(0, r) != DropStrategy.Drop.getID()) {
                return 1;
            }
        }
        return 2;
    }

    /**
     * Is the closed model in NORMAL USAGE, the domain of the Mitra-McKenna
     * PANACEA asymptotic expansion (J. ACM 33(3), 1986)?
     *
     * Normal usage asks that every queueing centre absorb the load the think
     * stations offer it: with rho_j0 = Ztot(j) the aggregate think demand of
     * chain j, r_ij = L_ij / rho_j0 and mu_i(Ntot) the saturation rate,
     * alpha_i = 1 - (sum_j N_j r_ij) / mu_i(Ntot) &gt; 0 at every centre i.
     * Outside it the {phi(n)} series DIVERGES, which is why Pfqn_panaceald
     * returns NaN there and Pfqn_ncld turns that NaN into a refusal. It is a
     * property of the demands and not of a declared construct, so it has no
     * feature-registry name.
     *
     * The rates are the ones Solver_ncld would build: 1 for an ordinary single
     * server, min(n,c) for a finite multiserver (the conversion runAnalyzer
     * performs on the "panald" arm), and the declared lldscaling row when the
     * model sets one. An infinite server is a think station and feeds Ztot.
     *
     * @param sn the network struct
     * @return true when the expansion applies
     */
    public static boolean ncIsNormalUsage(NetworkStruct sn) {
        if (snHasOpenClasses(sn)) {
            return true;
        }
        Ret.snGetDemands dem = snGetDemandsChain(sn);
        Matrix Lchain = dem.Dchain;
        Matrix Nchain = dem.Nchain;
        if (Lchain == null || Nchain == null) {
            return true;
        }
        int C = Nchain.getNumElements();
        double NtD = 0;
        for (int c = 0; c < C; c++) {
            double v = Nchain.get(c);
            if (!Utils.isInf(v) && !Double.isNaN(v)) {
                NtD += v;
            }
        }
        int Nt = (int) FastMath.round(NtD);
        if (Nt < 1) {
            return true;  // the empty network: G = 1, nothing to expand
        }
        int M = sn.nstations;

        double[][] mu = new double[M][Nt];
        boolean haveLld = sn.lldscaling != null && !sn.lldscaling.isEmpty();
        for (int i = 0; i < M; i++) {
            for (int n = 0; n < Nt; n++) {
                if (haveLld) {
                    int col = FastMath.min(n, sn.lldscaling.getNumCols() - 1);
                    mu[i][n] = sn.lldscaling.get(i, col);
                } else if (!Utils.isInf(sn.nservers.get(i)) && sn.nservers.get(i) > 1) {
                    mu[i][n] = FastMath.min(n + 1, sn.nservers.get(i));
                } else {
                    mu[i][n] = 1.0;
                }
            }
        }

        double[] Ztot = new double[C];
        for (int i = 0; i < M; i++) {
            if (Utils.isInf(sn.nservers.get(i))) {
                for (int c = 0; c < C; c++) {
                    Ztot[c] += Lchain.get(i, c);
                }
            }
        }
        for (int c = 0; c < C; c++) {
            if (Nchain.get(c) > 0 && Ztot[c] <= 0) {
                // no think station on the route of a populated chain: the
                // expansion parameter rho_j0 is undefined
                return false;
            }
        }

        for (int i = 0; i < M; i++) {
            if (Utils.isInf(sn.nservers.get(i))) {
                continue;
            }
            double muK = mu[i][Nt - 1];
            if (!(muK > 0) || Double.isNaN(muK) || Utils.isInf(muK)) {
                return false;
            }
            double lambda = 0;
            for (int c = 0; c < C; c++) {
                if (Ztot[c] > 0) {
                    lambda += Lchain.get(i, c) / Ztot[c] * Nchain.get(c);
                }
            }
            if (1.0 - lambda / muK <= 0) {
                return false;
            }
        }
        // Every centre cleared the test; a model with no queueing centre at all
        // reaches here too, and there the delay-only constant is exact.
        return true;
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
        return supportsModelMethod(method, true);
    }

    /**
     * The same gate, told WHICH QUESTION IS BEING ASKED. The report asks "should this
     * pair be offered" and the run asks "what does the reference do"; they differ only
     * for "mmint2"/"gleint", which the reference performs by name and answers with a
     * zero table. See {@link #ncMethodRefusal(NetworkStruct, String, SolverOptions,
     * boolean)}.
     *
     * <p>The distinction has to be plumbed because the shared runAnalyzerChecks gate
     * lives on the RUN path yet reaches this method through the same one-argument call
     * findSolver makes, so there is no other way to tell the two callers apart.</p>
     *
     * @param method    the concrete method name
     * @param forReport true when the caller is the report, false when it is the run
     * @return empty string if supported, else the offending reason
     */
    public String supportsModelMethod(String method, boolean forReport) {
        // Discrete-time (slotted) route: a finite buffer on a Bernoulli server
        // is the loss system of Daduna's corollary 2.8, which
        // Solver_nc_dt_analyzer solves exactly, so the coarse feature gate must
        // not fire. ncIsDtModel performs the real admissibility check and
        // reports a precise reason.
        if (Solver_nc_dt_analyzer.isSlotted(this.options)) {
            NetworkStruct sdt = this.sn != null ? this.sn : this.model.getStruct(false);
            Solver_nc_dt_analyzer.DtModel dt = Solver_nc_dt_analyzer.ncIsDtModel(sdt);
            if (!"none".equals(dt.kind)) {
                return "";
            }
            return "options.config.slotted is set but " + dt.reason;
        }
        if ("mem".equals(method)) {
            NetworkStruct s = this.sn != null ? this.sn : this.model.getStruct(false);
            // MEM now supports both open (Section 3.2) and closed (Section 3.3)
            // networks, so validate against the model's actual class composition.
            String r = memUnsupportedReason(s, model.hasOpenClasses(), model.hasClosedClasses());
            return r == null ? "" : r;
        }
        // The per-method envelope, not the flat one: "divdiff" drops
        // SchedStrategy_INF and the closed-population evaluators drop OpenClass,
        // which is how "no think time" and "closed only" are said in the registry.
        String featReason = FeatureSet.supportsReason(SolverNC.methodFeatureSet(method),
                ((Network) this.model).getUsedLangFeatures());
        if (!featReason.isEmpty()) {
            return featReason;
        }
        // A BINDING BUFFER IS A PER-MODEL RULE THE REGISTRY CANNOT NAME, and it
        // belongs here as well as in supports(Network). Only the latter carried
        // it, so Solver.supports refused the model as a whole while every nc
        // method name stayed in listValidMethods and findSolver -- 'nc.bk' among them
        // -- offering methods that cannot represent a buffer. Exempt for the
        // same three paths supports() exempts: MEM carries finite buffers
        // explicitly, the slotted route is Daduna's loss system, and the
        // single-station M/M/1/K with tail drop has the Qsys_mm1k_loss closed
        // form (the names it does not serve are refused by ncMethodRefusal).
        NetworkStruct scap0 = this.sn != null ? this.sn : ((Network) this.model).getStruct(false);
        if (!memFiniteBufferPath(scap0, this.options)
                && !Solver_nc_dt_analyzer.isSlotted(this.options)
                && !jline.api.sn.SnIsMm1kLoss.snIsMm1kLoss(scap0)) {
            NetworkStruct scap = this.sn != null ? this.sn : ((Network) this.model).getStruct(false);
            String capReason = NetworkSolver.bindingCapacityReason((Network) this.model, scap,
                                                                   "SolverNC");
            if (capReason != null) {
                return capReason;
            }
        }
        String exactReason = exactnessReason((Network) this.model, method);
        if (!exactReason.isEmpty()) {
            return exactReason;
        }
        // The structural per-method rules the feature registry cannot name: which
        // route a method has on THIS model, and whether it exists at all.
        // ncMethodRefusal is the single copy of them, asked here and by
        // runAnalyzer, so the report and the run cannot disagree.
        NetworkStruct sfull = this.sn != null ? this.sn : this.model.getStruct(false);
        return ncMethodRefusal(sfull, method, this.options, forReport);
    }

    /**
     * Product-form precondition of the normalizing-constant methods, the same
     * rule runAnalyzer enforces at solve time. Only "exact", "is" and
     * "panald" require it -- the other methods fall back to Seidmann's
     * comom on a non-product-form model -- and "is" on a pass-and-swap model
     * is exempt (Pfqn_pas_is). Product form has no registry feature name, so
     * the check cannot live in getMethodFeatureSet.
     *
     * @param model  the network model
     * @param method the concrete method name
     * @return empty string if the method may run, else the offending reason
     */
    public static String exactnessReason(Network model, String method) {
        if (method == null) {
            return "";
        }
        boolean needsProductForm = "exact".equalsIgnoreCase(method)
                || "is".equalsIgnoreCase(method) || "panald".equalsIgnoreCase(method);
        if (!needsProductForm || model.hasProductFormSolution()) {
            return "";
        }
        // A LOSS NETWORK (open, one DROP region holding a single Delay) IS product
        // form -- the truncated Poisson law the residue transform of
        // solver_nc_lossn_analyzer evaluates exactly under "exact" -- but
        // snHasBlocking reads any region as blocking, so hasProductFormSolution
        // says no. The shape is exempted here and at the "exact" arm of
        // runAnalyzer alike, through the one predicate both ask.
        if (ncLossnKind(model.getStruct(false)) == 2) {
            return "";
        }
        if ("is".equalsIgnoreCase(method)
                && (Solver_nc_pas_is.nc_is_pas_model(model.getStruct(false))
                    || jline.solvers.nc.handlers.Solver_nc_oi.nc_is_oi_model(model.getStruct(false)))) {
            return "";
        }
        return "method '" + method + "' requires a product-form solution; use "
                + "method 'comom' or SolverCTMC instead";
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
            if (tok.equals("mci") || tok.equals("imci") || tok.equals("ls") || tok.equals("sampling")
                    || tok.equals("is") || tok.equals("mcmc")) {
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
     * First-probe-time summary of the response time CDF, one scalar per station
     * and class.
     *
     * NOT A DISTRIBUTION. This was called getCdfRespT(AvgHandle...) and was
     * renamed because that name made it an OVERLOAD of the base
     * {@link NetworkSolver#getCdfRespT(AvgHandle)} rather than an override --
     * varargs erases to AvgHandle[], so one handle reached the base
     * implementation and zero handles reached this class, two implementations
     * behind one name selected by argument count. For the CDF itself, one
     * (T x 2) matrix of [F(t) t] per station and class as MATLAB's cell array
     * carries, call {@link #getCdfRespT()} and read its {@code cdfData}.
     *
     * @return first-probe-time CDF summary, indexed by station and class
     */
    public Matrix getCdfRespTFirstProbe() {
        NetworkStruct sn0 = getStruct();
        List<List<Matrix>> full = computeCdfRespT();
        if (full == null) return new Matrix(0, 0);
        Matrix summary = new Matrix((int) sn0.nstations, (int) sn0.nclasses);
        for (int i = 0; i < sn0.nstations; i++) {
            for (int j = 0; j < sn0.nclasses; j++) {
                Matrix cdf = full.get(i).get(j);
                if (cdf != null && cdf.getNumRows() > 0 && cdf.getNumCols() > 0) {
                    summary.set(i, j, Math.abs(cdf.get(0, 0)));
                }
            }
        }
        return summary;
    }

    /**
     * The response time CDF at every probe time, indexed [station][class], each
     * entry a (T x 2) matrix of [F(t) t] as MATLAB's cell array holds it, or
     * null where the station-class pair has no distribution.
     *
     * Shared by both getCdfRespT overloads so the work is done once and neither
     * can drift from the other; returns null when the model has no FCFS station,
     * which is the case the reference warns about.
     */
    private List<List<Matrix>> computeCdfRespT() {
        if (GlobalConstants.DummyMode) {
            return null;
        }

        long startTimeMillis = System.nanoTime();
        NetworkStruct sn = getStruct();
        
        // Get algorithm configuration
        String algorithm = options.method != null ? options.method : "exact";
        
        // [station][class], each entry a (T x 2) [F(t) t] matrix or null.
        // STATION-indexed, matching NetworkSolver:6629 and SolverJMT:635; this
        // method used to index the outer dimension by NODE, which every
        // cdfData consumer reads as a station.
        List<List<Matrix>> RD = null;

        try {
            // Get product form parameters
            Ret.snGetProductFormParams params = snGetProductFormParams(sn);
            Matrix D = params.D;  // Service demands
            Matrix N = params.N;  // Population vector  
            Matrix Z = params.Z;  // Think times
            Matrix S = params.S;  // Number of servers
            
            // THREE INDEX SPACES MEET HERE, and using one index for all three is
            // how this method came to be dead code. Named explicitly:
            //   STATION space  sn.stations / sn.sched / sn.rates, 0..nstations-1
            //   QUEUE space    the ROWS of D and S, one per Queue-type NODE in
            //                  node order (snGetProductFormParams builds them
            //                  from sn.nodetype == Queue)
            //   NODE space     the rows of the returned RD, 0..nnodes-1
            // pfqn_stdf applies ONE index to L, S and rates, so whatever is
            // passed as fcfsNodes must be in QUEUE space and rates must be too.

            // node -> QUEUE-space row, in the order snGetProductFormParams uses
            Map<Integer, Integer> queueRowOfNode = new HashMap<Integer, Integer>();
            int queueRows = 0;
            for (int nd = 0; nd < sn.nodetype.size(); nd++) {
                if (sn.nodetype.get(nd) == NodeType.Queue) {
                    queueRowOfNode.put(nd, queueRows++);
                }
            }

            List<Integer> fcfsQueueRows = new ArrayList<>();   // QUEUE space
            List<Integer> fcfsStationsList = new ArrayList<>(); // STATION space
            List<Integer> delayStationsList = new ArrayList<>();// STATION space

            for (int ist = 0; ist < sn.nstations; ist++) {
                // sn.sched is Map<Station,SchedStrategy>; get(int) on it silently
                // returns null, which left this list empty on EVERY model and made
                // the whole method return an empty matrix. Use the keyed lookup the
                // rest of the JAR uses (e.g. SnHasHomogeneousScheduling).
                SchedStrategy sched = sn.sched.get(sn.stations.get(ist));
                if (sched == SchedStrategy.FCFS) {
                    Integer qrow = queueRowOfNode.get((int) sn.stationToNode.get(ist));
                    if (qrow != null) {
                        fcfsQueueRows.add(qrow);
                        fcfsStationsList.add(ist);
                    }
                } else if (sched == SchedStrategy.INF) {
                    delayStationsList.add(ist);
                }
            }

            if (!fcfsQueueRows.isEmpty()) {
                // Calculate time horizon
                double totalPop = N.elementSum();

                // Rates in QUEUE space: one row per queueing station, NOT compacted
                // to the FCFS ones, so that the single pfqn_stdf index addresses
                // rates, L and S alike.
                Matrix fcfsRates = new Matrix(queueRows, sn.nclasses);
                for (int ist = 0; ist < sn.nstations; ist++) {
                    Integer qrow = queueRowOfNode.get((int) sn.stationToNode.get(ist));
                    if (qrow == null) continue;
                    for (int j = 0; j < sn.nclasses; j++) {
                        fcfsRates.set(qrow, j, sn.rates.get(ist, j));
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
                
                // fcfsNodes in QUEUE space, matching D, S and fcfsRates
                Matrix fcfsNodes = new Matrix(fcfsQueueRows.size(), 1);
                for (int i = 0; i < fcfsQueueRows.size(); i++) {
                    fcfsNodes.set(i, fcfsQueueRows.get(i).doubleValue());
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
                
                // Initialize result structure, STATION x CLASS
                RD = new ArrayList<List<Matrix>>();
                for (int i = 0; i < sn.nstations; i++) {
                    List<Matrix> row = new ArrayList<Matrix>();
                    for (int j = 0; j < sn.nclasses; j++) row.add(null);
                    RD.add(row);
                }
                
                // Process FCFS results. RDout is indexed in QUEUE space (pfqn_stdf
                // allocates Matrix[M][R] and writes RD[kIdx][r]), so it is read at
                // the queue row and STORED at the station row. The WHOLE (T x 2)
                // curve is kept, as MATLAB's RD{k,r} does.
                if (RDout != null && RDout.length > 0) {
                    for (int i = 0; i < fcfsQueueRows.size(); i++) {
                        int qrow = fcfsQueueRows.get(i);
                        int ist = fcfsStationsList.get(i);
                        for (int j = 0; j < sn.nclasses; j++) {
                            if (qrow < RDout.length && j < RDout[qrow].length && RDout[qrow][j] != null) {
                                Matrix cdf = RDout[qrow][j];
                                if (cdf.getNumRows() > 0 && cdf.getNumCols() > 0) {
                                    // real() in MATLAB: strip complex round-off
                                    Matrix cdfAbs = new Matrix(cdf.getNumRows(), cdf.getNumCols());
                                    for (int t = 0; t < cdf.getNumRows(); t++) {
                                        cdfAbs.set(t, 0, Math.abs(cdf.get(t, 0)));
                                        if (cdf.getNumCols() > 1) cdfAbs.set(t, 1, cdf.get(t, 1));
                                    }
                                    RD.get(ist).set(j, cdfAbs);
                                }
                            }
                        }
                    }
                }

                // Process delay results. sn.proc is Map<Station,Map<JobClass,..>>,
                // so it is keyed by the objects and NOT by an int -- the same defect
                // the sched scan had, and it would have thrown a NullPointerException
                // here the moment that one was fixed.
                for (int i = 0; i < delayStationsList.size(); i++) {
                    int ist = delayStationsList.get(i);
                    int nodeId = (int) sn.stationToNode.get(ist);
                    Map<JobClass, MatrixCell> stProc = sn.proc.get(sn.stations.get(ist));
                    if (stProc == null) continue;
                    for (int j = 0; j < sn.nclasses; j++) {
                        MatrixCell procCell = stProc.get(sn.jobclasses.get(j));
                        if (procCell != null && !procCell.isEmpty()) {
                            Matrix cdfResult = map_cdf(procCell, tset.transpose());
                            // map_cdf returns a 1 x T ROW (Map_cdf.java:31), not
                            // a column: reading it by rows keeps one point.
                            int nt = cdfResult.length();
                            if (nt > 0) {
                                // whole curve as [F(t) t], the layout MATLAB's
                                // getCdfRespT.m:37 builds for a delay station
                                Matrix cdf = new Matrix(nt, 2);
                                for (int t = 0; t < nt; t++) {
                                    cdf.set(t, 0, cdfResult.get(t));
                                    cdf.set(t, 1, tset.get(t));
                                }
                                RD.get(ist).set(j, cdf);
                            }
                        }
                    }
                }
                
                long endTimeMillis = System.nanoTime();
                double runtime = (endTimeMillis - startTimeMillis) / 1000000000.0;

                // setDistribResults takes the (M x K) summary; the full curves
                // travel back through the return value.
                Matrix summary = new Matrix((int) sn.nstations, (int) sn.nclasses);
                for (int i = 0; i < sn.nstations; i++) {
                    for (int j = 0; j < sn.nclasses; j++) {
                        Matrix cdf = RD.get(i).get(j);
                        if (cdf != null && cdf.getNumRows() > 0) summary.set(i, j, cdf.get(0, 0));
                    }
                }
                setDistribResults(summary, runtime);

            } else {
                line_warning(mfilename(new Object(){}), "getCdfRespT applies only to FCFS nodes.");
            }
            
        } catch (Exception e) {
            line_error(mfilename(new Object(){}), "Error in getCdfRespT: " + e.getMessage());
        }
        
        return RD;
    }

    /**
     * Get the response time CDF at FCFS and delay stations.
     *
     * {@code cdfData} is indexed [station][class], each entry a (T x 2) matrix
     * of [F(t) t], matching MATLAB's RD cell array and the layout SolverJMT
     * fills. Entries are null where the pair has no distribution. This used to
     * compute the curves and DISCARD them, returning the empty structure
     * initializeCdfData allocates.
     *
     * @return the response time distributions, station by class
     */
    public DistributionResult getCdfRespT() {
        NetworkStruct sn = getStruct();
        DistributionResult distResult =
            new DistributionResult((int) sn.nstations, (int) sn.nclasses, "response_time");
        List<List<Matrix>> full = computeCdfRespT();
        if (full == null) return distResult;
        distResult.cdfData = full;
        return distResult;
    }

    /**
     * Response time CDF, ignoring the handle argument as the reference does.
     *
     * THIS GENUINELY OVERRIDES {@link NetworkSolver#getCdfRespT(AvgHandle)}.
     * Before it existed, SolverNC declared only a varargs method, which erases
     * to AvgHandle[] and therefore did NOT override the base: a caller passing
     * one handle silently received the base class's EXPONENTIAL FIT of the mean
     * response time, while a caller passing none received this class's exact
     * law. Both returned plausible numbers, so nothing looked wrong. The two now
     * agree because both delegate to {@link #computeCdfRespT()}.
     *
     * @param R response time handles, unused -- the CDF is computed at every
     *          station and class, as MATLAB's getCdfRespT.m does
     * @return the response time distributions, station by class
     */
    @Override
    public DistributionResult getCdfRespT(AvgHandle R) {
        return getCdfRespT();
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

    /**
     * Mean duration of the busy period of order n for the subnetwork made of
     * the given stations, that is the time from the instant a job entering the
     * subnetwork finds n-1 jobs in it up to the next instant when fewer than n
     * remain.
     *
     * <p>H. Daduna, "Busy Periods for Subnetworks in Stochastic Networks: Mean
     * Value Analysis", J. ACM 35(3), 1988. The result is exact on the
     * single-chain product-form class and, by the insensitivity of Section 5 of
     * that paper, depends on the service processes only through their mean
     * rates.
     *
     * <p>The Java twin of MATLAB's {@code @SolverNC/getAvgBusyPeriod.m} and of
     * {@code solver_nc_busyp} in the C++ port, which this follows step for
     * step. It had no Java spelling until now, so {@code -a busyperiod} on the
     * CLI and every {@code lang="java"} caller had to fall back to SolverLDES,
     * i.e. to a simulation of a quantity there is a transform for.
     *
     * @param stations zero-based STATION indexes forming the subnetwork; a
     *                 non-empty proper subset of the stations
     * @param orders   busy period orders, 1 &lt;= n &lt;= population for a closed model
     * @return one mean duration per requested order
     */
    public double[] getAvgBusyPeriod(int[] stations, int[] orders) {
        if (GlobalConstants.DummyMode) {
            double[] dummy = new double[orders.length];
            Arrays.fill(dummy, Double.NaN);
            return dummy;
        }
        NetworkStruct snb = this.model.getStruct(false);
        if (snb.nchains > 1) {
            line_error(mfilename(new Object(){}),
                "The busy period of a subnetwork is defined for single-chain models only. "
                + "Section 5 of Daduna (1988) sketches the multichain extension, which is "
                + "not implemented.");
        }
        final int M = snb.nstations;
        final int K = snb.nclasses;
        Ret.snGetDemands dem = SnGetDemandsChain.snGetDemandsChain(snb);
        Pair<Matrix, Matrix> rtv = SnRtStations.snRtStations(snb);
        final Matrix rtst = rtv.getLeft();
        final Matrix Vst = rtv.getRight();

        // station-to-station routing of the chain: the class-level probabilities
        // weighted by the class visits, which is exact because it is a flow balance
        Matrix Pst = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            double vtot = 0.0;
            for (int r = 0; r < K; r++) {
                vtot += Vst.get(i, r);
            }
            for (int j = 0; j < M; j++) {
                double flow = 0.0;
                for (int r = 0; r < K; r++) {
                    for (int t = 0; t < K; t++) {
                        flow += Vst.get(i, r) * rtst.get(i * K + r, j * K + t);
                    }
                }
                Pst.set(i, j, vtot > 0 ? flow / vtot : 0.0);
            }
        }

        final Matrix STchain = dem.STchain;
        final Matrix Vchain = dem.Vchain;
        final Matrix lld = snb.lldscaling;
        final Matrix nservers = snb.nservers;
        // mu(j,k): the load-dependent rate of station j holding k jobs, in the
        // solver_ncld precedence -- the infinite server first, then the declared
        // scaling table, then the multiserver staircase. A table shorter than k
        // keeps its last entry.
        final Pfqn_busyp.RateFunction rateOf = new Pfqn_busyp.RateFunction() {
            @Override
            public double rate(int j, int k) {
                double scaling;
                double servers = nservers.get(j, 0);
                if (Double.isInfinite(servers)) {
                    scaling = k;
                } else if (lld != null && !lld.isEmpty() && lld.getNumCols() > 0) {
                    scaling = lld.get(j, Math.min(k, lld.getNumCols()) - 1);
                } else {
                    scaling = Math.min((double) k, servers);
                }
                return scaling / STchain.get(j, 0);
            }
        };

        double N = (Vchain == null || dem.Nchain == null || dem.Nchain.isEmpty())
                ? Double.POSITIVE_INFINITY : dem.Nchain.get(0);
        if (!Double.isFinite(N)) {
            // the Source is not a node of the Jackson network of the paper: its
            // outflow is the external stream gamma
            int source = -1;
            for (int i = 0; i < M; i++) {
                if (snb.nodetype.get((int) snb.stationToNode.get(i)) == NodeType.Source) {
                    source = i;
                    break;
                }
            }
            if (source < 0) {
                line_error(mfilename(new Object(){}), "An open model must own a Source station.");
            }
            for (int t = 0; t < stations.length; t++) {
                if (stations[t] == source) {
                    line_error(mfilename(new Object(){}),
                        "The Source cannot belong to the subnetwork.");
                }
            }
            double lambda = 0.0;
            for (int r = 0; r < K; r++) {
                double rate = snb.rates.get(source, r);
                if (!Double.isNaN(rate) && Double.isFinite(rate)) {
                    lambda += rate;
                }
            }
            final int[] keep = new int[M - 1];
            int[] remap = new int[M];
            int at = 0;
            for (int i = 0; i < M; i++) {
                if (i != source) {
                    remap[i] = at;
                    keep[at++] = i;
                }
            }
            Matrix alpha = new Matrix(1, keep.length);
            Matrix gamma = new Matrix(1, keep.length);
            Matrix P = new Matrix(keep.length, keep.length);
            for (int i = 0; i < keep.length; i++) {
                alpha.set(0, i, lambda * Vchain.get(keep[i], 0) / Vchain.get(source, 0));
                gamma.set(0, i, lambda * Pst.get(source, keep[i]));
                for (int j = 0; j < keep.length; j++) {
                    P.set(i, j, Pst.get(keep[i], keep[j]));
                }
            }
            int[] mapped = new int[stations.length];
            for (int t = 0; t < stations.length; t++) {
                mapped[t] = remap[stations[t]];
            }
            Pfqn_busyp.RateFunction keptRate = new Pfqn_busyp.RateFunction() {
                @Override
                public double rate(int j, int k) {
                    return rateOf.rate(keep[j], k);
                }
            };
            return Pfqn_busyp.pfqn_busyp(alpha, keptRate, P, Double.POSITIVE_INFINITY,
                    mapped, orders, gamma);
        }

        Matrix alpha = new Matrix(1, M);
        for (int i = 0; i < M; i++) {
            alpha.set(0, i, Vchain.get(i, 0));
        }
        return Pfqn_busyp.pfqn_busyp(alpha, rateOf, Pst, N, stations, orders, null,
                Pfqn_busyp.DEFAULT_TOL);
    }

    /** Single-order form of {@link #getAvgBusyPeriod(int[], int[])}. */
    public double getAvgBusyPeriod(int[] stations, int order) {
        return getAvgBusyPeriod(stations, new int[]{order})[0];
    }


    /**
     * Solves through a model transformation, with NC as its own inner solve.
     *
     * <p>The inner options carry {@code transform='none'} and a raised depth, so
     * a transformed submodel cannot re-enter the driver. An inner NC picking an
     * estimator is a KERNEL selection, not a second transformation, and is
     * correctly not cut.
     */
    private void runTransformAnalyzer() {
        long T0 = System.nanoTime();
        String token = jline.solvers.tr.TransformSolve.requested(options);
        SolverOptions sub = options.copy();
        sub.config.put("transform", token);
        jline.solvers.tr.TransformSolve.Result tr = jline.solvers.tr.TransformSolve.run(
                this.model, this.sn, sub, new jline.solvers.tr.TransformSolve.InnerSolve() {
                    @Override
                    public jline.solvers.tr.TransformSolve.Inner solve(jline.lang.Network submodel,
                                                                       SolverOptions opts) {
                        SolverNC inner = new SolverNC(submodel, opts);
                        try {
                            inner.runAnalyzer();
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                        double lG = (inner.result instanceof NCResult)
                                ? ((NCResult) inner.result).logNormConstAggr() : Double.NaN;
                        return new jline.solvers.tr.TransformSolve.Inner(
                                inner.getAvgQLen(), inner.getAvgUtil(), inner.getAvgRespT(),
                                inner.getAvgTput(), inner.getAvgSysTput(), lG,
                                inner.result == null ? "" : inner.result.method);
                    }
                });

        double runtime = (System.nanoTime() - T0) / 1000000000.0;
        String reported = options.method + "/" + tr.method;
        jline.solvers.AvgHandle TH = getAvgTputHandles();
        Matrix AN = jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput(sn, tr.T, TH);
        this.setAvgResults(tr.Q, tr.U, tr.R, tr.T, AN, new Matrix(0, 0),
                tr.C, tr.X, runtime, reported, tr.iter);
    }
}
