/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.GlobalConstants;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.io.Ret.ProbabilityResult;
import jline.solvers.mva.handlers.MVARunner;
import jline.solvers.mva.analyzers.Solver_mva_analyzer;
import jline.lang.state.ToMarginal;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.nodes.StatefulNode;
import jline.util.matrix.Matrix;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

import static jline.util.Maths.logBinomial;
import static jline.util.Maths.factln;
import static jline.api.sn.SnIsOpenModel.snIsOpenModel;
import static jline.api.sn.SnHasClassdepRouting.snHasClassdepRouting;
import static jline.io.InputOutput.line_debug;
import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

/**
 * SolverMVA implements Mean Value Analysis (MVA) for queueing networks.
 * MVA is an exact analytical method for computing performance measures
 * of closed queueing networks, particularly effective for networks with
 * product-form solutions.
 */
public class SolverMVA extends NetworkSolver {

    /**
     * MVA differentiates its own recursion: getSensitivityTable uses the analytic
     * branch (pfqn_sens) wherever the model is in scope.
     *
     * @return true
     */
    @Override
    public boolean supportsExactSensitivity() {
        return true;
    }

    /**
     * Cached MMT fork-join transformation, held here rather than on the per-call
     * MVARunner so that it survives across the outer iterations of SolverLN (which
     * keeps one solver per layer alive and re-solves it once per iteration). This
     * avoids re-running the serialisation deep copy of the layer every iteration.
     * <p>
     * The cache carries the NetworkStruct it was derived from; MVARunner reuses it
     * only while that struct is still the model's current one, and re-feeds its
     * rates from the base model before every reuse. It is deliberately not cleared
     * by reset(): correctness lives at the point of use, and clearing it there
     * would discard the transformation on every outer iteration.
     */
    private jline.io.Ret.FJApprox mmtCache = null;

    /**
     * Auxiliary-class arrival rates of the fork-join (MMT) fixed point, retained
     * across runAnalyzer calls so that an outer iteration (e.g. SolverLN) restarts
     * the MMT loop from the previous converged point instead of from
     * GlobalConstants.FineTol. Like options.init_sol it survives reset(); use
     * resetForkWarmStart to discard it.
     */
    private Matrix fjForkLambda = null;

    /**
     * Discards the retained MMT fixed point. Mirrors options.init_sol: the iterate
     * survives reset() and is invalidated explicitly by the caller when the chain
     * basis changes.
     */
    public void resetForkWarmStart() {
        this.fjForkLambda = null;
    }

    /**
     * Returns the retained MMT fixed-point iterate, or null if none is held.
     * Callers that drive this solver across an outer iteration without keeping the
     * instance alive (e.g. the MATLAB lang='java' dispatch, which rebuilds the
     * JLINE solver on every call) can carry it themselves via this accessor and
     * setForkWarmStart.
     */
    public Matrix getForkWarmStart() {
        return this.fjForkLambda;
    }

    /**
     * Seeds the MMT fixed point with a previously retained iterate. Conformance
     * with the current fork topology is checked at the point of use in MVARunner,
     * so a stale iterate is ignored rather than misapplied.
     */
    public void setForkWarmStart(Matrix fjForkLambda) {
        this.fjForkLambda = fjForkLambda;
    }

    /**
     * Creates a new SolverMVA instance with a specific method.
     *
     * @param model The network model to analyze
     * @param method The MVA method to use
     */
    public SolverMVA(Network model, String method) {
        super(model, "SolverMVA", SolverMVA.defaultOptions().method(method));
        this.sn = model.getStruct(false);
        this.result = new MVAResult();
    }

    /**
     * Creates a new SolverMVA instance with specific options.
     *
     * @param model The network model to analyze
     * @param options The solver options to use
     */
    public SolverMVA(Network model, SolverOptions options) {
        super(model, "SolverMVA", options);
        this.sn = model.getStruct(false);
        this.result = new MVAResult();
    }

    /**
     * Creates a new SolverMVA instance with variable arguments for options.
     *
     * @param model The network model to analyze
     * @param varargin Variable arguments for solver options
     */
    public SolverMVA(Network model, Object... varargin) {
        this(model, SolverMVA.defaultOptions());
        this.options = Solver.parseOptions(this.options, varargin);
    }

    /**
     * Creates a new SolverMVA instance with default options.
     *
     * @param model The network model to analyze
     */
    public SolverMVA(Network model) {
        super(model, "SolverMVA", SolverMVA.defaultOptions());
        this.sn = model.getStruct(false);
        this.result = new MVAResult();
    }

    /**
     * Returns the default solver options for the MVA solver.
     *
     * @return Default solver options with SolverType.MVA
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.MVA);
    }

    /**
     * Returns the feature set supported by the MVA solver
     *
     * @return - the feature set supported by the MVA solver
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink", "Source",
                "ClassSwitch", "Delay", "DelayStation", "Queue",
                "APH", "Coxian", "Erlang", "Exp", "HyperExp", "BMAP",
                "Pareto", "Weibull", "Lognormal", "Uniform", "Det",
                "StatelessClassSwitcher", "InfiniteServer", "SharedServer", "Buffer", "Dispatcher",
                "CacheClassSwitcher", "Cache", "CacheRetrieval",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_DPS", "SchedStrategy_FCFS", "SchedStrategy_SIRO", "SchedStrategy_HOL",
                "SchedStrategy_LCFS", "SchedStrategy_LCFSPR", "SchedStrategy_POLLING",
                "SchedStrategy_OI", "SchedStrategy_PAS",
                "Fork", "Forker", "Join", "Joiner",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_LRU",
                "ReplacementStrategy_HLRU",
                "MMAP", // marked MAP sources (cache LRU via cache_ttl_lrum_map)
                "ClosedClass", "SelfLoopingClass", "OpenClass", "Replayer",
                "LoadDependence", "ClassDependence", "JointDependence"
        });
        return featSupported;
    }

    /**
     * Per-method feature deltas applied to the base MVA envelope. QNA is a
     * two-moment open-network method, so it drops closed-class support. The
     * queueing-system and bounds methods are already structurally restricted by
     * listValidMethods and inherit the base envelope. RQNA (open only) natively
     * consumes the non-renewal MAP/MMPP/MMAP/RAP point processes through the
     * index-of-dispersion equations, so it declares them on top of the base
     * envelope; it dispatches from Solver_mva_analyzer to Solver_rqna, and the
     * default method also selects it for bursty single-class open networks.
     * Mirrors the MATLAB/Python SolverMVA.getMethodFeatureSet.
     *
     * @param method the concrete method name
     * @return the per-method FeatureSet
     */
    public static FeatureSet methodFeatureSet(String method) {
        FeatureSet featSupported = SolverMVA.getFeatureSet();
        if ("qna".equals(method)) {
            featSupported.setFalse(new String[]{"ClosedClass", "SelfLoopingClass"});
        } else if ("rqna".equals(method)) {
            featSupported.setTrue(new String[]{"MAP", "MMPP2", "MMAP", "RAP"});
            featSupported.setFalse(new String[]{"ClosedClass", "SelfLoopingClass"});
        }
        return featSupported;
    }

    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        return SolverMVA.methodFeatureSet(method);
    }

    /**
     * Returns the network structure used by this solver.
     * If the structure is not yet initialized, it will be created from the model.
     *
     * @return The network structure
     */
    public NetworkStruct getStruct() {
        if (this.sn == null)
            this.sn = this.model.getStruct(false);
        return this.sn;
    }

    /**
     * Sets the network structure for this solver.
     *
     * @param sn The network structure to use
     */
    public void setStruct(NetworkStruct sn) {
        this.sn = sn;
    }

    /**
     * Runs the MVA analyzer to solve the queueing network.
     * This method executes the Mean Value Analysis algorithm and stores
     * the results in the solver's result object.
     *
     * @throws IllegalAccessException if there are access issues during analysis
     */
    @Override
    public void runAnalyzer() throws IllegalAccessException {
        if (this.options == null)
            this.options = new SolverOptions(SolverType.MVA);
        // Propagate solver verbose level to global
        GlobalConstants.Verbose = options.verbose;
        if (this.sn == null)
            this.sn = this.model.getStruct(false);

        // see _kb/06-solver-catalog.md for rationale
        if (sn.nregions > 0) {
            throw new RuntimeException("This model uses a Finite Capacity Region (addRegion), "
                    + "which is not supported by SolverMVA (the region's aggregate job limit "
                    + "would be silently ignored). Use SolverJMT, or setCapacity for a "
                    + "single-station limit.");
        }

        // see _kb/06-solver-catalog.md for rationale
        String capacityReason = finiteCapacityReason(this.model, sn);
        if (capacityReason != null) {
            throw new RuntimeException(capacityReason);
        }

        if (sn.immfeed != null && sn.immfeed.elementSum() > 0) {
            line_warning(mfilename(new Object(){}), "SolverMVA does not handle immediate feedback (immfeed); the solver will treat self-loops as class-switching with re-queueing.");
        }

        line_debug(options.verbose, String.format("MVA solver starting: method=%s, nstations=%d, nclasses=%d",
            options.method, sn.nstations, sn.nclasses));

        String origMethod = options.method;
        // see _kb/06-solver-catalog.md for rationale
        MVARunner runner = new MVARunner(this.model, this.options, this.enableChecks, this.mmtCache,
                this.fjForkLambda);
        MVAResult ret = (MVAResult) runner.runAnalyzer(this.avgHandles);
        this.mmtCache = runner.getMmtCache();
        this.fjForkLambda = runner.getFjForkLambda();
        String resultMethod = ret.method;
        if (origMethod.equals("default") && !resultMethod.equals("default") && !resultMethod.startsWith("default/")) {
            resultMethod = "default/" + resultMethod;
        }
        line_debug(options.verbose, String.format("MVA solver completed: method=%s, iter=%d", resultMethod, ret.iter));
        this.setAvgResults(ret.QN, ret.UN, ret.RN, ret.TN, ret.AN, ret.WN, ret.CN, ret.XN, ret.runtime, resultMethod, ret.iter);
    }

    /**
     * Checks whether the given model is supported by the MVA solver
     *
     * @param model - the network model
     * @return - true if the model is supported, false otherwise
     */
    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverMVA.getFeatureSet();
        if (!FeatureSet.supports(featSupported, featUsed)) {
            return false;
        }
        // Registry inclusion cannot see finite capacity (there is no feature
        // name for it), so apply the structural gate as well.
        String reason = finiteCapacityReason(model, model.getStruct(false));
        if (reason != null) {
            line_warning(mfilename(new Object() {}), reason);
            return false;
        }
        return true;
    }

    /**
     * MVA-specific finite-capacity gate: Blocking-After-Service models are
     * exempt because MVA offers the Smith queue-decomposition method "sqd",
     * and Solver_mva_analyzer routes a BAS model to Solver_sqd under the
     * default method too, so the finite buffers ARE honoured on every MVA
     * path. Everything else defers to the shared product-form gate.
     *
     * @param model - the network model
     * @param sn    - the network structure
     * @return null if the model has no binding capacity, otherwise the reason
     */
    public static String finiteCapacityReason(Network model, NetworkStruct sn) {
        if (Solver_mva_analyzer.isBasModel(sn)) {
            return null;
        }
        return NetworkSolver.bindingCapacityReason(model, sn, "SolverMVA");
    }

    /**
     * Returns the logarithm of the normalizing constant for the aggregate state
     * probabilities
     *
     * @return the log of the normalizing constant
     */
    public ProbabilityResult getProbNormConstAggr() {
        if (this.result == null) {
            try {
                this.runAnalyzer();
            } catch (IllegalAccessException e) {
                throw new RuntimeException("Failed to run analyzer", e);
            }
        }

        MVAResult mvaResult = (MVAResult) this.result;
        if (mvaResult != null && !Double.isNaN(mvaResult.logNormConstAggr)) {
            return new ProbabilityResult(mvaResult.logNormConstAggr, true);
        } else {
            // If current method doesn't support exact calculations, try running with exact method
            SolverOptions exactOptions = this.options.copy();
            exactOptions.method = "exact";

            SolverMVA exactSolver = new SolverMVA(this.model, exactOptions);
            try {
                exactSolver.runAnalyzer();
                MVAResult exactResult = (MVAResult) exactSolver.result;
                if (exactResult != null && !Double.isNaN(exactResult.logNormConstAggr)) {
                    return new ProbabilityResult(exactResult.logNormConstAggr, true);
                }
            } catch (Exception e) {
                throw new RuntimeException("getProbNormConstAggr: exact method not supported for this model", e);
            }

            throw new RuntimeException("getProbNormConstAggr: solver must be run first or exact method not supported");
        }
    }

    /**
     * Get marginal state probabilities for a specific station
     *
     * @param ist station index (0-based)
     * @return ProbabilityResult with probability and log probability
     */
    public ProbabilityResult getProbAggr(int ist) {
        if (ist >= this.sn.nstations) {
            throw new RuntimeException("Station number exceeds the number of stations in the model.");
        }

        if (this.result == null) {
            try {
                this.runAnalyzer();
            } catch (IllegalAccessException e) {
                throw new RuntimeException("Failed to run analyzer", e);
            }
        }

        Matrix Q = this.getAvgQLen();
        Matrix N = this.sn.njobs;

        if (N.isFinite()) {
            switch (this.options.method) {
                case "exact":
                    throw new RuntimeException("Exact marginal state probabilities not available yet in SolverMVA.");
                default:
                    // Get state for this station  
                    Map<StatefulNode, Matrix> stateMap = this.sn.state;
                    int statefulIndex = (int) this.sn.stationToStateful.get(ist);

                    // Find the corresponding StatefulNode and get its state
                    Matrix state = null;
                    int currentIndex = 0;
                    for (Map.Entry<StatefulNode, Matrix> entry : stateMap.entrySet()) {
                        if (currentIndex == statefulIndex) {
                            state = entry.getValue();
                            break;
                        }
                        currentIndex++;
                    }

                    if (state == null) {
                        throw new RuntimeException("Could not find state for station " + ist);
                    }

                    jline.lang.state.State.StateMarginalStatistics margStats = ToMarginal.toMarginal(this.sn, ist, state, null, null, null, null, null);
                    Matrix nir = margStats.nir;

                    // Binomial approximation with mean fitted to queue-lengths
                    // Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997
                    double logPnir = 0.0;
                    for (int r = 0; r < nir.getNumCols(); r++) {
                        int nirVal = (int) nir.get(0, r);
                        int NVal = (int) N.get(r);
                        double QVal = Q.get(ist, r);

                        logPnir += logBinomial(NVal, nirVal);
                        if (QVal > 0 && NVal > 0) {
                            logPnir += nirVal * Math.log(QVal / NVal);
                            logPnir += (NVal - nirVal) * Math.log(1.0 - QVal / NVal);
                        }
                    }
                    double Pnir = Math.exp(logPnir);
                    ProbabilityResult result = new ProbabilityResult(Pnir);
                    result.logNormalizingConstant = logPnir;
                    return result;
            }
        } else {
            throw new RuntimeException("getProbAggr not yet implemented for models with open classes.");
        }
    }

    /**
     * Get joint system state probabilities
     *
     * @return ProbabilityResult with probability and log probability
     */
    public ProbabilityResult getProbSysAggr() {
        line_debug(options.verbose, "MVA: computing system aggregate probabilities");
        if (this.result == null) {
            try {
                this.runAnalyzer();
            } catch (IllegalAccessException e) {
                throw new RuntimeException("Failed to run analyzer", e);
            }
        }

        Matrix Q = this.getAvgQLen();
        Matrix N = this.sn.njobs;

        if (N.isFinite()) {
            switch (this.options.method) {
                case "exact":
                    throw new RuntimeException("Exact joint state probabilities not available yet in SolverMVA.");
                default:
                    Map<StatefulNode, Matrix> stateMap = this.sn.state;

                    // Binomial approximation with mean fitted to queue-lengths
                    // Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997
                    double logPn = 0.0;
                    for (int r = 0; r < N.length(); r++) {
                        logPn += factln((int) N.get(r));
                    }

                    for (int ist = 0; ist < this.sn.nstations; ist++) {
                        int statefulIndex = (int) this.sn.stationToStateful.get(ist);

                        // Find the corresponding StatefulNode and get its state
                        Matrix stateIst = null;
                        int currentIndex = 0;
                        for (Map.Entry<StatefulNode, Matrix> entry : stateMap.entrySet()) {
                            if (currentIndex == statefulIndex) {
                                stateIst = entry.getValue();
                                break;
                            }
                            currentIndex++;
                        }

                        if (stateIst == null) {
                            continue; // Skip if no state found
                        }

                        jline.lang.state.State.StateMarginalStatistics margStats = ToMarginal.toMarginal(this.sn, ist, stateIst, null, null, null, null, null);
                        Matrix nir = margStats.nir;

                        for (int r = 0; r < this.sn.nclasses; r++) {
                            int nirVal = (int) nir.get(0, r);
                            logPn -= factln(nirVal);
                            double QVal = Q.get(ist, r);
                            int NVal = (int) N.get(r);
                            if (QVal > 0 && NVal > 0) {
                                logPn += nirVal * Math.log(QVal / NVal);
                            }
                        }
                    }

                    double Pnir = Math.exp(logPn);
                    ProbabilityResult result = new ProbabilityResult(Pnir);
                    result.logNormalizingConstant = logPn;
                    return result;
            }
        } else {
            throw new RuntimeException("getProbSysAggr not yet implemented for models with open classes.");
        }
    }

    /**
     * Get marginalized state probabilities for a specific station and job class
     *
     * @param ist station index (0-based)
     * @param jobclass job class index (0-based)
     * @return ProbabilityResult with marginalized state probabilities
     */
    public ProbabilityResult getProbMarg(int ist, int jobclass) {
        if (ist >= this.sn.nstations) {
            throw new RuntimeException("Station number exceeds the number of stations in the model.");
        }
        if (jobclass >= this.sn.nclasses) {
            throw new RuntimeException("Job class index exceeds the number of classes in the model.");
        }

        if (this.result == null) {
            try {
                this.runAnalyzer();
            } catch (IllegalAccessException e) {
                throw new RuntimeException("Failed to run analyzer", e);
            }
        }

        Matrix N = this.sn.njobs;

        if (N.isFinite()) {
            switch (this.options.method) {
                case "exact":
                    throw new RuntimeException("Exact marginalized state probabilities not available yet in SolverMVA.");
                default:
                    // Use binomial approximation for marginalized probabilities 
                    // Similar to getProbAggr but for a single class
                    Matrix Q = this.getAvgQLen();
                    double qVal = Q.get(ist, jobclass);
                    int nVal = (int) N.get(jobclass);

                    // Create probability vector for this class at this station
                    Matrix Pmarg = new Matrix(1, nVal + 1);
                    for (int k = 0; k <= nVal; k++) {
                        // Binomial probability with mean qVal
                        // Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997
                        double logPk = logBinomial(nVal, k) + k * Math.log(qVal / nVal) + (nVal - k) * Math.log(1 - qVal / nVal);
                        Pmarg.set(0, k, Math.exp(logPk));
                    }

                    ProbabilityResult result = new ProbabilityResult(Pmarg);
                    return result;
            }
        } else {
            throw new RuntimeException("getProbMarg not yet implemented for models with open classes.");
        }
    }

    /**
     * Get marginalized state probabilities for a specific station and job class with state filter
     *
     * @param ist station index (0-based)
     * @param jobclass job class index (0-based)
     * @param state_m marginalized state vector to query (optional, null for all states)
     * @return ProbabilityResult with marginalized state probabilities
     */
    public ProbabilityResult getProbMarg(int ist, int jobclass, Matrix state_m) {
        ProbabilityResult result = getProbMarg(ist, jobclass);

        if (state_m != null && !state_m.isEmpty()) {
            // Filter results based on state_m
            Matrix filteredProb = new Matrix(1, state_m.length());
            for (int i = 0; i < state_m.length(); i++) {
                int stateIndex = (int) state_m.get(i);
                if (stateIndex < result.probability.getNumCols()) {
                    filteredProb.set(0, i, result.probability.get(0, stateIndex));
                }
            }
            return new ProbabilityResult(filteredProb);
        }

        return result;
    }

    /**
     * List all valid solution methods supported by this solver
     *
     * @return array of valid method names
     */
    public String[] listValidMethods() {
        // Base set of methods
        List<String> allMethods = new ArrayList<String>(Arrays.asList(
                "default", "mva", "exact", "amva",
                "sum", "esum",
                "qdlin", "amva.qdlin",
                "bs", "amva.bs",
                "sqni",
                "qd", "amva.qd",
                "qli", "amva.qli",
                "fli", "amva.fli",
                "lin", "egflin", "gflin", "amva.lin",
                "schmidt", "ab", "sqd"
        ));

        // see _kb/06-solver-catalog.md for rationale
        if (snIsOpenModel(this.sn)) {
            allMethods.addAll(4, Arrays.asList("qna", "rqna"));
        }

        // MVAC (exact mean value analysis by chain, pfqn_mvac): closed
        // single-server product-form networks only; rejects open/mixed and
        // multiserver at solve time.
        boolean anyOpen = false;
        for (int i = 0; i < this.sn.njobs.getNumCols(); i++) {
            if (Double.isInfinite(this.sn.njobs.get(i))) { anyOpen = true; break; }
        }
        if (!anyOpen) {
            allMethods.add("mvac");
        }

        // see _kb/06-solver-catalog.md for rationale
        if (!snIsOpenModel(this.sn) && !snHasClassdepRouting(this.sn)) {
            allMethods.addAll(Arrays.asList("marie", "amva.marie"));
        }

        // Add queueing system methods for open single-class 2-station models
        if (snIsOpenModel(this.sn) && this.sn.nstations == 2 && this.sn.nclasses == 1) {
            allMethods.addAll(Arrays.asList(
                    "mm1", "mmk", "mg1", "mgi1", "gm1", "gig1", "gim1", "gig1.kingman",
                    "gigk", "gigk.kingman_approx",
                    "gig1.gelenbe", "gig1.heyman", "gig1.kimura", "gig1.allen",
                    "gig1.kobayashi", "gig1.klb", "gig1.marchal"
            ));
        }

        return allMethods.toArray(new String[0]);
    }

    /**
     * Check if the network has multi-server stations
     *
     * @return true if any station has more than one server
     */
    private boolean hasMultiServerStations() {
        if (this.sn == null || this.sn.nservers == null) {
            return false;
        }
        
        for (int ist = 0; ist < this.sn.nstations; ist++) {
            if (this.sn.nservers.get(ist) > 1) {
                return true;
            }
        }
        
        return false;
    }
}
