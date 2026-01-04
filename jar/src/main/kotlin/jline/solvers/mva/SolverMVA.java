/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.io.Ret.ProbabilityResult;
import jline.solvers.mva.handlers.MVARunner;
import jline.lang.state.ToMarginal;
import jline.lang.nodes.StatefulNode;
import jline.util.matrix.Matrix;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

import static jline.util.Maths.logBinomial;
import static jline.util.Maths.factln;
import static jline.api.sn.SnIsOpenModelKt.snIsOpenModel;
import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * SolverMVA implements Mean Value Analysis (MVA) for queueing networks.
 * MVA is an exact analytical method for computing performance measures
 * of closed queueing networks, particularly effective for networks with
 * product-form solutions.
 */
public class SolverMVA extends NetworkSolver {

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
                "CacheClassSwitcher", "Cache",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_DPS", "SchedStrategy_FCFS", "SchedStrategy_SIRO", "SchedStrategy_HOL",
                "SchedStrategy_LCFS", "SchedStrategy_LCFSPR", "SchedStrategy_POLLING",
                "Fork", "Forker", "Join", "Joiner",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_LRU",
                "ClosedClass", "SelfLoopingClass", "OpenClass", "Replayer",
                "LoadDependence"
        });
        return featSupported;
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
        if (this.sn == null)
            this.sn = this.model.getStruct(false);
        if (this.options == null)
            this.options = new SolverOptions(SolverType.MVA);

        String origMethod = options.method;
        MVARunner runner = new MVARunner(this.model, this.options, this.enableChecks);
        MVAResult ret = (MVAResult) runner.runAnalyzer(this.avgHandles);
        String resultMethod = ret.method;
        if (origMethod.equals("default") && !resultMethod.equals("default") && !resultMethod.startsWith("default/")) {
            resultMethod = "default/" + resultMethod;
        }
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
        return FeatureSet.supports(featSupported, featUsed);
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
                "default", "mva", "exact", "amva", "qna",
                "qdlin", "amva.qdlin",
                "bs", "amva.bs",
                "sqni",
                "qd", "amva.qd",
                "qli", "amva.qli",
                "fli", "amva.fli",
                "lin", "egflin", "gflin", "amva.lin",
                "schmidt", "ab"
        ));

        // Add bounds for single-class closed models
        if (!snIsOpenModel(this.sn) && this.sn.nclasses == 1) {
            allMethods.addAll(Arrays.asList(
                    "aba.upper", "aba.lower", "bjb.upper", "bjb.lower",
                    "gb.upper", "gb.lower", "pb.upper", "pb.lower",
                    "sb.upper", "sb.lower"
            ));
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
