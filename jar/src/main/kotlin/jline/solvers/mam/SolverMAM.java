/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SolverType;
import jline.solvers.*;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.DistributionResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static jline.api.sn.SnGetArvRFromTputKt.snGetArvRFromTput;
import static jline.api.sn.SnGetResidTFromRespTKt.snGetResidTFromRespT;
import static jline.io.InputOutputKt.*;
import static jline.solvers.mam.analyzers.Solver_mam_analyzerKt.solver_mam_analyzer;
import static jline.solvers.mam.handlers.Solver_mam_passage_timeKt.solver_mam_passage_time;
import static jline.api.fj.FJValidationKt.isFJ;
import static jline.solvers.mam.SolverMAMFJKt.solverMAMFJ;
import static jline.solvers.mam.SolverMAMFJKt.interpolatePercentile;

import jline.api.fj.FJInfo;
import jline.lib.fjcodes.FJPercentileResult;
import jline.lib.butools.MMAPPH1FCFSKt;
import jline.lang.constant.NodeType;
import kotlin.Pair;
import java.util.HashMap;


/**
 * Solver for Matrix Analytic Methods (MAM) applied to queueing networks.
 * 
 * <p>SolverMAM implements matrix-analytic techniques for analyzing queueing networks
 * with Markovian arrival processes (MAP), phase-type service distributions, and
 * other non-exponential characteristics that go beyond product-form assumptions.</p>
 * 
 * <p>Key MAM solver capabilities:
 * <ul>
 *   <li>Markovian Arrival Process (MAP) modeling</li>
 *   <li>Phase-type (PH) service distribution analysis</li>
 *   <li>Matrix-geometric solution methods</li>
 *   <li>Quasi-Birth-Death (QBD) process analysis</li>
 *   <li>Non-product-form queueing network solutions</li>
 *   <li>Passage time distribution computation</li>
 * </ul>
 * </p>
 * 
 * <p>This solver is particularly useful for networks with correlated arrivals,
 * general service times, and complex dependency structures that cannot be
 * analyzed using traditional product-form methods.</p>
 * 
 * @see jline.api.mam
 * @see MAMResult
 * @see MAMOptions
 * @since 1.0
 */
public class SolverMAM extends NetworkSolver {

    private List<FJPercentileResult> percentileResults = null;

    public SolverMAM(Network model) {
        super(model, "SolverMAM", new SolverOptions(SolverType.MAM));
        this.result = new MAMResult();
    }

    public SolverMAM(Network model, String method) {
        super(model, "SolverMAM", Solver.defaultOptions().method(method));
        this.result = new MAMResult();
    }

    public SolverMAM(Network model, Object... varargin) {
        this(model, Solver.defaultOptions());
        Solver.parseOptions(this.options, varargin);
        this.result = new MAMResult();
    }


    public SolverMAM(Network model, SolverOptions options) {
        super(model, "SolverMAM", options);
        this.result = new MAMResult();
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.MAM);
    }

    /**
     * Returns the feature set supported by the MAM solver
     *
     * @return - the feature set supported by the MAM solver
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink", "Source",
                "Delay", "DelayStation", "Queue",
                "APH", "Coxian", "Erlang", "Exp", "HyperExp", "MMPP2", "MAP", "ME", "RAP",
                "Det", "Gamma", "Lognormal", "Pareto", "Uniform", "Weibull",
                "StatelessClassSwitcher", "InfiniteServer",
                "ClassSwitch",
                "SharedServer", "Buffer", "Dispatcher",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS", "SchedStrategy_HOL",
                "SchedStrategy_FCFS",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ClosedClass", "SelfLoopingClass",
                "OpenClass"
        });
        return featSupported;
    }

    public NetworkStruct getStruct() {
        return model.getStruct(true);
    }

    public List<String> listValidMethods(Network model) {
        return Arrays.asList("default", "dec.source", "dec.mmap", "dec.poisson", "mna", "inap", "inapplus", "exact", "ldqbd");
    }

    public List<String> listValidMethods() {
        return listValidMethods(null);
    }

    @Override
    public void runAnalyzer() {
        double start = System.nanoTime();
        NetworkStruct sn = getStruct();

        // Check if network has Fork-Join topology
        Pair<Boolean, FJInfo> fjCheck = isFJ(sn);
        if (fjCheck.getFirst()) {
            line_debug(options.verbose, String.format("Detected Fork-Join topology with K=%d parallel queues", fjCheck.getSecond().getK()));
            line_debug(options.verbose, "Using FJ_codes algorithm for percentile analysis");

            try {
                // Use FJ solver
                jline.solvers.mam.MAMFJResult fjResult = solverMAMFJ(
                    sn,
                    new double[]{0.50, 0.90, 0.95, 0.99},  // Default percentiles
                    100,     // C parameter
                    "NARE"   // T-matrix method
                );

                // Store percentile results
                this.percentileResults = fjResult.getPercentileResults();

                // Convert to standard MAM result
                double finish = System.nanoTime();
                SolverResult res = new SolverResult();
                res.QN = fjResult.getQN();
                res.UN = fjResult.getUN();
                res.RN = fjResult.getRN();
                res.TN = fjResult.getTN();
                res.XN = fjResult.getXN();
                res.runtime = (finish - start) / 1000000000.0;
                res.method = "fj/NARE";

                line_debug(options.verbose, String.format("FJ_codes solution completed in %.3f seconds", res.runtime));
                this.result = res;
                return;

            } catch (Exception e) {
                line_warning(mfilename(new Object() {}), "FJ_codes solver failed: " + e.getMessage());
                line_warning(mfilename(new Object() {}), "Falling back to standard MAM solver");
                // Fall through to standard MAM solver
            }
        }

        line_debug(options.verbose, String.format("MAM solver starting: method=%s, nstations=%d, nclasses=%d",
                options.method, sn.nstations, sn.nclasses));
        if (true) { // ~snHasMultipleClosedClasses(sn)
            line_debug(options.verbose, "Running MAM analysis, calling solver_mam_analyzer");
            MAMResult res = solver_mam_analyzer(sn, options);

            AvgHandle T = getAvgTputHandles();
            Matrix AN = snGetArvRFromTput(sn, res.TN, T);
            AvgHandle W = getAvgResidTHandles();
            Matrix WN = snGetResidTFromRespT(sn, res.RN, W);

            double finish = System.nanoTime();
            res.runtime = (finish - start) / 1000000000.0;
            String methodName;
            if (options.method.equals("default")) {
                methodName = "default/" + res.method;
            } else {
                methodName = options.method;
            }
            this.setAvgResults(res.QN, res.UN, res.RN, res.TN, AN, WN, res.CN, res.XN, res.runtime, methodName, res.iter);
        } else {
            line_warning(mfilename(new Object() {
            }), "SolverMAM supports at most a single closed class.");
        }
    }

    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverMAM.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * Returns cumulative distribution functions of response times at steady-state.
     * This method computes response time distributions using matrix-analytic methods.
     *
     * @param R response time handles (optional)
     * @return result containing CDFs for response times
     */
    @Override
    public DistributionResult getCdfRespT(AvgHandle R) {
        long startTime = System.nanoTime();

        // Use default handles if not provided
        if (R == null) {
            R = getAvgRespTHandles();
        }

        NetworkStruct sn = getStruct();

        // Get steady-state solution first
        getAvg();

        // Call solver_mam_passage_time with sn.proc, matching MATLAB implementation
        // MATLAB: RD = solver_mam_passage_time(sn, sn.proc, options)
        Map<Integer, MatrixCell> RD = solver_mam_passage_time(sn, sn.proc, options);

        double runtime = (System.nanoTime() - startTime) / 1000000000.0;

        // Convert RD to Matrix for setDistribResults (Java requirement)
        Matrix rdMatrix = new Matrix(sn.nstations, sn.nclasses);
        setDistribResults(rdMatrix, runtime);

        // Create DistributionResult for return value
        DistributionResult result = new DistributionResult(sn.nstations, sn.nclasses, "response_time");

        // Convert RD to DistributionResult format
        for (Integer stationIdx : RD.keySet()) {
            MatrixCell stationDistrs = RD.get(stationIdx);
            if (stationDistrs != null) {
                for (int k = 0; k < sn.nclasses; k++) {
                    Matrix cdfData = stationDistrs.get(k);
                    if (cdfData != null && !cdfData.isEmpty()) {
                        result.setCdf(stationIdx, k, cdfData);
                    }
                }
            }
        }

        return result;
    }

    /**
     * Returns cumulative distribution functions of response times at steady-state.
     * Uses default response time handles.
     *
     * @return result containing CDFs for response times
     */
    @Override
    public DistributionResult getCdfRespT() {
        return getCdfRespT(null);
    }

    /**
     * Returns cumulative distribution functions of passage times at steady-state.
     *
     * @param R response time handles (optional)
     * @return result containing CDFs for passage times
     */
    @Override
    public DistributionResult getCdfPassT(AvgHandle R) {
        NetworkStruct sn = getStruct();

        // Call solver_mam_passage_time with sn.proc, same as getCdfRespT
        Map<Integer, MatrixCell> RD = solver_mam_passage_time(sn, sn.proc, options);

        // Create DistributionResult for passage times
        DistributionResult result = new DistributionResult(sn.nstations, sn.nclasses, "passage_time");

        // Convert RD to DistributionResult format
        for (Integer stationIdx : RD.keySet()) {
            MatrixCell stationDistrs = RD.get(stationIdx);
            if (stationDistrs != null) {
                for (int k = 0; k < sn.nclasses; k++) {
                    Matrix cdfData = stationDistrs.get(k);
                    if (cdfData != null && !cdfData.isEmpty()) {
                        result.setCdf(stationIdx, k, cdfData);
                    }
                }
            }
        }

        return result;
    }

    /**
     * Returns cumulative distribution functions of passage times at steady-state.
     * Uses default response time handles.
     *
     * @return result containing CDFs for passage times
     */
    @Override
    public DistributionResult getCdfPassT() {
        return getCdfPassT(null);
    }

    /**
     * Returns cumulative distribution functions of passage times during transient analysis.
     *
     * @param R response time handles (optional)
     * @return result containing transient CDFs for passage times
     */
    @Override
    public DistributionResult getTranCdfPassT(AvgHandle R) {
        // MAM solver currently does not support transient passage time analysis
        // Return empty result
        NetworkStruct sn = getStruct();
        DistributionResult result = new DistributionResult(sn.nstations, sn.nclasses, "passage_time");
        return result;
    }

    /**
     * Returns cumulative distribution functions of passage times during transient analysis.
     * Uses default response time handles.
     *
     * @return result containing transient CDFs for passage times
     */
    @Override
    public DistributionResult getTranCdfPassT() {
        return getTranCdfPassT(null);
    }

    /**
     * Get marginal queue-length probability distribution for a job class.
     *
     * <p>Computes the probability distribution P(n jobs of class r) for n=0,1,...,N(r)
     * using MMAPPH1FCFS from BUTools.</p>
     *
     * <p><strong>Current limitations:</strong></p>
     * <ul>
     *   <li>Only supported for single queue models</li>
     *   <li>Requires Queue station with FCFS scheduling</li>
     * </ul>
     *
     * @param node Station/node index (0-based)
     * @param jobclass Job class index (0-based)
     * @param state_m Optional state levels to query (null for all states, 0-based indexing)
     * @return Probability result containing marginal probabilities
     * @throws IllegalArgumentException if station or class index is invalid
     * @throws UnsupportedOperationException if model structure is not supported
     */
    @Override
    public ProbabilityResult getProbMarg(int node, int jobclass, Matrix state_m) {
        NetworkStruct sn = getStruct();

        // Validate station index
        if (node >= sn.nstations) {
            throw new IllegalArgumentException("Station number exceeds the number of stations in the model.");
        }
        if (jobclass >= sn.nclasses) {
            throw new IllegalArgumentException("Job class index exceeds the number of classes in the model.");
        }

        // Check if this is a network (more than one queue station)
        int queueStations = 0;
        for (int i = 0; i < sn.nstations; i++) {
            int nodeIdx = (int) sn.stationToNode.get(i, 0);
            if (sn.nodetype.get(nodeIdx) == NodeType.Queue) {
                queueStations++;
            }
        }

        if (queueStations > 1) {
            throw new UnsupportedOperationException(
                "getProbMarg is currently only supported for single queue models in SolverMAM. " +
                "For networks, this method is not yet implemented.");
        }

        if (queueStations == 0) {
            throw new UnsupportedOperationException("Model does not contain any queue stations.");
        }

        // Ensure results are available
        if (result == null || ((SolverResult) result).QN == null) {
            runAnalyzer();
        }

        int K = sn.nclasses;
        Matrix N = sn.njobs.transpose();

        // Determine max queue length to compute
        int maxLevel = 100;
        if (options.cutoff != null && options.cutoff.length() > 0) {
            double cutoffVal = options.cutoff.get(0, 0);
            if (cutoffVal > 0 && Double.isFinite(cutoffVal)) {
                maxLevel = (int) cutoffVal;
            }
        }

        // For closed models, limit by population
        boolean isClosed = true;
        for (int k = 0; k < K; k++) {
            if (!Double.isFinite(N.get(k, 0))) {
                isClosed = false;
                break;
            }
        }

        if (isClosed) {
            int totalPop = 0;
            for (int k = 0; k < K; k++) {
                totalPop += (int) N.get(k, 0);
            }
            maxLevel = Math.min(maxLevel, totalPop + 1);
        }

        // Build arrival and service processes for MMAPPH1FCFS
        try {
            // Get throughput from results to approximate arrival rates
            SolverResult res = (SolverResult) result;
            double lambdaTotal = 0;
            for (int k = 0; k < K; k++) {
                lambdaTotal += res.TN.get(node, k);
            }

            if (lambdaTotal < 1e-10) {
                // No traffic - all probability at state 0
                Matrix Pmarg = new Matrix(1, maxLevel);
                Pmarg.set(0, 0, 1.0);

                ProbabilityResult probResult = new ProbabilityResult(filterByState(Pmarg, state_m));
                probResult.nodeIndex = node;
                return probResult;
            }

            // Build MMAP arrival process (simplified: use exponential arrivals based on throughput)
            MatrixCell D = new MatrixCell(K + 1);
            Matrix D0 = new Matrix(1, 1);
            D0.set(0, 0, -lambdaTotal);
            D.set(0, D0);  // D0
            for (int k = 0; k < K; k++) {
                Matrix Dk = new Matrix(1, 1);
                Dk.set(0, 0, res.TN.get(node, k));
                D.set(k + 1, Dk);
            }

            // Build service process parameters
            Map<Integer, Matrix> sigma = new HashMap<>();
            Map<Integer, Matrix> S = new HashMap<>();

            for (int k = 0; k < K; k++) {
                double rate = sn.rates.get(node, k);
                if (Double.isNaN(rate) || rate <= 0) {
                    rate = 1.0;  // Default rate
                }
                // Exponential service (1-phase PH)
                Matrix sigmaK = new Matrix(1, 1);
                sigmaK.set(0, 0, 1.0);
                sigma.put(k, sigmaK);
                Matrix Sk = new Matrix(1, 1);
                Sk.set(0, 0, -rate);
                S.put(k, Sk);
            }

            // Call MMAPPH1FCFS
            Map<String, java.util.Map<Integer, Matrix>> mmapResult = MMAPPH1FCFSKt.MMAPPH1FCFS(
                D, sigma, S,
                null,           // numOfQLMoms
                maxLevel,       // numOfQLProbs
                null,           // numOfSTMoms
                null,           // stDistr
                false,          // stDistrME
                false,          // stDistrPH
                1e-14,          // prec
                null            // classes
            );

            // Extract queue length distribution
            Matrix Pmarg;
            if (mmapResult.containsKey("ncDistr") && mmapResult.get("ncDistr").containsKey(jobclass)) {
                Pmarg = mmapResult.get("ncDistr").get(jobclass);
                // Normalize
                double sum = 0;
                for (int i = 0; i < Pmarg.length(); i++) {
                    sum += Math.abs(Pmarg.get(0, i));
                }
                if (sum > 0) {
                    for (int i = 0; i < Pmarg.length(); i++) {
                        Pmarg.set(0, i, Math.abs(Pmarg.get(0, i)) / sum);
                    }
                }
            } else {
                // Fallback: uniform distribution
                Pmarg = new Matrix(1, maxLevel);
                Pmarg.fill(1.0 / maxLevel);
            }

            ProbabilityResult probResult = new ProbabilityResult(filterByState(Pmarg, state_m));
            probResult.nodeIndex = node;
            return probResult;

        } catch (Exception e) {
            throw new RuntimeException("Failed to compute marginal probabilities: " + e.getMessage(), e);
        }
    }

    /**
     * Filter probability distribution by specific states if requested.
     */
    private Matrix filterByState(Matrix Pmarg, Matrix state_m) {
        if (state_m == null) {
            return Pmarg;
        }

        int numStates = (int) state_m.length();
        Matrix filtered = new Matrix(1, numStates);
        for (int i = 0; i < numStates; i++) {
            int stateIdx = (int) state_m.get(0, i);
            if (stateIdx < Pmarg.length()) {
                filtered.set(0, i, Pmarg.get(0, stateIdx));
            }
        }
        return filtered;
    }

    /**
     * Get marginal queue-length probability distribution for a job class (all states).
     *
     * <p><strong>NOT YET IMPLEMENTED</strong> - This method is a stub.</p>
     *
     * @param node Station/node index (0-based)
     * @param jobclass Job class index (0-based)
     * @return Probability result - currently throws UnsupportedOperationException
     * @throws UnsupportedOperationException Always thrown - method not yet implemented
     * @see #getProbMarg(int, int, Matrix)
     */
    @Override
    public ProbabilityResult getProbMarg(int node, int jobclass) {
        return getProbMarg(node, jobclass, null);
    }

    /**
     * Get response time percentiles from Fork-Join analysis
     *
     * <p>This method retrieves percentile values computed by the FJ_codes algorithm
     * for Fork-Join queueing systems. It automatically detects FJ topology and
     * computes percentiles using the algorithm from "Beyond the Mean in Fork-Join Queues"
     * (IFIP Performance 2015).</p>
     *
     * <p><strong>Requirements:</strong></p>
     * <ul>
     *   <li>Model must have valid Fork-Join topology: Source → Fork → K Queues → Join → Sink</li>
     *   <li>Solver must have been run first (runAnalyzer() called)</li>
     *   <li>Homogeneous service distributions across parallel queues</li>
     * </ul>
     *
     * @param percentiles Array of percentile levels (0-100 scale, e.g., {50, 90, 95, 99})
     * @return List of FJPercentileResult, one per job class, containing:
     *         - jobClass: class index
     *         - percentiles: requested percentile levels
     *         - values: computed percentile values
     *         - K: number of parallel queues
     *         - method: algorithm used (e.g., "FJ_NARE")
     * @throws IllegalStateException if model is not Fork-Join or solver not run
     */
    public List<FJPercentileResult> getPerctRespT(double[] percentiles) {
        if (percentileResults == null) {
            throw new IllegalStateException(
                "No percentile results available. " +
                "Ensure the model has a valid Fork-Join topology and run() has been called."
            );
        }

        // Interpolate stored results to requested percentiles
        List<FJPercentileResult> interpolated = new java.util.ArrayList<>();
        for (FJPercentileResult stored : percentileResults) {
            double[] values = new double[percentiles.length];
            for (int i = 0; i < percentiles.length; i++) {
                values[i] = interpolatePercentile(
                    stored.getPercentiles(),
                    stored.getRTp(),
                    percentiles[i]
                );
            }
            interpolated.add(new FJPercentileResult(
                stored.getK(),
                percentiles,
                values
            ));
        }
        return interpolated;
    }

    /**
     * Get response time percentiles using default values [50, 90, 95, 99]
     *
     * @return List of FJPercentileResult for default percentiles
     * @throws IllegalStateException if model is not Fork-Join or solver not run
     * @see #getPerctRespT(double[])
     */
    public List<FJPercentileResult> getPerctRespT() {
        return getPerctRespT(new double[]{50.0, 90.0, 95.0, 99.0});
    }

}



