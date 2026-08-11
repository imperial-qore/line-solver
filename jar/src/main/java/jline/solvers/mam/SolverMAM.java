/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.GlobalConstants;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.solvers.*;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.DistributionResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.api.sn.SnGetResidTFromRespT.snGetResidTFromRespT;
import static jline.io.InputOutput.*;
import static jline.solvers.mam.analyzers.Solver_mam_analyzer.solver_mam_analyzer;
import static jline.solvers.mam.handlers.Solver_mam_passage_time.solver_mam_passage_time;
import static jline.api.fj.FJValidation.isHomogeneous;
import static jline.solvers.mam.handlers.Solver_mam_fj.solver_mam_fj;
import static jline.solvers.mam.handlers.Solver_mam_fj.interpolatePercentile;
import static jline.solvers.mam.handlers.Solver_mam_ldqbd_transient.solver_mam_ldqbd_transient;
import static jline.api.sn.SnNonmarkovToPh.snNonmarkovToPh;
import jline.solvers.mam.handlers.TransientResult;
import jline.solvers.mam.handlers.Solver_mam_transient_qbd;

import jline.api.fj.FJInfo;
import jline.lib.fjcodes.MainFJ.FJPercentileResult;
import jline.lib.butools.MMAPPH1FCFS;
import jline.lang.constant.NodeType;
import jline.util.Pair;
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
        super(model, "SolverMAM", defaultOptions().method(method));
        this.result = new MAMResult();
    }

    public SolverMAM(Network model, Object... varargin) {
        this(model, defaultOptions());
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
     * Intermediate quantities of the matrix-analytic analysis of a single-queue
     * model, in addition to the mean performance measures returned by getAvg.
     *
     * <p>Mean values alone hide the objects the method is actually built on, so a
     * matrix-analytic result cannot be inspected, taught, or checked against a
     * published derivation. This accessor returns them.</p>
     *
     * <p>For a BMAP (or MAP) arrival stream feeding an exponential single server the
     * result is that of {@link jline.api.qsys.Qsys_bmapm1} and carries the
     * M/G/1-type quantities: the phase-process stationary vectors theta and alpha,
     * the randomized blocks A0, A1, B0 and Bk, the matrix G, the drift, the measured
     * decay rate and the level probabilities. For a retrial station the result is
     * that of {@link jline.api.qsys.Qsys_bmapphnn_retrial} and carries the
     * orbit-level stationary distribution together with the truncation level and its
     * residual.</p>
     *
     * @return the matrix-analytic internals
     */
    public MAMInternals getMAMResult() {
        NetworkStruct snl = this.model.getStruct(true);

        // A retrial station carries its own engine, whose result object already
        // exposes the orbit-level internals.
        jline.api.qsys.RetrialInfo retInfo = jline.api.qsys.Qsys_is_retrial.qsys_is_retrial(snl);
        if (retInfo.isRetrial()) {
            // The accessor reports internals; the engine's progress trace is noise here.
            SolverOptions opts = this.options;
            jline.VerboseLevel savedVerbose = opts.verbose;
            opts.verbose = jline.VerboseLevel.SILENT;
            MAMResult res;
            try {
                res = jline.solvers.mam.handlers.Solver_mam_retrial.solver_mam_retrial(snl, opts);
            } finally {
                opts.verbose = savedVerbose;
            }
            return new MAMInternals(null, res.retrialInternals);
        }

        // Otherwise: BMAP/MAP arrivals into a single exponential server.
        int sourceIdx = -1;
        int queueIdx = -1;
        for (int ist = 0; ist < snl.nstations; ist++) {
            int nodeIdx = (int) snl.stationToNode.get(ist);
            if (snl.nodetype.get(nodeIdx) == jline.lang.constant.NodeType.Source) {
                sourceIdx = ist;
            } else if (snl.nodetype.get(nodeIdx) == jline.lang.constant.NodeType.Queue) {
                if (queueIdx < 0) {
                    queueIdx = ist;
                } else {
                    throw new RuntimeException(
                            "getMAMResult exposes the matrix-analytic internals of a single-queue model only.");
                }
            }
        }
        if (sourceIdx < 0 || queueIdx < 0) {
            throw new RuntimeException("getMAMResult requires an open model with one Source and one Queue.");
        }
        if (snl.nclasses > 1) {
            throw new RuntimeException(
                    "getMAMResult exposes the matrix-analytic internals of a single-class model only.");
        }
        if (snl.nservers.get(queueIdx) != 1) {
            throw new RuntimeException("getMAMResult requires a single-server queue.");
        }

        Station sourceStation = snl.stations.get(sourceIdx);
        Station queueStation = snl.stations.get(queueIdx);
        JobClass jobclass = snl.jobclasses.get(0);
        MatrixCell arrivalProc = snl.proc.get(sourceStation).get(jobclass);
        if (arrivalProc == null || arrivalProc.size() < 2) {
            throw new RuntimeException(
                    "The arrival process has no Markovian (D0,D1,...) representation.");
        }
        MatrixCell serviceProc = snl.proc.get(queueStation).get(jobclass);
        if (serviceProc == null || serviceProc.size() < 1 || serviceProc.get(0).getNumRows() != 1) {
            throw new RuntimeException("getMAMResult exposes the M/G/1-type internals for exponential "
                    + "service only; the queue has a multi-phase service process.");
        }
        double mu = -serviceProc.get(0).get(0, 0);

        Matrix[] D = new Matrix[arrivalProc.size()];
        for (int k = 0; k < D.length; k++) {
            D[k] = arrivalProc.get(k);
        }
        return new MAMInternals(jline.api.qsys.Qsys_bmapm1.qsys_bmapm1(D, mu), null);
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
                "Fork", "Join", "Forker", "Joiner",  // Fork-Join support (via Solver_mam_fj)
                "Delay", "DelayStation", "Queue",
                "APH", "Coxian", "Erlang", "Exp", "HyperExp", "MMPP2", "MAP", "MMAP", "DMAP", "ME", "RAP",
                "PH", "BMAP",  // BMAP/PH/N/N retrial queue (Solver_mam_retrial)
                "Det", "Gamma", "Lognormal", "Pareto", "Uniform", "Weibull",
                "StatelessClassSwitcher", "InfiniteServer",
                "ClassSwitch",
                "SharedServer", "Buffer", "Dispatcher",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS", "SchedStrategy_HOL",
                "SchedStrategy_FCFS",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ClosedClass", "SelfLoopingClass",
                "OpenClass",
                "OpenSignal", "ClosedSignal", // G-network signals (Solver_mam_build_ag)
                "SignalType_NEGATIVE", "SignalType_CATASTROPHE",
                "SignalBatchRemoval", // AG reads sn.signalremdist
                "Retrial",
                // see _kb/06-solver-catalog.md for rationale
                "SetupDelayOff"
        });
        return featSupported;
    }

    public NetworkStruct getStruct() {
        return model.getStruct(true);
    }

    public List<String> listValidMethods(Network model) {
        return Arrays.asList("default", "dec.source", "dec.mmap", "dec.poisson", "mna", "inap", "inapplus", "inapinf", "exact", "ldqbd", "retrial");
    }

    public List<String> listValidMethods() {
        return listValidMethods(null);
    }

    @Override
    public void getTranAvg() {
        String savedMethod = options.method;
        options.method = "ldqbd";
        super.getTranAvg();
        options.method = savedMethod;
    }

    @Override
    public void runAnalyzer() {
        // Propagate solver verbose level to global
        if (this.options != null) {
            GlobalConstants.Verbose = options.verbose;
        }
        // see _kb/06-solver-catalog.md for rationale
        if (this.enableChecks) {
            String reason = this.supportsModelMethod(this.options.method);
            if (!reason.isEmpty()) {
                line_error(mfilename(new Object() {
                }), "This model contains features not supported by the solver. " + reason);
                return;
            }
        }

        // Finite Capacity Region: MAM does not enforce the aggregate per-region
        // job limit and would silently return the unconstrained answer.
        if (this.model.getStruct(false).nregions > 0) {
            throw new RuntimeException("This model uses a Finite Capacity Region (addRegion), "
                    + "which is not supported by SolverMAM (the region's aggregate job limit "
                    + "is not enforced). Use SolverCTMC, SolverJMT, SolverSSA or SolverLDES, "
                    + "or setCapacity for a single-station limit.");
        }

        double start = System.nanoTime();
        NetworkStruct sn = getStruct();

        // Check if transient analysis is requested
        boolean isTran = options.timespan != null && options.timespan.length >= 2
                && !Double.isInfinite(options.timespan[1]);
        if (isTran && options.method.equals("ldqbd")) {
            TransientResult tranRes;
            if (Solver_mam_transient_qbd.applicable(sn)) {
                // Correlated MAP arrival/service or non-Poisson arrival: use the
                // Laplace-domain transient QBD solver on the true MAP blocks.
                line_debug(options.verbose, "MAM: transient analysis via Laplace transient QBD");
                tranRes = Solver_mam_transient_qbd.solver_mam_transient_qbd(sn, options);
            } else {
                line_debug(options.verbose, "MAM: transient analysis via standard QBD");
                NetworkStruct snPh = snNonmarkovToPh(sn, options);
                tranRes = solver_mam_ldqbd_transient(snPh, options);
            }
            int M = sn.nstations;
            int K = sn.nclasses;
            Matrix[][] Qt = tranRes.getQt();
            Matrix[][] Ut = tranRes.getUt();
            Matrix[][] Tt = tranRes.getTt();
            Matrix[][] Rt = new Matrix[M][K];
            Matrix[][] Ct = new Matrix[1][K];
            Matrix[][] Xt = new Matrix[1][K];
            double finish = System.nanoTime();
            double runtime = (finish - start) / 1000000000.0;
            setTranAvgResults(Qt, Ut, Rt, Tt, Ct, Xt, runtime);
            return;
        }

        // see _kb/06-solver-catalog.md for rationale
        Pair<Boolean, FJInfo> fjCheck = isHomogeneous(sn);
        boolean fjRouted = "default".equals(options.method) || "dec.source".equals(options.method);
        if (fjCheck.getFirst() && fjRouted) {
            line_debug(options.verbose, String.format("Detected Fork-Join topology with K=%d parallel queues", fjCheck.getSecond().getK()));
            line_debug(options.verbose, "Using FJ_codes algorithm for percentile analysis");

            try {
                // Use FJ solver
                jline.solvers.mam.MAMFJResult fjResult = solver_mam_fj(
                    sn,
                    new double[]{0.50, 0.90, 0.95, 0.99},  // Default percentiles (0-1 probability scale)
                    0,       // C=0: auto-select QBD truncation from utilization (was hardcoded 100)
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
     * Method-aware gate for the genuine per-method restrictions of the MAM
     * analyzer (mirrors the guards in Solver_mam_analyzer and the MATLAB
     * SolverMAM.supportsModelMethod): 'mna' does not support mixed open/closed
     * models, and 'ldqbd' requires a single-class model. All other methods rely
     * on the coarse MAM feature set and the analyzer's topology-based routing
     * (a Fork-Join model is routed to the FJ solver, not rejected). Returns an
     * empty string when supported, otherwise the reason.
     */
    @Override
    public String supportsModelMethod(String method) {
        if ("mna".equals(method)
                && this.model.hasOpenClasses() && this.model.hasClosedClasses()) {
            return "The mna method does not support mixed open/closed models.";
        }
        if ("mna".equals(method)) {
            // see _kb/06-solver-catalog.md for rationale
            NetworkStruct snMna = getStruct();
            Matrix Vmna = snMna.visits.get(0);
            for (int c = 1; c < snMna.nchains; c++) {
                Vmna = Vmna.add(1.0, snMna.visits.get(c));
            }
            for (int k = 0; k < snMna.nclasses; k++) {
                // see _kb/06-solver-catalog.md for rationale
                if (Double.isInfinite(snMna.njobs.get(k))) {
                    continue;
                }
                int visited = -1;
                boolean multiple = false;
                for (int i = 0; i < snMna.nstations; i++) {
                    if (Vmna.get(i, k) > GlobalConstants.FineTol) {
                        if (visited >= 0) { multiple = true; break; }
                        visited = i;
                    }
                }
                if (!multiple && visited >= 0) {
                    SchedStrategy sc = snMna.sched.get(snMna.stations.get(visited));
                    if (sc != SchedStrategy.INF && sc != SchedStrategy.EXT) {
                        return "The mna method does not support self-looping classes (class "
                                + (k + 1) + " is confined to station " + (visited + 1)
                                + " with no inter-station flow to decompose). Use the dec.source method.";
                    }
                }
            }
        }
        if ("ldqbd".equals(method) && getStruct().nclasses != 1) {
            return "The ldqbd method requires a single-class model.";
        }
        if ("inap".equals(method) || "inapplus".equals(method)
                || "inapinf".equals(method) || "exact".equals(method)) {
            String rcatReason = rcatUnsupportedProcessReason(method);
            if (rcatReason != null) {
                return rcatReason;
            }
        }
        return supports((Network) this.model) ? "" : "Some features are not supported by the chosen solver.";
    }

    /**
     * Reject non-exponential processes for the RCAT methods, which model each
     * station-class by its mean rate only.
     *
     * <p>Solver_mam_build_ag builds one scalar birth-death chain per
     * (station,class) from sn.rates: the state is the queue length, with no
     * service-phase dimension. Every process is therefore collapsed to its mean
     * rate, so a non-exponential arrival or service process would be answered as
     * if it were exponential. Measured vs SolverCTMC on an open M/PH/1, INAP
     * returns the M/M/1 queue length for EVERY scv (1.000000 at rho=0.5,
     * 4.000000 at rho=0.8): 14% error at scv=0.5 and 43% at scv=4.0. Reject
     * rather than mis-answer.</p>
     *
     * @return the reason string, or null when every process is exponential
     */
    private String rcatUnsupportedProcessReason(String method) {
        NetworkStruct sn = getStruct();
        if (sn.procid == null) {
            return null;
        }
        for (int i = 0; i < sn.nstations; i++) {
            Station st = sn.stations.get(i);
            Map<JobClass, ProcessType> byClass = sn.procid.get(st);
            if (byClass == null) {
                continue;
            }
            boolean isSource = st instanceof Source;
            for (int r = 0; r < sn.nclasses; r++) {
                JobClass jc = sn.jobclasses.get(r);
                ProcessType pt = byClass.get(jc);
                if (pt == null || pt == ProcessType.EXP) {
                    continue;
                }
                // see _kb/06-solver-catalog.md for rationale
                if (!isSource && sn.issignal != null && sn.issignal.get(r, 0) > 0) {
                    continue;
                }
                // Only a process that is actually in use can mis-answer.
                double rate = sn.rates.get(i, r);
                if (!Double.isFinite(rate) || rate <= 0) {
                    continue;
                }
                return "The " + method + " method supports exponential processes only "
                        + "(RCAT models each station-class by its mean rate, with no "
                        + "service-phase dimension), but station " + (i + 1) + " class "
                        + (r + 1) + " is " + ProcessType.toText(pt) + ". Use the dec.source "
                        + "method for non-exponential models.";
            }
        }
        // see _kb/06-solver-catalog.md for rationale
        for (int i = 0; i < sn.nstations; i++) {
            double c = sn.nservers.get(i);
            if (Double.isFinite(c) && c > 1) {
                return "The " + method + " method supports single-server stations only "
                        + "(RCAT does not model sn.nservers, so a multiserver station is "
                        + "driven at rho = lambda/mu instead of lambda/(c*mu)), but station "
                        + (i + 1) + " has " + (int) c + " servers. Use the dec.source method "
                        + "for multiserver models.";
            }
        }
        return null;
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
            int nodeIdx = (int) sn.stationToNode.get(i);
            if (sn.nodetype.get(nodeIdx) == NodeType.Queue) {
                queueStations++;
            }
        }

        if (queueStations > 1) {
            throw new UnsupportedOperationException(
                "getProbMarg is not supported for networks with multiple queues in SolverMAM. " +
                "The MAM solver uses QBD (quasi-birth-death) analysis, which is fundamentally a single-queue method. " +
                "Use SolverCTMC or SolverSSA for marginal probabilities in networks with multiple queues.");
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
            Map<String, java.util.Map<Integer, Matrix>> mmapResult = MMAPPH1FCFS.MMAPPH1FCFS(
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
     * <p>Delegates to {@link #getProbMarg(int, int, Matrix)} with {@code state_m == null}
     * (all levels). Supported for single-queue models via QBD/MMAPPH1FCFS analysis; a
     * network with multiple queue stations is rejected there, matching the MATLAB
     * {@code @SolverMAM/getProbMarg.m} behaviour.</p>
     *
     * @param node Station/node index (0-based)
     * @param jobclass Job class index (0-based)
     * @return Probability result containing the marginal queue-length distribution
     * @throws UnsupportedOperationException if the model has more than one queue station
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
                // Both requested and stored percentile levels are on the 0-100
                // scale (MainFJ stores pers*100).
                values[i] = interpolatePercentile(
                    stored.percentiles,
                    stored.RTp,
                    percentiles[i]
                );
            }
            interpolated.add(new FJPercentileResult(
                stored.K,
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


    /**
     * Bundled third-party libraries used by the matrix-analytic solver:
     * MAMSolver for the M/G/1 and GI/M/1-type decompositions, Q-MAM for the
     * RCAT-based methods, and BUTools wherever a MAP/PH process is present.
     * Mirrors SolverMAM.getLibrariesUsed in MATLAB.
     */
    @Override
    public java.util.List<String> getLibrariesUsed(NetworkStruct sn, SolverOptions options) {
        java.util.List<String> libs = new java.util.ArrayList<String>();
        String method = (options == null || options.method == null) ? "" : options.method;
        if (method.equals("default") || method.equals("dec.source") || method.equals("dec.mmap")
                || method.equals("dec.poisson") || method.equals("dec.source.mmap")) {
            libs.add("MAMSolver");
        }
        if (method.equals("mna") || method.equals("inap") || method.equals("inapplus")
                || method.equals("inapinf")) {
            libs.add("Q-MAM");
        }
        if (sn != null && sn.proc != null && !sn.proc.isEmpty()) {
            libs.add("BUTools");
        }
        return libs;
    }
}
