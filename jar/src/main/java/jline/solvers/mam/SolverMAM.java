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
import jline.lang.constant.RoutingStrategy;
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

import jline.api.qsys.Qsys_is_retrial;
import jline.api.qsys.RetrialInfo;
import jline.api.sn.SnGetBufferSize;
import jline.api.sn.SnHasForkJoin;
import jline.api.sn.SnIsOpenModel;

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
                // Geometric and DiscreteUniform are the lattice laws of the
                // discrete-time path, solved by the Q-MAM discrete-time queues
                "Geometric", "DiscreteUniform",
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
                // "Retrial" is kept by the names that route to Solver_mam_retrial
                // only; methodFeatureSet withdraws it from the rest.
                "Retrial",
                // see _kb/06-solver-catalog.md for rationale
                "SetupDelayOff",
                // c-server stations (every analyzer reads sn.nservers; the slotted
                // path's single-server rule stays structural) and finite buffers as
                // LOSS buffers, see methodFeatureSet and bufferRefusal for the
                // closed-class refusal.
                "MultiServer", "FiniteCapacity"
        });
        return featSupported;
    }

    /**
     * Per-method feature deltas on the base MAM envelope. Only 'mna' resolves a
     * round-robin split (Npfqn_traffic_split_rr in Solver_mna_open); the closed
     * branch has no counterpart and is rejected in supportsModelMethod.
     * Mirrors the MATLAB/Python SolverMAM.getMethodFeatureSet.
     *
     * @param method the concrete method name
     * @return the per-method FeatureSet
     */
    public static FeatureSet methodFeatureSet(String method) {
        FeatureSet featSupported = SolverMAM.getFeatureSet();
        if ("mna".equals(method)) {
            featSupported.setTrue(new String[]{"RoutingStrategy_RROBIN"});
        }
        // dec.mmap (Solver_mam) is an OPEN-network departure-process fixed
        // point: it iterates on arrival streams a closed population does not
        // have, and its station ladder serves EXT, FCFS, HOL, FCFSPRPRIO and PS
        // only. Both restrictions are things the model HAS, so both belong here
        // rather than in supportsModelMethod; Solver_mam raises the matching
        // message when the method is named by hand. Until this delta existed
        // the gate offered dec.mmap on every closed model and the analyzer
        // answered with a table of zeros.
        if ("dec.mmap".equals(method)) {
            // Fork-join goes too: the sweep uses Solver_mam_traffic, the plain
            // traffic step, which has no synchronization -- which is why
            // Solver_mam_analyzer routes a fork-join topology from
            // 'default'/'dec.source' to Solver_mam_fj and never here. Left
            // declared, the gate offered dec.mmap on an open fork-join model and
            // the departure process it then built had no recurrent state.
            featSupported.setFalse(new String[]{
                    "ClosedClass", "SelfLoopingClass", "SchedStrategy_INF",
                    "Fork", "Join", "Forker", "Joiner"});
        }
        // Solver_mam_ldqbd is the only MAM algorithm that reads sn.lldscaling: it
        // applies the level-dependent factor in each departure block. 'default'
        // declares it because the single-class closed Delay+Queue shape routes
        // there, and supportsModelMethod refuses a load-dependent model outside
        // that shape -- a featset cannot see topology, and the dec.source
        // decomposition would otherwise solve every level at the nominal rate.
        if ("default".equals(method) || "ldqbd".equals(method)) {
            featSupported.setTrue(new String[]{"LoadDependence"});
        }
        // RETRIAL is served by Solver_mam_retrial alone, which "default" and
        // "dec.source" route to on the BMAP/PH/N/N shape and "retrial" names.
        // Every other algorithm reads no sn.retrial* field and would answer with
        // the refused jobs lost, so the base grant is withdrawn there;
        // supportsModelMethod refuses an orbit OUTSIDE that shape for the two
        // routing names, the same way it refuses reneging outside it.
        if (!"default".equals(method) && !"dec.source".equals(method)
                && !"retrial".equals(method)) {
            featSupported.setFalse(new String[]{"Retrial"});
        }
        // FINITECAPACITY: the open-network analyzers carry a finite buffer as a
        // LOSS buffer (Solver_mam_basic M/M/c/K and MMAP[K]/G/1/K, Solver_mam and
        // Solver_mna_open truncate-and-renormalize, Solver_mam_retrial's
        // bufferless N/N station), so the base envelope declares it. The two
        // chains that read no sn.cap withdraw it: the LD-QBD levels run to the
        // population or the cutoff, and the background chain to the state cap. A
        // buffer a CLOSED class can fill is refused for every method by
        // bufferRefusal (a closed job blocks, a loss formula does not).
        if ("ldqbd".equals(method) || "bgchain".equals(method)) {
            featSupported.setFalse(new String[]{"FiniteCapacity"});
        }
        return featSupported;
    }

    /**
     * Why no MAM method can answer a model whose finite buffer a CLOSED class
     * can fill, or "" when every binding buffer is reached by open classes only.
     *
     * <p>A MAM analyzer represents a finite buffer as a LOSS buffer:
     * Solver_mam_basic (default, dec.source, dec.poisson) solves an M/M/c/K or
     * an MMAP[K]/G/1/K, and Solver_mam (dec.mmap), Solver_mam_basic_mmap
     * (dec.source.mmap) and Solver_mna_open (mna) truncate and renormalize the
     * same way. That is the right model for an OPEN class, whose refused arrival
     * is lost. A closed job that finds no room BLOCKS instead -- LINE disables
     * the upstream departure and holds the job where it is -- and no MAM
     * analyzer blocks: the loss formulas answer a different system, and the
     * closed routes of "default" (Solver_mam_ldqbd, Solver_mam_bgchain,
     * Solver_mna_closed) read no sn.cap or sn.classcap at all. So the pair is
     * refused rather than answered.
     *
     * <p>ONE PREDICATE, TWO CALLERS: {@link #supportsModelMethod} reports it and
     * Solver_mam_analyzer raises it. Only a buffer that can BIND counts, which
     * is what {@link SnGetBufferSize} decides: refreshCapacity derives a finite
     * classcap (the chain population) at every station of every closed model, so
     * a plain finiteness test would refuse every closed model. A Cache builds
     * its own capped retrieval queues and is exempt, as in SnHasBlocking.
     * Mirrors MATLAB {@code mam_buffer_refusal}.
     *
     * @param sn the network structure
     * @return empty string when no closed class can fill a binding buffer
     */
    public static String bufferRefusal(NetworkStruct sn) {
        if (sn == null || sn.classcap == null || sn.nodetype == null) {
            return "";
        }
        for (int a = 0; a < sn.nodetype.size(); a++) {
            if (sn.nodetype.get(a) == NodeType.Cache) {
                return "";
            }
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            if (!Double.isFinite(SnGetBufferSize.snGetBufferSize(sn, ist))) {
                continue;
            }
            for (int r = 0; r < sn.nclasses && r < sn.classcap.getNumCols(); r++) {
                if (!Double.isFinite(sn.njobs.get(r)) || sn.classcap.get(ist, r) <= 0) {
                    continue;
                }
                return "Station " + sn.nodenames.get((int) sn.stationToNode.get(ist))
                        + " carries a finite capacity that binds for the closed class "
                        + sn.classnames.get(r) + ". SolverMAM represents a finite buffer as a "
                        + "LOSS buffer (M/M/c/K, MMAP[K]/G/1/K, truncate-and-renormalize), which "
                        + "is the open-class model: a closed job that finds no room blocks "
                        + "instead, and no MAM analyzer blocks, while the closed routes of the "
                        + "default method (ldqbd, bgchain, mna) read no capacity at all. Use "
                        + "SolverCTMC, SolverSSA or SolverLDES, or SolverMVA with method 'sqd' "
                        + "for blocking after service.";
            }
        }
        return "";
    }

    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        return SolverMAM.methodFeatureSet(method);
    }

    public NetworkStruct getStruct() {
        return model.getStruct(true);
    }

    public List<String> listValidMethods(Network model) {
        // SolverMAM.m verbatim, in its order, MINUS the RCAT names ("inap",
        // "inapplus", "inapinf", "exact") which moved to SolverAG;
        // supportsModelMethod redirects a caller that still asks for them.
        // 'dec.source.mmap' and 'retrial' are dispatched and so are advertised.
        return Arrays.asList("default", "dec.source", "dec.mmap", "dec.poisson", "mna", "ldqbd", "dec.source.mmap", "bgchain", "retrial");
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
        // The RCAT methods moved to SolverAG. Name them here rather than letting
        // listValidMethods report an unknown method, so a caller carrying an old
        // options.method is told where they went. The text comes from
        // unsupportedMethodReason, which checkDeclaredMethod also asks, so the
        // gate above the dispatcher and this one cannot drift into two answers.
        String moved = unsupportedMethodReason(method);
        if (!moved.isEmpty()) {
            return moved;
        }
        // G-network signals belong to the RCAT analyzer alone: no MAM algorithm
        // reads sn.issignal, so every one of them would solve the model with the
        // signals turned into ordinary customers and report that as the answer.
        if (hasSignalClass(getStruct())) {
            return "The " + method + " method does not support G-network signals: no MAM "
                    + "algorithm reads sn.issignal, so this method would solve the model with "
                    + "every signal turned into an ordinary customer. Use SolverAG, whose RCAT "
                    + "methods are the only ones that model signals.";
        }
        // A finite buffer a closed class can fill: no MAM analyzer blocks, so
        // every method refuses it (the open loss case is the one the envelope
        // declares). Same predicate Solver_mam_analyzer raises.
        String bufferWhy = bufferRefusal(getStruct());
        if (!bufferWhy.isEmpty()) {
            return bufferWhy;
        }
        // A retrial orbit reaches 'default' and 'dec.source' only through the
        // BMAP/PH/N/N shape their shared arm routes to Solver_mam_retrial on;
        // outside it the arm falls through to Solver_mam_basic, which reads no
        // sn.retrial* field and would answer with the refused jobs lost.
        if (("default".equals(method) || "dec.source".equals(method))
                && mamHasRetrial(getStruct())) {
            String retrialShape = mamRetrialRefusal(getStruct());
            if (!retrialShape.isEmpty()) {
                return "This model declares a retrial orbit outside the BMAP/PH/N/N bufferless "
                        + "shape (one open class, one bufferless Queue with a RETRIAL drop rule) "
                        + "that the " + method + " method solves through the retrial analyzer; "
                        + "dec.source would ignore the orbit. " + retrialShape + " Use SolverCTMC, "
                        + "SolverSSA, SolverJMT or SolverLDES.";
            }
        }
        // methodFeatureSet declares LoadDependence for 'default' because the
        // single-class closed Delay+Queue shape routes to Solver_mam_ldqbd, which
        // reads sn.lldscaling. Outside that shape 'default' falls through to the
        // dec.source decomposition, which never reads the field and would solve
        // every level at the nominal rate, so refuse by name here.
        if ("default".equals(method) && hasNonTrivialLoadDependence(getStruct())
                && !isClosedDelayQueue(getStruct())) {
            return "This model uses load-dependent service rates outside the single-class closed "
                    + "Delay+Queue shape that the default method routes to the level-dependent QBD "
                    + "(Solver_mam_ldqbd); the dec.source decomposition it would otherwise use does "
                    + "not read the scaling and would solve every level at the nominal rate. Call "
                    + "method 'ldqbd' directly, which refuses by name if the shape still does not fit.";
        }
        if ("mna".equals(method)
                && this.model.hasOpenClasses() && this.model.hasClosedClasses()) {
            return "The mna method does not support mixed open/closed models.";
        }
        if ("mna".equals(method)) {
            // the deterministic split is carried by the open traffic equations
            // only; Solver_mna_closed has no counterpart
            if (!this.model.hasOpenClasses() && hasRoundRobin(getStruct())) {
                return "The mna method supports round-robin routing in open models only.";
            }
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
        if ("retrial".equals(method)) {
            // A "must be present" rule, which a feature set cannot state: it
            // says which constructs are ACCEPTED, so it can refuse a model for
            // having something and never for lacking it. Solver_mam_retrial
            // needs the BMAP/PH/N/N bufferless retrial topology to analyze, and
            // mamRetrialRefusal is the same predicate Solver_mam_analyzer asks
            // before running it.
            String retrialWhy = mamRetrialRefusal(getStruct());
            if (!retrialWhy.isEmpty()) {
                return retrialWhy;
            }
        }
        if ("bgchain".equals(method)) {
            NetworkStruct snBg = getStruct();
            // The closed classes ARE the background chain, so a purely open model has
            // nothing to build it from. A purely CLOSED one is accepted: it is the
            // degenerate case where the chain answers alone, with no open work to
            // modulate it.
            if (SnIsOpenModel.snIsOpenModel(snBg)) {
                return "The bgchain method requires at least one closed class: the background chain IS "
                        + "the closed population vector, which a purely open model does not have. "
                        + "Use the dec.source method.";
            }
            // Priority disciplines need the per-class QBD of MMAPPH1PRPR, which has
            // no counterpart in the modulated level-dependent QBD this method builds.
            boolean prioSched = false;
            for (int i = 0; i < snBg.nstations; i++) {
                SchedStrategy sc = snBg.sched.get(snBg.stations.get(i));
                if (sc == SchedStrategy.HOL || sc == SchedStrategy.FCFSPRPRIO) {
                    prioSched = true;
                }
            }
            boolean mixedPrio = false;
            for (int k = 1; k < snBg.nclasses; k++) {
                if (snBg.classprio.get(0, k) != snBg.classprio.get(0, 0)) mixedPrio = true;
            }
            if (prioSched && mixedPrio) {
                return "The bgchain method does not support class priorities: it aggregates the open "
                        + "classes into one phase-type mixture per station, which cannot express a "
                        + "priority order. Use the dec.source method.";
            }
            if (SnHasForkJoin.snHasForkJoin(snBg)) {
                return "The bgchain method does not support fork-join: the background chain conserves "
                        + "the closed population per station, which a fork violates. Use the dec.source "
                        + "method.";
            }
        }
        return super.supportsModelMethod(method);
    }

    /** True for the RCAT family, which moved to SolverAG. Kept to redirect callers. */
    private static boolean isRcatMethod(String method) {
        return "inap".equals(method) || "inapplus".equals(method)
                || "inapinf".equals(method) || "exact".equals(method);
    }

    /**
     * The forwarding address for the RCAT names, which are SolverAG's now.
     *
     * <p>Asks nothing of the model, so {@code checkDeclaredMethod} can call it
     * before the struct is built; {@link #supportsModelMethod(String)} returns
     * the same string, so a caller gets one answer whichever gate it meets
     * first.
     */
    @Override
    protected String unsupportedMethodReason(String method) {
        if (!isRcatMethod(method)) {
            return "";
        }
        return "The " + method + " method moved to SolverAG: RCAT decomposes the model "
                + "into cooperating agents rather than decomposing traffic, and no MAM "
                + "algorithm shares its machinery. Use new SolverAG(model, \"" + method
                + "\").";
    }

    /**
     * Can the 'retrial' method answer this model? Empty string when it can.
     *
     * <p>THE RULE IS A "MUST BE PRESENT" ONE, which is why it cannot live in a
     * feature set: a FeatureSet says "I accept this construct", so it can refuse
     * a model for HAVING something and never for LACKING it.
     * Solver_mam_retrial needs the BMAP/PH/N/N bufferless retrial topology of
     * Dudin et al. (Mathematics 13(9), 2025) to analyze, and a model without one
     * is not a smaller retrial model, it is a different one.
     *
     * <p>ONE PREDICATE, TWO CALLERS: supportsModelMethod asks it so the method
     * is not offered by findSolver or SolverAUTO on a model it cannot answer,
     * and Solver_mam_analyzer's 'retrial' arm asks it so a caller naming the
     * method by hand gets the identical sentence.
     *
     * <p>Unlike MATLAB, this port has NO MAP/M/s+G reneging arm behind the
     * method, so a reneging patience does not make it applicable here.
     *
     * @param sn the model struct
     * @return empty string when applicable, otherwise the refusal
     */
    public static String mamRetrialRefusal(NetworkStruct sn) {
        RetrialInfo retInfo;
        try {
            retInfo = Qsys_is_retrial.qsys_is_retrial(sn);
        } catch (Exception e) {
            retInfo = new RetrialInfo();
        }
        if (retInfo.isRetrial()) {
            return "";
        }
        // qsys_is_retrial reports WHICH requirement the model missed (open
        // model, single class, a bufferless station, a retrial drop rule);
        // carrying it through is the difference between a bare no and a usable
        // answer. The wording is Solver_mam_retrial's own, so the gate and the
        // run agree.
        String detail = retInfo.getErrorMsg() == null ? "" : retInfo.getErrorMsg();
        return "No valid retrial configuration detected: " + detail;
    }

    /**
     * True when some station-class pair configures a retrial orbit
     * (Queue.setRetrial / setOrbit), read off sn.retrialProc, which the refresh
     * fills only for a configured delay that is not the Disabled placeholder.
     * retrialProc is the ambiguity-free test: retrialType is 0 both for "none"
     * and, in MATLAB, for an exponential delay (see CLAUDE.md).
     *
     * <p>{@link #mamRetrialRefusal} answers about the SHAPE Solver_mam_retrial
     * needs, not about whether an orbit EXISTS, and a model with an orbit
     * outside that shape falls through to Solver_mam_basic, which reads no
     * sn.retrial* field and answers with the refused jobs simply lost.
     * Mirrors MATLAB {@code mam_has_retrial}.
     *
     * @param sn the network structure
     * @return true when any station-class pair carries an orbit delay
     */
    public static boolean mamHasRetrial(NetworkStruct sn) {
        if (sn == null || sn.retrialProc == null) {
            return false;
        }
        for (Map<JobClass, MatrixCell> byClass : sn.retrialProc.values()) {
            if (byClass == null) {
                continue;
            }
            for (MatrixCell proc : byClass.values()) {
                if (proc != null) {
                    return true;
                }
            }
        }
        return false;
    }

    /** True when the model declares at least one G-network signal class. */
    private static boolean hasSignalClass(NetworkStruct sn) {
        if (sn.issignal == null || sn.issignal.isEmpty()) {
            return false;
        }
        for (int r = 0; r < sn.issignal.getNumRows(); r++) {
            if (sn.issignal.get(r, 0) > 0) {
                return true;
            }
        }
        return false;
    }

    /** True when some station declares a load-dependent scaling other than alpha == 1. */
    private static boolean hasNonTrivialLoadDependence(NetworkStruct sn) {
        if (sn.lldscaling == null || sn.lldscaling.isEmpty()) {
            return false;
        }
        for (int i = 0; i < sn.lldscaling.getNumRows(); i++) {
            for (int j = 0; j < sn.lldscaling.getNumCols(); j++) {
                if (sn.lldscaling.get(i, j) != 1.0) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * The exact regime Solver_mam_analyzer routes to Solver_mam_ldqbd on 'default':
     * one class, finite population, exactly two stations, one INF (Delay) and one
     * FCFS (Queue).
     */
    private static boolean isClosedDelayQueue(NetworkStruct sn) {
        if (sn.nclasses != 1 || sn.nstations != 2) {
            return false;
        }
        for (int k = 0; k < sn.njobs.length(); k++) {
            if (!Double.isFinite(sn.njobs.get(k))) {
                return false;
            }
        }
        int nDelay = 0, nQueue = 0;
        for (int i = 0; i < sn.nstations; i++) {
            SchedStrategy s = sn.sched.get(sn.stations.get(i));
            if (s == SchedStrategy.INF) {
                nDelay++;
            } else if (s == SchedStrategy.FCFS) {
                nQueue++;
            }
        }
        return nDelay == 1 && nQueue == 1;
    }

    /** True when any node dispatches round-robin for some class. */
    private static boolean hasRoundRobin(NetworkStruct sn) {
        if (sn.routing == null) {
            return false;
        }
        for (Map<JobClass, RoutingStrategy> perClass : sn.routing.values()) {
            for (RoutingStrategy rs : perClass.values()) {
                if (rs == RoutingStrategy.RROBIN) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * The process types the RCAT construction can give a phase dimension to.
     *
     * <p>After sn_nonmarkov_toph (which the RCAT methods run with phfit='ph' and
     * preserveDet=false) each of these holds a genuine (D0,D1) pair with
     * non-negative off-diagonal rates and a single arrival per epoch. The list is
     * an ALLOW-list on purpose: a process type nobody has checked against this
     * construction must be refused, not answered. Refused are the laws whose
     * matrices are not a generator (ME, RAP), those that are not
     * time-homogeneous (NHPP, MAPt, PHt), those that are not continuous-time
     * (DMAP), and those that arrive in batches (BMAP, MMAP), since a batch moves
     * the level by more than one.</p>
     */
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
    /**
     * Per-class phase-type service (sigma, S) at a station, read from sn.proc.
     *
     * <p>sn.proc holds the (D0, D1) MAP of the service process, so the PH pair
     * is sigma = map_pie(D0, D1) and S = D0. getProbMarg previously synthesized
     * a 1-phase exponential from sn.rates instead, which silently discarded
     * every non-exponential service law: an Erlang-2 station was analyzed as
     * M/M/1. Mirrors @SolverMAM/getProb.m, which reads sn.proc directly.</p>
     *
     * @param sn network struct
     * @param ist station index (0-based)
     * @param K number of classes
     * @param sigma out: per-class initial vectors
     * @param S out: per-class transient generators
     */
    private void mamServicePh(NetworkStruct sn, int ist, int K,
                              Map<Integer, Matrix> sigma, Map<Integer, Matrix> S) {
        Station station = sn.stations.get(ist);
        for (int k = 0; k < K; k++) {
            JobClass jobClass = sn.jobclasses.get(k);
            MatrixCell proc = (sn.proc != null && sn.proc.get(station) != null)
                    ? sn.proc.get(station).get(jobClass) : null;
            if (proc != null && proc.size() >= 2 && proc.get(0) != null
                    && proc.get(0).getNumRows() > 0) {
                Matrix D0 = proc.get(0);
                Matrix D1 = proc.get(1);
                sigma.put(k, jline.api.mam.Map_pie.map_pie(D0, D1));
                S.put(k, D0);
                continue;
            }
            // No MAP stored: fall back to the exponential implied by sn.rates.
            double rate = sn.rates.get(ist, k);
            if (Double.isNaN(rate) || rate <= 0) rate = 1.0;
            Matrix sigmaK = new Matrix(1, 1);
            sigmaK.set(0, 0, 1.0);
            sigma.put(k, sigmaK);
            Matrix Sk = new Matrix(1, 1);
            Sk.set(0, 0, -rate);
            S.put(k, Sk);
        }
    }

    /**
     * Joint (level, phase) state probability, port of @SolverMAM/getProb.m.
     *
     * <p>The level marginal comes from MMAPPH1FCFS over an aggregate MMAP built
     * from the class throughputs; the phase factor is the arrival-weighted
     * TIME-STATIONARY phase distribution map_prob, not map_pie. map_pie is the
     * embedded equilibrium at departure instants -- the phase a service STARTS
     * in -- so for an Erlang-2 it is [1 0] and the joint gave P(phase 2) = 0 at
     * every level for a server spending half its busy time in phase 2.</p>
     *
     * <p>Phases are taken independent of level, where an exact QBD would have
     * pi_n = pi_1 R^(n-1). That approximation is the reference's and is
     * reproduced rather than improved.</p>
     *
     * @param node node index (0-based)
     * @param state [level, phase] pair, or null for the whole matrix
     * @return probability of the state, or the (levels x phases) matrix
     */
    @Override
    public ProbabilityResult getProb(int node, Matrix state) {
        NetworkStruct sn = getStruct();
        if (node >= sn.nnodes) {
            throw new IllegalArgumentException(
                    "Node number exceeds the number of nodes in the model.");
        }
        int ist = (int) sn.nodeToStation.get(node);
        if (ist < 0) {
            throw new IllegalArgumentException("Specified node is not a station.");
        }
        int queueStations = 0;
        for (int i = 0; i < sn.nstations; i++) {
            if (sn.nodetype.get((int) sn.stationToNode.get(i)) == NodeType.Queue) queueStations++;
        }
        if (queueStations > 1) {
            throw new UnsupportedOperationException(
                    "getProb is not supported for networks with multiple queues in SolverMAM. "
                    + "The MAM solver uses QBD (quasi-birth-death) analysis, which is fundamentally "
                    + "a single-queue method. Use SolverCTMC or SolverSSA for state probabilities in "
                    + "networks with multiple queues.");
        }
        if (queueStations == 0) {
            throw new UnsupportedOperationException("Model does not contain any queue stations.");
        }
        if (result == null || ((SolverResult) result).QN == null) {
            runAnalyzer();
        }

        int K = sn.nclasses;
        SolverResult res = (SolverResult) result;
        double lambdaTotal = 0;
        for (int k = 0; k < K; k++) lambdaTotal += res.TN.get(ist, k);
        if (lambdaTotal < 1e-10) {
            throw new UnsupportedOperationException(
                    "getProb: the station carries no throughput, so no QBD exists.");
        }

        int totalPop = 0;
        for (int k = 0; k < K; k++) {
            double nk = sn.njobs.get(k);
            if (!Double.isInfinite(nk)) totalPop += (int) nk;
        }
        int maxLevel = totalPop > 0 ? totalPop + 1 : 100;

        MatrixCell D = new MatrixCell(K + 1);
        Matrix Dagg = new Matrix(1, 1);
        Dagg.set(0, 0, -lambdaTotal);
        D.set(0, Dagg);
        for (int k = 0; k < K; k++) {
            Matrix Dk = new Matrix(1, 1);
            Dk.set(0, 0, res.TN.get(ist, k));
            D.set(k + 1, Dk);
        }

        Map<Integer, Matrix> sigma = new HashMap<>();
        Map<Integer, Matrix> S = new HashMap<>();
        mamServicePh(sn, ist, K, sigma, S);

        Map<String, java.util.Map<Integer, Matrix>> mmapResult = MMAPPH1FCFS.MMAPPH1FCFS(
                D, sigma, S, null, maxLevel, null, null, false, false, 1e-14, null);

        Matrix pdistr = null;
        if (mmapResult.containsKey("ncDistr")) {
            for (int k = 0; k < K && pdistr == null; k++) {
                if (mmapResult.get("ncDistr").containsKey(k)) pdistr = mmapResult.get("ncDistr").get(k);
            }
        }
        if (pdistr == null) {
            throw new UnsupportedOperationException(
                    "getProb: MMAPPH1FCFS returned no queue-length distribution for this model.");
        }
        int levels = Math.max(pdistr.getNumRows(), pdistr.getNumCols());
        double[] p = new double[levels];
        double psum = 0;
        for (int n = 0; n < levels; n++) {
            double v = Math.abs(pdistr.getNumRows() == 1 ? pdistr.get(0, n) : pdistr.get(n, 0));
            p[n] = v;
            psum += v;
        }
        if (psum > 0) for (int n = 0; n < levels; n++) p[n] /= psum;

        int nPhases = 1;
        for (int k = 0; k < K; k++) nPhases = Math.max(nPhases, S.get(k).getNumRows());
        double[] avgPie = new double[nPhases];
        for (int k = 0; k < K; k++) {
            Matrix Sk = S.get(k);
            int m = Sk.getNumRows();
            Matrix D1k = new Matrix(m, m);
            Matrix sigK = sigma.get(k);
            for (int i = 0; i < m; i++) {
                double exit = 0;
                for (int j = 0; j < m; j++) exit -= Sk.get(i, j);
                for (int j = 0; j < m; j++) {
                    double sj = sigK.getNumRows() == 1 ? sigK.get(0, j) : sigK.get(j, 0);
                    D1k.set(i, j, exit * sj);
                }
            }
            Matrix piq = jline.api.mam.Map_prob.map_prob(Sk, D1k);
            for (int i = 0; i < m && i < nPhases; i++) {
                double v = piq.getNumRows() == 1 ? piq.get(0, i) : piq.get(i, 0);
                avgPie[i] += v * res.TN.get(ist, k) / lambdaTotal;
            }
        }
        double asum = 0;
        for (int i = 0; i < nPhases; i++) asum += avgPie[i];
        if (!(asum > 0) || Double.isNaN(asum)) {
            for (int i = 0; i < nPhases; i++) avgPie[i] = 1.0 / nPhases;
        } else {
            for (int i = 0; i < nPhases; i++) avgPie[i] /= asum;
        }

        if (state != null && state.getNumRows() * state.getNumCols() >= 2) {
            int level = (int) (state.getNumRows() == 1 ? state.get(0, 0) : state.get(0, 0));
            int phase = (int) (state.getNumRows() == 1 ? state.get(0, 1) : state.get(1, 0));
            if (level < 0 || phase < 1 || level >= levels || phase > nPhases) {
                return new ProbabilityResult(0.0);
            }
            double pv = p[level] * avgPie[phase - 1];
            return new ProbabilityResult(pv);
        }
        Matrix joint = new Matrix(levels, nPhases);
        for (int n = 0; n < levels; n++) {
            for (int i = 0; i < nPhases; i++) joint.set(n, i, p[n] * avgPie[i]);
        }
        ProbabilityResult out = new ProbabilityResult(joint);
        out.nodeIndex = node;
        return out;
    }

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

            // Read the real service law from sn.proc; synthesizing a 1-phase
            // exponential from sn.rates analyzed every Erlang/HyperExp/PH
            // station as M/M/1 and silently discarded its phase structure.
            mamServicePh(sn, node, K, sigma, S);

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
