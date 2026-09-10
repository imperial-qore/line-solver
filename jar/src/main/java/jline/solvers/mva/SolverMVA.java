/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.GlobalConstants;
import jline.api.sn.SnHasBurstyArrival;
import jline.api.sn.SnIsMm1kLoss;
import jline.api.sn.SnPatienceHandles;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.io.Ret.ProbabilityResult;
import jline.solvers.mva.handlers.MVARunner;
import jline.solvers.mva.handlers.Solver_mapqn;
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
import java.util.Collections;
import java.util.Map;

import static jline.util.Maths.logBinomial;
import static jline.util.Maths.factln;
import static jline.api.sn.SnIsOpenModel.snIsOpenModel;
import static jline.api.sn.SnOpenProbTerms.snOpenProbTerms;
import static jline.api.sn.SnHasClassdepRouting.snHasClassdepRouting;
import static jline.api.sn.SnHasSjn.snHasSjn;
import static jline.api.sn.SnHasOpenClasses.snHasOpenClasses;
import static jline.api.sn.SnHasProductForm.snHasProductForm;
import static jline.api.sn.SnHasProductFormNotHetFCFS.snHasProductFormNotHetFCFS;
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
                "SchedStrategy_INF", "SchedStrategy_PS", "SchedStrategy_FCFSPRPRIO",
                "SchedStrategy_DPS", "SchedStrategy_FCFS", "SchedStrategy_SIRO", "SchedStrategy_HOL",
                "SchedStrategy_LCFS", "SchedStrategy_LCFSPR", "SchedStrategy_POLLING",
                "SchedStrategy_OI", "SchedStrategy_PAS",
                // size-based M/G/1 disciplines, served by
                // Solver_mva_qsys_sizebased_analyzer (Wierman and
                // Harchol-Balter, SIGMETRICS 2003)
                "SchedStrategy_SRPT", "SchedStrategy_PSJF", "SchedStrategy_FB",
                "SchedStrategy_LRPT", "SchedStrategy_SETF",
                "SchedStrategy_SJF", // closed models only (Solver_mva_sjn_analyzer)
                "Fork", "Forker", "Join", "Joiner",
                // quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
                "JoinPartial",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_LRU",
                "ReplacementStrategy_HLRU",
                "MMAP", // marked MAP sources (cache LRU via cache_ttl_lrum_map)
                "ClosedClass", "SelfLoopingClass", "OpenClass", "Replayer",
                "LoadDependence", "ClassDependence", "JointDependence",
                // c-server stations: the exact recursion, every AMVA kernel,
                // qna/rqt and the M/M/k and G/G/k closed forms carry the count;
                // methodFeatureSet withdraws it from the single-server names.
                // FiniteCapacity is deliberately NOT here (see methodFeatureSet
                // and finiteCapacityReason).
                "MultiServer"
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
        method = baseMethod(method);
        if (isClosedPopulationMethod(method)) {
            // The closed-population AMVA family estimates the arrival-instant
            // queue length as a function of the population vector N and is handed
            // (L, N, Z) alone, so an open chain gives it nothing to recur on:
            // Solver_amva has no arm for any of these outside its closed
            // product-form branch, and falling through returned the qd-family
            // answer, or a table of zeros, under their name. Strict product form
            // has no registry name and is applied by closedPopulationReason.
            featSupported.setFalse(new String[]{"OpenClass"});
        }
        // Solver_mvald_analyzer serves a load-, class- or joint-dependent model
        // through "exact"/"mva" (load dependence only, it has no class- or
        // joint-dependent recursion) and through the default/amva/qd/lin/qdlin
        // arms, and refuses every other name by name. The queueing-system closed
        // forms are intercepted by MVARunner upstream of that analyzer and keep
        // the base envelope.
        List<String> ldRefused = new ArrayList<String>(CLOSED_POPULATION_METHODS);
        ldRefused.addAll(Arrays.asList("sum", "esum", "mvac", "qli", "fli",
                "gflin", "egflin", "qna", "rqna", "rqt"));
        if (ldRefused.contains(method)) {
            featSupported.setFalse(new String[]{"LoadDependence", "ClassDependence", "JointDependence"});
        } else if ("mva".equals(method) || "exact".equals(method)) {
            featSupported.setFalse(new String[]{"ClassDependence", "JointDependence"});
        }
        if (!"default".equals(method) && !"exact".equals(method)) {
            // An order-independent or pass-and-swap station is served by
            // Solver_mva_oi_analyzer alone, which MVARunner reaches only under
            // "default" or "exact"; every other name is refused there by name, so
            // it must not be advertised for such a model.
            featSupported.setFalse(new String[]{"SchedStrategy_OI", "SchedStrategy_PAS"});
        }
        if ("sum".equals(method) || "esum".equals(method)) {
            // Solver_mva_sum passes each station to the summation kernel as an
            // INF, PS, LCFS-PR, FCFS or SIRO centre and refuses the rest by name.
            featSupported.setFalse(NON_BCMP_SCHED_FEATURES);
        } else if ("mvac".equals(method)) {
            // pfqn_mvac recurs on the closed chains over single-server fixed-rate
            // (SSFR) and infinite-server centres; Solver_mva refuses every other
            // discipline by name.
            featSupported.setFalse(new String[]{"OpenClass"});
            featSupported.setFalse(NON_BCMP_SCHED_FEATURES);
        }
        if ("qna".equals(method)) {
            // round-robin dispatching enters as a deterministic traffic split
            // (Npfqn_traffic_split_rr), which the exact-MVA paths have no
            // counterpart for
            featSupported.setTrue(new String[]{"RoutingStrategy_RROBIN"});
            featSupported.setFalse(new String[]{"ClosedClass", "SelfLoopingClass"});
            // Solver_qna's station loop has an arm for INF, PS and FCFS and none
            // for anything else, so on a SIRO, LCFS-PR, HOL or priority station it
            // left that row of Q, U, R and T at zero and reported the table as a
            // solution.
            featSupported.setFalse(new String[]{"SchedStrategy_SIRO", "SchedStrategy_LCFSPR"});
            featSupported.setFalse(NON_BCMP_SCHED_FEATURES);
        } else if ("rqna".equals(method)) {
            featSupported.setTrue(new String[]{"MAP", "MMPP2", "MMAP", "RAP"});
            featSupported.setFalse(new String[]{"ClosedClass", "SelfLoopingClass"});
        } else if ("mapqn".equals(method)) {
            // the horizontal-cut MVA consumes a MAP service natively (a closed
            // delay + FCFS queue model, see Solver_mapqn.mapqnReason); declaring
            // MAP here is what keeps NetworkSolver.needsMapEnv from routing the
            // model through its random-environment image
            featSupported.setTrue(new String[]{"MAP", "MMPP2"});
            featSupported.setFalse(new String[]{"OpenClass", "Source", "Sink", "Fork", "Forker", "Join", "Joiner",
                    "JoinPartial", "ClassSwitch", "StatelessClassSwitcher", "Cache", "CacheClassSwitcher",
                    "CacheRetrieval", "LoadDependence", "ClassDependence", "JointDependence",
                    "SchedStrategy_PS", "SchedStrategy_SIRO", "SchedStrategy_LCFSPR", "SchedStrategy_SRPT",
                    "SchedStrategy_PSJF", "SchedStrategy_FB", "SchedStrategy_LRPT", "SchedStrategy_SETF",
                    "SchedStrategy_OI", "SchedStrategy_PAS"});
            featSupported.setFalse(NON_BCMP_SCHED_FEATURES);
        } else if ("rqt".equals(method)) {
            // robust queueing theory: single-class open networks, the primitives
            // entering the uncertainty sets are two moments
            featSupported.setFalse(new String[]{"ClosedClass", "SelfLoopingClass"});
        } else if ("erlanga".equals(method) || "mgisrgi".equals(method)) {
            // The only analytical methods in LINE that carry an abandonment
            // rate. Reneging stays OUT of the base MVA envelope: every other
            // method here would silently ignore the patience law and report the
            // no-abandonment answer.
            featSupported.setTrue(new String[]{"Reneging"});
            featSupported.setFalse(new String[]{"ClosedClass", "SelfLoopingClass"});
        }
        if ("rqna".equals(method) || "rqt".equals(method)) {
            // A Join is a synchronisation node, not a queue: it carries no service
            // process, so the index-of-dispersion curve these two read off every
            // station does not exist for it, and neither analyzer has a
            // synchronisation term to put in its place. QNA keeps Fork/Join -- its
            // station loop has an explicit Join arm.
            featSupported.setFalse(new String[]{"Fork", "Forker", "Join", "Joiner", "JoinPartial"});
        }
        // MULTISERVER (registry name since 2026-09-05). The single-server
        // recursions: AQL, QSA and Tay (closedPopulationReason), MVAC's SSFR
        // chain recursion (mvacReason), RQNA's GI/G/1 workload
        // (singleClassOpenReason), Kant's SJN recursion and the single-server
        // closed forms of the queueing-system analyzer, every M/G/1, G/M/1 and
        // G/G/1 name. Each predicate stays, wording the refusal for the run;
        // the delta is what makes it nameable. RQT, QNA, M/M/k, G/G/k and the
        // rest of the envelope carry a server count.
        if ("aql".equals(method) || "qsa".equals(method) || "tay".equals(method)
                || "mvac".equals(method) || "rqna".equals(method)
                || "sjn.mva".equals(method) || "sjn.amva".equals(method)
                || "mm1".equals(method) || "mg1".equals(method) || "mgi1".equals(method)
                || "gm1".equals(method) || "gim1".equals(method)
                || method.startsWith("gig1")) {
            featSupported.setFalse(new String[]{"MultiServer"});
        }
        // FINITECAPACITY (registry name since 2026-09-05) is NOT in the base
        // envelope: the product-form recursions solve a buffer away, which is
        // what finiteCapacityReason refuses. The names that honour one are
        // granted it here, and the structural predicate keeps the shape half of
        // each rule: "default" and "sqd" reach Solver_sqd, the one
        // Blocking-After-Service arm.
        if ("default".equals(method) || "sqd".equals(method)) {
            featSupported.setTrue(new String[]{"FiniteCapacity"});
        }
        return featSupported;
    }

    /**
     * The per-method envelope plus the one MVA rule that is judged on the MODEL
     * rather than on the name, which is why it cannot live in the static
     * {@link #methodFeatureSet} table.
     *
     * <p>A single-station M/M/1/K with tail drop keeps a closed form under every
     * name but "exact" (the moment-based {@code Qsys_mg1k_loss_mgs} branch of
     * {@code Solver_mva_qsys_analyzer}, exact at scv = 1 only), so the grant is
     * judged on {@link SnIsMm1kLoss}. {@code MVARunner} reaches that branch by
     * SHAPE rather than by name and {@link #finiteCapacityReason} carries the
     * matching exemption, so the report and the run agree on it.
     *
     * @param method the concrete method name
     * @return the per-method FeatureSet for this model
     */
    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        if (!(this.model instanceof Network)) {
            return SolverMVA.methodFeatureSet(method);
        }
        return SolverMVA.methodFeatureSet(((Network) this.model).getStruct(false), method);
    }

    /**
     * The same envelope, asked of a struct rather than of a solver instance, so
     * that the run gate in {@code MVARunner} and the report gate in
     * {@link #getMethodFeatureSet} apply ONE copy of the model-judged rule. The
     * runner reads the static table alone until this exists, and would refuse
     * the M/M/1/K shape under every name its own dispatch claims by shape.
     *
     * @param sn     the network structure, or null to get the name-only envelope
     * @param method the concrete method name
     * @return the per-method FeatureSet for this model
     */
    public static FeatureSet methodFeatureSet(NetworkStruct sn, String method) {
        FeatureSet featSupported = SolverMVA.methodFeatureSet(method);
        if (sn != null && !"exact".equals(baseMethod(method)) && SnIsMm1kLoss.snIsMm1kLoss(sn)) {
            featSupported.setTrue(new String[]{"FiniteCapacity"});
        }
        return featSupported;
    }

    /**
     * Feature-driven resolution of {@code options.method="default"}: a bursty
     * single-class open network has a non-renewal (MAP/MMPP) arrival process
     * whose autocorrelation a two-moment method cannot capture, so the default
     * dispatch of Solver_mva_analyzer selects RQNA. Naming that here is what
     * lets the feature gate admit the MAP family on the path that will actually
     * run: without it "default" carried the base envelope, the MAP/MMPP arrival
     * read as unsupported, and the model went to the random-environment
     * approximation instead of the exact-in-IDC analyzer. Mirrors the
     * MATLAB/Python SolverMVA.resolveMethod.
     *
     * @param options the solver options
     * @return the concrete method the analyzer will dispatch
     */
    @Override
    public String resolveMethod(SolverOptions options) {
        String method = (options == null) ? "default" : options.method;
        if (this.model == null) {
            return method;
        }
        return resolveMethodForStruct(this.getStruct(), method);
    }

    /**
     * The concrete method the MVA analyzer will dispatch, given a struct and a
     * requested name. Static and struct-only so that {@link
     * jline.solvers.mva.handlers.MVARunner} gates on exactly the name this
     * solver resolves: two copies of this rule drift, and a gate that resolves
     * differently from the analyzer either refuses a model the analyzer could
     * solve or admits one it cannot.
     *
     * @param sn     the network struct
     * @param method the requested method name
     * @return the concrete method name
     */
    public static String resolveMethodForStruct(NetworkStruct sn, String method) {
        if (sn == null || !"default".equals(method)) {
            return method;
        }
        boolean allOpen = true;
        for (int r = 0; r < sn.nclasses; r++) {
            if (!Double.isInfinite(sn.njobs.get(r))) {
                allOpen = false;
                break;
            }
        }
        if (sn.nclasses == 1 && allOpen && SnHasBurstyArrival.snHasBurstyArrival(sn)) {
            return "rqna";
        }
        if (sn.nclasses == 1 && allOpen && sn.nstations == 2) {
            // A single-class open station customers ABANDON is a different
            // model, not a correction to a G/G/k one: the resolution has to
            // happen here as well as in the analyzer, because the feature gate
            // runs on the resolved name and Reneging is admitted for these two
            // methods only.
            int qi = queueStationIndex(sn);
            if (qi >= 0) {
                SnPatienceHandles.Handles hpat = SnPatienceHandles.snPatienceHandles(sn, qi, 0);
                if (hpat != null) {
                    return hpat.isExponential ? "erlanga" : "mgisrgi";
                }
            }
        }
        return method;
    }

    /**
     * Station index of the single Queue node, or -1 when the model has none.
     *
     * @param sn the network struct
     * @return the station index of the Queue node
     */
    private static int queueStationIndex(NetworkStruct sn) {
        if (sn == null || sn.nodetype == null) {
            return -1;
        }
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Queue) {
                return (int) sn.nodeToStation.get(i);
            }
        }
        return -1;
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
        String capacityReason = finiteCapacityReason(this.model, sn, this.options.method);
        if (capacityReason != null) {
            throw new RuntimeException(capacityReason);
        }

        if (sn.immfeed != null && sn.immfeed.elementSum() > 0) {
            line_warning(mfilename(new Object(){}), "SolverMVA does not handle immediate feedback (immfeed); the solver will treat self-loops as class-switching with re-queueing.");
        }

        line_debug(options.verbose, String.format("MVA solver starting: method=%s, nstations=%d, nclasses=%d",
            options.method, sn.nstations, sn.nclasses));

        String origMethod = options.method;

        // MODEL TRANSFORMATION, opt-in through options.config.transform. The
        // strategy rewrites the model into subproblems, TransformSolve solves
        // each with an instance of THIS solver and maps the metrics back, so a
        // strategy written once serves MVA as well as CTMC. Mirrors the branch
        // MATLAB puts in the shared runAnalyzerPreamble.
        if (jline.solvers.tr.TransformSolve.isRequested(options)) {
            runTransformAnalyzer();
            return;
        }

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
        // setAvgResults copies only the metrics the base container declares, so the
        // MVA-specific fields have to be carried over explicitly: this.result is a
        // DIFFERENT MVAResult from the analyzer's. Without this the convergence flag
        // never reaches an API caller (the count alone cannot decide it on the
        // amvald route, where it aggregates the nested sweeps), and the normalizing
        // constant stays at the field default 0.0, which getProbNormConstAggr
        // accepts as computed because its guard is !isNaN.
        if (this.result instanceof MVAResult) {
            MVAResult mvaResult = (MVAResult) this.result;
            mvaResult.converged = ret.converged;
            mvaResult.logNormConstAggr = ret.logNormConstAggr;
        }
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
     * Method-level gate. The base gate answers from the feature set; method
     * "exact" additionally requires a product-form solution, which has no
     * registry feature name and so cannot live in getMethodFeatureSet.
     *
     * @param method the concrete method name
     * @return empty string if supported, else the offending reason
     */
    @Override
    public String supportsModelMethod(String method) {
        String reason = super.supportsModelMethod(method);
        if (!reason.isEmpty()) {
            return reason;
        }
        if (!(this.model instanceof Network)) {
            return "";
        }
        reason = exactnessReason((Network) this.model, method);
        if (!reason.isEmpty()) {
            return reason;
        }
        // A binding buffer is a per-model rule the registry can name only in
        // part: methodFeatureSet says WHICH names honour one, this says whether
        // THIS model has one the named route can carry. Same predicate
        // runAnalyzer stops on, so the report and the run cannot disagree.
        String capReason = finiteCapacityReason((Network) this.model,
                ((Network) this.model).getStruct(false), method);
        if (capReason != null) {
            return capReason;
        }
        // Product form, a class count and a server count have no registry feature
        // name, so these three rules cannot live in getMethodFeatureSet. Each is
        // the SAME predicate the analyzer raises on, so a row the report offers is
        // a row that runs.
        NetworkStruct s = ((Network) this.model).getStruct(false);
        reason = closedPopulationReason(s, method);
        if (!reason.isEmpty()) {
            return reason;
        }
        reason = singleClassOpenReason(s, method);
        if (!reason.isEmpty()) {
            return reason;
        }
        reason = mvacReason(s, method);
        if (!reason.isEmpty()) {
            return reason;
        }
        if ("mapqn".equals(baseMethod(method))) {
            reason = Solver_mapqn.mapqnReason(s);
            if (!reason.isEmpty()) {
                return reason;
            }
        }
        return schmidtExtReasonForStruct(s, method);
    }

    /**
     * {@link #schmidtExtReason} asked about a struct: it rebuilds the population
     * vector and per-station discipline that Solver_amva's schmidt-ext arm hands
     * the kernel, so the report and the run ask the same question of the same
     * numbers.
     *
     * <p>That arm recurs on the CHAIN populations under class switching and on the
     * class ones otherwise, which is the vector the kernel is handed; this asks
     * about the same one. All four codebases now agree on that basis.</p>
     *
     * @param sn     the network struct
     * @param method the requested method name
     * @return empty string if the method may run, else the offending reason
     */
    public static String schmidtExtReasonForStruct(NetworkStruct sn, String method) {
        if (sn == null || !"schmidt-ext".equals(baseMethod(method))) {
            return "";
        }
        List<Boolean> fcfs = new ArrayList<Boolean>();
        for (int ist = 0; ist < sn.nstations; ist++) {
            fcfs.add(sn.sched.get(sn.stations.get(ist)) == SchedStrategy.FCFS);
        }
        Matrix pop = sn.njobs;
        if (jline.api.sn.SnHasClassSwitching.snHasClassSwitching(sn)) {
            pop = jline.api.sn.SnGetProductFormChainParams.snGetProductFormChainParams(sn).N;
        }
        return schmidtExtReason(pop, fcfs, method);
    }

    /**
     * Refusal reason for the extended Schmidt method, empty when it may run.
     *
     * Schmidt's EXTENSION over plain Schmidt is an alpha correction applied at an
     * FCFS station, and the correction is computed from the network with ONE
     * class-r customer TAGGED, that is at population N - 1_r. A class holding no
     * customer has none to tag: the sub-problem is formed at a negative
     * population, whose state lattice prod(N+1) collapses to zero and the
     * recursion indexes an empty array. Plain "schmidt" forms no such
     * sub-problem, which is why the requirement is the -ext arm's alone.
     *
     * <p>THE TEST IS STATED AT THE FCFS STATION AND NOT AT A CLASS-DEPENDENT ONE,
     * because the four kernels differ on when they form the correction: MATLAB,
     * the C++ port and native python form it only where the station's demands
     * differ by class, while Pfqn_schmidt_amva forms it at EVERY FCFS station.
     * Stating the union is what keeps one rule safe for all four; the case it
     * costs -- an FCFS station whose demands are identical across classes, one of
     * them empty -- is one where the extension reduces to plain "schmidt", which
     * stays offered.</p>
     *
     * @param N      the population vector the kernel recurs on
     * @param fcfs   whether each station row is served FCFS
     * @param method the requested method name
     * @return empty string if the method may run, else the offending reason
     */
    public static String schmidtExtReason(Matrix N, List<Boolean> fcfs, String method) {
        if (!"schmidt-ext".equals(baseMethod(method)) || N == null || fcfs == null) {
            return "";
        }
        boolean anyFcfs = false;
        for (Boolean b : fcfs) {
            if (Boolean.TRUE.equals(b)) {
                anyFcfs = true;
                break;
            }
        }
        if (!anyFcfs) {
            return "";
        }
        for (int r = 0; r < N.getNumElements(); r++) {
            if (Double.isFinite(N.get(r)) && N.get(r) < 1.0) {
                return "the 'schmidt-ext' method corrects an FCFS station from the network "
                        + "with one customer of that class tagged, so it needs every class to "
                        + "hold at least one customer; class " + (r + 1) + " holds none. Use "
                        + "'schmidt' for the uncorrected recursion.";
            }
        }
        return "";
    }

    /**
     * The AMVA algorithms whose recursion is over a CLOSED population vector.
     *
     * Each approximates the arrival-instant queue length E[Q(N-1_r)] from E[Q(N)]
     * and is handed (L, N, Z) alone, with no arrival rate and no rate function, so
     * an open chain gives it nothing to recur on. Canonical spellings only; an
     * "amva." prefix is stripped before the list is consulted, exactly as
     * Solver_amva strips it before it selects an algorithm. ONE list, read by
     * methodFeatureSet (which drops OpenClass for these) and by Solver_amva (which
     * refuses them by name), so the gate and the run cannot drift apart.
     */
    public static final List<String> CLOSED_POPULATION_METHODS =
            Collections.unmodifiableList(Arrays.asList(
                    "bs", "aql", "qsa", "sqni", "tay", "scat", "lcp", "chow",
                    "pamb", "pami", "pamt", "clust", "dmlin", "ab", "schmidt", "schmidt-ext"));

    /**
     * Scheduling feature names OUTSIDE the BCMP set {INF, PS, FCFS, SIRO, LCFS-PR}
     * that the base MVA envelope declares. The chain algorithms that walk the
     * stations one by one (the summation method, MVAC, QNA) accept the BCMP set and
     * refuse the rest, so each drops these from its own envelope.
     */
    public static final String[] NON_BCMP_SCHED_FEATURES = new String[]{
            "SchedStrategy_HOL", "SchedStrategy_DPS", "SchedStrategy_FCFSPRPRIO",
            "SchedStrategy_LCFS", "SchedStrategy_POLLING", "SchedStrategy_SJF",
            "SchedStrategy_SRPT", "SchedStrategy_PSJF", "SchedStrategy_FB",
            "SchedStrategy_LRPT", "SchedStrategy_SETF",
            "SchedStrategy_OI", "SchedStrategy_PAS"};

    /**
     * The method name with any leading "amva." alias stripped.
     *
     * @param method the requested method name
     * @return the canonical spelling
     */
    public static String baseMethod(String method) {
        if (method == null) {
            return "";
        }
        String m = method.toLowerCase();
        return m.startsWith("amva.") ? m.substring(5) : m;
    }

    /**
     * Whether the name, alias stripped, is one of the closed-population algorithms.
     *
     * @param method the requested method name
     * @return true when it belongs to the family
     */
    public static boolean isClosedPopulationMethod(String method) {
        return CLOSED_POPULATION_METHODS.contains(baseMethod(method));
    }

    /**
     * Refusal reason for the closed-population AMVA family, empty when it may run.
     *
     * Open chains are also expressed in the registry (methodFeatureSet drops
     * OpenClass for these), which is what keeps them off the report; they are
     * repeated here because the analyzer must refuse by name with a sentence
     * rather than fall through and answer under a method nobody asked for. Strict
     * product form has no registry name at all, so this is its only home.
     *
     * @param sn     the network struct
     * @param method the requested method name
     * @return empty string if the method may run, else the offending reason
     */
    public static String closedPopulationReason(NetworkStruct sn, String method) {
        if (sn == null || !isClosedPopulationMethod(method)) {
            return "";
        }
        String base = baseMethod(method);
        if (snHasOpenClasses(sn)) {
            return "the '" + base + "' method approximates the arrival-instant queue length as a "
                    + "function of the closed population vector N, so it is defined for closed "
                    + "models only; use 'default', 'lin', 'qd' or 'qna' for a model with open classes";
        }
        // ab, schmidt and schmidt-ext ARE the class-dependent FCFS algorithms, so
        // heterogeneous FCFS service means are their subject matter rather than a
        // disqualification.
        boolean checkMeans = !("ab".equals(base) || "schmidt".equals(base) || "schmidt-ext".equals(base));
        if (snHasProductFormNotHetFCFS(sn, checkMeans)) {
            return "";
        }
        return "the '" + base + "' method is defined for strict product-form, load-independent "
                + "models; use 'default', 'lin' or 'qd' for this model";
    }

    /**
     * Refusal reason for RQNA and RQT, empty when the method may run.
     *
     * Both build one uncertainty set per flow out of the first two moments of a
     * SINGLE stream, so a multiclass model has no counterpart in their equations.
     * No registry feature names a class count, so the rule is structural.
     *
     * @param sn     the network struct
     * @param method the requested method name
     * @return empty string if the method may run, else the offending reason
     */
    public static String singleClassOpenReason(NetworkStruct sn, String method) {
        String base = baseMethod(method);
        if (sn == null || (!"rqna".equals(base) && !"rqt".equals(base))) {
            return "";
        }
        String label = base.toUpperCase();
        if (sn.nodetype != null) {
            for (NodeType nt : sn.nodetype) {
                if (nt == NodeType.Fork || nt == NodeType.Join) {
                    return label + " decomposes an open network into GI/G/1 queues and has no "
                            + "synchronisation term; a Join carries no service process for its "
                            + "index of dispersion to be read from. Use SolverMVA's 'default' "
                            + "method for a fork-join model.";
                }
            }
        }
        if (sn.nclasses == 1) {
            return "";
        }
        return label + " supports single-class open networks only. Use the 'qna' "
                + "method for multiclass models.";
    }

    /**
     * Refusal reason for MVAC, empty when the method may run.
     *
     * MVAC (Conway-de Souza e Silva-Lavenberg) is the exact chain recursion over
     * single-server fixed-rate (SSFR) queues and infinite-server centres of a
     * product-form network. Neither the server count nor product form has a
     * registry feature name, so both are structural; the scheduling restriction IS
     * nameable and lives in methodFeatureSet.
     *
     * @param sn     the network struct
     * @param method the requested method name
     * @return empty string if the method may run, else the offending reason
     */
    public static String mvacReason(NetworkStruct sn, String method) {
        if (sn == null || !"mvac".equals(baseMethod(method))) {
            return "";
        }
        if (!snHasProductForm(sn)) {
            return "MVAC requires a product-form model.";
        }
        int nq = 0;
        for (int ist = 0; ist < sn.nstations; ist++) {
            SchedStrategy ss = sn.sched.get(sn.stations.get(ist));
            if (ss == SchedStrategy.INF || ss == SchedStrategy.EXT) {
                continue;
            }
            nq++;
            // An infinite count is refused here too, as the reference does: a
            // station scheduled FCFS with infinitely many servers is not an IS
            // centre to MVAC, and infSET is built from the DISCIPLINE, not the count.
            double m = sn.nservers.get(ist);
            if (m != 1.0) {
                return "MVAC supports single-server (SSFR) queues only; use method 'exact' for "
                        + "multiserver stations.";
            }
        }
        if (nq == 0) {
            // Pfqn_mvac recurs on the queueing centres and is handed a zero-row
            // demand matrix without one. MATLAB and native python answer the
            // degenerate pure-delay network from the delay closed form instead of
            // refusing it; the JAR and the C++ port refuse, and this gate says so
            // rather than letting the report offer a row the analyzer throws on.
            return "MVAC recurs on the queueing centers and needs at least one; this model has "
                    + "only delay stations.";
        }
        return "";
    }

    /**
     * Refusal reason for QNA's station update, empty when every station has an arm.
     *
     * The update has an arm for INF, PS and FCFS and none for any other discipline,
     * so a SIRO, LCFS, LCFS-PR, HOL or priority station used to leave its whole row
     * of Q, U, R and T at zero and the table was returned as a solution. The
     * registry expresses this as well (methodFeatureSet drops the disciplines from
     * QNA's envelope); this is the analyzer's half of it.
     *
     * @param sn the network struct
     * @return empty string if QNA may run, else the offending reason
     */
    public static String qnaSchedulingReason(NetworkStruct sn) {
        if (sn == null) {
            return "";
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            Station st = sn.stations.get(ist);
            if (sn.nodetype != null && sn.nodetype.get((int) sn.stationToNode.get(ist)) == NodeType.Join) {
                continue; // a Join station carries no service and is skipped below
            }
            SchedStrategy ss = sn.sched.get(st);
            if (ss == SchedStrategy.EXT || ss == SchedStrategy.INF
                    || ss == SchedStrategy.PS || ss == SchedStrategy.FCFS) {
                continue;
            }
            return "QNA decomposes every station as a GI/G/m centre and has no arm for "
                    + SchedStrategy.toText(ss) + " scheduling. Use the 'default' or 'lin' methods.";
        }
        return "";
    }

    /**
     * Product-form precondition of method "exact", the same rule MVARunner
     * enforces at solve time. Order-independent and pass-and-swap stations are
     * exempt: Solver_mva_oi_analyzer is exact for them regardless of the
     * product-form test.
     *
     * @param model  the network model
     * @param method the concrete method name
     * @return empty string if the method may run, else the offending reason
     */
    public static String exactnessReason(Network model, String method) {
        if (!"exact".equals(method) || model.hasProductFormSolution()) {
            return "";
        }
        NetworkStruct s = model.getStruct(false);
        if (s.sched != null) {
            for (SchedStrategy sched : s.sched.values()) {
                if (sched == SchedStrategy.OI || sched == SchedStrategy.PAS) {
                    return "";
                }
            }
        }
        return "method 'exact' requires a product-form solution; use 'mva' for "
                + "the approximation based on the exact MVA algorithm";
    }

    /**
     * MVA-specific finite-capacity gate, asked of the METHOD as well as of the
     * model. Mirrors MATLAB {@code SolverMVA.supportsFiniteCapacity}.
     *
     * <p>A Blocking-After-Service model is exempt where the finite buffers ARE
     * honoured, which is Solver_sqd alone: "sqd" names it and the default ladder
     * of Solver_mva_analyzer selects it. Every other name (the exact recursion
     * under "mva", every AMVA kernel, sum, mvac) solves the buffers away, and so
     * does "default" on a load-, class- or joint-dependent model, whose ladder
     * (Solver_mvald_analyzer) has no sqd arm. Everything else defers to the
     * shared product-form gate.
     *
     * <p>THE SECOND EXEMPTION is a single-station M/M/1/K with tail drop, which
     * the moment-based (MacGregor Smith) {@code Qsys_mg1k_loss_mgs} branch of
     * Solver_mva_qsys_analyzer answers, MVARunner reaching it by SHAPE under
     * every name. That branch is exact at scv=1 only, so "exact" is NOT exempted
     * and must refuse; every other name is.
     *
     * @param model  - the network model
     * @param sn     - the network structure
     * @param method - the concrete method name, "default" when unknown
     * @return null if the model has no binding capacity, otherwise the reason
     */
    public static String finiteCapacityReason(Network model, NetworkStruct sn, String method) {
        String m = (method == null || method.isEmpty()) ? "default" : baseMethod(method);
        boolean noScaling = (sn.lldscaling == null || sn.lldscaling.getNumElements() == 0)
                && (sn.cdscaling == null || sn.cdscaling.isEmpty())
                && (sn.jdscaling == null || sn.jdscaling.isEmpty());
        if (Solver_mva_analyzer.isBasModel(sn) && noScaling
                && ("default".equals(m) || "sqd".equals(m))) {
            return null;
        }
        // Single-station M/M/1/K with tail drop is answered by the moment-based
        // (MacGregor Smith) Qsys_mg1k_loss_mgs branch of
        // Solver_mva_qsys_analyzer, exact only at scv=1. It is an approximation
        // in general, so "exact" is NOT exempted (it must refuse); every other
        // method is.
        if (!"exact".equals(m) && SnIsMm1kLoss.snIsMm1kLoss(sn)) {
            return null;
        }
        return NetworkSolver.bindingCapacityReason(model, sn, "SolverMVA");
    }

    /**
     * The solver-level form of {@link #finiteCapacityReason(Network,
     * NetworkStruct, String)}, answered as the default method, as MATLAB
     * answers it when supportsFiniteCapacity is called with no method.
     *
     * @param model - the network model
     * @param sn    - the network structure
     * @return null if the model has no binding capacity, otherwise the reason
     */
    public static String finiteCapacityReason(Network model, NetworkStruct sn) {
        return finiteCapacityReason(model, sn, "default");
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
        // Re-read the struct AFTER the analysis: the one captured at construction
        // predates the model's default initialization, and it is getAvgQLen that
        // triggers that initialization, so its state rows are otherwise not the
        // ones every metric above was computed for.
        this.sn = this.model.getStruct(false);
        Matrix N = this.sn.njobs;

        if (N.isFinite()) {
            switch (this.options.method) {
                case "exact":
                    throw new RuntimeException("Exact marginal state probabilities not available yet in SolverMVA.");
                default:
                    Matrix state = statefulState(ist);

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
            // Mixed or open model: the open classes' product form at this
            // station, times the closed classes' binomials.
            Matrix stateOpen = statefulState(ist);
            jline.lang.state.State.StateMarginalStatistics margStatsOpen =
                    ToMarginal.toMarginal(this.sn, ist, stateOpen, null, null, null, null, null);
            Matrix nirOpen = margStatsOpen.nir;
            double logPnir = snOpenProbTerms(this.sn, Q, this.getAvgUtil(), nirOpen, ist);
            for (int r = 0; r < this.sn.nclasses; r++) {
                if (Double.isInfinite(N.get(r))) {
                    continue;
                }
                int nirVal = (int) nirOpen.get(0, r);
                int NVal = (int) N.get(r);
                double QVal = Q.get(ist, r);
                logPnir += logBinomial(NVal, nirVal);
                if (QVal > 0 && NVal > 0) {
                    logPnir += nirVal * Math.log(QVal / NVal);
                    logPnir += (NVal - nirVal) * Math.log(1.0 - QVal / NVal);
                }
            }
            ProbabilityResult result = new ProbabilityResult(Math.exp(logPnir));
            result.logNormalizingConstant = logPnir;
            return result;
        }
    }

    /**
     * State row of the stateful node backing a station.
     *
     * @param ist station index (0-based)
     * @return the station's state row
     */
    private Matrix statefulState(int ist) {
        int statefulIndex = (int) this.sn.stationToStateful.get(ist);
        // Key the map by the node itself. Walking entrySet() and counting to
        // statefulIndex assumes the map iterates in stateful order, which a
        // HashMap does not: it returned another station's state row.
        Matrix state = this.sn.state.get(this.model.getStatefulNodes().get(statefulIndex));
        if (state == null) {
            throw new RuntimeException("Could not find state for station " + ist);
        }
        return state;
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
        // Re-read the struct AFTER the analysis: the one captured at construction
        // predates the model's default initialization, and it is getAvgQLen that
        // triggers that initialization, so its state rows are otherwise not the
        // ones every metric above was computed for.
        this.sn = this.model.getStruct(false);
        Matrix N = this.sn.njobs;

        if (N.isFinite()) {
            switch (this.options.method) {
                case "exact":
                    throw new RuntimeException("Exact joint state probabilities not available yet in SolverMVA.");
                default:
                    // Binomial approximation with mean fitted to queue-lengths
                    // Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997
                    double logPn = 0.0;
                    for (int r = 0; r < N.length(); r++) {
                        logPn += factln((int) N.get(r));
                    }

                    for (int ist = 0; ist < this.sn.nstations; ist++) {
                        Matrix stateIst = statefulState(ist);

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
            // Mixed or open model: the closed classes' multinomial over the
            // stations, plus each station's open-class product form.
            Matrix U = this.getAvgUtil();
            double logPn = 0.0;
            for (int r = 0; r < N.length(); r++) {
                if (!Double.isInfinite(N.get(r))) {
                    logPn += factln((int) N.get(r));
                }
            }
            for (int ist = 0; ist < this.sn.nstations; ist++) {
                jline.lang.state.State.StateMarginalStatistics margStats =
                        ToMarginal.toMarginal(this.sn, ist, statefulState(ist), null, null, null, null, null);
                Matrix nir = margStats.nir;
                logPn += snOpenProbTerms(this.sn, Q, U, nir, ist);
                for (int r = 0; r < this.sn.nclasses; r++) {
                    if (Double.isInfinite(N.get(r))) {
                        continue;
                    }
                    int nirVal = (int) nir.get(0, r);
                    logPn -= factln(nirVal);
                    if (Q.get(ist, r) > 0 && N.get(r) > 0) {
                        logPn += nirVal * Math.log(Q.get(ist, r) / N.get(r));
                    }
                }
            }
            ProbabilityResult result = new ProbabilityResult(Math.exp(logPn));
            result.logNormalizingConstant = logPn;
            return result;
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
                "aql", "amva.aql",
                "qsa", "amva.qsa",
                "sqni",
                "qd", "amva.qd",
                "qli", "amva.qli",
                "fli", "amva.fli",
                "lin", "egflin", "gflin", "amva.lin",
                "schmidt", "amva.schmidt", "schmidt-ext", "amva.schmidt-ext",
                "tay", "amva.tay",
                "scat", "amva.scat",
                "ab", "amva.ab",
                "lcp", "amva.lcp",
                "chow", "amva.chow",
                "pamb", "amva.pamb",
                "pami", "amva.pami",
                "pamt", "amva.pamt",
                "clust", "amva.clust",
                "dmlin", "amva.dmlin"
        ));

        // SQNI (pfqn_sqni) is a closed form for one queueing station with a
        // delay; listing it elsewhere named a method that cannot run.
        int nInfSqni = 0;
        for (Station st : this.sn.stations) {
            if (this.sn.sched.get(st) == SchedStrategy.INF) nInfSqni++;
        }
        if (this.sn.nstations != 2 || nInfSqni != 1) {
            allMethods.remove("sqni");
        }

        // AQL (pfqn_aql), QSA (pfqn_qsa) and Tay (pfqn_tay) reject multiserver stations at solve
        // time, so they are only advertised for single-server models.
        if (hasMultiServerStations()) {
            allMethods.remove("aql");
            allMethods.remove("amva.aql");
            allMethods.remove("qsa");
            allMethods.remove("amva.qsa");
            allMethods.remove("tay");
            allMethods.remove("amva.tay");
        }

        // see _kb/06-solver-catalog.md for rationale
        if (snIsOpenModel(this.sn)) {
            allMethods.addAll(4, Arrays.asList("qna", "rqna", "rqt"));
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

        // SJN (shortest-job-next, pfqn_mvasjn / pfqn_amvasjn): the conditional
        // waiting time equation is a population recursion, so the family runs on
        // a CLOSED model with an SJF station and nowhere else -- MVARunner
        // rejects an open one by name. Advertised only there, for the reason
        // 'sqni' is gated above: a name on the list is a name a caller is
        // invited to ask for. Withholding them entirely was worse than either,
        // because `checkDeclaredMethod` gates on this list, so `sjn.mva` was
        // dispatched by MVARunner and refused before it ever got there.
        if (!anyOpen && snHasSjn(this.sn)) {
            allMethods.addAll(Arrays.asList("sjn.mva", "sjn.amva"));
        }

        // SQD (Smith Queue Decomposition) is only valid for closed single-chain
        // Blocking-After-Service networks; Solver_sqd returns empty results on
        // anything else, so listing it unconditionally named a method that
        // cannot run. Mirrors SolverMVA.m and the C++ runner.
        if (Solver_mva_analyzer.isBasModel(this.sn)) {
            allMethods.add("sqd");
        }

        // see _kb/06-solver-catalog.md for rationale
        if (!snIsOpenModel(this.sn) && !snHasClassdepRouting(this.sn)) {
            allMethods.addAll(Arrays.asList("marie", "amva.marie"));
        }

        // amva.mapqn: the horizontal-cut MVA for one exponential delay and one FCFS
        // MAP queue (Mapqn_amva); offered only on that shape, which
        // Solver_mapqn.mapqnReason judges for the list, the report and the run
        if (Solver_mapqn.mapqnReason(this.sn).isEmpty()) {
            allMethods.add("amva.mapqn");
        }

        // priomva: preemptive-resume priority arm (Chandy-Lakshmi [ChaL83]),
        // offered only when a station actually uses FCFSPRPRIO. The arm itself
        // lives in Solver_amvald's forward step. Mirrors SolverMVA.m.
        if (hasPrsPrioStation()) {
            allMethods.addAll(Arrays.asList("priomva", "amva.priomva"));
        }

        // Add queueing system methods for open single-class 2-station models
        if (snIsOpenModel(this.sn) && this.sn.nstations == 2 && this.sn.nclasses == 1) {
            allMethods.addAll(Arrays.asList(
                    "mm1", "mmk", "mg1", "mgi1", "gm1", "gig1", "gim1", "gig1.kingman",
                    "gigk", "gigk.kingman_approx",
                    "gig1.gelenbe", "gig1.heyman", "gig1.kimura", "gig1.allen",
                    "gig1.kobayashi", "gig1.klb", "gig1.marchal",
                    "gigk.whitt", "qed", "gig1.extremal", "gigk.diffusion"
            ));
            // The two abandonment methods are listed only when the station
            // actually reneges: they have nothing to say about a queue nobody
            // walks away from, and listing them there would name a method that
            // cannot run.
            int qiList = queueStationIndex(this.sn);
            if (qiList >= 0 && SnPatienceHandles.snPatienceHandles(this.sn, qiList, 0) != null) {
                allMethods.addAll(Arrays.asList("erlanga", "mgisrgi"));
            }
        }

        return allMethods.toArray(new String[0]);
    }

    /**
     * Check whether any station schedules by preemptive-resume priority.
     *
     * @return true if some station is FCFSPRPRIO
     */
    private boolean hasPrsPrioStation() {
        if (this.sn == null || this.sn.sched == null) {
            return false;
        }
        for (SchedStrategy ss : this.sn.sched.values()) {
            if (ss == SchedStrategy.FCFSPRPRIO) {
                return true;
            }
        }
        return false;
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
            // Finite servers only. A delay station carries nservers = Inf, and
            // Inf > 1, so an unmasked test called EVERY model with a delay
            // multiserver and silently withheld aql/qsa/tay (and their amva.
            // aliases) from listValidMethods, six methods that run there. The
            // reference masks: sn_has_multi_server.m is
            // any(sn.nservers(isfinite(sn.nservers)) > 1), and the C++
            // has_multi_server() and Python both do the same.
            if (Double.isFinite(this.sn.nservers.get(ist)) && this.sn.nservers.get(ist) > 1) {
                return true;
            }
        }

        return false;
    }

    /**
     * Solves through a model transformation, with MVA as its own inner solve.
     *
     * <p>The inner options carry {@code transform='none'} and a raised depth, so
     * a transformed submodel cannot re-enter the driver.
     */
    private void runTransformAnalyzer() {
        long T0 = System.nanoTime();
        String methodName = jline.solvers.tr.TransformSolve.requested(options);
        SolverOptions sub = options.copy();
        sub.config.put("transform", methodName);
        jline.solvers.tr.TransformSolve.Result tr = jline.solvers.tr.TransformSolve.run(
                this.model, this.sn, sub, new jline.solvers.tr.TransformSolve.InnerSolve() {
                    @Override
                    public jline.solvers.tr.TransformSolve.Inner solve(jline.lang.Network submodel,
                                                                       SolverOptions opts) {
                        SolverMVA inner = new SolverMVA(submodel, opts);
                        try {
                            inner.runAnalyzer();
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                        return new jline.solvers.tr.TransformSolve.Inner(
                                inner.getAvgQLen(), inner.getAvgUtil(), inner.getAvgRespT(),
                                inner.getAvgTput(), inner.getAvgSysTput(), Double.NaN,
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
