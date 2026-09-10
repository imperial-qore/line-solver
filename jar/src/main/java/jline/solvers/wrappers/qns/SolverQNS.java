package jline.solvers.wrappers.qns;

import jline.GlobalConstants;
import static jline.GlobalConstants.Inf;
import static jline.api.sn.SnHasProductForm.snHasProductForm;
import static jline.api.sn.SnHasOpenClasses.snHasOpenClasses;
import static jline.io.InputOutput.line_debug;

import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkStruct;
import jline.io.QN2LQN;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.wrappers.lqns.SolverLQNS;
import jline.solvers.wrappers.qns.analyzers.Solver_qns_analyzer;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.solvers.AvgHandle;
import jline.util.matrix.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.List;

/**
 * SolverQNS class implements a queueing network solver that wraps the external qnsolver tool.
 * This solver provides various multiserver approximation methods for analyzing queueing networks.
 */
public class SolverQNS extends NetworkSolver {

    /**
     * Default constructor with network model
     */
    public SolverQNS(Network model) {
        super(model, "QNS");
        this.options = defaultOptions();
        this.result = new QNSResult();
    }

    /**
     * Constructor with network model and method string
     */
    public SolverQNS(Network model, String method) {
        super(model, "QNS");
        this.options = defaultOptions();
        this.options.method = method;
        this.result = new QNSResult();
    }

    /**
     * Constructor with network model and solver options
     */
    public SolverQNS(Network model, SolverOptions options) {
        super(model, "QNS");
        this.options = options;
        this.result = new QNSResult();
    }

    /**
     * Constructor with network model and variable arguments
     */
    public SolverQNS(Network model, Object... varargin) {
        super(model, "QNS");
        this.options = parseOptions(varargin);
        this.result = new QNSResult();
    }

    /**
     * Parse options from variable arguments
     */
    public static SolverOptions parseOptions(Object... varargin) {
        SolverOptions options = defaultOptions();

        for (int i = 0; i < varargin.length; i += 2) {
            if (i + 1 < varargin.length && varargin[i] instanceof String) {
                String key = (String) varargin[i];
                Object value = varargin[i + 1];

                switch (key.toLowerCase()) {
                    case "method":
                        if (value instanceof String) {
                            options.method = (String) value;
                        }
                        break;
                    case "multiserver":
                        if (value instanceof String) {
                            options.config.multiserver = (String) value;
                        }
                        break;
                    case "timespan":
                        if (value instanceof double[]) {
                            options.timespan = (double[]) value;
                        }
                        break;
                }
            }
        }

        return options;
    }

    /**
     * Get the default options for the QNS solver
     */
    public static SolverOptions defaultOptions() {
        SolverOptions options = new SolverOptions(SolverType.QNS);
        options.method = "default";
        options.config.multiserver = "default";
        options.timespan = new double[]{Inf, Inf};
        return options;
    }

    /**
     * Get the feature set supported by this solver
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();

        // Node types
        featSupported.setTrue("Source");
        featSupported.setTrue("Sink");
        featSupported.setTrue("Queue");
        featSupported.setTrue("Delay");
        featSupported.setTrue("DelayStation");
        featSupported.setTrue("Router");
        featSupported.setTrue("ClassSwitch");
        featSupported.setTrue("Fork");
        featSupported.setTrue("Join");
        featSupported.setTrue("Forker");
        featSupported.setTrue("Joiner");
        featSupported.setTrue("Logger");

        // Service distributions
        featSupported.setTrue("Exp");
        featSupported.setTrue("HyperExp");
        featSupported.setTrue("Coxian");
        featSupported.setTrue("Cox2");
        featSupported.setTrue("APH");
        featSupported.setTrue("Erlang");
        featSupported.setTrue("Det");
        featSupported.setTrue("Gamma");
        featSupported.setTrue("Lognormal");
        featSupported.setTrue("MAP");
        featSupported.setTrue("MMPP2");
        featSupported.setTrue("Normal");
        featSupported.setTrue("PH");
        featSupported.setTrue("Pareto");
        featSupported.setTrue("Weibull");
        featSupported.setTrue("Uniform");
        featSupported.setTrue("Trace");
        featSupported.setTrue("Replayer");

        // Server types
        featSupported.setTrue("InfiniteServer");
        featSupported.setTrue("SharedServer");
        featSupported.setTrue("Server");
        featSupported.setTrue("Buffer");
        featSupported.setTrue("Dispatcher");
        featSupported.setTrue("JobSink");
        featSupported.setTrue("RandomSource");
        featSupported.setTrue("ServiceTunnel");
        featSupported.setTrue("LogTunnel");
        featSupported.setTrue("StatelessClassSwitcher");

        // Petri Net elements
        featSupported.setTrue("Linkage");
        featSupported.setTrue("Enabling");
        featSupported.setTrue("Timing");
        featSupported.setTrue("Firing");
        featSupported.setTrue("Storage");
        featSupported.setTrue("Place");
        featSupported.setTrue("Transition");

        // Scheduling strategies
        featSupported.setTrue("SchedStrategy_INF");
        featSupported.setTrue("SchedStrategy_PS");
        featSupported.setTrue("SchedStrategy_DPS");
        featSupported.setTrue("SchedStrategy_FCFS");
        featSupported.setTrue("SchedStrategy_GPS");
        featSupported.setTrue("SchedStrategy_SIRO");
        featSupported.setTrue("SchedStrategy_HOL");
        featSupported.setTrue("SchedStrategy_LCFS");
        featSupported.setTrue("SchedStrategy_LCFSPR");
        featSupported.setTrue("SchedStrategy_SEPT");
        featSupported.setTrue("SchedStrategy_LEPT");
        featSupported.setTrue("SchedStrategy_SJF");
        featSupported.setTrue("SchedStrategy_LJF");
        featSupported.setTrue("SchedStrategy_EXT");

        // Routing strategies
        featSupported.setTrue("RoutingStrategy_PROB");
        featSupported.setTrue("RoutingStrategy_RAND");
        featSupported.setTrue("RoutingStrategy_RROBIN");
        featSupported.setTrue("RoutingStrategy_WRROBIN");
        featSupported.setTrue("RoutingStrategy_SQ");

        // Job classes
        featSupported.setTrue("OpenClass");
        featSupported.setTrue("ClosedClass");

        // c-server stations: the JMVA document carries the count as an
        // <ldstation> and the LQN as a host multiplicity; "suri" and "schmidt"
        // refuse one on the qnsolver path, which stays structural (it is a
        // condition on options.config.multiserver, not on the model).
        // FiniteCapacity is NOT declared: neither document has a buffer, which
        // is what the binding-capacity gate in supportsModelMethod refuses.
        featSupported.setTrue("MultiServer");

        return featSupported;
    }

    /**
     * Check if the solver supports the given model
     */
    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * Structural finite-capacity gate.
     *
     * <p>NOTHING under the QNS tree reads sn.cap or sn.classcap -- the model is
     * written out for {@code qnsolver}, whose MVA-family algorithms have no
     * representation of a finite buffer -- so a capped station was solved as an
     * unbounded one and the table reported the unconstrained answer under this
     * solver's name. There is no registry feature name for plain capacity, hence
     * the structural test; SolverMVA, SolverNC, SolverAG and SolverFluid gate
     * the same way through the same helper.
     *
     * <p>Without it SolverAUTO.listValidMethods offered all eight "qns" method names on
     * the BAS-blocking model of cqn_bas_blocking.
     *
     * @param method the concrete method name
     * @return empty string if supported, else the offending reason
     */
    @Override
    public String supportsModelMethod(String method) {
        // A BINDING FINITE BUFFER FIRST, ahead of the base gate. Nothing under
        // the QNS tree reads sn.cap or sn.classcap -- neither the JMVA document
        // qnsolver reads nor the LQN QN2LQN writes has a buffer -- and this
        // refusal names the station, the cap and the way out, where the feature
        // envelope can only say "(feature: FiniteCapacity)". Once FiniteCapacity
        // became a registry name on 2026-09-05 the base gate started answering
        // first and the useful sentence became unreachable.
        if (this.model != null) {
            String capReason = NetworkSolver.bindingCapacityReason(this.model,
                    this.model.getStruct(false), "SolverQNS");
            if (capReason != null && !capReason.isEmpty()) {
                return capReason;
            }
        }
        String reason = super.supportsModelMethod(method);
        if (reason != null && !reason.isEmpty()) {
            return reason;
        }
        if (this.model == null) {
            return "";
        }
        NetworkStruct snLocal = this.model.getStruct(false);
        String immfeed = qnsImmfeedRefusal(snLocal);
        if (!immfeed.isEmpty()) {
            return immfeed;
        }
        if (this.model.hasProductFormSolution() || this.model.hasOpenClasses()) {
            String ms = qnsMultiserverRefusal(snLocal, method);
            if (!ms.isEmpty()) {
                return ms;
            }
            String jmva = jline.solvers.wrappers.jmt.SolverJMT.jmtMethodRefusal(snLocal, method, null);
            if (jmva != null && !jmva.isEmpty()) {
                return jmva;
            }
        }
        return "";
    }

    /**
     * Why SolverQNS cannot serve a model with immediate feedback, or "" when the
     * model has none.
     *
     * <p>Immediate feedback (sn.immfeed) keeps a self-looping job on its server
     * instead of re-queueing it, and neither path of SolverQNS can state that:
     * the JMVA document qnsolver reads carries a mean demand and a visit count
     * per chain, and the LQN QN2LQN writes turns the routing into OR-fork
     * precedences of pseudo-activities on the reference task, where a repeated
     * visit is a new call. Either would answer for re-queueing under this
     * solver's name.
     *
     * <p>ONE PREDICATE, TWO CALLERS: {@link #supportsModelMethod} (the gate,
     * hence model.help and SolverAUTO) and {@link #runAnalyzer} (the run, for a
     * caller with enableChecks off). SolverJMT keeps its own wording in
     * jmtMethodRefusal. Mirrors matlab/src/solvers/wrappers/QNS/qns_immfeed_refusal.m.
     *
     * @param sn the network struct
     * @return the refusal, or "" when the model carries no immediate feedback
     */
    public static String qnsImmfeedRefusal(NetworkStruct sn) {
        if (sn == null || sn.immfeed == null || sn.immfeed.isEmpty()) {
            return "";
        }
        if (sn.immfeed.elementSum() <= 0) {
            return "";
        }
        return "SolverQNS does not support immediate feedback (sn.immfeed): neither the JMVA "
                + "document qnsolver reads nor the LQN QN2LQN writes can keep a self-looping job "
                + "on its server. Use SolverCTMC or SolverSSA, whose state space carries the "
                + "self-loop.";
    }

    /**
     * Whether qnsolver's own -m switch offers this multiserver approximation.
     *
     * <p>THE RULE IS INSIDE THE MULTISERVER BRANCH, and that is not a detail.
     * Without a multiserver station the reference emits no -m at all and answers
     * under the caller's method name, so refusing "suri" there would refuse a
     * model this solver does solve.
     *
     * <p>"qnsolver -m" accepts conway, reiser, rolia and zhou. "suri" and
     * "schmidt" are LQNS approximations, reachable only on the
     * non-product-form closed SolverLQNS branch, and qnsolver has no flag for
     * either. Mirrors matlab/src/solvers/wrappers/QNS/qns_multiserver_refusal.m
     * and the C++ is_qnsolver_multiserver.
     *
     * @param sn     the network struct
     * @param method the requested method name
     * @return the refusal, or "" when the pair is served
     */
    public static String qnsMultiserverRefusal(NetworkStruct sn, String method) {
        if (sn == null || method == null || method.isEmpty()) {
            return "";
        }
        boolean multiserver = false;
        if (sn.nservers != null) {
            for (int i = 0; i < sn.nservers.length() && !multiserver; i++) {
                double c = sn.nservers.get(i);
                if (c > 1 && !Double.isInfinite(c)) {
                    multiserver = true;
                }
            }
        }
        if (!multiserver) {
            // No multiserver station, so no -m flag is emitted and every method
            // name is served by the plain invocation.
            return "";
        }
        String ms = method.toLowerCase();
        if (ms.equals("default") || ms.equals("conway") || ms.equals("reiser")
                || ms.equals("rolia") || ms.equals("zhou")) {
            return "";
        }
        return "SolverQNS: the multiserver approximation '" + ms + "' is one LQNS offers and "
                + "qnsolver does not: 'qnsolver -m' accepts conway, reiser, rolia and zhou only; "
                + "suri and schmidt are available only on the non-product-form closed SolverLQNS "
                + "branch.";
    }

    /**
     * List valid methods for this solver
     */
    public String[] listValidMethods() {
        return new String[]{"default", "conway", "rolia", "zhou", "suri", "reiser", "schmidt"};
    }

    /**
     * Run the analyzer for the QNS solver
     */
    @Override
    public void runAnalyzer() throws IllegalAccessException, ParserConfigurationException, IOException {
        long startTime = System.nanoTime();

        if (this.model == null) {
            throw new RuntimeException("Model is not provided");
        }

        if (this.options == null) {
            this.options = defaultOptions();
        }
        // Propagate solver verbose level to global
        GlobalConstants.Verbose = options.verbose;

        if (this.sn == null) {
            this.sn = this.model.getStruct(false);
        }
        // The gate's own sentence for a caller running with enableChecks off:
        // neither path can keep a self-looping job on its server.
        String immfeedReason = qnsImmfeedRefusal(this.sn);
        if (!immfeedReason.isEmpty()) {
            throw new RuntimeException(immfeedReason);
        }
        jline.io.InputOutput.line_ack(options.verbose, "QNS");
        line_debug(options.verbose, String.format("QNS solver starting: method=%s, multiserver=%s, nstations=%d, nclasses=%d",
            options.method, options.config.multiserver, sn.nstations, sn.nclasses));

        // Map method to multiserver config (matches MATLAB lines 29-44)
        String method = this.options.method;
        switch (method) {
            case "conway":
                this.options.config.multiserver = "conway";
                break;
            case "rolia":
                this.options.config.multiserver = "rolia";
                break;
            case "zhou":
                this.options.config.multiserver = "zhou";
                break;
            case "suri":
                this.options.config.multiserver = "suri";
                break;
            case "reiser":
                this.options.config.multiserver = "reiser";
                break;
            case "schmidt":
                this.options.config.multiserver = "schmidt";
                break;
            case "default":
                this.options.config.multiserver = "rolia";
                break;
        }

        boolean isProductForm = snHasProductForm(sn);
        boolean isOpen = snHasOpenClasses(sn);

        if (isProductForm || isOpen) {
            // Product-form or open: use qnsolver directly (QN2LQN does not support Source/Sink)
            if (!Solver_qns_analyzer.isQNSolverAvailable()) {
                throw new RuntimeException("QNS solver requires the external 'qnsolver' tool for product-form and open networks. "
                        + "Obtain it from its authors at http://www.sce.carleton.ca/rads/lqns/; LINE ships no copy and "
                        + "runs none from a container image.");
            }
            line_debug(options.verbose, "QNS: product-form or open model, using qnsolver directly");
            Solver_qns_analyzer analyzer = new Solver_qns_analyzer(this);
            QNSResult result = analyzer.runAnalyzer();

            this.setAvgResults(result.QN, result.UN, result.RN, result.TN,
                    result.AN, result.WN, result.CN, result.XN,
                    result.runtime, result.method, result.iter);
        } else {
            // Non-product-form closed: convert to LQN and solve via SolverLQNS
            line_debug(options.verbose, "QNS: non-product-form closed model, converting to LQN and solving via SolverLQNS");
            LayeredNetwork lqnmodel = QN2LQN.convert(this.model);
            SolverOptions lqnsoptions = SolverLQNS.defaultOptions();
            lqnsoptions.verbose = this.options.verbose;

            // Map multiserver method to LQNS options
            String actualMethod = method;
            switch (method) {
                case "conway":
                    lqnsoptions.config.multiserver = "conway";
                    break;
                case "rolia":
                    lqnsoptions.config.multiserver = "rolia";
                    break;
                case "zhou":
                    lqnsoptions.config.multiserver = "zhou";
                    break;
                case "suri":
                    lqnsoptions.config.multiserver = "suri";
                    break;
                case "reiser":
                    lqnsoptions.config.multiserver = "reiser";
                    break;
                case "schmidt":
                    lqnsoptions.config.multiserver = "schmidt";
                    break;
                case "default":
                    lqnsoptions.config.multiserver = "rolia";
                    actualMethod = "rolia";
                    break;
            }

            LayeredNetworkAvgTable avgTable = new SolverLQNS(lqnmodel, lqnsoptions).getAvgTable();
            LayeredNetworkStruct lqn = lqnmodel.getStruct();

            int M = this.sn.nstations;
            int K = this.sn.nclasses;
            Matrix QN = new Matrix(M, K);
            Matrix UN = new Matrix(M, K);
            Matrix RN = new Matrix(M, K);
            Matrix TN = new Matrix(M, K);
            Matrix WN = new Matrix(M, K);

            List<Double> qlenList = avgTable.getQLen();
            List<Double> utilList = avgTable.getUtil();
            List<Double> respTList = avgTable.getRespT();
            List<Double> residTList = avgTable.getResidT();
            List<Double> tputList = avgTable.getTput();

            // see _kb/12-interfaces-and-docs.md (Wrappers: JAR subprocess-bridge notes: LQNS entry-phase utilization fallback)
            int entryOffset = lqn.eshift + sn.nchains;

            for (int r = 0; r < K; r++) {
                for (int i = 0; i < M; i++) {
                    // MATLAB: t = lqn.ashift + r + (i-1)*nclasses (1-based)
                    // Java: ashift is 0-based count of hosts+tasks+entries, list is 0-based
                    int t = lqn.ashift + r + i * K;
                    int e = entryOffset + r + i * K;
                    if (t < qlenList.size()) {
                        QN.set(i, r, qlenList.get(t));
                        // Use entry-level Util/Tput when activity-level values are NaN
                        double util = utilList.get(t);
                        double tput = tputList.get(t);
                        if (Double.isNaN(util) && e < utilList.size()) {
                            util = utilList.get(e);
                        }
                        if (Double.isNaN(tput) && e < tputList.size()) {
                            tput = tputList.get(e);
                        }
                        // see _kb/12-interfaces-and-docs.md (Wrappers: JAR subprocess-bridge notes: LQNS entry-phase utilization fallback)
                        // lqns sums the utilization over the host's servers, a Network station reports it per server
                        double nservers = this.sn.nservers.get(i);
                        if (!Double.isInfinite(nservers) && nservers > 0) {
                            util = util / nservers;
                        }
                        UN.set(i, r, util);
                        RN.set(i, r, respTList.get(t));
                        WN.set(i, r, residTList.get(t));
                        TN.set(i, r, tput);
                    }
                }
            }

            // Compute arrival rates from throughputs
            Matrix AN = new Matrix(M, K);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    AN.set(i, r, TN.get(i, r));
                }
            }

            double runtime = (System.nanoTime() - startTime) / 1_000_000_000.0;

            // Handle default method naming
            if ("default".equals(method)) {
                actualMethod = "default/" + actualMethod;
            }

            int C = this.sn.nchains;
            Matrix CN = new Matrix(1, C);
            Matrix XN = new Matrix(1, C);
            this.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN,
                    runtime, actualMethod, 0);
        }
    }

    /**
     * Check if the solver is available: a native {@code qnsolver} binary is on
     * the PATH. qnsolver ships with LQNS, whose licence forbids redistribution,
     * so LINE never runs it from a container image; use
     * {@code run-tests.sh --lqns-docker} to test a containerised build.
     */
    public static boolean isAvailable() {
        return Solver_qns_analyzer.hasNativeQNSolver();
    }

    // Probability methods - QNS solver does not support detailed state probability analysis

    @Override
    public ProbabilityResult getProbNormConstAggr() {
        throw new RuntimeException("getProbNormConstAggr not supported by QNS solver");
    }

    @Override
    public ProbabilityResult getProb(int node, Matrix state) {
        throw new RuntimeException("getProb not supported by QNS solver");
    }

    @Override
    public ProbabilityResult getProb(int node) {
        throw new RuntimeException("getProb not supported by QNS solver");
    }

    @Override
    public ProbabilityResult getProbSys() {
        throw new RuntimeException("getProbSys not supported by QNS solver");
    }

    @Override
    public ProbabilityResult getProbAggr(int node, Matrix state_a) {
        throw new RuntimeException("getProbAggr not supported by QNS solver");
    }

    @Override
    public ProbabilityResult getProbAggr(int node) {
        throw new RuntimeException("getProbAggr not supported by QNS solver");
    }

    @Override
    public ProbabilityResult getProbSysAggr() {
        throw new RuntimeException("getProbSysAggr not supported by QNS solver");
    }

    // Sampling methods - QNS solver does not support sampling

    @Override
    public SampleResult sample(int node, int numEvents) {
        throw new RuntimeException("sample not supported by QNS solver");
    }

    @Override
    public SampleResult sampleAggr(int node, int numEvents) {
        throw new RuntimeException("sampleAggr not supported by QNS solver");
    }

    @Override
    public SampleResult sampleSys(int numEvents) {
        throw new RuntimeException("sampleSys not supported by QNS solver");
    }

    @Override
    public SampleResult sampleSysAggr(int numEvents) {
        throw new RuntimeException("sampleSysAggr not supported by QNS solver");
    }

    // Distribution methods are NOT overridden: the reference SolverQNS declares
    // none, so getCdfRespT falls through to the NetworkSolver exponential
    // fallback and the transient/passage getters to the base refusals, exactly
    // as MATLAB's inheritance resolves them.
}