package jline.solvers.qns;

import static jline.GlobalConstants.Inf;
import static jline.api.sn.SnHasProductFormKt.snHasProductForm;
import static jline.api.sn.SnHasOpenClassesKt.snHasOpenClasses;

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
import jline.solvers.lqns.SolverLQNS;
import jline.solvers.qns.analyzers.Solver_qns_analyzer;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.solvers.AvgHandle;
import jline.util.matrix.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
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
        featSupported.setTrue("RoutingStrategy_KCHOICES");

        // Job classes
        featSupported.setTrue("OpenClass");
        featSupported.setTrue("ClosedClass");

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

        if (this.sn == null) {
            this.sn = this.model.getStruct(false);
        }

        if (this.options == null) {
            this.options = defaultOptions();
        }

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

        if (this.model.hasProductFormSolution() || this.model.hasOpenClasses()) {
            // Product-form or open: use qnsolver directly
            Solver_qns_analyzer analyzer = new Solver_qns_analyzer(this);
            QNSResult result = analyzer.runAnalyzer();

            this.setAvgResults(result.QN, result.UN, result.RN, result.TN,
                    result.AN, result.WN, result.CN, result.XN,
                    result.runtime, result.method, result.iter);
        } else {
            // Non-product-form closed: convert to LQN and solve via SolverLQNS
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

            for (int r = 0; r < K; r++) {
                for (int i = 0; i < M; i++) {
                    // MATLAB: t = lqn.ashift + r + (i-1)*nclasses (1-based)
                    // Java: ashift is 0-based count of hosts+tasks+entries, list is 0-based
                    int t = lqn.ashift + r + i * K;
                    if (t < qlenList.size()) {
                        QN.set(i, r, qlenList.get(t));
                        if (!Double.isInfinite(this.sn.nservers.get(i))) {
                            UN.set(i, r, utilList.get(t) / this.sn.nservers.get(i));
                        } else {
                            UN.set(i, r, utilList.get(t));
                        }
                        RN.set(i, r, respTList.get(t));
                        WN.set(i, r, residTList.get(t));
                        TN.set(i, r, tputList.get(t));
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
     * Check if the solver is available
     * Checks for the qnsolver command which is the actual executable used
     */
    public static boolean isAvailable() {
        try {
            String os = System.getProperty("os.name").toLowerCase();
            String command = "qnsolver -h";
            Process process;

            if (os.contains("win")) {
                process = Runtime.getRuntime().exec(new String[]{"cmd", "/c", command});
            } else {
                process = Runtime.getRuntime().exec(new String[]{"sh", "-c", command});
            }

            int exitCode = process.waitFor();

            // Check if the command exists and can be executed
            // exitCode 0 or 1 (help shows) typically means the command exists
            // exitCode 127 on Linux or 1 on Windows typically means command not found
            BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
            String line;
            while ((line = errorReader.readLine()) != null) {
                String lowerLine = line.toLowerCase();
                if (lowerLine.contains("command not found") ||
                        lowerLine.contains("not recognized") ||
                        lowerLine.contains("no such file")) {
                    errorReader.close();
                    return false;
                }
            }
            errorReader.close();

            // Also check stdout for command existence
            BufferedReader outReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            boolean hasOutput = outReader.readLine() != null;
            outReader.close();

            // If we got any output or a successful exit code, the command exists
            return hasOutput || exitCode == 0;

        } catch (IOException e) {
            return false;
        } catch (InterruptedException e) {
            return false;
        }
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

    // Distribution methods - QNS solver does not support distribution analysis

    @Override
    public DistributionResult getCdfRespT(AvgHandle R) {
        throw new RuntimeException("getCdfRespT not supported by QNS solver");
    }

    @Override
    public DistributionResult getCdfRespT() {
        throw new RuntimeException("getCdfRespT not supported by QNS solver");
    }

    @Override
    public DistributionResult getTranCdfRespT(AvgHandle R) {
        throw new RuntimeException("getTranCdfRespT not supported by QNS solver");
    }

    @Override
    public DistributionResult getTranCdfRespT() {
        throw new RuntimeException("getTranCdfRespT not supported by QNS solver");
    }

    @Override
    public DistributionResult getCdfPassT(AvgHandle R) {
        throw new RuntimeException("getCdfPassT not supported by QNS solver");
    }

    @Override
    public DistributionResult getCdfPassT() {
        throw new RuntimeException("getCdfPassT not supported by QNS solver");
    }

    @Override
    public DistributionResult getTranCdfPassT(AvgHandle R) {
        throw new RuntimeException("getTranCdfPassT not supported by QNS solver");
    }

    @Override
    public DistributionResult getTranCdfPassT() {
        throw new RuntimeException("getTranCdfPassT not supported by QNS solver");
    }
}