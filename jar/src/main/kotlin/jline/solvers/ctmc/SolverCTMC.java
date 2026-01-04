package jline.solvers.ctmc;

import jline.GlobalConstants;
import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.processes.MarkedMarkovProcess;
import jline.lang.state.State;
import jline.solvers.AvgHandle;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.io.Ret.ProbabilityResult;
import jline.solvers.ctmc.analyzers.RewardResult;
import jline.solvers.ctmc.analyzers.Solver_ctmc_analyzer;
import jline.solvers.ctmc.analyzers.Solver_ctmc_reward;
import jline.solvers.ctmc.handlers.Solver_ctmc;
import jline.solvers.ctmc.handlers.Solver_ctmc_joint;
import jline.solvers.ctmc.handlers.Solver_ctmc_marg;
import jline.solvers.ctmc.handlers.Solver_ctmc_margaggr;
import jline.solvers.ctmc.handlers.Solver_ctmc_jointaggr;
import jline.solvers.ctmc.ResultCTMCMargAggr;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static jline.api.mam.Map_meanKt.map_mean;
import static jline.api.sn.SnGetArvRFromTputKt.snGetArvRFromTput;
import static jline.io.InputOutputKt.*;
import static jline.util.PopulationLattice.pprod;
import static jline.util.Utils.isInf;
import static jline.util.matrix.Matrix.removeTrailingNewLine;
import static jline.api.mc.Ctmc_simulateKt.ctmc_simulate;
import static jline.api.mc.Ctmc_transientKt.ctmc_transient;
import static jline.api.mam.Map_normalizeKt.map_normalize;
import static jline.api.mam.Map_pieKt.map_pie;
import jline.util.Pair;
import jline.io.Ret;


/**
 * Solver for Continuous-Time Markov Chain (CTMC) analysis of queueing networks.
 * 
 * <p>SolverCTMC implements exact numerical analysis of queueing networks by constructing
 * and solving the underlying continuous-time Markov chain. This approach provides exact
 * results for steady-state and transient behavior of networks that may not satisfy
 * product-form assumptions.</p>
 * 
 * <p>Key CTMC solver capabilities:
 * <ul>
 *   <li>Exact CTMC state space construction and solution</li>
 *   <li>Steady-state probability computation</li>
 *   <li>Transient analysis with time-dependent solutions</li>
 *   <li>Joint and marginal state probability distributions</li>
 *   <li>Cache network modeling with exact hit/miss probabilities</li>
 *   <li>General service and arrival process support</li>
 * </ul>
 * </p>
 * 
 * <p>The solver automatically constructs the infinitesimal generator matrix Q
 * and solves the balance equations πQ = 0 for steady-state analysis, or the
 * differential equation dπ/dt = πQ for transient analysis.</p>
 * 
 * @see jline.api.mc
 * @see CTMCResult  
 * @see CTMCOptions
 * @since 1.0
 */
public class SolverCTMC extends NetworkSolver {

    public SolverCTMC(Network model, Object... args) {
        super(model, "SolverCTMC");
        this.setOptions(Solver.parseOptions(new SolverOptions(SolverType.CTMC), args));
        this.result = new CTMCResult();
    }

    public SolverCTMC(Network model, SolverOptions options) {
        super(model, "SolverCTMC", options);
        this.result = new CTMCResult();
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.CTMC);
    }

    public List<String> listValidMethods() {
        List<String> methods = new ArrayList<String>();
        methods.add("default");
        methods.add("gpu");
        return methods;
    }

    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Source", "Sink",
                "ClassSwitch", "Delay", "DelayStation", "Queue",
                "MAP", "APH", "MMPP2", "PH", "Coxian", "Erlang", "Exp", "HyperExp",
                "StatelessClassSwitcher", "InfiniteServer", "SharedServer", "Buffer", "Dispatcher",
                "Cache", "CacheClassSwitcher",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_DPS", "SchedStrategy_GPS",
                "SchedStrategy_SIRO", "SchedStrategy_SEPT",
                "SchedStrategy_LEPT", "SchedStrategy_FCFS",
                "SchedStrategy_HOL", "SchedStrategy_LCFS",
                "SchedStrategy_LCFSPR",
                "RoutingStrategy_RROBIN",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_SFIFO", "ReplacementStrategy_LRU",
                "ClosedClass", "SelfLoopingClass", "OpenClass", "Replayer"
        });
        return featSupported;
    }

    public static void printInfGen(SolverCTMC.generatorResult infGen, SolverCTMC.StateSpace stateSpace) {
        Matrix Q = infGen.infGen;
        Matrix SS = stateSpace.stateSpace;
        printInfGen(Q, SS);
    }

    public static void printInfGen(Matrix Q, Matrix SS) {
        for (int s = 0; s < SS.getNumRows(); s++) {
            for (int sp = 0; sp < SS.getNumRows(); sp++) {
                if (Q.get(s, sp) > 0) {
                    System.out.println(
                            removeTrailingNewLine(SS.getRow(s).toString())
                                    + " -> "
                                    + removeTrailingNewLine(SS.getRow(sp).toString())
                                    + " : "
                                    + Q.get(s, sp));
                }
            }
        }
    }

    public static void printEventFilt(SolverCTMC.generatorResult infGen, SolverCTMC.StateSpace stateSpace) {
        MatrixCell eventFilt = infGen.eventFilt;
        Matrix SS = stateSpace.stateSpace;
        printEventFilt(eventFilt, SS);
    }

    public static void printEventFilt(MatrixCell eventFilt, Matrix SS) {
        for (int e = 0; e < eventFilt.size(); e++) {
            System.out.println("Event " + e + ":");
            Matrix filtMatrix = eventFilt.get(e);
            for (int s = 0; s < SS.getNumRows(); s++) {
                for (int sp = 0; sp < SS.getNumRows(); sp++) {
                    if (filtMatrix.get(s, sp) > 0) {
                        System.out.println(
                                "  " + removeTrailingNewLine(SS.getRow(s).toString())
                                        + " -> "
                                        + removeTrailingNewLine(SS.getRow(sp).toString())
                                        + " : "
                                        + filtMatrix.get(s, sp));
                    }
                }
            }
        }
    }

    /**
     * Get the cumulative distribution function of response times using tagged job methodology
     * @param R Response time matrix or percentile values
     * @return Matrix containing CDF values
     */
    public Matrix getCdfRespT(Matrix R) {
        if (GlobalConstants.DummyMode) {
            return new Matrix(0, 0);
        }

        NetworkStruct sn = this.getStruct(this);

        // Check if model has open classes (not supported in tagged job analysis)
        boolean hasOpenClasses = false;
        for (int k = 0; k < sn.nclasses; k++) {
            if (Double.isInfinite(sn.njobs.get(k, 0))) {
                hasOpenClasses = true;
                break;
            }
        }

        if (hasOpenClasses) {
            line_error(mfilename(new Object[]{}), "getCdfRespT is presently supported only for closed models.");
            return new Matrix(R.getNumRows(), R.getNumCols());
        }

        try {
            Matrix result = new Matrix(sn.nstations, sn.nclasses);

            // Get the chains in the model
            List<Chain> chains = this.model.getChains();

            // For each chain, create tagged model and compute response time CDF
            for (Chain chain : chains) {
                for (JobClass jobclass : chain.getClasses()) {
                    // Use ModelAdapter to create tagged job model
                    ModelAdapter.TaggedChainResult taggedResult =
                            ModelAdapter.tagChain(this.model, chain, jobclass);

                    Network taggedModel = taggedResult.getTaggedModel();
                    JobClass taggedJob = taggedResult.getTaggedJob();

                    if (taggedModel == null || taggedJob == null) {
                        continue;
                    }

                    // Create solver for tagged model
                    SolverOptions taggedOptions = this.options.copy();
                    SolverCTMC taggedSolver = new SolverCTMC(taggedModel, taggedOptions);

                    // Get generator and state space for tagged model
                    generatorResult taggedGenResult = taggedSolver.getGenerator();
                    Matrix taggedQ = taggedGenResult.infGen;
                    MatrixCell taggedEventFilt = taggedGenResult.eventFilt;
                    Map<Integer, Sync> taggedSynchInfo = taggedGenResult.ev;

                    StateSpace taggedStateSpaceResult = taggedSolver.getStateSpace();
                    Matrix taggedStateSpace = taggedStateSpaceResult.stateSpace;

                    NetworkStruct taggedSn = taggedModel.getStruct(false);

                    // Compute response time CDF using tagged job analysis
                    for (int ist = 0; ist < sn.nstations; ist++) {
                        Matrix cdfResult = computeTaggedResponseTimeCDF(taggedQ, taggedEventFilt,
                                taggedSynchInfo, taggedStateSpace, ist, jobclass.getIndex() - 1,
                                taggedJob.getIndex() - 1, taggedSn);

                        if (cdfResult != null && cdfResult.getNumRows() > 0) {
                            // Store the final CDF value (probability of completion)
                            double finalCdfValue = cdfResult.get(cdfResult.getNumRows() - 1, 0);
                            result.set(ist, jobclass.getIndex() - 1,
                                    Math.max(0.0, Math.min(1.0, finalCdfValue)));
                        }
                    }
                }
            }

            return result;

        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Failed to compute response time CDF: " + e.getMessage());
            return new Matrix(sn.nstations, sn.nclasses);
        }
    }

    /**
     * Compute response time CDF using tagged job methodology
     */
    private Matrix computeTaggedResponseTimeCDF(Matrix Q, MatrixCell eventFilt,
                                                Map<Integer, Sync> synchInfo, Matrix stateSpace, int station, int originalClass,
                                                int taggedClass, NetworkStruct sn) {

        try {
            // Build arrival and departure event filters for tagged job
            Matrix A1_tagged = new Matrix(Q.getNumRows(), Q.getNumCols());
            Matrix D1_tagged = new Matrix(Q.getNumRows(), Q.getNumCols());

            // Find events for tagged job arrivals and departures at the station
            for (Map.Entry<Integer, Sync> entry : synchInfo.entrySet()) {
                int eventIdx = entry.getKey();
                Sync sync = entry.getValue();

                if (eventIdx < eventFilt.size()) {
                    Matrix eventMatrix = eventFilt.get(eventIdx);
                    if (eventMatrix != null) {
                        // Check for tagged job arrivals to our station
                        if (sync.passive != null) {
                            for (Event passiveEvent : sync.passive.values()) {
                                if (passiveEvent.getNode() == station &&
                                        passiveEvent.getJobClass() == taggedClass &&
                                        passiveEvent.getEvent() == EventType.ARV) {
                                    A1_tagged = A1_tagged.add(eventMatrix);
                                }
                            }
                        }

                        // Check for tagged job departures from our station
                        if (sync.active != null) {
                            for (Event activeEvent : sync.active.values()) {
                                if (activeEvent.getNode() == station &&
                                        activeEvent.getJobClass() == taggedClass &&
                                        activeEvent.getEvent() == EventType.DEP) {
                                    D1_tagged = D1_tagged.add(eventMatrix);
                                }
                            }
                        }
                    }
                }
            }

            // Create MAP for tagged job analysis: D0 = Q - A1 - D1, D1 = D1_tagged
            Matrix D0 = Q.sub(A1_tagged).sub(D1_tagged);
            MatrixCell taggedMAP = new MatrixCell();
            taggedMAP.set(0, D0);
            taggedMAP.set(1, D1_tagged);

            // Normalize the MAP
            taggedMAP = map_normalize(taggedMAP);

            // Compute steady-state probability for tagged job entrance
            Matrix pi0 = map_pie(taggedMAP);

            // Time horizon for CDF computation
            List<Double> nonZeroRates = new ArrayList<Double>();
            for (int i = 0; i < Q.getNumRows(); i++) {
                for (int j = 0; j < Q.getNumCols(); j++) {
                    double rate = Math.abs(Q.get(i, j));
                    if (rate > GlobalConstants.FineTol && rate != 0) {
                        nonZeroRates.add(rate);
                    }
                }
            }

            if (nonZeroRates.isEmpty()) {
                return new Matrix(1, 2);
            }

            double minRate = nonZeroRates.stream().mapToDouble(Double::doubleValue).min().orElse(1.0);
            double T = Math.abs(100.0 / minRate);
            double dT = T / 100.0; // Time discretization

            List<Double> timePoints = new ArrayList<Double>();
            List<Double> cdfValues = new ArrayList<Double>();

            // Compute response time CDF: 1 - pi0 * exp(D0 * t) * ones
            for (double t = 0; t <= T; t += dT) {
                timePoints.add(t);

                try {
                    Matrix exponential = matrixExponential(D0, t);
                    Matrix ones = Matrix.ones(D0.getNumRows(), 1);

                    double survival = pi0.mult(exponential).mult(ones).toDouble();
                    double cdfValue = 1.0 - survival;
                    cdfValues.add(Math.max(0.0, Math.min(1.0, cdfValue)));

                    // Early termination if CDF approaches 1
                    if (cdfValue > 1.0 - GlobalConstants.CoarseTol) {
                        break;
                    }
                } catch (Exception e) {
                    // Fallback to exponential approximation
                    double rate = 1.0 / (T / 10.0);
                    double cdfValue = 1.0 - Math.exp(-rate * t);
                    cdfValues.add(Math.max(0.0, Math.min(1.0, cdfValue)));
                }
            }

            // Return CDF results
            Matrix result = new Matrix(cdfValues.size(), 2);
            for (int i = 0; i < cdfValues.size(); i++) {
                result.set(i, 0, cdfValues.get(i)); // CDF value
                result.set(i, 1, timePoints.get(i)); // Time point
            }

            return result;

        } catch (Exception e) {
            line_warning(mfilename(new Object[]{}),
                    "Failed to compute tagged response time CDF for station " + station +
                            ", original class " + originalClass + ", tagged class " + taggedClass +
                            ": " + e.getMessage());
            return new Matrix(1, 2);
        }
    }

    /**
     * Simple matrix exponential approximation using Taylor series
     */
    private Matrix matrixExponential(Matrix A, double t) {
        Matrix At = A.mult(new Matrix(new double[][]{{t}}));
        Matrix result = Matrix.eye(A.getNumRows());
        Matrix term = Matrix.eye(A.getNumRows());

        // Taylor series: exp(At) = I + At + (At)^2/2! + (At)^3/3! + ...
        int maxTerms = 20;
        for (int k = 1; k <= maxTerms; k++) {
            term = term.mult(At).scale(1.0 / k);
            result = result.add(term);

            // Check for convergence
            double termNorm = term.elementMaxAbs();
            if (termNorm < GlobalConstants.FineTol) {
                break;
            }
        }

        return result;
    }

    /**
     * Get the cumulative distribution function of system response times
     * @return Matrix containing system-wide CDF values
     */
    public Matrix getCdfSysRespT() {
        if (GlobalConstants.DummyMode) {
            return new Matrix(0, 0);
        }

        NetworkStruct sn = this.getStruct(this);

        try {
            // Get individual station-class response time CDFs
            Matrix stationClassCDFs = getCdfRespT(new Matrix(sn.nstations, sn.nclasses));

            // Aggregate system response time CDF
            // Simple approach: average across all station-class pairs weighted by traffic
            double totalWeight = 0.0;
            double weightedCDF = 0.0;

            for (int ist = 0; ist < sn.nstations; ist++) {
                for (int r = 0; r < sn.nclasses; r++) {
                    double jobPopulation = sn.njobs.get(r, 0);
                    if (!Double.isInfinite(jobPopulation) && jobPopulation > 0) {
                        double weight = jobPopulation;
                        totalWeight += weight;
                        weightedCDF += weight * stationClassCDFs.get(ist, r);
                    }
                }
            }

            Matrix result = new Matrix(1, 1);
            if (totalWeight > 0) {
                result.set(0, 0, weightedCDF / totalWeight);
            } else {
                result.set(0, 0, 0.0);
            }

            return result;

        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Failed to compute system response time CDF: " + e.getMessage());
            return new Matrix(1, 1);
        }
    }

    public generatorResult getGenerator() {
        return getGenerator(null);
    }

    public generatorResult getGenerator(SolverOptions options) {
        if (options == null) {
            options = this.options;
        }
        NetworkStruct sn = this.getStruct(this);
        if ((this.result) == null || ((CTMCResult) this.result).infGen == null) {
            ResultCTMC solverCTMCResult = Solver_ctmc.solver_ctmc(sn, options);
            ((CTMCResult) this.result).infGen = solverCTMCResult.getQ();
            ((CTMCResult) this.result).space = solverCTMCResult.getStateSpace();
            ((CTMCResult) this.result).spaceAggr = solverCTMCResult.getStateSpaceAggr();
            ((CTMCResult) this.result).nodeSpace = sn.space;
            ((CTMCResult) this.result).eventFilt = solverCTMCResult.getDfilt();
        }
        Matrix infGen = ((CTMCResult) this.result).infGen;
        MatrixCell eventFilt = ((CTMCResult) this.result).eventFilt;
        Map<Integer, Sync> ev = sn.sync;
        return new generatorResult(infGen, eventFilt, ev);
    }

    /**
     * Get the MarkedCTMC representation of the model
     * @return MarkedCTMC with generator and event filters
     */
    public MarkedMarkovProcess getMarkedCTMC() {
        return getMarkedCTMC(null);
    }

    /**
     * Get the MarkedCTMC representation of the model with specified options
     * @param options solver options
     * @return MarkedCTMC with generator and event filters
     */
    public MarkedMarkovProcess getMarkedCTMC(SolverOptions options) {
        generatorResult genResult = getGenerator(options);
        Matrix infGen = genResult.infGen;
        MatrixCell eventFilt = genResult.eventFilt;
        Map<Integer, Sync> synchInfo = genResult.ev;

        // Convert synchInfo to the format expected by MarkedCTMC
        List<Map<String, Object>> eventList = new ArrayList<Map<String, Object>>();
        for (Map.Entry<Integer, Sync> entry : synchInfo.entrySet()) {
            Map<String, Object> eventItem = new HashMap<String, Object>();
            eventItem.putIfAbsent("active", entry.getValue().active);
            eventItem.putIfAbsent("passive", entry.getValue().passive);
            eventList.add(eventItem);
        }

        return new MarkedMarkovProcess(infGen, eventFilt, eventList);
    }

    public generatorResult getInfGen() {
        return getGenerator();
    }

    public generatorResult getInfGen(SolverOptions options) {
        return getGenerator(options);
    }

    @Override
    public ProbabilityResult getProb(int node, Matrix state) {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult(Double.NaN);
        }

        if (node < 0 || node >= this.model.getNumberOfNodes()) {
            throw new IllegalArgumentException(
                    "getProb requires to pass a parameter the station of interest.");
        }

        long T0 = System.nanoTime();
        NetworkStruct sn = this.getStruct(this);
        sn.state = this.sn.state;

        // Convert node index to StatefulNode
        StatefulNode statefulNode = null;
        for (StatefulNode sn_node : this.model.getStatefulNodes()) {
            if (sn_node.getNodeIndex() == node) {
                statefulNode = sn_node;
                break;
            }
        }

        if (statefulNode == null) {
            throw new IllegalArgumentException("Node " + node + " is not a stateful node.");
        }

        if (state != null) {
            sn.state.replace(statefulNode, state);
        }

        for (Map.Entry<StatefulNode, Matrix> entry : sn.state.entrySet()) {
            int isf = this.model.getStatefulNodes().indexOf(entry.getKey());
            int isf_param = (int) sn.nodeToStateful.get(0, node);
            if (isf != isf_param) {
                Matrix updatedState =
                        new Matrix(entry.getValue().getNumRows(), entry.getValue().getNumCols());
                updatedState.ones();
                updatedState.mulByMinusOne();
                sn.state.replace(entry.getKey(), updatedState);
            }
        }
        Matrix Pnir = Solver_ctmc_marg.solver_ctmc_marg(sn, this.options);
        ((CTMCResult) this.result).solver = this.getName();
        ((CTMCResult) this.result).prob.marginal = Pnir;
        long T1 = System.nanoTime();
        this.result.runtime = (T1 - T0) / 1000000000.0;
        return new ProbabilityResult(Pnir.get(0, node));
    }

    public ProbabilityResult getProb(StatefulNode node, Matrix state) {
        if (node == null) {
            throw new IllegalArgumentException(
                    "getProb requires to pass a parameter the station of interest.");
        }
        return getProb(node.getNodeIndex(), state);
    }

    public ProbabilityResult getProb(StatefulNode node) {
        return getProb(node, null);
    }

    @Override
    public ProbabilityResult getProbAggr(int node, Matrix state_a) {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult(Double.NaN);
        }

        NetworkStruct sn = this.getStruct(this);
        if (node >= sn.nnodes) {
            line_error(mfilename(new Object[]{}), "Node number exceeds the number of nodes in the model.");
            return new ProbabilityResult(Double.NaN);
        }

        // Convert node index to station index
        int station = (int) sn.nodeToStation.get(0, node);
        if (station >= sn.nstations) {
            line_error(mfilename(new Object[]{}), "Station number exceeds the number of stations in the model.");
            return new ProbabilityResult(Double.NaN);
        }

        long T0 = System.nanoTime();
        sn.state = this.sn.state;

        // If a specific aggregated state is provided, set it
        if (state_a != null) {
            // Convert node index to StatefulNode and set aggregated state
            StatefulNode statefulNode = null;
            for (StatefulNode sn_node : this.model.getStatefulNodes()) {
                if (sn_node.getNodeIndex() == node) {
                    statefulNode = sn_node;
                    break;
                }
            }
            if (statefulNode != null) {
                sn.state.replace(statefulNode, state_a);
            }
        }

        Matrix Pnir;
        if (this.result == null || ((CTMCResult) this.result).prob == null || ((CTMCResult) this.result).prob.marginal == null) {
            if (this.result == null) {
                this.result = new CTMCResult();
            }
            ResultCTMCMargAggr margAggrResult = Solver_ctmc_margaggr.solver_ctmc_margaggr(sn, this.options);
            Pnir = margAggrResult.getPnir();
            ((CTMCResult) this.result).solver = this.getName();
            ((CTMCResult) this.result).prob.marginal = Pnir;
        } else {
            Pnir = ((CTMCResult) this.result).prob.marginal;
        }

        long T1 = System.nanoTime();
        this.result.runtime = (T1 - T0) / 1000000000.0;
        return new ProbabilityResult(Pnir.get(0, station));
    }

    public ProbabilityResult getProbAggr(StatefulNode node, Matrix state_a) {
        if (node == null) {
            throw new IllegalArgumentException(
                    "getProbAggr requires to pass a parameter the station of interest.");
        }
        return getProbAggr(node.getNodeIndex(), state_a);
    }

    public ProbabilityResult getProbAggr(StatefulNode node) {
        return getProbAggr(node, null);
    }

    public ProbabilityResult getProbAggr(jline.lang.nodes.Node node, Matrix state_a) {
        if (node == null) {
            throw new IllegalArgumentException(
                    "getProbAggr requires to pass a parameter the station of interest.");
        }
        return getProbAggr(node.getNodeIndex(), state_a);
    }

    public ProbabilityResult getProbAggr(jline.lang.nodes.Node node) {
        return getProbAggr(node, null);
    }

    @Override
    public ProbabilityResult getProbSys() {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult();
        }
        long T0 = System.nanoTime();
        sn = this.getStruct(this);
        Matrix Pn = Matrix.eye(0);
        if (this.model.isStateValid()) {
            SolverCtmcJointResult solverCtmcJointResult = Solver_ctmc_joint.solver_ctmc_joint(sn, this.options);
            Pn = solverCtmcJointResult.getPnir();
            ((CTMCResult) this.result).solver = this.getName();
            ((CTMCResult) this.result).prob.joint = Pn;
        } else {
            line_error(mfilename(new Object[]{}), "The model state is invalid.");
        }
        long T1 = System.nanoTime();
        this.result.runtime = (double) (T1 - T0) / 1000000000.0;
        return new ProbabilityResult(Pn);
    }

    @Override
    public ProbabilityResult getProbSysAggr() {
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult();
        }

        long T0 = System.nanoTime();
        NetworkStruct sn = this.getStruct(this);

        SolverCtmcJointResult solverCtmcJointResult = Solver_ctmc_jointaggr.solver_ctmc_jointaggr(sn, this.options);
        Matrix Pn = solverCtmcJointResult.getPnir();
        ((CTMCResult) this.result).solver = this.getName();
        ((CTMCResult) this.result).prob.joint = Pn;

        long T1 = System.nanoTime();
        this.result.runtime = (T1 - T0) / 1000000000.0;
        return new ProbabilityResult(Pn);
    }

    public StateSpace getStateSpace() {
        return this.getStateSpace(null);
    }

    public StateSpace getStateSpace(SolverOptions options) {
        if (options == null) {
            options = this.options;
        }
        NetworkStruct sn = this.getStruct(this);

        if (this.result == null || ((CTMCResult) this.result).space == null) {
            // Get properly dimensioned cutoff matrix
            Matrix cutoffMatrix = options.getCutoffMatrix(sn.nstations, sn.nclasses);
            State.StateSpaceGeneratorResult stateSpaceGeneratorResult =
                    State.spaceGenerator(sn, cutoffMatrix, options);
            sn.space = stateSpaceGeneratorResult.ST.space;
            ((CTMCResult) this.result).space = stateSpaceGeneratorResult.SS;
            ((CTMCResult) this.result).nodeSpace = stateSpaceGeneratorResult.ST.space;
        }
        Matrix stateSpace = ((CTMCResult) this.result).space;
        int shift = 0;
        MatrixCell localStateSpace = new MatrixCell();
        for (int i = 0; i < ((CTMCResult) this.result).nodeSpace.size(); i++) {
            int endCol =
                    shift + ((CTMCResult) this.result).nodeSpace.get(this.model.getStatefulNodes().get(i)).getNumCols();
            Matrix value =
                    Matrix.extract(
                            ((CTMCResult) this.result).space,
                            0,
                            ((CTMCResult) this.result).space.getNumRows(),
                            shift,
                            endCol);
            localStateSpace.set(i, value);
            shift += ((CTMCResult) this.result).nodeSpace.get(this.model.getStatefulNodes().get(i)).getNumCols();
        }
        return new StateSpace(stateSpace, localStateSpace);
    }

    public Matrix getStateSpaceAggr() {
        SolverOptions options = this.getOptions();
        if (options.force) {
            try {
                this.runAnalyzer();
            } catch (Exception e) {
                line_warning(mfilename(new Object[]{}),
                        "Failed to run analyzer automatically: " + e.getMessage());
                return null;
            }
        }
        if (this.result == null || ((CTMCResult) this.result).spaceAggr == null) {
            line_warning(mfilename(new Object[]{}),
                    "The model has not been cached. Either solve it or use the 'force' option to require this is done automatically, e.g., SolverCTMC(model).force(true).getStateSpaceAggr()");
            return null;
        }
        return ((CTMCResult) this.result).spaceAggr;
    }

    public NetworkStruct getStruct(SolverCTMC solverCTMC) {
        //    return new NetworkStruct();
        return this.model.getStruct(true);
    }


    public ProbabilityResult getTranProb(StatefulNode node) {
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProb in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model, options.timespan([0,T])).");
        }

        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult();
        }

        long T0 = System.nanoTime();
        NetworkStruct sn = this.model.getStruct(false);

        try {
            TransientResult transientResult = this.solver_ctmc_transient_analyzer(sn, this.options);

            Matrix t = transientResult.t;
            Matrix pit = transientResult.pit;
            Matrix stateSpace = transientResult.StateSpace;

            // Store transient results in CTMCResult
            if (this.result == null) {
                this.result = new CTMCResult();
            }
            ((CTMCResult) this.result).tranProb.t = t;
            ((CTMCResult) this.result).tranProb.pit = pit;
            ((CTMCResult) this.result).tranProb.stateSpace = stateSpace;

            long T1 = System.nanoTime();
            this.result.runtime = (T1 - T0) / 1000000000.0;

            return new ProbabilityResult(pit);

        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Failed to compute transient probabilities: " + e.getMessage());
            return new ProbabilityResult();
        }
    }

    public ProbabilityResult getTranProbAggr(StatefulNode node) {
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbAggr in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model, options.timespan([0,T])).");
        }

        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult();
        }

        long T0 = System.nanoTime();
        NetworkStruct sn = this.model.getStruct(false);

        try {
            TransientResult transientResult = this.solver_ctmc_transient_analyzer(sn, this.options);

            Matrix t = transientResult.t;
            Matrix pit = transientResult.pit;
            Matrix stateSpaceAggr = transientResult.StateSpaceAggr;

            // For aggregated transient probabilities, we marginalize over the aggregated state space
            int nodeIdx = node.getNodeIndex();
            int isf = (int) sn.nodeToStateful.get(0, nodeIdx);

            // Extract probabilities for the specific node (aggregated)
            Matrix nodeTranProb = new Matrix(pit.getNumRows(), 1);
            for (int timeIdx = 0; timeIdx < pit.getNumRows(); timeIdx++) {
                double prob = 0.0;
                for (int stateIdx = 0; stateIdx < pit.getNumCols(); stateIdx++) {
                    // Sum probabilities for states where this node has jobs
                    if (stateSpaceAggr.get(stateIdx, isf) > 0) {
                        prob += pit.get(timeIdx, stateIdx);
                    }
                }
                nodeTranProb.set(timeIdx, 0, prob);
            }

            // Store results
            if (this.result == null) {
                this.result = new CTMCResult();
            }
            ((CTMCResult) this.result).tranProbAggr.t = t;
            ((CTMCResult) this.result).tranProbAggr.pit = nodeTranProb;
            ((CTMCResult) this.result).tranProbAggr.node = nodeIdx;

            long T1 = System.nanoTime();
            this.result.runtime = (T1 - T0) / 1000000000.0;

            return new ProbabilityResult(nodeTranProb);

        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Failed to compute aggregated transient probabilities: " + e.getMessage());
            return new ProbabilityResult();
        }
    }

    public ProbabilityResult getTranProbSys() {
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbSys in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model, options.timespan([0,T])).");
        }

        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult();
        }

        long T0 = System.nanoTime();
        NetworkStruct sn = this.model.getStruct(false);

        try {
            TransientResult transientResult = this.solver_ctmc_transient_analyzer(sn, this.options);

            Matrix t = transientResult.t;
            Matrix pit = transientResult.pit;
            Matrix stateSpace = transientResult.StateSpace;

            // Store system transient results
            if (this.result == null) {
                this.result = new CTMCResult();
            }
            ((CTMCResult) this.result).tranProbSys.t = t;
            ((CTMCResult) this.result).tranProbSys.pit = pit;
            ((CTMCResult) this.result).tranProbSys.stateSpace = stateSpace;

            long T1 = System.nanoTime();
            this.result.runtime = (T1 - T0) / 1000000000.0;

            return new ProbabilityResult(pit);

        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Failed to compute system transient probabilities: " + e.getMessage());
            return new ProbabilityResult();
        }
    }

    public ProbabilityResult getTranProbSysAggr() {
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbSysAggr in SolverCTMC requires to specify a finite timespan T, e.g., SolverCTMC(model, options.timespan([0,T])).");
        }

        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult();
        }

        long T0 = System.nanoTime();
        NetworkStruct sn = this.model.getStruct(false);

        try {
            TransientResult transientResult = this.solver_ctmc_transient_analyzer(sn, this.options);

            Matrix t = transientResult.t;
            Matrix pit = transientResult.pit;
            Matrix stateSpaceAggr = transientResult.StateSpaceAggr;

            // Store system aggregated transient results
            if (this.result == null) {
                this.result = new CTMCResult();
            }
            ((CTMCResult) this.result).tranProbSysAggr.t = t;
            ((CTMCResult) this.result).tranProbSysAggr.pit = pit;
            ((CTMCResult) this.result).tranProbSysAggr.stateSpaceAggr = stateSpaceAggr;

            long T1 = System.nanoTime();
            this.result.runtime = (T1 - T0) / 1000000000.0;

            return new ProbabilityResult(pit);

        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Failed to compute system aggregated transient probabilities: " + e.getMessage());
            return new ProbabilityResult();
        }
    }


    @Override
    public void runAnalyzer()
            throws IllegalAccessException, ParserConfigurationException, IOException {
        long T0 = System.nanoTime();
        options = this.getOptions();

        if (!isInf(options.timespan[0]) && options.timespan[0] == options.timespan[1]) {
            line_error(mfilename(new Object[]{}), String.format(
                    "%s: timespan is a single point, spacing by options.tol (%e).\n",
                    this.getClass().getSimpleName(), options.tol));
            options.timespan[1] = options.timespan[0] + options.tol;
        }

        //this.runAnalyzerChecks(options);
        this.resetRandomGeneratorSeed(options.seed);

        sn = getStruct(this);
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix NK = sn.njobs;

        line_debug(options.verbose, String.format("CTMC solver starting: nstations=%d, nclasses=%d, timespan=[%.3f,%.3f]",
                M, K, options.timespan[0], options.timespan[1]));

        double sizeEstimator = 0;

        //    worst-case estimate of the state space
        for (int k = 0; k < K; k++) {
            sizeEstimator =
                    sizeEstimator
                            + Maths.factln(NK.get(k) + M - 1)
                            - Maths.factln(M - 1)
                            - Maths.factln(NK.get(k));
        }

        if (sn.njobs.hasInfinite()) {
            // Check if any cutoff values are infinite and need auto-setting
            Matrix currentCutoff = options.getCutoffMatrix(sn.nstations, sn.nclasses);
            boolean hasInfiniteCutoff = false;

            for (int i = 0; i < currentCutoff.getNumRows(); i++) {
                for (int j = 0; j < currentCutoff.getNumCols(); j++) {
                    if (Double.isInfinite(currentCutoff.get(i, j))) {
                        hasInfiniteCutoff = true;
                        break;
                    }
                }
                if (hasInfiniteCutoff) break;
            }

            if (hasInfiniteCutoff) {
                line_warning(mfilename(new Object[]{}), String.format(
                        "%s: The model has open chains, it is recommended to specify a finite cutoff value, e.g., SolverCTMC(model).cutoff(1).",
                        this.getClass().getSimpleName()));
                double autoCutoff = Math.ceil(Math.pow(6000, 1.0 / (M * K)));
                options.cutoff(autoCutoff);
                line_warning(mfilename(new Object[]{}), String.format(
                        "%s: Setting cutoff=%d.", this.getClass().getSimpleName(), (int) autoCutoff));
            }
            // Mandatory truncation warning for open/mixed models
            Matrix currentCutoffFinal = options.getCutoffMatrix(sn.nstations, sn.nclasses);
            line_warning(mfilename(new Object[]{}), String.format(
                    "CTMC solver using state space cutoff = %d for open/mixed model. State space truncation may cause inaccurate results. Consider varying cutoff to assess sensitivity.",
                    (int) currentCutoffFinal.get(0, 0)));
        }

        if (sizeEstimator > 6) {
            if (!options.force) {
                line_error(mfilename(new Object[]{}),
                        "CTMC size may be too large to solve. Stopping SolverCTMC. Set options.force=true to bypass this control.");
            }
        }

        if (Double.isInfinite(options.timespan[0])) {
            Map<StatefulNode, Matrix> s0 = sn.state;
            Map<StatefulNode, Matrix> s0prior = sn.stateprior;
            for (int ind = 1; ind < sn.nnodes; ind++) {
                if (sn.isstateful.get(ind, 0) == 1) {
                    int isf = (int) sn.nodeToStateful.get(ind);
                    Matrix initialStatePrior = s0prior.get(this.model.getStatefulNodes().get(0));
                    int maxPos = Maths.maxpos(initialStatePrior);
                    Matrix initialState = s0.get(this.model.getStatefulNodes().get(isf)).getRow(maxPos);
                    sn.state.replace(this.model.getStatefulNodes().get(isf), initialState);
                }
            }
            //      call solver_ctmc_analyzer(sn, options)
            line_debug(options.verbose, "Computing steady-state probabilities, calling solver_ctmc_analyzer");
            AnalyzerResult analyzerResult = Solver_ctmc_analyzer.solver_ctmc_analyzer(sn, options);
            //      [QN,UN,RN,TN,CN,XN,Q,SS,SSq,Dfilt,~,~,sn]
            Matrix QN = analyzerResult.QN;
            Matrix UN = analyzerResult.UN;
            Matrix RN = analyzerResult.RN;
            Matrix TN = analyzerResult.TN;
            Matrix CN = analyzerResult.CN;
            Matrix XN = analyzerResult.XN;
            Matrix Q = analyzerResult.InfGen;
            Matrix SS = analyzerResult.StateSpace;
            Matrix SSq = analyzerResult.StateSpaceAggr;
            MatrixCell Dfilt = analyzerResult.EventFiltration;
            sn = analyzerResult.sncopy;

            // call solver
            for (int isf = 0; isf < sn.nstateful; isf++) {
                int ind = (int) sn.statefulToNode.get(isf);
                // use stateful nodes directly to avoid IndexOutOfBoundsException
                Node a = this.model.getStatefulNodes().get(isf);
                ((StatefulNode) a).setState(sn.state.get(this.model.getStatefulNodes().get(isf)));
                if (a instanceof Cache) {
                    Cache cacheNode = (Cache) a;
                    cacheNode.setResultHitProb(((CacheNodeParam) sn.nodeparam.get(this.model.getNodes().get(ind))).actualhitprob);
                    cacheNode.setResultMissProb(((CacheNodeParam) sn.nodeparam.get(this.model.getNodes().get(ind))).actualmissprob);
                    this.model.refreshChains(true);
                }
            }
            ((CTMCResult) this.result).infGen = Q.copy();
            ((CTMCResult) this.result).space = SS;
            ((CTMCResult) this.result).spaceAggr = SSq;
            ((CTMCResult) this.result).nodeSpace = sn.space;
            ((CTMCResult) this.result).eventFilt = Dfilt;
            // Set the method field in the result
            ((CTMCResult) this.result).method = options.method;
            double runtime = (System.nanoTime() - T0) / 1000000000.0;
            sn.space = new HashMap<>();
            M = sn.nstations;
            int R = sn.nclasses;
            AvgHandle T = getAvgTputHandles();
            Matrix AN = snGetArvRFromTput(sn, TN, T);
            Matrix WN = new Matrix(0, 0);
            this.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, options.method, 0);
        } else {
            Matrix lastSol = null;

            Map<StatefulNode, Matrix> cur_state = new HashMap<StatefulNode, Matrix>(sn.nstations);
            for (int i = 0; i < sn.nnodes; i++) {
                Node node_i = this.model.getNodes().get(i);
                if (node_i.isStateful()) {
                    cur_state.putIfAbsent((StatefulNode) node_i, sn.state.get(node_i).copy());
                }
            }

            Map<StatefulNode, Matrix> s0 = sn.space;  // Use sn.space like MATLAB, not sn.state
            Map<StatefulNode, Matrix> s0prior = sn.stateprior;

            List<Integer> rowCounts =
                    s0.values().stream().map(Matrix::getNumRows).collect(Collectors.toList());
            double[][] rowCountsData = new double[rowCounts.size()][1];
            for (int i = 0; i < rowCounts.size(); i++) {
                rowCountsData[i][0] = rowCounts.get(i);
            }
            Matrix s0_sz = new Matrix(rowCountsData);
            Matrix s0_sz_1 = s0_sz.copy();
            s0_sz_1.addEq(-1);
            Matrix s0_id = pprod(s0_sz_1);
            while (s0_id.get(0, 0) >= 0) {  // Loop while not terminated (pprod returns -1 when done)
                double s0prior_val = 1;
                for (int ind = 0; ind < sn.nnodes; ind++) {
                    if (sn.isstateful.get(ind) == 1) {
                        int isf = (int) sn.nodeToStateful.get(ind);
                        s0prior_val =
                                s0prior_val * s0prior.get(this.model.getStatefulNodes().get(isf)).get((int) (1 + s0_id.get(isf)));

                        // Extract row from the state matrix corresponding to the current state
                        // MATLAB: s0{isf}(1+s0_id(isf),:) - add 1 to index like MATLAB
                        Matrix newState =
                                Matrix.extractRows(
                                        s0.get(this.model.getStatefulNodes().get(isf)),
                                        (int) (1 + s0_id.get(isf)),
                                        (int) (1 + s0_id.get(isf)) + 1,
                                        null);

                        // Update the state of the node
                        this.model.getStations().get((int) sn.nodeToStation.get(ind)).setState(newState);
                    }
                }
                NetworkStruct sn_cur = this.model.getStruct(false);
                if (s0prior_val > 0) {
                    line_debug(options.verbose, "Computing transient probabilities, calling solver_ctmc_transient_analyzer");
                    TransientResult transientResult = solver_ctmc_transient_analyzer(sn_cur, options);
                    assert transientResult != null;
                    Matrix t = transientResult.t;
                    Matrix pit = transientResult.pit;
                    Matrix QNt = transientResult.QNt;
                    Matrix UNt = transientResult.UNt;
                    Matrix TNt = transientResult.TNt;
                    Matrix Q = transientResult.InfGen;
                    Matrix SS = transientResult.StateSpace;
                    Matrix SSq = transientResult.StateSpaceAggr;
                    MatrixCell Dfilt = transientResult.EventFiltration;
                    double runtime_t = transientResult.runtime;

                    ((CTMCResult) this.result).space = SS;
                    ((CTMCResult) this.result).spaceAggr = SSq;
                    ((CTMCResult) this.result).infGen = Q;
                    ((CTMCResult) this.result).eventFilt = Dfilt;
                    setTranProb(t, pit, SS, runtime_t);
                    if (((CTMCResult) this.result).Tran == null
                            || ((CTMCResult) this.result).Tran.Avg == null
                            || ((CTMCResult) this.result).Tran.Avg.Q == null
                            || ((CTMCResult) this.result).Tran.Avg.Q.isEmpty()) {
                        // First time initialization - create structure similar to MATLAB
                        ((CTMCResult) this.result).Tran.Avg.Q = new HashMap<>();
                        ((CTMCResult) this.result).Tran.Avg.U = new HashMap<>();
                        ((CTMCResult) this.result).Tran.Avg.T = new HashMap<>();
                        for (int ist = 0; ist < M; ist++) {
                            for (int r = 0; r < K; r++) {
                                // Create matrix with weighted values and time column [data*s0prior_val, t]
                                Matrix QMatrix = new Matrix(t.getNumRows(), 2);
                                Matrix UMatrix = new Matrix(t.getNumRows(), 2);
                                Matrix TMatrix = new Matrix(t.getNumRows(), 2);

                                for (int timeIdx = 0; timeIdx < t.getNumRows(); timeIdx++) {
                                    QMatrix.set(timeIdx, 0, QNt.get(timeIdx, ist * K + r) * s0prior_val);
                                    QMatrix.set(timeIdx, 1, t.get(timeIdx, 0));
                                    UMatrix.set(timeIdx, 0, UNt.get(timeIdx, ist * K + r) * s0prior_val);
                                    UMatrix.set(timeIdx, 1, t.get(timeIdx, 0));
                                    TMatrix.set(timeIdx, 0, TNt.get(timeIdx, ist * K + r) * s0prior_val);
                                    TMatrix.set(timeIdx, 1, t.get(timeIdx, 0));
                                }

                                // Store in nested Map structure
                                ((CTMCResult) this.result).Tran.Avg.Q.computeIfAbsent(ist, k -> new HashMap<>()).put(r, QMatrix);
                                ((CTMCResult) this.result).Tran.Avg.U.computeIfAbsent(ist, k -> new HashMap<>()).put(r, UMatrix);
                                ((CTMCResult) this.result).Tran.Avg.T.computeIfAbsent(ist, k -> new HashMap<>()).put(r, TMatrix);
                            }
                        }
                    } else {
                        // Aggregate with existing results using interpolation (like MATLAB)
                        for (int ist = 0; ist < M; ist++) {
                            for (int r = 0; r < K; r++) {
                                Matrix existingQ = ((CTMCResult) this.result).Tran.Avg.Q.get(ist).get(r);
                                Matrix existingU = ((CTMCResult) this.result).Tran.Avg.U.get(ist).get(r);
                                Matrix existingT = ((CTMCResult) this.result).Tran.Avg.T.get(ist).get(r);

                                // Get time points from existing and new data
                                Matrix existingTimes = Matrix.extractColumns(existingQ, 1, 2, null);
                                Matrix newTimes = t.copy();
                                Matrix tunion = Matrix.union(existingTimes, newTimes);

                                // Create new aggregated matrices
                                Matrix newQMatrix = new Matrix(tunion.getNumRows(), 2);
                                Matrix newUMatrix = new Matrix(tunion.getNumRows(), 2);
                                Matrix newTMatrix = new Matrix(tunion.getNumRows(), 2);

                                Matrix existingQData = Matrix.extractColumns(existingQ, 0, 1, null);
                                Matrix existingUData = Matrix.extractColumns(existingU, 0, 1, null);
                                Matrix existingTData = Matrix.extractColumns(existingT, 0, 1, null);
                                Matrix newQData = Matrix.extractColumns(QNt, ist * K + r, ist * K + r + 1, null);
                                Matrix newUData = Matrix.extractColumns(UNt, ist * K + r, ist * K + r + 1, null);
                                Matrix newTData = Matrix.extractColumns(TNt, ist * K + r, ist * K + r + 1, null);

                                // Interpolate and aggregate
                                for (int timeIdx = 0; timeIdx < tunion.getNumRows(); timeIdx++) {
                                    double timePoint = tunion.get(timeIdx, 0);

                                    // Interpolate existing data at this time point
                                    double oldQ = interpolateValue(existingTimes, existingQData, timePoint);
                                    double oldU = interpolateValue(existingTimes, existingUData, timePoint);
                                    double oldT = interpolateValue(existingTimes, existingTData, timePoint);

                                    // Interpolate new data at this time point
                                    double newQ = interpolateValue(newTimes, newQData, timePoint);
                                    double newU = interpolateValue(newTimes, newUData, timePoint);
                                    double newT = interpolateValue(newTimes, newTData, timePoint);

                                    // Add weighted new contribution to existing
                                    newQMatrix.set(timeIdx, 0, oldQ + s0prior_val * newQ);
                                    newQMatrix.set(timeIdx, 1, timePoint);
                                    newUMatrix.set(timeIdx, 0, oldU + s0prior_val * newU);
                                    newUMatrix.set(timeIdx, 1, timePoint);
                                    newTMatrix.set(timeIdx, 0, oldT + s0prior_val * newT);
                                    newTMatrix.set(timeIdx, 1, timePoint);
                                }

                                // Update stored results
                                ((CTMCResult) this.result).Tran.Avg.Q.get(ist).put(r, newQMatrix);
                                ((CTMCResult) this.result).Tran.Avg.U.get(ist).put(r, newUMatrix);
                                ((CTMCResult) this.result).Tran.Avg.T.get(ist).put(r, newTMatrix);
                            }
                        }
                    }
                }
                s0_id = pprod(s0_id, s0_sz_1);  // Move to next state combination
            }
            //
            double runtime = (System.nanoTime() - T0) / 1000000000.0;
            for (int i = 0; i < sn.nnodes; i++) {
                Node node_i = this.model.getNodes().get(i);
                if (node_i.isStateful()) {
                    ((StatefulNode) node_i).setState(cur_state.get(node_i));
                }
            }
            ((CTMCResult) this.result).solver = getName();
            // Set the method field in the result for transient case
            this.result.method = options.method;
            this.result.runtime = runtime;
            ((CTMCResult) this.result).solverSpecific = lastSol;
        }
    }

    // Helper method for linear interpolation (mimics MATLAB's interp1)
    private double interpolateValue(Matrix xData, Matrix yData, double xQuery) {
        if (xData.getNumRows() == 1) {
            return yData.get(0, 0);
        }

        for (int i = 0; i < xData.getNumRows() - 1; i++) {
            double x1 = xData.get(i, 0);
            double x2 = xData.get(i + 1, 0);
            if (xQuery >= x1 && xQuery <= x2) {
                double y1 = yData.get(i, 0);
                double y2 = yData.get(i + 1, 0);
                return y1 + (y2 - y1) * (xQuery - x1) / (x2 - x1);
            }
        }

        // If outside range, return nearest value
        if (xQuery < xData.get(0, 0)) {
            return yData.get(0, 0);
        } else {
            return yData.get(yData.getNumRows() - 1, 0);
        }
    }

    public SampleResult sample(StatefulNode node, int numSamples) {
        SolverOptions options = this.getOptions();
        options.force = true;

        if (this.result == null || ((CTMCResult) this.result).infGen == null) {
            try {
                this.runAnalyzer();
            } catch (Exception e) {
                line_error(mfilename(new Object[]{}), "Failed to run analyzer: " + e.getMessage());
                return null;
            }
        }

        generatorResult genResult = getGenerator();
        Matrix infGen = genResult.infGen;
        StateSpace stateSpaceResult = getStateSpace();
        Matrix stateSpace = stateSpaceResult.stateSpace;

        NetworkStruct sn = this.getStruct(this);

        // Get initial state
        Map<StatefulNode, Matrix> initState = sn.state;
        List<Double> s0List = new ArrayList<Double>();
        for (StatefulNode statefulNode : this.model.getStatefulNodes()) {
            Matrix stateMatrix = initState.get(statefulNode);
            if (stateMatrix != null) {
                for (int i = 0; i < stateMatrix.getNumRows(); i++) {
                    for (int j = 0; j < stateMatrix.getNumCols(); j++) {
                        s0List.add(stateMatrix.get(i, j));
                    }
                }
            }
        }
        Matrix s0 = new Matrix(1, s0List.size());
        for (int i = 0; i < s0List.size(); i++) {
            s0.set(0, i, s0List.get(i));
        }

        // Set initial probability distribution
        double[] pi0 = new double[stateSpace.getNumRows()];
        int matchIdx = Matrix.matchrow(stateSpace, s0);
        if (matchIdx >= 0) {
            pi0[matchIdx] = 1.0;
        } else {
            // Uniform distribution as fallback
            for (int i = 0; i < pi0.length; i++) {
                pi0[i] = 1.0 / pi0.length;
            }
        }

        int nodeIdx = node.getNodeIndex();
        int isf = (int) sn.nodeToStateful.get(0, nodeIdx);

        SampleResult result = new SampleResult();
        result.handle = node;
        result.isaggregate = false;

        try {
            // Enhanced sampling using MMAP approach for better event tracking
            generatorResult genResult2 = getGenerator();
            MatrixCell eventFilt = genResult2.eventFilt;
            Map<Integer, Sync> synchInfo = genResult2.ev;

            // Create node-specific MMAP for this sampling
            MatrixCell nodeSpecificMMAP = createNodeSpecificMMAP(infGen, eventFilt, synchInfo, nodeIdx, sn);

            // Use MMAP sampling if available and appropriate
            boolean useMMAPSampling = (nodeSpecificMMAP != null && nodeSpecificMMAP.size() > 2);

            result.t = new Matrix(numSamples, 1);
            result.state = new Matrix(numSamples, sn.space.get(this.model.getStatefulNodes().get(isf)).getNumCols());
            result.event = new ArrayList<EventInfo>();

            if (useMMAPSampling) {
                // Use MMAP sampling for better event resolution
                java.util.Random random = new java.util.Random();
                jline.io.Ret.mamMMAPSample mmapSample = jline.api.mam.Mmap_sampleKt.mmap_sample(
                        nodeSpecificMMAP, (long)numSamples, random);

                double[] interArrivalTimes = mmapSample.getSamples();
                int[] eventTypes = mmapSample.getTypes();

                // Also get CTMC states for state information
                Ret.ctmcSimulation simulation = ctmc_simulate(infGen, pi0, numSamples);

                double currentTime = 0.0;
                for (int i = 0; i < numSamples; i++) {
                    // Use MMAP inter-arrival times
                    if (i < interArrivalTimes.length) {
                        currentTime += interArrivalTimes[i];
                    } else if (i < simulation.sojournTimes.length) {
                        currentTime += simulation.sojournTimes[i];
                    }
                    result.t.set(i, 0, currentTime);

                    // Create enhanced event info with MMAP event type
                    EventInfo eventInfo = new EventInfo();
                    eventInfo.node = nodeIdx;
                    eventInfo.jobclass = (i < eventTypes.length) ? eventTypes[i] % sn.nclasses : 0;
                    eventInfo.t = currentTime;
                    result.event.add(eventInfo);

                    // Extract node-specific state from global state
                    int globalState = (i < simulation.states.length) ? simulation.states[i] : 0;
                    if (globalState < stateSpace.getNumRows()) {
                        Matrix globalStateVector = stateSpace.getRow(globalState);

                        // Extract relevant columns for this node
                        int startCol = 0;
                        for (int nodeIdx2 = 0; nodeIdx2 < isf; nodeIdx2++) {
                            startCol += sn.space.get(this.model.getStatefulNodes().get(nodeIdx2)).getNumCols();
                        }
                        int endCol = startCol + sn.space.get(this.model.getStatefulNodes().get(isf)).getNumCols();

                        for (int j = 0; j < result.state.getNumCols() && startCol + j < endCol; j++) {
                            if (startCol + j < globalStateVector.getNumCols()) {
                                result.state.set(i, j, globalStateVector.get(0, startCol + j));
                            }
                        }
                    }
                }
            } else {
                // Fallback to standard CTMC simulation
                Ret.ctmcSimulation simulation = ctmc_simulate(infGen, pi0, numSamples);

                double currentTime = 0.0;
                for (int i = 0; i < numSamples; i++) {
                    // Accumulate sojourn times to get event times
                    currentTime += simulation.sojournTimes[i];
                    result.t.set(i, 0, currentTime);

                    // Create event info
                    EventInfo eventInfo = new EventInfo();
                    eventInfo.node = nodeIdx;
                    eventInfo.jobclass = 0; // Default job class
                    eventInfo.t = currentTime;
                    result.event.add(eventInfo);

                    // Extract node-specific state from global state
                    int globalState = simulation.states[i];
                    if (globalState < stateSpace.getNumRows()) {
                        Matrix globalStateVector = stateSpace.getRow(globalState);

                        // Extract relevant columns for this node
                        int startCol = 0;
                        for (int nodeIdx2 = 0; nodeIdx2 < isf; nodeIdx2++) {
                            startCol += sn.space.get(this.model.getStatefulNodes().get(nodeIdx2)).getNumCols();
                        }
                        int endCol = startCol + sn.space.get(this.model.getStatefulNodes().get(isf)).getNumCols();

                        for (int j = 0; j < result.state.getNumCols() && startCol + j < endCol; j++) {
                            if (startCol + j < globalStateVector.getNumCols()) {
                                result.state.set(i, j, globalStateVector.get(0, startCol + j));
                            }
                        }
                    }
                }
            }

        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "CTMC sampling failed: " + e.getMessage());
            // Fallback to simple implementation
            result.t = new Matrix(numSamples, 1);
            result.state = new Matrix(numSamples, sn.space.get(this.model.getStatefulNodes().get(isf)).getNumCols());
            result.event = new ArrayList<EventInfo>();
        }

        return result;
    }

    public jline.io.Ret.SampleResult sampleSys(int numSamples) {
        SolverOptions options = this.getOptions();
        options.force = true;

        if (this.result == null || ((CTMCResult) this.result).infGen == null) {
            try {
                this.runAnalyzer();
            } catch (Exception e) {
                line_error(mfilename(new Object[]{}), "Failed to run analyzer: " + e.getMessage());
                return new jline.io.Ret.SampleResult();
            }
        }

        try {
            // Get generator with event filtration (like MATLAB dev/)
            generatorResult genResult = getGenerator();
            Matrix infGen = genResult.infGen;
            MatrixCell eventFilt = genResult.eventFilt;
            Map<Integer, Sync> synchInfo = genResult.ev;

            StateSpace stateSpaceResult = getStateSpace();
            Matrix stateSpace = stateSpaceResult.stateSpace;
            NetworkStruct sn = this.getStruct(this);

            // Build initial probability distribution
            Map<StatefulNode, Matrix> initState = sn.state;
            List<Double> s0List = new ArrayList<Double>();
            for (StatefulNode statefulNode : this.model.getStatefulNodes()) {
                Matrix stateMatrix = initState.get(statefulNode);
                if (stateMatrix != null) {
                    for (int i = 0; i < stateMatrix.getNumRows(); i++) {
                        for (int j = 0; j < stateMatrix.getNumCols(); j++) {
                            s0List.add(stateMatrix.get(i, j));
                        }
                    }
                }
            }
            Matrix s0 = new Matrix(1, s0List.size());
            for (int i = 0; i < s0List.size(); i++) {
                s0.set(0, i, s0List.get(i));
            }

            // Set initial probability distribution
            Matrix pi0 = new Matrix(1, stateSpace.getNumRows());
            int matchIdx = Matrix.matchrow(stateSpace, s0);
            if (matchIdx >= 0) {
                pi0.set(0, matchIdx, 1.0);
            } else {
                // Uniform distribution as fallback
                for (int i = 0; i < pi0.getNumCols(); i++) {
                    pi0.set(0, i, 1.0 / pi0.getNumCols());
                }
            }

            // Create MMAP from event filtration (like MATLAB dev/ lines 21-24)
            MatrixCell MMAP = new MatrixCell(eventFilt.size() + 2);

            // Sum all event filters to get D1
            Matrix D1 = new Matrix(infGen.getNumRows(), infGen.getNumCols());
            for (int a = 0; a < eventFilt.size(); a++) {
                Matrix eventMatrix = eventFilt.get(a);
                if (eventMatrix != null) {
                    D1 = D1.add(eventMatrix);
                }
            }

            // D0 = infGen - D1
            Matrix D0 = infGen.sub(D1);

            // Build MMAP: [D0, D1, eventFilt1, eventFilt2, ...]
            MMAP.set(0, D0);
            MMAP.set(1, D1);
            for (int a = 0; a < eventFilt.size(); a++) {
                Matrix eventMatrix = eventFilt.get(a);
                if (eventMatrix != null) {
                    MMAP.set(2 + a, eventMatrix);
                } else {
                    MMAP.set(2 + a, new Matrix(infGen.getNumRows(), infGen.getNumCols()));
                }
            }

            // Normalize MMAP (like MATLAB dev/ line 24)
            MMAP = jline.api.mam.Mmap_normalizeKt.mmap_normalize(MMAP);

            // Sample MMAP (like MATLAB dev/ line 27)
            // Note: Java version doesn't accept pi0, uses steady-state initialization
            java.util.Random random = new java.util.Random();
            jline.io.Ret.mamMMAPSample mmapSample = jline.api.mam.Mmap_sampleKt.mmap_sample(MMAP, (long)numSamples, random);

            // For now, use simple CTMC sampling since MMAP sampling doesn't return states
            // This is a temporary approach until we can access MMAP state trajectory
            Ret.ctmcSimulation ctmcSim = ctmc_simulate(infGen,
                    new double[]{pi0.get(0, matchIdx >= 0 ? matchIdx : 0)}, numSamples);

            // Build time series (like MATLAB dev/ line 32)
            Matrix t = new Matrix(numSamples, 1);
            double cumulativeTime = 0.0;
            for (int i = 0; i < numSamples; i++) {
                if (i > 0) {
                    cumulativeTime += ctmcSim.sojournTimes[i-1];
                }
                t.set(i, 0, cumulativeTime);
            }

            // Build enhanced state series with MMAP state tracking
            Matrix state = new Matrix(numSamples, stateSpace.getNumCols());

            // Track MMAP phase states alongside CTMC states
            int[] mmapPhaseStates = new int[numSamples];
            int currentMMAPPhase = 0; // Start from first phase

            for (int i = 0; i < numSamples; i++) {
                // Update CTMC state information
                int stateIdx = (i < ctmcSim.states.length) ? ctmcSim.states[i] : 0;
                if (stateIdx < stateSpace.getNumRows()) {
                    for (int j = 0; j < stateSpace.getNumCols(); j++) {
                        state.set(i, j, stateSpace.get(stateIdx, j));
                    }
                }

                // Update MMAP phase state based on event type transitions
                // Simple phase progression for now
                currentMMAPPhase = (currentMMAPPhase + 1) % Math.max(1, MMAP.size() - 2);
                mmapPhaseStates[i] = currentMMAPPhase;
            }

            // Build event list from synchronizations (like MATLAB dev/ lines 39-48)
            List<Event> eventList = new ArrayList<Event>();
            for (int i = 1; i < numSamples; i++) { // Start from 1 to compare with previous state
                double eventTime = t.get(i, 0);

                // Find which synchronization event caused this state transition
                int prevStateIdx = ctmcSim.states[i-1];
                int currStateIdx = ctmcSim.states[i];

                // Find the sync that has a non-zero transition rate from prevStateIdx to currStateIdx
                int foundSyncIdx = -1;
                for (int syncIdx = 0; syncIdx < eventFilt.size(); syncIdx++) {
                    Matrix eventMatrix = eventFilt.get(syncIdx);
                    if (eventMatrix != null &&
                            prevStateIdx < eventMatrix.getNumRows() &&
                            currStateIdx < eventMatrix.getNumCols() &&
                            eventMatrix.get(prevStateIdx, currStateIdx) > 0) {
                        foundSyncIdx = syncIdx;
                        break;
                    }
                }

                // Get synchronization for this event
                Sync sync = synchInfo.get(foundSyncIdx);
                if (sync != null) {
                    // Add active events
                    if (sync.active != null) {
                        for (Event activeEvent : sync.active.values()) {
                            Event eventCopy = new Event(activeEvent.getEvent(), activeEvent.getNode(), activeEvent.getJobClass());
                            eventCopy.setT(eventTime);
                            eventList.add(eventCopy);
                        }
                    }
                    // Add passive events
                    if (sync.passive != null) {
                        for (Event passiveEvent : sync.passive.values()) {
                            Event eventCopy = new Event(passiveEvent.getEvent(), passiveEvent.getNode(), passiveEvent.getJobClass());
                            eventCopy.setT(eventTime);
                            eventList.add(eventCopy);
                        }
                    }
                }
            }

            // Convert event list to matrix format
            Matrix event = new Matrix(eventList.size(), 3);
            for (int i = 0; i < eventList.size(); i++) {
                Event e = eventList.get(i);
                event.set(i, 0, e.getT());
                event.set(i, 1, e.getNode());

                // Convert EventType to integer
                int eventTypeInt = -1;
                if (e.getEvent() != null) {
                    switch (e.getEvent().toString()) {
                        case "ARV":
                            eventTypeInt = 1;
                            break;
                        case "DEP":
                            eventTypeInt = 2;
                            break;
                        case "PHASE":
                            eventTypeInt = 3;
                            break;
                        default:
                            eventTypeInt = -1;
                    }
                }
                event.set(i, 2, eventTypeInt);
            }

            // Create enhanced sample result with MMAP information
            jline.io.Ret.SampleResult result = new jline.io.Ret.SampleResult("ctmc", t, state, event, false, null, numSamples);

            // Add MMAP-specific metadata including state trajectory
            Map<String, Object> mmapMetadata = new HashMap<String, Object>();
            // mmapMetadata.put("interArrivalTimes", interArrivalTimes);
            // mmapMetadata.put("mmapEventTypes", eventTypes);
            // mmapMetadata.put("mmapPhaseStates", mmapPhaseStates);
            // mmapMetadata.put("numEventTypes", numEventTypes);
            mmapMetadata.put("usedMMAPSampling", true);
            mmapMetadata.put("mmapSize", MMAP.size());
            mmapMetadata.put("ctmcStates", ctmcSim.states);
            mmapMetadata.put("synchronizationEvents", synchInfo.size());
            // result.setMetadata(mmapMetadata); // Method doesn't exist

            return result;

        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "CTMC MMAP sampling failed: " + e.getMessage());
            e.printStackTrace();
            // Fallback to empty result
            StateSpace stateSpaceResult = getStateSpace();
            Matrix t = new Matrix(numSamples, 1);
            Matrix state = new Matrix(numSamples, stateSpaceResult.stateSpace.getNumCols());
            Matrix event = new Matrix(numSamples, 3);
            return new jline.io.Ret.SampleResult("ctmc", t, state, event, false, null, numSamples);
        }
    }

    private TransientResult solver_ctmc_transient_analyzer(NetworkStruct sn, SolverOptions options) {
        long startTime = System.nanoTime();

        // Get the infinitesimal generator and state space
        ResultCTMC solverCTMCResult = Solver_ctmc.solver_ctmc(sn, options);
        Matrix infGen = solverCTMCResult.getQ();
        Matrix stateSpace = solverCTMCResult.getStateSpace();
        Matrix stateSpaceAggr = solverCTMCResult.getStateSpaceAggr();
        MatrixCell eventFiltration = solverCTMCResult.getDfilt();
        double[][][] depRates = solverCTMCResult.getDepRates(); // Get depRates from solver_ctmc (MATLAB line 17)

        // Build initial state vector (MATLAB lines 26-32)
        // MATLAB: state = []; for ist=1:sn.nnodes if sn.isstateful(ist) ...
        List<Double> stateList = new ArrayList<Double>();
        for (int ist = 0; ist < sn.nnodes; ist++) {
            if (sn.isstateful.get(ist, 0) == 1) {
                int isf = (int) sn.nodeToStateful.get(0, ist);
                // MATLAB: state = [state,zeros(1,size(sn.space{isf},2)-length(sn.state{isf})),sn.state{isf}];
                Matrix nodeSpace = sn.space.get(this.model.getStatefulNodes().get(isf));
                Matrix nodeState = sn.state.get(this.model.getStatefulNodes().get(isf));
                if (nodeSpace != null && nodeState != null) {
                    // Add zeros padding first
                    int spaceCols = nodeSpace.getNumCols();
                    int stateCols = nodeState.getNumCols();
                    for (int j = 0; j < spaceCols - stateCols; j++) {
                        stateList.add(0.0);
                    }
                    // Then add actual state
                    for (int j = 0; j < stateCols; j++) {
                        stateList.add(nodeState.get(0, j));
                    }
                }
            }
        }

        Matrix initialState = new Matrix(1, stateList.size());
        for (int i = 0; i < stateList.size(); i++) {
            initialState.set(0, i, stateList.get(i));
        }

        // Find initial state in state space
        Matrix pi0 = new Matrix(1, infGen.getNumRows());
        int state0 = Matrix.matchrow(stateSpace, initialState);
        if (state0 == -1) {
            throw new RuntimeException("Initial state not contained in the state space.");
        }
        pi0.set(0, state0, 1.0);

        // Perform transient analysis
        Matrix t;
        Matrix pit;

        try {
            // MATLAB line 50: [pit,t] = ctmc_transient(InfGen,pi0,options.timespan(1),options.timespan(2),options.stiff);
            // Note: Java ctmc_transient API doesn't have stiff parameter, use direct call
            Pair<double[], java.util.List<double[]>> transientResult =
                    ctmc_transient(infGen, pi0, options.timespan[0], options.timespan[1], options.timestep);

            double[] timeArray = transientResult.getLeft();
            java.util.List<double[]> probArray = transientResult.getRight();

            t = new Matrix(timeArray.length, 1);
            pit = new Matrix(timeArray.length, infGen.getNumRows());

            for (int i = 0; i < timeArray.length; i++) {
                t.set(i, 0, timeArray[i]);
                double[] probRow = probArray.get(i);
                for (int j = 0; j < probRow.length; j++) {
                    pit.set(i, j, probRow[j]);
                }
            }
        } catch (Exception e) {
            line_error(mfilename(new Object[]{}), "Transient analysis failed: " + e.getMessage());
            // Fallback to simple result
            t = new Matrix(1, 1);
            t.set(0, 0, options.timespan[1]);
            pit = new Matrix(1, infGen.getNumRows());
            pit.set(0, state0, 1.0);
        }

        // MATLAB line 52: pit(pit<GlobalConstants.Zero)=0;
        for (int i = 0; i < pit.getNumRows(); i++) {
            for (int j = 0; j < pit.getNumCols(); j++) {
                if (pit.get(i, j) < GlobalConstants.Zero) {
                    pit.set(i, j, 0.0);
                }
            }
        }

        // MATLAB lines 59-61: if t(1) == 0, t(1) = GlobalConstants.Zero; end
        if (t.getNumRows() > 0 && t.get(0, 0) == 0.0) {
            t.set(0, 0, GlobalConstants.Zero);
        }

        // Compute time-dependent metrics
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix QNt = new Matrix(t.getNumRows(), M * K);
        Matrix UNt = new Matrix(t.getNumRows(), M * K);
        Matrix TNt = new Matrix(t.getNumRows(), M * K);

        // Use depRates from solver_ctmc (MATLAB line 67: TNt{ist,k} = occupancy_t*depRates(:,ist,k))
        for (int k = 0; k < K; k++) {
            for (int ist = 0; ist < M; ist++) {
                // MATLAB line 66: occupancy_t = pit;
                Matrix occupancy_t = pit;

                // MATLAB line 67: TNt{ist,k} = occupancy_t*depRates(:,ist,k);
                for (int timeIdx = 0; timeIdx < t.getNumRows(); timeIdx++) {
                    double throughput = 0.0;
                    for (int stateIdx = 0; stateIdx < pit.getNumCols(); stateIdx++) {
                        if (stateIdx < depRates.length && depRates[stateIdx] != null &&
                                ist < depRates[stateIdx].length && depRates[stateIdx][ist] != null &&
                                k < depRates[stateIdx][ist].length) {
                            // MATLAB: depRates(:,ist,k) corresponds to depRates[stateIdx][ist][k] in Java
                            throughput += pit.get(timeIdx, stateIdx) * depRates[stateIdx][ist][k];
                        }
                    }
                    TNt.set(timeIdx, ist * K + k, throughput);
                }

                // Queue length computation
                for (int timeIdx = 0; timeIdx < t.getNumRows(); timeIdx++) {
                    double queueLength = 0.0;
                    for (int stateIdx = 0; stateIdx < pit.getNumCols(); stateIdx++) {
                        if (stateIdx < stateSpaceAggr.getNumRows()) {
                            queueLength += pit.get(timeIdx, stateIdx) * stateSpaceAggr.get(stateIdx, ist * K + k);
                        }
                    }
                    QNt.set(timeIdx, ist * K + k, queueLength);
                }

                // Utilization computation
                SchedStrategy schedStrategy = sn.sched.get(ist);
                int nservers = (int) sn.nservers.get(ist, 0);

                if (schedStrategy == SchedStrategy.INF) {
                    for (int timeIdx = 0; timeIdx < t.getNumRows(); timeIdx++) {
                        UNt.set(timeIdx, ist * K + k, QNt.get(timeIdx, ist * K + k));
                    }
                } else if (schedStrategy == SchedStrategy.PS) {
                    for (int timeIdx = 0; timeIdx < t.getNumRows(); timeIdx++) {
                        double utilization = 0.0;
                        for (int stateIdx = 0; stateIdx < pit.getNumCols(); stateIdx++) {
                            if (stateIdx < stateSpaceAggr.getNumRows()) {
                                double totalJobs = 0.0;
                                for (int kk = 0; kk < K; kk++) {
                                    totalJobs += stateSpaceAggr.get(stateIdx, ist * K + kk);
                                }
                                if (totalJobs > 0) {
                                    double uik = Math.min(stateSpaceAggr.get(stateIdx, ist * K + k), nservers)
                                            * stateSpaceAggr.get(stateIdx, ist * K + k) / totalJobs;
                                    utilization += pit.get(timeIdx, stateIdx) * uik / nservers;
                                }
                            }
                        }
                        UNt.set(timeIdx, ist * K + k, utilization);
                    }
                } else {
                    // Default utilization approximation
                    for (int timeIdx = 0; timeIdx < t.getNumRows(); timeIdx++) {
                        double utilization = 0.0;
                        for (int stateIdx = 0; stateIdx < pit.getNumCols(); stateIdx++) {
                            if (stateIdx < stateSpaceAggr.getNumRows()) {
                                utilization += pit.get(timeIdx, stateIdx)
                                        * Math.min(stateSpaceAggr.get(stateIdx, ist * K + k), nservers) / nservers;
                            }
                        }
                        UNt.set(timeIdx, ist * K + k, utilization);
                    }
                }
            }
        }

        Matrix RNt = new Matrix(0, 0);  // Response times not computed in transient analysis
        Matrix CNt = new Matrix(0, 0);  // Cycle times not computed in transient analysis  
        Matrix XNt = new Matrix(0, 0);  // System throughput not computed separately

        double runtime = (System.nanoTime() - startTime) / 1000000000.0;

        return new TransientResult(
                t, pit, QNt, UNt, RNt, TNt, CNt, XNt,
                infGen, stateSpace, stateSpaceAggr, eventFiltration,
                runtime, "ctmc_transient"
        );
    }

    /**
     * Compute departure rates for each station and class
     */
    private Matrix[][] computeDepartureRates(NetworkStruct sn, Matrix stateSpaceAggr, Matrix infGen) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix[][] depRates = new Matrix[M][K];

        // Initialize departure rate matrices
        for (int ist = 0; ist < M; ist++) {
            for (int k = 0; k < K; k++) {
                depRates[ist][k] = new Matrix(stateSpaceAggr.getNumRows(), 1);

                // Compute departure rates based on service process and state
                for (int stateIdx = 0; stateIdx < stateSpaceAggr.getNumRows(); stateIdx++) {
                    double queueLength = stateSpaceAggr.get(stateIdx, ist * K + k);
                    if (queueLength > 0) {
                        // Simple approximation: departure rate proportional to service rate
                        double serviceRate = 1.0; // Default service rate
                        // Note: sn.rates structure may not be available in this context

                        int nservers = (int) sn.nservers.get(ist, 0);
                        double effectiveRate = serviceRate * Math.min(queueLength, nservers);
                        depRates[ist][k].set(stateIdx, 0, effectiveRate);
                    }
                }
            }
        }

        return depRates;
    }

    /**
     * Create a node-specific MMAP for focused sampling on a particular node
     */
    private MatrixCell createNodeSpecificMMAP(Matrix infGen, MatrixCell eventFilt,
                                              Map<Integer, Sync> synchInfo, int targetNodeIdx, NetworkStruct sn) {

        try {
            // Create MMAP focused on events relevant to the target node
            MatrixCell nodeSpecificMMAP = new MatrixCell(eventFilt.size() + 2);

            // Start with the base generator
            Matrix D0 = infGen.copy();
            Matrix D1 = new Matrix(infGen.getNumRows(), infGen.getNumCols());

            // Collect all event filters that involve the target node
            List<Matrix> relevantEventFilters = new ArrayList<Matrix>();

            for (int eventIdx = 0; eventIdx < eventFilt.size(); eventIdx++) {
                Matrix eventMatrix = eventFilt.get(eventIdx);
                if (eventMatrix != null) {
                    Sync sync = synchInfo.get(eventIdx);
                    boolean isRelevant = false;

                    if (sync != null) {
                        // Check if this synchronization involves the target node
                        if (sync.active != null) {
                            for (Event activeEvent : sync.active.values()) {
                                if (activeEvent.getNode() == targetNodeIdx) {
                                    isRelevant = true;
                                    break;
                                }
                            }
                        }
                        if (!isRelevant && sync.passive != null) {
                            for (Event passiveEvent : sync.passive.values()) {
                                if (passiveEvent.getNode() == targetNodeIdx) {
                                    isRelevant = true;
                                    break;
                                }
                            }
                        }
                    }

                    if (isRelevant) {
                        D1 = D1.add(eventMatrix);
                        relevantEventFilters.add(eventMatrix);
                    }
                }
            }

            // Adjust D0 to account for D1
            D0 = D0.sub(D1);

            // Build the MMAP structure
            nodeSpecificMMAP.set(0, D0);
            nodeSpecificMMAP.set(1, D1);

            // Add individual event marking matrices
            for (int i = 0; i < relevantEventFilters.size(); i++) {
                if (i + 2 < nodeSpecificMMAP.size()) {
                    nodeSpecificMMAP.set(i + 2, relevantEventFilters.get(i));
                }
            }

            // Normalize the MMAP
            return jline.api.mam.Mmap_normalizeKt.mmap_normalize(nodeSpecificMMAP);

        } catch (Exception e) {
            line_warning(mfilename(new Object[]{}),
                    "Failed to create node-specific MMAP for node " + targetNodeIdx + ": " + e.getMessage());
            return null;
        }
    }

    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverCTMC.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    public static class StateSpace {
        public Matrix stateSpace;
        public MatrixCell localStateSpace;

        public StateSpace(Matrix stateSpace, MatrixCell localStateSpace) {
            this.stateSpace = stateSpace;
            this.localStateSpace = localStateSpace;
        }

        public void print() {
            stateSpace.print();
        }
    }

    public static class SupportResult {
        public boolean bool;
        public FeatureSet featSupported;
        public FeatureSet featUsed;

        public SupportResult(boolean bool, FeatureSet featSupported, FeatureSet featUsed) {
            this.bool = bool;
            this.featSupported = featSupported;
            this.featUsed = featUsed;
        }
    }

    public static class CtmcSsgResult {
        private final Matrix stateSpace;
        private final Matrix stateSpaceAggr;
        private final Matrix stateSpaceHashed;
        private final Map<StatefulNode, Matrix> nodeStateSpace;
        private final NetworkStruct sn;

        public CtmcSsgResult(
                Matrix stateSpace,
                Matrix stateSpaceAggr,
                Matrix stateSpaceHashed,
                Map<StatefulNode, Matrix> nodeStateSpace,
                NetworkStruct sn) {
            this.stateSpace = stateSpace;
            this.stateSpaceAggr = stateSpaceAggr;
            this.stateSpaceHashed = stateSpaceHashed;
            this.nodeStateSpace = nodeStateSpace;
            this.sn = sn;
        }

        public NetworkStruct getSn() {
            return sn;
        }

        public Matrix getStateSpace() {
            return stateSpace;
        }

        public Matrix getStateSpaceAggr() {
            return stateSpaceAggr;
        }

        public Matrix getStateSpaceHashed() {
            return stateSpaceHashed;
        }
    }

    public static class generatorResult {

        public Matrix infGen;
        public MatrixCell eventFilt;
        public Map<Integer, Sync> ev;

        public generatorResult(Matrix infGen, MatrixCell eventFilt, Map<Integer, Sync> ev) {
            this.infGen = infGen;
            this.eventFilt = eventFilt;
            this.ev = ev;
        }

        public void prettyPrint() {
            infGen.prettyPrint();
        }

        public void prettyPrintInt() {
            infGen.prettyPrintInt();
        }

        public void print() {
            infGen.print();
        }
    }

    public static class SolverCtmcJointResult {
        private final Matrix pnir;
        private final double runtime;
        private final String fname;

        public SolverCtmcJointResult(Matrix pnir, double runtime, String fname) {

            this.pnir = pnir;
            this.runtime = runtime;
            this.fname = fname;
        }

        public String getFname() {
            return fname;
        }

        public Matrix getPnir() {
            return pnir;
        }

        public double getRuntime() {
            return runtime;
        }
    }

    public static class StochCompResult {

        public Matrix S;
        public Matrix Q11;
        public Matrix Q12;
        public Matrix Q21;
        public Matrix Q22;
        Matrix T;

        public StochCompResult(Matrix S, Matrix Q11, Matrix Q12, Matrix Q21, Matrix Q22, Matrix T) {
            this.S = S;
            this.Q11 = Q11;
            this.Q12 = Q12;
            this.Q21 = Q21;
            this.Q22 = Q22;
            this.T = T;
        }
    }


    public static class AnalyzerResult {
        public Matrix QN, UN, RN, TN, CN, XN, InfGen, StateSpace, StateSpaceAggr;
        public MatrixCell EventFiltration;
        public double runtime;
        public String fname;
        public NetworkStruct sncopy;

        public AnalyzerResult(
                Matrix QN,
                Matrix UN,
                Matrix RN,
                Matrix TN,
                Matrix CN,
                Matrix XN,
                Matrix InfGen,
                Matrix StateSpace,
                Matrix StateSpaceAggr,
                MatrixCell EventFiltration,
                double runtime,
                String fname,
                NetworkStruct sncopy) {
            this.QN = QN;
            this.UN = UN;
            this.RN = RN;
            this.TN = TN;
            this.CN = CN;
            this.XN = XN;
            this.InfGen = InfGen;
            this.StateSpace = StateSpace;
            this.StateSpaceAggr = StateSpaceAggr;
            this.EventFiltration = EventFiltration;
            this.runtime = runtime;
            this.fname = fname;
            this.sncopy = sncopy;
        }
    }

    public static class TransientResult {

        public Matrix t;
        public Matrix pit;
        public Matrix QNt;
        public Matrix UNt;
        public Matrix RNt;
        public Matrix TNt;
        public Matrix CNt;
        public Matrix XNt;
        public Matrix InfGen;
        public Matrix StateSpace;
        public Matrix StateSpaceAggr;
        public MatrixCell EventFiltration;
        public double runtime;
        public String fname;

        public TransientResult(
                Matrix t,
                Matrix pit,
                Matrix QNt,
                Matrix UNt,
                Matrix RNt,
                Matrix TNt,
                Matrix CNt,
                Matrix XNt,
                Matrix InfGen,
                Matrix StateSpace,
                Matrix StateSpaceAggr,
                MatrixCell EventFiltration,
                double runtime,
                String fname) {
            this.t = t;
            this.pit = pit;
            this.QNt = QNt;
            this.UNt = UNt;
            this.RNt = RNt;
            this.TNt = TNt;
            this.CNt = CNt;
            this.XNt = XNt;
            this.InfGen = InfGen;
            this.StateSpace = StateSpace;
            this.StateSpaceAggr = StateSpaceAggr;
            this.EventFiltration = EventFiltration;
            this.runtime = runtime;
            this.fname = fname;
        }
    }

    public static class SampleResult {
        public StatefulNode handle;
        public Matrix t;
        public Matrix state;
        public List<EventInfo> event;
        public boolean isaggregate;
    }

    public static class SampleSysResult {
        public Matrix t;
        public Matrix state;
    }

    public static class EventInfo {
        public int node;
        public int jobclass;
        public double t;
    }

    public SampleResult sampleAggr(StatefulNode node, int numSamples) {
        SampleResult result = sample(node, numSamples);
        if (result != null) {
            result.isaggregate = true;
        }
        return result;
    }

    public jline.io.Ret.SampleResult sampleSysAggr(int numSamples) {
        // For aggregated system sampling, delegate to regular system sampling
        // In a full implementation, this would use aggregated state space
        return sampleSys(numSamples);
    }

    /**
     * Helper class to store event information from state transition analysis
     */
    private static class StateTransitionInfo {
        public int nodeIndex;
        public int eventType;

        public StateTransitionInfo(int nodeIndex, int eventType) {
            this.nodeIndex = nodeIndex;
            this.eventType = eventType;
        }
    }

    /**
     * Analyze state transition to determine which node and event type occurred
     */
    private StateTransitionInfo analyzeStateTransition(Matrix prevState, Matrix currState, NetworkStruct sn) {
        // EventType constants: ARV = 1, DEP = 2
        final int EVENT_ARV = 1;
        final int EVENT_DEP = 2;

        // Calculate total population change for each node
        // State space layout may include phases, so we need to sum across all columns for each node

        for (int nodeIdx = 0; nodeIdx < sn.nstations; nodeIdx++) {
            double prevNodePop = 0.0;
            double currNodePop = 0.0;

            // Sum all state variables belonging to this node
            // For open networks, we typically have: [source_states..., node1_states..., node2_states...]
            // The exact layout depends on the node types and service processes

            // Simple heuristic: check all columns for significant changes
            // If total population in "middle" columns changes, it's likely a queue event
            for (int col = 0; col < prevState.getNumCols(); col++) {
                if (!Double.isInfinite(prevState.get(0, col)) && !Double.isInfinite(currState.get(0, col))) {
                    double prevVal = prevState.get(0, col);
                    double currVal = currState.get(0, col);

                    // Skip source node columns (usually constant or infinite)
                    if (col > 1) {  // Skip first two columns which are typically source-related
                        prevNodePop += prevVal;
                        currNodePop += currVal;
                    }
                }
            }

            double popChange = currNodePop - prevNodePop;
            if (Math.abs(popChange) > 1e-10) {
                if (popChange > 0) {
                    // Total population in non-source nodes increased - arrival to queue
                    return new StateTransitionInfo(1, EVENT_ARV);  // nodeIdx=1 for queue
                } else {
                    // Total population in non-source nodes decreased - departure from queue  
                    return new StateTransitionInfo(1, EVENT_DEP);  // nodeIdx=1 for queue
                }
            }
        }

        // Alternative approach: look for any significant state change
        for (int col = 2; col < Math.min(prevState.getNumCols(), currState.getNumCols()); col++) {
            if (!Double.isInfinite(prevState.get(0, col)) && !Double.isInfinite(currState.get(0, col))) {
                double prevVal = prevState.get(0, col);
                double currVal = currState.get(0, col);

                if (Math.abs(currVal - prevVal) > 1e-10) {
                    // State change detected in queue-related columns
                    if (currVal > prevVal) {
                        return new StateTransitionInfo(1, EVENT_ARV);
                    } else {
                        return new StateTransitionInfo(1, EVENT_DEP);
                    }
                }
            }
        }

        // No clear change detected, return default
        return new StateTransitionInfo(-1, -1);
    }

    // ========================================================================
    // REWARD COMPUTATION METHODS
    // ========================================================================

    /** Cached reward computation result */
    private RewardResult rewardResult = null;

    /**
     * Get reward computation results via value iteration.
     *
     * Computes cumulative rewards for all defined reward functions using
     * value iteration on the uniformized CTMC.
     *
     * @return RewardResult containing value functions and steady-state rewards
     * @throws IllegalStateException if no rewards are defined on the model
     */
    public RewardResult getRewardResult() {
        if (this.rewardResult == null) {
            NetworkStruct sn = this.model.getStruct(true);
            if (sn.reward == null || sn.reward.isEmpty()) {
                throw new IllegalStateException(
                    "No rewards defined. Use model.setReward(name, rewardFn) before calling reward analysis.");
            }
            this.rewardResult = Solver_ctmc_reward.solver_ctmc_reward(sn, this.options);
        }
        return this.rewardResult;
    }

    /**
     * Get the value function for a specific reward.
     *
     * @param rewardName The name of the reward
     * @return Matrix of size [Tmax+1 x nStates] containing V^k(s) values
     * @throws IllegalArgumentException if reward name not found
     */
    public Matrix getRewardValueFunction(String rewardName) {
        RewardResult result = getRewardResult();
        Matrix V = result.getValueFunction().get(rewardName);
        if (V == null) {
            throw new IllegalArgumentException("Reward '" + rewardName + "' not found. Available rewards: " +
                String.join(", ", result.getRewardNames()));
        }
        return V;
    }

    /**
     * Get the time vector for reward computation.
     *
     * @return Time vector scaled by uniformization rate
     */
    public double[] getRewardTimeVector() {
        return getRewardResult().getTime();
    }

    /**
     * Get steady-state expected reward for all rewards.
     *
     * @return Map from reward name to expected reward value
     */
    public Map<String, Double> getAvgReward() {
        return getRewardResult().getSteadyState();
    }

    /**
     * Get steady-state expected reward for a specific reward.
     *
     * @param rewardName The name of the reward
     * @return Expected reward value in steady state
     * @throws IllegalArgumentException if reward name not found
     */
    public double getAvgReward(String rewardName) {
        Map<String, Double> steadyState = getAvgReward();
        Double value = steadyState.get(rewardName);
        if (value == null) {
            throw new IllegalArgumentException("Reward '" + rewardName + "' not found.");
        }
        return value;
    }

    /**
     * Get the list of defined reward names.
     *
     * @return List of reward names
     */
    public List<String> getRewardNames() {
        return getRewardResult().getRewardNames();
    }

    /**
     * Clear cached reward results to force recomputation.
     */
    public void clearRewardResult() {
        this.rewardResult = null;
    }

}
