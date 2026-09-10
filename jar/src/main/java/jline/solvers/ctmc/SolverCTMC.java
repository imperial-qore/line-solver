package jline.solvers.ctmc;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.processes.MarkedMarkovProcess;
import jline.lang.processes.MarkovChain;
import jline.lang.processes.MarkovProcess;
import jline.solvers.fluid.handlers.FluidRateMultiplier;
import jline.lang.state.State;
import jline.lang.state.ToMarginal;
import jline.solvers.AvgHandle;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.io.Ret.ProbabilityResult;
import jline.solvers.ctmc.analyzers.RewardResult;
import jline.solvers.ctmc.analyzers.Solver_ctmc_analyzer;
import jline.solvers.ctmc.analyzers.Solver_ctmc_cftp_analyzer;
import jline.solvers.ctmc.analyzers.Solver_ctmc_mdd_analyzer;
import jline.solvers.ctmc.analyzers.Solver_ctmc_qrf_analyzer;
import jline.solvers.ctmc.analyzers.Solver_ctmc_reward;
import jline.solvers.ctmc.handlers.Solver_ctmc;
import jline.solvers.ctmc.handlers.Solver_ctmc_joint;
import jline.solvers.ctmc.handlers.Solver_ctmc_marg;
import jline.solvers.ctmc.handlers.Solver_ctmc_margaggr;
import jline.solvers.ctmc.handlers.Solver_ctmc_jointaggr;
import jline.solvers.ctmc.ResultCTMCMargAggr;
import jline.api.sn.SnNonmarkovToPh;
import jline.api.sym.SageRestEngine;
import jline.api.sym.SymEngine;
import jline.api.sym.SymEngines;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static jline.api.mam.Map_mean.map_mean;
import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.io.InputOutput.*;
import static jline.util.PopulationLattice.pprod;
import static jline.util.Utils.isInf;
import static jline.util.matrix.Matrix.removeTrailingNewLine;
import static jline.api.mc.Ctmc_simulate.ctmc_simulate;
import static jline.api.mc.Ctmc_makeinfgen.ctmc_makeinfgen;
import jline.api.mc.Ctmc_fau;
import static jline.api.mc.Ctmc_transient.ctmc_transient;
import static jline.api.mc.Ctmc_solve.ctmc_solve;
import static jline.api.mc.Ctmc_solve_reducible.ctmc_solve_reducible;
import static jline.api.mc.Dtmc_solve.dtmc_solve;
import static jline.api.mc.Dtmc_solve_reducible.dtmc_solve_reducible;
import static jline.api.mc.Dtmc_simulate.dtmc_simulate;
import static jline.api.mam.Map_normalize.map_normalize;
import static jline.api.mam.Map_pie.map_pie;
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

    /**
     * User-supplied CTMC in chain mode, null when the solver analyzes a Network.
     * In chain mode on a DTMC this holds its P-I image, which carries the same
     * stationary vector.
     */
    private MarkovProcess chainProcess;

    /** User-supplied DTMC in chain mode, null for a CTMC or for a Network. */
    private MarkovChain chainMatrix;

    /**
     * Solves a user-supplied CTMC directly, bypassing state-space generation.
     *
     * @param chain the continuous-time Markov chain to solve
     * @param args  solver options in key-value form
     */
    public SolverCTMC(MarkovProcess chain, Object... args) {
        this(chain, Solver.parseOptions(new SolverOptions(SolverType.CTMC), args));
    }

    /**
     * Solves a user-supplied CTMC directly, bypassing state-space generation.
     *
     * @param chain   the continuous-time Markov chain to solve
     * @param options solver options
     */
    public SolverCTMC(MarkovProcess chain, SolverOptions options) {
        super(null, "SolverCTMC", options);
        this.chainProcess = chain;
        this.result = new CTMCResult();
    }

    /**
     * Solves a user-supplied DTMC directly, bypassing state-space generation.
     *
     * @param chain the discrete-time Markov chain to solve
     * @param args  solver options in key-value form
     */
    public SolverCTMC(MarkovChain chain, Object... args) {
        this(chain, Solver.parseOptions(new SolverOptions(SolverType.CTMC), args));
    }

    /**
     * Solves a user-supplied DTMC directly, bypassing state-space generation.
     *
     * @param chain   the discrete-time Markov chain to solve
     * @param options solver options
     */
    public SolverCTMC(MarkovChain chain, SolverOptions options) {
        super(null, "SolverCTMC", options);
        this.chainMatrix = chain;
        this.chainProcess = chain.toCTMC();
        this.result = new CTMCResult();
    }

    /**
     * @return true when the solver was built from a MarkovProcess or a MarkovChain
     */
    public boolean isChainSolver() {
        return this.chainProcess != null;
    }

    /**
     * @return true when the solver was built from a MarkovChain (DTMC)
     */
    public boolean isDiscreteChain() {
        return this.chainMatrix != null;
    }

    /**
     * @return the transition matrix of the user-supplied DTMC (chain mode only)
     */
    public Matrix getTransMat() {
        if (!isDiscreteChain()) {
            throw new RuntimeException("getTransMat requires a SolverCTMC built from a MarkovChain.");
        }
        return this.chainMatrix.getTransMat();
    }

    /**
     * Guard for the entry points that need stations and classes.
     *
     * @param caller name of the calling method, used in the error message
     */
    private void assertNotChainModel(String caller) {
        if (isChainSolver()) {
            throw new RuntimeException(caller + " requires a Network model. This solver was built from a "
                    + (isDiscreteChain() ? "MarkovChain" : "MarkovProcess")
                    + ", which has no stations or classes: use getProbSys, getGenerator, getStateSpace, "
                    + "getTranProbSys or sampleSys instead.");
        }
    }

    /** Steady-state analysis of the user-supplied chain. */
    private void chainRunAnalyzer() {
        long T0 = System.nanoTime();
        Matrix infGen = this.chainProcess.getGenerator();
        int n = infGen.getNumRows();
        Matrix pi;
        if (isDiscreteChain()) {
            Matrix P = this.chainMatrix.getTransMat();
            pi = dtmc_solve(P);
            if (!isChainDistribution(pi, n)) {
                pi = dtmc_solve_reducible(P).getLeft();
            }
        } else {
            pi = ctmc_solve(infGen, this.options);
            if (!isChainDistribution(pi, n)) {
                pi = ctmc_solve_reducible(infGen).getLeft();
            }
        }
        CTMCResult res = (CTMCResult) this.result;
        res.solver = this.getName();
        res.infGen = infGen;
        res.space = chainStateSpace();
        res.pi = pi;
        res.prob.joint = pi;
        res.runtime = (System.nanoTime() - T0) / 1000000000.0;
    }

    /**
     * Rejects a solution the primary solver could not produce on a reducible
     * chain, so that the reducible fallback is used instead.
     */
    private static boolean isChainDistribution(Matrix pi, int n) {
        if (pi == null || pi.length() != n) {
            return false;
        }
        double total = 0;
        for (int i = 0; i < n; i++) {
            double p = pi.get(i);
            if (!Double.isFinite(p) || p < -GlobalConstants.FineTol) {
                return false;
            }
            total += p;
        }
        return Math.abs(total - 1) <= Math.sqrt(GlobalConstants.FineTol);
    }

    /**
     * Initial distribution of a chain-mode transient analysis or sample path:
     * options.init_sol when it matches the chain size, uniform otherwise. Built
     * fresh rather than reshaped, since reshaping a sparse matrix clears it.
     *
     * @param n number of states of the chain
     * @return the initial distribution as a row vector
     */
    private Matrix chainInitDistribution(int n) {
        Matrix src = this.options.init_sol;
        boolean given = src != null && src.length() == n;
        Matrix pi0 = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            pi0.set(0, i, given ? src.get(i) : 1.0 / n);
        }
        return pi0;
    }

    /** State space of the user-supplied chain, or the state indices. */
    private Matrix chainStateSpace() {
        Matrix space = isDiscreteChain() ? this.chainMatrix.getStateSpace() : this.chainProcess.getStateSpace();
        if (space != null && space.getNumRows() > 0) {
            return space;
        }
        int n = this.chainProcess.getGenerator().getNumRows();
        space = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            space.set(i, 0, i + 1);
        }
        return space;
    }

    /** Runs the chain analyzer unless its result is already cached. */
    private void chainEnsureAnalyzed() {
        if (((CTMCResult) this.result).pi == null) {
            chainRunAnalyzer();
        }
    }

    /**
     * Stationary probability of a single state of the user-supplied chain,
     * identified by its row in the chain state space or, when the chain carries
     * none, by its 1-based state index.
     *
     * @param state state row or state index
     * @return the stationary probability of that state
     */
    public double getProb(Matrix state) {
        if (!isChainSolver()) {
            throw new RuntimeException("getProb(Matrix) requires a SolverCTMC built from a MarkovProcess or a MarkovChain.");
        }
        chainEnsureAnalyzed();
        CTMCResult res = (CTMCResult) this.result;
        Matrix userSpace = isDiscreteChain() ? this.chainMatrix.getStateSpace() : this.chainProcess.getStateSpace();
        int idx;
        if (userSpace == null || userSpace.getNumRows() == 0) {
            if (state.length() != 1) {
                throw new RuntimeException("The chain carries no state space, so getProb requires a state index in 1.."
                        + res.pi.length() + ".");
            }
            idx = (int) state.get(0) - 1;
        } else {
            idx = Matrix.matchrow(userSpace, state);
            if (idx < 0) {
                throw new RuntimeException("The requested state is not in the chain state space.");
            }
        }
        if (idx < 0 || idx >= res.pi.length()) {
            throw new RuntimeException("The requested state index is out of range.");
        }
        return res.pi.get(idx);
    }

    /**
     * Transient distribution of the user-supplied chain over options.timespan.
     * A DTMC advances one step per unit of time, so the returned times are the
     * integer steps within the timespan.
     *
     * @return pair of the time points and the distribution at each of them
     */
    private Pair<Matrix, Matrix> chainTranProbSys() {
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbSys requires a finite timespan T, e.g., SolverCTMC(chain).timespan(0,T).");
        }
        int n = this.chainProcess.getGenerator().getNumRows();
        Matrix pi0 = chainInitDistribution(n);
        double t0 = Double.isFinite(this.options.timespan[0]) ? this.options.timespan[0] : 0;
        double t1 = this.options.timespan[1];
        if (isDiscreteChain()) {
            Matrix P = this.chainMatrix.getTransMat();
            int k0 = (int) Math.ceil(t0);
            int k1 = (int) Math.floor(t1);
            if (k1 < k0) {
                throw new RuntimeException("The timespan [" + t0 + "," + t1 + "] contains no integer step of the DTMC.");
            }
            Matrix t = new Matrix(k1 - k0 + 1, 1);
            Matrix pit = new Matrix(k1 - k0 + 1, n);
            Matrix pik = new Matrix(pi0);
            for (int k = 0; k < k0; k++) {
                pik = pik.mult(P);
            }
            for (int k = 0; k <= k1 - k0; k++) {
                t.set(k, 0, k0 + k);
                for (int i = 0; i < n; i++) {
                    pit.set(k, i, pik.get(0, i));
                }
                pik = pik.mult(P);
            }
            return new Pair<Matrix, Matrix>(t, pit);
        }
        Pair<double[], List<double[]>> tran = ctmc_transient(this.chainProcess.getGenerator(), pi0, t0, t1);
        double[] tvec = tran.getLeft();
        List<double[]> pivec = tran.getRight();
        Matrix t = new Matrix(tvec.length, 1);
        Matrix pit = new Matrix(tvec.length, n);
        for (int k = 0; k < tvec.length; k++) {
            t.set(k, 0, tvec[k]);
            for (int i = 0; i < n; i++) {
                pit.set(k, i, pivec.get(k)[i]);
            }
        }
        return new Pair<Matrix, Matrix>(t, pit);
    }

    /**
     * Sample path of the user-supplied chain, started from options.init_sol when
     * given and from the uniform distribution otherwise. A DTMC advances one
     * unit of time per step.
     *
     * @param numEvents number of transitions to sample
     * @return the sampled trajectory
     */
    private jline.io.Ret.SampleResult chainSampleSys(int numEvents) {
        chainEnsureAnalyzed();
        this.resetRandomGeneratorSeed(this.options.seed);
        Matrix space = ((CTMCResult) this.result).space;
        int n = this.chainProcess.getGenerator().getNumRows();
        Matrix pi0 = chainInitDistribution(n);
        Matrix t = new Matrix(numEvents, 1);
        Matrix state = new Matrix(numEvents, space.getNumCols());
        int[] sts;
        if (isDiscreteChain()) {
            sts = dtmc_simulate(this.chainMatrix.getTransMat(), pi0, numEvents);
            for (int k = 0; k < numEvents; k++) {
                t.set(k, 0, k);
            }
        } else {
            double[] pi0arr = new double[n];
            for (int i = 0; i < n; i++) {
                pi0arr[i] = pi0.get(i);
            }
            Ret.ctmcSimulation sim = ctmc_simulate(this.chainProcess.getGenerator(), pi0arr, numEvents,
                    RandomManager.getThreadRandomAsRandom());
            sts = sim.states;
            double clock = 0;
            for (int k = 0; k < numEvents; k++) {
                t.set(k, 0, clock);
                clock += sim.sojournTimes[k];
            }
        }
        for (int k = 0; k < numEvents && k < sts.length; k++) {
            for (int j = 0; j < space.getNumCols(); j++) {
                state.set(k, j, space.get(sts[k], j));
            }
        }
        return new jline.io.Ret.SampleResult("ctmc", t, state, new Matrix(numEvents, 3), false, null, numEvents);
    }

    @Override
    public jline.solvers.SolverResult getAvg() {
        assertNotChainModel("getAvg");
        return super.getAvg();
    }

    @Override
    public jline.solvers.NetworkAvgTable getAvgTable() {
        assertNotChainModel("getAvgTable");
        return super.getAvgTable();
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.CTMC);
    }

    /**
     * True when the worst-case CTMC state space of the model fits the host
     * memory budget. Same estimator and gate the analyzer runs, exposed so a
     * caller (e.g. SolverAUTO) can rank CTMC out before paying for state-space
     * generation. Mirrors MATLAB SolverCTMC.isStateSpaceTractable.
     *
     * @param model   the network under analysis
     * @param options solver options carrying cutoff, force and safety fraction
     * @return the gate decision; ok is true when CTMC is a viable candidate
     */
    public static MemoryGuard.GateResult isStateSpaceTractable(Network model, SolverOptions options) {
        if (options == null) {
            options = defaultOptions();
        }
        double logNstates;
        try {
            NetworkStruct sn = SnNonmarkovToPh.snNonmarkovToPh(model.getStruct(), options, false);
            logNstates = MemoryGuard.stateSpaceLogSize(sn, options);
        } catch (Exception e) {
            // An estimator failure must not be read as a refusal: the analyzer
            // runs its own gate and reports the real error.
            return new MemoryGuard.GateResult(true, e.getMessage());
        }
        return MemoryGuard.gate(logNstates, options.force, false,
                MemoryGuard.DEFAULT_SAFETY_FRACTION);
    }

    /**
     * Valid solution methods.
     *
     * <p>'exact' is an explicit alias for the default state-space path: it pins
     * the intent at the call site so an example or test cannot be re-baselined
     * by a later change of what 'default' selects, and it MUST stay
     * behaviourally identical to 'default'.</p>
     *
     * <p>'mdd' holds the reachable set in a decision diagram and solves K
     * coupled level-CTMCs instead of the |S|-state generator; it is exact on
     * product-form models and approximate otherwise, and is restricted to closed
     * single-class networks (Solver_ctmc_mdd_analyzer).</p>
     */
    /**
     * The forwarding address for the QRF reduction bounds, which are SolverBA's.
     *
     * <p>{@code runAnalyzer} carried this text and still does, for the
     * {@code enableChecks = false} path that skips the gate entirely -- but it
     * sits DOWNSTREAM of {@link jline.solvers.NetworkSolver#checkDeclaredMethod},
     * which had already reported the flat "the 'qrf.bas' method is unsupported
     * by this solver" and sent the caller looking for SolverBA on their own.
     * That gate asks this method first now, and both sites read the text from
     * here so they cannot drift into two answers.
     *
     * <p>Asks nothing of the model, which is what lets the name gate call it.
     *
     * @param method the requested method name
     * @return the forwarding address, or "" when the name is not a QRF one
     */
    @Override
    protected String unsupportedMethodReason(String method) {
        if (method == null || !method.startsWith("qrf")) {
            return "";
        }
        return "QRF bound method '" + method
                + "' has moved out of SolverCTMC into the dedicated SolverBA solver. "
                + "Use SolverBA(model, \"" + method + "\") (or aliases \"qr\"/\"lr\") instead.";
    }

    public List<String> listValidMethods() {
        List<String> methods = new ArrayList<String>();
        methods.add("default");
        methods.add("exact");
        methods.add("gpu");
        methods.add("mdd");
        methods.add("cftp");
        methods.add("cftp.approx");
        // qrf.* (Quadratic Reduction Framework) bounds moved to SolverBA.
        return methods;
    }

    /**
     * Per-method feature deltas applied to the base CTMC envelope.
     *
     * <p>Four of the six methods share it; "cftp"/"cftp.approx" and "mdd"
     * narrow it, because neither builds the explicit generator that carries the
     * rest of the envelope. Mirrors MATLAB SolverCTMC.getMethodFeatureSet.</p>
     *
     * @param method the concrete method name
     * @return the envelope of that method
     */
    @Override
    public FeatureSet getMethodFeatureSet(String method) {
        FeatureSet featSupported = getFeatureSet();
        if ("cftp".equals(method) || "cftp.approx".equals(method)) {
            // PERFECT SAMPLING FROM A BALANCE FUNCTION, not from a generator:
            // the sampler encodes the closed single-class product form of
            // Gordon-Newell and nothing else, so every construct outside it has
            // to leave the envelope. The class count and the station count have
            // no registry name and are checked structurally in
            // supportsModelMethod, against the same predicate the analyzer uses.
            //
            // The one-phase-per-station rule is deliberately NOT spelled as a
            // list of distribution names. The sampler refuses phases(i,0) > 1,
            // and a name is not a phase count: a one-phase Coxian passes and a
            // HyperExp does not, while Det/Gamma/Pareto only acquire their
            // phases in SnNonmarkovToPh.
            featSupported.setFalse(new String[]{
                    "OpenClass",
                    // Queue, Delay and Router are the only node kinds the sampler walks
                    "Source", "Sink", "RandomSource", "JobSink",
                    "ClassSwitch", "StatelessClassSwitcher",
                    "Cache", "CacheClassSwitcher", "CacheRetrieval",
                    "ReplacementStrategy_RR", "ReplacementStrategy_FIFO",
                    "ReplacementStrategy_SFIFO", "ReplacementStrategy_LRU",
                    "ReplacementStrategy_HLRU", "ReplacementStrategy_CLIMB",
                    "ReplacementStrategy_QLRU",
                    "Fork", "Join", "Forker", "Joiner",
                    "Place", "Transition", "Linkage", "Enabling", "Inhibiting",
                    "Timing", "Firing", "Storage",
                    // disciplines outside INF/PS/FCFS/SIRO/LCFSPR have no product form
                    "SchedStrategy_DPS", "SchedStrategy_GPS",
                    "SchedStrategy_SEPT", "SchedStrategy_LEPT",
                    "SchedStrategy_HOL", "SchedStrategy_LCFS",
                    "SchedStrategy_LCFSPRPRIO", "SchedStrategy_FCFSPRPRIO",
                    "SchedStrategy_FCFSPR", "SchedStrategy_LCFSPI", "SchedStrategy_FCFSPI",
                    "SchedStrategy_LCFSPIPRIO", "SchedStrategy_FCFSPIPRIO",
                    "SchedStrategy_PSPRIO", "SchedStrategy_DPSPRIO", "SchedStrategy_GPSPRIO",
                    "SchedStrategy_LPS", "SchedStrategy_PAS", "SchedStrategy_OI",
                    "SchedStrategy_POLLING",
                    "Region",
                    "LoadDependence", "ClassDependence", "JointDependence", "GlobalDependence",
                    // a state-dependent decision is not Markovian routing
                    "RoutingStrategy_RROBIN", "RoutingStrategy_WRROBIN",
                    "RoutingStrategy_JSQ", "RoutingStrategy_SQ", "RoutingStrategy_SDR",
                    // the balance function has no orbit, no abandonment and no
                    // outage; each is a State construct the sampler never walks
                    "Balking", "Reneging", "Retrial", "Breakdown",
                    // the Gordon-Newell balance function has no buffer:
                    // Solver_ctmc_cftp_analyzer.supportsReason refuses a finite one
                    "FiniteCapacity"});
        } else if ("mdd".equals(method)) {
            // The decision diagram holds the MARKING of a closed network; an
            // open stream makes it unbounded, so there is no finite diagram to
            // hold. The single-class rule is structural (no registry name for a
            // class count) and lives in supportsModelMethod. A stochastic Petri
            // net keeps the Place/Transition names: Spn_mdd reads the marking.
            //
            // A FORK-JOIN MODEL IS NEITHER of the two shapes it serves. The tag
            // augmentation a fork needs adds one auxiliary class per branch, so
            // the struct that reaches the analyzer is never single-class however
            // the model was written, and the level decomposition has no meaning
            // for a firing that does not conserve the per-chain population.
            //
            // THE DESCRIPTOR READS RATES, SERVERS, PHASES, DISCIPLINES AND THE
            // ROUTING CHAIN, AND NOTHING ELSE. An orbit, an abandonment, an
            // outage or a buffer would be dropped on the floor, so each leaves
            // the envelope here: the level decomposition reads no sn.cap or
            // sn.classcap and no sn.retrial*/balking*/breakdown field.
            featSupported.setFalse(new String[]{
                    "OpenClass", "Source", "Sink", "RandomSource", "JobSink",
                    "Fork", "Join", "Forker", "Joiner", "JoinPartial",
                    "Balking", "Reneging", "Retrial", "Breakdown",
                    "FiniteCapacity"});
        }
        return featSupported;
    }

    /**
     * The per-method rules the feature registry has no name for, asked of the
     * SAME predicates the analyzers use so that the report and the run cannot
     * answer differently.
     *
     * <p>Three of them: the class count and the station count that "cftp" and
     * "mdd" need (a class count is not a model feature), and the state-space
     * size that the explicit-generator methods need. The last one is why
     * "default"/"exact"/"gpu" were offered on models whose chain does not fit
     * memory -- the analyzer priced the state space and refused, and nothing
     * above it had asked.</p>
     *
     * <p>THE TWO STRUCTURAL PREDICATES ARE ASKED BEFORE THE FEATURE GATE, which
     * is the reverse of the usual order and deliberate: each is the analyzer's
     * own assert, so it refuses a strict superset of what the per-method
     * feature deltas refuse, and its wording names the offending station or
     * class count instead of a feature. Asking the feature gate first would
     * replace "the cftp method supports closed models only" with "(feature:
     * OpenClass)" on the very run the caller is about to make.</p>
     *
     * @param method the concrete method name
     * @return empty string if supported, else the offending reason
     */
    @Override
    public String supportsModelMethod(String method) {
        boolean isCftp = "cftp".equals(method) || "cftp.approx".equals(method);
        boolean isMdd = "mdd".equals(method);
        if (this.model != null) {
            NetworkStruct snGate = this.model.getStruct(false);
            if (isCftp || isMdd) {
                String shapeReason = isCftp
                        ? Solver_ctmc_cftp_analyzer.supportsReason(snGate, this.options)
                        : Solver_ctmc_mdd_analyzer.supportsReason(snGate);
                if (!shapeReason.isEmpty()) {
                    return shapeReason;
                }
            }
            // The fork-join model class, which EVERY method has to clear: the
            // tag augmentation runs before the state space, the decision diagram
            // and the sampler alike, so a model fjValidate refuses is refused
            // whichever name was asked for.
            String fjReason = ModelAdapter.fjSupportsReason(snGate);
            if (!fjReason.isEmpty()) {
                return fjReason;
            }
        }
        String reason = super.supportsModelMethod(method);
        if (!reason.isEmpty() || this.model == null || isCftp || isMdd) {
            return reason;
        }
        // The explicit state space is what the remaining methods enumerate, and
        // MemoryGuard.gate refuses it above the host budget. Asking the same
        // estimator here costs a combinatorial formula, not a state space, so
        // the report stays cheap.
        MemoryGuard.GateResult gate = isStateSpaceTractable(this.model, this.options);
        if (!gate.ok) {
            return gate.message;
        }
        return "";
    }

    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Source", "Sink",
                "ClassSwitch", "Delay", "DelayStation", "Queue", "Router",
                "MAP", "APH", "MMPP2", "MMAP", "PH", "Coxian", "Erlang", "Exp", "HyperExp", "ME",
                "Det", "Gamma", "Weibull", "Lognormal", "Pareto", "Uniform",
                "StatelessClassSwitcher", "InfiniteServer", "SharedServer", "Buffer", "Dispatcher",
                // Finite capacity regions: CTMC represents them exactly (the
                // blocked-job overflow buffer is part of the chain).
                "Region",
                "Cache", "CacheClassSwitcher", "CacheRetrieval",
                "Server", "JobSink", "RandomSource", "ServiceTunnel",
                "SchedStrategy_INF", "SchedStrategy_PS",
                "SchedStrategy_DPS", "SchedStrategy_GPS",
                "SchedStrategy_SIRO", "SchedStrategy_SEPT",
                "SchedStrategy_LEPT", "SchedStrategy_FCFS",
                "SchedStrategy_HOL", "SchedStrategy_LCFS",
                "SchedStrategy_LCFSPR", "SchedStrategy_LCFSPRPRIO", "SchedStrategy_FCFSPRPRIO",
                // the rest of the preempt family: AfterEventStation carries an arm
                // for all eight, so declaring three gated five reachable
                // disciplines off (matches SolverCTMC.m and solver_ctmc.py)
                "SchedStrategy_FCFSPR", "SchedStrategy_LCFSPI", "SchedStrategy_FCFSPI",
                "SchedStrategy_LCFSPIPRIO", "SchedStrategy_FCFSPIPRIO",
                "SchedStrategy_PSPRIO", "SchedStrategy_DPSPRIO", "SchedStrategy_GPSPRIO",
                "SchedStrategy_LPS",
                "SchedStrategy_PAS", "SchedStrategy_OI", "SchedStrategy_POLLING",
                "RoutingStrategy_RROBIN",
                "RoutingStrategy_WRROBIN",
                "RoutingStrategy_JSQ",
                "RoutingStrategy_SQ",
                "RoutingStrategy_SDR",
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_SFIFO", "ReplacementStrategy_LRU",
                "ReplacementStrategy_HLRU", "ReplacementStrategy_CLIMB", "ReplacementStrategy_QLRU",
                "ClosedClass", "SelfLoopingClass", "OpenClass", "Replayer",
                "OpenSignal", "ClosedSignal",
                "SignalType_NEGATIVE", "SignalType_CATASTROPHE", "SignalType_REPLY",
                "SignalBatchRemoval", "SignalRemovalPolicy",
                "Place", "Transition", "Linkage", "Enabling", "Inhibiting", "Timing", "Firing", "Storage",
                "Fork", "Join", "Forker", "Joiner",
                "Balking", "Reneging", "Retrial", "Breakdown",
                "LoadDependence", "ClassDependence", "JointDependence", "GlobalDependence",
                // c-server stations and binding buffers are both State constructs
                // (State.fromMarginal / afterEventStation): served by the explicit
                // generator, withdrawn from cftp and mdd in getMethodFeatureSet
                "MultiServer", "FiniteCapacity"
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
            line_error(mfilename(new Object(){}), "getCdfRespT is presently supported only for closed models.");
            return new Matrix(R.getNumRows(), R.getNumCols());
        }

        try {
            Matrix result = new Matrix(sn.nstations, sn.nclasses);
            // The curves computed on the way, so getCdfRespT(AvgHandle) can
            // return the distribution rather than the base class's fabricated
            // exponential of the mean.
            lastCdfCurves = new ArrayList<List<Matrix>>();
            for (int ist = 0; ist < sn.nstations; ist++) {
                List<Matrix> row = new ArrayList<Matrix>();
                for (int r = 0; r < sn.nclasses; r++) {
                    row.add(null);
                }
                lastCdfCurves.add(row);
            }
            List<List<Matrix>> lastCurves = lastCdfCurves;

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
                            // The CURVE is what a distribution getter owes its
                            // caller; this summary matrix keeps only its last
                            // point and is kept for the Matrix-typed overload's
                            // own callers. getCdfRespT(AvgHandle) returns the
                            // curves.
                            lastCurves.get(ist).set(jobclass.getIndex() - 1, cdfResult);
                            double finalCdfValue = cdfResult.get(cdfResult.getNumRows() - 1, 0);
                            result.set(ist, jobclass.getIndex() - 1,
                                    Math.max(0.0, Math.min(1.0, finalCdfValue)));
                        }
                    }
                }
            }

            return result;

        } catch (Exception e) {
            line_error(mfilename(new Object(){}), "Failed to compute response time CDF: " + e.getMessage());
            return new Matrix(sn.nstations, sn.nclasses);
        }
    }

    /**
     * Compute response time CDF using tagged job methodology
     */
    /** True when the matrix carries at least one nonzero rate. */
    private static boolean anyRate(Matrix M) {
        for (int i = 0; i < M.getNumRows(); i++) {
            for (int j = 0; j < M.getNumCols(); j++) {
                if (M.get(i, j) != 0.0) return true;
            }
        }
        return false;
    }

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

            // TWO MAPs, NOT ONE. The reference builds the ARRIVAL map to find
            // the state a tagged job sees when it joins the station, and a
            // separate DEPARTURE map whose D0 governs how long it then waits:
            //
            //     A   = map_normalize({Q - A1, A1});  pie = map_pie(A)
            //     D   = map_normalize({Q - D1, D1});  D0  = D{1}
            //     F(t) = 1 - pie exp(D0 t) 1
            //
            // The old code formed a single MAP with D0 = Q - A1 - D1 and took
            // pie from THAT, so the initial vector was the state seen at a
            // DEPARTURE and the sub-generator had the arrival transitions
            // removed as well. Both are wrong and neither is a tolerance.
            if (!anyRate(A1_tagged) || !anyRate(D1_tagged)) {
                // No tagged arrival or no tagged departure at this station: the
                // job never passes through it, so there is no law to report.
                return null;
            }

            MatrixCell arrivalMAP = new MatrixCell();
            arrivalMAP.set(0, Q.sub(A1_tagged));
            arrivalMAP.set(1, A1_tagged);
            arrivalMAP = map_normalize(arrivalMAP);
            Matrix pi0 = map_pie(arrivalMAP);

            MatrixCell departureMAP = new MatrixCell();
            departureMAP.set(0, Q.sub(D1_tagged));
            departureMAP.set(1, D1_tagged);
            departureMAP = map_normalize(departureMAP);
            Matrix D0 = departureMAP.get(0);

            // Time horizon: 100 events at the slowest rate in the chain, as the
            // reference chooses it.
            List<Double> nonZeroRates = new ArrayList<Double>();
            for (int i = 0; i < Q.getNumRows(); i++) {
                for (int j = 0; j < Q.getNumCols(); j++) {
                    double rate = Math.abs(Q.get(i, j));
                    if (rate > GlobalConstants.FineTol) {
                        nonZeroRates.add(rate);
                    }
                }
            }
            if (nonZeroRates.isEmpty()) {
                return null;
            }
            double minRate = Double.MAX_VALUE;
            for (double r : nonZeroRates) {
                minRate = Math.min(minRate, r);
            }
            final int intervals = 100000;
            double T = Math.abs(100.0 / minRate);
            double dT = T / intervals;

            // ONE matrix exponential, then propagate. The reference recomputes
            // expm(D0*t) at each of the 100001 grid points, which is the same
            // answer at a cost linear in the grid rather than constant.
            Matrix E = matrixExponential(D0, dT);
            int n = D0.getNumRows();
            double[] v = new double[n];
            for (int j = 0; j < n; j++) {
                v[j] = pi0.get(0, j);
            }

            List<Double> timePoints = new ArrayList<Double>();
            List<Double> cdfValues = new ArrayList<Double>();
            for (int k = 0; k <= intervals; k++) {
                if (k > 0) {
                    double[] w = new double[n];
                    for (int j = 0; j < n; j++) {
                        double acc = 0.0;
                        for (int i = 0; i < n; i++) {
                            acc += v[i] * E.get(i, j);
                        }
                        w[j] = acc;
                    }
                    v = w;
                }
                double survival = 0.0;
                for (int j = 0; j < n; j++) {
                    survival += v[j];
                }
                double cdfValue = 1.0 - survival;
                timePoints.add(k * dT);
                cdfValues.add(Math.max(0.0, Math.min(1.0, cdfValue)));
                if (cdfValue > 1.0 - GlobalConstants.CoarseTol) {
                    break;
                }
            }

            Matrix result = new Matrix(cdfValues.size(), 2);
            for (int i = 0; i < cdfValues.size(); i++) {
                result.set(i, 0, cdfValues.get(i)); // CDF value
                result.set(i, 1, timePoints.get(i)); // Time point
            }

            return result;

        } catch (Exception e) {
            line_warning(mfilename(new Object(){}),
                    "Failed to compute tagged response time CDF for station " + station +
                            ", original class " + originalClass + ", tagged class " + taggedClass +
                            ": " + e.getMessage());
            return new Matrix(1, 2);
        }
    }

    /**
     * exp(A t), by the scaling-and-squaring Pade routine on {@link Matrix}.
     *
     * <p>THIS USED TO BE A 20-TERM TAYLOR SERIES that never ran: it formed
     * {@code A.mult(new Matrix({{t}}))} to scale by t, which multiplies an
     * n-by-n matrix by a 1-by-1 one and throws MatrixDimensionException on
     * every call. The response-time CDF therefore always landed in its own
     * catch. Even with the scaling fixed, a truncated Taylor series is the
     * wrong tool for a sub-generator, whose norm grows with the horizon; the
     * Pade routine already in Matrix is what MATLAB's expm and the C++ port
     * both use.
     */
    private Matrix matrixExponential(Matrix A, double t) {
        return A.scale(t).expm();
    }

    /**
     * The per-(station, class) response-time CURVES computed by the most recent
     * {@link #getCdfRespT(Matrix)} call, each a (T x 2) matrix of [F(t) t].
     */
    private List<List<Matrix>> lastCdfCurves = null;

    /**
     * Response-time distribution, as the distribution and not as a summary.
     *
     * <p>THIS OVERRIDE IS THE POINT. {@code NetworkSolver.getCdfRespT(AvgHandle)}
     * fabricates an exponential law with the right mean, and the only other
     * method here takes a {@code Matrix}, so it OVERLOADS rather than overrides:
     * every caller holding a {@code NetworkSolver} silently received the
     * fabricated law even though the tagged-chain computation was available. The
     * curves are the ones MATLAB and C++ return, in their column order
     * [F(t) t].
     */
    @Override
    public Ret.DistributionResult getCdfRespT(AvgHandle R) {
        NetworkStruct sn = this.getStruct(this);
        getCdfRespT(new Matrix(sn.nstations, sn.nclasses));
        Ret.DistributionResult out =
                new Ret.DistributionResult(sn.nstations, sn.nclasses, "response_time");
        if (lastCdfCurves != null) {
            out.cdfData = lastCdfCurves;
        }
        return out;
    }

    /**
     * The SYSTEM response-time distribution: one law per CHAIN, not a number.
     *
     * <p>Port of matlab/src/solvers/CTMC/@SolverCTMC/getCdfSysRespT.m. Each
     * entry is a (T x 2) matrix of [F(t) t] for one chain, in chain order, as
     * MATLAB's {@code RD{1,c}} cell and the C++ {@code sysrespt} block are.
     *
     * <p>THE QUANTITY IS THE CYCLE TIME. The split is the tagged job's ARRIVAL
     * AT ITS OWN REFERENCE STATION, so a passage runs from one such arrival to
     * the next: the job's whole trip round the network, not its stay at one
     * station. That is why a single MAP suffices here where
     * {@link #getCdfRespT(Matrix)} needs two -- the arrival that starts the
     * passage and the one that ends it are the same event, so
     * {@code map_pie} and the sub-generator come from the same
     * {@code map_normalize({Q - D1, D1})}.
     *
     * <p>THIS REPLACED A WEIGHTED SCALAR AVERAGE. The old body called
     * {@code getCdfRespT} for its per-(station, class) summary values and
     * returned a single population-weighted mean of them in a 1x1 matrix. That
     * is not a distribution, and an average of per-station response-time values
     * is not the system response time even as a mean: the cycle time is their
     * SUM along the job's route, not their average.
     *
     * <p>Two constants differ from the per-station getter on purpose, matching
     * the reference: the grid is 10000 intervals rather than 100000, and the
     * truncation is at {@code 1 - FineTol} rather than {@code CoarseTol},
     * because a cycle time is longer and its tail matters more.
     */
    public List<Matrix> getCdfSysRespT() {
        NetworkStruct sn = this.getStruct(this);
        List<Matrix> RD = new ArrayList<Matrix>();
        if (GlobalConstants.DummyMode) {
            return RD;
        }

        for (int k = 0; k < sn.nclasses; k++) {
            if (Double.isInfinite(sn.njobs.get(k, 0))) {
                line_error(mfilename(new Object(){}),
                        "getCdfSysRespT is presently supported only for closed models.");
                return RD;
            }
        }

        int nchains = sn.chains.getNumElements();
        for (int c = 0; c < nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            int tagged_src = -1;
            for (int idx = 0; idx < inchain.length(); idx++) {
                int r = (int) inchain.get(idx);
                if (sn.njobs.get(r, 0) > 0) {
                    tagged_src = r;
                    break;
                }
            }
            if (tagged_src < 0) {
                RD.add(null);
                continue;
            }

            JobClass jobclass = this.model.getClassByIndex(tagged_src);
            Chain chain = this.model.getClassChain(jobclass);
            ModelAdapter.TaggedChainResult tr = ModelAdapter.tagChain(this.model, chain, jobclass);
            Network taggedModel = tr.getTaggedModel();
            JobClass taggedJob = tr.getTaggedJob();
            if (taggedModel == null || taggedJob == null) {
                RD.add(null);
                continue;
            }

            SolverOptions taggedOptions = this.options.copy();
            SolverCTMC taggedSolver = new SolverCTMC(taggedModel, taggedOptions);
            generatorResult g = taggedSolver.getGenerator();
            Matrix Q = g.infGen;
            MatrixCell filt = g.eventFilt;
            Map<Integer, Sync> ev = g.ev;
            NetworkStruct tsn = taggedModel.getStruct(false);

            int taggedClass = taggedJob.getIndex() - 1;
            // sn.refstat is a STATION index; the events carry NODE indices, so
            // the two must be mapped rather than compared directly. The
            // reference compares them as they stand, which coincides only when
            // every node is a station.
            int refStation = (int) tsn.refstat.get(taggedClass, 0);
            int refNode = (refStation >= 0 && refStation < tsn.stationToNode.length())
                    ? (int) tsn.stationToNode.get(refStation) : refStation;

            Matrix D1 = new Matrix(Q.getNumRows(), Q.getNumCols());
            for (Map.Entry<Integer, Sync> entry : ev.entrySet()) {
                int idx = entry.getKey();
                if (idx >= filt.size() || filt.get(idx) == null) {
                    continue;
                }
                Sync sync = entry.getValue();
                if (sync.passive == null) {
                    continue;
                }
                for (Event pas : sync.passive.values()) {
                    if (pas.getEvent() == EventType.ARV
                            && pas.getJobClass() == taggedClass
                            && pas.getNode() == refNode) {
                        D1 = D1.add(filt.get(idx));
                    }
                }
            }
            if (!anyRate(D1)) {
                RD.add(null);
                continue;
            }

            MatrixCell D = new MatrixCell();
            D.set(0, Q.sub(D1));
            D.set(1, D1);
            D = map_normalize(D);
            Matrix pie = map_pie(D);
            Matrix D0 = D.get(0);

            double minRate = Double.MAX_VALUE;
            for (int i = 0; i < Q.getNumRows(); i++) {
                for (int j = 0; j < Q.getNumCols(); j++) {
                    double rate = Math.abs(Q.get(i, j));
                    if (rate > GlobalConstants.FineTol) {
                        minRate = Math.min(minRate, rate);
                    }
                }
            }
            if (minRate == Double.MAX_VALUE) {
                RD.add(null);
                continue;
            }
            final int intervals = 10000;
            double T = Math.abs(100.0 / minRate);
            double dT = T / intervals;

            Matrix E = matrixExponential(D0, dT);
            int n = D0.getNumRows();
            double[] v = new double[n];
            for (int j = 0; j < n; j++) {
                v[j] = pie.get(0, j);
            }
            List<Double> tvals = new ArrayList<Double>();
            List<Double> Fvals = new ArrayList<Double>();
            for (int step = 0; step <= intervals; step++) {
                if (step > 0) {
                    double[] w = new double[n];
                    for (int j = 0; j < n; j++) {
                        double acc = 0.0;
                        for (int i = 0; i < n; i++) {
                            acc += v[i] * E.get(i, j);
                        }
                        w[j] = acc;
                    }
                    v = w;
                }
                double survival = 0.0;
                for (int j = 0; j < n; j++) {
                    survival += v[j];
                }
                double F = Math.max(0.0, Math.min(1.0, 1.0 - survival));
                tvals.add(step * dT);
                Fvals.add(F);
                if (F > 1.0 - GlobalConstants.FineTol) {
                    break;
                }
            }

            Matrix curve = new Matrix(Fvals.size(), 2);
            for (int i = 0; i < Fvals.size(); i++) {
                curve.set(i, 0, Fvals.get(i));
                curve.set(i, 1, tvals.get(i));
            }
            RD.add(curve);
        }
        return RD;
    }

    /**
     * Result of {@link #getCdfFirstPassT}: the [F(t), t] curve plus the grid,
     * density and resolved state sets, the {@code out} struct of the reference.
     */
    public static class FirstPassageResult {
        /** (n x 2) matrix, first column F(t), second column t. */
        public Matrix RD;
        public double[] tset;
        public double[] density;
        public int[] source;
        public int[] target;
        public double runtime;
    }

    /**
     * Distribution of the FIRST PASSAGE TIME from state set A into state set B,
     * on the CTMC underlying this model.
     *
     * <p>Mirrors {@code @SolverCTMC/getCdfFirstPassT.m}. A and B name states
     * either as 1-based ROW INDICES into the state space returned by
     * {@link #getStateSpace()}, or as matrices of state rows, which are resolved
     * against that space. An empty A starts from the conditional stationary law
     * on the complement of B.</p>
     *
     * <p>THIS IS NOT getCdfRespT. That getter times a tagged job between an
     * arrival at a station and its departure, through the event filtration; this
     * one times the chain between two sets of states the caller names, and
     * answers questions the filtration cannot express -- the writer cycle time
     * of a readers-writers model, the time to fill a buffer, the time to leave a
     * degraded region.</p>
     *
     * <p>Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time
     * Distributions in Large Markov Chains", 2002.</p>
     *
     * @param A source state set, or null/empty for the conditional stationary law
     * @param B target state set, which may not be empty
     * @return the [F(t), t] curve and its supporting data
     */
    public FirstPassageResult getCdfFirstPassT(Matrix A, Matrix B) {
        long startTime = System.nanoTime();
        generatorResult g = this.getGenerator();
        Matrix Q = g.infGen;
        int n = Q.getNumRows();

        int[] Bidx = resolveStateSet(B, n, "B");
        if (Bidx.length == 0) {
            throw new RuntimeException("The target state set B is empty: a first passage time "
                    + "into no state is undefined.");
        }
        int[] Aidx = resolveStateSet(A, n, "A");

        String method = "expm";
        if (this.options.config != null && this.options.config.passage_method != null
                && !this.options.config.passage_method.isEmpty()) {
            method = this.options.config.passage_method;
        }

        Matrix pi0 = null;
        if (Aidx.length > 0) {
            pi0 = new Matrix(1, n);
            for (int idx : Aidx) {
                pi0.set(0, idx, 1.0 / Aidx.length);
            }
        }

        // The horizon is chosen the way the response-time getter chooses it: 100
        // events at the slowest rate in the chain.
        double minRate = Double.POSITIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double v = Math.abs(Q.get(i, j));
                if (v > GlobalConstants.FineTol && v < minRate) {
                    minRate = v;
                }
            }
        }
        double thor = Math.abs(100.0 / minRate);
        double[] tset = new double[1000];
        for (int i = 0; i < tset.length; i++) {
            tset[i] = thor * i / (tset.length - 1);
        }

        jline.api.mc.PassageCurve curve =
                jline.api.mc.Ctmc_passage_time.ctmc_passage_time(Q, pi0, Bidx, tset, method, "euler");

        FirstPassageResult out = new FirstPassageResult();
        out.RD = new Matrix(tset.length, 2);
        for (int i = 0; i < tset.length; i++) {
            out.RD.set(i, 0, curve.F[i]);
            out.RD.set(i, 1, tset[i]);
        }
        out.tset = tset;
        out.density = curve.f;
        out.source = Aidx;
        out.target = Bidx;
        out.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        return out;
    }

    /** Result of {@link #getFirstPassTMoments}. */
    public static class FirstPassageMomentsResult {
        /** (1 x nmax) moments for a passage started uniformly in A. */
        public Matrix m;
        /** (nstates x nmax), one row per starting state; 0 on B, Inf where B is unreachable. */
        public Matrix mall;
        public int[] source;
        public int[] target;
        public double runtime;
    }

    /**
     * Moments of order 1..nmax of the first passage time from state set A into state set B.
     *
     * <p>Mirrors {@code @SolverCTMC/getFirstPassTMoments.m}. NO TRANSFORM INVERSION AND NO TIME
     * GRID ARE INVOLVED: the moments come from Eq. 3 of Harrison and Knottenbelt (2002), one
     * linear solve per order, so they are EXACT and are not limited by the horizon a CDF would
     * have to be truncated at. That is why this getter exists beside {@link #getCdfFirstPassT}:
     * the variance or the skewness of a passage time costs nmax solves here and a numerical
     * integration of a truncated curve there.</p>
     *
     * <p>A and B name states as in {@link #getCdfFirstPassT}.</p>
     *
     * @param A    source state set, or null/empty for the conditional stationary law
     * @param B    target state set, which may not be empty
     * @param nmax highest moment order
     * @return the moment vector, the per-source moments and the resolved state sets
     */
    public FirstPassageMomentsResult getFirstPassTMoments(Matrix A, Matrix B, int nmax) {
        long startTime = System.nanoTime();
        if (nmax < 1) {
            nmax = 3;
        }
        generatorResult g = this.getGenerator();
        Matrix Q = g.infGen;
        int n = Q.getNumRows();

        int[] Bidx = resolveStateSet(B, n, "B");
        if (Bidx.length == 0) {
            throw new RuntimeException("The target state set B is empty: a first passage time "
                    + "into no state is undefined.");
        }
        int[] Aidx = resolveStateSet(A, n, "A");

        Matrix pi0 = null;
        if (Aidx.length > 0) {
            pi0 = new Matrix(1, n);
            for (int idx : Aidx) {
                pi0.set(0, idx, 1.0 / Aidx.length);
            }
        }

        jline.api.mc.PassageMomentsResult pm =
                jline.api.mc.Ctmc_passage_moments.ctmc_passage_moments(Q, pi0, Bidx, nmax);

        FirstPassageMomentsResult out = new FirstPassageMomentsResult();
        out.mall = pm.mall;
        out.m = new Matrix(1, nmax);
        for (int k = 0; k < nmax; k++) {
            out.m.set(0, k, pm.m[k]);
        }
        out.source = Aidx;
        out.target = Bidx;
        out.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        return out;
    }

    /** Moments of order 1..3, the reference's default. */
    public FirstPassageMomentsResult getFirstPassTMoments(Matrix A, Matrix B) {
        return getFirstPassTMoments(A, B, 3);
    }

    /**
     * A state set given as 1-based row indices or as state rows, resolved to
     * 0-based row indices; an unrecognised row is an error rather than a silent
     * drop, since a passage into a state that is not in the space is not a slow
     * passage but an undefined one.
     */
    private int[] resolveStateSet(Matrix S, int n, String name) {
        if (S == null || S.isEmpty()) {
            return new int[0];
        }
        boolean isIndexVector = (S.getNumRows() == 1 || S.getNumCols() == 1);
        if (isIndexVector) {
            for (int i = 0; i < S.length(); i++) {
                double v = S.get(i);
                if (v != Math.rint(v) || v < 1 || v > n) {
                    isIndexVector = false;
                    break;
                }
            }
        }
        java.util.TreeSet<Integer> idx = new java.util.TreeSet<Integer>();
        if (isIndexVector) {
            for (int i = 0; i < S.length(); i++) {
                idx.add((int) S.get(i) - 1);
            }
        } else {
            Matrix space = this.getStateSpace().stateSpace;
            for (int i = 0; i < S.getNumRows(); i++) {
                Matrix row = new Matrix(1, S.getNumCols());
                for (int j = 0; j < S.getNumCols(); j++) {
                    row.set(0, j, S.get(i, j));
                }
                int r = Matrix.matchrow(space, row);
                if (r < 0) {
                    throw new RuntimeException("A state given in set " + name
                            + " is not in the state space.");
                }
                idx.add(r);
            }
        }
        int[] out = new int[idx.size()];
        int k = 0;
        for (Integer v : idx) {
            out[k++] = v;
        }
        return out;
    }

    /**
     * True if the model carries a Fork or a Join, so its chain cannot be
     * enumerated on the population lattice.
     */
    private static boolean isForkJoinModel(NetworkStruct sn) {
        for (jline.lang.constant.NodeType nt : sn.nodetype) {
            if (nt == jline.lang.constant.NodeType.Fork || nt == jline.lang.constant.NodeType.Join) {
                return true;
            }
        }
        return false;
    }

    /**
     * The tag augmentation a fork-join model needs before its chain can be
     * enumerated, as runAnalyzer applies it; null for a model without a Fork or
     * a Join.
     *
     * A fork firing does not conserve the per-chain population, so the
     * population lattice State.spaceGenerator walks produces an EMPTY local
     * space for the Join and hence a 0x0 chain. getGenerator and getStateSpace
     * used to return that empty chain rather than an answer or an error. See
     * BUGS.md BUG-88.
     */
    private jline.lang.ModelAdapter.FJTagResult fjAugment(NetworkStruct sn) {
        if (!isForkJoinModel(sn)) {
            return null;
        }
        return jline.lang.ModelAdapter.fjtag(this.model);
    }

    /**
     * Stateful nodes of the augmented copy the cached nodeSpace map is keyed by.
     * Every fjtag call mints a fresh model, so a second augmentation would not
     * reproduce these keys and the local blocks would come back null; the owners
     * are recorded once, where the map is filled.
     */
    private List<jline.lang.nodes.StatefulNode> fjBlockOwners = null;

    public generatorResult getGenerator() {
        return getGenerator(null);
    }

    /**
     * The asymptotic variance of the time-average of a reward along a sample
     * path of this model's CTMC.
     *
     * <p>WHAT IT IS FOR. A simulation estimate of a steady-state mean has a
     * standard error that shrinks like sqrt(sigma^2/t), where sigma^2 is NOT the
     * stationary variance of the reward but its ASYMPTOTIC variance, which also
     * carries the autocorrelation of the path. That number is what says how long
     * a run has to be, and {@link jline.api.sim.SimRunlength#sim_runlength} turns
     * it into a run length for a target precision. It cannot be guessed from the
     * stationary variance: on M/M/1 the two differ by a factor that blows up
     * like (1-rho)^-2.
     *
     * @param f one reward value per CTMC state, in the state order
     *          {@link #getGenerator()} returns
     * @return the map of sim_asymvar_ctmc: mean, variance, asymptoticVariance
     * @see jline.api.sim.SimRunlength
     */
    public Map<String, Double> getAsymptoticVariance(double[] f) {
        Matrix infGen = getGenerator().infGen;
        int n = infGen.getNumRows();
        if (f == null || f.length != n) {
            throw new RuntimeException("getAsymptoticVariance: the reward vector has "
                    + (f == null ? 0 : f.length) + " entries but the generator is " + n + "x" + n);
        }
        return jline.api.sim.SimRunlength.sim_asymvar_ctmc(infGen, f);
    }

    /**
     * The same, with the reward given as a function of the STATE ROW rather than
     * as a vector; the state space is the one the generator was built from.
     *
     * @param f the reward, applied to each row of the state space
     * @return the map of sim_asymvar_ctmc
     */
    public Map<String, Double> getAsymptoticVariance(
            java.util.function.Function<double[], Double> f) {
        Matrix infGen = getGenerator().infGen;
        Matrix space = getStateSpace().stateSpace;
        int n = infGen.getNumRows();
        if (space == null || space.getNumRows() != n) {
            throw new RuntimeException("getAsymptoticVariance: the state space has "
                    + (space == null ? 0 : space.getNumRows()) + " rows but the generator is " + n
                    + "x" + n + "; pass the reward as a vector instead");
        }
        double[] fvec = new double[n];
        for (int i = 0; i < n; i++) {
            double[] row = new double[space.getNumCols()];
            for (int j = 0; j < space.getNumCols(); j++) {
                row[j] = space.get(i, j);
            }
            fvec[i] = f.apply(row);
        }
        return jline.api.sim.SimRunlength.sim_asymvar_ctmc(infGen, fvec);
    }

    /** Alias of {@link #getGenerator()}. */
    public generatorResult generator() {
        return getGenerator(null);
    }

    public generatorResult getGenerator(SolverOptions options) {
        if (options == null) {
            options = this.options;
        }
        if (isChainSolver()) {
            // Chain mode: the generator is the user-supplied one (P-I for a
            // DTMC), and transitions carry no event labels to filter on.
            chainEnsureAnalyzed();
            return new generatorResult(((CTMCResult) this.result).infGen, new MatrixCell(),
                    new HashMap<Integer, Sync>());
        }
        NetworkStruct sn = this.getStruct(this);
        jline.lang.ModelAdapter.FJTagResult fjRet = fjAugment(sn);
        if (fjRet != null) {
            sn = fjRet.fjsn;
            this.fjBlockOwners = fjRet.fjmodel.getStatefulNodes();
        }
        if ((this.result) == null || ((CTMCResult) this.result).infGen == null) {
            ResultCTMC solverCTMCResult = Solver_ctmc.solver_ctmc(sn, options);
            if (solverCTMCResult.getQ() == null || solverCTMCResult.getQ().getNumRows() == 0) {
                throw new RuntimeException("SolverCTMC generated an EMPTY infinitesimal generator for this "
                        + "model. A 0x0 generator is a solver limitation, not a valid chain: every caller that "
                        + "loops over its states would silently do nothing.");
            }
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

    /**
     * Symbolic infinitesimal generator with each event filtration normalized by its
     * minimum positive rate and scaled by a symbolic variable x1, ..., xE.
     *
     * <p>Java has no symbolic algebra engine, but the symbolic generator is linear in
     * the event symbols, so it is represented exactly by one numeric coefficient matrix
     * per event. The result can be evaluated at any symbol assignment via
     * {@link symbolicGeneratorResult#evalInfGen(double[])} and inspected entry-wise via
     * {@link symbolicGeneratorResult#getSymbolicEntry(int, int)}. The MATLAB and Python
     * wrappers rebuild native symbolic objects from the coefficient matrices.</p>
     *
     * @return symbolic generator decomposition
     */
    public symbolicGeneratorResult getSymbolicGenerator() {
        return getSymbolicGenerator(false);
    }

    /**
     * Symbolic infinitesimal generator.
     *
     * @param invertSymbol if true, each event filtration is divided by its symbol
     *                     instead of multiplied
     * @return symbolic generator decomposition
     */
    public symbolicGeneratorResult getSymbolicGenerator(boolean invertSymbol) {
        generatorResult gen = getGenerator();
        StateSpace ss = getStateSpace();
        MatrixCell F = gen.eventFilt;
        int nEvents = F.size();
        MatrixCell eventFilt = new MatrixCell();
        MatrixCell infGenTerms = new MatrixCell();
        List<String> symbols = new ArrayList<String>();
        for (int e = 0; e < nEvents; e++) {
            Matrix Fe = F.get(e);
            double minF = Double.POSITIVE_INFINITY;
            for (int i = 0; i < Fe.getNumRows(); i++) {
                for (int j = 0; j < Fe.getNumCols(); j++) {
                    double val = Fe.get(i, j);
                    if (val > 0 && val < minF) {
                        minF = val;
                    }
                }
            }
            if (Double.isInfinite(minF)) {
                // no positive rates for this event: empty entry, as in MATLAB
                eventFilt.set(e, new Matrix(0, 0));
                infGenTerms.set(e, new Matrix(0, 0));
                symbols.add(null);
            } else {
                Matrix Fnorm = new Matrix(Fe.getNumRows(), Fe.getNumCols());
                for (int i = 0; i < Fe.getNumRows(); i++) {
                    for (int j = 0; j < Fe.getNumCols(); j++) {
                        double val = Fe.get(i, j);
                        if (val != 0) {
                            Fnorm.set(i, j, val / minF);
                        }
                    }
                }
                eventFilt.set(e, Fnorm);
                // ctmc_makeinfgen is linear, so the symbolic generator is the sum of
                // the per-event terms scaled by their symbols
                infGenTerms.set(e, ctmc_makeinfgen(Fnorm));
                symbols.add("x" + (e + 1));
            }
        }
        return new symbolicGeneratorResult(eventFilt, infGenTerms, symbols, invertSymbol,
                ss.stateSpace, ss.localStateSpace, gen.ev);
    }

    /**
     * Refuses a query whose answer is a per-state probability under an ME.
     *
     * A matrix-exponential service embeds in the generator with negative off-diagonal
     * entries, so the stationary vector is a SIGNED measure: only its aggregates over each
     * phase block are probabilities. Mean measures stay exact (they are linear in that
     * vector), but a per-state or transient answer is not a probability at all, and
     * uniformization -- a Poisson mixture of powers of I + Q/lambda -- diverges on a
     * signed generator. Such queries are refused rather than answered with a number that
     * looks like a probability.
     *
     * Mirrors MATLAB @SolverCTMC/assertPhaseTypeStates.m and the native Python
     * SolverCTMC._assert_phasetype_states.
     *
     * @param what name of the query, used in the error message
     */
    protected void assertPhaseTypeStates(String what) {
        NetworkStruct sn = this.model.getStruct();
        if (sn.isph == null) {
            return;
        }
        for (java.util.Map<JobClass, Boolean> row : sn.isph.values()) {
            for (Boolean v : row.values()) {
                if (v != null && !v) {
                    throw new RuntimeException(what + " is unavailable: the model has a "
                            + "matrix-exponential (ME) service or arrival process, so the "
                            + "stationary vector of the generator is a signed measure and "
                            + "per-state probabilities and uniformization-based transients "
                            + "do not exist. Mean measures (getAvg, getAvgTable) remain exact.");
                }
            }
        }
    }

    /**
     * Symbolic stationary distribution of the CTMC, as a function of the event
     * rate symbols x1, ..., xE.
     *
     * <p>The generator is built here from {@link #getSymbolicGenerator()}, which
     * needs no computer algebra because it is linear in the symbols; solving
     * pi Q = 0 over the rational function field does, and is delegated to the
     * backend resolved by {@link SymEngines} (SageMath in a container by
     * default, see {@code options.config.symbolic}).</p>
     *
     * <p>The returned expressions are not comparable with another codebase's by
     * text: symbol numbering follows event enumeration order and the printed
     * normal form depends on the engine version. Substitute rates and compare
     * numbers instead, as {@link SymEngine#eval} does.</p>
     *
     * @return the stationary distribution, one expression per state
     * @throws RuntimeException if no symbolic backend is available or the solve fails
     */
    public SymEngine.CTMCSolution getSymbolicSolution() {
        SymEngine engine = SymEngines.resolve(options.config.symbolic);
        if (engine == null) {
            throw new RuntimeException(
                    "No symbolic backend is available. Start one with "
                            + "'docker run -d -p 8080:8080 " + SymEngines.DOCKER_IMAGE
                            + "', point " + SymEngines.URL_ENV + " at a running service, or "
                            + "set options.config.symbolic to its URL.");
        }
        return getSymbolicSolution(engine);
    }

    /**
     * Symbolic stationary distribution, computed by a given backend.
     *
     * @param engine the computer algebra backend
     * @return the stationary distribution, one expression per state
     * @throws RuntimeException if the solve fails
     */
    public SymEngine.CTMCSolution getSymbolicSolution(SymEngine engine) {
        symbolicGeneratorResult sym = getSymbolicGenerator();
        if (engine instanceof SageRestEngine) {
            ((SageRestEngine) engine).setTimeoutSeconds(options.config.symbolic_timeout);
        }
        try {
            return engine.solveCTMC(sym.toExpressionMatrix(), sym.activeSymbols());
        } catch (IOException e) {
            throw new RuntimeException("Symbolic solve failed on backend "
                    + engine.name() + ": " + e.getMessage(), e);
        }
    }

    @Override
    public ProbabilityResult getProb(int node, Matrix state) {
        assertPhaseTypeStates("getProb");
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

        if (Pnir == null || Pnir.isEmpty()) {
            return new ProbabilityResult(Double.NaN);
        }

        // Convert node index to station index for Pnir access
        // Validate bounds before accessing nodeToStation
        if (sn.nodeToStation == null || node >= sn.nodeToStation.getNumCols()) {
            // Return full Pnir if we can't determine station index
            return new ProbabilityResult(Pnir);
        }
        int station = (int) sn.nodeToStation.get(0, node);
        if (station < 0 || station >= Pnir.getNumCols()) {
            // Fall back to returning the full Pnir row if station index is out of bounds
            return new ProbabilityResult(Pnir);
        }
        return new ProbabilityResult(Pnir.get(0, station));
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
        assertPhaseTypeStates("getProbAggr");
        if (GlobalConstants.DummyMode) {
            return new ProbabilityResult(Double.NaN);
        }

        NetworkStruct sn = this.getStruct(this);
        if (node >= sn.nnodes) {
            line_error(mfilename(new Object(){}), "Node number exceeds the number of nodes in the model.");
            return new ProbabilityResult(Double.NaN);
        }

        // Convert node index to station index
        int station = (int) sn.nodeToStation.get(0, node);
        if (station >= sn.nstations) {
            line_error(mfilename(new Object(){}), "Station number exceeds the number of stations in the model.");
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
        if (isChainSolver()) {
            // Chain mode: the stationary vector of the user-supplied chain.
            chainEnsureAnalyzed();
            return new ProbabilityResult(((CTMCResult) this.result).pi);
        }
        assertPhaseTypeStates("getProbSys");
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
            line_error(mfilename(new Object(){}), "The model state is invalid.");
        }
        long T1 = System.nanoTime();
        this.result.runtime = (double) (T1 - T0) / 1000000000.0;
        return new ProbabilityResult(Pn);
    }

    /**
     * getProbSysAggr for a caller-named system state, given as one row of
     * per-class job counts per STATION, station-major.
     *
     * The model interchange carries no per-station initial state, so a
     * delegated query -- the CLI's {@code -a prob-sys-aggr}, and through it
     * lang='java' -- would otherwise be answered at the JAR's own default
     * initialization: on statepr_sys_aggr_large that reported 0.0941, the
     * probability of all four jobs at Queue1, where the caller asked about all
     * four at Queue3 (0.000348). Naming the state here is what makes the
     * delegated getter answer the caller's question.
     *
     * @param sysState (nstations x nclasses) per-class counts, or null for the
     *                 model's own state
     */
    public ProbabilityResult getProbSysAggr(Matrix sysState) {
        if (sysState != null) {
            setSystemState(sysState);
        }
        return getProbSysAggr();
    }

    /** {@link #getProbSysAggr(Matrix)} for the joint (non-aggregated) getter. */
    public ProbabilityResult getProbSys(Matrix sysState) {
        if (sysState != null) {
            setSystemState(sysState);
        }
        return getProbSys();
    }

    /**
     * Place the caller's per-class counts on the MODEL, one row per station.
     *
     * Through {@link Network#initFromMarginal}, so a station with phases,
     * buffers or a class-switch slot gets the state its own space requires and
     * not a bare count vector. Writing into {@code sn.state} here instead does
     * not survive: every getter calls {@code getStruct(true)}, whose
     * {@code getState()} re-reads each StatefulNode's own state back over the
     * map, so the query would run at the default initialization again.
     */
    private void setSystemState(Matrix sysState) {
        this.model.initFromMarginal(sysState);
        this.sn = this.getStruct(this);
    }

    @Override
    public ProbabilityResult getProbSysAggr() {
        if (isChainSolver()) {
            // Chain mode: states carry no phase dimension to aggregate over.
            return getProbSys();
        }
        assertPhaseTypeStates("getProbSysAggr");
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

    /**
     * (stations x classes) rate at which a class-r job BEGINS or RESUMES
     * holding a server at station i, i.e. pi*F*e over the START filtration.
     *
     * <p>At a lossless station with no in-service abandonment
     * getStartRate == getAvgTput + getPreemptRate, because every job starts
     * service once per entry into a server and every preemption is followed by
     * exactly one later resume or restart. At a non-preemptive station this
     * collapses to startRate == throughput.</p>
     *
     * <p>An accessor, not a MetricType: it adds no getAvgTable column.</p>
     */
    public Matrix getStartRate() {
        if (this.result == null || ((CTMCResult) this.result).startRate == null) {
            runAnalyzerChecked();
        }
        Matrix startRate = ((CTMCResult) this.result).startRate;
        if (startRate == null) {
            throw new RuntimeException("This solver run produced no START rates.");
        }
        return startRate;
    }

    /**
     * (stations x classes) rate at which a class-r job HOLDING A SERVER at
     * station i is pushed back into the buffer. Identically zero at a
     * non-preemptive station. Preempt-resume and preempt-independent stations
     * report the SAME rate: which phase the displaced job resumes in is not a
     * property of how often it is displaced.
     */
    public Matrix getPreemptRate() {
        if (this.result == null || ((CTMCResult) this.result).preemptRate == null) {
            runAnalyzerChecked();
        }
        Matrix preemptRate = ((CTMCResult) this.result).preemptRate;
        if (preemptRate == null) {
            throw new RuntimeException("This solver run produced no PREEMPT rates.");
        }
        return preemptRate;
    }

    /**
     * Filtration of a DERIVED event type, indexed [station][class]: the (s,ns)
     * entry is the rate at which the transition s -&gt; ns carries one such event
     * at that station for that class.
     *
     * <p>EVENTTYPE must be {@link EventType#START} or {@link EventType#PREEMPT}.
     * The two are not synchronizations: they are tags on the ARV and DEP arcs
     * that cause them, so they are NOT part of the event filtration
     * getGenerator returns (which pairs one-to-one with sn.sync and is summed
     * as D1) and are kept here instead.</p>
     */
    public Matrix[][] getEventFiltration(EventType eventType) {
        if (eventType != EventType.START && eventType != EventType.PREEMPT) {
            throw new IllegalArgumentException("getEventFiltration serves the derived events only (START, PREEMPT); "
                    + eventType + " is a synchronization and its filtration is the one getGenerator returns.");
        }
        if (this.result == null || ((CTMCResult) this.result).startFilt == null) {
            runAnalyzerChecked();
        }
        Matrix[][] filt = (eventType == EventType.START)
                ? ((CTMCResult) this.result).startFilt
                : ((CTMCResult) this.result).preemptFilt;
        if (filt == null) {
            throw new RuntimeException("This model produced no derived event filtration.");
        }
        return filt;
    }

    /**
     * runAnalyzer with its checked exceptions wrapped: the accessors above are
     * plain getters and a caller of getStartRate has no separate recovery for a
     * parser or I/O failure of the analyzer.
     */
    private void runAnalyzerChecked() {
        try {
            this.runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public StateSpace getStateSpace() {
        return this.getStateSpace(null);
    }

    /** Alias of {@link #getStateSpace()}. */
    public StateSpace stateSpace() {
        return this.getStateSpace(null);
    }

    public StateSpace getStateSpace(SolverOptions options) {
        if (options == null) {
            options = this.options;
        }
        if (isChainSolver()) {
            // Chain mode: the state space is the one attached to the chain, or
            // the state indices when the user supplied none, and it is local to
            // the single component the chain represents.
            chainEnsureAnalyzed();
            Matrix space = ((CTMCResult) this.result).space;
            MatrixCell local = new MatrixCell();
            local.set(0, space);
            return new StateSpace(space, local);
        }
        NetworkStruct sn = this.getStruct(this);
        boolean isFJ = isForkJoinModel(sn);

        if (isFJ) {
            // Take the chain from getGenerator rather than enumerating here: the
            // fork-occupied and join-firable states are vanishing and are removed
            // by the stochastic complement in solver_ctmc, so enumerating
            // separately would return a state space with MORE rows than the
            // generator it is meant to label. One producer, one pairing -- and
            // one augmentation, since the nodeSpace map is keyed by the augmented
            // copy's own stateful nodes.
            if (this.result == null || ((CTMCResult) this.result).space == null) {
                getGenerator(options);
            }
        } else if (this.result == null || ((CTMCResult) this.result).space == null) {
            // Get properly dimensioned cutoff matrix
            Matrix cutoffMatrix = options.getCutoffMatrix(sn.nstations, sn.nclasses);
            State.StateSpaceGeneratorResult stateSpaceGeneratorResult =
                    State.spaceGenerator(sn, cutoffMatrix, options);
            sn.space = stateSpaceGeneratorResult.ST.space;
            sn.spaceHash = stateSpaceGeneratorResult.ST.spaceHash;
            if (stateSpaceGeneratorResult.SS == null || stateSpaceGeneratorResult.SS.getNumRows() == 0) {
                throw new RuntimeException("SolverCTMC generated an EMPTY state space for this model, so "
                        + "there is no chain to return. This is a solver limitation, not an empty model: check "
                        + "that every stateful node admits a local state (a Join outside a fork-join "
                        + "augmentation does not).");
            }
            ((CTMCResult) this.result).space = stateSpaceGeneratorResult.SS;
            ((CTMCResult) this.result).nodeSpace = stateSpaceGeneratorResult.ST.space;
        }
        Matrix stateSpace = ((CTMCResult) this.result).space;
        // The stateful nodes the local blocks belong to: the augmented copy owns
        // one more of them (the Fork) than the user's model does.
        List<jline.lang.nodes.StatefulNode> blockOwners =
                (isFJ && this.fjBlockOwners != null) ? this.fjBlockOwners : this.model.getStatefulNodes();
        int shift = 0;
        MatrixCell localStateSpace = new MatrixCell();
        for (int i = 0; i < ((CTMCResult) this.result).nodeSpace.size(); i++) {
            int endCol =
                    shift + ((CTMCResult) this.result).nodeSpace.get(blockOwners.get(i)).getNumCols();
            Matrix value =
                    Matrix.extract(
                            ((CTMCResult) this.result).space,
                            0,
                            ((CTMCResult) this.result).space.getNumRows(),
                            shift,
                            endCol);
            localStateSpace.set(i, value);
            shift += ((CTMCResult) this.result).nodeSpace.get(blockOwners.get(i)).getNumCols();
        }
        return new StateSpace(stateSpace, localStateSpace);
    }

    public Matrix getStateSpaceAggr() {
        if (isChainSolver()) {
            // Chain mode: no phases, so the aggregate space is the state space.
            chainEnsureAnalyzed();
            return ((CTMCResult) this.result).space;
        }
        SolverOptions options = this.getOptions();
        // The aggregate space is derivable from the model alone, exactly like the
        // space in getStateSpace, so build it on demand rather than returning null
        // when nothing has been cached. Matches MATLAB and native Python.
        if (options.force || this.result == null || ((CTMCResult) this.result).spaceAggr == null) {
            try {
                this.runAnalyzer();
            } catch (Exception e) {
                line_warning(mfilename(new Object(){}),
                        "Failed to run analyzer automatically: " + e.getMessage());
                return null;
            }
        }
        return ((CTMCResult) this.result).spaceAggr;
    }

    public NetworkStruct getStruct(SolverCTMC solverCTMC) {
        //    return new NetworkStruct();
        return this.model.getStruct(true);
    }


    public ProbabilityResult getTranProb(StatefulNode node) {
        assertPhaseTypeStates("getTranProb");
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
            line_error(mfilename(new Object(){}), "Failed to compute transient probabilities: " + e.getMessage());
            return new ProbabilityResult();
        }
    }

    public ProbabilityResult getTranProbAggr(StatefulNode node) {
        assertPhaseTypeStates("getTranProbAggr");
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
            line_error(mfilename(new Object(){}), "Failed to compute aggregated transient probabilities: " + e.getMessage());
            return new ProbabilityResult();
        }
    }

    public ProbabilityResult getTranProbSys() {
        if (isChainSolver()) {
            // Chain mode: integrate (CTMC) or iterate (DTMC) from options.init_sol,
            // defaulting to the uniform distribution.
            long TC0 = System.nanoTime();
            Pair<Matrix, Matrix> tran = chainTranProbSys();
            CTMCResult chainRes = (CTMCResult) this.result;
            chainRes.tranProbSys.t = tran.getLeft();
            chainRes.tranProbSys.pit = tran.getRight();
            chainRes.tranProbSys.stateSpace = chainStateSpace();
            chainRes.runtime = (System.nanoTime() - TC0) / 1000000000.0;
            return new ProbabilityResult(tran.getRight());
        }
        assertPhaseTypeStates("getTranProbSys");
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
            line_error(mfilename(new Object(){}), "Failed to compute system transient probabilities: " + e.getMessage());
            return new ProbabilityResult();
        }
    }

    public ProbabilityResult getTranProbSysAggr() {
        assertPhaseTypeStates("getTranProbSysAggr");
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
            line_error(mfilename(new Object(){}), "Failed to compute system aggregated transient probabilities: " + e.getMessage());
            return new ProbabilityResult();
        }
    }
    @Override
    public boolean supportsTransientAnalysis() {
        // Transient averages are available (uniformization of the generator over options.timespan).
        return true;
    }



    @Override
    public void runAnalyzer()
            throws IllegalAccessException, ParserConfigurationException, IOException {
        // see _kb/06-solver-catalog.md for rationale
        jline.util.LineTimeout.set(this.getOptions().timeout);
        try {
            runAnalyzerBody();
        } finally {
            jline.util.LineTimeout.clear();
        }
    }

    /** {@code options.config.chain_aggregation}, absent or false by default. */
    private boolean chainAggregationRequested() {
        Object v = options.config.get("chain_aggregation");
        return v instanceof Boolean && (Boolean) v;
    }

    /**
     * Solves the CHAIN-AGGREGATED model and maps its metrics back to the classes.
     *
     * <p>ModelAdapter.aggregateChains collapses every chain onto a single class,
     * class switching disappearing with it, and SnDeaggregateChainResults maps
     * chain-level metrics back through alpha, the per-station share of the
     * chain's visits each class carries. What is traded is exactness on a
     * non-product-form model: one aggregate service law replaces the per-class
     * ones. A caller who needs the exact multiclass answer leaves the flag off
     * and pays the state space.
     *
     * @param T0 the wall-clock marker of the enclosing analyzer
     */
    private void runChainAggregationAnalyzer(long T0)
            throws IllegalAccessException, ParserConfigurationException, IOException {
        // Driven by TransformSolve, so the aggregate is solved by an instance of
        // THIS solver rather than a hard-wired SolverCTMC. Clearing the flag
        // states that the aggregate must not be re-aggregated, rather than
        // relying on its nchains == nclasses guard to decline it.
        runTransformAnalyzer(T0, jline.solvers.tr.TransformMethod.CHAINS, "chainaggr");
    }

    /**
     * Runs whichever transformation METHOD NAME names, publishing under LABEL.
     *
     * <p>LABEL is separate from the token so the older
     * {@code options.config.chain_aggregation} entry keeps reporting
     * {@code /chainaggr} and stays in step with the MATLAB, python and C++
     * twins, while a user-supplied token reports itself.
     */
    private void runTransformAnalyzer(long T0, String token, String label)
            throws IllegalAccessException, ParserConfigurationException, IOException {
        SolverOptions sub = options.copy();
        sub.config.put("chain_aggregation", Boolean.FALSE);
        sub.config.put("transform", token);
        jline.solvers.tr.TransformSolve.Result tr = jline.solvers.tr.TransformSolve.run(
                this.model, sn, sub, new jline.solvers.tr.TransformSolve.InnerSolve() {
                    @Override
                    public jline.solvers.tr.TransformSolve.Inner solve(jline.lang.Network submodel,
                                                                       SolverOptions opts) {
                        SolverCTMC inner = new SolverCTMC(submodel, opts);
                        try {
                            inner.runAnalyzer();
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                        String innerMethod = (inner.result instanceof CTMCResult)
                                ? ((CTMCResult) inner.result).method : "";
                        return new jline.solvers.tr.TransformSolve.Inner(
                                inner.getAvgQLen(), inner.getAvgUtil(), inner.getAvgRespT(),
                                inner.getAvgTput(), inner.getAvgSysTput(), Double.NaN, innerMethod);
                    }
                });

        double runtime = (System.nanoTime() - T0) / 1000000000.0;
        String reported = options.method + "/" + label;
        ((CTMCResult) this.result).method = reported;
        AvgHandle TH = getAvgTputHandles();
        Matrix AN = snGetArvRFromTput(sn, tr.T, TH);
        this.setAvgResults(tr.Q, tr.U, tr.R, tr.T, AN, new Matrix(0, 0),
                tr.C, tr.X, runtime, reported, tr.iter);
    }

    /** {@code options.config.fes_stations}, the 0-based stations to collapse, or null. */
    private int[] fesStationsRequested() {
        Object v = options.config.get("fes_stations");
        if (v instanceof int[]) {
            int[] a = (int[]) v;
            return a.length == 0 ? null : a;
        }
        if (v instanceof Matrix) {
            Matrix m = (Matrix) v;
            if (m.length() == 0) return null;
            int[] a = new int[m.length()];
            for (int i = 0; i < a.length; i++) a[i] = (int) m.get(i);
            return a;
        }
        return null;
    }

    /**
     * Solves the model with a station subset replaced by a FLOW-EQUIVALENT SERVER.
     *
     * <p>ModelAdapter.aggregateFES has existed in all four codebases with no
     * solver consumer at all: it was exercised by examples and tests only, so
     * nothing in the solver stack depended on it. Flow-equivalent aggregation is
     * the standard route to HIERARCHICAL DECOMPOSITION -- a subnetwork is solved
     * in isolation and enters the outer chain as a single load-dependent station,
     * which is what makes an otherwise intractable state space tractable.</p>
     *
     * <p>The reduced model answers for the surviving stations directly. For a
     * collapsed station the answer is the Chandy-Herzog-Woo conditional sum
     * E[Q_i] = sum_n P(N_fes = n) * Q_i(n), with P read off the reduced chain's
     * stationary law and Q_i(n) from the isolated subnetwork. Throughput needs no
     * conditioning: flow is fixed by the routing and an exact reduction leaves
     * the chain throughput unchanged.</p>
     *
     * @param T0 the wall-clock marker of the enclosing analyzer
     */
    private void runFesAggregationAnalyzer(long T0)
            throws IllegalAccessException, ParserConfigurationException, IOException {
        int[] subset = fesStationsRequested();
        int M = sn.nstations;
        int K = sn.nclasses;
        if (subset.length < 2) {
            throw new RuntimeException("options.config.fes_stations must name at least two "
                    + "stations: collapsing one station into a flow-equivalent server saves nothing.");
        }
        if (subset.length >= M) {
            throw new RuntimeException("options.config.fes_stations names every station: there "
                    + "is no complement left to solve.");
        }
        java.util.List<jline.lang.nodes.Station> stationSubset =
                new java.util.ArrayList<jline.lang.nodes.Station>();
        for (int i : subset) {
            if (i < 0 || i >= M) {
                throw new RuntimeException("options.config.fes_stations must be 0-based station "
                        + "indices in 0.." + (M - 1) + ".");
            }
            stationSubset.add(this.model.getStations().get(i));
        }
        jline.api.fes.FESResult fes =
                jline.api.fes.FESAggregator.aggregateFES(this.model, stationSubset);
        jline.api.fes.FESDeaggInfo info = fes.getDeaggInfo();

        SolverOptions sub = options.copy();
        sub.config.put("fes_stations", new int[0]);
        SolverCTMC inner = new SolverCTMC(fes.getFesModel(), sub);
        inner.runAnalyzer();
        Matrix Qr = inner.getAvgQLen();
        Matrix Ur = inner.getAvgUtil();
        Matrix Tr = inner.getAvgTput();
        Matrix Xr = inner.getAvgSysTput();
        Matrix pi = ((CTMCResult) inner.result).pi;
        Matrix SSq = inner.getStateSpaceAggr();

        // P(N_fes = n): the aggregate state space carries K columns per stateful
        // node, so the FES's block is the one at its stateful index.
        NetworkStruct snRed = fes.getFesModel().getStruct(true);
        int fesIst = (int) snRed.nodeToStation.get(info.fesNodeIdx);
        int fesIsf = (int) snRed.nodeToStateful.get(info.fesNodeIdx);
        Matrix cutoffs = info.cutoffs;
        int tableSize = 1;
        for (int k = 0; k < K; k++) tableSize *= ((int) cutoffs.get(k) + 1);
        double[] Pn = new double[tableSize];
        for (int r = 0; r < SSq.getNumRows(); r++) {
            Matrix nvec = new Matrix(1, K);
            for (int k = 0; k < K; k++) nvec.set(0, k, SSq.get(r, fesIsf * K + k));
            Pn[jline.api.pfqn.ld.Ljd.ljd_linearize(nvec, cutoffs)] += pi.get(r);
        }

        jline.api.fes.FESAggregator.ConditionalMetrics cm =
                jline.api.fes.FESAggregator.computeConditionalMetrics(info.isolatedModel, cutoffs);

        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);
        for (int a = 0; a < info.complementIndices.length; a++) {
            int i = info.complementIndices[a];
            for (int k = 0; k < K; k++) {
                QN.set(i, k, Qr.get(a, k));
                UN.set(i, k, Ur.get(a, k));
                TN.set(i, k, Tr.get(a, k));
            }
        }
        int Msub = info.subsetIndices.length;
        double[][] Qsub = new double[Msub][K];
        double[][] Usub = new double[Msub][K];
        for (int idx = 0; idx < tableSize; idx++) {
            if (Pn[idx] <= 0) continue;
            Matrix q = cm.QN.get(idx);
            Matrix u = cm.UN.get(idx);
            for (int a = 0; a < Msub && a < q.getNumRows(); a++)
                for (int k = 0; k < K && k < q.getNumCols(); k++) {
                    Qsub[a][k] += Pn[idx] * q.get(a, k);
                    Usub[a][k] += Pn[idx] * u.get(a, k);
                }
        }
        for (int a = 0; a < Msub; a++) {
            int i = info.subsetIndices[a];
            for (int k = 0; k < K; k++) {
                QN.set(i, k, Qsub[a][k]);
                UN.set(i, k, Usub[a][k]);
                // Flow through a station is fixed by the routing, so it is the
                // FES's throughput scaled by the ratio of ORIGINAL visit ratios.
                TN.set(i, k, Tr.get(fesIst, k) * fesVisitRatio(snRed, i, fesIst, k));
            }
        }
        Matrix RN = new Matrix(M, K);
        Matrix CN = new Matrix(1, K);
        for (int i = 0; i < M; i++)
            for (int k = 0; k < K; k++) {
                double r = TN.get(i, k) > 0 ? QN.get(i, k) / TN.get(i, k) : 0.0;
                RN.set(i, k, r);
                CN.set(0, k, CN.get(0, k) + r);
            }

        double runtime = (System.nanoTime() - T0) / 1000000000.0;
        ((CTMCResult) this.result).method = options.method + "/fes";
        AvgHandle TH = getAvgTputHandles();
        Matrix AN = snGetArvRFromTput(sn, TN, TH);
        this.setAvgResults(QN, UN, RN, TN, AN, new Matrix(0, 0), CN, Xr, runtime,
                options.method + "/fes", 0);
    }

    /** Visits at ORIGINAL station {@code ist} per visit to the FES, for class k. */
    private double fesVisitRatio(NetworkStruct snRed, int ist, int fesIst, int k) {
        double out = 0.0;
        for (int c = 0; c < sn.nchains; c++) {
            Matrix V = sn.visits.get(c);
            Matrix Vr = snRed.visits.get(c);
            if (V == null || Vr == null) continue;
            int isf = (int) sn.stationToStateful.get(ist);
            int isfFes = (int) snRed.stationToStateful.get(fesIst);
            if (isf < V.getNumRows() && isfFes < Vr.getNumRows() && Vr.get(isfFes, k) > 0) {
                out += V.get(isf, k) / Vr.get(isfFes, k);
            }
        }
        return out;
    }

    private void runAnalyzerBody()
            throws IllegalAccessException, ParserConfigurationException, IOException {
        long T0 = System.nanoTime();
        options = this.getOptions();

        // Chain mode: the generator is user-supplied, so there is no state space
        // to generate and no performance metric to derive, only the stationary
        // vector.
        if (isChainSolver()) {
            chainRunAnalyzer();
            return;
        }

        if (!isInf(options.timespan[0]) && options.timespan[0] == options.timespan[1]) {
            line_error(mfilename(new Object(){}), String.format(
                    "%s: timespan is a single point, spacing by options.tol (%e).\n",
                    this.getClass().getSimpleName(), options.tol));
            options.timespan[1] = options.timespan[0] + options.tol;
        }

        // see _kb/06-solver-catalog.md for rationale
        if (this.enableChecks && !this.supports(this.model)) {
            throw new RuntimeException("This model contains features not supported by SolverCTMC.");
        }

        // see _kb/06-solver-catalog.md for rationale
        GlobalConstants.Verbose = options.verbose;
        this.resetRandomGeneratorSeed(options.seed);

        // see _kb/06-solver-catalog.md for rationale
        for (int ci = 0; ci < this.model.getNodes().size(); ci++) {
            if (this.model.getNodes().get(ci) instanceof Cache) {
                Cache cacheNode = (Cache) this.model.getNodes().get(ci);
                cacheNode.setResultHitProb(new Matrix(0, 0));
                cacheNode.setResultMissProb(new Matrix(0, 0));
                cacheNode.setResultDelayedHitProb(new Matrix(0, 0));
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        this.model.resetStruct();
        sn = getStruct(this);

        // Chain aggregation, opt-in through options.config.chain_aggregation. The
        // state space of a multiclass model grows with the per-class populations,
        // so collapsing every chain onto a single class is the standard way to
        // make an otherwise intractable model solvable. ModelAdapter.aggregateChains
        // builds the collapsed model and SnDeaggregateChainResults maps its metrics
        // back, both of which existed with no solver consumer until this branch.
        // EXACT on a product-form model, an approximation otherwise: one aggregate
        // service law, fitted to the alpha-weighted first two moments, replaces the
        // per-class ones.
        // A user-supplied transform method name runs whichever strategy it names. The
        // inner solve carries transform='none', so a transformed submodel
        // cannot re-enter the driver.
        String trToken = jline.solvers.tr.TransformMethod.canonical(options.config.get("transform"));
        if (!jline.solvers.tr.TransformMethod.NONE.equals(trToken)) {
            runTransformAnalyzer(T0, trToken, trToken);
            return;
        }

        if (chainAggregationRequested() && sn.nchains < sn.nclasses) {
            runChainAggregationAnalyzer(T0);
            return;
        }

        // Flow-equivalent server aggregation, opt-in through
        // options.config.fes_stations. ModelAdapter.aggregateFES collapses the
        // named station subset into one load-dependent station and the collapsed
        // stations' own metrics are recovered by conditioning on its population.
        // Exact when the subnetwork is product-form.
        if (fesStationsRequested() != null) {
            runFesAggregationAnalyzer(T0);
            return;
        }

        // The 'mdd' method never builds the explicit generator, so it returns
        // before the state-space path below and leaves the state space empty by
        // design.
        if (options.method != null && options.method.equalsIgnoreCase("mdd")) {
            line_debug(options.verbose, "CTMC: using MDD level aggregation");
            jline.solvers.ctmc.analyzers.Solver_ctmc_mdd_analyzer.MddResult mddRes =
                    jline.solvers.ctmc.analyzers.Solver_ctmc_mdd_analyzer.solver_ctmc_mdd(sn, options, this.model);
            double mddRuntime = (System.nanoTime() - T0) / 1000000000.0;
            ((CTMCResult) this.result).method = options.method;
            AvgHandle mddT = getAvgTputHandles();
            Matrix mddAN = snGetArvRFromTput(sn, mddRes.TN, mddT);
            this.setAvgResults(mddRes.QN, mddRes.UN, mddRes.RN, mddRes.TN, mddAN,
                    new Matrix(0, 0), mddRes.CN, mddRes.XN, mddRuntime, options.method, 0);
            return;
        }

        // Perfect sampling replaces enumeration: intercepted before the state space
        // is built, so the memory gate below never applies to it.
        if (options.method != null && options.method.toLowerCase().startsWith("cftp")) {
            line_debug(options.verbose, String.format("CTMC: using perfect sampling (%s), %d samples",
                    options.method, (int) options.samples));
            jline.solvers.ctmc.analyzers.Solver_ctmc_cftp_analyzer.CftpResult cftp =
                    jline.solvers.ctmc.analyzers.Solver_ctmc_cftp_analyzer.solver_ctmc_cftp(sn, options);
            double cftpRuntime = (System.nanoTime() - T0) / 1000000000.0;
            ((CTMCResult) this.result).space = cftp.spaceAggr;
            ((CTMCResult) this.result).spaceAggr = cftp.spaceAggr;
            ((CTMCResult) this.result).pi = cftp.pi;
            ((CTMCResult) this.result).cftpSamples = cftp.samples;
            ((CTMCResult) this.result).cftpHorizon = cftp.horizon;
            ((CTMCResult) this.result).method = options.method;
            AvgHandle cftpT = getAvgTputHandles();
            Matrix cftpAN = snGetArvRFromTput(sn, cftp.TN, cftpT);
            this.setAvgResults(cftp.QN, cftp.UN, cftp.RN, cftp.TN, cftpAN, new Matrix(0, 0),
                    cftp.CN, cftp.XN, cftpRuntime, options.method, 0);
            return;
        }

        // Native fork-join support: solve the tag-augmented copy exactly and
        // fold the auxiliary sibling classes back into the original classes
        boolean isFJ = false;
        for (jline.lang.constant.NodeType nt : sn.nodetype) {
            if (nt == jline.lang.constant.NodeType.Fork || nt == jline.lang.constant.NodeType.Join) {
                isFJ = true;
                break;
            }
        }
        jline.solvers.tr.FJTagTransform.Context fjctx = null;
        if (isFJ) {
            if (!Double.isInfinite(options.timespan[0])) {
                throw new RuntimeException("Transient analysis of fork-join models is not supported by SolverCTMC.");
            }
            if (options.method != null && options.method.startsWith("qrf")) {
                line_warning(mfilename(new Object(){}), "The qrf method does not support fork-join models, switching to the default method.");
                options.method = "default";
            }
            fjctx = jline.solvers.tr.FJTagTransform.expand(this.model, sn, options.verbose, "CTMC");
            sn = fjctx.fjsn;
        }

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
                line_warning(mfilename(new Object(){}), String.format(
                        "%s: The model has open chains, it is recommended to specify a finite cutoff value, e.g., SolverCTMC(model).cutoff(1).",
                        this.getClass().getSimpleName()));
                double autoCutoff = Math.ceil(Math.pow(6000, 1.0 / (M * K)));
                options.cutoff(autoCutoff);
                line_warning(mfilename(new Object(){}), String.format(
                        "%s: Setting cutoff=%d.", this.getClass().getSimpleName(), (int) autoCutoff));
            }
            // Mandatory truncation warning for open/mixed models
            Matrix currentCutoffFinal = options.getCutoffMatrix(sn.nstations, sn.nclasses);
            line_warning(mfilename(new Object(){}), String.format(
                    "CTMC solver using state space cutoff = %d for open/mixed model. State space truncation may cause inaccurate results. Consider varying cutoff to assess sensitivity.",
                    (int) currentCutoffFinal.get(0, 0)));
        }

        // see _kb/06-solver-catalog.md for rationale. The estimate is taken on
        // the PH-converted struct, as MATLAB and Python do: a non-Markovian
        // service becomes phases the raw struct does not carry.
        double logNstates = MemoryGuard.stateSpaceLogSize(
                SnNonmarkovToPh.snNonmarkovToPh(sn, options, false), options);
        line_debug(options.verbose, String.format("State space size estimate: exp(%f)", logNstates));
        MemoryGuard.GateResult gateRes = MemoryGuard.gate(logNstates, options.force,
                options.verbose != null && options.verbose == VerboseLevel.DEBUG,
                MemoryGuard.DEFAULT_SAFETY_FRACTION);
        if (!gateRes.ok) {
            line_error(mfilename(new Object(){}),
                    gateRes.message + " Stopping SolverCTMC.");
            return;
        }

        if (Double.isInfinite(options.timespan[0])) {
            // QRF (Quadratic/Linear Reduction Framework) LP-based bounds have
            // moved out of SolverCTMC into the dedicated SolverBA solver. The
            // text lives in unsupportedMethodReason, which checkDeclaredMethod
            // asks BEFORE it reports an unlisted method; this call is what
            // still refuses on the enableChecks = false path, which skips that
            // gate entirely.
            String movedQrf = unsupportedMethodReason(options.method);
            if (!movedQrf.isEmpty()) {
                throw new RuntimeException(movedQrf);
            }

            Map<StatefulNode, Matrix> s0 = sn.state;
            Map<StatefulNode, Matrix> s0prior = sn.stateprior;
            for (int ind = 1; ind < sn.nnodes; ind++) {
                if (sn.isstateful.get(ind, 0) == 1) {
                    int isf = (int) sn.nodeToStateful.get(ind);
                    // see _kb/06-solver-catalog.md for rationale
                    StatefulNode sfNode = sn.stateful.get(isf);
                    Matrix initialStatePrior = (s0prior != null) ? s0prior.get(sfNode) : null;
                    int maxPos = (initialStatePrior != null && !initialStatePrior.isEmpty()) ? Maths.maxpos(initialStatePrior) : 0;
                    Matrix stateMatrix = s0.get(sfNode);
                    if (stateMatrix != null && !stateMatrix.isEmpty()) {
                        Matrix initialState = stateMatrix.getRow(maxPos);
                        sn.state.replace(sfNode, initialState);
                    }
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

            // see _kb/06-solver-catalog.md for rationale
            Matrix rtOrig = sn.rt.copy();
            if (!isFJ) {
                // (skipped on fork-join models: the analyzed struct is the
                // tag-augmented copy, whose states do not fit the original model)
                for (int isf = 0; isf < sn.nstateful; isf++) {
                    int ind = (int) sn.statefulToNode.get(isf);
                    // use stateful nodes directly to avoid IndexOutOfBoundsException
                    Node a = this.model.getStatefulNodes().get(isf);
                    ((StatefulNode) a).setState(sn.state.get(this.model.getStatefulNodes().get(isf)));
                    if (a instanceof Cache) {
                        Cache cacheNode = (Cache) a;
                        CacheNodeParam _cnp = (CacheNodeParam) sn.nodeparam.get(this.model.getNodes().get(ind));
                        cacheNode.setResultHitProb(_cnp.actualhitprob);
                        cacheNode.setResultMissProb(_cnp.actualmissprob);
                        cacheNode.setResultResidT(_cnp.actualresidt);
                        computeCacheItemProb(sn, ind, isf, SS, analyzerResult.pi, cacheNode);
                        computeDelayedHitQLen(sn, ind, isf, SS, analyzerResult.pi, cacheNode);
                        applyDelayedHitSplit(sn, ind, isf, SS, analyzerResult.pi, Q, cacheNode, _cnp);
                        this.model.refreshChains(true);
                    }
                }
            }
            ((CTMCResult) this.result).infGen = Q.copy();
            ((CTMCResult) this.result).space = SS;
            ((CTMCResult) this.result).spaceAggr = SSq;
            ((CTMCResult) this.result).pi = analyzerResult.pi;
            ((CTMCResult) this.result).spaceWork = analyzerResult.StateSpaceWork;
            ((CTMCResult) this.result).spaceAggrWork = analyzerResult.StateSpaceAggrWork;
            ((CTMCResult) this.result).infGenWork = analyzerResult.InfGenWork;
            ((CTMCResult) this.result).nodeSpace = sn.space;
            ((CTMCResult) this.result).eventFilt = Dfilt;
            // Derived START/PREEMPT filtration and rates, in their own fields:
            // eventFilt is paired with sn.sync one-to-one and summed as D1, so
            // a filtration that rides on the same arcs must not join it.
            ((CTMCResult) this.result).startFilt = analyzerResult.startFilt;
            ((CTMCResult) this.result).preemptFilt = analyzerResult.preemptFilt;
            ((CTMCResult) this.result).startRate = analyzerResult.StartN;
            ((CTMCResult) this.result).preemptRate = analyzerResult.PreemptN;
            // Set the method field in the result
            ((CTMCResult) this.result).method = options.method;
            double runtime = (System.nanoTime() - T0) / 1000000000.0;
            sn.space = new HashMap<>();
            M = sn.nstations;
            int R = sn.nclasses;
            AvgHandle T = getAvgTputHandles();
            Matrix AN;
            if (isFJ) {
                // fold the auxiliary sibling classes back into the original
                // classes and report the Join per-sibling waiting time
                jline.solvers.tr.FJTagTransform.Lifted lifted =
                        jline.solvers.tr.FJTagTransform.lift(fjctx, QN, UN, RN, TN, CN, XN, T);
                QN = lifted.QN;
                UN = lifted.UN;
                RN = lifted.RN;
                TN = lifted.TN;
                CN = lifted.CN;
                XN = lifted.XN;
                AN = lifted.AN;
                // restore the original struct: downstream consumers (tables,
                // handles) index the folded matrices by the original classes
                sn = fjctx.snOrig;
            } else {
                Matrix rtRefreshed = sn.rt;
                sn.rt = rtOrig;
                jline.api.sn.SnPnAvgRates.snPnAvgRates(sn, QN, TN, null, RN);
                AN = snGetArvRFromTput(sn, TN, T);
                if (rtRefreshed != null) {
                    sn.rt = rtRefreshed;
                }
            }
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

            // see _kb/06-solver-catalog.md for rationale
            double[][] rowCountsData = new double[sn.nstateful][1];
            for (int isf = 0; isf < sn.nstateful; isf++) {
                StatefulNode node = this.model.getStatefulNodes().get(isf);
                rowCountsData[isf][0] = s0.get(node).getNumRows();
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
                        // s0_id is 0-based; MATLAB uses 1+s0_id for 1-based indexing, Java uses s0_id directly
                        s0prior_val =
                                s0prior_val * s0prior.get(this.model.getStatefulNodes().get(isf)).get((int) s0_id.get(isf));

                        // Extract row from the state matrix corresponding to the current state
                        Matrix newState =
                                Matrix.extractRows(
                                        s0.get(this.model.getStatefulNodes().get(isf)),
                                        (int) s0_id.get(isf),
                                        (int) s0_id.get(isf) + 1,
                                        null);

                        // see _kb/06-solver-catalog.md for rationale
                        this.model.getStatefulNodes().get(isf).setState(newState);
                    }
                }
                NetworkStruct sn_cur = this.model.getStruct(true);
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
                                // Compute sorted unique union of time values (MATLAB union semantics)
                                java.util.TreeSet<Double> timeSet = new java.util.TreeSet<>();
                                for (int ti = 0; ti < existingTimes.getNumRows(); ti++) {
                                    timeSet.add(existingTimes.get(ti, 0));
                                }
                                for (int ti = 0; ti < newTimes.getNumRows(); ti++) {
                                    timeSet.add(newTimes.get(ti, 0));
                                }
                                Matrix tunion = new Matrix(timeSet.size(), 1);
                                int tuIdx = 0;
                                for (Double tv : timeSet) {
                                    tunion.set(tuIdx++, 0, tv);
                                }

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

            // Populate result.QNt, result.UNt, result.TNt, result.t from Tran.Avg
            // so that SolverENV.finish() can read them (it uses these flat arrays)
            if (((CTMCResult) this.result).Tran != null
                    && ((CTMCResult) this.result).Tran.Avg != null
                    && ((CTMCResult) this.result).Tran.Avg.Q != null
                    && !((CTMCResult) this.result).Tran.Avg.Q.isEmpty()) {
                // Extract time points from first available entry
                Matrix firstEntry = ((CTMCResult) this.result).Tran.Avg.Q.values().iterator().next()
                        .values().iterator().next();
                int nTimePoints = firstEntry.getNumRows();
                this.result.t = new Matrix(nTimePoints, 1);
                for (int ti = 0; ti < nTimePoints; ti++) {
                    this.result.t.set(ti, 0, firstEntry.get(ti, 1));
                }

                this.result.QNt = new Matrix[M][K];
                this.result.UNt = new Matrix[M][K];
                this.result.TNt = new Matrix[M][K];
                for (int ist = 0; ist < M; ist++) {
                    for (int r = 0; r < K; r++) {
                        this.result.QNt[ist][r] = Matrix.extractColumns(
                                ((CTMCResult) this.result).Tran.Avg.Q.get(ist).get(r), 0, 1, null);
                        this.result.UNt[ist][r] = Matrix.extractColumns(
                                ((CTMCResult) this.result).Tran.Avg.U.get(ist).get(r), 0, 1, null);
                        this.result.TNt[ist][r] = Matrix.extractColumns(
                                ((CTMCResult) this.result).Tran.Avg.T.get(ist).get(r), 0, 1, null);
                    }
                }
            }
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

    public SampleResult sample(StatefulNode node, int numEvents) {
        assertPhaseTypeStates("sample");
        SolverOptions options = this.getOptions();
        options.force = true;

        if (this.result == null || ((CTMCResult) this.result).infGen == null) {
            try {
                this.runAnalyzer();
            } catch (Exception e) {
                line_error(mfilename(new Object(){}), "Failed to run analyzer: " + e.getMessage());
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

            result.t = new Matrix(numEvents, 1);
            result.state = new Matrix(numEvents, sn.space.get(this.model.getStatefulNodes().get(isf)).getNumCols());
            result.event = new ArrayList<EventInfo>();

            if (useMMAPSampling) {
                // Use MMAP sampling for better event resolution
                java.util.Random random = new java.util.Random();
                jline.io.Ret.mamMMAPSample mmapSample = jline.api.mam.Mmap_sample.mmap_sample(
                        nodeSpecificMMAP, (long)numEvents, random);

                double[] interArrivalTimes = mmapSample.getSamples();
                int[] eventTypes = mmapSample.getTypes();

                // Also get CTMC states for state information
                Ret.ctmcSimulation simulation = ctmc_simulate(infGen, pi0, numEvents);

                double currentTime = 0.0;
                for (int i = 0; i < numEvents; i++) {
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
                Ret.ctmcSimulation simulation = ctmc_simulate(infGen, pi0, numEvents);

                double currentTime = 0.0;
                for (int i = 0; i < numEvents; i++) {
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
            line_error(mfilename(new Object(){}), "CTMC sampling failed: " + e.getMessage());
            // Fallback to simple implementation
            result.t = new Matrix(numEvents, 1);
            result.state = new Matrix(numEvents, sn.space.get(this.model.getStatefulNodes().get(isf)).getNumCols());
            result.event = new ArrayList<EventInfo>();
        }

        return result;
    }

    public jline.io.Ret.SampleResult sampleSys(int numEvents) {
        if (isChainSolver()) {
            return chainSampleSys(numEvents);
        }
        assertPhaseTypeStates("sampleSys");
        SolverOptions options = this.getOptions();
        options.force = true;

        if (this.result == null || ((CTMCResult) this.result).infGen == null) {
            try {
                this.runAnalyzer();
            } catch (Exception e) {
                line_error(mfilename(new Object(){}), "Failed to run analyzer: " + e.getMessage());
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
            MatrixCell localStateSpace = stateSpaceResult.localStateSpace;
            NetworkStruct sn = this.getStruct(this);

            // see _kb/06-solver-catalog.md for rationale
            Map<StatefulNode, Matrix> initState = sn.state;
            List<Double> s0List = new ArrayList<Double>();
            List<StatefulNode> statefulNodes = this.model.getStatefulNodes();
            for (int isf = 0; isf < statefulNodes.size(); isf++) {
                StatefulNode statefulNode = statefulNodes.get(isf);
                Matrix stateMatrix = initState.get(statefulNode);
                int spaceCols = (localStateSpace != null && isf < localStateSpace.size()
                        && localStateSpace.get(isf) != null) ? localStateSpace.get(isf).getNumCols() : 0;
                int stateCols = (stateMatrix != null) ? stateMatrix.getNumCols() : 0;
                for (int j = 0; j < spaceCols - stateCols; j++) {
                    s0List.add(0.0);
                }
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
            if (matchIdx < 0) {
                // Retry with rounded state (e.g. fractional ENV-averaged states)
                Matrix roundedState = new Matrix(1, s0.getNumCols());
                for (int j = 0; j < s0.getNumCols(); j++) {
                    roundedState.set(0, j, Math.round(s0.get(0, j)));
                }
                matchIdx = Matrix.matchrow(stateSpace, roundedState);
            }
            if (matchIdx < 0) {
                throw new RuntimeException("Initial state not contained in the state space.");
            }
            pi0.set(0, matchIdx, 1.0);

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
            MMAP = jline.api.mam.Mmap_normalize.mmap_normalize(MMAP);

            // Sample MMAP (like MATLAB dev/ line 27)
            // [sjt, event, ~, ~, sts] = mmap_sample(MMAP, numEvents, pi0);
            java.util.Random random = new java.util.Random(options.seed);
            jline.io.Ret.mamMMAPSample mmapSample = jline.api.mam.Mmap_sample.mmap_sample(MMAP, (long)numEvents, random);

            double[] interArrivalTimes = mmapSample.getSamples();
            int[] eventTypes = mmapSample.getTypes();
            int[] mmapStates = mmapSample.getStates();

            // Build time series (like MATLAB dev/ line 32)
            // MATLAB: tranSysState.t = cumsum([0,sjt(1:end-1)']');
            Matrix t = new Matrix(numEvents, 1);
            double cumulativeTime = 0.0;
            for (int i = 0; i < numEvents; i++) {
                t.set(i, 0, cumulativeTime);
                if (i < numEvents - 1) {
                    cumulativeTime += interArrivalTimes[i];
                }
            }

            // Build state series from MMAP states (like MATLAB dev/ lines 33-36)
            // MATLAB: tranSysState.state{isf} = stateSpace(sts,(nst(isf):nst(isf+1)-1));
            Matrix state = new Matrix(numEvents, stateSpace.getNumCols());
            for (int i = 0; i < numEvents; i++) {
                int stateIdx = (mmapStates != null && i < mmapStates.length) ? mmapStates[i] : 0;
                if (stateIdx >= 0 && stateIdx < stateSpace.getNumRows()) {
                    for (int j = 0; j < stateSpace.getNumCols(); j++) {
                        state.set(i, j, stateSpace.get(stateIdx, j));
                    }
                }
            }

            // see _kb/06-solver-catalog.md for rationale
            List<Event> eventList = new ArrayList<Event>();
            for (int i = 0; i < numEvents; i++) {
                double eventTime = t.get(i, 0);
                int syncIdx = eventTypes[i];  // Direct index into synchInfo (0-indexed in Java)

                // Get synchronization for this event type
                Sync sync = synchInfo.get(syncIdx);
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

            // Convert event list to matrix format: [time, node, eventType, jobclass].
            // THE CODE TABLE LIVES IN Ret.SampleResult, which is what reads it back;
            // a switch here mapped everything but ARV/DEP/PHASE to -1, so those
            // events were sampled and then made unreadable.
            Matrix event = new Matrix(eventList.size(), 4);
            for (int i = 0; i < eventList.size(); i++) {
                Event e = eventList.get(i);
                event.set(i, 0, e.getT());
                event.set(i, 1, e.getNode());
                event.set(i, 2, e.getEvent() == null ? -1
                        : jline.io.Ret.SampleResult.codeOf(e.getEvent()));
                event.set(i, 3, e.getJobClass());
            }

            // Create sample result
            jline.io.Ret.SampleResult result = new jline.io.Ret.SampleResult("ctmc", t, state, event, false, null, numEvents);

            return result;

        } catch (Exception e) {
            line_error(mfilename(new Object(){}), "CTMC MMAP sampling failed: " + e.getMessage());
            e.printStackTrace();
            // Fallback to empty result
            StateSpace stateSpaceResult = getStateSpace();
            Matrix t = new Matrix(numEvents, 1);
            Matrix state = new Matrix(numEvents, stateSpaceResult.stateSpace.getNumCols());
            Matrix event = new Matrix(numEvents, 3);
            return new jline.io.Ret.SampleResult("ctmc", t, state, event, false, null, numEvents);
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
        Matrix carriedPi0 = solverCTMCResult.getPi0();
        if (state0 == -1 && carriedPi0 != null && carriedPi0.getNumCols() == infGen.getNumRows()) {
            // a vanishing initial state is absent from the complemented chain: the analyzer
            // starts from its first-entry distribution instead
            pi0 = carriedPi0.copy();
        } else {
            if (state0 == -1) {
                // Try rounding fractional states (e.g., from ENV solver weighted averages)
                Matrix roundedState = new Matrix(1, initialState.getNumCols());
                for (int j = 0; j < initialState.getNumCols(); j++) {
                    roundedState.set(0, j, Math.round(initialState.get(0, j)));
                }
                state0 = Matrix.matchrow(stateSpace, roundedState);
                if (state0 == -1) {
                    throw new RuntimeException("Initial state not contained in the state space.");
                }
            }
            pi0.set(0, state0, 1.0);
        }

        // Perform transient analysis
        Matrix t;
        Matrix pit;

        // see _kb/06-solver-catalog.md for rationale
        java.util.List<FluidRateMultiplier.RateEntry> rate_sched =
                (options.config != null) ? options.config.rate_sched : null;
        double[][][] mscale = null;

        if (rate_sched != null && !rate_sched.isEmpty()) {
            TimeVaryingTransient tv = ctmcTimeVarying(sn, options, infGen, pi0, rate_sched,
                    sn.nstations, sn.nclasses);
            t = tv.t;
            pit = tv.pit;
            mscale = tv.mscale;
        } else if (options.config != null && "fau".equalsIgnoreCase(options.config.transient_method)) {
            // see _kb/06-solver-catalog.md (CTMC section, transient methods)
            Pair<Matrix, Matrix> fau = ctmcFauTransient(infGen, pi0, options);
            t = fau.getLeft();
            pit = fau.getRight();
        } else {
            if (options.config != null && options.config.transient_method != null
                    && !"ode".equalsIgnoreCase(options.config.transient_method)) {
                line_error(mfilename(new Object(){}), "Unknown options.config.transient_method '"
                        + options.config.transient_method + "'; use 'ode' or 'fau'.");
            }
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
                line_error(mfilename(new Object(){}), "Transient analysis failed: " + e.getMessage());
                // Fallback to simple result
                t = new Matrix(1, 1);
                t.set(0, 0, options.timespan[1]);
                pit = new Matrix(1, infGen.getNumRows());
                for (int j = 0; j < infGen.getNumRows(); j++) {
                    pit.set(0, j, pi0.get(0, j));
                }
            }
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
                    // see _kb/06-solver-catalog.md for rationale
                    if (mscale != null) {
                        throughput *= mscale[ist][k][timeIdx];
                    }
                    TNt.set(timeIdx, ist * K + k, throughput);
                }

                // see _kb/06-solver-catalog.md for rationale
                int indSrc = (int) sn.stationToNode.get(ist);
                if (sn.nodetype.get(indSrc) == NodeType.Source) {
                    // QNt/UNt already zero-initialized for this station/class.
                    continue;
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
                SchedStrategy schedStrategy = sn.sched.get(sn.stations.get(ist));
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
     * Transient trajectory by FAST ADAPTIVE UNIFORMIZATION, selected with
     * {@code options.config.transient_method = "fau"} (see {@link jline.api.mc.Ctmc_fau}).
     *
     * <p>WHY IT IS NOT A METHOD NAME. {@code "fau"} is deliberately absent from the
     * solver's valid-method list: that list is enumerated by the sanity harness,
     * which then demands a recorded baseline per method, and this changes no
     * stationary answer at all -- it is the transient path only. The config key is
     * where this analyzer's other transient switches already live.</p>
     *
     * <p>MARCHED, NOT RESTARTED. pi(t_{k+1}) comes from pi(t_k) over the step
     * rather than from pi(0) over the whole horizon, which keeps the cost
     * proportional to the grid instead of quadratic in it. Every step removes a
     * little mass and none puts any back, so the per-step tolerance is
     * {@code fau_epsilon} divided by the number of steps and the total defect stays
     * below it; the accumulated defect is reported, not normalized away, the point
     * of the method being that its error is a measured quantity.</p>
     *
     * @return the time grid and the (ntimes x nstates) trajectory
     */
    private Pair<Matrix, Matrix> ctmcFauTransient(Matrix infGen, Matrix pi0, SolverOptions options) {
        double t0 = options.timespan[0];
        double t1 = options.timespan[1];
        if (Double.isInfinite(t1) || Double.isNaN(t1)) {
            line_error(mfilename(new Object(){}),
                    "transient_method 'fau' needs a finite horizon; options.timespan[1] is not finite.");
        }
        java.util.List<Double> grid = new java.util.ArrayList<Double>();
        if (options.timestep > 0) {
            for (double ti = t0; ti < t1; ti += options.timestep) {
                grid.add(ti);
            }
            grid.add(t1);
        } else {
            int ngrid = (options.config != null && options.config.fau_ngrid > 1)
                    ? options.config.fau_ngrid : 100;
            for (int i = 0; i < ngrid; i++) {
                grid.add(t0 + (t1 - t0) * i / (double) (ngrid - 1));
            }
        }
        int nt = grid.size();
        double epsilon = (options.config != null && options.config.fau_epsilon > 0)
                ? options.config.fau_epsilon : 1e-6;
        double delta = (options.config != null && options.config.fau_delta >= 0)
                ? options.config.fau_delta : 1e-12;
        double epsStep = epsilon / Math.max(1, nt - 1);

        int n = infGen.getNumRows();
        Matrix t = new Matrix(nt, 1);
        Matrix pit = new Matrix(nt, n);
        for (int j = 0; j < n; j++) {
            pit.set(0, j, pi0.get(0, j));
        }
        t.set(0, 0, grid.get(0));
        Matrix cur = pi0.copy();
        double defect = 0.0;
        int supportMax = 0;
        for (int k = 1; k < nt; k++) {
            double dt = grid.get(k) - grid.get(k - 1);
            Ctmc_fau.CtmcFauResult r = Ctmc_fau.ctmc_fau(cur, infGen, dt, epsStep, delta, -1);
            cur = r.pit;
            defect += r.errorBound;
            supportMax = Math.max(supportMax, r.supportMax);
            t.set(k, 0, grid.get(k));
            for (int j = 0; j < n; j++) {
                pit.set(k, j, cur.get(j));
            }
        }
        line_debug(options.verbose, "CTMC transient by FAU: " + nt + " grid points, support at most "
                + supportMax + " of " + n + " states, missing mass " + defect);
        if (defect > GlobalConstants.CoarseTol) {
            line_warning(mfilename(new Object(){}), "FAU transient discarded " + defect
                    + " of the probability mass over the horizon; tighten options.config.fau_epsilon"
                    + " or options.config.fau_delta.\n");
        }
        return new Pair<Matrix, Matrix>(t, pit);
    }

    /**
     * Result of the time-inhomogeneous transient propagation: the uniform time
     * grid, the state probability trajectory, and the per-(station,class)
     * throughput multiplier over that grid.
     */
    private static class TimeVaryingTransient {
        final Matrix t;
        final Matrix pit;
        final double[][][] mscale; // [station][class][timeIdx]

        TimeVaryingTransient(Matrix t, Matrix pit, double[][][] mscale) {
            this.t = t;
            this.pit = pit;
            this.mscale = mscale;
        }
    }

    /**
     * Integrates the time-inhomogeneous forward equation dpi/dt = pi Q(t) over a
     * uniform grid, with the generator frozen at each interval midpoint.
     *
     * <p>Mirrors MATLAB {@code local_ctmc_timevarying} in
     * {@code solver_ctmc_transient_analyzer.m}. The generator is
     * {@code Q(t) = Qbase + sum_sc (m_sc(t)-1) Qhat_sc}, where {@code Qhat_sc} is
     * the linear component of {@code Qbase} attributable to the scaled
     * (station,class) rate, extracted by a single probe rebuild (Q is linear in
     * {@code sn.rates}).</p>
     *
     * <p>The propagation deliberately uses the matrix exponential rather than
     * uniformization: an LN layer can carry a near-instantaneous reply-signal
     * sentinel rate (~1e9), making the generator stiff (q dt ~ 1e8);
     * uniformization then splits the step into more than 1e5 sub-segments and
     * leaks all probability mass to zero. {@code expm} of a valid generator is
     * exactly stochastic, so mass is conserved for any stiffness.</p>
     */
    private TimeVaryingTransient ctmcTimeVarying(NetworkStruct sn, SolverOptions options,
                                                 Matrix qbase, Matrix pi0,
                                                 java.util.List<FluidRateMultiplier.RateEntry> rate_sched,
                                                 int M, int K) {
        int ngrid = (options.config != null && options.config.ctmc_tv_ngrid > 1)
                ? options.config.ctmc_tv_ngrid : 100;
        double t0 = options.timespan[0];
        double t1 = options.timespan[1];
        Matrix t = new Matrix(ngrid, 1);
        for (int i = 0; i < ngrid; i++) {
            t.set(i, 0, t0 + (t1 - t0) * i / (double) (ngrid - 1));
        }
        int nS = qbase.getNumRows();
        int nsc = rate_sched.size();

        final double probe = 2.0;
        Matrix[] qhat = new Matrix[nsc];
        double[][] mtraj = new double[ngrid][nsc]; // [timeIdx][schedIdx]
        int[] scStation = new int[nsc];
        int[] scClass = new int[nsc];

        for (int s = 0; s < nsc; s++) {
            FluidRateMultiplier.RateEntry entry = rate_sched.get(s);
            int ist = entry.station;
            int r = entry.jobclass;
            scStation[s] = ist;
            scClass[s] = r;
            if (ist < 0 || ist >= M || r < 0 || r >= K) {
                line_error(mfilename(new Object(){}),
                        "rate_sched entry refers to a (station,class) outside the network.");
            }

            // see _kb/06-solver-catalog.md for rationale
            Station station = sn.stations.get(ist);
            JobClass jobClass = sn.jobclasses.get(r);
            double savedRate = sn.rates.get(ist, r);
            MatrixCell savedProc = (sn.proc != null && sn.proc.get(station) != null)
                    ? sn.proc.get(station).get(jobClass) : null;
            Matrix savedMu = (sn.mu != null && sn.mu.get(station) != null)
                    ? sn.mu.get(station).get(jobClass) : null;

            Matrix qp;
            try {
                sn.rates.set(ist, r, savedRate * probe);
                if (savedProc != null) {
                    MatrixCell scaledProc = new MatrixCell(savedProc.size());
                    for (int z = 0; z < savedProc.size(); z++) {
                        Matrix dz = savedProc.get(z);
                        scaledProc.set(z, dz != null ? dz.scale(probe) : null);
                    }
                    sn.proc.get(station).put(jobClass, scaledProc);
                }
                if (savedMu != null) {
                    sn.mu.get(station).put(jobClass, savedMu.scale(probe));
                }
                qp = Solver_ctmc.solver_ctmc(sn, options).getQ();
            } finally {
                sn.rates.set(ist, r, savedRate);
                if (savedProc != null) {
                    sn.proc.get(station).put(jobClass, savedProc);
                }
                if (savedMu != null) {
                    sn.mu.get(station).put(jobClass, savedMu);
                }
            }
            if (qp.getNumRows() != nS) {
                line_error(mfilename(new Object(){}),
                        "rate_sched probe changed the CTMC state-space size; cannot build time-varying generator.");
            }
            qhat[s] = qp.sub(qbase).scale(1.0 / (probe - 1.0));

            // multiplier m(t) = rate(t)/nominal (nominal defaults to sn.rates(ist,r))
            double nominal = (entry.nominal != null) ? entry.nominal.doubleValue() : savedRate;
            if (!(Math.abs(nominal) > 0.0)) {
                line_error(mfilename(new Object(){}),
                        "rate_sched entry has a zero nominal rate; cannot form the multiplier.");
            }
            for (int i = 0; i < ngrid; i++) {
                mtraj[i][s] = interpClamped(entry.tgrid, entry.rates, t.get(i, 0)) / nominal;
            }
        }

        // see _kb/06-solver-catalog.md for rationale
        Matrix pit = new Matrix(ngrid, nS);
        double[] cur = new double[nS];
        for (int j = 0; j < nS; j++) {
            cur[j] = pi0.get(0, j);
            pit.set(0, j, cur[j]);
        }
        for (int k = 0; k < ngrid - 1; k++) {
            double dt = t.get(k + 1, 0) - t.get(k, 0);
            Matrix qk = qbase.copy();
            for (int s = 0; s < nsc; s++) {
                double mk = 0.5 * (mtraj[k][s] + mtraj[k + 1][s]);
                if (mk != 1.0) {
                    qk = qk.add(1.0, qhat[s].scale(mk - 1.0));
                }
            }
            Matrix expQ = qk.scale(dt).expm();
            double[] next = new double[nS];
            for (int j = 0; j < nS; j++) {
                double acc = 0.0;
                for (int i = 0; i < nS; i++) {
                    if (cur[i] != 0.0) {
                        acc += cur[i] * expQ.get(i, j);
                    }
                }
                next[j] = acc;
                pit.set(k + 1, j, acc);
            }
            cur = next;
        }

        // Per-(station,class) throughput multiplier over time (1 where not scaled).
        double[][][] mscale = new double[M][K][ngrid];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                for (int p = 0; p < ngrid; p++) {
                    mscale[i][r][p] = 1.0;
                }
            }
        }
        for (int s = 0; s < nsc; s++) {
            for (int p = 0; p < ngrid; p++) {
                mscale[scStation[s]][scClass[s]][p] = mtraj[p][s];
            }
        }

        return new TimeVaryingTransient(t, pit, mscale);
    }

    /**
     * Clamped piecewise-linear interpolation of {@code (tg, val)} at time
     * {@code tt}; times outside the grid take the boundary value. Mirrors the
     * MATLAB {@code interp1(...,'linear')} on a clamped argument.
     */
    private static double interpClamped(double[] tg, double[] val, double tt) {
        int n = tg.length;
        if (n == 0) {
            return 0.0;
        }
        if (n == 1 || tt <= tg[0]) {
            return val[0];
        }
        if (tt >= tg[n - 1]) {
            return val[n - 1];
        }
        int lo = 0;
        int hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (tg[mid] <= tt) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double w = (tt - tg[lo]) / (tg[lo + 1] - tg[lo]);
        return (1.0 - w) * val[lo] + w * val[lo + 1];
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
            return jline.api.mam.Mmap_normalize.mmap_normalize(nodeSpecificMMAP);

        } catch (Exception e) {
            line_warning(mfilename(new Object(){}),
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

    public static class symbolicGeneratorResult {

        /** Normalized filtration matrix of each event (coefficient of its symbol); empty if the event has no positive rates */
        public MatrixCell eventFilt;
        /** Per-event contribution to the symbolic infinitesimal generator (diagonal included) */
        public MatrixCell infGenTerms;
        /** Symbol names x1..xE; null for events with no positive rates */
        public List<String> symbols;
        /** True if event filtrations are divided by their symbol instead of multiplied */
        public boolean invertSymbol;
        public Matrix stateSpace;
        public MatrixCell nodeStateSpace;
        public Map<Integer, Sync> syncInfo;

        public symbolicGeneratorResult(MatrixCell eventFilt, MatrixCell infGenTerms,
                                       List<String> symbols, boolean invertSymbol,
                                       Matrix stateSpace, MatrixCell nodeStateSpace,
                                       Map<Integer, Sync> syncInfo) {
            this.eventFilt = eventFilt;
            this.infGenTerms = infGenTerms;
            this.symbols = symbols;
            this.invertSymbol = invertSymbol;
            this.stateSpace = stateSpace;
            this.nodeStateSpace = nodeStateSpace;
            this.syncInfo = syncInfo;
        }

        /**
         * Evaluate the symbolic generator at the given symbol values.
         *
         * @param x value of symbol xe at index e-1, one per event
         * @return numeric infinitesimal generator
         */
        public Matrix evalInfGen(double[] x) {
            if (x.length != symbols.size()) {
                throw new IllegalArgumentException(
                        "Expected " + symbols.size() + " symbol values, got " + x.length);
            }
            int n = stateSpace.getNumRows();
            Matrix Q = new Matrix(n, n);
            for (int e = 0; e < symbols.size(); e++) {
                if (symbols.get(e) == null) {
                    continue;
                }
                double c = invertSymbol ? 1.0 / x[e] : x[e];
                Q = Q.add(c, infGenTerms.get(e));
            }
            return Q;
        }

        /**
         * Evaluate the symbolic generator at the given symbol assignment.
         *
         * @param assignment map from symbol name (e.g. "x1") to value
         * @return numeric infinitesimal generator
         */
        public Matrix evalInfGen(Map<String, Double> assignment) {
            double[] x = new double[symbols.size()];
            for (int e = 0; e < symbols.size(); e++) {
                String sym = symbols.get(e);
                if (sym == null) {
                    x[e] = 1.0;
                } else {
                    Double val = assignment.get(sym);
                    if (val == null) {
                        throw new IllegalArgumentException("No value assigned to symbol " + sym);
                    }
                    x[e] = val;
                }
            }
            return evalInfGen(x);
        }

        /**
         * Symbolic expression of generator entry (i,j), e.g. "2*x1 - 3*x2".
         *
         * @param i row index
         * @param j column index
         * @return expression string; "0" if the entry is zero
         */
        public String getSymbolicEntry(int i, int j) {
            StringBuilder sb = new StringBuilder();
            for (int e = 0; e < symbols.size(); e++) {
                if (symbols.get(e) == null) {
                    continue;
                }
                double c = infGenTerms.get(e).get(i, j);
                if (c == 0) {
                    continue;
                }
                if (sb.length() == 0) {
                    if (c < 0) {
                        sb.append("-");
                    }
                } else {
                    sb.append(c < 0 ? " - " : " + ");
                }
                double absC = Math.abs(c);
                if (invertSymbol) {
                    sb.append(formatCoeff(absC)).append("/").append(symbols.get(e));
                } else {
                    if (absC != 1.0) {
                        sb.append(formatCoeff(absC)).append("*");
                    }
                    sb.append(symbols.get(e));
                }
            }
            if (sb.length() == 0) {
                return "0";
            }
            return sb.toString();
        }

        private static String formatCoeff(double c) {
            if (c == Math.rint(c) && !Double.isInfinite(c)) {
                return Long.toString((long) c);
            }
            return Double.toString(c);
        }

        /**
         * The whole symbolic generator as expression strings, row major.
         *
         * <p>This is the wire form the computer algebra backend consumes, see
         * {@link jline.api.sym.SymEngine}.</p>
         *
         * @return n by n array of expressions, "0" where the entry is zero
         */
        public String[][] toExpressionMatrix() {
            int n = stateSpace.getNumRows();
            String[][] Q = new String[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    Q[i][j] = getSymbolicEntry(i, j);
                }
            }
            return Q;
        }

        /**
         * Symbols that actually occur in the generator.
         *
         * <p>An event with no positive rate contributes no symbol and is stored
         * as a null in {@link #symbols}; those nulls are dropped here so the
         * list matches the symbols the expressions mention.</p>
         *
         * @return the non-null symbol names, in event order
         */
        public List<String> activeSymbols() {
            List<String> active = new ArrayList<String>();
            for (int e = 0; e < symbols.size(); e++) {
                if (symbols.get(e) != null) {
                    active.add(symbols.get(e));
                }
            }
            return active;
        }

        /** Print the symbolic generator, one row of expressions per state. */
        public void prettyPrint() {
            int n = stateSpace.getNumRows();
            for (int i = 0; i < n; i++) {
                StringBuilder row = new StringBuilder();
                for (int j = 0; j < n; j++) {
                    if (j > 0) {
                        row.append("\t");
                    }
                    row.append(getSymbolicEntry(i, j));
                }
                System.out.println(row.toString());
            }
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
        public Matrix T;
        /** First-entry matrix (-Q22)^-1 Q21: row i is where vanishing state i lands. */
        public Matrix entry;
        public Object denseLUSolver;

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
        /**
         * Stationary distribution the analyzer solved for, and the generator and
         * state spaces it is indexed by. On a reducible chain the analyzer restricts
         * the generator to the component supporting the stationary distribution, so
         * pi is shorter than StateSpace; InfGenWork/StateSpaceWork/StateSpaceAggrWork
         * are the matching rows and must be consumed as one triple. On an irreducible
         * chain they are InfGen/StateSpace/StateSpaceAggr themselves.
         */
        public Matrix pi, StateSpaceWork, StateSpaceAggrWork, InfGenWork;
        /**
         * Derived START/PREEMPT rates, (stations x classes): how often per unit
         * time a class-r service starts at station i, and how often a class-r
         * job in service is pushed back into the buffer there. Annotations on
         * the arcs the generator already carries, so they are not MetricType
         * entries and add no getAvgTable column. At a lossless station with no
         * in-service abandonment StartN == TN + PreemptN.
         */
        public Matrix StartN, PreemptN;
        /** The filtrations those rates reduce, indexed [station][class]. */
        public Matrix[][] startFilt, preemptFilt;

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
            this(QN, UN, RN, TN, CN, XN, InfGen, StateSpace, StateSpaceAggr,
                    EventFiltration, runtime, fname, sncopy, null, null, null, null);
        }

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
                NetworkStruct sncopy,
                Matrix pi,
                Matrix StateSpaceWork,
                Matrix StateSpaceAggrWork,
                Matrix InfGenWork) {
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
            this.pi = pi;
            this.StateSpaceWork = StateSpaceWork;
            this.StateSpaceAggrWork = StateSpaceAggrWork;
            this.InfGenWork = InfGenWork;
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

    public SampleResult sampleAggr(StatefulNode node, int numEvents) {
        assertPhaseTypeStates("sampleAggr");
        SampleResult result = sample(node, numEvents);
        if (result == null) {
            return null;
        }
        // AGGREGATE means COLLAPSE, not relabel. The reference is
        // sampleAggr.m:85-86, whose only difference from sample.m:86 is the
        // State.toMarginal call on the sliced state; it keeps the SECOND
        // output, nir, the per-class job counts. Setting isaggregate=true on a
        // phase-detailed matrix (as this used to) returns the wrong width and
        // the wrong content under the aggregate label, with nothing to signal
        // it. Note toMarginal, NOT toMarginalAggr: that is what the reference
        // calls, and it is also the variant carrying the phasessz single-row
        // clamp. K/Ks are passed null so it slices sn.phasessz itself.
        NetworkStruct sn = this.getStruct(this);
        int nodeIdx = node.getNodeIndex();
        int R = sn.nclasses;
        Matrix phaseDetailed = result.state;
        Matrix collapsed = new Matrix(phaseDetailed.getNumRows(), R);
        for (int row = 0; row < phaseDetailed.getNumRows(); row++) {
            Matrix stateRow = new Matrix(1, phaseDetailed.getNumCols());
            Matrix.extract(phaseDetailed, row, row + 1, 0, phaseDetailed.getNumCols(), stateRow, 0, 0);
            State.StateMarginalStatistics st =
                    ToMarginal.toMarginal(sn, nodeIdx, stateRow, null, null, null, null, null);
            for (int r = 0; r < R; r++) {
                collapsed.set(row, r, st.nir.get(0, r));
            }
        }
        result.state = collapsed;
        result.isaggregate = true;
        return result;
    }

    public jline.io.Ret.SampleResult sampleSysAggr(int numEvents) {
        // Collapse each stateful node's block as sampleSysAggr.m does; the JAR keeps the
        // joint matrix shape rather than MATLAB's per-node cells, see _kb/07-cross-language-parity.md
        jline.io.Ret.SampleResult result = sampleSys(numEvents);
        if (result == null || !(result.state instanceof Matrix)) {
            return result;
        }
        NetworkStruct sn = this.getStruct(this);
        Matrix joint = (Matrix) result.state;
        int R = sn.nclasses;
        List<StatefulNode> sfNodes = this.model.getStatefulNodes();
        int nsf = sfNodes.size();
        Matrix collapsed = new Matrix(joint.getNumRows(), nsf * R);
        int colOffset = 0;
        for (int isf = 0; isf < nsf; isf++) {
            StatefulNode nd = sfNodes.get(isf);
            int width = sn.space.get(nd).getNumCols();
            int nodeIdx = nd.getNodeIndex();
            for (int row = 0; row < joint.getNumRows(); row++) {
                Matrix stateRow = new Matrix(1, width);
                Matrix.extract(joint, row, row + 1, colOffset, colOffset + width, stateRow, 0, 0);
                State.StateMarginalStatistics st =
                        ToMarginal.toMarginal(sn, nodeIdx, stateRow, null, null, null, null, null);
                for (int r = 0; r < R; r++) {
                    collapsed.set(row, isf * R + r, st.nir.get(0, r));
                }
            }
            colOffset += width;
        }
        result.state = collapsed;
        result.isAggregate = true;
        return result;
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

            // see _kb/06-solver-catalog.md for rationale

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

    /**
     * Run the reward analyzer and cache results.
     * Convenience wrapper calling solver_ctmc_reward and storing results.
     *
     * @return RewardResult containing value functions, time vector, names, and steady-state rewards
     */
    public RewardResult runRewardAnalyzer() {
        NetworkStruct sn = this.model.getStruct(true);
        if (sn.reward == null || sn.reward.isEmpty()) {
            throw new IllegalStateException(
                "No rewards defined. Use model.setReward(name, rewardFn) before calling reward analysis.");
        }
        this.rewardResult = Solver_ctmc_reward.solver_ctmc_reward(sn, this.options);
        return this.rewardResult;
    }

    /**
     * Get reward value function and state space, with optional filtering by reward name.
     * Alias matching MATLAB getReward() signature.
     *
     * @param rewardName Optional reward name to filter. If null, returns all rewards.
     * @return RewardResult containing value functions, time vector, names, state space
     */
    public RewardResult getReward(String rewardName) {
        RewardResult result = getRewardResult();
        if (rewardName == null) {
            return result;
        }
        // Filter to specific reward
        Matrix V = result.getValueFunction().get(rewardName);
        if (V == null) {
            throw new IllegalArgumentException("Reward '" + rewardName + "' not found. Available rewards: " +
                String.join(", ", result.getRewardNames()));
        }
        Map<String, Matrix> filteredV = new HashMap<String, Matrix>();
        filteredV.put(rewardName, V);
        Map<String, Double> filteredSS = new HashMap<String, Double>();
        filteredSS.put(rewardName, result.getSteadyState().get(rewardName));
        List<String> filteredNames = new ArrayList<String>();
        filteredNames.add(rewardName);
        return new RewardResult(filteredV, result.getTime(), filteredNames, result.getStateSpace(), filteredSS, result.getRuntime());
    }

    /**
     * Get reward value function and state space for all rewards.
     *
     * @return RewardResult containing all rewards
     */
    public RewardResult getReward() {
        return getReward(null);
    }

    /**
     * Get transient expected reward E[r(X(t))] over time.
     *
     * Computes transient expected rewards using:
     *   E[r(X(t))] = sum_s pi_t(s) * r(s)
     *
     * where pi_t is the transient probability distribution at time t.
     *
     * Requires a finite timespan set via SolverCTMC(model, options.timespan([0,T])).
     *
     * @param rewardName Optional reward name to filter. If null, returns all rewards.
     * @return Map from reward name to double[] of expected reward values at each time point.
     *         Use getRewardTimeVector() or the result's time field to get the corresponding time points.
     */
    public Map<String, double[]> getTranReward(String rewardName) {
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException(
                "getTranReward requires a finite timespan. Use SolverCTMC(model, options.timespan([0,T])).");
        }

        NetworkStruct sn = this.model.getStruct(true);
        if (sn.reward == null || sn.reward.isEmpty()) {
            throw new IllegalStateException(
                "No rewards defined. Use model.setReward(name, rewardFn) before calling getTranReward.");
        }

        // Get transient probabilities
        TransientResult transientResult = this.solver_ctmc_transient_analyzer(sn, this.options);
        Matrix t = transientResult.t;
        Matrix pit = transientResult.pit;
        Matrix stateSpaceAggr = transientResult.StateSpaceAggr;

        int nTimePoints = t.getNumRows();
        int nstates = pit.getNumCols();

        // Clamp negative probabilities to zero
        for (int i = 0; i < nTimePoints; i++) {
            for (int j = 0; j < nstates; j++) {
                if (pit.get(i, j) < 0) {
                    pit.set(i, j, 0.0);
                }
            }
        }

        // Build reward vectors
        Map<String, double[]> rewardVectors = new HashMap<String, double[]>();
        List<String> names = new ArrayList<String>();
        for (Map.Entry<String, jline.lang.reward.RewardFunction> entry : sn.reward.entrySet()) {
            String name = entry.getKey();
            jline.lang.reward.RewardFunction rewardFn = entry.getValue();
            double[] rv = new double[nstates];
            for (int s = 0; s < nstates; s++) {
                Matrix stateRow = stateSpaceAggr.getRow(s);
                rv[s] = rewardFn.compute(stateRow, sn);
            }
            rewardVectors.put(name, rv);
            names.add(name);
        }

        // Compute E[r(X(t))] = pit * R' for each time point
        Map<String, double[]> result = new HashMap<String, double[]>();
        for (String name : names) {
            if (rewardName != null && !rewardName.equals(name)) {
                continue;
            }
            double[] rv = rewardVectors.get(name);
            double[] tranReward = new double[nTimePoints];
            for (int ti = 0; ti < nTimePoints; ti++) {
                double sum = 0.0;
                for (int s = 0; s < nstates; s++) {
                    sum += pit.get(ti, s) * rv[s];
                }
                tranReward[ti] = sum;
            }
            result.put(name, tranReward);
        }

        return result;
    }

    /**
     * Get transient expected reward for all rewards.
     *
     * @return Map from reward name to transient expected reward time series
     */
    public Map<String, double[]> getTranReward() {
        return getTranReward(null);
    }


    /**
     * Time-stationary per-item occupancy of each cache list.
     *
     * <p>The cache-contents block of the local-variable vector holds the item index
     * resident in each cache position, so P(item i is held by list l) is a state reward
     * of the stationary distribution. This is the TIME-WEIGHTED occupancy, the CTMC
     * counterpart of the EMBEDDED (per-request) occupancy the NC/MVA cache algorithms
     * return; the two coincide only when requests see time averages (PASTA). Mirrors the
     * per-item block in the MATLAB solver_ctmc_analyzer.</p>
     *
     * @param sn network structure carrying the per-node state spaces
     * @param ind cache node index
     * @param isf cache stateful index
     * @param SS global state space
     * @param pi stationary distribution over SS
     * @param cacheNode cache node receiving the result
     */
    private void computeCacheItemProb(NetworkStruct sn, int ind, int isf, Matrix SS, Matrix pi, Cache cacheNode) {
        CacheNodeParam np = (CacheNodeParam) sn.nodeparam.get(this.model.getNodes().get(ind));
        if (np == null || SS == null || pi == null || sn.space == null || np.itemcap == null) {
            return;
        }
        int nitems = np.nitems;
        int h = np.itemcap.length();
        if (nitems <= 0 || h <= 0) {
            return;
        }
        int tcc = 0;
        for (int l = 0; l < h; l++) {
            tcc += (int) np.itemcap.get(l);
        }
        int colOff = 0;
        for (int k = 0; k < isf; k++) {
            Matrix sk = sn.space.get(sn.stateful.get(k));
            colOff += (sk == null) ? 0 : sk.getNumCols();
        }
        Matrix scache = sn.space.get(sn.stateful.get(isf));
        if (scache == null) {
            return;
        }
        // local row layout: [per-class server presence | contents | block A | block B]
        int tail = 0;
        if (np.retrievalSystemCapacity > 0) {
            tail = nitems + State.cacheRetrievalClassMap(sn, ind)[1].length;
        }
        int lvs = scache.getNumCols() - (tcc + tail);
        if (lvs < 0 || colOff + scache.getNumCols() > SS.getNumCols() || SS.getNumRows() != pi.length()) {
            return;
        }
        Matrix itemprob = new Matrix(nitems, h + 1);
        boolean[] present = new boolean[nitems];
        for (int row = 0; row < SS.getNumRows(); row++) {
            double w = pi.get(row);
            if (w == 0) {
                continue;
            }
            int off = colOff + lvs;
            for (int l = 0; l < h; l++) {
                int cap = (int) np.itemcap.get(l);
                java.util.Arrays.fill(present, false);
                for (int q = 0; q < cap; q++) {
                    int it = (int) SS.get(row, off + q);
                    if (it >= 1 && it <= nitems) {
                        present[it - 1] = true;
                    }
                }
                for (int i = 0; i < nitems; i++) {
                    if (present[i]) {
                        itemprob.set(i, l + 1, itemprob.get(i, l + 1) + w);
                    }
                }
                off += cap;
            }
        }
        for (int i = 0; i < nitems; i++) {
            double miss = 1.0;
            for (int l = 0; l < h; l++) {
                miss -= itemprob.get(i, l + 1);
            }
            itemprob.set(i, 0, miss);
        }
        cacheNode.setResultItemProb(itemprob);
    }

    /**
     * Exact delayed-hit queue length of a retrieval-system cache.
     *
     * <p>Block A of the cache local-variable vector marks the items being fetched and
     * block B counts, per retrieval class, the secondary requests merged onto those
     * fetches, so phi_i = P(a fetch of item i is in flight), d1_i = E[secondary requests
     * waiting on the fetch of item i] and dfull_i = d1_i + phi_i are state rewards of the
     * stationary distribution, hence exact. Mirrors the delayed-hit block in the MATLAB
     * solver_ctmc_analyzer.</p>
     *
     * @param sn network structure carrying the per-node state spaces
     * @param ind cache node index
     * @param isf cache stateful index
     * @param SS global state space
     * @param pi stationary distribution over SS
     * @param cacheNode cache node receiving the result
     */
    private void computeDelayedHitQLen(NetworkStruct sn, int ind, int isf, Matrix SS, Matrix pi, Cache cacheNode) {
        CacheNodeParam np = (CacheNodeParam) sn.nodeparam.get(this.model.getNodes().get(ind));
        if (np.retrievalSystemCapacity <= 0 || SS == null || pi == null || sn.space == null) {
            return;
        }
        int[][] rcMap = State.cacheRetrievalClassMap(sn, ind);
        int[] rcItems = rcMap[1];
        int nitems = np.nitems;
        int tcc = np.totalCacheCapacity;
        int colOff = 0;
        for (int k = 0; k < isf; k++) {
            Matrix sk = sn.space.get(sn.stateful.get(k));
            colOff += (sk == null) ? 0 : sk.getNumCols();
        }
        Matrix scache = sn.space.get(sn.stateful.get(isf));
        if (scache == null) {
            return;
        }
        int lvs = scache.getNumCols() - (tcc + nitems + rcItems.length);
        if (lvs < 0 || colOff + scache.getNumCols() > SS.getNumCols() || SS.getNumRows() != pi.length()) {
            return;
        }
        int a0 = colOff + lvs + tcc;
        int b0 = a0 + nitems;
        Matrix d1 = new Matrix(1, nitems);
        Matrix dfull = new Matrix(1, nitems);
        for (int i = 0; i < nitems; i++) {
            double phi = 0, d = 0;
            for (int row = 0; row < SS.getNumRows(); row++) {
                double w = pi.get(row);
                if (SS.get(row, a0 + i) != 0) {
                    phi += w;
                }
                for (int j = 0; j < rcItems.length; j++) {
                    if (rcItems[j] == i + 1) {
                        d += w * SS.get(row, b0 + j);
                    }
                }
            }
            d1.set(0, i, d);
            dfull.set(0, i, d + phi);
        }
        cacheNode.setResultDelayedHitQLen(d1, dfull);
    }

    /**
     * Exact delayed-hit rate per originating class, as a transition reward.
     *
     * <p>A fetch of item i completes on exactly the transitions that clear block A bit i,
     * and each such transition releases the block-B counts of item i as delayed hits. The
     * rate is therefore a reward over the GENERATOR, not over the states: the alternative
     * arrival-rate identity lambda_i*phi_i is only PASTA-exact.</p>
     *
     * @param sn network structure carrying the per-node state spaces
     * @param ind cache node index
     * @param isf cache stateful index
     * @param SS global state space
     * @param pi stationary distribution over SS
     * @param Q infinitesimal generator over SS
     * @return delayed-hit rate per originating job class
     */
    private Matrix delayedHitRate(NetworkStruct sn, int ind, int isf, Matrix SS, Matrix pi, Matrix Q) {
        // column nclasses (one past the classes) accumulates the total fetch-completion
        // rate, which equals the miss rate because a completed fetch yields exactly one
        // miss; it supplies the normalizer for the probability split.
        Matrix rate = new Matrix(1, sn.nclasses + 1);
        rate.zero();
        CacheNodeParam np = (CacheNodeParam) sn.nodeparam.get(this.model.getNodes().get(ind));
        if (np.retrievalSystemCapacity <= 0 || SS == null || pi == null || Q == null || sn.space == null) {
            return rate;
        }
        int[][] rcMap = State.cacheRetrievalClassMap(sn, ind);
        int[] rcItems = rcMap[1], rcOrig = rcMap[2];
        int nitems = np.nitems, tcc = np.totalCacheCapacity;
        int colOff = 0;
        for (int k = 0; k < isf; k++) {
            Matrix sk = sn.space.get(sn.stateful.get(k));
            colOff += (sk == null) ? 0 : sk.getNumCols();
        }
        Matrix scache = sn.space.get(sn.stateful.get(isf));
        if (scache == null) {
            return rate;
        }
        int lvs = scache.getNumCols() - (tcc + nitems + rcItems.length);
        if (lvs < 0 || SS.getNumRows() != pi.length()) {
            return rate;
        }
        int a0 = colOff + lvs + tcc;
        int b0 = a0 + nitems;
        for (int i = 0; i < nitems; i++) {
            for (int row = 0; row < SS.getNumRows(); row++) {
                if (SS.get(row, a0 + i) == 0) {
                    continue;
                }
                double cr = 0;
                for (int col = 0; col < Q.getNumCols(); col++) {
                    if (col == row) continue;
                    double q = Q.get(row, col);
                    if (q != 0 && SS.get(col, a0 + i) == 0) {
                        cr += q;
                    }
                }
                rate.set(0, sn.nclasses, rate.get(0, sn.nclasses) + pi.get(row) * cr);
            }
        }
        for (int j = 0; j < rcItems.length; j++) {
            int i = rcItems[j] - 1;
            for (int row = 0; row < SS.getNumRows(); row++) {
                if (SS.get(row, a0 + i) == 0 || SS.get(row, b0 + j) <= 0) {
                    continue;
                }
                double completionRate = 0;
                for (int col = 0; col < Q.getNumCols(); col++) {
                    if (col == row) continue;
                    double q = Q.get(row, col);
                    if (q != 0 && SS.get(col, a0 + i) == 0) {
                        completionRate += q;
                    }
                }
                if (completionRate != 0) {
                    rate.set(0, rcOrig[j], rate.get(0, rcOrig[j])
                            + pi.get(row) * SS.get(row, b0 + j) * completionRate);
                }
            }
        }
        return rate;
    }

    /**
     * Splits the cache hit-class rate into true hits and delayed hits.
     *
     * <p>Delayed hits depart in the hit class, so the hit-class rate is (true hits +
     * delayed hits). Applying the exact delayed rate makes hit + delayed + miss = 1,
     * matching the LDES/NC report.</p>
     *
     * @param sn network structure
     * @param ind cache node index
     * @param isf cache stateful index
     * @param SS global state space
     * @param pi stationary distribution
     * @param Q infinitesimal generator
     * @param cacheNode cache node receiving the result
     * @param cnp cache node parameters holding the hit/miss probabilities
     */
    private void applyDelayedHitSplit(NetworkStruct sn, int ind, int isf, Matrix SS, Matrix pi,
                                      Matrix Q, Cache cacheNode, CacheNodeParam cnp) {
        if (cnp.retrievalSystemCapacity <= 0 || cnp.actualhitprob == null || cnp.actualmissprob == null) {
            return;
        }
        Matrix rate = delayedHitRate(sn, ind, isf, SS, pi, Q);
        Matrix hp = cnp.actualhitprob.copy();
        Matrix mp = cnp.actualmissprob.copy();
        Matrix dp = new Matrix(hp.getNumRows(), hp.getNumCols());
        dp.zero();
        for (int k = 0; k < Math.min(cnp.hitclass.getNumCols(), hp.getNumCols()); k++) {
            int h = (int) cnp.hitclass.get(k);
            int m = (int) cnp.missclass.get(k);
            if (h < 0 || m < 0) continue;
            // hp/mp are already normalized probabilities; recover the rate normalizer
            // from missRate = denom * missProb, then express the delayed rate in the
            // same units.
            double hitProb = hp.get(k), missProb = mp.get(k);
            double missRate = rate.get(0, sn.nclasses);
            if (!(missProb > 0) || !(missRate > 0)) continue;
            double denom = missRate / missProb;
            double d = Math.min(rate.get(0, k) / denom, hitProb);
            hp.set(k, hitProb - d);
            dp.set(k, d);
        }
        cnp.actualhitprob = hp;
        cnp.actualmissprob = mp;
        cnp.actualdelayedhitprob = dp;
        cacheNode.setResultHitProb(hp);
        cacheNode.setResultMissProb(mp);
        cacheNode.setResultDelayedHitProb(dp);
    }
}
