package jline.solvers.auto;

import static jline.io.InputOutput.line_debug;
import static jline.io.InputOutput.line_warning;

import jline.GlobalConstants;
import jline.lang.Network;
import jline.lang.FeatureSet;
import jline.lang.NetworkStruct;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.ln.SolverLN;
import jline.VerboseLevel;
import jline.lang.nodes.Node;
import jline.solvers.NetworkAvgChainTable;
import jline.solvers.NetworkAvgNodeTable;
import jline.solvers.NetworkAvgSysTable;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.ag.SolverAG;
import jline.solvers.ba.SolverBA;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.env.SolverENV;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.JobClass;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Automatic solver selection for queueing network models.
 * <p>
 * This solver automatically selects the most appropriate solution method
 * based on model characteristics and requested performance metrics.
 */
public class SolverAUTO extends NetworkSolver {

    // JMT is not an AUTO candidate: LDES subsumes its feature set, so automatic
    // selection never dispatches to the external simulator.
    public static final int CANDIDATE_CTMC = 5;
    public static final int CANDIDATE_LDES = 6;
    public static final int CANDIDATE_FLUID = 3;
    public static final int CANDIDATE_MAM = 2;
    // Solver candidates constants
    public static final int CANDIDATE_MVA = 0;
    public static final int CANDIDATE_NC = 1;
    public static final int CANDIDATE_SSA = 4;
    // public static final String METHOD_AI = "ai";  // AI method not yet available
    // Method selection strategies
    public static final String METHOD_DEFAULT = "default";
    public static final String METHOD_HEURISTIC = "heur";
    public static final String METHOD_SIM = "sim";
    public static final String METHOD_EXACT = "exact";
    public static final String METHOD_FAST = "fast";
    public static final String METHOD_ACCURATE = "accurate";
    public static final String METHOD_BOUND = "bound";

    /**
     * Population at or below which an exact solver is preferred over an
     * approximation, when one is available.
     */
    public static final int EXACT_POPULATION_MAX = 5;

    // Solver instances
    private List<NetworkSolver> candidates;
    // Bound-analysis delegate for METHOD_BOUND, kept out of the candidate pool
    // because it returns bounds rather than point estimates.
    private NetworkSolver boundSolver;
    private Map<String, Integer> solverNameToId;
    private NetworkSolver selectedSolver;
    private final String selectionMethod;

    /**
     * The delegate a FAMILY method name pinned, or null when the method name was an intent.
     *
     * A family method name ("nc", "nc.comom", or a bare algorithm name a family
     * declares) asks for a named engine and bypasses the ranking entirely, which
     * is what the reference's selectionMode = tokenFamily does.
     */
    private NetworkSolver pinnedSolver;

    /** The submethod the family method name pinned, "default" for a bare family name. */
    private String pinnedMethod;

    // Options specific to AUTO solver
    private final AUTOptions autoOptions;

    /**
     * Constructor with model only
     */
    public SolverAUTO(Network model) {
        this(model, METHOD_DEFAULT);
    }

    /**
     * Constructor with model and method
     */
    public SolverAUTO(Network model, String method) {
        this(model, new AUTOptions(method));
    }

    /**
     * Constructor with model and options
     */
    public SolverAUTO(Network model, SolverOptions options) {
        super(model, "SolverAuto", options);

        if (options instanceof AUTOptions) {
            this.autoOptions = (AUTOptions) options;
        } else {
            this.autoOptions = new AUTOptions(options);
        }

        this.selectionMethod = autoOptions.selectionMethod;
        initializeCandidates();
        // A method name is either a selection intent or a request for a
        // specific algorithm. Resolving it HERE keeps the qualified form
        // 'family.submethod' intact and turns an unknown method name into an error
        // rather than into a silent fall-through to the heuristic, which is what
        // the switch in selectSolver used to do. Mirrors the MATLAB constructor.
        AutoToken tok = resolveMethodToken(this.selectionMethod);
        if ("family".equals(tok.kind)) {
            SolverOptions famOptions = this.options.copy();
            famOptions.method = tok.submethod;
            Object built = buildFamilySolver(tok.family, famOptions);
            this.pinnedSolver = (NetworkSolver) built;
            this.pinnedMethod = tok.submethod;
        }
    }

    /**
     * Constructor with model and varargs
     */
    public SolverAUTO(Network model, Object... varargin) {
        this(model, parseAUTOptions(varargin));
    }

    /**
     * Parse AUTOptions from varargs
     */
    private static AUTOptions parseAUTOptions(Object... varargin) {
        AUTOptions options = new AUTOptions();

        for (int i = 0; i < varargin.length - 1; i += 2) {
            if (varargin[i] instanceof String) {
                String key = (String) varargin[i];
                Object value = varargin[i + 1];

                switch (key.toLowerCase()) {
                    case "method":
                        options.selectionMethod = (String) value;
                        break;
                    case "force":
                        options.forceSolver = (String) value;
                        break;
                    case "verbose":
                        if (value instanceof VerboseLevel) {
                            options.verbose = (VerboseLevel) value;
                        } else if (value instanceof Integer) {
                            int level = (Integer) value;
                            if (level == 0) {
                                options.verbose = VerboseLevel.SILENT;
                            } else if (level == 1) {
                                options.verbose = VerboseLevel.STD;
                            } else {
                                options.verbose = VerboseLevel.DEBUG;
                            }
                        }
                        break;
                    default:
                        // Pass other options to base SolverOptions
                        break;
                }
            }
        }

        return options;
    }

    /**
     * Copy results from the selected solver to this solver
     */
    private void copyResultsFromSolver(NetworkSolver solver) {
        // Copy the result object
        if (solver.result != null) {
            this.result = solver.result;
        }

        // The base NetworkSolver class will handle result access through getAvg() etc.
        // No need to copy protected fields directly
    }

    /**
     * Ensure a solver is selected before delegating calls
     */
    private void ensureSolverSelected() {
        if (selectedSolver == null) {
            selectSolver();
        }
    }

    /**
     * Find a solver by its class type from the candidates
     */
    private NetworkSolver findSolverByType(Class<? extends NetworkSolver> solverClass) {
        for (NetworkSolver solver : candidates) {
            if (solverClass.isInstance(solver)) {
                return solver;
            }
        }
        return null;
    }

    /**
     * Get list of candidate solver names
     */
    public List<String> getCandidateSolverNames() {
        List<String> names = new ArrayList<String>();
        for (NetworkSolver solver : candidates) {
            names.add(solver.getName());
        }
        return names;
    }

    /**
     * The CTMC candidate, built on demand. State-space and generator accessors
     * are CTMC-only concepts, so they resolve here rather than through the
     * ranked selection.
     *
     * @return a SolverCTMC over the same model
     */
    private SolverCTMC ctmcSolver() {
        for (NetworkSolver solver : candidates) {
            if (solver instanceof SolverCTMC) {
                return (SolverCTMC) solver;
            }
        }
        return new SolverCTMC(model);
    }

    /**
     * Infinitesimal generator, always from SolverCTMC.
     *
     * @return the generator, event filtrations and synchronization info
     */
    public SolverCTMC.generatorResult getGenerator() {
        return ctmcSolver().getGenerator();
    }

    /**
     * Infinitesimal generator, always from SolverCTMC.
     *
     * @param options solver options for the CTMC solve
     * @return the generator, event filtrations and synchronization info
     */
    public SolverCTMC.generatorResult getGenerator(SolverOptions options) {
        return ctmcSolver().getGenerator(options);
    }

    /**
     * Symbolic infinitesimal generator, always from SolverCTMC.
     *
     * @return the symbolic generator decomposition
     */
    public SolverCTMC.symbolicGeneratorResult getSymbolicGenerator() {
        return ctmcSolver().getSymbolicGenerator();
    }

    /**
     * Symbolic infinitesimal generator, always from SolverCTMC.
     *
     * @param invertSymbol divide each event filtration by its symbol instead of
     *                     multiplying
     * @return the symbolic generator decomposition
     */
    public SolverCTMC.symbolicGeneratorResult getSymbolicGenerator(boolean invertSymbol) {
        return ctmcSolver().getSymbolicGenerator(invertSymbol);
    }

    /**
     * State space, always from SolverCTMC.
     *
     * @return the global and per-node state spaces
     */
    public SolverCTMC.StateSpace getStateSpace() {
        return ctmcSolver().getStateSpace();
    }

    /**
     * State space, always from SolverCTMC.
     *
     * @param options solver options for the state-space generation
     * @return the global and per-node state spaces
     */
    public SolverCTMC.StateSpace getStateSpace(SolverOptions options) {
        return ctmcSolver().getStateSpace(options);
    }

    /**
     * Model structure, as every other solver exposes it.
     *
     * @return the NetworkStruct of the model under analysis
     */
    public NetworkStruct getStruct() {
        if (this.sn == null) {
            this.sn = this.model.getStruct(false);
        }
        return this.sn;
    }

    /**
     * Every token the CONSTRUCTOR accepts, INDEPENDENT of the model: the
     * selection intents, the method families, and every family method in its
     * QUALIFIED form. The unqualified form is accepted too (resolveMethodToken
     * looks it up) and is left out here to keep the list unambiguous.
     *
     * <p>THIS is the list a method-NAME check must gate on.
     * {@link #listValidMethods()} narrows it to the model in hand, and gating a
     * name check on that would replace a rejection the delegate would have
     * EXPLAINED with a flat "the method is unsupported by this solver" -- the
     * same distinction SolverBA draws between its own listAllMethods and
     * listValidMethods.
     *
     * @return every AUTO token, regardless of the model
     */
    public String[] listAllMethods() {
        Set<String> methods = new TreeSet<String>(Arrays.asList(selectionIntents()));
        methods.add("auto");
        SolverOptions probeOptions = this.options.copy();
        probeOptions.verbose = VerboseLevel.SILENT;
        probeOptions.method = METHOD_DEFAULT;
        for (String fam : familyNames()) {
            List<String> declared;
            try {
                Object probe = buildFamilySolver(fam, probeOptions);
                declared = NetworkSolver.declaredAllMethods(probe);
                if (declared == null) {
                    declared = NetworkSolver.declaredValidMethods(probe);
                }
            } catch (Exception e) {
                // A family that cannot even be INSTANTIATED here ('ln', 'env',
                // 'lqns' and 'uq' do not solve a Network) contributes nothing.
                continue;
            }
            if (declared == null) {
                continue;
            }
            methods.add(fam);
            for (String d : declared) {
                methods.add(fam + "." + d);
            }
        }
        return methods.toArray(new String[0]);
    }

    /**
     * The method names of {@link #listAllMethods()} that THIS MODEL can actually run.
     *
     * <p>The gate is the one {@code chooseSolverRanked} already applies before
     * delegating, asked here of every candidate instead of the first feasible
     * one: a method whose {@code supportsModelMethod} refuses the model is not
     * offered, and a family with no method-level gate is judged by its flat
     * feature set through {@code supports(model)}. The method-level gate is
     * where the rules a feature set cannot express live (product form for
     * "exact", a binding finite buffer, NC "mem" applicability).
     *
     * <p>Without it this returned the token universe regardless of the model,
     * naming every SolverNC method on the two-station BAS-blocking model of
     * cqn_bas_blocking although SolverNC refuses that model method by method. A
     * caller enumerating the list was being invited to ask for an analysis no
     * candidate would perform.
     *
     * <p>A family whose every method is refused loses its bare method name too: "nc"
     * alone delegates to SolverNC, which is exactly the rejection the per-method
     * gate just returned.
     *
     * <p>This used to return the union of the intents and every candidate's
     * UNQUALIFIED method names -- a set the class then ignored, since any method name
     * but an intent fell through to the heuristic. It also cast each candidate's
     * answer to String[], which throws for the half of the solvers that return
     * List&lt;String&gt;, so those contributed nothing and the failure was
     * swallowed.
     *
     * @return the AUTO tokens this model can run
     */
    public String[] listValidMethods() {
        // IT IS THE RUNNABLE ROWS OF findSolver, projected onto their method
        // token. The narrowing used to be written out a second time here, and a
        // second copy of one gate is how two answers to one question start to
        // differ; findSolver owns it now, and this adds only the method names that
        // name no single method: the selection intents, which name a RANKING
        // rather than an algorithm, and each family's bare method name.
        //
        // A family whose every method is refused loses its bare method name too:
        // "nc" alone delegates to SolverNC, which is exactly the rejection the
        // per-method gate just returned. That falls out of the projection,
        // since a family with no runnable row contributes no row to name.
        Set<String> methods = new TreeSet<String>(Arrays.asList(selectionIntents()));
        methods.add("auto");
        for (SolverCandidate row : findSolver()) {
            methods.add(row.method);
            methods.add(row.solver);
        }
        return methods.toArray(new String[0]);
    }

    // ------------------------------------------------------------------
    // findSolver: which solvers and methods can analyze this model
    // Port of @@SolverAUTO/findSolver.m and its native python twin.
    // ------------------------------------------------------------------

    /**
     * The measure groups {@link #findSolver} reports on, in report order.
     *
     * <p>A group is a family of accessors that stand or fall together: a solver
     * that returns getCdfRespT returns getCdfPassT and getPerctRespT as well,
     * because all three read the same passage time, so listing the three
     * separately would say nothing extra.
     *
     * @return the group names
     */
    public static String[] metricGroups() {
        return new String[]{"avg", "tran", "cdf", "prob", "tranprob", "sample",
                "cache", "loss", "orbit", "moment", "sens"};
    }

    /**
     * The measure group an accessor belongs to, "" when the name is none.
     *
     * <p>A group name maps to itself, so findSolver("cdf") and
     * findSolver("getCdfRespT") ask the same question.
     *
     * <p>THIS IS NOT chooseSolverHeur's TABLE, although both are keyed by
     * accessor name. That one maps an accessor to a RANKING, i.e. which
     * candidate should be preferred; this one maps it to a CAPABILITY question,
     * i.e. which candidates can answer it at all. The two differ wherever a
     * family can serve a measure but is never the one AUTO would pick for it.
     *
     * @param name an accessor name or a group name
     * @return the group, or "" when the name belongs to none
     */
    public static String metricGroupOf(String name) {
        if (name == null || name.isEmpty()) {
            return "";
        }
        for (String g : metricGroups()) {
            if (g.equals(name)) {
                return g;
            }
        }
        if ("any".equalsIgnoreCase(name) || "all".equalsIgnoreCase(name)) {
            return "";
        }
        if (in(name, "getTranAvg", "getTranAvgVar", "tranAvg")) return "tran";
        if (in(name, "getCdfRespT", "getCdfPassT", "getPerctRespT", "getTranCdfPassT",
                "getTranCdfRespT", "getCdfSysRespT")) return "cdf";
        if (in(name, "getTranProb", "getTranProbSys", "getTranProbAggr",
                "getTranProbSysAggr")) return "tranprob";
        if (in(name, "getProb", "getProbAggr", "getProbSys", "getProbSysAggr",
                "getProbMarg", "getProbNormConstAggr")) return "prob";
        if (in(name, "sample", "sampleSys", "sampleAggr", "sampleSysAggr")) return "sample";
        if (in(name, "getAvgCacheTable", "getAvgCacheT", "getAvgItemTable", "getAvgItemT",
                "cacheAvgT", "itemAvgT", "aCaT", "aIT")) return "cache";
        if (in(name, "getAvgLossTable", "getAvgLossT", "getAvgRegionLossTable",
                "getAvgRegionLossT", "lossAvgT", "regionLossAvgT", "aLT", "aRLT")) return "loss";
        if (in(name, "getAvgOrbitTable", "getAvgOrbitT", "getAvgOrbit", "orbitAvgT",
                "aOT")) return "orbit";
        if (in(name, "getMomentTable", "getMomentChainTable", "getMomentStationTable",
                "getMomentT", "getMomentChainT", "getMomentStationT", "momentT",
                "momentChainT", "momentStationT", "mT", "mCT", "mST")) return "moment";
        if (in(name, "getSensitivityTable", "getSensitivityT", "sensitivityT", "sT",
                "getSensitivity", "getSensitivityRanking")) return "sens";
        // Everything else in the accessor surface is a mean measure: getAvg, its
        // chain, node and system forms, their handles and their short aliases.
        if (name.startsWith("getAvg") || name.startsWith("avg")
                || in(name, "getAvgSysRespT", "getAvgSysTput", "aT", "aNT", "aCT", "aST", "aNCT")) {
            return "avg";
        }
        return "";
    }

    private static boolean in(String name, String... names) {
        for (String n : names) {
            if (n.equals(name)) {
                return true;
            }
        }
        return false;
    }

    /**
     * The measure groups a method family can answer.
     *
     * <p>Every family answers "avg", which is what a solver is for; the rest is
     * the capability declaration this class owns.
     *
     * <p>SOURCES, so that a claim here can be checked rather than trusted:
     * "tran" is supportsTransientAnalysis, which FLD, CTMC, LDES and JMT
     * override to true and no one else does. "cdf", "prob", "tranprob" and
     * "sample" are the families that carry an implementation of the
     * corresponding accessor rather than inheriting the base refusal. The
     * remaining five groups are computed by NetworkSolver from a solver's own
     * results, so no per-solver method marks them: their lists are
     * chooseSolverHeur's rankings for the same accessors, which is where AUTO
     * already records who can serve them.
     *
     * <p>A family that gains or loses a measure must be edited here in the same
     * change, the way a solver that gains a feature is edited into its feature
     * set: an omission here does not fail, it silently hides the family from a
     * caller asking for that measure.
     *
     * @param family the method family
     * @return its measure groups
     */
    public static String[] familyMetrics(String family) {
        if ("mva".equals(family)) {
            return new String[]{"avg", "prob", "cache", "orbit", "moment", "sens"};
        }
        if ("nc".equals(family)) {
            return new String[]{"avg", "cdf", "prob", "cache", "moment", "sens"};
        }
        if ("ctmc".equals(family)) {
            return new String[]{"avg", "tran", "cdf", "prob", "tranprob", "sample",
                    "cache", "loss", "orbit", "moment"};
        }
        if ("fluid".equals(family)) {
            return new String[]{"avg", "tran", "cdf", "prob", "cache", "sens"};
        }
        if ("mam".equals(family)) {
            return new String[]{"avg", "cdf"};
        }
        if ("ag".equals(family)) {
            // The RCAT fixed point reports means only; the passage-time law it
            // answers is the base exponential fit, not its own.
            return new String[]{"avg"};
        }
        if ("ba".equals(family)) {
            // A bound brackets the mean measures and nothing else.
            return new String[]{"avg"};
        }
        if ("ssa".equals(family)) {
            return new String[]{"avg", "cdf", "prob", "sample", "loss"};
        }
        if ("ldes".equals(family)) {
            return new String[]{"avg", "tran", "cdf", "prob", "sample", "cache", "loss", "orbit"};
        }
        if ("jmt".equals(family)) {
            return new String[]{"avg", "tran", "cdf", "prob", "tranprob", "sample"};
        }
        if ("ln".equals(family)) {
            return new String[]{"avg", "tran", "cdf", "sens"};
        }
        if ("env".equals(family)) {
            return new String[]{"avg", "tran"};
        }
        return new String[]{"avg"};
    }

    /**
     * What KIND of answer a method returns: exact, approx, bound or simulation.
     *
     * <p>"simulation" is not decided here: isStochastic is the solver's own
     * isStochasticMethod, which already tokenizes qualified and
     * runtime-resolved names and is the only place that knowledge lives.
     *
     * <p>"exact" IS CLAIMED ONLY WHERE IT IS TRUE OF THIS MODEL, never of the
     * algorithm in the abstract. Exactness of a normalizing constant or of mean
     * value analysis is a property of the product-form model it is computed on,
     * and of the QBD shape for the matrix analytic methods, so both conditions
     * are passed in and a method that needs one reports "approx" without it.
     * The bias is deliberate: an under-claimed "approx" costs a user a better
     * method they could have had, an over-claimed "exact" costs them a wrong
     * number they trusted.
     *
     * <p>A CACHE IS THE THIRD CONDITION, and it was the over-claim the bias
     * above exists to prevent. snHasProductForm answers about the QUEUEING
     * network and knows nothing of a cache: the hit/miss split is a class switch
     * whose probabilities are not routing data but the output of a cache model,
     * so a network holding one reads as product form and "mva.exact" was
     * labelled exact on it. Measured on the tut06 shape with an LRU cache: exact
     * MVA returns QLen 0.2516 at the hit station where the CTMC returns 0.3022
     * and simulation 0.3023, a 17% error under a label that says there is none.
     * The analytic families are conditioned on it; SolverCTMC is NOT, because
     * its state space carries the cache contents and it is exact there, which is
     * what the two numbers above show.
     *
     * @param family        the method family
     * @param method        the unqualified method name
     * @param isStochastic  the solver's own isStochasticMethod verdict
     * @param isProductForm whether the model has a product-form solution
     * @param isQbdShape    whether the model is one queueing station fed by a Source
     * @param hasCache      whether the model holds a Cache node
     * @return one of the SolverCandidate CLASS_ constants
     */
    public static String methodClass(String family, String method, boolean isStochastic,
                                     boolean isProductForm, boolean isQbdShape,
                                     boolean hasCache) {
        if ("ba".equals(family)) {
            // Bounds are what SolverBA is for; every one of its methods returns
            // a bracket rather than an estimate.
            return SolverCandidate.CLASS_BOUND;
        }
        if (isStochastic) {
            return SolverCandidate.CLASS_SIMULATION;
        }
        if ("ctmc".equals(family)) {
            // The generator is solved as written, so every state-space route is
            // exact. "cftp.approx" says in its own name that it is not, and
            // "mdd" is exact on a product-form model and an approximation
            // otherwise.
            if ("cftp.approx".equals(method)) {
                return SolverCandidate.CLASS_APPROX;
            }
            if ("mdd".equals(method)) {
                return exactIf(isProductForm);
            }
            return SolverCandidate.CLASS_EXACT;
        }
        if ("nc".equals(family)) {
            // The normalizing-constant routes that evaluate G exactly rather
            // than expanding or estimating it. The asymptotic expansions
            // (pana, le, kt, bk, gm, ...) and the non-product-form
            // "morrison" are approximations by construction and are left out.
            if (in(method, "exact", "divdiff", "ca", "comom", "comomld", "rec", "ms",
                    "cub", "rgf")) {
                return exactIf(isProductForm && !hasCache);
            }
            return SolverCandidate.CLASS_APPROX;
        }
        if ("mva".equals(family)) {
            // Exact MVA; every "amva.*" arm is an approximation, and so are the
            // open-network QNA transforms.
            if (in(method, "exact", "mva")) {
                return exactIf(isProductForm && !hasCache);
            }
            return SolverCandidate.CLASS_APPROX;
        }
        if ("jmt".equals(family)) {
            // JMVA's exact algorithms. "jsim" and "replication" are simulation
            // and never reach here.
            if (in(method, "jmva.mva", "jmva.recal", "jmva.comom", "jmva.treeconv")) {
                return exactIf(isProductForm && !hasCache);
            }
            return SolverCandidate.CLASS_APPROX;
        }
        if ("mam".equals(family)) {
            // The QBD is solved exactly on the shape it is stated for, one
            // queueing station fed by a Source. Everything named "dec.*" is a
            // decomposition of a larger network into such queues and is
            // therefore an approximation of it.
            if (in(method, "default", "mna", "ldqbd", "bgchain", "retrial")) {
                return exactIf(isQbdShape);
            }
            return SolverCandidate.CLASS_APPROX;
        }
        if ("ag".equals(family)) {
            // Every RCAT arm estimates the reversed rate of each synchronising
            // action and iterates to a fixed point, an approximation by
            // construction; SolverAG's "exact" is a vestigial alias that warns
            // and runs "inap", so nothing here is claimed exact.
            return SolverCandidate.CLASS_APPROX;
        }
        return SolverCandidate.CLASS_APPROX;
    }

    private static String exactIf(boolean cond) {
        return cond ? SolverCandidate.CLASS_EXACT : SolverCandidate.CLASS_APPROX;
    }

    /**
     * Which solvers and solver methods can analyze this model: the runnable
     * pairs, for every measure.
     *
     * @return one row per runnable (family, method) pair
     */
    public List<SolverCandidate> findSolver() {
        return findSolver("", false);
    }

    /**
     * Which solvers and solver methods can analyze this model, and for the ones
     * that cannot, why not.
     *
     * <p>One row per (family, method) pair AUTO can be asked for; see
     * {@link SolverCandidate} for the columns.
     *
     * <p>THE GATE IS NOT A SECOND ONE. It is the gate chooseSolverRanked applies
     * before delegating, asked of every candidate instead of of the first
     * feasible one, which is exactly what {@link #listValidMethods} already did
     * -- that method is now the method column of the runnable rows, so the two
     * cannot disagree. What is new is that the REASON the gate produced is kept
     * rather than discarded, and that the answer carries the two facts a caller
     * needs in order to choose among the survivors: whether the method is exact
     * on this model, and which measures it can report.
     *
     * <p>WHY A REASON HAS TO BE RECONSTRUCTED for some rows. The base
     * supportsModelMethod returns a reason that names no feature when the solver
     * does not diverge per method, since it then falls back to supports(model),
     * which answers with a bare boolean. That is enough for a gate, which only
     * has to stop the run, and not enough for a report, whose whole content is
     * the explanation. So a refused row is re-asked against the solver's own
     * feature set, which is where the offending feature names are.
     *
     * @param metric  a measure group ("cdf") or the accessor that returns it
     *                ("getCdfRespT"); "" or "any" keeps every pair
     * @param showAll keep the refused pairs too, with the reason each was refused
     * @return the rows, in family order
     */
    public List<SolverCandidate> findSolver(String metric, boolean showAll) {
        String group = metricGroupOf(metric);
        if (group.isEmpty() && metric != null && !metric.isEmpty()
                && !"any".equalsIgnoreCase(metric) && !"all".equalsIgnoreCase(metric)) {
            throw new RuntimeException("'" + metric + "' names no measure. Pass a group ("
                    + String.join(", ", metricGroups())
                    + ") or the accessor that returns it, e.g. 'getCdfRespT'.");
        }

        // The three model properties an exactness claim can rest on, evaluated
        // once: a method whose exactness needs one of them reports "approx"
        // without it, see methodClass.
        boolean isProductForm = false;
        boolean isQbdShape = false;
        boolean hasCache = false;
        try {
            isProductForm = model.hasProductFormSolution();
            NetworkStruct sn = model.getStruct(false);
            int nsources = 0;
            for (NodeType t : sn.nodetype) {
                if (t == NodeType.Source) {
                    nsources++;
                }
            }
            boolean allOpen = sn.njobs.length() > 0;
            for (int r = 0; r < sn.njobs.length(); r++) {
                if (!Double.isInfinite(sn.njobs.get(r))) {
                    allOpen = false;
                }
            }
            isQbdShape = allOpen && (sn.nstations - nsources) == 1;
            for (NodeType t : sn.nodetype) {
                if (t == NodeType.Cache) {
                    hasCache = true;
                }
            }
        } catch (Exception e) {
            // A model whose struct cannot be refreshed here answers "approx"
            // everywhere, which is the safe direction: see methodClass.
        }

        SolverOptions probeOptions = this.options.copy();
        probeOptions.verbose = VerboseLevel.SILENT;
        probeOptions.method = METHOD_DEFAULT;

        // A REPORT MUST NOT PRINT. Asking a solver whether it supports the model
        // runs FeatureSet.supports, which warns naming the missing feature: a
        // side effect that is right on the solve path, where nobody asked to be
        // told, and wrong here, where every refused row would raise one and the
        // answer IS the table.
        VerboseLevel savedVerbose = GlobalConstants.getVerbose();
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
        List<SolverCandidate> rows = new ArrayList<SolverCandidate>();
        try {
            for (String family : familyNames()) {
                String[] groups = familyMetrics(family);
                if (!group.isEmpty() && !contains(groups, group)) {
                    continue;
                }
                Object probe;
                List<String> declared;
                try {
                    probe = buildFamilySolver(family, probeOptions);
                    declared = NetworkSolver.declaredValidMethods(probe);
                } catch (Exception e) {
                    // A family that cannot even be instantiated here (ln, env,
                    // lqns and uq take another model type) contributes nothing:
                    // there is no solver to report on and no gate to ask.
                    continue;
                }
                if (declared == null) {
                    continue;
                }
                String metricList = String.join(",", groups);
                // Once per family, not once per refused row: the flat feature
                // set does not vary with the method, and a family that refuses
                // every one of its forty methods would otherwise recompute the
                // same comparison forty times.
                String famFeatureReason = featureReason(probe);
                for (String name : declared) {
                    if (name.startsWith(family + ".")) {
                        // A spelling already qualified with its own family.
                        // SolverFluid declares both "dae" and "fluid.dae" so
                        // that its own gate takes either, and prefixing the
                        // family again yields "fluid.fluid.dae": a method name that
                        // does resolve, but that names the same method twice and
                        // would double every fluid row of this report.
                        continue;
                    }
                    if (isMethodAlias(family, name, declared)) {
                        // The same duplication under a different prefix.
                        // SolverMVA advertises every AMVA name twice, plain and
                        // "amva."-prefixed, and its dispatch strips the prefix,
                        // so the two spellings are one algorithm; that alone was
                        // 20 of the 49 mva rows of a report. The plain spelling
                        // is the one kept.
                        continue;
                    }
                    String reason = gateReason(probe, name);
                    boolean ok = reason.isEmpty();
                    if (!ok && isGenericReason(reason)) {
                        // A SPECIFIC GATE REASON WINS; the feature-set names
                        // only replace a generic one ("Some features are not
                        // supported by the chosen solver", which names none) or
                        // fill an empty one. Getting this the other way round
                        // buried real reasons under a list of every feature the
                        // model uses, for any solver whose feature set answers a
                        // different question than "what do I accept".
                        if (!famFeatureReason.isEmpty()) {
                            reason = famFeatureReason;
                        } else if (reason.isEmpty()) {
                            reason = "This solver refuses the model through its own "
                                    + "structural check.";
                        }
                    }
                    if (!ok && !showAll) {
                        continue;
                    }
                    rows.add(new SolverCandidate(family, family + "." + name, ok,
                            methodClass(family, name, stochasticVerdict(probe, name),
                                    isProductForm, isQbdShape, hasCache),
                            metricList, ok ? "" : reason));
                }
            }
        } finally {
            GlobalConstants.setVerbose(savedVerbose);
        }
        return rows;
    }

    /**
     * findSolver for a LayeredNetwork, which this class does not otherwise
     * solve: its constructor takes a Network, and buildFamilySolver says so.
     *
     * <p>THE FAMILY SET IS THE ONE THAT CAN BE BUILT OVER A LayeredNetwork, and
     * not a shorter hand-picked list. MATLAB's familyAcceptsModelClass admits
     * every family for a layered model and lets the ones that cannot be
     * constructed drop out, which is how "ldes" stays in -- it analyzes an LQN
     * natively -- while mva, nc, ctmc and the rest fall away because their
     * constructors take a Network. Naming only "ln" and "lqns" here would have
     * hidden the LDES row that MATLAB and native python both report.
     *
     * @param model   the layered model
     * @param metric  a measure group or the accessor that returns it; "" keeps all
     * @param showAll keep the refused pairs too
     * @return the rows, in family order
     */
    public static List<SolverCandidate> findSolverLayered(LayeredNetwork model, String metric,
                                                          boolean showAll) {
        String group = metricGroupOf(metric);
        if (group.isEmpty() && metric != null && !metric.isEmpty()
                && !"any".equalsIgnoreCase(metric) && !"all".equalsIgnoreCase(metric)) {
            throw new RuntimeException("'" + metric + "' names no measure. Pass a group ("
                    + String.join(", ", metricGroups())
                    + ") or the accessor that returns it, e.g. 'getCdfRespT'.");
        }
        SolverOptions probeOptions = new SolverOptions();
        probeOptions.verbose = VerboseLevel.SILENT;
        probeOptions.method = METHOD_DEFAULT;

        VerboseLevel savedVerbose = GlobalConstants.getVerbose();
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
        List<SolverCandidate> rows = new ArrayList<SolverCandidate>();
        try {
            String[] families = new String[]{"ldes", "ln", "lqns"};
            for (String family : families) {
                String[] groups = familyMetrics(family);
                if (!group.isEmpty() && !contains(groups, group)) {
                    continue;
                }
                Object probe;
                List<String> declared;
                try {
                    if ("ln".equals(family)) {
                        probe = new SolverLN(model, probeOptions);
                    } else if ("lqns".equals(family)) {
                        probe = new jline.solvers.wrappers.lqns.SolverLQNS(model, probeOptions);
                    } else {
                        probe = new jline.solvers.ldes.SolverLDES(model, probeOptions);
                    }
                    declared = NetworkSolver.declaredValidMethods(probe);
                } catch (Exception e) {
                    // SolverLQNS without the lqns binary is the case this
                    // catches: no solver to report on and no gate to ask.
                    continue;
                }
                if (declared == null) {
                    continue;
                }
                String metricList = String.join(",", groups);
                for (String name : declared) {
                    if (name.startsWith(family + ".")) {
                        continue;
                    }
                    String reason = gateReason(probe, name);
                    boolean ok = reason.isEmpty();
                    if (!ok && !showAll) {
                        continue;
                    }
                    // No product-form or QBD claim is made for a layered model:
                    // no layered solver is exact, which is what chooseSolverExact
                    // already records by ranking the NC layers first and calling
                    // them the closest available.
                    rows.add(new SolverCandidate(family, family + "." + name, ok,
                            methodClass(family, name, stochasticVerdict(probe, name), false, false,
                                    false),
                            metricList, ok ? "" : reason));
                }
            }
        } finally {
            GlobalConstants.setVerbose(savedVerbose);
        }
        return rows;
    }

    /** Array membership; the varargs {@code in} above cannot also take a String[]. */
    private static boolean contains(String[] hay, String needle) {
        for (String h : hay) {
            if (h.equals(needle)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Does the reason say only that SOME feature is unsupported, without naming
     * one? That is the sentence the shared base gate returns, and the one worth
     * replacing with the offending feature names. Empty counts: a gate that
     * answered with a bare boolean said nothing at all.
     */
    private static boolean isGenericReason(String reason) {
        if (reason == null || reason.isEmpty()) {
            return true;
        }
        return reason.toLowerCase().contains("features are not supported")
                && !reason.contains("(feature:");
    }

    /**
     * The method-level gate, asked without letting it raise: "" when the pair
     * runs, the reason otherwise. The base supportsModelMethod already falls
     * back to the flat supports(model) for a solver that does not diverge per
     * method, so this one call is both gates.
     */
    private static String gateReason(Object probe, String method) {
        if (!(probe instanceof Solver)) {
            return "";
        }
        try {
            String reason = ((Solver) probe).supportsModelMethod(method);
            return reason == null ? "" : reason;
        } catch (Exception e) {
            // A gate that raises has said something, and it is the only thing it
            // can say about this pair; reporting it beats swallowing it and
            // calling the pair runnable.
            String msg = e.getMessage();
            return (msg == null || msg.isEmpty()) ? e.toString() : msg;
        }
    }

    /**
     * isStochasticMethod asked without letting it raise; a solver that cannot
     * classify a name is taken at its class default, deterministic.
     */
    private static boolean stochasticVerdict(Object probe, String method) {
        if (!(probe instanceof Solver)) {
            return false;
        }
        try {
            return ((Solver) probe).isStochasticMethod(method);
        } catch (Exception e) {
            return false;
        }
    }

    /**
     * The offending feature names, or "" when the solver's feature set accepts
     * the model. The empty answer is meaningful and not a failure: it says the
     * refusal came from somewhere the feature set cannot see, so the caller
     * should keep whatever the gate itself said.
     *
     * <p>getFeatureSet is static on every solver and is not declared by any
     * shared supertype, so it is reached by reflection rather than by a cast.
     */
    private String featureReason(Object probe) {
        try {
            java.lang.reflect.Method m = probe.getClass().getMethod("getFeatureSet");
            Object fs = m.invoke(null);
            if (!(fs instanceof FeatureSet)) {
                return "";
            }
            // A SOLVER THAT DECLARES NOTHING HAS NO ENVELOPE, and a
            // missing-feature list against an empty set is not an explanation:
            // it is every feature the model uses. Such a refusal is structural,
            // so the caller keeps the gate's own words instead.
            FeatureSet declared = (FeatureSet) fs;
            boolean any = false;
            for (String feature : declared.featureNames()) {
                if (declared.inspectFeature(feature)) {
                    any = true;
                    break;
                }
            }
            if (!any) {
                return "";
            }
            return FeatureSet.supportsReason(declared, model.getUsedLangFeatures());
        } catch (Exception e) {
            return "";
        }
    }

    /**
     * Union of the feature sets of every solver AUTO can delegate to: a model
     * AUTO can analyze is one that at least one candidate supports.
     *
     * @return the union feature set over all candidate solvers
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        FeatureSet[] sets = new FeatureSet[]{
                SolverMVA.getFeatureSet(), SolverNC.getFeatureSet(),
                SolverMAM.getFeatureSet(), SolverFluid.getFeatureSet(),
                SolverSSA.getFeatureSet(),
                SolverCTMC.getFeatureSet(), SolverLDES.getFeatureSet()};
        for (FeatureSet fs : sets) {
            for (String feature : fs.featureNames()) {
                if (fs.inspectFeature(feature)) {
                    featSupported.setTrue(feature);
                }
            }
        }
        return featSupported;
    }

    /**
     * Get the name of the selected solver
     */
    public String getSelectedSolverName() {
        return selectedSolver != null ? selectedSolver.getName() : "none";
    }

    /**
     * If a delegate solver has been selected, classify by it (once run, it
     * knows the method it resolved at runtime). Before selection, the
     * delegate choice is unknown, so classify conservatively: true if any
     * candidate is stochastic.
     *
     * @return true if the (prospective) delegate returns stochastic estimates
     */
    @Override
    public boolean isStochastic() {
        if (selectedSolver != null) {
            return selectedSolver.isStochastic();
        }
        if (candidates != null) {
            for (NetworkSolver solver : candidates) {
                if (solver != null && solver.isStochastic()) {
                    return true;
                }
            }
        }
        return false;
    }

    // ------------------------------------------------------------------
    // Method-method name resolution.
    //
    // Port of @@SolverAUTO/resolveMethodToken.m, familyAlias, familyNames,
    // selectionIntents, familyDeclaringMethod and buildFamilySolver.
    //
    // WHAT WAS MISSING AND WHY IT MATTERED. This class used to switch on
    // selectionMethod alone, with a `default:` arm that fell through to the
    // heuristic, so every token that was not one of the six intents was
    // SILENTLY IGNORED: SolverAUTO(model, "nc.comom") ran the ranked heuristic
    // and answered with whatever it picked, under the caller's name. MATLAB,
    // native python and the C++ port all resolve the token first and error on
    // one no family owns. listValidMethods made it worse by advertising the
    // union of every candidate's unqualified method names, none of which were
    // honoured.
    // ------------------------------------------------------------------

    /** Tokens that state what the solution is for, rather than naming an algorithm. */
    public static String[] selectionIntents() {
        return new String[]{METHOD_DEFAULT, METHOD_HEURISTIC, METHOD_SIM, METHOD_EXACT,
                METHOD_FAST, METHOD_ACCURATE, METHOD_BOUND};
    }

    /**
     * Method families, in the order in which an unqualified method name is looked up.
     *
     * <p>"ag" SITS AFTER "mam", whose RCAT names it took over, and is a family
     * for the METHOD NAME and REPORT tables only: LINE(model, "ag.inap") and
     * model.help() reach SolverAG through it. It is deliberately NOT a candidate
     * of the automatic ranking (initializeCandidates builds that list by hand)
     * and NOT in the feature-set union of getFeatureSet, so a G-network is still
     * refused by "default" and has to be asked for by name; see
     * _kb/06-solver-catalog.md, "SolverAG owns the RCAT methods".
     */
    public static String[] familyNames() {
        return new String[]{"mva", "nc", "ctmc", "fluid", "mam", "ag", "ba", "ssa", "ldes", "jmt",
                "qns", "ln", "env", "lqns", "uq"};
    }

    /**
     * The prefixes under which a family advertises a SECOND SPELLING of a
     * method it already declares plainly.
     *
     * <p>THIS IS A DECLARATION, not a derivation, and belongs beside
     * familyMetrics and methodClass for the same reason: the knowledge lives in
     * the solver's own dispatch (SolverMVA strips a leading "amva." before
     * selecting an algorithm) and no accessor exposes it, so a family that gains
     * or loses an alias spelling must be edited into all four copies in the SAME
     * change. An omission does not fail; it puts the same algorithm in the
     * report twice.
     *
     * @param family the method family
     * @return the alias prefixes, empty when the family advertises none
     */
    public static String[] methodAliasPrefixes(String family) {
        if ("mva".equals(family)) {
            return new String[]{"amva."};
        }
        return new String[0];
    }

    /**
     * Is {@code name} a second spelling of another method this family declares?
     *
     * <p>The remainder has to be declared too, which is what keeps the rule from
     * eating a genuine method that merely starts with the prefix: it is an alias
     * only when the thing it aliases is there beside it.
     *
     * @param family   the method family
     * @param name     the declared method name under test
     * @param declared every method name the family declares on this model
     * @return true when the name only respells one of the others
     */
    public static boolean isMethodAlias(String family, String name, List<String> declared) {
        for (String prefix : methodAliasPrefixes(family)) {
            if (name.length() > prefix.length() && name.startsWith(prefix)
                    && declared.contains(name.substring(prefix.length()))) {
                return true;
            }
        }
        return false;
    }

    /**
     * Canonical family of a method name, or "" when the method name names none.
     *
     * @param name the token to classify
     * @return the family name, or "" when the method name names no family
     */
    public static String familyAlias(String name) {
        if (name == null || name.isEmpty()) {
            return "";
        }
        String n = name.toLowerCase();
        if (n.equals("mam")) return "mam";
        if (n.equals("ag")) return "ag";
        if (n.equals("mva")) return "mva";
        if (n.equals("nc")) return "nc";
        if (n.equals("fluid") || n.equals("fld")) return "fluid";
        if (n.equals("jmt")) return "jmt";
        if (n.equals("ssa")) return "ssa";
        if (n.equals("ctmc")) return "ctmc";
        if (n.equals("ldes") || n.equals("des")) return "ldes";
        if (n.equals("ba")) return "ba";
        if (n.equals("env")) return "env";
        if (n.equals("ln")) return "ln";
        if (n.equals("lqns") || n.equals("lqsim")) return "lqns";
        if (n.equals("qns")) return "qns";
        if (n.equals("uq")) return "uq";
        return "";
    }

    /** A method name split into a selection intent, or a family plus its submethod. */
    public static final class AutoToken {
        /** "intent" or "family". */
        public final String kind;
        /** The intent, when kind is "intent"; the family name otherwise. */
        public final String family;
        /** The method handed to the family; "default" for a bare family name. */
        public final String submethod;

        AutoToken(String kind, String family, String submethod) {
            this.kind = kind;
            this.family = family;
            this.submethod = submethod;
        }
    }

    /**
     * Find the family that declares an unqualified algorithm name, by asking each
     * family for its own method list. Keeping the question there avoids a second
     * copy of the name table here, which would drift.
     *
     * @param token the unqualified algorithm name
     * @return the family that declares it, or "" when none does
     */
    public String familyDeclaringMethod(String token) {
        SolverOptions probeOptions = new SolverOptions();
        probeOptions.verbose = VerboseLevel.SILENT;
        probeOptions.method = METHOD_DEFAULT;
        for (String fam : familyNames()) {
            try {
                Object probe = buildFamilySolver(fam, probeOptions);
                List<String> declared = NetworkSolver.declaredValidMethods(probe);
                if (declared == null) {
                    continue;
                }
                for (String d : declared) {
                    if (d.equalsIgnoreCase(token)) {
                        return fam;
                    }
                }
            } catch (Exception e) {
                // A family that cannot even be instantiated on this model cannot
                // own the token; the next one is asked instead.
            }
        }
        return "";
    }

    /**
     * Split a method name into a selection intent, or into a family and the
     * submethod handed to it. A qualified method name "family.submethod" KEEPS its
     * submethod: dropping it would silently downgrade a pinned method to the
     * family default.
     *
     * @param token the requested method name, may be null
     * @return the resolved token
     * @throws RuntimeException when no family owns the method name
     */
    public AutoToken resolveMethodToken(String token) {
        String tok = (token == null || token.isEmpty()) ? METHOD_DEFAULT : token;
        if ("auto".equalsIgnoreCase(tok)) {
            tok = METHOD_DEFAULT;
        }
        for (String intent : selectionIntents()) {
            if (intent.equalsIgnoreCase(tok)) {
                return new AutoToken("intent", intent.toLowerCase(), METHOD_DEFAULT);
            }
        }
        int dot = tok.indexOf('.');
        String head = dot < 0 ? tok : tok.substring(0, dot);
        String rest = dot < 0 ? "" : tok.substring(dot + 1);
        String fam = familyAlias(head);
        if (!fam.isEmpty()) {
            // A bare family name means its default method; for bounds the
            // default is the composite tightest-of-all family, not a single one.
            String sub = rest.isEmpty() ? ("ba".equals(fam) ? "auto" : METHOD_DEFAULT) : rest;
            return new AutoToken("family", fam, sub);
        }
        // Unqualified algorithm name, e.g. "comom" or "gb.upper". The family
        // that declares it owns it, so the name table stays in the families.
        fam = familyDeclaringMethod(tok);
        if (!fam.isEmpty()) {
            return new AutoToken("family", fam, tok);
        }
        StringBuilder fams = new StringBuilder();
        for (String f : familyNames()) {
            if (fams.length() > 0) fams.append(", ");
            fams.append(f);
        }
        throw new RuntimeException("Unrecognized method '" + tok
                + "'. Valid tokens are a selection intent (default, heur, sim, exact, fast, "
                + "accurate, bound), a method family (" + fams
                + "), or a qualified method name such as 'nc.comom'.");
    }

    /**
     * Instantiate the solver of a method family. {@code options.method} already
     * holds the submethod resolved by {@link #resolveMethodToken}, so a pinned
     * method reaches the family that runs it.
     *
     * @param family  the family name
     * @param options the options, carrying the resolved submethod
     * @return the family's solver
     */
    public Object buildFamilySolver(String family, SolverOptions options) {
        if ("mam".equals(family)) return new SolverMAM(model, options);
        if ("ag".equals(family)) return new SolverAG(model, options);
        if ("mva".equals(family)) return new SolverMVA(model, options);
        if ("nc".equals(family)) return new SolverNC(model, options);
        if ("fluid".equals(family)) return new SolverFluid(model, options);
        if ("jmt".equals(family)) return new jline.solvers.wrappers.jmt.SolverJMT(model, options);
        if ("ssa".equals(family)) return new SolverSSA(model, options);
        if ("ctmc".equals(family)) return new SolverCTMC(model, options);
        if ("ldes".equals(family)) return new SolverLDES(model, options);
        if ("ba".equals(family)) return new SolverBA(model, options);
        if ("qns".equals(family)) return new jline.solvers.wrappers.qns.SolverQNS(model, options);
        // 'ln', 'env', 'lqns' and 'uq' take a LayeredNetwork, an Environment or
        // an inner-solver factory rather than this Network, so they are not
        // constructible from here; the reference builds them in its own
        // buildFamilySolver, which this class reaches only for Network models.
        throw new RuntimeException("The '" + family + "' method family does not solve a Network "
                + "model from SolverAUTO; construct Solver" + family.toUpperCase()
                + " directly with its own model type.");
    }

    /**
     * Bound-analysis delegate: SolverBA with its own family selection ('auto').
     *
     * @return the lazily created SolverBA instance used by METHOD_BOUND
     */
    private NetworkSolver getBoundSolver() {
        if (boundSolver == null) {
            SolverOptions baOptions = SolverBA.defaultOptions();
            baOptions.method = "auto";
            baOptions.verbose = options.verbose;
            boundSolver = new SolverBA(model, baOptions);
        }
        return boundSolver;
    }

    /**
     * Initialize solver candidates based on model type
     */
    private void initializeCandidates() {
        candidates = new ArrayList<NetworkSolver>();
        solverNameToId = new HashMap<>();

        // Add all solver candidates
        candidates.add(new SolverMVA(model));
        solverNameToId.put("mva", CANDIDATE_MVA);

        candidates.add(new SolverNC(model));
        solverNameToId.put("nc", CANDIDATE_NC);

        candidates.add(new SolverMAM(model));
        solverNameToId.put("mam", CANDIDATE_MAM);

        candidates.add(new SolverFluid(model));
        solverNameToId.put("fluid", CANDIDATE_FLUID);

        candidates.add(new SolverSSA(model));
        solverNameToId.put("ssa", CANDIDATE_SSA);

        // CTMC solver integration
        candidates.add(new SolverCTMC(model));
        solverNameToId.put("ctmc", CANDIDATE_CTMC);

        // LDES solver integration
        candidates.add(new SolverLDES(model));
        solverNameToId.put("ldes", CANDIDATE_LDES);

        // Filter candidates by those that support the model
        List<NetworkSolver> supportedCandidates = new ArrayList<NetworkSolver>();
        for (NetworkSolver solver : candidates) {
            if (solver.supports(model)) {
                supportedCandidates.add(solver);
            }
        }
        candidates = supportedCandidates;

        if (candidates.isEmpty()) {
            throw new RuntimeException("No solver supports this model");
        }
    }

    /**
     * Override the main run method to ensure proper delegation
     */
    @Override
    public void runAnalyzer() {
        // Propagate solver verbose level to global
        if (this.options != null) {
            GlobalConstants.Verbose = options.verbose;
        }
        line_debug(options.verbose, String.format("AUTO solver starting: selectionMethod=%s", selectionMethod));
        // A FAMILY METHOD NAME BYPASSES THE RANKING. The caller named an engine, so
        // running the ranked heuristic instead would answer a different
        // algorithm under their name -- and there is no fallback either, for the
        // same reason an explicit 'bound' request must not degrade.
        if (pinnedSolver != null) {
            line_debug(options.verbose, String.format("AUTO: method token pinned solver=%s, method=%s",
                    pinnedSolver.getName(), pinnedMethod));
            selectedSolver = pinnedSolver;
            pinnedSolver.getAvg();
            copyResultsFromSolver(pinnedSolver);
            if (options.verbose != VerboseLevel.SILENT) {
                System.out.println("SolverAuto: Successfully used " + pinnedSolver.getName());
            }
            return;
        }
        ensureSolverSelected();
        line_debug(options.verbose, String.format("AUTO: selected solver=%s", selectedSolver.getName()));

        // Try the selected solver first by using getAvg() which will trigger runAnalyzer()
        try {
            selectedSolver.getAvg(); // This will call runAnalyzer() internally
            copyResultsFromSolver(selectedSolver);

            if (options.verbose != VerboseLevel.SILENT) {
                System.out.println("SolverAuto: Successfully used " + selectedSolver.getName());
            }
            return;
        } catch (Exception e) {
            if (options.verbose != VerboseLevel.SILENT) {
                System.out.println("SolverAuto: " + selectedSolver.getName() +
                        " failed: " + e.getMessage());
            }
            // An explicit 'bound' request must not silently degrade to a point
            // estimate from the candidate pool.
            if (METHOD_BOUND.equals(selectionMethod)) {
                throw new RuntimeException("SolverAuto: bound analysis failed: " + e.getMessage(), e);
            }
        }

        // Try other candidates if the selected solver fails
        for (NetworkSolver candidate : candidates) {
            if (candidate != selectedSolver) {
                try {
                    line_debug(options.verbose, String.format("AUTO: trying fallback solver=%s", candidate.getName()));
                    if (options.verbose != VerboseLevel.SILENT) {
                        System.out.println("SolverAuto: Trying " + candidate.getName());
                    }

                    candidate.getAvg(); // This will call runAnalyzer() internally
                    copyResultsFromSolver(candidate);
                    selectedSolver = candidate; // Update to successful solver

                    if (options.verbose != VerboseLevel.SILENT) {
                        System.out.println("SolverAuto: Successfully used " + candidate.getName());
                    }
                    return;
                } catch (Exception e) {
                    if (options.verbose == VerboseLevel.DEBUG) {
                        System.out.println("SolverAuto: " + candidate.getName() +
                                " also failed: " + e.getMessage());
                    }
                }
            }
        }

        throw new RuntimeException("SolverAuto: No solver could handle this model");
    }

    /**
     * Force selection of a specific solver
     */
    private void selectForcedSolver(String solverName) {
        Integer solverId = solverNameToId.get(solverName.toLowerCase());
        if (solverId != null && solverId < candidates.size()) {
            selectedSolver = candidates.get(solverId);
        } else {
            throw new RuntimeException("Unknown solver: " + solverName);
        }
    }

    /**
     * Select the most appropriate solver based on model characteristics
     */
    private void selectSolver() {
        if (autoOptions.forceSolver != null && !autoOptions.forceSolver.isEmpty()) {
            // Force specific solver
            line_debug(options.verbose, String.format("AUTO: forcing solver=%s", autoOptions.forceSolver));
            selectForcedSolver(autoOptions.forceSolver);
        } else {
            // Automatic selection based on method
            line_debug(options.verbose, String.format("AUTO: selecting solver using method=%s", selectionMethod));
            switch (selectionMethod) {
                case METHOD_DEFAULT:
                case METHOD_HEURISTIC:
                    selectSolverHeuristic();
                    break;
                // case METHOD_AI:  // AI method not yet available
                //     selectSolverAI();
                //     break;
                case METHOD_SIM:
                    selectedSolver = chooseSolverSim("getAvg");
                    if (selectedSolver == null) {
                        selectSolverHeuristic();
                    }
                    break;
                case METHOD_EXACT:
                    selectedSolver = chooseSolverExact("getAvg");
                    if (selectedSolver == null) {
                        throw new RuntimeException("SolverAuto: no exact solver supports this model; "
                                + "use method 'default' for the approximate heuristic");
                    }
                    break;
                case METHOD_FAST:
                    selectedSolver = chooseSolverRanked(new Class<?>[]{SolverMVA.class,
                            SolverNC.class, SolverFluid.class, SolverMAM.class});
                    if (selectedSolver == null) {
                        selectSolverHeuristic();
                    }
                    break;
                case METHOD_BOUND:
                    selectedSolver = getBoundSolver();
                    break;
                case METHOD_ACCURATE:
                    selectedSolver = chooseSolverRanked(new Class<?>[]{SolverFluid.class,
                            SolverMAM.class, SolverCTMC.class, SolverLDES.class});
                    if (selectedSolver == null) {
                        selectSolverHeuristic();
                    }
                    break;
                default:
                    selectSolverHeuristic();
            }
        }

        if (options.verbose != VerboseLevel.SILENT) {
            System.out.println("SolverAuto: Selected " + selectedSolver.getName());
        }
    }

    // /**
    //  * Select solver using AI-based method (not yet available)
    //  */
    // private void selectSolverAI() {
    //     // Placeholder - would require loading trained model
    //     // For now, fall back to heuristic
    //     selectSolverHeuristic();
    // }

    /**
     * First solver in the ranked order whose feature set accepts the model.
     *
     * @param order candidate solver classes, most preferred first
     * @return the first feasible solver, or null when none qualifies
     */
    private NetworkSolver chooseSolverRanked(Class<?>[] order) {
        return chooseSolverRanked(order, null);
    }

    /**
     * Ranked choice gated on a concrete method. With a non-null methodToken the
     * gate tightens from supports(model) to supportsModelMethod(methodToken),
     * which is where the rules a flat feature set cannot express live (product
     * form for 'exact', finite capacity, NC 'mem' applicability).
     *
     * @param order       ranked candidate classes
     * @param methodToken the method to gate on, or null for the coarse gate
     * @return the first feasible candidate, or null when none qualifies
     */
    private NetworkSolver chooseSolverRanked(Class<?>[] order, String methodToken) {
        for (int k = 0; k < order.length; k++) {
            for (NetworkSolver solver : candidates) {
                boolean supported = methodToken == null
                        ? solver.supports(model)
                        : solver.supportsModelMethod(methodToken).isEmpty();
                if (order[k].isInstance(solver) && supported) {
                    // Feature support is necessary but not sufficient for CTMC:
                    // the chain must also fit memory.
                    if (solver instanceof SolverCTMC
                            && !SolverCTMC.isStateSpaceTractable(model, solver.getOptions()).ok) {
                        continue;
                    }
                    // Pin the delegate's method: the gate above only ASKED
                    // whether the candidate can run methodToken, and delegate
                    // hands the solver its own options, so without this an
                    // exactness-gated choice would still run the solver's
                    // default (approximate) method.
                    solver.getOptions().method =
                            methodToken == null ? METHOD_DEFAULT : methodToken;
                    return solver;
                }
            }
        }
        return null;
    }

    /**
     * Cache node present: cache metrics invert the analytical order to NC first.
     *
     * @return true when the model contains a Cache node
     */
    private boolean hasCacheNode() {
        NetworkStruct s = getStruct();
        if (s.nodetype == null) {
            return false;
        }
        for (NodeType t : s.nodetype) {
            if (t == NodeType.Cache) {
                return true;
            }
        }
        return false;
    }

    /**
     * Autocorrelated arrival or service process: only MAM keeps the correlation.
     *
     * @return true when any process is a MAP or MMPP2
     */
    private boolean hasMAPProcess() {
        NetworkStruct s = getStruct();
        if (s.procid == null) {
            return false;
        }
        for (Map<JobClass, ProcessType> byClass : s.procid.values()) {
            for (ProcessType t : byClass.values()) {
                if (t == ProcessType.MAP || t == ProcessType.MMPP2) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Priority flavour, which ranks differently: preemptive priority is served
     * analytically only by MAM, and the PS family only by CTMC or a simulator.
     *
     * @return "preempt", "ps", "hol" or "none"
     */
    private String priorityKind() {
        NetworkStruct s = getStruct();
        if (s.sched == null) {
            return "none";
        }
        boolean hol = false;
        for (SchedStrategy sched : s.sched.values()) {
            if (sched == SchedStrategy.FCFSPRPRIO || sched == SchedStrategy.FCFSPIPRIO
                    || sched == SchedStrategy.LCFSPRPRIO || sched == SchedStrategy.LCFSPIPRIO) {
                return "preempt";
            }
            if (sched == SchedStrategy.PSPRIO || sched == SchedStrategy.DPSPRIO
                    || sched == SchedStrategy.GPSPRIO) {
                return "ps";
            }
            if (sched == SchedStrategy.HOL || sched == SchedStrategy.LCFSPRIO) {
                hol = true;
            }
        }
        return hol ? "hol" : "none";
    }

    /**
     * Ranked choice for mean-value metrics. Global order is MVA &gt; NC &gt; MAM,
     * inverted to NC &gt; MVA on cache models, with Fluid promoted for large
     * populations and MAM promoted on autocorrelated traffic. Mirrors MATLAB's
     * chooseAvgSolverHeur.m.
     *
     * @return the ranked solver for average metrics
     */
    private NetworkSolver chooseAvgSolverHeuristic() {
        ModelAnalyzer analyzer = new ModelAnalyzer(model);
        Class<?>[] order;
        String prio = priorityKind();
        // Small populations: an approximation buys nothing there, so take an
        // exact solver whenever one is available, preferring MVA over NC over
        // CTMC (inverted to NC first on caches). The "exact" token is what
        // makes this a claim rather than a preference: MVA and NC reject it
        // without a product-form solution and CTMC rejects it when the chain
        // does not fit memory.
        if (analyzer.getTotalJobs() > 0 && analyzer.getTotalJobs() <= EXACT_POPULATION_MAX) {
            Class<?>[] exactOrder = hasCacheNode()
                    ? new Class<?>[]{SolverNC.class, SolverMVA.class, SolverCTMC.class}
                    : new Class<?>[]{SolverMVA.class, SolverNC.class, SolverCTMC.class};
            NetworkSolver exact = chooseSolverRanked(exactOrder, METHOD_EXACT);
            if (exact != null) {
                return exact;
            }
        }
        if (hasCacheNode()) {
            order = new Class<?>[]{SolverNC.class, SolverMVA.class, SolverFluid.class,
                    SolverCTMC.class, SolverLDES.class};
        } else if (getStruct().nregions > 0) {
            if (analyzer.getTotalJobs() <= 10) {
                order = new Class<?>[]{SolverNC.class, SolverCTMC.class, SolverLDES.class};
            } else {
                order = new Class<?>[]{SolverNC.class, SolverLDES.class, SolverCTMC.class};
            }
        } else if ("preempt".equals(prio)) {
            order = new Class<?>[]{SolverMAM.class, SolverLDES.class, SolverCTMC.class, SolverSSA.class};
        } else if ("ps".equals(prio)) {
            order = new Class<?>[]{SolverCTMC.class, SolverLDES.class, SolverSSA.class};
        } else if (hasMAPProcess()) {
            order = new Class<?>[]{SolverMAM.class, SolverMVA.class, SolverFluid.class,
                    SolverLDES.class};
        } else if ("hol".equals(prio)) {
            order = new Class<?>[]{SolverMVA.class, SolverMAM.class, SolverFluid.class,
                    SolverCTMC.class, SolverLDES.class};
        } else if (analyzer.getAvgJobsPerChain() > 30) {
            order = new Class<?>[]{SolverFluid.class, SolverMVA.class, SolverNC.class};
        } else if (analyzer.getTotalJobs() > 0 && analyzer.getTotalJobs() <= EXACT_POPULATION_MAX) {
            // No exact solver was available at this population (tried above),
            // so keep the exact-leaning approximate order.
            order = new Class<?>[]{SolverNC.class, SolverMVA.class, SolverMAM.class};
        } else if (analyzer.hasHomogeneousScheduling(SchedStrategy.INF)) {
            order = new Class<?>[]{SolverMVA.class, SolverNC.class, SolverFluid.class};
        } else {
            order = new Class<?>[]{SolverMVA.class, SolverNC.class, SolverMAM.class,
                    SolverFluid.class, SolverLDES.class};
        }

        NetworkSolver solver = chooseSolverRanked(order);
        if (solver == null) {
            solver = chooseSolverRanked(new Class<?>[]{SolverMVA.class, SolverNC.class,
                    SolverMAM.class, SolverFluid.class, SolverLDES.class,
                    SolverCTMC.class, SolverSSA.class});
        }
        if (solver == null && !candidates.isEmpty()) {
            solver = candidates.get(0);
        }
        return solver;
    }

    /**
     * Metric-aware ranked choice. Mirrors MATLAB's chooseSolverHeur.m: LDES is
     * the only simulation candidate, MVA leads NC leads MAM for analytical
     * solvers (inverted on caches), and Fluid leads where a smooth answer is
     * wanted.
     *
     * @param methodName the accessor being delegated
     * @return the ranked solver for that accessor
     */
    private NetworkSolver chooseSolverHeuristic(String methodName) {
        ModelAnalyzer analyzer = new ModelAnalyzer(model);
        Class<?>[] order;

        switch (methodName) {
            // Average metrics - delegate to chooseAvgSolverHeuristic
            case "getAvgChainTable":
            case "getAvgTputTable":
            case "getAvgRespTTable":
            case "getAvgUtilTable":
            case "getAvgSysTable":
            case "getAvgNodeTable":
            case "getAvgNodeChainTable":
            case "getAvgTable":
            case "getAvg":
            case "getAvgChain":
            case "getAvgSys":
            case "getAvgNode":
            case "getAvgNodeChain":
            case "getAvgArvRChain":
            case "getAvgQLenChain":
            case "getAvgUtilChain":
            case "getAvgRespTChain":
            case "getAvgTputChain":
            case "getAvgSysRespT":
            case "getAvgSysTput":
            case "getAvgQLen":
            case "getAvgUtil":
            case "getAvgRespT":
            case "getAvgResidT":
            case "getAvgWaitT":
            case "getAvgTput":
            case "getAvgArvR":
            case "getAvgQLenTable":
            case "getAvgResidTChain":
            case "getAvgNodeQLenChain":
            case "getAvgNodeUtilChain":
            case "getAvgNodeRespTChain":
            case "getAvgNodeResidTChain":
            case "getAvgNodeTputChain":
            case "getAvgNodeArvRChain":
                return chooseAvgSolverHeuristic();

            case "getTranAvg":
            case "getTranCdfPassT":
            case "getTranCdfRespT":
                order = new Class<?>[]{SolverFluid.class, SolverLDES.class};
                break;

            case "getCdfRespT":
            case "getCdfPassT":
            case "getPerctRespT":
                if (analyzer.hasHomogeneousScheduling(SchedStrategy.FCFS) && analyzer.hasProductForm()) {
                    order = new Class<?>[]{SolverNC.class, SolverFluid.class, SolverLDES.class};
                } else {
                    order = new Class<?>[]{SolverFluid.class, SolverLDES.class};
                }
                break;

            case "getTranProb":
            case "getTranProbSys":
            case "getTranProbAggr":
            case "getTranProbSysAggr":
                order = new Class<?>[]{SolverCTMC.class};
                break;

            case "sample":
            case "sampleSys":
                order = new Class<?>[]{SolverSSA.class, SolverLDES.class};
                break;

            case "sampleAggr":
            case "sampleSysAggr":
                order = new Class<?>[]{SolverLDES.class, SolverSSA.class};
                break;

            case "getProb":
            case "getProbAggr":
            case "getProbSys":
            case "getProbSysAggr":
            case "getProbMarg":
            case "getProbNormConstAggr":
                if (analyzer.hasProductForm()) {
                    order = new Class<?>[]{SolverNC.class, SolverCTMC.class, SolverLDES.class};
                } else {
                    order = new Class<?>[]{SolverCTMC.class, SolverLDES.class};
                }
                break;

            case "getAvgCacheTable":
            case "getAvgItemTable":
                order = new Class<?>[]{SolverNC.class, SolverMVA.class, SolverFluid.class,
                        SolverCTMC.class, SolverLDES.class};
                break;

            case "getAvgLossTable":
            case "getAvgRegionLossTable":
                order = new Class<?>[]{SolverLDES.class, SolverCTMC.class, SolverSSA.class};
                break;

            case "getAvgOrbitTable":
            case "getAvgOrbit":
                order = new Class<?>[]{SolverMVA.class, SolverCTMC.class, SolverLDES.class};
                break;

            case "getMomentTable":
            case "getMomentChainTable":
            case "getMomentStationTable":
                order = new Class<?>[]{SolverMVA.class, SolverNC.class, SolverCTMC.class, SolverLDES.class};
                break;

            case "getSensitivityTable":
                order = new Class<?>[]{SolverFluid.class, SolverMVA.class, SolverNC.class};
                break;

            default:
                return chooseAvgSolverHeuristic();
        }

        NetworkSolver solver = chooseSolverRanked(order);
        if (solver == null) {
            solver = chooseAvgSolverHeuristic();
        }
        return solver;
    }

    /**
     * Ranked choice restricted to exact solvers. Exactness overrides the global
     * MVA &gt; NC order. Mirrors MATLAB's chooseSolverExact.m.
     *
     * @param methodName the accessor being delegated
     * @return the exact solver for that accessor, or null when none is feasible
     */
    private NetworkSolver chooseSolverExact(String methodName) {
        ModelAnalyzer analyzer = new ModelAnalyzer(model);
        Class<?>[] order;
        if ("sample".equals(methodName) || "sampleSys".equals(methodName)
                || "sampleAggr".equals(methodName) || "sampleSysAggr".equals(methodName)) {
            // A sample path is exact in distribution, not in the mean.
            order = new Class<?>[]{SolverSSA.class, SolverLDES.class};
        } else if (methodName.startsWith("getTranProb") || methodName.startsWith("getCdf")
                || methodName.startsWith("getTranCdf") || "getPerctRespT".equals(methodName)
                || "getTranAvg".equals(methodName)) {
            order = new Class<?>[]{SolverCTMC.class};
        } else if (methodName.startsWith("getProb")) {
            if (analyzer.hasProductForm()) {
                order = new Class<?>[]{SolverNC.class, SolverCTMC.class};
            } else {
                order = new Class<?>[]{SolverCTMC.class};
            }
        } else if (analyzer.hasProductForm() && !analyzer.hasMultiServer()) {
            order = new Class<?>[]{SolverNC.class, SolverCTMC.class};
        } else {
            order = new Class<?>[]{SolverCTMC.class, SolverNC.class};
        }
        // Gate on the method-level rule, not just the feature set: NC and MVA
        // reject 'exact' on a non-product-form model, which a flat feature set
        // cannot express.
        return chooseSolverRanked(order, "exact");
    }

    /**
     * Ranked choice restricted to simulators: LDES leads everywhere except
     * event-level sampling, where SSA is the native sample-path engine. Mirrors
     * MATLAB's chooseSolverSim.m.
     *
     * @param methodName the accessor being delegated
     * @return the simulator for that accessor, or null when none is feasible
     */
    private NetworkSolver chooseSolverSim(String methodName) {
        Class<?>[] order;
        if ("sample".equals(methodName) || "sampleSys".equals(methodName)) {
            order = new Class<?>[]{SolverSSA.class, SolverLDES.class};
        } else {
            order = new Class<?>[]{SolverLDES.class, SolverSSA.class};
        }
        return chooseSolverRanked(order);
    }

    /**
     * Select solver using heuristic rules (default for runAnalyzer)
     */
    private void selectSolverHeuristic() {
        selectedSolver = chooseAvgSolverHeuristic();
        if (selectedSolver == null && !candidates.isEmpty()) {
            selectedSolver = candidates.get(0);
        }
    }

    // ========== Delegation Methods ==========

    /**
     * Generic delegation method that tries multiple solvers until one succeeds
     */
    private Object[] delegate(String methodName, int numOutputs, Object... args) {
        List<NetworkSolver> proposedSolvers = new ArrayList<NetworkSolver>();

        // A FAMILY METHOD NAME BYPASSES THE RANKING AND THE FALLBACK SWEEP. The caller
        // named an engine, so answering from another one would report a
        // different algorithm under their name -- the same reason an explicit
        // 'bound' request must not degrade. Every accessor reaches this method
        // rather than runAnalyzer, so the pin has to be honoured here too.
        if (pinnedSolver != null) {
            selectedSolver = pinnedSolver;
            proposedSolvers.add(pinnedSolver);
            return delegateTo(proposedSolvers, methodName, numOutputs, args);
        }

        // An explicit 'bound' request must not silently degrade to a point
        // estimate, so the candidate pool is not consulted.
        if (METHOD_BOUND.equals(selectionMethod)) {
            proposedSolvers.add(getBoundSolver());
        } else if (candidates.size() > 1) {
            NetworkSolver chosenSolver = chooseSolver(methodName);
            if (chosenSolver != null && chosenSolver.supports(model)) {
                proposedSolvers.add(chosenSolver);
            }
            // Add all other candidates
            for (NetworkSolver candidate : candidates) {
                if (!proposedSolvers.contains(candidate)) {
                    proposedSolvers.add(candidate);
                }
            }
        } else {
            // Use the single solver
            proposedSolvers.addAll(candidates);
        }
        
        return delegateTo(proposedSolvers, methodName, numOutputs, args);
    }

    /**
     * Invoke {@code methodName} on the first solver of {@code proposedSolvers}
     * that can answer it.
     *
     * @param proposedSolvers the delegates to try, most preferred first
     * @param methodName      the accessor to invoke
     * @param numOutputs      how many values the caller expects back
     * @param args            the accessor's arguments
     * @return the delegate's answer
     */
    private Object[] delegateTo(List<NetworkSolver> proposedSolvers, String methodName,
                                int numOutputs, Object... args) {
        // Try each solver until one succeeds
        for (NetworkSolver solver : proposedSolvers) {
            try {
                java.lang.reflect.Method method = findMethod(solver.getClass(), methodName, args);
                if (method != null) {
                    Object result = method.invoke(solver, args);
                    selectedSolver = solver;
                    copyResultsFromSolver(solver);
                    
                    if (options.verbose != VerboseLevel.SILENT) {
                        System.out.println("Successful method execution completed by " + solver.getName());
                    }
                    
                    // Handle different return types
                    if (numOutputs == 1) {
                        return new Object[]{result};
                    } else {
                        // For multiple outputs, assume result is an array or return single result
                        if (result instanceof Object[]) {
                            return (Object[]) result;
                        } else {
                            return new Object[]{result};
                        }
                    }
                }
            } catch (Exception e) {
                if (e.getMessage() != null && e.getMessage().contains("Unrecognized method")) {
                    line_warning("SolverAUTO.invokeMethod", "Method unsupported by %s", solver.getName());
                } else {
                    line_warning("SolverAUTO.invokeMethod", "Error in %s: %s", solver.getName(), e.getMessage());
                }
            }
        }

        throw new RuntimeException("No solver could execute method: " + methodName);
    }
    
    /**
     * Find method by name and compatible parameter types
     */
    private java.lang.reflect.Method findMethod(Class<?> clazz, String methodName, Object... args) {
        java.lang.reflect.Method[] methods = clazz.getMethods();
        for (java.lang.reflect.Method method : methods) {
            if (method.getName().equals(methodName)) {
                Class<?>[] paramTypes = method.getParameterTypes();
                if (paramTypes.length == args.length) {
                    boolean compatible = true;
                    for (int i = 0; i < args.length; i++) {
                        if (args[i] != null && !paramTypes[i].isAssignableFrom(args[i].getClass())) {
                            // Check for primitive type compatibility
                            if (!isPrimitiveCompatible(paramTypes[i], args[i].getClass())) {
                                compatible = false;
                                break;
                            }
                        }
                    }
                    if (compatible) {
                        return method;
                    }
                }
            }
        }
        return null;
    }
    
    /**
     * Check if primitive types are compatible
     */
    private boolean isPrimitiveCompatible(Class<?> paramType, Class<?> argType) {
        if (paramType == int.class && argType == Integer.class) return true;
        if (paramType == double.class && argType == Double.class) return true;
        if (paramType == boolean.class && argType == Boolean.class) return true;
        if (paramType == long.class && argType == Long.class) return true;
        return false;
    }
    
    /**
     * Choose solver for a specific method based on selection method.
     * Matches MATLAB's chooseSolver.m
     */
    private NetworkSolver chooseSolver(String methodName) {
        // If a solver is forced, always use it
        if (autoOptions.forceSolver != null && !autoOptions.forceSolver.isEmpty()) {
            Integer solverId = solverNameToId.get(autoOptions.forceSolver.toLowerCase());
            if (solverId != null) {
                // Find the solver by its type based on the ID
                switch (solverId) {
                    case CANDIDATE_MVA:
                        return findSolverByType(SolverMVA.class);
                    case CANDIDATE_NC:
                        return findSolverByType(SolverNC.class);
                    case CANDIDATE_MAM:
                        return findSolverByType(SolverMAM.class);
                    case CANDIDATE_FLUID:
                        return findSolverByType(SolverFluid.class);
                    case CANDIDATE_SSA:
                        return findSolverByType(SolverSSA.class);
                    case CANDIDATE_CTMC:
                        return findSolverByType(SolverCTMC.class);
                }
            }
        }

        // Otherwise use selection method
        switch (selectionMethod) {
            // case METHOD_AI:  // AI method not yet available
            //     return chooseSolverHeuristic(methodName);
            case METHOD_BOUND:
                return getBoundSolver();
            case METHOD_EXACT: {
                NetworkSolver exact = chooseSolverExact(methodName);
                if (exact == null) {
                    throw new RuntimeException("SolverAuto: no exact solver supports this model for "
                            + methodName + "; use method 'default' for the approximate heuristic");
                }
                return exact;
            }
            case METHOD_SIM: {
                NetworkSolver sim = chooseSolverSim(methodName);
                return sim != null ? sim : chooseSolverHeuristic(methodName);
            }
            case METHOD_FAST: {
                NetworkSolver fast = chooseSolverRanked(new Class<?>[]{SolverMVA.class,
                        SolverNC.class, SolverFluid.class, SolverMAM.class});
                return fast != null ? fast : chooseSolverHeuristic(methodName);
            }
            case METHOD_ACCURATE: {
                NetworkSolver acc = chooseSolverRanked(new Class<?>[]{SolverFluid.class,
                        SolverMAM.class, SolverCTMC.class, SolverLDES.class});
                return acc != null ? acc : chooseSolverHeuristic(methodName);
            }
            case METHOD_HEURISTIC:
            case METHOD_DEFAULT:
            default:
                return chooseSolverHeuristic(methodName);
        }
    }

    @Override
    public boolean supports(Network model) {
        // AUTO solver supports any model that at least one candidate supports
        return !candidates.isEmpty();
    }

    // ========== Basic Network Analysis Methods ==========

    public NetworkAvgChainTable getAvgChainTable() {
        Object[] results = delegate("getAvgChainTable", 1);
        if (results.length > 0 && results[0] instanceof NetworkAvgChainTable) {
            return (NetworkAvgChainTable) results[0];
        }
        // Fallback to empty table if delegation fails
        return new NetworkAvgChainTable(new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    public NetworkAvgSysTable getAvgSysTable() {
        Object[] results = delegate("getAvgSysTable", 1);
        if (results.length > 0 && results[0] instanceof NetworkAvgSysTable) {
            return (NetworkAvgSysTable) results[0];
        }
        // Fallback to empty table if delegation fails
        return new NetworkAvgSysTable(new ArrayList<Double>(), new ArrayList<Double>(), this.options);
    }

    public NetworkAvgNodeTable getAvgNodeTable() {
        Object[] results = delegate("getAvgNodeTable", 1);
        if (results.length > 0 && results[0] instanceof NetworkAvgNodeTable) {
            return (NetworkAvgNodeTable) results[0];
        }
        // Fallback to empty table if delegation fails
        return new NetworkAvgNodeTable(new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    public NetworkAvgTable getAvgTable() {
        Object[] results = delegate("getAvgTable", 1);
        if (results.length > 0 && results[0] instanceof NetworkAvgTable) {
            return (NetworkAvgTable) results[0];
        }
        // Fallback to empty table if delegation fails
        return new NetworkAvgTable(new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    @Override
    public SolverResult getAvg() {
        Object[] results = delegate("getAvg", 1);
        SolverResult result = new SolverResult();
        if (results.length > 0) {
            result.QN = (Matrix) results[0];
        }
        return result;
    }

    public SolverResult getAvgChain() {
        Object[] results = delegate("getAvgChain", 4);
        SolverResult result = new SolverResult();
        result.QN = (Matrix) results[0];
        result.UN = (Matrix) results[1];
        result.RN = (Matrix) results[2];
        result.TN = (Matrix) results[3];
        return result;
    }

    public void getAvgSys() {
        delegate("getAvgSys", 2);
    }

    public SolverResult getAvgNode() {
        Object[] results = delegate("getAvgNode", 6);
        SolverResult result = new SolverResult();
        result.QN = (Matrix) results[0];
        result.UN = (Matrix) results[1];
        result.RN = (Matrix) results[2];
        result.TN = (Matrix) results[3];
        result.AN = (Matrix) results[4];
        result.WN = (Matrix) results[5];
        return result;
    }

    // ========== Chain Methods ==========

    public Matrix getAvgArvRChain() {
        Object[] results = delegate("getAvgArvRChain", 1);
        return (Matrix) results[0];
    }

    public Matrix getAvgQLenChain() {
        Object[] results = delegate("getAvgQLenChain", 1);
        return (Matrix) results[0];
    }

    public Matrix getAvgUtilChain() {
        Object[] results = delegate("getAvgUtilChain", 1);
        return (Matrix) results[0];
    }

    public Matrix getAvgRespTChain() {
        Object[] results = delegate("getAvgRespTChain", 1);
        return (Matrix) results[0];
    }

    public Matrix getAvgTputChain() {
        Object[] results = delegate("getAvgTputChain", 1);
        return (Matrix) results[0];
    }

    // ========== System Methods ==========

    public Matrix getAvgSysRespT() {
        Object[] results = delegate("getAvgSysRespT", 1);
        return (Matrix) results[0];
    }

    public Matrix getAvgSysTput() {
        Object[] results = delegate("getAvgSysTput", 1);
        return (Matrix) results[0];
    }

    // ========== Transient Analysis Methods ==========

    public void getTranAvg() {
        delegate("getTranAvg", 3);
    }

    public DistributionResult getTranCdfPassT() {
        Object[] results = delegate("getTranCdfPassT", 1);
        return (DistributionResult) results[0];
    }

    public DistributionResult getTranCdfRespT() {
        Object[] results = delegate("getTranCdfRespT", 1);
        return (DistributionResult) results[0];
    }

    // ========== Probability Methods ==========

    public Matrix[] getTranProb(Node node) {
        Object[] results = delegate("getTranProb", 2, node);
        return new Matrix[]{(Matrix) results[0], (Matrix) results[1]};
    }

    public Matrix[] getTranProbAggr(Node node) {
        Object[] results = delegate("getTranProbAggr", 2, node);
        return new Matrix[]{(Matrix) results[0], (Matrix) results[1]};
    }

    public Matrix[] getTranProbSys() {
        Object[] results = delegate("getTranProbSys", 2);
        return new Matrix[]{(Matrix) results[0], (Matrix) results[1]};
    }

    public Matrix[] getTranProbSysAggr() {
        Object[] results = delegate("getTranProbSysAggr", 2);
        return new Matrix[]{(Matrix) results[0], (Matrix) results[1]};
    }

    // ========== Sampling Methods ==========

    public SampleResult sample(Node node, int numEvents) {
        Object[] results = delegate("sample", 1, node, numEvents);
        return (SampleResult) results[0];
    }

    public SampleResult sampleAggr(Node node, int numEvents) {
        Object[] results = delegate("sampleAggr", 1, node, numEvents);
        return (SampleResult) results[0];
    }

    public SampleResult sampleSys(int numEvents) {
        Object[] results = delegate("sampleSys", 1, numEvents);
        return (SampleResult) results[0];
    }

    public SampleResult sampleSysAggr(int numEvents) {
        Object[] results = delegate("sampleSysAggr", 1, numEvents);
        return (SampleResult) results[0];
    }

    // ========== Distribution Methods ==========

    public DistributionResult getCdfRespT() {
        Object[] results = delegate("getCdfRespT", 1);
        return (DistributionResult) results[0];
    }

    /**
     * System response time CDF per chain, always from SolverCTMC -- like the
     * state-space accessors, it is a property of the CTMC representation, so it
     * resolves on the CTMC candidate rather than the ranked selection, as the
     * reference {@code @SolverAUTO/getCdfSysRespT} does.
     */
    public java.util.List<Matrix> getCdfSysRespT() {
        return ctmcSolver().getCdfSysRespT();
    }

    /**
     * A SolverFluid over the same model, reusing the fluid candidate when one
     * exists -- the {@code fldSolver()} helper of the reference.
     */
    private SolverFluid fluidSolver() {
        for (NetworkSolver solver : candidates) {
            if (solver instanceof SolverFluid) {
                return (SolverFluid) solver;
            }
        }
        return new SolverFluid(model);
    }

    /** Backward-compatible passage-time CDF name, resolved on the fluid family. */
    public DistributionResult getCdfPT() {
        return fluidSolver().getCdfPT();
    }

    /** Age-of-Information CDFs, a property of the ODE representation, so fluid-only. */
    public Matrix[] getCdfAoI() {
        return fluidSolver().getCdfAoI();
    }

    /** Age-of-Information CDFs on caller-supplied time points, fluid-only. */
    public Matrix[] getCdfAoI(Matrix tValues) {
        return fluidSolver().getCdfAoI(tValues);
    }

    public double getProb(Node node, Matrix state) {
        Object[] results = delegate("getProb", 1, node, state);
        return (Double) results[0];
    }

    public double getProbAggr(Node node, Matrix state_a) {
        Object[] results = delegate("getProbAggr", 1, node, state_a);
        return (Double) results[0];
    }

    public ProbabilityResult getProbSys() {
        Object[] results = delegate("getProbSys", 1);
        return new ProbabilityResult((Double) results[0]);
    }

    public ProbabilityResult getProbSysAggr() {
        Object[] results = delegate("getProbSysAggr", 1);
        return new ProbabilityResult((Double) results[0]);
    }

    public ProbabilityResult getProbNormConstAggr() {
        Object[] results = delegate("getProbNormConstAggr", 1);
        return new ProbabilityResult((Double) results[0]);
    }

    // ========== Basic Metric Methods ==========

    /**
     * Get average queue lengths at steady-state
     */
    public Matrix getAvgQLen() {
        Object[] results = delegate("getAvgQLen", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average utilizations at steady-state
     */
    public Matrix getAvgUtil() {
        Object[] results = delegate("getAvgUtil", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average response times at steady-state
     */
    public Matrix getAvgRespT() {
        Object[] results = delegate("getAvgRespT", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average residence times at steady-state
     */
    public Matrix getAvgResidT() {
        Object[] results = delegate("getAvgResidT", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average waiting times (queue time excluding service)
     */
    public Matrix getAvgWaitT() {
        Object[] results = delegate("getAvgWaitT", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average throughputs at steady-state
     */
    public Matrix getAvgTput() {
        Object[] results = delegate("getAvgTput", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average arrival rates at steady-state
     */
    public Matrix getAvgArvR() {
        Object[] results = delegate("getAvgArvR", 1);
        return (Matrix) results[0];
    }

    // ========== Additional Table Methods ==========

    /**
     * Get average queue length table
     */
    public NetworkAvgTable getAvgQLenTable() {
        Object[] results = delegate("getAvgQLenTable", 1);
        if (results.length > 0 && results[0] instanceof NetworkAvgTable) {
            return (NetworkAvgTable) results[0];
        }
        return new NetworkAvgTable(new ArrayList<Double>(), new ArrayList<Double>(),
            new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    /**
     * Get average utilization table
     */
    public NetworkAvgTable getAvgUtilTable() {
        Object[] results = delegate("getAvgUtilTable", 1);
        if (results.length > 0 && results[0] instanceof NetworkAvgTable) {
            return (NetworkAvgTable) results[0];
        }
        return new NetworkAvgTable(new ArrayList<Double>(), new ArrayList<Double>(),
            new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    /**
     * Get average response time table
     */
    public NetworkAvgTable getAvgRespTTable() {
        Object[] results = delegate("getAvgRespTTable", 1);
        if (results.length > 0 && results[0] instanceof NetworkAvgTable) {
            return (NetworkAvgTable) results[0];
        }
        return new NetworkAvgTable(new ArrayList<Double>(), new ArrayList<Double>(),
            new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    /**
     * Get average throughput table
     */
    public NetworkAvgTable getAvgTputTable() {
        Object[] results = delegate("getAvgTputTable", 1);
        if (results.length > 0 && results[0] instanceof NetworkAvgTable) {
            return (NetworkAvgTable) results[0];
        }
        return new NetworkAvgTable(new ArrayList<Double>(), new ArrayList<Double>(),
            new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    // ========== Additional Chain Methods ==========

    /**
     * Get average residence time by chain
     */
    public Matrix getAvgResidTChain() {
        Object[] results = delegate("getAvgResidTChain", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average node queue length by chain
     */
    public Matrix getAvgNodeQLenChain() {
        Object[] results = delegate("getAvgNodeQLenChain", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average node utilization by chain
     */
    public Matrix getAvgNodeUtilChain() {
        Object[] results = delegate("getAvgNodeUtilChain", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average node response time by chain
     */
    public Matrix getAvgNodeRespTChain() {
        Object[] results = delegate("getAvgNodeRespTChain", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average node residence time by chain
     */
    public Matrix getAvgNodeResidTChain() {
        Object[] results = delegate("getAvgNodeResidTChain", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average node throughput by chain
     */
    public Matrix getAvgNodeTputChain() {
        Object[] results = delegate("getAvgNodeTputChain", 1);
        return (Matrix) results[0];
    }

    /**
     * Get average node arrival rate by chain
     */
    public Matrix getAvgNodeArvRChain() {
        Object[] results = delegate("getAvgNodeArvRChain", 1);
        return (Matrix) results[0];
    }

    // ========== Additional Distribution Methods ==========

    /**
     * Get CDF of passage times at steady-state
     */
    public DistributionResult getCdfPassT() {
        Object[] results = delegate("getCdfPassT", 1);
        return (DistributionResult) results[0];
    }

    /**
     * Get response time percentiles
     */
    public Matrix getPerctRespT(double[] percentiles) {
        Object[] results = delegate("getPerctRespT", 1, percentiles);
        return (Matrix) results[0];
    }

    /**
     * Get marginalized state probability
     */
    public double getProbMarg(Node node, int jobclass, Matrix state_m) {
        Object[] results = delegate("getProbMarg", 1, node, jobclass, state_m);
        return (Double) results[0];
    }

    // ========== Aliases ==========
    // These provide shorter method names without the 'get' prefix

    /** Alias for getAvgTable() */
    public NetworkAvgTable avgTable() { return getAvgTable(); }

    /** Alias for getAvgSysTable() */
    public NetworkAvgSysTable avgSysTable() { return getAvgSysTable(); }

    /** Alias for getAvgNodeTable() */
    public NetworkAvgNodeTable avgNodeTable() { return getAvgNodeTable(); }

    /** Alias for getAvgChainTable() */
    public NetworkAvgChainTable avgChainTable() { return getAvgChainTable(); }

    /** Alias for getAvg() */
    public SolverResult avg() { return getAvg(); }

    /** Alias for getAvgChain() */
    public SolverResult avgChain() { return getAvgChain(); }

    /** Alias for getAvgNode() */
    public SolverResult avgNode() { return getAvgNode(); }

    /** Alias for getAvgQLen() */
    public Matrix avgQLen() { return getAvgQLen(); }

    /** Alias for getAvgUtil() */
    public Matrix avgUtil() { return getAvgUtil(); }

    /** Alias for getAvgRespT() */
    public Matrix avgRespT() { return getAvgRespT(); }

    /** Alias for getAvgResidT() */
    public Matrix avgResidT() { return getAvgResidT(); }

    /** Alias for getAvgWaitT() */
    public Matrix avgWaitT() { return getAvgWaitT(); }

    /** Alias for getAvgTput() */
    public Matrix avgTput() { return getAvgTput(); }

    /** Alias for getAvgArvR() */
    public Matrix avgArvR() { return getAvgArvR(); }

    /** Alias for getAvgQLenChain() */
    public Matrix avgQLenChain() { return getAvgQLenChain(); }

    /** Alias for getAvgUtilChain() */
    public Matrix avgUtilChain() { return getAvgUtilChain(); }

    /** Alias for getAvgRespTChain() */
    public Matrix avgRespTChain() { return getAvgRespTChain(); }

    /** Alias for getAvgResidTChain() */
    public Matrix avgResidTChain() { return getAvgResidTChain(); }

    /** Alias for getAvgTputChain() */
    public Matrix avgTputChain() { return getAvgTputChain(); }

    /** Alias for getAvgArvRChain() */
    public Matrix avgArvRChain() { return getAvgArvRChain(); }

    /** Alias for getAvgSysRespT() */
    public Matrix avgSysRespT() { return getAvgSysRespT(); }

    /** Alias for getAvgSysTput() */
    public Matrix avgSysTput() { return getAvgSysTput(); }

    /** Alias for getCdfRespT() */
    public DistributionResult cdfRespT() { return getCdfRespT(); }

    /** Alias for getCdfPassT() */
    public DistributionResult cdfPassT() { return getCdfPassT(); }

    /** Alias for getPerctRespT() */
    public Matrix perctRespT(double[] percentiles) { return getPerctRespT(percentiles); }

    /** Alias for getProb() */
    public double prob(Node node, Matrix state) { return getProb(node, state); }

    /** Alias for getProbAggr() */
    public double probAggr(Node node, Matrix state_a) { return getProbAggr(node, state_a); }

    /** Alias for getProbSys() */
    public ProbabilityResult probSys() { return getProbSys(); }

    /** Alias for getProbSysAggr() */
    public ProbabilityResult probSysAggr() { return getProbSysAggr(); }

    /** Alias for getProbNormConstAggr() */
    public ProbabilityResult probNormConstAggr() { return getProbNormConstAggr(); }

    /** Alias for getProbMarg() */
    public double probMarg(Node node, int jobclass, Matrix state_m) { return getProbMarg(node, jobclass, state_m); }

    /** Alias for sample() */
    public SampleResult sample(Node node) { return sample(node, 1000); }

    /** Alias for sampleSys() */
    public SampleResult sampleSys() { return sampleSys(1000); }
}