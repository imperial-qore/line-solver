package jline.solvers.auto;

import static jline.io.InputOutputKt.line_warning;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.VerboseLevel;
import jline.lang.nodes.Node;
import jline.solvers.NetworkAvgChainTable;
import jline.solvers.NetworkAvgNodeTable;
import jline.solvers.NetworkAvgSysTable;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.des.SolverDES;
import jline.solvers.env.SolverENV;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;
import jline.lang.constant.SchedStrategy;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Automatic solver selection for queueing network models.
 * <p>
 * This solver automatically selects the most appropriate solution method
 * based on model characteristics and requested performance metrics.
 */
public class SolverAUTO extends NetworkSolver {

    public static final int CANDIDATE_CTMC = 6;
    public static final int CANDIDATE_DES = 7;
    public static final int CANDIDATE_FLUID = 3;
    public static final int CANDIDATE_JMT = 4;
    public static final int CANDIDATE_MAM = 2;
    // Solver candidates constants
    public static final int CANDIDATE_MVA = 0;
    public static final int CANDIDATE_NC = 1;
    public static final int CANDIDATE_SSA = 5;
    // public static final String METHOD_AI = "ai";  // AI method not yet available
    // Method selection strategies
    public static final String METHOD_DEFAULT = "default";
    public static final String METHOD_HEURISTIC = "heur";
    public static final String METHOD_SIM = "sim";
    public static final String METHOD_EXACT = "exact";
    public static final String METHOD_FAST = "fast";
    public static final String METHOD_ACCURATE = "accurate";

    // Solver instances
    private List<NetworkSolver> candidates;
    private Map<String, Integer> solverNameToId;
    private NetworkSolver selectedSolver;
    private final String selectionMethod;

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
     * Get the name of the selected solver
     */
    public String getSelectedSolverName() {
        return selectedSolver != null ? selectedSolver.getName() : "none";
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

        candidates.add(new SolverJMT(model));
        solverNameToId.put("jmt", CANDIDATE_JMT);

        candidates.add(new SolverSSA(model));
        solverNameToId.put("ssa", CANDIDATE_SSA);

        // CTMC solver integration
        candidates.add(new SolverCTMC(model));
        solverNameToId.put("ctmc", CANDIDATE_CTMC);

        // DES solver integration
        candidates.add(new SolverDES(model));
        solverNameToId.put("des", CANDIDATE_DES);

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
        ensureSolverSelected();

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
        }

        // Try other candidates if the selected solver fails
        for (NetworkSolver candidate : candidates) {
            if (candidate != selectedSolver) {
                try {
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
            selectForcedSolver(autoOptions.forceSolver);
        } else {
            // Automatic selection based on method
            switch (selectionMethod) {
                case METHOD_HEURISTIC:
                case METHOD_DEFAULT:
                    selectSolverHeuristic();
                    break;
                // case METHOD_AI:  // AI method not yet available
                //     selectSolverAI();
                //     break;
                case METHOD_SIM:
                    // Best simulator - JMT
                    selectedSolver = findSolverByType(SolverJMT.class);
                    if (selectedSolver == null) {
                        selectedSolver = findSolverByType(SolverSSA.class);
                    }
                    if (selectedSolver == null) {
                        selectSolverHeuristic();
                    }
                    break;
                case METHOD_EXACT:
                    // Best exact method - NC for closed networks, CTMC otherwise
                    NetworkStruct sn = model.getStruct(true);
                    boolean isClosed = true;
                    for (int i = 0; i < sn.nclasses; i++) {
                        if (Double.isInfinite(sn.njobs.get(i))) {
                            isClosed = false;
                            break;
                        }
                    }
                    if (isClosed) {
                        selectedSolver = findSolverByType(SolverNC.class);
                    } else {
                        selectedSolver = findSolverByType(SolverCTMC.class);
                    }
                    if (selectedSolver == null) {
                        selectSolverHeuristic();
                    }
                    break;
                case METHOD_FAST:
                    // Fast approximate method - MVA
                    selectedSolver = findSolverByType(SolverMVA.class);
                    if (selectedSolver == null) {
                        selectSolverHeuristic();
                    }
                    break;
                case METHOD_ACCURATE:
                    // Accurate approximate method - Fluid
                    selectedSolver = findSolverByType(SolverFluid.class);
                    if (selectedSolver == null) {
                        selectedSolver = findSolverByType(SolverMAM.class);
                    }
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
     * Select solver using heuristic rules for average metrics.
     * Matches MATLAB's chooseAvgSolverHeur.m
     */
    private NetworkSolver chooseAvgSolverHeuristic() {
        ModelAnalyzer analyzer = new ModelAnalyzer(model);
        NetworkSolver solver = null;

        // For single chain networks - use NC for exact O(1) solution
        if (analyzer.hasSingleChain()) {
            solver = findSolverByType(SolverNC.class);
            if (solver != null) return solver;
        }

        // For multi-chain networks
        if (analyzer.hasMultiChain()) {
            // All infinite servers - use MVA
            if (analyzer.hasHomogeneousScheduling(SchedStrategy.INF)) {
                solver = findSolverByType(SolverMVA.class);
                if (solver != null) return solver;
            }

            // Product form without multi-servers - use NC
            if (analyzer.hasProductForm() && !analyzer.hasMultiServer()) {
                solver = findSolverByType(SolverNC.class);
                if (solver != null) return solver;
            }

            // FCFS without multi-servers
            if (analyzer.hasHomogeneousScheduling(SchedStrategy.FCFS) && !analyzer.hasMultiServer()) {
                double avgJobsPerChain = analyzer.getAvgJobsPerChain();
                int totalJobs = analyzer.getTotalJobs();

                if (avgJobsPerChain > 30) {
                    // Heavy load - use FLUID (likely fluid regime)
                    solver = findSolverByType(SolverFluid.class);
                    if (solver != null) return solver;
                } else if (avgJobsPerChain > 10) {
                    // Medium/heavy load - use MVA
                    solver = findSolverByType(SolverMVA.class);
                    if (solver != null) return solver;
                } else if (totalJobs < 5) {
                    // Light load - use NC (avoids AMVA errors)
                    solver = findSolverByType(SolverNC.class);
                    if (solver != null) return solver;
                }
            }

            // PS/PSPRIO with multi-servers - use MVA
            if ((analyzer.hasHomogeneousScheduling(SchedStrategy.PS) ||
                 analyzer.hasHomogeneousScheduling(SchedStrategy.PSPRIO)) && analyzer.hasMultiServer()) {
                solver = findSolverByType(SolverMVA.class);
                if (solver != null) return solver;
            }

            // FCFS with multi-servers
            if (analyzer.hasHomogeneousScheduling(SchedStrategy.FCFS) && analyzer.hasMultiServer()) {
                int totalJobs = analyzer.getTotalJobs();

                if (totalJobs < 5) {
                    // Light load - use NC
                    solver = findSolverByType(SolverNC.class);
                    if (solver != null) return solver;
                } else {
                    // Use MVA
                    solver = findSolverByType(SolverMVA.class);
                    if (solver != null) return solver;
                }
            }
        }

        // Default fallback - prefer MVA
        solver = findSolverByType(SolverMVA.class);
        if (solver != null) return solver;

        // Fallback to first available
        if (!candidates.isEmpty()) {
            return candidates.get(0);
        }
        return null;
    }

    /**
     * Select solver using heuristic rules based on the method being called.
     * Matches MATLAB's chooseSolverHeur.m
     */
    private NetworkSolver chooseSolverHeuristic(String methodName) {
        ModelAnalyzer analyzer = new ModelAnalyzer(model);
        NetworkSolver solver = null;

        switch (methodName) {
            // Average metrics - delegate to chooseAvgSolverHeuristic
            case "getAvgChainTable":
            case "getAvgTputTable":
            case "getAvgRespTTable":
            case "getAvgUtilTable":
            case "getAvgSysTable":
            case "getAvgNodeTable":
            case "getAvgTable":
            case "getAvg":
            case "getAvgChain":
            case "getAvgSys":
            case "getAvgNode":
            case "getAvgArvRChain":
            case "getAvgQLenChain":
            case "getAvgUtilChain":
            case "getAvgRespTChain":
            case "getAvgTputChain":
            case "getAvgSysRespT":
            case "getAvgSysTput":
            // Basic metric methods
            case "getAvgQLen":
            case "getAvgUtil":
            case "getAvgRespT":
            case "getAvgResidT":
            case "getAvgWaitT":
            case "getAvgTput":
            case "getAvgArvR":
            // Additional table methods
            case "getAvgQLenTable":
            // Additional chain methods
            case "getAvgResidTChain":
            case "getAvgNodeQLenChain":
            case "getAvgNodeUtilChain":
            case "getAvgNodeRespTChain":
            case "getAvgNodeResidTChain":
            case "getAvgNodeTputChain":
            case "getAvgNodeArvRChain":
            // Percentiles
            case "getPerctRespT":
                return chooseAvgSolverHeuristic();

            // Transient average - Fluid if supported, else JMT
            case "getTranAvg":
                solver = findSolverByType(SolverFluid.class);
                if (solver != null && solver.supports(model)) {
                    return solver;
                }
                solver = findSolverByType(SolverJMT.class);
                if (solver != null) return solver;
                break;

            // CDF response time - NC if FCFS+product form, else Fluid, else JMT
            case "getCdfRespT":
                if (analyzer.hasHomogeneousScheduling(SchedStrategy.FCFS) &&
                    analyzer.hasProductForm()) {
                    solver = findSolverByType(SolverNC.class);
                    if (solver != null && solver.supports(model)) {
                        return solver;
                    }
                }
                solver = findSolverByType(SolverFluid.class);
                if (solver != null && solver.supports(model)) {
                    return solver;
                }
                solver = findSolverByType(SolverJMT.class);
                if (solver != null) return solver;
                break;

            // Transient CDF pass/response time - Fluid if supported, else JMT
            case "getTranCdfPassT":
            case "getTranCdfRespT":
                solver = findSolverByType(SolverFluid.class);
                if (solver != null && solver.supports(model)) {
                    return solver;
                }
                solver = findSolverByType(SolverJMT.class);
                if (solver != null) return solver;
                break;

            // Transient probability - use CTMC
            case "getTranProb":
            case "getTranProbSys":
            case "getTranProbAggr":
            case "getTranProbSysAggr":
                solver = findSolverByType(SolverCTMC.class);
                if (solver != null) return solver;
                break;

            // Sampling - use SSA
            case "sample":
            case "sampleSys":
                solver = findSolverByType(SolverSSA.class);
                if (solver != null) return solver;
                break;

            // Aggregate sampling - use JMT
            case "sampleAggr":
            case "sampleSysAggr":
                solver = findSolverByType(SolverJMT.class);
                if (solver != null) return solver;
                break;

            // Probability - NC if product form, else JMT
            case "getProb":
            case "getProbAggr":
            case "getProbSys":
            case "getProbSysAggr":
            case "getProbNormConstAggr":
                if (analyzer.hasProductForm()) {
                    solver = findSolverByType(SolverNC.class);
                    if (solver != null && solver.supports(model)) {
                        return solver;
                    }
                }
                solver = findSolverByType(SolverJMT.class);
                if (solver != null) return solver;
                break;

            default:
                // For unknown methods, use average heuristic as default
                return chooseAvgSolverHeuristic();
        }

        // Fallback to first available solver
        if (!candidates.isEmpty()) {
            return candidates.get(0);
        }
        return null;
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
        
        // Add chosen solver first if we have multiple candidates
        if (candidates.size() > 1) {
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
                    case CANDIDATE_JMT:
                        return findSolverByType(SolverJMT.class);
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
        // Convert matrix to appropriate format - simplified implementation
        return new NetworkAvgChainTable(new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>(), new ArrayList<Double>());
    }

    public NetworkAvgSysTable getAvgSysTable() {
        Object[] results = delegate("getAvgSysTable", 1);
        // Convert matrix to appropriate format - simplified implementation
        return new NetworkAvgSysTable(new ArrayList<Double>(), new ArrayList<Double>(), this.options);
    }

    public NetworkAvgNodeTable getAvgNodeTable() {
        Object[] results = delegate("getAvgNodeTable", 1);
        // Convert matrix to appropriate format - simplified implementation
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

    public SampleResult sample(Node node, int numSamples) {
        Object[] results = delegate("sample", 1, node, numSamples);
        return (SampleResult) results[0];
    }

    public SampleResult sampleAggr(Node node, int numSamples) {
        Object[] results = delegate("sampleAggr", 1, node, numSamples);
        return (SampleResult) results[0];
    }

    public SampleResult sampleSys(int numSamples) {
        Object[] results = delegate("sampleSys", 1, numSamples);
        return (SampleResult) results[0];
    }

    public SampleResult sampleSysAggr(int numSamples) {
        Object[] results = delegate("sampleSysAggr", 1, numSamples);
        return (SampleResult) results[0];
    }

    // ========== Distribution Methods ==========

    public DistributionResult getCdfRespT() {
        Object[] results = delegate("getCdfRespT", 1);
        return (DistributionResult) results[0];
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

    // ========== Kotlin-style Aliases ==========
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