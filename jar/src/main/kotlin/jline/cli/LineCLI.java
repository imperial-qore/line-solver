/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.cli;

import static jline.io.InputOutputKt.line_warning;

import jline.lang.Model;
import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.lang.nodes.Node;
import jline.solvers.AvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ctmc.analyzers.RewardResult;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.ln.SolverLN;
import jline.solvers.lqns.SolverLQNS;
import jline.solvers.mva.SolverMVA;
import jline.solvers.mam.SolverMAM;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.ssa.SampleNodeState;
import jline.solvers.ssa.SampleSysState;
import jline.io.M2M;
import jline.util.matrix.Matrix;
import org.apache.commons.io.FilenameUtils;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import jline.util.RandomManager;
import java.util.*;
import java.util.Scanner;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * The LineCLI class provides a command-line interface for configuring and running the LINE Solver.
 * It supports various options for input and output formats, solvers, analysis types, and more.
 * This class includes a static method to parse command-line arguments and set the appropriate
 * configuration options.
 */
public class LineCLI {
    /**
     * Default constructor for the LineCLI class.
     */
    public LineCLI() {
    }

    /**
     * Prints detailed help information for all runtime options.
     */
    private static void printDetailedHelp() {
        System.out.println("====================================================================");
        System.out.println("LINE Solver - Command Line Interface");
        System.out.println("Copyright (c) 2012-2026, QORE Lab, Imperial College London");
        System.out.printf("Version %s. All rights reserved.%n", new Model("").getVersion());
        System.out.println("====================================================================");
        System.out.println();
        System.out.println("USAGE:");
        System.out.println("  java -jar jline.jar [OPTIONS]");
        System.out.println("  cat model.jsimg | java -jar jline.jar [OPTIONS]");
        System.out.println();
        System.out.println("RUNTIME OPTIONS:");
        System.out.println();
        
        System.out.println("-p, --port <port>");
        System.out.println("    Run LINE solver in server mode on the specified port.");
        System.out.println("    The server will listen for WebSocket connections and process");
        System.out.println("    model solving requests. Type 'q' to quit server mode.");
        System.out.println("    Default: 5863");
        System.out.println();
        
        System.out.println("-f, --file <filepath>");
        System.out.println("    Specify input model file path. If not provided, the solver");
        System.out.println("    will read from standard input (stdin).");
        System.out.println("    Example: -f /path/to/model.jsimg");
        System.out.println();
        
        System.out.println("-i, --input <format>");
        System.out.println("    Input file format. Supported formats:");
        System.out.println("    • jsim   - JSIM format (default)");
        System.out.println("    • jsimg  - JSIM graphics format");
        System.out.println("    • jsimw  - JSIM workspace format");
        System.out.println("    • lqnx   - LQN XML format");
        System.out.println("    • xml    - XML format");
        System.out.println("    Default: jsim");
        System.out.println();
        
        System.out.println("-o, --output <format>");
        System.out.println("    Output format for results. Supported formats:");
        System.out.println("    • readable - Human-readable text format (default)");
        System.out.println("    • json     - JSON format");
        System.out.println("    Default: readable");
        System.out.println();
        
        System.out.println("-s, --solver <solver>");
        System.out.println("    Solver algorithm to use. Available solvers:");
        System.out.println("    • auto   - Automatic solver selection based on input format");
        System.out.println("    • mva    - Mean Value Analysis (default)");
        System.out.println("    • ctmc   - Continuous-Time Markov Chain");
        System.out.println("    • fluid  - Fluid/Mean-Field ODE Solver");
        System.out.println("    • jmt    - Java Modelling Tools simulation");
        System.out.println("    • mam    - Matrix Analytic Methods");
        System.out.println("    • nc     - Normalizing Constant Analyzer");
        System.out.println("    • ssa    - Stochastic Simulation Algorithm");
        System.out.println("    • ln     - Layered Network solver");
        System.out.println("    • lqns   - LQN solver (for LQN models only)");
        System.out.println("    Default: mva");
        System.out.println();
        
        System.out.println("-a, --analysis <type>");
        System.out.println("    Analysis type(s) to perform. Can be comma-separated for multiple.");
        System.out.println("    Basic:");
        System.out.println("    • all         - Both avg and sys metrics (default)");
        System.out.println("    • avg         - Average performance metrics");
        System.out.println("    • sys         - System-level metrics");
        System.out.println("    • stage       - Stage-based metrics (multi-stage service)");
        System.out.println("    • chain       - Chain-level averages");
        System.out.println("    • node        - Node-level averages");
        System.out.println("    • nodechain   - Node-chain level averages");
        System.out.println("    Distribution:");
        System.out.println("    • cdf-respt   - Response time CDF");
        System.out.println("    • cdf-passt   - Passage time CDF");
        System.out.println("    • perct-respt - Response time percentiles (MAM solver)");
        System.out.println("    Transient:");
        System.out.println("    • tran-avg       - Transient average metrics");
        System.out.println("    • tran-cdf-respt - Transient response time CDF");
        System.out.println("    • tran-cdf-passt - Transient passage time CDF");
        System.out.println("    Probability:");
        System.out.println("    • prob         - State probability at node (requires -n)");
        System.out.println("    • prob-aggr    - Aggregated state probability (requires -n)");
        System.out.println("    • prob-marg    - Marginal state probability (requires -n, -c)");
        System.out.println("    • prob-sys     - System state probability");
        System.out.println("    • prob-sys-aggr - Aggregated system state probability");
        System.out.println("    Sampling (SSA solver):");
        System.out.println("    • sample         - Sample node state trajectory (requires -n)");
        System.out.println("    • sample-aggr    - Sample aggregated node state (requires -n)");
        System.out.println("    • sample-sys     - Sample system state trajectory");
        System.out.println("    • sample-sys-aggr - Sample aggregated system state");
        System.out.println("    Reward (CTMC solver):");
        System.out.println("    • reward       - Compute reward metrics");
        System.out.println("    • reward-steady - Steady-state reward");
        System.out.println("    • reward-value  - Reward value function (requires --reward-name)");
        System.out.println("    Default: all");
        System.out.println();
        
        System.out.println("-d, --seed <number>");
        System.out.println("    Random number seed for stochastic solvers (JMT, SSA).");
        System.out.println("    Use the same seed for reproducible results.");
        System.out.println("    Default: random number between 1 and 100000");
        System.out.println();
        
        System.out.println("-v, --verbosity <level>");
        System.out.println("    Verbosity level for solver output. Options:");
        System.out.println("    • normal - Standard output (default)");
        System.out.println("    • silent - Minimal output");
        System.out.println("    Default: normal");
        System.out.println();
        
        System.out.println("-n, --node <index>");
        System.out.println("    Node index for prob/sample analysis (0-based).");
        System.out.println("    Required for: prob, prob-aggr, prob-marg, sample, sample-aggr");
        System.out.println();

        System.out.println("-c, --class <index>");
        System.out.println("    Job class index for prob-marg analysis (0-based).");
        System.out.println("    Required for: prob-marg");
        System.out.println();

        System.out.println("--state <values>");
        System.out.println("    State vector for prob analysis (comma-separated integers).");
        System.out.println("    Example: --state 1,0,2");
        System.out.println();

        System.out.println("--events <number>");
        System.out.println("    Number of events for sample analysis.");
        System.out.println("    Default: 1000");
        System.out.println();

        System.out.println("--percentiles <values>");
        System.out.println("    Percentile values for perct-respt (comma-separated).");
        System.out.println("    Example: --percentiles 50,90,95,99");
        System.out.println("    Default: 50,90,95,99");
        System.out.println();

        System.out.println("--reward-name <name>");
        System.out.println("    Built-in reward function name for reward-value analysis.");
        System.out.println("    Valid names: QLen, Tput, Util, RespT, WaitT, ArvR, ResidT");
        System.out.println();

        System.out.println("-m, --maxreq <number>");
        System.out.println("    Maximum number of requests to process in server mode.");
        System.out.println("    Server will quit after processing this many requests.");
        System.out.println("    Currently not implemented.");
        System.out.println();

        System.out.println("-h, --help");
        System.out.println("    Display this detailed help information.");
        System.out.println();
        
        System.out.println("-V, --version");
        System.out.println("    Display version information.");
        System.out.println();
        
        System.out.println("SOLVER COMPATIBILITY:");
        System.out.println("┌─────────────────┬─────────────────────────────────────────────────┐");
        System.out.println("│ Input Format    │ Compatible Solvers                              │");
        System.out.println("├─────────────────┼─────────────────────────────────────────────────┤");
        System.out.println("│ jsim/jsimg/jsimw│ ctmc, fluid, jmt, mva, nc, ssa                  │");
        System.out.println("│ lqnx/xml        │ ln, lqns, mva (as ln.mva), nc (as ln.comom)     │");
        System.out.println("└─────────────────┴─────────────────────────────────────────────────┘");
        System.out.println();
        
        System.out.println("EXAMPLES:");
        System.out.println("  # Basic usage with file input");
        System.out.println("  java -jar jline.jar -f model.jsimg -s mva -a avg");
        System.out.println();
        System.out.println("  # Multiple analysis types (comma-separated)");
        System.out.println("  java -jar jline.jar -f model.jsimg -s mva -a avg,stage,chain");
        System.out.println();
        System.out.println("  # CDF analysis");
        System.out.println("  java -jar jline.jar -f model.jsimg -s jmt -a cdf-respt");
        System.out.println();
        System.out.println("  # Percentile analysis with custom values");
        System.out.println("  java -jar jline.jar -f model.jsimg -s mam -a perct-respt --percentiles 50,90,95,99");
        System.out.println();
        System.out.println("  # Sampling with SSA solver");
        System.out.println("  java -jar jline.jar -f model.jsimg -s ssa -a sample -n 1 --events 5000");
        System.out.println();
        System.out.println("  # Probability analysis at specific node");
        System.out.println("  java -jar jline.jar -f model.jsimg -s ctmc -a prob -n 1");
        System.out.println();
        System.out.println("  # Reward analysis");
        System.out.println("  java -jar jline.jar -f model.jsimg -s ctmc -a reward");
        System.out.println();
        System.out.println("  # Server mode");
        System.out.println("  java -jar jline.jar -p 8080");
        System.out.println();
        System.out.println("  # LQN model with layered network solver");
        System.out.println("  java -jar jline.jar -f model.lqnx -i lqnx -s ln -o readable");
        System.out.println();
        System.out.println("  # Docker usage");
        System.out.println("  cat model.jsimg | docker run -i --rm line-solver -i jsimg -s mva -a sys");
        System.out.println();
    }

    /**
     * Validates that a parameter value is provided and not empty.
     */
    private static boolean validateParameter(String paramName, String value) {
        if (value == null || value.trim().isEmpty()) {
            System.err.println("Error: Parameter " + paramName + " requires a value.");
            return false;
        }
        return true;
    }

    /**
     * Validates input format parameter.
     */
    private static boolean validateInputFormat(String format) {
        String[] validFormats = {"jsim", "jsimg", "jsimw", "lqnx", "xml"};
        for (String validFormat : validFormats) {
            if (validFormat.equals(format)) {
                return true;
            }
        }
        System.err.println("Error: Invalid input format '" + format + "'.");
        System.err.println("Valid formats: jsim, jsimg, jsimw, lqnx, xml");
        return false;
    }

    /**
     * Validates output format parameter.
     */
    private static boolean validateOutputFormat(String format) {
        String[] validFormats = {"json", "readable"};
        for (String validFormat : validFormats) {
            if (validFormat.equals(format)) {
                return true;
            }
        }
        System.err.println("Error: Invalid output format '" + format + "'.");
        System.err.println("Valid formats: json, readable");
        return false;
    }

    /**
     * Validates solver parameter.
     */
    private static boolean validateSolver(String solver) {
        String[] validSolvers = {"auto", "ctmc", "fluid", "jmt", "mam", "mva", "nc", "ssa", "ln", "lqns"};
        for (String validSolver : validSolvers) {
            if (validSolver.equals(solver)) {
                return true;
            }
        }
        System.err.println("Error: Invalid solver '" + solver + "'.");
        System.err.println("Valid solvers: auto, ctmc, fluid, jmt, mam, mva, nc, ssa, ln, lqns");
        return false;
    }

    /**
     * Selects an appropriate solver based on input format when 'auto' is specified.
     * @param inputFormat The input file format
     * @return The selected solver name
     */
    private static String autoSelectSolver(String inputFormat) {
        if (inputFormat.equals("lqnx") || inputFormat.equals("xml")) {
            return "ln";
        }
        // For JMT formats (jsim, jsimg, jsimw), use MVA as the default analytical solver
        return "mva";
    }

    /**
     * All valid analysis types.
     */
    private static final Set<String> VALID_ANALYSIS_TYPES = new HashSet<>(Arrays.asList(
        // Basic
        "all", "avg", "sys", "stage", "chain", "node", "nodechain",
        // Distribution
        "cdf-respt", "cdf-passt", "perct-respt",
        // Transient
        "tran-avg", "tran-cdf-respt", "tran-cdf-passt",
        // Probability
        "prob", "prob-aggr", "prob-marg", "prob-sys", "prob-sys-aggr",
        // Sampling
        "sample", "sample-aggr", "sample-sys", "sample-sys-aggr",
        // Reward
        "reward", "reward-steady", "reward-value"
    ));

    /**
     * Analysis types that require specific solvers.
     */
    private static final Map<String, Set<String>> ANALYSIS_SOLVER_COMPAT = new HashMap<>();
    static {
        ANALYSIS_SOLVER_COMPAT.put("sample", new HashSet<>(Arrays.asList("ssa")));
        ANALYSIS_SOLVER_COMPAT.put("sample-aggr", new HashSet<>(Arrays.asList("ssa")));
        ANALYSIS_SOLVER_COMPAT.put("sample-sys", new HashSet<>(Arrays.asList("ssa")));
        ANALYSIS_SOLVER_COMPAT.put("sample-sys-aggr", new HashSet<>(Arrays.asList("ssa")));
        ANALYSIS_SOLVER_COMPAT.put("reward", new HashSet<>(Arrays.asList("ctmc")));
        ANALYSIS_SOLVER_COMPAT.put("reward-steady", new HashSet<>(Arrays.asList("ctmc")));
        ANALYSIS_SOLVER_COMPAT.put("reward-value", new HashSet<>(Arrays.asList("ctmc")));
        ANALYSIS_SOLVER_COMPAT.put("perct-respt", new HashSet<>(Arrays.asList("mam")));
        ANALYSIS_SOLVER_COMPAT.put("prob", new HashSet<>(Arrays.asList("ctmc", "ssa")));
        ANALYSIS_SOLVER_COMPAT.put("prob-aggr", new HashSet<>(Arrays.asList("ctmc", "ssa")));
        ANALYSIS_SOLVER_COMPAT.put("prob-marg", new HashSet<>(Arrays.asList("ctmc", "ssa")));
        ANALYSIS_SOLVER_COMPAT.put("prob-sys", new HashSet<>(Arrays.asList("ctmc", "ssa")));
        ANALYSIS_SOLVER_COMPAT.put("prob-sys-aggr", new HashSet<>(Arrays.asList("ctmc", "ssa")));
    }

    /**
     * Analysis types that require node index.
     */
    private static final Set<String> ANALYSIS_REQUIRES_NODE = new HashSet<>(Arrays.asList(
        "prob", "prob-aggr", "prob-marg", "sample", "sample-aggr"
    ));

    /**
     * Analysis types that require class index.
     */
    private static final Set<String> ANALYSIS_REQUIRES_CLASS = new HashSet<>(Arrays.asList(
        "prob-marg"
    ));

    /**
     * Valid built-in reward names.
     */
    private static final Set<String> VALID_REWARD_NAMES = new HashSet<>(Arrays.asList(
        "QLen", "Tput", "Util", "RespT", "WaitT", "ArvR", "ResidT"
    ));

    /**
     * Validates analysis type parameter. Supports comma-separated values.
     */
    private static boolean validateAnalysis(String analysis) {
        String[] types = analysis.split(",");
        for (String type : types) {
            String trimmed = type.trim();
            if (!VALID_ANALYSIS_TYPES.contains(trimmed)) {
                System.err.println("Error: Invalid analysis type '" + trimmed + "'.");
                System.err.println("Valid types: " + String.join(", ", VALID_ANALYSIS_TYPES));
                return false;
            }
        }
        return true;
    }

    /**
     * Validates that analysis types are compatible with the chosen solver.
     */
    private static boolean validateAnalysisSolverCompat(String analysis, String solver) {
        // 'auto' solver will be resolved later - skip strict compatibility check
        // but warn if using analysis types that require specific solvers
        if (solver.equals("auto")) {
            String[] types = analysis.split(",");
            for (String type : types) {
                String trimmed = type.trim();
                Set<String> requiredSolvers = ANALYSIS_SOLVER_COMPAT.get(trimmed);
                if (requiredSolvers != null) {
                    System.err.println("Warning: Analysis type '" + trimmed + "' requires solver: " +
                        String.join(" or ", requiredSolvers) + ". Auto-selection may not choose a compatible solver.");
                }
            }
            return true;
        }
        String[] types = analysis.split(",");
        for (String type : types) {
            String trimmed = type.trim();
            Set<String> requiredSolvers = ANALYSIS_SOLVER_COMPAT.get(trimmed);
            if (requiredSolvers != null && !requiredSolvers.contains(solver)) {
                System.err.println("Error: Analysis type '" + trimmed + "' requires solver: " +
                    String.join(" or ", requiredSolvers) + ", but '" + solver + "' was specified.");
                return false;
            }
        }
        return true;
    }

    /**
     * Validates that required parameters are provided for analysis types.
     */
    private static boolean validateAnalysisParams(String analysis, Integer nodeIndex, Integer classIndex, String rewardName) {
        String[] types = analysis.split(",");
        for (String type : types) {
            String trimmed = type.trim();
            if (ANALYSIS_REQUIRES_NODE.contains(trimmed) && nodeIndex == null) {
                System.err.println("Error: Analysis type '" + trimmed + "' requires -n/--node parameter.");
                return false;
            }
            if (ANALYSIS_REQUIRES_CLASS.contains(trimmed) && classIndex == null) {
                System.err.println("Error: Analysis type '" + trimmed + "' requires -c/--class parameter.");
                return false;
            }
            if (trimmed.equals("reward-value") && (rewardName == null || rewardName.isEmpty())) {
                System.err.println("Error: Analysis type 'reward-value' requires --reward-name parameter.");
                return false;
            }
        }
        return true;
    }

    /**
     * Validates node index parameter.
     */
    private static boolean validateNodeIndex(String nodeStr) {
        try {
            int node = Integer.parseInt(nodeStr);
            if (node < 0) {
                System.err.println("Error: Node index must be non-negative.");
                return false;
            }
            return true;
        } catch (NumberFormatException e) {
            System.err.println("Error: Node index must be a valid integer.");
            return false;
        }
    }

    /**
     * Validates class index parameter.
     */
    private static boolean validateClassIndex(String classStr) {
        try {
            int classIdx = Integer.parseInt(classStr);
            if (classIdx < 0) {
                System.err.println("Error: Class index must be non-negative.");
                return false;
            }
            return true;
        } catch (NumberFormatException e) {
            System.err.println("Error: Class index must be a valid integer.");
            return false;
        }
    }

    /**
     * Validates events count parameter.
     */
    private static boolean validateEvents(String eventsStr) {
        try {
            int events = Integer.parseInt(eventsStr);
            if (events <= 0) {
                System.err.println("Error: Events count must be positive.");
                return false;
            }
            return true;
        } catch (NumberFormatException e) {
            System.err.println("Error: Events count must be a valid integer.");
            return false;
        }
    }

    /**
     * Validates percentiles parameter.
     */
    private static boolean validatePercentiles(String percentilesStr) {
        try {
            String[] parts = percentilesStr.split(",");
            for (String part : parts) {
                double p = Double.parseDouble(part.trim());
                if (p < 0 || p > 100) {
                    System.err.println("Error: Percentile values must be between 0 and 100.");
                    return false;
                }
            }
            return true;
        } catch (NumberFormatException e) {
            System.err.println("Error: Percentiles must be comma-separated numbers.");
            return false;
        }
    }

    /**
     * Validates reward name parameter.
     */
    private static boolean validateRewardName(String rewardName) {
        if (!VALID_REWARD_NAMES.contains(rewardName)) {
            System.err.println("Error: Invalid reward name '" + rewardName + "'.");
            System.err.println("Valid names: " + String.join(", ", VALID_REWARD_NAMES));
            return false;
        }
        return true;
    }

    /**
     * Parses state vector from comma-separated string.
     */
    private static Matrix parseState(String stateStr) {
        String[] parts = stateStr.split(",");
        Matrix state = new Matrix(1, parts.length);
        for (int i = 0; i < parts.length; i++) {
            state.set(0, i, Integer.parseInt(parts[i].trim()));
        }
        return state;
    }

    /**
     * Parses percentiles from comma-separated string.
     */
    private static double[] parsePercentiles(String percentilesStr) {
        String[] parts = percentilesStr.split(",");
        double[] percentiles = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            percentiles[i] = Double.parseDouble(parts[i].trim());
        }
        return percentiles;
    }

    /**
     * Validates verbosity level parameter.
     */
    private static boolean validateVerbosity(String verbosity) {
        String[] validLevels = {"normal", "silent"};
        for (String validLevel : validLevels) {
            if (validLevel.equals(verbosity)) {
                return true;
            }
        }
        System.err.println("Error: Invalid verbosity level '" + verbosity + "'.");
        System.err.println("Valid levels: normal, silent");
        return false;
    }

    /**
     * Validates port number parameter.
     */
    private static boolean validatePort(String portStr) {
        try {
            int port = Integer.parseInt(portStr);
            if (port < 1 || port > 65535) {
                System.err.println("Error: Port number must be between 1 and 65535.");
                return false;
            }
            return true;
        } catch (NumberFormatException e) {
            System.err.println("Error: Port must be a valid integer.");
            return false;
        }
    }

    /**
     * Validates seed parameter.
     */
    private static boolean validateSeed(String seedStr) {
        try {
            Integer.parseInt(seedStr);
            return true;
        } catch (NumberFormatException e) {
            System.err.println("Error: Seed must be a valid integer.");
            return false;
        }
    }

    /**
     * Validates solver compatibility with input format.
     */
    private static boolean validateSolverCompatibility(String inputFormat, String solver) {
        // 'auto' is always compatible - it will be resolved to an appropriate solver later
        if (solver.equals("auto")) {
            return true;
        }
        if (inputFormat.equals("lqnx") || inputFormat.equals("xml")) {
            String[] validLqnSolvers = {"ln", "lqns", "mva", "nc"};
            for (String validSolver : validLqnSolvers) {
                if (validSolver.equals(solver)) {
                    return true;
                }
            }
            System.err.println("Error: Solver '" + solver + "' is not compatible with input format '" + inputFormat + "'.");
            System.err.println("Valid solvers for LQN/XML formats: ln, lqns, mva, nc");
            return false;
        } else {
            String[] validJsimSolvers = {"ctmc", "fluid", "jmt", "mam", "mva", "nc", "ssa"};
            for (String validSolver : validJsimSolvers) {
                if (validSolver.equals(solver)) {
                    return true;
                }
            }
            if (solver.equals("ln") || solver.equals("lqns")) {
                System.err.println("Error: Solver '" + solver + "' is not compatible with input format '" + inputFormat + "'.");
                System.err.println("Valid solvers for JSIM formats: ctmc, fluid, jmt, mam, mva, nc, ssa");
                return false;
            }
        }
        return true;
    }

    // Temporary, only needed if return table in string type
    private static String consoleOutputToString(AvgTable avgTable) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(baos);
        PrintStream old = System.out;

        System.setOut(ps);
        avgTable.print();
        System.setOut(old);

        return baos.toString();
    }

    /**
     * Parses the command-line arguments provided to configure the LINE Solver.
     * This method processes various options such as input and output file formats,
     * solver selection, analysis type, server mode, and others.
     *
     * @param varargin an array of strings representing the command-line arguments.
     * @return a string, currently unused, but can be extended to return status or configuration details.
     * @throws IOException if an input or output exception occurs.
     */
    public static String parseArgs(String[] varargin) throws IOException {
        String ret = null;
        String inputext = "jsim";
        String solver = "mva";
        String analysis = "all";
        String outputext = "readable";
        String file = null;
        String verbosity = "normal";
        int randomSeed = 1 + RandomManager.nextInt((int) 1e5);
        boolean serverMode = false;
        int serverPort = 5863;

        // New parameters for extended analysis
        Integer nodeIndex = null;
        Integer classIndex = null;
        String stateStr = null;
        int numEvents = 1000;
        String percentilesStr = "50,90,95,99";
        String rewardName = null;

        // validate argument count - show help if no arguments
        if (varargin.length == 0) {
            printDetailedHelp();
            return null;
        }

        // initialise user options
        for (int v = 0; v < varargin.length; v += 2) {
            // Check if we have a value for parameters that require one
            if (v + 1 >= varargin.length && !varargin[v].equals("-h") && !varargin[v].equals("--help") &&
                !varargin[v].equals("-V") && !varargin[v].equals("--version")) {
                System.err.println("Error: Parameter " + varargin[v] + " requires a value.");
                System.err.println("Use -h or --help for usage information.");
                return null;
            }

            switch (varargin[v]) {
                case "-p":
                case "--port":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validatePort(varargin[v + 1])) {
                        return null;
                    }
                    serverMode = true;
                    serverPort = Integer.parseInt(varargin[v + 1]);
                    break;
                case "-m":
                case "--maxreq":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        return null;
                    }
                    // Maximum requests parameter is not currently implemented
                    line_warning("LineCLI", "--maxreq parameter is not currently implemented.");
                    break;
                case "-s":
                case "--solver":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateSolver(varargin[v + 1])) {
                        return null;
                    }
                    solver = varargin[v + 1];
                    break;
                case "-a":
                case "--analysis":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateAnalysis(varargin[v + 1])) {
                        return null;
                    }
                    analysis = varargin[v + 1];
                    break;
                case "-f":
                case "--file":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        return null;
                    }
                    file = varargin[v + 1];
                    break;
                case "-v":
                case "--verbosity":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateVerbosity(varargin[v + 1])) {
                        return null;
                    }
                    verbosity = varargin[v + 1];
                    break;
                case "-i":
                case "--input":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateInputFormat(varargin[v + 1])) {
                        return null;
                    }
                    inputext = varargin[v + 1];
                    break;
                case "-o":
                case "--output":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateOutputFormat(varargin[v + 1])) {
                        return null;
                    }
                    outputext = varargin[v + 1];
                    break;
                case "-d":
                case "--seed":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateSeed(varargin[v + 1])) {
                        return null;
                    }
                    randomSeed = Integer.parseInt(varargin[v + 1]);
                    break;
                case "-n":
                case "--node":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateNodeIndex(varargin[v + 1])) {
                        return null;
                    }
                    nodeIndex = Integer.parseInt(varargin[v + 1]);
                    break;
                case "-c":
                case "--class":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateClassIndex(varargin[v + 1])) {
                        return null;
                    }
                    classIndex = Integer.parseInt(varargin[v + 1]);
                    break;
                case "--state":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        return null;
                    }
                    stateStr = varargin[v + 1];
                    break;
                case "--events":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateEvents(varargin[v + 1])) {
                        return null;
                    }
                    numEvents = Integer.parseInt(varargin[v + 1]);
                    break;
                case "--percentiles":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validatePercentiles(varargin[v + 1])) {
                        return null;
                    }
                    percentilesStr = varargin[v + 1];
                    break;
                case "--reward-name":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateRewardName(varargin[v + 1])) {
                        return null;
                    }
                    rewardName = varargin[v + 1];
                    break;
                case "-h":
                case "--help":
                    printDetailedHelp();
                    return ret;
                case "-V":
                case "--version":
                    System.out.printf("%s%n", new Model("").getVersion());
                    return ret;
                default:
                    System.err.println("Error: Unknown parameter: " + varargin[v]);
                    System.err.println("Use -h or --help for usage information.");
                    return null;
            }
        }

        // Validate solver compatibility with input format
        if (!validateSolverCompatibility(inputext, solver)) {
            return null;
        }

        // Validate analysis-solver compatibility
        if (!validateAnalysisSolverCompat(analysis, solver)) {
            return null;
        }

        // Validate required parameters for analysis types
        if (!validateAnalysisParams(analysis, nodeIndex, classIndex, rewardName)) {
            return null;
        }

        Scanner scanner = new Scanner(System.in);

        if (serverMode) {
            LineWebSocketServer server = new LineWebSocketServer(serverPort);
            server.start(); // start LINE server mode

            while (true) {
                String cmd = scanner.nextLine();
                if (cmd.equalsIgnoreCase("q")) {
                    System.out.println("Shutting down. Please hold on, it may take several seconds.");
                    try {
                        server.stop();
                        return "Server closed.";
                    } catch (InterruptedException e) {
                        System.out.println("Unable to shut down the server. Please try again later.");
                    }
                }
            }
        }

        Path modelfile = Files.createTempFile("linetmp", "." + inputext);
        if (file == null) {
            String filecontent = scanner.nextLine();
            while (!filecontent.isEmpty()) {
                try {
                    Files.write(modelfile, filecontent.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
                    filecontent = scanner.nextLine();
                } catch (Exception e) {
                    break;
                }
            }
        } else {
            Files.copy(Paths.get(file), modelfile, StandardCopyOption.REPLACE_EXISTING);
        }
        String fileext = FilenameUtils.getExtension(modelfile.toString());
        String name = FilenameUtils.getBaseName(modelfile.toString());

        // initialise solver options
        SolverOptions solverOptions = new SolverOptions();
        solverOptions.seed(randomSeed);
        solverOptions.verbose(verbosity.equals("normal"));

        // Resolve 'auto' solver based on input format
        if (solver.equals("auto")) {
            solver = autoSelectSolver(inputext);
        }

        // choose solver
        Model model = null;
        Solver solverObj = null;
        switch (fileext) {
            case "jsimg":
            case "jsimw":
            case "jsim":
                model = new M2M().JSIM2LINE(modelfile.toString());
                switch (solver) {
                    case "ctmc":
                        solverOptions.force(true);
                        solverObj = new SolverCTMC((Network) model, solverOptions);
                        break;
                    case "fluid":
                        solverObj = new SolverFluid((Network) model, solverOptions);
                        break;
                    case "jmt":
                        solverObj = new SolverJMT((Network) model, solverOptions);
                        break;
                    case "mva":
                        solverObj = new SolverMVA((Network) model, solverOptions);
                        break;
                    case "mam":
                        solverObj = new SolverMAM((Network) model, solverOptions);
                        break;
                    case "nc":
                        solverObj = new SolverNC((Network) model, solverOptions);
                        break;
                    case "ssa":
                        solverObj = new SolverSSA((Network) model, solverOptions);
                        break;
                    default:
                        line_error(mfilename(new Object[]{}), "Unknown solver type: " + solver);
                }
                break;
            case "lqnx":
            case "xml":
                model = new M2M().LQN2LINE(modelfile.toString(), name);
                switch (solver) {
                    case "lqns":
                        solverObj = new SolverLQNS((LayeredNetwork) model, solverOptions);
                        break;
                    case "nc":
                    case "ln":
                    case "ln.comom":
                        solverObj = new SolverLN((LayeredNetwork) model,
                            (net) -> new SolverNC(net, solverOptions), solverOptions);
                        break;
                    case "mva":
                    case "ln.mva":
                        solverObj = new SolverLN((LayeredNetwork) model, solverOptions);
                        break;
                    default:
                        line_error(mfilename(new Object[]{}), "Unknown solver type: " + solver);
                }
                break;
        }

        // Parse analysis parameters
        Matrix stateMatrix = (stateStr != null) ? parseState(stateStr) : null;
        double[] percentiles = parsePercentiles(percentilesStr);

        // Execute multi-analysis
        Map<String, Object> analysisResults = new LinkedHashMap<>();
        String[] analysisTypes = analysis.split(",");

        for (String analysisType : analysisTypes) {
            String type = analysisType.trim();
            try {
                // Pass Network model only for NetworkSolver cases; LayeredNetwork doesn't cast to Network
                Network networkModel = (model instanceof Network) ? (Network) model : null;
                Object result = executeAnalysis(solverObj, networkModel, type,
                    nodeIndex, classIndex, stateMatrix, numEvents, percentiles, rewardName);
                if (result != null) {
                    analysisResults.put(type, result);
                }
            } catch (Exception e) {
                System.err.println("Error executing analysis '" + type + "': " + e.getMessage());
                if (verbosity.equals("normal")) {
                    e.printStackTrace();
                }
            }
        }

        // Format output
        switch (outputext) {
            case "json":
                ret = formatResultsAsJSON(analysisResults);
                break;
            case "readable":
                ret = formatResultsAsReadable(analysisResults);
                break;
        }

        return ret;
    }

    /**
     * Execute a single analysis type on the solver.
     */
    private static Object executeAnalysis(Solver solverObj, Network model, String analysisType,
            Integer nodeIndex, Integer classIndex, Matrix stateMatrix,
            int numEvents, double[] percentiles, String rewardName) throws Exception {

        switch (analysisType) {
            // Basic analysis types
            case "all":
                Map<String, Object> allResults = new LinkedHashMap<>();
                if (solverObj instanceof NetworkSolver) {
                    allResults.put("avg", ((NetworkSolver) solverObj).getAvgTable());
                    allResults.put("sys", ((NetworkSolver) solverObj).getAvgSysTable());
                } else if (solverObj instanceof SolverLN) {
                    allResults.put("avg", ((SolverLN) solverObj).getAvgTable());
                } else if (solverObj instanceof SolverLQNS) {
                    allResults.put("avg", ((SolverLQNS) solverObj).getAvgTable());
                }
                return allResults;

            case "avg":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgTable();
                } else if (solverObj instanceof SolverLN) {
                    return ((SolverLN) solverObj).getAvgTable();
                } else if (solverObj instanceof SolverLQNS) {
                    return ((SolverLQNS) solverObj).getAvgTable();
                }
                break;

            case "sys":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgSysTable();
                }
                break;

            case "stage":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getStageTable();
                }
                break;

            case "chain":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgChainTable();
                }
                break;

            case "node":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgNodeTable();
                }
                break;

            case "nodechain":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgNodeChainTable();
                }
                break;

            // Distribution analysis types
            case "cdf-respt":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getCdfRespT();
                }
                break;

            case "cdf-passt":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getCdfPassT();
                }
                break;

            case "perct-respt":
                if (solverObj instanceof SolverMAM) {
                    return ((SolverMAM) solverObj).getPerctRespT(percentiles);
                }
                break;

            // Transient analysis types
            case "tran-avg":
                if (solverObj instanceof NetworkSolver) {
                    ((NetworkSolver) solverObj).getTranAvg();
                    return "Transient analysis completed";
                }
                break;

            case "tran-cdf-respt":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getTranCdfRespT();
                }
                break;

            case "tran-cdf-passt":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getTranCdfPassT();
                }
                break;

            // Probability analysis types
            case "prob":
                if (solverObj instanceof NetworkSolver && nodeIndex != null) {
                    if (stateMatrix != null) {
                        return ((NetworkSolver) solverObj).getProb(nodeIndex, stateMatrix);
                    } else {
                        return ((NetworkSolver) solverObj).getProb(nodeIndex);
                    }
                }
                break;

            case "prob-aggr":
                if (solverObj instanceof NetworkSolver && nodeIndex != null) {
                    if (stateMatrix != null) {
                        return ((NetworkSolver) solverObj).getProbAggr(nodeIndex, stateMatrix);
                    } else {
                        return ((NetworkSolver) solverObj).getProbAggr(nodeIndex);
                    }
                }
                break;

            case "prob-marg":
                if (solverObj instanceof NetworkSolver && nodeIndex != null && classIndex != null) {
                    if (stateMatrix != null) {
                        return ((NetworkSolver) solverObj).getProbMarg(nodeIndex, classIndex, stateMatrix);
                    } else {
                        return ((NetworkSolver) solverObj).getProbMarg(nodeIndex, classIndex);
                    }
                }
                break;

            case "prob-sys":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getProbSys();
                }
                break;

            case "prob-sys-aggr":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getProbSysAggr();
                }
                break;

            // Sampling analysis types (SSA only)
            case "sample":
                if (solverObj instanceof SolverSSA && nodeIndex != null && model != null) {
                    Node node = model.getStatefulNodes().get(nodeIndex);
                    return ((SolverSSA) solverObj).sample(node, numEvents);
                }
                break;

            case "sample-aggr":
                if (solverObj instanceof SolverSSA && nodeIndex != null && model != null) {
                    Node node = model.getStatefulNodes().get(nodeIndex);
                    return ((SolverSSA) solverObj).sampleAggr(node, numEvents);
                }
                break;

            case "sample-sys":
                if (solverObj instanceof SolverSSA) {
                    return ((SolverSSA) solverObj).sampleSys(numEvents);
                }
                break;

            case "sample-sys-aggr":
                if (solverObj instanceof SolverSSA) {
                    return ((SolverSSA) solverObj).sampleSysAggr(numEvents);
                }
                break;

            // Reward analysis types (CTMC only)
            case "reward":
                if (solverObj instanceof SolverCTMC) {
                    return ((SolverCTMC) solverObj).getRewardResult();
                }
                break;

            case "reward-steady":
                if (solverObj instanceof SolverCTMC) {
                    return ((SolverCTMC) solverObj).getAvgReward();
                }
                break;

            case "reward-value":
                if (solverObj instanceof SolverCTMC && rewardName != null) {
                    return ((SolverCTMC) solverObj).getRewardValueFunction(rewardName);
                }
                break;
        }

        return null;
    }

    /**
     * Format analysis results as JSON.
     */
    private static String formatResultsAsJSON(Map<String, Object> results) {
        StringBuilder json = new StringBuilder();
        json.append("{");

        boolean first = true;
        for (Map.Entry<String, Object> entry : results.entrySet()) {
            if (!first) {
                json.append(",");
            }
            first = false;

            json.append("\"").append(entry.getKey()).append("\":");
            json.append(objectToJSON(entry.getValue()));
        }

        json.append("}");
        return json.toString();
    }

    /**
     * Convert an object to JSON representation.
     */
    private static String objectToJSON(Object obj) {
        if (obj == null) {
            return "null";
        }

        if (obj instanceof Map) {
            StringBuilder json = new StringBuilder();
            json.append("{");
            boolean first = true;
            for (Map.Entry<?, ?> entry : ((Map<?, ?>) obj).entrySet()) {
                if (!first) {
                    json.append(",");
                }
                first = false;
                json.append("\"").append(entry.getKey()).append("\":");
                json.append(objectToJSON(entry.getValue()));
            }
            json.append("}");
            return json.toString();
        }

        if (obj instanceof List) {
            StringBuilder json = new StringBuilder();
            json.append("[");
            boolean first = true;
            for (Object item : (List<?>) obj) {
                if (!first) {
                    json.append(",");
                }
                first = false;
                json.append(objectToJSON(item));
            }
            json.append("]");
            return json.toString();
        }

        if (obj instanceof AvgTable) {
            return avgTableToJSON((AvgTable) obj);
        }

        if (obj instanceof SampleNodeState) {
            return sampleNodeStateToJSON((SampleNodeState) obj);
        }

        if (obj instanceof SampleSysState) {
            return sampleSysStateToJSON((SampleSysState) obj);
        }

        if (obj instanceof Matrix) {
            return matrixToJSON((Matrix) obj);
        }

        if (obj instanceof double[]) {
            return doubleArrayToJSON((double[]) obj);
        }

        if (obj instanceof Number) {
            return obj.toString();
        }

        if (obj instanceof String) {
            return "\"" + escapeJSONString((String) obj) + "\"";
        }

        // Default: convert to string representation
        return "\"" + escapeJSONString(obj.toString()) + "\"";
    }

    /**
     * Convert SampleNodeState to JSON.
     */
    private static String sampleNodeStateToJSON(SampleNodeState result) {
        StringBuilder json = new StringBuilder();
        json.append("{\"type\":\"SampleNodeState\"");
        json.append(",\"isaggregate\":").append(result.isaggregate);
        if (result.t != null) {
            json.append(",\"t\":").append(matrixToJSON(result.t));
        }
        if (result.state != null) {
            json.append(",\"state\":").append(matrixToJSON(result.state));
        }
        json.append("}");
        return json.toString();
    }

    /**
     * Convert SampleSysState to JSON.
     */
    private static String sampleSysStateToJSON(SampleSysState result) {
        StringBuilder json = new StringBuilder();
        json.append("{\"type\":\"SampleSysState\"");
        json.append(",\"isaggregate\":").append(result.isaggregate);
        if (result.t != null) {
            json.append(",\"t\":").append(matrixToJSON(result.t));
        }
        if (result.state != null) {
            json.append(",\"state\":[");
            boolean first = true;
            for (Matrix state : result.state) {
                if (!first) json.append(",");
                first = false;
                json.append(matrixToJSON(state));
            }
            json.append("]");
        }
        json.append("}");
        return json.toString();
    }

    /**
     * Convert Matrix to JSON array.
     */
    private static String matrixToJSON(Matrix matrix) {
        if (matrix == null) {
            return "null";
        }
        StringBuilder json = new StringBuilder();
        int rows = matrix.getNumRows();
        int cols = matrix.getNumCols();

        if (rows == 1) {
            // Return as 1D array
            json.append("[");
            for (int j = 0; j < cols; j++) {
                if (j > 0) json.append(",");
                json.append(matrix.get(0, j));
            }
            json.append("]");
        } else {
            // Return as 2D array
            json.append("[");
            for (int i = 0; i < rows; i++) {
                if (i > 0) json.append(",");
                json.append("[");
                for (int j = 0; j < cols; j++) {
                    if (j > 0) json.append(",");
                    json.append(matrix.get(i, j));
                }
                json.append("]");
            }
            json.append("]");
        }
        return json.toString();
    }

    /**
     * Convert double array to JSON.
     */
    private static String doubleArrayToJSON(double[] arr) {
        StringBuilder json = new StringBuilder();
        json.append("[");
        for (int i = 0; i < arr.length; i++) {
            if (i > 0) json.append(",");
            json.append(arr[i]);
        }
        json.append("]");
        return json.toString();
    }

    /**
     * Format analysis results as readable text.
     */
    private static String formatResultsAsReadable(Map<String, Object> results) {
        StringBuilder sb = new StringBuilder();

        for (Map.Entry<String, Object> entry : results.entrySet()) {
            if (sb.length() > 0) {
                sb.append(System.lineSeparator()).append(System.lineSeparator());
            }

            String type = entry.getKey();
            Object value = entry.getValue();

            sb.append("=== ").append(type.toUpperCase()).append(" ===");
            sb.append(System.lineSeparator());

            if (value instanceof Map) {
                for (Map.Entry<?, ?> subEntry : ((Map<?, ?>) value).entrySet()) {
                    sb.append("--- ").append(subEntry.getKey()).append(" ---");
                    sb.append(System.lineSeparator());
                    sb.append(objectToReadable(subEntry.getValue()));
                    sb.append(System.lineSeparator());
                }
            } else {
                sb.append(objectToReadable(value));
            }
        }

        return sb.toString();
    }

    /**
     * Convert object to readable string.
     */
    private static String objectToReadable(Object obj) {
        if (obj == null) {
            return "(no result)";
        }

        if (obj instanceof AvgTable) {
            return consoleOutputToString((AvgTable) obj);
        }

        if (obj instanceof Matrix) {
            return ((Matrix) obj).toString();
        }

        if (obj instanceof SampleNodeState) {
            SampleNodeState state = (SampleNodeState) obj;
            StringBuilder sb = new StringBuilder();
            sb.append("Sample Node State (aggregate: ").append(state.isaggregate).append(")");
            sb.append(System.lineSeparator());
            if (state.t != null) {
                sb.append("Time points: ").append(state.t.length());
            }
            return sb.toString();
        }

        if (obj instanceof SampleSysState) {
            SampleSysState state = (SampleSysState) obj;
            StringBuilder sb = new StringBuilder();
            sb.append("Sample System State (aggregate: ").append(state.isaggregate).append(")");
            sb.append(System.lineSeparator());
            if (state.t != null) {
                sb.append("Time points: ").append(state.t.length());
            }
            return sb.toString();
        }

        if (obj instanceof RewardResult) {
            RewardResult result = (RewardResult) obj;
            StringBuilder sb = new StringBuilder();
            sb.append("Reward Result");
            if (result.getRewardNames() != null) {
                sb.append(System.lineSeparator());
                sb.append("Rewards: ").append(String.join(", ", result.getRewardNames()));
            }
            return sb.toString();
        }

        return obj.toString();
    }

    /**
     * Convert AvgTable objects to JSON format
     * Matches the MATLAB implementation using jsonencode(output)
     * @param avgTable Average performance metrics table
     * @param avgSysTable Average system performance metrics table  
     * @return JSON string representation
     */
    private static String convertToJSON(AvgTable avgTable, AvgTable avgSysTable) {
        StringBuilder json = new StringBuilder();
        json.append("[");
        
        boolean hasContent = false;
        
        if (avgTable != null) {
            json.append(avgTableToJSON(avgTable));
            hasContent = true;
        }
        
        if (avgSysTable != null) {
            if (hasContent) {
                json.append(",");
            }
            json.append(avgTableToJSON(avgSysTable));
            hasContent = true;
        }
        
        json.append("]");
        return json.toString();
    }
    
    /**
     * Convert a single AvgTable to JSON object
     * @param table AvgTable to convert
     * @return JSON object string
     */
    private static String avgTableToJSON(AvgTable table) {
        if (table == null) {
            return "null";
        }
        
        StringBuilder json = new StringBuilder();
        json.append("{");
        
        // Add table metadata
        json.append("\"type\":\"AvgTable\",");
        
        // Convert table data - this is a simplified conversion
        // In a production system, you might want to use a proper JSON library
        try {
            String tableString = consoleOutputToString(table);
            // Escape the string for JSON
            String escapedString = escapeJSONString(tableString);
            json.append("\"data\":\"").append(escapedString).append("\"");
        } catch (Exception e) {
            json.append("\"data\":\"Error converting table: ").append(escapeJSONString(e.getMessage())).append("\"");
        }
        
        json.append("}");
        return json.toString();
    }
    
    /**
     * Escape a string for JSON format
     * @param str String to escape
     * @return Escaped string
     */
    private static String escapeJSONString(String str) {
        if (str == null) {
            return "";
        }
        
        return str.replace("\\", "\\\\")
                 .replace("\"", "\\\"")
                 .replace("\n", "\\n")
                 .replace("\r", "\\r")
                 .replace("\t", "\\t");
    }
    
    /**
     * Main entry point for the LINE CLI
     * @param args Command line arguments
     */
    public static void main(String[] args) {
        try {
            String result = parseArgs(args);
            if (result != null) {
                System.out.println(result);
            }
        } catch (IOException e) {
            System.err.println("Error: " + e.getMessage());
            System.exit(1);
        } catch (Exception e) {
            System.err.println("Unexpected error: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }
    }
}
