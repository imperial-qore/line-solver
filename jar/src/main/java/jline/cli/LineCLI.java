/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.cli;

import static jline.io.InputOutput.line_warning;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.Environment;
import jline.lang.Model;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.layered.LayeredNetwork;
import jline.lang.nodes.Node;
import jline.lang.constant.SolverType;
import jline.solvers.AvgTable;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkAvgNodeTable;
import jline.solvers.NetworkAvgChainTable;
import jline.solvers.NetworkAvgNodeChainTable;
import jline.solvers.NetworkAvgSysTable;
import jline.solvers.NetworkAvgCacheTable;
import jline.solvers.NetworkAvgItemTable;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.lib.fjcodes.MainFJ.FJPercentileResult;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.ag.SolverAG;
import jline.solvers.ba.SolverBA;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.env.ENV;
import jline.solvers.env.SolverENV;
import jline.solvers.ctmc.analyzers.RewardResult;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.auto.SolverAUTO;
import jline.solvers.ln.SolverFactory;
import jline.solvers.ln.SolverLN;
import jline.solvers.wrappers.lqns.SolverLQNS;
import jline.solvers.mva.SolverMVA;
import jline.solvers.mam.SolverMAM;
import jline.solvers.nc.SolverNC;
import jline.solvers.wrappers.qns.SolverQNS;
import jline.solvers.uq.SolverUQ;
import jline.solvers.NetworkAvgOrbitTable;
import jline.solvers.NetworkLossTable;
import jline.solvers.NetworkSensitivityTable;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.ssa.SampleNodeState;
import jline.solvers.ssa.SampleSysState;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.LDESResult;
import jline.io.LineModelIO;
import jline.io.PnmlIO;
import jline.io.LDESResultIO;
import jline.io.M2M;
import jline.util.matrix.Matrix;
import jline.solvers.QrfParams;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.apache.commons.io.FilenameUtils;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.Scanner;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * The LineCLI class provides a command-line interface for configuring and running the LINE Solver.
 * It supports various options for input and output formats, solvers, analysis types, and more.
 * This class includes a static method to parse command-line arguments and set the appropriate
 * configuration options.
 */
public class LineCLI {
    /** Set by parseArgs when it rejects the command line, as opposed to serving help. */
    private static boolean usageError = false;

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
        System.out.println("Type model.help() for the solvers that support a model, "
                + "solver.libraries() for third-party dependencies, solver.citations() for references.");
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
        System.out.println("      jsim   - JSIM format (default)");
        System.out.println("      jsimg  - JSIM graphics format");
        System.out.println("      jsimw  - JSIM workspace format");
        System.out.println("      lqnx   - LQN XML format");
        System.out.println("      xml    - XML format");
        System.out.println("      json   - LINE portable model (Network, LayeredNetwork or Environment)");
        System.out.println("      pnml   - PNML place/transition net (ISO/IEC 15909-2)");
        System.out.println("    Default: jsim");
        System.out.println();
        
        System.out.println("-o, --output <format>");
        System.out.println("    Output format for results. Supported formats:");
        System.out.println("      readable - Human-readable text format (default)");
        System.out.println("      json     - JSON format");
        System.out.println("    Default: readable");
        System.out.println();
        
        System.out.println("-s, --solver <solver>");
        System.out.println("    Solver algorithm to use. Available solvers:");
        System.out.println("      auto   - Automatic solver selection");
        System.out.println("      mva    - Mean Value Analysis (default)");
        System.out.println("      ag     - Agent-based RCAT/INAP analyzer");
        System.out.println("      ba     - Bound analysis");
        System.out.println("      ctmc   - Continuous-Time Markov Chain");
        System.out.println("      fld    - Fluid/Mean-Field ODE solver (alias: fluid)");
        System.out.println("      jmt    - Java Modelling Tools simulation");
        System.out.println("      mam    - Matrix Analytic Methods");
        System.out.println("      nc     - Normalizing Constant analyzer");
        System.out.println("      ssa    - Stochastic Simulation Algorithm");
        System.out.println("      ldes   - LDES discrete-event simulator (alias: des)");
        System.out.println("      qns    - External qnsolver binary (alias: qnsolver)");
        System.out.println("      uq     - Uncertainty quantification over a Prior (needs --uq-solver)");
        System.out.println("      env    - Random environment (Environment json models only)");
        System.out.println("      ln     - Layered Network solver (MVA layers)");
        System.out.println("      ln.mva / ln.nc / ln.comom - Layered solver, layer engine named");
        System.out.println("      lqns   - LQN solver (LQN models only)");
        System.out.println("    Default: mva");
        System.out.println();
        
        System.out.println("-a, --analysis <type>");
        System.out.println("    Analysis type(s) to perform. Comma-separated for multiple.");
        System.out.println("    Basic:");
        System.out.println("      all         - Both avg and sys metrics (default)");
        System.out.println("      avg         - Average performance metrics");
        System.out.println("      sys         - System-level metrics");
        System.out.println("      stage       - Stage-based metrics (multi-stage service)");
        System.out.println("      chain       - Chain-level averages");
        System.out.println("      node        - Node-level averages");
        System.out.println("      nodechain   - Node-chain level averages");
        System.out.println("      cache       - Cache hit/miss metrics");
        System.out.println("      item        - Per-item cache metrics");
        System.out.println("      orbit       - Retrial orbit metrics");
        System.out.println("      loss        - Class-level loss metrics");
        System.out.println("      region-loss - Finite-capacity region loss metrics");
        System.out.println("      deadline    - Deadline-miss metrics (EDD/EDF)");
        System.out.println("      normconst   - Log normalizing constant (nc, mva)");
        System.out.println("      busyperiod  - Subnetwork busy period (nc, ldes; see --busyperiod*)");
        System.out.println("      sens        - Sensitivity to the service demands");
        System.out.println("    Distribution:");
        System.out.println("      cdf-respt   - Response time CDF");
        System.out.println("      cdf-passt   - Passage time CDF");
        System.out.println("      perct-respt - Response time percentiles (MAM solver)");
        System.out.println("    Transient:");
        System.out.println("      tran-avg       - Transient average metrics");
        System.out.println("      tran-cdf-respt - Transient response time CDF");
        System.out.println("      tran-cdf-passt - Transient passage time CDF");
        System.out.println("    Probability:");
        System.out.println("      prob         - State probability at node (requires -n)");
        System.out.println("      prob-aggr    - Aggregated state probability (requires -n)");
        System.out.println("      prob-marg    - Marginal state probability (requires -n, -c)");
        System.out.println("      prob-sys     - System state probability");
        System.out.println("      prob-sys-aggr - Aggregated system state probability");
        System.out.println("      prob-sys-marg - System marginal probability (requires --state)");
        System.out.println("    Sampling (SSA solver):");
        System.out.println("      sample         - Sample node state trajectory (requires -n)");
        System.out.println("      sample-aggr    - Sample aggregated node state (requires -n)");
        System.out.println("      sample-sys     - Sample system state trajectory");
        System.out.println("      sample-sys-aggr - Sample aggregated system state");
        System.out.println("    Reward (CTMC solver):");
        System.out.println("      reward       - Compute reward metrics");
        System.out.println("      reward-steady - Steady-state reward");
        System.out.println("      reward-value  - Reward value function (requires --reward-name)");
        System.out.println("    Solver-internal:");
        System.out.println("      generator   - CTMC infinitesimal generator and state space");
        System.out.println("      statevec    - Fluid ODE state vector");
        System.out.println("      moments     - Second-order moment-closure report (fld)");
        System.out.println("      interval    - Design-point envelope (uq)");
        System.out.println("    Default: all");
        System.out.println();
        
        System.out.println("-d, --seed <number>");
        System.out.println("    Random number seed for stochastic solvers (JMT, SSA).");
        System.out.println("    Use the same seed for reproducible results.");
        System.out.println("    Default: random number between 1 and 100000");
        System.out.println();
        
        System.out.println("--timespan <t0,t1>   (alias --tspan)");
        System.out.println("    Time span of the transient analyses. Without it getTranAvg");
        System.out.println("    substitutes 30/minRate, which is a different question.");
        System.out.println();

        System.out.println("--uq-solver <solver>");
        System.out.println("    Engine SolverUQ runs at each design point. Required by -s uq:");
        System.out.println("    UQ computes nothing itself, it expands the model's Prior.");
        System.out.println();

        System.out.println("--busyperiod <n[,n...]>  --busyperiod-subnet <i[,i...]>");
        System.out.println("    Orders and 0-based station indexes of -a busyperiod. The");
        System.out.println("    subnetwork has no default; the orders default to 1.");
        System.out.println();

        System.out.println("--sens-method <name>  --sens-scheme <name>  --sens-step <h>");
        System.out.println("    Differentiation of -a sens: analytic where the solver has one,");
        System.out.println("    else the named finite-difference scheme at step h.");
        System.out.println();

        System.out.println("--qrf-params <json>  --qrf-alpha <json>  --level <n>");
        System.out.println("    Parameterisation of the -s ba QRF reduction bounds. Inline JSON");
        System.out.println("    or a path; queue indexes stay 1-based, as in the MATLAB option.");
        System.out.println();

        System.out.println("-v, --verbosity <level>");
        System.out.println("    Verbosity level for solver output. Options:");
        System.out.println("      normal - Standard output (default; alias: standard)");
        System.out.println("      silent - Minimal output");
        System.out.println("      debug  - Turns on the solver console: a running");
        System.out.println("               progress log of every solver run");
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
        System.out.println("    State for a prob analysis: one node's per-class counts as");
        System.out.println("    comma-separated integers, or, for prob-sys/prob-sys-aggr, one");
        System.out.println("    such row per station separated by ';'.");
        System.out.println("    Example: --state 1,0,2     --state 0,0;1,2");
        System.out.println();

        System.out.println("--events <number>");
        System.out.println("    Number of events for sample analysis.");
        System.out.println("    Default: 1000");
        System.out.println();

        System.out.println("--samples <number>");
        System.out.println("    Number of simulation samples / Monte Carlo draws (SSA, NC).");
        System.out.println("    Default: 10000");
        System.out.println();

        System.out.println("--warmupfrac <fraction>");
        System.out.println("    SSA warmup discard: fraction of the trajectory dropped before");
        System.out.println("    computing steady-state averages and CI batch means (0 = disabled).");
        System.out.println("    Default: 0");
        System.out.println();

        System.out.println("--cutoff <number|matrix>");
        System.out.println("    State-space cutoff for CTMC/SSA on open/mixed models: one number,");
        System.out.println("    or a per-(station,class) matrix as ';'-separated rows of");
        System.out.println("    ','-separated cells.");
        System.out.println("    Example: --cutoff 1     --cutoff 1,1,0;3,3,0;0,0,3");
        System.out.println();

        System.out.println("--timespan <T0,T1>");
        System.out.println("    Time span of a transient analysis (tran-avg, tran-cdf-*).");
        System.out.println("    Example: --timespan 0,4");
        System.out.println();

        System.out.println("--timestep <number>");
        System.out.println("    Fixed output step of a transient analysis (default: adaptive).");
        System.out.println("    Example: --timestep 0.1");
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

        System.out.println("--method <name>");
        System.out.println("    Solution method for the selected solver (e.g. exact, amva, comom).");
        System.out.println("    Default: solver-specific default method.");
        System.out.println();

        System.out.println("--tol <number>");
        System.out.println("    General numerical tolerance. Default: solver-specific.");
        System.out.println();

        System.out.println("--iter_tol <number>");
        System.out.println("    Iteration convergence tolerance. Default: solver-specific.");
        System.out.println();

        System.out.println("--iter_max <number>");
        System.out.println("--multiserver <rule>");
        System.out.println("    Maximum number of iterations. Default: solver-specific.");
        System.out.println();

        System.out.println("--map-env <auto|off>");
        System.out.println("    Random-environment fallback for non-renewal (MAP/MMPP) processes a");
        System.out.println("    solver cannot consume natively. 'off' rejects such a model instead.");
        System.out.println("    Default: auto.");
        System.out.println();

        System.out.println("--map-env-method <auto|meanfield|dec|avg>");
        System.out.println("    Recombination of the environment stages. Default: auto.");
        System.out.println();

        System.out.println("--pstar <p[,p...]>");
        System.out.println("    Fluid p-norm smoothing exponent, one value or one per station.");
        System.out.println("    Default: unset, i.e. the unsmoothed min() drift.");
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
        System.out.println("  jsim/jsimg/jsimw  ctmc, fld, jmt, ldes, mam, mva, nc, ssa");
        System.out.println("  lqnx/xml          ln, ln.mva, ln.nc, ln.comom, lqns, mva, nc");
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
        String[] validFormats = {"jsim", "jsimg", "jsimw", "lqnx", "xml", "json", "pnml"};
        for (String validFormat : validFormats) {
            if (validFormat.equals(format)) {
                return true;
            }
        }
        System.err.println("Error: Invalid input format '" + format + "'.");
        System.err.println("Valid formats: jsim, jsimg, jsimw, lqnx, xml, json, pnml");
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
        // `ag` and `ba` ARE ENGINES OF THIS JAR and were missing from this list
        // alone: SolverAG (RCAT/INAP) and SolverBA (bounds) both exist, both have
        // a SolverType, and `solverTypeOf` already named BA -- so the only thing
        // between them and a caller was an argument error. The C++ line-cli
        // carries both, so a cross-codebase row that asked either of them of this
        // CLI reported "solver AG missing from output" rather than a refusal.
        // `qns` and `uq` ARE ENGINES OF THIS JAR and were missing from this list
        // alone, the same defect `ag` and `ba` had: SolverQNS (the qnsolver
        // wrapper) and SolverUQ (the prior-expanding design-of-models wrapper)
        // both exist and `solverTypeOf` already named QNS. The root `line-cli.py`
        // advertised `qns` and forwarded it here verbatim, so every such call
        // died on an argument error rather than running the wrapper.
        String[] validSolvers = {"ag", "auto", "ba", "ctmc", "env", "fld", "jmt", "ldes", "mam", "mva", "nc", "qns", "ssa", "uq", "ln", "ln.mva", "ln.nc", "ln.comom", "lqns"};
        String canonical = canonicalizeSolver(solver);
        for (String validSolver : validSolvers) {
            if (validSolver.equals(canonical)) {
                return true;
            }
        }
        System.err.println("Error: Invalid solver '" + solver + "'.");
        System.err.println("Valid solvers: ag, auto, ba, ctmc, env, fld, jmt, ldes, mam, mva, nc, qns, ssa, uq, ln, ln.mva, ln.nc, ln.comom, lqns");
        return false;
    }

    /**
     * Resolves a solver alias to its canonical method name. `fluid` is an alias of
     * `fld`, `qnsolver` of `qns` and `des` of `ldes` -- the last two are the
     * spellings the root `line-cli.py` wrapper offers, and accepting them here
     * is what makes that wrapper's own vocabulary reach a solver.
     */
    private static String canonicalizeSolver(String solver) {
        if ("fluid".equals(solver)) {
            return "fld";
        }
        if ("qnsolver".equals(solver)) {
            return "qns";
        }
        if ("des".equals(solver)) {
            return "ldes";
        }
        return solver;
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
        // Cache
        "cache", "item",
        // Station-class tables the reference publishes beside the AvgTable and
        // that had no token here, so a caller could reach them only in process:
        // the retrial orbit, the two loss views (class-level and finite-capacity
        // region) and the deadline-miss table of the EDD/EDF disciplines.
        "orbit", "loss", "region-loss", "deadline",
        // Distribution
        "cdf-respt", "cdf-passt", "perct-respt",
        // Transient
        "tran-avg", "tran-cdf-respt", "tran-cdf-passt",
        // Probability
        "prob", "prob-aggr", "prob-marg", "prob-sys", "prob-sys-aggr", "prob-sys-marg",
        // Normalizing constant, the quantity SolverNC exists to compute
        "normconst",
        // Busy periods of a subnetwork (Daduna 1988), the SolverNC/SolverLDES report
        "busyperiod",
        // Sampling
        "sample", "sample-aggr", "sample-sys", "sample-sys-aggr",
        // Reward
        "reward", "reward-steady", "reward-value",
        // Sensitivity of the mean metrics to the service demands
        "sens",
        // The SolverUQ views beside its prior-weighted -a avg
        "interval",
        // Solver-internal structures
        "generator", "statevec", "moments"
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
        // fld: SolverFluid answers getProbAggr from the moment closure's joint
        // normal, which only it holds -- a delegating caller that cannot ask for
        // it here falls back to the first-order binomial and reports an
        // approximation under the name of the closure.
        ANALYSIS_SOLVER_COMPAT.put("prob-aggr", new HashSet<>(Arrays.asList("ctmc", "ssa", "fld")));
        ANALYSIS_SOLVER_COMPAT.put("prob-marg", new HashSet<>(Arrays.asList("ctmc", "ssa", "mva", "nc", "mam")));
        ANALYSIS_SOLVER_COMPAT.put("prob-sys", new HashSet<>(Arrays.asList("ctmc", "ssa")));
        ANALYSIS_SOLVER_COMPAT.put("prob-sys-aggr", new HashSet<>(Arrays.asList("ctmc", "ssa")));
        // The joint marginal over the whole system. It is the same transform
        // family as prob-marg and served by the same engines, so it takes that
        // list rather than the two-simulator one prob-sys carries.
        ANALYSIS_SOLVER_COMPAT.put("prob-sys-marg", new HashSet<>(Arrays.asList("ctmc", "ssa", "mva", "nc", "mam")));
        // log G(N) is the normalizing constant of the product form, so it is
        // defined for the two solvers that compute one.
        ANALYSIS_SOLVER_COMPAT.put("normconst", new HashSet<>(Arrays.asList("nc", "mva")));
        // Daduna's subnetwork busy period is a normalizing-constant transform in
        // SolverNC and a sample-path measurement in SolverLDES.
        ANALYSIS_SOLVER_COMPAT.put("busyperiod", new HashSet<>(Arrays.asList("nc", "ldes")));
        ANALYSIS_SOLVER_COMPAT.put("generator", new HashSet<>(Arrays.asList("ctmc")));
        ANALYSIS_SOLVER_COMPAT.put("statevec", new HashSet<>(Arrays.asList("fld")));
        ANALYSIS_SOLVER_COMPAT.put("moments", new HashSet<>(Arrays.asList("fld")));
        // The interval report is SolverUQ's own: it is the support of the
        // design, which no single-model solver has.
        ANALYSIS_SOLVER_COMPAT.put("interval", new HashSet<>(Arrays.asList("uq")));
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
    /**
     * A comma-separated list of non-negative integers, as `--busyperiod` and
     * `--busyperiod-subnet` take it.
     *
     * @param spec the argument text
     * @return the parsed values, or null when any token is not an integer or
     *         the list is empty -- the caller turns that into a usage error
     */
    private static int[] parseIntList(String spec) {
        if (spec == null) {
            return null;
        }
        String[] parts = spec.split(",");
        List<Integer> vals = new ArrayList<Integer>();
        for (String part : parts) {
            String t = part.trim();
            if (t.isEmpty()) {
                continue;
            }
            try {
                vals.add(Integer.valueOf(Integer.parseInt(t)));
            } catch (NumberFormatException e) {
                return null;
            }
        }
        if (vals.isEmpty()) {
            return null;
        }
        int[] out = new int[vals.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = vals.get(i).intValue();
        }
        return out;
    }

    /**
     * Reads the argument of a JSON-valued flag: the text itself when it parses
     * as JSON, otherwise the contents of the file it names. Both spellings are
     * accepted by the C++ CLI's `--qrf-params`, so both are accepted here.
     *
     * @param spec inline JSON or a path to a file holding it
     * @return the JSON text, or null when a named file cannot be read
     */
    private static String readJsonArg(String spec) {
        String trimmed = spec.trim();
        if (trimmed.startsWith("{") || trimmed.startsWith("[")) {
            return trimmed;
        }
        try {
            return new String(Files.readAllBytes(Paths.get(trimmed)), StandardCharsets.UTF_8);
        } catch (IOException e) {
            System.err.println("Error: cannot read '" + trimmed + "': " + e.getMessage());
            return null;
        }
    }

    /**
     * The QRF blocking document of `-s ba --qrf-params`, in the field layout
     * {@link jline.solvers.QrfParams} and MATLAB's {@code options.config.qrf_params}
     * both use. THE QUEUE INDEXES IN f, MM AND MM1 STAY 1-BASED on the wire, as
     * they are in the MATLAB option, and the analyzer shifts them itself.
     *
     * @param spec inline JSON or a path to it
     * @return the parsed parameters, or null on a malformed document
     */
    private static QrfParams parseQrfParams(String spec) {
        String text = readJsonArg(spec);
        if (text == null) {
            return null;
        }
        try {
            JsonObject j = new JsonParser().parse(text).getAsJsonObject();
            String[] required = {"f", "MR", "BB", "MM", "MM1", "ZZ"};
            for (String key : required) {
                if (!j.has(key)) {
                    System.err.println("Error: --qrf-params is missing the required field '"
                        + key + "'.");
                    return null;
                }
            }
            QrfParams qp = new QrfParams();
            qp.f = j.get("f").getAsInt();
            qp.MR = j.get("MR").getAsInt();
            qp.BB = jsonTable(j.get("BB"));
            qp.MM = jsonTable(j.get("MM"));
            qp.MM1 = jsonTable(j.get("MM1"));
            qp.ZZ = jsonIntVector(j.get("ZZ"));
            // ZM is DERIVED from ZZ; a supplied one is ignored, as it is in the
            // C++ CLI and in the MATLAB option.
            int zm = 0;
            for (int i = 0; i < qp.ZZ.length; i++) {
                zm = Math.max(zm, qp.ZZ[i]);
            }
            qp.ZM = zm;
            if (j.has("F")) {
                qp.F = jsonIntVector(j.get("F"));
            }
            return qp;
        } catch (RuntimeException e) {
            System.err.println("Error: --qrf-params is not a valid QRF document: " + e.getMessage());
            return null;
        }
    }

    /** A JSON array of arrays as a Matrix; a flat array is read as one row. */
    private static Matrix jsonTable(JsonElement el) {
        JsonArray outer = el.getAsJsonArray();
        if (outer.size() > 0 && outer.get(0).isJsonArray()) {
            int rows = outer.size();
            int cols = outer.get(0).getAsJsonArray().size();
            Matrix m = new Matrix(rows, cols);
            for (int i = 0; i < rows; i++) {
                JsonArray row = outer.get(i).getAsJsonArray();
                for (int jx = 0; jx < row.size(); jx++) {
                    m.set(i, jx, row.get(jx).getAsDouble());
                }
            }
            return m;
        }
        Matrix m = new Matrix(1, outer.size());
        for (int i = 0; i < outer.size(); i++) {
            m.set(0, i, outer.get(i).getAsDouble());
        }
        return m;
    }

    /** A JSON array of integers. */
    private static int[] jsonIntVector(JsonElement el) {
        JsonArray arr = el.getAsJsonArray();
        int[] out = new int[arr.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = arr.get(i).getAsInt();
        }
        return out;
    }

    /** A JSON array of arrays as a double[][]; a flat array is read as one row. */
    private static double[][] parseDoubleTable(String spec) {
        String text = readJsonArg(spec);
        if (text == null) {
            return null;
        }
        try {
            JsonArray outer = new JsonParser().parse(text).getAsJsonArray();
            if (outer.size() > 0 && outer.get(0).isJsonArray()) {
                double[][] out = new double[outer.size()][];
                for (int i = 0; i < outer.size(); i++) {
                    JsonArray row = outer.get(i).getAsJsonArray();
                    out[i] = new double[row.size()];
                    for (int jx = 0; jx < row.size(); jx++) {
                        out[i][jx] = row.get(jx).getAsDouble();
                    }
                }
                return out;
            }
            double[][] out = new double[1][outer.size()];
            for (int i = 0; i < outer.size(); i++) {
                out[0][i] = outer.get(i).getAsDouble();
            }
            return out;
        } catch (RuntimeException e) {
            return null;
        }
    }

    private static Matrix parseState(String stateStr) {
        // ONE ROW PER STATEFUL UNIT. A node-level query (-a prob, prob-aggr,
        // prob-marg) names one row of per-class counts and is written without
        // ';'; a SYSTEM query (-a prob-sys, prob-sys-aggr) names one row per
        // station, station-major, and the rows are ';'-separated. Both spell a
        // cell the same way, so the single-row form is unchanged.
        String[] rows = stateStr.split(";");
        String[][] cells = new String[rows.length][];
        int ncols = -1;
        for (int i = 0; i < rows.length; i++) {
            cells[i] = rows[i].split(",");
            if (ncols < 0) {
                ncols = cells[i].length;
            } else if (cells[i].length != ncols) {
                throw new IllegalArgumentException(
                        "--state rows must all name the same number of classes");
            }
        }
        Matrix state = new Matrix(rows.length, ncols);
        for (int i = 0; i < rows.length; i++) {
            for (int j = 0; j < ncols; j++) {
                state.set(i, j, Integer.parseInt(cells[i][j].trim()));
            }
        }
        return state;
    }

    /**
     * Parses a per-(station,class) cutoff matrix: ';'-separated rows of
     * ','-separated cells, e.g. {@code 1,1,0;3,3,0;0,0,3}. Returns null when
     * the rows are ragged or a cell is not a number, so the caller can report
     * the spelling rather than solve a matrix it guessed at.
     */
    private static Matrix parseCutoffMatrix(String cutoffStr) {
        String[] rows = cutoffStr.split(";");
        String[][] cells = new String[rows.length][];
        int ncols = -1;
        for (int i = 0; i < rows.length; i++) {
            cells[i] = rows[i].split(",");
            if (ncols < 0) {
                ncols = cells[i].length;
            } else if (cells[i].length != ncols) {
                return null;
            }
        }
        if (ncols <= 0) {
            return null;
        }
        Matrix out = new Matrix(rows.length, ncols);
        for (int i = 0; i < rows.length; i++) {
            for (int j = 0; j < ncols; j++) {
                try {
                    out.set(i, j, Double.parseDouble(cells[i][j].trim()));
                } catch (NumberFormatException e) {
                    return null;
                }
            }
        }
        return out;
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
        // `standard` is what the C++ line-cli's help calls this level, and a
        // script written against that help used to die here on an argument
        // error. Accepted as a synonym of `normal` rather than renamed, so both
        // CLIs answer to both spellings and neither vocabulary is the odd one.
        String[] validLevels = {"normal", "standard", "silent", "debug", "verbose"};
        for (String validLevel : validLevels) {
            if (validLevel.equals(verbosity)) {
                return true;
            }
        }
        System.err.println("Error: Invalid verbosity level '" + verbosity + "'.");
        System.err.println("Valid levels: normal (alias standard), silent, debug");
        return false;
    }

    /** Canonical spelling of a verbosity level, folding the C++ CLI's synonyms. */
    private static String canonicalizeVerbosity(String verbosity) {
        if ("standard".equals(verbosity)) {
            return "normal";
        }
        if ("verbose".equals(verbosity)) {
            return "debug";
        }
        return verbosity;
    }

    /**
     * The {@link VerboseLevel} named by the {@code -v} token.
     *
     * <p>{@code debug} is what switches the SOLVER CONSOLE on, the running
     * progress log of {@link jline.io.LineConsole}: the console has no flag of
     * its own, it IS this level.</p>
     */
    private static VerboseLevel verboseLevelOf(String verbosity) {
        if ("silent".equals(verbosity)) {
            return VerboseLevel.SILENT;
        }
        if ("debug".equals(verbosity)) {
            return VerboseLevel.DEBUG;
        }
        return VerboseLevel.STD;
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
     * Maps a CLI solver method name onto the {@link SolverType} whose per-solver
     * defaults the options object must carry.
     *
     * The options object has to be built with this type, not with the no-arg
     * {@code new SolverOptions()}: the no-arg constructor supplies only the
     * generic defaults, so every solver-specific default is silently lost on
     * the CLI path (and therefore on the Python-native lang="java" dispatch,
     * which drives the CLI). For CTMC that meant losing
     * {@code config.hide_immediate=true}, {@code config.state_space_gen="full"}
     * and {@code cutoff=10}, so a model with immediate transitions (any
     * Fork/Join) was solved without eliminating them by stochastic
     * complementation: on the closed fork-join model of test_fj_driver_nc the
     * queue lengths came out 1.5e-8 away from the exact answer that MATLAB,
     * native Python and the in-process JLINE dispatch of the same JAR all agree
     * on to 1e-15.
     *
     * A method name with no matching type maps to null, which SolverOptions treats as
     * "generic defaults only" -- the previous behaviour, kept as the fallback.
     *
     * @param solver solver method name as given to -s
     * @return the matching SolverType, or null when the method name has none
     */
    private static SolverType solverTypeOf(String solver) {
        if (solver == null) {
            return null;
        }
        // Layered method names (ln, ln.mva, ln.nc, ln.comom) select the layer solver whose
        // options the LN solver forwards, so they take that solver's defaults.
        // Bare 'ln' is MVA layers, as in buildLayeredSolver; keep the two in step
        if (solver.equals("ln.comom") || solver.equals("ln.nc")) {
            return SolverType.NC;
        }
        if (solver.equals("ln") || solver.equals("ln.mva")) {
            return SolverType.MVA;
        }
        if (solver.equals("mva")) {
            return SolverType.MVA;
        }
        if (solver.equals("nc")) {
            return SolverType.NC;
        }
        if (solver.equals("ctmc")) {
            return SolverType.CTMC;
        }
        if (solver.equals("fld")) {
            return SolverType.FLUID;
        }
        if (solver.equals("mam")) {
            return SolverType.MAM;
        }
        if (solver.equals("ssa")) {
            return SolverType.SSA;
        }
        if (solver.equals("jmt")) {
            return SolverType.JMT;
        }
        if (solver.equals("ldes")) {
            return SolverType.LDES;
        }
        if (solver.equals("lqns")) {
            return SolverType.LQNS;
        }
        if (solver.equals("qns")) {
            return SolverType.QNS;
        }
        if (solver.equals("env")) {
            return SolverType.ENV;
        }
        if (solver.equals("ba")) {
            return SolverType.BA;
        }
        if (solver.equals("ag")) {
            return SolverType.AG;
        }
        if (solver.equals("auto")) {
            return SolverType.AUTO;
        }
        return null;
    }

    /**
     * Builds an analytical/simulation solver for a Network model, shared by the
     * JSIM and line-model JSON input branches.
     *
     * @param solver        solver token (mva, nc, ctmc, fld, mam, ssa, jmt, ldes, qns, ag, ba, uq, auto)
     * @param model         the Network model
     * @param solverOptions common solver options (seed, verbosity)
     * @param randomSeed    seed forwarded to the LDES options
     * @param verbosity     verbosity string ("normal"/"silent")
     * @param uqSolver      engine SolverUQ runs at each design point, from
     *                      `--uq-solver`; read only by the `uq` method name
     * @return the constructed solver, or null if the method name is unknown
     */
    private static Solver buildNetworkSolver(String solver, Network model,
            SolverOptions solverOptions, int randomSeed, String verbosity,
            final String uqSolver) {
        switch (solver) {
            // RCAT/INAP. The options carry `--method` already, which is what
            // picks between inap, inapplus, inapinf and exact; `default` resolves
            // to inap inside the analyzer, as it does in every codebase.
            case "ag":
                return new SolverAG(model, solverOptions);
            // Bound analysis. Closed-form, so it converges nothing and reads
            // none of the iteration knobs the ladder above validates.
            case "ba":
                return new SolverBA(model, solverOptions);
            case "ctmc":
                solverOptions.force(true);
                return new SolverCTMC(model, solverOptions);
            case "fld":
                return new SolverFluid(model, solverOptions);
            case "jmt":
                return new SolverJMT(model, solverOptions);
            case "mva":
                return new SolverMVA(model, solverOptions);
            case "mam":
                return new SolverMAM(model, solverOptions);
            case "nc":
                return new SolverNC(model, solverOptions);
            case "ssa":
                return new SolverSSA(model, solverOptions);
            case "ldes":
                LDESOptions ldesOpts = new LDESOptions();
                ldesOpts.seed(randomSeed);
                ldesOpts.verbose(verboseLevelOf(verbosity));
                return new SolverLDES(model, ldesOpts);
            // The external `qnsolver` binary. A wrapper, so the numeric knobs the
            // ladder validates describe a child process it runs rather than an
            // algorithm it implements; the options still carry seed and verbosity.
            case "qns":
                return new SolverQNS(model, solverOptions);
            // SolverUQ COMPUTES NOTHING ITSELF: it expands the model's Prior into
            // a design and runs another solver at each point, so it needs one
            // named. Defaulting the inner engine would attribute the numbers to a
            // solver the caller never chose, which is why `--uq-solver` is
            // required rather than assumed -- the same contract as the C++ CLI.
            case "uq":
                if (uqSolver == null || uqSolver.isEmpty()) {
                    line_error(mfilename(new Object(){}),
                        "-s uq needs --uq-solver: UQ expands the Prior and runs another solver "
                        + "at each design point (the CLI spelling of UQ(model, @SolverMVA))");
                    return null;
                }
                final SolverOptions uqInnerOptions = solverOptions;
                final int uqSeed = randomSeed;
                final String uqVerbosity = verbosity;
                return new SolverUQ(model, new SolverUQ.SolverFactory() {
                    @Override
                    public NetworkSolver create(Network net) {
                        Solver inner = buildNetworkSolver(uqSolver, net, uqInnerOptions,
                                uqSeed, uqVerbosity, null);
                        if (!(inner instanceof NetworkSolver)) {
                            line_error(mfilename(new Object(){}),
                                "--uq-solver '" + uqSolver + "' is not a Network solver");
                            return null;
                        }
                        return (NetworkSolver) inner;
                    }
                }, solverOptions);
            case "auto":
                return new SolverAUTO(model, solverOptions);
            default:
                // A KNOWN TOKEN FOR A DIFFERENT MODEL CLASS IS NOT AN UNKNOWN
                // TOKEN. `ln`, `ln.*`, `lqns` and `env` all reach this arm when
                // the document turned out to hold a flat Network, and reporting
                // "Unknown solver type: ln" named the wrong problem: the method name
                // is perfectly good and the MODEL is not the one it solves,
                // which is what the caller has to fix.
                if (solver.equals("ln") || solver.startsWith("ln.") || solver.equals("lqns")) {
                    line_error(mfilename(new Object(){}),
                        "Solver '" + solver + "' solves a LayeredNetwork and this document holds "
                        + "a Network; use a flat-network solver (mva, nc, ctmc, fld, mam, ssa, "
                        + "ldes, jmt, qns, ag, ba, uq) or supply an .lqnx / layered .json model.");
                    return null;
                }
                if (solver.equals("env")) {
                    line_error(mfilename(new Object(){}),
                        "Solver 'env' solves an Environment and this document holds a Network; "
                        + "an Environment model.json carries the stage networks and the "
                        + "transition process, which a flat Network does not.");
                    return null;
                }
                line_error(mfilename(new Object(){}), "Unknown solver type: " + solver);
                return null;
        }
    }

    /**
     * Builds a layered-network solver for a LayeredNetwork model, shared by the
     * LQN (lqnx/xml) input branch and the line-model JSON input branch (used by
     * the Python-native lang="java" dispatch of SolverLN/SolverLQNS).
     *
     * @param solver        solver token (ln, ln.mva, ln.nc, ln.comom, lqns, mva, nc, auto)
     * @param model         the LayeredNetwork model
     * @param solverOptions LAYER solver options (the defaults of the layer solver the
     *                      token names, plus the caller's overrides)
     * @param lnOptions     SolverLN's own options (LN defaults -- relaxation, iter_max,
     *                      iter_tol, interlocking -- plus the same caller overrides).
     *                      The fixed point reads these and NOT the layer object: see the
     *                      note where the two are built
     * @return the constructed solver, or null if the method name is unknown
     */
    private static Solver buildLayeredSolver(String solver, LayeredNetwork model,
            SolverOptions solverOptions, SolverOptions lnOptions) {
        switch (solver) {
            case "lqns":
                return new SolverLQNS(model, solverOptions);
            case "nc":
            case "ln.comom":
            case "ln.nc":
            case "auto":
                // NC layers, reached only by name; 'auto' resolves here for LQN
                return new SolverLN(model, (net) -> new SolverNC(net, solverOptions), lnOptions);
            case "mva":
            case "ln":
            case "ln.mva":
                // see _kb/12-interfaces-and-docs.md: bare 'ln' is SolverLN's own default
                return new SolverLN(model, lnOptions);
            default:
                line_error(mfilename(new Object(){}), "Unknown solver type: " + solver);
                return null;
        }
    }

    /**
     * Builds a random-environment solver for an Environment model, reached from
     * the line-model JSON input branch.
     *
     * <p>THE STAGE SOLVER IS NAMED, NOT GUESSED. SolverENV couples its stages
     * through their TRANSIENT means, which both the fluid analyzer and the
     * enumerated CTMC produce -- and the two are different models rather than
     * two routes to one answer, since a chain holds whole jobs (SolverENV's own
     * `roundMarginalForDiscreteSolver` runs for exactly the non-fluid case).
     * `model.json` carries the stage Networks and the transition process and NOT
     * the ensemble's solver choice, so leaving the factory fixed answered a
     * different model in silence: on renv_threestages_repairmen, whose stages
     * are SolverCTMC, that returned Queue1 throughput 1.5577 against the 1.3333
     * the CTMC stages give. `--stage-solver` states it.
     *
     * @param model         the Environment model
     * @param solverOptions common solver options (seed, verbosity, method, tol)
     * @param stageSolver   `fluid` (the default) or `ctmc`; null means fluid
     * @return the constructed solver
     */
    private static Solver buildEnvSolver(Environment model, SolverOptions solverOptions,
            final String stageSolver) {
        // THE STAGE HORIZON IS PART OF THE COUPLING, so it cannot be left to
        // each stage's own guess. SolverFluid with no timespan picks one from
        // the model it was handed ("End time of transient analysis unspecified,
        // setting the timespan option to [0, 37.5]"), and on renv_node_breakdown
        // that horizon is short against an Exp(0.1) breakdown clock: the blend
        // then runs its full 100 iterations without converging and reports
        // Server QLen 44.95 at throughput 1.85, above the 0.8 the source
        // admits, which is not an approximation but an answer refuted by flow
        // balance. Given a horizon it converges in 52. The default matches
        // `EnvOptions::timespan_end` in the C++ engine, so both CLIs answer the
        // same question when the caller states nothing; `--timespan` overrides.
        if (solverOptions.timespan == null || solverOptions.timespan.length < 2
                || !Double.isFinite(solverOptions.timespan[1])) {
            solverOptions.timespan = new double[]{0, 100};
        }
        final boolean ctmcStages = "ctmc".equals(stageSolver);
        return new ENV(model, new SolverFactory() {
            @Override
            public NetworkSolver at(Network net) {
                if (ctmcStages) {
                    return new SolverCTMC(net, solverOptions);
                }
                return new SolverFluid(net, solverOptions);
            }
        }, solverOptions);
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
            String[] validLqnSolvers = {"ln", "ln.mva", "ln.nc", "ln.comom", "lqns", "mva", "nc"};
            for (String validSolver : validLqnSolvers) {
                if (validSolver.equals(solver)) {
                    return true;
                }
            }
            System.err.println("Error: Solver '" + solver + "' is not compatible with input format '" + inputFormat + "'.");
            System.err.println("Valid solvers for LQN/XML formats: ln, ln.mva, ln.nc, ln.comom, lqns, mva, nc");
            return false;
        } else if (inputFormat.equals("json")) {
            // see _kb/12-interfaces-and-docs.md for LineCLI's model-type/dispatch rationale
            String[] validJsonSolvers = {"ag", "ba", "ctmc", "env", "fld", "jmt", "ldes", "mam", "mva",
                "nc", "qns", "ssa", "uq", "ln", "ln.mva", "ln.nc", "ln.comom", "lqns"};
            for (String validSolver : validJsonSolvers) {
                if (validSolver.equals(solver)) {
                    return true;
                }
            }
            System.err.println("Error: Solver '" + solver + "' is not compatible with input format 'json'.");
            System.err.println("Valid solvers for json: ag, ba, ctmc, env, fld, jmt, ldes, mam, mva, nc, qns, ssa, uq, ln, ln.mva, ln.nc, ln.comom, lqns");
            return false;
        } else if (inputFormat.equals("pnml")) {
            // A PNML document is a place/transition net, so only the solvers whose
            // feature set declares Transition can answer for it; the product-form
            // solvers cannot.
            String[] validPnmlSolvers = {"ctmc", "jmt", "ldes", "ssa"};
            for (String validSolver : validPnmlSolvers) {
                if (validSolver.equals(solver)) {
                    return true;
                }
            }
            System.err.println("Error: Solver '" + solver + "' is not compatible with input format 'pnml'.");
            System.err.println("Valid solvers for pnml: ctmc, jmt, ldes, ssa");
            return false;
        } else {
            // A JSIM document carries a flat Network, so every Network engine may
            // answer for it -- `ag`, `ba`, `qns` and `uq` included. They were
            // absent for no reason but this list: the json branch above already
            // admits them and both branches build the model through the same
            // `buildNetworkSolver`.
            String[] validJsimSolvers = {"ag", "ba", "ctmc", "fld", "jmt", "ldes", "mam", "mva", "nc", "qns", "ssa", "uq"};
            for (String validSolver : validJsimSolvers) {
                if (validSolver.equals(solver)) {
                    return true;
                }
            }
            if (solver.equals("ln") || solver.equals("lqns") || solver.equals("env")
                    || solver.startsWith("ln.")) {
                System.err.println("Error: Solver '" + solver + "' is not compatible with input format '" + inputFormat + "'.");
                System.err.println("Valid solvers for JSIM formats: ag, ba, ctmc, fld, jmt, ldes, mam, mva, nc, qns, ssa, uq");
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
        // see _kb/12-interfaces-and-docs.md (null-return ambiguity: usage error vs help/version)
        usageError = false;
        String ret = null;
        String inputext = "jsim";
        String solver = "mva";
        String analysis = "all";
        String outputext = "readable";
        String file = null;
        String verbosity = "normal";
        // LINE's own default stream: a CLI that draws a fresh seed per invocation
        // answers a different sample path every run, so a simulated solve is not
        // reproducible. The C++ line-cli and the MATLAB JMTIO default alike to 23000.
        int randomSeed = 23000;
        boolean serverMode = false;
        int serverPort = 5863;

        // New parameters for extended analysis
        Integer nodeIndex = null;
        Integer classIndex = null;
        String stateStr = null;
        int numEvents = 1000;
        String percentilesStr = "50,90,95,99";
        String rewardName = null;
        double cutoffVal = Double.NaN; // scalar state-space cutoff (NaN = not set)
        // --stage-solver: which solver runs each STAGE of an Environment. It is
        // not --layer-solver: a layer of an LQN is solved in steady state and a
        // stage of a random environment transiently, so their admissible sets
        // are different sets for different reasons. null = not given, i.e. the
        // coupling's own default (fluid).
        String stageSolver = null;
        Matrix cutoffMat = null; // per-(station,class) state-space cutoff (null = not set)
        double[] timespanVal = null; // transient analysis time span (null = not set)
        double timestepVal = Double.NaN; // fixed transient output step (NaN = adaptive)
        int samplesVal = -1; // simulation samples / Monte Carlo draws (-1 = not set)
        double warmupfracVal = Double.NaN; // SSA warmup discard fraction (NaN = not set)
        String multiserverVal = null; // AMVA multiserver rule (null = leave the default)
        String mapEnvVal = null;      // non-renewal random-environment gate (null = leave the default)
        String mapEnvMethodVal = null; // environment recombination (null = leave the default)
        String pstarVal = null; // fluid p-norm smoothing exponent(s) (null = unsmoothed min())
        String methodVal = null; // solution method (null = solver default)
        double tolVal = Double.NaN; // general tolerance (NaN = not set)
        double iterTolVal = Double.NaN; // iteration convergence tolerance (NaN = not set)
        int iterMaxVal = -1; // maximum iterations (-1 = not set)
        // SolverUQ computes nothing itself: it runs this engine at each design
        // point. Required by -s uq rather than defaulted, so the numbers are
        // never attributed to a solver the caller did not name.
        String uqSolverVal = null;
        String sensMethodVal = "auto";   // sensitivity differentiation (auto = analytic where available)
        String sensSchemeVal = "forward"; // finite-difference scheme
        double sensStepVal = Double.NaN;  // finite-difference step (NaN = the solver's own)
        // The subnetwork and the orders of -a busyperiod. The flag names are
        // LdesCLI's, so one spelling drives the transform and the simulation.
        int[] busyOrdersVal = null;
        int[] busySubnetVal = null;
        String qrfParamsVal = null; // -s ba QRF blocking document (inline JSON or a path)
        String qrfAlphaVal = null;  // -s ba QRF load-dependent scaling (inline JSON or a path)
        int levelVal = -1;          // -s ba hierarchy level / iteration count (-1 = not set)

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
                usageError = true;
                return null;
            }

            switch (varargin[v]) {
                case "-p":
                case "--port":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validatePort(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    serverMode = true;
                    serverPort = Integer.parseInt(varargin[v + 1]);
                    break;
                case "-m":
                case "--maxreq":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    // Maximum requests parameter is not currently implemented
                    line_warning("LineCLI", "--maxreq parameter is not currently implemented.");
                    break;
                case "-s":
                case "--solver":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateSolver(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    solver = canonicalizeSolver(varargin[v + 1]);
                    break;
                case "-a":
                case "--analysis":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateAnalysis(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    analysis = varargin[v + 1];
                    break;
                case "-f":
                case "--file":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    file = varargin[v + 1];
                    break;
                case "-v":
                case "--verbosity":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateVerbosity(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    verbosity = canonicalizeVerbosity(varargin[v + 1]);
                    // Switched on for the WHOLE process, not just for the
                    // solver's own options, so that the model compile narrates
                    // too: LineConsole.writes() falls back to the session level
                    // when no run is open, and reading the model happens before
                    // any solver exists.
                    GlobalConstants.setVerbose(verboseLevelOf(verbosity));
                    break;
                case "-i":
                case "--input":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateInputFormat(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    inputext = varargin[v + 1];
                    break;
                case "-o":
                case "--output":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateOutputFormat(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    outputext = varargin[v + 1];
                    break;
                case "-d":
                case "--seed":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateSeed(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    randomSeed = Integer.parseInt(varargin[v + 1]);
                    break;
                case "-n":
                case "--node":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateNodeIndex(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    nodeIndex = Integer.parseInt(varargin[v + 1]);
                    break;
                case "-c":
                case "--class":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateClassIndex(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    classIndex = Integer.parseInt(varargin[v + 1]);
                    break;
                case "--state":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    stateStr = varargin[v + 1];
                    break;
                case "--events":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateEvents(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    numEvents = Integer.parseInt(varargin[v + 1]);
                    break;
                case "--stage-solver":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    stageSolver = varargin[v + 1].trim().toLowerCase();
                    if (!stageSolver.equals("fluid") && !stageSolver.equals("fld")
                            && !stageSolver.equals("ctmc")) {
                        System.err.println("Error: --stage-solver '" + stageSolver
                                + "' is not available: the environment coupling needs a TRANSIENT "
                                + "stage solve, and only the fluid analyzer and the enumerated "
                                + "CTMC provide one.");
                        usageError = true;
                        return null;
                    }
                    break;
                case "--cutoff":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    // A CUTOFF IS NOT ALWAYS A SCALAR. options.cutoff is a
                    // (station x class) matrix in every codebase, and an
                    // example that pins one -- oqn_cs_routing writes
                    // [1,1,0;3,3,0;0,0,3] -- is asking for a DIFFERENT
                    // truncation per cell. Forwarded as nothing, the JAR solved
                    // its own default 10 and reported the answer under the
                    // caller's name; forwarded as one number it would solve a
                    // third model. Rows are separated by ';' and cells by ','.
                    if (varargin[v + 1].indexOf(';') >= 0 || varargin[v + 1].indexOf(',') >= 0) {
                        cutoffMat = parseCutoffMatrix(varargin[v + 1]);
                        if (cutoffMat == null) {
                            System.err.println("Error: --cutoff takes a number or a "
                                    + "';'-separated list of ','-separated rows.");
                            usageError = true;
                            return null;
                        }
                    } else {
                        cutoffVal = Double.parseDouble(varargin[v + 1]);
                    }
                    break;
                // `--tspan` is the C++ line-cli's spelling of this flag. Both
                // CLIs now take both, so a script moving between them does not
                // die on an argument error over one name for one horizon.
                case "--tspan":
                case "--timespan":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    // BOTH SEPARATORS, and a bare end time. This CLI has always
                    // written the horizon `t0,t1` and the C++ line-cli `t0:t1`;
                    // accepting only one made a command line naming a horizon
                    // unportable between them even after the flag names were
                    // reconciled. A single value is the END of the horizon,
                    // starting at 0, which is what a transient from the initial
                    // state means -- the C++ port's reading of a bare value.
                    String tsArg = varargin[v + 1].trim();
                    String[] tsParts = tsArg.contains(":") ? tsArg.split(":") : tsArg.split(",");
                    if (tsParts.length > 2 || tsParts.length < 1) {
                        System.err.println("Error: --timespan takes T1, T0,T1 or T0:T1.");
                        usageError = true;
                        return null;
                    }
                    try {
                        if (tsParts.length == 1) {
                            timespanVal = new double[]{0.0, Double.parseDouble(tsParts[0].trim())};
                        } else {
                            timespanVal = new double[]{Double.parseDouble(tsParts[0].trim()),
                                    Double.parseDouble(tsParts[1].trim())};
                        }
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid timespan values.");
                        usageError = true;
                        return null;
                    }
                    break;
                case "--timestep":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    timestepVal = Double.parseDouble(varargin[v + 1]);
                    break;
                case "--samples":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateEvents(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    samplesVal = Integer.parseInt(varargin[v + 1]);
                    break;
                case "--warmupfrac":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    try {
                        warmupfracVal = Double.parseDouble(varargin[v + 1]);
                    } catch (NumberFormatException e) {
                        usageError = true;
                        return null;
                    }
                    break;
                case "--multiserver":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    multiserverVal = varargin[v + 1];
                    break;
                case "--map-env":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    mapEnvVal = varargin[v + 1];
                    break;
                case "--map-env-method":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    mapEnvMethodVal = varargin[v + 1];
                    break;
                case "--pstar":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    pstarVal = varargin[v + 1];
                    break;
                case "--percentiles":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validatePercentiles(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    percentilesStr = varargin[v + 1];
                    break;
                case "--method":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    methodVal = varargin[v + 1];
                    break;
                case "--tol":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    try {
                        tolVal = Double.parseDouble(varargin[v + 1]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: --tol requires a numeric value.");
                        usageError = true;
                        return null;
                    }
                    break;
                case "--iter_tol":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    try {
                        iterTolVal = Double.parseDouble(varargin[v + 1]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: --iter_tol requires a numeric value.");
                        usageError = true;
                        return null;
                    }
                    break;
                case "--iter_max":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    try {
                        iterMaxVal = Integer.parseInt(varargin[v + 1]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: --iter_max requires an integer value.");
                        usageError = true;
                        return null;
                    }
                    break;
                case "--reward-name":
                    if (!validateParameter(varargin[v], varargin[v + 1]) || !validateRewardName(varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    rewardName = varargin[v + 1];
                    break;
                case "--uq-solver":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    uqSolverVal = canonicalizeSolver(varargin[v + 1]);
                    break;
                case "--sens-method":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    sensMethodVal = varargin[v + 1];
                    break;
                case "--sens-scheme":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    sensSchemeVal = varargin[v + 1];
                    break;
                case "--sens-step":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    try {
                        sensStepVal = Double.parseDouble(varargin[v + 1]);
                    } catch (NumberFormatException e) {
                        usageError = true;
                        return null;
                    }
                    break;
                case "--busyperiod":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    busyOrdersVal = parseIntList(varargin[v + 1]);
                    if (busyOrdersVal == null) {
                        System.err.println("Error: --busyperiod takes a comma-separated list of "
                            + "positive integer orders, e.g. --busyperiod 1,2,3.");
                        usageError = true;
                        return null;
                    }
                    break;
                case "--busyperiod-subnet":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    busySubnetVal = parseIntList(varargin[v + 1]);
                    if (busySubnetVal == null) {
                        System.err.println("Error: --busyperiod-subnet takes a comma-separated list "
                            + "of 0-based station indexes, e.g. --busyperiod-subnet 1,2.");
                        usageError = true;
                        return null;
                    }
                    break;
                case "--qrf-params":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    qrfParamsVal = varargin[v + 1];
                    break;
                case "--qrf-alpha":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    qrfAlphaVal = varargin[v + 1];
                    break;
                case "--level":
                    if (!validateParameter(varargin[v], varargin[v + 1])) {
                        usageError = true;
                        return null;
                    }
                    try {
                        levelVal = Integer.parseInt(varargin[v + 1]);
                    } catch (NumberFormatException e) {
                        usageError = true;
                        return null;
                    }
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
                    usageError = true;
                    return null;
            }
        }

        // Validate solver compatibility with input format
        if (!validateSolverCompatibility(inputext, solver)) {
            usageError = true;
            return null;
        }

        // Validate analysis-solver compatibility
        if (!validateAnalysisSolverCompat(analysis, solver)) {
            usageError = true;
            return null;
        }

        // Validate required parameters for analysis types
        if (!validateAnalysisParams(analysis, nodeIndex, classIndex, rewardName)) {
            usageError = true;
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

        // see _kb/12-interfaces-and-docs.md for LineCLI's per-solver-defaults + option-forwarding rationale
        SolverOptions solverOptions = new SolverOptions(solverTypeOf(solver));
        // A LAYERED method name builds TWO solvers -- SolverLN, and the layer solver whose
        // defaults solverTypeOf names -- and one options object cannot carry both
        // sets. The object handed to SolverLN used to carry the LAYER defaults, so
        // the fixed point ran with config.relax "none" (the generic default) instead
        // of the LN "fixed"/0.5, and with the layer solver's iter_max/iter_tol. On
        // model_C2_L2_T4_P2_1_doesnotconverge_v3 that is the difference between
        // converging in 53 iterations onto the MATLAB answer and oscillating for
        // 1000 iterations without converging, reporting c_task_1 QLen 59.8 against a
        // population of 22. The caller's own overrides are applied to BOTH objects,
        // so every flag still reaches whichever of the two solvers reads it.
        SolverOptions lnOptions = new SolverOptions(SolverType.LN);
        SolverOptions[] cliOptions = new SolverOptions[]{solverOptions, lnOptions};
        for (SolverOptions opts : cliOptions) {
            opts.seed(randomSeed);
            opts.verbose(verboseLevelOf(verbosity));
            // see _kb/12-interfaces-and-docs.md for LineCLI's per-solver-defaults + option-forwarding rationale
            if (cutoffMat != null) {
                opts.cutoff(cutoffMat);
            } else if (!Double.isNaN(cutoffVal)) {
                opts.cutoff(cutoffVal);
            }
            // without this a transient analysis falls back to the [Inf,Inf] default and
            // getTranAvg substitutes 30/minRate, which is not the caller's time span
            if (timespanVal != null) {
                opts.timespan = timespanVal;
            }
            // null leaves the output grid to the adaptive ODE solver, as MATLAB does
            if (!Double.isNaN(timestepVal) && timestepVal > 0) {
                opts.timestep = Double.valueOf(timestepVal);
            }
            // see _kb/12-interfaces-and-docs.md for LineCLI's per-solver-defaults + option-forwarding rationale
            if (samplesVal > 0) {
                opts.samples(samplesVal);
            }
            // SSA warmup discard fraction forwarded from the caller (mean-estimate
            // trajectory discard plus CI batch-means transient discard).
            if (!Double.isNaN(warmupfracVal) && warmupfracVal > 0) {
                opts.config.warmupfrac = warmupfracVal;
            }
            // The AMVA multiserver rule decides which algorithm serves a multiserver
            // model (Seidmann transform, solver_amvald, conway, krzesinski), so a solve
            // that does not receive it silently answers under the DEFAULT rule while
            // reporting the caller's. Forwarded like any other option.
            if (multiserverVal != null && !multiserverVal.isEmpty()) {
                opts.config.multiserver = multiserverVal;
            }
            // The non-renewal random-environment gate decides whether a MAP/MMPP model
            // this solver cannot consume natively is APPROXIMATED through its
            // environment image or REJECTED. A delegated solve that never receives it
            // answers where the caller asked for a refusal, which is the one outcome
            // map_env="off" exists to produce.
            if (mapEnvVal != null && !mapEnvVal.isEmpty()) {
                opts.config.map_env = mapEnvVal;
            }
            if (mapEnvMethodVal != null && !mapEnvMethodVal.isEmpty()) {
                opts.config.map_env_method = mapEnvMethodVal;
            }
            // Fluid p-norm smoothing exponent (Ruuskanen et al., PEVA 151 (2021), eq.
            // (26)). It selects the DRIFT the matrix method integrates, so a solve that
            // never receives it returns the unsmoothed mean-field fixed point (lambda
            // for an M/M/1) while the caller believes it asked for the smoothed one.
            if (pstarVal != null && !pstarVal.isEmpty()) {
                opts.config.pstar = new ArrayList<Double>();
                for (String p : pstarVal.split(",")) {
                    p = p.trim();
                    if (p.isEmpty()) {
                        continue;
                    }
                    try {
                        opts.config.pstar.add(Double.parseDouble(p));
                    } catch (NumberFormatException e) {
                        usageError = true;
                        return null;
                    }
                }
            }
            // see _kb/12-interfaces-and-docs.md for LineCLI's per-solver-defaults + option-forwarding rationale
            if (methodVal != null) {
                opts.method(methodVal);
            }
            if (!Double.isNaN(tolVal)) {
                opts.tol = tolVal;
            }
            if (!Double.isNaN(iterTolVal)) {
                opts.iter_tol = iterTolVal;
            }
            if (iterMaxVal > 0) {
                opts.iter_max = iterMaxVal;
            }
            // The QRF reduction bounds are the only SolverBA arms that take a
            // parameterisation, and `-s ba` reached this CLI (2026-08-27)
            // without a way to state one: `qrf.bas.mem` requires F, and a
            // caller who cannot supply the blocking document gets the derived
            // structure or a refusal, never the one it meant. The C++ line-cli
            // has carried all three flags since 2026-08-02 -- these are their
            // Java spelling, reading the same document key for key.
            if (qrfParamsVal != null) {
                opts.qrfParams = parseQrfParams(qrfParamsVal);
                if (opts.qrfParams == null) {
                    usageError = true;
                    return null;
                }
            }
            if (qrfAlphaVal != null) {
                opts.qrfAlpha = parseDoubleTable(qrfAlphaVal);
                if (opts.qrfAlpha == null) {
                    System.err.println("Error: --qrf-alpha takes a JSON array of arrays, "
                        + "inline or as a path.");
                    usageError = true;
                    return null;
                }
            }
            // The hierarchy level of pbh/cbh/sib and the iteration count of
            // pbk/bjbk. Without it a caller setting options.level silently got
            // the default 2.
            if (levelVal > 0) {
                opts.level = levelVal;
            }
        }

        // see _kb/12-interfaces-and-docs.md for LineCLI's 'auto' solver resolution rationale
        if (solver.equals("auto") && (inputext.equals("lqnx") || inputext.equals("xml"))) {
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
                solverObj = buildNetworkSolver(solver, (Network) model, solverOptions, randomSeed, verbosity, uqSolverVal);
                break;
            case "json":
                // see _kb/12-interfaces-and-docs.md for LineCLI's model-type/dispatch rationale
                Object loadedModel = LineModelIO.load(modelfile.toString());
                if (loadedModel instanceof LayeredNetwork) {
                    model = (Model) loadedModel;
                    solverObj = buildLayeredSolver(solver, (LayeredNetwork) model, solverOptions, lnOptions);
                } else if (loadedModel instanceof Network) {
                    model = (Model) loadedModel;
                    solverObj = buildNetworkSolver(solver, (Network) model, solverOptions, randomSeed, verbosity, uqSolverVal);
                } else if (loadedModel instanceof Environment) {
                    // An ENVIRONMENT is its own model class, not a Network with
                    // extra fields, and `-s env` is the only token that solves
                    // one: SolverENV couples the stages, so asking any Network
                    // solver for it would have to pick a stage and answer about
                    // a model the document does not describe.
                    model = (Model) loadedModel;
                    if (!solver.equals("env")) {
                        line_error(mfilename(new Object(){}),
                            "An Environment model is solved by '-s env'; '" + solver
                            + "' solves a Network, and this document holds the coupling "
                            + "rather than any one stage.");
                    }
                    solverObj = buildEnvSolver((Environment) model, solverOptions, stageSolver);
                } else {
                    line_error(mfilename(new Object(){}),
                        "JSON input supports Network, LayeredNetwork and Environment models only.");
                }
                break;
            case "lqnx":
            case "xml":
                model = new M2M().LQN2LINE(modelfile.toString(), name);
                solverObj = buildLayeredSolver(solver, (LayeredNetwork) model, solverOptions, lnOptions);
                break;
            case "pnml":
                model = PnmlIO.load(modelfile.toString());
                solverObj = buildNetworkSolver(solver, (Network) model, solverOptions, randomSeed, verbosity, uqSolverVal);
                break;
        }

        // Parse analysis parameters
        Matrix stateMatrix = (stateStr != null) ? parseState(stateStr) : null;
        double[] percentiles = parsePercentiles(percentilesStr);

        // Execute multi-analysis
        Map<String, Object> analysisResults = new LinkedHashMap<>();
        String[] analysisTypes = analysis.split(",");
        Exception firstError = null;
        String firstErrorType = null;

        for (String analysisType : analysisTypes) {
            String type = analysisType.trim();
            try {
                // Pass Network model only for NetworkSolver cases; LayeredNetwork doesn't cast to Network
                Network networkModel = (model instanceof Network) ? (Network) model : null;
                Object result = executeAnalysis(solverObj, networkModel, type,
                    nodeIndex, classIndex, stateMatrix, numEvents, percentiles, rewardName,
                    sensMethodVal, sensSchemeVal, sensStepVal, busySubnetVal, busyOrdersVal);
                if (result != null) {
                    analysisResults.put(type, result);
                }
            } catch (Exception e) {
                System.err.println("Error executing analysis '" + type + "': " + e.getMessage());
                if (!verbosity.equals("silent")) {
                    e.printStackTrace();
                }
                if (firstError == null) {
                    firstError = e;
                    firstErrorType = type;
                }
            }
        }

        // see _kb/12-interfaces-and-docs.md (an all-failed run must not report as an empty success)
        if (analysisResults.isEmpty() && firstError != null) {
            throw new RuntimeException("analysis '" + firstErrorType + "' failed: "
                    + firstError.getMessage(), firstError);
        }

        // see _kb/12-interfaces-and-docs.md (SolverLN convergence flag exposed for lang="java" callers)
        if (solverObj instanceof SolverLN && !analysisResults.isEmpty()) {
            analysisResults.put("hasconverged", ((SolverLN) solverObj).hasconverged);
        }

        // see _kb/12-interfaces-and-docs.md (post-dispatch method name, iteration
        // count and convergence flag exposed for lang="java" callers). The method
        // name is what a delegating caller reads to tell WHICH algorithm answered:
        // SolverNC alone dispatches an FCR loss network to lossn.exact, erlangfp
        // or mci from the same 'default' request, and without this key the caller
        // cannot distinguish an exact transform from an approximation. Without the
        // count and flag it cannot raise the non-convergence warning its native
        // path raises, so the same model reports "possibly not converged" under
        // lang='python' and stays silent under lang='java' -- the silence being
        // exactly the defect the warning exists for.
        if (solverObj instanceof NetworkSolver && !analysisResults.isEmpty()) {
            jline.solvers.SolverResult netRes = ((NetworkSolver) solverObj).result;
            if (netRes != null) {
                if (netRes.method != null) {
                    analysisResults.put("method", netRes.method);
                }
                // The log normalizing constant is a first-class NC result, not a
                // diagnostic: it is what a caller compares against a reference
                // g(N), and reading it off a delegated solve is otherwise
                // impossible.
                if (netRes instanceof jline.solvers.nc.NCResult) {
                    double lg = ((jline.solvers.nc.NCResult) netRes).logNormConstAggr();
                    if (!Double.isNaN(lg)) {
                        analysisResults.put("lG", Double.valueOf(lg));
                    }
                }
                analysisResults.put("iter", Integer.valueOf(netRes.iter));
                if (netRes instanceof jline.solvers.mva.MVAResult) {
                    Boolean conv = ((jline.solvers.mva.MVAResult) netRes).converged;
                    if (conv != null) {
                        analysisResults.put("converged", conv);
                    }
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
            int numEvents, double[] percentiles, String rewardName,
            String sensMethod, String sensScheme, double sensStep,
            int[] busySubnet, int[] busyOrders) throws Exception {

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
                } else if (solverObj instanceof SolverENV) {
                    allResults.put("avg", ((SolverENV) solverObj).getAvgTable());
                }
                return allResults;

            case "avg":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgTable();
                } else if (solverObj instanceof SolverLN) {
                    return ((SolverLN) solverObj).getAvgTable();
                } else if (solverObj instanceof SolverLQNS) {
                    return ((SolverLQNS) solverObj).getAvgTable();
                } else if (solverObj instanceof SolverENV) {
                    return ((SolverENV) solverObj).getAvgTable();
                } else if (solverObj instanceof SolverUQ) {
                    // SolverUQ is an EnsembleSolver, not a NetworkSolver, so it
                    // needs its own arm: -a avg is the PRIOR-WEIGHTED
                    // expectation over the design, not one model's table.
                    return ((SolverUQ) solverObj).getAvgTable();
                }
                break;

            // The support-only range over the design points, SolverUQ's second
            // report beside the weighted mean.
            case "interval":
                if (solverObj instanceof SolverUQ) {
                    return intervalToMap(((SolverUQ) solverObj).getInterval());
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

            // Cache analysis types
            case "cache":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgCacheTable();
                }
                break;

            case "item":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgItemTable();
                }
                break;

            // The station-class tables the reference publishes beside the
            // AvgTable. Each was reachable in process and from the C++ CLI and
            // had no token here, so a delegating caller could not read the
            // retrial orbit, either loss view, or the deadline-miss counts at
            // all.
            case "orbit":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgOrbitTable();
                }
                break;

            case "loss":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgLossTable();
                }
                break;

            case "region-loss":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getAvgRegionLossTable();
                }
                break;

            case "deadline":
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getDeadlineTable();
                }
                break;

            // Sensitivity of the mean metrics to the service demands. `auto`
            // takes the analytic derivative where the solver publishes one and
            // falls back to the named finite-difference scheme otherwise.
            case "sens":
                if (solverObj instanceof NetworkSolver) {
                    NetworkSolver sensSolver = (NetworkSolver) solverObj;
                    if (Double.isNaN(sensStep)) {
                        return sensSolver.getSensitivityTable();
                    }
                    return sensSolver.getSensitivityTable(sensMethod, sensStep, sensScheme);
                }
                break;

            // log G(N). The analyzers park it in the package-private prob block
            // and leave the flat NCResult.lG at zero, so this is the only route
            // to the quantity SolverNC exists to compute.
            case "normconst":
                if (solverObj instanceof NetworkSolver) {
                    NetworkSolver ncSolver = (NetworkSolver) solverObj;
                    ncSolver.getAvg();
                    return ncSolver.getProbNormConstAggr();
                }
                break;

            // Mean busy period of a subnetwork, Daduna (J. ACM 35(3), 1988).
            // SolverNC evaluates the transform and SolverLDES measures it on the
            // sample path; both answer the one -a token, with the flag names
            // LdesCLI already used.
            case "busyperiod": {
                if (busySubnet == null || busySubnet.length == 0) {
                    throw new IllegalArgumentException(
                        "-a busyperiod needs --busyperiod-subnet: the busy period is defined "
                        + "for a named subnetwork, which no default can choose.");
                }
                int[] orders = (busyOrders == null || busyOrders.length == 0)
                        ? new int[]{1} : busyOrders;
                Map<String, Object> bp = new LinkedHashMap<String, Object>();
                // as Lists, not int[]: the JSON writer renders an unrecognised
                // object with String.valueOf, so a raw array reaches the caller
                // as "[I@a74868d" rather than as the subnetwork it names
                List<Integer> subnetOut = new ArrayList<Integer>();
                for (int t = 0; t < busySubnet.length; t++) {
                    subnetOut.add(Integer.valueOf(busySubnet[t]));
                }
                List<Integer> ordersOut = new ArrayList<Integer>();
                for (int t = 0; t < orders.length; t++) {
                    ordersOut.add(Integer.valueOf(orders[t]));
                }
                bp.put("subnet", subnetOut);
                bp.put("orders", ordersOut);
                double[] durations;
                if (solverObj instanceof SolverNC) {
                    durations = ((SolverNC) solverObj).getAvgBusyPeriod(busySubnet, orders);
                } else if (solverObj instanceof SolverLDES) {
                    SolverLDES ldesSolver = (SolverLDES) solverObj;
                    durations = new double[orders.length];
                    for (int t = 0; t < orders.length; t++) {
                        // the LDES report is per class; -1 is its all-classes total
                        durations[t] = ldesSolver.getAvgBusyPeriod(busySubnet, -1, orders[t]);
                    }
                } else {
                    break;
                }
                // a List, for the same reason subnet and orders are: the JSON
                // writer and the readable formatter both render an unrecognised
                // object with String.valueOf, so a raw double[] reaches the
                // caller as "[D@2a70a3d8" instead of as the durations
                List<Double> bOut = new ArrayList<Double>();
                for (int t = 0; t < durations.length; t++) {
                    bOut.add(Double.valueOf(durations[t]));
                }
                bp.put("b", bOut);
                return bp;
            }

            // Solver-internal structures (CTMC generator / state space, fluid state vector)
            case "generator":
                if (solverObj instanceof SolverCTMC) {
                    SolverCTMC ctmcSolver = (SolverCTMC) solverObj;
                    // see _kb/12-interfaces-and-docs.md for the CTMC chain/generator ordering rationale
                    ctmcSolver.getAvg();
                    jline.solvers.SolverResult cres = ctmcSolver.getResults();
                    jline.solvers.ctmc.CTMCResult ctmcRes =
                            (cres instanceof jline.solvers.ctmc.CTMCResult)
                                    ? (jline.solvers.ctmc.CTMCResult) cres : null;
                    Matrix infGen = ctmcRes != null ? ctmcRes.infGen : null;
                    Matrix space = ctmcRes != null ? ctmcRes.space : null;
                    jline.util.matrix.MatrixCell eventFiltCell =
                            ctmcRes != null ? ctmcRes.eventFilt : null;
                    if (infGen == null || space == null) {
                        // see _kb/12-interfaces-and-docs.md for the CTMC chain/generator ordering rationale
                        SolverCTMC.generatorResult gr = ctmcSolver.getGenerator();
                        SolverCTMC.StateSpace ss = ctmcSolver.getStateSpace();
                        infGen = gr.infGen;
                        space = ss.stateSpace;
                        eventFiltCell = gr.eventFilt;
                        cres = ctmcSolver.getResults();
                        ctmcRes = (cres instanceof jline.solvers.ctmc.CTMCResult)
                                ? (jline.solvers.ctmc.CTMCResult) cres : null;
                    }
                    java.util.Map<String, Object> gen = new LinkedHashMap<String, Object>();
                    gen.put("infgen", infGen);
                    gen.put("space", space);
                    gen.put("space_aggr", ctmcRes != null ? ctmcRes.spaceAggr : null);
                    gen.put("pi", ctmcRes != null ? ctmcRes.pi : null);
                    // The work triple: the generator and state spaces pi is indexed
                    // by. A consumer adopts all three together or none of them.
                    gen.put("infgen_work", ctmcRes != null ? ctmcRes.infGenWork : null);
                    gen.put("space_work", ctmcRes != null ? ctmcRes.spaceWork : null);
                    gen.put("space_aggr_work", ctmcRes != null ? ctmcRes.spaceAggrWork : null);
                    java.util.List<Object> eventFilt = new java.util.ArrayList<Object>();
                    if (eventFiltCell != null) {
                        for (int i = 0; i < eventFiltCell.size(); i++) {
                            eventFilt.add(eventFiltCell.get(i));
                        }
                    }
                    gen.put("eventFilt", eventFilt);
                    // Widths of the per-stateful-node blocks of a state row, in
                    // stateful order. A client has to cut `space` into node blocks
                    // to answer getStateSpace or any per-station probability, and
                    // the widths are the only thing that says where the cuts go: a
                    // 2-server FCFS station carries a buffer column AND a service
                    // column for one class, so guessing one column per (station,
                    // class) silently drops the jobs in service.
                    jline.lang.Network gmodel = ctmcSolver.getModel();
                    jline.lang.NetworkStruct gsn = gmodel.getStruct();
                    java.util.List<Integer> spaceWidths = new java.util.ArrayList<Integer>();
                    for (int isf = 0; isf < gsn.nstateful; isf++) {
                        Matrix blk = gsn.space == null ? null
                                : gsn.space.get(gmodel.getStatefulNodes().get(isf));
                        spaceWidths.add(Integer.valueOf(blk == null ? 0 : blk.getNumCols()));
                    }
                    gen.put("spaceWidths", spaceWidths);
                    return gen;
                }
                break;

            case "statevec":
                if (solverObj instanceof SolverFluid) {
                    SolverFluid fluidSolver = (SolverFluid) solverObj;
                    fluidSolver.getAvg();  // ensure the ODE solve has produced odeStateVec
                    jline.solvers.SolverResult fres = fluidSolver.getResults();
                    java.util.Map<String, Object> sv = new LinkedHashMap<String, Object>();
                    sv.put("xvec", (fres instanceof jline.solvers.fluid.FluidResult)
                            ? ((jline.solvers.fluid.FluidResult) fres).odeStateVec : null);
                    return sv;
                }
                break;

            // Second-order report of the moment-closure methods. It exists only
            // here: the covariance is the whole content of "minnormal", and a
            // delegating caller that cannot read it has no way to tell a closure
            // solve from a first-order one except by the numbers.
            case "moments":
                if (solverObj instanceof SolverFluid) {
                    SolverFluid momSolver = (SolverFluid) solverObj;
                    momSolver.getAvg();
                    jline.solvers.SolverResult mres = momSolver.getResults();
                    if (!(mres instanceof jline.solvers.fluid.FluidResult)) {
                        return null;
                    }
                    return momentsToMap((jline.solvers.fluid.FluidResult) mres);
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
                    // Ensure the fork-join analysis has run so percentile results exist.
                    ((SolverMAM) solverObj).getAvg();
                    return ((SolverMAM) solverObj).getPerctRespT(percentiles);
                }
                break;

            // Transient analysis types
            case "tran-avg":
                if (solverObj instanceof NetworkSolver) {
                    NetworkSolver ns = (NetworkSolver) solverObj;
                    ns.getTranAvg();
                    jline.solvers.SolverResult tr = ns.getResults();
                    java.util.Map<String, Object> tran = new java.util.LinkedHashMap<String, Object>();
                    tran.put("t", tr != null ? tr.t : null);
                    tran.put("QNt", matrixGridToList(tr != null ? tr.QNt : null));
                    tran.put("UNt", matrixGridToList(tr != null ? tr.UNt : null));
                    tran.put("TNt", matrixGridToList(tr != null ? tr.TNt : null));
                    return tran;
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
                if (solverObj instanceof SolverCTMC && stateMatrix != null) {
                    // A system state is one row of per-class counts per station;
                    // the model interchange carries none, so without it the
                    // query is answered at the default initialization.
                    return ((SolverCTMC) solverObj).getProbSys(stateMatrix);
                }
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getProbSys();
                }
                break;

            case "prob-sys-aggr":
                if (solverObj instanceof SolverCTMC && stateMatrix != null) {
                    return ((SolverCTMC) solverObj).getProbSysAggr(stateMatrix);
                }
                if (solverObj instanceof NetworkSolver) {
                    return ((NetworkSolver) solverObj).getProbSysAggr();
                }
                break;

            // The joint marginal over the whole system, the getProbSysMarg the
            // JAR CLI had no name for while the C++ one served it as -a sysmarg.
            case "prob-sys-marg":
                if (solverObj instanceof NetworkSolver) {
                    if (stateMatrix == null) {
                        throw new IllegalArgumentException(
                            "-a prob-sys-marg needs --state: the marginal is evaluated at a "
                            + "stated population vector.");
                    }
                    return ((NetworkSolver) solverObj).getProbSysMarg(stateMatrix);
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
    /**
     * Second-order report of the moment-closure methods, in the field layout of
     * MATLAB {@code @SolverFLD/getMoments} and its python twin. Index blocks
     * travel as lists of integers rather than {@code int[][]}, which the JSON
     * writer would stringify.
     */
    /**
     * SolverUQ's interval report as a structured map: the lower and upper
     * envelope of each mean metric over the design points, plus how the range
     * was obtained. `exact` says whether it is the true support (an MVA
     * monotonicity argument) or a sampled estimate, which the numbers alone do
     * not reveal -- so a delegating caller reading only the bounds would report
     * a Monte Carlo range as a guarantee.
     *
     * @param iv the interval, never null
     * @return the map published under the `interval` key
     */
    private static Map<String, Object> intervalToMap(SolverUQ.Interval iv) {
        Map<String, Object> out = new LinkedHashMap<String, Object>();
        out.put("QLen", boundsPair(iv.Qlo, iv.Qup));
        out.put("Util", boundsPair(iv.Ulo, iv.Uup));
        out.put("RespT", boundsPair(iv.Rlo, iv.Rup));
        out.put("Tput", boundsPair(iv.Tlo, iv.Tup));
        out.put("WaitT", boundsPair(iv.Wlo, iv.Wup));
        out.put("X", iv.X);
        out.put("Rtot", iv.Rtot);
        out.put("exact", Boolean.valueOf(iv.exact));
        out.put("method", iv.method);
        out.put("reason", iv.reason);
        return out;
    }

    /** One metric's lower/upper envelope. */
    private static Map<String, Object> boundsPair(Matrix lo, Matrix up) {
        Map<String, Object> out = new LinkedHashMap<String, Object>();
        out.put("lower", lo);
        out.put("upper", up);
        return out;
    }

    private static Map<String, Object> momentsToMap(jline.solvers.fluid.FluidResult fr) {
        Map<String, Object> out = new LinkedHashMap<String, Object>();
        out.put("Sigma", fr.momentSigma);
        out.put("QVar", fr.momentQVar);
        out.put("sigma2", fr.momentSigma2);
        // The variance the DRIFT was closed at, which is not sigma2: a delegating
        // caller that measures a passage time has to re-integrate on the same
        // drift, and without this it falls back to the first-order min() and
        // reports the empty-station service law (cdf_respt_open_twoclasses read
        // mean 1.0007 for a response time whose closure mean is 1.5010).
        out.put("sigma2Drift", fr.momentSigma2Drift);
        out.put("outerIters", Integer.valueOf(fr.momentOuterIters));
        out.put("stationBlock", indexBlocksToList(fr.momentStationBlock));
        out.put("classBlock", indexBlocksToList(fr.momentClassBlock));
        if (fr.momentCacheSigma != null) {
            List<Object> caches = new ArrayList<Object>();
            for (int c = 0; c < fr.momentCacheSigma.length; c++) {
                Matrix sigma = fr.momentCacheSigma[c];
                Matrix pi0 = (fr.momentCachePi0 != null && c < fr.momentCachePi0.length)
                        ? fr.momentCachePi0[c] : null;
                Map<String, Object> entry = new LinkedHashMap<String, Object>();
                entry.put("node", Integer.valueOf(
                        (fr.momentCacheNode != null && c < fr.momentCacheNode.length)
                                ? fr.momentCacheNode[c] : -1));
                entry.put("pi0", pi0);
                entry.put("Sigma", sigma);
                // the occupancy variance is the diagonal of the item block, which
                // the caller would otherwise have to know how to cut out
                int n = (pi0 == null) ? 0 : pi0.length();
                Matrix pi0Var = new Matrix(n, 1, n);
                for (int i = 0; i < n && sigma != null && i < sigma.getNumRows(); i++) {
                    pi0Var.set(i, 0, Math.max(0.0, sigma.get(i, i)));
                }
                entry.put("pi0Var", pi0Var);
                if (fr.momentCacheMissProbVar != null && c < fr.momentCacheMissProbVar.getNumRows()) {
                    int K = fr.momentCacheMissProbVar.getNumCols();
                    Matrix mv = new Matrix(K, 1, K);
                    for (int r = 0; r < K; r++) {
                        mv.set(r, 0, fr.momentCacheMissProbVar.get(c, r));
                    }
                    entry.put("missProbVar", mv);
                }
                caches.add(entry);
            }
            out.put("cache", caches);
        }
        return out;
    }

    /** {@code int[][]} as a list of integer lists, or null. */
    private static List<Object> indexBlocksToList(int[][] blocks) {
        if (blocks == null) {
            return null;
        }
        List<Object> out = new ArrayList<Object>();
        for (int i = 0; i < blocks.length; i++) {
            List<Object> row = new ArrayList<Object>();
            if (blocks[i] != null) {
                for (int j = 0; j < blocks[i].length; j++) {
                    row.add(Integer.valueOf(blocks[i][j]));
                }
            }
            out.add(row);
        }
        return out;
    }

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

        // Cache subclasses must be checked before AvgTable (they extend it).
        if (obj instanceof NetworkAvgCacheTable) {
            return cacheTableToJSON((NetworkAvgCacheTable) obj);
        }

        if (obj instanceof NetworkAvgItemTable) {
            return itemTableToJSON((NetworkAvgItemTable) obj);
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

        if (obj instanceof DistributionResult) {
            return distributionResultToJSON((DistributionResult) obj);
        }

        if (obj instanceof ProbabilityResult) {
            return probabilityResultToJSON((ProbabilityResult) obj);
        }

        if (obj instanceof FJPercentileResult) {
            return fjPercentileResultToJSON((FJPercentileResult) obj);
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

        // Emit an unquoted literal: the default branch below would stringify it,
        // and "false" is truthy to most consumers.
        if (obj instanceof Boolean) {
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
        json.append(",\"event\":").append(eventListToJSON(result.event));
        json.append("}");
        return json.toString();
    }

    /**
     * Serialize a list of simulation Events (arrival/departure traces) as JSON,
     * so the Python-native SSA sampler can rebuild its SamplePath under lang="java".
     */
    private static String eventListToJSON(java.util.List<jline.lang.Event> events) {
        StringBuilder json = new StringBuilder();
        json.append("[");
        if (events != null) {
            for (int i = 0; i < events.size(); i++) {
                jline.lang.Event e = events.get(i);
                if (i > 0) json.append(",");
                json.append("{\"t\":").append(e.getT());
                json.append(",\"node\":").append(e.getNode());
                json.append(",\"jobclass\":").append(e.getJobClass());
                json.append(",\"event\":\"").append(e.getEvent()).append("\"}");
            }
        }
        json.append("]");
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

        // Every probability getter answers with one of these, so without an arm
        // here `-a normconst`, `-a prob*` and `-a prob-sys-marg` rendered as
        // `jline.io.Ret$ProbabilityResult@7161d8d1` in readable mode while the
        // JSON output carried the numbers in full.
        if (obj instanceof jline.io.Ret.ProbabilityResult) {
            jline.io.Ret.ProbabilityResult pr = (jline.io.Ret.ProbabilityResult) obj;
            StringBuilder sb = new StringBuilder();
            if (!Double.isNaN(pr.logNormalizingConstant)) {
                sb.append("logNormConstAggr: ").append(pr.logNormalizingConstant);
                sb.append(System.lineSeparator());
            }
            if (pr.nodeIndex != null) {
                sb.append("node: ").append(pr.nodeIndex).append(System.lineSeparator());
            }
            if (pr.state != null && !pr.state.isEmpty()) {
                sb.append("state: ").append(pr.state.toString()).append(System.lineSeparator());
            }
            sb.append("aggregate: ").append(pr.isAggregated).append(System.lineSeparator());
            if (pr.probability != null && !pr.probability.isEmpty()) {
                sb.append("probability: ").append(pr.probability.toString());
            }
            return sb.toString();
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

        // The VIEW names itself, as line-cli's envelope already does: a consumer
        // that reads the columns has to know which two label columns it got.
        json.append("\"type\":\"").append(avgTableJSONType(table)).append("\",");

        // see _kb/12-interfaces-and-docs.md (every AvgTable view carries columns)
        if (table instanceof NetworkAvgTable) {
            NetworkAvgTable nt = (NetworkAvgTable) table;
            appendStringArray(json, "Station", nt.getStationNames());
            appendStringArray(json, "JobClass", nt.getClassNames());
            appendDoubleArray(json, "QLen", nt.getQLen());
            appendDoubleArray(json, "Util", nt.getUtil());
            appendDoubleArray(json, "RespT", nt.getRespT());
            appendDoubleArray(json, "ResidT", nt.getResidT());
            appendDoubleArray(json, "ArvR", nt.getArvR());
            appendDoubleArray(json, "Tput", nt.getTput());
        } else if (table instanceof NetworkAvgNodeTable) {
            NetworkAvgNodeTable nt = (NetworkAvgNodeTable) table;
            appendStringArray(json, "Node", nt.getNodeNames());
            appendStringArray(json, "JobClass", nt.getClassNames());
            appendDoubleArray(json, "QLen", nt.getQLen());
            appendDoubleArray(json, "Util", nt.getUtil());
            appendDoubleArray(json, "RespT", nt.getRespT());
            appendDoubleArray(json, "ResidT", nt.getResidT());
            appendDoubleArray(json, "ArvR", nt.getArvR());
            appendDoubleArray(json, "Tput", nt.getTput());
        } else if (table instanceof NetworkAvgChainTable) {
            NetworkAvgChainTable ct = (NetworkAvgChainTable) table;
            appendStringArray(json, "Station", ct.getStationNames());
            appendStringArray(json, "Chain", ct.getChainNames());
            appendStringArray(json, "JobClasses", ct.getInChainNames());
            appendDoubleArray(json, "QLen", ct.getQLen());
            appendDoubleArray(json, "Util", ct.getUtil());
            appendDoubleArray(json, "RespT", ct.getRespT());
            appendDoubleArray(json, "ResidT", ct.getResidT());
            appendDoubleArray(json, "ArvR", ct.getArvR());
            appendDoubleArray(json, "Tput", ct.getTput());
        } else if (table instanceof NetworkAvgNodeChainTable) {
            NetworkAvgNodeChainTable ct = (NetworkAvgNodeChainTable) table;
            appendStringArray(json, "Node", ct.getNodeNames());
            appendStringArray(json, "Chain", ct.getChainNames());
            appendStringArray(json, "JobClasses", ct.getInChainNames());
            appendDoubleArray(json, "QLen", ct.getQLen());
            appendDoubleArray(json, "Util", ct.getUtil());
            appendDoubleArray(json, "RespT", ct.getRespT());
            appendDoubleArray(json, "ResidT", ct.getResidT());
            appendDoubleArray(json, "ArvR", ct.getArvR());
            appendDoubleArray(json, "Tput", ct.getTput());
        } else if (table instanceof NetworkAvgSysTable) {
            NetworkAvgSysTable st = (NetworkAvgSysTable) table;
            appendStringArray(json, "Chain", st.getChainNames());
            appendStringArray(json, "JobClasses", st.getInChainNames());
            appendDoubleArray(json, "SysRespT", st.getSysRespT());
            appendDoubleArray(json, "SysTput", st.getSysTput());
        } else if (table instanceof LayeredNetworkAvgTable) {
            // see _kb/12-interfaces-and-docs.md (structured JSON columns are NetworkAvgTable-only)
            LayeredNetworkAvgTable lt = (LayeredNetworkAvgTable) table;
            appendStringArray(json, "Node", lt.getNodeNames());
            appendStringArray(json, "NodeType", lt.getNodeTypes());
            appendDoubleArray(json, "QLen", lt.getQLen());
            appendDoubleArray(json, "Util", lt.getUtil());
            appendDoubleArray(json, "RespT", lt.getRespT());
            appendDoubleArray(json, "ResidT", lt.getResidT());
            appendDoubleArray(json, "ArvR", lt.getArvR());
            appendDoubleArray(json, "Tput", lt.getTput());
        }

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
     * The `type` an AvgTable view announces in JSON, matching line-cli's names.
     * @param table the table being serialized
     * @return the view name (AvgTable, AvgNodeTable, AvgChainTable, ...)
     */
    private static String avgTableJSONType(AvgTable table) {
        if (table instanceof NetworkAvgNodeTable) {
            return "AvgNodeTable";
        }
        if (table instanceof NetworkAvgChainTable) {
            return "AvgChainTable";
        }
        if (table instanceof NetworkAvgNodeChainTable) {
            return "AvgNodeChainTable";
        }
        if (table instanceof NetworkAvgSysTable) {
            return "AvgSysTable";
        }
        return "AvgTable";
    }

    /**
     * Converts a Matrix[][] grid (station x class transient trajectories) into a
     * nested List so objectToJSON serializes each cell via matrixToJSON.
     */
    private static List<List<Object>> matrixGridToList(Matrix[][] grid) {
        List<List<Object>> out = new ArrayList<List<Object>>();
        if (grid == null) {
            return out;
        }
        for (int i = 0; i < grid.length; i++) {
            List<Object> row = new ArrayList<Object>();
            for (int k = 0; k < grid[i].length; k++) {
                row.add(grid[i][k]);
            }
            out.add(row);
        }
        return out;
    }

    /**
     * Appends a JSON array of strings under the given key, with a trailing comma.
     */
    private static void appendStringArray(StringBuilder json, String key, List<String> values) {
        json.append("\"").append(key).append("\":[");
        if (values != null) {
            for (int i = 0; i < values.size(); i++) {
                if (i > 0) json.append(",");
                json.append("\"").append(escapeJSONString(values.get(i))).append("\"");
            }
        }
        json.append("],");
    }

    /**
     * Appends a JSON array of doubles under the given key, with a trailing comma.
     *
     * An INFINITY rides as the string "Infinity" / "-Infinity", the same spelling
     * the model direction already accepts for every numeric field, because JSON
     * has no infinite literal and a null is how this wire says "no value" -- which
     * a host reads as zero. The two are different answers: a queue at its
     * stability boundary has an infinite mean waiting time, and reporting 0 there
     * inverts the conclusion. A NaN keeps riding as null, which is what the
     * readable table's nan means. line-cli emits the identical form.
     */
    private static void appendDoubleArray(StringBuilder json, String key, List<Double> values) {
        json.append("\"").append(key).append("\":[");
        if (values != null) {
            for (int i = 0; i < values.size(); i++) {
                if (i > 0) json.append(",");
                Double v = values.get(i);
                if (v == null || v.isNaN()) {
                    json.append("null");
                } else if (v.isInfinite()) {
                    json.append(v.doubleValue() > 0 ? "\"Infinity\"" : "\"-Infinity\"");
                } else {
                    json.append(v.doubleValue());
                }
            }
        }
        json.append("],");
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
     * Convert a DistributionResult (response/passage time CDFs) to structured JSON.
     * Emits per-(station,class) probability and time arrays. Each cdfData matrix is
     * (npoints x 2) with column 0 = probability and column 1 = time.
     */
    private static String distributionResultToJSON(DistributionResult dr) {
        if (dr == null) {
            return "null";
        }
        StringBuilder json = new StringBuilder();
        json.append("{\"type\":\"DistributionResult\"");
        json.append(",\"numStations\":").append(dr.numStations);
        json.append(",\"numClasses\":").append(dr.numClasses);
        json.append(",\"isTransient\":").append(dr.isTransient);
        json.append(",\"cdf\":[");
        for (int i = 0; i < dr.cdfData.size(); i++) {
            if (i > 0) json.append(",");
            json.append("[");
            java.util.List<Matrix> row = dr.cdfData.get(i);
            for (int j = 0; j < row.size(); j++) {
                if (j > 0) json.append(",");
                Matrix cdf = row.get(j);
                json.append("{\"p\":[");
                if (cdf != null && cdf.getNumRows() > 0) {
                    for (int q = 0; q < cdf.getNumRows(); q++) {
                        if (q > 0) json.append(",");
                        json.append(cdf.get(q, 0));
                    }
                }
                json.append("],\"t\":[");
                if (cdf != null && cdf.getNumRows() > 0 && cdf.getNumCols() > 1) {
                    for (int q = 0; q < cdf.getNumRows(); q++) {
                        if (q > 0) json.append(",");
                        json.append(cdf.get(q, 1));
                    }
                }
                json.append("]}");
            }
            json.append("]");
        }
        json.append("]}");
        return json.toString();
    }

    /**
     * Convert a ProbabilityResult to structured JSON: the scalar probability, the
     * full probability matrix (for aggregated/marginal vectors), the log
     * normalizing constant, and the aggregation flag.
     */
    private static String probabilityResultToJSON(ProbabilityResult pr) {
        if (pr == null) {
            return "null";
        }
        StringBuilder json = new StringBuilder();
        json.append("{\"type\":\"ProbabilityResult\"");
        double scalar = 0.0;
        try {
            scalar = pr.getScalarProbability();
        } catch (Exception e) {
            scalar = 0.0;
        }
        json.append(",\"scalar\":").append(scalar);
        json.append(",\"logNormalizingConstant\":").append(pr.logNormalizingConstant);
        json.append(",\"isAggregated\":").append(pr.isAggregated);
        if (pr.probability != null && pr.probability.getNumRows() * pr.probability.getNumCols() > 0) {
            json.append(",\"probability\":").append(matrixToJSON(pr.probability));
        }
        json.append("}");
        return json.toString();
    }

    /**
     * Convert a fork-join FJPercentileResult to JSON (percentile levels 0..100 and
     * their response-time values).
     */
    private static String fjPercentileResultToJSON(FJPercentileResult r) {
        if (r == null) {
            return "null";
        }
        StringBuilder json = new StringBuilder();
        json.append("{\"type\":\"FJPercentileResult\",\"K\":").append(r.K);
        json.append(",\"percentiles\":");
        json.append(doubleArrayToJSON(r.percentiles));
        json.append(",\"RTp\":");
        json.append(doubleArrayToJSON(r.RTp));
        json.append("}");
        return json.toString();
    }

    /**
     * Convert a NetworkAvgCacheTable to structured JSON with the same columns as
     * the native getAvgCacheTable output.
     */
    private static String cacheTableToJSON(NetworkAvgCacheTable t) {
        if (t == null) {
            return "null";
        }
        StringBuilder json = new StringBuilder();
        json.append("{\"type\":\"NetworkAvgCacheTable\",");
        appendStringArray(json, "Node", t.getNodeNames());
        appendStringArray(json, "JobClass", t.getClassNames());
        appendDoubleArray(json, "List", t.getList());
        appendDoubleArray(json, "ListCap", t.getListCap());
        appendDoubleArray(json, "Items", t.getItems());
        appendDoubleArray(json, "HitProb", t.getHitProb());
        appendDoubleArray(json, "DelayedHitProb", t.getDelayedHitProb());
        appendDoubleArray(json, "MissProb", t.getMissProb());
        appendDoubleArray(json, "HitRate", t.getHitRate());
        appendDoubleArray(json, "DelayedHitRate", t.getDelayedHitRate());
        appendDoubleArray(json, "MissRate", t.getMissRate());
        appendDoubleArray(json, "ArvR", t.getArvR());
        appendDoubleArray(json, "ResidT", t.getResidT());
        // strip the trailing comma left by the last append
        if (json.charAt(json.length() - 1) == ',') {
            json.deleteCharAt(json.length() - 1);
        }
        json.append("}");
        return json.toString();
    }

    /**
     * Convert a NetworkAvgItemTable to structured JSON with the same columns as
     * the native getAvgItemTable output.
     */
    private static String itemTableToJSON(NetworkAvgItemTable t) {
        if (t == null) {
            return "null";
        }
        StringBuilder json = new StringBuilder();
        json.append("{\"type\":\"NetworkAvgItemTable\",");
        appendStringArray(json, "Node", t.getNodeNames());
        appendDoubleArray(json, "Item", t.getItem());
        appendDoubleArray(json, "List", t.getList());
        appendDoubleArray(json, "ListCap", t.getListCap());
        appendDoubleArray(json, "Prob", t.getProb());
        appendDoubleArray(json, "DelayedHitQLen", t.getDelayedHitQLen());
        appendDoubleArray(json, "DelayedHitQLenFull", t.getDelayedHitQLenFull());
        if (json.charAt(json.length() - 1) == ',') {
            json.deleteCharAt(json.length() - 1);
        }
        json.append("}");
        return json.toString();
    }

    // ========================================================================
    // SOLVE subcommand: java -jar ldes.jar solve model.json -o result.json [options]
    // ========================================================================

    /**
     * Handles the 'solve' subcommand for LDES simulation.
     *
     * <p>Delegates to {@link LdesCLI#handleSolveCommand}, which is the single
     * implementation of this subcommand. This class previously carried its own
     * copy, which drifted: the copy never gained the {@code -e}/{@code --maxevents}
     * and {@code --maxtime} options, so whether they worked depended on which main
     * class the bundle had been built with. Delegating keeps the two entry points
     * (this class and {@link LdesCLI}) behaviorally identical by construction.
     *
     * @param args command-line arguments after "solve"
     * @return 0 on success, non-zero on error
     */
    private static int handleSolveCommand(String[] args) {
        return LdesCLI.handleSolveCommand(args);
    }

    /**
     * Main entry point for the LINE CLI
     * @param args Command line arguments
     */
    public static void main(String[] args) {
        try {
            // Check for 'solve' subcommand (used by Python native LDES solver)
            if (args.length > 0 && "solve".equals(args[0])) {
                String[] solveArgs = new String[args.length - 1];
                System.arraycopy(args, 1, solveArgs, 0, solveArgs.length);
                int exitCode = handleSolveCommand(solveArgs);
                System.exit(exitCode);
                return;
            }

            String result = parseArgs(args);
            if (result != null) {
                System.out.println(result);
            } else if (usageError) {
                System.exit(1);
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
