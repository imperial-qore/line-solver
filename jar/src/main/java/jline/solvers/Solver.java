/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.GlobalConstants;
import jline.lang.Model;
import jline.lang.Network;
import jline.lang.FeatureSet;
import jline.VerboseLevel;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.*;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Abstract base class for model solution algorithms and analysis tools.
 * <p>
 * This class provides the fundamental infrastructure for solving queueing models
 * using various analytical and simulation algorithms. It manages solver configuration,
 * result storage, random number generation, and validation of solver options.
 * <p>
 * Concrete implementations must provide the {@link #runAnalyzer()} method to perform
 * the actual model solution.
 *
 * @see SolverOptions
 * @see SolverResult
 * @see NetworkSolver
 */
public abstract class Solver {

    /**
     * The model to be solved
     */
    public Model model;

    /**
     * Name identifier for this solver instance
     */
    public String name;

    /**
     * Configuration options for the solver
     */
    public SolverOptions options;

    /**
     * Results from the most recent solver execution
     */
    public SolverResult result;

    /**
     * Flag controlling whether to perform validation checks
     */
    public boolean enableChecks;

    /**
     * Random number generator for stochastic algorithms
     */
    public Random random;

    /**
     * Constructs a solver with the specified name and options.
     *
     * @param name    the solver name identifier
     * @param options configuration options for the solver
     */
    protected Solver(String name, SolverOptions options) {
        this(null, name, options);
    }

    /**
     * Constructs a solver with the specified model, name, and options.
     *
     * @param model   the model to be solved
     * @param name    the solver name identifier
     * @param options configuration options for the solver
     */
    protected Solver(Model model, String name, SolverOptions options) {
        this.model = model;
        this.name = name;
        this.setOptions(options.copy());
        //this.result = new SolverResult();
        this.enableChecks = true;
        // Set thread-local seed via ThreadLocalRandom which delegates to RandomManager
        // This maintains backward compatibility with existing test expectations
        RandomManager.getThreadRandom().setSeed(options.seed);
    }

    /**
     * Constructs a solver with the specified name using default options.
     *
     * @param name the solver name identifier
     */
    protected Solver(String name) {
        this(name, defaultOptions());
    }

    /**
     * Returns a new SolverOptions instance with default settings.
     *
     * @return default solver options
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(null);
    }

    /**
     * Returns lists of valid options and methods supported by solvers.
     *
     * @return map containing "allOptions" and "allMethods" lists
     */
    public static Map<String, List<String>> listValidOptions() {

        List<String> allOptions = Arrays.asList(
                "cache", "cutoff", "force", "init_sol", "iter_max", "iter_tol", "lang", "tol",
                "keep", "method", "odesolvers", "samples", "seed", "stiff", "timespan", "verbose", "config.multiserver"
        );

        List<String> allMethods = Arrays.asList(
                "cache", "cutoff", "force", "init_sol", "iter_max", "iter_tol", "lang", "tol",
                "keep", "method", "odesolvers", "samples", "seed", "stiff", "timespan", "verbose", "config.multiserver",
                "default", "exact", "auto", "heur", "tree", "sim", "fast", "accurate", "bound", // SolverAUTO selection method names
                "ctmc", "ctmc.gpu", "gpu", "mva", "mva.exact", "mva.amva", "mva.qna", "sqrt", "mva.sqrt",
                "amva", "amva.bs", "amva.qd", "bs", "qd", "amva.qli", "qli", "amva.fli", "fli", "amva.aql", "aql", "amva.qdaql", "qdaql", "amva.lin", "lin", "amva.qdlin", "qdlin",
                "nrm", "ssa", "ssa.serial.hash", "ssa.para.hash", "ssa.parallel.hash", "ssa.serial", "ssa.para", "ssa",
                "ssa.parallel", "serial.hash", "serial", "para", "parallel", "para.hash", "parallel.hash",
                "jmt", "jsim", "jmva", "jmva.amva", "jmva.mva", "jmva.recal", "jmva.mom", "jmva.comom", "jmva.chow", "jmva.bs", "jmva.aql", "jmva.lin", "jmva.dmlin", "jmva.ls",
                "jmt.jsim", "jmt.jmva", "jmt.jmva.mva", "jmt.jmva.amva", "jmt.jmva.recal", "jmt.jmva.comom", "jmt.jmva.chow", "jmt.jmva.bs", "jmt.jmva.aql", "jmt.jmva.lin", "jmt.jmva.dmlin", "jmt.jmva.ls",
                "brute", "ca", "comom", "comomrm", "comomld", "gm", "propfair", "recal", "kt", "bkt", "lekt", "rd", "nr.probit", "nr.logit", "nc.brute", "nc.ca", "nc.comom", "nc.comomld", "nc.gm", "nc.propfair", "nc.recal", "nc.kt", "nc.bkt", "nc.lekt", "nc.rd", "nc.nr.probit", "nc.nr.logit",
                "fluid", "matrix", "softmin", "statedep", "closing", "fluid.softmin", "fluid.statedep", "fluid.closing", "fluid.matrix",
                "nc", "nc.exact", "nc.imci", "ls", "nc.ls", "nc.cub", "cub", "le", "nc.le", "ble", "nc.ble", "aghq", "nc.aghq", "mcmc", "nc.mcmc", "nc.pana", "pana", "nc.panald", "panald", "nc.mmint2", "mmint2", "nc.gleint", "gleint", "mam", "dec.source", "dec.mmap",
                "mmk", "gigk", "gigk.kingman_approx",
                "mm1", "mg1", "gm1", "gig1", "gim1", "gig1.kingman", "gig1.gelenbe", "gig1.heyman", "gig1.kimura", "gig1.allen", "gig1.kobayashi", "gig1.klb", "gig1.marchal",
                "aba.upper", "aba.lower", "gb.upper", "gb.lower", "sb.upper", "sb.lower", "bjb.upper", "bjb.lower", "pb.upper", "pb.lower"
        );

        Map<String, List<String>> lists = new HashMap<>();
        lists.put("allOptions", allOptions);
        lists.put("allMethods", allMethods);

        return lists;
    }

    /**
     * Parses option parameters into a SolverOptions data structure.
     *
     * @param varargin variable arguments in key-value pairs
     * @return parsed solver options
     * @throws IllegalArgumentException if arguments are invalid
     */
    public static SolverOptions parseOptions(Object... varargin) throws IllegalArgumentException {
        SolverOptions options = new SolverOptions();
        return parseOptions(options, varargin);
    }

    /**
     * Parses option parameters into an existing SolverOptions instance.
     *
     * @param options  existing options object to modify
     * @param varargin variable arguments in key-value pairs
     * @return modified solver options
     * @throws IllegalArgumentException if arguments are invalid
     */
    public static SolverOptions parseOptions(SolverOptions options, Object... varargin) throws IllegalArgumentException {
        // Handle empty arguments
        if (varargin == null || varargin.length == 0) {
            return options;
        }
        
        // If first argument is SolverOptions, use it as base
        if (varargin.length == 1 && varargin[0] instanceof SolverOptions) {
            return ((SolverOptions) varargin[0]).copy();
        }
        
        // Handle single numeric argument as cutoff value (for backward compatibility)
        if (varargin.length == 1 && (varargin[0] instanceof Number)) {
            if (varargin[0] instanceof Double) {
                options.cutoff((Double) varargin[0]);
            } else if (varargin[0] instanceof Integer) {
                options.cutoff((Integer) varargin[0]);
            }
            return options;
        }
        
        // Parse key-value pairs, plus the positional method name
        for (int i = 0; i < varargin.length; i++) {
            if (varargin[i] instanceof String) {
                String key = (String) varargin[i];

                if (isOptionKey(key)) {
                    if (i + 1 >= varargin.length) {
                        // A dangling option name used to be dropped in silence,
                        // leaving the solver on the default for that option.
                        line_error(mfilename(new Object() {
                        }), String.format("Option '%s' was given without a value.", key));
                    }
                    parseOptionPair(options, key, varargin[i + 1]);
                    i++;   // skip the value
                    continue;
                }

                // An unrecognized name followed by a non-string is still a
                // key-value pair whose key parseOptionPair does not implement
                // (e.g. "warmup", 1e4): consume both and let it fall through
                // there, as before. Reading such a key as a method instead would
                // break every caller that passes one.
                if (i + 1 < varargin.length && !(varargin[i + 1] instanceof String)) {
                    parseOptionPair(options, key, varargin[i + 1]);
                    i++;
                    continue;
                }

                // What is left is the METHOD, which is how MATLAB and python
                // spell SolverCTMC(model,'mdd'). Dropping it in silence (the
                // previous behaviour, which only exempted the literal "exact")
                // ran the DEFAULT method under the name of another one, so every
                // comparison of that method against the default agreed trivially.
                options.method(key);
            }
        }

        return options;
    }

    /**
     * Names {@link #parseOptionPair} recognizes as an option, as opposed to a
     * method name given positionally.
     *
     * <p>Kept in step with the switch in parseOptionPair by hand: a key missing
     * here is read as a method name and its value is then parsed as a second
     * positional argument, which is why the two lists must agree.</p>
     */
    private static final Set<String> OPTION_KEYS = new HashSet<String>(Arrays.asList(
            "cache", "compress", "container", "cutoff", "eventcache", "force", "fork_join",
            "hide_immediate", "highvar", "init_sol", "interlocking", "iter_max", "iter_tol",
            "keep", "lang", "merge", "method", "multiserver", "np_priority", "pstar",
            "remote", "remote_endpoint", "resturl", "rest_url", "samples", "seed",
            "space_max", "state_space_gen", "stiff", "timespan", "tol", "verbose"));

    /** True when the string names an option rather than a method. */
    private static boolean isOptionKey(String key) {
        if (key == null) {
            return false;
        }
        return key.startsWith("config.") || OPTION_KEYS.contains(key.toLowerCase());
    }
    
    /**
     * Parses a single key-value pair for solver options.
     *
     * @param options the options object to modify
     * @param key     the option name
     * @param value   the option value
     */
    private static void parseOptionPair(SolverOptions options, String key, Object value) {
        // Handle config.* options with dot notation
        if (key.startsWith("config.")) {
            String configKey = key.substring(7); // Remove "config." prefix
            parseConfigOption(options, configKey, value);
            return;
        }
        
        // Handle standard options
        switch (key.toLowerCase()) {
            case "cache":
                if (value instanceof Boolean) {
                    options.cache = (Boolean) value;
                }
                break;
            case "cutoff":
                if (value instanceof Double) {
                    options.cutoff((Double) value);
                } else if (value instanceof Integer) {
                    options.cutoff((Integer) value);
                } else if (value instanceof String) {
                    options.cutoff(Double.parseDouble((String) value));
                } else if (value instanceof Matrix) {
                    options.cutoff((Matrix) value);
                } else if (value instanceof double[][]) {
                    options.cutoff(new Matrix((double[][]) value));
                } else if (value instanceof int[][]) {
                    int[][] src = (int[][]) value;
                    double[][] dst = new double[src.length][];
                    for (int i = 0; i < src.length; i++) {
                        dst[i] = new double[src[i].length];
                        for (int j = 0; j < src[i].length; j++) {
                            dst[i][j] = src[i][j];
                        }
                    }
                    options.cutoff(new Matrix(dst));
                } else {
                    // A silently dropped cutoff leaves the solver on its default
                    // scalar cutoff, which changes the truncated state space and
                    // hence the results, so refuse the value instead.
                    line_error(mfilename(new Object() {
                    }), String.format("Unsupported type for option 'cutoff': %s.",
                            value == null ? "null" : value.getClass().getName()));
                }
                break;
            case "force":
                if (value instanceof Boolean) {
                    options.force = (Boolean) value;
                }
                break;
            case "hide_immediate":
                if (value instanceof Boolean) {
                    options.hide_immediate = (Boolean) value;
                }
                break;
            case "init_sol":
                if (value instanceof Matrix) {
                    options.init_sol = (Matrix) value;
                }
                break;
            case "iter_max":
                if (value instanceof Integer) {
                    options.iter_max = (Integer) value;
                } else if (value instanceof Double) {
                    options.iter_max = ((Double) value).intValue();
                } else if (value instanceof String) {
                    options.iter_max = Integer.parseInt((String) value);
                }
                break;
            case "iter_tol":
                if (value instanceof Double) {
                    options.iter_tol = (Double) value;
                } else if (value instanceof Integer) {
                    options.iter_tol = ((Integer) value).doubleValue();
                } else if (value instanceof String) {
                    options.iter_tol = Double.parseDouble((String) value);
                }
                break;
            case "keep":
                if (value instanceof Boolean) {
                    options.keep = (Boolean) value;
                }
                break;
            case "lang":
                if (value instanceof String) {
                    options.lang = (String) value;
                }
                break;
            case "method":
                if (value instanceof String) {
                    options.method = (String) value;
                }
                break;
            case "remote":
                if (value instanceof Boolean) {
                    options.remote = (Boolean) value;
                }
                break;
            case "remote_endpoint":
                if (value instanceof String) {
                    options.remote_endpoint = (String) value;
                }
                break;
            case "restUrl":
            case "rest_url":
                if (value instanceof String) {
                    options.restUrl = (String) value;
                }
                break;
            case "container":
                if (value instanceof String) {
                    options.container = (String) value;
                }
                break;
            case "samples":
                if (value instanceof Integer) {
                    options.samples = (Integer) value;
                } else if (value instanceof Double) {
                    options.samples = ((Double) value).intValue();
                } else if (value instanceof String) {
                    options.samples = Integer.parseInt((String) value);
                }
                break;
            case "seed":
                if (value instanceof Integer) {
                    options.seed = (Integer) value;
                } else if (value instanceof Double) {
                    options.seed = ((Double) value).intValue();
                } else if (value instanceof String) {
                    options.seed = Integer.parseInt((String) value);
                }
                break;
            case "stiff":
                if (value instanceof Boolean) {
                    options.stiff = (Boolean) value;
                }
                break;
            case "timespan":
                if (value instanceof double[]) {
                    options.timespan = (double[]) value;
                } else if (value instanceof Matrix) {
                    Matrix m = (Matrix) value;
                    if (m.getNumCols() == 2 && m.getNumRows() == 1) {
                        options.timespan = new double[]{m.get(0, 0), m.get(0, 1)};
                    }
                }
                break;
            case "tol":
                if (value instanceof Double) {
                    options.tol = (Double) value;
                } else if (value instanceof Integer) {
                    options.tol = ((Integer) value).doubleValue();
                } else if (value instanceof String) {
                    options.tol = Double.parseDouble((String) value);
                }
                break;
            case "warmup":
                // an absolute count of samples to discard; resolved against
                // options.samples when the run is configured, since the two may
                // be given in either order
                if (value instanceof Number) {
                    options.config.warmup = ((Number) value).doubleValue();
                } else if (value instanceof String) {
                    options.config.warmup = Double.parseDouble((String) value);
                } else {
                    line_error(mfilename(new Object() {
                    }), String.format("Unsupported type for option 'warmup': %s.",
                            value == null ? "null" : value.getClass().getName()));
                }
                break;
            case "verbose":
                if (value instanceof Boolean) {
                    options.verbose((Boolean) value);
                } else if (value instanceof VerboseLevel) {
                    options.verbose = (VerboseLevel) value;
                } else if (value instanceof String) {
                    String strVal = ((String) value).toLowerCase();
                    switch (strVal) {
                        case "silent":
                        case "false":
                            options.verbose = VerboseLevel.SILENT;
                            break;
                        case "std":
                        case "standard":
                        case "true":
                            options.verbose = VerboseLevel.STD;
                            break;
                        case "debug":
                            options.verbose = VerboseLevel.DEBUG;
                            break;
                    }
                }
                break;
            default:
                // Ignore unknown options (consistent with MATLAB behavior)
                break;
        }
    }
    
    /**
     * Parses configuration sub-options (config.*).
     *
     * @param options   the options object to modify
     * @param configKey the configuration key
     * @param value     the configuration value
     */
    private static void parseConfigOption(SolverOptions options, String configKey, Object value) {
        if (options.config == null) {
            options.config = new SolverOptions.Config();
        }
        
        switch (configKey.toLowerCase()) {
            case "highvar":
                if (value instanceof String) {
                    options.config.highvar = (String) value;
                }
                break;
            case "multiserver":
                if (value instanceof String) {
                    options.config.multiserver = (String) value;
                }
                break;
            case "np_priority":
                if (value instanceof String) {
                    options.config.np_priority = (String) value;
                }
                break;
            case "fork_join":
                if (value instanceof String) {
                    options.config.fork_join = (String) value;
                }
                break;
            case "merge":
                if (value instanceof String) {
                    options.config.merge = (String) value;
                }
                break;
            case "compress":
                if (value instanceof String) {
                    options.config.compress = (String) value;
                }
                break;
            case "space_max":
                if (value instanceof Integer) {
                    options.config.space_max = (Integer) value;
                } else if (value instanceof Double) {
                    options.config.space_max = ((Double) value).intValue();
                } else if (value instanceof String) {
                    options.config.space_max = Integer.parseInt((String) value);
                }
                break;
            case "interlocking":
                if (value instanceof Boolean) {
                    options.config.interlocking = (Boolean) value;
                }
                break;
            case "eventcache":
                if (value instanceof Boolean) {
                    options.config.eventcache = (Boolean) value;
                }
                break;
            case "hide_immediate":
                if (value instanceof Boolean) {
                    options.config.hide_immediate = (Boolean) value;
                }
                break;
            case "state_space_gen":
                if (value instanceof String) {
                    options.config.state_space_gen = (String) value;
                }
                break;
            case "pstar":
                if (value instanceof List) {
                    @SuppressWarnings("unchecked")
                    List<Double> pstarList = (List<Double>) value;
                    options.config.pstar = pstarList;
                } else if (value instanceof double[]) {
                    double[] arr = (double[]) value;
                    options.config.pstar = new ArrayList<>();
                    for (double d : arr) {
                        options.config.pstar.add(d);
                    }
                }
                break;
            default:
                // Ignore unknown config options
                break;
        }
    }

    /**
     * Returns the name identifier of this solver.
     *
     * @return the solver name
     */
    public String getName() {
        return name;
    }

    /**
     * Returns the current solver options.
     *
     * @return the solver options
     */
    public SolverOptions getOptions() {
        return options;
    }

    /**
     * Sets new solver options.
     *
     * @param options the new solver options to set
     */
    public void setOptions(SolverOptions options) {
        this.options = options;
    }

    /**
     * Returns the results from the most recent solver execution.
     *
     * @return the solver results
     */
    public SolverResult getResults() {
        return result;
    }

    /**
     * Checks if the solver has computed results.
     *
     * @return true if results are available, false otherwise
     */
    public boolean hasResults() {
        if (result != null && result.QN == null) {
            return false;
        } else {
            return !result.QN.isEmpty();
        }
    }

    /**
     * Checks whether the solver, with its currently configured method, returns
     * stochastic estimates, i.e. results that depend on the random seed, as in
     * simulation or Monte Carlo integration.
     *
     * <p>A solver run with method "default" may resolve the actual method only
     * at runtime; once results are available the classification therefore uses
     * the method recorded in {@code result.method} (e.g. "default/imci" when
     * the NC default path resolved to Monte Carlo integration).</p>
     *
     * @return true if the solver returns stochastic estimates
     */
    public boolean isStochastic() {
        String method = (options != null) ? options.method : null;
        if (result != null && result.method != null && !result.method.isEmpty()) {
            method = result.method;
        }
        return isStochasticMethod(method);
    }

    /**
     * Classifies a (possibly runtime-resolved) method name of this solver as
     * stochastic. Deterministic by default; subclasses with simulation-based
     * or sampling-based methods override this.
     *
     * @param method the method name to classify
     * @return true if the method returns stochastic estimates
     */
    public boolean isStochasticMethod(String method) {
        return false;
    }

    /**
     * Checks if Java runtime is available for solver execution.
     * Always returns true in this Java implementation.
     *
     * @return true indicating Java is available
     */
    public boolean isJavaAvailable() {
        return true;
    }

    /**
     * Checks if the specified option name is valid for this solver.
     *
     * @param optName the option name to validate
     * @return true if the option is valid, false otherwise
     */
    public boolean isValidOption(String optName) {
        Map<String, List<String>> options = listValidOptions();
        List<String> allOpts = options.get("allOpt");
        return allOpts.contains(optName);
    }

    /**
     * Clears previously stored results and resets the random number generator.
     */
    public void reset() {
        this.result.reset();
        resetRandomGeneratorSeed(this.options.seed);
    }

    /**
     * Assigns a new seed to the random number generator.
     * This sets the master seed for all random number generation in the system.
     *
     * @param seed the seed value for random number generation
     */
    public void resetRandomGeneratorSeed(long seed) {
        // Set thread-local seed via ThreadLocalRandom which delegates to RandomManager
        // This maintains backward compatibility with existing test expectations
        RandomManager.getThreadRandom().setSeed((int) seed);
    }

    /**
     * Executes the solver algorithm to analyze the model.
     * This abstract method must be implemented by concrete solver classes.
     *
     * @throws IllegalAccessException       if access to required resources is denied
     * @throws ParserConfigurationException if XML parsing configuration fails
     * @throws IOException                  if I/O operations fail
     */
    public abstract void runAnalyzer() throws IllegalAccessException, ParserConfigurationException, IOException;

    /**
     * Tests whether a wall-clock time budget has been exceeded. Used as a
     * cooperative checkpoint inside iterative solver loops.
     *
     * @param startNanos     value of System.nanoTime() captured at solver launch
     * @param timeoutSeconds wall-clock budget in seconds (Double.POSITIVE_INFINITY = no budget)
     * @return true if (now - startNanos) exceeds timeoutSeconds
     */
    public static boolean timeExceeded(long startNanos, double timeoutSeconds) {
        if (Double.isInfinite(timeoutSeconds) || timeoutSeconds <= 0) {
            return false;
        }
        return (System.nanoTime() - startNanos) / 1e9 > timeoutSeconds;
    }

    /**
     * Performs validation checks before running the analyzer.
     * Verifies model compatibility and method validity.
     *
     * @param options the solver options to validate
     * @throws RuntimeException if validation fails
     */
    public void runAnalyzerChecks(SolverOptions options) {
        // Propagate solver verbose level to global so that model-level
        // messages (e.g., priority info in refreshStruct) respect it
        if (options != null) {
            GlobalConstants.Verbose = options.verbose;
        }
        if (!this.enableChecks) {
            return;
        }
        List<String> allMethods = listValidOptions().get("allMethods");
        if (options != null && !allMethods.contains(options.method)) {
            line_error(mfilename(new Object() {
            }), "The " + options.method + " method is unsupported by this solver.");
            return;
        }
        // Method-aware feature gate: resolve the concrete method that will run
        // (for options.method='default' this may map to a specific method) and
        // gate against that method's per-method feature set rather than the
        // solver-level union.
        String method = resolveMethod(options);
        String reason = supportsModelMethod(method);
        if (!reason.isEmpty()) {
            String prefix = (options != null && method.equals(options.method))
                    ? "This model contains features not supported by the solver. "
                    : "This model contains features not supported by the solver's '" + method + "' method. ";
            line_error(mfilename(new Object() {
            }), prefix + reason);
        }
    }

    /**
     * Resolve the concrete method that will run. Base behavior returns
     * options.method unchanged; solvers that perform feature-driven selection
     * for options.method='default' override this (typically via selectMethod).
     *
     * @param options the solver options
     * @return the concrete method name
     */
    public String resolveMethod(SolverOptions options) {
        return options == null ? "default" : options.method;
    }

    /**
     * Per-method feature set, or null to signal "this solver does not diverge
     * per method" (the coarse supports(model) is then used, preserving any
     * structural checks it carries). Divergent solvers (e.g. MVA, MAM, NC)
     * override this to return the base envelope with per-method deltas applied.
     *
     * @param method the concrete method name
     * @return the per-method FeatureSet, or null
     */
    public FeatureSet getMethodFeatureSet(String method) {
        return null;
    }

    /**
     * Does this solver produce transient averages, i.e. does getTranAvg return
     * trajectories on a finite options.timespan? Declared false here and
     * overridden by the solvers that populate result.Tran (Fluid, CTMC, LDES,
     * JMT). It is a capability claim, not a state test: it must answer before
     * any run has taken place, because the MAP/MMPP random-environment fallback
     * uses it to decide whether the environment stages can be coupled by the
     * mean-field analyzer (which needs getTranAvg) or only by the two
     * steady-state limits.
     *
     * @return true if the solver can return transient averages
     */
    public boolean supportsTransientAnalysis() {
        return false;
    }

    /**
     * Fine, method-aware gate. Returns an empty string when the model fits the
     * concrete METHOD, else a human-readable reason. Base behavior derives the
     * answer from getMethodFeatureSet(method); when that is null the solver's
     * own supports(model) is used. Solvers with non-feature-set structural
     * per-method rules override this.
     *
     * @param method the concrete method name
     * @return empty string if supported, else the offending reason
     */
    public String supportsModelMethod(String method) {
        FeatureSet fs = getMethodFeatureSet(method);
        if (fs == null) {
            return supports((Network) this.model) ? "" : "Some features are not supported by the chosen solver.";
        }
        return FeatureSet.supportsReason(fs, ((Network) this.model).getUsedLangFeatures());
    }

    /**
     * Feature-driven method selection: returns the first method in
     * preferenceList whose per-method feature set covers the model, or the last
     * entry when none fully covers (the gate then reports the precise features).
     *
     * @param preferenceList ordered candidate methods
     * @return the selected method name
     */
    public String selectMethod(List<String> preferenceList) {
        String chosen = preferenceList.get(preferenceList.size() - 1);
        for (String cand : preferenceList) {
            if (supportsModelMethod(cand).isEmpty()) {
                return cand;
            }
        }
        return chosen;
    }

    /**
     * Enables or disables validation checks during solver execution.
     *
     * @param bool true to enable checks, false to disable
     */
    public void setChecks(boolean bool) {
        enableChecks = bool;
    }

    /**
     * Checks if this solver supports the given network model.
     * Default implementation returns true; subclasses should override
     * to provide specific feature validation.
     *
     * @param model the network model to check
     * @return true if the model is supported, false otherwise
     */
    public boolean supports(Network model) {
        return true;
    }

    // NOTE: the following LINE methods have not been migrated to JLINE
    // - isValidOption() - all options are always available as part of SolverOptions class
    // - supports() - static method at specific Solver level rather than abstract within Solver class
}
