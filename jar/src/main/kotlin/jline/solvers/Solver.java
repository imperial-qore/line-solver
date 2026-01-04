/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.lang.Model;
import jline.lang.Network;
import jline.VerboseLevel;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.*;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

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
                "default", "exact", "auto", "ctmc", "ctmc.gpu", "gpu", "mva", "mva.exact", "mva.amva", "mva.qna", "sqrt", "mva.sqrt",
                "amva", "amva.bs", "amva.qd", "bs", "qd", "amva.qli", "qli", "amva.fli", "fli", "amva.aql", "aql", "amva.qdaql", "qdaql", "amva.lin", "lin", "amva.qdlin", "qdlin",
                "nrm", "ssa", "ssa.serial.hash", "ssa.para.hash", "ssa.parallel.hash", "ssa.serial", "ssa.para", "ssa", "taussa", "tauleap",
                "ssa.parallel", "serial.hash", "serial", "para", "parallel", "para.hash", "parallel.hash",
                "jmt", "jsim", "jmva", "jmva.amva", "jmva.mva", "jmva.recal", "jmva.mom", "jmva.comom", "jmva.chow", "jmva.bs", "jmva.aql", "jmva.lin", "jmva.dmlin", "jmva.ls",
                "jmt.jsim", "jmt.jmva", "jmt.jmva.mva", "jmt.jmva.amva", "jmt.jmva.recal", "jmt.jmva.comom", "jmt.jmva.chow", "jmt.jmva.bs", "jmt.jmva.aql", "jmt.jmva.lin", "jmt.jmva.dmlin", "jmt.jmva.ls",
                "brute", "ca", "comomrm", "comomld", "gm", "mom", "propfair", "recal", "kt", "rd", "nr.probit", "nr.logit", "nc.brute", "nc.ca", "nc.comom", "nc.comomld", "nc.gm", "nc.mom", "nc.propfair", "nc.recal", "nc.kt", "nc.rd", "nc.nr.probit", "nc.nr.logit",
                "fluid", "matrix", "softmin", "statedep", "closing", "fluid.softmin", "fluid.statedep", "fluid.closing", "fluid.matrix",
                "nc", "nc.exact", "nc.imci", "ls", "nc.ls", "nc.cub", "cub", "le", "nc.le", "nc.panacea", "panacea", "nc.mmint2", "mmint2", "nc.gleint", "gleint", "mam", "dec.source", "dec.mmap",
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
        
        // Parse key-value pairs
        for (int i = 0; i < varargin.length; i++) {
            if (varargin[i] instanceof String) {
                String key = (String) varargin[i];
                
                // Handle single keywords without values
                if (key.equals("exact")) {
                    options.method("exact");
                    continue;
                }
                
                // Handle key-value pairs
                if (i + 1 < varargin.length) {
                    Object value = varargin[i + 1];
                    parseOptionPair(options, key, value);
                    i++; // Skip the value
                }
            }
        }
        
        return options;
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
     * Performs validation checks before running the analyzer.
     * Verifies model compatibility and method validity.
     *
     * @param options the solver options to validate
     * @throws RuntimeException if validation fails
     */
    public void runAnalyzerChecks(SolverOptions options) {
        List<String> allMethods = listValidOptions().get("allMethods");

        if (this.enableChecks && !supports((Network) this.model)) {
            line_error(mfilename(new Object() {
            }), "This model contains features not supported by the solver.");
        }
        if (this.enableChecks && !allMethods.contains(options.method)) {
            line_error(mfilename(new Object() {
            }), "The " + options.method + " method is unsupported by this solver.");
        }
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
