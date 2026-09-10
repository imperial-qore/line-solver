package jline.solvers.auto;

import jline.solvers.SolverOptions;

/**
 * Options specific to the AUTO solver
 */
public class AUTOptions extends SolverOptions {

    /**
     * Method for solver selection. Either a SELECTION INTENT -- "default" or
     * "heur" (heuristic), "sim", "exact", "fast", "accurate", "bound" -- or a
     * method FAMILY that pins the delegate outright ("nc", "mva", ...,
     * optionally qualified as "nc.comom").
     *
     * <p>"ai" and "nn" are NOT among them. Both were documented here while
     * {@code SolverAUTO.METHOD_AI} sat commented out as "not yet available", so
     * this javadoc invited a caller to ask for a classifier that does not exist;
     * {@code SolverAUTO.resolveMethodToken} now rejects an unknown token by name
     * and lists what it will take.
     */
    public String selectionMethod = "default";

    /**
     * Force a specific solver (overrides automatic selection)
     * Valid values: "mva", "nc", "mam", "fluid", "jmt", "ssa", "ctmc"
     */
    public String forceSolver = null;

    /**
     * Default constructor
     */
    public AUTOptions() {
        super();
    }

    /**
     * Constructor with method
     */
    public AUTOptions(String selectionMethod) {
        super();
        this.selectionMethod = selectionMethod;
    }

    /**
     * Copy constructor from base options
     */
    public AUTOptions(SolverOptions options) {
        super();
        // Copy base options
        this.verbose = options.verbose;
        this.samples = options.samples;
        this.seed = options.seed;
        this.timespan = options.timespan;
        this.iter_max = options.iter_max;
        this.iter_tol = options.iter_tol;
        this.tol = options.tol;
        this.cutoff = options.cutoff;
        this.init_sol = options.init_sol;
        this.remote = options.remote;
        this.remote_endpoint = options.remote_endpoint;
        this.cache = options.cache;
        this.keep = options.keep;
        this.force = options.force;
        this.hide_immediate = options.hide_immediate;
        this.lang = options.lang;
        this.method = options.method;
        this.stiff = options.stiff;
        this.config = options.config;
        this.odesolvers = options.odesolvers;
    }

    public SolverOptions copy() {
        AUTOptions copy = new AUTOptions();

        // Copy base options using super's clone
        SolverOptions baseClone = super.copy();
        copy.verbose = baseClone.verbose;
        copy.samples = baseClone.samples;
        copy.seed = baseClone.seed;
        copy.timespan = baseClone.timespan;
        copy.iter_max = baseClone.iter_max;
        copy.iter_tol = baseClone.iter_tol;
        copy.tol = baseClone.tol;
        copy.cutoff = baseClone.cutoff;
        copy.init_sol = baseClone.init_sol;
        copy.remote = baseClone.remote;
        copy.remote_endpoint = baseClone.remote_endpoint;
        copy.cache = baseClone.cache;
        copy.keep = baseClone.keep;
        copy.force = baseClone.force;
        copy.hide_immediate = baseClone.hide_immediate;
        copy.lang = baseClone.lang;
        copy.method = baseClone.method;
        copy.stiff = baseClone.stiff;
        copy.config = baseClone.config;
        copy.odesolvers = baseClone.odesolvers;

        // Copy AUTO-specific options
        copy.selectionMethod = this.selectionMethod;
        copy.forceSolver = this.forceSolver;

        return copy;
    }
}