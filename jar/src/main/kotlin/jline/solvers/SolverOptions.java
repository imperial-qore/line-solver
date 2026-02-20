/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers;

import jline.GlobalConstants;
import static jline.GlobalConstants.Inf;
import jline.lang.constant.SolverType;
import jline.VerboseLevel;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.List;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

/**
 * Configuration options for queueing network solvers.
 * <p>
 * This class contains all configurable parameters that control solver behavior,
 * including convergence criteria, numerical tolerances, algorithmic choices,
 * and output preferences. Options can be customized per solver type or
 * globally across all solvers.
 * <p>
 * Many options are solver-specific and may be ignored by solvers that don't
 * support them. The class provides builder-style setter methods for convenient
 * configuration chaining.
 *
 * @see Solver
 * @see SolverType
 * @see VerboseLevel
 */
public class SolverOptions {

    /**
     * Enable result caching to avoid recomputation
     */
    public boolean cache;

    /**
     * Cutoff threshold for numerical computations.
     * Can be either a scalar (uniform cutoff for all stations and classes)
     * or a matrix (station × class specific cutoffs).
     * Matrix dimensions should be nstations × nclasses.
     */
    public Matrix cutoff;

    /**
     * Advanced configuration options
     */
    public Config config;

    /**
     * Force solver execution even when validation fails
     */
    public boolean force;

    /**
     * Hide immediate transitions from output when possible
     */
    public boolean hide_immediate;

    /**
     * Initial solution for iterative solvers
     */
    public Matrix init_sol;

    /**
     * Maximum number of iterations for iterative algorithms
     */
    public int iter_max;

    /**
     * Tolerance for iteration convergence
     */
    public double iter_tol;

    /**
     * General numerical tolerance
     */
    public double tol;

    /**
     * Keep intermediate results and temporary files
     */
    public boolean keep;

    /**
     * Target language for execution ("java", "matlab", etc.)
     */
    public String lang;

    /**
     * Solution algorithm/method to use
     */
    public String method;

    /**
     * Enable remote solver execution
     */
    public boolean remote;

    /**
     * Remote solver endpoint address
     */
    public String remote_endpoint;

    /**
     * ODE solver configurations for fluid analysis
     */
    public ODESolvers odesolvers;

    /**
     * Number of samples/operations for simulation-based methods.
     *
     * <p>Semantics by solver:
     * <ul>
     *   <li><b>DES</b>: Maximum service completion events (default: 1,000,000)
     *       - Simulation stops when this many completions are processed
     *       - 20% of events used for warmup (MSER-5 for automatic truncation)
     *   </li>
     *   <li><b>SSA</b>: Maximum random number generation operations (default: 10,000)
     *       - Approximately 2 RNG ops per reaction event (Gillespie algorithm)
     *   </li>
     *   <li><b>JMT</b>: Samples per performance metric, minimum 5000 (default: 10,000)
     *       - Statistical accuracy parameter, not simulation duration
     *   </li>
     *   <li><b>NC</b>: Monte Carlo samples for sampling-based methods (default: 100,000)
     *   </li>
     *   <li><b>LQNS</b>: Simulation length parameter for lqsim tool
     *   </li>
     *   <li><b>Analytical solvers</b> (MVA, MAM, CTMC, Fluid, etc.): Not used
     *   </li>
     * </ul>
     */
    public int samples;

    /**
     * Random number generator seed for reproducibility
     */
    public int seed;

    /**
     * Use stiff ODE solvers for numerical integration
     */
    public boolean stiff;

    /**
     * Time interval for transient analysis [start, end]
     */
    public double[] timespan;

    /**
     * Fixed timestep for transient analysis (null for adaptive stepping)
     */
    public Double timestep;

    /**
     * Verbosity level for solver output
     */
    public VerboseLevel verbose;

    /**
     * Number of value iterations for CTMC reward computation.
     * Used by SolverCTMC for computing cumulative rewards via value iteration.
     * Default is 1000 iterations.
     */
    public Integer rewardIterations;

    /**
     * Confidence interval level for simulation-based solvers.
     * When set to a value between 0 and 1 (e.g., 0.95 for 95% confidence),
     * the solver will compute and return confidence interval bounds.
     * When set to 0 or negative, confidence intervals are disabled.
     * Default is 0 (disabled).
     */
    public double confint;

    /**
     * Creates SolverOptions with default settings.
     * Initializes all parameters to sensible defaults suitable for most models.
     */
    public SolverOptions() {

        // Solver Default Options
        this.cache = true;
        this.cutoff = null; // Will be initialized when network dimensions are known
        this.config = new Config();
        this.config.highvar = "default";
        this.config.multiserver = "default";
        this.config.np_priority = "default";
        this.config.fork_join = "default";
        this.config.eventcache = true;
        this.config.hide_immediate = false; // Default for config
        this.config.pstar = new ArrayList<>();
        this.config.variates = "none"; // Variance reduction disabled by default
        this.config.nonmkv = "bernstein"; // Method for non-Markovian distribution conversion
        this.config.nonmkvorder = 20; // Order (number of phases) for non-Markovian approximation
        this.config.da = "courtois"; // CTMC decomposition/aggregation method for Env solver
        this.config.da_iter = 10; // Number of iterations for kms/takahashi
        this.config.relax = "none"; // Default is no relaxation (will be overridden for LN solver)
        this.config.relax_factor = 0.1; // Relaxation factor when enabled
        this.config.relax_min = 0.1; // Minimum relaxation factor for adaptive mode
        this.config.relax_history = 5; // Error history window for adaptive mode
        this.force = false;
        this.hide_immediate = true; // Hide immediate transitions if possible
        this.init_sol = new Matrix(0, 0);
        this.iter_max = 100;
        this.iter_tol = 0.0001; // Convergence tolerance to stop iterations
        this.tol = 0.0001; // Tolerance for all other uses
        this.keep = true;
        this.lang = "java";
        this.method = "default";
        this.remote = false;
        this.remote_endpoint = "127.0.0.1";

        this.odesolvers = new ODESolvers();
        this.odesolvers.odeminstep = 0.001;
        //this.odeMinStep = 0.00000001;
        this.odesolvers.odemaxstep = Inf;
        this.odesolvers.fastODESolver = new ClassicalRungeKuttaIntegrator(this.odesolvers.odemaxstep);
        this.odesolvers.accurateODESolver =
                new DormandPrince54Integrator(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol);
        this.odesolvers.fastStiffODESolver =
                new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 3, 3);
        this.odesolvers.accurateStiffODESolver =
                new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 12, 5);

        this.samples = 10000;
        this.seed = RandomManager.generateRandomSeed();
        this.stiff = true;
        this.timespan = new double[2];
        this.timespan[0] = Inf;
        this.timespan[1] = Inf;
        this.timestep = null;
        this.verbose = GlobalConstants.getInstance().getVerbose();
        this.confint = 0; // disabled by default
    }

    /**
     * Creates SolverOptions with defaults customized for a specific solver type.
     *
     * @param solverType the type of solver to configure defaults for, or null for generic defaults
     */
    public SolverOptions(SolverType solverType) {
        this();
        if (solverType == null) {
            return;
        }

        // Solver-specific Defaults
        switch (solverType) {
            case CTMC:
                this.timespan = new double[]{Inf, Inf};
                this.verbose = VerboseLevel.STD;
                this.config.hide_immediate = true;
                this.config.state_space_gen = "full";
                this.rewardIterations = 1000; // Number of value iterations for reward computation
                break;
            case ENV:
                this.iter_max = 100;
                this.verbose = VerboseLevel.SILENT;
                break;
            case FLUID:
                this.config.highvar = "default";
                this.config.hide_immediate = true; // Eliminate immediate transitions by default
                this.config.put("immediate_tol", 1e7); // Threshold for detecting immediate transitions
                this.iter_max = 200;
                this.stiff = true;
                this.timespan[0] = 0;
                // Reduce min step to handle Immediate transition rates (~1e8)
                // MATLAB's ode15s has effectively no minimum step constraint
                this.odesolvers.odeminstep = 1e-14;
                this.odesolvers.accurateODESolver =
                    new DormandPrince54Integrator(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol);
                this.odesolvers.fastStiffODESolver =
                    new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 3, 3);
                this.odesolvers.accurateStiffODESolver =
                    new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 12, 5);
                break;
            case LN:
                this.config.interlocking = false; // do not change, true unstable as of 2.0.38
                this.config.multiserver = "default";
                // Under-relaxation options for convergence improvement
                this.config.relax = "auto"; // 'auto' | 'fixed' | 'adaptive' | 'none'
                this.config.relax_factor = 0.1; // Relaxation factor when enabled (0 < omega <= 1)
                this.config.relax_min = 0.1; // Minimum relaxation factor for adaptive mode
                this.config.relax_history = 5; // Error history window for adaptive mode
                // MOL (Method of Layers) options for hierarchical iteration
                this.config.mol_task_inner_max = 50; // Max task inner iterations per host outer iteration
                this.config.mol_task_inner_tol = 1e-4; // Task layer convergence tolerance
                this.config.mol_host_outer_tol = 1e-4; // Host layer convergence tolerance
                this.config.mol_min_steps = 2; // Minimum outer iterations before checking convergence
                this.timespan = new double[]{Inf, Inf};
                this.keep = true;
                this.iter_max = 200; // More iterations for difficult LQN models
                this.iter_tol = 5e-3; // Convergence tolerance (looser than default for LQN models)
                this.tol = 1e-4;
                break;
            case LQNS:
                this.timespan = new double[]{Inf, Inf};
                this.keep = true;
                this.verbose = VerboseLevel.SILENT;  // MATLAB has "false" which maps to SILENT
                this.config.multiserver = "rolia";
                break;
            case MAM:
                this.iter_max = 100;
                this.timespan = new double[]{Inf, Inf};
                // num_cdf_pts uses global default of 200
                break;
            case MVA:
                this.iter_max = 1000;
                this.iter_tol = 1e-6;
                break;
            case NC:
                this.samples = 100000;
                this.timespan = new double[]{Inf, Inf};
                this.config.highvar = "interp";
                break;
            case SSA:
                this.timespan = new double[]{0, Inf};
                this.verbose = VerboseLevel.STD; // MATLAB has "true" 
                this.config.state_space_gen = "none";
                // eventcache is set based on lang in MATLAB, but Java always uses eventcache = true
                this.config.eventcache = true;
                break;
            case JMT:
                // JMT uses default settings (as per MATLAB implementation)
                break;
            case QNS:
                this.config.multiserver = "rolia";
                break;
            default: // Global options unless overridden by a solver
        }
    }

    /**
     * Creates a deep copy of this SolverOptions instance.
     * All fields including nested objects are properly cloned.
     *
     * @return a deep copy of this options object
     */
    public SolverOptions copy() {
        SolverOptions cloned = new SolverOptions();

        cloned.cache = this.cache;
        cloned.cutoff = this.cutoff != null ? this.cutoff.copy() : null;

        // Cloning config with null checks
        if (this.config != null) {
            cloned.config = new Config();
            cloned.config.highvar = this.config.highvar;
            cloned.config.multiserver = this.config.multiserver;
            cloned.config.np_priority = this.config.np_priority;
            cloned.config.pstar = this.config.pstar != null ? new ArrayList<>(this.config.pstar) : null;
            cloned.config.variates = this.config.variates;
            cloned.config.fork_join = this.config.fork_join;
            cloned.config.merge = this.config.merge;
            cloned.config.compress = this.config.compress;
            cloned.config.space_max = this.config.space_max;
            cloned.config.interlocking = this.config.interlocking;
            cloned.config.eventcache = this.config.eventcache;
            cloned.config.hide_immediate = this.config.hide_immediate;
            cloned.config.state_space_gen = this.config.state_space_gen;
            cloned.config.nonmkv = this.config.nonmkv;
            cloned.config.nonmkvorder = this.config.nonmkvorder;
            cloned.config.da = this.config.da;
            cloned.config.da_iter = this.config.da_iter;
            cloned.config.relax = this.config.relax;
            cloned.config.relax_factor = this.config.relax_factor;
            cloned.config.relax_min = this.config.relax_min;
            cloned.config.relax_history = this.config.relax_history;
            cloned.config.mol_task_inner_max = this.config.mol_task_inner_max;
            cloned.config.mol_task_inner_tol = this.config.mol_task_inner_tol;
            cloned.config.mol_host_outer_tol = this.config.mol_host_outer_tol;
            cloned.config.mol_min_steps = this.config.mol_min_steps;
            cloned.config.num_cdf_pts = this.config.num_cdf_pts;
            cloned.config.remote = this.config.remote;
            cloned.config.remote_url = this.config.remote_url;
        }

        cloned.force = this.force;
        cloned.hide_immediate = this.hide_immediate;
        cloned.init_sol = this.init_sol != null ? this.init_sol : new Matrix(0, 0);
        cloned.iter_max = this.iter_max;
        cloned.iter_tol = this.iter_tol;
        cloned.tol = this.tol;
        cloned.keep = this.keep;
        cloned.lang = this.lang;
        cloned.method = this.method;
        cloned.remote = this.remote;
        cloned.remote_endpoint = this.remote_endpoint;

        // Cloning ODESolvers - create new instances because ODE solvers have mutable internal state
        if (this.odesolvers != null) {
            cloned.odesolvers = new ODESolvers();
            cloned.odesolvers.odeminstep = this.odesolvers.odeminstep;
            cloned.odesolvers.odemaxstep = this.odesolvers.odemaxstep;
            // Create fresh ODE solver instances to avoid shared mutable state
            cloned.odesolvers.fastODESolver = new ClassicalRungeKuttaIntegrator(cloned.odesolvers.odemaxstep);
            cloned.odesolvers.accurateODESolver = new DormandPrince54Integrator(
                    cloned.odesolvers.odeminstep, cloned.odesolvers.odemaxstep, this.tol, this.tol);
            cloned.odesolvers.fastStiffODESolver = new LSODA(
                    cloned.odesolvers.odeminstep, cloned.odesolvers.odemaxstep, this.tol, this.tol, 3, 3);
            cloned.odesolvers.accurateStiffODESolver = new LSODA(
                    cloned.odesolvers.odeminstep, cloned.odesolvers.odemaxstep, this.tol, this.tol, 12, 5);
        }

        cloned.samples = this.samples;
        cloned.seed = this.seed;
        cloned.stiff = this.stiff;
        cloned.timespan = this.timespan != null ? this.timespan.clone() : new double[]{Inf, Inf};
        cloned.timestep = this.timestep;
        cloned.verbose = this.verbose;
        cloned.confint = this.confint;

        return cloned;
    }

    /**
     * Sets the numerical cutoff threshold as a scalar value (builder pattern).
     * This creates a uniform cutoff matrix where all stations and classes use the same cutoff.
     * The matrix will be properly dimensioned when the solver runs.
     *
     * @param s cutoff value to apply uniformly
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions cutoff(int s) {
        return cutoff((double) s);
    }

    /**
     * Sets the numerical cutoff threshold as a scalar value (builder pattern).
     * This creates a uniform cutoff matrix where all stations and classes use the same cutoff.
     * The matrix will be properly dimensioned when the solver runs.
     *
     * @param s cutoff value to apply uniformly
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions cutoff(double s) {
        // Create a 1x1 matrix to indicate scalar cutoff - will be expanded later
        this.cutoff = new Matrix(1, 1);
        this.cutoff.set(0, 0, s);
        return this;
    }

    /**
     * Sets the numerical cutoff threshold as a matrix (builder pattern).
     * The matrix should have dimensions nstations × nclasses where each element
     * specifies the cutoff for a specific station-class combination.
     *
     * @param cutoffMatrix Matrix of cutoff values with dimensions nstations × nclasses
     * @return this SolverOptions instance for method chaining
     * @throws IllegalArgumentException if the matrix is null or has invalid dimensions
     */
    public SolverOptions cutoff(Matrix cutoffMatrix) {
        if (cutoffMatrix == null) {
            throw new IllegalArgumentException("Cutoff matrix cannot be null");
        }
        this.cutoff = cutoffMatrix.copy();
        return this;
    }

    /**
     * Ensures the cutoff is properly dimensioned for the given network structure.
     * If cutoff is null, initializes it to POSITIVE_INFINITY for all stations and classes.
     * If cutoff is a 1x1 matrix (scalar), expands it to full dimensions.
     * If cutoff is already properly dimensioned, leaves it unchanged.
     *
     * @param nstations Number of stations in the network
     * @param nclasses Number of job classes in the network
     * @return The properly dimensioned cutoff matrix
     */
    public Matrix getCutoffMatrix(int nstations, int nclasses) {
        if (this.cutoff == null) {
            // Initialize with infinite cutoffs
            Matrix cutoffMatrix = new Matrix(nstations, nclasses);
            cutoffMatrix.fill(Inf);
            return cutoffMatrix;
        } else if (this.cutoff.getNumRows() == 1 && this.cutoff.getNumCols() == 1) {
            // Expand scalar cutoff to full matrix
            Matrix cutoffMatrix = new Matrix(nstations, nclasses);
            double scalarValue = this.cutoff.get(0, 0);
            cutoffMatrix.fill(scalarValue);
            return cutoffMatrix;
        } else if (this.cutoff.getNumRows() == nstations && this.cutoff.getNumCols() == nclasses) {
            // Already properly dimensioned
            return this.cutoff.copy();
        } else {
            // Incorrect dimensions
            throw new IllegalArgumentException(
                String.format("Cutoff matrix has dimensions %dx%d, but network requires %dx%d (nstations × nclasses)",
                    this.cutoff.getNumRows(), this.cutoff.getNumCols(), nstations, nclasses)
            );
        }
    }

    /**
     * Sets whether to keep intermediate results and temporary files (builder pattern).
     *
     * @param s true to keep files, false to clean up
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions keep(boolean s) {
        this.keep = s;
        return this;
    }

    /**
     * Sets the solution method/algorithm (builder pattern).
     *
     * @param s method name (e.g., "mva", "ctmc", "ssa", "fluid")
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions method(String s) {
        this.method = s;
        return this;
    }

    /**
     * Sets the number of samples for simulation methods (builder pattern).
     *
     * @param s number of samples to generate
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions samples(int s) {
        this.samples = s;
        return this;
    }

    /**
     * Sets the random number generator seed (builder pattern).
     *
     * @param s seed value for reproducible random number generation
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions seed(int s) {
        this.seed = s;
        return this;
    }

    /**
     * Sets the maximum step size for ODE solvers and updates all integrators.
     *
     * @param odeMaxStep maximum step size for numerical integration
     */
    public void setODEMaxStep(double odeMaxStep) {
        this.odesolvers.odemaxstep = odeMaxStep;
        this.odesolvers.fastODESolver = new ClassicalRungeKuttaIntegrator(this.odesolvers.odemaxstep);
        this.odesolvers.accurateODESolver =
                new DormandPrince54Integrator(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol);
        this.odesolvers.fastStiffODESolver =
                new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 3, 3);
        this.odesolvers.accurateStiffODESolver =
                new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 12, 5);
    }

    /**
     * Sets the minimum step size for ODE solvers and updates all integrators.
     *
     * @param odeMinStep minimum step size for numerical integration
     */
    public void setODEMinStep(double odeMinStep) {
        this.odesolvers.odeminstep = odeMinStep;
        this.odesolvers.fastODESolver = new ClassicalRungeKuttaIntegrator(this.odesolvers.odeminstep);
        this.odesolvers.accurateODESolver =
                new DormandPrince54Integrator(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol);
        this.odesolvers.fastStiffODESolver =
                new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 3, 3);
        this.odesolvers.accurateStiffODESolver =
                new LSODA(this.odesolvers.odeminstep, this.odesolvers.odemaxstep, tol, tol, 12, 5);
    }

    /**
     * Sets the verbosity level for solver output (builder pattern).
     *
     * @param s verbosity level
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions verbose(VerboseLevel s) {
        this.verbose = s;
        return this;
    }

    /**
     * Sets the verbosity level using a boolean flag (builder pattern).
     *
     * @param s true for standard output, false for silent mode
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions verbose(boolean s) {
        if (s) {
            this.verbose = VerboseLevel.STD;
        } else {
            this.verbose = VerboseLevel.SILENT;
        }
        return this;
    }

    /**
     * Forces solver execution even when validation fails (builder pattern).
     *
     * @param force true to force execution despite validation failures
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions force(boolean force) {
        this.force = force;
        return this;
    }

    /**
     * Sets the confidence interval level for simulation-based solvers (builder pattern).
     * A value between 0 and 1 enables CI computation (e.g., 0.95 for 95% confidence).
     * A value of 0 or negative disables CI computation.
     *
     * @param level confidence level (0-1) or 0 to disable
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions confint(double level) {
        this.confint = level;
        return this;
    }

    /**
     * Advanced configuration options for specialized solver features.
     * These options control solver-specific behavior and may not be
     * applicable to all solver types.
     */
    public static class Config {

        /**
         * High-variance class handling strategy
         */
        public String highvar;

        /**
         * Multi-server scheduling strategy
         */
        public String multiserver;

        /**
         * Non-preemptive priority handling
         */
        public String np_priority;

        /**
         * P-norm smoothing parameters for fluid solvers
         */
        public List<Double> pstar;

        /**
         * Variance reduction technique for DES simulation.
         * <p>Available options:</p>
         * <ul>
         *   <li><b>"none"</b>: No variance reduction (standard simulation)</li>
         *   <li><b>"antithetic"</b>: Antithetic variates using synchronized 1-U method.
         *       Generates paired samples with negative correlation by using antithetic
         *       variates to reduce variance through negative correlation.</li>
         *   <li><b>"control"</b>: Control variates using mean-based correction.
         *       Applies post-hoc corrections based on deviation of sampled means
         *       from known theoretical means (E[arrival]=1/λ, E[service]=1/μ).</li>
         *   <li><b>"both"</b>: Combined antithetic and control variates</li>
         * </ul>
         * <p>Default: "none"</p>
         */
        public String variates;

        /**
         * Alpha parameter for Env solver Courtois decomposition (default: 0.015)
         * Controls the coupling threshold for environment state grouping
         */
        public Double env_alpha;

        /**
         * Fork-join handling strategy
         */
        public String fork_join;

        /**
         * State merging strategy
         */
        public String merge;

        /**
         * State space compression method
         */
        public String compress;

        /**
         * Maximum state space size
         */
        public int space_max;

        /**
         * Enable interlocking optimization
         */
        public boolean interlocking;

        /**
         * Enable event caching for SSA
         */
        public boolean eventcache;

        /**
         * Hide immediate transitions from analysis
         */
        public boolean hide_immediate;

        /**
         * State space generation strategy
         */
        public String state_space_gen;

        /**
         * Method for non-Markovian distribution conversion.
         * <p>Available options:</p>
         * <ul>
         *   <li><b>"none"</b>: No conversion, keep distributions as-is</li>
         *   <li><b>"bernstein"</b>: Convert using Bernstein polynomial approximation to phase-type</li>
         * </ul>
         * <p>Default: "bernstein"</p>
         */
        public String nonmkv = "bernstein";

        /**
         * Order (number of phases) for non-Markovian distribution approximation.
         * Higher values provide more accurate approximations but increase computational cost.
         * <p>Default: 20</p>
         */
        public int nonmkvorder = 20;

        /**
         * CTMC decomposition/aggregation method for Env solver.
         * <p>Available options:</p>
         * <ul>
         *   <li><b>"courtois"</b>: Courtois decomposition (default)</li>
         *   <li><b>"kms"</b>: Koury-McAllister-Stewart method</li>
         *   <li><b>"takahashi"</b>: Takahashi's method</li>
         *   <li><b>"multi"</b>: Multigrid method</li>
         * </ul>
         * <p>Default: "courtois"</p>
         */
        public String da;

        /**
         * Number of iterations for iterative decomposition/aggregation methods (kms, takahashi).
         * <p>Default: 10</p>
         */
        public int da_iter;

        /**
         * Under-relaxation mode for SolverLN convergence improvement.
         * <p>Available options:</p>
         * <ul>
         *   <li><b>"auto"</b>: Start without relaxation, enable when oscillation detected (default)</li>
         *   <li><b>"fixed"</b>: Always use relax_factor</li>
         *   <li><b>"adaptive"</b>: Adjust omega based on error trajectory</li>
         *   <li><b>"none"</b>: Disable relaxation</li>
         * </ul>
         * <p>Default: "auto"</p>
         */
        public String relax;

        /**
         * Relaxation factor (omega) when relaxation is enabled.
         * Value should be between 0 and 1, where lower values provide more damping.
         * <p>Default: 0.1</p>
         */
        public double relax_factor;

        /**
         * Minimum relaxation factor for adaptive mode.
         * <p>Default: 0.1</p>
         */
        public double relax_min;

        /**
         * Error history window size for adaptive relaxation mode.
         * <p>Default: 5</p>
         */
        public int relax_history;

// MOL (Method of Layers) options for hierarchical iteration
        /**
         * Maximum task inner iterations per host outer iteration.
         * <p>Default: 50</p>
         */
        public int mol_task_inner_max;

        /**
         * Task layer convergence tolerance (utilization delta).
         * <p>Default: 1e-4</p>
         */
        public double mol_task_inner_tol;

        /**
         * Host layer convergence tolerance (utilization delta).
         * <p>Default: 1e-4</p>
         */
        public double mol_host_outer_tol;

        /**
         * Minimum outer iterations before checking host convergence.
         * <p>Default: 2</p>
         */
        public int mol_min_steps;

        /**
         * Number of points for CDF (Cumulative Distribution Function) computation.
         * Used by SolverFluid and SolverMAM for response time distribution analysis.
         * <p>Default: 200 (100 for MAM solver)</p>
         */
        public int num_cdf_pts = 200;

        /**
         * Enable remote execution via REST API (LQNS-specific).
         * When true, solver uses HTTP to communicate with lqns-rest server.
         * <p>Default: false</p>
         */
        public boolean remote = false;

        /**
         * URL of lqns-rest server for remote execution (LQNS-specific).
         * <p>Default: "http://localhost:8080"</p>
         */
        public String remote_url = "http://localhost:8080";

        /**
         * Additional configuration parameters stored as key-value pairs
         */
        private java.util.Map<String, Object> additionalParams = new java.util.HashMap<>();

        /**
         * Store a configuration parameter
         */
        public void put(String key, Object value) {
            additionalParams.put(key, value);
        }

        /**
         * Retrieve a configuration parameter
         */
        public Object get(String key) {
            return additionalParams.get(key);
        }

        /**
         * Check if a configuration parameter exists
         */
        public boolean containsKey(String key) {
            return additionalParams.containsKey(key);
        }
    }

    /**
     * Configuration for ordinary differential equation solvers used in fluid analysis.
     * Contains different integrators optimized for various problem characteristics.
     */
    public static class ODESolvers {
        /**
         * Minimum step size for ODE integration
         */
        public double odeminstep;

        /**
         * Maximum step size for ODE integration
         */
        public double odemaxstep;

        /**
         * Fast integrator for non-stiff problems
         */
        public FirstOrderIntegrator fastODESolver;

        /**
         * Accurate integrator for non-stiff problems
         */
        public FirstOrderIntegrator accurateODESolver;

        /**
         * Fast integrator for stiff problems
         */
        public LSODA fastStiffODESolver;

        /**
         * Accurate integrator for stiff problems
         */
        public LSODA accurateStiffODESolver;
    }
}
