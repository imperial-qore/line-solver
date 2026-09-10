/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers;

import jline.GlobalConstants;
import jline.lang.constant.SolverType;
import jline.VerboseLevel;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.solvers.fluid.LSODAExt;
import jline.solvers.fluid.handlers.FluidRateMultiplier;
import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
import org.apache.commons.math3.ode.nonstiff.BogackiShampine23Integrator;
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
     * Environment sojourn mode for the SolverENV state-vector / mean-field
     * analyzers: null/"stochastic" (default) averages each stage's transient
     * over the random holding-time CDF; "deterministic" evaluates it at the
     * fixed mean stage duration.
     */
    public String sojourn;

    /**
     * Enable remote solver execution
     */
    public boolean remote;

    /**
     * Remote solver endpoint address
     */
    public String remote_endpoint;

    /**
     * Base URL of a solver REST server, e.g. "http://localhost:8080". Empty or
     * null means the solver runs its tool locally. Used by SolverLDES (through
     * LDESOptions.restUrl) and by SolverJMT, whose backend dispatch is
     * described in {@link jline.solvers.wrappers.jmt.JmtBackend}.
     */
    public String restUrl;

    /**
     * Docker image to dispatch an external tool through, overriding the
     * default candidates. Only consulted when no local binary is available.
     */
    public String container;

    /**
     * ODE solver configurations for fluid analysis
     */
    public ODESolvers odesolvers;

    /**
     * Statistical sample budget for simulation- and sampling-based methods.
     *
     * <p>Semantics by solver:
     * <ul>
     *   <li><b>LDES</b>: Deprecated alias of {@link #events} - maximum service
     *       completion events (default: 200,000, see LDESOptions.DEFAULT_SAMPLES);
     *       MSER-5 (or a fixed warmupfrac) governs warmup truncation
     *   </li>
     *   <li><b>SSA</b>: Deprecated alias of {@link #events} - number of fired
     *       transitions, one collected state sample per event (default: 10,000)
     *   </li>
     *   <li><b>JMT</b>: Samples collected per performance metric, minimum 5000
     *       (default: 10,000); JMT runs with disableStatisticStop, so this
     *       fixes the run length rather than adapting it to CI width
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
     * Event budget for discrete-event simulation solvers (SSA: fired
     * transitions; LDES: service completion events). Default -1 (unset): the
     * solvers fall back to {@link #samples}, which remains accepted as a
     * deprecated alias for the event budget. When positive, this value
     * overrides samples in SolverSSA and SolverLDES.
     */
    public int events;

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
     * Instants the caller wants the transient trajectory AT, increasing, or null.
     *
     * <p>When set, the fluid integrator adds them to its recorded grid from its
     * own continuous extension, instead of reporting only the steps it happened
     * to take. Port of MATLAB's {@code options.tranpoints}
     * (solver_fluid_iteration.m) and of the native-python twin;
     * {@link jline.solvers.env.SolverENV} sets it to the sojourn quadrature
     * grid, because reading that grid off a LINEAR interpolation of the step
     * grid is a first-order error the dense output does not have. Null leaves
     * every other caller unchanged.</p>
     */
    public double[] tranpoints = null;

    /**
     * Fixed timestep for transient analysis (null for adaptive stepping)
     */
    public Double timestep;

    /**
     * Verbosity level for solver output
     */
    public VerboseLevel verbose;

    /*
     * The solver console (jline.io.LineConsole) has NO option of its own: it IS
     * VerboseLevel.DEBUG, so `verbose = VerboseLevel.DEBUG` is what asks for it.
     */

    /**
     * Wall-clock time budget in seconds from solver launch to end.
     * When the elapsed wall-clock time exceeds this value, the solver stops at
     * the next cooperative checkpoint and returns either an interim solution (if
     * one is available) or an empty result with a warning. Default is
     * Double.POSITIVE_INFINITY (no time budget). This is wall-clock time and
     * must not be confused with timespan, which is simulated/model time.
     */
    public double timeout;

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
     * Load-dependent service rate scaling factors for QRF methods.
     * Dimensions: alpha[i][n] where i is station index and n is population level.
     * When null, defaults to ones(M, N) in the QRF solver.
     */
    public double[][] qrfAlpha;

    /**
     * Blocking parameters for the QRF bounds qrf.bas / qrf.rsrd. Null means the
     * caller supplied none, and those two methods then REJECT the model rather
     * than assume no blocking; see {@link QrfParams}.
     */
    public QrfParams qrfParams;

    /**
     * Hierarchy level for SolverBA bound hierarchies (pbh/cbh MVA-step depth,
     * pbk/bjbk iteration count k). Default 2.
     */
    public int level = 2;

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
        this.config.fj_warmstart = true;
        this.config.map_env = "auto";
        this.config.map_env_method = "auto";
        this.config.map_env_maxstages = 64;
        this.config.eventcache = true;
        this.config.hide_immediate = false; // Default for config
        this.config.pstar = new ArrayList<>();
        this.config.variates = "none"; // Variance reduction disabled by default
        this.config.nonmkv = "bernstein"; // Method for non-Markovian distribution conversion
        this.config.nonmkvorder = 20; // Order (number of phases) for non-Markovian approximation
        this.config.phfit = "cme"; // Surrogate family for concrete distributions
        this.config.da = "courtois"; // CTMC decomposition/aggregation method for Env solver
        this.config.da_iter = 10; // Number of iterations for kms/takahashi
        this.config.relax = "none"; // Default is no relaxation (will be overridden for LN solver)
        this.config.relax_factor = 0.1; // Relaxation factor when enabled
        this.config.relax_min = 0.1; // Minimum relaxation factor for adaptive mode
        this.config.relax_history = 5; // Error history window for adaptive mode
        // Stochastic iteration options (SolverLN), used when layer solvers are
        // simulation-based or Monte Carlo based and thus return noisy estimates
        this.config.stochiter = "auto"; // 'auto' | 'rm' | 'crn' | 'off'
        this.config.stochiter_alpha = 0.6; // Robbins-Monro step decay exponent, in (0.5,1]
        this.config.stochiter_a0 = 1.0; // Robbins-Monro initial step after burn-in
        this.config.stochiter_burnin = 5; // Picard burn-in iterations before step decay starts
        this.config.stochiter_conseq = 3; // consecutive sub-tolerance iterations required to stop
        this.force = false;
        this.hide_immediate = true; // Hide immediate transitions if possible
        this.init_sol = new Matrix(0, 0);
        this.iter_max = 1000;
        this.iter_tol = 0.0001; // Convergence tolerance to stop iterations
        this.tol = 0.0001; // Tolerance for all other uses
        this.keep = true;
        this.lang = "java";
        this.method = "default";
        this.remote = false;
        this.remote_endpoint = "127.0.0.1";
        this.restUrl = null;
        this.container = null;

        this.odesolvers = new ODESolvers();
        this.odesolvers.setDefaults(0.001, Inf, tol);

        this.samples = 10000;
        this.events = -1; // unset: DES solvers fall back to samples
        this.seed = RandomManager.generateRandomSeed();
        this.stiff = true;
        this.timespan = new double[2];
        this.timespan[0] = Inf;
        this.timespan[1] = Inf;
        this.timestep = null;
        this.verbose = GlobalConstants.getVerbose();
        this.confint = 0; // disabled by default
        this.timeout = Inf; // no wall-clock time budget by default
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
                // Finite per-class state-space cutoff for open/mixed models. Must stay
                // in step with MATLAB (SolverOptions.m, case 'CTMC') and native Python
                // (solver_ctmc.py), which both default to 10: the cutoff selects the
                // truncated chain, so an unequal default is an unequal model and the
                // three codebases disagree numerically. Leaving it null defers to the
                // ceil(6000^(1/(M*K))) auto-cutoff in runAnalyzer, which MATLAB reaches
                // only when the user explicitly asks for cutoff=Inf.
                this.cutoff(10.0);
                this.timespan = new double[]{Inf, Inf};
                // verbose inherits the global default set in this(); like MVA/NC/MAM
                // it must honor GlobalConstants.setVerbose(SILENT). MATLAB's CTMC case
                // likewise leaves the base STD default untouched.
                this.config.hide_immediate = true;
                this.config.state_space_gen = "full";
                this.rewardIterations = 1000; // Number of value iterations for reward computation
                break;
            case ENV:
                this.iter_max = 100;
                // STD like every other solver, as MATLAB and python now are:
                // at SILENT this ensemble printed no banner of its own and
                // could not narrate on the console either.
                this.verbose = VerboseLevel.STD;
                break;
            case FLUID:
                this.config.highvar = "default";
                // TRUE: a coordinate whose exit rate is GlobalConstants.Immediate is LINE's
                // stand-in for infinity, not a fast rate, and every fluid route now
                // stochastic-complements it out of the event set rather than integrating it
                // (ImmediateElimination, reached through FluidHideImmediate). It read false while
                // MATLAB's stiff slot was @ode15s, which coped; under LSODA -- which this port has
                // always used -- a single LN task layer carrying one InfRate phase took 145 s here
                // for the answer the reduced system returns in 0.35 s.
                this.config.hide_immediate = true;
                // No immediate_tol default: the in-code threshold is GlobalConstants.Immediate*(1-1e-2),
                // the same rule MATLAB uses, so that only the InfRate sentinel qualifies and a
                // genuinely fast rate the user wrote does not. Set the key to override it.
                this.iter_max = 200;
                this.stiff = true;
                this.timespan[0] = 0;
                // Reduce min step to handle Immediate transition rates (~1e8)
                // MATLAB's ode15s has effectively no minimum step constraint
                this.odesolvers.setDefaults(1e-14, this.odesolvers.odemaxstep, tol);
                break;
            case LN:
                this.config.interlocking = true;
                this.config.layering = "srvn";
                this.config.multiserver = "default";
                // Under-relaxation options for convergence improvement
                this.config.relax = "fixed"; // 'auto' | 'fixed' | 'adaptive' | 'none'
                this.config.relax_factor = 0.5; // Relaxation factor (0 < omega <= 1)
                this.config.relax_min = 0.1; // Minimum relaxation factor for adaptive mode
                this.config.relax_history = 5; // Error history window for adaptive mode
                this.timespan = new double[]{Inf, Inf};
                this.keep = true;
                this.iter_max = 200; // More iterations for difficult LQN models
                this.iter_tol = 5e-3; // Convergence tolerance (looser than default for LQN models)
                this.tol = 1e-4;
                break;
            case LQNS:
                this.timespan = new double[]{Inf, Inf};
                this.keep = true;
                // STD like every other solver, matching MATLAB's flipped default
                this.verbose = VerboseLevel.STD;
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
                // verbose inherits the global default set in this() so that
                // GlobalConstants.setVerbose(SILENT) silences SSA like every other
                // solver; the base default is STD, matching MATLAB's "true".
                this.config.state_space_gen = "none";
                // eventcache is set based on lang in MATLAB, but Java always uses eventcache = true
                this.config.eventcache = true;
                // config.warmupfrac stays null (disabled), matching MATLAB's 0.0
                break;
            case JMT:
                // JMT removes its temporary working directory after each solve
                // unless the user opts in, matching MATLAB @SolverJMT/runAnalyzer.m
                // (options.keep defaults to false for JMT).
                this.keep = false;
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
            // The additionalParams map carries every config key set through
            // Config.put -- the ones with no declared field, such as the MAM
            // 'bgaggr'/'bgstates_max'/'qbdphases_max' and 'etaqa_trunc'. It was
            // NOT copied here, and Solver's own constructor calls copy(), so a
            // put() key never reached the analyzer: the solver silently ran its
            // default while reporting the requested method, which is
            // indistinguishable from the option having no effect at all.
            cloned.config.additionalParams = new java.util.HashMap<>(this.config.additionalParams);
            // NOTE: the rest of this clone is field-by-field, so a declared
            // Config field added without a line here is silently dropped on
            // every copy().
            cloned.config.slotted = this.config.slotted;
            cloned.config.slotlength = this.config.slotlength;
            cloned.config.highvar = this.config.highvar;
            cloned.config.fluid_earlystop = this.config.fluid_earlystop;
            cloned.config.warmupfrac = this.config.warmupfrac;
            cloned.config.multiserver = this.config.multiserver;
            cloned.config.aghq_nodes = this.config.aghq_nodes;
            cloned.config.mcmc_batches = this.config.mcmc_batches;
            cloned.config.mcmc_burnin = this.config.mcmc_burnin;
            cloned.config.np_priority = this.config.np_priority;
            cloned.config.pstar = this.config.pstar != null ? new ArrayList<>(this.config.pstar) : null;
            cloned.config.variates = this.config.variates;
            cloned.config.fork_join = this.config.fork_join;
            cloned.config.sjn_lattice_max = this.config.sjn_lattice_max;
            cloned.config.sjn_ns = this.config.sjn_ns;
            cloned.config.sjn_lfactor = this.config.sjn_lfactor;
            cloned.config.sjn_umax = this.config.sjn_umax;
            cloned.config.fj_warmstart = this.config.fj_warmstart;
            cloned.config.map_env = this.config.map_env;
            cloned.config.map_env_method = this.config.map_env_method;
            cloned.config.map_env_maxstages = this.config.map_env_maxstages;
            cloned.config.merge = this.config.merge;
            cloned.config.warmup = this.config.warmup;
            cloned.config.compress = this.config.compress;
            cloned.config.space_max = this.config.space_max;
            cloned.config.interlocking = this.config.interlocking;
            cloned.config.interlock = this.config.interlock != null ? this.config.interlock.copy() : null;
            cloned.config.interlock_chain = this.config.interlock_chain != null ? this.config.interlock_chain.copy() : null;
            cloned.config.layering = this.config.layering;
            cloned.config.layer_init = this.config.layer_init;
            cloned.config.eventcache = this.config.eventcache;
            cloned.config.hide_immediate = this.config.hide_immediate;
            cloned.config.state_space_gen = this.config.state_space_gen;
            cloned.config.nonmkv = this.config.nonmkv;
            cloned.config.nonmkvorder = this.config.nonmkvorder;
            cloned.config.phfit = this.config.phfit;
            cloned.config.da = this.config.da;
            cloned.config.da_iter = this.config.da_iter;
            cloned.config.relax = this.config.relax;
            cloned.config.relax_factor = this.config.relax_factor;
            cloned.config.relax_min = this.config.relax_min;
            cloned.config.relax_history = this.config.relax_history;
            cloned.config.stochiter = this.config.stochiter;
            cloned.config.stochiter_alpha = this.config.stochiter_alpha;
            cloned.config.stochiter_a0 = this.config.stochiter_a0;
            cloned.config.stochiter_burnin = this.config.stochiter_burnin;
            cloned.config.stochiter_conseq = this.config.stochiter_conseq;
            cloned.config.num_cdf_pts = this.config.num_cdf_pts;
            cloned.config.remote = this.config.remote;
            cloned.config.remote_url = this.config.remote_url;
            cloned.config.symbolic = this.config.symbolic;
            cloned.config.symbolic_timeout = this.config.symbolic_timeout;
            cloned.config.tbi_tol = this.config.tbi_tol;
            cloned.config.tbi_iter_max = this.config.tbi_iter_max;
            cloned.config.gmres_restart = this.config.gmres_restart;
            cloned.config.tbi_cellsize = this.config.tbi_cellsize;
            cloned.config.tbi_cells = this.config.tbi_cells;
            cloned.config.rate_traj_tgrid = this.config.rate_traj_tgrid != null
                    ? this.config.rate_traj_tgrid.clone() : null;
            cloned.config.rate_traj_mmat = this.config.rate_traj_mmat != null
                    ? this.config.rate_traj_mmat.copy() : null;
            cloned.config.nhpp_sched = this.config.nhpp_sched != null
                    ? new ArrayList<>(this.config.nhpp_sched) : null;
            cloned.config.rate_sched = this.config.rate_sched != null
                    ? new ArrayList<>(this.config.rate_sched) : null;
            cloned.config.ctmc_tv_ngrid = this.config.ctmc_tv_ngrid;
            cloned.config.transient_method = this.config.transient_method;
            cloned.config.fau_epsilon = this.config.fau_epsilon;
            cloned.config.fau_delta = this.config.fau_delta;
            cloned.config.fau_ngrid = this.config.fau_ngrid;
            cloned.config.ln_transient = this.config.ln_transient;
            cloned.config.ln_transient_channels = this.config.ln_transient_channels;
            cloned.config.ln_transient_iter_max = this.config.ln_transient_iter_max;
            cloned.config.ln_transient_tol = this.config.ln_transient_tol;
            cloned.config.orbit_maxlevel = this.config.orbit_maxlevel;
            cloned.config.orbit_tailtol = this.config.orbit_tailtol;
            cloned.config.kp_init_sol = this.config.kp_init_sol != null
                    ? this.config.kp_init_sol.clone() : null;
            cloned.config.init_cov = this.config.init_cov != null
                    ? new Matrix(this.config.init_cov) : null;
            cloned.config.moment_sigma2 = this.config.moment_sigma2 != null
                    ? this.config.moment_sigma2.clone() : null;
            cloned.config.moment_cov = this.config.moment_cov != null
                    ? this.config.moment_cov.clone() : null;
            cloned.config.moment_maxstate = this.config.moment_maxstate;
            cloned.config.dae_maxstate = this.config.dae_maxstate;
            cloned.config.dae_maxcov = this.config.dae_maxcov;
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
        cloned.level = this.level;
        cloned.remote = this.remote;
        cloned.remote_endpoint = this.remote_endpoint;
        cloned.restUrl = this.restUrl;
        cloned.container = this.container;

        // Fresh instances of the DEFAULT integrators (they hold mutable state), but a
        // caller-assigned integrator is kept: rebuilding it would discard the choice.
        if (this.odesolvers != null) {
            cloned.odesolvers = this.odesolvers.copy(this.tol);
        }

        cloned.samples = this.samples;
        cloned.events = this.events;
        cloned.seed = this.seed;
        cloned.stiff = this.stiff;
        cloned.timespan = this.timespan != null ? this.timespan.clone() : new double[]{Inf, Inf};
        cloned.tranpoints = this.tranpoints != null ? this.tranpoints.clone() : null;
        cloned.timestep = this.timestep;
        cloned.verbose = this.verbose;
        cloned.confint = this.confint;
        cloned.timeout = this.timeout;
        if (this.qrfAlpha != null) {
            cloned.qrfAlpha = new double[this.qrfAlpha.length][];
            for (int i = 0; i < this.qrfAlpha.length; i++) {
                cloned.qrfAlpha[i] = this.qrfAlpha[i].clone();
            }
        }
        if (this.qrfParams != null) {
            cloned.qrfParams = this.qrfParams.copy();
        }

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
     * Sets the event budget for discrete-event simulation solvers (builder pattern).
     * Overrides samples in SolverSSA (fired transitions) and SolverLDES
     * (service completion events).
     *
     * @param s number of events to simulate
     * @return this SolverOptions instance for method chaining
     */
    public SolverOptions events(int s) {
        this.events = s;
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
        this.odesolvers.setDefaults(this.odesolvers.odeminstep, odeMaxStep, tol);
    }

    /**
     * Sets the minimum step size for ODE solvers and updates all integrators.
     *
     * @param odeMinStep minimum step size for numerical integration
     */
    public void setODEMinStep(double odeMinStep) {
        this.odesolvers.setDefaults(odeMinStep, this.odesolvers.odemaxstep, tol);
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
     * Set the maximum number of internal steps for an LSODA solver instance
     * via reflection (the library does not expose a public setter).
     */
    public static void setLsodaMaxSteps(LSODA solver, int maxSteps) {
        try {
            java.lang.reflect.Field mxstepField = LSODA.class.getDeclaredField("mxstep");
            mxstepField.setAccessible(true);
            mxstepField.setInt(solver, maxSteps);
        } catch (Exception ignored) {
        }
    }

    /**
     * Advanced configuration options for specialized solver features.
     * These options control solver-specific behavior and may not be
     * applicable to all solver types.
     */
    public static class Config {

        /**
         * Run the solver on a discrete time scale (slot lattice) rather than on
         * the continuous time axis. SolverLDES simulates the lattice directly;
         * SolverNC routes to the discrete-time product-form analyzer
         * (Solver_nc_dt_analyzer) and refuses a model outside it, rather than
         * falling back to a continuous-time approximation. Default: false.
         */
        public boolean slotted = false;

        /**
         * Slot length in model time units, used when {@link #slotted} is true;
         * null keeps unit slots. Metrics are computed per slot and rescaled by
         * this factor, so the caller reads back the units it built the model in.
         */
        public Double slotlength;

        /**
         * High-variance class handling strategy
         */
        public String highvar;

        /**
         * SolverCTMC.getCdfFirstPassT: "expm" (default) or "lt", selecting the
         * dense matrix-exponential route or the Laplace-transform inversion of
         * ctmc_passage_time. Null means the default.
         */
        public String passage_method;

        /**
         * Fluid: stop the ODE window loop once the geometric tail of the
         * iteration says the state is at a fixed point, instead of always
         * running iter_max windows. Null means the default, which is on.
         */
        public Boolean fluid_earlystop;

        /**
         * SSA warmup discard fraction. When set and positive the serial SSA
         * engine discards this fraction of the collected trajectory from the
         * mean estimates (mirroring MATLAB solver_ssa.m) and the batch-means
         * confidence intervals use it as their transient discard; null or 0
         * disables the mean discard and the CI falls back to its legacy 10%.
         */
        public Double warmupfrac;

        /**
         * Warmup length as an ABSOLUTE number of samples to discard; null leaves
         * the transient filter on {@link #warmupfrac}.
         *
         * <p>The simulators only take a FRACTION of the run, so a count is
         * meaningful only against {@code options.samples} and is converted to
         * {@code warmupfrac = warmup / samples} when the run is configured, not
         * when it is parsed: the two may be given in either order. A count at
         * or beyond the sample budget would discard the whole run and is
         * refused rather than clamped.</p>
         */
        public Double warmup;

        /**
         * Convergence tolerance of the MDD level iteration (SolverCTMC method
         * 'mdd'); null keeps the Mdd_mcd default of 1e-12.
         *
         * <p>Kept separate from {@code iter_tol} on purpose: that is the
         * solver-level fixed-point tolerance, sized for AMVA outer loops, while
         * the level iteration is an inner numerical solve whose result is
         * checked against the population invariant at 1e-6. Feeding iter_tol
         * here stops the iteration short of the fixed point and trips that
         * guard.</p>
         */
        public Double mdd_tol;

        /**
         * Maximum sweeps of the MDD level iteration; null keeps the Mdd_mcd
         * default of 500.
         */
        public Integer mdd_maxiter;

        /**
         * Number of non-overlapping batches the Chen-O'Cinneide MCMC estimator
         * (Pfqn_mcmc) splits its run into for the batch-means confidence
         * intervals; null keeps the 30 of Schmeiser (1982), the count used in
         * the tables of the paper.
         */
        /**
         * options.config.aghq_nodes: nodes per simplex direction of the adaptive
         * Gauss-Hermite rule. q = 1 reproduces pfqn_le; the rule costs q^(M-1)
         * evaluations, so the default stays small.
         */
        public Integer aghq_nodes;
        public Integer mcmc_batches;

        /**
         * Warm-up fraction the MCMC estimator discards before it starts
         * accumulating; null keeps 0.1. The paper discards none and ignores the
         * initialization bias.
         */
        public Double mcmc_burnin;

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
         * RQNA (Robust Queueing Network Analyzer) near-immediate feedback
         * elimination toggle. Null (default) applies elimination (Whitt-You
         * Algorithm 2 / Section 4.2); set FALSE for plain RQNA (Algorithm 1).
         */
        public Boolean rqna_feedback_elim;

        /**
         * RQNA alpha correction terms toggle (eq. 34). Null (default) = enabled.
         */
        public Boolean rqna_alpha;

        /**
         * RQNA beta correction terms toggle (eqs. 38-39). Null (default) = enabled.
         */
        public Boolean rqna_beta;

        /**
         * RQT (Robust Queueing Theory) service adaptation regime of Table 1:
         * "independent" (default, service distribution unknown), "normal" or
         * "pareto".
         */
        public String rqt_regime;

        /**
         * RQT toggle replacing the closed-form bound of Theorem 3 by the exact
         * worst case over the uncertainty sets. Null (default) = closed form.
         */
        public Boolean rqt_exact;

        /**
         * RQT tail coefficient in (1,2] of the external arrival processes.
         * Null (default) = 2, the finite-variance regime.
         */
        public Double rqt_alpha_a;

        /**
         * RQT tail coefficient in (1,2] of the service processes.
         * Null (default) = 2, the finite-variance regime.
         */
        public Double rqt_alpha_s;

        /**
         * Variance reduction technique for LDES simulation.
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
         * Lattice size above which the shortest-job-next dispatch prefers the Schweitzer fixed
         * point over the population recursion, under method 'default'. Null keeps the built-in
         * threshold. Mirrors MATLAB/Python options.config.sjn_lattice_max.
         */
        public Double sjn_lattice_max;

        /**
         * Number of grid subdivisions of the job size axis at a shortest-job-next station, even.
         * Null keeps the built-in value. Mirrors MATLAB/Python options.config.sjn_ns.
         */
        public Integer sjn_ns;

        /**
         * Grid extent at a shortest-job-next station, in units of the largest mean service time.
         * Null keeps the built-in value. Mirrors MATLAB/Python options.config.sjn_lfactor.
         */
        public Double sjn_lfactor;

        /**
         * Utilization cap at a shortest-job-next station, strictly below one. Null keeps the
         * built-in value. Mirrors MATLAB/Python options.config.sjn_umax.
         */
        public Double sjn_umax;

        /**
         * Resume the fork-join (MMT) fixed point from the iterate retained by the
         * previous runAnalyzer call on the same solver, instead of restarting from
         * GlobalConstants.FineTol. Only has an effect under an outer iteration such
         * as SolverLN, which re-solves each layer once per outer iteration.
         */
        public boolean fj_warmstart;

        /**
         * Random-environment fallback for MAP/MMPP/MMAP models on solvers that
         * cannot consume a non-renewal process: "auto" approximates the model
         * through its environment image, "off" rejects it as before.
         */
        public String map_env;

        /**
         * Environment recombination used by that fallback: "auto", "meanfield"
         * (transient-capable stage solvers only), "dec" (slow-environment
         * limit) or "avg" (fast-environment limit).
         */
        public String map_env_method;

        /** Cap on the number of environment stages, i.e. the product of the phase orders. */
        public int map_env_maxstages;

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
         * Interlock matrix of Franks (1999), Eq. (4.7), CLASS-indexed: interlock(r,s) is the
         * share of the class-s queue that a class-r arrival must not see, because that work was
         * itself caused by the class-r request. Null for every model but the layers of SolverLN.
         */
        public Matrix interlock;

        /**
         * The same matrix aggregated to the chain basis, set by the AMVA handler for the
         * iteration it is about to run. Never set from outside the MVA solvers.
         */
        public Matrix interlock_chain;

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
         * Family used when a concrete distribution is replaced by a Markovian surrogate:
         * "cme" fits a concentrated matrix exponential plus an exponential tail, "ph"
         * keeps the Erlang/Bernstein phase-type. At a budget of nonmkvorder phases the ME
         * reaches an SCV of O(1/n^2) where the Erlang stops at 1/n, and it matches the
         * first two moments exactly. SSA, Fluid and JMT must use "ph".
         */
        public String phfit = "cme";
        /**
         * `options.config.runLengthPlan`: ask a simulation solver how long its
         * run SHOULD have been for a target relative precision.
         *
         * Null leaves the plan uncomputed. A Double is the target relative
         * precision at the run's own confidence level; a
         * {@code Map<String,Double>} with keys {@code relprecision} and
         * {@code confidence} sets both. The plan lands in
         * {@code result.runLengthPlan} and is computed by
         * {@link jline.api.sim.SimRunlength#sim_runlength_plan}.
         */
        public Object runLengthPlan = null;

        /**
         * Whether to preserve deterministic distributions during non-Markovian
         * approximation (true keeps Det, false approximates with PH).
         */
        public boolean preserveDet = false;

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
         * SolverLN layer initialization strategy. "bound" initializes the layer
         * throughputs from the Majumdar-Woodside robust box bounds before the
         * first solve; null/"none" uses the default (zero) initialization.
         */
        public String layer_init;

        /**
         * SolverLN layering strategy. "srvn" (default) builds one submodel per
         * server; "flat"/"squashed" builds a single submodel holding every
         * processor and task. See _kb/06-solver-catalog.md (LN section).
         */
        public String layering;

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

        /**
         * Stochastic iteration mode for SolverLN when one or more layer
         * solvers return noisy estimates (simulation or Monte Carlo).
         * <p>Available options:</p>
         * <ul>
         *   <li><b>"auto"</b>: enable Robbins-Monro mode when a stochastic layer solver is detected (default)</li>
         *   <li><b>"rm"</b>: Robbins-Monro step decay with Polyak-Ruppert averaging</li>
         *   <li><b>"crn"</b>: common random numbers (pinned per-layer seeds), standard convergence test</li>
         *   <li><b>"off"</b>: deterministic Picard iteration</li>
         * </ul>
         */
        public String stochiter;

        /**
         * Robbins-Monro step decay exponent alpha in (0.5,1].
         * <p>Default: 0.6</p>
         */
        public double stochiter_alpha;

        /**
         * Robbins-Monro initial step a0 applied after burn-in.
         * <p>Default: 1.0</p>
         */
        public double stochiter_a0;

        /**
         * Picard burn-in iterations before the Robbins-Monro step decay starts.
         * <p>Default: 5</p>
         */
        public int stochiter_burnin;

        /**
         * Consecutive sub-tolerance iterations of the averaged-iterate drift
         * required to stop.
         * <p>Default: 3</p>
         */
        public int stochiter_conseq;

        /**
         * Number of points for CDF (Cumulative Distribution Function) computation.
         * Used by SolverFluid and SolverMAM for response time distribution analysis.
         * <p>Default: 200 (100 for MAM solver)</p>
         */
        public int num_cdf_pts = 200;

        /**
         * AoI preemption/replacement probability override for Age of Information analysis.
         * Set to a value between 0 and 1 to override automatic detection.
         * NaN means auto-detect from scheduling strategy.
         * <p>Default: NaN (auto-detect)</p>
         */
        public double aoi_preemption = Double.NaN;

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
         * Computer algebra backend for symbolic analysis: "auto" searches for a
         * line-sage-rest service and starts one from a local image if needed,
         * "none" disables it, and any other value is either a service URL or a
         * Docker image name. See {@link jline.api.sym.SymEngines}.
         * <p>Default: "auto"</p>
         */
        public String symbolic = "auto";

        /**
         * Per-request timeout of the symbolic backend, in seconds. Symbolic
         * solves grow superpolynomially in the number of states, so this bounds
         * both the client wait and the server side computation.
         * <p>Default: 300</p>
         */
        public int symbolic_timeout = 300;

        /**
         * Trajectory-based iteration (TBI) fluid method: absolute sup-norm
         * tolerance on the state trajectory used to stop the Jacobi waveform
         * relaxation sweeps within a time segment.
         * <p>Default: 1e-3</p>
         */
        public double tbi_tol = 1e-3;

        /**
         * Trajectory-based iteration (TBI) fluid method: maximum number of
         * waveform relaxation sweeps per time segment.
         * <p>Default: 50</p>
         */
        public int tbi_iter_max = 50;

        /**
         * Krylov subspace dimension between restarts for the GMRES path in
         * ctmc_solve. Nonpositive leaves the kernel default of min(n, 50).
         * <p>Default: 0</p>
         */
        public int gmres_restart = 0;

        /**
         * Trajectory-based iteration (TBI) fluid method: target number of
         * stations per cell used by the greedy agglomerative partition when
         * tbi_cells is not provided.
         * <p>Default: 5</p>
         */
        public int tbi_cellsize = 5;

        /**
         * Trajectory-based iteration (TBI) fluid method: explicit station
         * partition. Each int[] is a cell of zero-based station indices; the
         * cells must be a disjoint cover of 0..nstations-1. Null (default)
         * triggers the greedy agglomerative partition.
         */
        public java.util.List<int[]> tbi_cells = null;

        /**
         * Fluid closing ODE: time grid of the caller-supplied event rate
         * multiplier trajectory. Paired with {@link #rate_traj_mmat}; both must
         * be set for the channel to be active. Mirrors the MATLAB
         * {@code options.config.rate_traj} cell {tgrid, Mmat}.
         * <p>Default: null (no time-varying multiplier)</p>
         */
        public double[] rate_traj_tgrid = null;

        /**
         * Fluid closing ODE: event rate multiplier matrix of the caller-supplied
         * trajectory, of size numEvents x numel(rate_traj_tgrid). Used by the
         * coupled LN layer transient to inject time-varying inter-layer demand.
         * <p>Default: null</p>
         */
        public Matrix rate_traj_mmat = null;

        /**
         * Fluid closing ODE: non-homogeneous Poisson (NHPP) source intensities
         * to track over the transient, one entry per (station, class). Injected
         * by {@code SolverFluid.getTranAvg}; steady-state analysis leaves this
         * null, the NHPP steady state being its time average. Mirrors the
         * MATLAB {@code options.config.nhpp_sched}.
         * <p>Default: null</p>
         */
        public java.util.List<FluidRateMultiplier.NhppEntry> nhpp_sched = null;

        /**
         * Fluid closing ODE: explicit per-(station,class) rate trajectories.
         * Mirrors the MATLAB {@code options.config.rate_sched}.
         * <p>Default: null</p>
         */
        public java.util.List<FluidRateMultiplier.RateEntry> rate_sched = null;

        /**
         * Fluid 'kp': initial state vector in the KO-PENDER layout -- one offset
         * counter walking the stations in order, an arrival-phase block at each EXT
         * station-class and a service-phase block at every other, with no mass
         * returning to the source. Mirrors the MATLAB
         * {@code options.config.kp_init_sol}.
         *
         * <p>NOT {@link SolverOptions#init_sol}: that one is laid out for the
         * CLOSING state vector, and consuming it here silently zeroes the source
         * phase mass and with it the whole network. A wrong-sized seed is
         * refused rather than ignored.</p>
         * <p>Default: null</p>
         */
        public double[] kp_init_sol = null;

        /**
         * Fluid 'kp': initial covariance Sigma(0), dim-by-dim in the same layout
         * as {@link #kp_init_sol}. A caller that carries a DISTRIBUTION across a
         * handoff supplies the second moment beside the mean, so the next stage
         * does not restart from a point mass it never had. Null keeps the
         * default diag(theta) - theta theta' of the initial arrival phase.
         * <p>Default: null</p>
         */
        public Matrix init_cov = null;

        /**
         * Fluid moment closure: per-station population variance closing the
         * {@code E[min(X,c)]} term of the closing drift. Null or all-zero
         * selects the first-order closure, which is what every method other than
         * {@code minnormal} uses, so the legacy code path stays bit-identical.
         * Set by the outer fixed point of the moment-closure driver.
         * <p>Default: null</p>
         */
        public double[] moment_sigma2 = null;

        /**
         * Fluid moment closure: per-station coordinate covariance blocks closing
         * the capacity-share RATIO of the PS/FCFS/DPS/GPS branches, which is a
         * separate closure from the min(). One entry per station, null where the
         * station needs none. Paired with {@link #moment_sigma2}.
         * <p>Default: null</p>
         */
        public Matrix[] moment_cov = null;

        /**
         * Fluid moment closure: largest phase-resolved state the covariance
         * (Lyapunov) equation is attempted on. The solve is cubic in this
         * dimension and the covariance is dense, so the driver refuses rather
         * than silently crawling above it.
         * <p>Default: 200</p>
         */
        public int moment_maxstate = 200;

        /**
         * Fluid 'dae': largest phase-resolved state the SIMULTANEOUS closure
         * solve is attempted on. Lower than {@link #moment_maxstate} because
         * this route takes a finite-difference Jacobian over the unknowns rather
         * than solving one Lyapunov equation, so its cost is quartic and not
         * cubic. Above it the driver refuses and names
         * {@code options.method='minnormal'}, which computes the same closure by
         * successive substitution.
         * <p>Default: 100</p>
         */
        public int dae_maxstate = 100;

        /**
         * Fluid 'dae' transient: largest covariance dimension integrated
         * ALONGSIDE the mean. The covariance adds nc^2 differential states and
         * the Jacobian is formed by finite differences over all of them, so the
         * cost grows as nc^4. Above this the mean is still integrated as a DAE
         * -- population conservation stays an algebraic equation -- but the
         * variance is held at its stationary value, which is what 'minnormal'
         * does for the whole of its transient anyway.
         * <p>Default: 25</p>
         */
        public int dae_maxcov = 25;

        /**
         * CTMC time-varying transient: number of points of the uniform time grid
         * over which the non-homogeneous generator is propagated (one matrix
         * exponential per interval, multiplier frozen at the interval midpoint).
         * Only read when {@link #rate_sched} is set. Mirrors the MATLAB
         * {@code options.config.ctmc_tv_ngrid}.
         * <p>Default: 100</p>
         */
        public int ctmc_tv_ngrid = 100;

        /**
         * CTMC transient method: {@code "ode"} (default) integrates the forward
         * equation, {@code "fau"} marches fast adaptive uniformization
         * ({@link jline.api.mc.Ctmc_fau}) over the output grid. It is a config
         * key rather than a {@code method} name because it changes no stationary
         * answer, and because a new entry in the solver's valid-method list is
         * enumerated by the sanity harness, which then wants a baseline per
         * method. Mirrors the MATLAB {@code options.config.transient_method}.
         * <p>Default: "ode"</p>
         */
        public String transient_method = "ode";

        /**
         * CTMC transient by {@code "fau"}: total probability mass the whole
         * horizon may discard. It is divided by the number of grid steps, each
         * step removing mass and none putting any back, so the accumulated
         * defect stays below this. Mirrors {@code options.config.fau_epsilon}.
         * <p>Default: 1e-6</p>
         */
        public double fau_epsilon = 1e-6;

        /**
         * CTMC transient by {@code "fau"}: occupancy below which a state is
         * dropped from the support. Mirrors {@code options.config.fau_delta}.
         * <p>Default: 1e-12</p>
         */
        public double fau_delta = 1e-12;

        /**
         * CTMC transient by {@code "fau"}: number of points of the uniform
         * output grid, used when {@code options.timestep} is unset. Mirrors
         * {@code options.config.fau_ngrid}.
         * <p>Default: 100</p>
         */
        public int fau_ngrid = 100;

        /**
         * SolverLN transient coupling mode: {@code "coupled"} (default) runs the
         * waveform-relaxation coupled layered transient, {@code "decoupled"}
         * freezes the inter-layer demands at the converged fixed point. Mirrors
         * the MATLAB {@code options.config.ln_transient}.
         */
        public String ln_transient = "coupled";

        /**
         * SolverLN coupled transient: which inter-layer coupling channels are
         * injected, {@code "both"} (default), {@code "thinkt"} (client-delay
         * only) or {@code "callservt"} (synchronous-call service only). Mirrors
         * the MATLAB {@code options.config.ln_transient_channels}.
         */
        public String ln_transient_channels = "both";

        /**
         * SolverLN coupled transient: maximum number of waveform-relaxation
         * iterations. Mirrors {@code options.config.ln_transient_iter_max}.
         * <p>Default: 20</p>
         */
        public int ln_transient_iter_max = 20;

        /**
         * SolverLN coupled transient: sup-norm tolerance on the queue-length
         * trajectory change between relaxation iterations. Mirrors
         * {@code options.config.ln_transient_tol}.
         * <p>Default: 1e-2</p>
         */
        public double ln_transient_tol = 1e-2;

        /**
         * Fixed orbit truncation level for the matrix-analytic retrial solver.
         * Non-positive (default) selects the adaptive, residual-driven level.
         * <p>Default: 0</p>
         */
        public int orbit_maxlevel = 0;

        /**
         * Relative orbit-truncation error target of the adaptive refinement in
         * the matrix-analytic retrial solver.
         * <p>Default: 1e-6</p>
         */
        public double orbit_tailtol = 1e-6;

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

        // The instances handed out by the last setDefaults, so copy() can tell a
        // default apart from one the caller assigned and only rebuild the former.
        private FirstOrderIntegrator defaultFastODESolver;
        private FirstOrderIntegrator defaultAccurateODESolver;
        private LSODA defaultFastStiffODESolver;
        private LSODA defaultAccurateStiffODESolver;

        /**
         * Installs the built-in integrators, mirroring MATLAB {@code SolverOptions.m}:
         * {@code fastOdeSolver}/{@code accurateOdeSolver} are the non-stiff pair
         * (@ode23 / @ode113 there, their adaptive Runge-Kutta and Adams-Moulton
         * counterparts here) and the stiff pair is LSODA, as @ode15s is there.
         *
         * @param minStep minimum step size
         * @param maxStep maximum step size
         * @param tol     absolute and relative tolerance
         */
        public void setDefaults(double minStep, double maxStep, double tol) {
            this.odeminstep = minStep;
            this.odemaxstep = maxStep;
            // @ode23 twin: the same Bogacki-Shampine 3(2) pair, supplied here
            // because commons-math3 ships no order-3 embedded Runge-Kutta method
            this.fastODESolver = new BogackiShampine23Integrator(minStep, maxStep, tol, tol);
            // @ode113 twin: Adams-Moulton PECE. ode113 varies its order up to 13;
            // the commons-math3 implementation is fixed-order, and order 5 is where
            // its Nordsieck start-up stays stable at the tolerances used here.
            this.accurateODESolver = new AdamsMoultonIntegrator(5, minStep, maxStep, tol, tol);
            this.fastStiffODESolver = newStiff(minStep, maxStep, tol, 3, 3);
            this.accurateStiffODESolver = newStiff(minStep, maxStep, tol, 12, 5);
            this.defaultFastODESolver = this.fastODESolver;
            this.defaultAccurateODESolver = this.accurateODESolver;
            this.defaultFastStiffODESolver = this.fastStiffODESolver;
            this.defaultAccurateStiffODESolver = this.accurateStiffODESolver;
        }

        /**
         * The stiff integrator to use over {@code [t0,t1]}.
         *
         * <p>MATLAB's {@code odeset} leaves {@code MaxStep} at {@code 0.1*|tf-t0|}
         * unless the caller pins it, and every MATLAB fluid ODE call but
         * {@code solver_fluid_kp} takes that default. The JAR left the step
         * unbounded, and an unbounded step lets LSODA's Nordsieck interpolation
         * drift off a fixed point by ~1e-6 instead of settling on it, so a closed
         * model's throughputs stopped balancing flow. The bound is applied only to
         * the built-in integrators and only while {@code odemaxstep} is infinite,
         * so an integrator or a maximum step the caller supplied is left alone.</p>
         *
         * @param t0   window start
         * @param t1   window end
         * @param tol  absolute and relative tolerance
         * @param fast true for the coarse-tolerance integrator
         * @return the integrator to drive over this window
         */
        public LSODA stiffIntegratorFor(double t0, double t1, double tol, boolean fast) {
            LSODA chosen = fast ? this.fastStiffODESolver : this.accurateStiffODESolver;
            LSODA builtin = fast ? this.defaultFastStiffODESolver : this.defaultAccurateStiffODESolver;
            double maxStep = defaultMaxStep(t0, t1);
            if (chosen != builtin || maxStep <= 0) {
                return chosen;
            }
            return fast ? newStiff(this.odeminstep, maxStep, tol, 3, 3)
                    : newStiff(this.odeminstep, maxStep, tol, 12, 5);
        }

        /**
         * The stiff slot's integrator, pinned to the BDF half.
         *
         * <p>MATLAB fills this slot with {@code @ode15s} and the native Python with
         * scipy's BDF, both A-stable. LSODA is an Adams/BDF auto-switcher that STARTS
         * on Adams, and its Adams half loses its stability bound at a fixed point,
         * where the corrector converges before the eigenvalue estimate behind that
         * bound is taken: the step then grows to MaxStep and the solution wanders
         * around the fixed point at the amplitude the error test tolerates instead of
         * settling on it. On Delay(1) -> PS(0.8), N=4 at tol=1e-4 that is 2.5e-5 of
         * queue length after 200 windows, against 1e-12 in both references, and it
         * grows with the horizon. Pinning the BDF half makes this slot the ode15s
         * twin it is documented to be; it also took FEWER steps on that model.</p>
         */
        private LSODA newStiff(double minStep, double maxStep, double tol, int maxOrderN, int maxOrderS) {
            LSODAExt s = new LSODAExt(minStep, maxStep, tol, tol, maxOrderN, maxOrderS, 10000000);
            s.setForceStiff(true);
            return s;
        }

        /**
         * The non-stiff integrator to use over {@code [t0,t1]}, under the same
         * {@code odeset} MaxStep default as {@link #stiffIntegratorFor}.
         *
         * @param t0   window start
         * @param t1   window end
         * @param tol  absolute and relative tolerance
         * @param fast true for the coarse-tolerance integrator
         * @return the integrator to drive over this window
         */
        public FirstOrderIntegrator integratorFor(double t0, double t1, double tol, boolean fast) {
            FirstOrderIntegrator chosen = fast ? this.fastODESolver : this.accurateODESolver;
            FirstOrderIntegrator builtin = fast ? this.defaultFastODESolver : this.defaultAccurateODESolver;
            double maxStep = defaultMaxStep(t0, t1);
            if (chosen != builtin || maxStep <= 0) {
                return chosen;
            }
            return fast ? new BogackiShampine23Integrator(this.odeminstep, maxStep, tol, tol)
                    : new AdamsMoultonIntegrator(5, this.odeminstep, maxStep, tol, tol);
        }

        /** MATLAB odeset's MaxStep default, or 0 when the caller already pinned one. */
        private double defaultMaxStep(double t0, double t1) {
            if (!Double.isInfinite(this.odemaxstep)) {
                return 0;
            }
            double width = Math.abs(t1 - t0);
            return (width > 0 && !Double.isInfinite(width)) ? 0.1 * width : 0;
        }

        /**
         * Returns a copy for a cloned SolverOptions. Integrators carry mutable state
         * (step handlers, and LSODA's whole working set), so each copy gets its own
         * instance of every DEFAULT integrator -- but an integrator the caller
         * assigned is carried over as-is, because there is no way to rebuild it and
         * silently replacing it would discard the caller's choice.
         *
         * @param tol tolerance for the rebuilt default integrators
         * @return the copy
         */
        public ODESolvers copy(double tol) {
            ODESolvers cloned = new ODESolvers();
            cloned.setDefaults(this.odeminstep, this.odemaxstep, tol);
            if (this.fastODESolver != this.defaultFastODESolver) {
                cloned.fastODESolver = this.fastODESolver;
            }
            if (this.accurateODESolver != this.defaultAccurateODESolver) {
                cloned.accurateODESolver = this.accurateODESolver;
            }
            if (this.fastStiffODESolver != this.defaultFastStiffODESolver) {
                cloned.fastStiffODESolver = this.fastStiffODESolver;
            }
            if (this.accurateStiffODESolver != this.defaultAccurateStiffODESolver) {
                cloned.accurateStiffODESolver = this.accurateStiffODESolver;
            }
            return cloned;
        }
    }
}
