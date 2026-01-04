/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.des;

import jline.VerboseLevel;
import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;

/**
 * Configuration options for Discrete Event Simulation (DES) solver.
 *
 * <p>DESOptions extends the base solver options to provide DES-specific
 * configuration parameters. DES uses the SSJ library for discrete-event simulation
 * that generates sample paths of the queueing network state evolution over time.</p>
 *
 * <p>Key DES characteristics:
 * <ul>
 *   <li>Discrete event simulation approach using SSJ library</li>
 *   <li>Handles multiclass Jackson queueing networks</li>
 *   <li>Supports steady-state analysis</li>
 *   <li>Provides statistical estimates with confidence intervals</li>
 *   <li>Uses event-count based stopping (samples = max service completions)</li>
 *   <li>Default: 200,000 service completion events</li>
 * </ul>
 * </p>
 *
 * <p><b>Default Configuration (based on benchmark analysis):</b>
 * <pre>
 * samples     = 200,000   (good accuracy/speed balance, ~5% mean error)
 * tranfilter  = "mser5"   (adaptive warmup detection)
 * warmupfrac  = 0.20      (20% warmup for fixed filter)
 * cimethod    = "obm"     (overlapping batch means)
 * obmoverlap  = 0.50      (50% overlap)
 * cnvgon      = false     (fixed event count)
 * </pre>
 * For heavy-tailed workloads (high SCV), consider increasing samples to 500,000+.
 * </p>
 *
 * @see SolverDES
 * @see SolverOptions
 * @since 1.0
 */
public class DESOptions extends SolverOptions {

    /** Default number of service completion events for DES solver (200,000) */
    public static final int DEFAULT_SAMPLES = 200000;

    /** Default convergence tolerance (5% relative precision) */
    public static final double DEFAULT_CNVG_TOL = 0.05;

    /** Default minimum number of batches before checking convergence */
    public static final int DEFAULT_CNVG_BATCH = 20;

    /**
     * Enable convergence-based stopping.
     * When enabled, simulation stops when confidence interval half-width relative
     * to mean falls below the convergence tolerance for all metrics.
     * Default: false (disabled to maintain backward compatibility with existing tests)
     */
    public boolean cnvgon = false;

    /**
     * Convergence tolerance - stop when (CI half-width / mean) is below this threshold.
     * A value of 0.05 means 5% relative precision.
     * Default: 0.05
     */
    public double cnvgtol = DEFAULT_CNVG_TOL;

    /**
     * Minimum number of batches required before checking for convergence.
     * More batches provide more reliable confidence interval estimates.
     * Default: 20
     */
    public int cnvgbatch = DEFAULT_CNVG_BATCH;

    /**
     * Number of events between convergence checks.
     * A value of 0 means auto-calculate (samples / 50).
     * Default: 0 (auto)
     */
    public int cnvgchk = 0;

    // ==================== Transient Detection Options ====================

    /** Default MSER batch size (MSER-5 standard) */
    public static final int DEFAULT_MSER_BATCH = 5;

    /** Default warmup fraction when using fixed method */
    public static final double DEFAULT_WARMUP_FRAC = 0.2;

    /**
     * Transient filter method for warmup detection and removal.
     * Valid values: "mser5" (automatic MSER-5), "fixed" (fixed fraction), "none" (no filtering).
     * Default: "mser5"
     */
    public String tranfilter = "mser5";

    /**
     * Batch size for MSER transient detection algorithm.
     * Only used when tranfilter = "mser5".
     * Default: 5 (standard MSER-5)
     */
    public int mserbatch = DEFAULT_MSER_BATCH;

    /**
     * Warmup fraction when using fixed transient filter.
     * Fraction of total events discarded as warmup (0.0 to 1.0).
     * Only used when tranfilter = "fixed".
     * Default: 0.2 (20% warmup)
     */
    public double warmupfrac = DEFAULT_WARMUP_FRAC;

    // ==================== Confidence Interval Options ====================

    /** Default OBM overlap fraction (50% overlap) */
    public static final double DEFAULT_OBM_OVERLAP = 0.5;

    /** Default minimum batch size for CI computation */
    public static final int DEFAULT_CI_MIN_BATCH = 10;

    /** Default minimum observations required for CI computation */
    public static final int DEFAULT_CI_MIN_OBS = 100;

    /**
     * Confidence interval computation method.
     * Valid values: "obm" (overlapping batch means), "bm" (batch means), "none" (no CI).
     * Default: "obm"
     */
    public String cimethod = "obm";

    /**
     * Overlap fraction for OBM confidence intervals (0.0 to 1.0).
     * 0.5 = 50% overlap (standard), 0.0 = non-overlapping (same as bm).
     * Only used when cimethod = "obm".
     * Default: 0.5
     */
    public double obmoverlap = DEFAULT_OBM_OVERLAP;

    /**
     * Minimum batch size for confidence interval computation.
     * Actual batch size is max(ciminbatch, sqrt(n)) where n is sample count.
     * Default: 10
     */
    public int ciminbatch = DEFAULT_CI_MIN_BATCH;

    /**
     * Minimum number of post-warmup observations required before computing CIs.
     * If fewer observations available, CI is set to 0 (unavailable).
     * Default: 100
     */
    public int ciminobs = DEFAULT_CI_MIN_OBS;

    /**
     * Enable detailed BAS (Blocking After Service) tracing for debugging.
     * Default: false
     */
    public boolean basTrace = false;

    // ==================== Parallel Replication Options ====================

    /** Default number of replications (1 = single run, backward compatible) */
    public static final int DEFAULT_REPLICATIONS = 1;

    /**
     * Number of independent replications to run.
     * When > 1, enables parallel replication mode where each replication
     * runs independently with its own RNG stream, and results are aggregated
     * using cross-replication variance for confidence intervals.
     * Default: 1 (single long run with batch means, backward compatible)
     */
    public int replications = DEFAULT_REPLICATIONS;

    /**
     * Number of parallel threads for replication execution.
     * Default: Runtime.getRuntime().availableProcessors() / 2
     * Only used when replications > 1.
     */
    public int numThreads = (int) Math.ceil(Runtime.getRuntime().availableProcessors() / 2.0);

    /**
     * Constructs a new DESOptions instance with default DES solver configuration.
     * Sets the maximum service completion events to 200,000 by default.
     * Verbose output is disabled by default for DES.
     */
    public DESOptions() {
        super(SolverType.DES);
        this.samples = DEFAULT_SAMPLES;
        this.verbose = VerboseLevel.SILENT;
    }

    /**
     * Sets whether convergence-based stopping is enabled.
     * @param enabled true to enable convergence checking
     * @return this options instance for method chaining
     */
    public DESOptions setCnvgon(boolean enabled) {
        this.cnvgon = enabled;
        return this;
    }

    /**
     * Sets the convergence tolerance.
     * @param tolerance the relative precision threshold (e.g., 0.05 for 5%)
     * @return this options instance for method chaining
     */
    public DESOptions setCnvgtol(double tolerance) {
        this.cnvgtol = tolerance;
        return this;
    }

    /**
     * Sets the minimum number of batches before checking convergence.
     * @param minBatches the minimum number of batches required
     * @return this options instance for method chaining
     */
    public DESOptions setCnvgbatch(int minBatches) {
        this.cnvgbatch = minBatches;
        return this;
    }

    /**
     * Sets the number of events between convergence checks.
     * @param interval the number of events (0 for auto-calculation)
     * @return this options instance for method chaining
     */
    public DESOptions setCnvgchk(int interval) {
        this.cnvgchk = interval;
        return this;
    }

    // ==================== Transient Detection Setters ====================

    /**
     * Sets the transient filter method.
     * @param method the transient filter method ("mser5", "fixed", or "none")
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if method is not valid
     */
    public DESOptions setTranfilter(String method) {
        if (!"mser5".equals(method) && !"fixed".equals(method) && !"none".equals(method)) {
            throw new IllegalArgumentException("tranfilter must be 'mser5', 'fixed', or 'none'");
        }
        this.tranfilter = method;
        return this;
    }

    /**
     * Sets the MSER batch size.
     * @param batchSize the batch size for MSER algorithm (must be >= 1)
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if batchSize is less than 1
     */
    public DESOptions setMserbatch(int batchSize) {
        if (batchSize < 1) {
            throw new IllegalArgumentException("mserbatch must be at least 1");
        }
        this.mserbatch = batchSize;
        return this;
    }

    /**
     * Sets the warmup fraction for fixed transient filter.
     * @param fraction the warmup fraction (0.0 to 1.0, exclusive)
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if fraction is not in [0, 1)
     */
    public DESOptions setWarmupfrac(double fraction) {
        if (fraction < 0.0 || fraction >= 1.0) {
            throw new IllegalArgumentException("warmupfrac must be in [0, 1)");
        }
        this.warmupfrac = fraction;
        return this;
    }

    // ==================== Confidence Interval Setters ====================

    /**
     * Sets the confidence interval computation method.
     * @param method the CI method ("obm", "bm", or "none")
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if method is not valid
     */
    public DESOptions setCimethod(String method) {
        if (!"obm".equals(method) && !"bm".equals(method) && !"none".equals(method)) {
            throw new IllegalArgumentException("cimethod must be 'obm', 'bm', or 'none'");
        }
        this.cimethod = method;
        return this;
    }

    /**
     * Sets the OBM overlap fraction.
     * @param overlap the overlap fraction (0.0 to 1.0)
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if overlap is not in [0, 1]
     */
    public DESOptions setObmoverlap(double overlap) {
        if (overlap < 0.0 || overlap > 1.0) {
            throw new IllegalArgumentException("obmoverlap must be in [0, 1]");
        }
        this.obmoverlap = overlap;
        return this;
    }

    /**
     * Sets the minimum batch size for CI computation.
     * @param minSize the minimum batch size (must be >= 2)
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if minSize is less than 2
     */
    public DESOptions setCiminbatch(int minSize) {
        if (minSize < 2) {
            throw new IllegalArgumentException("ciminbatch must be at least 2");
        }
        this.ciminbatch = minSize;
        return this;
    }

    /**
     * Sets the minimum observations required for CI computation.
     * @param minObs the minimum number of observations (must be >= 10)
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if minObs is less than 10
     */
    public DESOptions setCiminobs(int minObs) {
        if (minObs < 10) {
            throw new IllegalArgumentException("ciminobs must be at least 10");
        }
        this.ciminobs = minObs;
        return this;
    }

    // ==================== Parallel Replication Setters ====================

    /**
     * Sets the number of independent replications to run.
     * When > 1, enables parallel replication mode with cross-replication CI.
     * @param numReps number of replications (must be >= 1)
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if numReps is less than 1
     */
    public DESOptions setReplications(int numReps) {
        if (numReps < 1) {
            throw new IllegalArgumentException("replications must be at least 1");
        }
        this.replications = numReps;
        return this;
    }

    /**
     * Sets the number of parallel threads for replication execution.
     * @param threads number of threads (must be >= 1)
     * @return this options instance for method chaining
     * @throws IllegalArgumentException if threads is less than 1
     */
    public DESOptions setNumThreads(int threads) {
        if (threads < 1) {
            throw new IllegalArgumentException("numThreads must be at least 1");
        }
        this.numThreads = threads;
        return this;
    }

    /**
     * Creates a deep copy of this DESOptions object.
     * Overrides the base class copy() to preserve DES-specific fields.
     *
     * @return a deep copy of this DESOptions object
     */
    @Override
    public DESOptions copy() {
        DESOptions cloned = new DESOptions();

        // Copy base class fields
        cloned.cache = this.cache;
        cloned.cutoff = this.cutoff != null ? this.cutoff.copy() : null;

        if (this.config != null) {
            cloned.config = new Config();
            cloned.config.highvar = this.config.highvar;
            cloned.config.multiserver = this.config.multiserver;
            cloned.config.np_priority = this.config.np_priority;
            cloned.config.pstar = this.config.pstar != null ? new java.util.ArrayList<>(this.config.pstar) : null;
            cloned.config.fork_join = this.config.fork_join;
            cloned.config.merge = this.config.merge;
            cloned.config.compress = this.config.compress;
            cloned.config.space_max = this.config.space_max;
            cloned.config.interlocking = this.config.interlocking;
            cloned.config.eventcache = this.config.eventcache;
            cloned.config.hide_immediate = this.config.hide_immediate;
            cloned.config.state_space_gen = this.config.state_space_gen;
        }

        cloned.force = this.force;
        cloned.hide_immediate = this.hide_immediate;
        cloned.init_sol = this.init_sol != null ? this.init_sol : new jline.util.matrix.Matrix(0, 0);
        cloned.iter_max = this.iter_max;
        cloned.iter_tol = this.iter_tol;
        cloned.tol = this.tol;
        cloned.keep = this.keep;
        cloned.lang = this.lang;
        cloned.method = this.method;
        cloned.remote = this.remote;
        cloned.remote_endpoint = this.remote_endpoint;

        if (this.odesolvers != null) {
            cloned.odesolvers = new ODESolvers();
            cloned.odesolvers.odeminstep = this.odesolvers.odeminstep;
            cloned.odesolvers.odemaxstep = this.odesolvers.odemaxstep;
            cloned.odesolvers.fastODESolver = this.odesolvers.fastODESolver;
            cloned.odesolvers.accurateODESolver = this.odesolvers.accurateODESolver;
            cloned.odesolvers.fastStiffODESolver = this.odesolvers.fastStiffODESolver;
            cloned.odesolvers.accurateStiffODESolver = this.odesolvers.accurateStiffODESolver;
        }

        cloned.samples = this.samples;
        cloned.seed = this.seed;
        cloned.stiff = this.stiff;
        cloned.timespan = this.timespan != null ? this.timespan.clone() : new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        cloned.timestep = this.timestep;
        cloned.verbose = this.verbose;
        cloned.confint = this.confint;

        // Copy DES-specific fields
        cloned.cnvgon = this.cnvgon;
        cloned.cnvgtol = this.cnvgtol;
        cloned.cnvgbatch = this.cnvgbatch;
        cloned.cnvgchk = this.cnvgchk;

        // Copy transient detection options
        cloned.tranfilter = this.tranfilter;
        cloned.mserbatch = this.mserbatch;
        cloned.warmupfrac = this.warmupfrac;

        // Copy confidence interval options
        cloned.cimethod = this.cimethod;
        cloned.obmoverlap = this.obmoverlap;
        cloned.ciminbatch = this.ciminbatch;
        cloned.ciminobs = this.ciminobs;

        // Copy parallel replication options
        cloned.replications = this.replications;
        cloned.numThreads = this.numThreads;

        return cloned;
    }
}
