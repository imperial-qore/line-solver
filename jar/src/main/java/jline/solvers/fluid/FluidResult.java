/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid;

import jline.api.aoi.AoiMfqResult;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * Result container for Fluid solver analysis of queueing networks.
 *
 * <p>FluidResult extends the base SolverResult to include fluid-specific
 * performance metrics and analysis results. Fluid solvers use continuous
 * approximations and differential equations to analyze queueing networks
 * with general service and arrival processes.</p>
 *
 * <p>Fluid analysis results include:
 * <ul>
 *   <li>ODE-based steady-state and transient solutions</li>
 *   <li>Phase-type distribution approximations</li>
 *   <li>Response time and passage time distributions</li>
 *   <li>Aggregated state probabilities</li>
 *   <li>Age of Information (AoI) metrics and distributions</li>
 * </ul>
 * </p>
 *
 * @see jline.solvers.fluid.SolverFluid
 * @see SolverResult
 * @since 1.0
 */
public class FluidResult extends SolverResult {

    /** ODE state vector solution from fluid approximation */
    public Matrix odeStateVec;

    /** Modified network structure after phase-type expansion (if applied) */
    public jline.lang.NetworkStruct snFinal;

    /** Distribution matrices for response time and passage time CDFs [stations][classes] */
    public Matrix[][] distribC;

    /** Runtime required for distribution computation */
    public double distribRuntime;

    /** Aggregate state probability */
    public double Pnir;

    /** Logarithm of aggregate state probability */
    public double logPnir;

    /** Age of Information results from MFQ analysis (null if not an AoI model) */
    public AoiMfqResult aoiResults;

    /**
     * Transient queue-length VARIANCE per station and class, [stations][classes], each a
     * column over the time grid in {@code t}. Only the "kp" method fills this: it integrates
     * the covariance of the Ko-Pender diffusion limit alongside the fluid mean, and every
     * other fluid method carries a first moment only, leaving this null.
     */
    public Matrix[][] QVart;

    /**
     * Full state covariance Sigma(t) of the "kp" method, one dim-by-dim matrix per time
     * point, so cross-station and cross-class covariances survive rather than only the
     * per-block totals in {@link #QVart}. Null for every other method.
     */
    public Matrix[] Sigmat;

    /** Cache hit probabilities [nodes x classes] from cacheqn analysis */
    public Matrix hitProb;

    /** Cache miss probabilities [nodes x classes] from cacheqn analysis */
    public Matrix missProb;

    /**
     * State-level stationary covariance of the moment-closure methods, i.e. the
     * solution of the linear noise Lyapunov equation at the fixed point. Null for
     * every first-order method, which computes no second moment at all. Mirrors
     * the MATLAB {@code SolverFLD.getMoments} report.
     */
    public Matrix momentSigma;

    /** Per-(station,class) queue-length variance of the moment-closure methods. */
    public Matrix momentQVar;

    /** Per-station population variance of the moment-closure methods. */
    public Matrix momentSigma2;

    /**
     * The same variance as it entered the DRIFT, i.e. zero at the delay stations.
     * The passage-time ODE closes its capacity term at this one so that the
     * response-time distribution is measured on the drift that produced the mean.
     */
    public Matrix momentSigma2Drift;

    /** Outer (mean, covariance) iterations performed by the moment-closure methods. */
    public int momentOuterIters;

    /**
     * State coordinates of each (station,class) in {@link #momentSigma}, flattened
     * as {@code i*nclasses+k}. Null for every first-order method. Needed to read a
     * per-class population off the covariance, whose coordinates are service
     * phases rather than classes.
     */
    public int[][] momentClassBlock;

    /**
     * State coordinates of each station in {@link #momentSigma}: every class and
     * service phase it serves. Null for every first-order method.
     */
    public int[][] momentStationBlock;

    /**
     * Per-cache stationary covariance of the item occupancy under the linear
     * noise approximation, one entry per cache that has a drift-based fluid
     * model. Null unless the moment closure ran on a cache model.
     */
    public Matrix[] momentCacheSigma;

    /** Per-cache per-item miss probability, aligned with {@link #momentCacheSigma}. */
    public Matrix[] momentCachePi0;

    /** Per-cache, per-class variance of the miss probability that class sees. */
    public Matrix momentCacheMissProbVar;

    /** Node index of each entry in {@link #momentCacheSigma}. */
    public int[] momentCacheNode;

    // ---- the 'dae' route ------------------------------------------------
    // What the SIMULTANEOUS solve did, which the substitution route has no
    // counterpart for: it stops at outer_max with no indication either way,
    // which is the failure mode the DAE form replaces.

    /** Infinity norm of the residual the simultaneous closure solve stopped at. */
    /**
     * What the Petri route computes that the station table has no column for: the marking and its
     * variance, the per-mode firing flows, the invariants and which capacities bound. Null on
     * every model that is not a stochastic Petri net.
     */
    public jline.solvers.fluid.petri.PetriSolver.PetriReport petri;
    public double daeResidual;

    /** Whether that residual reached {@code options.tol}. */
    public boolean daeConverged;

    /**
     * Largest violation of the population-conservation rows at the reported
     * point. On this route conservation is an EQUATION, so this is a solver
     * tolerance rather than the accumulated integrator drift every other fluid
     * method leaves behind.
     */
    public double daeConservation;

    /** Human-readable origin of each finite capacity constraint. */
    public String[] daeCapacityLabel;

    /** Right-hand side of each capacity constraint. */
    public Matrix daeCapacityB;

    /** Value each constrained quantity actually took at the fixed point. */
    public Matrix daeCapacityValue;

    /** Which capacity constraints bound at the reported point. */
    public int[] daeCapacityActive;

    /**
     * Mass held OUTSIDE a capped region, per staging coordinate.
     * REPORTED HERE AND NOT FOLDED INTO QN, because a blocked job is at no
     * station -- the same choice LDES makes, whose station queues likewise sum
     * to less than N.
     */
    public Matrix daeStaging;

    /** Region each staging coordinate belongs to. */
    public int[] daeStagingRegion;

    /** Class each staging coordinate carries. */
    public int[] daeStagingClass;

    /** Total blocked mass; station queues sum to N minus this. */
    public double daeBlocked;

    /** Rate each waiting queue drains at, one per binding cap. */
    public Matrix daeDrain;

    /** True where a cap stages the blocked job rather than holding or losing it. */
    public boolean[] daeCapacityStaged;

    /** Which region each cap came from, -1 for a station buffer. */
    public int[] daeCapacityRegion;

    /** Which station each cap came from, -1 for a region cap. */
    public int[] daeCapacityStation;

    /**
     * Every time the trajectory made a cap start or stop binding, as
     * {t, row, kind} with kind 0 = release, 1 = activate, 2 = a crossing the cap
     * could not hold. Empty for a steady-state solve.
     */
    public java.util.List<double[]> daeCapacitySwitches;
}
