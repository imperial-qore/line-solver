/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag;

import java.util.List;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * Result of the agent-based (RCAT) solver.
 *
 * <p>Beside the mean measures every solver reports, this carries the three
 * objects the decomposition is actually about: the converged reversed rates,
 * one scalar per synchronizing action, and the per-agent generators and
 * stationary vectors they induce. They are the whole coupling between agents,
 * which is why a distributed sweep exchanges only {@link #actionRates}.</p>
 */
public class AGResult extends SolverResult {

    /** Sweeps taken by the reversed-rate fixed point. */
    public int iter;

    /** Converged reversed rates, one per synchronizing action. */
    public Matrix actionRates;

    /** Stationary vector of every agent at the fixed point. */
    public List<Matrix> equilibrium;

    /** Generator of every agent at the fixed point. */
    public List<Matrix> generators;
}
