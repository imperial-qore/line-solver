/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid;

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
}
