/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.handlers;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;
import org.apache.commons.math3.util.FastMath;

/**
 * Step handler for capturing transient data during fluid solver ODE integration.
 * 
 * <p>This class implements the Apache Commons Math StepHandler interface to collect
 * time series data during numerical integration of ordinary differential equations
 * in fluid queueing network analysis. It captures state vectors at each integration
 * step for transient performance analysis.</p>
 * 
 * <p>Key capabilities:
 * <ul>
 *   <li>Dynamic storage of time points and state vectors</li>
 *   <li>Automatic array resizing for long simulations</li>
 *   <li>Non-negative value enforcement for physical constraints</li>
 *   <li>Memory optimization through array shrinking</li>
 * </ul>
 * </p>
 * 
 * @see jline.solvers.fluid.SolverFluid
 * @see StepHandler
 * @since 1.0
 */
public class TransientDataHandler implements StepHandler {

    /** Time vector storing integration time points */
    public Matrix tVec;
    
    /** State vector matrix storing system state at each time point [time x dimensions] */
    public Matrix xVec;
    
    /** Maximum number of steps before array resizing */
    private int maxSteps;
    
    /** Current step count */
    private int stepCount;

    /**
     * Constructs a transient data handler for the specified number of state dimensions.
     * 
     * @param numDimensions number of state variables in the ODE system
     */
    public TransientDataHandler(int numDimensions) {

        maxSteps = 1000; // Increased for better CDF resolution
        tVec = new Matrix(maxSteps, 1);
        xVec = new Matrix(maxSteps, numDimensions);
        stepCount = 0;
    }

    /**
     * Handles each integration step by storing the time and state vector data.
     * 
     * <p>This method is called by the ODE integrator at each step. It captures
     * the current time and interpolated state vector, enforcing non-negative
     * constraints and managing dynamic array resizing as needed.</p>
     * 
     * @param interpolator provides access to current time and interpolated state
     * @param isLast indicates if this is the final integration step
     */
    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast) {

        double currentTime = interpolator.getCurrentTime();
        tVec.set(stepCount, 0, currentTime);
        

        double[] x = interpolator.getInterpolatedState();
        for (int i = 0; i < x.length; i++) {
            xVec.set(stepCount, i, FastMath.max(0, x[i])); // Equivalent of NonNegative odeset in LINE
        }
        stepCount++;

        if (isLast) {
            // Shrinking arrays to minimum number of rows
            tVec.shrinkNumRows(stepCount);
            xVec.shrinkNumRows(stepCount);
            return;
        }

        // Dynamic resizing of JLineMatrix objects per Dynamic Array logic
        if (stepCount == maxSteps) {
            maxSteps *= 2;
            tVec.expandMatrix(maxSteps, tVec.getNumCols(), maxSteps * tVec.getNumCols());
            xVec.expandMatrix(maxSteps, xVec.getNumCols(), maxSteps * xVec.getNumCols());
        }
    }

    @Override
    public void init(double t0, double[] x0, double t) {

        tVec.set(stepCount, 0, t0);
        for (int i = 0; i < x0.length; i++) {
            xVec.set(stepCount, i, x0[i]);
        }
        stepCount++;
    }
}
