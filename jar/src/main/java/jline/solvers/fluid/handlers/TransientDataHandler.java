/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.handlers;

import jline.util.matrix.Matrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.ops.DConvertMatrixStruct;
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
 * <p>Uses dense DMatrixRMaj internally for O(1) element access during step handling,
 * converting to Matrix objects only when the integration completes.</p>
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

    /** Dense time storage for O(1) writes during integration */
    private DMatrixRMaj denseT;

    /** Dense state storage for O(1) writes during integration */
    private DMatrixRMaj denseX;

    /** Maximum number of steps before array resizing */
    private int maxSteps;

    /** Current step count */
    private int stepCount;

    /** Number of state dimensions */
    private int numDimensions;

    /**
     * Constructs a transient data handler for the specified number of state dimensions.
     *
     * @param numDimensions number of state variables in the ODE system
     */
    public TransientDataHandler(int numDimensions) {
        this.numDimensions = numDimensions;
        maxSteps = 1000;
        denseT = new DMatrixRMaj(maxSteps, 1);
        denseX = new DMatrixRMaj(maxSteps, numDimensions);
        stepCount = 0;
    }

    /**
     * Handles each integration step by storing the time and state vector data.
     *
     * @param interpolator provides access to current time and interpolated state
     * @param isLast indicates if this is the final integration step
     */
    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast) {

        double currentTime = interpolator.getCurrentTime();
        denseT.set(stepCount, 0, currentTime);

        double[] x = interpolator.getInterpolatedState();
        for (int i = 0; i < x.length; i++) {
            denseX.set(stepCount, i, FastMath.max(0, x[i]));
        }
        stepCount++;

        if (isLast) {
            // Trim dense arrays to actual step count and convert to sparse Matrix
            // for downstream compatibility (many Matrix methods assume sparse CSC)
            DMatrixRMaj trimmedT = new DMatrixRMaj(stepCount, 1);
            DMatrixRMaj trimmedX = new DMatrixRMaj(stepCount, numDimensions);
            for (int i = 0; i < stepCount; i++) {
                trimmedT.set(i, 0, denseT.get(i, 0));
                for (int j = 0; j < numDimensions; j++) {
                    trimmedX.set(i, j, denseX.get(i, j));
                }
            }
            // O(n) conversion from dense to sparse CSC via DConvertMatrixStruct
            tVec = new Matrix(DConvertMatrixStruct.convert(trimmedT, (DMatrixSparseCSC) null, 0.0));
            xVec = new Matrix(DConvertMatrixStruct.convert(trimmedX, (DMatrixSparseCSC) null, 0.0));
            return;
        }

        // Dynamic resizing of dense arrays
        if (stepCount == maxSteps) {
            maxSteps *= 2;
            DMatrixRMaj newT = new DMatrixRMaj(maxSteps, 1);
            DMatrixRMaj newX = new DMatrixRMaj(maxSteps, numDimensions);
            for (int i = 0; i < stepCount; i++) {
                newT.set(i, 0, denseT.get(i, 0));
                for (int j = 0; j < numDimensions; j++) {
                    newX.set(i, j, denseX.get(i, j));
                }
            }
            denseT = newT;
            denseX = newX;
        }
    }

    @Override
    public void init(double t0, double[] x0, double t) {

        denseT.set(stepCount, 0, t0);
        for (int i = 0; i < x0.length; i++) {
            denseX.set(stepCount, i, x0[i]);
        }
        stepCount++;
    }
}
