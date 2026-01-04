/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.handlers;

import jline.lang.Network;
import jline.solvers.fluid.SolverFluid;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.util.FastMath;

public class PStarOptimisationFunction implements MultivariateFunction {

    private final Matrix QNSSA;
    private final Network model;
    private final boolean stiff;
    private int iteration;

    public PStarOptimisationFunction(Matrix QNSSA, Network model, boolean stiff) {
        this.QNSSA = QNSSA;
        this.model = model;
        this.stiff = stiff;
        this.iteration = 0;
    }

    private double computeError(Matrix QNFluid) {
        // Compute Error (L2 Norm - Euclidean)
        double tmpErrorValue = 0;
        for (int i = 0; i < QNSSA.getNumRows(); i++) {
            for (int j = 0; j < QNSSA.getNumCols(); j++) {
                tmpErrorValue += FastMath.pow(QNFluid.get(i, j) - QNSSA.get(i, j), 2);
            }
        }
        return tmpErrorValue;
    }

    @Override
    public double value(double[] doubles) {

        System.out.format("Iteration %d\n", this.iteration);
        for (int i = 0; i < doubles.length; i++) {
            System.out.format("pStar for Station %d: %f\n", i, doubles[i]);
        }

        SolverFluid solverFluid = new SolverFluid(this.model);
        solverFluid.options.stiff = this.stiff;

        for (int i = 0; i < doubles.length; i++) {
            solverFluid.options.config.pstar.add(i, doubles[i]);
        }
        solverFluid.runAnalyzer();
        Matrix QNFluid = solverFluid.result.QN;
        double errorValue = computeError(QNFluid);

        System.out.format("Error Value: %f\n\n", errorValue);

        iteration++;
        return errorValue;
    }
}
