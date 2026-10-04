/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * The bare solution algorithm of the custom-solver template.
 *
 * <p>This is the file to fill in: it receives the {@link NetworkStruct} and the
 * options and must return the six average-metric matrices. The template returns
 * zeros and says so, exactly as the MATLAB twin
 * {@code examples/advanced/example_custom_solver/solver_custom.m} does.</p>
 */
public class SolverCustomAlgorithm {

    /**
     * The solution algorithm proper: QN, UN, RN, TN, CN and XN for the model.
     *
     * @param sn      the model structure
     * @param options the solver options
     * @return a result carrying the six average-metric matrices
     */
    public static SolverResult solver_custom(NetworkStruct sn, SolverOptions options) {
        int M = sn.nstations;                 // number of stations
        int K = sn.nclasses;                  // number of classes
        Matrix N = sn.njobs;                  // job populations
        Matrix rates = sn.rates;              // arrival and service rates
        java.util.Map<Integer, Matrix> V = sn.visits;   // visits

        SolverResult result = new SolverResult();
        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        result.CN = new Matrix(1, K);
        result.XN = new Matrix(1, K);

        System.out.println("The solution algorithm needs to be implemented in "
                + "SolverCustomAlgorithm.java: returning with no result.");
        return result;
    }
}
