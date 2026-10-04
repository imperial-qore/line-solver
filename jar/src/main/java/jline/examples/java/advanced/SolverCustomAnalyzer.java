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
 * Everything around the solution algorithm of the custom-solver template.
 *
 * <p>Any activity prior to or after launching the solution algorithm belongs
 * here: the template times the call and scrubs the NaNs the algorithm may leave
 * behind, exactly as the MATLAB twin
 * {@code examples/advanced/example_custom_solver/solver_custom_analyzer.m}
 * does.</p>
 */
public class SolverCustomAnalyzer {

    /**
     * Wraps {@link SolverCustomAlgorithm#solver_custom} with the bookkeeping a
     * solver owes its caller.
     *
     * @param sn      the model structure
     * @param options the solver options
     * @return the result, with {@code runtime} set
     */
    public static SolverResult solver_custom_analyzer(NetworkStruct sn, SolverOptions options) {
        long tstart = System.nanoTime();

        System.out.println("Any activity prior or after launching the solution algorithm "
                + "needs to be implemented in SolverCustomAnalyzer.java.");
        SolverResult result = SolverCustomAlgorithm.solver_custom(sn, options);

        zeroNaNs(result.QN);
        zeroNaNs(result.UN);
        zeroNaNs(result.RN);
        zeroNaNs(result.TN);
        zeroNaNs(result.CN);
        zeroNaNs(result.XN);

        result.runtime = (System.nanoTime() - tstart) / 1e9;
        return result;
    }

    private static void zeroNaNs(Matrix m) {
        if (m == null) {
            return;
        }
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                if (Double.isNaN(m.get(i, j))) {
                    m.set(i, j, 0.0);
                }
            }
        }
    }
}
