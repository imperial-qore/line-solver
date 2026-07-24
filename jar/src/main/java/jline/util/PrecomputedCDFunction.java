/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * A pre-computed class dependence function beta_i(n) that stores function values
 * for all possible state combinations.
 *
 * <p>This class is used to convert MATLAB class-dependence handles to Java by
 * pre-computing all possible function values in MATLAB and storing them in a
 * lookup table (see {@link PrecomputedTableFunction}). When the function is
 * called in Java, it looks up the pre-computed value.</p>
 *
 * <p>A chain-independent handle returns the dimensionless scaling as a scalar,
 * which MATLAB broadcasts over the classes by elementwise multiplication in
 * State.afterEventStation. Such a scalar is returned here as a 1x1 Matrix, which
 * the JAR state machinery ({@code AfterEventStation.cdScalar}) broadcasts over
 * the classes in the same way.</p>
 *
 * <p>A chain-specific handle instead returns the length-R row vector
 * [beta_1(n), ..., beta_R(n)] (Sauer's mu_{r,i}(n)), as the flow-equivalent-server
 * aggregation produces; see {@link jline.api.fes.FesBetaFunction}. Such a state is
 * tabulated whole via {@link PrecomputedTableFunction#addValue(Matrix, double[])}
 * and returned here as a 1xR Matrix. Collapsing it to a scalar would silently
 * drop the per-class resolution and yield the unscaled network.</p>
 */
public class PrecomputedCDFunction extends PrecomputedTableFunction
        implements SerializableFunction<Matrix, Matrix>, Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * Creates a new pre-computed class dependence function.
     *
     * @param numClasses   the number of classes in the model
     * @param defaultValue the default value to return when a state is not found
     */
    public PrecomputedCDFunction(int numClasses, double defaultValue) {
        super(numClasses, defaultValue);
    }

    /**
     * Creates a new pre-computed class dependence function with default value of 1.0,
     * i.e. the neutral scaling.
     *
     * @param numClasses the number of classes in the model
     */
    public PrecomputedCDFunction(int numClasses) {
        this(numClasses, 1.0);
    }

    /**
     * Applies the class dependence function to the given state.
     *
     * @param ni the state vector as a Matrix (1 x numClasses or numClasses x 1)
     * @return the pre-computed scaling: a 1xR Matrix when the state was tabulated
     *         per class, otherwise a 1x1 Matrix broadcast over the classes
     */
    @Override
    public Matrix apply(Matrix ni) {
        double[] perClass = lookupVector(ni);
        if (perClass != null) {
            Matrix v = new Matrix(1, perClass.length);
            for (int r = 0; r < perClass.length; r++) {
                v.set(0, r, perClass[r]);
            }
            return v;
        }
        Matrix v = new Matrix(1, 1);
        v.set(0, 0, lookup(ni));
        return v;
    }
}
