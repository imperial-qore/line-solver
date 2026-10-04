/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * A pre-computed scalar rate function that stores function values for all
 * possible state combinations.
 *
 * <p>This class is used to convert the MATLAB mu(c) handle of a pass-and-swap
 * (PAS) or order-independent (OI) queue to Java, by pre-computing the rate over
 * every ordered class sequence in MATLAB and storing it in a lookup table (see
 * {@link PrecomputedTableFunction}). The state c is the ordered list of 0-based
 * class indices, keyed as "c1,c2,...,cp".</p>
 *
 * <p>Unlike the dimensionless class dependence of {@link PrecomputedCDFunction},
 * mu(c) is a total service RATE and is consumed as a scalar
 * ({@code QueueNodeParam.svcRateFun}), hence the Double return type.</p>
 */
public class PrecomputedRateFunction extends PrecomputedTableFunction
        implements SerializableFunction<Matrix, Double>, Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * Creates a new pre-computed rate function.
     *
     * @param numClasses   the number of classes in the model
     * @param defaultValue the default rate to return when a state is not found
     */
    public PrecomputedRateFunction(int numClasses, double defaultValue) {
        super(numClasses, defaultValue);
    }

    /**
     * Creates a new pre-computed rate function defaulting to a zero rate on
     * states absent from the table.
     *
     * @param numClasses the number of classes in the model
     */
    public PrecomputedRateFunction(int numClasses) {
        this(numClasses, 0.0);
    }

    /**
     * Applies the rate function to the given ordered class sequence.
     *
     * @param c the ordered list of 0-based class indices as a row Matrix
     * @return the pre-computed total service rate, or the default if not found
     */
    @Override
    public Double apply(Matrix c) {
        return lookup(c);
    }
}
