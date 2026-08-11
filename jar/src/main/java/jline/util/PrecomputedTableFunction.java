/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * A pre-computed lookup table over per-class state vectors.
 *
 * <p>This is the shared storage used to convert MATLAB function handles to Java:
 * the handle is evaluated in MATLAB over every reachable state and the resulting
 * values are stored here, so that a call in Java is a table lookup. The state is
 * keyed by the string "n1,n2,...,nR", where n_r is the entry of the state vector
 * for class r.</p>
 *
 * <p>Subclasses bind the table to a particular functional contract, since the two
 * consumers disagree on the return type: {@link PrecomputedCDFunction} yields the
 * dimensionless class-dependence beta as a Matrix, whereas
 * {@link PrecomputedRateFunction} yields a scalar service rate mu(c) as a Double.</p>
 */
public abstract class PrecomputedTableFunction implements Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * Map from state key (e.g., "1,2,3") to the tabulated value.
     */
    private final Map<String, Double> valueMap;

    /**
     * Map from state key to a tabulated per-class value, for handles that return
     * one entry per class rather than a single shared scalar (see
     * {@link #addValue(Matrix, double[])}). Kept separate from valueMap so that
     * the scalar contract of {@link PrecomputedRateFunction} is untouched.
     */
    private final Map<String, double[]> vectorMap;

    /**
     * Value returned when the state is absent from the table.
     */
    private final double defaultValue;

    /**
     * Number of classes in the model.
     */
    private final int numClasses;

    /**
     * Creates a new pre-computed table.
     *
     * @param numClasses   the number of classes in the model
     * @param defaultValue the value to return when a state is not found
     */
    protected PrecomputedTableFunction(int numClasses, double defaultValue) {
        this.numClasses = numClasses;
        this.defaultValue = defaultValue;
        this.valueMap = new HashMap<String, Double>();
        this.vectorMap = new HashMap<String, double[]>();
    }

    /**
     * Adds a pre-computed value for a specific state.
     *
     * @param state the state vector as a Matrix (1 x numClasses)
     * @param value the function value for this state
     */
    public void addValue(Matrix state, double value) {
        valueMap.put(stateToKey(state), value);
    }

    /**
     * Adds a pre-computed value for a specific state given as an array.
     *
     * @param state the state as an array of integers
     * @param value the function value for this state
     */
    public void addValue(int[] state, double value) {
        valueMap.put(arrayToKey(state), value);
    }

    /**
     * Adds a pre-computed per-class value for a specific state, for a handle that
     * returns one entry per class (as the flow-equivalent-server aggregation
     * does, beta_r(n) = X_r(n)|n|/n_r). Storing it as a scalar would collapse the
     * per-class resolution, so it is kept whole here and returned intact by
     * {@link PrecomputedCDFunction#apply(Matrix)}.
     *
     * @param state the state vector as a Matrix (1 x numClasses)
     * @param value the per-class function values for this state
     */
    public void addValue(Matrix state, double[] value) {
        vectorMap.put(stateToKey(state), value.clone());
    }

    /**
     * Looks up the tabulated value for a state.
     *
     * @param state the state vector (1 x numClasses or numClasses x 1)
     * @return the pre-computed value, or the default value if not found
     */
    protected double lookup(Matrix state) {
        return valueMap.getOrDefault(stateToKey(state), defaultValue);
    }

    /**
     * Looks up the tabulated per-class value for a state.
     *
     * @param state the state vector (1 x numClasses or numClasses x 1)
     * @return the per-class values, or null when the state was tabulated as a
     *         scalar (or not at all)
     */
    protected double[] lookupVector(Matrix state) {
        return vectorMap.get(stateToKey(state));
    }

    /**
     * Converts a Matrix state to a string key.
     */
    private String stateToKey(Matrix state) {
        StringBuilder sb = new StringBuilder();
        int len = Math.max(state.getNumRows(), state.getNumCols());
        for (int i = 0; i < len; i++) {
            if (i > 0) {
                sb.append(",");
            }
            // Handle both row and column vectors
            double val;
            if (state.getNumRows() == 1) {
                val = state.get(0, i);
            } else {
                val = state.get(i, 0);
            }
            sb.append((int) val);
        }
        return sb.toString();
    }

    /**
     * Converts an int array state to a string key.
     */
    private String arrayToKey(int[] state) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < state.length; i++) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(state[i]);
        }
        return sb.toString();
    }

    /**
     * Gets the number of pre-computed values stored.
     *
     * @return the number of pre-computed values
     */
    public int size() {
        return valueMap.size() + vectorMap.size();
    }

    /**
     * Gets the number of classes.
     *
     * @return the number of classes
     */
    public int getNumClasses() {
        return numClasses;
    }

    @Override
    public String toString() {
        return getClass().getSimpleName() + "[numClasses=" + numClasses + ", entries=" + size() + "]";
    }
}
