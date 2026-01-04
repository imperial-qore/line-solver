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
 * A pre-computed class dependence function that stores function values for all possible
 * state combinations.
 *
 * <p>This class is used to convert MATLAB function handles to Java by pre-computing
 * all possible function values in MATLAB and storing them in a lookup table. When
 * the function is called in Java, it looks up the pre-computed value.</p>
 *
 * <p>The state is represented as a string key of the form "n1,n2,...,nR" where n_r
 * is the number of jobs of class r at the station.</p>
 */
public class PrecomputedCDFunction implements SerializableFunction<Matrix, Double>, Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * Map from state key (e.g., "1,2,3") to scaling factor.
     */
    private final Map<String, Double> valueMap;

    /**
     * Default value to return when the state is not found.
     */
    private final double defaultValue;

    /**
     * Number of classes in the model.
     */
    private final int numClasses;

    /**
     * Creates a new pre-computed class dependence function.
     *
     * @param numClasses the number of classes in the model
     * @param defaultValue the default value to return when a state is not found
     */
    public PrecomputedCDFunction(int numClasses, double defaultValue) {
        this.numClasses = numClasses;
        this.defaultValue = defaultValue;
        this.valueMap = new HashMap<String, Double>();
    }

    /**
     * Creates a new pre-computed class dependence function with default value of 1.0.
     *
     * @param numClasses the number of classes in the model
     */
    public PrecomputedCDFunction(int numClasses) {
        this(numClasses, 1.0);
    }

    /**
     * Adds a pre-computed value for a specific state.
     *
     * @param state the state vector as a Matrix (1 x numClasses)
     * @param value the function value for this state
     */
    public void addValue(Matrix state, double value) {
        String key = stateToKey(state);
        valueMap.put(key, value);
    }

    /**
     * Adds a pre-computed value for a specific state given as an array.
     *
     * @param state the state as an array of integers
     * @param value the function value for this state
     */
    public void addValue(int[] state, double value) {
        String key = arrayToKey(state);
        valueMap.put(key, value);
    }

    /**
     * Applies the function to the given state.
     *
     * @param ni the state vector as a Matrix (1 x numClasses or numClasses x 1)
     * @return the pre-computed function value, or the default value if not found
     */
    @Override
    public Double apply(Matrix ni) {
        String key = stateToKey(ni);
        return valueMap.getOrDefault(key, defaultValue);
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
        return valueMap.size();
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
        return "PrecomputedCDFunction[numClasses=" + numClasses + ", entries=" + valueMap.size() + "]";
    }
}
