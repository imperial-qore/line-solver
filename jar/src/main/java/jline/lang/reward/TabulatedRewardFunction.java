/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.reward;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * A reward function backed by a lookup table over aggregated state vectors.
 *
 * <p>This class is used by language wrappers (e.g., the MATLAB JLINE bridge) to
 * marshal opaque reward function handles: the handle values are pre-computed in
 * the source language over the enumerable domain of aggregated states and
 * stored in a lookup table, which the CTMC reward analyzer then queries per
 * row of stateSpaceAggr.</p>
 *
 * <p>The state is encoded as a string key of the form "n11,n12,...,nMK" where
 * n_{i,k} is the number of class-k jobs at station i (the same encoding as
 * {@link jline.util.PrecomputedCDFunction}).</p>
 */
public class TabulatedRewardFunction implements RewardFunction, Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * Map from aggregated state key (e.g., "1,2,3") to reward value.
     */
    private final Map<String, Double> valueMap;

    public TabulatedRewardFunction() {
        this.valueMap = new HashMap<String, Double>();
    }

    /**
     * Adds a pre-computed reward value for a specific aggregated state.
     *
     * @param state the aggregated state vector (1 x M*K, station-major)
     * @param value the reward value for this state
     */
    public void addValue(Matrix state, double value) {
        valueMap.put(stateToKey(state), value);
    }

    /**
     * Gets the number of tabulated states.
     *
     * @return the number of pre-computed reward values
     */
    public int size() {
        return valueMap.size();
    }

    @Override
    public double compute(Matrix state, NetworkStruct sn) {
        String key = stateToKey(state);
        Double value = valueMap.get(key);
        if (value == null) {
            throw new RuntimeException("TabulatedRewardFunction: no tabulated reward value for aggregated state ["
                    + key + "]; the pre-computed domain does not cover the CTMC state space "
                    + "(set finite class capacities or use the native solver of the source language).");
        }
        return value;
    }

    /**
     * Converts a Matrix state to a string key (comma-separated integer counts).
     */
    private static String stateToKey(Matrix state) {
        StringBuilder sb = new StringBuilder();
        int len = Math.max(state.getNumRows(), state.getNumCols());
        for (int i = 0; i < len; i++) {
            if (i > 0) {
                sb.append(",");
            }
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

    @Override
    public String toString() {
        return "TabulatedRewardFunction[entries=" + valueMap.size() + "]";
    }
}
