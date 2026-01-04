/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.reward;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * Functional interface for defining reward functions on CTMC states.
 *
 * A reward function maps a state vector and network structure to a scalar reward value.
 * This interface is used for CTMC reward computation via value iteration.
 *
 * Example usage:
 * <pre>
 * // Queue length reward
 * RewardFunction qlenReward = (state, sn) -> state.get(0, 1);
 *
 * // Throughput reward
 * RewardFunction tputReward = (state, sn) -> {
 *     double n = state.get(0, 1);
 *     return n > 0 ? n * sn.rates.get(1, 0) / n : 0;
 * };
 * </pre>
 */
@FunctionalInterface
public interface RewardFunction extends Serializable {

    /**
     * Compute the reward value for a given state.
     *
     * @param state The state vector (row from stateSpaceAggr matrix)
     *              Format: [n_{1,1}, n_{1,2}, ..., n_{M,K}] where n_{i,k} is
     *              the number of class-k jobs at station i
     * @param sn The NetworkStruct providing access to rates, parameters, etc.
     * @return The reward value for this state
     */
    double compute(Matrix state, NetworkStruct sn);
}
