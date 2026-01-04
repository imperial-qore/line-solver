/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Constants for specifying balking strategies (customer refusal to join queue).
 *
 * Balking occurs when a customer arrives at a queue but decides not to join.
 * The decision can be based on:
 * - QUEUE_LENGTH: Balking probability depends on current queue length
 * - EXPECTED_WAIT: Balking probability depends on expected waiting time
 * - COMBINED: Both conditions are evaluated (OR logic - balk if either triggers)
 *
 * ID values match MATLAB BalkingStrategy constants.
 */
public enum BalkingStrategy implements Serializable {
    /**
     * Balking based on queue length.
     * Uses thresholds mapping queue length ranges to balking probabilities.
     * Example: balk with 50% probability when queue length is between 5-10.
     */
    QUEUE_LENGTH(1),

    /**
     * Balking based on expected waiting time.
     * Customer estimates waiting time and balks if it exceeds their patience threshold.
     * Requires knowledge of service rate for estimation.
     */
    EXPECTED_WAIT(2),

    /**
     * Combined balking strategy.
     * Customer balks if EITHER the queue-length threshold OR expected-wait
     * condition triggers. Provides flexibility for complex balking behaviors.
     */
    COMBINED(3);

    private final int id;

    BalkingStrategy(int id) {
        this.id = id;
    }

    /**
     * Get the numeric ID for this balking strategy.
     * IDs match MATLAB BalkingStrategy constants.
     *
     * @return numeric ID
     */
    public int getID() {
        return id;
    }

    /**
     * Convert numeric ID to BalkingStrategy enum value.
     *
     * @param id numeric balking strategy ID
     * @return corresponding BalkingStrategy enum value
     * @throws IllegalArgumentException if ID is not recognized
     */
    public static BalkingStrategy fromID(int id) {
        for (BalkingStrategy strategy : values()) {
            if (strategy.id == id) {
                return strategy;
            }
        }
        throw new IllegalArgumentException("Unrecognized balking strategy ID: " + id);
    }

    /**
     * Convert BalkingStrategy to human-readable text.
     *
     * @param strategy the balking strategy
     * @return text description
     */
    public static String toText(BalkingStrategy strategy) {
        switch (strategy) {
            case QUEUE_LENGTH:
                return "queue length";
            case EXPECTED_WAIT:
                return "expected wait";
            case COMBINED:
                return "combined";
            default:
                return strategy.name().toLowerCase();
        }
    }
}
