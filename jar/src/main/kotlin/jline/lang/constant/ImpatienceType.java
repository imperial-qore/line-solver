/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Constants for specifying customer impatience types.
 *
 * Impatience refers to customers leaving a queueing system without receiving service.
 * The main types are:
 * - RENEGING: Customer abandons after joining the queue (timer-based abandonment)
 * - BALKING: Customer refuses to join based on queue state (queue-length or wait-time based)
 * - RETRIAL: Customer moves to orbit and retries entry after a delay
 *
 * ID values match MATLAB ImpatienceType constants.
 */
public enum ImpatienceType implements Serializable {
    /**
     * Reneging - customer abandons the queue after waiting too long.
     * This is the standard "impatience" where a job joins the queue but leaves
     * before service if its patience time expires.
     * Currently supported by JMT solver.
     */
    RENEGING(1),

    /**
     * Balking - customer refuses to join the queue based on queue state.
     * This occurs when a customer arrives but decides not to join based on
     * factors like queue length or expected waiting time.
     */
    BALKING(2),

    /**
     * Retrial - customer moves to orbit and retries queue entry after a delay.
     * This occurs when a customer is rejected due to capacity constraints
     * and enters an "orbit" buffer, periodically attempting to re-enter.
     * Supports both unlimited retries and configurable maximum attempts.
     */
    RETRIAL(3);

    private final int id;

    ImpatienceType(int id) {
        this.id = id;
    }

    /**
     * Get the numeric ID for this impatience type.
     * IDs match MATLAB ImpatienceType constants.
     *
     * @return numeric ID
     */
    public int getID() {
        return id;
    }

    /**
     * Convert numeric ID to ImpatienceType enum value.
     *
     * @param id numeric impatience type ID
     * @return corresponding ImpatienceType enum value
     * @throws IllegalArgumentException if ID is not recognized
     */
    public static ImpatienceType fromID(int id) {
        for (ImpatienceType type : values()) {
            if (type.id == id) {
                return type;
            }
        }
        throw new IllegalArgumentException("Unrecognized impatience type ID: " + id);
    }

    /**
     * Convert ImpatienceType to human-readable text.
     *
     * @param type the impatience type
     * @return text description
     */
    public static String toText(ImpatienceType type) {
        switch (type) {
            case RENEGING:
                return "reneging";
            case BALKING:
                return "balking";
            case RETRIAL:
                return "retrial";
            default:
                return type.name().toLowerCase();
        }
    }
}
