/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Enumeration of polling service types for polling scheduling strategy.
 * <p>
 * These constants specify how a polling server handles jobs at each queue
 * during a polling cycle:
 * <ul>
 * <li>GATED - Serves only jobs present at the start of the polling cycle</li>
 * <li>EXHAUSTIVE - Serves all jobs until the queue is empty</li>
 * <li>KLIMITED - Serves at most K jobs per polling cycle</li>
 * </ul>
 */
public enum PollingType {
    GATED(0),
    EXHAUSTIVE(1),
    KLIMITED(2);

    private final int value;

    PollingType(int value) {
        this.value = value;
    }

    /**
     * Gets the numeric value of this polling type.
     *
     * @return the numeric value
     */
    public int getValue() {
        return value;
    }

    /**
     * Converts a numeric value to a PollingType enum value.
     *
     * @param value the numeric value
     * @return the corresponding PollingType enum value
     * @throws RuntimeException if the value is not recognized
     */
    public static PollingType fromValue(int value) {
        switch (value) {
            case 0:
                return GATED;
            case 1:
                return EXHAUSTIVE;
            case 2:
                return KLIMITED;
            default:
                throw new RuntimeException("Unable to return a PollingType for value: " + value);
        }
    }

    /**
     * Converts a PollingType enum value to its text representation.
     *
     * @param type the polling type to convert
     * @return the text representation of the type
     */
    public static String toText(PollingType type) {
        switch (type) {
            case GATED:
                return "Gated";
            case EXHAUSTIVE:
                return "Exhaustive";
            case KLIMITED:
                return "K-Limited";
            default:
                return "";
        }
    }

    /**
     * Gets the ID (numeric value) of the polling type.
     *
     * @param type the polling type
     * @return the ID value
     */
    public static int toId(PollingType type) {
        return type.getValue();
    }
}