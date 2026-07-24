/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Retrial policies for an orbiting population.
 *
 * <p>The policy fixes how the aggregate rate at which the orbit attempts to
 * re-enter the station depends on the orbit size n:</p>
 *
 * <ul>
 *   <li>{@code LINEAR} - each orbiting customer retries independently at rate
 *       nu, so the aggregate retrial rate is n*nu. This is the classical
 *       retrial queue of Falin and Templeton.</li>
 *   <li>{@code CONSTANT} - the orbit as a whole retries at rate nu whenever it
 *       is non-empty, independently of n. This models a single retrial
 *       controller shared by the orbit rather than per-customer timers.</li>
 * </ul>
 *
 * <p>Declared as a constant class rather than an enum so that the numeric
 * values match the MATLAB RetrialPolicy constants (LINEAR = 1, CONSTANT = 2)
 * carried in sn.retrialPolicy across the codebases.</p>
 */
public final class RetrialPolicy {

    /** Each orbiting customer retries at its own rate: aggregate rate n*nu. */
    public static final int LINEAR = 1;

    /** The orbit retries as a whole at rate nu whenever it is non-empty. */
    public static final int CONSTANT = 2;

    private RetrialPolicy() {
    }

    /**
     * Text representation of a retrial policy.
     *
     * @param type the policy constant
     * @return the policy name
     */
    public static String toText(int type) {
        switch (type) {
            case LINEAR:
                return "linear";
            case CONSTANT:
                return "constant";
            default:
                throw new IllegalArgumentException("Unrecognized retrial policy type: " + type);
        }
    }

    /**
     * Parses a retrial policy from its text representation.
     *
     * @param text the policy name
     * @return the policy constant
     */
    public static int fromText(String text) {
        if (text == null) {
            throw new IllegalArgumentException("Unrecognized retrial policy: null");
        }
        String t = text.toLowerCase();
        if (t.equals("linear") || t.equals("per-customer") || t.equals("falin")) {
            return LINEAR;
        }
        if (t.equals("constant") || t.equals("fixed") || t.equals("controller")) {
            return CONSTANT;
        }
        throw new IllegalArgumentException("Unrecognized retrial policy: " + text);
    }
}
