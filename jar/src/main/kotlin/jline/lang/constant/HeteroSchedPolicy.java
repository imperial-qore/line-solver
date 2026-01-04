/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Enumeration of scheduling policies for heterogeneous multiserver queues.
 * <p>
 * These policies determine how jobs are assigned to servers when multiple
 * server types are available and a job's class is compatible with more than
 * one server type.
 * <p>
 * The policies are based on JMT's heterogeneous server scheduling:
 * <ul>
 * <li>ORDER - Assign to first available compatible server type (in definition order)</li>
 * <li>ALIS - Assign Longest Idle Server</li>
 * <li>ALFS - Assign Longest Free Server (with fairness sorting)</li>
 * <li>FAIRNESS - Fair distribution across compatible server types</li>
 * <li>FSF - Fastest Server First (based on expected service time)</li>
 * <li>RAIS - Random Available Idle Server</li>
 * </ul>
 */
public enum HeteroSchedPolicy {
    /**
     * Assign to first available compatible server type in definition order.
     * This is the default policy. Servers are checked in the order they were
     * added to the queue, and the first compatible server with availability is used.
     */
    ORDER,

    /**
     * Assign Longest Idle Server.
     * Servers cycle through in order, with fully busy servers moving to the back
     * of the list. Provides round-robin behavior among compatible servers.
     */
    ALIS,

    /**
     * Assign Longest Free Server with fairness sorting.
     * Server types are sorted by compatibility coverage (servers with exclusive
     * classes first), then all servers cycle together to ensure fairness.
     */
    ALFS,

    /**
     * Fair distribution across compatible server types.
     * Simple round-robin among compatible servers, with used servers
     * moving to the back of the list.
     */
    FAIRNESS,

    /**
     * Fastest Server First.
     * Always selects the compatible server type with the fastest expected
     * service time for the job's class. Optimizes for lowest latency.
     */
    FSF,

    /**
     * Random Available Idle Server.
     * Randomly selects among compatible server types that have available capacity.
     */
    RAIS;

    /**
     * Converts a string representation to a HeteroSchedPolicy enum value.
     *
     * @param text the string representation (case-insensitive)
     * @return the corresponding HeteroSchedPolicy enum value
     * @throws IllegalArgumentException if the string is not recognized
     */
    public static HeteroSchedPolicy fromText(String text) {
        if (text == null) {
            throw new IllegalArgumentException("Policy text cannot be null");
        }
        switch (text.toUpperCase()) {
            case "ORDER":
                return ORDER;
            case "ALIS":
                return ALIS;
            case "ALFS":
                return ALFS;
            case "FAIRNESS":
                return FAIRNESS;
            case "FSF":
                return FSF;
            case "RAIS":
                return RAIS;
            default:
                throw new IllegalArgumentException(
                    "Unknown heterogeneous scheduling policy: " + text +
                    ". Valid values are: ORDER, ALIS, ALFS, FAIRNESS, FSF, RAIS");
        }
    }

    /**
     * Converts this policy to its text representation.
     *
     * @return the text representation of this policy
     */
    public String toText() {
        return this.name();
    }

    /**
     * Converts a HeteroSchedPolicy to its text representation.
     *
     * @param policy the policy to convert
     * @return the text representation
     */
    public static String toText(HeteroSchedPolicy policy) {
        return policy.name();
    }
}
