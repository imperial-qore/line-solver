/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Enumeration of removal policies for negative signals in G-networks.
 *
 * <p>When a negative signal (or catastrophe) arrives at a queue and removes
 * positive customers, the removal policy determines which customers are
 * selected for removal.
 *
 * <p>Policies:
 * <ul>
 *   <li>{@link #RANDOM} - Uniform random selection from all jobs at the station</li>
 *   <li>{@link #FCFS} - Remove oldest jobs first (first-come-first-served order)</li>
 *   <li>{@link #LCFS} - Remove newest jobs first (last-come-first-served order)</li>
 * </ul>
 *
 * @see jline.lang.Signal
 * @see jline.lang.constant.SignalType#CATASTROPHE
 */
public enum RemovalPolicy implements Serializable {
    RANDOM(0),
    FCFS(1),
    LCFS(2);

    private final int id;

    RemovalPolicy(int id) {
        this.id = id;
    }

    public int getID() {
        return id;
    }

    public static RemovalPolicy fromID(int id) {
        for (RemovalPolicy policy : values()) {
            if (policy.id == id) {
                return policy;
            }
        }
        return null;
    }

    public static String toText(RemovalPolicy policy) {
        if (policy == null) {
            return null;
        }
        switch (policy) {
            case RANDOM:
                return "random";
            case FCFS:
                return "fcfs";
            case LCFS:
                return "lcfs";
            default:
                return policy.name();
        }
    }

    public static RemovalPolicy fromText(String text) {
        if (text == null || text.isEmpty()) {
            return null;
        }
        switch (text.toLowerCase()) {
            case "random":
                return RANDOM;
            case "fcfs":
                return FCFS;
            case "lcfs":
                return LCFS;
            default:
                return null;
        }
    }
}
