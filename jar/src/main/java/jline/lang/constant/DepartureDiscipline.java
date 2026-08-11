/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Departure disciplines for the depository of a queueing place (QPN semantics).
 *
 * A queueing place serves tokens in its embedded queue and, on service completion,
 * moves them to a depository from which they become available to the output transitions.
 * The departure discipline governs the order in which depository tokens become available.
 *
 * - Normal: tokens become available immediately upon service completion (standard QPN).
 * - FIFO: tokens become available in their order of arrival to the depository (Spinner et al.).
 */
public enum DepartureDiscipline implements Serializable {
    Normal(0),
    FIFO(1);

    private final int id;

    DepartureDiscipline(int id) {
        this.id = id;
    }

    public int getID() {
        return id;
    }

    public static DepartureDiscipline fromID(int id) {
        for (DepartureDiscipline discipline : values()) {
            if (discipline.id == id) {
                return discipline;
            }
        }
        return Normal;
    }

    public static String toText(DepartureDiscipline discipline) {
        switch (discipline) {
            case Normal:
                return "normal";
            case FIFO:
                return "fifo";
            default:
                return discipline.name();
        }
    }
}
