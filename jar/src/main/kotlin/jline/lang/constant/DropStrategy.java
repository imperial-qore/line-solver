/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Constants for specifying drop strategies at stations when capacity is exceeded.
 * ID values match MATLAB DropStrategy constants.
 *
 * Strategies:
 * - WaitingQueue: Job waits in queue (default, infinite capacity)
 * - Drop: Job is rejected and leaves the system
 * - BlockingAfterService: Job blocks server after service until destination available
 * - BlockingBeforeService: Job blocks server before service until destination available
 * - ReServiceOnRejection: Job is re-serviced when rejected
 * - Retrial: Job moves to orbit and retries after a delay (unlimited attempts)
 * - RetrialWithLimit: Job moves to orbit with maximum retry attempts, dropped after limit
 */
public enum DropStrategy implements Serializable {
    WaitingQueue(-1),
    Drop(1),
    BlockingAfterService(2),
    BlockingBeforeService(3),
    ReServiceOnRejection(4),
    Retrial(5),
    RetrialWithLimit(6);

    private final int id;

    DropStrategy(int id) {
        this.id = id;
    }

    public int getID() {
        return id;
    }

    public static DropStrategy fromID(int id) {
        for (DropStrategy strategy : values()) {
            if (strategy.id == id) {
                return strategy;
            }
        }
        return Drop;  // Default
    }

    public static String toText(DropStrategy strategy) {
        switch (strategy) {
            case WaitingQueue:
                return "waiting queue";
            case Drop:
                return "drop";
            case BlockingAfterService:
                return "BAS blocking";
            case BlockingBeforeService:
                return "BBS blocking";
            case ReServiceOnRejection:
                return "re-service on rejection";
            case Retrial:
                return "retrial";
            case RetrialWithLimit:
                return "retrial with limit";
            default:
                return strategy.name();
        }
    }
}
