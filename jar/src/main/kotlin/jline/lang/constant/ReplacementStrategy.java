/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

/**
 * Constants for specifying a cache replacement strategy
 */
public enum ReplacementStrategy {
    RR,
    FIFO,
    SFIFO,
    LRU;
    //HLRU;

    public static ReplacementStrategy fromText(String text) {
        if (text == null || text.trim().isEmpty()) {
            throw new IllegalArgumentException("Replacement strategy text cannot be null or empty");
        }

        String normalized = text.trim().toUpperCase();
        switch (normalized) {
            case "RR":
                return RR;
            case "FIFO":
                return FIFO;
            case "SFIFO":
                return SFIFO;
            case "LRU":
                return LRU;
            default:
                throw new IllegalArgumentException("Unknown replacement strategy: " + text);
        }
    }

    public static String toText(ReplacementStrategy r) {
        switch (r) {
            case RR:
                return "rr";
            case FIFO:
                return "fifo";
            case SFIFO:
                return "strict-fifo";
            case LRU:
                return "lru";
//            case HLRU:
//                return "hlru";
            default:
                throw new RuntimeException("Unrecognized replacement strategy");
        }
    }

    public static String toFeature(ReplacementStrategy r) {
        switch (r) {
            case RR:
                return "ReplacementStrategy_RR";
            case FIFO:
                return "ReplacementStrategy_FIFO";
            case SFIFO:
                return "ReplacementStrategy_SFIFO";
            case LRU:
                return "ReplacementStrategy_LRU";
//            case HLRU:
//                return "ReplacementStrategy_HLRU";
            default:
                throw new RuntimeException("Unrecognized replacement strategy");
        }
    }
}
