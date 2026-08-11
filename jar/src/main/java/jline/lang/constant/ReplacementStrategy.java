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
    SFIFO,   // strict FIFO: FIFO eviction, no reinsertion/promotion on hit (not segmented FIFO)
    LRU,
    HLRU,    // hierarchical/k-LRU: h lists, LRU discipline, promote i->i+1 on hit
    CLIMB,   // move-up-one-position on hit (transposition rule)
    QLRU;    // q-LRU: LRU discipline with probabilistic admission q on a miss

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
            case "HLRU":
                return HLRU;
            case "CLIMB":
                return CLIMB;
            case "QLRU":
                return QLRU;
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
            case HLRU:
                return "hlru";
            case CLIMB:
                return "climb";
            case QLRU:
                return "qlru";
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
            case HLRU:
                return "ReplacementStrategy_HLRU";
            case CLIMB:
                return "ReplacementStrategy_CLIMB";
            case QLRU:
                return "ReplacementStrategy_QLRU";
            default:
                throw new RuntimeException("Unrecognized replacement strategy");
        }
    }
}
